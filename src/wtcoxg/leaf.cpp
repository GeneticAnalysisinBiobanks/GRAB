// leaf.cpp — LEAF full implementation

#include "wtcoxg/leaf.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "wtcoxg/wtcoxg.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// ======================================================================
// loadAndMatchRefAf — load one plink2 .afreq and match vs PlinkData
// ======================================================================

std::vector<PopMatchedAF> loadAndMatchRefAf(
    const GenoMeta &plinkData,
    const std::string &afreqFile
) {

    bool isNumeric = false;
    auto records = loadRefAfFile(afreqFile, &isNumeric);

    if (isNumeric) {
        // Numeric fallback: rows in .bim order, AF is col 5 freq directly
        const auto &markers = plinkData.markerInfo();
        if (records.size() != markers.size())throw std::runtime_error(
                      "ref-af numeric fallback: row count (" + std::to_string(records.size()) +
                      ") != bim marker count (" + std::to_string(markers.size()) + ")"
        );
        std::vector<PopMatchedAF> matched;
        matched.reserve(markers.size());
        for (size_t i = 0; i < markers.size(); ++i) {
            PopMatchedAF pm;
            pm.genoIndex = markers[i].genoIndex;
            pm.af = records[i].alt_freq;
            pm.obs_ct = records[i].obs_ct;
            matched.push_back(pm);
        }
        return matched;
    }

    // Build lookup: "chrom:id" → index
    auto makeKey = [](const std::string &chr, const std::string &id) -> std::string {
        std::string k;
        k.reserve(chr.size() + 1 + id.size());
        k += chr;
        k += ':';
        k += id;
        return k;
    };

    std::unordered_map<std::string, size_t> refMap;
    refMap.reserve(records.size());
    for (size_t i = 0; i < records.size(); ++i)
        refMap.emplace(makeKey(records[i].chrom, records[i].id), i);

    std::vector<PopMatchedAF> matched;
    matched.reserve(plinkData.markerInfo().size());
    for (const auto &mi : plinkData.markerInfo()) {
        auto key = makeKey(mi.chrom, mi.id);
        auto it = refMap.find(key);
        if (it == refMap.end()) continue;

        const auto &rec = records[it->second];
        PopMatchedAF pm;
        pm.genoIndex = mi.genoIndex;

        if (rec.alt_allele == mi.ref && rec.ref_allele == mi.alt) {
            pm.af = rec.alt_freq;
        } else if (rec.ref_allele == mi.ref && rec.alt_allele == mi.alt) {
            pm.af = 1.0 - rec.alt_freq;
        } else {
            continue;
        }
        pm.obs_ct = rec.obs_ct;
        matched.push_back(pm);
    }
    return matched;
}

// ======================================================================
// kmeansCluster — K-means with K-means++ init and nstart restarts
//
// Distance computation uses the identity
//   ||x_i − c_j||² = ||x_i||² − 2 x_i·c_j + ||c_j||²
// so the n×k distance matrix is obtained via a single BLAS dgemm
// (Eigen maps to BLAS when available), making the assignment step
// O(n·k·p) with excellent vectorisation and cache behaviour.
// ======================================================================

Eigen::VectorXi kmeansCluster(
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    int k,
    int nstart,
    uint64_t seed,
    int nthreads
) {

    const Eigen::Index n = X.rows();
    const Eigen::Index p = X.cols();

    if (k <= 0 || k > n) throw std::runtime_error("kmeansCluster: k must be in [1, n]");

    // Pre-compute ||x_i||² (read-only, shared across threads).
    Eigen::VectorXd xSqNorm = X.rowwise().squaredNorm(); // (n)

    const uint64_t baseSeed = (seed != 0) ? seed : static_cast<uint64_t>(std::random_device{}());
    const int maxIter = 300;

    // Per-thread state for parallel restarts.
    struct ThreadState {
        Eigen::VectorXi bestLabels;
        double bestInertia = std::numeric_limits<double>::infinity();
    };

    const int nWorkers = std::min(nthreads, nstart);
    std::vector<ThreadState> threadState(nWorkers);
    for (auto &ts : threadState)
        ts.bestLabels = Eigen::VectorXi::Zero(n);

    std::atomic<int> nextRestart{0};

    auto worker = [&](int workerIdx) {
        // Per-thread RNG (deterministic: baseSeed + restart index).
        // Per-thread scratch buffers.
        Eigen::MatrixXd centers(k, p);
        Eigen::MatrixXd XC(n, k);
        Eigen::VectorXd cSqNorm(k);
        Eigen::VectorXi labels(n);
        Eigen::VectorXd minDist(n);
        Eigen::VectorXd counts(k);

        auto &ts = threadState[workerIdx];

        for (;;) {
            int restart = nextRestart.fetch_add(1, std::memory_order_relaxed);
            if (restart >= nstart) break;

            std::mt19937_64 rng(baseSeed + static_cast<uint64_t>(restart));

            // ── K-means++ initialisation ──────────────────────────────
            {
                std::uniform_int_distribution<Eigen::Index> uniIdx(0, n - 1);
                centers.row(0) = X.row(uniIdx(rng));
                minDist = (X.rowwise() - centers.row(0)).rowwise().squaredNorm();

                for (int c = 1; c < k; ++c) {
                    std::discrete_distribution<Eigen::Index> dd(minDist.data(), minDist.data() + n);
                    Eigen::Index pick = dd(rng);
                    centers.row(c) = X.row(pick);
                    Eigen::VectorXd newDist = (X.rowwise() - centers.row(c)).rowwise().squaredNorm();
                    minDist = minDist.cwiseMin(newDist);
                }
            }

            // ── Lloyd iterations ──────────────────────────────────────
            double inertia = 0.0;

            for (int iter = 0; iter < maxIter; ++iter) {
                cSqNorm = centers.rowwise().squaredNorm();
                XC.noalias() = X * centers.transpose(); // BLAS dgemm

                inertia = 0.0;
                bool changed = false;
                for (Eigen::Index i = 0; i < n; ++i) {
                    double best = std::numeric_limits<double>::infinity();
                    int bestJ = 0;
                    const double xi2 = xSqNorm[i];
                    for (int j = 0; j < k; ++j) {
                        double d2 = xi2 - 2.0 * XC(i, j) + cSqNorm[j];
                        if (d2 < best) {
                            best = d2;
                            bestJ = j;
                        }
                    }
                    if (labels[i] != bestJ) {
                        labels[i] = bestJ;
                        changed = true;
                    }
                    inertia += best;
                }

                centers.setZero();
                counts.setZero();
                for (Eigen::Index i = 0; i < n; ++i) {
                    int c = labels[i];
                    centers.row(c) += X.row(i);
                    counts[c] += 1.0;
                }
                for (int j = 0; j < k; ++j) {
                    if (counts[j] > 0.0) centers.row(j) /= counts[j];
                }

                if (!changed) break;
            }

            if (inertia < ts.bestInertia) {
                ts.bestInertia = inertia;
                ts.bestLabels = labels;
            }
        }
    };

    // Launch workers.
    std::vector<std::thread> threads;
    threads.reserve(nWorkers - 1);
    for (int t = 1; t < nWorkers; ++t)
        threads.emplace_back(worker, t);
    worker(0); // main thread participates
    for (auto &th : threads) th.join();

    // Reduce: pick global best across threads.
    int bestIdx = 0;
    for (int t = 1; t < nWorkers; ++t)
        if (threadState[t].bestInertia < threadState[bestIdx].bestInertia)
            bestIdx = t;

    Eigen::VectorXi bestLabels = std::move(threadState[bestIdx].bestLabels);
    bestLabels.array() += 1; // 0-based → 1-based (R convention)
    return bestLabels;
}

// ======================================================================
// Summix — exact active-set enumeration with KKT / Lagrange multiplier
//
// For each non-empty subset S ⊆ {0..K-1} (K ≤ 8 → max 255 subsets),
// solve the equality-constrained LS:
//   min ||D_S * p_S - o||²   s.t.  1'p_S = 1
// via the (|S|+1) × (|S|+1) KKT system:
//   [2*D_S'D_S  1] [p_S  ]   [2*D_S'o]
//   [1'         0] [lambda] = [1      ]
// Keep the solution with all p_S ≥ 0 and lowest objective.
// ======================================================================

Eigen::VectorXd summixEstimate(
    const Eigen::VectorXd &observedAF,
    const Eigen::MatrixXd &refAF
) {

    const int K = static_cast<int>(refAF.cols());
    if (K == 0) return Eigen::VectorXd();
    if (K > 8) throw std::runtime_error("summix: nPop > 8 not supported (got " + std::to_string(K) + ")");

    // Filter out rows with NaN
    std::vector<Eigen::Index> validRows;
    validRows.reserve(observedAF.size());
    for (Eigen::Index i = 0; i < observedAF.size(); ++i) {
        if (std::isnan(observedAF[i])) continue;
        bool ok = true;
        for (int p = 0; p < K; ++p)
            if (std::isnan(refAF(i, p))) {
                ok = false;
                break;
            }
        if (ok) validRows.push_back(i);
    }

    if (validRows.empty()) return Eigen::VectorXd::Constant(K, 1.0 / K);

    const int nValid = static_cast<int>(validRows.size());
    Eigen::VectorXd obs(nValid);
    Eigen::MatrixXd D(nValid, K);
    for (int i = 0; i < nValid; ++i) {
        obs[i] = observedAF[validRows[i]];
        D.row(i) = refAF.row(validRows[i]);
    }

    // Precompute full D'D and D'o
    Eigen::MatrixXd DtD = D.transpose() * D;   // K×K
    Eigen::VectorXd Dto = D.transpose() * obs; // K
    double oTo = obs.squaredNorm();

    Eigen::VectorXd bestP = Eigen::VectorXd::Constant(K, 1.0 / K);
    double bestObj = std::numeric_limits<double>::infinity();

    // Enumerate all non-empty subsets of {0..K-1}
    const int nSubsets = (1 << K);
    for (int mask = 1; mask < nSubsets; ++mask) {
        // Count and collect active indices
        int s = 0;
        int idx[8];
        for (int j = 0; j < K; ++j)
            if (mask & (1 << j)) idx[s++] = j;

        // Build the (s+1) × (s+1) KKT system
        // [2*A   e] [p_S  ]   [2*b]
        // [e'   0] [lambda] = [1  ]
        // where A = D_S'D_S (s×s), b = D_S'o (s), e = ones(s)
        const int dim = s + 1;
        Eigen::MatrixXd KKT = Eigen::MatrixXd::Zero(dim, dim);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dim);

        for (int i = 0; i < s; ++i) {
            for (int j = 0; j < s; ++j)
                KKT(i, j) = 2.0 * DtD(idx[i], idx[j]);
            KKT(i, s) = 1.0; // e column
            KKT(s, i) = 1.0; // e' row
            rhs[i] = 2.0 * Dto[idx[i]];
        }
        rhs[s] = 1.0; // sum constraint

        // Solve
        Eigen::VectorXd sol = KKT.partialPivLu().solve(rhs);

        // Check all p_S >= 0
        bool feasible = true;
        for (int i = 0; i < s; ++i) {
            if (sol[i] < -1e-12) {
                feasible = false;
                break;
            }
        }
        if (!feasible) continue;

        // Clamp tiny negatives to 0
        for (int i = 0; i < s; ++i)
            if (sol[i] < 0.0) sol[i] = 0.0;

        // Compute objective: ||D_S p_S - o||²  = p_S' A p_S - 2 b' p_S + o'o
        Eigen::VectorXd pS = sol.head(s);
        double obj = 0.0;
        for (int i = 0; i < s; ++i)
            for (int j = 0; j < s; ++j)
                obj += pS[i] * DtD(idx[i], idx[j]) * pS[j];
        for (int i = 0; i < s; ++i)
            obj -= 2.0 * Dto[idx[i]] * pS[i];
        obj += oTo;

        if (obj < bestObj) {
            bestObj = obj;
            bestP.setZero();
            for (int i = 0; i < s; ++i)
                bestP[idx[i]] = pS[i];
        }
    }

    return bestP;
}

// ======================================================================
// LEAFMethod — MethodBase implementation
// ======================================================================

LEAFMethod::LEAFMethod(
    std::vector<std::unique_ptr<WtCoxGMethod> > clusterMethods,
    std::vector<std::vector<uint32_t> > clusterIndices
)
    : m_nCluster(static_cast<int>(clusterMethods.size())),
      m_clusterMethods(std::move(clusterMethods)),
      m_clusterIndices(std::move(clusterIndices))
{
    m_clusterGVec.resize(m_nCluster);
    for (int c = 0; c < m_nCluster; ++c)
        m_clusterGVec[c].resize(m_clusterIndices[c].size());
}

std::unique_ptr<MethodBase> LEAFMethod::clone() const {
    std::vector<std::unique_ptr<WtCoxGMethod> > cloned;
    cloned.reserve(m_nCluster);
    for (const auto &m : m_clusterMethods) {
        auto p = m->clone(); // returns unique_ptr<MethodBase>
        cloned.push_back(std::unique_ptr<WtCoxGMethod>(static_cast<WtCoxGMethod *>(p.release())));
    }
    return std::make_unique<LEAFMethod>(std::move(cloned), m_clusterIndices);
}

int LEAFMethod::resultSize() const {
    return 2 + 3 * m_nCluster;
}

std::string LEAFMethod::getHeaderColumns() const {
    std::ostringstream oss;
    oss << "\tmeta.p_ext\tmeta.p_noext";
    for (int i = 1; i <= m_nCluster; ++i)
        oss << "\tcl" << i << ".p_ext\tcl" << i << ".p_noext\tcl" << i << ".p_batch";
    return oss.str();
}

void LEAFMethod::prepareChunk(const std::vector<uint64_t> &gIndices) {
    for (auto &m : m_clusterMethods)
        m->prepareChunk(gIndices);
}

void LEAFMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double /*altFreq*/,
    int markerInChunkIdx,
    std::vector<double> &result
) {

    std::vector<double> pExt(m_nCluster), pNoext(m_nCluster);
    std::vector<double> sExt(m_nCluster), sNoext(m_nCluster);

    for (int c = 0; c < m_nCluster; ++c) {
        // Gather cluster genotypes from full GVec
        const auto &idx = m_clusterIndices[c];
        Eigen::VectorXd &gClu = m_clusterGVec[c];
        for (size_t k = 0; k < idx.size(); ++k)
            gClu[static_cast<Eigen::Index>(k)] = GVec[idx[k]];

        auto dr = m_clusterMethods[c]->computeDual(gClu, markerInChunkIdx);
        pExt[c] = dr.p_ext;
        pNoext[c] = dr.p_noext;
        sExt[c] = dr.score_ext;
        sNoext[c] = dr.score_noext;
    }

    // Fixed-effects meta-analysis: pool scores across clusters
    auto metaP = [](const std::vector<double> &scores, const std::vector<double> &pvals) -> double {
        double sumScore = 0.0, sumVar = 0.0;
        for (size_t c = 0; c < scores.size(); ++c) {
            if (std::isnan(scores[c]) || std::isnan(pvals[c]) || pvals[c] <= 0.0 ||
                pvals[c] >= 1.0) continue;
            double chisq = math::qchisq(pvals[c], 1.0, false, false);
            if (chisq < 1e-30) chisq = 1e-30;
            double var = (scores[c] * scores[c]) / chisq;
            if (std::isnan(var)) continue;
            sumScore += scores[c];
            sumVar += var;
        }
        if (sumVar <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        double z = sumScore / std::sqrt(sumVar);
        return std::erfc(std::fabs(z) / std::sqrt(2.0));              // two-sided p-value
    };

    result.push_back(metaP(sExt, pExt));
    result.push_back(metaP(sNoext, pNoext));
    for (int c = 0; c < m_nCluster; ++c) {
        result.push_back(pExt[c]);
        result.push_back(pNoext[c]);
        result.push_back(m_clusterMethods[c]->chunkRefInfoAt(markerInChunkIdx).pvalue_bat);
    }
}

// ======================================================================
// runLEAFPheno — --pheno path: K-means + per-cluster regression
//
// Follows the R reference code (LEAF.md):
//   1. Load phenotype + covariates
//   2. K-means on PCs → nClusters clusters
//   3. Per-cluster: compute calRegrWeight, fit glm/coxph, get residuals
//   4. Build per-cluster SubjectData with disjoint subject masks
//   5. Load PLINK with union mask, load ref-af files
//   6. Per-cluster: summix + batch-effect testing
//   7. Build LEAFMethod → markerEngine
// ======================================================================

#include "wtcoxg/regression.hpp"

void runLEAFPheno(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &pcColNames,
    int nClusters,
    uint64_t seed,
    const GenoSpec &geno,
    const std::vector<std::string> &refAfFiles,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile
) {

    const bool isSurv = (phenoNames.size() >= 2);
    const int nCluster = nClusters;

    // ---- 1. Load phenotype / covariate data ----
    infoMsg("Loading phenotype file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    const uint32_t nFamTotal = static_cast<uint32_t>(famIIDs.size());
    SubjectData sdFull(famIIDs); // copy — need famIIDs later
    sdFull.loadPhenoFile(phenoFile);
    if (!covarFile.empty()) {
        infoMsg("Loading covariate file: %s", covarFile.c_str());
        // Need both covarNames and pcColNames from the covar file
        std::vector<std::string> needed = covarNames;
        needed.insert(needed.end(), pcColNames.begin(), pcColNames.end());
        sdFull.loadCovar(covarFile, needed);
    }
    sdFull.setKeepRemove(keepFile, removeFile);
    if (!spgrmGrabFile.empty() || !spgrmGctaFile.empty())sdFull.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile,
                                                                                                          spgrmGctaFile,
                                                                                                          sdFull.famIIDs
                                                                                                              ()));
    sdFull.setGenoLabel(geno.flagLabel());
    sdFull.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sdFull.finalize();
    // Drop subjects with NA in the selected phenotype column(s)
    if (isSurv) {
        sdFull.dropNaInColumns({phenoNames[0], phenoNames[1]});
    } else {
        sdFull.dropNaInColumns({phenoNames[0]});
    }
    infoMsg("  %u subjects loaded", sdFull.nUsed());

    const Eigen::Index N = static_cast<Eigen::Index>(sdFull.nUsed());

    // ---- Extract response ----
    Eigen::VectorXd indicator;
    Eigen::VectorXd survTime;

    if (isSurv) {
        survTime = sdFull.getColumn(phenoNames[0]);
        indicator = sdFull.getColumn(phenoNames[1]);
        infoMsg("  Survival phenotype: time=%s, event=%s", phenoNames[0].c_str(), phenoNames[1].c_str());
    } else {
        indicator = sdFull.getColumn(phenoNames[0]);
        // Validate binary: all values must be 0 or 1
        for (Eigen::Index i = 0; i < indicator.size(); ++i) {
            double v = indicator[i];
            if (v != 0.0 && v != 1.0)
                throw std::runtime_error(
                          "Phenotype '" + phenoNames[0] + "' is not binary"
                          " (found value " + std::to_string(v) + ")."
                          " For survival phenotypes use --pheno-name TIME:EVENT.");
        }
        infoMsg("  Binary phenotype: %s", phenoNames[0].c_str());
    }

    // ---- Build design matrices ----
    // Logistic: [1 | covariates] (intercept needed)
    // Cox PH:   [covariates]     (no intercept — absorbed into baseline hazard)
    Eigen::MatrixXd covarMat;
    if (!covarNames.empty()) {
        covarMat = sdFull.getColumns(covarNames);
    } else if (sdFull.hasCovar()) {
        covarMat = sdFull.covar();
    }
    const int nCov = static_cast<int>(covarMat.cols());
    if (nCov > 0) infoMsg("  %d covariate(s)", nCov);
    Eigen::MatrixXd logisticDesign; // only built if needed
    if (!isSurv) {
        logisticDesign.resize(N, 1 + nCov);
        logisticDesign.col(0).setOnes();
        if (nCov > 0) logisticDesign.rightCols(nCov) = covarMat;
    }

    // ---- 2. K-means clustering on PCs ----
    Eigen::MatrixXd PCs = sdFull.getColumns(pcColNames);
    infoMsg("  %d PCs for K-means clustering", static_cast<int>(PCs.cols()));
    infoMsg("K-means clustering into %d clusters...", nCluster);
    Eigen::VectorXi clusterLabels = kmeansCluster(PCs, nCluster, /*nstart=*/ 25, seed, nthread);

    // ---- 3. Per-cluster regression → resid/weight/indicator ----
    //
    // Following R reference (LEAF.md):
    //   Indicator_i <- indicator[idx_i]
    //   weight_i    <- calRegrWeight(Indicator_i, RefPrevalence)
    //   residual_i  <- glm.fit(designMat[idx_i,], Indicator_i,
    //                          weights=weight_i, family=binomial(),
    //                          intercept=TRUE)$residuals
    //
    // Note: R residuals(glm(...), type="response") returns y − μ̂.
    // Our regression::logisticResiduals returns the same response residuals.

    // Map cluster labels → member indices (into full used-subject array)
    std::vector<std::vector<uint32_t> > clusterMemberIdx(nCluster);
    for (Eigen::Index i = 0; i < N; ++i)
        clusterMemberIdx[clusterLabels[i] - 1].push_back(static_cast<uint32_t>(i));

    // Per-cluster residuals/weights/indicators (indexed by cluster member order)
    struct ClusterRWI {
        Eigen::VectorXd resid, weight, ind;
    };

    std::vector<ClusterRWI> clusterRWI(nCluster);

    infoMsg("Per-cluster regression (%d clusters)...", nCluster);
    for (int c = 0; c < nCluster; ++c) {
        const auto &members = clusterMemberIdx[c];
        const Eigen::Index Nc = static_cast<Eigen::Index>(members.size());

        // Subset cluster data
        Eigen::VectorXd cInd(Nc), cTime(Nc);
        for (Eigen::Index j = 0; j < Nc; ++j) {
            auto i = members[j];
            cInd[j] = indicator[i];
            if (isSurv) cTime[j] = survTime[i];
        }

        // calRegrWeight per cluster (R: weight_i <- calRegrWeight(Indicator_i, ...))
        Eigen::VectorXd cWeight = regression::calRegrWeight(refPrevalence, cInd);

        // Subset design matrix and fit regression per cluster
        Eigen::VectorXd cResid;
        if (isSurv) {
            // Cox PH: no intercept
            Eigen::MatrixXd cDesign(Nc, covarMat.cols());
            for (Eigen::Index j = 0; j < Nc; ++j)
                cDesign.row(j) = covarMat.row(members[j]);
            cResid = regression::coxResiduals(cTime, cInd, cDesign, cWeight);
        } else {
            // Logistic: with intercept
            Eigen::MatrixXd cDesign(Nc, logisticDesign.cols());
            for (Eigen::Index j = 0; j < Nc; ++j)
                cDesign.row(j) = logisticDesign.row(members[j]);
            cResid = regression::logisticResiduals(cInd, cDesign, cWeight);
        }
        infoMsg("  Cluster %d: %d subjects, resid computed", c + 1, static_cast<int>(Nc));

        clusterRWI[c] = {std::move(cResid), std::move(cWeight), std::move(cInd)};
    }

    // ---- 4. Build per-cluster SubjectData with disjoint masks ----
    const size_t nMaskWords = (nFamTotal + 63) / 64;

    // Build per-cluster bitmasks from fullMask + clusterLabels
    std::vector<std::vector<uint64_t> > clusterMasks(nCluster);
    for (int c = 0; c < nCluster; ++c)
        clusterMasks[c].assign(nMaskWords, 0);

    {
        const auto &fullMask = sdFull.usedMask();
        uint32_t usedIdx = 0;
        for (uint32_t f = 0; f < nFamTotal; ++f) {
            if (fullMask[f / 64] & (1ULL << (f % 64))) {
                int cl = clusterLabels[usedIdx] - 1; // 0-based
                clusterMasks[cl][f / 64] |= 1ULL << (f % 64);
                ++usedIdx;
            }
        }
    }

    // Count per-cluster subjects from masks
    std::vector<uint32_t> clusterN(nCluster, 0);
    for (int c = 0; c < nCluster; ++c)
        for (size_t w = 0; w < nMaskWords; ++w)
            clusterN[c] += static_cast<uint32_t>(__builtin_popcountll(clusterMasks[c][w]));

    // Create per-cluster SubjectData with proper disjoint masks
    std::vector<SubjectData> clusterSD;
    clusterSD.reserve(nCluster);
    for (int c = 0; c < nCluster; ++c) {
        auto fc = famIIDs; // copy
        clusterSD.emplace_back(std::move(fc));
        clusterSD.back().initFromMask(
            clusterMasks[c],
            clusterN[c],
            std::move(clusterRWI[c].resid),
            std::move(clusterRWI[c].weight),
            std::move(clusterRWI[c].ind)
        );
    }

    // ---- Build union bitmask and cluster indices (same as runLEAF) ----
    std::vector<uint64_t> unionMask(nMaskWords, 0);
    for (int c = 0; c < nCluster; ++c) {
        const auto &cmask = clusterSD[c].usedMask();
        for (size_t w = 0; w < nMaskWords; ++w)
            unionMask[w] |= cmask[w];
    }
    uint32_t nUnion = 0;
    for (size_t w = 0; w < nMaskWords; ++w)
        nUnion += static_cast<uint32_t>(__builtin_popcountll(unionMask[w]));

    std::vector<std::vector<uint32_t> > clusterIndices(nCluster);
    for (int c = 0; c < nCluster; ++c)
        clusterIndices[c].reserve(clusterSD[c].nUsed());

    uint32_t denseIdx = 0;
    for (size_t w = 0; w < nMaskWords; ++w) {
        uint64_t bits = unionMask[w];
        while (bits) {
            int bit = __builtin_ctzll(bits);
            uint64_t bitVal = uint64_t(1) << bit;
            for (int c = 0; c < nCluster; ++c) {
                if (clusterSD[c].usedMask()[w] & bitVal) {
                    clusterIndices[c].push_back(denseIdx);
                    break;
                }
            }
            ++denseIdx;
            bits &= bits - 1;
        }
    }
    infoMsg("  Total subjects (union): %u", nUnion);

    // ---- 5. Load genotype data with union bitmask ----
    auto genoData = makeGenoData(geno, unionMask, nFamTotal, nUnion, nSnpPerChunk);
    infoMsg("  %u subjects matched, %u markers", genoData->nSubjUsed(), genoData->nMarkers());

    // ---- Load per-population ref-af files (plink2 .afreq) ----
    const int nPop = static_cast<int>(refAfFiles.size());
    infoMsg("Loading %d reference AF files...", nPop);
    std::vector<std::vector<PopMatchedAF> > popMatched(nPop);
    for (int p = 0; p < nPop; ++p) {
        infoMsg("  Pop %d: %s", p + 1, refAfFiles[p].c_str());
        popMatched[p] = loadAndMatchRefAf(*genoData, refAfFiles[p]);
        infoMsg("    %zu markers matched", popMatched[p].size());
    }

    // ---- 6. Build intersection + summix + batch-effect (same as runLEAF) ----

    // Build intersection of markers present in all populations
    struct MultiPopEntry {
        double af[8];
        double obs_ct[8];
    };

    std::unordered_map<uint64_t, MultiPopEntry> allPopMap;
    if (nPop > 0) {
        for (const auto &pm : popMatched[0]) {
            MultiPopEntry e{};
            e.af[0] = pm.af;
            e.obs_ct[0] = pm.obs_ct;
            allPopMap.emplace(pm.genoIndex, e);
        }
    }
    for (int p = 1; p < nPop; ++p) {
        std::unordered_map<uint64_t, MultiPopEntry> next;
        for (const auto &pm : popMatched[p]) {
            auto it = allPopMap.find(pm.genoIndex);
            if (it == allPopMap.end()) continue;
            auto entry = it->second;
            entry.af[p] = pm.af;
            entry.obs_ct[p] = pm.obs_ct;
            next.emplace(pm.genoIndex, entry);
        }
        allPopMap = std::move(next);
    }

    // Convert to sorted vector for deterministic order
    struct MatchedMultiPop {
        uint64_t genoIndex;
        double popAF[8];
        double popObsCt[8];
    };

    std::vector<MatchedMultiPop> matchedMulti;
    matchedMulti.reserve(allPopMap.size());
    for (auto &[gi, e] : allPopMap) {
        MatchedMultiPop mm{};
        mm.genoIndex = gi;
        for (int p = 0; p < nPop; ++p) {
            mm.popAF[p] = e.af[p];
            mm.popObsCt[p] = e.obs_ct[p];
        }
        matchedMulti.push_back(mm);
    }
    std::sort(
        matchedMulti.begin(),
        matchedMulti.end(),
        [](const auto &a, const auto &b) {
        return a.genoIndex < b.genoIndex;
    }
    );
    const size_t nMatched = matchedMulti.size();
    infoMsg("  %zu markers present in all %d populations", nMatched, nPop);

    // Per-cluster summix + AF synthesis (scan genotypes for per-cluster AF)
    infoMsg("Per-cluster ancestry estimation and AF synthesis...");

    auto cursor = genoData->makeCursor();

    const uint32_t nFull = genoData->nSubjUsed();
    Eigen::VectorXd fullGVec(nFull);

    if (!matchedMulti.empty()) cursor->beginSequentialBlock(matchedMulti.front().genoIndex);

    struct ClMarkerStat {
        double intAF;
        double mu0, mu1, n0, n1, mu_int;
    };

    std::vector<std::vector<ClMarkerStat> > clStats(nCluster, std::vector<ClMarkerStat>(nMatched));

    infoMsg("  Scanning %zu matched markers for per-cluster allele frequencies...", nMatched);
    for (size_t m = 0; m < nMatched; ++m) {
        cursor->getGenotypesSimple(matchedMulti[m].genoIndex, fullGVec);

        for (int c = 0; c < nCluster; ++c) {
            const auto &idx = clusterIndices[c];
            const auto &ind = clusterSD[c].indicator();
            double sum0 = 0, sum1 = 0, cnt0 = 0, cnt1 = 0;
            for (size_t k = 0; k < idx.size(); ++k) {
                double g = fullGVec[idx[k]];
                if (std::isnan(g)) continue;
                if (ind[static_cast<Eigen::Index>(k)] == 1.0) {
                    sum1 += g;
                    cnt1 += 1.0;
                } else {
                    sum0 += g;
                    cnt0 += 1.0;
                }
            }
            double total = cnt0 + cnt1;
            auto &st = clStats[c][m];
            st.intAF = (total > 0) ? (sum0 + sum1) / (2.0 * total) : std::numeric_limits<double>::quiet_NaN();
            st.mu0 = (cnt0 > 0) ? sum0 / (2.0 * cnt0) : 0.0;
            st.mu1 = (cnt1 > 0) ? sum1 / (2.0 * cnt1) : 0.0;
            st.n0 = cnt0;
            st.n1 = cnt1;
            st.mu_int = st.intAF;
        }
    }

    // Per cluster: summix → AF_ref, obs_ct → MatchedMarkerInfo
    std::vector<std::vector<MatchedMarkerInfo> > clMatchedInfo(nCluster);

    for (int c = 0; c < nCluster; ++c) {
        infoMsg("  Cluster %d: running summix...", c + 1);

        Eigen::VectorXd obsAF(nMatched);
        Eigen::MatrixXd refMat(nMatched, nPop);
        for (size_t m = 0; m < nMatched; ++m) {
            obsAF[m] = clStats[c][m].intAF;
            for (int p = 0; p < nPop; ++p)
                refMat(m, p) = matchedMulti[m].popAF[p];
        }

        Eigen::VectorXd proportions = summixEstimate(obsAF, refMat);
        for (int p = 0; p < nPop; ++p)
            infoMsg("    pop%d: %.4f", p + 1, proportions[p]);

        clMatchedInfo[c].resize(nMatched);
        for (size_t m = 0; m < nMatched; ++m) {
            auto &mi = clMatchedInfo[c][m];
            mi.genoIndex = matchedMulti[m].genoIndex;

            double af_ref = 0.0;
            for (int p = 0; p < nPop; ++p)
                af_ref += proportions[p] * matchedMulti[m].popAF[p];
            mi.AF_ref = af_ref;

            double denom = 0.0;
            for (int p = 0; p < nPop; ++p) {
                double oc = matchedMulti[m].popObsCt[p];
                if (oc > 0 && proportions[p] > 0) denom += (proportions[p] * proportions[p]) / oc;
            }
            mi.obs_ct = (denom > 0) ? 1.0 / denom : 0.0;

            mi.mu0 = clStats[c][m].mu0;
            mi.mu1 = clStats[c][m].mu1;
            mi.n0 = clStats[c][m].n0;
            mi.n1 = clStats[c][m].n1;
            mi.mu_int = clStats[c][m].mu_int;
        }
    }

    // Per-cluster batch-effect testing
    infoMsg("Per-cluster batch-effect testing...");

    std::vector<std::shared_ptr<std::unordered_map<uint64_t, WtCoxGRefInfo> > > clRefMaps(nCluster);

    for (int c = 0; c < nCluster; ++c) {
        infoMsg("  Cluster %d: batch-effect testing (%zu markers)...", c + 1, clMatchedInfo[c].size());

        std::unique_ptr<SparseGRM> grm;
        if (!spgrmGctaFile.empty() || !spgrmGrabFile.empty()) {
            grm = std::make_unique<SparseGRM>(
                SparseGRM::load(spgrmGrabFile, spgrmGctaFile, clusterSD[c].usedIIDs(), clusterSD[c].famIIDs())
            );
            infoMsg("    Sparse GRM: %u subjects, %zu non-zeros", grm->nSubjects(), grm->nnz());
        }

        clRefMaps[c] = testBatchEffects(
            clMatchedInfo[c],
            clusterSD[c].residuals(),
            clusterSD[c].weights(),
            clusterSD[c].indicator(),
            grm.get(),
            refPrevalence,
            cutoff
        );
        infoMsg("    %zu markers retained", clRefMaps[c]->size());
    }

    // ---- 7. Build LEAFMethod and run marker engine ----
    infoMsg("Building LEAF method (%d clusters)...", nCluster);

    std::vector<std::unique_ptr<WtCoxGMethod> > clMethods;
    clMethods.reserve(nCluster);
    for (int c = 0; c < nCluster; ++c) {
        clMethods.push_back(
            std::make_unique<WtCoxGMethod>(
                clusterSD[c].residuals(),
                clusterSD[c].weights(),
                cutoff,
                spaCutoff,
                clRefMaps[c]
            )
        );
    }

    auto method = std::make_unique<LEAFMethod>(std::move(clMethods), clusterIndices);

    // Build PhenoTask (single trait)
    std::string traitName = isSurv ? phenoNames[0] + "_" + phenoNames[1] : phenoNames[0];
    std::vector<PhenoTask> tasks(1);
    tasks[0].phenoName = traitName;
    tasks[0].method = std::move(method);
    // Identity mapping: all subjects in union = all subjects in phenotype (K=1)
    tasks[0].unionToLocal.resize(genoData->nSubjUsed());
    std::iota(tasks[0].unionToLocal.begin(), tasks[0].unionToLocal.end(), 0u);
    tasks[0].nUsed = genoData->nSubjUsed();

    infoMsg("Running LEAF marker engine (%d thread(s))...", nthread);
    multiPhenoEngine(
        *genoData,
        tasks,
        outPrefix,
        "LEAF",
        compression,
        compressionLevel,
        nthread,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        hweCutoff
    );
}

// ======================================================================
// runLEAF — multi-phenotype entry point
//
// Phases:
//   A  Shared data loading, K-means, masks, genoData, ref-AF, GRM
//   B  Parallel null-model regression  (min(T, P) threads, K fits each)
//   C  Shared genotype scan, summix, AF synthesis
//   D  Parallel batch-effect testing   (min(T, P×K) threads)
//   E  Single multiPhenoEngine call    (T-thread chunk parallel)
// ======================================================================

static std::vector<std::string> parsePhenoSpecLEAF(const std::string &spec) {
    std::vector<std::string> cols;
    auto colon = spec.find(':');
    if (colon != std::string::npos) {
        cols.push_back(spec.substr(0, colon));
        cols.push_back(spec.substr(colon + 1));
    } else {
        cols.push_back(spec);
    }
    return cols;
}

void runLEAF(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    const std::vector<std::string> &phenoSpecs,
    const std::vector<std::string> &pcColNames,
    int nClusters,
    uint64_t seed,
    const GenoSpec &geno,
    const std::vector<std::string> &refAfFiles,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const int P = static_cast<int>(phenoSpecs.size());
    const int K = nClusters;
    const int nPop = static_cast<int>(refAfFiles.size());

    // ── Parse each spec into column names ───────────────────────────
    std::vector<std::vector<std::string> > phenoCols(P);
    std::vector<std::string> allPhenoCols;
    for (int p = 0; p < P; ++p) {
        phenoCols[p] = parsePhenoSpecLEAF(phenoSpecs[p]);
        for (const auto &col : phenoCols[p])
            allPhenoCols.push_back(col);
    }

    // ── Phase A: shared data loading ────────────────────────────────
    infoMsg("LEAF: Loading data (%d phenotypes, %d clusters)", P, K);
    auto famIIDs = parseGenoIIDs(geno);
    const uint32_t nFamTotal = static_cast<uint32_t>(famIIDs.size());
    SubjectData sdFull(famIIDs); // copy — need famIIDs later
    sdFull.loadPhenoFile(phenoFile);
    if (!covarFile.empty()) {
        std::vector<std::string> needed = covarNames;
        needed.insert(needed.end(), pcColNames.begin(), pcColNames.end());
        sdFull.loadCovar(covarFile, needed);
    }
    sdFull.setKeepRemove(keepFile, removeFile);
    if (!spgrmGrabFile.empty() || !spgrmGctaFile.empty())
        sdFull.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sdFull.famIIDs()));
    sdFull.setGenoLabel(geno.flagLabel());
    sdFull.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sdFull.finalize();
    sdFull.dropNaInColumns(allPhenoCols);

    const Eigen::Index N = static_cast<Eigen::Index>(sdFull.nUsed());
    infoMsg("  %lld subjects after intersection", static_cast<long long>(N));

    // Shared covariate matrix
    Eigen::MatrixXd covarMat;
    if (!covarNames.empty()) {
        covarMat = sdFull.getColumns(covarNames);
    } else if (sdFull.hasCovar()) {
        covarMat = sdFull.covar();
    }
    const int nCov = static_cast<int>(covarMat.cols());
    Eigen::MatrixXd logisticDesign;
    // Will be built per-cluster, not globally

    // K-means clustering on PCs (phenotype-independent)
    Eigen::MatrixXd PCs = sdFull.getColumns(pcColNames);
    infoMsg("  %d PCs for K-means clustering", static_cast<int>(PCs.cols()));
    infoMsg("K-means clustering into %d clusters...", K);
    Eigen::VectorXi clusterLabels = kmeansCluster(PCs, K, /*nstart=*/ 25, seed, nthreads);

    // Cluster member indices (into full used-subject array)
    std::vector<std::vector<uint32_t> > clusterMemberIdx(K);
    for (Eigen::Index i = 0; i < N; ++i)
        clusterMemberIdx[clusterLabels[i] - 1].push_back(static_cast<uint32_t>(i));

    // Build per-cluster bitmasks from fullMask + clusterLabels
    const size_t nMaskWords = (nFamTotal + 63) / 64;
    std::vector<std::vector<uint64_t> > clusterMasks(K);
    for (int c = 0; c < K; ++c)
        clusterMasks[c].assign(nMaskWords, 0);
    {
        const auto &fullMask = sdFull.usedMask();
        uint32_t usedIdx = 0;
        for (uint32_t f = 0; f < nFamTotal; ++f) {
            if (fullMask[f / 64] & (1ULL << (f % 64))) {
                int cl = clusterLabels[usedIdx] - 1;
                clusterMasks[cl][f / 64] |= 1ULL << (f % 64);
                ++usedIdx;
            }
        }
    }

    std::vector<uint32_t> clusterN(K, 0);
    for (int c = 0; c < K; ++c)
        for (size_t w = 0; w < nMaskWords; ++w)
            clusterN[c] += static_cast<uint32_t>(__builtin_popcountll(clusterMasks[c][w]));

    // Union mask and cluster→union index mapping
    std::vector<uint64_t> unionMask(nMaskWords, 0);
    for (int c = 0; c < K; ++c)
        for (size_t w = 0; w < nMaskWords; ++w)
            unionMask[w] |= clusterMasks[c][w];
    uint32_t nUnion = 0;
    for (size_t w = 0; w < nMaskWords; ++w)
        nUnion += static_cast<uint32_t>(__builtin_popcountll(unionMask[w]));

    std::vector<std::vector<uint32_t> > clusterIndices(K);
    for (int c = 0; c < K; ++c)
        clusterIndices[c].reserve(clusterN[c]);
    {
        uint32_t denseIdx = 0;
        for (size_t w = 0; w < nMaskWords; ++w) {
            uint64_t bits = unionMask[w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                uint64_t bitVal = uint64_t(1) << bit;
                for (int c = 0; c < K; ++c) {
                    if (clusterMasks[c][w] & bitVal) {
                        clusterIndices[c].push_back(denseIdx);
                        break;
                    }
                }
                ++denseIdx;
                bits &= bits - 1;
            }
        }
    }
    infoMsg("  Total subjects (union): %u", nUnion);

    // Shared genotype data with union mask
    auto genoData = makeGenoData(geno, unionMask, nFamTotal, nUnion, nSnpPerChunk);
    infoMsg("  %u subjects matched, %u markers", genoData->nSubjUsed(), genoData->nMarkers());

    // Shared ref-AF loading
    infoMsg("Loading %d reference AF files...", nPop);
    std::vector<std::vector<PopMatchedAF> > popMatched(nPop);
    for (int r = 0; r < nPop; ++r) {
        infoMsg("  Pop %d: %s", r + 1, refAfFiles[r].c_str());
        popMatched[r] = loadAndMatchRefAf(*genoData, refAfFiles[r]);
        infoMsg("    %zu markers matched", popMatched[r].size());
    }

    // Build intersection of markers across all populations
    struct MultiPopEntry {
        double af[8];
        double obs_ct[8];
    };

    std::unordered_map<uint64_t, MultiPopEntry> allPopMap;
    if (nPop > 0) {
        for (const auto &pm : popMatched[0]) {
            MultiPopEntry e{};
            e.af[0] = pm.af;
            e.obs_ct[0] = pm.obs_ct;
            allPopMap.emplace(pm.genoIndex, e);
        }
    }
    for (int r = 1; r < nPop; ++r) {
        std::unordered_map<uint64_t, MultiPopEntry> next;
        for (const auto &pm : popMatched[r]) {
            auto it = allPopMap.find(pm.genoIndex);
            if (it == allPopMap.end()) continue;
            auto entry = it->second;
            entry.af[r] = pm.af;
            entry.obs_ct[r] = pm.obs_ct;
            next.emplace(pm.genoIndex, entry);
        }
        allPopMap = std::move(next);
    }

    struct MatchedMultiPop {
        uint64_t genoIndex;
        double popAF[8];
        double popObsCt[8];
    };

    std::vector<MatchedMultiPop> matchedMulti;
    matchedMulti.reserve(allPopMap.size());
    for (auto &[gi, e] : allPopMap) {
        MatchedMultiPop mm{};
        mm.genoIndex = gi;
        for (int r = 0; r < nPop; ++r) {
            mm.popAF[r] = e.af[r];
            mm.popObsCt[r] = e.obs_ct[r];
        }
        matchedMulti.push_back(mm);
    }
    std::sort(matchedMulti.begin(), matchedMulti.end(),
              [](const auto &a, const auto &b) {
        return a.genoIndex < b.genoIndex;
    });
    const size_t nMatched = matchedMulti.size();
    infoMsg("  %zu markers present in all %d populations", nMatched, nPop);

    // Build per-cluster usedIIDs for GRM loading (phenotype-independent)
    std::vector<std::vector<std::string> > clusterUsedIIDs(K);
    for (int c = 0; c < K; ++c) {
        clusterUsedIIDs[c].reserve(clusterN[c]);
        for (uint32_t f = 0; f < nFamTotal; ++f)
            if (clusterMasks[c][f / 64] & (1ULL << (f % 64)))
                clusterUsedIIDs[c].push_back(famIIDs[f]);
    }

    // Load K GRMs once (phenotype-independent)
    std::vector<std::unique_ptr<SparseGRM> > clGRMs(K);
    if (!spgrmGctaFile.empty() || !spgrmGrabFile.empty()) {
        infoMsg("Loading %d per-cluster sparse GRMs...", K);
        for (int c = 0; c < K; ++c) {
            clGRMs[c] = std::make_unique<SparseGRM>(
                SparseGRM::load(spgrmGrabFile, spgrmGctaFile, clusterUsedIIDs[c], famIIDs));
            infoMsg("  Cluster %d: %u subjects, %zu non-zeros",
                    c + 1, clGRMs[c]->nSubjects(), clGRMs[c]->nnz());
        }
    }

    // ── Phase B: parallel null-model regression ─────────────────────
    // Per-phenotype, per-cluster: residuals, weights, indicator
    struct ClusterRWI {
        Eigen::VectorXd resid, weight, ind;
    };

    // clRWI[p][c]
    std::vector<std::vector<ClusterRWI> > clRWI(P, std::vector<ClusterRWI>(K));
    std::vector<std::string> traitNames(P);

    {
        const int totalTasks = P * K;
        const int nWorkers = std::min(nthreads, totalTasks);
        infoMsg("LEAF: Fitting %d × %d null models with %d threads", P, K, nWorkers);
        std::atomic<int> nextTask{0};
        std::vector<std::string> regrErrors(totalTasks);

        // Pre-compute traitNames (outside worker to avoid races).
        for (int p = 0; p < P; ++p) {
            const auto &cols = phenoCols[p];
            const bool isSurv = (cols.size() >= 2);
            traitNames[p] = isSurv ? cols[0] + "_" + cols[1] : cols[0];
        }

        auto regrWorker = [&]() {
            for (;;) {
                int t = nextTask.fetch_add(1, std::memory_order_relaxed);
                if (t >= totalTasks) break;
                int p = t / K;
                int c = t % K;
                try {
                    const auto &cols = phenoCols[p];
                    const bool isSurv = (cols.size() >= 2);

                    // Extract full indicator/survTime
                    Eigen::VectorXd fullInd, fullTime;
                    if (isSurv) {
                        fullTime = sdFull.getColumn(cols[0]);
                        fullInd = sdFull.getColumn(cols[1]);
                    } else {
                        fullInd = sdFull.getColumn(cols[0]);
                        // Validate binary: all values must be 0 or 1
                        for (Eigen::Index i = 0; i < fullInd.size(); ++i) {
                            double v = fullInd[i];
                            if (v != 0.0 && v != 1.0)
                                throw std::runtime_error(
                                          "Phenotype '" + cols[0] + "' is not binary"
                                          " (found value " + std::to_string(v) + ")."
                                          " For survival phenotypes use --pheno-name TIME:EVENT.");
                        }
                    }

                    const auto &members = clusterMemberIdx[c];
                    const Eigen::Index Nc = static_cast<Eigen::Index>(members.size());

                    // Subset for cluster
                    Eigen::VectorXd cInd(Nc), cTime(Nc);
                    for (Eigen::Index j = 0; j < Nc; ++j) {
                        cInd[j] = fullInd[members[j]];
                        if (isSurv) cTime[j] = fullTime[members[j]];
                    }

                    Eigen::VectorXd cWeight = regression::calRegrWeight(refPrevalence, cInd);

                    Eigen::VectorXd cResid;
                    if (isSurv) {
                        Eigen::MatrixXd cCov(Nc, covarMat.cols());
                        for (Eigen::Index j = 0; j < Nc; ++j)
                            cCov.row(j) = covarMat.row(members[j]);
                        cResid = regression::coxResiduals(cTime, cInd, cCov, cWeight);
                    } else {
                        Eigen::MatrixXd cDesign(Nc, 1 + nCov);
                        cDesign.col(0).setOnes();
                        if (nCov > 0) {
                            for (Eigen::Index j = 0; j < Nc; ++j)
                                cDesign.row(j).tail(nCov) = covarMat.row(members[j]);
                        }
                        cResid = regression::logisticResiduals(cInd, cDesign, cWeight);
                    }
                    clRWI[p][c] = {std::move(cResid), std::move(cWeight), std::move(cInd)};
                    infoMsg("  [%s cl%d] null model done", traitNames[p].c_str(), c + 1);
                } catch (const std::exception &ex) {
                    regrErrors[t] = ex.what();
                }
            }
        };

        std::vector<std::thread> threads;
        threads.reserve(nWorkers - 1);
        for (int t = 0; t < nWorkers - 1; ++t)
            threads.emplace_back(regrWorker);
        regrWorker();
        for (auto &th : threads) th.join();

        for (int t = 0; t < totalTasks; ++t)
            if (!regrErrors[t].empty()) {
                int p = t / K, c = t % K;
                throw std::runtime_error(
                          "LEAF null model failed for '" + traitNames[p] +
                          "' cluster " + std::to_string(c + 1) + ": " + regrErrors[t]);
            }
    }

    // ── Phase C: shared genotype scan + summix + AF synthesis ───────
    // One I/O pass: compute intAF (shared) + per-phenotype mu0/mu1
    infoMsg("LEAF: Scanning %zu matched markers for %d clusters × %d phenotypes",
            nMatched, K, P);

    // intAF: [K][nMatched] — phenotype-independent
    std::vector<std::vector<double> > clIntAF(K, std::vector<double>(nMatched));
    // Per-phenotype per-cluster scan: mu0/mu1/n0/n1/mu_int
    struct MarkerPhStat {
        double mu0, mu1, n0, n1, mu_int;
    };

    // phClStats[p][c][m]
    std::vector<std::vector<std::vector<MarkerPhStat> > > phClStats(
        P, std::vector<std::vector<MarkerPhStat> >(K, std::vector<MarkerPhStat>(nMatched)));

    {
        const uint32_t nFull = genoData->nSubjUsed();
        auto cursor = genoData->makeCursor();
        if (!matchedMulti.empty())
            cursor->beginSequentialBlock(matchedMulti.front().genoIndex);
        Eigen::VectorXd fullGVec(nFull);

        for (size_t m = 0; m < nMatched; ++m) {
            cursor->getGenotypesSimple(matchedMulti[m].genoIndex, fullGVec);

            for (int c = 0; c < K; ++c) {
                const auto &idx = clusterIndices[c];
                const size_t Nc = idx.size();

                // One pass over cluster members: compute phenotype-independent total
                // and per-phenotype case/control sums
                double totalSum = 0, totalCnt = 0;
                std::vector<double> sum0(P, 0), sum1(P, 0), cnt0(P, 0), cnt1(P, 0);

                for (size_t k = 0; k < Nc; ++k) {
                    double g = fullGVec[idx[k]];
                    if (std::isnan(g)) continue;
                    totalSum += g;
                    totalCnt += 1.0;
                    for (int p = 0; p < P; ++p) {
                        if (clRWI[p][c].ind[static_cast<Eigen::Index>(k)] == 1.0) {
                            sum1[p] += g;
                            cnt1[p] += 1.0;
                        } else {
                            sum0[p] += g;
                            cnt0[p] += 1.0;
                        }
                    }
                }

                clIntAF[c][m] = (totalCnt > 0)
                    ? (totalSum) / (2.0 * totalCnt)
                    : std::numeric_limits<double>::quiet_NaN();

                for (int p = 0; p < P; ++p) {
                    auto &st = phClStats[p][c][m];
                    st.mu0 = (cnt0[p] > 0) ? sum0[p] / (2.0 * cnt0[p]) : 0.0;
                    st.mu1 = (cnt1[p] > 0) ? sum1[p] / (2.0 * cnt1[p]) : 0.0;
                    st.n0 = cnt0[p];
                    st.n1 = cnt1[p];
                    st.mu_int = clIntAF[c][m];
                }
            }
        }
    }

    // Per-cluster summix (phenotype-independent — uses only intAF)
    // K calls are independent; parallelize with min(T, K) workers.
    std::vector<Eigen::VectorXd> clProportions(K);
    {
        const int nWorkers = std::min(nthreads, K);
        infoMsg("LEAF: Running summix for %d clusters with %d threads", K, nWorkers);
        std::atomic<int> nextCl{0};

        // Build refMat once (shared read-only across clusters).
        Eigen::MatrixXd refMat(nMatched, nPop);
        for (size_t m = 0; m < nMatched; ++m)
            for (int r = 0; r < nPop; ++r)
                refMat(m, r) = matchedMulti[m].popAF[r];

        auto summixWorker = [&]() {
            for (;;) {
                int c = nextCl.fetch_add(1, std::memory_order_relaxed);
                if (c >= K) break;
                Eigen::VectorXd obsAF(nMatched);
                for (size_t m = 0; m < nMatched; ++m)
                    obsAF[m] = clIntAF[c][m];
                clProportions[c] = summixEstimate(obsAF, refMat);
                for (int r = 0; r < nPop; ++r)
                    infoMsg("  Cluster %d pop%d: %.4f", c + 1, r + 1, clProportions[c][r]);
            }
        };

        std::vector<std::thread> threads;
        threads.reserve(nWorkers - 1);
        for (int t = 0; t < nWorkers - 1; ++t)
            threads.emplace_back(summixWorker);
        summixWorker();
        for (auto &th : threads) th.join();
    }

    // Per-cluster AF_ref and obs_ct synthesis (shared across phenotypes)
    // Then per-phenotype × per-cluster: build MatchedMarkerInfo with phenotype-specific mu0/mu1
    // clMatchedInfo[p][c] — vector of MatchedMarkerInfo
    std::vector<std::vector<std::vector<MatchedMarkerInfo> > > clMatchedInfo(
        P, std::vector<std::vector<MatchedMarkerInfo> >(K));

    for (int c = 0; c < K; ++c) {
        // Synthesize shared AF_ref and obs_ct for cluster c
        std::vector<double> sharedAFRef(nMatched);
        std::vector<double> sharedObsCt(nMatched);
        for (size_t m = 0; m < nMatched; ++m) {
            double af_ref = 0.0;
            for (int r = 0; r < nPop; ++r)
                af_ref += clProportions[c][r] * matchedMulti[m].popAF[r];
            sharedAFRef[m] = af_ref;

            double denom = 0.0;
            for (int r = 0; r < nPop; ++r) {
                double oc = matchedMulti[m].popObsCt[r];
                if (oc > 0 && clProportions[c][r] > 0)
                    denom += (clProportions[c][r] * clProportions[c][r]) / oc;
            }
            sharedObsCt[m] = (denom > 0) ? 1.0 / denom : 0.0;
        }

        for (int p = 0; p < P; ++p) {
            clMatchedInfo[p][c].resize(nMatched);
            for (size_t m = 0; m < nMatched; ++m) {
                auto &mi = clMatchedInfo[p][c][m];
                mi.genoIndex = matchedMulti[m].genoIndex;
                mi.AF_ref = sharedAFRef[m];
                mi.obs_ct = sharedObsCt[m];
                mi.mu0 = phClStats[p][c][m].mu0;
                mi.mu1 = phClStats[p][c][m].mu1;
                mi.n0 = phClStats[p][c][m].n0;
                mi.n1 = phClStats[p][c][m].n1;
                mi.mu_int = phClStats[p][c][m].mu_int;
            }
        }
    }

    // ── Phase D: parallel batch-effect testing ──────────────────────
    // P × K tasks, parallelized with min(T, P*K) workers
    // clRefMaps[p][c]
    std::vector<std::vector<std::shared_ptr<std::unordered_map<uint64_t, WtCoxGRefInfo> > > > clRefMaps(
        P, std::vector<std::shared_ptr<std::unordered_map<uint64_t, WtCoxGRefInfo> > >(K));
    {
        const int totalTasks = P * K;
        const int nWorkers = std::min(nthreads, totalTasks);
        infoMsg("LEAF: Batch-effect testing for %d phenotypes × %d clusters with %d threads",
                P, K, nWorkers);
        std::atomic<int> nextTask{0};
        std::vector<std::string> batchErrors(totalTasks);

        auto batchWorker = [&]() {
            for (;;) {
                int t = nextTask.fetch_add(1, std::memory_order_relaxed);
                if (t >= totalTasks) break;
                int p = t / K;
                int c = t % K;
                try {
                    clRefMaps[p][c] = testBatchEffects(
                        clMatchedInfo[p][c],
                        clRWI[p][c].resid,
                        clRWI[p][c].weight,
                        clRWI[p][c].ind,
                        clGRMs[c].get(),
                        refPrevalence,
                        cutoff);
                    infoMsg("  [%s cl%d] %zu markers retained",
                            traitNames[p].c_str(), c + 1, clRefMaps[p][c]->size());
                } catch (const std::exception &ex) {
                    batchErrors[t] = ex.what();
                }
            }
        };

        std::vector<std::thread> threads;
        threads.reserve(nWorkers - 1);
        for (int t = 0; t < nWorkers - 1; ++t)
            threads.emplace_back(batchWorker);
        batchWorker();
        for (auto &th : threads) th.join();

        for (int t = 0; t < totalTasks; ++t)
            if (!batchErrors[t].empty()) {
                int p = t / K, c = t % K;
                throw std::runtime_error(
                          "LEAF batch-effect failed for '" + traitNames[p] +
                          "' cluster " + std::to_string(c + 1) + ": " + batchErrors[t]);
            }
    }

    // ── Phase E: multi-phenotype marker engine ──────────────────────
    std::vector<PhenoTask> tasks(P);
    for (int p = 0; p < P; ++p) {
        std::vector<std::unique_ptr<WtCoxGMethod> > clMethods;
        clMethods.reserve(K);
        for (int c = 0; c < K; ++c) {
            clMethods.push_back(std::make_unique<WtCoxGMethod>(
                                    clRWI[p][c].resid, clRWI[p][c].weight,
                                    cutoff, spaCutoff, clRefMaps[p][c]));
        }
        auto method = std::make_unique<LEAFMethod>(std::move(clMethods), clusterIndices);
        tasks[p].phenoName = traitNames[p];
        tasks[p].method = std::move(method);
        tasks[p].unionToLocal.resize(genoData->nSubjUsed());
        std::iota(tasks[p].unionToLocal.begin(), tasks[p].unionToLocal.end(), 0u);
        tasks[p].nUsed = genoData->nSubjUsed();
    }

    infoMsg("LEAF: Starting multi-phenotype association (%d phenotypes, %d clusters, %d threads)",
            P, K, nthreads);
    multiPhenoEngine(
        *genoData,
        tasks,
        outPrefix,
        "LEAF",
        compression,
        compressionLevel,
        nthreads,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        hweCutoff
    );
}
