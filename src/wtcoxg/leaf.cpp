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
// writeKmeansClusterTsv — emit cluster labels for downstream inspection
//
// Writes <outPrefix>.KmeansCluster.tsv with header "#IID\tcluster".
// Rows follow the order of usedIIDs (= the order of clusterLabels).
// ======================================================================

static void writeKmeansClusterTsv(
    const std::string &outPrefix,
    const std::vector<std::string> &iids,
    const Eigen::VectorXi &clusterLabels
) {
    if (static_cast<Eigen::Index>(iids.size()) != clusterLabels.size())
        throw std::runtime_error(
            "writeKmeansClusterTsv: iids/clusterLabels size mismatch (" +
            std::to_string(iids.size()) + " vs " +
            std::to_string(clusterLabels.size()) + ")"
        );
    const std::string path = outPrefix + ".KmeansCluster.tsv";
    std::ofstream ofs(path);
    if (!ofs)
        throw std::runtime_error("writeKmeansClusterTsv: cannot open " + path);
    ofs << "#IID\tcluster\n";
    for (Eigen::Index i = 0; i < clusterLabels.size(); ++i)
        ofs << iids[static_cast<size_t>(i)] << '\t' << clusterLabels[i] << '\n';
    if (!ofs)
        throw std::runtime_error("writeKmeansClusterTsv: write failed for " + path);
    infoMsg("  Wrote cluster labels: %s", path.c_str());
}

// ======================================================================
// parseLeafClusterFile — read pre-computed cluster labels
//
// File layout (header required):
//   <IID-col>   <cluster-col>   [other columns ignored]
// The IID column is detected by name (`#IID`, `IID`, case-sensitive in
// that order); the cluster column by name (`cluster`, `Cluster`,
// `CLUSTER`).  Cluster values must be integers; the inferred K is
// max(cluster).  If `K_inout > 0`, the inferred K must match.
// ======================================================================

Eigen::VectorXi parseLeafClusterFile(
    const std::string &path,
    const std::vector<std::string> &usedIIDs,
    int &K_inout
) {
    std::ifstream ifs(path);
    if (!ifs)
        throw std::runtime_error("parseLeafClusterFile: cannot open " + path);

    std::string line;
    if (!std::getline(ifs, line))
        throw std::runtime_error("parseLeafClusterFile: empty file " + path);
    if (!line.empty() && line.back() == '\r') line.pop_back();

    // Split header on whitespace
    auto splitWS = [](const std::string &s) {
        std::vector<std::string> out;
        std::istringstream iss(s);
        std::string tok;
        while (iss >> tok) out.push_back(tok);
        return out;
    };
    std::vector<std::string> header = splitWS(line);
    int iidCol = -1, clusterCol = -1;
    for (int i = 0; i < (int)header.size(); ++i) {
        const std::string &h = header[i];
        if (iidCol < 0 && (h == "#IID" || h == "IID")) iidCol = i;
        if (clusterCol < 0 && (h == "cluster" || h == "Cluster" || h == "CLUSTER"))
            clusterCol = i;
    }
    if (iidCol < 0)
        throw std::runtime_error("parseLeafClusterFile: header must contain"
                                 " an IID column (`#IID` or `IID`) in " + path);
    if (clusterCol < 0)
        throw std::runtime_error("parseLeafClusterFile: header must contain"
                                 " a cluster column (`cluster`, `Cluster`,"
                                 " or `CLUSTER`) in " + path);

    std::unordered_map<std::string, int> iidToCluster;
    iidToCluster.reserve(usedIIDs.size() * 2);
    int observedMax = 0;
    size_t lineno = 1;
    while (std::getline(ifs, line)) {
        ++lineno;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> cols = splitWS(line);
        if ((int)cols.size() <= std::max(iidCol, clusterCol)) {
            throw std::runtime_error("parseLeafClusterFile: " + path +
                                     " line " + std::to_string(lineno) +
                                     " has too few columns");
        }
        const std::string &iid = cols[iidCol];
        const std::string &cStr = cols[clusterCol];
        if (iid.empty() || cStr.empty()) continue;
        // Parse as integer; reject non-integer values (e.g. "1.5")
        char *endp = nullptr;
        long c = std::strtol(cStr.c_str(), &endp, 10);
        if (endp == cStr.c_str() || (endp && *endp != '\0')) {
            throw std::runtime_error("parseLeafClusterFile: " + path +
                                     " line " + std::to_string(lineno) +
                                     ": cluster value '" + cStr +
                                     "' is not an integer");
        }
        if (c < 1) {
            throw std::runtime_error("parseLeafClusterFile: " + path +
                                     " line " + std::to_string(lineno) +
                                     ": cluster value " + std::to_string(c) +
                                     " is < 1");
        }
        auto [it, inserted] = iidToCluster.emplace(iid, (int)c);
        if (!inserted) {
            throw std::runtime_error("parseLeafClusterFile: " + path +
                                     " line " + std::to_string(lineno) +
                                     ": duplicate IID '" + iid + "'");
        }
        if ((int)c > observedMax) observedMax = (int)c;
    }
    if (observedMax < 2)
        throw std::runtime_error("parseLeafClusterFile: " + path +
                                 " produced fewer than 2 distinct clusters");

    // Cross-check inferred K against caller's hint
    const int K_inferred = observedMax;
    if (K_inout > 0 && K_inout != K_inferred) {
        throw std::runtime_error("parseLeafClusterFile: --leaf-nclusters=" +
                                 std::to_string(K_inout) +
                                 " but the cluster file has max(cluster)=" +
                                 std::to_string(K_inferred));
    }
    K_inout = K_inferred;

    // Build the dense N-vector aligned to usedIIDs
    Eigen::VectorXi labels(usedIIDs.size());
    std::vector<std::string> missing;
    for (size_t i = 0; i < usedIIDs.size(); ++i) {
        auto it = iidToCluster.find(usedIIDs[i]);
        if (it == iidToCluster.end()) {
            if (missing.size() < 5) missing.push_back(usedIIDs[i]);
            else if (missing.size() == 5) missing.push_back("…");
            labels[i] = 0; // sentinel; will throw below
        } else {
            labels[i] = it->second;
        }
    }
    if (!missing.empty()) {
        std::string ex = missing.front();
        for (size_t i = 1; i < missing.size(); ++i) ex += ", " + missing[i];
        throw std::runtime_error("parseLeafClusterFile: " + path +
                                 " is missing entries for " +
                                 std::to_string(missing.size()) +
                                 "+ subjects in the analysis set: " + ex);
    }
    return labels;
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

// Inner solver: constrained LS  min ‖D_active · p − o‖²
//   s.t. p ≥ 0, Σ p = 1, with the active reference set restricted to the
//   `allowed` columns of refAF.  Returns a length-K vector with 0 in the
//   non-allowed positions.  This is the body of Summix's `summix_calc`
//   step; the caller (`summixEstimate`) wraps it for the 1% drop pass.
static Eigen::VectorXd summixCalcConstrainedLS(
    const Eigen::MatrixXd &DtD_full,
    const Eigen::VectorXd &Dto_full,
    double oTo,
    int K,
    const std::vector<int> &allowed
) {
    const int Ka = static_cast<int>(allowed.size());
    Eigen::VectorXd bestP = Eigen::VectorXd::Zero(K);
    if (Ka == 0) return bestP;

    double bestObj = std::numeric_limits<double>::infinity();
    bool foundFeasible = false;

    const int nSubsets = (1 << Ka);
    for (int mask = 1; mask < nSubsets; ++mask) {
        int s = 0;
        int idx[8];
        for (int j = 0; j < Ka; ++j)
            if (mask & (1 << j)) idx[s++] = allowed[j];

        const int dim = s + 1;
        Eigen::MatrixXd KKT = Eigen::MatrixXd::Zero(dim, dim);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dim);
        for (int i = 0; i < s; ++i) {
            for (int j = 0; j < s; ++j)
                KKT(i, j) = 2.0 * DtD_full(idx[i], idx[j]);
            KKT(i, s) = 1.0;
            KKT(s, i) = 1.0;
            rhs[i] = 2.0 * Dto_full[idx[i]];
        }
        rhs[s] = 1.0;

        Eigen::VectorXd sol = KKT.partialPivLu().solve(rhs);

        bool feasible = true;
        for (int i = 0; i < s; ++i)
            if (sol[i] < -1e-12) { feasible = false; break; }
        if (!feasible) continue;

        for (int i = 0; i < s; ++i)
            if (sol[i] < 0.0) sol[i] = 0.0;

        Eigen::VectorXd pS = sol.head(s);
        double obj = 0.0;
        for (int i = 0; i < s; ++i)
            for (int j = 0; j < s; ++j)
                obj += pS[i] * DtD_full(idx[i], idx[j]) * pS[j];
        for (int i = 0; i < s; ++i)
            obj -= 2.0 * Dto_full[idx[i]] * pS[i];
        obj += oTo;

        if (obj < bestObj) {
            bestObj = obj;
            foundFeasible = true;
            bestP.setZero();
            for (int i = 0; i < s; ++i)
                bestP[idx[i]] = pS[i];
        }
    }

    if (!foundFeasible) {
        // Degenerate: fall back to uniform over the allowed set.
        for (int k : allowed) bestP[k] = 1.0 / static_cast<double>(Ka);
    }
    return bestP;
}

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

    Eigen::MatrixXd DtD = D.transpose() * D;
    Eigen::VectorXd Dto = D.transpose() * obs;
    double oTo = obs.squaredNorm();

    // Pass 1: solve the constrained LS over all K references.  Matches
    // Summix's first `summix_calc` call inside `summix()`.
    std::vector<int> allAllowed(K);
    std::iota(allAllowed.begin(), allAllowed.end(), 0);
    Eigen::VectorXd pPass1 = summixCalcConstrainedLS(DtD, Dto, oTo, K, allAllowed);

    // Pass 2: Summix's `<1%` reference-removal step.  Drop any reference
    // whose pass-1 proportion is below 0.01 and refit on the remainder.
    // (Summix.R lines 317-329: `if(sum(globTest[1,5:ncol(globTest)] < 0.01) > 0)
    //  ... reference <- reference[-toremove] ; sum_res <- summix_calc(...)`).
    std::vector<int> kept;
    kept.reserve(K);
    for (int k = 0; k < K; ++k)
        if (pPass1[k] >= 0.01) kept.push_back(k);

    if (kept.empty() || static_cast<int>(kept.size()) == K) {
        // No reference fell below 1% — or all did (degenerate).  Use pass-1.
        return pPass1;
    }
    return summixCalcConstrainedLS(DtD, Dto, oTo, K, kept);
}

// ======================================================================
// LEAFMethod — MethodBase implementation
// ======================================================================

LEAFMethod::LEAFMethod(
    std::vector<std::unique_ptr<WtCoxGMethod> > clusterMethods,
    std::vector<std::vector<uint32_t> > clusterIndices
) {
    auto g = std::make_shared<ClusterGeometry>();
    g->nCluster = static_cast<int>(clusterMethods.size());
    g->clusterIndices = std::move(clusterIndices);
    g->clusterN.resize(g->nCluster);
    for (int c = 0; c < g->nCluster; ++c)
        g->clusterN[c] = static_cast<int>(g->clusterIndices[c].size());

    m_nCluster = g->nCluster;
    m_clusterMethods = std::move(clusterMethods);
    m_clusterGVec.resize(m_nCluster);
    for (int c = 0; c < m_nCluster; ++c)
        m_clusterGVec[c].resize(g->clusterIndices[c].size());
    m_geom = std::move(g);
}

LEAFMethod::LEAFMethod(
    std::vector<std::unique_ptr<WtCoxGMethod> > clusterMethods,
    std::shared_ptr<const ClusterGeometry> geom
)
    : m_nCluster(geom->nCluster),
      m_clusterMethods(std::move(clusterMethods)),
      m_geom(std::move(geom))
{
    m_clusterGVec.resize(m_nCluster);
    for (int c = 0; c < m_nCluster; ++c)
        m_clusterGVec[c].resize(m_geom->clusterIndices[c].size());
}

std::unique_ptr<MethodBase> LEAFMethod::clone() const {
    std::vector<std::unique_ptr<WtCoxGMethod> > cloned;
    cloned.reserve(m_nCluster);
    for (const auto &m : m_clusterMethods) {
        auto p = m->clone(); // returns unique_ptr<MethodBase>
        cloned.push_back(std::unique_ptr<WtCoxGMethod>(static_cast<WtCoxGMethod *>(p.release())));
    }
    // Share the cluster geometry by shared_ptr copy — no per-worker
    // duplication of clusterIndices / clusterN.
    return std::make_unique<LEAFMethod>(std::move(cloned), m_geom);
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
        const auto &idx = m_geom->clusterIndices[c];
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

// ── Fused-GEMM hooks ────────────────────────────────────────────────
//
// LEAFMethod is itself a phenotype (the engine sees one task per LEAF
// trait), but internally each marker is tested K times — once per
// cluster — and the K cluster results are then meta-analysed.  The
// fused GEMM batches all K · B inner products into a single matmul:
//   AugResid: union_N × (K residual + K mask) cols
//   GEMM     : (2K × union_N) · (union_N × B)  →  (2K × B)
// giving per-marker (R_c.dot(g_c), gSum_c) for every cluster c.

void LEAFMethod::fillUnionResiduals(
    Eigen::Ref<Eigen::MatrixXd> dest,
    const std::vector<uint32_t> &unionToLocal
) const {
    // unionToLocal is the LEAF-level union mapping (= LEAF's pheno mask
    // applied to the engine's union, which here is identical because
    // there is exactly one LEAF task).  Each entry maps union index →
    // LEAFMethod's pheno-local index (= union index for LEAF).
    //
    // m_geom->clusterIndices[c] holds union-level indices of cluster c's
    // members.  Use those directly for the scatter.
    const int K = m_nCluster;
    for (int c = 0; c < K; ++c) {
        const auto &idx = m_geom->clusterIndices[c];
        const Eigen::VectorXd &R = m_clusterMethods[c]->residuals();
        for (size_t k = 0; k < idx.size(); ++k) {
            const uint32_t u = idx[k];
            // dest is indexed by engine union; LEAF's union equals the
            // engine union because LEAF is the only fuseable phenotype.
            if (unionToLocal[u] != UINT32_MAX) {
                dest(u, c)     = R[k];
                dest(u, K + c) = 1.0;
            }
        }
    }
}

void LEAFMethod::fillResidualSums(double *dest) const {
    const int K = m_nCluster;
    for (int c = 0; c < K; ++c) {
        dest[c]     = m_clusterMethods[c]->residuals().sum();
        dest[K + c] = static_cast<double>(m_geom->clusterN[c]);
    }
}

void LEAFMethod::processScoreBatch(
    const Eigen::Ref<const Eigen::MatrixXd> &scores,
    const double *gSums,
    const double *gSumSqs,
    uint32_t nUsed,
    const std::vector<double> &altFreqs,
    const std::vector<int> &chunkIdxs,
    std::vector<std::vector<double> > &results
) {
    (void)gSums; (void)gSumSqs; (void)nUsed; (void)altFreqs;

    const int B = static_cast<int>(scores.cols());
    const int K = m_nCluster;
    results.resize(B);

    std::vector<double> pExt(K), pNoext(K), sExt(K), sNoext(K);

    auto metaP = [](const std::vector<double> &scoresArr, const std::vector<double> &pvals) -> double {
        double sumScore = 0.0, sumVar = 0.0;
        for (size_t c = 0; c < scoresArr.size(); ++c) {
            if (std::isnan(scoresArr[c]) || std::isnan(pvals[c]) ||
                pvals[c] <= 0.0 || pvals[c] >= 1.0) continue;
            double chisq = math::qchisq(pvals[c], 1.0, false, false);
            if (chisq < 1e-30) chisq = 1e-30;
            double var = (scoresArr[c] * scoresArr[c]) / chisq;
            if (std::isnan(var)) continue;
            sumScore += scoresArr[c];
            sumVar += var;
        }
        if (sumVar <= 0.0) return std::numeric_limits<double>::quiet_NaN();
        double z = sumScore / std::sqrt(sumVar);
        return std::erfc(std::fabs(z) / std::sqrt(2.0));
    };

    for (int b = 0; b < B; ++b) {
        const int chunkIdx = chunkIdxs[b];

        for (int c = 0; c < K; ++c) {
            const double R_dot_g = scores(c, b);
            const double gSum    = scores(K + c, b);
            auto dr = m_clusterMethods[c]->computeDualFromScalars(
                R_dot_g, gSum, m_geom->clusterN[c], chunkIdx);
            pExt[c]   = dr.p_ext;
            pNoext[c] = dr.p_noext;
            sExt[c]   = dr.score_ext;
            sNoext[c] = dr.score_noext;
        }

        auto &r = results[b];
        r.clear();
        r.reserve(2 + 3 * K);
        r.push_back(metaP(sExt, pExt));
        r.push_back(metaP(sNoext, pNoext));
        for (int c = 0; c < K; ++c) {
            r.push_back(pExt[c]);
            r.push_back(pNoext[c]);
            r.push_back(m_clusterMethods[c]->chunkRefInfoAt(chunkIdx).pvalue_bat);
        }
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

#include "util/null_model.hpp"
#include "util/regression.hpp"
#include "wtcoxg/regression.hpp"

// LEAF.R's weight definition:  weight = ifelse(Event==1, 1, (1-prev)/prev)
// — a fixed (1-prev)/prev for every control, regardless of sample
// case/control ratio.  This differs from regression::calRegrWeight()
// (used by other methods), which divides the sample case/control ratio
// by the population ratio.  Use LEAF.R's form so that residuals and w1
// downstream of testBatchEffects match the reference R pipeline exactly.
static Eigen::VectorXd leafRegrWeight(
    double prevalence,
    const Eigen::Ref<const Eigen::VectorXd> &indicator
) {
    const double wCtrl = (1.0 - prevalence) / prevalence;
    Eigen::VectorXd w(indicator.size());
    for (Eigen::Index i = 0; i < indicator.size(); ++i)
        w[i] = (indicator[i] > 0.5) ? 1.0 : wCtrl;
    return w;
}

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
    double outlierRatio,
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
    writeKmeansClusterTsv(outPrefix, sdFull.usedIIDs(), clusterLabels);

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

    // Fit ONE global null model on all subjects (matching LEAF.R, where
    //   obj.wglm = glm(Event ~ SEX+AGE+PC1..4, weight = weight, family = binomial)
    // is computed on the full pheno table; per-cluster residuals are then
    // a slice of obj.wglm$y - obj.wglm$fitted.values).  This shares the
    // covariate coefficients across clusters, in contrast to a per-cluster
    // refit which would absorb ancestry effects into each cluster's local
    // baseline.
    infoMsg("Fitting global null model on %lld subjects...", static_cast<long long>(N));
    Eigen::VectorXd fullWeight = leafRegrWeight(refPrevalence, indicator);
    Eigen::VectorXd fullResid;
    if (isSurv) {
        fullResid = regression::coxResiduals(survTime, indicator, covarMat, fullWeight);
    } else {
        fullResid = regression::logisticResiduals(indicator, logisticDesign, fullWeight);
    }

    for (int c = 0; c < nCluster; ++c) {
        const auto &members = clusterMemberIdx[c];
        const Eigen::Index Nc = static_cast<Eigen::Index>(members.size());
        Eigen::VectorXd cResid(Nc), cWeight(Nc), cInd(Nc);
        for (Eigen::Index j = 0; j < Nc; ++j) {
            const auto i = members[j];
            cResid[j] = fullResid[i];
            cWeight[j] = fullWeight[i];
            cInd[j] = indicator[i];
        }
        infoMsg("  Cluster %d: %d subjects (global residuals subset)", c + 1, static_cast<int>(Nc));
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

    // Per-cluster batch-effect testing.  Pass the GLOBAL weight sum (across
    // all clusters) as the w1 normaliser so that var_ratio_w0 mirrors
    // LEAF.R's `weight1 = weight/(2*sum(pheno$weight))` semantics.
    infoMsg("Per-cluster batch-effect testing...");

    std::vector<std::shared_ptr<std::unordered_map<uint64_t, WtCoxGRefInfo> > > clRefMaps(nCluster);

    const double globalSumWeight = fullWeight.sum();

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
            cutoff,
            globalSumWeight
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
                outlierRatio,
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

// PhenoSpec parsing delegated to nullmodel::parsePhenoSpecAuto.  Both
// binary ("COL") and survival ("TIME:EVENT") forms are inferred per token,
// matching the WtCoxG / SPACox / SPAGRM / SPAmix syntax.

void runLEAF(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    const std::vector<nullmodel::PhenoSpec> &parsedSpecs,
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
    double outlierRatio,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &clusterFile
) {
    const int P = static_cast<int>(parsedSpecs.size());
    const int nPop = static_cast<int>(refAfFiles.size());
    int K = nClusters; // may be 0 when --leaf-cluster-file infers it

    // ── Collect the union of phenotype columns the regression will touch ──
    std::vector<std::string> allPhenoCols;
    for (int p = 0; p < P; ++p) {
        if (nullmodel::isCoxSpec(parsedSpecs[p])) {
            allPhenoCols.push_back(parsedSpecs[p].timeColumn);
            allPhenoCols.push_back(parsedSpecs[p].eventColumn);
        } else {
            allPhenoCols.push_back(parsedSpecs[p].yColumn);
        }
    }

    // ── Phase A: shared data loading ────────────────────────────────
    if (K > 0) infoMsg("LEAF: Loading data (%d phenotypes, %d clusters)", P, K);
    else       infoMsg("LEAF: Loading data (%d phenotypes, K inferred from --leaf-cluster-file)", P);
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

    // Cluster labels: either read from --leaf-cluster-file or run K-means on PCs.
    Eigen::VectorXi clusterLabels;
    if (!clusterFile.empty()) {
        infoMsg("LEAF: Reading cluster labels from %s", clusterFile.c_str());
        clusterLabels = parseLeafClusterFile(clusterFile, sdFull.usedIIDs(), K);
        infoMsg("  %d clusters, %lld subjects assigned", K,
                static_cast<long long>(clusterLabels.size()));
    } else {
        Eigen::MatrixXd PCs = sdFull.getColumns(pcColNames);
        infoMsg("  %d PCs for K-means clustering", static_cast<int>(PCs.cols()));
        infoMsg("K-means clustering into %d clusters...", K);
        clusterLabels = kmeansCluster(PCs, K, /*nstart=*/ 25, seed, nthreads);
        writeKmeansClusterTsv(outPrefix, sdFull.usedIIDs(), clusterLabels);
    }

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
// Per-phenotype global weight sum (filled inside Phase B; consumed by
// Phase D as the w1 normaliser, matching LEAF.R's global `weight1`).
    std::vector<double> globalSumWeightPerPheno(P, 0.0);

    {
        const int totalTasks = P * K;
        const int nWorkers = std::min(nthreads, totalTasks);
        infoMsg("LEAF: Fitting %d × %d null models with %d threads", P, K, nWorkers);
        std::atomic<int> nextTask{0};
        std::vector<std::string> regrErrors(totalTasks);

        // Pre-compute traitNames and per-phenotype indicator/time vectors with
        // inference + recode applied once (instead of P × K times inside the
        // worker).  This also keeps inference log lines unique per phenotype.
        std::vector<Eigen::VectorXd> fullIndPre(P);
        std::vector<Eigen::VectorXd> fullTimePre(P);
        std::vector<char> isSurvPre(P);
        const std::vector<std::string> unionIIDsForInfer = sdFull.usedIIDs();
        for (int p = 0; p < P; ++p) {
            const auto &spec = parsedSpecs[p];
            traitNames[p] = spec.name;
            const bool isSurv = nullmodel::isCoxSpec(spec);
            isSurvPre[p] = isSurv ? 1 : 0;
            if (isSurv) {
                fullTimePre[p] = sdFull.getColumn(spec.timeColumn);
                fullIndPre[p] = sdFull.getColumn(spec.eventColumn);
                nullmodel::inferCoxSurvival(
                    fullTimePre[p], fullIndPre[p],
                    spec.timeColumn, spec.eventColumn,
                    unionIIDsForInfer);
            } else {
                fullIndPre[p] = sdFull.getColumn(spec.yColumn);
                auto info = nullmodel::inferModelFromColumn(
                    fullIndPre[p], spec.yColumn, unionIIDsForInfer);
                if (info.model != nullmodel::RegressionModel::Logistic)
                    throw std::runtime_error(
                        "LEAF requires a binary phenotype for '" +
                        spec.yColumn + "' but inference returned " +
                        nullmodel::regressionModelName(info.model) +
                        " (use --pheno-name TIME:EVENT for survival)");
                if (info.needRecode) {
                    double v0 = info.sortedDistinct[0];
                    double v1 = info.sortedDistinct[1];
                    for (Eigen::Index i = 0; i < fullIndPre[p].size(); ++i) {
                        double v = fullIndPre[p][i];
                        if (std::isnan(v)) continue;
                        fullIndPre[p][i] = (v == v0) ? 0.0 : 1.0;
                    }
                    infoMsg("  Recoded binary column '%s': {%g, %g} -> {0, 1}",
                            spec.yColumn.c_str(), v0, v1);
                }
            }
        }

        // LEAF.R fits a SINGLE global glm/coxph per phenotype on all
        // 42 640 subjects, then takes per-cluster slices of the resulting
        // residual vector.  Parallelise across phenotypes (P tasks); the
        // K cluster subsets are then a cheap pointer-scatter inside each
        // worker.  This matches LEAF.R and replaces the prior per-cluster
        // refit (which absorbed ancestry into each cluster's local baseline).
        std::vector<Eigen::VectorXd> fullResidPre(P);
        std::vector<Eigen::VectorXd> fullWeightPre(P);
        Eigen::MatrixXd logisticDesignFull;
        {
            const bool anyLogistic = std::any_of(
                isSurvPre.begin(), isSurvPre.end(),
                [](char s) { return s == 0; });
            if (anyLogistic) {
                logisticDesignFull.resize(N, 1 + nCov);
                logisticDesignFull.col(0).setOnes();
                if (nCov > 0) logisticDesignFull.rightCols(nCov) = covarMat;
            }
        }

        const int nGlobalWorkers = std::min(nthreads, P);
        infoMsg("LEAF: Fitting %d global null model(s) with %d threads", P, nGlobalWorkers);
        std::atomic<int> nextGlobal{0};
        std::vector<std::string> globalErrors(P);
        auto globalWorker = [&]() {
            for (;;) {
                int p = nextGlobal.fetch_add(1, std::memory_order_relaxed);
                if (p >= P) break;
                try {
                    const bool isSurv = (isSurvPre[p] != 0);
                    fullWeightPre[p] = leafRegrWeight(refPrevalence, fullIndPre[p]);
                    if (isSurv) {
                        fullResidPre[p] = regression::coxResiduals(
                            fullTimePre[p], fullIndPre[p], covarMat, fullWeightPre[p]);
                    } else {
                        fullResidPre[p] = regression::logisticResiduals(
                            fullIndPre[p], logisticDesignFull, fullWeightPre[p]);
                    }
                    infoMsg("  [%s] global null model done", traitNames[p].c_str());
                } catch (const std::exception &ex) {
                    globalErrors[p] = ex.what();
                }
            }
        };
        {
            std::vector<std::thread> ths;
            ths.reserve(nGlobalWorkers - 1);
            for (int t = 0; t < nGlobalWorkers - 1; ++t)
                ths.emplace_back(globalWorker);
            globalWorker();
            for (auto &th : ths) th.join();
        }
        for (int p = 0; p < P; ++p)
            if (!globalErrors[p].empty())
                throw std::runtime_error(
                          "LEAF global null model failed for '" + traitNames[p] +
                          "': " + globalErrors[p]);

        // Slice per-cluster residuals/weights/indicator from the global fit.
        for (int p = 0; p < P; ++p) {
            const Eigen::VectorXd &fullResid = fullResidPre[p];
            const Eigen::VectorXd &fullWeight = fullWeightPre[p];
            const Eigen::VectorXd &fullInd = fullIndPre[p];
            globalSumWeightPerPheno[p] = fullWeight.sum();
            for (int c = 0; c < K; ++c) {
                const auto &members = clusterMemberIdx[c];
                const Eigen::Index Nc = static_cast<Eigen::Index>(members.size());
                Eigen::VectorXd cResid(Nc), cWeight(Nc), cInd(Nc);
                for (Eigen::Index j = 0; j < Nc; ++j) {
                    const auto i = members[j];
                    cResid[j]  = fullResid[i];
                    cWeight[j] = fullWeight[i];
                    cInd[j]    = fullInd[i];
                }
                clRWI[p][c] = {std::move(cResid), std::move(cWeight), std::move(cInd)};
            }
        }
        // Silence unused-variable warnings from the prior worker scaffolding.
        (void)regrErrors;
        (void)totalTasks;
        (void)nWorkers;
        (void)nextTask;
    }

    // ── Phase C: per-phenotype × per-cluster scan + summix + AF synthesis ─
    // One I/O pass computes per-(p,c) mu0/mu1/n0/n1.  Summix observed AF is
    // the equal-blend `0.5*(mu0+mu1)` (matches LEAF.R `mu.target`); the
    // folded form is stored as `mu_int` for MAF-group binning (matches
    // LEAF.R `mu.int = ifelse(mu.target>0.5, 1-mu.target, mu.target)`).
    // Because the half-blend depends on the phenotype's case/control split,
    // summix proportions and the synthesised AF_ref / obs_ct are per (p,c).
    infoMsg("LEAF: Scanning %zu matched markers for %d clusters × %d phenotypes",
            nMatched, K, P);

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

                std::vector<double> sum0(P, 0), sum1(P, 0), cnt0(P, 0), cnt1(P, 0);

                for (size_t k = 0; k < Nc; ++k) {
                    double g = fullGVec[idx[k]];
                    if (std::isnan(g)) continue;
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

                for (int p = 0; p < P; ++p) {
                    auto &st = phClStats[p][c][m];
                    st.mu0 = (cnt0[p] > 0) ? sum0[p] / (2.0 * cnt0[p]) : 0.0;
                    st.mu1 = (cnt1[p] > 0) ? sum1[p] / (2.0 * cnt1[p]) : 0.0;
                    st.n0 = cnt0[p];
                    st.n1 = cnt1[p];
                    if (cnt0[p] > 0 && cnt1[p] > 0) {
                        double mu_target = 0.5 * (st.mu0 + st.mu1);
                        st.mu_int = (mu_target > 0.5) ? (1.0 - mu_target) : mu_target;
                    } else {
                        st.mu_int = std::numeric_limits<double>::quiet_NaN();
                    }
                }
            }
        }
    }

    // Per-phenotype × per-cluster summix.  P*K tasks parallelised with min(T, P*K).
    std::vector<std::vector<Eigen::VectorXd> > clProportions(
        P, std::vector<Eigen::VectorXd>(K));
    {
        const int totalTasks = P * K;
        const int nWorkers = std::min(nthreads, totalTasks);
        infoMsg("LEAF: Running summix for %d phenotypes × %d clusters with %d threads",
                P, K, nWorkers);
        std::atomic<int> nextTask{0};

        // Shared refMat (nMatched × nPop)
        Eigen::MatrixXd refMat(nMatched, nPop);
        for (size_t m = 0; m < nMatched; ++m)
            for (int r = 0; r < nPop; ++r)
                refMat(m, r) = matchedMulti[m].popAF[r];

        auto summixWorker = [&]() {
            for (;;) {
                int t = nextTask.fetch_add(1, std::memory_order_relaxed);
                if (t >= totalTasks) break;
                int p = t / K;
                int c = t % K;
                Eigen::VectorXd obsAF(nMatched);
                for (size_t m = 0; m < nMatched; ++m) {
                    const auto &st = phClStats[p][c][m];
                    // Markers with no cases or no controls in the cluster
                    // get NaN; summixEstimate filters them out.
                    obsAF[m] = (st.n0 > 0 && st.n1 > 0)
                        ? 0.5 * (st.mu0 + st.mu1)
                        : std::numeric_limits<double>::quiet_NaN();
                }
                clProportions[p][c] = summixEstimate(obsAF, refMat);
            }
        };

        std::vector<std::thread> threads;
        threads.reserve(nWorkers - 1);
        for (int t = 0; t < nWorkers - 1; ++t)
            threads.emplace_back(summixWorker);
        summixWorker();
        for (auto &th : threads) th.join();

        for (int p = 0; p < P; ++p)
            for (int c = 0; c < K; ++c)
                for (int r = 0; r < nPop; ++r)
                    infoMsg("  [%s cl%d] pop%d: %.4f",
                            traitNames[p].c_str(), c + 1, r + 1,
                            clProportions[p][c][r]);
    }

    // Per-(p,c) AF_ref and obs_ct synthesis from per-(p,c) proportions.
    std::vector<std::vector<std::vector<MatchedMarkerInfo> > > clMatchedInfo(
        P, std::vector<std::vector<MatchedMarkerInfo> >(K));
    for (int p = 0; p < P; ++p) {
        for (int c = 0; c < K; ++c) {
            clMatchedInfo[p][c].resize(nMatched);
            const auto &props = clProportions[p][c];
            for (size_t m = 0; m < nMatched; ++m) {
                double af_ref = 0.0;
                for (int r = 0; r < nPop; ++r)
                    af_ref += props[r] * matchedMulti[m].popAF[r];

                double denom = 0.0;
                for (int r = 0; r < nPop; ++r) {
                    double oc = matchedMulti[m].popObsCt[r];
                    if (oc > 0 && props[r] > 0)
                        denom += (props[r] * props[r]) / oc;
                }
                double obs_ct_loc = (denom > 0) ? 1.0 / denom : 0.0;

                auto &mi = clMatchedInfo[p][c][m];
                mi.genoIndex = matchedMulti[m].genoIndex;
                mi.AF_ref    = af_ref;
                mi.obs_ct    = obs_ct_loc;
                mi.mu0       = phClStats[p][c][m].mu0;
                mi.mu1       = phClStats[p][c][m].mu1;
                mi.n0        = phClStats[p][c][m].n0;
                mi.n1        = phClStats[p][c][m].n1;
                mi.mu_int    = phClStats[p][c][m].mu_int;
            }
        }
    }

    // ── Phase D: parallel batch-effect testing ──────────────────────
    // P × K tasks, parallelized with min(T, P*K) workers.  Pass each
    // phenotype's global weight sum so that w1 = weight / (2 · globalSum)
    // matches LEAF.R's `weight1` (defined on the full pheno table).
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
                        cutoff,
                        globalSumWeightPerPheno[p]);
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
                                    cutoff, spaCutoff, outlierRatio, clRefMaps[p][c]));
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
