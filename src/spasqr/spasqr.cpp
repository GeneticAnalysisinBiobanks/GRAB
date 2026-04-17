// spasqr.cpp — SPAsqr: SPA-squared multi-tau marker association (pure C++17 / Eigen / Boost)

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "spasqr/spasqr.hpp"
#include "engine/loco.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "spagrm/grm_null.hpp"
#include "spagrm/spagrm.hpp"
#include "spasqr/conquer.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <Eigen/Dense>

// ══════════════════════════════════════════════════════════════════════
// SPAsqrMethod — MethodBase adapter wrapping ntaus SPAGRMClass instances
// ══════════════════════════════════════════════════════════════════════

namespace {

class SPAsqrMethod : public MethodBase {
  public:
// Constructor with tau labels (used by pheno path):
//   labels like "tau0.1", "tau0.3", ...
//   Output order: P_CCT  P_tau0.1 ... P_tau0.9  Z_tau0.1 ... Z_tau0.9
    SPAsqrMethod(
        int ntaus,
        std::vector<SPAGRMClass> spagrm_vec,
        std::vector<std::string> tauLabels
    )
        : m_ntaus(ntaus),
          m_spagrm_vec(std::move(spagrm_vec)),
          m_tauLabels(std::move(tauLabels)),
          m_hasLabels(true)
    {
        buildResidMat();
    }

    std::unique_ptr<MethodBase> clone() const override {
        // SPAGRMClass copies share read-only data via shared_ptr.
        // m_methodShared (residMat + residSums) is also shared, not deep-copied.
        return std::make_unique<SPAsqrMethod>(*this);
    }

    int resultSize() const override {
        return 2 * m_ntaus + 1;
    }

    std::string getHeaderColumns() const override {
        std::ostringstream oss;
        if (m_hasLabels) {
            // Pheno path: P_CCT first, then P_tau{val}..., then Z_tau{val}...
            oss << "\tP_CCT";
            for (int i = 0; i < m_ntaus; ++i)
                oss << "\tP_" << m_tauLabels[i];
            for (int i = 0; i < m_ntaus; ++i)
                oss << "\tZ_" << m_tauLabels[i];
        } else {
            // Legacy: P_CCT  P_tau1..P_tauK  Z_tau1..Z_tauK
            oss << "\tP_CCT";
            for (int i = 1; i <= m_ntaus; ++i)
                oss << "\tP_tau" << i;
            for (int i = 1; i <= m_ntaus; ++i)
                oss << "\tZ_tau" << i;
        }
        return oss.str();
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int /*markerInChunkIdx*/,
        std::vector<double> &result
    ) override {
        result.clear();
        result.reserve(2 * m_ntaus + 1);

        const double gMean = GVec.mean();
        Eigen::VectorXd scores = m_methodShared->residMat.transpose() * GVec;
        for (int i = 0; i < m_ntaus; ++i)
            scores[i] -= gMean * m_methodShared->residSums[i];

        processOneMarker(scores.data(), altFreq, result);
    }

    // ── Batched analysis: B markers at once ────────────────────────────
    void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        std::vector<std::vector<double> > &results
    ) override {
        const int B = static_cast<int>(GBatch.cols());
        results.resize(B);

        Eigen::MatrixXd scoreMatrix;
        scoreMatrix.noalias() = m_methodShared->residMat.transpose() * GBatch;

        const Eigen::VectorXd gMeans = GBatch.colwise().mean();
        scoreMatrix.noalias() -= m_methodShared->residSums * gMeans.transpose();

        for (int b = 0; b < B; ++b)
            processOneMarker(scoreMatrix.col(b).data(), altFreqs[b], results[b]);
    }

    int preferredBatchSize() const override {
        return std::min(std::max(4, 2 * m_ntaus), 16);
    }

    // ── Fused union-level GEMM interface ───────────────────────────────
    bool supportsFusedGemm() const override {
        return true;
    }

    int fusedGemmColumns() const override {
        return m_ntaus;
    }

    void fillUnionResiduals(
        Eigen::Ref<Eigen::MatrixXd> dest,
        const std::vector<uint32_t> &unionToLocal
    ) const override {
        // dest is pre-zeroed, N_union × ntaus.
        // Scatter N_p-level residMat rows into union-level rows.
        const auto &residMat = m_methodShared->residMat;
        const uint32_t nUnion = static_cast<uint32_t>(unionToLocal.size());
        for (uint32_t i = 0; i < nUnion; ++i) {
            const uint32_t li = unionToLocal[i];
            if (li != UINT32_MAX)
                dest.row(i) = residMat.row(li);
        }
    }

    void fillResidualSums(double *dest) const override {
        const auto &rs = m_methodShared->residSums;
        std::copy(rs.data(), rs.data() + m_ntaus, dest);
    }

    void processScoreBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &scores,
        const double *gSums,
        uint32_t nUsed,
        const std::vector<double> &altFreqs,
        std::vector<std::vector<double> > &results
    ) override {
        const int B = static_cast<int>(scores.cols());
        results.resize(B);

        // Center: score[t][b] -= (gSums[b] / nUsed) * residSums[t]
        // Create a temporary to apply centering.
        Eigen::MatrixXd centered = scores;
        const double invN = 1.0 / static_cast<double>(nUsed);
        for (int b = 0; b < B; ++b) {
            const double gMean = gSums[b] * invN;
            for (int t = 0; t < m_ntaus; ++t)
                centered(t, b) -= gMean * m_methodShared->residSums[t];
        }

        for (int b = 0; b < B; ++b)
            processOneMarker(centered.col(b).data(), altFreqs[b], results[b]);
    }

  private:
    int m_ntaus;
    std::vector<SPAGRMClass> m_spagrm_vec;
    std::vector<std::string> m_tauLabels;
    bool m_hasLabels;

    // Shared read-only data: residual matrix and sums for fused dot products.
    struct SharedMethodData {
        Eigen::MatrixXd residMat;   // N × ntaus
        Eigen::VectorXd residSums;  // ntaus
    };

    std::shared_ptr<const SharedMethodData> m_methodShared;

    // ── Internal: SPA + CCT for one marker given centered scores ──────
    // centeredScores points to ntaus doubles (pre-centered by gMean×residSums).
    void processOneMarker(
        const double *centeredScores,
        double altFreq,
        std::vector<double> &result
    ) {
        result.clear();
        result.reserve(2 * m_ntaus + 1);

        std::vector<double> zScores(m_ntaus);
        std::vector<double> pvals(m_ntaus);

        for (int i = 0; i < m_ntaus; ++i) {
            double z;
            double p = m_spagrm_vec[i].getMarkerPvalFromScore(centeredScores[i], altFreq, z);
            zScores[i] = z;
            pvals[i] = p;
        }

        // CCT (Cauchy combination test) p-value
        std::vector<double> valid_p;
        valid_p.reserve(m_ntaus);
        for (int i = 0; i < m_ntaus; ++i)
            if (!std::isnan(pvals[i])) valid_p.push_back(pvals[i]);

        double pCCT = std::numeric_limits<double>::quiet_NaN();
        if (!valid_p.empty()) {
            bool hasZero = false;
            double tStat = 0.0;
            for (double p : valid_p) {
                if (p <= 0.0) {
                    hasZero = true;
                    break;
                }
                double pc = (p >= 1.0) ? 0.999 : p;
                tStat += std::tan((0.5 - pc) * M_PI);
            }
            if (hasZero) {
                pCCT = 0.0;
            } else {
                tStat /= static_cast<double>(valid_p.size());
                pCCT = (tStat > 1e15) ? (1.0 / tStat) / M_PI : 0.5 - std::atan(tStat) / M_PI;
            }
        }

        result.push_back(pCCT);
        for (int i = 0; i < m_ntaus; ++i)
            result.push_back(pvals[i]);
        for (int i = 0; i < m_ntaus; ++i)
            result.push_back(zScores[i]);
    }

    void buildResidMat() {
        const Eigen::Index N = m_spagrm_vec[0].resid().size();
        auto sd = std::make_shared<SharedMethodData>();
        sd->residMat.resize(N, m_ntaus);
        sd->residSums.resize(m_ntaus);
        for (int t = 0; t < m_ntaus; ++t) {
            sd->residMat.col(t) = m_spagrm_vec[t].resid();
            sd->residSums[t] = m_spagrm_vec[t].residSum();
        }
        m_methodShared = std::move(sd);
    }

};

// ══════════════════════════════════════════════════════════════════════
// Outlier detection (IQR-based, per column)
// ══════════════════════════════════════════════════════════════════════

// Returns an N × K boolean matrix (as std::vector<std::vector<bool>>).
// outlierIqrRatio  = multiplier for IQR (default 1.5)
// outlierAbsBound  = absolute clamp for cutoffs (default 0.55)
struct OutlierInfo {
    // per-column: indices of outlier subjects
    std::vector<std::vector<int> > outlierIdx;
    // per-column: boolean mask (size N)
    std::vector<std::vector<bool> > isOutlier;
};

OutlierInfo detectOutliers(
    const Eigen::MatrixXd &ResidMat,
    double outlierIqrRatio,
    double outlierAbsBound
) {
    const Eigen::Index N = ResidMat.rows();
    const Eigen::Index K = ResidMat.cols();

    OutlierInfo info;
    info.outlierIdx.resize(K);
    info.isOutlier.resize(K);

    // Scratch for sorting
    std::vector<double> scratch(N);

    for (Eigen::Index col = 0; col < K; ++col) {
        // Copy column for sorting
        for (Eigen::Index i = 0; i < N; ++i)
            scratch[i] = ResidMat(i, col);

        std::sort(scratch.begin(), scratch.end());

        // Q1, Q3 via linear interpolation (same as R type=7)
        auto quantile = [&](double prob) -> double {
            const double idx = prob * (N - 1);
            const Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
            const Eigen::Index hi = std::min(lo + 1, N - 1);
            const double frac = idx - lo;
            return scratch[lo] * (1.0 - frac) + scratch[hi] * frac;
        };

        const double Q1 = quantile(0.25);
        const double Q3 = quantile(0.75);
        const double IQR = Q3 - Q1;

        double cutLo = Q1 - outlierIqrRatio * IQR;
        double cutHi = Q3 + outlierIqrRatio * IQR;
        cutLo = std::max(cutLo, -outlierAbsBound);
        cutHi = std::min(cutHi, outlierAbsBound);

        info.isOutlier[col].resize(N, false);
        int nOutlier = 0;
        for (Eigen::Index i = 0; i < N; ++i) {
            const double v = ResidMat(i, col);
            if (v < cutLo || v > cutHi) {
                info.isOutlier[col][i] = true;
                info.outlierIdx[col].push_back(static_cast<int>(i));
                ++nOutlier;
            }
        }
    }
    return info;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// Shared pipeline: outlier detection → GRM → SPAGRM → marker engine
// ══════════════════════════════════════════════════════════════════════

// GRMEntry is now declared in spasqr.hpp

// Load GRM entries from disk and convert to flat GRMEntry vector.
std::vector<GRMEntry> loadGrmEntries(
    const std::vector<std::string> &subjOrder,
    const std::vector<std::string> &famIIDs,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile
) {
    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, subjOrder, famIIDs);
    infoMsg("Sparse GRM: %zu entries used", grm.nnz());
    std::vector<GRMEntry> entries;
    entries.reserve(grm.nnz());
    for (const auto &e : grm.entries())
        entries.push_back({e.row, e.col, e.value, (e.row == e.col) ? 1.0 : 2.0});
    return entries;
}

// Build SPAsqrMethod from a pre-computed residual matrix and pre-loaded GRM entries.
std::unique_ptr<MethodBase> buildSPAsqrMethod(
    Eigen::MatrixXd &ResidMat,
    const std::vector<GRMEntry> &grmEntries,
    uint32_t nUsed,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    double minMafCutoff,
    double minMacCutoff,
    std::vector<std::string> tauLabels,
    std::vector<double> *outlierRatiosOut
)
{
    const Eigen::Index K = ResidMat.cols();

    // ── 3. Outlier detection ───────────────────────────────────────────
    OutlierInfo outlierInfo = detectOutliers(ResidMat, outlierIqrRatio, outlierAbsBound);

    // Populate outlier ratios for caller to format
    if (outlierRatiosOut) {
        outlierRatiosOut->resize(K);
        const Eigen::Index N = ResidMat.rows();
        for (Eigen::Index c = 0; c < K; ++c)
            (*outlierRatiosOut)[c] = static_cast<double>(outlierInfo.outlierIdx[c].size()) / N;
    }

    // ── 4. Compute per-column variance terms using pre-loaded GRM ────
    const int ntaus = static_cast<int>(K);
    std::vector<double> R_GRM_R_vec(ntaus, 0.0);
    std::vector<double> R_GRM_R_nonOutlier_vec(ntaus, 0.0);
    std::vector<double> sum_R_nonOutlier_vec(ntaus, 0.0);
    std::vector<Eigen::VectorXd> resid_outliers_lst(ntaus);

    for (int col = 0; col < ntaus; ++col) {
        const auto &isOut = outlierInfo.isOutlier[col];

        // R_GRM_R and R_GRM_R_nonOutlier
        double rgrm_r = 0.0;
        double rgrm_r_no = 0.0;
        for (const auto &e : grmEntries) {
            const double contrib = e.factor * e.value * ResidMat(e.row, col) * ResidMat(e.col, col);
            rgrm_r += contrib;
            if (!isOut[e.row] && !isOut[e.col]) rgrm_r_no += contrib;
        }
        R_GRM_R_vec[col] = rgrm_r;
        R_GRM_R_nonOutlier_vec[col] = rgrm_r_no;

        // sum_R_nonOutlier and resid outlier values
        double sumNO = 0.0;
        std::vector<double> outVals;
        outVals.reserve(outlierInfo.outlierIdx[col].size());
        for (uint32_t i = 0; i < nUsed; ++i) {
            if (!isOut[i])sumNO += ResidMat(i, col);
            else outVals.push_back(ResidMat(i, col));
        }
        sum_R_nonOutlier_vec[col] = sumNO;
        resid_outliers_lst[col] =
            Eigen::Map<Eigen::VectorXd>(outVals.data(), static_cast<Eigen::Index>(outVals.size()));
    }

    // ── 6. Build SPAGRMClass instances (one per tau) ───────────────────
    const std::vector<double> MAF_interval = nsGRMNull::buildMafInterval(minMafCutoff, minMacCutoff, nUsed);
    const double zeta = 0.01;
    const double tol = 1e-6;

    std::vector<SPAGRMClass> spagrm_vec;
    spagrm_vec.reserve(ntaus);

    for (int col = 0; col < ntaus; ++col) {
        // SPAsqr: no families — empty FamilyData
        nsSPAGRM::FamilyData fam;
        fam.resid_unrelated_outliers = std::move(resid_outliers_lst[col]);

        spagrm_vec.emplace_back(
            ResidMat.col(col),                     // resid
            sum_R_nonOutlier_vec[col],
            R_GRM_R_nonOutlier_vec[col],
            0.0,                     // R_GRM_R_TwoSubjOutlier = 0 (no families)
            R_GRM_R_vec[col],
            MAF_interval,
            std::move(fam),
            spaCutoff,
            zeta,
            tol
        );
    }

    // ── 7. Build method ──────────────────────────────────────────────
    std::unique_ptr<SPAsqrMethod> method = std::make_unique<SPAsqrMethod>(
        ntaus,
        std::move(spagrm_vec),
        std::move(tauLabels)
    );

    return method;
}

// ══════════════════════════════════════════════════════════════════════
// Re-index GRM entries from union-dense space to pheno-dense space
// ══════════════════════════════════════════════════════════════════════

static std::vector<GRMEntry> reindexGrm(
    const std::vector<GRMEntry> &unionGrm,
    const std::vector<uint32_t> &unionToLocal, // union index → pheno index (UINT32_MAX = absent)
    uint32_t nUnion
) {
    std::vector<GRMEntry> out;
    out.reserve(unionGrm.size());
    for (const auto &e : unionGrm) {
        if (e.row >= nUnion || e.col >= nUnion) continue;
        uint32_t lr = unionToLocal[e.row];
        uint32_t lc = unionToLocal[e.col];
        if (lr == UINT32_MAX || lc == UINT32_MAX) continue;
        out.push_back({lr, lc, e.value, e.factor});
    }
    return out;
}

// ══════════════════════════════════════════════════════════════════════
// runSPAsqr — multi-phenotype entry point with parallel conquer fits
// ══════════════════════════════════════════════════════════════════════

void runSPAsqr(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<double> &taus,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    double spasqrTol,
    double spasqrH,
    double spasqrHScale,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const int K = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());

    // ── 1. Load phenotype/covariate data (union mask) ───────────────
    // Union = subjects with genotype ∩ GRM ∩ keep/remove.  Per-phenotype
    // NA filtering is deferred — each phenotype uses its own non-missing
    // subset of the union.
    infoMsg("SPAsqr: Loading phenotype and covariate data (%d phenotypes, %d taus)", K, ntaus);
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile, phenoNames);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();

    const uint32_t nUnion = sd.nUsed();
    const Eigen::Index N = static_cast<Eigen::Index>(nUnion);

    // ── 2. Extract union-space covariates ───────────────────────────
    Eigen::MatrixXd unionX = covarNames.empty()
        ? (sd.hasCovar() ? Eigen::MatrixXd(sd.covar()) : Eigen::MatrixXd(N, 0))
        : sd.getColumns(covarNames);
    const int nCov = static_cast<int>(unionX.cols());

    // ── 3. Per-phenotype: build non-missing mask, extract Y/X ───────
    struct PhenoWork {
        std::vector<uint32_t> unionToLocal; // size nUnion; UINT32_MAX = absent
        uint32_t nk;                        // non-missing count
        Eigen::VectorXd Y;                  // nk
        Eigen::MatrixXd X;                  // nk × nCov
        double h;                           // bandwidth
        Eigen::MatrixXd ResidMat;           // nk × ntaus (filled by conquer)
    };

    std::vector<PhenoWork> pw(K);

    for (int k = 0; k < K; ++k) {
        Eigen::VectorXd fullY = sd.getColumn(phenoNames[k]);
        pw[k].unionToLocal.resize(nUnion, UINT32_MAX);
        uint32_t localIdx = 0;
        for (uint32_t i = 0; i < nUnion; ++i) {
            if (!std::isnan(fullY[i])) {
                pw[k].unionToLocal[i] = localIdx++;
            }
        }
        pw[k].nk = localIdx;
        if (pw[k].nk == 0)
            throw std::runtime_error("SPAsqr: phenotype '" + phenoNames[k] + "' has no non-missing subjects");

        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        pw[k].Y.resize(Nk);
        pw[k].X.resize(Nk, nCov);
        for (uint32_t i = 0; i < nUnion; ++i) {
            uint32_t li = pw[k].unionToLocal[i];
            if (li == UINT32_MAX) continue;
            pw[k].Y[li] = fullY[i];
            if (nCov > 0) pw[k].X.row(li) = unionX.row(i);
        }

        // Per-phenotype bandwidth
        if (spasqrH >= 0.0) {
            pw[k].h = spasqrH;
        } else {
            std::vector<double> ysorted(Nk);
            Eigen::VectorXd::Map(ysorted.data(), Nk) = pw[k].Y;
            std::sort(ysorted.begin(), ysorted.end());
            auto quantile = [&](double prob) -> double {
                double idx = prob * (Nk - 1);
                Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
                Eigen::Index hi = std::min(lo + 1, Nk - 1);
                double frac = idx - lo;
                return ysorted[lo] * (1.0 - frac) + ysorted[hi] * frac;
            };
            double iqr = quantile(0.75) - quantile(0.25);
            double scale = (spasqrHScale >= 0.0) ? spasqrHScale : 3.0;
            pw[k].h = iqr / scale;
            if (pw[k].h <= 0.0)
                pw[k].h = std::max(std::pow((std::log(Nk) + nCov) / static_cast<double>(Nk), 0.4), 0.05);
        }

        pw[k].ResidMat.resize(Nk, ntaus);
    }

    // Log N/bandwidth table
    {
        infoMsg("Sample size and smooth bandwidth per phenotype:");
        size_t nameW = 4;
        for (int k = 0; k < K; ++k)
            nameW = std::max(nameW, phenoNames[k].size());
        nameW += 2;
        char row[256];
        std::snprintf(row, sizeof(row), "    %-*s %10s %12s\n", static_cast<int>(nameW), "", "N", "bandwidth");
        fprintf(stderr, "%s", row);
        for (int k = 0; k < K; ++k) {
            std::snprintf(row, sizeof(row), "    %-*s %10u %12.6f\n",
                          static_cast<int>(nameW), phenoNames[k].c_str(), pw[k].nk, pw[k].h);
            fprintf(stderr, "%s", row);
        }
    }

    // ── 4. Parallel conquer fits: K × ntaus ─────────────────────────
    const int totalFits = K * ntaus;
    const int nWorkers = std::min(nthreads, totalFits);
    infoMsg("SPAsqr: Running %d conquer fits with %d threads", totalFits, nWorkers);

    std::atomic<int> nextFit{0};
    std::vector<std::string> fitErrors(totalFits);

    auto fitWorker = [&]() {
        for (;;) {
            int idx = nextFit.fetch_add(1, std::memory_order_relaxed);
            if (idx >= totalFits) break;
            int k = idx / ntaus;
            int t = idx % ntaus;
            const double h1 = 1.0 / pw[k].h;

            try {
                Eigen::VectorXd resid;
                Eigen::VectorXd beta = conquer::smqrGauss(pw[k].X, pw[k].Y, taus[t], pw[k].h, &resid, spasqrTol);
                infoMsg("[%s] tau=%.4f intercept=%.6f", phenoNames[k].c_str(), taus[t], beta(0));

                const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
                for (Eigen::Index i = 0; i < Nk; ++i)
                    pw[k].ResidMat(i, t) = taus[t] - math::pnorm(-resid(i) * h1);
            } catch (const std::exception &ex) {
                fitErrors[idx] = ex.what();
            }
        }
    };

    {
        std::vector<std::thread> threads;
        threads.reserve(nWorkers - 1);
        for (int t = 0; t < nWorkers - 1; ++t)
            threads.emplace_back(fitWorker);
        fitWorker();
        for (auto &th : threads)
            th.join();
    }

    for (int idx = 0; idx < totalFits; ++idx) {
        if (!fitErrors[idx].empty()) {
            int k = idx / ntaus;
            int t = idx % ntaus;
            throw std::runtime_error("SPAsqr: conquer failed for phenotype '" +
                                     phenoNames[k] + "' tau=" + std::to_string(taus[t]) + ": " + fitErrors[idx]);
        }
    }

    // ── 5. Build tau labels ────────────────────────────────────────────
    std::vector<std::string> tauLabels;
    tauLabels.reserve(ntaus);
    for (double tau : taus) {
        std::ostringstream oss;
        oss << "tau" << tau;
        tauLabels.push_back(oss.str());
    }

    // ── 6. Load genotype data and GRM once (shared, union space) ────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    std::vector<GRMEntry> unionGrm = loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), spgrmGrabFile, spgrmGctaFile);

    // ── 7. Build per-phenotype SPAsqrMethod + PhenoTask ─────────────
    std::vector<PhenoTask> tasks(K);
    std::vector<std::vector<double> > allOutlierRatios(K);

    for (int k = 0; k < K; ++k) {
        infoMsg("[%s] Building SPAsqr method (%d taus, %u subjects)",
                phenoNames[k].c_str(), ntaus, pw[k].nk);

        // Re-index GRM to pheno-dense space
        auto phenoGrm = reindexGrm(unionGrm, pw[k].unionToLocal, nUnion);

        auto method = buildSPAsqrMethod(
            pw[k].ResidMat,
            phenoGrm,
            pw[k].nk,
            spaCutoff,
            outlierIqrRatio,
            outlierAbsBound,
            minMafCutoff,
            minMacCutoff,
            tauLabels,
            &allOutlierRatios[k]
        );

        tasks[k].phenoName = phenoNames[k];
        tasks[k].method = std::move(method);
        tasks[k].unionToLocal = pw[k].unionToLocal;
        tasks[k].nUsed = pw[k].nk;
    }

    // Print outlier ratio table line-by-line
    {
        size_t nameW = 4;
        for (int k = 0; k < K; ++k)
            nameW = std::max(nameW, phenoNames[k].size());
        nameW += 2;

        infoMsg("Outlier ratios (IQR=%.2f, bound=%.2f):", outlierIqrRatio, outlierAbsBound);

        std::ostringstream hdr;
        hdr << "    " << std::setw(static_cast<int>(nameW)) << std::left << "";
        for (const auto &tl : tauLabels)
            hdr << std::setw(10) << std::right << tl;
        fprintf(stderr, "%s\n", hdr.str().c_str());

        for (int k = 0; k < K; ++k) {
            std::ostringstream row;
            row << "    " << std::setw(static_cast<int>(nameW)) << std::left << phenoNames[k];
            for (double r : allOutlierRatios[k])
                row << std::setw(10) << std::right << std::fixed << std::setprecision(4) << r;
            fprintf(stderr, "%s\n", row.str().c_str());
        }
    }

    // Free per-phenotype work data
    pw.clear();

    // ── 8. Run multi-phenotype engine ───────────────────────────────
    infoMsg("SPAsqr: Starting multi-phenotype association (%d phenotypes, %d taus, %d threads)", K, ntaus, nthreads);
    multiPhenoEngine(
        *genoData, tasks, outPrefix, "SPAsqr", compression, compressionLevel,
        nthreads, missingCutoff, minMafCutoff, minMacCutoff, hweCutoff
    );
}

// ══════════════════════════════════════════════════════════════════════
// runSPAsqrLoco — LOCO entry point
//
// Matches the serial workflow (2.serial_loco.sh):
//   - conquer fits are chromosome-independent (same Y, same X)
//   - per-phenotype bandwidth h
//   - SPAsqrMethod built once per phenotype, cloned per chromosome
//   - locoEngine iterates chromosomes, testing each chromosome's markers
// ══════════════════════════════════════════════════════════════════════

void runSPAsqrLoco(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<double> &taus,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &predListFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    double spasqrTol,
    double spasqrH,
    double spasqrHScale,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const int K = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());

    // ── 1. Load phenotype/covariate data (union mask) ───────────────
    // Union = subjects with genotype ∩ GRM ∩ keep/remove.  Per-phenotype
    // NA filtering is deferred — each phenotype uses its own non-missing
    // subset of the union.
    infoMsg("SPAsqr-LOCO: Loading phenotype and covariate data (%d phenotypes, %d taus)", K, ntaus);
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile, phenoNames);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();

    const uint32_t nUnion = sd.nUsed();
    const Eigen::Index N = static_cast<Eigen::Index>(nUnion);

    // ── 2. Extract union-space covariates ───────────────────────────
    Eigen::MatrixXd unionX = covarNames.empty()
        ? (sd.hasCovar() ? Eigen::MatrixXd(sd.covar()) : Eigen::MatrixXd(N, 0))
        : sd.getColumns(covarNames);
    const int nCov = static_cast<int>(unionX.cols());

    // ── 3. Per-phenotype: build non-missing mask, extract Y/baseX ──
    struct PhenoWork {
        std::vector<uint32_t> unionToLocal; // size nUnion; UINT32_MAX = absent
        uint32_t nk;                        // non-missing count
        Eigen::VectorXd Y;                  // nk
        Eigen::MatrixXd baseX;              // nk × nCov
        double h;                           // bandwidth
    };

    std::vector<PhenoWork> pw(K);

    for (int k = 0; k < K; ++k) {
        Eigen::VectorXd fullY = sd.getColumn(phenoNames[k]);
        pw[k].unionToLocal.resize(nUnion, UINT32_MAX);
        uint32_t localIdx = 0;
        for (uint32_t i = 0; i < nUnion; ++i) {
            if (!std::isnan(fullY[i])) {
                pw[k].unionToLocal[i] = localIdx++;
            }
        }
        pw[k].nk = localIdx;
        if (pw[k].nk == 0)
            throw std::runtime_error("SPAsqr-LOCO: phenotype '" + phenoNames[k] + "' has no non-missing subjects");

        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        pw[k].Y.resize(Nk);
        pw[k].baseX.resize(Nk, nCov);
        for (uint32_t i = 0; i < nUnion; ++i) {
            uint32_t li = pw[k].unionToLocal[i];
            if (li == UINT32_MAX) continue;
            pw[k].Y[li] = fullY[i];
            if (nCov > 0) pw[k].baseX.row(li) = unionX.row(i);
        }

        // Per-phenotype bandwidth
        if (spasqrH >= 0.0) {
            pw[k].h = spasqrH;
        } else {
            std::vector<double> ysorted(Nk);
            Eigen::VectorXd::Map(ysorted.data(), Nk) = pw[k].Y;
            std::sort(ysorted.begin(), ysorted.end());
            auto quantile = [&](double prob) -> double {
                double idx = prob * (Nk - 1);
                Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
                Eigen::Index hi = std::min(lo + 1, Nk - 1);
                double frac = idx - lo;
                return ysorted[lo] * (1.0 - frac) + ysorted[hi] * frac;
            };
            double iqr = quantile(0.75) - quantile(0.25);
            double scale = (spasqrHScale >= 0.0) ? spasqrHScale : 3.0;
            pw[k].h = iqr / scale;
            if (pw[k].h <= 0.0)
                pw[k].h = std::max(std::pow((std::log(Nk) + nCov) / static_cast<double>(Nk), 0.4), 0.05);
        }
    }

    // Log N/bandwidth table
    {
        infoMsg("Sample size and smooth bandwidth per phenotype:");
        size_t nameW = 4;
        for (int k = 0; k < K; ++k)
            nameW = std::max(nameW, phenoNames[k].size());
        nameW += 2;
        char row[256];
        std::snprintf(row, sizeof(row), "    %-*s %10s %12s\n", static_cast<int>(nameW), "", "N", "bandwidth");
        fprintf(stderr, "%s", row);
        for (int k = 0; k < K; ++k) {
            std::snprintf(row, sizeof(row), "    %-*s %10u %12.6f\n",
                          static_cast<int>(nameW), phenoNames[k].c_str(), pw[k].nk, pw[k].h);
            fprintf(stderr, "%s", row);
        }
    }

    // ── 4. Load LOCO predictions (union space) ─────────────────────
    LocoData loco = LocoData::load(predListFile, phenoNames, sd.usedIIDs(), sd.famIIDs());
    auto locoChroms = loco.availableChromosomes();
    infoMsg("LOCO: %zu chromosomes available across all phenotypes", locoChroms.size());

    // ── 5. Build tau labels ─────────────────────────────────────────
    std::vector<std::string> tauLabels;
    tauLabels.reserve(ntaus);
    for (double tau : taus) {
        std::ostringstream oss;
        oss << "tau" << tau;
        tauLabels.push_back(oss.str());
    }

    // ── 6. Load genotype data and GRM once (union space) ────────────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    std::vector<GRMEntry> unionGrm = loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), spgrmGrabFile, spgrmGctaFile);

    // Pre-compute per-phenotype re-indexed GRM (shared across chromosomes)
    std::vector<std::vector<GRMEntry> > phenoGrms(K);
    for (int k = 0; k < K; ++k)
        phenoGrms[k] = reindexGrm(unionGrm, pw[k].unionToLocal, nUnion);

    // ── 7. Build LocoTaskBuilder callback ──────────────────────────
    // For each chromosome, augment base covariates with the LOCO column
    // (mapped to pheno-dense space), run K × ntaus conquer fits, and
    // build K SPAsqrMethods.
    auto buildTasks = [&](const std::string &chr, std::vector<PhenoTask> &tasks) {
        tasks.resize(K);

        // Augment base covariates with LOCO column per phenotype (pheno-dense)
        std::vector<Eigen::MatrixXd> X_augs(K);
        for (int k = 0; k < K; ++k) {
            const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
            X_augs[k].resize(Nk, nCov + 1);
            X_augs[k].leftCols(nCov) = pw[k].baseX;
            // Extract LOCO scores from union to pheno-dense space
            const auto &locoVec = loco.scores.at(phenoNames[k]).at(chr);
            for (uint32_t i = 0; i < nUnion; ++i) {
                uint32_t li = pw[k].unionToLocal[i];
                if (li != UINT32_MAX)
                    X_augs[k](li, nCov) = locoVec[i];
            }
        }

        // Parallel conquer fits: K × ntaus
        std::vector<Eigen::MatrixXd> ResidMats(K);
        for (int k = 0; k < K; ++k)
            ResidMats[k].resize(static_cast<Eigen::Index>(pw[k].nk), ntaus);

        const int totalFits = K * ntaus;
        const int nWorkers = std::min(nthreads, totalFits);
        infoMsg("SPAsqr-LOCO chr%s: Running %d conquer fits with %d threads",
                chr.c_str(), totalFits, nWorkers);

        std::atomic<int> nextFit{0};
        std::vector<std::string> fitErrors(totalFits);

        auto fitWorker = [&]() {
            for (;;) {
                int idx = nextFit.fetch_add(1, std::memory_order_relaxed);
                if (idx >= totalFits) break;
                int k = idx / ntaus;
                int t = idx % ntaus;
                const double h1 = 1.0 / pw[k].h;

                try {
                    Eigen::VectorXd resid;
                    conquer::smqrGauss(X_augs[k], pw[k].Y, taus[t], pw[k].h, &resid, spasqrTol);

                    const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
                    for (Eigen::Index i = 0; i < Nk; ++i)
                        ResidMats[k](i, t) = taus[t] - math::pnorm(-resid(i) * h1);
                } catch (const std::exception &ex) {
                    fitErrors[idx] = ex.what();
                }
            }
        };

        {
            std::vector<std::thread> threads;
            threads.reserve(nWorkers - 1);
            for (int t = 0; t < nWorkers - 1; ++t)
                threads.emplace_back(fitWorker);
            fitWorker();
            for (auto &th : threads)
                th.join();
        }

        for (int idx = 0; idx < totalFits; ++idx) {
            if (!fitErrors[idx].empty()) {
                int k = idx / ntaus;
                int t = idx % ntaus;
                throw std::runtime_error("SPAsqr-LOCO chr" + chr + ": conquer failed for phenotype '" +
                                         phenoNames[k] + "' tau=" + std::to_string(taus[t]) + ": " + fitErrors[idx]);
            }
        }

        // Build SPAsqrMethod for each phenotype (pheno-dense space)
        std::vector<std::vector<double> > allOutlierRatios(K);
        for (int k = 0; k < K; ++k) {
            auto method = buildSPAsqrMethod(
                ResidMats[k],
                phenoGrms[k],
                pw[k].nk,
                spaCutoff,
                outlierIqrRatio,
                outlierAbsBound,
                minMafCutoff,
                minMacCutoff,
                tauLabels,
                &allOutlierRatios[k]
            );

            tasks[k].phenoName = phenoNames[k];
            tasks[k].method = std::move(method);
            tasks[k].unionToLocal = pw[k].unionToLocal;
            tasks[k].nUsed = pw[k].nk;
        }

        // Print outlier ratio table for this chromosome line-by-line
        {
            size_t nameW = 4;
            for (int k = 0; k < K; ++k)
                nameW = std::max(nameW, phenoNames[k].size());
            nameW += 2;

            infoMsg("chr%s outlier ratios (IQR=%.2f, bound=%.2f):",
                    chr.c_str(), outlierIqrRatio, outlierAbsBound);

            std::ostringstream hdr;
            hdr << "    " << std::setw(static_cast<int>(nameW)) << std::left << "";
            for (const auto &tl : tauLabels)
                hdr << std::setw(10) << std::right << tl;
            fprintf(stderr, "%s\n", hdr.str().c_str());

            for (int k = 0; k < K; ++k) {
                std::ostringstream row;
                row << "    " << std::setw(static_cast<int>(nameW)) << std::left << phenoNames[k];
                for (double r : allOutlierRatios[k])
                    row << std::setw(10) << std::right << std::fixed << std::setprecision(4) << r;
                fprintf(stderr, "%s\n", row.str().c_str());
            }
        }
    };

    // ── 8. Run LOCO engine ─────────────────────────────────────────
    const int nChroms = static_cast<int>(locoChroms.size());
    infoMsg("SPAsqr-LOCO: Starting LOCO association (%d phenotypes, %d taus, %d chroms, %d threads)",
            K, ntaus, nChroms, nthreads);
    locoEngine(
        *genoData,
        locoChroms,
        phenoNames,
        buildTasks,
        outPrefix,
        "SPAsqr",
        compression,
        compressionLevel,
        nthreads,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        hweCutoff
    );
}
