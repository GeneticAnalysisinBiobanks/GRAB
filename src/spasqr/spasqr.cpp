// spasqr.cpp — SPAsqr: SPA-squared multi-tau marker association (pure C++17 / Eigen / Boost)

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "spasqr/spasqr.hpp"
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
//   labels like "Tau0.1", "Tau0.3", ...
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
    }

    std::unique_ptr<MethodBase> clone() const override {
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

        std::vector<double> zScores(m_ntaus);
        std::vector<double> pvals(m_ntaus);

        for (int i = 0; i < m_ntaus; ++i) {
            double z;
            double p = m_spagrm_vec[i].getMarkerPval(GVec, altFreq, z);
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

        if (m_hasLabels) {
            // P_CCT first, then P_tau{val}..., then Z_tau{val}...
            result.push_back(pCCT);
            for (int i = 0; i < m_ntaus; ++i)
                result.push_back(pvals[i]);
            for (int i = 0; i < m_ntaus; ++i)
                result.push_back(zScores[i]);
        } else {
            // Legacy: P_CCT  P_tau1..P_tauK  Z_tau1..Z_tauK
            result.push_back(pCCT);
            for (int i = 0; i < m_ntaus; ++i)
                result.push_back(pvals[i]);
            for (int i = 0; i < m_ntaus; ++i)
                result.push_back(zScores[i]);
        }
    }

  private:
    int m_ntaus;
    std::vector<SPAGRMClass> m_spagrm_vec;
    std::vector<std::string> m_tauLabels;
    bool m_hasLabels;
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

        infoMsg(
            "  Column %d: outlier ratio = %.4f (%d / %lld)",
            static_cast<int>(col + 1),
            static_cast<double>(nOutlier) / N,
            nOutlier,
            static_cast<long long>(N)
        );
    }
    return info;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// Shared pipeline: outlier detection → GRM → SPAGRM → marker engine
// ══════════════════════════════════════════════════════════════════════

struct GRMEntry {
    uint32_t row, col;
    double value;
    double factor; // 1 for diagonal, 2 for off-diagonal
};

// Load GRM entries from disk and convert to flat GRMEntry vector.
static std::vector<GRMEntry> loadGrmEntries(
    const std::vector<std::string> &subjOrder,
    const std::vector<std::string> &famIIDs,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile
) {
    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, subjOrder, famIIDs);
    infoMsg("Sparse GRM: %zu entries after filtering", grm.nnz());
    std::vector<GRMEntry> entries;
    entries.reserve(grm.nnz());
    for (const auto &e : grm.entries())
        entries.push_back({e.row, e.col, e.value, (e.row == e.col) ? 1.0 : 2.0});
    return entries;
}

// Build SPAsqrMethod from a pre-computed residual matrix and pre-loaded GRM entries.
static std::unique_ptr<SPAsqrMethod> buildSPAsqrMethod(
    Eigen::MatrixXd &ResidMat,
    const std::vector<GRMEntry> &grmEntries,
    uint32_t nUsed,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    double minMafCutoff,
    double minMacCutoff,
    std::vector<std::string> tauLabels
)
{
    const Eigen::Index K = ResidMat.cols();

    // ── 3. Outlier detection ───────────────────────────────────────────
    infoMsg("Detecting outliers (IQR ratio=%.2f, abs bound=%.2f)", outlierIqrRatio, outlierAbsBound);
    OutlierInfo outlierInfo = detectOutliers(ResidMat, outlierIqrRatio, outlierAbsBound);

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

    infoMsg("SPAsqr method built (%d taus)", ntaus);
    return method;
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
    infoMsg("SPAsqr: Loading phenotype and covariate data (%d phenotypes, %d taus)", K, ntaus);
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();
    // Drop subjects with NA in ANY phenotype column (complete cases across traits).
    sd.dropNaInColumns(phenoNames);

    const Eigen::Index N = static_cast<Eigen::Index>(sd.nUsed());
    infoMsg("Subjects after intersection: %lld", static_cast<long long>(N));

    // ── 2. Extract covariates (shared across all phenotypes) ────────
    Eigen::MatrixXd X = covarNames.empty() ? Eigen::MatrixXd(N, 0) : sd.getColumns(covarNames);
    const int p = static_cast<int>(X.cols());

    // ── 3. Compute bandwidth h (shared) ─────────────────────────────
    // Use the first phenotype's IQR for bandwidth; all phenotypes share h.
    Eigen::VectorXd Y0 = sd.getColumn(phenoNames[0]);
    double h;
    if (spasqrH >= 0.0) {
        h = spasqrH;
    } else {
        std::vector<double> ysorted(N);
        Eigen::VectorXd::Map(ysorted.data(), N) = Y0;
        std::sort(ysorted.begin(), ysorted.end());
        auto quantile = [&](double prob) -> double {
            double idx = prob * (N - 1);
            Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
            Eigen::Index hi = std::min(lo + 1, N - 1);
            double frac = idx - lo;
            return ysorted[lo] * (1.0 - frac) + ysorted[hi] * frac;
        };
        double iqr = quantile(0.75) - quantile(0.25);
        double scale = (spasqrHScale >= 0.0) ? spasqrHScale : 3.0;
        h = iqr / scale;
        if (h <= 0.0) h = std::max(std::pow((std::log(N) + p) / static_cast<double>(N), 0.4), 0.05);
    }
    infoMsg("Bandwidth h = %.6f", h);

    // ── 4. Parallel conquer fits: K × ntaus independent regressions ─
    // Build flat work items: (pheno_idx, tau_idx).
    // Each worker runs conquer::smqrGauss and writes one column of the
    // per-phenotype ResidMat.
    std::vector<Eigen::VectorXd> phenoY(K);
    for (int k = 0; k < K; ++k)
        phenoY[k] = sd.getColumn(phenoNames[k]);

    // Per-phenotype residual matrices
    std::vector<Eigen::MatrixXd> ResidMats(K);
    for (int k = 0; k < K; ++k)
        ResidMats[k].resize(N, ntaus);

    const int totalFits = K * ntaus;
    const int nWorkers = std::min(nthreads, totalFits);
    infoMsg("SPAsqr: Running %d conquer fits with %d threads", totalFits, nWorkers);

    std::atomic<int> nextFit{0};
    std::vector<std::string> fitErrors(totalFits);

    auto fitWorker = [&]() {
        const double h1 = 1.0 / h;
        for (;;) {
            int idx = nextFit.fetch_add(1, std::memory_order_relaxed);
            if (idx >= totalFits) break;
            int k = idx / ntaus;
            int t = idx % ntaus;

            try {
                infoMsg("[%s] conquer tau=%.4f ...", phenoNames[k].c_str(), taus[t]);
                Eigen::VectorXd resid;
                Eigen::VectorXd beta = conquer::smqrGauss(X, phenoY[k], taus[t], h, &resid, spasqrTol);
                infoMsg("[%s] tau=%.4f intercept=%.6f", phenoNames[k].c_str(), taus[t], beta(0));

                for (Eigen::Index i = 0; i < N; ++i) ResidMats[k](i, t) = taus[t] - math::pnorm(-resid(i) * h1);
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
        fitWorker(); // caller participates
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
        oss << "Tau" << tau;
        tauLabels.push_back(oss.str());
    }

    // ── 6. Load genotype data and GRM once (shared) ────────────────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    std::vector<GRMEntry> grmEntries = loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), spgrmGrabFile, spgrmGctaFile);

    // ── 7. Build per-phenotype SPAsqrMethod + PhenoTask ─────────────
    std::vector<PhenoTask> tasks(K);

    for (int k = 0; k < K; ++k) {
        infoMsg("[%s] Building SPAsqr method (%d taus)", phenoNames[k].c_str(), ntaus);
        auto method = buildSPAsqrMethod(
            ResidMats[k],
            grmEntries,
            genoData->nSubjUsed(),
            spaCutoff,
            outlierIqrRatio,
            outlierAbsBound,
            minMafCutoff,
            minMacCutoff,
            tauLabels // copy; last iteration could move but let's keep it simple
        );

        tasks[k].phenoName = phenoNames[k];
        tasks[k].method = std::move(method);
        // Identity mapping: all subjects are in the union (complete cases)
        tasks[k].unionToLocal.resize(genoData->nSubjUsed());
        std::iota(tasks[k].unionToLocal.begin(), tasks[k].unionToLocal.end(), 0u);
        tasks[k].nUsed = genoData->nSubjUsed();
    }

    // ── 8. Run multi-phenotype engine ───────────────────────────────
    infoMsg("SPAsqr: Starting multi-phenotype association (%d phenotypes, %d taus, %d threads)", K, ntaus, nthreads);
    multiPhenoEngine(
        *genoData,
        tasks,
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
