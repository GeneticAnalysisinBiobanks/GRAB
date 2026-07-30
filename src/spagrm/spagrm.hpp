// spagrm.hpp — Saddlepoint approximation with GRM (pure C++17 / Eigen / Boost)
//
// Translated from mtSPAGRM.h / mtSPAGRM.cpp (RcppArmadillo) to Eigen.
// Provides the nsSPAGRM namespace and the SPAGRMClass marker-level evaluator
// used by SPAGRM and SAGELD.
//
// spa_unify Stage 4.  The root finder is the shared, dynamically bracketed
// safeguarded Newton in util/spa.hpp and the tail is spa::bnTail /
// spa::bnTailLog; SPAGRM's own three-term-class CGF is tier 3 and lives in
// spagrm_cgf.hpp.  The private Family-B Newton iteration (nsSPAGRM::fastGetRoot),
// the private copy of the Barndorff-Nielsen root (nsSPAGRM::getProbSpa), the
// MGF0/MGF1/MGF2/temp workspace (nsSPAGRM::mgf, MgfWorkspace) and the dead
// SPAGRMClass::getMarkerPval are all gone.  See
// dev-notes/methods/spa_unify/02_design.md.
#pragma once

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "engine/marker.hpp"          // MethodBase
#include "geno_factory/geno_data.hpp" // GenoSpec
#include "spagrm/spagrm_cgf.hpp"      // tier-3 CGF + ThreeSubjTable
#include "util/math_helper.hpp"       // math::zFromPval
#include "util/spa.hpp"               // spa::TwoSided, spa::Status

namespace nsSPAGRM {

// Per-marker updated data for three-or-more-subject families.  The type now
// belongs to the CGF that consumes it; the name is retained.
using UpdatedThreeSubj = spagrm_cgf::ThreeSubjTable;

// Family data for one tau (or one residual column).
struct FamilyData {
    Eigen::VectorXd resid_unrelated_outliers;
    std::vector<std::array<double, 2> > twoSubj_resid;
    std::vector<std::vector<double> > twoSubj_rho;
    std::vector<std::vector<double> > threeSubj_standS;
    std::vector<Eigen::MatrixXd> threeSubj_CLT;
};

// Two-sided saddlepoint p-value for one marker.
//
// `absScore` is |Score_adj|, the variance-ratio-rescaled score; both tails are
// solved (upper at +absScore, lower at -absScore) with spa::solveSaddlepoint
// over the tier-3 CGF, evaluated with spa::bnTail / spa::bnTailLog, and
// assembled by spa::combineTails.  `initZeta` is the upper-tail initial
// abscissa; the lower tail starts at -initZeta (see the convention note in
// spagrm.cpp).  `tol` is the relative residual tolerance and reaches BOTH
// tails, which is the D7 repair.
spa::TwoSided twoSidedSpa(
    const spagrm_cgf::Context &cgf,
    double absScore,
    double initZeta,
    double tol
);

} // namespace nsSPAGRM

// ══════════════════════════════════════════════════════════════════════
// SPAGRMClass — per-column (per-tau) marker evaluator
// ══════════════════════════════════════════════════════════════════════

class SPAGRMClass {
  public:
    SPAGRMClass(
        Eigen::VectorXd resid,
        double sum_R_nonOutlier,
        double R_GRM_R_nonOutlier,
        double R_GRM_R_TwoSubjOutlier,
        double R_GRM_R,
        std::vector<double> MAF_interval,
        nsSPAGRM::FamilyData fam,
        double SPA_Cutoff,
        double tol
    );

    // Deep copy (including scratch state) for per-thread isolation.
    SPAGRMClass(const SPAGRMClass &o);
    SPAGRMClass &operator=(const SPAGRMClass &) = delete;

    // Caller pre-computes Score = GVec.dot(m_resid) - gMean * m_resid_sum, via
    // a fused multi-phenotype matrix multiply where one is available.
    //
    // `outScoreVar`, when non-null, receives the nominal variance
    // Var(S) = 2·MAF·(1−MAF)·Rᵀ Φ R that drives the normal-approximation
    // z-score; callers reporting BETA / SE consume it to form
    // BETA = S / Var(S) and SE = 1 / sqrt(Var(S)).  Set to zero for
    // monomorphic or degenerate markers (Var(S) ≤ 0).
    //
    // Returns P, −log10(P) and the spa::Status of the worse of the two tails.
    // P is NaN on every status other than Converged and NormalBranch, so a
    // saddlepoint failure reports NA plus a named reason rather than an
    // ordinary-looking finite number (spa_unify L2, D5, D6).
    spa::TwoSided getMarkerPvalFromScore(
        double Score,
        double altFreq,
        double &zScore,
        double *outScoreVar = nullptr
    );

    const Eigen::VectorXd &resid() const {
        return m_shared->resid;
    }

    double residSum() const {
        return m_shared->resid_sum;
    }

    // Read-only data shared across clones via shared_ptr.
    struct SharedData {
        Eigen::VectorXd resid;
        Eigen::VectorXd resid_unrelated_outliers;
        double sum_unrelated_outliers2;
        double sum_R_nonOutlier;
        double R_GRM_R_nonOutlier;
        double R_GRM_R_TwoSubjOutlier;
        double R_GRM_R;
        double resid_sum;
        std::vector<double> MAF_interval;
        std::vector<std::array<double, 2> > TwoSubj_resid_list;
        std::vector<std::vector<double> > TwoSubj_rho_list;
        std::vector<std::vector<double> > ThreeSubj_standS_list;
        std::vector<Eigen::MatrixXd> ThreeSubj_CLT_list;
        double SPA_Cutoff;
        double tol;
    };

    const std::shared_ptr<const SharedData> &sharedData() const {
        return m_shared;
    }

  private:
    std::shared_ptr<const SharedData> m_shared;

    // Per-thread mutable scratch (rebuilt on copy).  arr_prob is refreshed per
    // marker; m_cgfScratch holds one class-3 family's tilted weights so the
    // mean-centred K'' needs no second exponential pass.  The MGF0 / MGF1 /
    // MGF2 / temp vectors that used to live here — 8 × nOutlier doubles per
    // clone — are gone with the fused CGF (spa_unify P6).
    std::vector<nsSPAGRM::UpdatedThreeSubj> m_threeSubj_scratch;
    std::vector<double> m_cgfScratch;

    void rebuildScratch();

};

// ══════════════════════════════════════════════════════════════════════
// SPAGRMMethod — MethodBase adapter wrapping a single SPAGRMClass
// ══════════════════════════════════════════════════════════════════════

class SPAGRMMethod : public MethodBase {
  public:
    explicit SPAGRMMethod(SPAGRMClass spagrm)
        : m_spagrm(std::move(spagrm))
    {
    }

    std::unique_ptr<MethodBase> clone() const override {
        return std::make_unique<SPAGRMMethod>(*this);
    }

    int resultSize() const override {
        return 7;
    }

// P, LOG10P, Z, Z_Norm, BETA, SE, SPA_STATUS.
//
// LOG10P is −log10(P), computed through spa::bnTailLog / spa::combineTailsLog
// so that it stays meaningful past the point where the linear-scale tail
// underflows (Φ(−38.5) flushes to zero, i.e. p ≈ 1e-316).
//
// SPA_STATUS carries the spa::Status of the saddlepoint that produced P.  It is
// emitted as the integer enumerator rather than the token spelled by
// spa::statusName because MethodBase hands the engine a std::vector<double> and
// every result cell is formatted by numToChars; a string column would require a
// new hook in the MethodBase contract, which dev-notes/methods/spa_unify/
// 02_design.md places out of scope for the per-method migration stages.  The
// encoding is static_cast<uint8_t>(spa::Status), the same one Stage 3
// established for SPACox:
//
//     0 OK (converged)     3 GUARD_CURV     6 NORMAL (|Z| ≤ --spa-z-threshold,
//     1 MAXITER            4 GUARD_W          saddlepoint never attempted)
//     2 GUARD_TEMP         5 NONFINITE
//
// P and LOG10P are NA for every status other than 0 and 6.
    std::string getHeaderColumns() const override {
        return "\tP\tLOG10P\tZ\tZ_Norm\tBETA\tSE\tSPA_STATUS";
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int /*markerInChunkIdx*/,
        std::vector<double> &result
    ) override {
        result.clear();

        // Compute Score externally so the BETA = S / Var(S) reduction below
        // has Score in hand; getMarkerPvalFromScore returns Var(S) via
        // outScoreVar so the SPA-adjusted p-value and the score-test BETA
        // share the same nominal variance.
        const double gMean = GVec.mean();
        const double Score =
            GVec.dot(m_spagrm.resid()) - gMean * m_spagrm.residSum();
        double z, scoreVar;
        const spa::TwoSided ts =
            m_spagrm.getMarkerPvalFromScore(Score, altFreq, z, &scoreVar);
        pushResult(ts, Score, scoreVar, result);
    }

    // Batch analysis: fuse B dot products into one matrix-vector multiply.
    // Score[b] = GBatch.col(b).dot(resid) - gMean[b] * resid_sum
    // becomes: scores = GBatch^T * resid  (B×1, one Eigen matvec)
    void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        const std::vector<int> & /*chunkIdxs*/,
        std::vector<std::vector<double> > &results
    ) override {
        const int B = static_cast<int>(GBatch.cols());
        results.resize(B);

        // Fused: B dot products in one matrix-vector multiply.
        // GBatch is N × B, resid is N × 1 → scores is B × 1
        const Eigen::VectorXd &resid = m_spagrm.resid();
        const double resid_sum = m_spagrm.residSum();

        Eigen::VectorXd scores;
        scores.noalias() = GBatch.transpose() * resid;

        // gMeans (length B) and mean adjustment
        const Eigen::VectorXd gMeans = GBatch.colwise().mean();
        scores.array() -= gMeans.array() * resid_sum;

        // Per-marker SPA (not batchable)
        for (int b = 0; b < B; ++b) {
            double z, scoreVar;
            const spa::TwoSided ts = m_spagrm.getMarkerPvalFromScore(
                scores[b], altFreqs[b], z, &scoreVar);
            pushResult(ts, scores[b], scoreVar, results[b]);
        }
    }

    int preferredBatchSize() const override {
        // Saddlepoint-dominated (fused GEMM ~0.8% of runtime); 16 is the engine
        // floor and B>16 gives no measured gain, so request the floor exactly.
        return 16;
    }

    // ── Fused union-level GEMM interface ───────────────────────────────
    bool supportsFusedGemm() const override {
        return true;
    }

    int fusedGemmColumns() const override {
        return 1;
    }

    void fillUnionResiduals(
        Eigen::Ref<Eigen::MatrixXd> dest,
        const std::vector<uint32_t> &unionToLocal
    ) const override {
        // dest is pre-zeroed, N_union × 1.
        const Eigen::VectorXd &r = m_spagrm.resid();
        const uint32_t nUnion = static_cast<uint32_t>(unionToLocal.size());
        for (uint32_t i = 0; i < nUnion; ++i) {
            const uint32_t li = unionToLocal[i];
            if (li != UINT32_MAX)
                dest(i, 0) = r[li];
        }
    }

    void fillResidualSums(double *dest) const override {
        dest[0] = m_spagrm.residSum();
    }

    void processScoreBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &scores,
        const double *gSums,
        const double *gSumSqs,
        uint32_t nUsed,
        const std::vector<double> &altFreqs,
        const std::vector<int> & /*chunkIdxs*/,
        std::vector<std::vector<double> > &results
    ) override {
        (void)gSumSqs;
        const int B = static_cast<int>(scores.cols());
        results.resize(B);
        const double residSum = m_spagrm.residSum();
        const double invN = 1.0 / static_cast<double>(nUsed);

        for (int b = 0; b < B; ++b) {
            const double gMean = gSums[b] * invN;
            const double centeredScore = scores(0, b) - gMean * residSum;
            double z, scoreVar;
            const spa::TwoSided ts = m_spagrm.getMarkerPvalFromScore(
                centeredScore, altFreqs[b], z, &scoreVar);
            pushResult(ts, centeredScore, scoreVar, results[b]);
        }
    }

  private:
    // Emit (P, LOG10P, Z, Z_Norm, BETA, SE, SPA_STATUS) using the SPAGRM
    // score-test reduction
    //   Z_Norm = Score / sqrt(Var(S))           (raw normal-approx score z)
    //   Z      = sign(Z_Norm) · Phi^{-1}(1−p/2)  (z consistent with p, so
    //            2·pnorm(−|Z|) == p even after SPA recalibrates p)
    //   BETA   = Score / Var(S)
    //   SE     = 1 / sqrt(Var(S))
    // (so Z_Norm = BETA × SE).  Var(S) ≤ 0 marks monomorphic / degenerate
    // markers; Z, Z_Norm, BETA and SE become NaN there so downstream code
    // recognises them as missing.
    static void pushResult(
        const spa::TwoSided &ts,
        double score,
        double scoreVar,
        std::vector<double> &out
    ) {
        out.clear();
        out.reserve(7);
        out.push_back(ts.p);          // P
        out.push_back(ts.negLog10p);  // LOG10P
        if (scoreVar > 0.0) {
            const double sd = std::sqrt(scoreVar);
            const double zNorm = score / sd;
            out.push_back(math::zFromPval(ts.p, zNorm)); // Z (p-consistent)
            out.push_back(zNorm);                        // Z_Norm (raw score z)
            out.push_back(score / scoreVar);             // BETA
            out.push_back(1.0 / sd);                     // SE
        } else {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            out.push_back(nan); // Z
            out.push_back(nan); // Z_Norm
            out.push_back(nan); // BETA
            out.push_back(nan); // SE
        }
        out.push_back(static_cast<double>(static_cast<uint8_t>(ts.status)));
    }

  private:
    SPAGRMClass m_spagrm;
};

// ══════════════════════════════════════════════════════════════════════
// runSPAGRM — full workflow entry point
// ══════════════════════════════════════════════════════════════════════

void runSPAGRM(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    bool controlOutlier,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    const std::string &regressionModelStr = {},
    const std::string &phenoNameSpec = {},
    const std::string &covarFile = {},
    const std::vector<std::string> &covarNames = {},
    bool saveResid = false,
    bool longitudinal = false                             // --longitudinal: fit Y ~ X + (1|IID), use R_G
);

// LOCO variant of runSPAGRM.  For each chromosome, the null model is refit with
// that chromosome's LOCO PGS appended as one estimated covariate column, and
// the SPAGRM null model (R_GRM_R, outlier/family/Chow-Liu structures) is rebuilt
// from the refreshed residuals.  The sparse-GRM topology and IBD are residual-
// independent and are computed once.  Requires the in-process fit path
// (--pheno-name); precomputed residuals and --longitudinal are rejected in
// dispatch.  The retrospective GRM variance (R^T Phi R) is compatible with LOCO
// (the GRM models genotype covariance among relatives; LOCO adjusts the mean).
void runSPAGRMLoco(
    const std::string &phenoFile,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &predListFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    bool controlOutlier,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames
);
