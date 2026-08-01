// spacox.hpp — SPACox: Saddlepoint Approximation for Cox Model
//
// Port of ref_code/src/mtSPACox into the pure C++17 / Eigen framework.
//
// Workflow:
//   1. Pre-compute empirical CGF interpolation tables from residuals
//   2. Per-marker: score test → normal approx or SPA tail probability
//      Stage 1: unadjusted genotype (bucketed weights, see spacox_cgf.hpp)
//      Stage 2: covariate-adjusted if p < pVal_covaAdj_cutoff
//   3. Output: [P, LOG10P, Z, Z_Norm, BETA, SE, SPA_STATUS]
//
// The saddlepoint root finder and the Barndorff-Nielsen tail are the shared
// ones in util/spa.hpp; the empirical CGF is SPACox's own and lives in
// spacox_cgf.hpp.  See dev-notes/methods/spa_unify/02_design.md.
#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "spacox/spacox_cgf.hpp"
#include "util/spa.hpp"

// ======================================================================
// DesignMatrix — covariate projection (spacox-only)
// ======================================================================

class DesignMatrix {
  public:
    explicit DesignMatrix(const Eigen::MatrixXd &X);

    int nRows() const {
        return static_cast<int>(m_X.rows());
    }

    int nCols() const {
        return static_cast<int>(m_X.cols());
    }

    const Eigen::MatrixXd &X() const {
        return m_X;
    }

    const Eigen::MatrixXd &tX() const {
        return m_tX;
    }

    const Eigen::MatrixXd &XinvXX() const {
        return m_XinvXX;
    }

    void adjustGenotype(
        const double *G,
        const uint32_t *nzIdx,
        int nNz,
        Eigen::Ref<Eigen::VectorXd> adjG
    ) const;

  private:
    Eigen::MatrixXd m_X;      // N × p
    Eigen::MatrixXd m_tX;     // p × N
    Eigen::MatrixXd m_XinvXX; // N × p
};

// ======================================================================
// SPACoxMethod — MethodBase implementation
// ======================================================================

class SPACoxMethod : public MethodBase {
  public:
// Construct from pre-computed CGF table, residuals, and design matrix.
// All large const objects are taken by const reference (shared, read-only).
    SPACoxMethod(
        const Eigen::VectorXd &residuals,
        double varResid,
        const CumulantTable &cumul,
        const DesignMatrix &design,
        double pvalCovAdjCut,
        double spaCutoff
    );

// ---- MethodBase interface ----
    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 7;
    }

// P, LOG10P, Z, Z_Norm, BETA, SE, SPA_STATUS.
//
// LOG10P is -log10(P), computed through spa::bnTailLog / spa::combineTailsLog
// so that it stays meaningful past the point where the linear-scale tail
// underflows (Phi(-38.5) flushes to zero, i.e. p ~ 1e-316).
//
// SPA_STATUS carries the spa::Status of the saddlepoint that produced P.  It
// is emitted as the integer enumerator rather than the token spelled by
// spa::statusName because MethodBase hands the engine a std::vector<double>
// and every result cell is formatted by numToChars; a string column would
// require a new hook in the MethodBase contract, which
// dev-notes/methods/spa_unify/02_design.md places out of scope for the
// per-method migration stages.  The mapping is spa::Status's own:
//
//     0 SPA_OK               3 FALLBACK_MAXITER      7 NA_POST_FAIL
//     1 NORMAL               4 FALLBACK_GUARD_TEMP   8 NA_NO_TEST
//     2 SPA_W_SINGULAR       5 FALLBACK_GUARD_CURV
//                            6 FALLBACK_NONFINITE
//
// The order is a contract (log10p_unify decision D4): SPA_STATUS <= 2 means
// LOG10P is trustworthy, 3..6 that the saddlepoint failed and LOG10P is the
// substituted two-sided normal tail -log10(2*Phi(-|Z_Norm|)) -- which is
// anti-conservative at low MAC, so filter with SPA_STATUS <= 2 before judging
// significance -- and >= 7 that LOG10P is NA.
    std::string getHeaderColumns() const override {
        return "\tP\tLOG10P\tZ\tZ_Norm\tBETA\tSE\tSPA_STATUS";
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result
    ) override;

// Batched per-marker analysis.  Folds the per-marker dot product
// S_b = GBatch.col(b)ᵀ · m_resid into a single Eigen GEMV across all
// B markers in the batch (BLAS-2), then delegates to the same
// scalar SPA core used by getResultVec for the per-marker work.
    void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        const std::vector<int> &chunkIdxs,
        std::vector<std::vector<double> > &results
    ) override;

    int preferredBatchSize() const override {
        return 16;
    }

  private:
// ---- Two-sided saddlepoint p-value ----
//
// Both tails are solved with spa::solveSaddlepoint and evaluated with
// spa::bnTail / spa::bnTailLog, then assembled by spa::assemble.  The
// only SPACox-specific part is the CGF callable, which is why these take one.
//
// `getProbSpaBucketed` is stage 1 (weights collapsed into buckets, P1);
// `getProbSpaDense` is stage 2 (covariate-adjusted weights, no bucket
// structure).  Both return P, -log10(P) and the status of the worse tail.
// `zNorm` is the raw score z; spa::assemble reports its normal tail, under a
// status naming the substitution, when the saddlepoint fails (decision D5).
// SPACox standardizes its weights, so |zNorm| equals `absZ` — the two are
// still passed separately so the call site, not a comment, is what has to
// hold.
    spa::TwoSided getProbSpaBucketed(
        const spacox_cgf::GenoWeights &gw,
        double absZ,
        double zNorm
    ) const;

    spa::TwoSided getProbSpaDense(
        const double *w,
        int n,
        double absZ,
        double zNorm
    ) const;

// ---- Per-marker score test ----
//
// `S` is the pre-computed inner product GVec · m_resid; passing it in
// lets getResultBatch fold all B such dot products into a single GEMV
// and avoids recomputing the per-marker dot inside the SPA core.
//
// `outScoreVar` returns the score variance used to scale `S`: the
// stage-1 value `varResid · Σ (G − 2·MAF)²` or, when the covariate-
// adjusted SPA branch fires, the stage-2 value
// `varResid · ‖adjGVec‖²`.  Reported BETA = S / Var(S),
// SE = 1 / sqrt(Var(S)) consume this variance so that both stages
// share the same nominal Fisher information.  Set to zero for
// degenerate markers where Var(S) ≤ 0.
    spa::TwoSided getMarkerPvalCore(
        const Eigen::Ref<const Eigen::VectorXd> &GVec,
        double altFreq,
        double S,
        double &zScore,
        double &outScoreVar
    );

// Push P, LOG10P, Z, Z_Norm, BETA, SE, SPA_STATUS in header order.
    static void pushResult(
        std::vector<double> &out,
        const spa::TwoSided &ts,
        double S,
        double scoreVar
    );

// Read-only shared data (references are stable because the owner outlives all clones).
// SPA scratch (adjGNorm / adjGVec / nzSet) is held in a translation-unit-local
// thread_local struct inside spacox.cpp, so K phenotype-clones running in the
// same worker thread share one scratch — the SPACoxMethod object itself only
// stores references and small scalars and is therefore cheap to clone.
    const Eigen::VectorXd &m_resid;
    double m_varResid;
    const CumulantTable &m_cumul;
    const DesignMatrix &m_design;
    int m_N;
    double m_pvalCovAdjCut;
    double m_spaCutoff;
};

// ======================================================================
// Orchestration
// ======================================================================

void runSPACox(
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &covarNames,             // empty = no covariates
    const std::string &phenoFile,                            // pheno file (for residuals + covar columns)
    const std::string &covarFile,                          // covar file (for covar columns)
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double pvalCovAdjCut,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    const std::string &regressionModelStr = {},            // empty → residual-passthrough path
    const std::string &phenoNameSpec = {},                 // raw --pheno-name string (Cox: TIME:EVENT,...)
    bool saveResid = false,
    bool longitudinal = false                             // --longitudinal: fit Y ~ X + (1|IID), use R_G
);

// LOCO variant of runSPACox.  For each chromosome, the null model is refit with
// that chromosome's LOCO PGS appended as one estimated covariate column (the
// regression estimates its coefficient — no offset mechanism), then SPACox
// score tests run on that chromosome's markers via the generic locoEngine.
// Requires the in-process fit path (--pheno-name); precomputed residuals and
// --longitudinal are rejected upstream in dispatch.  The pred.list is keyed on
// the canonical phenotype/spec name (e.g. "Time_Event" for a TIME:EVENT spec).
void runSPACoxLoco(
    const std::vector<std::string> &covarNames,
    const std::string &phenoFile,
    const std::string &covarFile,
    const GenoSpec &geno,
    const std::string &predListFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double pvalCovAdjCut,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec
);
