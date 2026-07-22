// spagxe.hpp — SPAGxE_CCT gene–environment interaction test (pure C++17 / Eigen)
//
// Retrospective, residual-based, saddlepoint (SPA) G×E score test
// (Ma, Zhao, Zhang & Bi, Nat. Commun. 2025).  A single genotype-independent
// null model is fit once (trait ~ X + E); per variant the marginal genetic
// effect is screened and the G×E interaction is tested by a retrospective SPA.
//
//   --method spagxe   base SPAGxE (unrelated)  +  SPAGxE+ (via --sp-grm-*)
//                     +  SPAGxE_CCT (Branch-B Wald + Cauchy combination).
//
// Design references (git-ignored, in the feat/sadge tree):
//   dev-notes/methods/SPAGxE_claude            model + equations + R conflicts
//   dev-notes/methods/SPAGxE_claude_plan       GRAB2 reuse / design map
//   dev-notes/methods/SPAGxE_claude_impl_plan  phase-by-phase build plan
//
// Base vs +: a sparse GRM Φ enters only through the score variance, as a
// retrospective quadratic form  xᵀΦx  evaluated by SparseGRM::spaVariance
// (2·Σstored − Σx²).  The base (unrelated) path is exactly the Φ = identity
// special case (spaVariance → Σx²), so a single code path serves both:
//   marginal   Var(S_G)      = 2q(1−q)·RᵀΦR
//   Branch A   λ_GRM         = RᵀΦR_E / RᵀΦR ,  w = (E − λ_GRM)R
//              Var(S_GxE)     = 2q(1−q)·wᵀΦw
//   Branch B   Var(Ŝ_GxE)    = 2q(1−q)·(E∘R̃)ᵀΦ(E∘R̃)   (per marker)
// The SPA is applied to the independence CGF with a SAIGE-style variance ratio
// √(indepVar / grmVar) rescaling the statistic (paper eq. 14–17; the non-mix
// convention, reflecting about the fitted mean).  In the significant-marginal
// Branch B the paper's SPAGxE+ runs no Wald test.
#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "engine/marker.hpp"          // MethodBase
#include "geno_factory/geno_data.hpp" // GenoSpec
#include "spagxe/spagxe_wald.hpp"     // spagxe_wald::WaldData (Branch-B Wald leg)
#include "util/outlier.hpp"           // OutlierData

class SparseGRM; // held by shared_ptr; full type only needed in the .cpp

// ══════════════════════════════════════════════════════════════════════
// SPAGxEMethod — MethodBase adapter (one clone per worker thread)
// ══════════════════════════════════════════════════════════════════════
//
// Output schema (mirrors SAGELD's multi-env layout, src/spagrm/sageld.cpp):
//   col 0..3        :  P_G  Z_G  BETA_G  SE_G        marginal genetic block
//                                                    (always normal-approx)
//   col 4 + 6·e+0   :  P_Gx<Ee>       final G×E p-value.  Branch A (and the
//                                     SPAGxE+ GRM / residual-mode paths): the
//                                     SPA/normal score p.  Branch B, base path:
//                                     CCT(p_spa, p_wald).
//   col 4 + 6·e+1   :  P_Wald_Gx<Ee>  Branch-B Wald p of the G:E coefficient
//                                     (NaN in Branch A, the GRM path, and
//                                     residual mode — no Wald there).
//   col 4 + 6·e+2   :  Z_Gx<Ee>       p-consistent z of P_Gx<Ee>: 2·pnorm(−|Z|)==P
//   col 4 + 6·e+3   :  Z_Norm_Gx<Ee>  raw score z (Score/√Var, SPA-uncalibrated)
//   col 4 + 6·e+4   :  BETA_Gx<Ee>    score-derived interaction effect
//   col 4 + 6·e+5   :  SE_Gx<Ee>      score-derived interaction SE
// resultSize() = 4 + 6·nEnv.
class SPAGxEMethod : public MethodBase {
  public:
    // resid    — the null-model residual R (per-phenotype dense, length N).
    // envNames — environment column names, one 5-wide output block each.
    // envVecs  — the environment vectors E_e (aligned to resid), one per env.
    // grm      — per-phenotype sparse GRM for the (+) variant; null → base.
    // wald     — per-phenotype Branch-B Wald data (raw phenotype + covariates).
    //            null, or trait==None, or a non-null grm ⇒ no Wald leg (the
    //            SPAGxE+ GRM path and residual mode stay SPA-only).
    // marginalCutoff — ε, the Branch A / Branch B routing threshold.
    // spaCutoff      — r (--spa-z-threshold), the |z| normal↔SPA switch.
    // outlierRatio   — IQR multiplier for the SPA outlier partition.
    SPAGxEMethod(
        Eigen::VectorXd resid,
        std::vector<std::string> envNames,
        std::vector<Eigen::VectorXd> envVecs,
        std::shared_ptr<const SparseGRM> grm,
        std::shared_ptr<const spagxe_wald::WaldData> wald,
        double marginalCutoff,
        double spaCutoff,
        double outlierRatio
    );

    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 4 + 6 * static_cast<int>(m_envNames.size());
    }

    std::string getHeaderColumns() const override;

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result
    ) override;

    // Batched: one GEMM  scores = GBatchᵀ · [R | w_1 | … | w_nEnv]  supplies the
    // marginal and every Branch-A G×E score; per-marker Branch B (rare) uses the
    // genotype column directly.
    void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        const std::vector<int> &chunkIdxs,
        std::vector<std::vector<double> > &results
    ) override;

  private:
    // Retrospective quadratic form  xᵀΦx: SparseGRM::spaVariance when a GRM is
    // present (2·Σstored − Σx²), Σx² otherwise (Φ = identity → base path).
    double phiQuad(const Eigen::VectorXd &x) const;

    // Per-environment precomputed Branch-A machinery (marker-independent).
    struct EnvData {
        std::string name;
        Eigen::VectorXd E;         // environment vector (Branch B needs it raw)
        double lambda;             // λ = RᵀΦR_E / RᵀΦR
        Eigen::VectorXd w;         // Branch-A weight (E − λ)∘R
        OutlierData wOutlier;      // IQR partition of w (SPA outlier split)
        double wSum;               // Σ w  (≈ 0 by residual orthogonality)
        double wSq;                // Σ w²         (independence quadratic form)
        double wPhiQuad;           // wᵀΦw         (GRM quadratic form; == wSq base)
    };

    // Evaluate one marker given the pre-multiplied raw scores
    //   rawScores[0]      = Σ G_i R_i             (marginal)
    //   rawScores[1 + e]  = Σ G_i w_{e,i}         (Branch-A G×E, env e)
    void evalMarker(
        const Eigen::Ref<const Eigen::VectorXd> &GVec,
        double altFreq,
        const Eigen::VectorXd &rawScores,
        std::vector<double> &out
    );

    Eigen::VectorXd m_resid;                    // R (null-model residual)
    double m_residSum;                          // Σ R
    double m_R_phiQuad;                         // RᵀΦR (marginal variance factor)
    Eigen::MatrixXd m_W;                        // N × (1+nEnv): [R | w_1 | … ]
    std::vector<std::string> m_envNames;        // env column names (header order)
    std::vector<EnvData> m_envs;                // per-env Branch-A state
    std::shared_ptr<const SparseGRM> m_grm;     // (+) variant; null → base
    std::shared_ptr<const spagxe_wald::WaldData> m_wald; // Branch-B Wald; null → SPA-only
    double m_marginalCutoff;                    // ε
    double m_spaCutoff;                         // r
    double m_outlierRatio;                      // IQR multiplier (Branch B repart.)
};

// ══════════════════════════════════════════════════════════════════════
// runSPAGxE — full workflow entry point
// ══════════════════════════════════════════════════════════════════════
//
// Mirrors runSPAGRM (src/spagrm/spagrm.cpp): load SubjectData, fit or load the
// genotype-independent residual, extract each --envir-name column, build one
// PhenoTask per phenotype, and drive multiPhenoEngine.  When a sparse GRM is
// supplied (--sp-grm-*, optional) the (+) variance correction engages; absent,
// the base unrelated path runs.
void runSPAGxE(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double marginalCutoff,
    double spaCutoff,
    double outlierIqrRatio,
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
    bool saveResid = false
);
