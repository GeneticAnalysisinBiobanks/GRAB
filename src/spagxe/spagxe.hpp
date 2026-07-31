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
// retrospective quadratic form  xᵀΦx  evaluated by SparseGRM::quadForm.
// The base (unrelated) path is exactly the Φ = identity
// special case (phiQuad → Σx²), so a single code path serves both:
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
// SPAGxEmix — per-individual allele frequency (admixture)
// ══════════════════════════════════════════════════════════════════════
//
// SPAGxEmix (Ma et al. 2025, `SPAGxEmix_CCT`) replaces the single scalar MAF q
// of the base test with a per-individual ALT frequency q̂_i estimated from the
// SNP-derived PCs, exactly the cascade GRAB2's SPAmix already uses
// (src/spamix/indiv_af.hpp: MAC ≤ 20 → uniform; OLS of g on [1|PCs]; logistic
// fallback on the significant PCs).  The retrospective genotype law becomes
// G_i ~ Bin(2, q̂_i), so the marginal screen is re-centred by 2·Σ R_i q̂_i, the
// λ-orthogonalisation weight is variance-weighted and recomputed PER MARKER
// (λ = Σ 2q̂_i(1−q̂_i) E_i R_i² / Σ 2q̂_i(1−q̂_i) R_i²), and the SPA uses the
// per-individual afVec (the SPAmix / SPAmixPlus kernel), NOT the scalar-q path.
// A sparse GRM (SPAGxEmix+) is out of scope: mix stays on the diagonal SPAmix
// variance, never a GRM quadratic form.
//
// AFData holds the [1|PCs] design and the OLS matrices needed by computeAFVec
// (indiv_af.hpp).  One per phenotype, shared read-only across worker clones via
// shared_ptr<const>; null → the base / SPAGxE+ scalar-q path.
struct AFData {
    Eigen::MatrixXd onePlusPCs;        // N × (1+nPC) = [1 | PCs]
    Eigen::MatrixXd XtX_inv_Xt;        // (1+nPC) × N = (XᵀX)⁻¹Xᵀ  (OLS solve)
    Eigen::VectorXd sqrt_XtX_inv_diag; // (1+nPC) = sqrt(diag((XᵀX)⁻¹))  (PC SE)
    int nPC = 0;
};

// ══════════════════════════════════════════════════════════════════════
// SPAGxEMethod — MethodBase adapter (one clone per worker thread)
// ══════════════════════════════════════════════════════════════════════
//
// Output schema (mirrors SAGELD's multi-env layout, src/spagrm/sageld.cpp):
//   col 0..3        :  P_G  Z_G  BETA_G  SE_G        marginal genetic block
//                                                    (always normal-approx)
//   col 4 + 8·e+0   :  P_Gx<Ee>       final G×E p-value.  Branch A (and the
//                                     SPAGxE+ GRM / residual-mode paths): the
//                                     SPA/normal score p.  Branch B, base path:
//                                     CCT(p_spa, p_wald).
//   col 4 + 8·e+1   :  LOG10P_Gx<Ee>  −log10(P_Gx<Ee>)
//   col 4 + 8·e+2   :  P_Wald_Gx<Ee>  Branch-B Wald p of the G:E coefficient
//                                     (NaN in Branch A, the GRM path, and
//                                     residual mode — no Wald there).
//   col 4 + 8·e+3   :  Z_Gx<Ee>       p-consistent z of P_Gx<Ee>: 2·pnorm(−|Z|)==P
//   col 4 + 8·e+4   :  Z_Norm_Gx<Ee>  raw score z (Score/√Var, SPA-uncalibrated)
//   col 4 + 8·e+5   :  BETA_Gx<Ee>    score-derived interaction effect
//   col 4 + 8·e+6   :  SE_Gx<Ee>      score-derived interaction SE
//   col 4 + 8·e+7   :  SPA_STATUS_Gx<Ee>  saddlepoint outcome, as
//                                     static_cast<uint8_t>(spa::Status)
// resultSize() = 4 + 8·nEnv.
//
// spa_unify L2 / L3 (Stage 5).  LOG10P_Gx comes from the log-domain tail
// assembly wherever the reported p IS the saddlepoint p, so it survives past
// the point at which the linear-scale sum of the two tails underflows to zero.
// In Branch B with a Wald leg the reported p is CCT(p_spa, p_wald) and
// math::cauchyCombine has no log-domain form, so there LOG10P_Gx is only
// −log10 of the linear combination — the same limitation SPAsqr records for
// its P_CCT column.
//
// The marginal block deliberately gains neither column.  It is always the
// normal approximation, never a saddlepoint, so SPA_STATUS_G would be a
// constant; and Z_G is exact, so −log10(P_G) is recoverable from it to full
// precision without a stored column, which is not true of the G×E block (Z_Gx
// is derived FROM P_Gx through math::zFromPval).
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
    // afData   — per-individual-AF design (SPAGxEmix).  Non-null ⇒ mix mode:
    //            q̂_i is estimated per marker via computeAFVec, the marginal
    //            screen / λ / SPA all use afVec, and grm MUST be null (mix has
    //            no GRM path).  Null ⇒ the base / SPAGxE+ scalar-q path.
    SPAGxEMethod(
        Eigen::VectorXd resid,
        std::vector<std::string> envNames,
        std::vector<Eigen::VectorXd> envVecs,
        std::shared_ptr<const SparseGRM> grm,
        std::shared_ptr<const spagxe_wald::WaldData> wald,
        double marginalCutoff,
        double spaCutoff,
        double outlierRatio,
        std::shared_ptr<const AFData> afData = nullptr
    );

    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 4 + 8 * static_cast<int>(m_envNames.size());
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
    // Retrospective quadratic form  xᵀΦx: SparseGRM::quadForm when a GRM is
    // present, Σx² otherwise (Φ = identity → base path).
    double phiQuad(const Eigen::VectorXd &x) const;

    // Per-environment precomputed Branch-A machinery.  In the base / SPAGxE+
    // scalar-q path λ and the weight w = (E−λ)∘R are marker-independent and
    // precomputed here.  In mix mode λ is per-marker (variance-weighted on the
    // per-individual q̂_i), so only E and the naive weight E∘R (the m_W column)
    // are meaningful; lambda/w/wOutlier/wSum/wSq/wPhiQuad are left unset.
    struct EnvData {
        std::string name;
        Eigen::VectorXd E;         // environment vector (Branch B needs it raw)
        double lambda;             // λ = RᵀΦR_E / RᵀΦR       (base only)
        Eigen::VectorXd w;         // Branch-A weight (E − λ)∘R (base only)
        OutlierData wOutlier;      // IQR partition of w (SPA outlier split; base)
        double wSum;               // Σ w  (≈ 0 by residual orthogonality; base)
        double wSq;                // Σ w²         (independence quadratic form; base)
        double wPhiQuad;           // wᵀΦw         (GRM quadratic form; == wSq base)
    };

    // Evaluate one marker given the pre-multiplied raw scores.  The m_W columns
    // differ by mode, so the raw-score interpretation follows suit:
    //   base:  rawScores[0]=Σ G_i R_i,  rawScores[1+e]=Σ G_i w_{e,i}  (Branch-A score)
    //   mix:   rawScores[0]=Σ G_i R_i,  rawScores[1+e]=Σ G_i E_{e,i} R_i (naive S2;
    //          the Branch-A score S_GxE = S2 − λ_e·S_G is formed per marker).
    void evalMarker(
        const Eigen::Ref<const Eigen::VectorXd> &GVec,
        double altFreq,
        const Eigen::VectorXd &rawScores,
        std::vector<double> &out
    );

    // Per-individual-AF marker evaluation (SPAGxEmix; m_afData != null).
    void evalMarkerMix(
        const Eigen::Ref<const Eigen::VectorXd> &GVec,
        double altFreq,
        const Eigen::VectorXd &rawScores,
        std::vector<double> &out
    );

    Eigen::VectorXd m_resid;                    // R (null-model residual)
    Eigen::VectorXd m_resid2;                   // R∘R (marginal / λ variance factor)
    double m_residSum;                          // Σ R
    double m_R_phiQuad;                         // RᵀΦR (marginal variance factor; base)
    Eigen::MatrixXd m_W;                        // N × (1+nEnv): base [R|w_e], mix [R|E_e∘R]
    std::vector<std::string> m_envNames;        // env column names (header order)
    std::vector<EnvData> m_envs;                // per-env Branch-A state
    std::shared_ptr<const SparseGRM> m_grm;     // (+) variant; null → base / mix
    std::shared_ptr<const spagxe_wald::WaldData> m_wald; // Branch-B Wald; null → SPA-only
    std::shared_ptr<const AFData> m_afData;     // SPAGxEmix per-individual AF; null → scalar-q
    double m_marginalCutoff;                    // ε
    double m_spaCutoff;                         // r
    double m_outlierRatio;                      // IQR multiplier (Branch B repart.)

    // Per-thread scratch for the mix path (freshly sized per clone; reused
    // across markers).  m_afVec = q̂_i, m_wVec = 2 q̂_i(1−q̂_i).
    Eigen::VectorXd m_afVec;
    Eigen::VectorXd m_wVec;
    Eigen::VectorXd m_wScratch; // per-marker Branch-A weight (E_e − λ_e)∘R
    // Outlier-position gather of m_afVec, handed to spa_cgf::binomIndiv.  A
    // per-clone buffer rather than a per-marker std::vector: the predecessor
    // allocated one per marker per environment inside spaScorePvalMix.
    std::vector<double> m_afOutlier;
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

// ══════════════════════════════════════════════════════════════════════
// runSPAGxEmix — SPAGxEmix workflow entry point (per-individual AF)
// ══════════════════════════════════════════════════════════════════════
//
// A separate --method spagxemix.  Identical to runSPAGxE except that it builds a
// per-phenotype [1|PCs] AF design (from the --pc-cols columns, which must also
// appear in --covar-name so they adjust the null model) and threads it into each
// SPAGxEMethod, switching the marginal screen / λ / SPA to the per-individual
// q̂_i (SPAGxEmix_CCT).  A sparse GRM is NOT accepted (SPAGxEmix+ is out of
// scope).  The Branch-B Wald + CCT leg is reused unchanged (the interaction
// model trait ~ X + E + g + g:E is identical).
void runSPAGxEmix(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::vector<std::string> &pcColNames,
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
