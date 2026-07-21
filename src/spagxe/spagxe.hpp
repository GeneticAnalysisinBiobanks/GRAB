// spagxe.hpp — SPAGxE_CCT gene–environment interaction test (pure C++17 / Eigen)
//
// Retrospective, residual-based, saddlepoint (SPA) G×E score test
// (Ma, Zhao, Zhang & Bi, Nat. Commun. 2025).  A single genotype-independent
// null model is fit once (trait ~ X + E); per variant the marginal genetic
// effect is screened and the G×E interaction is tested by a retrospective SPA.
//
//   --method spagxe   base SPAGxE (unrelated)  +  SPAGxE+ (via --sp-grm-*)
//                     +  SPAGxE_CCT (Wald + CCT, added in a later phase).
//
// Design references (git-ignored, in the feat/sadge tree):
//   dev-notes/methods/SPAGxE_claude            model + equations + R conflicts
//   dev-notes/methods/SPAGxE_claude_plan       GRAB2 reuse / design map
//   dev-notes/methods/SPAGxE_claude_impl_plan  phase-by-phase build plan
//
// Phase 1: base SPAGxE, SPA-only, unrelated.  Branch A (λ-orthogonalised G×E
// score) and Branch B (genotype-adjusted residual) both go through the SPA /
// normal hybrid; no Wald / CCT yet (that is Phase 3), and no sparse-GRM
// variance (Phase 2).
#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "engine/marker.hpp"          // MethodBase
#include "geno_factory/geno_data.hpp" // GenoSpec
#include "util/outlier.hpp"           // OutlierData

// ══════════════════════════════════════════════════════════════════════
// SPAGxEMethod — MethodBase adapter (one clone per worker thread)
// ══════════════════════════════════════════════════════════════════════
//
// Output schema (mirrors SAGELD's multi-env layout, src/spagrm/sageld.cpp):
//   col 0..3        :  P_G  Z_G  BETA_G  SE_G        marginal genetic block
//                                                    (always normal-approx)
//   col 4 + 5·e+0   :  P_Gx<Ee>       G×E p-value (SPA/normal; CCT in Phase 3)
//   col 4 + 5·e+1   :  Z_Gx<Ee>       p-consistent z:  2·pnorm(−|Z|) == P
//   col 4 + 5·e+2   :  Z_Norm_Gx<Ee>  raw score z (Score/√Var, SPA-uncalibrated)
//   col 4 + 5·e+3   :  BETA_Gx<Ee>
//   col 4 + 5·e+4   :  SE_Gx<Ee>
// resultSize() = 4 + 5·nEnv.
class SPAGxEMethod : public MethodBase {
  public:
    // resid    — the null-model residual R (union-dense per phenotype, length N).
    // envNames — environment column names, one 5-wide output block each.
    // envVecs  — the environment vectors E_e (aligned to resid), one per env.
    // marginalCutoff — ε, the Branch A / Branch B routing threshold.
    // spaCutoff      — r (--spa-z-threshold), the |z| normal↔SPA switch.
    // outlierRatio   — IQR multiplier for the SPA outlier partition.
    SPAGxEMethod(
        Eigen::VectorXd resid,
        std::vector<std::string> envNames,
        std::vector<Eigen::VectorXd> envVecs,
        double marginalCutoff,
        double spaCutoff,
        double outlierRatio
    );

    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 4 + 5 * static_cast<int>(m_envNames.size());
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
    // Per-environment precomputed Branch-A machinery (marker-independent).
    struct EnvData {
        std::string name;
        Eigen::VectorXd E;         // environment vector (Branch B needs it raw)
        double lambda;             // λ = Σ E R² / Σ R²
        Eigen::VectorXd w;         // Branch-A weight (E − λ)∘R
        OutlierData wOutlier;      // IQR partition of w (SPA outlier split)
        double wSum;               // Σ w  (≈ 0 by residual orthogonality)
        double wSq;                // Σ w²
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

    Eigen::VectorXd m_resid;               // R (null-model residual)
    double m_residSum;                     // Σ R
    double m_residSq;                      // Σ R²
    Eigen::MatrixXd m_W;                   // N × (1 + nEnv): [R | w_1 | … | w_nEnv]
    std::vector<std::string> m_envNames;   // env column names (header order)
    std::vector<EnvData> m_envs;           // per-env Branch-A state
    double m_marginalCutoff;               // ε
    double m_spaCutoff;                    // r
    double m_outlierRatio;                 // IQR multiplier (Branch B re-partition)
};

// ══════════════════════════════════════════════════════════════════════
// runSPAGxE — full workflow entry point
// ══════════════════════════════════════════════════════════════════════
//
// Mirrors runSPAGRM (src/spagrm/spagrm.cpp): load SubjectData, fit or load the
// genotype-independent residual, extract each --envir-name column, build one
// PhenoTask per phenotype, and drive multiPhenoEngine.  When a sparse GRM is
// supplied (--sp-grm-*, optional) the (+) variance correction engages (Phase 2);
// absent, the base unrelated path runs.
void runSPAGxE(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
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
