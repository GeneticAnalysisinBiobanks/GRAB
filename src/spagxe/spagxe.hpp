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
// Phase 0 scaffold: CLI registration + MethodBase interface only.  The base
// SPA-only marker logic (Branch A / Branch B) lands in Phase 1; the sparse-GRM
// (+) variant in Phase 2.
#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "engine/marker.hpp"          // MethodBase
#include "geno_factory/geno_data.hpp" // GenoSpec

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
    // resid    — the null-model residual R (union-dense, length N).
    // envNames — environment column names, one 5-wide output block each.
    // marginalCutoff — ε, the Branch A / Branch B routing threshold.
    // spaCutoff      — r (--spa-z-threshold), the |z| normal↔SPA switch.
    SPAGxEMethod(
        Eigen::VectorXd resid,
        std::vector<std::string> envNames,
        double marginalCutoff,
        double spaCutoff
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

  private:
    Eigen::VectorXd m_resid;               // R (null-model residual)
    std::vector<std::string> m_envNames;   // environment column names
    double m_marginalCutoff;               // ε
    double m_spaCutoff;                    // r
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
