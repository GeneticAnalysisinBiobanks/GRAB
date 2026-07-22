// spagxe_wald.hpp — Branch-B Wald leg of SPAGxE_CCT (self-contained, per marker)
//
// SPAGxE_CCT (Ma et al., Nat. Commun. 2025) combines, for the small fraction of
// variants with a significant marginal genetic effect (Branch B, p_marg ≤ ε), a
// retrospective saddlepoint G×E p-value with a prospective Wald p-value of the
// interaction coefficient in the FULL model
//
//     trait ~ [covar] + g + g:E          (model doc §4.2, R/SPAGxECCT.R:620)
//
// and reports  P = CCT(p_spa, p_wald)  (Cauchy combination; math::cauchyCombine).
//
// No GRAB2 residual fitter exposes coefficient covariance, so the per-marker
// Wald refit is net-new for every trait type.  This header declares a
// self-contained fitter (implemented in spagxe_wald.cpp) that mirrors the
// numeric patterns of src/util/regression.cpp but returns the Wald p-value of
// the G:E coefficient instead of a residual.  The shared residual fitters are
// left untouched (guard rails: do not modify shared/validated code).
//
// Scope: fires only for the base (non-GRM) path.  SPAGxE+ (sparse GRM) keeps no
// Wald test (paper: per-marker mixed-model refit is too expensive), and residual
// mode (--resid-name) has no raw phenotype, so both leave the trait as `None`.
#pragma once

#include <Eigen/Dense>

namespace spagxe_wald {

// Trait family of the null model, resolved per phenotype in runSPAGxE.
// `None` disables the Wald leg (residual mode, or the SPAGxE+ GRM path).
enum class TraitType { None, Linear, Logistic, Cox, Ordinal };

// Per-phenotype raw data for the full-interaction-model refit.  All
// vectors/matrices are in the phenotype's dense subject order (aligned to the
// method's residual, environment, and genotype vectors).  `covar` already
// contains the environment column(s) and the PCs as main effects — exactly the
// covariate design of the genotype-independent null model, minus the intercept
// (which Linear/Logistic add internally; Cox/Ordinal supply it through the
// baseline hazard / thresholds).  Held behind a shared_ptr<const>, one per
// phenotype, shared across worker-thread clones.
struct WaldData {
    TraitType trait = TraitType::None;
    Eigen::MatrixXd covar;  // (n × c) main-effect covariates: PCs + env(s); no intercept
    Eigen::VectorXd y;      // (n) Linear (raw) / Logistic ({0,1}) / Ordinal ({0..J−1}) response
    Eigen::VectorXd time;   // (n) Cox survival time (unused otherwise)
    Eigen::VectorXd event;  // (n) Cox event indicator {0,1} (unused otherwise)
};

// Two-sided Wald p-value of the G:E interaction coefficient in the full model
//   trait ~ [covar] + g + g:E ,
// where the interaction regressor is the element-wise product g∘E.  Uses a
// Student-t reference for Linear (matching lm) and a normal reference for
// Logistic / Cox / Ordinal (matching glm / coxph / clm).  Returns NaN on a
// singular information matrix, non-convergence, or a degenerate design; the
// caller then falls back to the SPA p through the NaN-skipping Cauchy
// combination.  `trait == None` also returns NaN.
double waldInteractionPval(
    const WaldData &wd,
    const Eigen::Ref<const Eigen::VectorXd> &g,
    const Eigen::Ref<const Eigen::VectorXd> &E);

}  // namespace spagxe_wald
