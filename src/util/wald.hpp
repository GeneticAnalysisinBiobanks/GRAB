// wald.hpp — Standard-model Wald tests for the last coefficient of a design.
//
// Four self-contained maximum-likelihood / least-squares fitters that return
// the MAGNITUDE  L = −log10 P  of the two-sided Wald p-value of the LAST
// coefficient of a supplied design.
// They mirror the numeric patterns of src/util/regression.cpp (OLS / logistic
// IRLS / Breslow-Cox Newton / proportional-odds Fisher scoring), but where the
// regression.cpp fitters return residuals, these retain the coefficient
// covariance [information]⁻¹ and extract its (last, last) element — so the
// coefficient under test must be the last column of the supplied design.
//
// These were relocated out of src/spagxe/spagxe_wald.cpp so that any method
// needing a per-marker Wald refit can share them without depending on SPAGxE.
// The G×E-interaction dispatcher (design assembly M = [covar | g | g∘E],
// complete-case filtering, and the Cauchy fall-back to the SPA p) stays in
// src/spagxe/spagxe_wald.{hpp,cpp}.
//
// The four returned a linear p until log10p_unify Stage 7.  Both references
// they invert — the Student-t two-sided tail and the normal two-sided tail —
// underflow to exactly zero well inside the range a genome-wide scan reaches
// (2·Φ(−|z|) is zero for |z| ≳ 38.6), and the only consumer, the SPAGxE
// Branch-B Cauchy combination, is itself taken over magnitudes since Stage 5.
// Returning L removes the last place in that path where a p-value had a
// ceiling: `math::ptLog` and `math::normalTwoSidedLog` evaluate the same two
// tails in the log domain (decision D1, 01_numerics §3.3).
//
// A singular information matrix, non-convergence, or a degenerate design
// returns NaN (the caller decides how to fold that in — for the output-column
// contract that is `8 NA_NO_TEST` with LOG10P = NA).  References for each
// leg (R equivalents): Linear = lm (Student-t), Logistic = glm (normal),
// Cox = coxph with Breslow ties (normal), Ordinal = ordinal::clm (normal).
#pragma once

#include <Eigen/Dense>

namespace wald {

// Covariance estimator for the ordinal (proportional-odds) leg.
//   Observed — observed information −∂²ℓ/∂θ∂θᵀ at the MLE (analytic Hessian),
//              matching R's ordinal::clm; the default and the more accurate
//              finite-sample estimator (Efron & Hinkley, 1978).
//   BHHH     — outer-product-of-gradients information Σ sᵢsᵢᵀ; asymptotically
//              equivalent but noisier in finite samples.  Retained as an option.
enum class OrdinalInfo { Observed, BHHH };

// ── Linear: OLS, Var(β̂) = σ̂²(ZᵀZ)⁻¹, Student-t with n−k df ─────────────────
// Z includes the intercept column; the tested coefficient is Z's last column.
// Returns −log10 of the two-sided Pr(>|t|), or NaN on a singular / degenerate
// design.  The tail is `math::ptLog`, not log(2·pt(−|t|)).
double lastCoefLinearLog10P(const Eigen::MatrixXd &Z, const Eigen::VectorXd &y);

// ── Logistic: IRLS to convergence, Var(β̂) = (ZᵀWZ)⁻¹, normal reference ─────
// Z includes the intercept column; the tested coefficient is Z's last column.
// Returns −log10 of the two-sided Pr(>|z|), or NaN.
double lastCoefLogisticLog10P(const Eigen::MatrixXd &Z, const Eigen::VectorXd &y);

// ── Cox: Breslow partial-likelihood Newton, Var(β̂) = H⁻¹, normal reference ──
// X has no intercept (baseline hazard); the tested coefficient is X's last
// column.  Ties handled by the Breslow approximation (as regression::coxResiduals).
// Returns −log10 of the two-sided Pr(>|z|), or NaN.
double lastCoefCoxLog10P(
    const Eigen::VectorXd &time,
    const Eigen::VectorXd &event,
    const Eigen::MatrixXd &X);

// ── Ordinal: proportional-odds Fisher scoring, normal reference ─────────────
// X has no intercept (thresholds are the intercepts); the tested coefficient is
// X's last column.  y is integer-coded 0..J−1; the parameter order is
// θ = [ε(J−1) | β(p)], so the tested coefficient is the final θ element.
// `info` selects the covariance estimator (default: observed, clm-matching).
// Returns −log10 of the two-sided Pr(>|z|), or NaN.
double lastCoefOrdinalLog10P(
    const Eigen::VectorXi &y,
    const Eigen::MatrixXd &X,
    OrdinalInfo info = OrdinalInfo::Observed);

}  // namespace wald
