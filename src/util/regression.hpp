// regression.hpp — Null-model fitting and residual construction
//
// Provides four upstream null-model fits whose residuals are designed to be
// consumed directly by SPACox / SPAGRM / SPAmix (and any other downstream
// score-test machinery in GRAB):
//
//   linearResiduals       — weighted ordinary least squares, returns y − Xβ̂
//   logisticResiduals     — weighted logistic IRLS, returns y − μ̂
//   coxResiduals          — weighted Breslow partial likelihood,
//                           returns martingale residuals δ − Λ̂₀(t) exp(Xβ̂)
//   cumulativeLogitFit    — fixed-effects cumulative-logit (proportional-odds)
//                           fit; returns β̂, ε̂ and the per-subject surrogate
//                           residual rᵢ = F̂(Yᵢ* | Xᵢ) − 1/2  where Yᵢ* is
//                           sampled from the truncated latent logistic on
//                           (ε̂_{Yᵢ−1}, ε̂_{Yᵢ}]  (Liu & Zheng, 2018, JASA)
//
// Residual-mean conventions:
//   linearResiduals    — Σ wᵢ rᵢ = 0 by construction when X contains an
//                        intercept column (XᵀW r = 0).
//   logisticResiduals  — Σ wᵢ rᵢ = 0 by IRLS convergence when X contains an
//                        intercept column.  Σ rᵢ is exactly zero only when
//                        weights are constant.
//   coxResiduals       — Σ wᵢ rᵢ = 0 by construction of the Breslow baseline
//                        hazard, independent of X.
//   cumulativeLogitFit — under H₀ the surrogate residual satisfies r ⊥⊥ X with
//                        E[r] = 0 by the probability integral transform on the
//                        latent logistic.  In a single finite sample the
//                        arithmetic mean Σᵢ rᵢ / n is O(1/√n) but not exactly
//                        zero; cumulativeLogitFit subtracts the algebraic mean
//                        before returning, restoring Σᵢ rᵢ = 0.  See §4 of
//                        docs/methods/residuals.md.
//
// All four functions accept a pre-formed design matrix and drop rows that
// contain any NaN before fitting.  calRegrWeight (case-control sampling-
// weight correction used by WtCoxG / LEAF) is intentionally kept under
// src/wtcoxg/regression.hpp because its semantics are tied to WtCoxG.
#pragma once

#include <Eigen/Dense>

namespace regression {

// ──────────────────────────────────────────────────────────────────────
// § 1  Weighted linear regression  (OLS via weighted normal equations)
//
// Fits   y = X β + ε,   ε ~ (0, σ² / w),  by solving (XᵀW X) β = XᵀW y.
// Returns response residuals  y − Xβ̂  for the rows that survive NaN
// filtering.
//
// X should include an intercept column (typically column 0 of all ones);
// otherwise the residuals will not be mean-zero.
//
// y:        (n) response
// X:        (n × p) design matrix
// weights:  (n) case weights (all > 0)
// tol, maxIter: accepted for API symmetry with logisticResiduals/coxResiduals;
//               linearResiduals solves the system in a single LDLᵀ step.
//
// Rows where any of y/weights/X contain NaN are dropped.
// Returns residuals for the kept rows (length ≤ n).
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd linearResiduals(
    const Eigen::Ref<const Eigen::VectorXd> &y,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &weights,
    double tol = 1e-9,
    int maxIter = 50
);

// ──────────────────────────────────────────────────────────────────────
// § 2  Weighted Cox regression  (Breslow partial likelihood, Newton-Raphson)
//
// Fits   h(t|X) = h₀(t) exp(Xβ)   by maximising the weighted Breslow partial
// log-likelihood via Newton-Raphson, then returns the martingale residuals
//     rᵢ = δᵢ − Λ̂₀(tᵢ) exp(Xᵢ β̂)
// weighted by wᵢ (equivalent to survival::coxph(...)$residuals in R).
//
// Convergence (matching R coxph.control logic):
//     |ℓ − ℓₒₗd| / |ℓₒₗd| < tol     OR    |Δℓ| < √tol
//
// time:     (n) survival times
// event:    (n) event indicators (1 = event, 0 = censored)
// X:        (n × p) covariate matrix (no intercept; Cox has none)
// weights:  (n) case weights (all > 0)
//
// Rows where any of time/event/weights/X contain NaN are dropped.
// Returns martingale residuals for the kept rows (length ≤ n).
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd coxResiduals(
    const Eigen::Ref<const Eigen::VectorXd> &time,
    const Eigen::Ref<const Eigen::VectorXd> &event,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &weights,
    double tol = 1e-10,
    int maxIter = 40
);

// ──────────────────────────────────────────────────────────────────────
// § 3  Weighted logistic regression  (IRLS)
//
// Fits   logit P(y = 1 | X) = X β   by Iteratively Reweighted Least Squares,
// then returns response residuals  (y − μ̂)
// (equivalent to residuals(glm(..., family=binomial()), type="response") in R).
//
// Convergence (matching R glm.control logic):
//     |dev − devₒₗd| / (|dev| + 0.1) < tol
//
// y:        (n) binary response (0 / 1)
// X:        (n × p) design matrix including an intercept column of ones
// weights:  (n) case weights (all > 0)
//
// Rows where any of y/weights/X contain NaN are dropped.
// Returns response residuals for the kept rows (length ≤ n).
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd logisticResiduals(
    const Eigen::Ref<const Eigen::VectorXd> &y,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &weights,
    double tol = 1e-9,
    int maxIter = 50
);

// ──────────────────────────────────────────────────────────────────────
// § 4  Fixed-effects cumulative-logit / proportional-odds fit
//
// Fits the cumulative-link logit model
//     logit P(Yᵢ ≤ j | Xᵢ) = εⱼ − Xᵢ β   for j = 0, …, J−2
// by Fisher scoring with score-outer-product approximation to the Hessian
// (mirrors R's MASS::polr(method = "logistic") / ordinal::clm(link = "logit")).
//
// This is the fixed-effects form of a proportional-odds model; a mixed-model
// extension would add a kinship-scaled random effect bᵢ to the linear
// predictor.  The fitted β̂ / ε̂ are exposed in the raw parameterization (ε
// is NOT shifted to ε₀ = 0); downstream consumers that require an
// identifiability convention (for example, ε₀ ≡ 0 with the intercept of β
// absorbing the shift) should apply the shift themselves.
//
// Surrogate residual returned in `residuals`:
//
//     νᵢⱼ = logistic(εⱼ − Xᵢβ̂),       ν_{i,−1} = 0,   ν_{i,J−1} = 1
//     Yᵢ* | Yᵢ = j, Xᵢ ~ Logistic(Xᵢβ̂, 1)  truncated to (εⱼ−1, εⱼ]
//     Fᵢ(Yᵢ* | Xᵢ) | Yᵢ = j, Xᵢ ~ Uniform(ν_{i,j−1}, ν_{i,j})
//     rᵢ  = Fᵢ(Yᵢ* | Xᵢ) − 1/2  ~ Uniform(ν_{i,Yᵢ−1} − 1/2,  ν_{i,Yᵢ} − 1/2)
//
// Implementation: rᵢ = ν_{i,Yᵢ−1} + Uᵢ · (ν_{i,Yᵢ} − ν_{i,Yᵢ−1}) − 1/2,
// with Uᵢ ~ Uniform(0, 1) generated by a local std::mt19937 seeded from the
// `seed` argument (`seed = 0` → std::random_device).  After sampling the
// algebraic mean is subtracted to restore Σᵢ rᵢ = 0.
//
// Under H₀ (proportional-odds model is correctly specified, no genetic
// effect) the marginal distribution of rᵢ is Uniform(−1/2, 1/2) with
// rᵢ ⊥⊥ Xᵢ.  Consequently the score statistic S = Σᵢ Gᵢrᵢ is asymptotically
// the sum of independent random variables under H₀, which is the
// configuration the saddlepoint approximation downstream is designed for.
//
// y:        (n) ordinal response, integer-coded in {0, 1, …, J−1} contiguously.
//           J is inferred from y.maxCoeff() + 1.
// X:        (n × p) design matrix.  Inclusion of an intercept column is
//           accepted (the parameterization is overspecified by 1 but the
//           Newton step is regularised with a small ridge on the Hessian).
// tol, maxIter: convergence threshold on ‖Δθ‖ and iteration cap.
// seed:     RNG seed for surrogate-residual sampling.  seed = 0 falls back
//           to std::random_device (non-reproducible).  Pass a non-zero seed
//           for reproducible p-values when ordinal phenotypes are present.
//
// Rows of X that contain any NaN are dropped, with the corresponding
// entries of y removed.
// ──────────────────────────────────────────────────────────────────────

struct CumulativeLogitFitResult {
    Eigen::VectorXd beta;      // (p)   fixed-effect coefficients
    Eigen::VectorXd eps;       // (J−1) cumulative-link thresholds
    Eigen::VectorXd residuals; // (n_kept) per-subject surrogate residual, mean-zero
};

CumulativeLogitFitResult cumulativeLogitFit(
    const Eigen::Ref<const Eigen::VectorXi> &y,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    double tol = 1e-7,
    int maxIter = 50,
    uint64_t seed = 0
);

} // namespace regression
