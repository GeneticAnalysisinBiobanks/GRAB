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
//                           fit; returns β̂, ε̂ and the per-subject analytic
//                           score residual rᵢ = ν̂_{i,Yᵢ−1} + ν̂_{i,Yᵢ} − 1
//                           (POLMM; cf. develop POLMM.cpp).  Deterministic:
//                           given (Y, β̂, ε̂) it is uniquely determined.
//
// Residual-mean conventions:
//   linearResiduals    — Σ wᵢ rᵢ = 0 by construction when X contains an
//                        intercept column (XᵀW r = 0).
//   logisticResiduals  — Σ wᵢ rᵢ = 0 by IRLS convergence when X contains an
//                        intercept column.  Σ rᵢ is exactly zero only when
//                        weights are constant.
//   coxResiduals       — Σ wᵢ rᵢ = 0 by construction of the Breslow baseline
//                        hazard, independent of X.
//   cumulativeLogitFit — the analytic score residual satisfies Xᵀr = 0 (the β
//                        estimating equations) and Σᵢ rᵢ = 0 (the threshold
//                        estimating equations) exactly at the MLE, with no RNG
//                        and no post-hoc centering.  See §4 of
//                        dev-notes/methods/residuals.md.
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
// by Fisher scoring with score-outer-product (BHHH) approximation to the
// Hessian (mirrors R's MASS::polr(method = "logistic") /
// ordinal::clm(link = "logit")).
//
// This is the fixed-effects form of a proportional-odds model; a mixed-model
// extension would add a kinship-scaled random effect bᵢ to the linear
// predictor.  The fitted β̂ / ε̂ are exposed in the raw parameterization (ε
// is NOT shifted to ε₀ = 0); downstream consumers that require an
// identifiability convention should apply the shift themselves.
//
// Analytic score residual returned in `residuals` (POLMM; cf. develop
// POLMM.cpp:75-80):
//
//     νᵢⱼ = logistic(εⱼ − Xᵢβ̂),       ν_{i,−1} = 0,   ν_{i,J−1} = 1
//     rᵢ  = ν_{i,Yᵢ−1} + ν_{i,Yᵢ} − 1  =  P(Y < Yᵢ) − P(Y > Yᵢ)
//
// rᵢ is the efficient score of the per-subject log-likelihood with respect to
// a shift in the linear predictor, i.e. rᵢ = (f_{Yᵢ−1} − f_{Yᵢ}) / μ_{i,Yᵢ}
// with f = ν(1−ν) and μ_{i,j} the category probability; the two forms are
// algebraically identical.  It is DETERMINISTIC: given (Y, β̂, ε̂) it is
// uniquely determined (no RNG).
//
// At the MLE the residual is orthogonal to the design EXACTLY (to optimizer
// tolerance):
//     Xᵀr = 0   — the β estimating equations are Σᵢ Xᵢ rᵢ = 0;
//     Σᵢ rᵢ = 0 — the threshold estimating equations telescope (ε act as the
//                 per-category intercepts), so this holds even though X
//                 carries no intercept column.
// Consequently the score statistic S = Σᵢ Gᵢrᵢ for a genetic effect is, under
// H₀, free of the covariate-explained part of G — the configuration the
// saddlepoint approximation downstream is designed for.  No post-hoc
// centering is performed.
//
// y:        (n) ordinal response, integer-coded in {0, 1, …, J−1} contiguously.
//           J is inferred from y.maxCoeff() + 1.
// X:        (n × p) covariate matrix.  Pass the covariates ALONE (NO intercept
//           column): the thresholds ε supply the per-category intercepts, so
//           an intercept column would be confounded with ε.
// tol, maxIter: convergence threshold on ‖Δθ‖ and iteration cap.
//
// Rows of X that contain any NaN are dropped, with the corresponding
// entries of y removed.
// ──────────────────────────────────────────────────────────────────────

struct CumulativeLogitFitResult {
    Eigen::VectorXd beta;      // (p)   fixed-effect coefficients
    Eigen::VectorXd eps;       // (J−1) cumulative-link thresholds
    Eigen::VectorXd residuals; // (n_kept) per-subject analytic score residual
};

CumulativeLogitFitResult cumulativeLogitFit(
    const Eigen::Ref<const Eigen::VectorXi> &y,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    double tol = 1e-7,
    int maxIter = 50
);

} // namespace regression
