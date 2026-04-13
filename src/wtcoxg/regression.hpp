// regression.hpp — Weighted Cox and logistic regression residuals
//
// Provides:
//   calRegrWeight        — case-control sampling weights from prevalence
//   coxResiduals         — weighted Cox PH → martingale residuals
//   logisticResiduals    — weighted logistic IRLS → response residuals
//
// Both regression functions accept a pre-formed design matrix (with
// intercept column) and drop rows containing any NaN before fitting.
#pragma once

#include <Eigen/Dense>
#include <string>

namespace regression {

// ──────────────────────────────────────────────────────────────────────
// calRegrWeight — case-control sampling-weight correction
//
// Equivalent to R:
//   sample_ratio <- sum(Indicator) / sum(1 - Indicator)
//   population_ratio <- RefPrevalence / (1 - RefPrevalence)
//   weight[control] <- sample_ratio / population_ratio
//
// prevalence: disease prevalence in reference population (0,1)
// indicator:  (n) case/control indicator (1 = case, 0 = control)
// Returns: weight vector (n)
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd calRegrWeight(
    double prevalence,
    const Eigen::Ref<const Eigen::VectorXd> &indicator
);

// ──────────────────────────────────────────────────────────────────────
// § 1  Weighted Cox regression  (Breslow partial likelihood, Newton-Raphson)
//
// Fits the model   h(t|X) = h0(t) exp(X β)   by maximising the weighted
// Breslow partial log-likelihood via Newton-Raphson, then returns the
// martingale residuals   r_i = δ_i − Λ̂_0(t_i) exp(X_i β̂)
// weighted by w_i (equivalent to survival::coxph(...)$residuals in R).
//
// Convergence: relative log-partial-likelihood change + absolute fallback
// (matching R coxph.control logic):
//   |loglik - loglikOld| / |loglikOld| < tol  OR  |loglik - loglikOld| < √tol
//
// time:     (n) survival times
// event:    (n) event indicators (1 = event, 0 = censored)
// X:        (n × p) covariate matrix (NO intercept — Cox has none)
// weights:  (n) case weights (all > 0)
//
// Rows where any of time/event/weight/X contain NaN are dropped.
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
// § 2  Weighted logistic regression  (IRLS)
//
// Fits   logit P(y=1|X) = X β   by Iteratively Reweighted Least Squares,
// then returns response residuals  (y − μ̂)
// (equivalent to residuals(glm(..., family=binomial()), type="response") in R).
//
// Convergence: relative deviance change (matching R glm.control logic)
//   |dev - devOld| / (|dev| + 0.1) < tol
//
// y:        (n) binary response (0/1)
// X:        (n × p) design matrix including intercept column of 1s
// weights:  (n) case weights (all > 0)
//
// Rows where any of y/weight/X contain NaN are dropped.
// Returns response residuals for the kept rows (length ≤ n).
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd logisticResiduals(
    const Eigen::Ref<const Eigen::VectorXd> &y,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &weights,
    double tol = 1e-9,
    int maxIter = 50
);

} // namespace regression
