// regression.cpp — Weighted Cox and logistic regression residuals
#include "regression.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace regression {

// ──────────────────────────────────────────────────────────────────────
// calRegrWeight — case-control sampling-weight correction
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd calRegrWeight(
    double prevalence,
    const Eigen::Ref<const Eigen::VectorXd> &indicator
) {
    const Eigen::Index N = indicator.size();
    double nCase = 0;
    for (Eigen::Index i = 0; i < N; ++i)
        if (indicator[i] > 0.5) ++nCase;
    double nCtrl = N - nCase;

    if (nCase < 1 || nCtrl < 1) throw std::runtime_error("calRegrWeight: need at least 1 case and 1 control");
    if (prevalence <= 0.0 || prevalence >= 1.0) throw std::runtime_error("calRegrWeight: prevalence must be in (0, 1)");

    // R equivalent:
    //   sample_ratio <- sum(Indicator) / sum(1 - Indicator)
    //   population_ratio <- RefPrevalence / (1 - RefPrevalence)
    //   weight[control] <- sample_ratio / population_ratio
    double sample_ratio = nCase / nCtrl;
    double population_ratio = prevalence / (1.0 - prevalence);
    double wCtrl = sample_ratio / population_ratio;

    Eigen::VectorXd w = Eigen::VectorXd::Ones(N);
    for (Eigen::Index i = 0; i < N; ++i)
        if (indicator[i] < 0.5) w[i] = wCtrl;
    return w;
}

// ──────────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────────

namespace {

// Return indices of rows where no element is NaN across all vectors/matrix.
std::vector<Eigen::Index> completeRows(
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    std::initializer_list<const Eigen::Ref<const Eigen::VectorXd> *> vecs
) {
    const Eigen::Index n = X.rows();
    std::vector<Eigen::Index> keep;
    keep.reserve(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        bool ok = true;
        for (auto *v : vecs)
            if (std::isnan((*v)[i])) {
                ok = false;
                break;
            }
        if (!ok) continue;
        for (Eigen::Index j = 0; j < X.cols(); ++j)
            if (std::isnan(X(i, j))) {
                ok = false;
                break;
            }
        if (ok) keep.push_back(i);
    }
    return keep;
}

// Subset rows by index vector.
Eigen::MatrixXd subsetRows(
    const Eigen::Ref<const Eigen::MatrixXd> &M,
    const std::vector<Eigen::Index> &idx
) {
    Eigen::MatrixXd out(static_cast<Eigen::Index>(idx.size()), M.cols());
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(idx.size()); ++i)
        out.row(i) = M.row(idx[i]);
    return out;
}

Eigen::VectorXd subsetVec(
    const Eigen::Ref<const Eigen::VectorXd> &v,
    const std::vector<Eigen::Index> &idx
) {
    Eigen::VectorXd out(static_cast<Eigen::Index>(idx.size()));
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(idx.size()); ++i)
        out[i] = v[idx[i]];
    return out;
}

// Validate post-NaN-drop data before Newton/IRLS.
//
// Checks (in order):
//   1. n > 0              — not empty after NaN removal
//   2. n > p              — overdetermined system (enough subjects)
//   3. all w[i] > 0       — positive weights
//   4. no constant column — zero-variance covariate means collinearity
//
// startCol: first column to check for zero variance.
//   Pass 0 for Cox (no intercept), 1 for logistic (col 0 = all-ones intercept).
void validateSubset(
    const std::string &fname,
    Eigen::Index n,
    Eigen::Index p,
    const Eigen::VectorXd &w,
    const Eigen::MatrixXd &X,
    Eigen::Index startCol = 0
) {
    if (n == 0) throw std::runtime_error(fname + ": no complete cases after NaN removal");
    if (n <= p)
        throw std::runtime_error(fname + ": n (" + std::to_string(n) + ") must exceed p (" + std::to_string(p) +
                                 ") — too few subjects for the number of covariates");
    for (Eigen::Index i = 0; i < n; ++i) {
        if (!(w[i] > 0.0))
            throw std::runtime_error(fname + ": non-positive weight (" + std::to_string(w[i]) + ") at row " +
                                     std::to_string(i));
    }
    for (Eigen::Index j = startCol; j < p; ++j) {
        if (X.col(j).maxCoeff() - X.col(j).minCoeff() < 1e-14)
            throw std::runtime_error(fname + ": covariate column " + std::to_string(j) +
                                     " has zero variance"
                                     " (constant, duplicate, or perfectly collinear column)");
    }
}

} // anonymous namespace

// ──────────────────────────────────────────────────────────────────────
// § 1  Weighted Cox regression — Breslow partial likelihood
//
// Algorithm overview (single-pass cumulative sums per Newton step):
//   1. Sort subjects by descending survival time.
//   2. Compute θ_i = exp(X_i β), weighted risk quantities via
//      reverse-cumulative sums over the sorted order.
//   3. Score U = X^T W (δ − dΛ exp(Xβ))  and
//      Hessian H are accumulated in one forward pass over events.
//   4. Newton update: β ← β + H^{-1} U.
//   5. After convergence, martingale residuals = δ_i − Ĥ(t_i) θ_i.
//
// Convergence: relative log-partial-likelihood change + absolute fallback
// (matching R coxph.control):
//   |loglik - loglikOld| / |loglikOld| < tol  OR  |Δloglik| < √tol
//
// Breslow estimator for cumulative hazard at observed event times:
//   dΛ_j = (w_j δ_j) / S0_j ,  where S0_j = Σ_{l in R_j} w_l θ_l
//
// We avoid explicit risk-set enumeration by sorting and cumsum.
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd coxResiduals(
    const Eigen::Ref<const Eigen::VectorXd> &time,
    const Eigen::Ref<const Eigen::VectorXd> &event,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &weights,
    double tol,
    int maxIter
) {

    // ── Drop rows with NaN ─────────────────────────────────────────────
    auto keep = completeRows(X, {&time, &event, &weights});
    const Eigen::Index n = static_cast<Eigen::Index>(keep.size());
    const Eigen::Index p = X.cols();

    Eigen::VectorXd t_ = subsetVec(time, keep);
    Eigen::VectorXd d_ = subsetVec(event, keep);
    Eigen::MatrixXd X_ = subsetRows(X, keep);
    Eigen::VectorXd w_ = subsetVec(weights, keep);

    // ── Input validation ──────────────────────────────────────────────
    validateSubset("coxResiduals", n, p, w_, X_, /*startCol=*/ 0);
    if (d_.sum() < 1.0) throw std::runtime_error("coxResiduals: no observed events (all censored)");

    // ── Sort by descending time (stable, so tied times keep order) ─────
    std::vector<Eigen::Index> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(), [&](Eigen::Index a, Eigen::Index b) {
        return t_[a] > t_[b];
    });

    Eigen::VectorXd ts(n), ds(n), ws(n);
    Eigen::MatrixXd Xs(n, p);
    for (Eigen::Index i = 0; i < n; ++i) {
        ts[i] = t_[ord[i]];
        ds[i] = d_[ord[i]];
        ws[i] = w_[ord[i]];
        Xs.row(i) = X_.row(ord[i]);
    }

    // ── Newton-Raphson on Breslow partial log-likelihood ───────────────
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd eta(n), theta(n); // linear predictor, exp(eta)

    // Scratch for cumulative sums (allocated once)
    Eigen::VectorXd S0(n);    // cumsum of w*theta
    Eigen::MatrixXd S1(n, p); // cumsum of w*theta*x
    // S2 (p×p per subject) not stored; we accumulate Hessian directly.

    double loglik = -std::numeric_limits<double>::infinity();
    const double sqrtTol = std::sqrt(tol);

    for (int iter = 0; iter < maxIter; ++iter) {
        eta.noalias() = Xs * beta;
        for (Eigen::Index i = 0; i < n; ++i) {
            eta[i] = std::clamp(eta[i], -500.0, 500.0);
            theta[i] = std::exp(eta[i]);
        }

        // ── Forward cumulative sums (descending-time order) ──────────────
        S0[0] = ws[0] * theta[0];
        S1.row(0) = ws[0] * theta[0] * Xs.row(0);
        for (Eigen::Index i = 1; i < n; ++i) {
            double wt = ws[i] * theta[i];
            S0[i] = S0[i - 1] + wt;
            S1.row(i) = S1.row(i - 1) + wt * Xs.row(i);
        }

        // ── Compute weighted Breslow partial log-likelihood ──────────────
        // ℓ(β) = Σ_{events j} w_j [η_j − log S0(t_j)]
        double loglikOld = loglik;
        loglik = 0.0;

        // ── Score and Hessian ────────────────────────────────────────────
        Eigen::VectorXd U = Eigen::VectorXd::Zero(p);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(p, p);

        Eigen::MatrixXd S2_cum = Eigen::MatrixXd::Zero(p, p);
        Eigen::Index s2_built_to = -1;

        auto advanceS2 = [&](Eigen::Index target) {
            for (Eigen::Index l = s2_built_to + 1; l <= target; ++l) {
                double wt = ws[l] * theta[l];
                S2_cum.selfadjointView<Eigen::Lower>().rankUpdate(Xs.row(l).transpose(), wt);
            }
            s2_built_to = target;
        };

        Eigen::Index i = 0;
        while (i < n) {
            if (ds[i] == 0.0) {
                ++i;
                continue;
            }

            // Tie group: events sharing the same time.
            Eigen::Index j = i;
            while (j + 1 < n && ts[j + 1] == ts[i] && ds[j + 1] != 0.0)
                ++j;
            // Extend to include censored subjects at the same time.
            Eigen::Index k = j;
            while (k + 1 < n && ts[k + 1] == ts[i])
                ++k;

            double s0 = S0[k];
            double logS0 = std::log(s0);
            Eigen::VectorXd xbar = S1.row(k).transpose() / s0;

            advanceS2(k);
            Eigen::MatrixXd S2_k = S2_cum.selfadjointView<Eigen::Lower>();

            double wdSum = 0.0;
            for (Eigen::Index m = i; m <= j; ++m) {
                if (ds[m] == 0.0) continue;
                wdSum += ws[m];
                // Log-likelihood contribution: w_m * (η_m − log S0)
                loglik += ws[m] * (eta[m] - logS0);
                U.noalias() += ws[m] * (Xs.row(m).transpose() - xbar);
            }

            H.noalias() += wdSum * (S2_k / s0 - xbar * xbar.transpose());

            i = k + 1;
        }

        // ── Convergence check (R coxph.control style) ────────────────────
        // Relative change OR absolute fallback.
        if (iter > 0) {
            double absChange = std::fabs(loglik - loglikOld);
            double relChange = (loglikOld != 0.0) ? absChange / std::fabs(loglikOld) : absChange;
            if (relChange < tol || absChange < sqrtTol) break;
        }

        // ── Newton step ──────────────────────────────────────────────────
        if (p > 0) {
            auto llt = H.selfadjointView<Eigen::Lower>().llt();
            if (llt.info() == Eigen::NumericalIssue)
                throw std::runtime_error("coxResiduals: information matrix is singular at iteration " +
                                         std::to_string(iter) + " (collinearity or complete separation)");
            beta += llt.solve(U);
        }
    }

    // ── Martingale residuals ───────────────────────────────────────────
    // r_i = δ_i − Λ̂_0(t_i) exp(X_i β̂)
    //
    // Breslow cumulative hazard: Λ̂_0(t) = Σ_{t_j ≤ t, δ_j=1} dΛ_j
    //   dΛ_j = (Σ_{m in tie group} w_m δ_m) / S0(t_j)
    //
    // Since subjects are sorted by descending time, we scan right-to-left
    // (ascending time) to accumulate Λ̂_0.

    eta.noalias() = Xs * beta;
    for (Eigen::Index i = 0; i < n; ++i) {
        eta[i] = std::clamp(eta[i], -500.0, 500.0);
        theta[i] = std::exp(eta[i]);
    }

    // Recompute S0 for final beta.
    S0[0] = ws[0] * theta[0];
    for (Eigen::Index i = 1; i < n; ++i)
        S0[i] = S0[i - 1] + ws[i] * theta[i];

    // Assign cumulative hazard Λ̂_0(t_i) to each subject.
    Eigen::VectorXd Lambda0(n); // Λ̂_0 at each subject's time
    Lambda0.setZero();

    // Scan from right (smallest time) to left (largest time).
    double cumHaz = 0.0;
    Eigen::Index i2 = n - 1;
    while (i2 >= 0) {
        // Find the start of this tie group (all with same time).
        Eigen::Index start = i2;
        while (start > 0 && ts[start - 1] == ts[i2])
            --start;

        // Risk set is 0..i2 (all with time ≥ ts[i2] in descending sort =
        // all subjects from index 0 to the last with this time).
        // Actually, for the risk set at time ts[i2], we need everyone with
        // time >= ts[i2].  In descending-sorted order, that's indices
        // 0..k where k is the last index with ts[k] == ts[i2].
        // Since we're scanning from right (smallest time) and ts is
        // descending, i2 is actually the first (rightmost) index with
        // this time, and start would be the leftmost.  All subjects
        // 0..i2 have time >= ts[i2].

        // Sum of weighted events in this tie group.
        double dLambda = 0.0;
        for (Eigen::Index m = start; m <= i2; ++m)
            dLambda += ws[m] * ds[m];
        if (dLambda > 0.0) dLambda /= S0[i2];

        cumHaz += dLambda;

        // Assign to all subjects in this tie group.
        for (Eigen::Index m = start; m <= i2; ++m)
            Lambda0[m] = cumHaz;

        i2 = start - 1;
    }

    // Martingale residuals: r_i = δ_i − Λ̂_0(t_i) θ_i
    Eigen::VectorXd resid(n);
    for (Eigen::Index i = 0; i < n; ++i)
        resid[i] = ds[i] - Lambda0[i] * theta[i];

    // ── Un-sort: map back to the order of kept rows ────────────────────
    Eigen::VectorXd result(n);
    for (Eigen::Index i = 0; i < n; ++i)
        result[ord[i]] = resid[i];

    return result;
}

// ──────────────────────────────────────────────────────────────────────
// § 2  Weighted logistic regression — IRLS
//
// Standard IRLS for logistic regression with case weights.
// Returns response residuals: (y − μ̂)
// which matches residuals(glm(..., family=binomial()), type="response") in R.
//
// Convergence: relative deviance change (matching R glm.control):
//   |dev - devOld| / (|dev| + 0.1) < tol
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd logisticResiduals(
    const Eigen::Ref<const Eigen::VectorXd> &y,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &weights,
    double tol,
    int maxIter
) {

    // ── Drop rows with NaN ─────────────────────────────────────────────
    auto keep = completeRows(X, {&y, &weights});
    const Eigen::Index p = X.cols();
    const Eigen::Index nk = static_cast<Eigen::Index>(keep.size());

    Eigen::VectorXd ys = subsetVec(y, keep);
    Eigen::MatrixXd Xs = subsetRows(X, keep);
    Eigen::VectorXd ws = subsetVec(weights, keep);

    // ── Input validation ──────────────────────────────────────────────
    // startCol=1: skip col 0 (all-ones intercept), check covariate cols only.
    validateSubset("logisticResiduals", nk, p, ws, Xs, /*startCol=*/ 1);
    {
        double nCase = ys.sum();
        if (nCase < 1.0) throw std::runtime_error("logisticResiduals: no cases (y=1) after NaN removal");
        if (static_cast<double>(nk) - nCase < 1.0)
            throw std::runtime_error("logisticResiduals: no controls (y=0) after NaN removal");
    }

    // Helper: compute weighted binomial deviance
    //   D = -2 Σ w_i [y_i log(μ_i) + (1−y_i) log(1−μ_i)]
    auto binomialDeviance = [&](const Eigen::ArrayXd &mu) -> double {
        double dev = 0.0;
        for (Eigen::Index i = 0; i < nk; ++i) {
            double m = mu[i];
            if (ys[i] > 0.5)
                dev -= 2.0 * ws[i] * std::log(m);
            else
                dev -= 2.0 * ws[i] * std::log(1.0 - m);
        }
        return dev;
    };

    // ── IRLS ───────────────────────────────────────────────────────────
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    double deviance = std::numeric_limits<double>::infinity();

    for (int iter = 0; iter < maxIter; ++iter) {
        Eigen::VectorXd eta = Xs * beta;
        Eigen::ArrayXd mu = 1.0 / (1.0 + (-eta.array()).exp());
        mu = mu.max(1e-10).min(1.0 - 1e-10);

        // IRLS variance function: V(μ) = μ(1−μ)
        Eigen::ArrayXd varmu = mu * (1.0 - mu);

        // Combined weight: case weight × variance
        Eigen::ArrayXd W = ws.array() * varmu;

        // Working response: z = η + (y − μ) / V(μ)
        Eigen::VectorXd z = eta.array() + (ys.array() - mu) / varmu;

        // Weighted normal equations: (X^T W X) β = X^T W z
        Eigen::MatrixXd XtW = Xs.transpose() * W.matrix().asDiagonal();
        Eigen::MatrixXd XtWXs = XtW * Xs;
        auto llt = XtWXs.selfadjointView<Eigen::Lower>().llt();
        if (llt.info() == Eigen::NumericalIssue)
            throw std::runtime_error("logisticResiduals: weighted information matrix is singular at iteration " +
                                     std::to_string(iter) + " (collinearity or complete separation)");
        Eigen::VectorXd betaNew = llt.solve(XtW * z);

        // ── Convergence check (R glm.control style) ──────────────────────
        // |dev - devOld| / (|dev| + 0.1) < tol
        double devOld = deviance;
        deviance = binomialDeviance(mu);
        if (iter > 0) {
            double relChange = std::fabs(deviance - devOld) / (std::fabs(deviance) + 0.1);
            if (relChange < tol) {
                beta = betaNew;
                break;
            }
        }
        beta = betaNew;
    }

    // ── Response residuals: y − μ̂ ──────────────────────────────────────
    Eigen::ArrayXd eta = (Xs * beta).array();
    Eigen::ArrayXd mu = 1.0 / (1.0 + (-eta).exp());
    mu = mu.max(1e-10).min(1.0 - 1e-10);

    return (ys.array() - mu).matrix();
}

} // namespace regression
