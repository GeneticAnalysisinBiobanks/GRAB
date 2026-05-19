// regression.cpp — Null-model fitting and residual construction
//
// Implements regression::{linearResiduals, logisticResiduals, coxResiduals,
// cumulativeLogitFit}.  Internal helpers (NaN row filtering, row subsetting,
// post-filter validation, logistic CDF) live in the anonymous namespace
// below.  See util/regression.hpp for the public API documentation.
#include "util/regression.hpp"

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// ──────────────────────────────────────────────────────────────────────
// Logistic helpers — numerically stable sigmoid and its derivative.
// ──────────────────────────────────────────────────────────────────────

inline double logistic(double x) {
    if (x >= 0.0) {
        double ez = std::exp(-x);
        return 1.0 / (1.0 + ez);
    }
    double ez = std::exp(x);
    return ez / (1.0 + ez);
}

inline double dlogistic(double x) {
    double p = logistic(x);
    return p * (1.0 - p);
}

inline double clampProb(double p) {
    constexpr double lo = 1e-10;
    constexpr double hi = 1.0 - 1e-10;
    return std::max(lo, std::min(hi, p));
}

// ──────────────────────────────────────────────────────────────────────
// Row filtering / subsetting helpers (moved from src/wtcoxg/regression.cpp).
// ──────────────────────────────────────────────────────────────────────

// Return indices of rows where no element is NaN across X and the supplied
// vectors.
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

Eigen::VectorXi subsetVecXi(
    const Eigen::Ref<const Eigen::VectorXi> &v,
    const std::vector<Eigen::Index> &idx
) {
    Eigen::VectorXi out(static_cast<Eigen::Index>(idx.size()));
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
//   Pass 0 for Cox (no intercept), 1 for logistic/linear (col 0 = intercept).
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

namespace regression {

// ══════════════════════════════════════════════════════════════════════
// § 1  Weighted linear regression  (OLS via weighted normal equations)
// ══════════════════════════════════════════════════════════════════════

Eigen::VectorXd linearResiduals(
    const Eigen::Ref<const Eigen::VectorXd> &y,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &weights,
    double /*tol*/,
    int /*maxIter*/
) {

    // ── Drop rows with NaN ─────────────────────────────────────────────
    auto keep = completeRows(X, {&y, &weights});
    const Eigen::Index nk = static_cast<Eigen::Index>(keep.size());
    const Eigen::Index p = X.cols();

    Eigen::VectorXd ys = subsetVec(y, keep);
    Eigen::MatrixXd Xs = subsetRows(X, keep);
    Eigen::VectorXd ws = subsetVec(weights, keep);

    // ── Input validation ──────────────────────────────────────────────
    // startCol=1: skip col 0 (intercept; allowed to be constant).
    validateSubset("linearResiduals", nk, p, ws, Xs, /*startCol=*/ 1);

    // ── Weighted normal equations:  (XᵀW X) β = XᵀW y ──────────────────
    Eigen::MatrixXd XtW = Xs.transpose() * ws.asDiagonal();
    Eigen::MatrixXd XtWX = XtW * Xs;
    auto llt = XtWX.selfadjointView<Eigen::Lower>().llt();
    if (llt.info() == Eigen::NumericalIssue)
        throw std::runtime_error("linearResiduals: weighted information matrix is singular "
                                 "(collinearity or duplicated columns)");
    Eigen::VectorXd beta = llt.solve(XtW * ys);

    // ── Residuals: y − Xβ̂ ──────────────────────────────────────────────
    return (ys - Xs * beta);
}

// ══════════════════════════════════════════════════════════════════════
// § 2  Weighted Cox regression — Breslow partial likelihood
//
// Algorithm overview (single-pass cumulative sums per Newton step):
//   1. Sort subjects by descending survival time.
//   2. Compute θᵢ = exp(Xᵢ β), weighted risk quantities via reverse-
//      cumulative sums over the sorted order.
//   3. Score U = Xᵀ W (δ − dΛ exp(Xβ)) and Hessian H accumulate in one
//      forward pass over events.
//   4. Newton update: β ← β + H⁻¹ U.
//   5. After convergence, martingale residuals = δᵢ − Λ̂₀(tᵢ) θᵢ.
//
// Breslow estimator for cumulative hazard at observed event times:
//     dΛⱼ = (wⱼ δⱼ) / S0ⱼ ,   S0ⱼ = Σ_{l ∈ Rⱼ} wₗ θₗ
// ══════════════════════════════════════════════════════════════════════

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
    Eigen::VectorXd Lambda0(n);
    Lambda0.setZero();

    double cumHaz = 0.0;
    Eigen::Index i2 = n - 1;
    while (i2 >= 0) {
        Eigen::Index start = i2;
        while (start > 0 && ts[start - 1] == ts[i2])
            --start;

        double dLambda = 0.0;
        for (Eigen::Index m = start; m <= i2; ++m)
            dLambda += ws[m] * ds[m];
        if (dLambda > 0.0) dLambda /= S0[i2];

        cumHaz += dLambda;

        for (Eigen::Index m = start; m <= i2; ++m)
            Lambda0[m] = cumHaz;

        i2 = start - 1;
    }

    Eigen::VectorXd resid(n);
    for (Eigen::Index i = 0; i < n; ++i)
        resid[i] = ds[i] - Lambda0[i] * theta[i];

    // ── Un-sort: map back to the order of kept rows ────────────────────
    Eigen::VectorXd result(n);
    for (Eigen::Index i = 0; i < n; ++i)
        result[ord[i]] = resid[i];

    return result;
}

// ══════════════════════════════════════════════════════════════════════
// § 3  Weighted logistic regression — IRLS
// ══════════════════════════════════════════════════════════════════════

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

    // Helper: weighted binomial deviance.
    //   D = −2 Σ wᵢ [yᵢ log(μᵢ) + (1 − yᵢ) log(1 − μᵢ)]
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

        Eigen::ArrayXd varmu = mu * (1.0 - mu);
        Eigen::ArrayXd W = ws.array() * varmu;
        Eigen::VectorXd z = eta.array() + (ys.array() - mu) / varmu;

        Eigen::MatrixXd XtW = Xs.transpose() * W.matrix().asDiagonal();
        Eigen::MatrixXd XtWXs = XtW * Xs;
        auto llt = XtWXs.selfadjointView<Eigen::Lower>().llt();
        if (llt.info() == Eigen::NumericalIssue)
            throw std::runtime_error("logisticResiduals: weighted information matrix is singular at iteration " +
                                     std::to_string(iter) + " (collinearity or complete separation)");
        Eigen::VectorXd betaNew = llt.solve(XtW * z);

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

// ══════════════════════════════════════════════════════════════════════
// § 4  Fixed-effects cumulative-logit fit (proportional-odds model)
//
// Newton / Fisher-scoring fit of   logit P(Y ≤ j | X) = εⱼ − Xβ.
// After convergence the function reconstructs μ̂ᵢⱼ and iR̂ᵢⱼ and aggregates
// the per-subject working residual.  See header for the residual definition
// and the mean-zero convention.
// ══════════════════════════════════════════════════════════════════════

CumulativeLogitFitResult cumulativeLogitFit(
    const Eigen::Ref<const Eigen::VectorXi> &y,
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    double tol,
    int maxIter
) {

    if (y.size() != X.rows())
        throw std::runtime_error("cumulativeLogitFit: y and X must have the same number of rows");

    // ── Drop rows whose X entries contain NaN (y is integer-coded) ─────
    auto keep = completeRows(X, {});
    const Eigen::Index nk = static_cast<Eigen::Index>(keep.size());
    const Eigen::Index p = X.cols();

    Eigen::VectorXi ys = subsetVecXi(y, keep);
    Eigen::MatrixXd Xs = subsetRows(X, keep);

    if (nk == 0) throw std::runtime_error("cumulativeLogitFit: no complete cases after NaN removal");
    if (nk <= p)
        throw std::runtime_error("cumulativeLogitFit: n (" + std::to_string(nk) + ") must exceed p (" +
                                 std::to_string(p) + ")");

    const int n = static_cast<int>(nk);
    const int yMin = ys.minCoeff();
    const int yMax = ys.maxCoeff();
    if (yMin < 0) throw std::runtime_error("cumulativeLogitFit: y must be non-negative integers");
    if (yMax < 1)
        throw std::runtime_error("cumulativeLogitFit: y must contain at least two distinct categories");
    const int J = yMax + 1;
    const int Jm1 = J - 1;

    // ── Initial values  (empirical cumulative log-odds for ε; β = 0) ───
    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd eps(Jm1);

    Eigen::VectorXd cumProp = Eigen::VectorXd::Zero(J);
    for (int i = 0; i < n; ++i) cumProp(ys(i)) += 1.0;
    cumProp /= static_cast<double>(n);
    double cs = 0.0;
    for (int j = 0; j < Jm1; ++j) {
        cs += cumProp(j);
        double q = std::min(0.999, std::max(0.001, cs));
        eps(j) = std::log(q / (1.0 - q));
    }

    const int nTheta = Jm1 + static_cast<int>(p);

    for (int it = 0; it < maxIter; ++it) {
        Eigen::VectorXd eta = Xs * beta;
        Eigen::VectorXd grad = Eigen::VectorXd::Zero(nTheta);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nTheta, nTheta);

        for (int i = 0; i < n; ++i) {
            int yi = ys(i);
            double F_yi = (yi < Jm1) ? logistic(eps(yi) - eta(i)) : 1.0;
            double F_yim1 = (yi > 0) ? logistic(eps(yi - 1) - eta(i)) : 0.0;
            double p_yi = clampProb(F_yi - F_yim1);
            double f_yi = (yi < Jm1) ? dlogistic(eps(yi) - eta(i)) : 0.0;
            double f_yim1 = (yi > 0) ? dlogistic(eps(yi - 1) - eta(i)) : 0.0;
            double dLogPdEta = -(f_yi - f_yim1) / p_yi;

            Eigen::VectorXd score = Eigen::VectorXd::Zero(nTheta);
            if (yi < Jm1) score(yi) = f_yi / p_yi;
            if (yi > 0) score(yi - 1) -= f_yim1 / p_yi;
            score.tail(p) = dLogPdEta * Xs.row(i).transpose();
            grad += score;
            H += score * score.transpose();
        }
        H.diagonal().array() += 1e-6;
        Eigen::VectorXd delta = H.ldlt().solve(grad);
        for (int j = 0; j < Jm1; ++j) eps(j) += delta(j);
        beta += delta.tail(p);
        if (delta.norm() < tol) break;
    }
    for (int j = 1; j < Jm1; ++j)
        eps(j) = std::max(eps(j), eps(j - 1) + 0.01);

    // ── Working residual reconstruction ────────────────────────────────
    //   νᵢⱼ = logistic(εⱼ − ηᵢ),     ν_{i,-1} = 0,  ν_{i,J-1} = 1
    //   μᵢⱼ = νᵢⱼ − ν_{i,j-1}
    //   mᵢⱼ = νᵢⱼ + ν_{i,j-1} − 1
    //   iRᵢⱼ = 1 / (mᵢⱼ − m_{i,J-1})   for j = 0,…,J-2
    //   rᵢ   = Σⱼ (1{Yᵢ = j} − μᵢⱼ) / iRᵢⱼ
    Eigen::VectorXd eta = Xs * beta;
    Eigen::VectorXd resid(n);
    for (int i = 0; i < n; ++i) {
        // Forward sweep over j = 0..Jm1-1 to build νⱼ and mⱼ; track m_{i,J-1}
        // for the iR denominator.
        double nu_prev = 0.0;
        double m_top = 0.0;            // m_{i,J-1} = ν_{i,J-1} + ν_{i,J-2} − 1
        double last_nu = 0.0;
        std::vector<double> mu(Jm1), m(Jm1);
        for (int j = 0; j < Jm1; ++j) {
            const double nu = logistic(eps(j) - eta(i));
            mu[j] = std::max(nu - nu_prev, 1e-20);
            m[j] = nu + nu_prev - 1.0;
            last_nu = nu;
            nu_prev = nu;
        }
        // ν_{i,J-1} = 1
        m_top = 1.0 + last_nu - 1.0; // = last_nu

        double r_i = 0.0;
        const int yi = ys(i);
        for (int j = 0; j < Jm1; ++j) {
            double denom = m[j] - m_top;
            if (std::abs(denom) < 1e-20) denom = std::copysign(1e-20, denom == 0.0 ? -1.0 : denom);
            const double iR = 1.0 / denom;
            const double y_ij = (yi == j) ? 1.0 : 0.0;
            r_i += (y_ij - mu[j]) / iR;
        }
        resid(i) = r_i;
    }

    // Mean-zero centering: see header §4 commentary.  This restores
    // Σᵢ rᵢ = 0 but is not equivalent to the working-space projection
    // performed by a mixed-model proportional-odds fit.
    resid.array() -= resid.mean();

    CumulativeLogitFitResult out;
    out.beta = std::move(beta);
    out.eps = std::move(eps);
    out.residuals = std::move(resid);
    return out;
}

} // namespace regression
