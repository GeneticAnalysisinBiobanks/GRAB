// wald.cpp — Standard-model Wald tests for the last coefficient of a design.
//
// Fits a standard model and returns  L = −log10 P  of the two-sided Wald
// p-value of the LAST coefficient of the supplied design.  The numeric patterns mirror
// src/util/regression.cpp (OLS / IRLS / Breslow-Cox Newton / proportional-odds
// Fisher scoring), but here the fitters retain the coefficient covariance —
// [information]⁻¹ — instead of a residual, and extract the (last, last) element
// of the inverse information via an SPD solve against the last basis vector.
//
// References for each leg (R equivalents):
//   Linear      lm  , summary(...)$coefficients[last,"Pr(>|t|)"]   (Student-t)
//   Logistic    glm , summary(...)$coefficients[last,"Pr(>|z|)"]   (normal)
//   Cox         coxph, summary(...)$coefficients[last,"Pr(>|z|)"]  (normal)
//   Ordinal     clm , summary(...)$coefficients[last,"Pr(>|z|)"]   (normal)
//
// Ordinal uses the observed information (analytic −∂²ℓ/∂θ∂θᵀ at the MLE) by
// default, matching R's ordinal::clm; the BHHH outer-product information is
// retained as an option (OrdinalInfo::BHHH) for diagnostics.  One documented
// deviation from R: Cox ties use the Breslow partial likelihood (as
// regression::coxResiduals does), whereas R coxph defaults to Efron; the two
// differ only when event times tie.
//
// A singular information matrix, non-convergence, or a degenerate design
// returns NaN.
//
// The four fitters returned a linear p until log10p_unify Stage 7; both tails
// they invert underflow to exactly zero inside the range a genome-wide scan
// reaches, so the magnitude is now what leaves the file (decision D1).  The
// two tails are evaluated by `math::ptLog` (Student-t, DLMF 8.17.7 below
// Boost's underflow — 01_numerics §3.3) and `math::normalTwoSidedLog`
// (normal); neither has a floor, so neither imposes one on L.

#include "util/wald.hpp"

#include "util/math_helper.hpp" // math::kLn10, math::ptLog, math::normalTwoSidedLog

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace wald {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Numerically stable logistic and its derivative (mirrors regression.cpp).
inline double logistic(double x) {
    if (x >= 0.0) {
        const double e = std::exp(-x);
        return 1.0 / (1.0 + e);
    }
    const double e = std::exp(x);
    return e / (1.0 + e);
}

inline double dlogistic(double x) {
    const double p = logistic(x);
    return p * (1.0 - p);
}

inline double clampProb(double p) {
    return std::max(1e-10, std::min(1.0 - 1e-10, p));
}

// [M⁻¹]_{k−1,k−1} for an SPD information matrix M (k×k), via the Cholesky solve
// M x = e_{k−1}; x_{k−1} equals the requested inverse-diagonal element.  Returns
// a negative sentinel when M is not positive definite (singular / collinear).
double invDiagLast(const Eigen::MatrixXd &M) {
    const Eigen::Index k = M.rows();
    Eigen::LLT<Eigen::MatrixXd> llt(M);
    if (llt.info() != Eigen::Success) return -1.0;
    Eigen::VectorXd e = Eigen::VectorXd::Zero(k);
    e[k - 1] = 1.0;
    const double d = llt.solve(e)[k - 1];
    return std::isfinite(d) ? d : -1.0;
}

// L = −log10 P from ln P.  The sign of zero is normalized: at the tested
// coefficient exactly zero both tails return ln P = +0.0, and −(+0.0)/ln10 is
// −0.0, which `plink2::dtoa_g` prints as "-0".  Nothing else is adjusted — a
// genuinely negative L would be a C2 violation (04_validation §2) and must
// reach the reader rather than be clamped away.
inline double negLog10FromLn(double lnP) noexcept {
    double L = -lnP / math::kLn10;
    if (L == 0.0) L = 0.0;
    return L;
}

// L of the two-sided normal tail 2Φ(−|z|), shared by the three fitters whose
// reference distribution is the normal (logistic / Cox / ordinal).
inline double normalLog10P(double z) noexcept {
    return negLog10FromLn(math::normalTwoSidedLog(z));
}

} // namespace

// ── Linear: OLS, Var(β̂) = σ̂²(ZᵀZ)⁻¹, Student-t with n−k df ─────────────────
// Z includes the intercept column; the tested coefficient is Z's last column.
double lastCoefLinearLog10P(const Eigen::MatrixXd &Z, const Eigen::VectorXd &y) {
    const Eigen::Index n = Z.rows(), k = Z.cols();
    if (n <= k) return kNaN;
    const Eigen::MatrixXd ZtZ = Z.transpose() * Z;
    Eigen::LLT<Eigen::MatrixXd> llt(ZtZ);
    if (llt.info() != Eigen::Success) return kNaN;
    const Eigen::VectorXd beta = llt.solve(Z.transpose() * y);
    const double rss = (y - Z * beta).squaredNorm();
    const double df = static_cast<double>(n - k);
    if (!(df > 0.0)) return kNaN;
    const double sigma2 = rss / df;
    const double invLast = invDiagLast(ZtZ);
    if (!(invLast > 0.0)) return kNaN;
    const double se = std::sqrt(sigma2 * invLast);
    if (!(se > 0.0) || !std::isfinite(se)) return kNaN;
    const double t = beta[k - 1] / se;
    // −log10 of the two-sided Pr(>|t|).  math::ptLog is the same quantity
    // 2·pt(−|t|, df) evaluated as ln I_x(df/2, 1/2), x = df/(df+t²): Boost's
    // regularized incomplete beta wherever that returns a normal double, and
    // the DLMF 8.17.7 hypergeometric series below it.
    return negLog10FromLn(math::ptLog(t, df));
}

// ── Logistic: IRLS to convergence, Var(β̂) = (ZᵀWZ)⁻¹, normal reference ─────
// Z includes the intercept column; the tested coefficient is Z's last column.
double lastCoefLogisticLog10P(const Eigen::MatrixXd &Z, const Eigen::VectorXd &y) {
    const Eigen::Index n = Z.rows(), k = Z.cols();
    if (n <= k) return kNaN;
    const double nCase = y.sum();
    if (nCase < 1.0 || static_cast<double>(n) - nCase < 1.0) return kNaN;

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(k);
    double devOld = std::numeric_limits<double>::infinity();
    for (int it = 0; it < 50; ++it) {
        const Eigen::VectorXd eta = Z * beta;
        Eigen::ArrayXd mu = 1.0 / (1.0 + (-eta.array()).exp());
        mu = mu.max(1e-10).min(1.0 - 1e-10);
        const Eigen::ArrayXd w = mu * (1.0 - mu);
        const Eigen::VectorXd z = eta.array() + (y.array() - mu) / w;
        const Eigen::MatrixXd ZtW = Z.transpose() * w.matrix().asDiagonal();
        const Eigen::MatrixXd info = ZtW * Z;
        Eigen::LLT<Eigen::MatrixXd> llt(info);
        if (llt.info() != Eigen::Success) return kNaN;
        const Eigen::VectorXd betaNew = llt.solve(ZtW * z);
        double dev = 0.0;
        for (Eigen::Index i = 0; i < n; ++i)
            dev -= 2.0 * (y[i] > 0.5 ? std::log(mu[i]) : std::log(1.0 - mu[i]));
        beta = betaNew;
        if (it > 0 && std::fabs(dev - devOld) / (std::fabs(dev) + 0.1) < 1e-9) break;
        devOld = dev;
    }

    // Observed information at the converged β.
    Eigen::VectorXd eta = Z * beta;
    Eigen::ArrayXd mu = 1.0 / (1.0 + (-eta.array()).exp());
    mu = mu.max(1e-10).min(1.0 - 1e-10);
    const Eigen::ArrayXd w = mu * (1.0 - mu);
    const Eigen::MatrixXd info = Z.transpose() * w.matrix().asDiagonal() * Z;
    const double invLast = invDiagLast(info);
    if (!(invLast > 0.0)) return kNaN;
    const double se = std::sqrt(invLast);
    if (!(se > 0.0) || !std::isfinite(se)) return kNaN;
    const double zstat = beta[k - 1] / se;
    return normalLog10P(zstat);
}

// ── Cox: Breslow partial-likelihood Newton, Var(β̂) = H⁻¹, normal reference ──
// X has no intercept (baseline hazard); the tested coefficient is X's last column.
double lastCoefCoxLog10P(
    const Eigen::VectorXd &time,
    const Eigen::VectorXd &event,
    const Eigen::MatrixXd &X
) {
    const Eigen::Index n = X.rows(), p = X.cols();
    if (n <= p) return kNaN;
    if (event.sum() < 1.0) return kNaN;

    // Sort by descending time (stable), so the risk set of an event is a prefix.
    std::vector<Eigen::Index> ord(n);
    std::iota(ord.begin(), ord.end(), Eigen::Index{0});
    std::stable_sort(ord.begin(), ord.end(),
                     [&](Eigen::Index a, Eigen::Index b) { return time[a] > time[b]; });
    Eigen::VectorXd ts(n), ds(n);
    Eigen::MatrixXd Xs(n, p);
    for (Eigen::Index i = 0; i < n; ++i) {
        ts[i] = time[ord[i]];
        ds[i] = event[ord[i]];
        Xs.row(i) = X.row(ord[i]);
    }

    // Partial log-likelihood, score U, and observed information H at β (Breslow).
    auto computeUH = [&](const Eigen::VectorXd &beta, Eigen::VectorXd &U,
                         Eigen::MatrixXd &H) {
        Eigen::VectorXd eta = (Xs * beta).array().min(500.0).max(-500.0).matrix();
        Eigen::VectorXd theta = eta.array().exp().matrix();
        Eigen::VectorXd S0(n);
        Eigen::MatrixXd S1(n, p);
        S0[0] = theta[0];
        S1.row(0) = theta[0] * Xs.row(0);
        for (Eigen::Index i = 1; i < n; ++i) {
            S0[i] = S0[i - 1] + theta[i];
            S1.row(i) = S1.row(i - 1) + theta[i] * Xs.row(i);
        }
        U = Eigen::VectorXd::Zero(p);
        H = Eigen::MatrixXd::Zero(p, p);
        Eigen::MatrixXd S2 = Eigen::MatrixXd::Zero(p, p);
        Eigen::Index built = -1;
        auto advance = [&](Eigen::Index target) {
            for (Eigen::Index l = built + 1; l <= target; ++l)
                S2.selfadjointView<Eigen::Lower>().rankUpdate(Xs.row(l).transpose(),
                                                              theta[l]);
            built = target;
        };
        Eigen::Index i = 0;
        while (i < n) {
            if (ds[i] == 0.0) { ++i; continue; }
            // All observations sharing this death time form the tie group
            // [i..k]; the risk set is the prefix 0..k (descending-time order).
            // The event loop must range over the WHOLE group and skip censored
            // rows — a censored observation may sit BETWEEN two tied events in
            // the stable sort, so bounding by the consecutive-event prefix would
            // drop the later event from the Breslow score, Hessian, and count.
            Eigen::Index k = i;
            while (k + 1 < n && ts[k + 1] == ts[i]) ++k;
            const double s0 = S0[k];
            const Eigen::VectorXd xbar = S1.row(k).transpose() / s0;
            advance(k);
            const Eigen::MatrixXd S2k = S2.selfadjointView<Eigen::Lower>();
            double wd = 0.0;
            for (Eigen::Index m = i; m <= k; ++m) {
                if (ds[m] == 0.0) continue; // censored rows: in the risk set, not deaths
                wd += 1.0;
                U.noalias() += (Xs.row(m).transpose() - xbar);
            }
            H.noalias() += wd * (S2k / s0 - xbar * xbar.transpose());
            i = k + 1;
        }
    };

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd U(p);
    Eigen::MatrixXd H(p, p);
    for (int it = 0; it < 40; ++it) {
        computeUH(beta, U, H);
        Eigen::LLT<Eigen::MatrixXd> llt(H);
        if (llt.info() != Eigen::Success) return kNaN;
        const Eigen::VectorXd step = llt.solve(U);
        beta += step;
        if (step.norm() < 1e-8) break;
    }
    computeUH(beta, U, H); // information at the converged β
    const double invLast = invDiagLast(H);
    if (!(invLast > 0.0)) return kNaN;
    const double se = std::sqrt(invLast);
    if (!(se > 0.0) || !std::isfinite(se)) return kNaN;
    const double zstat = beta[p - 1] / se;
    return normalLog10P(zstat);
}

// ── Ordinal: proportional-odds Fisher scoring, normal reference ─────────────
// X has no intercept (thresholds are the intercepts); tested coefficient is X's
// last column.  y is integer-coded 0..J−1.  Parameter order θ = [ε(J−1) | β(p)],
// so the tested coefficient is the final θ element.  The covariance uses the
// observed information (analytic −∂²ℓ/∂θ∂θᵀ at the MLE, matching clm) by default,
// or the BHHH outer-product information when `info == OrdinalInfo::BHHH`.
double lastCoefOrdinalLog10P(const Eigen::VectorXi &y, const Eigen::MatrixXd &X,
                           OrdinalInfo info) {
    const Eigen::Index n = X.rows(), p = X.cols();
    if (n <= p) return kNaN;
    const int yMin = y.minCoeff();
    const int yMax = y.maxCoeff();
    if (yMin < 0 || yMax < 1) return kNaN;
    const int J = yMax + 1;
    const int Jm1 = J - 1;
    const int nTheta = Jm1 + static_cast<int>(p);

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p);
    Eigen::VectorXd eps(Jm1);
    {
        Eigen::VectorXd cumProp = Eigen::VectorXd::Zero(J);
        for (Eigen::Index i = 0; i < n; ++i) cumProp(y(i)) += 1.0;
        cumProp /= static_cast<double>(n);
        double cs = 0.0;
        for (int j = 0; j < Jm1; ++j) {
            cs += cumProp(j);
            const double q = std::min(0.999, std::max(0.001, cs));
            eps(j) = std::log(q / (1.0 - q));
        }
    }

    // Accumulate score gradient and BHHH information (Σ score·scoreᵀ) at θ.
    auto accumulate = [&](Eigen::VectorXd &grad, Eigen::MatrixXd &H) {
        const Eigen::VectorXd eta = X * beta;
        grad = Eigen::VectorXd::Zero(nTheta);
        H = Eigen::MatrixXd::Zero(nTheta, nTheta);
        for (Eigen::Index i = 0; i < n; ++i) {
            const int yi = y(i);
            const double F_hi = (yi < Jm1) ? logistic(eps(yi) - eta(i)) : 1.0;
            const double F_lo = (yi > 0) ? logistic(eps(yi - 1) - eta(i)) : 0.0;
            const double p_yi = clampProb(F_hi - F_lo);
            const double f_hi = (yi < Jm1) ? dlogistic(eps(yi) - eta(i)) : 0.0;
            const double f_lo = (yi > 0) ? dlogistic(eps(yi - 1) - eta(i)) : 0.0;
            const double dEta = -(f_hi - f_lo) / p_yi;
            Eigen::VectorXd score = Eigen::VectorXd::Zero(nTheta);
            if (yi < Jm1) score(yi) = f_hi / p_yi;
            if (yi > 0) score(yi - 1) -= f_lo / p_yi;
            score.tail(p) = dEta * X.row(i).transpose();
            grad += score;
            H += score * score.transpose();
        }
    };

    // Observed information −∂²ℓ/∂θ∂θᵀ at the current (ε, β), the analytic Hessian
    // that clm uses.  Per subject the log-likelihood is log μ, μ = F(a) − F(b),
    // upper cut a = ε_{yi} − η (F(a)=1 at yi=J−1), lower cut b = ε_{yi−1} − η
    // (F(b)=0 at yi=0); F = logistic, f = F(1−F), f′ = f(1−2F).  The blocks below
    // are the negated second derivatives of that log-likelihood (direct
    // differentiation of the score assembled in `accumulate`).
    auto observedInfo = [&]() -> Eigen::MatrixXd {
        const Eigen::VectorXd eta = X * beta;
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nTheta, nTheta);
        for (Eigen::Index i = 0; i < n; ++i) {
            const int yi = y(i);
            const bool hasHi = (yi < Jm1); // upper threshold ε_{yi} present
            const bool hasLo = (yi > 0);   // lower threshold ε_{yi−1} present
            const double Fa = hasHi ? logistic(eps(yi) - eta(i)) : 1.0;
            const double Fb = hasLo ? logistic(eps(yi - 1) - eta(i)) : 0.0;
            const double fa = hasHi ? Fa * (1.0 - Fa) : 0.0;
            const double fb = hasLo ? Fb * (1.0 - Fb) : 0.0;
            const double fpa = hasHi ? fa * (1.0 - 2.0 * Fa) : 0.0; // f′(a)
            const double fpb = hasLo ? fb * (1.0 - 2.0 * Fb) : 0.0; // f′(b)
            const double mu = clampProb(Fa - Fb);
            const double mu2 = mu * mu;
            const double D = fa - fb;
            const int t1 = yi;     // index of ε_{yi}    (valid iff hasHi)
            const int t0 = yi - 1; // index of ε_{yi−1}  (valid iff hasLo)
            if (hasHi) H(t1, t1) += (fa * fa - fpa * mu) / mu2;
            if (hasLo) H(t0, t0) += (fb * fb + fpb * mu) / mu2;
            if (hasHi && hasLo) {
                const double v = -(fa * fb) / mu2;
                H(t1, t0) += v;
                H(t0, t1) += v;
            }
            if (hasHi) { // ε_{yi} × β block (and its transpose)
                const double c = (fpa * mu - fa * D) / mu2;
                H.row(t1).segment(Jm1, p) += c * X.row(i);
                H.col(t1).segment(Jm1, p) += c * X.row(i).transpose();
            }
            if (hasLo) { // ε_{yi−1} × β block (and its transpose)
                const double c = (-fpb * mu + fb * D) / mu2;
                H.row(t0).segment(Jm1, p) += c * X.row(i);
                H.col(t0).segment(Jm1, p) += c * X.row(i).transpose();
            }
            const double cbb = (D * D - (fpa - fpb) * mu) / mu2; // β × β block
            H.block(Jm1, Jm1, p, p) += cbb * (X.row(i).transpose() * X.row(i));
        }
        return H;
    };

    // Total log-likelihood at the current (ε, β) — drives the line search.
    auto logLik = [&]() {
        const Eigen::VectorXd eta = X * beta;
        double ll = 0.0;
        for (Eigen::Index i = 0; i < n; ++i) {
            const int yi = y(i);
            const double Fhi = (yi < Jm1) ? logistic(eps(yi) - eta(i)) : 1.0;
            const double Flo = (yi > 0) ? logistic(eps(yi - 1) - eta(i)) : 0.0;
            ll += std::log(clampProb(Fhi - Flo));
        }
        return ll;
    };

    // Newton-Raphson on the observed Hessian, with backtracking line search and
    // gradient-based convergence — reaches the same MLE as clm.  A BHHH-only step
    // (Fisher-scoring) can stall on an ill-conditioned marker and stop short of
    // the optimum, so the step matrix is the observed information (with a BHHH
    // fall-back when it is not positive definite).  If the fit fails to converge,
    // return NaN so the caller can degrade gracefully rather than emit a wrong p.
    Eigen::VectorXd grad;
    Eigen::MatrixXd Hbhhh;
    double ll = logLik();
    bool converged = false;
    for (int it = 0; it < 100; ++it) {
        accumulate(grad, Hbhhh);
        if (grad.norm() < 1e-9) { converged = true; break; }
        const Eigen::MatrixXd Hobs = observedInfo();
        Eigen::LLT<Eigen::MatrixXd> llt(Hobs);
        Eigen::VectorXd delta;
        if (llt.info() == Eigen::Success)
            delta = llt.solve(grad);       // Newton step (observed Hessian)
        else
            delta = Hbhhh.ldlt().solve(grad); // BHHH fall-back when Hobs not PD
        if (!delta.allFinite()) return kNaN;
        const Eigen::VectorXd eps0 = eps, beta0 = beta;
        double step = 1.0;
        for (int h = 0; h < 40; ++h) { // backtrack until the log-lik increases
            eps = eps0 + step * delta.head(Jm1);
            beta = beta0 + step * delta.tail(p);
            const double llNew = logLik();
            if (std::isfinite(llNew) && llNew >= ll - 1e-12) { ll = llNew; break; }
            step *= 0.5;
        }
    }
    if (!converged) {
        accumulate(grad, Hbhhh);
        if (grad.norm() > 1e-6) return kNaN; // did not reach the MLE
    }
    // Covariance from the requested information matrix at the MLE.
    Eigen::MatrixXd infoMat;
    if (info == OrdinalInfo::BHHH)
        accumulate(grad, infoMat);      // Σ score·scoreᵀ
    else
        infoMat = observedInfo();       // −∂²ℓ/∂θ∂θᵀ (clm-matching, default)
    const double invLast = invDiagLast(infoMat);
    if (!(invLast > 0.0)) return kNaN;
    const double se = std::sqrt(invLast);
    if (!(se > 0.0) || !std::isfinite(se)) return kNaN;
    const double zstat = beta[p - 1] / se; // tested coefficient is β's last element
    return normalLog10P(zstat);
}

} // namespace wald
