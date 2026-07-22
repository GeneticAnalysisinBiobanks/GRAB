// spagxe_wald.cpp — Branch-B Wald leg of SPAGxE_CCT (self-contained, per marker).
//
// Fits the full interaction model  trait ~ [covar] + g + g:E  per marker and
// returns the two-sided Wald p-value of the G:E coefficient.  The numeric
// patterns mirror src/util/regression.cpp (OLS / IRLS / Breslow-Cox Newton /
// proportional-odds Fisher scoring) but here the fitters retain the coefficient
// covariance — [information]⁻¹ — instead of a residual.  Every trait extracts
// the (last, last) element of the inverse information (the interaction column is
// appended last) via an SPD solve against the last basis vector.
//
// References for each leg (R/SPAGxECCT.R:620–677, SPAGxE_CCT_one_SNP):
//   Linear      lm  , summary(...)$coefficients["E:g","Pr(>|t|)"]   (Student-t)
//   Logistic    glm , summary(...)$coefficients["E:g","Pr(>|z|)"]   (normal)
//   Cox         coxph, summary(...)$coefficients["E:g","Pr(>|z|)"]  (normal)
//   Ordinal     clm , summary(...)$coefficients["E:g","Pr(>|z|)"]   (normal)
//
// Ordinal uses the observed information (analytic −∂²ℓ/∂θ∂θᵀ at the MLE) by
// default, matching R's ordinal::clm; the BHHH outer-product information is
// retained as an option (OrdinalInfo::BHHH) for diagnostics.  One documented
// deviation from R remains (tolerance parity, not byte-exact — plan rule 5):
//   • Cox ties: this fitter uses the Breslow partial likelihood (as
//     regression::coxResiduals does), whereas R coxph defaults to Efron; the two
//     differ only when event times tie.
//
// A singular information matrix, non-convergence, or a degenerate design returns
// NaN, which the caller folds into the Cauchy combination as a fall-back to the
// SPA p (never as p = 0 — GRAB2 convention; model doc §7 B4/B7).

#include "spagxe/spagxe_wald.hpp"

#include "util/math_helper.hpp" // math::pnorm, math::pt

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace spagxe_wald {

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

// ── Linear: OLS, Var(β̂) = σ̂²(ZᵀZ)⁻¹, Student-t with n−k df ─────────────────
// Z includes the intercept column; the interaction is Z's last column.
double waldLinear(const Eigen::MatrixXd &Z, const Eigen::VectorXd &y) {
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
    return 2.0 * math::pt(-std::abs(t), df); // two-sided Pr(>|t|)
}

// ── Logistic: IRLS to convergence, Var(β̂) = (ZᵀWZ)⁻¹, normal reference ─────
// Z includes the intercept column; the interaction is Z's last column.
double waldLogistic(const Eigen::MatrixXd &Z, const Eigen::VectorXd &y) {
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
    return 2.0 * math::pnorm(-std::abs(zstat));
}

// ── Cox: Breslow partial-likelihood Newton, Var(β̂) = H⁻¹, normal reference ──
// X has no intercept (baseline hazard); the interaction is X's last column.
double waldCox(
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
    return 2.0 * math::pnorm(-std::abs(zstat));
}

// ── Ordinal: proportional-odds Fisher scoring, normal reference ─────────────
// X has no intercept (thresholds are the intercepts); interaction is X's last
// column.  y is integer-coded 0..J−1.  Parameter order θ = [ε(J−1) | β(p)], so
// the interaction coefficient is the final θ element.  The covariance uses the
// observed information (analytic −∂²ℓ/∂θ∂θᵀ at the MLE, matching clm) by default,
// or the BHHH outer-product information when `info == OrdinalInfo::BHHH`.
double waldOrdinal(const Eigen::VectorXi &y, const Eigen::MatrixXd &X,
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
    // return NaN so the CCT degrades to the SPA p rather than emitting a wrong p.
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
        if (grad.norm() > 1e-6) return kNaN; // did not reach the MLE → SPA fall-back
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
    const double zstat = beta[p - 1] / se; // interaction is β's last element
    return 2.0 * math::pnorm(-std::abs(zstat));
}

} // namespace

// ══════════════════════════════════════════════════════════════════════════
// Dispatcher — assemble  M = [covar | g | g∘E], drop NaN rows, fit per trait.
// ══════════════════════════════════════════════════════════════════════════

double waldInteractionPval(
    const WaldData &wd,
    const Eigen::Ref<const Eigen::VectorXd> &g,
    const Eigen::Ref<const Eigen::VectorXd> &E,
    OrdinalInfo ordInfo
) {
    if (wd.trait == TraitType::None) return kNaN;
    const Eigen::Index n = wd.covar.rows();
    const Eigen::Index c = wd.covar.cols();
    if (g.size() != n || E.size() != n) return kNaN;

    // M = [covar | g | g∘E]; the interaction g∘E is the last column.
    Eigen::MatrixXd M(n, c + 2);
    if (c > 0) M.leftCols(c) = wd.covar;
    M.col(c) = g;
    M.col(c + 1) = g.array() * E.array();

    const bool isCox = (wd.trait == TraitType::Cox);

    // Complete-case rows: no NaN in M or in the trait's response(s).
    std::vector<Eigen::Index> keep;
    keep.reserve(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i) {
        bool ok = M.row(i).allFinite();
        if (ok) {
            if (isCox)
                ok = std::isfinite(wd.time[i]) && std::isfinite(wd.event[i]);
            else
                ok = std::isfinite(wd.y[i]);
        }
        if (ok) keep.push_back(i);
    }
    const Eigen::Index m = static_cast<Eigen::Index>(keep.size());
    if (m == 0) return kNaN;

    Eigen::MatrixXd Mf(m, M.cols());
    for (Eigen::Index i = 0; i < m; ++i) Mf.row(i) = M.row(keep[i]);

    switch (wd.trait) {
    case TraitType::Linear:
    case TraitType::Logistic: {
        // Prepend the intercept: Z = [1 | covar | g | g∘E].
        Eigen::MatrixXd Z(m, M.cols() + 1);
        Z.col(0).setOnes();
        Z.rightCols(M.cols()) = Mf;
        Eigen::VectorXd yf(m);
        for (Eigen::Index i = 0; i < m; ++i) yf[i] = wd.y[keep[i]];
        return (wd.trait == TraitType::Linear) ? waldLinear(Z, yf)
                                               : waldLogistic(Z, yf);
    }
    case TraitType::Cox: {
        Eigen::VectorXd tf(m), ef(m);
        for (Eigen::Index i = 0; i < m; ++i) {
            tf[i] = wd.time[keep[i]];
            ef[i] = wd.event[keep[i]];
        }
        return waldCox(tf, ef, Mf);
    }
    case TraitType::Ordinal: {
        Eigen::VectorXi yf(m);
        for (Eigen::Index i = 0; i < m; ++i)
            yf[i] = static_cast<int>(std::lround(wd.y[keep[i]]));
        return waldOrdinal(yf, Mf, ordInfo);
    }
    default:
        return kNaN;
    }
}

} // namespace spagxe_wald
