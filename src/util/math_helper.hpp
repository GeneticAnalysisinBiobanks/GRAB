// math_helper.hpp — Shared math/statistics utilities (pure C++17 / Eigen / Boost)
//
// Provides:
//   § 1  Distribution wrappers    pnorm, qnorm, qchisq, pt (Boost)
//   § 2  Bivariate normal CDF     bvnCdf, pmvnorm2d (Genz quadrature)
//   § 3  Brent root-finding       findRootBrent
//   § 4  Diploid genotype MGF     mG0/mG1/mG2, kG0/kG1/kG2 (scalar)
//   § 5  Logistic regression      logisticRegressionBeta, logisticRegression (IRLS, Eigen)
#pragma once

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <Eigen/Dense>
#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>

namespace math {

// ──────────────────────────────────────────────────────────────────────
// § 1  Distribution wrappers
// ──────────────────────────────────────────────────────────────────────

// Normal CDF: P(X ≤ x) or complementary tail, optionally on log scale.
inline double pnorm(
    double x,
    double mean = 0.0,
    double sd = 1.0,
    bool lower_tail = true,
    bool log_p = false
) {
    if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite(x)) {
        double r = ((x > 0) == lower_tail) ? 1.0 : 0.0;
        return log_p ? std::log(r + 1e-300) : r;
    }
    boost::math::normal_distribution<double> dist(mean, sd);
    double result = lower_tail ? boost::math::cdf(dist, x)
                               : boost::math::cdf(boost::math::complement(dist, x));
    if (log_p) result = std::log(result);
    return result;
}

// Normal quantile (inverse CDF).
inline double qnorm(
    double p,
    double mean = 0.0,
    double sd = 1.0,
    bool lower_tail = true,
    bool log_p = false
) {
    if (log_p) p = std::exp(p);
    if (!lower_tail) p = 1.0 - p;
    p = std::clamp(p, 1e-300, 1.0 - 1e-15);
    boost::math::normal_distribution<double> dist(mean, sd);
    return boost::math::quantile(dist, p);
}

// Chi-squared quantile.
inline double qchisq(
    double p,
    double df,
    bool lower_tail = true,
    bool log_p = false
) {
    if (log_p) p = std::exp(p);
    if (!lower_tail) p = 1.0 - p;
    p = std::clamp(p, 1e-300, 1.0 - 1e-15);
    boost::math::chi_squared_distribution<double> dist(df);
    return boost::math::quantile(dist, p);
}

// Student-t CDF (two-tailed p-value helper).
inline double pt(
    double t,
    double df,
    bool lower_tail = true
) {
    if (std::isnan(t)) return std::numeric_limits<double>::quiet_NaN();
    boost::math::students_t_distribution<double> dist(df);
    return lower_tail ? boost::math::cdf(dist, t) : boost::math::cdf(boost::math::complement(dist, t));
}

// ──────────────────────────────────────────────────────────────────────
// § 2  Bivariate normal CDF  (Genz 2004, 6-point quadrature)
// ──────────────────────────────────────────────────────────────────────

// P(X ≤ dh, Y ≤ dk) where (X,Y) ~ BVN(0, 0, 1, 1, r).
double bvnCdf(
    double dh,
    double dk,
    double r
);

// P(lo1 < X ≤ hi1, lo2 < Y ≤ hi2) where (X,Y) ~ BVN(0, [var1, cov12; cov12, var2]).
double pmvnorm2d(
    double lo1,
    double hi1,
    double lo2,
    double hi2,
    double var1,
    double cov12,
    double var2
);

// ──────────────────────────────────────────────────────────────────────
// § 3  Brent root-finding
// ──────────────────────────────────────────────────────────────────────

// Find x in [a,b] s.t. f(x) ≈ 0 using Brent's method.
// Requires f(a) and f(b) to have opposite signs.
template <typename F> double findRootBrent(
    F &&f,
    double a,
    double b,
    double tol = 1e-8
) {
    double fa = f(a);
    double fb = f(b);
    if (fa * fb > 0.0) throw std::runtime_error("findRootBrent: root not bracketed");

    if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
    }

    double c = a, fc = fa;
    const int max_iter = std::clamp(static_cast<int>(20.0 / tol), 15, 50);
    tol = std::max(tol, 1e-9);
    double last_improvement = std::abs(b - a);
    int stagnant = 0;

    for (int iter = 0; iter < max_iter; ++iter) {
        double gap = std::abs(b - a);
        if (gap < tol || std::abs(fb) < tol * 10.0) return b;

        if (gap >= 0.95 * last_improvement) {
            if (++stagnant > 3) break;
        } else {
            stagnant = 0;
        }
        last_improvement = gap;

        double s;
        if (fa != fc && fb != fc && std::abs(fa - fc) > 1e-15 && std::abs(fb - fc) > 1e-15) {
            double d1 = (fa - fb) * (fa - fc);
            double d2 = (fb - fa) * (fb - fc);
            double d3 = (fc - fa) * (fc - fb);
            if (std::abs(d1) > 1e-12 && std::abs(d2) > 1e-12 && std::abs(d3) > 1e-12)
                s = a * fb * fc / d1 + b * fa * fc / d2 + c * fa * fb / d3;
            else
                s = b - fb * (b - a) / (fb - fa);
        } else {
            s = b - fb * (b - a) / (fb - fa);
        }

        if (s <= std::min(a, b) || s >= std::max(a, b)) s = (a + b) / 2.0;

        double fs = f(s);
        c = b;
        fc = fb;
        if (fa * fs < 0) {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if (std::abs(fa) < std::abs(fb)) {
            std::swap(a, b);
            std::swap(fa, fb);
        }
    }
    return b;
}

// ──────────────────────────────────────────────────────────────────────
// § 4  Scalar diploid genotype MGF/cumulant
//
// Genotype G ~ Binomial(2, MAF).  Moment-generating function M(t):
//   M(t) = [(1-p) + p·e^t]^2   where p = MAF
// Cumulant generating function K(t) = log M(t).
// ──────────────────────────────────────────────────────────────────────

// MGF: M_G^(k)(t) for k = 0, 1, 2.
inline double mG0(
    double t,
    double MAF
) {
    double a = 1.0 - MAF + MAF * std::exp(t);
    return a * a;
}

inline double mG1(
    double t,
    double MAF
) {
    double e = MAF * std::exp(t);
    return 2.0 * e * (1.0 - MAF + e);
}

inline double mG2(
    double t,
    double MAF
) {
    double e = MAF * std::exp(t);
    return 2.0 * e * e + 2.0 * e * (1.0 - MAF + e);
}

// CGF: K^(k)(t) = d^k/dt^k log M(t).
inline double kG0(
    double t,
    double MAF
) {
    return std::log(mG0(t, MAF));
}

inline double kG1(
    double t,
    double MAF
) {
    return mG1(t, MAF) / mG0(t, MAF);
}

inline double kG2(
    double t,
    double MAF
) {
    double m0 = mG0(t, MAF);
    double m1 = mG1(t, MAF);
    double m2 = mG2(t, MAF);
    return (m0 * m2 - m1 * m1) / (m0 * m0);
}

// Fused: compute K0, K1, K2 with a single exp() call per subject.
// ~2x fewer floating-point ops than calling kG0/kG1/kG2 separately.
inline void kG012(
    double t,
    double MAF,
    double &K0,
    double &K1,
    double &K2
) {
    const double e = MAF * std::exp(t);
    const double a = 1.0 - MAF + e; // (1-p) + p*e^t
    const double m0 = a * a;
    const double m1 = 2.0 * e * a;
    const double m2 = 2.0 * e * (e + a);
    K0 = std::log(m0);
    K1 = m1 / m0;
    K2 = (m0 * m2 - m1 * m1) / (m0 * m0);
}

// ──────────────────────────────────────────────────────────────────────
// § 5  Logistic regression (IRLS)  — Eigen implementation
//
// Iteratively Reweighted Least Squares for logistic regression.
// Returns coefficient vector β (including intercept as first element).
// ──────────────────────────────────────────────────────────────────────

// Solve logistic(y | X) → β  using IRLS.
// X: (n × p) covariate matrix (intercept NOT included; added internally).
// y: (n × 1) binary outcome in {0, 1}.
Eigen::VectorXd logisticRegressionBeta(
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &y,
    double tol = 1e-6,
    int maxIter = 100
);

// Logistic regression → predicted allele frequency transform:
//   mu = sigmoid(X_new * beta),  return 1 - sqrt(1 - mu).
Eigen::VectorXd logisticRegression(
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &y
);

// ──────────────────────────────────────────────────────────────────────
// § 6  Nelder-Mead simplex optimiser (n-dimensional, unconstrained)
// ──────────────────────────────────────────────────────────────────────

struct OptimResult {
    std::vector<double> par; // best parameters found
    double value;            // objective at par
    int niter;               // iterations used
};

// Minimise f(x) starting from `init` using the Nelder-Mead simplex method.
// f:       objective  R^n → R
// init:    starting point  (length n)
// tol:     convergence tolerance on the simplex diameter
// maxIter: iteration cap
OptimResult nelderMead(
    std::function<double(const std::vector<double> &)> f,
    const std::vector<double> &init,
    double tol = 1e-8,
    int maxIter = 500
);

// ──────────────────────────────────────────────────────────────────────
// § 7  Bounded 1-D minimiser  (Brent's method for minima)
// ──────────────────────────────────────────────────────────────────────

// Minimise f(x) on [lo, hi] using Brent's parabolic / golden-section method.
template <typename F> double brentMin(
    F &&f,
    double lo,
    double hi,
    double tol = 1e-6,
    int maxIter = 200
) {
    constexpr double golden = 0.3819660112501051; // (3 - sqrt(5)) / 2
    double a = lo, b = hi;
    double x = a + golden * (b - a);
    double w = x, v = x;
    double fx = f(x), fw = fx, fv = fx;
    double d = 0.0, e = 0.0;

    for (int iter = 0; iter < maxIter; ++iter) {
        double midpt = 0.5 * (a + b);
        double tol1 = tol * std::abs(x) + 1e-10;
        double tol2 = 2.0 * tol1;
        if (std::abs(x - midpt) <= tol2 - 0.5 * (b - a)) return x;

        double u;
        if (std::abs(e) > tol1) {
            // Parabolic interpolation
            double r = (x - w) * (fx - fv);
            double q = (x - v) * (fx - fw);
            double p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0)
                p = -p;
            else
                q = -q;
            if (std::abs(p) < std::abs(0.5 * q * e) && p > q * (a - x) && p < q * (b - x)) {
                e = d;
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2) d = (x < midpt) ? tol1 : -tol1;
            } else {
                e = (x < midpt ? b : a) - x;
                d = golden * e;
            }
        } else {
            e = (x < midpt ? b : a) - x;
            d = golden * e;
        }
        u = x + (std::abs(d) >= tol1 ? d : (d > 0 ? tol1 : -tol1));
        double fu = f(u);

        if (fu <= fx) {
            if (u < x)
                b = x;
            else
                a = x;
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        } else {
            if (u < x)
                a = u;
            else
                b = u;
            if (fu <= fw || w == x) {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }
    return x;
}

} // namespace math
