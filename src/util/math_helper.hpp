// math_helper.hpp — Shared math/statistics utilities (pure C++17 / Eigen / Boost)
//
// Provides:
//   § 1  Distribution wrappers    pnorm, pnormLog, qnorm, qchisq, pt (Boost)
//   § 1b Log-domain inversions    zFromNegLog10P, chisq1FromNegLog10P, ptLog
//   § 2  Bivariate normal probability  pmvnorm2dHalfRect / pmvnorm2dHalfRectLog
//   § 3  Brent root-finding       findRootBrent
//   § 4  Diploid genotype MGF     mG0/mG1/mG2, kG0/kG1/kG2 (scalar)
//   § 5  Logistic regression      logisticRegressionBeta, logisticRegression (IRLS, Eigen)
//   § 7b Cauchy combination       cauchyCombine / cauchyCombineLog10
#pragma once

// Some standard library headers transitively pull in <cmath> before we
// reach the line below.  On MinGW, M_PI is only exposed when
// _USE_MATH_DEFINES is set *before* <cmath> is first processed; once
// <cmath>'s include guard fires, defining the macro later has no
// effect.  Therefore we both set the opt-in AND fall back to a
// hand-rolled definition after <cmath>, so any include order leaves
// M_PI defined.
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <Eigen/Dense>
#include <algorithm>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/policies/policy.hpp>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

namespace math {

// Boost's default policy promotes `double` arguments to `long double`
// internally on x86-64 (BOOST_MATH_PROMOTE_DOUBLE_POLICY resolves to true
// wherever the long-double libm functions exist), so the normal CDF and
// quantile are evaluated in 80-bit x87 arithmetic.  The additional
// precision is discarded at the double-precision return type; the only
// observable effect is roughly a fourfold slowdown on the scalar x87 FPU
// relative to the SSE2 double path.  The policy below disables that
// promotion, so pnorm/qnorm compute directly in double.  The results
// agree with the promoted path to within a few ULP.  Across
// examples/baseline.sh this is byte-identical for SPACox, SPAmix, SPAGRM,
// SAGELD, and most SPAsqr columns; for the heaviest consumers (WtCoxG and
// the SPAsqr fit-time residual transform) the few-ULP difference can shift
// the sixth significant figure of a printed p-value, and it can permute
// LEAF's kmeans cluster labels (the combined meta_P result stays
// byte-identical).  Every such difference is numerically identical up to
// GRAB's output precision.
using NoPromote = boost::math::policies::policy<
    boost::math::policies::promote_double<false> >;

// ──────────────────────────────────────────────────────────────────────
// § 1  Distribution wrappers
// ──────────────────────────────────────────────────────────────────────

// Normal CDF: P(X ≤ x) or the complementary tail.
//
// There is deliberately no `log_p` flag.  The only way to serve one from this
// body is to evaluate the CDF on the linear scale and take its logarithm, and
// that inherits the linear scale's underflow: Φ(−x) becomes denormal at
// x ≈ 38.0 and flushes to zero at x ≈ 38.5, so the flag would return −∞ for
// everything below p ≈ 1e-316 while advertising log-domain accuracy
// (spa_unify N2).  Callers that want the logarithm call `pnormLog` below,
// which is a genuine log-domain evaluation.
inline double pnorm(
    double x,
    double mean = 0.0,
    double sd = 1.0,
    bool lower_tail = true
) {
    if (std::isnan(x)) return std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite(x)) return ((x > 0) == lower_tail) ? 1.0 : 0.0;
    boost::math::normal_distribution<double, NoPromote> dist(mean, sd);
    return lower_tail ? boost::math::cdf(dist, x)
                      : boost::math::cdf(boost::math::complement(dist, x));
}

// log(√(2π)), the normalizing constant of the Mills-ratio expansion below.
inline constexpr double kLogSqrt2Pi = 0.91893853320467274178;

// Standardized abscissa below which the linear-scale normal CDF has lost all
// relative accuracy, so the asymptotic expansion is used instead.  Φ(−37) is
// 5.7e-300 — comfortably normal — so the two branches overlap in a regime
// where both are accurate, and the join is continuous to full double
// precision rather than to whatever the denormal range can still express.
inline constexpr double kPnormLogAsymptote = -37.0;

// log Φ((x − mean)/sd), or log of the complementary tail.  A genuine
// log-domain evaluation: the result stays finite and accurate far past the
// point where the probability itself underflows to zero (spa_unify N2, L3).
//
// Above the asymptote the linear CDF is a normal double whose few-ULP
// RELATIVE error becomes a few-ULP ABSOLUTE error in the logarithm, so
// deferring to `pnorm` there is both accurate and exactly what the removed
// `log_p` flag produced — the change is confined to the region where the old
// spelling was returning log(0).
//
// Below the asymptote the Mills-ratio asymptotic expansion is used:
//
//     Φ(x) = φ(x)/(−x) · (1 − 1/x² + 3/x⁴ − 15/x⁶ + 105/x⁸ − …)
//     log Φ(x) = −x²/2 − log(−x) − log(√(2π))
//                + log1p(−1/x² + 3/x⁴ − 15/x⁶ + …)
//
// At |x| = 37 the six retained terms leave a truncation error of order 1e-17
// in the log1p argument, hence 1e-17 absolute in a logarithm whose magnitude
// is ≈ 690 — a relative accuracy near 1e-20.  The series is asymptotic rather
// than convergent, which is exactly why it is used only far out in the tail.
inline double pnormLog(
    double x,
    double mean = 0.0,
    double sd = 1.0,
    bool lower_tail = true
) {
    if (std::isnan(x) || std::isnan(mean) || std::isnan(sd))
        return std::numeric_limits<double>::quiet_NaN();
    if (!std::isfinite(x))
        return ((x > 0) == lower_tail)
                   ? 0.0
                   : -std::numeric_limits<double>::infinity();

    // Standardized abscissa of the tail being asked for: in every case the
    // answer is log Φ(t), so there is one branch rather than two.
    const double t = lower_tail ? (x - mean) / sd : (mean - x) / sd;
    if (std::isnan(t)) return std::numeric_limits<double>::quiet_NaN();
    if (!(t <= kPnormLogAsymptote)) return std::log(pnorm(x, mean, sd, lower_tail));

    const double y = 1.0 / (t * t);
    const double series =
        y * (-1.0 + y * (3.0 + y * (-15.0 + y * (105.0 + y * (-945.0 + y * 10395.0)))));
    return -0.5 * t * t - std::log(-t) - kLogSqrt2Pi + std::log1p(series);
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
    p = std::clamp(p, 1e-300, 1.0 - 1e-15);
    boost::math::normal_distribution<double, NoPromote> dist(mean, sd);
    // Use Boost's complement for the upper tail so very small p (e.g.
    // 1e-300) does not collapse to 1.0 via subtractive cancellation.
    if (lower_tail) return boost::math::quantile(dist, p);
    return boost::math::quantile(boost::math::complement(dist, p));
}

// Two-sided p-value → signed z such that 2*pnorm(-|z|) == p, with the sign
// taken from zNormForSign.  Used to emit the "Z" column that is consistent
// with the (possibly SPA-recalibrated) p-value, alongside the raw normal
// score z "Z_Norm".  qnorm(p/2, lower_tail=false) routes through Boost's
// complement, so a tiny p (e.g. 1e-300) keeps full precision instead of
// collapsing via 1 - p/2.  Degenerate markers propagate NaN from either arg.
inline double zFromPval(double p, double zNormForSign) {
    if (std::isnan(p) || std::isnan(zNormForSign))
        return std::numeric_limits<double>::quiet_NaN();
    const double sign = (zNormForSign >= 0.0) ? 1.0 : -1.0;
    return sign * qnorm(0.5 * p, 0.0, 1.0, /*lower_tail=*/false);
}

// Chi-squared quantile.
inline double qchisq(
    double p,
    double df,
    bool lower_tail = true,
    bool log_p = false
) {
    if (log_p) p = std::exp(p);
    p = std::clamp(p, 1e-300, 1.0 - 1e-15);
    boost::math::chi_squared_distribution<double> dist(df);
    // Use Boost's complement for the upper tail so very small p (e.g.
    // 1e-300) does not collapse to 1.0 via subtractive cancellation.
    if (lower_tail) return boost::math::quantile(dist, p);
    return boost::math::quantile(boost::math::complement(dist, p));
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
// § 1b  Log-domain p-value inversion and evaluation
//
// The log10p_unify project (dev-notes/methods/log10p_unify/) makes
// L = −log10(P) the sole p-value representation.  Every quantity that used
// to be derived from a linear P — the Z column, the df = 1 chi-squared
// weight of a meta-analysis, the Cauchy combination, the Wald leg's
// Student-t tail — is derived from L here instead, so that it keeps
// meaning past P ≈ 1e-308 where the linear representation stops existing.
//
// These functions are ADDED by Stage 1 and are not yet called; the linear
// spellings above (qnorm, qchisq, pt, cauchyCombine) remain in place until
// Stages 4, 5 and 7 switch their call sites over.
// ──────────────────────────────────────────────────────────────────────

// ln 10.  One spelling for the tier.  (spa::detail::kLn10 in src/util/spa.hpp
// carries the same value; the saddlepoint tier should forward to this one when
// it is next edited — see 02_inventory.md §2.2.)
inline constexpr double kLn10 = 2.30258509299404568402;

// Two-sided p-value → signed z, taking the p-value as L = −log10(P).
//
// This is the log-domain inverse of `2·Φ(−|z|)`, i.e. it solves
//
//     ln 2 + ln Φ(−z) = −L · ln 10        for z ≥ 0,
//
// and applies the sign of `zNormForSign`.  It replaces the composition
// `qnorm(0.5·P, upper)` used by `zFromPval`, whose argument is clamped at
// 1e-300 and which therefore SATURATES: every marker with
// L ≥ −log10(2e-300) = 299.698970 came back as |Z| = 37.0470962993612 while
// the adjacent LOG10P column kept rising.  L = 300 is already inside that
// clamp, so the saturation is not an exotic corner (01_numerics §3.1).
//
// Method.  For L below `kZNewtonMinL` the linear argument 0.5·10^(−L) is a
// perfectly ordinary double and Boost's quantile is used unchanged — there is
// no underflow to avoid there, and the Newton criterion below degenerates as
// L → 0 because it is relative to L.  Above it, Newton is run on the
// LOGARITHM of the tail, whose derivative is available in closed form and
// free of underflow:
//
//     f(z)  = ln 2 + ln Φ(−z) + L·ln 10
//     f'(z) = −φ(z)/Φ(−z) = −exp( ln φ(z) − ln Φ(−z) )
//
// started from the analytic inversion of the Mills-ratio asymptote,
// z₀ = √(A − ln A + ln(2/π)) with A = 2·L·ln 10, which is accurate enough
// that two to three iterations converge.  The convergence criterion is
// RELATIVE, |f| ≤ 4·ε·|L·ln 10|, and deliberately contains no absolute
// constant: the same dimensional-soundness requirement the saddlepoint solver
// is held to (CLAUDE.md, src/util/spa.hpp; 04_validation.md §2 invariant 3).
//
// L ≤ 0 (P ≥ 1) returns ±0; L = +∞ returns ±∞; a NaN in either argument
// propagates.
double zFromNegLog10P(double negLog10p, double zNormForSign);

// Upper-tail quantile of the χ²₁ distribution, taking the p-value as
// L = −log10(P).
//
// Analytic, not a numerical inversion: the df = 1 chi-squared upper tail IS
// the two-sided normal tail, P(χ²₁ > q) = 2Φ(−√q), so the quantile is the
// square of the two-sided normal quantile.  This removes the 1e-300 clamp
// inside `qchisq`, which saturated at q = 1373.87 and thereby bounded both
// the per-cluster weight in `metaPvalueScorePool` and the variance WtCoxG
// recovers by inverting a saddlepoint p-value (01_numerics §3.2).
inline double chisq1FromNegLog10P(double negLog10p) {
    const double z = zFromNegLog10P(negLog10p, 1.0);
    return z * z;
}

// Natural logarithm of the TWO-SIDED Student-t tail, P(|T_df| > |t|).
//
// Equal to ln I_x(df/2, 1/2) with x = df/(df + t²).  Boost's regularized
// incomplete beta is used wherever it returns a normal double; below that the
// same quantity is evaluated from the hypergeometric representation
// DLMF 8.17.7,
//
//     I_x(a,b) = x^a (1−x)^b / (a·B(a,b)) · ₂F₁(a+b, 1; a+1; x),
//
// whose series Σ_n [(a+b)_n/(a+1)_n]·xⁿ has all-positive terms and a one
// multiply-divide recurrence.  Note that the Euler form
// I_x(a,b) = x^a/(a·B(a,b))·₂F₁(a, 1−b; a+1; x) is equally correct but has NO
// (1−x)^b factor; combining that ₂F₁ with this prefactor is wrong by exactly
// that factor and only becomes visible at large df, where x → 1
// (01_numerics §3.3).  The stopping rule is relative and the term cap is
// large, because the series ratio tends to x and needs a few hundred terms
// when df is in the thousands.
double ptLog(double t, double df);

// ──────────────────────────────────────────────────────────────────────
// § 2  Bivariate normal probability over a half-infinite rectangle
// ──────────────────────────────────────────────────────────────────────
//
// Compute  P(X ≤ s_hi,  sb_lo ≤ Y ≤ sb_hi)
// where (X, Y) ~ BVN(0, [var1, cov12; cov12, var2]).
// Either or both Y bounds may be ±∞.
//
// Implementation strategy (avoids the subtractive cancellation of the
// 4-corner inclusion-exclusion formula):
//   - both Y bounds infinite  →  marginal Φ(s_hi / √var1)
//   - one Y bound infinite    →  single bivariate-normal CDF evaluation
//                                (Genz 2004)
//   - both Y bounds finite    →  direct 1-D integration of the conditional
//                                tail probability via 20-point Gauss-
//                                Legendre quadrature on [sb_lo, sb_hi]
//                                (no subtraction; ≈ 10⁻¹³ relative error
//                                for |ρ| < 0.925), with the |ρ| → 1
//                                degeneracy handed back to bvnCdf
double pmvnorm2dHalfRect(
    double s_hi,
    double sb_lo,
    double sb_hi,
    double var1,
    double cov12,
    double var2
);

// Natural logarithm of the same probability.
//
// Same arguments, same conventions, same |ρ| > 1 rejection (NaN).  What
// changes is that the answer survives past the point where the probability
// itself underflows: WtCoxG's conditional branch takes −log10 of a ratio of
// two such probabilities, and when the numerator flushes to zero the present
// linear route reports LOG10P_EXT = +Inf rather than a magnitude
// (CLAUDE.md "known gaps"; 01_numerics §3.4).
//
// Both the finite-rectangle case and the half-infinite case are reduced to
// ONE integral, the conditional tail
//
//     P(X ≤ s_hi, a ≤ Y ≤ b) = ∫_a^b φ(u)·Φ((h − ρu)/√(1−ρ²)) du,
//
// whose integrand is positive, so its logarithm is a log-sum-exp over the
// quadrature nodes with no subtraction anywhere.  The half-infinite case is
// the same integral with a = −∞; it is NOT routed through `bvnCdf`, whose
// |r| < 0.925 branch carries an asin(r) factor that is negative for r < 0 and
// therefore is not a positive sum a log-sum-exp could be taken over.
//
// The quadrature is the same 20-point Gauss-Legendre rule the linear routine
// uses, applied on adaptively sized panels rather than once over the whole
// interval, and truncated where the integrand has fallen `kLnDropTrunc`
// below its peak.  Both the panel rule and the truncation act on a LOG
// PROBABILITY and are therefore dimensionless, as invariant 3 of
// 04_validation.md §2 requires.  The panels matter: with a single panel the
// answer is wrong by a factor of e^0.27 as soon as the peak sits on the
// boundary of the integration range with a steep gradient, which is the
// normal situation in the deep tail.
//
// Returns −∞, never +∞, when the probability is genuinely zero or when the
// degenerate |ρ| ≥ 1 − 1e-12 branch underflows; the caller reports that as
// NA with the status code that names it, never as a +Inf LOG10P.
double pmvnorm2dHalfRectLog(
    double s_hi,
    double sb_lo,
    double sb_hi,
    double var1,
    double cov12,
    double var2
);

// Standard bivariate normal CDF, Φ₂(dh, dk; r) = P(X ≤ dh, Y ≤ dk) for
// (X, Y) ~ BVN(0, 0, 1, 1, r).  Genz (2004): Plackett's identity under a
// 6-point Gauss rule for |r| < 0.925, and the asymptotic expansion about the
// degenerate line under a 20-point rule for |r| ≥ 0.925.  Exact at |r| = 1.
//
// Declared here rather than kept file-local so that tests/wtcoxg_cgf_test.cpp
// can pin it against an independent quadrature: the |r| ≥ 0.925 branch was
// mis-transcribed and returned probabilities outside [0, 1] (see the comment
// on the definition in math_helper.cpp), and a defect of that kind must stay
// under test rather than under review.
double bvnCdf(double dh, double dk, double r);

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

// Inverse rank normal transform with Blom plotting position.
//   p_i = (rank_i - 3/8) / (N + 1/4),  Y_int[i] = qnorm(p_i)
// Tied values get the average of their ranks (R `ties.method="average"`).
// Input is assumed to be NaN-free (caller filters); throws on N == 0.
Eigen::VectorXd inverseRankNormal(const Eigen::Ref<const Eigen::VectorXd> &y);

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

// Minimise f(x) starting from `init` using BFGS quasi-Newton with a
// central-difference numerical gradient and Armijo backtracking line
// search.  Designed for smooth, unconstrained problems in low dimension.
//
// Convergence test: |Δf| ≤ reltol · (|f| + reltol)   (same form as R's
// optim(method="BFGS")).
//
// f:       objective  R^n → R
// init:    starting point  (length n)
// reltol:  relative-function convergence tolerance
// maxIter: iteration cap
OptimResult bfgs(
    std::function<double(const std::vector<double> &)> f,
    const std::vector<double> &init,
    double reltol = 1e-10,
    int maxIter = 200
);

// ──────────────────────────────────────────────────────────────────────
// § 7b  Cauchy combination test (CCT)
//
// Combines n independent or dependent p-values into a single p-value via
//   T = (1/n) Σ tan( (0.5 − pᵢ) · π )
// and returns the upper-tail probability of the standard Cauchy at T.
// NaN entries in pvals[] are skipped; if all entries are NaN the result
// is NaN.  Any pᵢ ≤ 0 short-circuits the combined p-value to 0; pᵢ ≥ 1
// is clamped to 0.999.  The tail formulas use 1/(πT) for |T| > 1e15 and
// 1/(πpᵢ) for pᵢ < 1e−15 to avoid overflow in tan().
//
// Reference: Liu & Xie (2020), "Cauchy combination test: a powerful test
// with analytic p-value calculation under arbitrary dependency
// structures", J. Amer. Statist. Assoc. 115 (529), 393–402.
// ──────────────────────────────────────────────────────────────────────

inline double cauchyCombine(
    const double *pvals,
    int n
) {
    if (n <= 0) return std::numeric_limits<double>::quiet_NaN();

    int nValid = 0;
    bool hasZero = false;
    double tStat = 0.0;
    for (int i = 0; i < n; ++i) {
        const double p = pvals[i];
        if (std::isnan(p)) continue;
        ++nValid;
        if (p <= 0.0) { hasZero = true; break; }
        const double pc = (p >= 1.0) ? 0.999 : p;
        tStat += (pc < 1e-15) ? (1.0 / (pc * M_PI))
                              : std::tan((0.5 - pc) * M_PI);
    }

    if (nValid == 0) return std::numeric_limits<double>::quiet_NaN();
    if (hasZero) return 0.0;

    tStat /= static_cast<double>(nValid);
    return (tStat > 1e15) ? (1.0 / tStat) / M_PI
                          : 0.5 - std::atan(tStat) / M_PI;
}

// The same test carried entirely on the −log10 scale: the inputs are
// Lᵢ = −log10(pᵢ) and the result is L_CCT = −log10(p_CCT).
//
// This is not a cosmetic change of units.  The Cauchy statistic's own terms
// are cot(π pᵢ) → 1/(π pᵢ) = 10^Lᵢ/π, so the STATISTIC overflows as soon as
// any pᵢ falls below about 1e-308 — `cauchyCombine` above returns exactly 0
// from `tStat = +Inf` for an input of 1e-320, and 0 is not a p-value.  The
// statistic is therefore never formed: it is carried as T = 10^M·A with
// M = max_i Lᵢ, and the upper tail is inverted from ln T
// (01_numerics §2.2).
//
// Structure worth stating, because it is easy to get backwards: in the tail
//
//     L_CCT ≈ log10( Σ wᵢ · 10^(Lᵢ) )
//
// — a log-sum-exp acting on the Lᵢ and dominated by the LARGEST Lᵢ, which
// for equal weights is the Bonferroni form L_max − log10(n).  It is NOT
// −log10(Σ wᵢ·10^(−Lᵢ)), which is the log of the mean p and is dominated by
// the largest p (01_numerics §2.4).
//
// Semantics carried over from `cauchyCombine` unchanged: NaN entries are
// skipped and n is the count of the survivors; all-NaN returns NaN; a
// non-positive p (here L = +∞) short-circuits to L_CCT = +∞.  The clamp
// p ≥ 1 → 0.999 becomes Lᵢ ← max(Lᵢ, −log10(0.999)), which additionally
// regularizes p ∈ (0.999, 1) where cot(π p) diverges.
double cauchyCombineLog10(
    const double *negLog10p,
    int n
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
