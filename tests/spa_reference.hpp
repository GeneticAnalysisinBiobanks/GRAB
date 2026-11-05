// spa_reference.hpp — Reference implementations of the diploid binomial CGF.
//
// The saddlepoint machinery in GRAB evaluates, for a subject carrying
// genotype G ~ Binomial(2, p) weighted by residual r, the cumulant
// generating function of r·G:
//
//     M(t) = E[e^{t·r·G}] = Σ_{g=0}^{2} C(2,g)·p^g·(1−p)^{2−g}·e^{t·r·g}
//          = ((1−p) + p·e^{t·r})²
//     K(t) = log M(t)
//
// This header provides three evaluations of (K, K′, K″), used together to
// establish the Stage 2 change:
//
//   1. `cgfEnumerate`  — the definition, by explicit summation over the three
//      genotype states, in long double.  Independent of any algebraic
//      simplification; validates that the closed forms below are in fact the
//      binomial CGF.
//
//   2. `cgfDiffForm`   — the algebra used in production today
//      (src/util/math_helper.hpp:297-312, src/spasqr/spasqr.cpp:76-95,
//      src/wtcoxg/wtcoxg.cpp:68-83).  K″ is formed as the difference
//      MGF2/MGF0 − (MGF1/MGF0)².
//
//   3. `cgfStableForm` — the cancellation-free algebra Stage 2 adopts.
//
// The reduction, with λ = e^{t·r}, u = 1 − p, e = p·λ, α = u + e:
//
//     MGF0 = α²                       K   = 2·log α
//     MGF1 = 2·r·e·α                  K′  = 2·r·e/α
//     MGF2 = 2·r²·e·(e + α)           K″  = 2·r²·e·u/α²
//
// The K″ identity is the point of the exercise.  Writing π = e/α, the
// production form computes 2r²π(1+π) − 4r²π², whose two terms both approach
// 4r² as t·r → +∞ while their difference approaches zero.  The stable form
// contains no subtraction at all: e and u are both positive, so every
// operation is a multiplication or a division and the relative error is
// bounded by a small multiple of the unit roundoff over the whole domain.
//
// Note that the naive "fix" K″ = 2r²·π·(1−π) is NOT sufficient: 1−π cancels
// exactly where the original expression does.  The correct substitution is
// 1 − π = u/α, which is what `cgfStableForm` uses.

#pragma once

#include <cmath>
#include <limits>

namespace sparef {

struct Cgf {
    double K0;   // K(t)
    double K1;   // K′(t)
    double K2;   // K″(t)
};

struct CgfL {
    long double K0;
    long double K1;
    long double K2;
};

// ──────────────────────────────────────────────────────────────────────
// 1. Definition — explicit enumeration over genotype states, long double.
//
// Valid only where e^{t·r·2} does not overflow long double (|t·r| < ~5600);
// the callers restrict to the benign region where this path is used.
// ──────────────────────────────────────────────────────────────────────
inline CgfL cgfEnumerate(long double t, long double r, long double p) {
    const long double u = 1.0L - p;
    const long double w[3] = {u * u, 2.0L * p * u, p * p};   // Binomial(2,p) pmf
    const long double g[3] = {0.0L, 1.0L, 2.0L};

    long double m0 = 0.0L, m1 = 0.0L, m2 = 0.0L;
    for (int i = 0; i < 3; ++i) {
        const long double rg = r * g[i];
        const long double ex = w[i] * std::exp(t * rg);
        m0 += ex;                  // E[e^{trG}]
        m1 += rg * ex;             // E[rG·e^{trG}]
        m2 += rg * rg * ex;        // E[(rG)²·e^{trG}]
    }
    const long double d1 = m1 / m0;
    return CgfL{std::log(m0), d1, m2 / m0 - d1 * d1};
}

// ──────────────────────────────────────────────────────────────────────
// 2. Production algebra (double) — reproduces math_helper.hpp:297-312.
//
// Kept bit-for-bit as written there, including the per-subject r² folding
// performed by the callers, so that the characterization test measures the
// error actually present in the shipped binary.
// ──────────────────────────────────────────────────────────────────────
inline Cgf cgfDiffForm(double t, double r, double p) {
    const double e  = p * std::exp(t * r);
    const double a  = 1.0 - p + e;          // (1−p) + p·e^{tr}
    const double m0 = a * a;
    const double m1 = 2.0 * e * a;
    const double m2 = 2.0 * e * (e + a);

    const double K0 = std::log(m0);
    const double K1 = m1 / m0;
    const double K2 = (m0 * m2 - m1 * m1) / (m0 * m0);

    // Callers fold the residual weight in afterwards: K′ ← r·K′, K″ ← r²·K″.
    return Cgf{K0, r * K1, r * r * K2};
}

// ──────────────────────────────────────────────────────────────────────
// 3. Stable algebra (double) — what Stage 2 adopts.
//
// One exp, one reciprocal, no subtraction.  `logNeeded` skips the logarithm
// on the Newton-iteration path, where only K′ and K″ are consumed
// (src/spasqr/spasqr.cpp:292-293, src/wtcoxg/wtcoxg.cpp:464-465,
// src/spamix/common.cpp:81-84).
// ──────────────────────────────────────────────────────────────────────
inline Cgf cgfStableForm(double t, double r, double p, bool logNeeded = true) {
    const double u    = 1.0 - p;            // hoisted out of the loop in the kernels
    const double e    = p * std::exp(t * r);
    const double a    = u + e;
    const double inva = 1.0 / a;

    const double K0 = logNeeded ? 2.0 * std::log(a) : 0.0;
    const double K1 = 2.0 * r * e * inva;
    const double K2 = 2.0 * r * r * e * u * inva * inva;
    return Cgf{K0, K1, K2};
}

// ──────────────────────────────────────────────────────────────────────
// 3b. Stable algebra with an accurate K near t = 0.
//
// `cgfStableForm` above fixes K″ but leaves a second, independent precision
// loss in K itself.  Writing alpha = 1 + delta with delta = p·(lambda − 1),
// the saddlepoint solver visits |t·r| values small enough that delta is far
// below 1; rounding alpha to a double then discards delta's low bits, and
// log(alpha) inherits a relative error of order eps/delta.  At delta ~ 2.5e-4
// that is a relative error near 1e-13 in K — a thousand times worse than the
// K′ and K″ this same routine returns.
//
// The remedy is to never materialize 1 + delta:
//
//     delta = p·expm1(t·r),        alpha = 1 + delta
//     K     = 2·log1p(delta)
//
// expm1 and log1p cost the same as exp and log, so this is free.  It matters
// because K enters w = sgn(zeta)·sqrt(2·(zeta·s − K(zeta))) directly, and w is
// the dominant term of the modified signed root.
//
// Note that expm1 must be used ONLY for delta, never to recover lambda.
// Reconstructing lambda as 1 + expm1(t·r) cancels catastrophically when t·r is
// very negative: expm1 returns a value near −1 with absolute accuracy ~eps, so
// the sum has absolute accuracy ~eps while lambda itself is exponentially
// small, destroying its relative accuracy.  Since K′ and K″ need lambda
// relatively accurate, they keep plain exp — which is harmless for them,
// because alpha = u + e is a sum of two positive quantities and never cancels.
//
// The resulting division of labour matches the cost structure exactly.  The
// SIMD layer (src/util/simd_math.hpp) provides vectorized exp and log but no
// expm1 or log1p, and it does not need them: the Newton loop consumes K′ and
// K″ alone and stays on the vectorized exp path.  The extra expm1 is paid only
// on the single terminal evaluation at the converged root, where K is finally
// required.
// ──────────────────────────────────────────────────────────────────────
inline Cgf cgfStableFormL1P(double t, double r, double p, bool logNeeded = true) {
    const double u    = 1.0 - p;
    const double tr   = t * r;
    const double lam  = std::exp(tr);
    const double e    = p * lam;
    const double a    = u + e;              // sum of positives — never cancels
    const double inva = 1.0 / a;

    // alpha − 1 = p·(lambda − 1), formed without materializing 1 + delta.
    const double K0 = logNeeded ? 2.0 * std::log1p(p * std::expm1(tr)) : 0.0;
    const double K1 = 2.0 * r * e * inva;
    const double K2 = 2.0 * r * r * e * u * inva * inva;
    return Cgf{K0, K1, K2};
}

// ──────────────────────────────────────────────────────────────────────
// Reference: the stable algebra evaluated in long double.
//
// Because the stable form has no cancellation, its long-double evaluation is
// accurate to within a few units in the last place of the long-double format
// (64-bit significand on x86-64), which is ~2048x finer than double.  It is
// therefore an adequate reference for measuring the double-precision error of
// both forms.  `cgfEnumerate` above independently confirms that this algebra
// is the binomial CGF.
// ──────────────────────────────────────────────────────────────────────
inline CgfL cgfReference(long double t, long double r, long double p) {
    const long double u    = 1.0L - p;
    const long double tr   = t * r;
    const long double e    = p * std::exp(tr);
    const long double a    = u + e;
    const long double inva = 1.0L / a;
    // log1p in long double keeps K accurate near t = 0, so the reference does
    // not itself become the limiting factor when the double forms are measured.
    return CgfL{2.0L * std::log1p(p * std::expm1(tr)),
                2.0L * r * e * inva,
                2.0L * r * r * e * u * inva * inva};
}

// True when the platform's long double carries meaningfully more precision
// than double.  On x86-64 Linux/macOS this is the 80-bit x87 format
// (LDBL_MANT_DIG == 64); on arm64 macOS long double is an alias for double,
// in which case the reference is no finer than the value under test and the
// precision assertions must be skipped rather than reported as failures.
inline bool longDoubleIsWider() {
    return sizeof(long double) > sizeof(double) &&
           std::numeric_limits<long double>::digits >
               std::numeric_limits<double>::digits;
}

} // namespace sparef
