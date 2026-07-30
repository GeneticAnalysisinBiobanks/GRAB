// spa_cgf.hpp — Cancellation-free binomial cumulant generating functions
//
// Stage 2 of the saddlepoint rework.  Three shared CGF variants serve the four
// binomial saddlepoint methods; the SPAGRM family blocks and the SPACox
// empirical table are deliberately NOT expressed here (they are tier-3,
// per-method code).
//
//   binomUniform    G_i ~ Binomial(2, p),   one p for all subjects
//                   consumers: SPAsqr, WtCoxG, LEAF
//   binomIndiv      G_i ~ Binomial(2, p_i), per-individual allele frequency
//                   consumers: SPAmix, SPAGxE, SPAGxEmix
//   binomHapcount   G_i ~ Binomial(h_i, q), per-individual trial count
//                   consumers: SPAmixLocalPlus
//
// ══════════════════════════════════════════════════════════════════════
// The algebra
// ══════════════════════════════════════════════════════════════════════
//
// A subject carries genotype G ~ Binomial(h, p) and enters the statistic with
// residual weight r, so its contribution to the CGF of the score is the CGF of
// r·G.  With
//
//     lambda = e^{t·r},   u = 1 - p,   e = p·lambda,   alpha = u + e
//
// the moment generating function is, by the binomial theorem,
//
//     M(t) = E[e^{t·r·G}] = ((1-p) + p·e^{t·r})^h = alpha^h
//
// and therefore
//
//     K(t)   = h·log(alpha)
//     K'(t)  = h·r·e/alpha             = h·r·pi              , pi := e/alpha
//     K''(t) = h·r^2·e·u/alpha^2       = h·r^2·pi·(1-pi)
//
// ── Derivation of K'' (the point of this file) ────────────────────────
//
//     d(alpha)/dt = p·r·lambda = r·e
//
//     K'  = h · (d alpha/dt)/alpha = h·r·e/alpha
//
//     K'' = h·r · d/dt (e/alpha)
//         = h·r · [ (de/dt)·alpha - e·(d alpha/dt) ] / alpha^2
//         = h·r · [ r·e·alpha - e·r·e ] / alpha^2
//         = h·r^2·e·(alpha - e) / alpha^2
//         = h·r^2·e·u / alpha^2                       since alpha - e = u
//
// The last substitution is the whole point.  `alpha - e` is a difference of two
// positive quantities that coincide to machine precision once t·r grows: for
// t·r >~ 36 with p = 0.5, alpha and e agree in every representable bit, so the
// difference is pure rounding noise.  Replacing it with the loop invariant
// u = 1 - p — which is exact, computed once, and independent of t — removes the
// cancellation entirely.  Every operation in the resulting expression is a
// multiplication or a division of non-negative quantities.
//
// Setting h == 2 recovers the diploid form K'' = 2·r^2·e·u/alpha^2 used by
// binomUniform and binomIndiv, so all three variants share one derivation.
//
// Note that the superficially similar "fix" K'' = h·r^2·pi·(1-pi) with `1-pi`
// evaluated as a subtraction is NOT a fix: 1-pi cancels in exactly the regime
// where alpha-e did.  This file always forms 1-pi as the quotient u/alpha, never
// as a difference.  See 01_findings.md D1.
//
// ── Why K'' is regrouped as h·r^2·pi·(1-pi) rather than h·r^2·e·u/alpha^2 ──
//
// The two are algebraically identical; the grouping below evaluates
//
//     pi  = e * inva ,  omp = u * inva ,  K'' = h*r*r * pi * omp
//
// instead of the left-to-right `h*r*r*e*u*inva*inva`.  The reason is overflow,
// not accuracy.  At the top of the clamped domain e reaches ~3e307, so the
// partial product h*r*r*e overflows to +inf for |r| >~ 2, after which the two
// subsequent multiplications by inva (~3e-308) cannot recover it and K'' is
// reported as +inf instead of ~0.  Because pi and omp are both confined to
// [0, 1] by construction, the regrouped form cannot overflow anywhere in the
// domain, and it makes the absence of cancellation visually evident.
//
// ══════════════════════════════════════════════════════════════════════
// The k12 / kFull split
// ══════════════════════════════════════════════════════════════════════
//
// The Newton loop consumes K' and K'' alone; K is needed only once per tail, at
// the converged root, where it enters w = sgn(zeta)·sqrt(2·(zeta·s - K)).  The
// production kernels compute a logarithm on every iteration and discard it
// (01_findings.md P2), and `avx512_log_pd` costs one `vdivpd` on top of its
// polynomial — the divides, not the polynomials, dominate this kernel.  Hence:
//
//   k12    1 exp, 1 reciprocal, no logarithm.       Newton loop.
//   kFull  k12 + K.                                 Terminal evaluation only.
//
// `expm1` and `log1p` appear ONLY in the K branch, i.e. only in kFull.  This is
// deliberate and load-bearing; see the N1 note on `logAlpha` below.  The SIMD
// layer therefore needs no vectorized expm1/log1p and
// `src/util/simd_math.hpp` is not modified by this stage.
//
// ══════════════════════════════════════════════════════════════════════
// Interface shape
// ══════════════════════════════════════════════════════════════════════
//
// The production entry points are reductions over a subject set: they return
// the SUM over the supplied subjects, because that is what every consumer needs.
// Each caller then adds its own analytic non-outlier block and subtracts the
// score in its own association order — those terms differ per method (SPAsqr's
// variance-ratio rescaling, WtCoxG's deliberately inconsistent (K, K'') pair,
// SPAmixLocalPlus's hapcount-weighted mean) and are tier-3 by the boundary rule.
//
// The per-subject scalar cores are also exposed, inline, because they are the
// shared loop body and because point-by-point comparison against a long-double
// reference is how the reductions are validated.
//
// ══════════════════════════════════════════════════════════════════════
// SIMD
// ══════════════════════════════════════════════════════════════════════
//
// Each of the three k12 reductions ships the mandated scalar + AVX2 + AVX-512
// triple with a `simdLevel()` dispatch site, following
// `src/spasqr/spasqr.cpp`.  The AVX-512 variants use a masked tail; the AVX2
// variants a scalar tail.  Inactive AVX-512 tail lanes load r = 0, and both K'
// and K'' carry a factor r, so masked lanes contribute exactly zero.
//
// The kFull reductions dispatch their (K', K'') half through the same three-tier
// k12 kernels and accumulate K in a scalar pass.  K cannot be vectorized without
// a vectorized expm1/log1p, which this stage is forbidden to add and, per N1,
// does not need: kFull runs twice per marker rather than twice per Newton
// iteration, so its scalar half is not on the hot path.

#pragma once

#include <cmath>

namespace spa_cgf {

// ══════════════════════════════════════════════════════════════════════
// Return types
// ══════════════════════════════════════════════════════════════════════

// Newton-loop payload.  For the reductions, K1 and K2 are sums over subjects.
struct Cgf12 {
    double K1;   // sum_i K'_i(t)
    double K2;   // sum_i K''_i(t)
};

// Terminal payload.  For the reductions, all three are sums over subjects.
struct Cgf012 {
    double K0;   // sum_i K_i(t)
    double K1;   // sum_i K'_i(t)
    double K2;   // sum_i K''_i(t)
};

// ══════════════════════════════════════════════════════════════════════
// Domain clamp
// ══════════════════════════════════════════════════════════════════════
//
// t·r is clamped to +/- kTrClamp before `exp`, in every tier.
//
// Two reasons, both about agreement rather than about accuracy:
//
//   1. Unclamped, `std::exp(800)` returns +inf, giving alpha = +inf,
//      inva = 0 and K' = h·r·(+inf)·0 = NaN.  The vectorized `avx2_exp_pd` /
//      `avx512_exp_pd` clamp their argument internally to +/-709, so they would
//      return a finite result where the scalar tier returns NaN.  Clamping
//      explicitly, at a bound strictly inside the vectorized kernels' own, makes
//      all three tiers agree by construction and removes the NaN.
//
//   2. Inside the bound the clamp is inert, and at the bound the kernel is still
//      exact: at t·r = 708 with p = 0.3 it returns K'' = 1.5435e-307, which
//      agrees with a long-double evaluation to every printed digit.
//
//      PAST the bound the clamp does change the value, and deliberately so.  It
//      returns the derivatives evaluated AT the bound rather than their true,
//      smaller magnitudes: K'' saturates at ~1.5e-307·r^2 whereas the true K''
//      keeps decaying and underflows to zero near t·r = 745.  The alternative at
//      those arguments is not a more accurate number but NaN — unclamped,
//      exp(800) = +inf makes alpha = +inf, inva = 0 and pi = inf·0 = NaN — and a
//      Newton loop that divides by K'' is far better served by a finite positive
//      saturation than by NaN.  In absolute terms the substitution is below
//      1e-307 in a quantity the solver compares against O(1) variances.
//
//      The saturation also keeps K'' STRICTLY positive in the far tails, which
//      the unclamped algebra does not: at t·r = -800 an unclamped exp underflows
//      to exactly zero, giving e = 0 and K'' = 0, and a zero curvature is a
//      division by zero in the Newton step.
//
// exp(708) = 3.0e307 and exp(-708) = 3.3e-308 are both finite and normal, so no
// intermediate of the clamped kernel is subnormal.
//
// K is NOT clamped: `logAlpha` grows linearly in t·r and is overflow-free by
// construction, so the terminal K stays correct past the bound even though the
// derivatives have saturated. That asymmetry is intentional and correct.
//
// A non-finite t·r is outside the contract of these kernels.
inline constexpr double kTrClamp = 708.0;

inline double clampTr(double x) noexcept {
    return (x > kTrClamp) ? kTrClamp : ((x < -kTrClamp) ? -kTrClamp : x);
}

// ══════════════════════════════════════════════════════════════════════
// log(alpha), accurate at t = 0 and overflow-free at both extremes
// ══════════════════════════════════════════════════════════════════════
//
// Returns log(alpha) = log((1-p) + p·e^x).  Multiply by h to obtain K.
//
// ── N1: why expm1/log1p, and where they must NOT be used ─────────────
//
// Near x = 0, alpha = 1 + delta with delta = p·(lambda - 1) far below 1.
// Rounding alpha to a double discards delta's low bits, so log(alpha) inherits a
// relative error of order eps/delta — measured at up to 3.05e-7, nine orders of
// magnitude worse than the K' and K'' the same routine returns.  Forming delta
// directly as p·expm1(x) and taking log1p(delta) never materializes 1 + delta
// and restores full precision at identical cost.
//
// The critical caveat: expm1 is used ONLY for delta, never to recover lambda.
// Reconstructing lambda as 1 + expm1(x) cancels catastrophically for x very
// negative — expm1 returns a value near -1 with absolute accuracy ~eps, so the
// sum carries absolute accuracy ~eps while lambda itself is exponentially
// small.  K' and K'' consequently keep plain `exp`, which is harmless for them
// because alpha = u + e is a sum of two non-negative quantities and never
// cancels.
//
// ── The two branches ─────────────────────────────────────────────────
//
// The log1p form is used over the whole region the saddlepoint solver visits,
// |x| <= kLogAlphaReflect = 700, on both sides of zero.  It is correct
// throughout: at x -> -inf, expm1 saturates at -1 and the result is
// log1p(-p) = log(u), which is the exact limit.
//
// Beyond x = 700 the form cannot be used because expm1(x) overflows at
// x ~ 709.8.  There the exponential is factored out instead:
//
//     log(u + p·e^x) = x + log(p + u·e^{-x})
//
// The parenthesis is a sum of two non-negative quantities, so it never cancels,
// and e^{-x} merely underflows towards zero, leaving x + log(p) — the exact
// limit.  This branch is deliberately NOT written as x + log1p(u·expm1(-x)):
// that spelling reintroduces 1 - u implicitly and loses p entirely whenever
// u rounds to exactly 1.0, i.e. for every p below 2^-53, returning -infinity
// where the true value is finite.
//
// Splitting at 700 rather than at 0 matters for accuracy, not only for
// overflow.  For small positive x the reflected spelling would evaluate
// x + log(...) as a difference of two nearly equal quantities and lose about
// four digits, defeating N1 on the t > 0 side; the log1p branch keeps full
// precision there.
//
// ── The two endpoints ────────────────────────────────────────────────
//
// p == 1 and p == 0 are handled first, and must be: with p == 1 the x <= 0
// branch degenerates to log1p(expm1(x)), which returns -inf once expm1(x) rounds
// to exactly -1, and with p == 0 the x > 0 branch degenerates to
// x + log1p(expm1(-x)), which returns -inf for the same reason.  Both are
// exactly representable instead:
//
//     p == 1  =>  alpha = e^x,  log(alpha) = x     (G == h almost surely)
//     p == 0  =>  alpha = 1,    log(alpha) = 0     (G == 0 almost surely)
//
// Preconditions: p in [0, 1].  Outside that range alpha may vanish by
// cancellation and the result is unspecified.
// Above this argument expm1 overflows and the factored branch takes over.
// expm1 overflows at x ~ 709.78; p·expm1(700) <= 1.02e304 is comfortably finite.
inline constexpr double kLogAlphaReflect = 700.0;

inline double logAlpha(double x, double p) noexcept {
    const double u = 1.0 - p;
    if (u == 0.0) return x;              // p == 1: alpha = e^x exactly
    if (p == 0.0) return 0.0;            // p == 0: alpha = 1 exactly
    if (x <= kLogAlphaReflect) return std::log1p(p * std::expm1(x));
    return x + std::log(p + u * std::exp(-x));
}

// ══════════════════════════════════════════════════════════════════════
// Per-subject scalar cores
// ══════════════════════════════════════════════════════════════════════
//
// One core covers all three variants: binomUniform is h == 2 with p shared,
// binomIndiv is h == 2 with p per-subject, binomHapcount is p == q shared with
// h per-subject.
//
// ── The alpha == 0 guard ─────────────────────────────────────────────
//
// alpha = u + e with u >= 0 and e >= 0, so alpha == 0 requires both to vanish.
// u == 0 means p == 1 exactly (for any double p < 1, Sterbenz gives
// 1 - p >= 2^-53 > 0), and e == 0 then means lambda underflowed.  The exact law
// in that corner is G == h almost surely, hence pi = 1 and pi·(1-pi) = 0.
// Substituting (e, alpha) <- (1, 1) while leaving u == 0 reproduces both values
// exactly, and is the only degenerate input this core can receive.
//
// This is the corner that `kG012Local` currently mishandles: it sets
// K0 = -infinity when its `base <= 1e-15`, which makes zeta·s - K0 = +inf, passes
// the `temp <= 0` guard, and yields w = inf and a p-value of
// numeric_limits<double>::min() — a fabricated genome-wide-significant hit
// (01_findings.md D3).  Here the values are finite and correct instead.
//
// h == 0 needs no special case: K', K'' and K are all proportional to h, and
// every factor they multiply is finite, so all three return exactly zero.
inline Cgf12 subjectK12(double t, double r, double p, double h) noexcept {
    const double u  = 1.0 - p;
    const double tr = clampTr(t * r);
    const double e  = p * std::exp(tr);

    double a  = u + e;                   // sum of non-negatives — never cancels
    double en = e;
    if (a == 0.0) { a = 1.0; en = 1.0; } // p == 1 and lambda underflowed: pi = 1

    const double inva = 1.0 / a;
    const double pi   = en * inva;       // pi   = e/alpha  in [0, 1]
    const double omp  = u  * inva;       // 1-pi = u/alpha  in [0, 1] — a quotient,
                                         //        never the subtraction 1 - pi
    const double hr = h * r;
    return Cgf12{hr * pi, hr * r * pi * omp};
}

// Terminal core: adds K = h·log(alpha).  See `logAlpha` for the N1 treatment.
inline Cgf012 subjectKFull(double t, double r, double p, double h) noexcept {
    const Cgf12 d = subjectK12(t, r, p, h);
    return Cgf012{h * logAlpha(t * r, p), d.K1, d.K2};
}

// ══════════════════════════════════════════════════════════════════════
// Production entry points — sums over the supplied subject set
// ══════════════════════════════════════════════════════════════════════
//
// Each dispatches to the widest SIMD tier the host supports.  `n == 0` returns
// all zeros.  Preconditions: `resid` has n entries; allele frequencies lie in
// [0, 1]; hapcounts are non-negative.

// binomUniform — G_i ~ Binomial(2, p), residual weight resid[i].
//
// The p == 1 and p == 0 endpoints are resolved in closed form before the loop is
// entered, so the vector body may assume 0 < p < 1 and therefore
// alpha >= u > 0.  That removes the per-lane alpha == 0 guard from the hot path
// of the variant that serves SPAsqr, WtCoxG and LEAF.
Cgf12  binomUniformK12(double t, const double *resid, int n, double p) noexcept;
Cgf012 binomUniformKFull(double t, const double *resid, int n, double p) noexcept;

// binomIndiv — G_i ~ Binomial(2, af[i]), residual weight resid[i].
//
// af varies per lane, so the alpha == 0 guard is applied per lane here; it
// cannot be hoisted.
Cgf12  binomIndivK12(double t, const double *resid, const double *af, int n) noexcept;
Cgf012 binomIndivKFull(double t, const double *resid, const double *af, int n) noexcept;

// binomHapcount — G_i ~ Binomial(hap[i], q), residual weight resid[i].
//
// q is shared, so its endpoints are hoisted exactly as in binomUniform; hap
// varies per lane but needs no guard, every cumulant being proportional to it.
Cgf12  binomHapcountK12(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
Cgf012 binomHapcountKFull(
    double t, const double *resid, const double *hap, int n, double q) noexcept;

// ══════════════════════════════════════════════════════════════════════
// Tier-specific variants — exposed for cross-tier equality testing
// ══════════════════════════════════════════════════════════════════════
//
// Production code should call the dispatched entry points above.  These are
// declared so `tests/spa_cgf_test.cpp` can compare the tiers against each other
// and against the long-double reference on every tier the host supports.  A
// variant whose ISA the host lacks must never be called: guard with
// `simdLevelValue()`, following `tests/lanc_simd_test.cpp`.
namespace tier {

Cgf12 binomUniformK12_scalar(double t, const double *resid, int n, double p) noexcept;
Cgf12 binomIndivK12_scalar(
    double t, const double *resid, const double *af, int n) noexcept;
Cgf12 binomHapcountK12_scalar(
    double t, const double *resid, const double *hap, int n, double q) noexcept;

#if defined(__x86_64__) || defined(_M_X64)
Cgf12 binomUniformK12_avx2(double t, const double *resid, int n, double p) noexcept;
Cgf12 binomUniformK12_avx512(double t, const double *resid, int n, double p) noexcept;
Cgf12 binomIndivK12_avx2(
    double t, const double *resid, const double *af, int n) noexcept;
Cgf12 binomIndivK12_avx512(
    double t, const double *resid, const double *af, int n) noexcept;
Cgf12 binomHapcountK12_avx2(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
Cgf12 binomHapcountK12_avx512(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
#endif

}  // namespace tier

// Runtime SIMD tier actually selected: 0 = scalar, 1 = AVX2, 2 = AVX-512.
// Mirrors `simdLevel()` without exposing the enum to test translation units.
int simdLevelValue();

}  // namespace spa_cgf
