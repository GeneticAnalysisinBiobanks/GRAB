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
// ══════════════════════════════════════════════════════════════════════
// The terminal K: why it is vectorized, and what that costs numerically
// ══════════════════════════════════════════════════════════════════════
//
// Stage 2 accumulated K in a scalar loop calling `logAlpha` once per subject —
// one `expm1` plus one `log1p` each — on the argument that kFull runs twice per
// marker rather than twice per Newton iteration and so "is not on the hot path".
// Stage 4 measured that argument to be wrong.  At n = 651 outliers on AVX-512:
//
//     k12,   dispatched                    0.75 ns/subject      490 ns/call
//     kFull, scalar expm1/log1p pass       9.19 ns/subject     6033 ns/call
//
// One kFull cost 12x one k12 and 4.7x an entire old hand-written kernel
// evaluation.  Against a Newton schedule of three to four k12 calls per tail,
// the two kFull calls per marker were the LARGER half of all CGF time, not a
// rounding error.  The terminal K is therefore vectorized too, through the
// existing `avx2_log_pd` / `avx512_log_pd`, and `logAlphaFast` below replaces
// the log1p/expm1 spelling in the production path.
//
// ── The numerical argument, and its measurement ──────────────────────
//
// N1 is framed as a RELATIVE error in K, and on that metric
// `log1p(p*expm1(x))` is far better than `log(u + p*e^x)`: near x = 0 the latter
// materializes alpha = 1 + delta and rounding alpha discards delta's low bits.
// But K is consumed ONLY through
//
//     w = sgn(zeta) * sqrt(2 * (zeta*s - K))
//
// and every SPA site gates entry on |z| > spaCutoff (default 2.0), so
// zeta*s - K is bounded away from zero and w is of order 2 or more.  What
// propagates into w is therefore the ABSOLUTE error in K:
// dw = -dK/w, so |dw| <= |dK|/2.  For small delta, log(alpha) ~ delta and
// rounding alpha = 1 + delta to a double costs an absolute error of order eps —
// the same order as the log1p spelling costs.  The expensive spelling buys
// relative accuracy in a quantity nothing consumes relatively.
//
// The scope of that argument, and how Stage 8 kept it true.  The premise is "w
// is of order 2 or more", which the |z| > spaCutoff entry gate delivers — the
// measurement above records min |w| = 1.779 over the production grid — and the
// absolute-error argument for the cheap K holds only while it does.  Perturbing
// K alone moves w by dK/w and r* by roughly dK/w^3, since d/dw[log(v/w)/w] is
// -(1 + log(v/w))/w^2 and log(v/w) -> 0 as w -> 0.  At |w| = 1e-4 that turns an
// absolute error of 1e-11 in K — well inside the budget this section defends —
// into an error of order 10 in r*.  Stage 8 therefore did NOT set
// `spa::kWSingularity` to the 1e-4 the plan proposed: at n = 2e5 the measured
// |dK| of 2.1e-11 puts the abscissa where dK/|w|^3 reaches one at 2.8e-4, ABOVE
// that threshold, so 1e-4 would have left the catastrophic region partly
// unguarded.  The constant is 1e-3, chosen to sit above the catastrophe onset
// for cohorts up to n ~ 1e7 and a decade below the |w| ~ 1e-2 that the
// CLI-enforced minimum spaCutoff of 0.01 makes reachable, and below it the tail
// degrades to Phi(+/-w) with r* never formed.  Nothing in this file is licence
// to evaluate r* at small |w|.  See the kWSingularity comment in util/spa.hpp.
//
// Measured against a long-double reference over the domain the solver visits
// (MAF 5e-4 to 0.5, residual scale 0.2 to 8, |z| 2 to 12, both tails, the
// analytic Gaussian block scaled over six orders of magnitude so that |zeta|
// ranges from 1e-5 to 50, and both mixed-sign and one-sign residual sets):
//
//     spelling of K                    |dK|abs    |dK|rel    |dw|     |dLOG10P|
//     ------------------------------   --------   --------   ------   ---------
//     n = 651
//     log1p(p*expm1(x))   [Stage 2]    2.41e-13   1.62e-14   2.6e-14   1.07e-13
//     log(u + p*e^x)      scalar       2.70e-13   8.61e-14   4.3e-14   1.07e-13
//     log(u + p*e^x)      avx512       4.23e-13   8.79e-14   4.3e-14   1.99e-13
//     log(u + p*e^x)      avx2         4.25e-13   8.85e-14   4.2e-14   1.85e-13
//     n = 5000
//     log1p(p*expm1(x))   [Stage 2]    1.96e-12   4.06e-13   2.4e-13   8.88e-13
//     log(u + p*e^x)      avx512       7.35e-13   3.36e-11   2.9e-13   3.27e-13
//
// The absolute error and the error in w are the same order for every spelling —
// within a factor of two, and at n = 5000 the vectorized spelling is the more
// accurate of the two in absolute terms, because there the error of the
// REDUCTION dominates the error of any individual term.  The induced error in
// -log10(p) never exceeds 1e-12, ten orders of magnitude below the tightest
// budget any stage of this rework was held to.  The relative-error column is the
// one place the log1p spelling remains clearly better, by two orders of
// magnitude at n = 5000; `binom*KFullExact` below keeps that spelling available.
//
// ── What `simd_math.hpp` costs here, and what it does not ────────────
//
// `avx512_log_pd` and `avx2_log_pd` were characterized for this change at
// 1.48 ULP worst case over alpha in [1e-16, 3e307], and 3.7e-16 absolute over
// alpha in [0.05, 20] — so the header's "~1 ULP" claim holds for `log`.  It does
// NOT hold for `exp`, which Stage 4 measured at up to 38.9 ULP, and that is the
// real limit on the vectorized K's relative accuracy near t = 0: delta must be
// recovered as p*(lambda - 1), whose relative error is eps*lambda/(lambda - 1)
// however accurate lambda itself is.  Restoring N1's relative accuracy in a
// vectorized kernel would need a vectorized `expm1`, which would have to live in
// `simd_math.hpp` — validated shared infrastructure this stage may not modify.
// The consequence is recorded rather than worked around: the production kernels
// are relatively accurate in K to the level the vendored vector exp allows, and
// `binom*KFullExact` is the entry point that is relatively accurate in K
// unconditionally.
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
// Each of the three k12 reductions, and each of the three K reductions, ships
// the mandated scalar + AVX2 + AVX-512 triple with a `simdLevel()` dispatch
// site, following `src/spasqr/spasqr.cpp`.  The AVX-512 variants use a masked
// tail; the AVX2 variants a scalar tail.
//
// Masked lanes contribute exactly zero in both halves, but for different
// reasons, and both are load-bearing:
//
//   * k12 — inactive lanes load r = 0, and both K' and K'' carry a factor r.
//   * K   — inactive lanes load r = 0 and, where they are vectors, af = 0 and
//           hap = 0.  With r = 0 the argument is t*r = 0, alpha = u + p = 1
//           EXACTLY for every double p in (0, 1) (the rounding error of
//           fl(1 - p) is at most 2^-54, and 1 + eta rounds to 1.0 for
//           |eta| <= 2^-54), and log(1) = 0 exactly in all three tiers.  For
//           binomHapcount the factor hap = 0 zeroes the lane a second time.
//
// The kFull reductions dispatch their (K', K'') half through the same k12
// kernels rather than recomputing them, so kFull's K' and K'' are bit-identical
// to k12's at the same abscissa.  That identity is asserted by
// tests/spa_cgf_test.cpp and matters to the solver: the terminal K'' is the one
// the tail kernel divides by, and it must be the same number the Newton loop
// converged against.

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
// log(alpha), the vectorizable spelling — the production path
// ══════════════════════════════════════════════════════════════════════
//
// Returns log(alpha) = log((1-p) + p·e^x) formed the way the SIMD tiers form it,
// so the scalar tier differs from them only by the vendored vector exp and log
// and not by the algebra.  See the terminal-K section at the top of this header
// for the measurement that justifies preferring this to `logAlpha`, and
// `binom*KFullExact` for the entry point that keeps `logAlpha`.
//
// Three properties are relied on and are not obvious:
//
//   1. THE CLAMP IS UNDONE FOR K.  `alpha` is built from the clamped argument,
//      because that is what k12 does and what keeps every intermediate normal.
//      Past the clamp the derivatives saturate deliberately (see kTrClamp) but K
//      does not: log(u + p·e^x) -> x + log(p) as x -> +inf, i.e. K keeps growing
//      linearly.  Adding back `x - clamp(x)` recovers that exactly:
//      log(alpha_clamped) + (x - 708) = 708 + log(p + u·e^-708) + x - 708, which
//      differs from the true x + log(p + u·e^-x) by under 1e-307.
//
//   2. THE CORRECTION IS DIRECTIONAL, AND VANISHES AT p == 0.  For x < -708 with
//      u > 0 the true value is log(u + p·e^x) -> log(u), which
//      log(alpha_clamped) already gives to full precision since
//      p·e^-708 <= 3.3e-308 is negligible beside u >= 2^-53.  Adding the deficit
//      there would be wrong.  When u == 0 exactly (p == 1) the true value is x
//      itself, and then the deficit is exactly what recovers it.  And at p == 0
//      alpha is identically 1 and log(alpha) is identically 0, so NO correction
//      applies however large x is.  Hence: the full excess when u == 0, its
//      positive part when 0 < p < 1, and nothing at all when p == 0.
//
//   3. alpha IS NEVER ZERO AND NEVER SUBNORMAL.  alpha = u + p·e^{clamp(x)} with
//      u >= 0 and, for p > 0, p·e^{clamp(x)} >= p·3.3e-308.  u == 0 forces p == 1
//      (for any double p < 1, Sterbenz gives 1 - p >= 2^-53), and then
//      alpha = e^{clamp(x)} >= 3.3e-308, which is normal.  This matters because
//      `avx512_log_pd` decomposes its argument by masking the exponent field and
//      would return a wrong answer, not an infinity, on a subnormal.
inline double logAlphaFast(double x, double p) noexcept {
    const double u  = 1.0 - p;
    const double xc = clampTr(x);
    const double a  = u + p * std::exp(xc);
    const double excess = x - xc;              // 0 strictly inside the clamp
    double add = 0.0;                          // see property 2 above
    if (u == 0.0)                    add = excess;
    else if (p > 0.0 && excess > 0.0) add = excess;
    return std::log(a) + add;
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
// The exact-relative-accuracy terminal entry points
// ══════════════════════════════════════════════════════════════════════
//
// Identical to `binom*KFull` in K' and K'' — the same dispatched k12 kernels —
// but accumulating K through `logAlpha`, i.e. N1's log1p(p*expm1(x)) spelling,
// in a scalar loop.  Relatively accurate in K unconditionally, and about six
// times slower.
//
// Production uses `binom*KFull`.  These exist because the relative-accuracy
// property is worth keeping proved and reachable: `tests/spa_cgf_test.cpp` pins
// N1 on them, and any future consumer that needs K itself rather than
// zeta*s - K (a moment-generating-function evaluation, say, or a diagnostic that
// reports K directly) has an entry point that delivers it.
Cgf012 binomUniformKFullExact(
    double t, const double *resid, int n, double p) noexcept;
Cgf012 binomIndivKFullExact(
    double t, const double *resid, const double *af, int n) noexcept;
Cgf012 binomHapcountKFullExact(
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

// The K half of kFull, as a bare reduction.  binomUniform and binomIndiv return
// sum_i log(alpha_i) with the common factor h == 2 left to the caller, which is
// bit-identical to folding it in because multiplication by two is exact;
// binomHapcount returns sum_i hap_i * log(alpha_i), h varying per subject.
//
// Preconditions match the k12 kernels: binomUniform and binomHapcount require
// 0 < p < 1, their endpoints being resolved by the public entry point.
double binomUniformK0_scalar(double t, const double *resid, int n, double p) noexcept;
double binomIndivK0_scalar(
    double t, const double *resid, const double *af, int n) noexcept;
double binomHapcountK0_scalar(
    double t, const double *resid, const double *hap, int n, double q) noexcept;

#if defined(__x86_64__) || defined(_M_X64)
double binomUniformK0_avx2(double t, const double *resid, int n, double p) noexcept;
double binomUniformK0_avx512(double t, const double *resid, int n, double p) noexcept;
double binomIndivK0_avx2(
    double t, const double *resid, const double *af, int n) noexcept;
double binomIndivK0_avx512(
    double t, const double *resid, const double *af, int n) noexcept;
double binomHapcountK0_avx2(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
double binomHapcountK0_avx512(
    double t, const double *resid, const double *hap, int n, double q) noexcept;
#endif

}  // namespace tier

// Runtime SIMD tier actually selected: 0 = scalar, 1 = AVX2, 2 = AVX-512.
// Mirrors `simdLevel()` without exposing the enum to test translation units.
int simdLevelValue();

}  // namespace spa_cgf
