// spamixlocalp_cgf.hpp — the saddlepoint CGF of SPAmixLocalPlus (tier 3)
//
// Stage 7 of the saddlepoint rework, and the only stage that is a REPAIR
// rather than a migration: SPAmixLocalPlus is under active development, does
// not appear in examples/baseline.sh or examples/tutorial.sh, and none of its
// output is pinned.  It also carried the worst saddlepoint defects in the
// repository (01_findings.md D3, D5, N3, P2, P4), so the old code is deleted
// outright rather than translated.
//
// The statistic is the ancestry-specific score
//
//     S_k = sum_i R_i * G_{i,k}
//
// where G_{i,k} is subject i's count of alternate alleles carried on the
// haplotypes that RFMix assigned to ancestry k, and h_{i,k} in {0, 1, 2} is
// how many haplotypes that is.  Under the null,
//
//     G_{i,k} ~ Binomial(h_{i,k}, q_k)   independently over i,
//
// with ONE allele frequency q_k shared by every subject and a TRIAL COUNT that
// varies per subject.  That is `spa_cgf::binomHapcount`, the tier-2 variant
// built for this method and for no other.  As in SPAmix and WtCoxG the
// subjects are split by the 1.5*IQR rule: the outliers keep the exact binomial
// CGF and the non-outliers contribute a Gaussian block of matching mean and
// variance,
//
//     K_gauss(t) = mean*t + var*t^2/2,  K'_gauss = mean + var*t,  K''_gauss = var.
//
// ══════════════════════════════════════════════════════════════════════
// What this replaces
// ══════════════════════════════════════════════════════════════════════
//
// `kG012Local`, `findRoot`, `RootResult` and `lugannamiRicePval` in
// spamixlocalp.cpp — "Family C" of the three incompatible root-finder families
// 01_findings.md identified.  Point by point:
//
//   * `kG012Local` (spamixlocalp.cpp:876-898) was the one binomial kernel in
//     the repository already free of D1's cancellation — it formed
//     K'' = h*p*e^t*(1-p)/base^2 rather than a difference — but it computed
//     `K0 = h*log(base)` on every Newton iteration and the loop discarded it
//     (P2, :890 computed, :929-930 used), it materialized `base = 1 + q(e^t-1)`
//     rather than using log1p (N1), and it had no SIMD of any kind.  It is
//     replaced by `spa_cgf::binomHapcount`, which supplies all three: the
//     k12/kFull split so the logarithm runs once per tail instead of once per
//     iteration, `logAlpha`'s log1p spelling in the `KFullExact` entry point,
//     and the mandated scalar + AVX2 + AVX-512 triple.
//
//   * `findRoot` (:906-951) restarted a damped Newton iteration from each of
//     {0, -1, 1, -2, 2} until one of them reported convergence.  Its tolerance
//     was an absolute, unscaled `tol = 0.001` on K1total (D7's criticism,
//     which does not list this file), it tested |dt| < tol and returned
//     `conv = true` WITHOUT re-checking K1total — so a stalled iteration was
//     reported as converged — and it maintained no bracket, so nothing
//     guaranteed the point it returned was near a root at all.  The shared
//     solver's criterion is relative and sqrt(K'')-scaled, it applies the step
//     before testing, and it returns only points inside a bracket in which the
//     residual has been observed to change sign.
//
//   * `lugannamiRicePval` (:954-998, the identifier misspells Lugannani)
//     re-derived K and K'' in a second O(nOutlier) scalar pass and then
//     evaluated the Barndorff-Nielsen tail with the defects catalogued under
//     D3 below.
//
// ══════════════════════════════════════════════════════════════════════
// D3 — what each of the three parts turned out to be
// ══════════════════════════════════════════════════════════════════════
//
// D3(a) — MIN_P_VALUE.  `static constexpr double MIN_P_VALUE =
// std::numeric_limits<double>::min()` (:41) is 2.2250738585072014e-308, and it
// was returned from four places: the |w| or |v| below 1e-12 test (:987), the
// non-finite r* test (:990), the `pNorm <= 0` clamp (:1016) and the
// `pval <= 0` clamp (:1059).  At :987 the answer being substituted for is
// approximately ONE HALF — w is the signed square root of 2(zeta*s - K), so
// |w| ~ 0 means zeta ~ 0 means s is at the CGF mean — and 2.2e-308 is a hit at
// any genome-wide threshold ever proposed.  It is the only anti-conservative
// guard in the repository; everything else fails toward 1 or toward NaN.  The
// replacement is `spa::Status::SpaWSingular` with P reported as NA (L2), and the
// clamps at :1016 and :1059 are gone with it: `spa::combineTailsLog` clamps at
// p = 1 and has no lower clamp at all, so a p that genuinely underflows keeps
// its magnitude in LOG10P and appears as an honest 0 in P (L3).
//
//   Reachability, since the audit and its re-verification disagreed.  The
//   audit argued the sibling branch at :894 (`base <= 1e-15` yields
//   K0 = -infinity, hence zeta*s - K0 = +inf, hence w = inf and p =
//   MIN_P_VALUE) was live.  The re-verification showed the argument as stated
//   is wrong — the double below 1.0 is 1 - 1.11e-16, which satisfies
//   (1 - q) <= 1e-15 without q being 1.0 — but that the branch is nevertheless
//   dead here for a different reason: q = dosSum/hapSum with both sums built
//   from integral hapcounts (spamixlocalp.cpp:1496, lanc_io.cpp:923-928), so q
//   is either exactly 1.0 or at most 1 - 1/hapSum, and the interval
//   (1 - 1e-15, 1) is unreachable; and q == 1.0 exactly gives qTerm = 0, hence
//   varS = 0, hence z = 0, which returns on the normal branch before the
//   saddlepoint is entered.  Two routes the audit missed are real, though, and
//   both are closed by construction here rather than by argument:
//   `if (base > 1e-15)` is FALSE when base is NaN, so a NaN residual became
//   K0 = -infinity instead of propagating; and a subject with h == 0 was given
//   K0 = -infinity when the correct value is 0.  `spa_cgf::subjectK12` has no
//   `base`-conditioned branch at all: h == 0 returns exactly zero in all three
//   cumulants because every one of them carries h as a factor, and q == 0 and
//   q == 1 are resolved in closed form at the entry point.  This is asserted
//   by `hapcount_degenerate_cases` in tests/spa_cgf_test.cpp.
//
// D3(b) — the half p-value.  Each tail was added only `if (root.converged)`
// (:1049-1054) and NaN was produced only when BOTH failed (:1056), so a marker
// whose upper tail converged and whose lower tail did not was reported at
// approximately half its correct p-value, with nothing in the output to say
// so.  Halving a p-value is a factor-of-two anti-conservative error applied
// silently and selectively to the hardest markers.  `spa::combineTailsLog`
// returns NaN whenever either tail failed, and names the worse of the two
// statuses in SPA_STATUS.
//
// D3(c) — the Boost policy.  Five sites (:981, :994, :996, :1015 and the
// `norm` object at :980/:992/:1014) called `boost::math::cdf` on a
// default-policy `normal_distribution<>`.  Two consequences beyond the ULP
// difference.  First, the default policy has `promote_double = true`, so Phi
// was evaluated in 80-bit x87 while every other saddlepoint site in the
// repository uses `math::NoPromote` (math_helper.hpp:39-56); the two therefore
// disagreed in the last bits for no reason.  Second, the default
// `domain_error` action is `throw_on_error`, and there is no try/catch
// anywhere in spamixlocalp.cpp while the workers are launched bare
// (`workers.emplace_back(workerFn)`), so a non-finite argument reaching :1015
// would call std::terminate.  A NaN cannot in fact reach it — a NaN residual
// makes varS NaN and `(varS > 0.0) ? ... : 0.0` maps that to z = 0 — but an
// infinite residual can: varS = +inf > 0 and z = (NaN - inf)/inf = NaN.
// Everything here routes through `math::pnorm`, which is `NoPromote` and
// returns NaN for NaN rather than throwing.
//
// D3, a fourth part the audit does not list.  The `temp <= 0 || K2 <= 0`
// fallback (:977-982) returned `Phi_upper(|s| / sqrt(K2))` — which omits the
// mean shift, since K'(0) = q*sum_i h_i r_i is the CGF mean and is not zero —
// and IGNORED its `upperTail` argument, so both tail calls returned the same
// one-sided number and :1049-1054 added it twice.  That is the behaviour at
// exactly zeta = 0, where temp = 0 satisfies `temp <= 0`.  There is no
// fallback here: a non-positive temp or curvature is GUARD_TEMP or GUARD_CURV
// and P is NA.
//
// ══════════════════════════════════════════════════════════════════════
// P4 — the O(N) loop that the batched GEMM had already done
// ══════════════════════════════════════════════════════════════════════
//
// The Gaussian block's mean and variance were accumulated by looping over
// `outlier.posNonOutlier` (:1037-1040), which is ~99 % of the cohort, once per
// (marker x ancestry x phenotype) that entered the SPA branch.  Both totals
// were already available: the caller passes `sMean = q * (H^T R)` and
// `varDiag = q(1-q) * (H^T R^2)`, computed for every (marker, ancestry,
// phenotype) at once by the three fused GEMMs at spamixlocalp.cpp:1549-1551.
// Since `detectOutliers` produces a complete partition, the non-outlier totals
// are the all-subject totals minus the outlier totals, and the outlier totals
// are sums over the gather loop that already runs.  O(N) becomes O(nOutlier).
//
// The subtraction is a genuine cancellation, not a free rewrite: outliers are
// by construction the subjects with the largest |R|, so their share of
// sum_i h_i R_i^2 is far above their share of the sample, and the result will
// not be bit-identical to the loop it replaces.  Nothing pins this method's
// output, and the two forms agree to the working precision of a sum that the
// GEMM computes in a different association order anyway.  `varNorm` is floored
// at zero: it is a variance, the exact value is non-negative because the
// outliers are a subset, and only rounding can drive the difference below zero.
//
// The finding's third P4 sub-claim — that a Boost normal CDF is computed
// unconditionally and discarded — is wrong, as its re-verification records:
// `pNorm` at :1015 IS the returned p-value on the ordinary non-SPA path
// (`return {pNorm, pNorm}` at :1018-1020).  It is wasted only on the rare SPA
// branch.  It is retained, as `spa::normalBranch`, for exactly that reason.
//
// ══════════════════════════════════════════════════════════════════════
// The two tails, and the variance-ratio rescaling
// ══════════════════════════════════════════════════════════════════════
//
// The tails are placed by reflection about the fitted mean, which is what the
// old code's `max(sNew, 2*sMeanNew - sNew)` / `min(...)` spelled out:
//
//     upper = P(S >= sMean + |dev|)      lower = P(S <= sMean - |dev|)
//
// Both s values are rescaled by sqrt(varDiag / varS) first, so that a score
// whose variance was inflated by the phi (relatedness) block is compared
// against the independence CGF this header evaluates.  The CGF itself is NOT
// rescaled, so when varS != varDiag the reflection point is not the CGF's own
// mean K'(0) = sMean.  That is the same construction SPAmixPlus and SPAGxE+
// use (see spamix_cgf.hpp), it is pre-existing, and it is preserved: changing
// it would change the statistic rather than the way the statistic is computed,
// which is outside this stage.  It is also inert on an unrelated cohort, where
// the phi off-diagonal block is empty and varS == varDiag exactly.
//
// For the same reason `spa::SolveOpts::scoreSign` is left at zero: sgn(s)
// carries no information when s is reflected about a non-zero, rescaled mean.
// Nothing is lost — K' is non-decreasing, so the root is unique and the
// solver brackets it.
//
// The initial abscissa is the first-order saddlepoint |dev| / K''(0) =
// |dev| / varDiag, capped at 1.2 exactly as SPAsqr, SPAGRM and SPAmix cap
// theirs.  The cap matters: uncapped, a marker with a very small independence
// variance produces an enormous initial abscissa at which spa_cgf's t*r clamp
// has already saturated K'' to ~1e-307, and the solver's first bracket probe
// (the Newton step from there) overflows.

#pragma once

#include <cmath>
#include <cstdint>

#include "util/spa.hpp"
#include "util/spa_cgf.hpp"

namespace spamixlocalp_cgf {

// ══════════════════════════════════════════════════════════════════════
// Context
// ══════════════════════════════════════════════════════════════════════
//
// `resid` and `hap` are the outlier subset gathered into contiguous storage,
// in the same order; `q` is the ancestry's allele frequency, shared by every
// subject.  `mean` and `var` are the Gaussian block over the non-outliers,
// obtained by subtraction (P4).
struct Context {
    const double *resid = nullptr;   // R_i over the outliers, nOutlier entries
    const double *hap = nullptr;     // h_{i,k} over the same subjects
    int nOutlier = 0;
    double q = 0.0;                  // ancestry allele frequency, in (0, 1)

    double mean = 0.0;               // q * sum_non h_i R_i
    double var = 0.0;                // q(1-q) * sum_non h_i R_i^2
};

// Newton-loop payload: {K'(t) - s, K''(t)}.  No logarithm (P2).
inline spa::K12 k12(double t, const Context &c, double s) noexcept {
    const spa_cgf::Cgf12 d =
        spa_cgf::binomHapcountK12(t, c.resid, c.hap, c.nOutlier, c.q);
    return spa::K12{d.K1 + c.mean + c.var * t - s, d.K2 + c.var};
}

// Terminal payload: {K(t), K'(t) - s, K''(t)}.  Runs once per tail.
//
// The Gaussian block's K is written t * (mean + 0.5 * var * t) so that a large
// |t| multiplies once rather than squaring an already-large quantity, and so
// that the two terms are combined before the binomial K is added.
inline spa::K012 kFull(double t, const Context &c, double s) noexcept {
    const spa_cgf::Cgf012 d =
        spa_cgf::binomHapcountKFull(t, c.resid, c.hap, c.nOutlier, c.q);
    return spa::K012{d.K0 + t * (c.mean + 0.5 * c.var * t),
                     d.K1 + c.mean + c.var * t - s,
                     d.K2 + c.var};
}

// ══════════════════════════════════════════════════════════════════════
// Two-sided saddlepoint p-value
// ══════════════════════════════════════════════════════════════════════
//
// `sMean` is the reflection point and `absDev` the (non-negative) deviation
// from it; `indepVar` is K''(0), used only to size the initial abscissa.
//
// Failure semantics are `spa::assemble`'s: a failure in EITHER tail discards
// the saddlepoint for both, and the row reports the two-sided normal tail at
// `zNorm` under a status naming the substitution (decision D5).  That is
// D3(b)'s repair, now with a named estimator behind it instead of NA.
inline spa::Result twoSidedSpa(
    const Context &c,
    double sMean,
    double absDev,
    double indepVar,
    double zNorm,
    double rtol = 1e-6
) noexcept {
    double zeta0 = 0.0;
    if (indepVar > 0.0) {
        const double e = absDev / indepVar;
        if (std::isfinite(e) && e > 0.0) zeta0 = (e > 1.2) ? 1.2 : e;
    }

    double logp[2];
    spa::Status st[2];

    for (int k = 0; k < 2; ++k) {
        const bool lowerTail = (k == 1);
        const double s = lowerTail ? (sMean - absDev) : (sMean + absDev);

        spa::SolveOpts opt;
        opt.init = lowerTail ? -zeta0 : zeta0;
        opt.rtol = rtol;
        // opt.scoreSign deliberately left at 0; see the header note.

        const spa::Saddle sd = spa::solveSaddlepoint(
            s,
            [&](double t) { return k12(t, c, s); },
            [&](double t) { return kFull(t, c, s); },
            opt);

        spa::Status stLog = spa::Status::SpaOk;
        logp[k] = spa::bnTailLog(sd.zeta, s, sd.K0, sd.K2, lowerTail, stLog);
        st[k] = spa::worseStatus(sd.status, stLog);
    }

    // The two-sided assembly and the D5 fallback when either tail failed —
    // one call, in spa.hpp.  One tail evaluation per side: the linear sibling
    // of `bnTailLog` that used to run beside it is gone with the P column's
    // parallel assembly (log10p_unify Stage 3).
    return spa::assemble(logp[0], logp[1], st[0], st[1], zNorm);
}

// The SPA_STATUS output column: 0 SPA_OK, 1 NORMAL, 2 SPA_W_SINGULAR,
// 3..6 the FALLBACK_* codes, 7 NA_POST_FAIL, 8 NA_NO_TEST.  SPAmixLocalPlus formats its own output
// rather than going through the marker engine, so it COULD emit the token
// `spa::statusName` returns; it emits the numeric enumerator instead so that
// the column means the same thing in every method's output, which is the
// encoding Stages 3 to 6 established for SPACox, SPAGRM, SPAsqr, SPAmix,
// SPAGxE, WtCoxG and LEAF.
inline double statusCode(spa::Status s) noexcept {
    return static_cast<double>(static_cast<uint8_t>(s));
}

}  // namespace spamixlocalp_cgf
