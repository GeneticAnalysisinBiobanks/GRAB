// spamixlocalp_cgf_test.cpp — the SPAmixLocalPlus CGF (spa_unify Stage 7).
//
// src/localplus/spamixlocalp_cgf.hpp is header-only; it needs util/spa_cgf.o for
// the dispatched binomial kernels and nothing else from src/.
//
// This stage is a repair rather than a migration — SPAmixLocalPlus is unpinned,
// appearing in neither examples/baseline.sh nor examples/tutorial.sh — so there
// is no regression baseline that could catch a mistake here.  These tests are
// the whole of the automated evidence; the rest is the null-calibration run
// described in 03_stages.md.
//
// What distinguishes this CGF from the other three is the TRIAL COUNT: the
// genotype is G_i ~ Binomial(h_i, q) with h_i in {0, 1, 2} varying per subject
// and one q shared by all of them.  Every failure mode the old `kG012Local`
// had lives in that h, so most of these tests are about it.
//
// In order, the tests establish:
//
//   1. THE TRIAL COUNT IS INSIDE THE KERNEL EXACTLY ONCE.  Sum_i K'_i and
//      Sum_i K''_i match a long-double reference over a mixture of h in
//      {0, 1, 2} and residuals spanning four orders of magnitude, so neither a
//      missing nor a duplicated factor of h or of r survives.  With h == 2
//      everywhere this CGF is indistinguishable from binomUniform's, which is
//      exactly why the fixture must not use a constant h.
//
//   2. h == 0 CONTRIBUTES EXACTLY NOTHING.  01_findings.md D3, as re-verified,
//      records that `kG012Local` assigned K0 = -infinity to any subject whose
//      `base` failed the `> 1e-15` test, and that a subject with h == 0 is one
//      such: its cumulants are all zero, so -infinity was simply wrong, and it
//      then propagated into zeta*s - K0 = +infinity and out through the
//      `temp <= 0` guard.  Here the h == 0 subjects must be removable from the
//      subject set without changing the answer.
//
//   3. A NaN RESIDUAL PROPAGATES.  The same guard was written on `base` rather
//      than on its argument, and `base > 1e-15` is FALSE for NaN, so a NaN
//      residual was silently converted into K0 = -infinity.  It must now come
//      out as NaN.
//
//   4. THE DERIVATIVE BOOKKEEPING IS RIGHT, Gaussian block included.
//
//   5. kFull AND k12 AGREE BIT-FOR-BIT on K' and K''.
//
//   6. THE SCORE IS SUBTRACTED IN K1 ONLY.
//
//   7. THE SOLVER IS EXACT ON A GAUSSIAN CGF, end to end through both tails
//      and the two-sided assembly, against the closed form.
//
//   8. D3(b) — A FAILURE IN ONE TAIL IS NOT HALF A P-VALUE.  The predecessor
//      added each tail only if that tail's root had converged and produced NaN
//      only when both failed, so a one-sided failure was reported at
//      approximately half the correct p-value with no diagnostic.  Constructed
//      here as a CGF whose upper tail is solvable and whose lower tail is not.
//
//   9. THE q ENDPOINTS ARE EXACT AND FINITE.  q == 1 is the corner that drove
//      `kG012Local`'s `base <= 1e-15` branch, and q == 0 is its mirror.

#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <vector>

#include "localplus/spamixlocalp_cgf.hpp"
#include "tinytest.hpp"

namespace {

using spamixlocalp_cgf::Context;
using ld = long double;

// ──────────────────────────────────────────────────────────────────────
// Long-double reference, written from the mathematics
// ──────────────────────────────────────────────────────────────────────
//
// One subject contributes the CGF of r*G with G ~ Binomial(h, q):
//     K   = h*log((1-q) + q*e^{t r})
//     K'  = h*r*q*e^{t r} / alpha
//     K'' = h*r^2*q*(1-q)*e^{t r} / alpha^2
// Both the h and the r factors are written out explicitly so that the test
// states the convention rather than inheriting it from the code under test.
struct Ref {
    ld K0 = 0.0L, K1 = 0.0L, K2 = 0.0L;
};

Ref refSubject(ld t, ld r, ld q, ld h) {
    const ld u = 1.0L - q;
    const ld lam = std::exp(t * r);
    const ld e = q * lam;
    const ld a = u + e;
    Ref o;
    o.K0 = h * std::log(a);
    o.K1 = h * r * e / a;
    o.K2 = h * r * r * e * u / (a * a);
    return o;
}

Ref refSum(double t, const std::vector<double> &r, const std::vector<double> &h,
           double q, double mean, double var) {
    Ref o;
    for (size_t i = 0; i < r.size(); ++i) {
        const Ref s = refSubject(static_cast<ld>(t), static_cast<ld>(r[i]),
                                 static_cast<ld>(q), static_cast<ld>(h[i]));
        o.K0 += s.K0;
        o.K1 += s.K1;
        o.K2 += s.K2;
    }
    const ld tl = static_cast<ld>(t);
    o.K0 += static_cast<ld>(mean) * tl + 0.5L * static_cast<ld>(var) * tl * tl;
    o.K1 += static_cast<ld>(mean) + static_cast<ld>(var) * tl;
    o.K2 += static_cast<ld>(var);
    return o;
}

double relErr(double got, ld want) {
    const ld w = std::fabs(want);
    if (w == 0.0L) return std::fabs(got);
    return static_cast<double>(std::fabs(static_cast<ld>(got) - want) / w);
}

// Residuals spanning four orders of magnitude and hapcounts covering all three
// admissible values.  A constant h would make this CGF indistinguishable from
// binomUniform's and would hide a missing factor of h entirely.
struct Fixture {
    std::vector<double> resid, hap;
    double q = 0.0, mean = 0.0, var = 0.0;
};

Fixture makeFixture(int n, uint32_t seed, double q = 0.23) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uExp(-2.0, 2.0);   // |r| in 1e-2..1e2
    std::uniform_real_distribution<double> uSign(-1.0, 1.0);
    std::uniform_int_distribution<int> uHap(0, 2);
    Fixture f;
    f.resid.resize(static_cast<size_t>(n));
    f.hap.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double mag = std::pow(10.0, uExp(rng));
        f.resid[static_cast<size_t>(i)] = (uSign(rng) < 0.0) ? -mag : mag;
        f.hap[static_cast<size_t>(i)] = static_cast<double>(uHap(rng));
    }
    f.q = q;
    f.mean = 3.5;
    f.var = 12.25;
    return f;
}

Context ctxOf(const Fixture &f) {
    Context c;
    c.resid = f.resid.data();
    c.hap = f.hap.data();
    c.nOutlier = static_cast<int>(f.resid.size());
    c.q = f.q;
    c.mean = f.mean;
    c.var = f.var;
    return c;
}

// ══════════════════════════════════════════════════════════════════════
// 1  The h and r conventions
// ══════════════════════════════════════════════════════════════════════

TEST(hapcount_and_residual_factors_are_inside_the_kernel_exactly_once) {
    const Fixture f = makeFixture(97, 20260731u);
    const Context c = ctxOf(f);

    // t scaled so that t*r stays inside the clamp for the largest |r| ~ 100.
    const double ts[] = {-3.0, -0.4, -1e-7, 0.0, 1e-7, 0.4, 3.0};
    double worst0 = 0.0, worst1 = 0.0, worst2 = 0.0;
    for (double t : ts) {
        const double tt = t / 100.0;
        const Ref want = refSum(tt, f.resid, f.hap, f.q, f.mean, f.var);
        const spa::K12 got12 = spamixlocalp_cgf::k12(tt, c, 0.0);
        const spa::K012 got012 = spamixlocalp_cgf::kFull(tt, c, 0.0);
        worst1 = std::fmax(worst1, relErr(got12.k1, want.K1));
        worst2 = std::fmax(worst2, relErr(got12.k2, want.K2));
        // K is judged in the ABSOLUTE sense in which the tail consumes it; see
        // util/spa_cgf.hpp's terminal-K section for why the production spelling
        // has no relative accuracy near t = 0.
        worst0 = std::fmax(
            worst0, std::fabs(got012.k0 - static_cast<double>(want.K0)));
    }
    // A missing factor of h would show up at the 50 % level, a spurious r^2 at
    // the 1e4 level — not at 1e-15.
    CHECK(worst1 < 1e-12);
    CHECK(worst2 < 1e-12);
    CHECK(worst0 < 1e-11);
    std::printf("    worst err vs long double: K %.3e (abs)  "
                "K' %.3e (rel)  K'' %.3e (rel)\n",
                worst0, worst1, worst2);
}

TEST(the_fixture_could_detect_a_missing_h_or_a_spurious_r_squared) {
    // Guards the guard: the fixture must actually vary h and must keep |r| far
    // from unity, or test 1 proves nothing.
    const Fixture f = makeFixture(97, 20260731u);
    bool saw0 = false, saw1 = false, saw2 = false;
    for (double h : f.hap) {
        if (h == 0.0) saw0 = true;
        if (h == 1.0) saw1 = true;
        if (h == 2.0) saw2 = true;
    }
    CHECK(saw0 && saw1 && saw2);
    double lo = 1e300, hi = 0.0;
    for (double r : f.resid) {
        lo = std::fmin(lo, std::fabs(r));
        hi = std::fmax(hi, std::fabs(r));
    }
    CHECK(lo < 0.05);
    CHECK(hi > 20.0);
}

// ══════════════════════════════════════════════════════════════════════
// 2  h == 0 contributes exactly nothing (D3)
// ══════════════════════════════════════════════════════════════════════

TEST(zero_hapcount_subjects_are_removable) {
    const Fixture f = makeFixture(64, 5150u, 0.4);
    std::vector<double> rSub, hSub;
    for (size_t i = 0; i < f.hap.size(); ++i) {
        if (f.hap[i] != 0.0) {
            rSub.push_back(f.resid[i]);
            hSub.push_back(f.hap[i]);
        }
    }
    CHECK(rSub.size() < f.hap.size());   // the fixture does contain h == 0

    Context full = ctxOf(f);
    Context sub = full;
    sub.resid = rSub.data();
    sub.hap = hSub.data();
    sub.nOutlier = static_cast<int>(rSub.size());

    // Includes t*r far enough negative that the exponential underflows, which
    // is the corner in which `kG012Local` produced -infinity rather than 0.
    double worst = 0.0;
    for (double t : {-4.0, -0.05, 0.0, 0.05, 4.0}) {
        const spa::K012 a = spamixlocalp_cgf::kFull(t, full, 0.0);
        const spa::K012 b = spamixlocalp_cgf::kFull(t, sub, 0.0);
        CHECK(std::isfinite(a.k0) && std::isfinite(a.k1) && std::isfinite(a.k2));
        // Not bit-identical: dropping lanes changes how many terms land in the
        // vector body versus the tail, hence the reduction order.  The criterion
        // is the reduction-order budget, as in tests/spa_cgf_test.cpp.
        const double sc0 = std::fabs(a.k0) + std::fabs(b.k0) + 1e-300;
        const double sc1 = std::fabs(a.k1) + std::fabs(b.k1) + 1e-300;
        const double sc2 = std::fabs(a.k2) + std::fabs(b.k2) + 1e-300;
        worst = std::fmax(worst, std::fabs(a.k0 - b.k0) / sc0);
        worst = std::fmax(worst, std::fabs(a.k1 - b.k1) / sc1);
        worst = std::fmax(worst, std::fabs(a.k2 - b.k2) / sc2);
    }
    CHECK(worst < 1e-13);
    std::printf("    h == 0 subjects perturbed the cumulants by %.3e (relative)\n",
                worst);
}

// ══════════════════════════════════════════════════════════════════════
// 3  A NaN residual propagates instead of becoming -infinity (D3)
// ══════════════════════════════════════════════════════════════════════

TEST(a_nan_residual_propagates_rather_than_becoming_minus_infinity) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> r{0.7, nan, -1.1}, h{2.0, 1.0, 2.0};
    Context c;
    c.resid = r.data();
    c.hap = h.data();
    c.nOutlier = 3;
    c.q = 0.3;

    const spa::K012 k = spamixlocalp_cgf::kFull(0.2, c, 0.0);
    // K' and K'' carry the residual as an outright factor, so a NaN residual
    // propagates through them in every tier.
    CHECK(std::isnan(k.k1));
    CHECK(std::isnan(k.k2));
    // K does NOT carry r as a factor — r enters only inside log(alpha) — and
    // it is exactly there that the predecessor substituted -infinity, because
    // its guard `if (base > 1e-15)` is FALSE for a NaN base.  Whatever K comes
    // out as, it must not be that.
    //
    // K is deliberately NOT asserted to be NaN.  util/spa_cgf.hpp states that a
    // non-finite t*r is outside the contract of these kernels, and it is: with
    // AVX-512 selected, `avx512_log_pd` decomposes its argument by masking the
    // exponent field and returns a finite value (~708 + log(q)) for a NaN
    // input, so binomHapcountK0_avx512 launders the NaN while the scalar and
    // AVX2 tiers return it.  That divergence is in util/simd_math.hpp, which is
    // shared validated infrastructure this stage does not modify, and it is
    // harmless here precisely because K' and K'' above are NaN in every tier:
    // the solver's residual is non-finite, the retreat toward the origin finds
    // nothing usable, and the result is a failure status rather than a number.
    CHECK(k.k0 != -std::numeric_limits<double>::infinity());

    // Which is the property that actually matters: the tail must report the
    // condition rather than turning it into a probability.  With no z to fall
    // back to (decision D5 does not manufacture one) the row is NA; with a z
    // it would be the named normal substitution, never the saddlepoint value.
    const spa::TwoSided ts = spamixlocalp_cgf::twoSidedSpa(
        c, 0.0, 1.0, 1.0, std::numeric_limits<double>::quiet_NaN());
    CHECK(std::isnan(ts.p));
    CHECK(std::isnan(ts.negLog10p));
    CHECK(ts.status == spa::Status::NaNoTest);
    std::printf("    NaN residual: K'' = %s, K = %.6g, status = %s\n",
                std::isnan(k.k2) ? "nan" : "finite", k.k0,
                spa::statusName(ts.status));
}

// ══════════════════════════════════════════════════════════════════════
// 4  Derivative bookkeeping, Gaussian block included
// ══════════════════════════════════════════════════════════════════════

TEST(k1_and_k2_are_the_derivatives_of_k0) {
    Fixture f = makeFixture(41, 907u, 0.31);
    for (double &r : f.resid) r *= 0.05;   // keep the central differences tame
    const Context c = ctxOf(f);
    const double s = 2.75;

    double worst1 = 0.0, worst2 = 0.0;
    for (double t : {-0.6, -0.1, 0.1, 0.6}) {
        const double hh = 1e-5;
        const double k0p = spamixlocalp_cgf::kFull(t + hh, c, s).k0;
        const double k0m = spamixlocalp_cgf::kFull(t - hh, c, s).k0;
        const double k1p = spamixlocalp_cgf::kFull(t + hh, c, s).k1;
        const double k1m = spamixlocalp_cgf::kFull(t - hh, c, s).k1;
        // kFull.k1 is K' - s, so the central difference of K0 equals k1 + s.
        const double fd1 = (k0p - k0m) / (2.0 * hh);
        const double fd2 = (k1p - k1m) / (2.0 * hh);
        const spa::K012 at = spamixlocalp_cgf::kFull(t, c, s);
        worst1 = std::fmax(worst1,
                           std::fabs(fd1 - (at.k1 + s)) / std::fabs(at.k1 + s));
        worst2 = std::fmax(worst2, std::fabs(fd2 - at.k2) / std::fabs(at.k2));
    }
    CHECK(worst1 < 1e-7);
    CHECK(worst2 < 1e-7);
    std::printf("    central-difference agreement  K' %.3e  K'' %.3e\n",
                worst1, worst2);
}

// ══════════════════════════════════════════════════════════════════════
// 5  kFull and k12 agree bit-for-bit on the derivatives
// ══════════════════════════════════════════════════════════════════════

TEST(kfull_and_k12_agree_bitwise_on_the_derivatives) {
    const Fixture f = makeFixture(53, 313u, 0.17);
    const Context c = ctxOf(f);
    for (double t : {-0.02, -1e-12, 0.0, 1e-12, 0.02, 0.007}) {
        const spa::K12 a = spamixlocalp_cgf::k12(t, c, 0.9);
        const spa::K012 b = spamixlocalp_cgf::kFull(t, c, 0.9);
        CHECK(a.k1 == b.k1);
        CHECK(a.k2 == b.k2);
    }
}

// ══════════════════════════════════════════════════════════════════════
// 6  The score enters K1 and nothing else
// ══════════════════════════════════════════════════════════════════════

TEST(the_score_is_subtracted_in_k1_only) {
    const Fixture f = makeFixture(33, 77u, 0.44);
    const Context c = ctxOf(f);
    const double t = 0.013;
    const spa::K012 a = spamixlocalp_cgf::kFull(t, c, 0.0);
    const spa::K012 b = spamixlocalp_cgf::kFull(t, c, 5.5);
    CHECK(a.k0 == b.k0);
    CHECK(a.k2 == b.k2);
    CHECK(std::fabs((a.k1 - b.k1) - 5.5) < 1e-12 * 5.5);
}

// ══════════════════════════════════════════════════════════════════════
// 7  The solver is exact on a purely Gaussian CGF
// ══════════════════════════════════════════════════════════════════════

TEST(twosided_reproduces_the_normal_tail_when_there_are_no_outliers) {
    Context c;
    c.resid = nullptr;
    c.hap = nullptr;
    c.nOutlier = 0;
    c.q = 0.0;
    c.mean = -4.25;          // K'(0)
    c.var = 6.5 * 6.5;       // K''(0) = sigma^2

    const double sigma = 6.5;
    double worstP = 0.0, worstL = 0.0;
    for (double z : {2.5, 3.0, 5.0, 8.0, 12.0, 20.0}) {
        const double absDev = z * sigma;
        const spa::TwoSided ts =
            spamixlocalp_cgf::twoSidedSpa(c, c.mean, absDev, c.var, z);
        CHECK(ts.status == spa::Status::SpaOk);
        const double want = spa::normalTwoSided(z);
        worstP = std::fmax(worstP, std::fabs(ts.p - want) / want);
        const double wantL = -spa::normalTwoSidedLog(z) / std::log(10.0);
        worstL = std::fmax(worstL, std::fabs(ts.negLog10p - wantL) / wantL);
    }
    // For a Gaussian CGF v == w, so log(v/w) == 0 and r* == w exactly; the only
    // error is the solver's residual tolerance and the two square roots.
    CHECK(worstP < 1e-9);
    CHECK(worstL < 1e-11);
    std::printf("    Gaussian CGF: worst rel err  P %.3e   LOG10P %.3e\n",
                worstP, worstL);
}

// ══════════════════════════════════════════════════════════════════════
// 8  D3(b) — one failed tail is NA, not half a p-value
// ══════════════════════════════════════════════════════════════════════

TEST(a_one_sided_failure_is_na_and_not_half_a_pvalue) {
    // A CGF with no Gaussian block and a single outlier carrying h = 2:
    // K'(t) = h*r*q*e^{tr}/alpha is confined to (0, h*r) for r > 0, so a score
    // outside that open interval has no saddlepoint at all.
    //
    // The reflection point is deliberately NOT K'(0): SPAmixLocalPlus reflects
    // about the VARIANCE-RESCALED mean sMean*sqrt(varDiag/varS), which is not
    // the CGF's own mean whenever the phi block is non-empty, so an asymmetric
    // pair of targets is the production shape and not a contrivance.  With
    // sup K' = 3 and the reflection point at 2, the upper target 3.5 lies above
    // the supremum and is unreachable while the lower target 0.5 is interior.
    std::vector<double> r{1.5}, h{2.0};
    Context c;
    c.resid = r.data();
    c.hap = h.data();
    c.nOutlier = 1;
    c.q = 0.25;
    c.mean = 0.0;
    c.var = 0.0;

    const double sup = 2.0 * 1.5;                    // sup K' = h*r = 3
    const double centre = 2.0;                       // inside (0, 3)
    const double absDev = 1.5;                       // upper 3.5 > sup, lower 0.5
    CHECK(centre + absDev > sup);
    CHECK(centre - absDev > 0.0 && centre - absDev < sup);

    // One tail fails.  The whole two-sided p falls back — it is NOT the
    // surviving tail, and it is not a sum of one saddlepoint tail and one
    // normal tail either, which would be neither quantity.
    const double zBad = 2.5;
    const spa::TwoSided ts =
        spamixlocalp_cgf::twoSidedSpa(c, centre, absDev, 1.0, zBad);
    CHECK(spa::statusIsFallback(ts.status));
    CHECK(ts.p == spa::normalTwoSided(zBad));
    CHECK(ts.negLog10p == -spa::normalTwoSidedLog(zBad) / std::log(10.0));

    // The predecessor would have returned the lower tail alone, a number
    // between 0 and 1 with nothing to mark it as half an answer.  Confirm that
    // the surviving tail really is a plausible-looking probability, so that the
    // test is about the reporting and not about an unreachable branch.
    spa::Status stLower = spa::Status::SpaOk;
    const double sLower = centre - absDev;
    const spa::Saddle sd = spa::solveSaddlepoint(
        sLower,
        [&](double t) { return spamixlocalp_cgf::k12(t, c, sLower); },
        [&](double t) { return spamixlocalp_cgf::kFull(t, c, sLower); },
        spa::SolveOpts{});
    const double pLower = spa::bnTail(sd.zeta, sLower, sd.K0, sd.K2, true, stLower);
    CHECK(sd.status == spa::Status::SpaOk);
    CHECK(stLower == spa::Status::SpaOk);
    CHECK(pLower > 0.0 && pLower < 1.0);
    std::printf("    one-sided failure: combined = normal fallback (%s); "
                "the tail the predecessor would have reported alone = %.6g\n",
                spa::statusName(ts.status), pLower);
}

// ══════════════════════════════════════════════════════════════════════
// 9  The q endpoints are exact and finite
// ══════════════════════════════════════════════════════════════════════

TEST(saturated_and_vanishing_allele_frequency_are_exact) {
    std::vector<double> r{0.7, -1.3, 2.1, 0.4}, h{2.0, 1.0, 2.0, 0.0};
    Context c;
    c.resid = r.data();
    c.hap = h.data();
    c.nOutlier = 4;
    c.mean = 0.0;
    c.var = 0.0;

    double sumHR = 0.0;
    for (int i = 0; i < 4; ++i) sumHR += h[static_cast<size_t>(i)] * r[static_cast<size_t>(i)];

    // q == 1: G_i == h_i almost surely.  K = t*sum(h r), K' = sum(h r), K'' = 0.
    c.q = 1.0;
    for (double t : {-400.0, -1.0, 0.0, 1.0, 400.0}) {
        const spa::K012 k = spamixlocalp_cgf::kFull(t, c, 0.0);
        CHECK(std::isfinite(k.k0) && std::isfinite(k.k1) && std::isfinite(k.k2));
        CHECK_REL(k.k1, sumHR, 1e-15);
        CHECK(k.k2 == 0.0);
        CHECK_REL(k.k0, t * sumHR, 1e-15);
    }
    // Zero curvature everywhere: no root.  With no z supplied there is
    // nothing to substitute, so the row is NA with NA_NO_TEST.
    const spa::TwoSided tsOne = spamixlocalp_cgf::twoSidedSpa(
        c, sumHR, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN());
    CHECK(std::isnan(tsOne.p));
    CHECK(tsOne.status == spa::Status::NaNoTest);

    // q == 0: G_i == 0 almost surely, every cumulant exactly zero.
    c.q = 0.0;
    for (double t : {-400.0, -1.0, 0.0, 1.0, 400.0}) {
        const spa::K012 k = spamixlocalp_cgf::kFull(t, c, 0.0);
        CHECK(k.k0 == 0.0 && k.k1 == 0.0 && k.k2 == 0.0);
    }
    std::printf("    q endpoints: q=1 status %s, q=0 cumulants exactly zero\n",
                spa::statusName(tsOne.status));
}

TEST(status_code_matches_the_documented_integer_encoding) {
    CHECK(spamixlocalp_cgf::statusCode(spa::Status::SpaOk) == 0.0);
    CHECK(spamixlocalp_cgf::statusCode(spa::Status::Normal) == 1.0);
    CHECK(spamixlocalp_cgf::statusCode(spa::Status::SpaWSingular) == 2.0);
    CHECK(spamixlocalp_cgf::statusCode(spa::Status::FallbackMaxIter) == 3.0);
    CHECK(spamixlocalp_cgf::statusCode(spa::Status::FallbackGuardTemp) == 4.0);
    CHECK(spamixlocalp_cgf::statusCode(spa::Status::FallbackGuardCurv) == 5.0);
    CHECK(spamixlocalp_cgf::statusCode(spa::Status::FallbackNonFinite) == 6.0);
    CHECK(spamixlocalp_cgf::statusCode(spa::Status::NaPostFail) == 7.0);
    CHECK(spamixlocalp_cgf::statusCode(spa::Status::NaNoTest) == 8.0);
}

}  // namespace

TINYTEST_MAIN
