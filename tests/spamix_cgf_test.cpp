// spamix_cgf_test.cpp — the SPAmix / SPAGxE / SPAGxEmix CGF (spa_unify Stage 5).
//
// src/spamix/spamix_cgf.hpp is header-only; it needs util/spa_cgf.o for the
// dispatched binomial kernels and nothing else from src/.
//
// The migration this file guards is unusual in one specific way, and that is
// what most of the tests are about.  `math::kG012`, the kernel being replaced,
// differentiates with respect to x = t*r and leaves the chain-rule factors to
// its caller:
//
//     K'  += resid[i]         * k1
//     K'' += resid[i]*resid[i] * k2
//
// whereas `spa_cgf` differentiates with respect to t and takes the residual as
// a kernel argument, so its K' and K'' already carry r and r^2.  Porting
// 01_findings.md D1's cancellation-free rewrite into kG012 while leaving the
// old call sites alone would have multiplied K'' by a spurious r^2 — a factor
// that is invisible when every |r| is near 1 and catastrophic otherwise.  Test
// 1 is that factor, checked against a long-double reference written from the
// mathematics rather than from either kernel.
//
// In order, the tests establish:
//
//   1. THE RESIDUAL CONVENTION IS RIGHT.  Sum_i K'_i and Sum_i K''_i match a
//      long-double reference over residuals spanning four orders of magnitude,
//      so neither an extra nor a missing power of r survives.  Run with
//      deliberately non-unit residuals: at |r| = 1 a spurious r^2 is invisible.
//
//   2. THE UNIFORM AND PER-INDIVIDUAL PATHS AGREE where they must.  SPAGxE base
//      is switched from the per-individual kernel (fed a constant vector) to
//      the uniform one; with af[i] identically q the two must produce the same
//      cumulants.
//
//   3. THE DERIVATIVE BOOKKEEPING IS RIGHT, including the Gaussian block.
//      K'(t) is the central difference of K(t) and K''(t) that of K'(t), with
//      the non-outlier mean and variance populated, so a mis-scaled
//      mean*t + var*t^2/2 term shows up here.
//
//   4. kFull AND k12 AGREE BIT-FOR-BIT on K' and K''.  The tail is evaluated
//      with kFull's cumulants and the root was found with k12's; if they
//      drifted, the two would describe different problems.
//
//   5. THE SCORE IS SUBTRACTED IN K1 ONLY.  s must not leak into K0 or K2.
//
//   6. THE SOLVER IS EXACT ON A GAUSSIAN CGF.  With no outliers the CGF is
//      exactly mu*t + sigma^2*t^2/2, for which the Barndorff-Nielsen tail is
//      exact: zeta = (s-mu)/sigma^2 gives v == w and r* == w.  twoSidedSpa
//      must therefore reproduce 2*Phi(-|dev|/sigma) to near machine precision.
//      This exercises the bracketing, the Newton loop, both tails and the
//      two-sided assembly end to end against a closed form.
//
//   7. A SATURATED ALLELE FREQUENCY IS REPORTED, NOT SILENTLY NUMBERED.  When
//      every af[i] is exactly 1 the genotype is deterministic, K'' vanishes and
//      no test exists.  This is the corner three markers of examples/1kg reach
//      through a completely separated logistic AF model, and the result must be
//      NaN with a naming status rather than a finite-looking p-value.

#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <vector>

#include "spamix/spamix_cgf.hpp"
#include "tinytest.hpp"

namespace {

using spamix_cgf::Context;
using ld = long double;

// ──────────────────────────────────────────────────────────────────────
// Long-double reference, written from the mathematics
// ──────────────────────────────────────────────────────────────────────
//
// One subject contributes the CGF of r*G with G ~ Binomial(2, p):
//     K   = 2*log((1-p) + p*e^{t r})
//     K'  = 2*r*p*e^{t r} / alpha
//     K'' = 2*r^2*p*(1-p)*e^{t r} / alpha^2
// The residual factors are written out explicitly here so that the test states
// the convention rather than inheriting it.
struct Ref {
    ld K0 = 0.0L, K1 = 0.0L, K2 = 0.0L;
};

Ref refSubject(ld t, ld r, ld p) {
    const ld u = 1.0L - p;
    const ld lam = std::exp(t * r);
    const ld e = p * lam;
    const ld a = u + e;
    Ref o;
    o.K0 = 2.0L * std::log(a);
    o.K1 = 2.0L * r * e / a;
    o.K2 = 2.0L * r * r * e * u / (a * a);
    return o;
}

Ref refSum(double t, const std::vector<double> &r, const std::vector<double> &p,
           double mean, double var) {
    Ref o;
    for (size_t i = 0; i < r.size(); ++i) {
        const Ref s = refSubject(static_cast<ld>(t), static_cast<ld>(r[i]),
                                 static_cast<ld>(p[i]));
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

// A subject set with residuals spanning four orders of magnitude, so that a
// spurious or missing power of r cannot hide.
struct Fixture {
    std::vector<double> resid, af;
    double mean = 0.0, var = 0.0;
};

Fixture makeFixture(int n, uint32_t seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<double> uExp(-2.0, 2.0);   // |r| in 1e-2..1e2
    std::uniform_real_distribution<double> uAf(0.001, 0.999);
    std::uniform_real_distribution<double> uSign(-1.0, 1.0);
    Fixture f;
    f.resid.resize(static_cast<size_t>(n));
    f.af.resize(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double mag = std::pow(10.0, uExp(rng));
        f.resid[static_cast<size_t>(i)] = (uSign(rng) < 0.0) ? -mag : mag;
        f.af[static_cast<size_t>(i)] = uAf(rng);
    }
    f.mean = 3.5;
    f.var = 12.25;
    return f;
}

Context ctxOf(const Fixture &f, bool perIndiv, double q = 0.0) {
    Context c;
    c.resid = f.resid.data();
    c.af = perIndiv ? f.af.data() : nullptr;
    c.nOutlier = static_cast<int>(f.resid.size());
    c.q = q;
    c.mean = f.mean;
    c.var = f.var;
    return c;
}

// ══════════════════════════════════════════════════════════════════════
// 1  The residual-factor convention
// ══════════════════════════════════════════════════════════════════════

TEST(residual_factor_is_inside_the_kernel_exactly_once) {
    const Fixture f = makeFixture(97, 20260730u);
    const Context c = ctxOf(f, /*perIndiv=*/true);

    // t chosen so that t*r stays inside the clamp for the largest |r| ~ 100.
    const double ts[] = {-3.0, -0.4, -1e-7, 0.0, 1e-7, 0.4, 3.0};
    double worst1 = 0.0, worst2 = 0.0, worst0 = 0.0;
    for (double t : ts) {
        const double tt = t / 100.0;
        const Ref want = refSum(tt, f.resid, f.af, f.mean, f.var);
        const spa::K12 got12 = spamix_cgf::k12(tt, c, 0.0);
        const spa::K012 got012 = spamix_cgf::kFull(tt, c, 0.0);
        worst1 = std::fmax(worst1, relErr(got12.k1, want.K1));
        worst2 = std::fmax(worst2, relErr(got12.k2, want.K2));
        // K is judged in the ABSOLUTE sense in which the tail consumes it.
        // util/spa_cgf.hpp's terminal-K section is explicit that the production
        // spelling has no relative accuracy near t = 0 — alpha = 1 + delta is
        // rounded before the logarithm, costing an absolute error of order eps
        // per subject whatever delta is — and the tail reads K only through
        // w = sgn(zeta)*sqrt(2*(zeta*s - K)), where dw = -dK/w with |w| of
        // order 2.  A relative criterion here would test the fixture's
        // cancellation, not the kernel: at t*r ~ 1e-9 the per-subject terms are
        // 2*p*r*t and their mixed-sign sum is six orders below the terms that
        // formed it.
        worst0 = std::fmax(
            worst0, std::fabs(got012.k0 - static_cast<double>(want.K0)));
    }
    // A spurious r^2 would show up as an error of order 1e4 in K'', not 1e-15.
    CHECK(worst1 < 1e-12);
    CHECK(worst2 < 1e-12);
    CHECK(worst0 < 1e-11);
    std::printf("    worst err vs long double: K %.3e (abs)  "
                "K' %.3e (rel)  K'' %.3e (rel)\n",
                worst0, worst1, worst2);
}

TEST(a_spurious_residual_squared_would_be_detected) {
    // Guards the guard: confirm the fixture's residuals really are far enough
    // from unity that the test above could see an extra r^2 factor.
    const Fixture f = makeFixture(97, 20260730u);
    double lo = 1e300, hi = 0.0;
    for (double r : f.resid) {
        lo = std::fmin(lo, std::fabs(r));
        hi = std::fmax(hi, std::fabs(r));
    }
    CHECK(lo < 0.05);
    CHECK(hi > 20.0);
}

// ══════════════════════════════════════════════════════════════════════
// 2  Uniform q reproduces the constant-vector per-individual call
// ══════════════════════════════════════════════════════════════════════

TEST(uniform_and_per_individual_agree_when_af_is_constant) {
    Fixture f = makeFixture(64, 7u);
    const double q = 0.37;
    for (double &p : f.af) p = q;

    const Context cIndiv = ctxOf(f, /*perIndiv=*/true);
    const Context cUnif = ctxOf(f, /*perIndiv=*/false, q);

    double worst = 0.0;
    for (double t : {-0.05, -1e-9, 0.0, 1e-9, 0.05, 0.31}) {
        const spa::K012 a = spamix_cgf::kFull(t, cIndiv, 1.25);
        const spa::K012 b = spamix_cgf::kFull(t, cUnif, 1.25);
        worst = std::fmax(worst, std::fabs(a.k0 - b.k0) / (std::fabs(a.k0) + 1e-300));
        worst = std::fmax(worst, std::fabs(a.k1 - b.k1) / (std::fabs(a.k1) + 1e-300));
        worst = std::fmax(worst, std::fabs(a.k2 - b.k2) / (std::fabs(a.k2) + 1e-300));
    }
    // Same algebra, same lane order; only the hoisting of `1 - p` differs.
    CHECK(worst < 1e-14);
}

// ══════════════════════════════════════════════════════════════════════
// 3  Derivative bookkeeping, Gaussian block included
// ══════════════════════════════════════════════════════════════════════

TEST(k1_and_k2_are_the_derivatives_of_k0) {
    Fixture f = makeFixture(41, 99u);
    for (double &r : f.resid) r *= 0.05;   // keep the central differences tame
    const Context c = ctxOf(f, /*perIndiv=*/true);
    const double s = 2.75;

    double worst1 = 0.0, worst2 = 0.0;
    for (double t : {-0.6, -0.1, 0.1, 0.6}) {
        const double h = 1e-5;
        const double k0p = spamix_cgf::kFull(t + h, c, s).k0;
        const double k0m = spamix_cgf::kFull(t - h, c, s).k0;
        const double k1p = spamix_cgf::kFull(t + h, c, s).k1;
        const double k1m = spamix_cgf::kFull(t - h, c, s).k1;
        // kFull.k1 is K' - s, so the central difference of K0 must equal k1 + s.
        const double fd1 = (k0p - k0m) / (2.0 * h);
        const double fd2 = (k1p - k1m) / (2.0 * h);
        const spa::K012 at = spamix_cgf::kFull(t, c, s);
        worst1 = std::fmax(worst1, std::fabs(fd1 - (at.k1 + s)) / std::fabs(at.k1 + s));
        worst2 = std::fmax(worst2, std::fabs(fd2 - at.k2) / std::fabs(at.k2));
    }
    CHECK(worst1 < 1e-7);
    CHECK(worst2 < 1e-7);
}

// ══════════════════════════════════════════════════════════════════════
// 4  kFull and k12 agree bit-for-bit on the derivatives
// ══════════════════════════════════════════════════════════════════════

TEST(kfull_and_k12_agree_bitwise_on_the_derivatives) {
    const Fixture f = makeFixture(53, 4242u);
    const Context cI = ctxOf(f, true);
    const Context cU = ctxOf(f, false, 0.21);
    for (double t : {-0.02, -1e-12, 0.0, 1e-12, 0.02, 0.007}) {
        for (const Context *c : {&cI, &cU}) {
            const spa::K12 a = spamix_cgf::k12(t, *c, 0.9);
            const spa::K012 b = spamix_cgf::kFull(t, *c, 0.9);
            CHECK(a.k1 == b.k1);
            CHECK(a.k2 == b.k2);
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
// 5  The score enters K1 and nothing else
// ══════════════════════════════════════════════════════════════════════

TEST(the_score_is_subtracted_in_k1_only) {
    const Fixture f = makeFixture(33, 11u);
    const Context c = ctxOf(f, true);
    const double t = 0.013;
    const spa::K012 a = spamix_cgf::kFull(t, c, 0.0);
    const spa::K012 b = spamix_cgf::kFull(t, c, 5.5);
    CHECK(a.k0 == b.k0);
    CHECK(a.k2 == b.k2);
    CHECK(std::fabs((a.k1 - b.k1) - 5.5) < 1e-12 * 5.5);
}

// ══════════════════════════════════════════════════════════════════════
// 6  The solver is exact on a purely Gaussian CGF
// ══════════════════════════════════════════════════════════════════════

TEST(twosided_reproduces_the_normal_tail_when_there_are_no_outliers) {
    Context c;
    c.resid = nullptr;
    c.af = nullptr;
    c.nOutlier = 0;
    c.q = 0.0;
    c.mean = -4.25;          // K'(0)
    c.var = 6.5 * 6.5;       // K''(0) = sigma^2

    const double sigma = 6.5;
    double worstP = 0.0, worstL = 0.0;
    for (double z : {2.5, 3.0, 5.0, 8.0, 12.0, 20.0}) {
        const double absDev = z * sigma;
        const spa::Result ts = spamix_cgf::twoSidedSpa(c, c.mean, absDev, c.var, z);
        CHECK(ts.status == spa::Status::SpaOk);
        // `spa::normalTwoSided` is gone with the linear tail path
        // (log10p_unify Stage 3); this is the expression it evaluated.
        const double want =
            2.0 * math::pnorm(z, 0.0, 1.0, /*lower_tail=*/false);
        worstP = std::fmax(
            worstP, std::fabs(spa::pFromNegLog10P(ts.negLog10p) - want) / want);
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

TEST(twosided_agrees_with_the_normal_tail_when_the_outlier_block_is_negligible) {
    // One outlier carrying a vanishing residual: the CGF is Gaussian to within
    // that subject's contribution, so the tail must still be close to normal.
    std::vector<double> r{1e-6}, p{0.4};
    Context c;
    c.resid = r.data();
    c.af = p.data();
    c.nOutlier = 1;
    c.mean = 0.0;
    c.var = 4.0;
    const spa::Result ts = spamix_cgf::twoSidedSpa(c, 0.0, 3.0 * 2.0, c.var, 3.0);
    CHECK(ts.status == spa::Status::SpaOk);
    CHECK_REL(spa::pFromNegLog10P(ts.negLog10p),
              2.0 * math::pnorm(3.0, 0.0, 1.0, /*lower_tail=*/false), 1e-5);
}

// ══════════════════════════════════════════════════════════════════════
// 7  A saturated allele frequency is reported, not silently numbered
// ══════════════════════════════════════════════════════════════════════

TEST(saturated_allele_frequency_yields_zero_curvature_not_a_wrong_number) {
    // af == 1 exactly: G == 2 almost surely, so K' = 2*r and K'' = 0.
    std::vector<double> r{0.7, -1.3, 2.1}, p{1.0, 1.0, 1.0};
    Context c;
    c.resid = r.data();
    c.af = p.data();
    c.nOutlier = 3;
    c.mean = 0.0;
    c.var = 0.0;

    const spa::K012 k = spamix_cgf::kFull(0.05, c, 0.0);
    const double sumR = r[0] + r[1] + r[2];
    CHECK(std::fabs(k.k1 - 2.0 * sumR) < 1e-12 * std::fabs(2.0 * sumR));
    CHECK(k.k2 == 0.0);
    CHECK(std::fabs(k.k0 - 2.0 * 0.05 * sumR) < 1e-12);

    // With zero curvature everywhere there is no root, so the saddlepoint
    // fails.  Under decision D5 the row is not silently numbered and not
    // discarded either: it reports the two-sided normal tail at the caller's
    // own z, under a FALLBACK_* status naming which guard fired.  The
    // saddlepoint value itself never reaches the output.
    const double zBad = 4.0;
    const spa::Result ts = spamix_cgf::twoSidedSpa(c, 2.0 * sumR, 1.0, 0.0, zBad);
    CHECK(spa::statusIsFallback(ts.status));
    CHECK(ts.negLog10p == -spa::normalTwoSidedLog(zBad) / std::log(10.0));
    std::printf("    saturated af: status = %s\n", spa::statusName(ts.status));

    // ...and when the caller has no z either, there is nothing to fall back
    // to and the row stays NA.
    const spa::Result tsNa = spamix_cgf::twoSidedSpa(
        c, 2.0 * sumR, 1.0, 0.0, std::numeric_limits<double>::quiet_NaN());
    CHECK(std::isnan(tsNa.negLog10p));
    CHECK(tsNa.status == spa::Status::NaNoTest);
}

TEST(status_code_matches_the_documented_integer_encoding) {
    CHECK(spamix_cgf::statusCode(spa::Status::SpaOk) == 0.0);
    CHECK(spamix_cgf::statusCode(spa::Status::Normal) == 1.0);
    CHECK(spamix_cgf::statusCode(spa::Status::SpaWSingular) == 2.0);
    CHECK(spamix_cgf::statusCode(spa::Status::FallbackMaxIter) == 3.0);
    CHECK(spamix_cgf::statusCode(spa::Status::FallbackGuardTemp) == 4.0);
    CHECK(spamix_cgf::statusCode(spa::Status::FallbackGuardCurv) == 5.0);
    CHECK(spamix_cgf::statusCode(spa::Status::FallbackNonFinite) == 6.0);
    CHECK(spamix_cgf::statusCode(spa::Status::NaPostFail) == 7.0);
    CHECK(spamix_cgf::statusCode(spa::Status::NaNoTest) == 8.0);
}

}  // namespace

TINYTEST_MAIN
