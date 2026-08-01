// wtcoxg_cgf_test.cpp — the WtCoxG / LEAF saddlepoint and the bivariate
// normal integrator it feeds (spa_unify Stage 6).
//
// src/wtcoxg/wtcoxg_cgf.hpp is header-only over util/spa_cgf.o; the bivariate
// tests additionally need util/math_helper.o.  Nothing else from src/ is
// linked, so this file exercises the migrated kernel rather than the method.
//
// In order, the tests establish:
//
//   1. THE CUMULANTS ARE RIGHT.  K, K' and K'' match a long-double reference
//      written from the mathematics, over residuals spanning four orders of
//      magnitude so that a missing or spurious power of r cannot hide.
//
//   2. D4 IS PRESERVED.  WtCoxG hands the kernel a K'' that is NOT the second
//      derivative of its K: K carries var_n_K01 (the batch-effect Gaussian
//      variance) and K'' carries var_n_K2 (the finite-reference-panel
//      correction), and the two differ whenever obs_ct > 0.  This is
//      documented as intentional, and 02_design.md rules it a modelling
//      question rather than a defect.  The test pins the asymmetry in both
//      directions: kFull must report the mismatched pair unchanged, and
//      bnTail must consume it without "correcting" either member.  A future
//      change that quietly makes the pair consistent will fail here, which is
//      the point.
//
//   3. THE DERIVATIVE BOOKKEEPING IS RIGHT.  K'(t) is the central difference
//      of K(t) and K''(t) that of K'(t) — with varK01 == varK2 so that the
//      relation actually holds — so a mis-scaled Gaussian block shows up.
//
//   4. kFull AND k12 AGREE BIT-FOR-BIT on K' and K''.
//
//   5. THE SOLVER IS EXACT ON A GAUSSIAN CGF.  With no outliers the CGF is
//      exactly mean*t + var*t^2/2, for which the Barndorff-Nielsen tail is
//      exact, so twoSidedSpa must reproduce 2*Phi(-|dev|/sigma).
//
//   6. D2 IS REPAIRED.  A degenerate input must produce NaN with a naming
//      status, never the 1.0 that `std::min(1.0, pval1 + pval2)` returned for
//      a NaN tail.  This is tested by construction because the bundled
//      fixture never reaches a failing saddlepoint: over all 2318 SPA-branch
//      entries of examples/baseline.sh, zero tails returned NaN, so the
//      repair cannot be demonstrated from the pinned output alone.
//
//   7. math::bvnCdf IS CORRECT AT |r| >= 0.925.  That branch was
//      mis-transcribed from Genz (2004) and returned joint probabilities
//      larger than the marginals bounding them, and values outside [0, 1];
//      it is the origin of the 166 markers with P_EXT > 1 in the pinned
//      baseline.  Pinned here against an independent Simpson reduction of
//      Phi2(h,k;r) = int_{-inf}^{h} phi(x) Phi((k-rx)/sqrt(1-r^2)) dx.
//
//   8. pmvnorm2dHalfRect NEVER RETURNS A JOINT ABOVE ITS OWN MARGINALS, and
//      reports rather than clamps a covariance matrix that is not positive
//      semi-definite.  The second half is what stops WtCoxG's conditional
//      p-value from exceeding 1 or collapsing to a fabricated 1e-13.
//
//   9. THE IMMATERIAL-LEG RULE IS EXACTLY WHAT IT CLAIMS TO BE.  Item 8 makes
//      a mixture leg unavailable whenever sigma2 drives |rho| past 1, which on
//      a degenerate batch fit is most markers.  wtcoxg_cond::conditionalP
//      keeps such a marker only when the interval that leg's unknown
//      contribution spans is narrower than kImmaterialLeg times the reported
//      value; the tests pin the arithmetic of that interval, both sides of the
//      threshold, that the surviving leg's status is what gets reported, and
//      that a leg which does matter still produces NaN rather than the 1.0
//      that `std::min(1.0, ...)` used to manufacture.

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "tinytest.hpp"
#include "util/math_helper.hpp"
#include "wtcoxg/conditional_p.hpp"
#include "wtcoxg/wtcoxg_cgf.hpp"

namespace {

using wtcoxg_cgf::Context;
using ld = long double;

// ──────────────────────────────────────────────────────────────────────
// Long-double reference, written from the mathematics
// ──────────────────────────────────────────────────────────────────────
//
// One subject contributes the CGF of r*G with G ~ Binomial(2, p):
//     K   = 2*log((1-p) + p*e^{t r})
//     K'  = 2*r*p*e^{t r} / alpha
//     K'' = 2*r^2*p*(1-p)*e^{t r} / alpha^2
struct Ref {
    ld K0 = 0.0L, K1 = 0.0L, K2 = 0.0L;
};

Ref refSum(double t, const std::vector<double> &r, double p,
           double mean, double varK01, double varK2) {
    Ref o;
    for (double ri : r) {
        const ld tl = static_cast<ld>(t), rl = static_cast<ld>(ri),
                 pl = static_cast<ld>(p);
        const ld u = 1.0L - pl;
        const ld e = pl * std::exp(tl * rl);
        const ld a = u + e;
        o.K0 += 2.0L * std::log(a);
        o.K1 += 2.0L * rl * e / a;
        o.K2 += 2.0L * rl * rl * e * u / (a * a);
    }
    const ld tl = static_cast<ld>(t);
    o.K0 += static_cast<ld>(mean) * tl + 0.5L * static_cast<ld>(varK01) * tl * tl;
    o.K1 += static_cast<ld>(mean) + static_cast<ld>(varK01) * tl;
    o.K2 += static_cast<ld>(varK2);
    return o;
}

double relErr(double got, ld want) {
    const ld w = std::fabs(want);
    if (w == 0.0L) return std::fabs(got);
    return static_cast<double>(std::fabs(static_cast<ld>(got) - want) / w);
}

// Residuals spanning four orders of magnitude, so a spurious r or r^2 shows.
const std::vector<double> kResid = {
    -1.7e-2, 3.1e-1, -8.4e-1, 2.05, -3.6, 1.1e-3, 7.7e-1, -2.2,
    5.9e-2, -4.4e-1, 1.35, -6.6e-3, 9.1e-1, -1.9, 2.7e-1, 4.8e-2,
    -3.3e-1, 1.02, -7.5e-2, 6.3e-1
};

Context makeCtx(double p, double mean, double varK01, double varK2) {
    Context c;
    c.resid = kResid.data();
    c.nOutlier = static_cast<int>(kResid.size());
    c.maf = p;
    c.mean = mean;
    c.varK01 = varK01;
    c.varK2 = varK2;
    return c;
}

double Phi(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }

// Independent reference for Phi2(h, k; r), by composite Simpson on
//   int_{-inf}^{h} phi(x) * Phi((k - r x) / sqrt(1 - r^2)) dx.
// Deliberately unrelated to Genz's expansion.
double bvnRef(double h, double k, double r) {
    if (std::fabs(r) >= 1.0)
        return (r > 0.0) ? std::fmin(Phi(h), Phi(k))
                         : std::fmax(0.0, Phi(h) + Phi(k) - 1.0);
    const double s = std::sqrt((1.0 - r) * (1.0 + r));
    const double lo = -14.0;
    const double hi = std::fmin(h, 14.0);
    if (hi <= lo) return 0.0;
    const long N = 400000;
    const double dx = (hi - lo) / N;
    double sum = 0.0;
    for (long i = 0; i <= N; ++i) {
        const double x = lo + i * dx;
        const double f = std::exp(-0.5 * x * x) / std::sqrt(2.0 * M_PI) *
                         Phi((k - r * x) / s);
        const double w = (i == 0 || i == N) ? 1.0 : ((i % 2) ? 4.0 : 2.0);
        sum += w * f;
    }
    return sum * dx / 3.0;
}

// ──────────────────────────────────────────────────────────────────────
// 1  Cumulants against the long-double reference
// ──────────────────────────────────────────────────────────────────────

TEST(cumulants_match_the_long_double_reference) {
    const Context c = makeCtx(0.23, 1.7, 0.9, 0.9);
    double worst0 = 0.0, worst1 = 0.0, worst2 = 0.0;
    for (double t : {-2.0, -0.7, -0.05, 0.0, 0.05, 0.7, 2.0, 6.0}) {
        const Ref want = refSum(t, kResid, c.maf, c.mean, c.varK01, c.varK2);
        const spa::K012 got = wtcoxg_cgf::kFull(t, c, 0.0);
        worst0 = std::fmax(worst0, relErr(got.k0, want.K0));
        worst1 = std::fmax(worst1, relErr(got.k1, want.K1));
        worst2 = std::fmax(worst2, relErr(got.k2, want.K2));
    }
    std::printf("    worst rel err  K %.3e  K' %.3e  K'' %.3e\n",
                worst0, worst1, worst2);
    CHECK(worst0 < 1e-13);
    CHECK(worst1 < 1e-13);
    CHECK(worst2 < 1e-13);
}

// ──────────────────────────────────────────────────────────────────────
// 2  D4 — the (K, K'') pair stays inconsistent, on purpose
// ──────────────────────────────────────────────────────────────────────

TEST(d4_the_mismatched_variance_pair_is_carried_through_unchanged) {
    // var_n_K01 != var_n_K2 is what WtCoxG's Branch B supplies whenever the
    // external reference panel is retained (obs_ct > 0).
    const double v01 = 0.90, v2 = 1.45;
    const Context c = makeCtx(0.23, 1.7, v01, v2);

    // kFull must report K built from v01 and K'' built from v2.
    const Ref withV01 = refSum(0.4, kResid, c.maf, c.mean, v01, v01);
    const Ref withV2  = refSum(0.4, kResid, c.maf, c.mean, v2,  v2);
    const spa::K012 got = wtcoxg_cgf::kFull(0.4, c, 0.0);
    CHECK(relErr(got.k0, withV01.K0) < 1e-13);   // K   uses varK01
    CHECK(relErr(got.k2, withV2.K2)  < 1e-13);   // K'' uses varK2
    CHECK(std::fabs(got.k2 - static_cast<double>(withV01.K2)) > 0.5 * (v2 - v01));

    // k12's residual uses varK01 while its Jacobian uses varK2, exactly as the
    // predecessor's Newton step did.
    const spa::K12 d = wtcoxg_cgf::k12(0.4, c, 0.0);
    CHECK(d.k1 == got.k1);
    CHECK(d.k2 == got.k2);

    // And spa::bnTail must consume the pair as two opaque scalars: feeding it
    // the same zeta and s with K2 replaced by the "consistent" value must move
    // the answer, i.e. the kernel is not silently recomputing one from the
    // other.  s is taken at the saddlepoint, s = K'(zeta), so that
    // zeta*s - K(zeta) >= 0 by convexity and no guard fires.
    const double sAtRoot = got.k1;   // kFull was called with s = 0
    spa::Status stA = spa::Status::SpaOk, stB = spa::Status::SpaOk;
    const double pA = spa::bnTail(0.4, sAtRoot, got.k0, got.k2, false, stA);
    const double pB = spa::bnTail(0.4, sAtRoot, got.k0,
                                  static_cast<double>(withV01.K2), false, stB);
    CHECK(stA == spa::Status::SpaOk);
    CHECK(stB == spa::Status::SpaOk);
    CHECK(pA != pB);
}

// ──────────────────────────────────────────────────────────────────────
// 3  Derivative bookkeeping
// ──────────────────────────────────────────────────────────────────────

TEST(k1_and_k2_are_the_derivatives_of_k0) {
    // varK01 == varK2 so that K'' really is d2K/dt2 (see D4 above).
    const Context c = makeCtx(0.31, -0.8, 1.25, 1.25);
    const double h = 1e-4;
    double worst1 = 0.0, worst2 = 0.0;
    for (double t : {-1.5, -0.3, 0.0, 0.3, 1.5}) {
        const double kp = wtcoxg_cgf::kFull(t + h, c, 0.0).k0;
        const double km = wtcoxg_cgf::kFull(t - h, c, 0.0).k0;
        const double k00 = wtcoxg_cgf::kFull(t, c, 0.0).k0;
        const double d1 = (kp - km) / (2.0 * h);
        const double d2 = (kp - 2.0 * k00 + km) / (h * h);
        const spa::K012 got = wtcoxg_cgf::kFull(t, c, 0.0);
        worst1 = std::fmax(worst1, std::fabs(d1 - got.k1) / std::fabs(got.k1));
        worst2 = std::fmax(worst2, std::fabs(d2 - got.k2) / std::fabs(got.k2));
    }
    std::printf("    central-difference agreement  K' %.3e  K'' %.3e\n", worst1, worst2);
    CHECK(worst1 < 1e-8);
    CHECK(worst2 < 1e-6);
}

// ──────────────────────────────────────────────────────────────────────
// 4  kFull and k12 agree
// ──────────────────────────────────────────────────────────────────────

TEST(kfull_and_k12_agree_bit_for_bit) {
    const Context c = makeCtx(0.17, 2.4, 0.7, 1.1);
    for (double t : {-3.0, -0.5, 0.0, 0.5, 3.0, 9.0}) {
        const spa::K012 f = wtcoxg_cgf::kFull(t, c, 1.75);
        const spa::K12  g = wtcoxg_cgf::k12(t, c, 1.75);
        CHECK(f.k1 == g.k1);
        CHECK(f.k2 == g.k2);
    }
}

// ──────────────────────────────────────────────────────────────────────
// 5  The solver is exact on a Gaussian CGF
// ──────────────────────────────────────────────────────────────────────

TEST(gaussian_cgf_reproduces_the_closed_form_two_sided_p) {
    // No outliers: K(t) = mean*t + var*t^2/2 exactly, and K'(0) = mean.  The
    // WtCoxG convention places the tails at +/-|dev| about ZERO, so the closed
    // form is Phi((-dev - mean)/sd) + 1 - Phi((dev - mean)/sd); with mean = 0
    // that is 2*Phi(-|dev|/sd), which is the case the method actually reaches
    // (K'(0) is identically zero for its CGF).
    Context c;
    c.resid = nullptr;
    c.nOutlier = 0;
    c.maf = 0.25;
    c.mean = 0.0;
    c.varK01 = 4.0;
    c.varK2 = 4.0;
    const double sd = 2.0;
    for (double dev : {1.0, 3.0, 7.0, 15.0, 40.0}) {
        const spa::TwoSided ts = wtcoxg_cgf::twoSidedSpa(c, dev, c.varK01, dev / sd);
        const double want = 2.0 * Phi(-dev / sd);
        CHECK(ts.status == spa::Status::SpaOk);
        if (want > 0.0) {
            const double rel = std::fabs(ts.p - want) / want;
            std::printf("    dev=%5.1f  p=%.12e  closed form=%.12e  rel=%.2e\n",
                        dev, ts.p, want, rel);
            CHECK(rel < 1e-9);
        }
        // LOG10P must agree with -log10(p) wherever the linear scale survives.
        if (ts.p > 1e-300) {
            const double lg = -std::log10(ts.p);
            CHECK(std::fabs(ts.negLog10p - lg) < 1e-9 * std::fmax(1.0, std::fabs(lg)));
        }
    }
}

// ──────────────────────────────────────────────────────────────────────
// 6  D2 — a failure is NaN with a status, never a masked 1.0
// ──────────────────────────────────────────────────────────────────────

TEST(d2_a_degenerate_saddlepoint_is_named_never_a_masked_one) {
    // p == 1 makes the genotype deterministic: alpha == p*lambda, K'' vanishes
    // identically, and no saddlepoint exists.  The predecessor's
    // `std::min(1.0, pval1 + pval2)` returned 1.0 here, i.e. a perfectly null
    // marker, because std::min(1.0, NaN) is 1.0.
    //
    // Under decision D5 the row is neither 1.0 nor NA: it is the two-sided
    // normal tail at the caller's z, under a FALLBACK_* status that names the
    // guard.  What is forbidden is the SILENT substitution, and the status
    // column is what makes this one loud.
    Context c = makeCtx(1.0, 0.0, 0.0, 0.0);
    const double zBad = 3.0;
    const spa::TwoSided ts = wtcoxg_cgf::twoSidedSpa(c, 3.0, 0.0, zBad);
    CHECK(spa::statusIsFallback(ts.status));
    CHECK(ts.p == spa::normalTwoSided(zBad));
    CHECK(ts.negLog10p == -spa::normalTwoSidedLog(zBad) / std::log(10.0));
    std::printf("    deterministic genotype: p=%g status=%s\n",
                ts.p, spa::statusName(ts.status));

    // With no z either, there is nothing to fall back to and the row is NA.
    const spa::TwoSided tsNa = wtcoxg_cgf::twoSidedSpa(
        c, 3.0, 0.0, std::numeric_limits<double>::quiet_NaN());
    CHECK(std::isnan(tsNa.p));
    CHECK(std::isnan(tsNa.negLog10p));
    CHECK(tsNa.status == spa::Status::NaNoTest);

    // And the idiom itself: this is what the old code did.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    CHECK(std::min(1.0, nan) == 1.0);              // the defect, demonstrated
    CHECK(std::isnan(spa::combineTails(0.5, nan, spa::Status::SpaOk,
                                       spa::Status::SpaOk).p));  // the repair
}

TEST(status_code_matches_the_documented_integer_encoding) {
    CHECK(wtcoxg_cgf::statusCode(spa::Status::SpaOk) == 0.0);
    CHECK(wtcoxg_cgf::statusCode(spa::Status::Normal) == 1.0);
    CHECK(wtcoxg_cgf::statusCode(spa::Status::SpaWSingular) == 2.0);
    CHECK(wtcoxg_cgf::statusCode(spa::Status::FallbackMaxIter) == 3.0);
    CHECK(wtcoxg_cgf::statusCode(spa::Status::FallbackGuardTemp) == 4.0);
    CHECK(wtcoxg_cgf::statusCode(spa::Status::FallbackGuardCurv) == 5.0);
    CHECK(wtcoxg_cgf::statusCode(spa::Status::FallbackNonFinite) == 6.0);
    CHECK(wtcoxg_cgf::statusCode(spa::Status::NaPostFail) == 7.0);
    CHECK(wtcoxg_cgf::statusCode(spa::Status::NaNoTest) == 8.0);
}

// ──────────────────────────────────────────────────────────────────────
// 7  math::bvnCdf over the whole correlation range
// ──────────────────────────────────────────────────────────────────────

// The two arms are pinned separately because they are different algorithms
// with different accuracy.  |r| < 0.925 is Plackett's identity under Genz's
// 6-point rule, good to ~1e-7 — this is the arm every call in the bundled
// fixture that is not degenerate lands in, and it is left exactly as it was.
// |r| >= 0.925 is the asymptotic expansion under the 20-point rule, good to
// ~1e-13, and it is the arm that was broken.
TEST(bvncdf_matches_an_independent_quadrature) {
    const double hs[] = {-4.0, -2.0, -1.0, -0.005, 0.0, 0.4, 1.5, 3.0};
    const double ks[] = {-3.0, -2.0, -0.392, 0.0, 0.392, 1.0, 2.5};
    const double lowR[]  = {-0.9, -0.5, -0.1, 0.0, 0.1, 0.5, 0.9, 0.9249};
    const double highR[] = {-1.0, -0.99999, -0.999, -0.99, -0.96, -0.93, -0.925,
                            0.925, 0.93, 0.96, 0.99, 0.999, 0.99999, 1.0};
    int outOfRange = 0;
    double wLow = 0.0, wHigh = 0.0, wh = 0, wk = 0, wr = 0;
    for (double r : lowR) for (double h : hs) for (double k : ks) {
        const double got = math::bvnCdf(h, k, r);
        if (!(got >= 0.0 && got <= 1.0)) ++outOfRange;
        wLow = std::fmax(wLow, std::fabs(got - bvnRef(h, k, r)));
    }
    for (double r : highR) for (double h : hs) for (double k : ks) {
        const double got = math::bvnCdf(h, k, r);
        if (!(got >= 0.0 && got <= 1.0)) ++outOfRange;
        const double e = std::fabs(got - bvnRef(h, k, r));
        if (e > wHigh) { wHigh = e; wh = h; wk = k; wr = r; }
    }
    std::printf("    worst |bvnCdf - reference|:  |r| < 0.925: %.4e   "
                "|r| >= 0.925: %.4e at (h,k,r) = (%g, %g, %g)\n",
                wLow, wHigh, wh, wk, wr);
    std::printf("    returns outside [0,1]: %d\n", outOfRange);
    CHECK(outOfRange == 0);
    // FLAGGED FOR THE MAINTAINER, deliberately not changed in this stage.
    // The measured worst error of the |r| < 0.925 arm is 1.20e-06 absolute,
    // because it uses Genz's 3-node (6-point) rule at every |r| where Genz
    // himself selects 6 / 12 / 20 nodes for |r| < 0.3 / < 0.75 / else.  That
    // is roughly 1e-6 of absolute error in p0 and p1, hence ~1e-5 in the
    // conditional p-value once divided by a p_deno of order 0.3 — immaterial
    // for a p near 1, but a large RELATIVE error for a p near 1e-5.  Raising
    // the node count would move every Branch-A and Branch-B marker, which is a
    // re-baselining decision rather than a migration one, so this stage
    // reproduces the arm exactly and records the number here instead.
    CHECK(wLow < 2e-6);
    CHECK(wHigh < 1e-11);
}

TEST(bvncdf_is_exact_at_the_degenerate_correlations) {
    // At |r| = 1 the pair is a single random variable and the CDF is closed
    // form.  This is the case WtCoxG's Branch B reaches after the covariance
    // clamp, and the case the old transcription got wrong by adding a spurious
    // Phi(h) + Phi(k) - 1.
    for (double h : {-2.0, -0.005, 0.0, 1.3}) for (double k : {-1.0, -0.392, 0.392, 2.0}) {
        CHECK_NEAR(math::bvnCdf(h, k, 1.0), std::fmin(Phi(h), Phi(k)), 1e-14);
        CHECK_NEAR(math::bvnCdf(h, k, -1.0), std::fmax(0.0, Phi(h) + Phi(k) - 1.0), 1e-14);
    }
}

// ──────────────────────────────────────────────────────────────────────
// 8  pmvnorm2dHalfRect: bounded by its marginals, and honest about |rho| > 1
// ──────────────────────────────────────────────────────────────────────

TEST(pmvnorm_joint_never_exceeds_either_marginal) {
    // The joint P(X <= s, a <= Y <= b) cannot exceed P(X <= s) nor
    // P(a <= Y <= b).  Violating the second is exactly how the old bvnCdf let
    // WtCoxG's conditional p-value reach 2.99: the numerator came back larger
    // than the denominator that normalises it.
    const double var1 = 17.694, var2 = 3.5031e-5;
    const double sd1 = std::sqrt(var1), sd2 = std::sqrt(var2);
    int checked = 0;
    for (double rho : {-1.0, -0.99, -0.95, -0.9, -0.5, 0.0, 0.5, 0.9, 0.95, 0.99, 1.0}) {
        const double cov = rho * sd1 * sd2;
        for (double hz : {-3.0, -0.0051, 0.0, 0.8}) {
            for (double az : {0.05, 0.3923, 1.2}) {
                const double p = math::pmvnorm2dHalfRect(
                    hz * sd1, -az * sd2, az * sd2, var1, cov, var2);
                const double mX = Phi(hz);
                const double mY = Phi(az) - Phi(-az);
                CHECK(p >= -1e-12);
                CHECK(p <= mX + 1e-12);
                CHECK(p <= mY + 1e-12);
                ++checked;
            }
        }
    }
    std::printf("    %d (rho, bound) combinations bounded by both marginals\n", checked);
}

TEST(pmvnorm_reports_a_covariance_matrix_that_is_not_psd) {
    const double var1 = 17.694, var2 = 3.5031e-5;
    const double sd1 = std::sqrt(var1), sd2 = std::sqrt(var2);

    // Round-off above 1 is still absorbed: a few ULP must not turn a marker NA.
    const double covJust = (1.0 + 1e-13) * sd1 * sd2;
    const double pJust = math::pmvnorm2dHalfRect(-0.005 * sd1, -0.39 * sd2,
                                                 0.39 * sd2, var1, covJust, var2);
    CHECK(std::isfinite(pJust));

    // A genuinely indefinite pair is reported, not clamped.  1.0729 is the
    // median excess measured over the 177 such calls in examples/baseline.sh.
    for (double rho : {1.00144, 1.0729, 2.71564, -1.0729}) {
        const double cov = rho * sd1 * sd2;
        const double p = math::pmvnorm2dHalfRect(-0.005 * sd1, -0.39 * sd2,
                                                 0.39 * sd2, var1, cov, var2);
        CHECK(std::isnan(p));
    }
}

// ──────────────────────────────────────────────────────────────────────
// 9  The immaterial-leg rule
// ──────────────────────────────────────────────────────────────────────

using wtcoxg_cond::ConditionalP;
using wtcoxg_cond::conditionalP;
using wtcoxg_cond::kImmaterialLeg;
using spa::Status;

// Both legs present: the assembly is the original ratio, unchanged.
TEST(conditionalp_keeps_the_two_leg_formula_when_both_legs_are_present) {
    const double w1 = 0.3, m1 = 0.42, p1 = 0.11;
    const double w0 = 0.7, m0 = 0.9,  p0 = 0.25;
    const double pd = w1 * m1 + w0 * m0;
    const ConditionalP c = conditionalP(w1, m1, p1, Status::SpaOk,
                                        w0, m0, p0, Status::Normal, pd);
    CHECK_NEAR(c.p, 2.0 * (w1 * p1 + w0 * p0) / pd, 0.0);
    CHECK_NEAR(c.negLog10p, -std::log10(c.p), 0.0);
    // Converged outranks NormalBranch at equal severity, per spa::worseStatus.
    CHECK(c.status == Status::SpaOk);
}

// The configuration that motivated the rule, taken from the cohort: the
// sigma2 > 0 leg is weighted by TPR = 1.87e-26 and, because sigma2 = 4.02e+16
// spreads S_bat 1e10 times wider than the acceptance interval, its
// conditioning mass m1 is a further 5.5e-12.  The leg's whole contribution is
// bounded by 1.1e-37 of the answer; the marker must survive, and the reported
// value must be the mixture with that leg's joint probability set to zero.
TEST(conditionalp_keeps_a_marker_whose_missing_leg_cannot_move_the_answer) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double w1 = 1.8724e-26, m1 = 5.456e-12;
    const double w0 = 1.0 - w1,   m0 = 0.9, p0 = 0.25;
    const double pd = w1 * m1 + w0 * m0;

    const ConditionalP c = conditionalP(w1, m1, nan, Status::Normal,
                                        w0, m0, p0, Status::SpaOk, pd);
    CHECK(std::isfinite(c.p));
    CHECK_NEAR(c.p, 2.0 * w0 * p0 / pd, 0.0);
    // The status reported is the SURVIVING leg's, so the invariant that P is a
    // number exactly for statuses 0, 4 and 6 is preserved.
    CHECK(c.status == Status::SpaOk);
    CHECK(spa::statusIsUsable(c.status));

    // And the interval really is below the printed resolution: six significant
    // digits cannot separate the two ends.
    const double width = w1 * m1 / pd;
    CHECK(width / c.p < 1e-30);
}

// The same rule must refuse a leg that does carry weight.  These numbers are
// the other population measured on the same cohort: TPR = 8.55e-4 against a
// sigma2 of 3.95e-4, so the leg holds 4.4e-4 of the conditioning mass — three
// orders of magnitude above the printed resolution.
TEST(conditionalp_reports_na_when_the_missing_leg_does_move_the_answer) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double w1 = 8.54726e-4, m1 = 0.4636;
    const double w0 = 1.0 - w1,   m0 = 0.9, p0 = 0.25;
    const double pd = w1 * m1 + w0 * m0;

    const ConditionalP c = conditionalP(w1, m1, nan, Status::Normal,
                                        w0, m0, p0, Status::SpaOk, pd);
    CHECK(std::isnan(c.p));
    CHECK(std::isnan(c.negLog10p));
    // D2: the answer is NaN, never the 1.0 that std::min(1.0, NaN) returned.
    CHECK(c.p != 1.0);
    // The leg's own saddlepoint succeeded (NORMAL), so the loss happened in
    // the bivariate integral — downstream of the saddlepoint.  That is
    // NA_POST_FAIL: the substituted normal tail would answer a question about
    // a quantity that never failed.  Reporting NORMAL here would advertise a
    // usable LOG10P.
    CHECK(c.status == Status::NaPostFail);
    CHECK(!spa::statusIsUsable(c.status));
}

// The threshold is on the RELATIVE width, so the same missing leg is
// immaterial against a p-value of one half and material against a p-value in
// the tail.  This is the property that stops the rule from quietly flattening
// a genuine association down to the leg's mass.
TEST(conditionalp_threshold_is_relative_not_absolute) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double w1 = 1e-12, m1 = 1.0;
    const double w0 = 1.0,   m0 = 1.0;
    const double pd = w1 * m1 + w0 * m0;

    // p = 0.5: an interval of width 1e-12 is 2e-12 of it, well inside 1e-8.
    const ConditionalP big = conditionalP(w1, m1, nan, Status::SpaOk,
                                          w0, m0, 0.25, Status::SpaOk, pd);
    CHECK(std::isfinite(big.p));
    CHECK(w1 * m1 / pd <= kImmaterialLeg * big.p);

    // p = 2e-8: the same 1e-12 interval is now 5e-5 of the answer, and moves
    // -log10(P) in its fifth decimal.  The marker has no determined value.
    const ConditionalP tail = conditionalP(w1, m1, nan, Status::SpaOk,
                                           w0, m0, 1e-8, Status::SpaOk, pd);
    CHECK(std::isnan(tail.p));
}

// Straddle the threshold from both sides with everything else held fixed.
TEST(conditionalp_threshold_fires_where_it_says_it_does) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    // w0 = m0 = 1 and p0 = 1/2 give p_deno = 1 + eps, p = 1/(1+eps) and a
    // relative interval width of exactly eps.
    for (double eps : {0.9 * kImmaterialLeg, 1.1 * kImmaterialLeg}) {
        const double pd = eps + 1.0;
        const ConditionalP c = conditionalP(eps, 1.0, nan, Status::SpaOk,
                                            1.0, 1.0, 0.5, Status::SpaOk, pd);
        const bool immaterial = eps < kImmaterialLeg;
        CHECK(std::isfinite(c.p) == immaterial);
        if (immaterial) CHECK_NEAR(c.p, 1.0 / pd, 1e-15);
    }
}

// A leg of weight exactly zero contributes exactly zero: no threshold is
// consulted, and the marker is reported from the other leg alone.
TEST(conditionalp_drops_a_zero_weight_leg_without_a_test) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const ConditionalP c = conditionalP(0.0, 0.5, nan, Status::FallbackNonFinite,
                                        1.0, 0.9, 0.25, Status::SpaOk, 0.9);
    CHECK(std::isfinite(c.p));
    CHECK_NEAR(c.p, 2.0 * 0.25 / 0.9, 0.0);
    CHECK(c.status == Status::SpaOk);
}

// Neither leg available: there is nothing to report, and a failure status must
// come out even when both legs' saddlepoints had succeeded.
TEST(conditionalp_reports_na_when_neither_leg_is_available) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const ConditionalP c = conditionalP(0.3, 0.5, nan, Status::Normal,
                                        0.7, 0.9, nan, Status::SpaOk, 0.78);
    CHECK(std::isnan(c.p));
    CHECK(c.status == Status::NaPostFail);

    // A leg whose OWN saddlepoint had already left it without an answer keeps
    // that reason.  A fallback code does not: the saddlepoint there produced a
    // number, so the loss is downstream and NA_POST_FAIL is the accurate one.
    const ConditionalP g = conditionalP(0.3, 0.5, nan, Status::NaNoTest,
                                        0.7, 0.9, nan, Status::SpaOk, 0.78);
    CHECK(std::isnan(g.p));
    CHECK(g.status == Status::NaNoTest);            // a named reason survives
}

// D2 sweep.  A NaN leg must leave either a NaN or a value the immaterial-leg
// rule certified — never the 1.0 that `std::min(1.0, pval1 + pval2)` returned.
// q0 is the surviving leg's CONDITIONAL tail, which the mixture bounds by 1/2;
// it is kept strictly below 1/2 here so that the honest answer is strictly
// below 1 and a returned 1.0 could only have been manufactured.
TEST(conditionalp_never_manufactures_one_from_a_nan_leg) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    int finite = 0, na = 0;
    for (double w1 : {0.0, 1e-30, 1e-12, 1e-6, 1e-3, 0.5, 1.0}) {
        for (double m1 : {1e-20, 1e-6, 0.1, 0.9, 1.0}) {
            for (double q0 : {1e-12, 1e-6, 0.05, 0.3, 0.499}) {
                const double w0 = 1.0 - w1, m0 = 0.9, p0 = q0 * m0;
                const double pd = w1 * m1 + w0 * m0;
                if (!(pd > 0.0)) continue;
                const ConditionalP c =
                    conditionalP(w1, m1, nan, Status::Normal,
                                 w0, m0, p0, Status::SpaOk, pd);
                CHECK(c.p != 1.0);
                if (std::isnan(c.p)) {
                    ++na;
                    CHECK(!spa::statusIsUsable(c.status));
                } else {
                    ++finite;
                    // The reported value is the mixture with the missing leg's
                    // joint probability set to zero, exactly.
                    CHECK_NEAR(c.p, 2.0 * w0 * p0 / pd, 0.0);
                    // Certified: the dropped leg's whole interval is below the
                    // reported value times the threshold.
                    CHECK(w1 * m1 / pd <= kImmaterialLeg * c.p);
                    CHECK(spa::statusIsUsable(c.status));
                    CHECK(c.p >= 0.0 && c.p <= 1.0);
                }
            }
        }
    }
    std::printf("    %d certified-immaterial, %d NA, 0 fabricated 1.0\n", finite, na);
    CHECK(finite > 0);
    CHECK(na > 0);
}

}  // namespace

TINYTEST_MAIN
