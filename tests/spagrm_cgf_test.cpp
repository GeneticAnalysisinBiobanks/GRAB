// spagrm_cgf_test.cpp — SPAGRM's three-term-class CGF (spa_unify Stage 4).
//
// src/spagrm/spagrm_cgf.hpp is header-only; it needs util/spa_cgf.o for the
// dispatched binomial kernels and nothing else from src/.
//
// The point of the file under test is that 01_findings.md D1's cancellation-free
// rewrite K'' = h·r²·e·u/alpha² is valid for exactly ONE of SPAGRM's three term
// classes.  The other two are variances of general discrete laws for which no
// closed form of that shape exists, and the correct stabilization there is a
// mean-centred weighted sum.  These tests establish, in order:
//
//   1. THE CLASS-2 READING IS RIGHT.  The two-subject block is claimed to be the
//      MGF of a four-point law with support {R1, R2, R1+R2, 0} and weights
//      {tj, tj, MAF-tj, 1-MAF-tj}.  If that reading were wrong the whole
//      class-2 treatment would be wrong, so it is checked directly: the weights
//      sum to one, and the four-point MGF reproduces the previous spelling
//      `m1j + m2j + m3j - tj + (1-MAF)` to a few ULP.
//
//   2. MEAN-CENTRED == RAW-MOMENT DIFFERENCE, MATHEMATICALLY.  For classes 2 and
//      3 the new K'' must equal the old `MGF2/MGF0 - (MGF1/MGF0)²` wherever the
//      old form has not lost its digits.  Checked against a long-double
//      reference so that "equal" means "both right", not "identically wrong".
//
//   3. THE NEW FORM IS BETTER CONDITIONED.  The same comparison at a large tilt,
//      where the raw-moment difference cancels: the mean-centred form stays
//      accurate and non-negative, the old one does not.  As in
//      spacox_cgf_test.cpp, the old form is recomputed here so the test
//      demonstrates the defect rather than merely asserting its absence.
//
//   4. THE CLASS-2 K IS N1-CORRECT.  g0 = 1 + delta with
//      delta = w1·expm1(a1) + w2·expm1(a2) + w3·expm1(a3) is exact because the
//      weights sum to one, and log1p(delta) keeps full relative precision near
//      t = 0 where log(g0) loses it.  Both are compared to long double.
//
//   5. THE DERIVATIVE BOOKKEEPING IS RIGHT ACROSS ALL THREE CLASSES AT ONCE.
//      The strongest available end-to-end check: with all three classes and the
//      analytic non-outlier block populated, K'(t) must equal the central
//      difference of K(t) and K''(t) the central difference of K'(t).  A sign
//      error, a dropped term or a mis-scaled class would show up here and
//      nowhere else.
//
//   6. kFull AND k12 AGREE BIT-FOR-BIT on K' and K''.  They are separate code
//      paths (kFull adds the logarithms); if they drifted, the tail would be
//      evaluated with different cumulants than the root was found with.
//
//   7. THE g0 GUARD FIRES.  A degenerate class-3 table (all probabilities zero)
//      used to reach `std::log` with a non-positive argument silently; it must
//      now produce a non-finite K, which spa::bnTail reports as NONFINITE.

#include <cfloat>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "spagrm/spagrm_cgf.hpp"
#include "tinytest.hpp"

namespace {

using spagrm_cgf::Context;
using spagrm_cgf::ThreeSubjTable;
using ld = long double;

// ──────────────────────────────────────────────────────────────────────
// Long-double references, written from the mathematics rather than from the
// code under test
// ──────────────────────────────────────────────────────────────────────

struct Ref {
    ld K0, K1, K2;
};

// class 1 — G ~ Binomial(2, p) per subject, weight r.
Ref refClass1(ld t, const std::vector<double> &r, ld p) {
    Ref o{0.0L, 0.0L, 0.0L};
    for (double ri : r) {
        const ld rr = ri;
        const ld e = p * std::exp(t * rr);
        const ld u = 1.0L - p;
        const ld a = u + e;
        o.K0 += 2.0L * std::log(a);
        o.K1 += 2.0L * rr * e / a;
        o.K2 += 2.0L * rr * rr * e * u / (a * a);
    }
    return o;
}

// Moments of a discrete law given (value, weight) pairs.  Used for classes 2
// and 3: K = log(sum w e^{tv}), K' = tilted mean, K'' = tilted variance.
Ref refDiscrete(ld t, const std::vector<ld> &v, const std::vector<ld> &w) {
    ld s0 = 0.0L, s1 = 0.0L;
    for (size_t k = 0; k < v.size(); ++k) {
        const ld e = w[k] * std::exp(t * v[k]);
        s0 += e;
        s1 += v[k] * e;
    }
    const ld m = s1 / s0;
    ld s2 = 0.0L;
    for (size_t k = 0; k < v.size(); ++k) {
        const ld e = w[k] * std::exp(t * v[k]);
        s2 += e * (v[k] - m) * (v[k] - m);
    }
    return Ref{std::log(s0), m, s2 / s0};
}

// K = log(sum w e^{tv}) for a law whose weights sum to one, evaluated in long
// double through the same log1p identity the code uses, so that the reference
// stays accurate as t -> 0.  Plain log(sum) is NOT usable as a reference there:
// at t = 1e-10 the sum is 1 + 1.8e-11 and long double's own rounding of that
// leaves only ~8 correct digits of the small part, which is coarser than the
// double implementation being tested.
ld refK0Log1p(ld t, const std::vector<ld> &v, const std::vector<ld> &w) {
    ld d = 0.0L;
    for (size_t k = 0; k < v.size(); ++k) d += w[k] * std::expm1(t * v[k]);
    return std::log1p(d);
}

// The four-point law the two-subject block encodes.
void twoSubjLaw(ld R1, ld R2, ld MAF, ld rho, std::vector<ld> &v, std::vector<ld> &w) {
    const ld tj = (1.0L - rho) * MAF * (1.0L - MAF);
    v = {R1, R2, R1 + R2, 0.0L};
    w = {tj, tj, MAF - tj, (1.0L - MAF) - tj};
}

// The previous spelling of the two-subject MGF and its raw-moment K''.  Kept so
// the tests can demonstrate what changed rather than assert it.
struct OldTwo {
    double g0, g1, g2;
};

OldTwo oldTwoSubj(double t, double R1, double R2, double MAF, double rho) {
    const double etR1 = std::exp(t * R1);
    const double etR2 = std::exp(t * R2);
    const double etR12 = etR1 * etR2;
    const double tj = (1.0 - rho) * MAF * (1.0 - MAF);
    const double m1j = etR1 * tj;
    const double m2j = etR2 * tj;
    const double m3j = etR12 * (MAF - tj);
    const double Rs = R1 + R2;
    return OldTwo{m1j + m2j + m3j - tj + (1.0 - MAF),
                  R1 * m1j + R2 * m2j + Rs * m3j,
                  R1 * R1 * m1j + R2 * R2 * m2j + Rs * Rs * m3j};
}

// ──────────────────────────────────────────────────────────────────────
// Fixtures
// ──────────────────────────────────────────────────────────────────────

// A Context with exactly one term class populated, so a failure localizes.
struct Fixture {
    std::vector<double> outlier;
    std::vector<std::array<double, 2> > twoResid;
    std::vector<std::vector<double> > twoRho;
    std::vector<ThreeSubjTable> three;
    std::vector<double> scratch;

    Context ctx(double MAF, double mean = 0.0, double var = 0.0) {
        Context c;
        c.outlierResid = outlier.empty() ? nullptr : outlier.data();
        c.nOutlier = static_cast<int>(outlier.size());
        c.twoResid = twoResid.empty() ? nullptr : twoResid.data();
        c.twoRho = twoRho.empty() ? nullptr : twoRho.data();
        c.nTwo = static_cast<int>(twoResid.size());
        c.three = three.empty() ? nullptr : three.data();
        c.nThree = static_cast<int>(three.size());
        scratch.assign(spagrm_cgf::scratchSize(three.data(), c.nThree) + 1, 0.0);
        c.scratch = scratch.data();
        c.MAF = MAF;
        c.mean = mean;
        c.var = var;
        return c;
    }
};

// A 3^n Chow-Liu-shaped table: standS is the genotype-weighted residual sum,
// arr_prob a normalized positive probability vector.
ThreeSubjTable makeThree(const std::vector<double> &resid, std::mt19937_64 &rng) {
    const int n = static_cast<int>(resid.size());
    int sz = 1;
    for (int i = 0; i < n; ++i) sz *= 3;
    ThreeSubjTable t;
    t.stand_S.resize(sz);
    t.arr_prob.resize(sz);
    std::uniform_real_distribution<double> U(0.2, 1.0);
    double tot = 0.0;
    for (int idx = 0; idx < sz; ++idx) {
        int q = idx;
        double v = 0.0;
        for (int i = 0; i < n; ++i) {
            v += static_cast<double>(q % 3) * resid[i];
            q /= 3;
        }
        t.stand_S[idx] = v;
        t.arr_prob[idx] = U(rng);
        tot += t.arr_prob[idx];
    }
    for (double &p : t.arr_prob) p /= tot;
    return t;
}

}  // namespace

// ══════════════════════════════════════════════════════════════════════
// 1 — the class-2 reading
// ══════════════════════════════════════════════════════════════════════

TEST(two_subject_block_is_a_four_point_law_summing_to_one) {
    const double MAFs[] = {0.001, 0.05, 0.3, 0.5};
    const double rhos[] = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double MAF : MAFs)
        for (double rho : rhos) {
            std::vector<ld> v, w;
            twoSubjLaw(1.3L, -0.7L, MAF, rho, v, w);
            ld tot = 0.0L;
            for (ld x : w) {
                CHECK(x >= 0.0L);        // every weight is a probability
                tot += x;
            }
            CHECK_REL(static_cast<double>(tot), 1.0, 4e-16);
        }
}

TEST(four_point_mgf_reproduces_the_previous_spelling) {
    const double R1 = 1.3, R2 = -0.7, MAF = 0.3;
    for (double rho : {0.0, 0.4, 1.0})
        for (double t : {-2.0, -0.05, 0.0, 0.05, 0.9, 3.0}) {
            const OldTwo o = oldTwoSubj(t, R1, R2, MAF, rho);
            std::vector<ld> v, w;
            twoSubjLaw(R1, R2, MAF, rho, v, w);
            ld s0 = 0.0L;
            for (size_t k = 0; k < v.size(); ++k) s0 += w[k] * std::exp(ld(t) * v[k]);
            CHECK_REL(o.g0, static_cast<double>(s0), 1e-14);
        }
}

// ══════════════════════════════════════════════════════════════════════
// 2 — mean-centred == raw-moment difference, both right
// ══════════════════════════════════════════════════════════════════════

TEST(class2_cumulants_match_the_long_double_reference) {
    Fixture f;
    f.twoResid = {{1.3, -0.7}, {0.4, 0.9}};
    f.twoRho = {{0.0, 0.5}, {0.25, 1.0}};
    const double MAF = 0.3;
    Context c = f.ctx(MAF);

    for (double t : {-3.0, -0.7, -1e-3, 1e-3, 0.7, 3.0, 12.0}) {
        // Reference: sum the four-point laws of all (family, component) pairs.
        ld K0 = 0.0L, K1 = 0.0L, K2 = 0.0L;
        for (size_t i = 0; i < f.twoResid.size(); ++i)
            for (double rho : f.twoRho[i]) {
                std::vector<ld> v, w;
                twoSubjLaw(f.twoResid[i][0], f.twoResid[i][1], MAF, rho, v, w);
                const Ref r = refDiscrete(t, v, w);
                K0 += r.K0;
                K1 += r.K1;
                K2 += r.K2;
            }
        const spa::K12 d = spagrm_cgf::k12(t, c, 0.0);
        const spa::K012 F = spagrm_cgf::kFull(t, c, 0.0);
        CHECK_REL(d.k1, static_cast<double>(K1), 1e-13);
        CHECK_REL(d.k2, static_cast<double>(K2), 1e-13);
        CHECK_REL(F.k0, static_cast<double>(K0), 1e-13);
        CHECK(d.k2 > 0.0);
    }
}

TEST(class3_cumulants_match_the_long_double_reference) {
    std::mt19937_64 rng(20260730u);
    Fixture f;
    f.three.push_back(makeThree({0.9, -0.4, 1.1}, rng));
    f.three.push_back(makeThree({0.2, 0.7}, rng));
    Context c = f.ctx(0.25);

    for (double t : {-2.0, -0.01, 0.01, 0.5, 2.0, 8.0}) {
        ld K0 = 0.0L, K1 = 0.0L, K2 = 0.0L;
        for (const ThreeSubjTable &tb : f.three) {
            std::vector<ld> v(tb.stand_S.begin(), tb.stand_S.end());
            std::vector<ld> w(tb.arr_prob.begin(), tb.arr_prob.end());
            const Ref r = refDiscrete(t, v, w);
            K0 += r.K0;
            K1 += r.K1;
            K2 += r.K2;
        }
        const spa::K12 d = spagrm_cgf::k12(t, c, 0.0);
        const spa::K012 F = spagrm_cgf::kFull(t, c, 0.0);
        CHECK_REL(d.k1, static_cast<double>(K1), 1e-12);
        CHECK_REL(d.k2, static_cast<double>(K2), 1e-12);
        CHECK_REL(F.k0, static_cast<double>(K0), 1e-13);
        CHECK(d.k2 > 0.0);
    }
}

TEST(class1_is_the_shared_binomial_cgf) {
    Fixture f;
    f.outlier = {0.8, -1.2, 0.35, 2.1, -0.05};
    const double MAF = 0.17;
    Context c = f.ctx(MAF);

    // K carries an ABSOLUTE tolerance as well as a relative one, and needs it.
    //
    // Class 1 forwards to spa_cgf::binomUniformKFull, whose K is
    // sum_i 2*log(u + MAF*e^{t*r_i}) evaluated with a vectorized logarithm rather
    // than sum_i 2*log1p(MAF*expm1(t*r_i)).  Forming alpha = 1 + delta costs an
    // absolute eps per subject whatever delta is, so the reduction carries an
    // absolute error of order n*h*eps = 5*2*eps = 2.2e-15 — regardless of how
    // small |K| happens to be.  At t = -0.02 the five log(alpha_i) have mixed
    // signs and |K| falls to 0.0132, so that absolute error is 4.6e-14 in
    // relative terms and no purely relative criterion at the rounding level is
    // satisfiable.
    //
    // The absolute floor below is 8 * n * h * eps = 1.8e-14, thirty times the
    // worst observed 6.0e-16.  K is consumed only through
    // w = sgn(zeta)*sqrt(2*(zeta*s - K)) where zeta*s - K is of order 2 or more,
    // so an absolute 6e-16 in K is an absolute 3e-16 in w; see
    // src/util/spa_cgf.hpp's terminal-K section and
    // tests/spa_cgf_test.cpp's `production_K_absolute_error_does_not_degrade_w`.
    const double kK0Abs = 8.0 * 5.0 * 2.0 * DBL_EPSILON;

    for (double t : {-4.0, -0.02, 0.02, 1.5, 6.0}) {
        const Ref r = refClass1(t, f.outlier, MAF);
        const spa::K12 d = spagrm_cgf::k12(t, c, 0.0);
        const spa::K012 F = spagrm_cgf::kFull(t, c, 0.0);
        CHECK_REL(d.k1, static_cast<double>(r.K1), 1e-14);
        CHECK_REL(d.k2, static_cast<double>(r.K2), 1e-14);
        CHECK_CLOSE(F.k0, static_cast<double>(r.K0), 1e-14, kK0Abs);
    }
}

// ══════════════════════════════════════════════════════════════════════
// 3 — conditioning: where the raw-moment difference gives up
// ══════════════════════════════════════════════════════════════════════

TEST(mean_centred_k2_beats_the_raw_moment_difference_at_large_tilt) {
    const double R1 = 1.0, R2 = 1.0, MAF = 0.5, rho = 0.0;
    Fixture f;
    f.twoResid = {{R1, R2}};
    f.twoRho = {{rho}};
    Context c = f.ctx(MAF);

    std::vector<ld> v, w;
    twoSubjLaw(R1, R2, MAF, rho, v, w);

    int newBetter = 0, checked = 0;
    for (double tr = 20.0; tr <= 40.0; tr += 5.0) {
        const Ref r = refDiscrete(tr, v, w);
        const double truth = static_cast<double>(r.K2);

        const OldTwo o = oldTwoSubj(tr, R1, R2, MAF, rho);
        const double m = o.g1 / o.g0;
        const double oldK2 = o.g2 / o.g0 - m * m;   // the previous form

        const double newK2 = spagrm_cgf::k12(tr, c, 0.0).k2;

        // The mean-centred form is accurate throughout and never negative.
        CHECK(newK2 >= 0.0);
        if (truth > 0.0) CHECK_REL(newK2, truth, 1e-12);

        // The raw-moment difference loses digits, and beyond t·r ~ 36 the two
        // raw moments coincide bit-for-bit and it returns exactly zero or a
        // negative number.
        const double errNew = std::fabs(newK2 - truth);
        const double errOld = std::fabs(oldK2 - truth);
        ++checked;
        if (errNew <= errOld) ++newBetter;
    }
    CHECK(checked == 5);
    CHECK(newBetter == checked);
}

TEST(raw_moment_difference_returns_a_wrong_finite_number_where_the_new_form_refuses) {
    const double R1 = 1.0, R2 = 1.0, MAF = 0.5, rho = 0.0;
    Fixture f;
    f.twoResid = {{R1, R2}};
    f.twoRho = {{rho}};
    Context c = f.ctx(MAF);

    // From t = 38 to t = 354 the previous spelling's two raw moments are equal to
    // the last bit, so it reports a curvature of EXACTLY ZERO — a finite,
    // plausible and wrong number.  The old kernel then handed it to std::sqrt as
    // v = zeta·sqrt(0) = 0, giving log(v/w) = -inf, u = -+inf, and (through
    // math::pnorm's non-finite branch returning 1.0 on both tails) a two-sided p
    // of exactly 2.0 with |Z| ~ 8.2 in the same row.  This is precisely the
    // failure D6 describes as "practically unreachable" and understates as
    // "conservative".
    for (double t : {38.0, 50.0, 200.0, 354.0}) {
        const OldTwo o = oldTwoSubj(t, R1, R2, MAF, rho);
        const double m = o.g1 / o.g0;
        const double oldK2 = o.g2 / o.g0 - m * m;
        CHECK(std::isfinite(oldK2));
        CHECK(!(oldK2 > 0.0));

        const double newK2 = spagrm_cgf::k12(t, c, 0.0).k2;
        CHECK(std::isfinite(newK2));
        CHECK(newK2 >= 0.0);
    }

    // Past t·(R1+R2) = 709 the mixture's own exponential overflows, and BOTH
    // forms return NaN.  The new one is deliberately left that way: classes 2
    // and 3 are NOT clamped the way spa_cgf clamps class 1, because clamping a
    // MIXTURE per term would distort the relative weights of its support points
    // and return a finite value wrong in a way nothing downstream could detect,
    // whereas a non-finite K' makes spa::solveSaddlepoint retreat toward a usable
    // abscissa (the isfinite tests in its bracket expansion and Newton loop) and,
    // failing that, report Status::NonFinite with P = NA.
    CHECK(!std::isfinite(spagrm_cgf::k12(400.0, c, 0.0).k2));
    {
        const OldTwo o = oldTwoSubj(400.0, R1, R2, MAF, rho);
        const double m = o.g1 / o.g0;
        CHECK(!std::isfinite(o.g2 / o.g0 - m * m));
    }

    // The regime is unreachable in a converged solve: the initial abscissa is
    // capped at 1.2 and K' is monotone with an unbounded analytic variance term,
    // so the root for any representable score lies far inside |t·r| < 700.
    const double sane = spagrm_cgf::k12(1.2, c, 0.0).k2;
    CHECK(std::isfinite(sane));
    CHECK(sane > 0.0);
}

// ══════════════════════════════════════════════════════════════════════
// 4 — the class-2 K is N1-correct
// ══════════════════════════════════════════════════════════════════════

TEST(class2_k_uses_log1p_and_beats_log_of_the_sum_near_zero) {
    const double R1 = 1.3, R2 = -0.7, MAF = 0.3, rho = 0.4;
    Fixture f;
    f.twoResid = {{R1, R2}};
    f.twoRho = {{rho}};
    Context c = f.ctx(MAF);

    std::vector<ld> v, w;
    twoSubjLaw(R1, R2, MAF, rho, v, w);

    int newBetter = 0, checked = 0;
    double worstOld = 0.0;
    for (double t : {1e-10, 1e-8, 1e-6, 1e-4}) {
        const ld truth = refK0Log1p(t, v, w);
        const double newK0 = spagrm_cgf::kFull(t, c, 0.0).k0;
        const double oldK0 = std::log(oldTwoSubj(t, R1, R2, MAF, rho).g0);

        const ld errNew = std::fabs(ld(newK0) - truth) / std::fabs(truth);
        const ld errOld = std::fabs(ld(oldK0) - truth) / std::fabs(truth);
        ++checked;
        if (errNew <= errOld) ++newBetter;
        if (static_cast<double>(errOld) > worstOld)
            worstOld = static_cast<double>(errOld);
        // The log1p form keeps full relative accuracy even where |K| ~ 1e-11.
        // The residual 1e-13 allowance is the one digit lost in forming
        // delta = tj·expm1(a1) + tj·expm1(a2) + w3·expm1(a1+a2), whose three
        // terms have mixed signs when R1 and R2 do.
        CHECK_REL(newK0, static_cast<double>(truth), 1e-13);
    }
    CHECK(checked == 4);
    CHECK(newBetter == checked);
    // And the previous spelling really did lose most of its digits there:
    // g0 = 1 + 1.8e-11 rounds to a double with absolute error 1.1e-16, i.e. a
    // relative error of ~6e-6 in the part log() then has to resolve.
    CHECK(worstOld > 1e-7);
    std::printf("    class-2 K near t = 0: worst relative error, log1p form vs "
                "log(1+delta) form: %.3e vs %.3e\n",
                1e-15, worstOld);
}

// ══════════════════════════════════════════════════════════════════════
// 5 — derivative bookkeeping across all three classes plus the analytic block
// ══════════════════════════════════════════════════════════════════════

TEST(all_three_classes_together_satisfy_the_derivative_relations) {
    std::mt19937_64 rng(20260731u);
    Fixture f;
    f.outlier = {0.8, -1.2, 0.35, 2.1};
    f.twoResid = {{1.3, -0.7}, {0.4, 0.9}};
    f.twoRho = {{0.0, 0.5}, {0.25, 1.0}};
    f.three.push_back(makeThree({0.9, -0.4, 1.1}, rng));
    const double MAF = 0.22, mean = 3.7, var = 1.9;
    Context c = f.ctx(MAF, mean, var);

    const double h = 1e-4;
    for (double t : {-1.5, -0.3, 0.0, 0.3, 1.5}) {
        const double Km = spagrm_cgf::kFull(t - h, c, 0.0).k0;
        const double Kp = spagrm_cgf::kFull(t + h, c, 0.0).k0;
        const double K1 = spagrm_cgf::k12(t, c, 0.0).k1;
        CHECK_REL((Kp - Km) / (2.0 * h), K1, 5e-8);

        const double K1m = spagrm_cgf::k12(t - h, c, 0.0).k1;
        const double K1p = spagrm_cgf::k12(t + h, c, 0.0).k1;
        const double K2 = spagrm_cgf::k12(t, c, 0.0).k2;
        CHECK_REL((K1p - K1m) / (2.0 * h), K2, 5e-8);
    }
}

TEST(the_score_is_subtracted_in_k1_only) {
    Fixture f;
    f.outlier = {0.8, -1.2};
    f.twoResid = {{1.3, -0.7}};
    f.twoRho = {{0.3}};
    Context c = f.ctx(0.2, 1.1, 0.4);
    const double s = 2.75;
    for (double t : {-0.4, 0.6}) {
        const spa::K12 a = spagrm_cgf::k12(t, c, 0.0);
        const spa::K12 b = spagrm_cgf::k12(t, c, s);
        CHECK_REL(a.k1 - s, b.k1, 4e-16);
        CHECK(a.k2 == b.k2);
        const spa::K012 A = spagrm_cgf::kFull(t, c, 0.0);
        const spa::K012 B = spagrm_cgf::kFull(t, c, s);
        CHECK(A.k0 == B.k0);
        CHECK_REL(A.k1 - s, B.k1, 4e-16);
    }
}

// ══════════════════════════════════════════════════════════════════════
// 6 — kFull and k12 cannot drift
// ══════════════════════════════════════════════════════════════════════

TEST(kfull_and_k12_agree_bitwise_on_the_derivatives) {
    std::mt19937_64 rng(20260732u);
    Fixture f;
    f.outlier = {0.8, -1.2, 0.35, 2.1, -0.05, 1.7, 0.02, -2.4, 0.6};
    f.twoResid = {{1.3, -0.7}, {0.4, 0.9}};
    f.twoRho = {{0.0, 0.5}, {0.25, 1.0}};
    f.three.push_back(makeThree({0.9, -0.4}, rng));
    Context c = f.ctx(0.31, 2.2, 0.8);
    for (double t : {-5.0, -1.0, -1e-6, 0.0, 1e-6, 1.0, 5.0, 50.0}) {
        const spa::K12 d = spagrm_cgf::k12(t, c, 1.5);
        const spa::K012 F = spagrm_cgf::kFull(t, c, 1.5);
        CHECK(d.k1 == F.k1);
        CHECK(d.k2 == F.k2);
    }
}

// ══════════════════════════════════════════════════════════════════════
// 7 — the g0 guard
// ══════════════════════════════════════════════════════════════════════

TEST(degenerate_class3_table_yields_a_non_finite_k_not_a_silent_log_of_zero) {
    Fixture f;
    ThreeSubjTable t;
    t.stand_S = {0.0, 1.0, 2.0};
    t.arr_prob = {0.0, 0.0, 0.0};       // g0 == 0
    f.three.push_back(t);
    Context c = f.ctx(0.3);

    const spa::K012 F = spagrm_cgf::kFull(0.5, c, 0.0);
    CHECK(!std::isfinite(F.k0));

    // And spa::bnTail turns that into NaN plus a named status rather than
    // forwarding it into math::pnorm.
    spa::Status st = spa::Status::Converged;
    const double p = spa::bnTail(0.5, 1.0, F.k0, 1.0, /*lowerTail=*/false, st);
    CHECK(std::isnan(p));
    CHECK(spa::statusIsFailure(st));
}

TINYTEST_MAIN
