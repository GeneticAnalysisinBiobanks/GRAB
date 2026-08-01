// log10p_test.cpp — The log-domain distribution tier of the log10p_unify
// project.
//
// ══════════════════════════════════════════════════════════════════════
// What this suite is for
// ══════════════════════════════════════════════════════════════════════
//
// dev-notes/methods/log10p_unify/ makes L = -log10(P) the sole p-value
// representation.  Stage 1 adds the five functions that make that possible to
// src/util/math_helper.{hpp,cpp} without changing a single call site, so the
// binary's output is byte-identical and the ONLY evidence that the new code is
// right is this file.  Every assertion below is therefore a cross-check
// against something that is not the code under test: the linear spelling it
// will replace (wherever that spelling still works), a long-double evaluation,
// an analytic special case, or the measured table in 01_numerics.md.
//
// Sections, in the order the functions appear in math_helper.hpp:
//
//   § 0  the -log10 convention itself
//   § 1  zFromNegLog10P
//   § 2  chisq1FromNegLog10P
//   § 3  cauchyCombineLog10
//   § 4  ptLog
//   § 5  pmvnorm2dHalfRectLog
//
// Stage 9 appends its HweLnP cross-check as a new section immediately before
// TINYTEST_MAIN at the end of the file; nothing here has to move for it.
//
// Build:  make test   (or: make build/tests/log10p_test && ./build/tests/log10p_test)

#include "tinytest.hpp"

#include "util/math_helper.hpp"
#include "util/spa.hpp"

#include <boost/math/special_functions/beta.hpp>

#include <cmath>
#include <limits>
#include <vector>

namespace {

constexpr double kLn10 = 2.30258509299404568402;
constexpr double kLn2  = 0.69314718055994530942;

const double kInf = std::numeric_limits<double>::infinity();
const double kNaN = std::numeric_limits<double>::quiet_NaN();

// ── long-double references ─────────────────────────────────────────────
//
// long double on x86-64 carries a 64-bit significand and an exponent range to
// 1e4932, so a probability that is subnormal or zero in double is an ordinary
// number here.  That is what makes these usable as references in exactly the
// regime the functions under test exist for.

long double lnPhiL(long double x) {
    if (x > -40.0L) return std::log(0.5L * std::erfc(-x / std::sqrt(2.0L)));
    const long double y = 1.0L / (x * x);
    const long double s =
        y * (-1.0L + y * (3.0L + y * (-15.0L + y * (105.0L + y * (-945.0L + y * 10395.0L)))));
    return -0.5L * x * x - std::log(-x)
           - 0.918938533204672741780329736405617639L + std::log1p(s);
}

// -log10 of the standard Cauchy upper tail at T, evaluated in long double.
long double cauchyUpperNegLog10L(long double T) {
    const long double pi = 3.14159265358979323846264338327950288L;
    const long double tail = (T > 0.0L) ? std::atan(1.0L / T) / pi
                                        : 0.5L - std::atan(T) / pi;
    return -std::log10(tail);
}

// The whole Cauchy combination in long double, formed the direct way: the
// statistic T = (1/n)·Σ cot(π pᵢ) is built explicitly, which double cannot do
// past pᵢ ≈ 1e-308 but long double can carry to pᵢ ≈ 1e-4900.
long double cauchyCombineNegLog10RefL(const double *L, int n) {
    const long double pi = 3.14159265358979323846264338327950288L;
    const long double clampL = -std::log10(0.999L);
    long double T = 0.0L;
    int nValid = 0;
    for (int i = 0; i < n; ++i) {
        if (std::isnan(L[i])) continue;
        ++nValid;
        const long double Li = std::fmax(static_cast<long double>(L[i]), clampL);
        const long double p  = std::pow(10.0L, -Li);
        T += std::cos(pi * p) / std::sin(pi * p);
    }
    if (nValid == 0) return std::numeric_limits<long double>::quiet_NaN();
    return cauchyUpperNegLog10L(T / static_cast<long double>(nValid));
}

// ln of the conditional tail integral that pmvnorm2dHalfRectLog evaluates,
// by composite Simpson in long double on a window located by scanning.  A
// different rule, a different precision and a different way of finding the
// peak than the routine under test uses.
long double lnCondTailRefL(long double h, long double rho, long double lo, long double hi) {
    const long double c = 1.0L / std::sqrt((1.0L - rho) * (1.0L + rho));
    const auto M = [&](long double u) {
        return -0.5L * u * u - 0.918938533204672741780329736405617639L
               + lnPhiL(c * (h - rho * u));
    };
    // Peak by golden section over an expanding bracket.
    long double a = -1.0L, b = 1.0L;
    for (int i = 0; i < 300; ++i) {
        if (M(a) < M(0.5L * a) && M(b) < M(0.5L * b)) break;
        a *= 2.0L; b *= 2.0L;
    }
    const long double gr = 0.6180339887498948482045868343656381L;
    long double x1 = b - gr * (b - a), x2 = a + gr * (b - a);
    long double f1 = M(x1), f2 = M(x2);
    for (int i = 0; i < 400; ++i) {
        if (f1 < f2) { a = x1; x1 = x2; f1 = f2; x2 = a + gr * (b - a); f2 = M(x2); }
        else         { b = x2; x2 = x1; f2 = f1; x1 = b - gr * (b - a); f1 = M(x1); }
        if (b - a < 1e-19L * (std::fabs(a) + std::fabs(b) + 1e-30L)) break;
    }
    const long double uPeak = std::fmin(std::fmax(0.5L * (a + b), lo), hi);
    const long double mPeak = M(uPeak);
    const auto drop = [&](int dir) {
        long double s0 = 0.0L, s1 = 1e-14L;
        for (int i = 0; i < 400 && s1 < 1e6L; ++i) {
            if (M(uPeak + dir * s1) < mPeak - 200.0L) break;
            s0 = s1; s1 *= 2.0L;
        }
        for (int i = 0; i < 400; ++i) {
            const long double mid = 0.5L * (s0 + s1);
            if (mid == s0 || mid == s1) break;
            if (M(uPeak + dir * mid) < mPeak - 200.0L) s1 = mid; else s0 = mid;
        }
        return s1;
    };
    const long double A = std::fmax(lo, uPeak - drop(-1));
    const long double B = std::fmin(hi, uPeak + drop(+1));
    if (!(B > A)) return -std::numeric_limits<long double>::infinity();

    const int N = 100000;                    // even, so Simpson is well posed
    const long double step = (B - A) / static_cast<long double>(N);
    long double sum = 0.0L;
    for (int i = 0; i <= N; ++i) {
        const long double u = A + step * static_cast<long double>(i);
        const long double w = (i == 0 || i == N) ? 1.0L : ((i % 2) ? 4.0L : 2.0L);
        sum += w * std::exp(M(u) - mPeak);
    }
    return mPeak + std::log(sum * step / 3.0L);
}

} // namespace

// ══════════════════════════════════════════════════════════════════════
// § 0  The convention
// ══════════════════════════════════════════════════════════════════════

// The convention decision D2 fixes for every p-value column this project
// produces: LOG10P is -log10(P), not log10(P) and not ln(P).
TEST(log10p_is_minus_log10_of_p) {
    CHECK_NEAR(-std::log10(1e-8), 8.0, 1e-15);
    CHECK_NEAR(-std::log10(1.0), 0.0, 0.0);
    // A p-value of 1 is LOG10P = 0, and no valid p-value yields a negative
    // LOG10P.  This is invariant C2 of 04_validation.md §2, checked here on the
    // convention and by tests/regress.py on every emitted column.
    CHECK(-std::log10(1.0) >= 0.0);
}

// ══════════════════════════════════════════════════════════════════════
// § 1  zFromNegLog10P
// ══════════════════════════════════════════════════════════════════════

// The acceptance point of the whole project on the Z column.  01_numerics §3.1
// tabulates the correct z against what math::zFromPval returns today: the two
// agree up to L = 299.698970 = -log10(2e-300), which is where qnorm's argument
// clamp engages, and from there on zFromPval is frozen at 37.0470962993612
// while the true z keeps rising.  L = 300 is already past the clamp, so this
// is not an exotic corner: any marker a GWAS would call the top hit is in it.
TEST(zFromNegLog10P_absolute_pins) {
    CHECK_REL(math::zFromNegLog10P(299.69897, 1.0), 37.0470962990919, 1e-14);
    CHECK_REL(math::zFromNegLog10P(300.0,     1.0), 37.0657878807721, 1e-14);
    CHECK_REL(math::zFromNegLog10P(400.0,     1.0), 42.8264064911712, 1e-14);

    // ... and the same three, showing that the linear route is stuck.
    CHECK_REL(math::zFromPval(std::pow(10.0, -299.69897), 1.0), 37.0470962990919, 1e-14);
    CHECK_REL(math::zFromPval(std::pow(10.0, -300.0),     1.0), 37.0470962993612, 1e-14);
    CHECK_REL(math::zFromPval(std::pow(10.0, -400.0),     1.0), 37.0470962993612, 1e-14);
    // The clamp is a plateau, not a rounding: the two saturated values above
    // are bit-identical to each other although the p-values differ by 100
    // orders of magnitude.
    CHECK(math::zFromPval(std::pow(10.0, -300.0), 1.0) ==
          math::zFromPval(std::pow(10.0, -400.0), 1.0));
    // The replacement is strictly increasing across the same span.
    CHECK(math::zFromNegLog10P(300.0, 1.0) < math::zFromNegLog10P(301.0, 1.0));
    CHECK(math::zFromNegLog10P(301.0, 1.0) < math::zFromNegLog10P(400.0, 1.0));
}

// Round trip against the tail it inverts.  spa::normalTwoSidedLog(z) is
// ln(2*Phi(-|z|)) as the saddlepoint tier's normal branch reports it, so this
// closes the loop the Z column will actually travel: LOG10P is produced from a
// z, and Stage 4 will produce a z back from that LOG10P.
TEST(zFromNegLog10P_round_trip_against_normalTwoSidedLog) {
    for (double lz = -1.0; lz <= 6.0001; lz += 0.05) {
        const double z = std::pow(10.0, lz);
        const double L = -spa::normalTwoSidedLog(z) / kLn10;
        CHECK_REL(math::zFromNegLog10P(L, 1.0), z, 1e-14);
        CHECK_REL(math::zFromNegLog10P(L, -1.0), -z, 1e-14);
    }
    // And far past anything a finite sample can mean, where the point is only
    // that nothing saturates or overflows: 01_numerics §1.4 puts the
    // representational ceiling at L ~ 3.9e307.
    for (double L : {1e3, 1e5, 1e7, 1e100, 3.8e307}) {
        const double z = math::zFromNegLog10P(L, 1.0);
        CHECK(std::isfinite(z));
        CHECK_REL(-spa::normalTwoSidedLog(z) / kLn10, L, 1e-14);
    }
}

// Where the linear route is valid the two must agree; the sweep runs up to
// L = 299.6, the last L at which 0.5*P is still above qnorm's 1e-300 clamp.
TEST(zFromNegLog10P_matches_zFromPval_below_the_clamp) {
    for (double L = 0.005; L < 299.6; L *= 1.03) {
        const double p = std::pow(10.0, -L);
        CHECK_REL(math::zFromNegLog10P(L, 1.0), math::zFromPval(p, 1.0), 1e-13);
    }
    // The Boost-quantile / Newton seam at L = 1 is continuous: both sides are
    // accurate there, which is the same overlap argument that fixes pnormLog's
    // own branch at |t| = 37.
    for (double L = 0.90; L <= 1.1001; L += 0.005) {
        const double p = std::pow(10.0, -L);
        CHECK_REL(math::zFromNegLog10P(L, 1.0), math::zFromPval(p, 1.0), 1e-14);
    }
}

TEST(zFromNegLog10P_edge_cases) {
    CHECK(std::isnan(math::zFromNegLog10P(kNaN, 1.0)));
    CHECK(std::isnan(math::zFromNegLog10P(5.0, kNaN)));
    // P >= 1 is z = 0, with the sign carried through exactly as zFromPval did.
    CHECK(math::zFromNegLog10P(0.0, 1.0) == 0.0);
    CHECK(math::zFromNegLog10P(-1.0, 1.0) == 0.0);
    CHECK(math::zFromNegLog10P(0.0, -1.0) == 0.0);
    // P = 0 is the only input for which an infinite z is the right answer.
    CHECK(math::zFromNegLog10P(kInf, 1.0) == kInf);
    CHECK(math::zFromNegLog10P(kInf, -1.0) == -kInf);
    // Sign comes from the score, magnitude from the p-value.
    CHECK(math::zFromNegLog10P(8.0, -0.001) < 0.0);
    CHECK(math::zFromNegLog10P(8.0, 0.0) > 0.0);
    CHECK_REL(math::zFromNegLog10P(8.0, -3.0), -math::zFromNegLog10P(8.0, 3.0), 0.0);
}

// ══════════════════════════════════════════════════════════════════════
// § 2  chisq1FromNegLog10P
// ══════════════════════════════════════════════════════════════════════

// P(chisq_1 > q) = 2*Phi(-sqrt(q)) is an identity, not an approximation, so
// the df = 1 quantile is exactly the square of the two-sided normal quantile
// and needs no numerical inversion of its own.
TEST(chisq1FromNegLog10P_is_the_square_of_the_normal_quantile) {
    for (double L = 0.01; L < 1e6; L *= 1.7) {
        const double z = math::zFromNegLog10P(L, 1.0);
        CHECK_REL(math::chisq1FromNegLog10P(L), z * z, 0.0);
    }
}

// Against the Boost quantile it replaces, over the range where that one is
// not clamped.
TEST(chisq1FromNegLog10P_matches_qchisq_below_the_clamp) {
    for (double L = 0.01; L < 299.0; L *= 1.05) {
        const double p = std::pow(10.0, -L);
        CHECK_REL(math::chisq1FromNegLog10P(L), math::qchisq(p, 1.0, false, false), 1e-12);
    }
}

// The saturation this removes: qchisq clamps its argument at 1e-300 and so
// tops out at 1373.87, which is the ceiling CLAUDE.md records on the
// per-cluster weight in metaPvalueScorePool and on the variance WtCoxG
// recovers by inverting a saddlepoint p-value.
TEST(chisq1FromNegLog10P_removes_the_1373_ceiling) {
    CHECK_REL(math::qchisq(std::pow(10.0, -300.0), 1.0, false, false), 1373.87, 1e-4);
    CHECK(math::qchisq(std::pow(10.0, -300.0), 1.0, false, false) ==
          math::qchisq(std::pow(10.0, -900.0), 1.0, false, false));
    CHECK(math::chisq1FromNegLog10P(300.0) > 1373.87);
    CHECK(math::chisq1FromNegLog10P(900.0) > 3.0 * math::chisq1FromNegLog10P(300.0));
    // 2*Phi(-sqrt(q)) = 10^(-L) read back, in the log domain.
    for (double L : {300.0, 900.0, 5000.0}) {
        const double q = math::chisq1FromNegLog10P(L);
        CHECK_REL(-spa::normalTwoSidedLog(std::sqrt(q)) / kLn10, L, 1e-14);
    }
}

// ══════════════════════════════════════════════════════════════════════
// § 3  cauchyCombineLog10
// ══════════════════════════════════════════════════════════════════════

// The eight cases measured in 01_numerics §2.3.  Cases 4-6 are the point of
// the exercise: the linear cauchyCombine returns exactly 0 there, so its
// -log10 is +Inf, which is invariant C1's forbidden value.
TEST(cauchyCombineLog10_the_eight_reference_cases) {
    const double c1[] = {0.3, 1.0, 2.0};
    const double c2[] = {8.0, 1.0, 0.5};
    const double c3[] = {300.0, 1.0, 0.5};
    const double c4[] = {320.0, 1.0, 0.5};
    const double c5[] = {400.0, 1.0, 0.5};
    const double c6[] = {400.0, 400.0};
    const double c7[] = {0.1, 0.05, 0.2};
    const double c8[] = {8.0, kNaN};

    CHECK_REL(math::cauchyCombineLog10(c1, 3), 1.56385183,     1e-8);
    CHECK_REL(math::cauchyCombineLog10(c2, 3), 7.5228788,      1e-7);
    CHECK_REL(math::cauchyCombineLog10(c3, 3), 299.522879,     1e-8);
    CHECK_REL(math::cauchyCombineLog10(c4, 3), 319.522879,     1e-8);
    CHECK_REL(math::cauchyCombineLog10(c5, 3), 399.522879,     1e-8);
    CHECK_REL(math::cauchyCombineLog10(c6, 2), 400.0,          1e-14);
    CHECK_REL(math::cauchyCombineLog10(c7, 3), 0.0887143509,   1e-9);
    CHECK_REL(math::cauchyCombineLog10(c8, 2), 8.0,            1e-14);

    // What the linear form does with the same inputs.
    const double p3[] = {1e-300, 0.1, std::pow(10.0, -0.5)};
    const double p4[] = {1e-320, 0.1, std::pow(10.0, -0.5)};
    CHECK_REL(-std::log10(math::cauchyCombine(p3, 3)), 299.522879, 1e-6);
    CHECK(math::cauchyCombine(p4, 3) == 0.0);              // and -log10 is +Inf
    CHECK(std::isfinite(math::cauchyCombineLog10(c4, 3))); // whereas this is not
}

// A single p-value combined with itself alone is that p-value: the Cauchy
// upper tail at cot(pi*p) is p exactly, because atan(cot(x)) = pi/2 - x on
// (0, pi).  That makes the whole tail-inversion path testable against an
// identity rather than against a table, and it sweeps lnT continuously.
TEST(cauchyCombineLog10_single_input_is_the_identity) {
    const double clampL = -std::log10(0.999);
    for (double L = clampL; L < 1e4; L *= 1.01) {
        const double got = math::cauchyCombineLog10(&L, 1);
        // Absolute floor because L -> 0 means p -> 1, where -log10(p) has
        // absolute and not relative accuracy; the floor is one ulp of ln(1).
        CHECK_CLOSE(got, L, 1e-14, 2e-16);
    }
    // Below the clamp everything collapses onto -log10(0.999), which is the
    // L-scale spelling of cauchyCombine's own p >= 1 -> 0.999 rule.
    const double tiny = 1e-9;
    CHECK_REL(math::cauchyCombineLog10(&tiny, 1), clampL, 1e-13);
    double zero = 0.0;
    CHECK_REL(math::cauchyCombineLog10(&zero, 1), clampL, 1e-13);
}

// The atan inversion, swept over the lnT values 01_numerics §2.2 calls out,
// including both sides of the lnT = 350 asymptotic switch and the region just
// above lnT = 1 where the truncated log(atan u)/u series an earlier draft
// proposed is wrong by 3.5e-5 -- seven orders of magnitude outside this
// tolerance.
TEST(cauchyCombineLog10_atan_inversion_against_long_double) {
    const long double lnPi = 1.144729885849400174143427351353058712L;
    for (double lnT : {1.0, 1.5, 2.0, 5.0, 20.0, 100.0, 349.0, 351.0, 700.0}) {
        // Choose the single input whose statistic is (very nearly) exp(lnT).
        const double L = (lnT + static_cast<double>(lnPi)) / kLn10;
        const long double p = std::pow(10.0L, -static_cast<long double>(L));
        const long double pi = 3.14159265358979323846264338327950288L;
        const long double T = std::cos(pi * p) / std::sin(pi * p);
        const long double ref = cauchyUpperNegLog10L(T);
        CHECK_REL(math::cauchyCombineLog10(&L, 1), static_cast<double>(ref), 1e-14);
    }
}

// Against a long-double evaluation of the same statistic over the whole
// three-term family, and against the linear implementation over the narrow
// range where the linear implementation is itself accurate.
//
// The linear form's accuracy is the point worth recording.  `cauchyCombine`
// finishes with `0.5 - atan(T)/pi`, and for small p the statistic T grows like
// 10^L, so that subtraction cancels: at L_max = 6.6 it has already lost
// relative accuracy to 2e-11, at 13.8 to 2e-4, and by 1e15 it is replaced by a
// reciprocal asymptote.  The log-domain form inverts atan(1/T) instead and
// never forms the difference, so the two disagree by exactly the cancellation
// -- with the long-double reference siding with the log-domain form.
TEST(cauchyCombineLog10_matches_long_double_and_uncancelled_linear) {
    for (double Lmax = 0.5; Lmax < 290.0; Lmax *= 1.15) {
        const double Ls[] = {Lmax, 1.0, 0.5};
        CHECK_REL(math::cauchyCombineLog10(Ls, 3),
                  static_cast<double>(cauchyCombineNegLog10RefL(Ls, 3)), 1e-13);
        if (Lmax < 3.0) {
            const double ps[] = {std::pow(10.0, -Lmax), 0.1, std::pow(10.0, -0.5)};
            CHECK_REL(math::cauchyCombineLog10(Ls, 3),
                      -std::log10(math::cauchyCombine(ps, 3)), 1e-11);
        }
    }
    {
        // The cancellation, measured.  Both are computing the same quantity.
        const double Ls[] = {13.8354668408, 1.0, 0.5};
        const double ps[] = {std::pow(10.0, -13.8354668408), 0.1, std::pow(10.0, -0.5)};
        const double ref = static_cast<double>(cauchyCombineNegLog10RefL(Ls, 3));
        CHECK_REL(math::cauchyCombineLog10(Ls, 3), ref, 1e-13);
        CHECK(std::fabs(-std::log10(math::cauchyCombine(ps, 3)) - ref) / ref > 1e-5);
    }
    // 01_numerics §2.4: with equal weights the tail is L_max - log10(n), and
    // it is dominated by the LARGEST L (smallest p).  The opposite spelling,
    // -log10(sum w_i 10^(-L_i)), would give 0.858 for the first row below.
    for (double Lmax : {20.0, 50.0, 100.0, 300.0, 1000.0, 100000.0}) {
        const double Ls[] = {Lmax, 1.0, 0.5};
        CHECK_NEAR(math::cauchyCombineLog10(Ls, 3), Lmax - std::log10(3.0), 5e-7);
    }
    const double two[] = {400.0, 400.0};
    CHECK_REL(math::cauchyCombineLog10(two, 2), 400.0, 1e-14);
}

TEST(cauchyCombineLog10_edge_cases) {
    const double allNaN[] = {kNaN, kNaN};
    CHECK(std::isnan(math::cauchyCombineLog10(allNaN, 2)));
    CHECK(std::isnan(math::cauchyCombineLog10(nullptr, 0)));
    // L = +Inf is p = 0, which short-circuits the combination to p = 0 exactly
    // as cauchyCombine's own `hasZero` branch does.
    const double withInf[] = {kInf, 1.0};
    CHECK(math::cauchyCombineLog10(withInf, 2) == kInf);
    // NaN entries are skipped and n is the count of the survivors, so the
    // padding does not dilute the result.
    const double padded[] = {8.0, kNaN, kNaN};
    CHECK_REL(math::cauchyCombineLog10(padded, 3), 8.0, 1e-14);
    // Never negative (invariant C2) and never -Inf, whatever the mixture.
    const double mixed[] = {0.0, 0.0, 0.0, 700.0};
    CHECK(math::cauchyCombineLog10(mixed, 4) >= 0.0);
    CHECK(std::isfinite(math::cauchyCombineLog10(mixed, 4)));
}

// ══════════════════════════════════════════════════════════════════════
// § 4  ptLog
// ══════════════════════════════════════════════════════════════════════

// The grid 01_numerics §3.3 specifies, against Boost's regularized incomplete
// beta wherever that returns a normal double.  Two of the twelve cells --
// df = 10000 with |t| = 50 and 200 -- are below Boost's underflow and are
// checked against the analytic asymptote instead.
TEST(ptLog_matches_boost_ibeta_on_the_reference_grid) {
    for (double df : {5.0, 100.0, 10000.0}) {
        for (double t : {2.0, 20.0, 50.0, 200.0}) {
            const double a = 0.5 * df, b = 0.5, x = df / (df + t * t);
            const double ib = boost::math::ibeta(a, b, x);
            const double got = math::ptLog(t, df);
            if (ib > 1e-300) {
                CHECK_REL(got, std::log(ib), 1e-13);
            } else {
                CHECK(std::isfinite(got));
                CHECK(got < -700.0);
            }
            CHECK_REL(math::ptLog(-t, df), got, 0.0);   // two-sided: even in t
        }
    }
}

// The seam between the Boost branch and the DLMF 8.17.7 series, scanned on
// both sides at df = 10000 where it falls near |t| = 38.4.  This is the test
// that separates the two 2F1 forms: mixing the Euler numerator with the
// (1-x)^b prefactor of 8.17.7 is wrong by exactly (1-x)^(1/2), which here is
// 0.36 -- a discrepancy of 1.0 in ln p, eleven orders outside the tolerance.
// At small df the same defect is invisible, because x is then near 0.
TEST(ptLog_seam_and_series_branch_agree_with_boost) {
    bool sawBoostBranch = false, sawSeriesBranch = false;
    for (double t = 30.0; t <= 38.9; t += 0.05) {
        const double df = 10000.0;
        const double a = 0.5 * df, b = 0.5, x = df / (df + t * t);
        const double ib = boost::math::ibeta(a, b, x);
        if (!(ib > 1e-305)) continue;
        CHECK_REL(math::ptLog(t, df), std::log(ib), 1e-13);
        // lnLead is the switch's discriminant; -690 is where math_helper puts
        // the branch.  Both sides of it are covered by this scan.
        const double lnLead = a * std::log(x) + b * std::log1p(-x)
                              - std::log(a)
                              - (std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
        if (lnLead >= -690.0) sawBoostBranch = true; else sawSeriesBranch = true;
    }
    CHECK(sawBoostBranch);
    CHECK(sawSeriesBranch);
}

// Against the linear two-sided tail it replaces, wherever that one has not
// underflowed.
TEST(ptLog_matches_two_times_pt_where_linear_survives) {
    for (double df : {1.0, 2.0, 5.0, 30.0, 100.0, 1000.0, 10000.0}) {
        for (double t = 0.25; t < 40.0; t *= 1.35) {
            const double lin = 2.0 * math::pt(-t, df, true);
            if (!(lin > 1e-300)) continue;
            CHECK_REL(math::ptLog(t, df), std::log(lin), 1e-11);
        }
    }
}

TEST(ptLog_deep_tail_and_edge_cases) {
    // Monotone in |t| and finite far past where the linear tail is zero.
    double prev = 1.0;
    for (double t : {50.0, 100.0, 200.0, 1000.0, 10000.0}) {
        const double v = math::ptLog(t, 10000.0);
        CHECK(std::isfinite(v));
        CHECK(v < prev);
        prev = v;
    }
    // -log10 p, the column this feeds, stays an ordinary number.
    CHECK_NEAR(-math::ptLog(50.0, 10000.0) / kLn10, 486.2988243, 1e-6);
    CHECK_NEAR(-math::ptLog(200.0, 10000.0) / kLn10, 3496.899648, 1e-5);

    CHECK(math::ptLog(0.0, 5.0) == 0.0);             // P(|T| > 0) = 1
    CHECK(math::ptLog(kInf, 5.0) == -kInf);
    CHECK(std::isnan(math::ptLog(kNaN, 5.0)));
    CHECK(std::isnan(math::ptLog(2.0, kNaN)));
    CHECK(std::isnan(math::ptLog(2.0, 0.0)));
    CHECK(std::isnan(math::ptLog(2.0, -1.0)));
    // Never positive: it is the logarithm of a probability.
    for (double df : {1.0, 7.0, 500.0})
        for (double t = 0.0; t < 5.0; t += 0.25)
            CHECK(math::ptLog(t, df) <= 0.0);
}

// ══════════════════════════════════════════════════════════════════════
// § 5  pmvnorm2dHalfRectLog
// ══════════════════════════════════════════════════════════════════════

// Both Y bounds infinite reduces to the marginal, exactly as the linear
// routine does.
TEST(pmvnorm2dHalfRectLog_marginal_case) {
    for (double s : {2.0, 0.0, -1.5, -10.0, -50.0}) {
        for (double rho : {-0.7, 0.0, 0.4}) {
            const double var1 = 2.0, var2 = 3.0;
            const double cov = rho * std::sqrt(var1 * var2);
            CHECK_REL(math::pmvnorm2dHalfRectLog(s, -kInf, kInf, var1, cov, var2),
                      math::pnormLog(s / std::sqrt(var1)), 0.0);
        }
    }
}

// Finite rectangle, moderate |rho|, probability well above underflow: here the
// linear routine performs the very same quadrature on the linear scale, and
// the two must agree to round-off.
TEST(pmvnorm2dHalfRectLog_matches_the_linear_routine_on_finite_rectangles) {
    const double var1 = 1.7, var2 = 0.8;
    for (double rho : {-0.9, -0.5, -0.2, 0.0, 0.3, 0.6, 0.9}) {
        const double cov = rho * std::sqrt(var1 * var2);
        for (double s : {1.0, 0.0, -1.0, -3.0, -6.0}) {
            for (double lo : {-2.0, -0.5, 0.5}) {
                const double hi = lo + 1.0;
                const double lin = math::pmvnorm2dHalfRect(s, lo, hi, var1, cov, var2);
                if (!(lin > 1e-280)) continue;
                CHECK_REL(std::exp(math::pmvnorm2dHalfRectLog(s, lo, hi, var1, cov, var2)),
                          lin, 1e-11);
            }
        }
    }
}

// The half-infinite case is a different story and the difference is a finding,
// not a tolerance to be widened.  The linear routine answers it through
// `bvnCdf`, whose |r| < 0.925 arm is Plackett's identity under the SIX-point
// rule Genz specifies for |r| < 0.3; math_helper.cpp's own comment records
// that it agrees with an independent Simpson reduction "to ~1e-7", and that is
// an ABSOLUTE accuracy, so it says nothing at all about a probability below
// 1e-7.  Measured on the grid below at (rho, s) = (-0.9, -6): the linear
// routine returns 1.729e-10 where the true value is 5.215e-24.
//
// Three checks, in increasing strength:
//   * against the linear routine, at the absolute accuracy the linear routine
//     actually has;
//   * against the complement identity P(Y <= k) + P(Y > k) = Phi(h), which
//     needs no external reference and holds to 4e-15 here;
//   * against the long-double Simpson reference, in the test below.
TEST(pmvnorm2dHalfRectLog_half_infinite_is_the_accurate_one) {
    const double var1 = 1.7, var2 = 0.8;
    for (double rho : {-0.9, -0.5, -0.2, 0.0, 0.3, 0.6, 0.9}) {
        const double cov = rho * std::sqrt(var1 * var2);
        for (double s : {1.0, 0.0, -1.0, -3.0, -6.0}) {
            const double lgLo = math::pmvnorm2dHalfRectLog(s, -kInf, 0.4, var1, cov, var2);
            const double lgHi = math::pmvnorm2dHalfRectLog(s, 0.4, kInf, var1, cov, var2);

            // Absolute agreement with bvnCdf, which is all bvnCdf offers.
            CHECK_NEAR(std::exp(lgLo),
                       math::pmvnorm2dHalfRect(s, -kInf, 0.4, var1, cov, var2), 1e-6);
            CHECK_NEAR(std::exp(lgHi),
                       math::pmvnorm2dHalfRect(s, 0.4, kInf, var1, cov, var2), 1e-6);

            // The two half-infinite pieces partition the marginal exactly.
            const double marginal = std::exp(math::pnormLog(s / std::sqrt(var1)));
            CHECK_REL(std::exp(lgLo) + std::exp(lgHi), marginal, 1e-13);

            // Against the long-double reference, relatively, everywhere.
            const long double ref =
                lnCondTailRefL(s / std::sqrt(var1), rho,
                               -std::numeric_limits<long double>::infinity(),
                               0.4 / std::sqrt(var2));
            CHECK_REL(lgLo, static_cast<double>(ref), 1e-12);
        }
    }
    // The specific cell in which the linear route is wrong by fourteen orders
    // of magnitude, pinned so that a later change to bvnCdf is noticed here.
    {
        const double cov = -0.9 * std::sqrt(1.7 * 0.8);
        const double lg = math::pmvnorm2dHalfRectLog(-6.0, -kInf, 0.4, 1.7, cov, 0.8);
        const long double ref =
            lnCondTailRefL(-6.0 / std::sqrt(1.7), -0.9,
                           -std::numeric_limits<long double>::infinity(),
                           0.4 / std::sqrt(0.8));
        CHECK_REL(lg, static_cast<double>(ref), 1e-12);
        CHECK_REL(std::exp(lg), 5.2146929281199e-24, 1e-10);
        CHECK(math::pmvnorm2dHalfRect(-6.0, -kInf, 0.4, 1.7, cov, 0.8) > 1e-11);
    }
}

// Branch (a) against branch (b): a finite lower bound placed far enough below
// the effective support must give the same answer as -inf.  In this
// implementation the two are the SAME integral with the same truncation rule,
// which is why they agree exactly rather than merely to 1e-12; the independent
// evidence is the long-double Simpson check below, not this one.
TEST(pmvnorm2dHalfRectLog_finite_and_infinite_lower_bound_agree) {
    const double var1 = 1.0, var2 = 1.0;
    for (double rho : {-0.95, -0.6, 0.0, 0.35, 0.85}) {
        const double cov = rho;
        for (double s : {2.0, -1.0, -8.0, -30.0, -70.0}) {
            for (double hi : {1.0, -2.0, -25.0}) {
                const double a = math::pmvnorm2dHalfRectLog(s, hi - 400.0, hi, var1, cov, var2);
                const double b = math::pmvnorm2dHalfRectLog(s, -kInf, hi, var1, cov, var2);
                CHECK_REL(a, b, 1e-12);
            }
        }
    }
}

// The independent reference: composite Simpson in long double, on a window
// located by golden section rather than by the routine's own bisection.  The
// cases below straddle the linear routine's underflow, so the last of them
// cannot be checked any other way.
TEST(pmvnorm2dHalfRectLog_against_long_double_simpson) {
    struct Case { double h, rho, k; };
    const Case cases[] = {
        {-1.0, 0.5, -1.0}, {-3.0, 0.3, -2.0}, {-8.0, -0.4, -6.0},
        {2.0, 0.6, 1.0},   {-20.0, 0.7, -20.0}, {-40.0, 0.9, -35.0},
        {-60.0, -0.8, -50.0}, {-100.0, 0.85, -90.0}, {-30.0, 0.0, -30.0},
    };
    for (const Case &c : cases) {
        const double got = math::pmvnorm2dHalfRectLog(c.h, -kInf, c.k, 1.0, c.rho, 1.0);
        const long double ref = lnCondTailRefL(c.h, c.rho,
                                               -std::numeric_limits<long double>::infinity(),
                                               c.k);
        CHECK_REL(got, static_cast<double>(ref), 1e-12);
    }
    // The same, on a finite rectangle.
    struct FCase { double h, rho, a, b; };
    const FCase fcases[] = {
        {-1.0, 0.5, -1.5, -0.5}, {-2.0, -0.3, 0.0, 1.0}, {1.0, 0.7, -0.4, 0.6},
        {-5.0, 0.2, -2.0, -1.0}, {-1.0, 0.5, -4.0, 3.0}, {-10.0, 0.6, -12.0, -8.0},
        {-40.0, -0.5, -30.0, -20.0},
    };
    for (const FCase &f : fcases) {
        const double got = math::pmvnorm2dHalfRectLog(f.h, f.a, f.b, 1.0, f.rho, 1.0);
        const long double ref = lnCondTailRefL(f.h, f.rho, f.a, f.b);
        CHECK_REL(got, static_cast<double>(ref), 1e-12);
    }
}

// Two analytic identities, chosen because they hold for EVERY rho and so walk
// the routine into the degenerate regime without needing a reference that
// itself degrades there.
//
//   * the orthant probability  P(X <= 0, Y <= 0) = 1/4 + asin(rho)/(2*pi);
//   * rho = 0 factorizes,  P = Phi(h)*Phi(k).
//
// The orthant sweep is the test that the Phi-argument step cap earns: at
// rho = 1 - 1e-9 the conditional CDF is a cliff of width 1/(c|rho|) = 4.5e-5
// sitting exactly on the peak, and a panel chosen from the gradient and
// curvature at its near edge alone steps straight over it.
TEST(pmvnorm2dHalfRectLog_analytic_identities_including_degenerate_rho) {
    const double pi = 3.14159265358979323846;
    for (double rho : {-0.999999999, -0.999999, -0.99, -0.9, -0.5, 0.0,
                       0.5, 0.9, 0.99, 0.999999, 0.999999999}) {
        const double exact = 0.25 + std::asin(rho) / (2.0 * pi);
        CHECK_REL(std::exp(math::pmvnorm2dHalfRectLog(0.0, -kInf, 0.0, 1.0, rho, 1.0)),
                  exact, 1e-11);
        // Both Y bounds finite, over a rectangle that is the same orthant
        // minus a strip whose probability is known the same way.
        const double lower = 0.25 + std::asin(rho) / (2.0 * pi);   // Y <= 0
        const double both  = std::exp(math::pmvnorm2dHalfRectLog(0.0, -20.0, 0.0, 1.0, rho, 1.0));
        CHECK_REL(both, lower, 1e-11);
    }
    // rho = 0 factorizes exactly, including where both factors have underflowed.
    for (double h : {1.0, -2.0, -25.0, -60.0})
        for (double k : {0.5, -3.0, -40.0})
            CHECK_REL(math::pmvnorm2dHalfRectLog(h, -kInf, k, 1.0, 0.0, 1.0),
                      math::pnormLog(h) + math::pnormLog(k), 1e-12);
    // And the rho -> 1 limit itself, P -> Phi(min(h,k)).  The residual is the
    // genuine O(sqrt(1-rho^2)) correction, not a quadrature error, so the
    // tolerance is written in those terms rather than as a constant.
    for (double h : {2.0, 0.0, -1.0, -3.0, -8.0, -20.0}) {
        for (double k : {3.0, 0.0, -2.0}) {
            for (double rho : {1.0 - 1e-6, 1.0 - 1e-9, 1.0 - 1e-13}) {
                const double got = math::pmvnorm2dHalfRectLog(h, -kInf, k, 1.0, rho, 1.0);
                const double slack = 1.0 * std::sqrt((1.0 - rho) * (1.0 + rho)) + 1e-9;
                CHECK_NEAR(got, math::pnormLog(std::fmin(h, k)), slack);
            }
        }
    }
}

// Invariant C1: this routine may return -Inf, never +Inf, and it reports an
// indefinite covariance as NaN rather than clamping it into a different
// problem -- the same contract the linear routine carries.
TEST(pmvnorm2dHalfRectLog_failure_modes) {
    CHECK(math::pmvnorm2dHalfRectLog(0.0, -1.0, 1.0, -1.0, 0.0, 1.0) == -kInf);
    CHECK(math::pmvnorm2dHalfRectLog(0.0, -1.0, 1.0, 1.0, 0.0, 0.0) == -kInf);
    CHECK(math::pmvnorm2dHalfRectLog(0.0, 1.0, -1.0, 1.0, 0.0, 1.0) == -kInf);  // empty
    CHECK(std::isnan(math::pmvnorm2dHalfRectLog(0.0, -1.0, 1.0, 1.0, 1.5, 1.0)));
    CHECK(std::isnan(math::pmvnorm2dHalfRect(0.0, -1.0, 1.0, 1.0, 1.5, 1.0)));

    // Nothing anywhere in a broad sweep comes back +Inf or positive.
    for (double h = -80.0; h <= 5.0; h += 2.5) {
        for (double rho : {-0.999999, -0.9, -0.3, 0.0, 0.3, 0.9, 0.999999}) {
            for (double k : {3.0, -1.0, -18.0, -60.0}) {
                const double v = math::pmvnorm2dHalfRectLog(h, -kInf, k, 1.0, rho, 1.0);
                CHECK(!(v > 0.0));
                CHECK(v != kInf);
                const double w = math::pmvnorm2dHalfRectLog(h, k - 3.0, k, 1.0, rho, 1.0);
                CHECK(!(w > 0.0));
                CHECK(w != kInf);
            }
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
// Stage 9 appends its HweLnP section here.
// ══════════════════════════════════════════════════════════════════════

TINYTEST_MAIN
