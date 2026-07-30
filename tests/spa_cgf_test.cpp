// spa_cgf_test.cpp — Correctness of the binomial CGF evaluation.
//
// ══════════════════════════════════════════════════════════════════════
// Two layers, added by two stages
// ══════════════════════════════════════════════════════════════════════
//
// PART I (Stage 0) characterizes the defect using the reference copies in
// tests/spa_reference.hpp.  Those tests are retained verbatim: `cgfDiffForm` is
// a frozen copy of the OLD production algebra, kept so the before/after
// comparison stays measurable after the production code is replaced, and the
// characterization assertions fail if anyone reintroduces the differencing form.
//
// PART II (Stage 2) tests the REAL kernels in src/util/spa_cgf.{hpp,cpp}:
//
//   * all three variants (binomUniform, binomIndiv, binomHapcount), both entry
//     points (k12, kFull), against a long-double reference;
//   * the reductions, not only the per-subject algebra, so the accumulation is
//     covered as well as the closed form;
//   * cross-tier equality — scalar vs AVX2 vs AVX-512 vs the dispatched entry
//     point — on every tier the host supports, guarded by simdLevelValue() so an
//     unsupported tier is never invoked (SIGILL);
//   * K'' strictly positive and finite over the grid where the production form
//     produces 20 non-positive or non-finite values out of 288;
//   * the degenerate inputs the production `kG012Local` mishandles: hapcount
//     h = 0, and allele frequency at exactly 0 or exactly 1.
//
// Every random draw is seeded from a fixed constant, so the suite is
// reproducible run to run.
//
// ── Original Stage 0 header follows ──────────────────────────────────
//
// This file establishes, before any production code is touched:
//
//   (a) that the closed-form algebra in tests/spa_reference.hpp really is the
//       CGF of r·G with G ~ Binomial(2, p), by checking it against explicit
//       enumeration over the three genotype states;
//
//   (b) that the two closed forms agree to within rounding wherever the
//       saddlepoint argument is benign, so the Stage 2 substitution is an
//       algebraic identity and not a change of model;
//
//   (c) that the production form loses precision catastrophically as t·r
//       grows, while the stable form does not — the defect Stage 2 fixes.
//
// Point (c) is written as an assertion rather than a printout so that it
// becomes a regression test in the opposite direction once Stage 2 lands: if
// someone reintroduces the differencing algebra, `stable_form_has_no_
// cancellation` fails.
//
// Build:  make test        (or: make build/tests/spa_cgf_test && ./build/tests/spa_cgf_test)

#include "tinytest.hpp"
#include "spa_reference.hpp"

#include "util/spa_cgf.hpp"
// For the terminal-K argument: w and the Barndorff-Nielsen tail are what K is
// consumed through, so the test that justifies the cheaper spelling has to
// evaluate them rather than reason about them.
#include "util/spa.hpp"

#include <cfloat>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using sparef::Cgf;
using sparef::CgfL;

namespace {

// Representative minor-allele frequencies, spanning the range GRAB admits:
// the rare end is where alpha is dominated by (1-p), the common end where the
// two terms are comparable.
const std::vector<double> kMafs = {1e-6, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.25, 0.5};

// Residual weights.  GRAB centres residuals, so both signs occur; the
// magnitude spread reflects what the IQR outlier filter admits.
const std::vector<double> kResids = {-8.0, -3.5, -1.0, -0.25, 0.25, 1.0, 3.5, 8.0};

// Saddlepoint arguments in the benign region, where |t*r| stays small enough
// that no term of either algebra is near cancellation or overflow.
const std::vector<double> kTBenign = {-0.5, -0.2, -0.05, -1e-3, 0.0, 1e-3, 0.05, 0.2, 0.5};

} // namespace

// ──────────────────────────────────────────────────────────────────────
// (a) The closed form is the binomial CGF.
// ──────────────────────────────────────────────────────────────────────

TEST(closed_form_matches_genotype_enumeration) {
    for (double p : kMafs) {
        for (double r : kResids) {
            for (double t : kTBenign) {
                const CgfL ref = sparef::cgfReference(t, r, p);
                const CgfL enu = sparef::cgfEnumerate(t, r, p);

                // K vanishes at t = 0 and is O(1e-8) near it, so the criterion
                // must be mixed absolute/relative; a pure relative bound is
                // unsatisfiable at the zero.  K′ and K″ have no zero in this
                // domain (r != 0, p > 0) and take a relative bound.
                //
                // The enumeration path forms K″ as E[X²] − E[X]², so it carries
                // the very cancellation the stable form avoids; in long double
                // at these benign arguments that costs a few digits, which the
                // K″ tolerance accommodates.
                CHECK_CLOSE(static_cast<double>(ref.K0), static_cast<double>(enu.K0),
                            1e-14, 1e-18);
                CHECK_REL(static_cast<double>(ref.K1), static_cast<double>(enu.K1), 1e-16);
                CHECK_REL(static_cast<double>(ref.K2), static_cast<double>(enu.K2), 1e-12);
            }
        }
    }
}

// K(0) = 0, K'(0) = E[rG] = 2pr, K''(0) = Var(rG) = 2p(1-p)r².
// These are the values the saddlepoint machinery relies on when it folds the
// non-outlier subjects into a Gaussian block, so an error here would
// mis-centre every method.
TEST(cumulants_at_zero_match_moments) {
    for (double p : kMafs) {
        for (double r : kResids) {
            const Cgf s = sparef::cgfStableForm(0.0, r, p);
            CHECK_NEAR(s.K0, 0.0, 0.0);                       // exactly zero
            CHECK_REL(s.K1, 2.0 * p * r, 1e-15);
            CHECK_REL(s.K2, 2.0 * p * (1.0 - p) * r * r, 1e-15);
        }
    }
}

// ──────────────────────────────────────────────────────────────────────
// (b) The two closed forms are the same algebra.
// ──────────────────────────────────────────────────────────────────────

// Rather than assert that the two forms agree with each other — which would
// only bound their difference, not their error — measure both against the
// long-double reference.  This establishes the substitution is an algebraic
// identity (both converge on the same value) and simultaneously establishes
// which of the two is the more accurate double-precision evaluation.
TEST(production_and_stable_forms_agree_where_benign) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    double worstProdK2 = 0.0, worstStableK2 = 0.0;
    int stableWorseCount = 0, nPoints = 0;

    for (double p : kMafs) {
        for (double r : kResids) {
            for (double t : kTBenign) {
                const CgfL ref = sparef::cgfReference(t, r, p);
                const Cgf  d   = sparef::cgfDiffForm(t, r, p);
                const Cgf  s   = sparef::cgfStableForm(t, r, p);

                const double refK0 = static_cast<double>(ref.K0);
                const double refK1 = static_cast<double>(ref.K1);
                const double refK2 = static_cast<double>(ref.K2);

                // K is deliberately not asserted here.  Neither of these two
                // forms evaluates it accurately near t = 0 — both materialize
                // alpha = 1 + delta before taking the logarithm — and that is a
                // separate defect with its own test, `log1p_form_fixes_K_near_
                // zero`.  Asserting a loose bound on K here would only obscure
                // which form is responsible for which error.
                (void)refK0;
                CHECK_REL(d.K1, refK1, 1e-14);
                CHECK_REL(s.K1, refK1, 1e-15);

                const double ed = tinytest::relDiff(d.K2, refK2);
                const double es = tinytest::relDiff(s.K2, refK2);
                worstProdK2   = std::fmax(worstProdK2, ed);
                worstStableK2 = std::fmax(worstStableK2, es);
                ++nPoints;

                // The stable form must never be materially worse.  A small
                // slack absorbs points where both are already at rounding
                // level and the ordering is arbitrary.
                if (es > ed + 1e-15) ++stableWorseCount;

                // The stable form must be at rounding level everywhere in the
                // benign region; the production form is not held to this.
                CHECK_MSG(es <= 1e-14, "stable K'' exceeded rounding level");
            }
        }
    }

    std::fprintf(stderr,
                 "      over %d benign points: worst K'' relative error — "
                 "production %.3e, stable %.3e\n",
                 nPoints, worstProdK2, worstStableK2);
    CHECK_MSG(stableWorseCount == 0,
              "stable form was materially less accurate at one or more points");
}

// K itself, as opposed to K″, loses relative precision near t = 0 in every
// form that materializes alpha = 1 + delta before taking the logarithm.  K
// enters w = sgn(zeta)·sqrt(2·(zeta·s − K(zeta))) directly, so this error
// propagates into the dominant term of the modified signed root.  The
// expm1/log1p formulation removes it at identical cost.
TEST(log1p_form_fixes_K_near_zero) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    struct Row { double tr; double relPlain; double relL1p; };
    std::vector<Row> table;

    const double p = 0.3;
    const double r = 1.0;

    for (double tr : {1e-8, 1e-6, 1e-4, 1e-2, 1.0}) {
        const CgfL ref = sparef::cgfReference(
            static_cast<long double>(tr), 1.0L, static_cast<long double>(p));
        const double refK0 = static_cast<double>(ref.K0);

        const double plain = sparef::cgfStableForm(tr, r, p).K0;
        const double l1p   = sparef::cgfStableFormL1P(tr, r, p).K0;

        const double ep = tinytest::relDiff(plain, refK0);
        const double el = tinytest::relDiff(l1p, refK0);
        table.push_back(Row{tr, ep, el});

        // The log1p form must stay at rounding level all the way to t = 0.
        CHECK_MSG(el <= 1e-15, "log1p form failed to keep K at rounding level");
        CHECK(el <= ep + 1e-16);
    }

    std::fprintf(stderr, "\n      K relative error vs long-double reference (MAF=0.3, r=1)\n");
    std::fprintf(stderr, "      %10s  %14s  %14s\n", "t*r", "2*log(alpha)", "2*log1p(delta)");
    for (const Row &row : table)
        std::fprintf(stderr, "      %10.0e  %14.3e  %14.3e\n",
                     row.tr, row.relPlain, row.relL1p);
    std::fprintf(stderr, "\n");

    // Characterization: by t*r = 1e-6 the plain form has lost at least six
    // decimal digits of K, and by 1e-8 nearly nine.  This assertion pins the
    // motivation for adopting expm1/log1p in Stage 2; it fails if someone
    // reverts to the plain logarithm.
    CHECK_MSG(table[1].relPlain > 1e-10,
              "expected 2*log(alpha) to have lost >=6 digits of K at t*r=1e-6");
    CHECK_MSG(table[0].relPlain > 1e-10,
              "expected 2*log(alpha) to have lost >=6 digits of K at t*r=1e-8");
}

// Adopting the log1p formulation must be a precision change, not a model
// change.  Two claims, over the full parameter grid:
//
//   1. K' and K'' are bit-identical between the two forms — the log1p variant
//      touches only the K branch, so nothing the Newton loop consumes moves.
//   2. K from the log1p form is at rounding level against the reference
//      everywhere, and never worse than the plain form.
TEST(log1p_form_is_a_pure_precision_improvement) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    double worstPlain = 0.0, worstL1p = 0.0;

    for (double p : kMafs) {
        for (double r : kResids) {
            for (double t : {-2.0, -0.5, -1e-3, 0.0, 1e-3, 0.5, 2.0}) {
                const Cgf a = sparef::cgfStableForm(t, r, p);
                const Cgf b = sparef::cgfStableFormL1P(t, r, p);

                // (1) bit-identical, not merely close.
                CHECK_MSG(a.K1 == b.K1, "K' moved when only K should have");
                CHECK_MSG(a.K2 == b.K2, "K'' moved when only K should have");

                // (2) accuracy of K against the reference.
                const double refK0 = static_cast<double>(
                    sparef::cgfReference(t, r, p).K0);
                const double ep = tinytest::relDiff(a.K0, refK0);
                const double el = tinytest::relDiff(b.K0, refK0);
                worstPlain = std::fmax(worstPlain, ep);
                worstL1p   = std::fmax(worstL1p, el);

                CHECK_MSG(el <= 1e-15, "log1p K exceeded rounding level");
                // "Never materially worse."  Where both forms are already at
                // rounding level the ordering between them is arbitrary, so
                // the comparison is only meaningful once the plain form's
                // error exceeds a few units in the last place.
                CHECK_MSG(el <= ep || el <= 4e-16,
                          "log1p K was materially worse than the plain form");
            }
        }
    }

    std::fprintf(stderr,
                 "      worst K relative error over the grid — "
                 "2*log(alpha) %.3e, 2*log1p(delta) %.3e\n",
                 worstPlain, worstL1p);
}

// ──────────────────────────────────────────────────────────────────────
// (c) The production form cancels; the stable form does not.
// ──────────────────────────────────────────────────────────────────────

// Accuracy of each form against the long-double reference, as t*r grows.
// The saddlepoint solver drives t*r to large positive values for markers in
// the extreme tail, which is precisely where the reported p-value matters.
TEST(stable_form_has_no_cancellation) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr,
                     "      SKIP: long double is not wider than double on this "
                     "platform; no reference available.\n");
        return;
    }

    struct Row { double tr; double relDiffErr; double relStableErr; };
    std::vector<Row> table;

    const double p = 0.5;   // worst case: alpha's two terms are comparable
    const double r = 1.0;

    for (double tr : {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0}) {
        const CgfL ref = sparef::cgfReference(
            static_cast<long double>(tr), 1.0L, static_cast<long double>(p));
        const Cgf d = sparef::cgfDiffForm(tr, r, p);
        const Cgf s = sparef::cgfStableForm(tr, r, p);

        const double refK2 = static_cast<double>(ref.K2);
        const double ed = tinytest::relDiff(d.K2, refK2);
        const double es = tinytest::relDiff(s.K2, refK2);
        table.push_back(Row{tr, ed, es});

        // The stable form must stay at rounding level across the whole range.
        CHECK_MSG(es <= 1e-14,
                  "stable form lost precision where it must not");

        // The stable form must be at least as accurate as the differencing
        // form everywhere, and strictly better once cancellation bites.
        CHECK(es <= ed || ed <= 1e-15);
    }

    std::fprintf(stderr, "\n      K'' relative error vs long-double reference (MAF=0.5, r=1)\n");
    std::fprintf(stderr, "      %8s  %14s  %14s\n", "t*r", "production", "stable");
    for (const Row &row : table)
        std::fprintf(stderr, "      %8.1f  %14.3e  %14.3e\n",
                     row.tr, row.relDiffErr, row.relStableErr);
    std::fprintf(stderr, "\n");

    // Characterization of the defect: by t*r = 30 the production form has lost
    // at least six decimal digits.  Once Stage 2 replaces the production
    // algebra this assertion still holds, because cgfDiffForm in this test
    // file is a frozen copy of the OLD algebra kept for exactly this purpose.
    const double errAt30 = table[table.size() - 2].relDiffErr;
    CHECK_MSG(errAt30 > 1e-7,
              "expected the differencing form to have lost >=6 digits at t*r=30");
}

// The saddlepoint solver divides by K'' every Newton iteration.  A K'' that
// underflows to zero or turns negative through cancellation produces an
// infinite or NaN step.  The stable form is positive by construction (a
// product of positive quantities), which this test asserts over a wide
// domain including the region where the production form fails.
TEST(stable_form_curvature_is_strictly_positive) {
    int nProductionNonPositive = 0;

    for (double p : kMafs) {
        for (double r : {-8.0, -1.0, 1.0, 8.0}) {
            for (double tr : {-60.0, -40.0, -20.0, -5.0, 0.0, 5.0, 20.0, 40.0, 60.0}) {
                const double t = tr / r;
                const Cgf s = sparef::cgfStableForm(t, r, p);
                CHECK_MSG(s.K2 > 0.0 && std::isfinite(s.K2),
                          "stable K'' must be finite and strictly positive");

                const Cgf d = sparef::cgfDiffForm(t, r, p);
                if (!(d.K2 > 0.0) || !std::isfinite(d.K2)) ++nProductionNonPositive;
            }
        }
    }

    std::fprintf(stderr,
                 "      production form yielded non-positive or non-finite K'' "
                 "in %d of %zu sampled points\n",
                 nProductionNonPositive, kMafs.size() * 4 * 9);
}

// ══════════════════════════════════════════════════════════════════════
// PART II — the real kernels in src/util/spa_cgf.{hpp,cpp}
// ══════════════════════════════════════════════════════════════════════

namespace {

// ── A long-double reference with a variable trial count ──────────────
//
// sparef::cgfReference fixes h = 2.  binomHapcount needs h to vary, so this
// generalizes it.  The generalization is anchored to the audited reference by
// `hapcount_reference_generalizes_the_diploid_reference` below, which pins the
// two together at h = 2 rather than letting this stand on its own authority.
struct RefL { long double K0, K1, K2; };

RefL refHap(long double t, long double r, long double p, long double h) {
    const long double u    = 1.0L - p;
    const long double tr   = t * r;
    const long double e    = p * std::exp(tr);
    const long double a    = u + e;
    const long double inva = 1.0L / a;
    return RefL{h * std::log1p(p * std::expm1(tr)),
                h * r * e * inva,
                h * r * r * e * u * inva * inva};
}

// ── Grids ────────────────────────────────────────────────────────────
//
// Residual weights restricted to powers of two for the accuracy tests.  The
// grids below drive the kernel with t = tr / r; when r is a power of two that
// division is exact and t * r reproduces tr bit-for-bit, in double and in long
// double alike.  Without that the reference and the kernel would evaluate at
// arguments differing by ~eps*|tr|, which at tr = 40 propagates to a ~4e-15
// relative difference in lambda and the measurement would be reporting argument
// rounding rather than the algebra under test.
const std::vector<double> kResidsP2 = {-8.0, -1.0, -0.25, 0.25, 1.0, 8.0};

// Saddlepoint arguments spanning the near-zero region (where N1 bites) and the
// extreme tail (where D1 bites).
const std::vector<double> kTrWide =
    {-40.0, -20.0, -5.0, -1.0, -1e-3, 0.0, 1e-3, 1.0, 5.0, 20.0, 40.0};

// Trial counts SPAmixLocalPlus supplies: a subject carries 0, 1 or 2 haplotypes
// of the ancestry under test.
const std::vector<double> kHaps = {0.0, 1.0, 2.0};

// ── Tier enumeration ─────────────────────────────────────────────────
//
// A tier is added to the list only when the host's ISA supports it, so an
// unsupported variant is never called and cannot raise SIGILL.  Pattern from
// tests/lanc_simd_test.cpp.  Level: 0 = scalar, 1 = AVX2, 2 = AVX-512.
const int kLevel = spa_cgf::simdLevelValue();

using UniformFn   = spa_cgf::Cgf12 (*)(double, const double *, int, double) noexcept;
using IndivFn     = spa_cgf::Cgf12 (*)(double, const double *, const double *, int) noexcept;
using HapcountFn  = spa_cgf::Cgf12 (*)(
    double, const double *, const double *, int, double) noexcept;

template <class Fn>
struct Tier { const char *name; Fn fn; };

std::vector<Tier<UniformFn>> uniformTiers() {
    std::vector<Tier<UniformFn>> v;
    v.push_back({"scalar", spa_cgf::tier::binomUniformK12_scalar});
#if defined(__x86_64__) || defined(_M_X64)
    if (kLevel >= 1) v.push_back({"avx2", spa_cgf::tier::binomUniformK12_avx2});
    if (kLevel >= 2) v.push_back({"avx512", spa_cgf::tier::binomUniformK12_avx512});
#endif
    v.push_back({"dispatch", spa_cgf::binomUniformK12});
    return v;
}

std::vector<Tier<IndivFn>> indivTiers() {
    std::vector<Tier<IndivFn>> v;
    v.push_back({"scalar", spa_cgf::tier::binomIndivK12_scalar});
#if defined(__x86_64__) || defined(_M_X64)
    if (kLevel >= 1) v.push_back({"avx2", spa_cgf::tier::binomIndivK12_avx2});
    if (kLevel >= 2) v.push_back({"avx512", spa_cgf::tier::binomIndivK12_avx512});
#endif
    v.push_back({"dispatch", spa_cgf::binomIndivK12});
    return v;
}

std::vector<Tier<HapcountFn>> hapcountTiers() {
    std::vector<Tier<HapcountFn>> v;
    v.push_back({"scalar", spa_cgf::tier::binomHapcountK12_scalar});
#if defined(__x86_64__) || defined(_M_X64)
    if (kLevel >= 1) v.push_back({"avx2", spa_cgf::tier::binomHapcountK12_avx2});
    if (kLevel >= 2) v.push_back({"avx512", spa_cgf::tier::binomHapcountK12_avx512});
#endif
    v.push_back({"dispatch", spa_cgf::binomHapcountK12});
    return v;
}

// The same enumeration for the terminal-K reductions.  These return the bare
// sum: sum_i log(alpha_i) for the two diploid variants, with the common factor
// h == 2 left to the caller, and sum_i hap_i*log(alpha_i) for binomHapcount.
using UniformK0Fn  = double (*)(double, const double *, int, double) noexcept;
using IndivK0Fn    = double (*)(double, const double *, const double *, int) noexcept;
using HapcountK0Fn = double (*)(
    double, const double *, const double *, int, double) noexcept;

std::vector<Tier<UniformK0Fn>> uniformK0Tiers() {
    std::vector<Tier<UniformK0Fn>> v;
    v.push_back({"scalar", spa_cgf::tier::binomUniformK0_scalar});
#if defined(__x86_64__) || defined(_M_X64)
    if (kLevel >= 1) v.push_back({"avx2", spa_cgf::tier::binomUniformK0_avx2});
    if (kLevel >= 2) v.push_back({"avx512", spa_cgf::tier::binomUniformK0_avx512});
#endif
    return v;
}

std::vector<Tier<IndivK0Fn>> indivK0Tiers() {
    std::vector<Tier<IndivK0Fn>> v;
    v.push_back({"scalar", spa_cgf::tier::binomIndivK0_scalar});
#if defined(__x86_64__) || defined(_M_X64)
    if (kLevel >= 1) v.push_back({"avx2", spa_cgf::tier::binomIndivK0_avx2});
    if (kLevel >= 2) v.push_back({"avx512", spa_cgf::tier::binomIndivK0_avx512});
#endif
    return v;
}

std::vector<Tier<HapcountK0Fn>> hapcountK0Tiers() {
    std::vector<Tier<HapcountK0Fn>> v;
    v.push_back({"scalar", spa_cgf::tier::binomHapcountK0_scalar});
#if defined(__x86_64__) || defined(_M_X64)
    if (kLevel >= 1) v.push_back({"avx2", spa_cgf::tier::binomHapcountK0_avx2});
    if (kLevel >= 2) v.push_back({"avx512", spa_cgf::tier::binomHapcountK0_avx512});
#endif
    return v;
}

// ══════════════════════════════════════════════════════════════════════
// The ULP budgets, and where each number comes from
// ══════════════════════════════════════════════════════════════════════
//
// Two budgets, because the two things being bounded are different: the ALGEBRA
// in src/util/spa_cgf.hpp, and the VENDORED EXP in src/util/simd_math.hpp that
// the vector tiers call.  Collapsing them into one number would hide which of
// the two any future regression came from.
//
// ── kAlgebraUlp = 4 — the scalar tier against the long-double reference ──
//
// The scalar tier calls glibc's `exp`, measured at 0.8 ULP over this grid and
// effectively correctly rounded.  What this budget bounds is therefore this
// file's own arithmetic: one rounding each for e, alpha, inva, pi, omp and the
// closing products.
//
// Measured worst over the full grid, all three variants: K' 1.0 ULP,
// K'' 2.0 ULP, K 1.0 ULP.  4 ULP is a 2x headroom.
//
// This is the assertion that actually validates the cancellation-free identity.
// The differencing form it replaces reaches 3.6e-1 — fourteen orders of
// magnitude away — so this bound is not remotely at risk of passing a wrong
// K'' identity.
//
// ── kVectorUlp = kCrossTierUlp = 64 — anything touching the vector exp ──
//
// The additional error at the vector tiers is entirely the vendored exp, not
// the algebra.  `avx2_exp_pd` and `avx512_exp_pd` are Cephes-style degree-11
// minimax kernels whose header claims "~1 ULP accuracy over the full double
// range"; measured against a long-double reference they reach 38.9 ULP over
// x in [-350, 350], worst at x = 349, and 17.7 ULP at x = 35.  Both tiers share
// the same polynomial and reduction, so both carry the same error.
//
// The attribution is direct: substituting the vector exp for std::exp and
// changing nothing else moves the per-subject K'' from 2.0 ULP to 10.0 ULP.
// simd_math.hpp is shared, validated infrastructure that Stage 2 is explicitly
// not permitted to modify, so the budget accommodates it rather than fixing it.
//
// A relative perturbation d on lambda reaches pi = e/alpha with gain
// (1 - pi) <= 1 and pi*(1-pi) with gain |1 - 2*pi| <= 1, so a summand inherits
// the exp's error roughly one-for-one and can reach ~39 ULP on its own.
//
// Reduction order contributes on top of that: the scalar tier accumulates
// sequentially, AVX2 carries 4 partial sums closed by a pairwise horizontal
// add, and AVX-512 carries 8 closed by _mm512_reduce_add_pd's log-depth tree.
// Two orderings of n terms can differ by up to ~n*eps and these tests run n up
// to 512.  Bit-identity across tiers is therefore not merely unachieved, it is
// not the correct expectation.
//
// Measured worst: cross-tier 10.0 ULP, reduction against the long-double sum
// 10.4 ULP, AVX2 vector body 13.8 ULP.  64 ULP is a 4.6x headroom over the
// worst observed and still eight or more orders tighter than an algebra error.
//
// ── Why the L1 scale rather than the value of the sum ──
//
// K1's summands h*r_i*pi_i carry the sign of r_i, and GRAB centres its
// residuals: at t = 0 every pi_i is equal, so the sum collapses to
// 2*p*sum(r_i) ~ 0 while the summands themselves are O(1).  A purely relative
// criterion on such a cancelling sum is unsatisfiable by construction — two
// tiers can agree to the last bit on every summand and still differ by 100% in
// the cancelled total.  Comparing against sum_i |term_i| is the correct
// criterion for a mixed-sign reduction.  K2's summands h*r_i^2*pi_i*(1-pi_i)
// are all non-negative, so there the L1 scale equals the sum and the criterion
// degenerates to the ordinary relative one at no loss.
constexpr double kAlgebraUlp   = 4.0;
constexpr double kVectorUlp    = 64.0;
constexpr double kCrossTierUlp = 64.0;

// ── The budget for the PRODUCTION terminal K, and why it is not relative ──
//
// The production kFull forms K as sum_i h_i*log(u + p_i*e^{t*r_i}) with a
// vectorized logarithm, rather than sum_i h_i*log1p(p_i*expm1(t*r_i)).  The two
// spellings differ in RELATIVE accuracy near t = 0, where alpha = 1 + delta and
// rounding alpha discards delta's low bits, but K is consumed only through
// w = sgn(zeta)*sqrt(2*(zeta*s - K)) at abscissae where zeta*s - K is bounded
// away from zero by the |z| > spaCutoff entry gate, so what propagates is the
// ABSOLUTE error.  See the terminal-K section of src/util/spa_cgf.hpp for the
// measurement, and `production_K_absolute_error_does_not_degrade_w` below for the
// assertion that carries it.
//
// The right scale for that absolute error is therefore NOT |K| but
//
//     sum_i h_i * max(1, |log alpha_i|)
//
// one unit of binomial weight per subject, floored at one.  The floor is the
// point: alpha is within rounding distance of 1 near t = 0, so forming it costs
// an absolute eps whatever log(alpha) itself happens to be, and dividing that by
// a vanishing |K| would be measuring the conditioning of the problem rather than
// the quality of the kernel.  Above |log alpha| = 1 the scale becomes the
// ordinary L1 scale and the criterion degenerates to the relative one.
//
// ── Where the number 64 comes from, and why it is not 512 ──
//
// This budget was 512, on the reasoning that the vendored `avx512_exp_pd` was
// measured at up to 38.9 ULP, that its error enters alpha multiplied by
// p*lambda/alpha which the grid drives to 1, and that 512 "leaves room for the
// reduction on top of that".  Two things are wrong with that as a justification.
//
// First, the reduction is not on top of it.  Every assertion that uses this
// budget is an n = 1 assertion — the II.b per-subject grids, the near-zero grid
// and the symmetric grid.  The multi-subject reductions are asserted against
// kCrossTierUlp, not against this constant.  There is no reduction for the
// headroom to accommodate.
//
// Second, the 38.9 ULP figure does not survive the scale.  The exp's worst case
// is at x = 349, and there alpha ~ p*e^x so |log alpha| ~ 349, which is exactly
// the regime in which k0Scale stops being the floor of 1 and becomes the L1
// scale.  A relative error of 39 ULP in e^x is 39*eps absolute in log(alpha),
// which is 39/|log alpha| ULP of the scale — under one ULP where the exp is at
// its worst.  Near t = 0, where the scale is at its floor and 39 ULP would cost
// 39 ULP of it, the exp is at its most accurate.  The two ends do not coincide,
// which is why the measured worst over EVERY site in this file is 3.50 ULP
// (reported by `terminal_K_budget_headroom` at the end of the file).
//
// 512 was therefore a 146x headroom in a file whose stated convention is 2x to
// 4.6x, justified by a mechanism no assertion exercises.  The budget is now 64,
// which is `kVectorUlp` — deliberately the same number, because this quantity
// touches exactly the same vendored vector exp that kVectorUlp was set for, and
// two budgets bounding one error source should not be two numbers.  Against the
// measured 3.50 that is 18x, still above the convention; against the a-priori
// per-subject bound of 38.9 ULP that the vendored exp could in principle deliver
// if a future grid point put the worst case and the scale floor together, it is
// 1.6x, and that is the constraint the budget genuinely has to cover.
constexpr double kK0ProdUlp = 64.0;

// A ULP budget expressed as a relative tolerance, for CHECK_REL / CHECK_CLOSE.
constexpr double relTol(double ulp) { return ulp * DBL_EPSILON; }

// True when this tier's arithmetic goes through the vendored vector exp.  At
// n = 1 the AVX2 kernel takes its scalar tail and so agrees with the scalar tier
// bit-for-bit; only AVX-512, whose tail is masked rather than scalar, exercises
// the vector exp at n = 1.
bool usesVectorExp(const char *tierName, int n) {
    const std::string s(tierName);
    if (s == "scalar") return false;
    if (s == "avx2") return n >= 4;
    return true;   // avx512 and dispatch
}

// L1 scale of the per-subject summands, in double, via the shared scalar core.
// K0's summands h_i*log(alpha_i) are mixed-sign too (negative for t < 0), so it
// needs the same treatment as K1.
struct L1 { double K0, K1, K2; };

L1 l1Scale(
    double t,
    const std::vector<double> &r,
    const std::vector<double> &p,
    const std::vector<double> &h
) {
    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (size_t i = 0; i < r.size(); ++i) {
        const spa_cgf::Cgf012 c = spa_cgf::subjectKFull(t, r[i], p[i], h[i]);
        s0 += std::fabs(c.K0);
        s1 += std::fabs(c.K1);
        s2 += std::fabs(c.K2);
    }
    return L1{s0, s1, s2};
}

// The scale for K's ABSOLUTE error: one unit of binomial weight per subject,
// floored at one.  See kK0ProdUlp for the motivation.
double k0Scale(
    double t,
    const std::vector<double> &r,
    const std::vector<double> &p,
    const std::vector<double> &h
) {
    double s = 0.0;
    for (size_t i = 0; i < r.size(); ++i)
        s += h[i] * std::fmax(1.0, std::fabs(spa_cgf::logAlpha(t * r[i], p[i])));
    return s;
}

double k0Scale1(double t, double r, double p, double h) {
    return k0Scale(t, std::vector<double>{r}, std::vector<double>{p},
                   std::vector<double>{h});
}

// Worst deviation seen at ANY kK0ProdUlp assertion site in this file, fed by
// `k0ProdUlp` below and reported by `terminal_K_budget_headroom` at the end.
// A budget whose headroom is not printed is a budget nobody can tighten, and
// this one was carrying 146x when the file's stated convention is 2x to 4.6x.
double gWorstK0Prod = 0.0;
const char *gWorstK0ProdWhere = "(none)";

// Difference of a and b expressed in ULP of `scale`.  Exact equality is 0;
// a zero scale with unequal values is infinite, which fails any finite budget.
double ulpDiff(double a, double b, double scale) {
    if (a == b) return 0.0;
    if (std::isnan(a) || std::isnan(b))
        return std::numeric_limits<double>::infinity();
    if (!(scale > 0.0)) return std::numeric_limits<double>::infinity();
    return std::fabs(a - b) / (DBL_EPSILON * scale);
}

// The same quantity, recorded against the file-wide worst so the headroom on
// kK0ProdUlp is measured rather than assumed.  Every production-K absolute
// assertion in this file goes through here.
double k0ProdUlp(double got, double ref, double scale, const char *where) {
    const double u = ulpDiff(got, ref, scale);
    if (u > gWorstK0Prod) { gWorstK0Prod = u; gWorstK0ProdWhere = where; }
    return u;
}

// Deterministic subject sets.  A fixed seed is required: the suite must be
// reproducible between runs and between stages.
constexpr uint64_t kSeed = 20260730u;

struct Subjects {
    std::vector<double> r, af, hap;
};

Subjects makeSubjects(int n) {
    std::mt19937_64 rng(kSeed);
    std::normal_distribution<double> nd(0.0, 1.5);
    std::uniform_real_distribution<double> ud(1e-4, 0.5);
    std::uniform_int_distribution<int> hd(0, 2);

    Subjects s;
    s.r.reserve(static_cast<size_t>(n));
    s.af.reserve(static_cast<size_t>(n));
    s.hap.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        s.r.push_back(nd(rng));
        s.af.push_back(ud(rng));
        s.hap.push_back(static_cast<double>(hd(rng)));
    }
    return s;
}

// Subject counts exercising the empty set, every partial vector, exact lane
// multiples for both widths, and multiples +/- 1.
const std::vector<int> kNGrid =
    {0, 1, 2, 3, 4, 5, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 255,
     256, 257, 512};

} // namespace

// ──────────────────────────────────────────────────────────────────────
// II.a  The hapcount reference is the diploid reference generalized
// ──────────────────────────────────────────────────────────────────────

TEST(hapcount_reference_generalizes_the_diploid_reference) {
    for (double p : kMafs) {
        for (double r : kResidsP2) {
            for (double tr : kTrWide) {
                const double t = tr / r;
                const CgfL a = sparef::cgfReference(t, r, p);
                const RefL b = refHap(t, r, p, 2.0L);
                // Same expression with h substituted for the literal 2, so at
                // h = 2 the agreement must be exact, not merely close.
                CHECK_MSG(a.K0 == b.K0, "hapcount reference K diverged at h = 2");
                CHECK_MSG(a.K1 == b.K1, "hapcount reference K' diverged at h = 2");
                CHECK_MSG(a.K2 == b.K2, "hapcount reference K'' diverged at h = 2");
            }
        }
    }
}

// ──────────────────────────────────────────────────────────────────────
// II.b  Per-subject accuracy of the real kernels vs the long-double reference
// ──────────────────────────────────────────────────────────────────────
//
// K' and K'' are driven through the production reduction entry points with
// n = 1 rather than through the inline scalar core, so for those two the shipped
// code path is what is measured.
//
// K is NOT one number and must not be reported as one.  Two entry points produce
// it, only one of them ships, and they are held to two different criteria:
//
//   * `binom*KFull` — the production terminal K, dispatched, log(u + p*e^x)
//     through the vendored vector logarithm.  It is called by every SPA site.
//     The criterion it is held to is the ABSOLUTE error, expressed in ULP of
//     sum_i h_i*max(1, |log alpha_i|), because K reaches the p-value only
//     through w = sgn(zeta)*sqrt(2*(zeta*s - K)) with dw = -dK/w.  See
//     kK0ProdUlp above for the scale and `production_K_absolute_error_does_not_
//     degrade_w` for the argument executed end to end.
//   * `binom*KFullExact` — the N1 spelling, log1p(p*expm1(x)), scalar
//     throughout.  It has NO production caller anywhere in src/; it exists so
//     that relative accuracy in K stays reachable and pinned.  The criterion it
//     is held to is the RELATIVE error.
//
// A previous revision of this file printed a single line labelled "K (any) ...
// 1.0" beneath the heading "worst error in ULP", under a section comment
// claiming the shipped path was what was measured.  The 1.0 was `KFullExact`'s
// relative error.  The shipped path's worst RELATIVE error is 3.95e-08 on this
// grid — 1.8e+08 ULP, eight orders away from the printed figure — and 1.1e-01 on
// the finer grid of `real_kernel_K_accuracy_is_symmetric_about_zero`, which
// reaches t*r = 1e-9 where K itself is ~1e-13 and an absolute eps swamps it
// entirely.  The shipped terminal K has, near t = 0, no relative accuracy at
// all; that is the whole content of finding N1 and it is precisely why the
// absolute criterion is the one this file asserts on.  The old number was not
// wrong for what it measured, but the label was wrong for what it claimed, and a
// reader of `make test` output would have taken the shipped terminal K to be
// correctly rounded.  Both entry points are therefore reported below, each
// against its own criterion, with the production path's relative error printed
// alongside and marked as a record rather than a bound.

// Accumulates, per tier, the worst relative error seen against the reference.
// K' and K'' only: kFull reuses k12's derivatives bit-for-bit, which the tests
// below assert separately, so one set of tier numbers covers both entry points.
struct TierErr { double k1 = 0.0, k2 = 0.0; };

// The terminal K, worst over the grid.  Three numbers because there are three
// distinct claims, and collapsing them is what produced the misleading report.
struct K0Err {
    double prodAbsUlp = 0.0;  // production kFull, ULP of the absolute k0Scale
    double prodRel    = 0.0;  // production kFull, relative — recorded, not a bound
    double exactRel   = 0.0;  // kFullExact, relative — the N1 criterion
};

void reportTierErrs(
    const char *variant,
    const std::vector<std::string> &names,
    const std::vector<TierErr> &errs,
    const K0Err &k0
) {
    std::fprintf(stderr,
                 "      %s: K' and K'' vs long-double reference, worst RELATIVE "
                 "error in ULP\n",
                 variant);
    for (size_t i = 0; i < names.size(); ++i)
        std::fprintf(stderr, "        %-9s K' %6.1f   K'' %6.1f\n",
                     names[i].c_str(), errs[i].k1 / DBL_EPSILON,
                     errs[i].k2 / DBL_EPSILON);
    std::fprintf(stderr,
                 "      %s: terminal K vs long-double reference, worst error\n",
                 variant);
    std::fprintf(stderr,
                 "        kFull      SHIPPED, dispatched   ABSOLUTE %8.1f ULP of "
                 "sum h*max(1,|log alpha|)   (budget %.0f)\n",
                 k0.prodAbsUlp, kK0ProdUlp);
    std::fprintf(stderr,
                 "        kFull      SHIPPED, dispatched   relative %8.2e "
                 "= %8.2e ULP   (recorded, not a criterion)\n",
                 k0.prodRel, k0.prodRel / DBL_EPSILON);
    std::fprintf(stderr,
                 "        kFullExact no production caller  relative %8.1f ULP"
                 "                             (budget %.0f)\n",
                 k0.exactRel / DBL_EPSILON, kAlgebraUlp);
}

TEST(real_binom_uniform_matches_long_double_reference) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    const auto tiers = uniformTiers();
    std::vector<TierErr> errs(tiers.size());
    std::vector<std::string> names;
    for (const auto &tt : tiers) names.push_back(tt.name);
    K0Err k0;

    for (double p : kMafs) {
        for (double r : kResidsP2) {
            for (double tr : kTrWide) {
                const double t = tr / r;
                const double rv[1] = {r};
                const CgfL ref = sparef::cgfReference(t, r, p);
                const double R1 = static_cast<double>(ref.K1);
                const double R2 = static_cast<double>(ref.K2);

                for (size_t k = 0; k < tiers.size(); ++k) {
                    const spa_cgf::Cgf12 d = tiers[k].fn(t, rv, 1, p);
                    const double budget =
                        usesVectorExp(tiers[k].name, 1) ? kVectorUlp : kAlgebraUlp;
                    errs[k].k1 = std::fmax(errs[k].k1, tinytest::relDiff(d.K1, R1));
                    errs[k].k2 = std::fmax(errs[k].k2, tinytest::relDiff(d.K2, R2));
                    CHECK_REL(d.K1, R1, relTol(budget));
                    CHECK_REL(d.K2, R2, relTol(budget));
                }

                // kFull must reuse k12's derivatives verbatim: they come from
                // the same dispatched kernel, so any difference is a bug.  Both
                // terminal entry points are checked, since both are documented
                // to differ from k12 in K alone.
                const spa_cgf::Cgf12  d = spa_cgf::binomUniformK12(t, rv, 1, p);
                const spa_cgf::Cgf012 f = spa_cgf::binomUniformKFull(t, rv, 1, p);
                const spa_cgf::Cgf012 x =
                    spa_cgf::binomUniformKFullExact(t, rv, 1, p);
                CHECK_MSG(f.K1 == d.K1, "kFull K' differs from k12 K'");
                CHECK_MSG(f.K2 == d.K2, "kFull K'' differs from k12 K''");
                CHECK_MSG(x.K1 == d.K1, "kFullExact K' differs from k12 K'");
                CHECK_MSG(x.K2 == d.K2, "kFullExact K'' differs from k12 K''");

                const double R0 = static_cast<double>(ref.K0);
                // The exact entry point keeps N1's relative accuracy in K.  It
                // is scalar in every tier, so it carries the algebra budget
                // regardless of the dispatched tier, and it passes through zero
                // at t = 0 where only a mixed criterion is satisfiable.
                CHECK_CLOSE(x.K0, R0, relTol(kAlgebraUlp), 1e-300);
                // The production entry point is held to the absolute criterion
                // K is actually consumed under.
                const double uProd =
                    k0ProdUlp(f.K0, R0, k0Scale1(t, r, p, 2.0), "binomUniform n=1");
                CHECK_MSG(uProd <= kK0ProdUlp, "production K off the absolute scale");
                k0.prodAbsUlp = std::fmax(k0.prodAbsUlp, uProd);
                if (R0 != 0.0) {
                    k0.exactRel = std::fmax(k0.exactRel, tinytest::relDiff(x.K0, R0));
                    k0.prodRel  = std::fmax(k0.prodRel,  tinytest::relDiff(f.K0, R0));
                }
            }
        }
    }

    reportTierErrs("binomUniform", names, errs, k0);
}

TEST(real_binom_indiv_matches_long_double_reference) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    const auto tiers = indivTiers();
    std::vector<TierErr> errs(tiers.size());
    std::vector<std::string> names;
    for (const auto &tt : tiers) names.push_back(tt.name);
    K0Err k0;

    for (double p : kMafs) {
        for (double r : kResidsP2) {
            for (double tr : kTrWide) {
                const double t = tr / r;
                const double rv[1] = {r};
                const double pv[1] = {p};
                const CgfL ref = sparef::cgfReference(t, r, p);
                const double R1 = static_cast<double>(ref.K1);
                const double R2 = static_cast<double>(ref.K2);

                for (size_t k = 0; k < tiers.size(); ++k) {
                    const spa_cgf::Cgf12 d = tiers[k].fn(t, rv, pv, 1);
                    const double budget =
                        usesVectorExp(tiers[k].name, 1) ? kVectorUlp : kAlgebraUlp;
                    errs[k].k1 = std::fmax(errs[k].k1, tinytest::relDiff(d.K1, R1));
                    errs[k].k2 = std::fmax(errs[k].k2, tinytest::relDiff(d.K2, R2));
                    CHECK_REL(d.K1, R1, relTol(budget));
                    CHECK_REL(d.K2, R2, relTol(budget));
                }

                const spa_cgf::Cgf12  d = spa_cgf::binomIndivK12(t, rv, pv, 1);
                const spa_cgf::Cgf012 f = spa_cgf::binomIndivKFull(t, rv, pv, 1);
                const spa_cgf::Cgf012 x =
                    spa_cgf::binomIndivKFullExact(t, rv, pv, 1);
                CHECK_MSG(f.K1 == d.K1, "kFull K' differs from k12 K'");
                CHECK_MSG(f.K2 == d.K2, "kFull K'' differs from k12 K''");
                CHECK_MSG(x.K1 == d.K1, "kFullExact K' differs from k12 K'");
                CHECK_MSG(x.K2 == d.K2, "kFullExact K'' differs from k12 K''");

                const double R0 = static_cast<double>(ref.K0);
                CHECK_CLOSE(x.K0, R0, relTol(kAlgebraUlp), 1e-300);
                const double uProd =
                    k0ProdUlp(f.K0, R0, k0Scale1(t, r, p, 2.0), "binomIndiv n=1");
                CHECK_MSG(uProd <= kK0ProdUlp, "production K off the absolute scale");
                k0.prodAbsUlp = std::fmax(k0.prodAbsUlp, uProd);
                if (R0 != 0.0) {
                    k0.exactRel = std::fmax(k0.exactRel, tinytest::relDiff(x.K0, R0));
                    k0.prodRel  = std::fmax(k0.prodRel,  tinytest::relDiff(f.K0, R0));
                }
            }
        }
    }

    reportTierErrs("binomIndiv", names, errs, k0);
}

TEST(real_binom_hapcount_matches_long_double_reference) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    const auto tiers = hapcountTiers();
    std::vector<TierErr> errs(tiers.size());
    std::vector<std::string> names;
    for (const auto &tt : tiers) names.push_back(tt.name);
    K0Err k0;

    for (double q : kMafs) {
        for (double r : kResidsP2) {
            for (double tr : kTrWide) {
                for (double h : kHaps) {
                    const double t = tr / r;
                    const double rv[1] = {r};
                    const double hv[1] = {h};
                    const RefL ref = refHap(t, r, q, h);
                    const double R0 = static_cast<double>(ref.K0);
                    const double R1 = static_cast<double>(ref.K1);
                    const double R2 = static_cast<double>(ref.K2);

                    for (size_t k = 0; k < tiers.size(); ++k) {
                        const spa_cgf::Cgf12 d = tiers[k].fn(t, rv, hv, 1, q);
                        // h = 0 makes every cumulant identically zero; a
                        // relative criterion is vacuous there, exactness is not.
                        if (h == 0.0) {
                            CHECK_MSG(d.K1 == 0.0 && d.K2 == 0.0,
                                      "h = 0 must give exactly zero derivatives");
                            continue;
                        }
                        const double budget =
                            usesVectorExp(tiers[k].name, 1) ? kVectorUlp : kAlgebraUlp;
                        errs[k].k1 = std::fmax(errs[k].k1, tinytest::relDiff(d.K1, R1));
                        errs[k].k2 = std::fmax(errs[k].k2, tinytest::relDiff(d.K2, R2));
                        CHECK_REL(d.K1, R1, relTol(budget));
                        CHECK_REL(d.K2, R2, relTol(budget));
                    }

                    const spa_cgf::Cgf12  d = spa_cgf::binomHapcountK12(t, rv, hv, 1, q);
                    const spa_cgf::Cgf012 f = spa_cgf::binomHapcountKFull(t, rv, hv, 1, q);
                    const spa_cgf::Cgf012 x =
                        spa_cgf::binomHapcountKFullExact(t, rv, hv, 1, q);
                    CHECK_MSG(f.K1 == d.K1, "kFull K' differs from k12 K'");
                    CHECK_MSG(f.K2 == d.K2, "kFull K'' differs from k12 K''");
                    CHECK_MSG(x.K1 == d.K1, "kFullExact K' differs from k12 K'");
                    CHECK_MSG(x.K2 == d.K2, "kFullExact K'' differs from k12 K''");

                    if (h == 0.0) {
                        CHECK_MSG(f.K0 == 0.0, "h = 0 must give exactly zero K");
                        CHECK_MSG(x.K0 == 0.0, "h = 0 must give exactly zero K");
                        continue;
                    }
                    CHECK_CLOSE(x.K0, R0, relTol(kAlgebraUlp), 1e-300);
                    const double uProd =
                        k0ProdUlp(f.K0, R0, k0Scale1(t, r, q, h), "binomHapcount n=1");
                    CHECK_MSG(uProd <= kK0ProdUlp,
                              "production K off the absolute scale");
                    k0.prodAbsUlp = std::fmax(k0.prodAbsUlp, uProd);
                    if (R0 != 0.0) {
                        k0.exactRel =
                            std::fmax(k0.exactRel, tinytest::relDiff(x.K0, R0));
                        k0.prodRel =
                            std::fmax(k0.prodRel, tinytest::relDiff(f.K0, R0));
                    }
                }
            }
        }
    }

    reportTierErrs("binomHapcount", names, errs, k0);
}

// ──────────────────────────────────────────────────────────────────────
// II.c  The reductions, not only the closed form
// ──────────────────────────────────────────────────────────────────────
//
// Accumulating n per-subject contributions is where a masked tail, a scalar
// tail or a horizontal reduction can go wrong, and none of that is reachable
// from an n = 1 test.  The comparison is against the long-double sum of the
// per-subject references, under the L1-scaled criterion motivated above.

TEST(real_kernel_reductions_match_long_double_sums) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    double worstU = 0.0, worstI = 0.0, worstH = 0.0;

    for (int n : kNGrid) {
        const Subjects s = makeSubjects(n);
        const std::vector<double> twos(static_cast<size_t>(n), 2.0);

        for (double t : {-0.7, -0.05, 0.0, 0.05, 0.7, 3.0}) {
            // ── binomUniform ──
            {
                const double p = 0.23;
                const std::vector<double> pv(static_cast<size_t>(n), p);
                long double R1 = 0.0L, R2 = 0.0L, R0 = 0.0L;
                for (int i = 0; i < n; ++i) {
                    const RefL x = refHap(t, s.r[i], p, 2.0L);
                    R0 += x.K0; R1 += x.K1; R2 += x.K2;
                }
                const spa_cgf::Cgf012 f =
                    spa_cgf::binomUniformKFull(t, s.r.data(), n, p);
                const L1 sc = l1Scale(t, s.r, pv, twos);

                const double u1 = ulpDiff(f.K1, static_cast<double>(R1), sc.K1);
                const double u2 = ulpDiff(f.K2, static_cast<double>(R2), sc.K2);
                worstU = std::fmax(worstU, std::fmax(u1, u2));
                CHECK_MSG(u1 <= kCrossTierUlp, "binomUniform reduction K' off scale");
                CHECK_MSG(u2 <= kCrossTierUlp, "binomUniform reduction K'' off scale");
                if (n > 0)
                    CHECK_MSG(ulpDiff(f.K0, static_cast<double>(R0), sc.K0)
                                  <= kCrossTierUlp,
                              "reduction K0 off the L1 scale");
            }

            // ── binomIndiv ──
            {
                long double R1 = 0.0L, R2 = 0.0L, R0 = 0.0L;
                for (int i = 0; i < n; ++i) {
                    const RefL x = refHap(t, s.r[i], s.af[i], 2.0L);
                    R0 += x.K0; R1 += x.K1; R2 += x.K2;
                }
                const spa_cgf::Cgf012 f =
                    spa_cgf::binomIndivKFull(t, s.r.data(), s.af.data(), n);
                const L1 sc = l1Scale(t, s.r, s.af, twos);

                const double u1 = ulpDiff(f.K1, static_cast<double>(R1), sc.K1);
                const double u2 = ulpDiff(f.K2, static_cast<double>(R2), sc.K2);
                worstI = std::fmax(worstI, std::fmax(u1, u2));
                CHECK_MSG(u1 <= kCrossTierUlp, "binomIndiv reduction K' off scale");
                CHECK_MSG(u2 <= kCrossTierUlp, "binomIndiv reduction K'' off scale");
                if (n > 0)
                    CHECK_MSG(ulpDiff(f.K0, static_cast<double>(R0), sc.K0)
                                  <= kCrossTierUlp,
                              "reduction K0 off the L1 scale");
            }

            // ── binomHapcount ──
            {
                const double q = 0.41;
                const std::vector<double> qv(static_cast<size_t>(n), q);
                long double R1 = 0.0L, R2 = 0.0L, R0 = 0.0L;
                for (int i = 0; i < n; ++i) {
                    const RefL x = refHap(t, s.r[i], q, s.hap[i]);
                    R0 += x.K0; R1 += x.K1; R2 += x.K2;
                }
                const spa_cgf::Cgf012 f =
                    spa_cgf::binomHapcountKFull(t, s.r.data(), s.hap.data(), n, q);
                const L1 sc = l1Scale(t, s.r, qv, s.hap);

                const double u1 = ulpDiff(f.K1, static_cast<double>(R1), sc.K1);
                const double u2 = ulpDiff(f.K2, static_cast<double>(R2), sc.K2);
                worstH = std::fmax(worstH, std::fmax(u1, u2));
                CHECK_MSG(u1 <= kCrossTierUlp, "binomHapcount reduction K' off scale");
                CHECK_MSG(u2 <= kCrossTierUlp, "binomHapcount reduction K'' off scale");
                if (n > 0)
                    CHECK_MSG(ulpDiff(f.K0, static_cast<double>(R0), sc.K0)
                                  <= kCrossTierUlp,
                              "reduction K0 off the L1 scale");
            }
        }
    }

    std::fprintf(stderr,
                 "      reduction vs long-double sum, worst deviation in ULP of the "
                 "L1 scale (budget %.0f)\n"
                 "        binomUniform %.2f, binomIndiv %.2f, binomHapcount %.2f\n",
                 kCrossTierUlp, worstU, worstI, worstH);
}

// ──────────────────────────────────────────────────────────────────────
// II.d  Cross-tier equality
// ──────────────────────────────────────────────────────────────────────
//
// Every tier the host supports, plus the dispatched entry point, against the
// scalar tier as reference.  Bit-identity is not the expectation and would be
// the wrong assertion; see the kCrossTierUlp commentary.

TEST(cross_tier_equality_binom_uniform) {
    double worst = 0.0;
    int nComparisons = 0;
    const auto tiers = uniformTiers();

    for (int n : kNGrid) {
        const Subjects s = makeSubjects(n);
        const std::vector<double> twos(static_cast<size_t>(n), 2.0);

        for (double p : {1e-4, 0.05, 0.23, 0.5, 0.97}) {
            const std::vector<double> pv(static_cast<size_t>(n), p);
            for (double t : {-3.0, -0.7, -0.05, 0.0, 0.05, 0.7, 3.0}) {
                const spa_cgf::Cgf12 ref =
                    spa_cgf::tier::binomUniformK12_scalar(t, s.r.data(), n, p);
                const L1 sc = l1Scale(t, s.r, pv, twos);

                for (const auto &tt : tiers) {
                    const spa_cgf::Cgf12 got = tt.fn(t, s.r.data(), n, p);
                    const double u1 = ulpDiff(got.K1, ref.K1, sc.K1);
                    const double u2 = ulpDiff(got.K2, ref.K2, sc.K2);
                    worst = std::fmax(worst, std::fmax(u1, u2));
                    ++nComparisons;
                    CHECK_MSG(u1 <= kCrossTierUlp, std::string("K' tier ") + tt.name);
                    CHECK_MSG(u2 <= kCrossTierUlp, std::string("K'' tier ") + tt.name);
                }
            }
        }
    }

    std::fprintf(stderr,
                 "      binomUniform: %d comparisons over %zu tiers, worst %.2f ULP "
                 "of L1 scale (budget %.0f)\n",
                 nComparisons, tiers.size(), worst, kCrossTierUlp);
}

TEST(cross_tier_equality_binom_indiv) {
    double worst = 0.0;
    int nComparisons = 0;
    const auto tiers = indivTiers();

    for (int n : kNGrid) {
        const Subjects s = makeSubjects(n);
        const std::vector<double> twos(static_cast<size_t>(n), 2.0);

        for (double t : {-3.0, -0.7, -0.05, 0.0, 0.05, 0.7, 3.0}) {
            const spa_cgf::Cgf12 ref =
                spa_cgf::tier::binomIndivK12_scalar(t, s.r.data(), s.af.data(), n);
            const L1 sc = l1Scale(t, s.r, s.af, twos);

            for (const auto &tt : tiers) {
                const spa_cgf::Cgf12 got = tt.fn(t, s.r.data(), s.af.data(), n);
                const double u1 = ulpDiff(got.K1, ref.K1, sc.K1);
                const double u2 = ulpDiff(got.K2, ref.K2, sc.K2);
                worst = std::fmax(worst, std::fmax(u1, u2));
                ++nComparisons;
                CHECK_MSG(u1 <= kCrossTierUlp, std::string("K' tier ") + tt.name);
                CHECK_MSG(u2 <= kCrossTierUlp, std::string("K'' tier ") + tt.name);
            }
        }
    }

    std::fprintf(stderr,
                 "      binomIndiv:   %d comparisons over %zu tiers, worst %.2f ULP "
                 "of L1 scale (budget %.0f)\n",
                 nComparisons, tiers.size(), worst, kCrossTierUlp);
}

TEST(cross_tier_equality_binom_hapcount) {
    double worst = 0.0;
    int nComparisons = 0;
    const auto tiers = hapcountTiers();

    for (int n : kNGrid) {
        const Subjects s = makeSubjects(n);

        for (double q : {1e-4, 0.05, 0.41, 0.5, 0.97}) {
            const std::vector<double> qv(static_cast<size_t>(n), q);
            for (double t : {-3.0, -0.7, -0.05, 0.0, 0.05, 0.7, 3.0}) {
                const spa_cgf::Cgf12 ref = spa_cgf::tier::binomHapcountK12_scalar(
                    t, s.r.data(), s.hap.data(), n, q);
                const L1 sc = l1Scale(t, s.r, qv, s.hap);

                for (const auto &tt : tiers) {
                    const spa_cgf::Cgf12 got =
                        tt.fn(t, s.r.data(), s.hap.data(), n, q);
                    const double u1 = ulpDiff(got.K1, ref.K1, sc.K1);
                    const double u2 = ulpDiff(got.K2, ref.K2, sc.K2);
                    worst = std::fmax(worst, std::fmax(u1, u2));
                    ++nComparisons;
                    CHECK_MSG(u1 <= kCrossTierUlp, std::string("K' tier ") + tt.name);
                    CHECK_MSG(u2 <= kCrossTierUlp, std::string("K'' tier ") + tt.name);
                }
            }
        }
    }

    std::fprintf(stderr,
                 "      binomHapcount: %d comparisons over %zu tiers, worst %.2f ULP "
                 "of L1 scale (budget %.0f)\n",
                 nComparisons, tiers.size(), worst, kCrossTierUlp);
}

// ──────────────────────────────────────────────────────────────────────
// II.d'  Cross-tier equality of the terminal K
// ──────────────────────────────────────────────────────────────────────
//
// CLAUDE.md requires the full scalar + AVX2 + AVX-512 triple for every SIMD
// kernel; the terminal K is now one, so it gets the same cross-tier treatment
// the k12 kernels get.  The scale is the absolute one (k0Scale) rather than the
// L1 one, for the reason set out at kK0ProdUlp: near t = 0 the L1 scale vanishes
// while the per-lane error does not.
TEST(cross_tier_equality_terminal_K) {
    double worst = 0.0;
    int nComparisons = 0;

    const auto uT = uniformK0Tiers();
    const auto iT = indivK0Tiers();
    const auto hT = hapcountK0Tiers();

    for (int n : kNGrid) {
        const Subjects s = makeSubjects(n);
        const std::vector<double> twos(static_cast<size_t>(n), 2.0);

        for (double t : {-3.0, -0.7, -0.05, 0.0, 0.05, 0.7, 3.0}) {
            // ── binomUniform ──
            for (double p : {1e-4, 0.05, 0.23, 0.5, 0.97}) {
                const std::vector<double> pv(static_cast<size_t>(n), p);
                const double ref =
                    spa_cgf::tier::binomUniformK0_scalar(t, s.r.data(), n, p);
                const double sc = k0Scale(t, s.r, pv, twos);
                for (const auto &tt : uT) {
                    const double got = tt.fn(t, s.r.data(), n, p);
                    const double u = ulpDiff(2.0 * got, 2.0 * ref, sc);
                    worst = std::fmax(worst, u);
                    ++nComparisons;
                    CHECK_MSG(u <= kCrossTierUlp,
                              std::string("binomUniform K tier ") + tt.name);
                }
            }

            // ── binomIndiv ──
            {
                const double ref = spa_cgf::tier::binomIndivK0_scalar(
                    t, s.r.data(), s.af.data(), n);
                const double sc = k0Scale(t, s.r, s.af, twos);
                for (const auto &tt : iT) {
                    const double got = tt.fn(t, s.r.data(), s.af.data(), n);
                    const double u = ulpDiff(2.0 * got, 2.0 * ref, sc);
                    worst = std::fmax(worst, u);
                    ++nComparisons;
                    CHECK_MSG(u <= kCrossTierUlp,
                              std::string("binomIndiv K tier ") + tt.name);
                }
            }

            // ── binomHapcount ──
            for (double q : {1e-4, 0.41, 0.97}) {
                const std::vector<double> qv(static_cast<size_t>(n), q);
                const double ref = spa_cgf::tier::binomHapcountK0_scalar(
                    t, s.r.data(), s.hap.data(), n, q);
                const double sc = k0Scale(t, s.r, qv, s.hap);
                for (const auto &tt : hT) {
                    const double got = tt.fn(t, s.r.data(), s.hap.data(), n, q);
                    const double u = ulpDiff(got, ref, sc);
                    worst = std::fmax(worst, u);
                    ++nComparisons;
                    CHECK_MSG(u <= kCrossTierUlp,
                              std::string("binomHapcount K tier ") + tt.name);
                }
            }
        }
    }

    // The triple must actually have been exercised on this host.
    CHECK_MSG(uT.size() == static_cast<size_t>(kLevel) + 1,
              "expected one K tier per supported ISA level");
    std::fprintf(stderr,
                 "      terminal K: %d comparisons over %zu tiers, worst %.2f ULP "
                 "of the absolute scale (budget %.0f)\n",
                 nComparisons, uT.size(), worst, kCrossTierUlp);
}

// ──────────────────────────────────────────────────────────────────────
// II.d''  The numerical argument for the cheaper terminal K
// ──────────────────────────────────────────────────────────────────────
//
// Stage 4 replaced sum_i h*log1p(p*expm1(t*r_i)) with sum_i h*log(u + p*e^{t*r_i})
// and vectorized it, on the argument that K is consumed only through
// w = sgn(zeta)*sqrt(2*(zeta*s - K)) at abscissae where the |z| > spaCutoff entry
// gate keeps zeta*s - K bounded away from zero, so the ABSOLUTE error in K is
// what propagates and the two spellings are equivalent on that metric.
//
// This test is that argument, executed.  For each population it locates the
// saddlepoint in long double, evaluates K there both ways and in long double,
// and compares the absolute error in K, the error in w, and the error in
// -log10(p) through the Barndorff-Nielsen tail.  The claim is not that the
// production spelling is as good relatively — it is not, and
// `real_kernel_K_accuracy_near_zero` tabulates by how much — but that the
// difference does not reach w.
TEST(production_K_absolute_error_does_not_degrade_w) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    std::mt19937_64 rng(20260731u);
    std::normal_distribution<double> nd(0.0, 1.0);

    double maxAbsExact = 0.0, maxAbsProd = 0.0;
    double maxWExact = 0.0, maxWProd = 0.0;
    double maxGExact = 0.0, maxGProd = 0.0;
    double minAbsW = std::numeric_limits<double>::infinity();
    int nCases = 0;

    for (double p : {5e-4, 0.01, 0.1, 0.5}) {
        for (double rsc : {0.5, 2.0, 6.0}) {
            for (double z : {2.0, 3.0, 6.0, 10.0}) {
                for (int sgn = -1; sgn <= 1; sgn += 2) {
                    const int n = 240;
                    std::vector<double> r(static_cast<size_t>(n));
                    double sumR2 = 0.0;
                    for (int i = 0; i < n; ++i) {
                        r[i] = nd(rng) * rsc;
                        sumR2 += r[i] * r[i];
                    }
                    // The analytic non-outlier Gaussian block every production
                    // site adds, three times the outlier block's own variance.
                    const double var = 2.0 * p * (1.0 - p) * sumR2 * 3.0;

                    // K'(t) and K''(t) of the whole CGF, in long double.
                    auto k12L = [&](long double t, long double &K1, long double &K2) {
                        K1 = static_cast<long double>(var) * t;
                        K2 = static_cast<long double>(var);
                        const long double pl = p, ul = 1.0L - pl;
                        for (int i = 0; i < n; ++i) {
                            const long double ri = r[i];
                            const long double e = pl * std::expl(t * ri);
                            const long double a = ul + e;
                            K1 += 2.0L * ri * e / a;
                            K2 += 2.0L * ri * ri * e * ul / (a * a);
                        }
                    };

                    long double K1, K2;
                    k12L(0.0L, K1, K2);
                    const double s =
                        static_cast<double>(K1 + static_cast<long double>(sgn) *
                                                     static_cast<long double>(z) *
                                                     std::sqrtl(K2));

                    // Bisect for the root of K'(t) = s in long double.
                    long double lo = -600.0L, hi = 600.0L;
                    for (int it = 0; it < 200; ++it) {
                        const long double mid = 0.5L * (lo + hi);
                        k12L(mid, K1, K2);
                        if (K1 - static_cast<long double>(s) < 0.0L) lo = mid;
                        else hi = mid;
                    }
                    const long double zetaL = 0.5L * (lo + hi);
                    const double zeta = static_cast<double>(zetaL);
                    k12L(zetaL, K1, K2);
                    const double K2d = static_cast<double>(K2);
                    if (!(K2d > 0.0)) continue;

                    // K at the root: long-double reference, and both entry points.
                    long double KrefOut = static_cast<long double>(var) * zetaL *
                                          zetaL * 0.5L;
                    for (int i = 0; i < n; ++i)
                        KrefOut += 2.0L * std::log1pl(
                            static_cast<long double>(p) *
                            std::expm1l(zetaL * static_cast<long double>(r[i])));
                    const double gaussK = static_cast<double>(
                        static_cast<long double>(var) * zetaL * zetaL * 0.5L);

                    const double kProd =
                        spa_cgf::binomUniformKFull(zeta, r.data(), n, p).K0 + gaussK;
                    const double kExact =
                        spa_cgf::binomUniformKFullExact(zeta, r.data(), n, p).K0 +
                        gaussK;

                    const long double tempRef =
                        zetaL * static_cast<long double>(s) - KrefOut;
                    if (!(tempRef > 0.0L)) continue;
                    ++nCases;

                    const double wRef = static_cast<double>(
                        (zetaL > 0 ? 1.0L : -1.0L) * std::sqrtl(2.0L * tempRef));
                    minAbsW = std::fmin(minAbsW, std::fabs(wRef));

                    auto wOf = [&](double K0) {
                        const double tmp = zeta * s - K0;
                        return spa::signOf(zeta) * std::sqrt(2.0 * tmp);
                    };
                    auto gOf = [&](double K0) {
                        spa::Status st;
                        return -spa::bnTailLog(zeta, s, K0, K2d, sgn < 0, st) /
                               2.30258509299404568402;
                    };
                    const double gRef = gOf(static_cast<double>(KrefOut));

                    maxAbsExact = std::fmax(
                        maxAbsExact,
                        std::fabs(static_cast<double>(static_cast<long double>(kExact) -
                                                      KrefOut)));
                    maxAbsProd = std::fmax(
                        maxAbsProd,
                        std::fabs(static_cast<double>(static_cast<long double>(kProd) -
                                                      KrefOut)));
                    maxWExact = std::fmax(maxWExact, std::fabs(wOf(kExact) - wRef));
                    maxWProd  = std::fmax(maxWProd,  std::fabs(wOf(kProd)  - wRef));
                    maxGExact = std::fmax(maxGExact, std::fabs(gOf(kExact) - gRef));
                    maxGProd  = std::fmax(maxGProd,  std::fabs(gOf(kProd)  - gRef));
                }
            }
        }
    }

    CHECK_MSG(nCases >= 90, "only " + std::to_string(nCases) + " usable cases");

    std::fprintf(stderr,
                 "\n      terminal K at the saddlepoint, n = 240, %d cases; "
                 "min |w| reached %.3f\n"
                 "      %-24s %12s %12s %14s\n"
                 "      %-24s %12.3e %12.3e %14.3e\n"
                 "      %-24s %12.3e %12.3e %14.3e\n\n",
                 nCases, minAbsW,
                 "spelling", "max |dK|", "max |dw|", "max |dLOG10P|",
                 "log1p(p*expm1(x))", maxAbsExact, maxWExact, maxGExact,
                 "log(u + p*e^x), vector", maxAbsProd, maxWProd, maxGProd);

    // (a) the entry gate really does keep w away from zero, which is the
    //     premise the whole argument rests on.
    CHECK_MSG(minAbsW > 1.5, "min |w| " + std::to_string(minAbsW));

    // (b) the production spelling's absolute error in K is within a small
    //     constant factor of the exact spelling's.  Measured worst ratio on this
    //     grid is 9.0 (1.48e-13 against 1.64e-14); the bound is 20.  Compare the
    //     factor separating their RELATIVE errors, which reaches 1/p ~ 2e3 at
    //     p = 5e-4: it is precisely because the two metrics disagree by two
    //     orders of magnitude that the absolute one has to be the one asserted.
    char buf[256];
    std::snprintf(buf, sizeof(buf), "production |dK| %.4g vs exact %.4g",
                  maxAbsProd, maxAbsExact);
    CHECK_MSG(maxAbsProd <= 20.0 * maxAbsExact + 1e-15, std::string(buf));

    // (c) what actually matters: the error reaching -log10(p) is negligible on
    //     any scale a GWAS result is read at.  The stage budget was 1e-3, and
    //     the measured worst is 6.8e-14.
    std::snprintf(buf, sizeof(buf), "production |dLOG10P| %.4g, |dw| %.4g",
                  maxGProd, maxWProd);
    CHECK_MSG(maxGProd <= 1e-9, std::string(buf));
    CHECK_MSG(maxWProd <= 1e-9, std::string(buf));
}

// The domain clamp exists so the tiers agree past the point where a scalar
// std::exp returns +inf while the vectorized kernels clamp internally.  Without
// it the scalar tier yields h*r*inf*0 = NaN and the vector tiers a finite value.
TEST(clamp_makes_tiers_agree_past_the_exp_domain) {
    const auto tiers = uniformTiers();
    const std::vector<double> r = {-40.0, -25.0, 25.0, 40.0, 1.0, -1.0, 12.0, -12.0};
    const std::vector<double> twos(r.size(), 2.0);

    int nChecked = 0;
    for (double p : {1e-4, 0.3, 0.5}) {
        const std::vector<double> pv(r.size(), p);
        // t*r reaches +/- 1200 here, far outside exp's domain.
        for (double t : {-30.0, -20.0, 20.0, 30.0}) {
            const spa_cgf::Cgf12 ref = spa_cgf::tier::binomUniformK12_scalar(
                t, r.data(), static_cast<int>(r.size()), p);
            const L1 sc = l1Scale(t, r, pv, twos);

            CHECK_MSG(std::isfinite(ref.K1) && std::isfinite(ref.K2),
                      "clamped scalar kernel must stay finite past exp's domain");
            CHECK_MSG(ref.K2 >= 0.0, "K'' must not go negative under the clamp");

            for (const auto &tt : tiers) {
                const spa_cgf::Cgf12 got =
                    tt.fn(t, r.data(), static_cast<int>(r.size()), p);
                CHECK_MSG(std::isfinite(got.K1) && std::isfinite(got.K2),
                          std::string("non-finite past exp domain, tier ") + tt.name);
                CHECK_MSG(ulpDiff(got.K1, ref.K1, sc.K1) <= kCrossTierUlp,
                          std::string("K' clamp disagreement, tier ") + tt.name);
                CHECK_MSG(ulpDiff(got.K2, ref.K2, sc.K2) <= kCrossTierUlp,
                          std::string("K'' clamp disagreement, tier ") + tt.name);
                ++nChecked;
            }
        }
    }
    std::fprintf(stderr,
                 "      %d tier evaluations at |t*r| up to 1200: all finite and in "
                 "agreement\n", nChecked);
}

// ──────────────────────────────────────────────────────────────────────
// II.e  K'' strictly positive and finite where the production form fails
// ──────────────────────────────────────────────────────────────────────
//
// The saddlepoint solver divides by K'' on every Newton iteration, so a K''
// that underflows to zero or turns negative through cancellation produces an
// infinite or NaN step.  The stable algebra is positive by construction: it is
// a product of non-negative factors, and critically 1-pi is formed as the
// quotient u/alpha rather than as a subtraction, so it stays strictly positive
// exactly where `1 - pi` would round to zero.
//
// The grid is the one Part I measures the production form on, where it yields 20
// non-positive or non-finite values out of 288.

TEST(real_kernel_curvature_is_strictly_positive) {
    const std::vector<double> rs = {-8.0, -1.0, 1.0, 8.0};
    const std::vector<double> trs =
        {-60.0, -40.0, -20.0, -5.0, 0.0, 5.0, 20.0, 40.0, 60.0};

    int nProductionBad = 0, nUniformBad = 0, nIndivBad = 0, nHapBad = 0, nPoints = 0;

    for (double p : kMafs) {
        for (double r : rs) {
            for (double tr : trs) {
                const double t = tr / r;
                const double rv[1] = {r};
                const double pv[1] = {p};
                const double hv[1] = {2.0};
                ++nPoints;

                const spa_cgf::Cgf12 u = spa_cgf::binomUniformK12(t, rv, 1, p);
                const spa_cgf::Cgf12 i = spa_cgf::binomIndivK12(t, rv, pv, 1);
                const spa_cgf::Cgf12 h = spa_cgf::binomHapcountK12(t, rv, hv, 1, p);

                if (!(u.K2 > 0.0) || !std::isfinite(u.K2)) ++nUniformBad;
                if (!(i.K2 > 0.0) || !std::isfinite(i.K2)) ++nIndivBad;
                if (!(h.K2 > 0.0) || !std::isfinite(h.K2)) ++nHapBad;

                CHECK_MSG(u.K2 > 0.0 && std::isfinite(u.K2),
                          "binomUniform K'' must be finite and strictly positive");
                CHECK_MSG(i.K2 > 0.0 && std::isfinite(i.K2),
                          "binomIndiv K'' must be finite and strictly positive");
                CHECK_MSG(h.K2 > 0.0 && std::isfinite(h.K2),
                          "binomHapcount K'' must be finite and strictly positive");

                // The three variants describe the same law at h = 2, so they
                // must also agree numerically with each other.
                CHECK_REL(u.K2, i.K2, relTol(kAlgebraUlp));
                CHECK_REL(u.K2, h.K2, relTol(kAlgebraUlp));

                const Cgf d = sparef::cgfDiffForm(t, r, p);
                if (!(d.K2 > 0.0) || !std::isfinite(d.K2)) ++nProductionBad;
            }
        }
    }

    std::fprintf(stderr,
                 "\n      non-positive or non-finite K'' over %d sampled points\n"
                 "        production (differencing form)  %d\n"
                 "        binomUniform                   %d\n"
                 "        binomIndiv                     %d\n"
                 "        binomHapcount                  %d\n\n",
                 nPoints, nProductionBad, nUniformBad, nIndivBad, nHapBad);

    CHECK_MSG(nUniformBad == 0 && nIndivBad == 0 && nHapBad == 0,
              "a real kernel produced a non-positive or non-finite K''");
    // Pins the motivation.  Fails if cgfDiffForm is ever "repaired" in place,
    // which would silently destroy the before/after comparison.
    CHECK_MSG(nProductionBad > 0,
              "expected the frozen differencing form to still fail on this grid");
}

// D1 quantifies the differencing form's failure as a 36% relative error at
// t*r = 35.  It is considerably worse than that further out, in two escalating
// steps that the audit does not record, and both are visible on the same p = 0.3
// ray:
//
//   * by t*r = 100 the K'' it returns is NEGATIVE (-2.66e-17).  That is a sign
//     error, not a loss of precision: the cancellation has overwhelmed the
//     quantity entirely.  A negative curvature sends the Newton step the wrong
//     way and trips every `K2 <= 0` guard that exists.
//
//   * by t*r = 356 it is NaN, because the form materializes MGF0 = alpha^2 and
//     alpha^2 overflows once alpha > sqrt(DBL_MAX) = 1.34e154, i.e. at
//     t*r > log(1.34e154 / p).  MGF1 and MGF2 overflow with it, and K'' becomes
//     inf/inf - inf/inf.
//
// The stable form never squares alpha — it multiplies by inva twice, and both pi
// and 1-pi are confined to [0, 1] — so it stays finite and strictly positive
// across the whole range.
TEST(production_form_goes_negative_then_nan_in_the_far_tail) {
    const double p = 0.3, r = 1.0;
    const double rv[1] = {r};

    int nNeg = 0, nNaN = 0;
    for (double tr : {35.0, 100.0, 300.0, 356.0, 400.0, 800.0}) {
        const Cgf old = sparef::cgfDiffForm(tr, r, p);
        const spa_cgf::Cgf12 got = spa_cgf::binomUniformK12(tr, rv, 1, p);

        if (std::isnan(old.K2)) ++nNaN;
        else if (old.K2 < 0.0) ++nNeg;

        // The real kernel must be finite and strictly positive throughout.
        CHECK_MSG(std::isfinite(got.K2) && got.K2 > 0.0,
                  "stable K'' must stay finite and positive in the far tail");
        CHECK_MSG(std::isfinite(got.K1),
                  "stable K' must stay finite in the far tail");
    }

    std::fprintf(stderr,
                 "      far tail (p 0.3, r 1, t*r in [35, 800]): production form gave "
                 "%d negative and %d NaN K''\n", nNeg, nNaN);

    // Characterization: pins both escalation steps so a future "repair" of
    // cgfDiffForm cannot silently erase the evidence.
    CHECK_MSG(nNeg >= 1, "expected the differencing form to return a negative K''");
    CHECK_MSG(nNaN >= 1, "expected the differencing form to return a NaN K''");
}

// ──────────────────────────────────────────────────────────────────────
// II.f  The two accuracy tables, before and after
// ──────────────────────────────────────────────────────────────────────

// D1: K'' as t*r grows.  This is the regime the solver drives markers in the
// extreme tail into, i.e. precisely the markers a scan is run to find.
TEST(real_kernel_K2_accuracy_across_large_tr) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    const double p = 0.5;   // worst case: alpha's two terms stay comparable
    const double r = 1.0;

    struct Row { double tr, before, afterScalar, afterDispatch; };
    std::vector<Row> table;

    for (double tr : {5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0}) {
        const CgfL ref = sparef::cgfReference(tr, r, p);
        const double refK2 = static_cast<double>(ref.K2);

        const double rv[1] = {r};
        const spa_cgf::Cgf12 sc = spa_cgf::tier::binomUniformK12_scalar(tr, rv, 1, p);
        const spa_cgf::Cgf12 dp = spa_cgf::binomUniformK12(tr, rv, 1, p);
        const Cgf old = sparef::cgfDiffForm(tr, r, p);

        table.push_back(Row{tr,
                            tinytest::relDiff(old.K2, refK2),
                            tinytest::relDiff(sc.K2, refK2),
                            tinytest::relDiff(dp.K2, refK2)});

        // The scalar tier is held to the algebra budget: this is the assertion
        // that the cancellation-free identity is correct.  The dispatched tier
        // additionally carries the vendored vector exp's error.
        CHECK_REL(sc.K2, refK2, relTol(kAlgebraUlp));
        CHECK_REL(dp.K2, refK2, relTol(kVectorUlp));
        CHECK_MSG(sc.K2 > 0.0 && dp.K2 > 0.0,
                  "real kernel K'' went non-positive at large t*r");
    }

    std::fprintf(stderr,
                 "\n      K'' relative error vs long-double reference (MAF 0.5, r 1)\n");
    std::fprintf(stderr, "      %8s  %14s  %14s  %14s\n",
                 "t*r", "before (prod)", "after (scalar)", "after (dispatch)");
    for (const Row &row : table)
        std::fprintf(stderr, "      %8.1f  %14.3e  %14.3e  %14.3e\n",
                     row.tr, row.before, row.afterScalar, row.afterDispatch);
    std::fprintf(stderr,
                 "      (the dispatched column's excess over the scalar column is the\n"
                 "       vendored avx512_exp_pd, measured separately at up to 38.9 ULP)\n\n");
}

// N1: K near t = 0.  Every form that materializes alpha = 1 + delta before
// taking the logarithm loses relative precision here, and K enters
// w = sgn(zeta)*sqrt(2*(zeta*s - K)) directly.
//
// The strict relative assertion is made on `binomUniformKFullExact`, the entry
// point that promises it.  The production entry point's relative error is
// tabulated alongside for the record — it is intermediate between the frozen
// production form and the exact one, and Stage 4's rework accepts it because the
// absolute error, which is what w consumes, is unchanged; see
// `production_K_absolute_error_does_not_degrade_w`.
TEST(real_kernel_K_accuracy_near_zero) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    const double p = 0.3;
    const double r = 1.0;

    struct Row { double tr, before, exact, prod, absProd; };
    std::vector<Row> table;

    for (double tr : {1e-8, 1e-6, 1e-4, 1e-2, 1.0}) {
        const CgfL ref = sparef::cgfReference(tr, r, p);
        const double refK0 = static_cast<double>(ref.K0);

        const double rv[1] = {r};
        const spa_cgf::Cgf012 exact = spa_cgf::binomUniformKFullExact(tr, rv, 1, p);
        const spa_cgf::Cgf012 prod  = spa_cgf::binomUniformKFull(tr, rv, 1, p);
        // The shipped pre-rework production K is log(MGF0) = log(alpha^2), per
        // math_helper.hpp's kG012; cgfDiffForm is its frozen copy.
        const Cgf old = sparef::cgfDiffForm(tr, r, p);

        table.push_back(Row{tr,
                            tinytest::relDiff(old.K0, refK0),
                            tinytest::relDiff(exact.K0, refK0),
                            tinytest::relDiff(prod.K0, refK0),
                            std::fabs(prod.K0 - refK0)});

        CHECK_REL(exact.K0, refK0, relTol(kAlgebraUlp));
        CHECK_MSG(k0ProdUlp(prod.K0, refK0, k0Scale1(tr, r, p, 2.0),
                            "near-zero grid") <= kK0ProdUlp,
                  "production K off the absolute scale");
    }

    std::fprintf(stderr,
                 "\n      K error vs long-double reference (MAF 0.3, r 1)\n");
    std::fprintf(stderr, "      %10s  %14s  %14s  %14s  %14s\n", "t*r",
                 "rel: frozen", "rel: exact", "rel: prod", "abs: prod");
    for (const Row &row : table)
        std::fprintf(stderr, "      %10.0e  %14.3e  %14.3e  %14.3e  %14.3e\n",
                     row.tr, row.before, row.exact, row.prod, row.absProd);
    std::fprintf(stderr, "\n");

    CHECK_MSG(table[0].before > 1e-10,
              "expected the frozen production form to have lost >=6 digits at t*r=1e-8");
    // The exact entry point must be strictly better than the production one
    // somewhere, or it is not earning its place in the interface.
    CHECK_MSG(table[0].exact < table[0].prod,
              "kFullExact is no more accurate in K than kFull at t*r = 1e-8");
}

// The log1p form must hold up on the positive side of zero too.  This is why
// logAlpha splits at 700 rather than at 0: the overflow-safe reflected spelling
// would evaluate x + log(...) as a cancelling difference for small positive x
// and lose about four digits there.
TEST(real_kernel_K_accuracy_is_symmetric_about_zero) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }

    double worst = 0.0, worstProd = 0.0, worstProdRel = 0.0;
    for (double p : kMafs) {
        for (double mag : {1e-9, 1e-7, 1e-5, 1e-3, 1e-1, 1.0, 10.0}) {
            for (double sgn : {-1.0, 1.0}) {
                const double tr = sgn * mag;
                const double rv[1] = {1.0};
                const double refK0 =
                    static_cast<double>(sparef::cgfReference(tr, 1.0, p).K0);
                const spa_cgf::Cgf012 got =
                    spa_cgf::binomUniformKFullExact(tr, rv, 1, p);
                const double e = tinytest::relDiff(got.K0, refK0);
                worst = std::fmax(worst, e);
                CHECK_MSG(e <= relTol(kAlgebraUlp),
                          "K left rounding level on one side of zero");

                // The production spelling has no reflection branch at all — it
                // is log(u + p*e^x) throughout — so its absolute error must be
                // symmetric in the sign of x too.
                const spa_cgf::Cgf012 pr =
                    spa_cgf::binomUniformKFull(tr, rv, 1, p);
                const double u = k0ProdUlp(pr.K0, refK0, k0Scale1(tr, 1.0, p, 2.0),
                                           "symmetric grid");
                worstProd = std::fmax(worstProd, u);
                CHECK_MSG(u <= kK0ProdUlp,
                          "production K off the absolute scale on one side of zero");
                // Recorded, not asserted: the shipped spelling's RELATIVE error
                // near t = 0 is the quantity N1 objects to, and printing it
                // beside the absolute one is what keeps the two claims apart.
                if (refK0 != 0.0)
                    worstProdRel =
                        std::fmax(worstProdRel, tinytest::relDiff(pr.K0, refK0));
            }
        }
    }
    std::fprintf(stderr,
                 "      worst K over both signs of t*r:\n"
                 "        kFullExact  relative %.3e   (budget %.3e)\n"
                 "        kFull       ABSOLUTE %.1f ULP of the k0Scale "
                 "(budget %.0f)\n"
                 "        kFull       relative %.3e   (recorded, not a criterion)\n",
                 worst, relTol(kAlgebraUlp), worstProd, kK0ProdUlp, worstProdRel);
}

// ──────────────────────────────────────────────────────────────────────
// II.g  Degenerate inputs
// ──────────────────────────────────────────────────────────────────────

// D3: `kG012Local` sets K0 = -infinity whenever its `base <= 1e-15`, which makes
// zeta*s - K0 = +inf, passes the `temp <= 0` guard, and yields w = inf and a
// p-value of numeric_limits<double>::min() — a fabricated genome-wide
// significant hit.  `base` is (1-q) + q*e^t, so that branch is reached exactly
// when q is at or adjacent to 1 and t*r is far negative.  These are the inputs;
// the kernel must be finite and equal to the analytic law on every one.
TEST(hapcount_degenerate_cases) {
    const std::vector<double> r = {-9.0, -2.0, -0.5, 0.5, 2.0, 9.0, 3.25, -7.75};
    const int n = static_cast<int>(r.size());

    // (1) q = 1: G_i == hap_i almost surely, so r_i*G_i is deterministic.
    //     K = t*sum(h_i r_i), K' = sum(h_i r_i), K'' = 0 exactly.
    {
        const std::vector<double> h = {0.0, 1.0, 2.0, 2.0, 1.0, 0.0, 2.0, 1.0};
        double sumHR = 0.0;
        for (int i = 0; i < n; ++i) sumHR += h[i] * r[i];

        // Includes t = -400: with |t*r| up to 3600 the exponential underflows to
        // exactly zero, which is the corner that makes alpha vanish.
        for (double t : {-400.0, -5.0, -1.0, 0.0, 1.0, 5.0, 400.0}) {
            const spa_cgf::Cgf012 f =
                spa_cgf::binomHapcountKFull(t, r.data(), h.data(), n, 1.0);
            CHECK_MSG(std::isfinite(f.K0) && std::isfinite(f.K1) && std::isfinite(f.K2),
                      "q = 1 must give finite cumulants, not -inf or NaN");
            CHECK_REL(f.K1, sumHR, 1e-15);
            CHECK_MSG(f.K2 == 0.0, "q = 1 is deterministic: K'' must be exactly 0");
            CHECK_REL(f.K0, t * sumHR, 1e-15);
        }
    }

    // (2) q = 0: G_i == 0 almost surely, so every cumulant vanishes exactly.
    {
        const std::vector<double> h = {0.0, 1.0, 2.0, 2.0, 1.0, 0.0, 2.0, 1.0};
        for (double t : {-400.0, -1.0, 0.0, 1.0, 400.0}) {
            const spa_cgf::Cgf012 f =
                spa_cgf::binomHapcountKFull(t, r.data(), h.data(), n, 0.0);
            CHECK_MSG(f.K0 == 0.0 && f.K1 == 0.0 && f.K2 == 0.0,
                      "q = 0 must give exactly zero cumulants");
        }
    }

    // (3) h = 0 for every subject: G_i == 0, all cumulants exactly zero, for any
    //     q.  Note this must hold without a special case in the kernel — every
    //     cumulant carries h as a factor, and no factor it multiplies is
    //     infinite, so the product is exactly zero rather than NaN.
    {
        const std::vector<double> h(r.size(), 0.0);
        for (double q : {0.0, 1e-6, 0.5, 1.0 - 1e-12, 1.0}) {
            for (double t : {-400.0, -1.0, 0.0, 1.0, 400.0}) {
                const spa_cgf::Cgf012 f =
                    spa_cgf::binomHapcountKFull(t, r.data(), h.data(), n, q);
                CHECK_MSG(f.K0 == 0.0 && f.K1 == 0.0 && f.K2 == 0.0,
                          "h = 0 everywhere must give exactly zero cumulants");
            }
        }
    }

    // (4) h = 0 for some subjects: those must contribute exactly nothing, so the
    //     reduction must agree with the reduction over the surviving subset.
    //
    //     None of the three is bit-identical, and requiring that would be a wrong
    //     assertion rather than a stronger one.  Each h = 0 lane does contribute
    //     exactly 0.0, but dropping the zero lanes changes n, which changes how
    //     many surviving terms land in the vector body versus the tail and how
    //     _mm512_reduce_add_pd pairs them.  The surviving addends are therefore
    //     summed in a different ORDER, and floating-point addition is not
    //     associative.  The correct criterion is the reduction-order budget.
    //
    //     STAGE 4.  K was previously asserted bit-identical here, and correctly
    //     so at the time: it was accumulated by a sequential scalar loop, into
    //     which inserting an exact zero addend cannot change anything.  Now that
    //     it too is a vector reduction it falls under the same reduction-order
    //     argument as K' and K''.  `binomHapcountKFullExact` is still the
    //     sequential scalar loop, so the bit-identity is asserted on it.
    {
        std::vector<double> hMixed = {0.0, 1.0, 0.0, 2.0, 1.0, 0.0, 2.0, 0.0};
        std::vector<double> rSub, hSub;
        for (int i = 0; i < n; ++i)
            if (hMixed[i] != 0.0) { rSub.push_back(r[i]); hSub.push_back(hMixed[i]); }

        for (double q : {1e-6, 0.2, 0.5, 0.9}) {
            for (double t : {-2.0, -0.3, 0.0, 0.3, 2.0}) {
                const std::vector<double> qv(r.size(), q);
                const spa_cgf::Cgf012 a = spa_cgf::binomHapcountKFull(
                    t, r.data(), hMixed.data(), n, q);
                const spa_cgf::Cgf012 b = spa_cgf::binomHapcountKFull(
                    t, rSub.data(), hSub.data(), static_cast<int>(rSub.size()), q);
                const L1 sc = l1Scale(t, r, qv, hMixed);

                const spa_cgf::Cgf012 ax = spa_cgf::binomHapcountKFullExact(
                    t, r.data(), hMixed.data(), n, q);
                const spa_cgf::Cgf012 bx = spa_cgf::binomHapcountKFullExact(
                    t, rSub.data(), hSub.data(), static_cast<int>(rSub.size()), q);
                CHECK_MSG(ax.K0 == bx.K0,
                          "h = 0 subjects perturbed the sequential K");
                CHECK_MSG(ulpDiff(a.K0, b.K0, k0Scale(t, r, qv, hMixed))
                              <= kCrossTierUlp,
                          "h = 0 subjects moved K beyond the reduction-order budget");
                CHECK_MSG(ulpDiff(a.K1, b.K1, sc.K1) <= kCrossTierUlp,
                          "h = 0 subjects moved K' beyond the reduction-order budget");
                CHECK_MSG(ulpDiff(a.K2, b.K2, sc.K2) <= kCrossTierUlp,
                          "h = 0 subjects moved K'' beyond the reduction-order budget");
            }
        }
    }

    // (5) q just inside 1, where kG012Local's `base <= 1e-15` branch fires while
    //     the law is still non-degenerate.  Nothing may be infinite here.
    {
        const std::vector<double> h = {1.0, 2.0, 2.0, 1.0, 2.0, 1.0, 2.0, 1.0};
        for (double q : {1.0 - 1e-16, 1.0 - 1e-15, 1.0 - 1e-12, 1.0 - 1e-6}) {
            for (double t : {-400.0, -60.0, -5.0, 0.0, 5.0, 400.0}) {
                const spa_cgf::Cgf012 f =
                    spa_cgf::binomHapcountKFull(t, r.data(), h.data(), n, q);
                CHECK_MSG(std::isfinite(f.K0), "K went non-finite near q = 1");
                CHECK_MSG(std::isfinite(f.K1), "K' went non-finite near q = 1");
                CHECK_MSG(std::isfinite(f.K2) && f.K2 >= 0.0,
                          "K'' went non-finite or negative near q = 1");
            }
        }
    }

    std::fprintf(stderr,
                 "      hapcount degenerate contract: h = 0 and q in {0, 1} are exact, "
                 "finite, and free of the -inf K0 branch\n");
}

// The same corner reached through binomIndiv, where the allele frequency is
// per-subject and the guard therefore cannot be hoisted out of the loop.
TEST(indiv_degenerate_allele_frequency_cases) {
    const std::vector<double> r = {-9.0, -2.0, -0.5, 0.5, 2.0, 9.0, 3.25, -7.75};
    const int n = static_cast<int>(r.size());

    // af = 1 for every subject: G_i == 2 almost surely.  K = 2*t*sum(r),
    // K' = 2*sum(r), K'' = 0.  At t = -400 the exponential underflows and alpha
    // vanishes, which is exactly what the per-lane guard is for.
    {
        const std::vector<double> af(r.size(), 1.0);
        double sumR = 0.0;
        for (double x : r) sumR += x;

        for (double t : {-400.0, -5.0, 0.0, 5.0, 400.0}) {
            const spa_cgf::Cgf012 f =
                spa_cgf::binomIndivKFull(t, r.data(), af.data(), n);
            CHECK_MSG(std::isfinite(f.K0) && std::isfinite(f.K1) && std::isfinite(f.K2),
                      "af = 1 must give finite cumulants");
            CHECK_REL(f.K1, 2.0 * sumR, 1e-15);
            CHECK_MSG(f.K2 == 0.0, "af = 1 is deterministic: K'' must be exactly 0");
            CHECK_REL(f.K0, 2.0 * t * sumR, 1e-15);
        }
    }

    // af = 0 for every subject: G_i == 0, all cumulants exactly zero.
    {
        const std::vector<double> af(r.size(), 0.0);
        for (double t : {-400.0, 0.0, 400.0}) {
            const spa_cgf::Cgf012 f =
                spa_cgf::binomIndivKFull(t, r.data(), af.data(), n);
            CHECK_MSG(f.K0 == 0.0 && f.K1 == 0.0 && f.K2 == 0.0,
                      "af = 0 must give exactly zero cumulants");
        }
    }

    // Mixed: interior, exactly 0 and exactly 1 in the same vector, so the guard
    // must act per lane.  Every tier must handle it, and agree.
    {
        const std::vector<double> af = {0.0, 1.0, 0.5, 1e-6, 1.0, 0.0, 0.3, 1.0 - 1e-16};
        const std::vector<double> twos(r.size(), 2.0);
        const auto tiers = indivTiers();

        for (double t : {-400.0, -60.0, -1.0, 0.0, 1.0, 60.0, 400.0}) {
            const spa_cgf::Cgf12 ref =
                spa_cgf::tier::binomIndivK12_scalar(t, r.data(), af.data(), n);
            const L1 sc = l1Scale(t, r, af, twos);
            CHECK_MSG(std::isfinite(ref.K1) && std::isfinite(ref.K2),
                      "mixed af produced a non-finite scalar result");
            CHECK_MSG(ref.K2 >= 0.0, "mixed af produced a negative K''");

            for (const auto &tt : tiers) {
                const spa_cgf::Cgf12 got = tt.fn(t, r.data(), af.data(), n);
                CHECK_MSG(std::isfinite(got.K1) && std::isfinite(got.K2),
                          std::string("mixed af non-finite, tier ") + tt.name);
                CHECK_MSG(ulpDiff(got.K1, ref.K1, sc.K1) <= kCrossTierUlp,
                          std::string("mixed af K' disagreement, tier ") + tt.name);
                CHECK_MSG(ulpDiff(got.K2, ref.K2, sc.K2) <= kCrossTierUlp,
                          std::string("mixed af K'' disagreement, tier ") + tt.name);
            }
        }
    }

    std::fprintf(stderr,
                 "      indiv degenerate contract: af in {0, 1} exact per lane, guard "
                 "verified on every host tier\n");
}

// K(0) = 0, K'(0) = E[rG] = h*p*r, K''(0) = Var(rG) = h*p*(1-p)*r^2.  These are
// the values every consumer folds its non-outlier block around, so an error here
// would mis-centre the whole method.
TEST(real_kernel_cumulants_at_zero_match_moments) {
    for (double p : kMafs) {
        for (double r : kResidsP2) {
            const double rv[1] = {r};
            const double pv[1] = {p};

            const spa_cgf::Cgf012 u = spa_cgf::binomUniformKFull(0.0, rv, 1, p);
            CHECK_MSG(u.K0 == 0.0, "K(0) must be exactly zero");
            CHECK_REL(u.K1, 2.0 * p * r, 1e-15);
            CHECK_REL(u.K2, 2.0 * p * (1.0 - p) * r * r, 1e-15);

            const spa_cgf::Cgf012 i = spa_cgf::binomIndivKFull(0.0, rv, pv, 1);
            CHECK_MSG(i.K0 == 0.0, "K(0) must be exactly zero");
            CHECK_REL(i.K1, 2.0 * p * r, 1e-15);
            CHECK_REL(i.K2, 2.0 * p * (1.0 - p) * r * r, 1e-15);

            for (double h : kHaps) {
                const double hv[1] = {h};
                const spa_cgf::Cgf012 g =
                    spa_cgf::binomHapcountKFull(0.0, rv, hv, 1, p);
                CHECK_MSG(g.K0 == 0.0, "K(0) must be exactly zero");
                CHECK_REL(g.K1, h * p * r, 1e-15);
                CHECK_REL(g.K2, h * p * (1.0 - p) * r * r, 1e-15);
            }
        }
    }
}

// Reports which tiers this host exercised, so a green run on a machine without
// AVX-512 cannot be mistaken for evidence that the AVX-512 kernels were tested.
//
// Also asserts that the dispatch actually SELECTED the widest supported tier,
// rather than merely that such a tier exists and passes when called directly.
// The check is behavioural rather than a function-pointer comparison, because
// the public entry point is a wrapper around the pointer (it resolves the p == 0
// and p == 1 endpoints first).  For an interior p the wrapper forwards verbatim,
// so a bit-identity requirement against the expected tier is exact: two
// different tiers cannot agree bit-for-bit here, since they use different exp
// implementations and different reduction orders.
TEST(report_simd_tier) {
    const char *name = (kLevel >= 2) ? "AVX-512" : ((kLevel >= 1) ? "AVX2" : "scalar");
    std::fprintf(stderr,
                 "      simdLevelValue() = %d (%s); tiers compared: scalar%s%s + dispatch\n",
                 kLevel, name,
                 (kLevel >= 1) ? " + AVX2" : "",
                 (kLevel >= 2) ? " + AVX512" : "");
#if defined(__x86_64__) || defined(_M_X64)
    CHECK_MSG(kLevel >= 1,
              "x86-64 host reported no AVX2: the AVX2 tier went untested");
#endif

    // n = 100 and p = 0.23: interior p, and enough subjects that a full vector
    // body plus a tail runs at both widths.
    const Subjects s = makeSubjects(100);
    const int n = 100;
    const double p = 0.23, q = 0.23, t = 0.37;

    const spa_cgf::Cgf12 du = spa_cgf::binomUniformK12(t, s.r.data(), n, p);
    const spa_cgf::Cgf12 di = spa_cgf::binomIndivK12(t, s.r.data(), s.af.data(), n);
    const spa_cgf::Cgf12 dh =
        spa_cgf::binomHapcountK12(t, s.r.data(), s.hap.data(), n, q);

    spa_cgf::Cgf12 eu, ei, eh;
#if defined(__x86_64__) || defined(_M_X64)
    if (kLevel >= 2) {
        eu = spa_cgf::tier::binomUniformK12_avx512(t, s.r.data(), n, p);
        ei = spa_cgf::tier::binomIndivK12_avx512(t, s.r.data(), s.af.data(), n);
        eh = spa_cgf::tier::binomHapcountK12_avx512(t, s.r.data(), s.hap.data(), n, q);
    } else if (kLevel >= 1) {
        eu = spa_cgf::tier::binomUniformK12_avx2(t, s.r.data(), n, p);
        ei = spa_cgf::tier::binomIndivK12_avx2(t, s.r.data(), s.af.data(), n);
        eh = spa_cgf::tier::binomHapcountK12_avx2(t, s.r.data(), s.hap.data(), n, q);
    } else
#endif
    {
        eu = spa_cgf::tier::binomUniformK12_scalar(t, s.r.data(), n, p);
        ei = spa_cgf::tier::binomIndivK12_scalar(t, s.r.data(), s.af.data(), n);
        eh = spa_cgf::tier::binomHapcountK12_scalar(t, s.r.data(), s.hap.data(), n, q);
    }

    CHECK_MSG(du.K1 == eu.K1 && du.K2 == eu.K2,
              "binomUniform dispatch did not select the widest supported tier");
    CHECK_MSG(di.K1 == ei.K1 && di.K2 == ei.K2,
              "binomIndiv dispatch did not select the widest supported tier");
    CHECK_MSG(dh.K1 == eh.K1 && dh.K2 == eh.K2,
              "binomHapcount dispatch did not select the widest supported tier");

    std::fprintf(stderr,
                 "      dispatch verified bit-identical to the %s tier for all three "
                 "variants\n", name);
}

// ──────────────────────────────────────────────────────────────────────
// III  The budgets, reported rather than asserted in silence
// ──────────────────────────────────────────────────────────────────────
//
// Registration order is definition order within a translation unit, so this
// case runs after every other and sees the file-wide worst.  It exists because
// kK0ProdUlp sat at 512 against a worst observed 3.5 — a 146x headroom in a
// file whose stated convention is 2x to 4.6x — and nothing in `make test`
// output made that visible.  The budget has since been brought to 64; this
// test is what will make the next drift visible, in BOTH directions:
// a headroom that grows past 40x means the budget has stopped bounding
// anything, and one that falls below 2x means it is about to start failing on
// an unrelated host.
TEST(terminal_K_budget_headroom) {
    if (!sparef::longDoubleIsWider()) {
        std::fprintf(stderr, "      SKIP: no wider reference type on this platform.\n");
        return;
    }
    const double headroom =
        (gWorstK0Prod > 0.0) ? kK0ProdUlp / gWorstK0Prod
                             : std::numeric_limits<double>::infinity();
    std::fprintf(stderr,
                 "      kK0ProdUlp %.0f, worst observed %.2f ULP at [%s], "
                 "headroom %.1fx\n",
                 kK0ProdUlp, gWorstK0Prod, gWorstK0ProdWhere, headroom);
    CHECK_MSG(gWorstK0Prod > 0.0,
              "no production-K assertion ran; the accumulator is not wired up");
    CHECK_MSG(gWorstK0Prod <= kK0ProdUlp, "budget exceeded");
    CHECK_MSG(headroom <= 40.0,
              "kK0ProdUlp headroom is " + std::to_string(headroom) +
                  "x; the budget has stopped bounding anything, tighten it");
}

TINYTEST_MAIN
