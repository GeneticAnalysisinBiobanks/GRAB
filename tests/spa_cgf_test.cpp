// spa_cgf_test.cpp — Correctness of the diploid binomial CGF evaluation.
//
// Stage 0 of the saddlepoint rework.  This file establishes, before any
// production code is touched:
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

#include <cstdio>
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

TINYTEST_MAIN
