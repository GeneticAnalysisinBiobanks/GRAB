// log10p_test.cpp — The log-domain distribution tier of the log10p_unify
// project.
//
// ══════════════════════════════════════════════════════════════════════
// Stage 0 — skeleton only
// ══════════════════════════════════════════════════════════════════════
//
// This file is created empty of substance on purpose.  Stage 0 of
// dev-notes/methods/log10p_unify/03_stages.md establishes the baseline and the
// tooling and changes no production behavior, so it has nothing to assert yet;
// what it does establish is that the binary is built, run and reported by
// `make test`, so that Stage 1 adds assertions to a suite that is already
// wired rather than discovering the wiring and the mathematics at once.
//
// Stage 1 fills this file with the cross-checks listed in 03_stages.md §Stage 1
// for the five functions it adds to src/util/math_helper.{hpp,cpp}:
//
//   * zFromNegLog10P        — round-trip against normalTwoSidedLog over
//                             z in [0.1, 1e6], plus the absolute pins of
//                             01_numerics §3.1, of which
//                             zFromNegLog10P(300, +1) == 37.0657878807721 is
//                             the project's central acceptance point (the
//                             present Z column saturates at 37.047096);
//   * chisq1FromNegLog10P   — the identity chisq1FromNegLog10P(L) ==
//                             zFromNegLog10P(L, +1)^2;
//   * cauchyCombineLog10    — the eight cases of 01_numerics §2.3, including
//                             {400, 400} -> 400 and {8, NaN} -> 8, and the
//                             accuracy of the atan inversion for lnT slightly
//                             above 1, where the three-term series is wrong by
//                             3.5e-5 and the direct log(atan(exp(-lnT))) form
//                             is not;
//   * ptLog                 — agreement with boost::math::ibeta wherever the
//                             latter does not underflow, over
//                             df in {5, 100, 10000} x |t| in {2, 20, 50, 200},
//                             which is the grid that separates the two 2F1
//                             forms of 01_numerics §3.3;
//   * pmvnorm2dHalfRectLog  — branch (a) against branch (b) in their overlap
//                             region, and both against a long-double Simpson
//                             reference.
//
// Stage 9 adds the HweLnP cross-check of decision D7.
//
// Build:  make test   (or: make build/tests/log10p_test && ./build/tests/log10p_test)

#include "tinytest.hpp"

#include <cmath>

// The convention decision D2 fixes for every p-value column this project
// produces: LOG10P is -log10(P), not log10(P) and not ln(P).  The assertion is
// deliberately trivial — it exists so that the Stage 0 binary asserts the
// meaning of the column it is named after rather than asserting nothing at
// all, and so that `make test` reports a passing case rather than an empty
// suite.
TEST(log10p_is_minus_log10_of_p) {
    CHECK_NEAR(-std::log10(1e-8), 8.0, 1e-15);
    CHECK_NEAR(-std::log10(1.0), 0.0, 0.0);
    // A p-value of 1 is LOG10P = 0, and no valid p-value yields a negative
    // LOG10P.  This is invariant C2 of 04_validation.md §2, checked here on the
    // convention and by tests/regress.py on every emitted column.
    CHECK(-std::log10(1.0) >= 0.0);
}

TINYTEST_MAIN
