// linear_p.hpp — the linear p-value, for tests only.
//
// `spa::pFromNegLog10P` was deleted from the shipped tree by log10p_unify
// Stage 8, together with the `P` column it existed to produce: decision D1
// makes `LOG10P` = −log10(P) the sole p-value representation, and nothing in
// `src/` re-forms a linear p from it any more.
//
// The tests still need the conversion, for a reason the production code does
// not have: a large part of the suite pins the magnitude by comparing it
// against an independently computed LINEAR reference — Boost's `pnorm`, an
// analytic special case, a `long double` evaluation — and the comparison has
// to happen in one representation or the other.  Doing it in the linear
// representation is the correct choice wherever the reference itself is a
// linear quantity that does not underflow, because it leaves the reference
// untouched; where the reference does underflow, the surrounding test compares
// magnitudes directly and does not call this at all.
//
// One definition, in one header, for the same reason the shipped helper had
// one: a conversion re-spelled in six test files would drift between them and
// the drift would be invisible, since each file would still pass against its
// own spelling.
//
// The body is `spa::pFromNegLog10P`'s, verbatim, so a test that pinned the
// round trip before Stage 8 pins exactly the same numbers after it:
//
//   * `std::pow(10, -L)` and not `std::exp(-L*ln10)` — the latter rounds the
//     product `L*ln10` first, and an error of 1 ulp in an exponent of that
//     magnitude is a RELATIVE error of the same size in the result, so at
//     L = 300 it would cost ~7e-14 where `pow` costs ~1e-16.
//   * NaN propagates: a p-value that does not exist must not become a number.
//   * L <= 0 returns exactly 1 rather than a `pow()` of a negative zero.

#ifndef GRAB_TESTS_LINEAR_P_HPP
#define GRAB_TESTS_LINEAR_P_HPP

#include <cmath>
#include <limits>

namespace tref {

inline double pFromNegLog10P(double negLog10p) noexcept {
    if (std::isnan(negLog10p)) return std::numeric_limits<double>::quiet_NaN();
    if (!(negLog10p > 0.0)) return 1.0;
    return std::pow(10.0, -negLog10p);
}

}  // namespace tref

#endif  // GRAB_TESTS_LINEAR_P_HPP
