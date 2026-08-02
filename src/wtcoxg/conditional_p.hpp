// conditional_p.hpp — WtCoxG / LEAF's two-component conditional p-value, and
// the rule that decides whether a mixture leg which could not be computed
// still leaves the marker with an answer.
//
// Header-only.  Depends on util/spa.hpp for spa::Status and on
// util/math_helper.hpp for the two log-domain combiners; it performs no
// saddlepoint work of its own, which is why tests/wtcoxg_cgf_test.cpp can
// exercise it without linking any method object.

#pragma once

#include "util/math_helper.hpp"
#include "util/spa.hpp"

#include <cmath>
#include <limits>

namespace wtcoxg_cond {

// ── The conditional p-value of the two batch-effect branches ─────────
//
//     p_con = 2 * (w1*p1 + w0*p0) / p_deno,     p_deno = w1*m1 + w0*m0
//
// where p0 and p1 are joint probabilities P(S' <= -|s|, S_bat in B) under
// sigma2 = 0 and sigma2 > 0, w1 = TPR and w0 = 1 - TPR are the mixture
// weights, and m0, m1 are the matching S_bat marginals P(S_bat in B).  Every
// numerator term is bounded by its denominator term, and the conditioning
// event B is symmetric under (S', S_bat) -> (-S', -S_bat) while the joint law
// is centred, so the conditional tail is at most one half and p_con is at most
// 1 — in exact arithmetic.
//
// ── Every argument is a natural logarithm ────────────────────────────
//
// log10p_unify Stage 6.  The expression above is evaluated as
//
//     ln p_con = ln2 + logaddexp(ln w1 + ln p1, ln w0 + ln p0) - ln p_deno
//
// with no exponential taken anywhere between the bivariate integrals and the
// reported -log10(P).  It has to be: the ratio underflows.  On the 50 000 x
// 2 000 witness cohort of dev-notes/methods/log10p_unify (WtCoxG,
// --prevalence 0.3804) seven markers had a conditional p below DBL_MIN, and
// the linear assembly reported each of them as LOG10P_EXT = +Inf with
// SPA_STATUS = 0 SPA_OK — a saddlepoint failure that never happened, attached
// to a p-value that is not a number, on the six markers of the fixture with
// the strongest association.  Nothing in the -log10 of a linear ratio can
// distinguish that from a converged result, which is why the ratio is now
// never formed.
//
// The inputs arrive in the log domain from where they are already computed
// there: ln p0 and ln p1 from math::pmvnorm2dHalfRectLog, and ln m0 / ln m1
// from the two math::pnormLog calls the caller makes to build p_deno.  ln w1
// and ln w0 are logs of TPR and 1 - TPR.
//
// A note on the clamp at the end.  It guards rounding, not a defect: the
// adaptive quadrature and the asymptotic expansion at |rho| = 1 are each good
// to ~1e-13, and a ratio of two such quantities can still land a few ULP above
// 1 when the true value is 1.  (It was NOT at most 1 in a much older tree: 166
// of the 3000 markers of the bundled fixture reported P_EXT above 1, up to
// 2.98664, because math::bvnCdf's |rho| >= 0.925 branch was mis-transcribed
// from Genz (2004) and returned a joint probability larger than the marginal
// that bounds it.  With that repaired the excess is gone.)  It is written as
// an explicit comparison rather than std::min(1.0, .) precisely because
// std::min(1.0, NaN) returns 1.0 — the D2 idiom the saddlepoint unification
// removed — so a NaN must be tested for, not minimised against.
//
// ── When one leg of the mixture cannot be computed ───────────────────
//
// Divide the ratio through by p_deno and the reported quantity is
//
//     p_con = 2 * (omega * q1 + (1 - omega) * q0),
//     omega = w1*m1 / p_deno,   q_k = p_k / m_k,
//
// i.e. twice a posterior-weighted average of the two components' CONDITIONAL
// tail probabilities, weighted by omega = P(batch effect | S_bat in B), the
// posterior probability that the batch-effect component is the operative one
// given the outcome of the batch test.  Two facts about q_k follow from the
// definitions alone, with no distributional assumption beyond the one the
// formula already makes:
//
//   * q_k >= 0, because p_k is a probability; and
//   * q_k <= 1/2, because p_k is the joint probability of {S' <= -|s|}
//     intersected with the very event m_k measures, and because B is
//     symmetric under (S', S_bat) -> (-S', -S_bat) while the joint law is
//     centred, so the conditional law of S' given S_bat in B is symmetric
//     about zero and its tail at -|s| <= 0 is at most 1/2.
//
// So a leg that cannot be computed does not leave the answer unknown: it
// leaves it confined to an interval
//
//     p_con in [ 2*(1-omega)*q0 , 2*(1-omega)*q0 + omega ]
//
// of width exactly omega = (that leg's own term in p_deno) / p_deno.  The
// interval, not the leg's mere existence, is what decides whether the marker
// still has an answer:
//
//   * If the interval is narrower than the precision at which the answer is
//     reported, no admissible value of the missing leg can change a single
//     printed digit.  Discarding the marker then destroys a usable p-value in
//     exchange for nothing.  Report the interval's lower end — the mixture
//     with the missing leg's joint probability set to zero — and the SURVIVING
//     leg's status, which is the status of the saddlepoint that actually
//     produced the reported number.
//   * Otherwise the missing leg does move the answer, decision L2 applies
//     unchanged, and the marker is NA with a status.
//
// WHAT COUNTS AS A LEG THAT CANNOT BE COMPUTED changed with the log domain,
// and it is the one semantic change of Stage 6 that is not a change of
// representation.  The test is that ln p_k be FINITE.  NaN means the same
// thing it did (math::pmvnorm2dHalfRectLog was handed a triple that is not a
// covariance matrix), but -infinity is new: it is what the degenerate
// |rho| >= 1 - 1e-12 branch returns when the linear expansion it still uses
// underflows to zero.  Such a leg is treated as MISSING rather than as a leg
// whose probability is zero, because that is what it is — the integral was not
// evaluated, only bounded below by the smallest positive double.  Calling it
// zero would understate p_con by an unbounded factor, i.e. overstate the
// significance, which is the direction an error must never be allowed to take.
// The interval rule then decides, from a bound that holds regardless
// (q_k <= 1/2), whether the marker keeps an answer.  The linear predecessor
// treated an underflowed leg as exactly zero and said nothing.
//
// The threshold is a RELATIVE width: the interval must be narrower than
// kImmaterialLeg times the value being reported.  `numToChars` formats every
// result cell through plink2's dtoa_g, which emits six significant digits, so
// a relative perturbation of 1e-8 is a hundredth of a unit in the last printed
// digit of P and moves -log10(P) by at most 4.3e-9.  In the log domain the
// same test is a difference of logarithms, ln width - ln p_con <= ln(1e-8);
// the threshold's VALUE is unchanged, because the measurement it rests on
// (below) is a statement about a ratio and a ratio does not know which
// representation it is held in.
//
// The bound is stated relative to the reported value, not as an absolute
// probability, precisely so that the rule cannot quietly flatten a genuinely
// small p-value: a marker deep enough in the tail for the missing leg to
// matter against its own p-value is reported NA, not silently raised to the
// leg's mass.
//
// The constant is not delicate.  Measured over 800 000 markers of the 5000 x
// 1e6 synthetic null cohort with an external reference (LEAF, two clusters,
// two phenotypes, --ref-af), the relative widths separate into two
// populations with eight empty decades between them: 777 903 in
// [1.2e-37, 3.2e-12] and 22 097 in [4.3e-4, 7.5].  Any threshold strictly
// inside (3.2e-12, 4.3e-4) partitions the markers identically; 1e-8 is near
// the middle of that gap on the log scale.
//
// This subsumes, and is strictly weaker than, a threshold on the mixture
// weight itself.  On that cohort TPR reaches 1.9e-26, but m1 — the sigma2 > 0
// component's conditioning mass — is a further 1e-11 smaller, because
// sigma2 = 4.0e+16 spreads S_bat over a range 1e10 times wider than the
// acceptance interval [lb, ub].  Both factors belong in the bound and both are
// available without any saddlepoint: m1 is exactly the quantity Branch B's two
// pnormLog calls already compute for p_deno.
constexpr double kImmaterialLeg = 1e-8;

// -log10 of the conditional p-value, and the status that says what produced
// it.  There is no linear `p` field: it was filled here from the same ratio
// that underflowed, and a caller that still emits a P column derives it once,
// through spa::pFromNegLog10P, exactly as every other method does since
// log10p_unify Stage 3.
struct ConditionalP {
    double negLog10p;
    spa::Status status;
};

// The status to attribute to a leg that produced no usable joint probability.
//
// If the leg's own saddlepoint already left it with no answer, keep that
// reason.  Otherwise the leg was lost DOWNSTREAM of the saddlepoint:
// math::pmvnorm2dHalfRectLog returns NaN when the (var_S, cov, var_Sbat)
// triple it is handed is not a covariance matrix, which is what happens when
// sigma2 drives |rho| past 1, and -infinity when the degenerate branch
// underflows.  Both are NA_POST_FAIL, and separating them from the
// saddlepoint's own failures is the point of the re-partition: on the bundled
// fixture all 177 such rows have |Z_Norm_EXT| far below the SPA cutoff, so the
// saddlepoint was never attempted and the normal tail at that Z would answer a
// question nobody asked.  Those rows stay NA (decision D5).
inline spa::Status legFailureStatus(spa::Status st) {
    return spa::statusIsNA(st) ? st : spa::Status::NaPostFail;
}

// lnW1/lnM1/lnP1 describe the sigma2 > 0 component, lnW0/lnM0/lnP0 the
// sigma2 = 0 one; all six are natural logarithms of the quantities named in
// the header comment, and lnDeno is ln(w1*m1 + w0*m0).
//
// ON THE STATUS OF A LEG WHOSE SADDLEPOINT FELL BACK.  Statuses 3..6 mean the
// leg's own p-value is the substituted normal tail (decision D5).  Such a
// status is propagated to the marker unchanged, and the reader should know
// what it then means on THIS column: not that LOG10P_EXT is itself a normal
// tail — it is the conditional probability, assembled from the bivariate
// integrals as always — but that one of the variances those integrals were
// built from was recovered by inverting a p-value that the saddlepoint did not
// produce.  The D4 ordering still says exactly the right thing
// (SPA_STATUS <= 2 selects the rows with no substitution anywhere in them,
// 3..6 the degraded ones, >= 7 the absent ones), and C3 is untouched since
// 3..6 requires only that LOG10P be a finite number.  Demoting such a row to
// NA instead would contradict D5, which reserves NA for the case where there
// is nothing to fall back to; inventing a tenth status code would contradict
// D4, which is locked.
inline ConditionalP conditionalP(
    double lnW1,
    double lnM1,
    double lnP1,
    spa::Status st1,
    double lnW0,
    double lnM0,
    double lnP0,
    spa::Status st0,
    double lnDeno
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    const bool ok1 = std::isfinite(lnP1);
    const bool ok0 = std::isfinite(lnP0);

    double lnNumer;        // ln of every available leg's contribution
    double lnMissMass;     // ln of the unavailable leg's term in p_deno
    spa::Status status;    // status of the leg(s) that produced `lnNumer`
    spa::Status lost;      // status to report if the marker has no answer
    if (ok1 && ok0) {
        lnNumer = math::logAddExp(lnW1 + lnP1, lnW0 + lnP0);
        lnMissMass = -inf;
        status = spa::worseStatus(st0, st1);
        lost = status;
    } else if (ok0) {
        lnNumer = lnW0 + lnP0;
        lnMissMass = lnW1 + lnM1;
        status = st0;
        lost = legFailureStatus(st1);
    } else if (ok1) {
        lnNumer = lnW1 + lnP1;
        lnMissMass = lnW0 + lnM0;
        status = st1;
        lost = legFailureStatus(st0);
    } else {
        return {nan, legFailureStatus(spa::worseStatus(st0, st1))};
    }

    // An unusable denominator is a failure of the conditional assembly, not of
    // any saddlepoint: NA_POST_FAIL, never a normal tail.  A numerator that is
    // -infinity lands here too — both mixture weights zero — and so would a
    // conditional probability small enough to underflow, if the log domain
    // could produce one.  What cannot come out of this line is +Inf.
    const double lnP = math::kLn2 + lnNumer - lnDeno;
    if (!std::isfinite(lnP)) return {nan, spa::Status::NaPostFail};

    // The immaterial-leg rule.  `lnWidth` is the whole interval the missing
    // leg spans, relative to p_deno; with both legs present it is -infinity
    // and the test passes without consulting anything.  The comparison is
    // written negated so that a NaN takes the conservative branch: a leg whose
    // mass is not even a number is never certified immaterial.
    const double lnWidth = lnMissMass - lnDeno;
    if (!(lnWidth <= std::log(kImmaterialLeg) + lnP)) return {nan, lost};

    double neg = -std::fmin(lnP, 0.0) / math::kLn10;   // p clamped at 1
    if (neg == 0.0) neg = 0.0;   // -log10(1) is -0.0; normalise the sign
    return {neg, status};
}

}  // namespace wtcoxg_cond
