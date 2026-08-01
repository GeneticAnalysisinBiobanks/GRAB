// conditional_p.hpp — WtCoxG / LEAF's two-component conditional p-value, and
// the rule that decides whether a mixture leg which could not be computed
// still leaves the marker with an answer.
//
// Header-only.  Depends on util/spa.hpp for spa::Status alone; it performs no
// saddlepoint work of its own, which is why tests/wtcoxg_cgf_test.cpp can
// exercise it without linking any method object.

#pragma once

#include "util/spa.hpp"

#include <cmath>
#include <limits>

namespace wtcoxg_cond {

// ── The conditional p-value of the two batch-effect branches ─────────
//
//     p_con = 2 * (TPR*p1 + (1-TPR)*p0) / p_deno
//
// where p0 and p1 are joint probabilities P(S' <= -|s|, S_bat in B) under
// sigma2 = 0 and sigma2 > 0 and p_deno is the matching mixture of the S_bat
// marginals P(S_bat in B).  Every numerator term is bounded by its
// denominator term, and the conditioning event B is symmetric under
// (S', S_bat) -> (-S', -S_bat) while the joint law is centred, so the
// conditional tail is at most one half and p_con is at most 1 — in exact
// arithmetic.
//
// It was NOT at most 1 in practice: 166 of the 3000 markers of the bundled
// fixture reported P_EXT above 1, up to 2.98664, which is not a probability
// and is the only column anywhere in examples_output outside [0, 1].  The
// cause was not this expression but math::bvnCdf, whose |rho| >= 0.925 branch
// was mis-transcribed from Genz (2004) and returned a joint probability
// larger than the marginal that bounds it; see the comment on its definition
// in util/math_helper.cpp.  With that repaired the excess is gone.
//
// The clamp below therefore guards rounding, not a defect: 20-point
// Gauss-Legendre on the rectangle and the asymptotic expansion at |rho| = 1
// are each good to ~1e-13, and a ratio of two such quantities can still land a
// few ULP above 1 when the true value is 1.  It is written as an explicit
// two-sided clamp rather than std::min(1.0, p_con) precisely because
// std::min(1.0, NaN) returns 1.0 — the D2 idiom the saddlepoint unification
// removed — so a NaN must be tested for, not minimised against.
//
// ── When one leg of the mixture cannot be computed ───────────────────
//
// Divide the ratio through by p_deno and the reported quantity is
//
//     p_con = 2 * (omega * q1 + (1 - omega) * q0),
//     omega = TPR*m1 / p_deno,   q_k = p_k / m_k,
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
// The threshold is a RELATIVE width: the interval must be narrower than
// kImmaterialLeg times the value being reported.  `numToChars` formats every
// result cell through plink2's dtoa_g, which emits six significant digits, so
// a relative perturbation of 1e-8 is a hundredth of a unit in the last printed
// digit of P and moves -log10(P) by at most 4.3e-9.
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

struct ConditionalP {
    double p;
    double negLog10p;
    spa::Status status;
};

// The status to attribute to a leg that produced no usable joint probability.
//
// If the leg's own saddlepoint already left it with no answer, keep that
// reason.  Otherwise the leg was lost DOWNSTREAM of the saddlepoint:
// math::pmvnorm2dHalfRect returns NaN when the (var_S, cov, var_Sbat) triple
// it is handed is not a covariance matrix, which is what happens when sigma2
// drives |rho| past 1.  That is NA_POST_FAIL, and separating it from the
// saddlepoint's own failures is the point of the re-partition: on the bundled
// fixture all 177 such rows have |Z_Norm_EXT| far below the SPA cutoff, so the
// saddlepoint was never attempted and the normal tail at that Z would answer a
// question nobody asked.  Those rows stay NA (decision D5).
inline spa::Status legFailureStatus(spa::Status st) {
    return spa::statusIsNA(st) ? st : spa::Status::NaPostFail;
}

// w1/m1/p1 describe the sigma2 > 0 component, w0/m0/p0 the sigma2 = 0 one:
// w is the mixture weight (TPR and 1 - TPR), m the conditioning mass
// P(S_bat in B | component), p the joint probability, so that
// p_deno = w1*m1 + w0*m0 and the numerator is w1*p1 + w0*p0.
inline ConditionalP conditionalP(
    double w1,
    double m1,
    double p1,
    spa::Status st1,
    double w0,
    double m0,
    double p0,
    spa::Status st0,
    double p_deno
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const bool ok1 = std::isfinite(p1);
    const bool ok0 = std::isfinite(p0);

    double numer;             // every available leg's contribution
    double missingHalfMass;   // (1/2) * the unavailable leg's term in p_deno
    spa::Status status;       // status of the leg(s) that produced `numer`
    spa::Status lost;         // status to report if the marker has no answer
    if (ok1 && ok0) {
        numer = w1 * p1 + w0 * p0;
        missingHalfMass = 0.0;
        status = spa::worseStatus(st0, st1);
        lost = status;
    } else if (ok0) {
        numer = w0 * p0;
        missingHalfMass = 0.5 * w1 * m1;
        status = st0;
        lost = legFailureStatus(st1);
    } else if (ok1) {
        numer = w1 * p1;
        missingHalfMass = 0.5 * w0 * m0;
        status = st1;
        lost = legFailureStatus(st0);
    } else {
        return {nan, nan, legFailureStatus(spa::worseStatus(st0, st1))};
    }

    // An unusable denominator is a failure of the conditional assembly, not
    // of any saddlepoint: NA_POST_FAIL, never a normal tail.
    double p = 2.0 * numer / p_deno;
    if (!std::isfinite(p)) return {nan, nan, spa::Status::NaPostFail};

    // The immaterial-leg rule.  `width` is the whole interval the missing leg
    // spans.  Both comparisons are written negated so that a NaN takes the
    // conservative branch: a leg whose mass is not even a number is never
    // certified immaterial.
    if (!(missingHalfMass <= 0.0)) {
        const double width = 2.0 * missingHalfMass / p_deno;
        if (!(width <= kImmaterialLeg * p)) return {nan, nan, lost};
    }

    if (p > 1.0) p = 1.0;
    if (p < 0.0) p = 0.0;
    double neg = -std::log10(p);
    if (neg == 0.0) neg = 0.0;   // -log10(1) is -0.0; normalise the sign
    return {p, neg, status};
}

}  // namespace wtcoxg_cond
