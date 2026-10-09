// meta_pvalue.hpp — Fixed-effects inverse-variance meta-analysis of
// per-cluster score statistics.
//
// The combining formula is the score-pooled meta:
//
//     Var_c    = S_c² / chisq1FromNegLog10P(L_c)
//     Z_meta   = Σ S_c / sqrt(Σ Var_c)
//     L_meta   = −log10( 2 · Φ(−|Z_meta|) )      (two-sided)
//
// where L_c = −log10(p_c) is the per-cluster p-value in the sole
// representation the tree carries (log10p_unify decision D1).
//
// The df = 1 chi-squared upper-tail quantile is analytic — P(χ²₁ > q) = 2Φ(−√q),
// so q is the square of the two-sided normal quantile — and is evaluated by
// `math::chisq1FromNegLog10P`, which has no argument clamp.  Until
// log10p_unify Stage 4 this line read `math::qchisq(p_c, 1, lower.tail = FALSE)`
// on a linear p_c clamped from below at META_P_FLOOR = 1e-300, because qchisq
// clamps its own argument there and saturates at q = 1373.87.  Both the clamp
// and the constant are gone: they existed only to keep the saturated branch
// away from the ±∞ its own clamp would otherwise have produced, and the
// analytic inversion has neither the saturation nor the need for the guard.
//
// The remaining bound is at the OTHER end and is not a p-value floor:
//
//   p_c → 1  (L_c → 0) :  q → 0 ⇒ Var_c = S²/0 = +∞ ⇒ the cluster poisons
//              the sum.  L_c is held at or above META_L_FLOOR, keeping q above
//              ≈ 1.57e-30, so such a cluster contributes a vanishingly small
//              weight instead of an infinite variance.  This is a degeneracy
//              of the recovered-variance formula at a NULL result, not a floor
//              on a significant one, and it moves no reported magnitude.
//
// The opposite degeneracy — L_c = +∞, i.e. a cluster p-value that underflowed
// to exactly zero — now yields Var_c = 0 rather than a clamped variance: the
// natural limit of the formula, in which that cluster's score is treated as
// exact.  If it is the only contributing cluster the pool has no variance at
// all and reports NA with NA_POST_FAIL; otherwise it enters the numerator with
// no matching denominator term, and the pooled magnitude is then bounded only
// by the OTHER clusters' variances — measured at L_meta = 3.5e4 on one
// bench_data/witness marker, which is not a number any reader should trust.
//
// That input no longer arises.  Its only producer was WtCoxG's conditional
// p-value taking −log10 of a linear ratio, and log10p_unify Stage 6 rewrote
// `src/wtcoxg/conditional_p.hpp` in the log domain: `wtcoxg_cond::conditionalP`
// assembles the mixture through `math::logAddExp` / `math::logSubExp` over
// `math::pmvnorm2dHalfRectLog`, returns a magnitude directly, and a denominator
// that is not usable comes back as NA with NA_POST_FAIL rather than as an
// infinite L.  The paragraph above is kept because the LIMIT it describes is
// still the behaviour of this function for an infinite L; what changed is that
// no caller in the tree supplies one.  It is still deliberately not
// intercepted here.
//
// Clusters with NaN score or NaN L are skipped (their per-cluster computation
// failed upstream — that is a distinct signal from a finite boundary p-value
// and must not be folded into the bound).

#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "util/math_helper.hpp"
#include "util/spa.hpp"

namespace math {

// The L-scale spelling of the former META_P_CEIL = 1 − 1e-15.  A cluster
// reporting p ≥ 1 − 1e-15 is reporting |z| ≈ 0; holding L here keeps
// chisq1FromNegLog10P above ≈ 1.57e-30 so the recovered variance stays finite
// and the cluster contributes near zero.  Its former partner META_P_FLOOR was
// deleted in log10p_unify Stage 4 (04_validation.md invariant C5): there is no
// lower bound on a p-value anywhere on the output path.
//
// The value is −log10 of the DOUBLE nearest 1 − 1e-15, which is
// 0.999999999999999 and not 1 − 1e-15 exactly: doubles are spaced 1.11e-16
// apart there, so the literal the predecessor clamped to already carried an
// 8e-4 relative error against the real number it was written as.  Taking the
// image of the double rather than the asymptotic 1e-15/ln 10 = 4.342944819e-16
// keeps this bound at exactly the same place it stood before the stage.
constexpr double META_L_FLOOR = 4.339473599489794e-16;

// Pooled result: the magnitude −log10(P) and the status of the POOLING.
//
// The linear pooled p was dropped in log10p_unify Stage 4 together with the
// linear inversion it was paired with, and the meta_P_EXT / meta_P_NOEXT
// columns LEAF derived from it went in Stage 8 (decision D1).
//
// The status is the worst spa::Status among the clusters that actually
// contributed.  Contributing requires a finite score and a finite L.  Under
// decision D5 a cluster whose saddlepoint failed has a finite L again — the
// substituted normal tail — so it contributes, and the pooled status carries
// its 3..6 code up to the meta row.  That is the intended behaviour: the
// substitution stays visible in `meta_SPA_STATUS` rather than being laundered
// by pooling.  A cluster with no test at all still has a NaN L and is still
// skipped.
//
// When the pooled variance is not positive there is no pooled statistic, and
// the status is NA_POST_FAIL — decision D4 lists "a meta pool with sum Var <= 0"
// under 7 by name, and 7 is what CLAUDE.md's status table has always said.
// This covers the empty pool too, and the post-audit repair is exactly that
// unification: Stage 2 wrote NA_NO_TEST here, which contradicted D4, the CLAUDE.md
// table and the LEAF `--help` block all at once, and made `meta_SPA_STATUS`
// claim that the MARKER has no statistic when what failed is the POOLING of
// the per-cluster ones.  Code 8 is reserved for "no statistic exists in this
// stratum" — no informative subject, a monomorphic stratum, Var(S) <= 0, a
// non-finite Z — which are properties of a cluster, not of the pool over them.
// Either way LOG10P is NA and there is no fallback: both codes are >= 7, so the
// C3 invariant is unaffected and no reader who filters on SPA_STATUS <= 2 or
// >= 7 sees any change.
struct MetaPooled {
    double negLog10p;
    spa::Status status;
};

inline MetaPooled metaPvalueScorePool(
    const std::vector<double> &scores,
    const std::vector<double> &negLog10ps,
    const std::vector<spa::Status> &statuses
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    double sumScore = 0.0, sumVar = 0.0;
    // Placeholder only: `any` guards every read of `st`, and a pool with no
    // contributor exits through the sum Var <= 0 branch below.
    spa::Status st = spa::Status::NaNoTest;
    bool any = false;
    const std::size_t K = scores.size();
    for (std::size_t c = 0; c < K; ++c) {
        if (std::isnan(scores[c]) || std::isnan(negLog10ps[c])) continue;
        double L = negLog10ps[c];
        if (L < META_L_FLOOR) L = META_L_FLOOR;
        const double chisq = math::chisq1FromNegLog10P(L);
        const double var   = (scores[c] * scores[c]) / chisq;
        sumScore += scores[c];
        sumVar   += var;
        st  = any ? spa::worseStatus(st, statuses[c]) : statuses[c];
        any = true;
    }
    if (sumVar <= 0.0) return MetaPooled{nan, spa::Status::NaPostFail};
    const double z = sumScore / std::sqrt(sumVar);
    return MetaPooled{spa::normalBranch(z).negLog10p, st};
}

}  // namespace math
