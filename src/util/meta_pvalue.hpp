// meta_pvalue.hpp — Fixed-effects inverse-variance meta-analysis of
// per-cluster score statistics, with symmetric p-value clamping.
//
// The combining formula is the score-pooled meta:
//
//     Var_c    = S_c² / qchisq(p_c, 1, lower.tail = FALSE)
//     Z_meta   = Σ S_c / sqrt(Σ Var_c)
//     P_meta   = 2 · (1 − Φ(|Z_meta|))         (two-sided)
//
// Per-cluster p_c is clamped to [P_FLOOR, P_CEIL] before back-recovering
// Var_c.  This avoids two distinct numerical degeneracies:
//
//   p_c → 0 :  qchisq returns +∞ ⇒ Var_c = 0 ⇒ infinite weight, and
//              naively the sum of weights overflows.  Clamping at
//              P_FLOOR (= 1e-300) caps chisq at ≈ 1380, giving a very
//              large but finite weight that lets the cluster dominate
//              the meta as it should.
//   p_c → 1 :  qchisq returns 0 ⇒ Var_c = S²/0 = +∞ ⇒ the cluster
//              poisons the sum.  Clamping at P_CEIL (= 1 − 1e-15) keeps
//              chisq above ≈ 1.57e-30, giving a vanishingly small
//              weight that correctly contributes near zero.
//
// Clusters with NaN score or NaN p-value are skipped (their per-cluster
// computation failed upstream — that is a distinct signal from a finite
// boundary p-value and must not be folded into the clamp).

#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "util/math_helper.hpp"
#include "util/spa.hpp"

namespace math {

// Clamp bounds.  Asymmetric because in GWAS the lower tail (very small
// p) is the scientifically meaningful one, whereas the upper tail (p ≈ 1)
// is essentially a numerical artefact of the SPA rounding |z| ≈ 0.
//
//   P_FLOOR = 1e-300  ⇒  qchisq(P_FLOOR, 1, lower=F) ≈ 1380
//   P_CEIL  = 1 − 1e-15 ⇒  qchisq(P_CEIL,  1, lower=F) ≈ 1.57e-30
//
// Both bounds keep chisq strictly positive and finite, so the post-hoc
// "chisq < 1e-30 → 1e-30" guard formerly present in LEAFMethod::metaP
// becomes unnecessary.
constexpr double META_P_FLOOR = 1e-300;
constexpr double META_P_CEIL  = 1.0 - 1e-15;

// Pooled result: the p-value, its −log10, and the status of the POOLING.
//
// spa_unify Stage 6 (L3) adds the log-domain magnitude, so a pooled p below
// the linear underflow floor is still reported; and the status, which is the
// worst spa::Status among the clusters that actually contributed.  Contributing
// requires a finite score and a finite p, and a cluster whose saddlepoint
// failed now has a NaN p (D2) and is skipped — so a contributing cluster
// always carries one of the three non-failure statuses (Converged, GuardW or
// NormalBranch), and the pooled status is 0, 4 or 6 whenever `p` is a number,
// NonFinite when it is NA.  That is the same "P is NA for every status other
// than 0, 4 and 6" invariant the six migrated methods carry.
struct MetaPooled {
    double p;
    double negLog10p;
    spa::Status status;
};

inline MetaPooled metaPvalueScorePool(
    const std::vector<double> &scores,
    const std::vector<double> &pvals,
    const std::vector<spa::Status> &statuses
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    double sumScore = 0.0, sumVar = 0.0;
    spa::Status st = spa::Status::NonFinite;
    bool any = false;
    const std::size_t K = scores.size();
    for (std::size_t c = 0; c < K; ++c) {
        if (std::isnan(scores[c]) || std::isnan(pvals[c])) continue;
        double p = pvals[c];
        if (p < META_P_FLOOR) p = META_P_FLOOR;
        else if (p > META_P_CEIL) p = META_P_CEIL;
        const double chisq = math::qchisq(p, 1.0, false, false);
        const double var   = (scores[c] * scores[c]) / chisq;
        sumScore += scores[c];
        sumVar   += var;
        st  = any ? spa::worseStatus(st, statuses[c]) : statuses[c];
        any = true;
    }
    if (sumVar <= 0.0) return MetaPooled{nan, nan, spa::Status::NonFinite};
    const double z = sumScore / std::sqrt(sumVar);
    // The linear-scale p is left exactly as it was — erfc(|z|/sqrt(2)), not
    // 2*pnorm(-|z|) — so that meta_P does not move by a ULP under a change
    // whose purpose is to add a column.  The magnitude comes from the log
    // domain, where erfc has already underflowed for |z| beyond ~38.
    return MetaPooled{std::erfc(std::fabs(z) / std::sqrt(2.0)),
                      spa::normalBranch(z).negLog10p, st};
}

}  // namespace math
