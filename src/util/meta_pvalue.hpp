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

inline double metaPvalueScorePool(
    const std::vector<double> &scores,
    const std::vector<double> &pvals
) {
    double sumScore = 0.0, sumVar = 0.0;
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
    }
    if (sumVar <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    const double z = sumScore / std::sqrt(sumVar);
    return std::erfc(std::fabs(z) / std::sqrt(2.0));
}

}  // namespace math
