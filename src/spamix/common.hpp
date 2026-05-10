// common.hpp — Shared SPA helpers for SPAmix and SPAmixPlus
//
// OutlierData / detectOutliers live in util/outlier.hpp and are reused here.
// spa::getProbSpaG: Lugannani-Rice tail probability with outlier/non-outlier split.
#pragma once

#include "util/outlier.hpp"

// ======================================================================
// SPA helpers — Lugannani-Rice tail probability with outlier/non-outlier split
// ======================================================================

namespace spa {

// Saddlepoint tail probability for the score statistic.
//   mafOutlier / residOutlier: per-individual MAF and residuals at outlier positions
//   s:            score statistic value to evaluate
//   lowerTail:    true → P(S ≤ s), false → P(S > s)
//   mean_nonOutlier / var_nonOutlier: normal-approx terms for non-outlier part
double getProbSpaG(
    const double *mafOutlier,
    const double *residOutlier,
    int nOutlier,
    double s,
    bool lowerTail,
    double mean_nonOutlier,
    double var_nonOutlier
);

} // namespace spa
