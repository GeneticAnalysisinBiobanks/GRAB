// common.hpp — Shared types and SPA helpers for SPAmix and SPAmixPlus
//
// OutlierData:    IQR-based outlier detection for residuals
// spa namespace:  Saddlepoint approximation tail probability
#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <vector>

// ======================================================================
// Outlier data — precomputed from residuals
// ======================================================================

struct OutlierData {
    std::vector<uint32_t> posOutlier;
    std::vector<uint32_t> posNonOutlier;
    Eigen::VectorXd residOutlier;     // resid[posOutlier]
    Eigen::VectorXd residNonOutlier;  // resid[posNonOutlier]
    Eigen::VectorXd resid2NonOutlier; // resid[posNonOutlier]^2
};

// Detect outlier residuals using IQR method.
// outlierRatio shrinks by 0.8x until at least one outlier is found.
OutlierData detectOutliers(
    const Eigen::VectorXd &resid,
    double outlierRatio
);

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
