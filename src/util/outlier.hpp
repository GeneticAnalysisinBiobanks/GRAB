// outlier.hpp — IQR-based residual outlier detection (shared by spamix / wtcoxg).
//
// Splits a residual vector into outlier and non-outlier subsets using the
// 1.5×IQR rule (or any other multiplier).  The non-outlier sub-vector lets
// downstream SPA implementations replace an O(N) empirical CGF sum with a
// closed-form Gaussian approximation of equal mean and variance.
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
