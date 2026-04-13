// conquer.hpp — Smoothed quantile regression (Gaussian kernel, low-dimensional)
//
// Pure C++17 / Eigen port of the conquer R package (He et al., 2021).
// Only the Gaussian kernel variant (`smqrGauss`) is implemented here —
// that is the only kernel used by the SPAsqr workflow.
//
// Reference:
//   Xuming He, Xiaoou Pan, Kean Ming Tan, Wen-Xin Zhou (2021).
//   "conquer: Convolution-Type Smoothed Quantile Regression."
//   R package version 1.3.3, https://CRAN.R-project.org/package=conquer
#pragma once

#include <Eigen/Dense>
#include <vector>

namespace conquer {

// Run smoothed quantile regression for a single tau.
//
//   X  : n × p  design matrix (raw covariates, NOT including intercept)
//   Y  : n × 1  response vector
//   tau: quantile level in (0,1)
//   h  : bandwidth (> 0).  Typical choice: IQR(Y) / 3.
//
// Returns (p+1) coefficient vector:  beta(0) = intercept,  beta(1..p) = slopes.
// Also stores the final n-vector of residuals in `residOut` if non-null.
Eigen::VectorXd smqrGauss(
    const Eigen::MatrixXd &X,
    const Eigen::VectorXd &Y,
    double tau,
    double h,
    Eigen::VectorXd *residOut = nullptr,
    double tol = 1e-7,
    int iteMax = 5000,
    double stepMax = 100.0
);

} // namespace conquer
