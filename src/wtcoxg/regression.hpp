// regression.hpp — WtCoxG / LEAF case-control sampling-weight correction
//
// Provides:
//   calRegrWeight        — case-control sampling weights from prevalence
//
// The generic regression residual functions (logisticResiduals, coxResiduals,
// linearResiduals, cumulativeLogitFit) used to live here but have been moved
// to src/util/regression.hpp so that SPACox / SPAGRM / SPAmix can consume them
// without depending on WtCoxG.  calRegrWeight remains in this header because
// its semantics (control-row reweighting from disease prevalence) are tied to
// the case-control sampling design that WtCoxG / LEAF assume.
#pragma once

#include <Eigen/Dense>

#include "util/regression.hpp"

namespace regression {

// ──────────────────────────────────────────────────────────────────────
// calRegrWeight — case-control sampling-weight correction
//
// Equivalent to R:
//   sample_ratio <- sum(Indicator) / sum(1 - Indicator)
//   population_ratio <- RefPrevalence / (1 - RefPrevalence)
//   weight[control] <- sample_ratio / population_ratio
//
// prevalence: disease prevalence in reference population (0,1)
// indicator:  (n) case/control indicator (1 = case, 0 = control)
// Returns: weight vector (n)
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd calRegrWeight(
    double prevalence,
    const Eigen::Ref<const Eigen::VectorXd> &indicator
);

} // namespace regression
