// regression.cpp — WtCoxG / LEAF case-control sampling-weight correction
//
// The logistic / Cox residual fits and their internal helpers
// (completeRows / subsetRows / subsetVec / validateSubset) live in
// src/util/regression.cpp.  This file only carries calRegrWeight, whose
// semantics are specific to the WtCoxG case-control design.
#include "wtcoxg/regression.hpp"

#include <stdexcept>

namespace regression {

// ──────────────────────────────────────────────────────────────────────
// calRegrWeight — case-control sampling-weight correction
// ──────────────────────────────────────────────────────────────────────

Eigen::VectorXd calRegrWeight(
    double prevalence,
    const Eigen::Ref<const Eigen::VectorXd> &indicator
) {
    const Eigen::Index N = indicator.size();
    double nCase = 0;
    for (Eigen::Index i = 0; i < N; ++i)
        if (indicator[i] > 0.5) ++nCase;
    double nCtrl = N - nCase;

    if (nCase < 1 || nCtrl < 1) throw std::runtime_error("calRegrWeight: need at least 1 case and 1 control");
    if (prevalence <= 0.0 || prevalence >= 1.0) throw std::runtime_error("calRegrWeight: prevalence must be in (0, 1)");

    // R equivalent:
    //   sample_ratio <- sum(Indicator) / sum(1 - Indicator)
    //   population_ratio <- RefPrevalence / (1 - RefPrevalence)
    //   weight[control] <- sample_ratio / population_ratio
    double sample_ratio = nCase / nCtrl;
    double population_ratio = prevalence / (1.0 - prevalence);
    double wCtrl = sample_ratio / population_ratio;

    Eigen::VectorXd w = Eigen::VectorXd::Ones(N);
    for (Eigen::Index i = 0; i < N; ++i)
        if (indicator[i] < 0.5) w[i] = wCtrl;
    return w;
}

} // namespace regression
