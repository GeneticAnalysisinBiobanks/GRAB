// outlier.cpp — IQR-based residual outlier detection (shared by spamix / wtcoxg).

#include "util/outlier.hpp"

#include <algorithm>
#include <vector>

namespace {

// R-style type-7 quantile on a pre-sorted vector.
double quantile7(
    const std::vector<double> &sorted,
    double p
) {
    const int n = static_cast<int>(sorted.size());
    double h = (n - 1) * p;
    int j = static_cast<int>(h);
    double g = h - j;
    if (j + 1 >= n) return sorted[n - 1];
    return sorted[j] + g * (sorted[j + 1] - sorted[j]);
}

} // anonymous namespace

OutlierData detectOutliers(
    const Eigen::VectorXd &resid,
    double outlierRatio
) {
    const int N = static_cast<int>(resid.size());

    std::vector<double> sorted(resid.data(), resid.data() + N);
    std::sort(sorted.begin(), sorted.end());

    double q25 = quantile7(sorted, 0.25);
    double q75 = quantile7(sorted, 0.75);
    double iqr = q75 - q25;

    auto findOutliers = [&](double r) {
        double lo = q25 - r * iqr;
        double hi = q75 + r * iqr;
        std::vector<uint32_t> pos;
        for (int i = 0; i < N; ++i)
            if (resid[i] < lo || resid[i] > hi) pos.push_back(static_cast<uint32_t>(i));
        return pos;
    };

    auto posOutlier = findOutliers(outlierRatio);
    while (posOutlier.empty()) {
        outlierRatio *= 0.8;
        posOutlier = findOutliers(outlierRatio);
    }

    std::vector<bool> isOutlier(N, false);
    for (uint32_t idx : posOutlier)
        isOutlier[idx] = true;
    std::vector<uint32_t> posNonOutlier;
    posNonOutlier.reserve(static_cast<size_t>(N) - posOutlier.size());
    for (int i = 0; i < N; ++i)
        if (!isOutlier[i]) posNonOutlier.push_back(static_cast<uint32_t>(i));

    OutlierData od;
    od.posOutlier = std::move(posOutlier);
    od.posNonOutlier = std::move(posNonOutlier);

    const int nOut = static_cast<int>(od.posOutlier.size());
    const int nNon = static_cast<int>(od.posNonOutlier.size());
    od.residOutlier.resize(nOut);
    od.residNonOutlier.resize(nNon);
    od.resid2NonOutlier.resize(nNon);

    for (int i = 0; i < nOut; ++i)
        od.residOutlier[i] = resid[od.posOutlier[i]];
    for (int i = 0; i < nNon; ++i) {
        double r = resid[od.posNonOutlier[i]];
        od.residNonOutlier[i] = r;
        od.resid2NonOutlier[i] = r * r;
    }

    return od;
}
