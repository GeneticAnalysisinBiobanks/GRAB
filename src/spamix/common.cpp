// common.cpp — Shared SPA helpers and outlier detection for SPAmix / SPAmixPlus

#include "spamix/common.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "util/math_helper.hpp"

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

// Evaluate outlier CGF K0 and K2 at saddlepoint parameter t.
inline std::pair<double, double> evalOutlierK0K2(
    double t,
    const double *resid,
    const double *maf,
    int n
) {
    double K0 = 0.0, K2 = 0.0;
    for (int i = 0; i < n; ++i) {
        const double tr = t * resid[i];
        double k0, k1, k2;
        math::kG012(tr, maf[i], k0, k1, k2);
        K0 += k0;
        K2 += resid[i] * resid[i] * k2;
    }
    return {K0, K2};
}

// Evaluate outlier K1_adj (= K1_outlier − s) and K2 at saddlepoint parameter t.
inline std::pair<double, double> evalOutlierK1adjK2(
    double t,
    double s,
    const double *resid,
    const double *maf,
    int n
) {
    double K1_adj = -s, K2 = 0.0;
    for (int i = 0; i < n; ++i) {
        const double tr = t * resid[i];
        double k0, k1, k2;
        math::kG012(tr, maf[i], k0, k1, k2);
        K1_adj += resid[i] * k1;
        K2 += resid[i] * resid[i] * k2;
    }
    return {K1_adj, K2};
}

// Newton-Raphson root-finding for K'_total(t) = s.
struct RootResult {
    double root;
    bool converge;
    double K2;
};

RootResult fastGetRootK1(
    double s,
    const double *mafOutlier,
    const double *residOutlier,
    int nOutlier,
    double mean_nonOutlier,
    double var_nonOutlier
) {

    double x = 0.0, oldX;
    double K1 = 0.0, K2 = 0.0, oldK1;
    double diffX = std::numeric_limits<double>::infinity(), oldDiffX;
    bool converge = true;
    constexpr double tol = 0.001;
    constexpr int maxiter = 100;

    for (int iter = 0; iter < maxiter; ++iter) {
        oldX = x;
        oldDiffX = diffX;
        oldK1 = K1;

        auto [k1_adj, k2_val] = evalOutlierK1adjK2(x, s, residOutlier, mafOutlier, nOutlier);

        K1 = k1_adj + mean_nonOutlier + var_nonOutlier * x;
        K2 = k2_val + var_nonOutlier;

        diffX = -K1 / K2;

        if (!std::isfinite(K1)) {
            x = std::numeric_limits<double>::infinity();
            K2 = 0.0;
            break;
        }

        if (iter > 0 && ((K1 > 0) != (oldK1 > 0))) {
            while (std::abs(diffX) > std::abs(oldDiffX) - tol)
                diffX *= 0.5;
        }

        if (std::abs(diffX) < tol) break;
        x = oldX + diffX;

        if (iter == maxiter - 1) converge = false;
    }

    return {x, converge, K2};
}

} // anonymous namespace

// ======================================================================
// Outlier detection (IQR-based, R-default quantile)
// ======================================================================

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

    // Non-outlier indices
    std::vector<bool> isOutlier(N, false);
    for (uint32_t idx : posOutlier)
        isOutlier[idx] = true;
    std::vector<uint32_t> posNonOutlier;
    posNonOutlier.reserve(static_cast<size_t>(N) - posOutlier.size());
    for (int i = 0; i < N; ++i)
        if (!isOutlier[i]) posNonOutlier.push_back(static_cast<uint32_t>(i));

    // Build sub-vectors
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

// ======================================================================
// spa::getProbSpaG — Lugannani-Rice saddlepoint tail probability
// ======================================================================

double spa::getProbSpaG(
    const double *mafOutlier,
    const double *residOutlier,
    int nOutlier,
    double s,
    bool lowerTail,
    double mean_nonOutlier,
    double var_nonOutlier
) {

    auto rootRes = fastGetRootK1(s, mafOutlier, residOutlier, nOutlier, mean_nonOutlier, var_nonOutlier);
    double zeta = rootRes.root;

    auto [k0_outlier, k2_outlier] = evalOutlierK0K2(zeta, residOutlier, mafOutlier, nOutlier);
    double k0_total = k0_outlier + mean_nonOutlier * zeta + 0.5 * var_nonOutlier * zeta * zeta;
    double k2_total = k2_outlier + var_nonOutlier;

    double temp1 = zeta * s - k0_total;

    if (!std::isfinite(zeta) || temp1 < 0.0 || k2_total <= 0.0) return std::numeric_limits<double>::quiet_NaN();

    double signZ = (zeta > 0.0) ? 1.0 : ((zeta < 0.0) ? -1.0 : 0.0);
    double w = signZ * std::sqrt(2.0 * temp1);
    double v = zeta * std::sqrt(k2_total);

    if (w == 0.0 || v == 0.0 || (v / w) <= 0.0) return std::numeric_limits<double>::quiet_NaN();

    double signFactor = lowerTail ? 1.0 : -1.0;
    return math::pnorm(signFactor * (w + std::log(v / w) / w));
}
