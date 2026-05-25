// common.cpp — Shared SPA helpers for SPAmix / SPAmixPlus
// Outlier detection lives in util/outlier.cpp.

#include "spamix/common.hpp"

#include <cmath>
#include <limits>
#include <utility>

#include "util/math_helper.hpp"

namespace {

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
