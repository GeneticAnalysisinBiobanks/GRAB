// empirical_cgf.cpp — table construction for the empirical residual CGF.
//
// See empirical_cgf.hpp for the model and for why this differs from
// src/spacox/spacox.cpp's buildCumulantTable.
#include "util/empirical_cgf.hpp"

#include <algorithm>
#include <cmath>

namespace ecgf {
namespace {

Eigen::VectorXd computeSlopes(const Eigen::VectorXd &x, const Eigen::VectorXd &y) {
    const int n = static_cast<int>(x.size());
    Eigen::VectorXd s(n - 1);
    for (int i = 0; i < n - 1; ++i)
        s[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    return s;
}

} // namespace

EmpiricalCgf buildEmpiricalCgf(
    const double *resid,
    int n,
    int nBins,
    int L,
    double rangeMax
) {
    EmpiricalCgf ct;
    if (n <= 0 || L < 2 || nBins < 1) return ct;

    // ── 1. Bin the residuals ──────────────────────────────────────────
    double rlo = resid[0], rhi = resid[0], rsum = 0.0;
    for (int i = 0; i < n; ++i) {
        const double r = resid[i];
        if (r < rlo) rlo = r;
        if (r > rhi) rhi = r;
        rsum += r;
    }
    ct.rMin = rlo;
    ct.rMax = rhi;

    const double span = rhi - rlo;
    std::vector<double> binCnt(static_cast<size_t>(nBins), 0.0);
    std::vector<double> binSum(static_cast<size_t>(nBins), 0.0);
    if (span <= 0.0) {
        binCnt[0] = static_cast<double>(n);
        binSum[0] = rsum;
    } else {
        const double inv = static_cast<double>(nBins) / span;
        for (int i = 0; i < n; ++i) {
            int k = static_cast<int>((resid[i] - rlo) * inv);
            if (k < 0) k = 0;
            else if (k >= nBins) k = nBins - 1;
            binCnt[static_cast<size_t>(k)] += 1.0;
            binSum[static_cast<size_t>(k)] += resid[i];
        }
    }

    // Occupied bins only; each contributes through its own mean, which makes
    // the binning second-order accurate.
    std::vector<double> rb, wb;
    rb.reserve(static_cast<size_t>(nBins));
    wb.reserve(static_cast<size_t>(nBins));
    const double invN = 1.0 / static_cast<double>(n);
    for (int k = 0; k < nBins; ++k) {
        const double c = binCnt[static_cast<size_t>(k)];
        if (c == 0.0) continue;
        const double r = binSum[static_cast<size_t>(k)] / c;
        const double w = c * invN;
        rb.push_back(r);
        wb.push_back(w);
    }
    const int nb = static_cast<int>(rb.size());

    // K0''(0) = population variance of the residuals.  Taken from the raw
    // residuals rather than the bins so it is exact — it is the quantity the
    // variance ratio is built from.
    const double rbar = rsum * invN;
    double var = 0.0;
    for (int i = 0; i < n; ++i) {
        const double d = resid[i] - rbar;
        var += d * d;
    }
    ct.varResid = var * invN;

    // ── 2. Cauchy-quantile grid, rescaled so max|x| = rangeMax ────────
    Eigen::VectorXd xGrid(L);
    double gridScale;
    {
        double maxAbs = 0.0;
        for (int k = 0; k < L; ++k) {
            const double p = static_cast<double>(k + 1) / static_cast<double>(L + 1);
            xGrid[k] = std::tan(M_PI * (p - 0.5));
            const double a = std::fabs(xGrid[k]);
            if (a > maxAbs) maxAbs = a;
        }
        gridScale = rangeMax / maxAbs;
        xGrid *= gridScale;
    }

    // ── 3. Reduce over bins at every grid point ───────────────────────
    // With the exponent shifted by m = max_b r_b t for numerical stability:
    //   S0 = sum_b w_b e^(r_b t - m),  S1 = sum_b w_b r_b e^(...),  S2 = sum_b w_b r_b^2 e^(...)
    // The shift cancels out of K1 and K2 and adds m to K0, so any residual
    // scale is safe (SPACox's unshifted loop overflows once |r t| > ~709).
    Eigen::VectorXd yK0(L), yK1(L), yK2(L);
    const Eigen::Map<const Eigen::ArrayXd> rbA(rb.data(), nb), wbA(wb.data(), nb);
    Eigen::ArrayXd e(nb);
    for (int i = 0; i < L; ++i) {
        const double t = xGrid[i];
        const double m = (t >= 0.0) ? ct.rMax * t : ct.rMin * t;
        e = (rbA * t - m).exp() * wbA;
        const double s0 = e.sum();
        const double s1 = (e * rbA).sum();
        const double s2 = (e * rbA * rbA).sum();
        if (!(s0 > 0.0)) {   // unreachable after the shift; keep it total
            yK0[i] = m;
            yK1[i] = (t >= 0.0) ? ct.rMax : ct.rMin;
            yK2[i] = 0.0;
            continue;
        }
        const double k1 = s1 / s0;
        yK0[i] = m + std::log(s0);
        yK1[i] = k1;
        yK2[i] = s2 / s0 - k1 * k1;
    }

    ct.xGrid   = std::move(xGrid);
    ct.nGrid   = L;
    ct.invScale = 1.0 / gridScale;
    ct.Lp1     = static_cast<double>(L + 1);
    ct.slopeK0 = computeSlopes(ct.xGrid, yK0);
    ct.slopeK1 = computeSlopes(ct.xGrid, yK1);
    ct.slopeK2 = computeSlopes(ct.xGrid, yK2);
    ct.yK0     = std::move(yK0);
    ct.yK1     = std::move(yK1);
    ct.yK2     = std::move(yK2);
    return ct;
}

} // namespace ecgf
