// empirical_cgf.hpp — Empirical CGF of a residual vector, tabulated on a
// Cauchy-quantile grid.
//
// This is the SPACox construction (Bi et al.): instead of modelling the
// genotype as Binomial(2, mu) and holding the residuals fixed, the residual
// is treated as a draw from the empirical residual distribution and the
// genotype is held fixed.  For a centred genotype vector Gtilde the score
// S = sum_i R_i * Gtilde_i has CGF
//
//     K(t)   = sum_i K0(t * Gtilde_i),      K0(v) = log( mean_j exp(R_j v) )
//     K'(t)  = sum_i Gtilde_i    * K1(t * Gtilde_i),   K1 = M1/M0
//     K''(t) = sum_i Gtilde_i^2  * K2(t * Gtilde_i),   K2 = (M0*M2 - M1^2)/M0^2
//
// K0/K1/K2 depend only on the residuals, so they are tabulated once and
// interpolated per marker.  Storage is on the *log* scale (K0 = log M0) and
// interpolation is piecewise linear with constant extrapolation outside the
// grid, matching SPACox_Null_Model()'s approxfun(..., rule = 2).
//
// One departure from src/spacox/spacox.cpp's buildCumulantTable, driven by
// SPAsqr needing one table per chromosome AND per quantile (22 x ntaus) where
// SPACox needs one per phenotype: the table is built from a *binned* residual
// histogram, turning the O(L*N) double loop into O(N) binning + O(L*B).  Each
// bin contributes through its own mean, which makes the approximation second
// order: the relative error on M0 is about (1/2) * s_b^2 * v^2 with s_b the
// within-bin sd (bin width / sqrt(12)).  With the default 4096 bins over the
// residual range and the |v| <~ 10 region the saddlepoint actually visits,
// that is below 1e-7.  The exact residual min/max are kept separately so the
// K1 asymptotes (and hence the support bound) stay exact.
//
// Evaluation follows SPACox: a marker below the MAC cutoff has few carriers
// (subjects with a non-zero genotype), so K and its derivatives are summed
// over the carriers one by one, while the non-carriers all share the same
// centred value -gMean and enter as a single closed-form term n0 * K0(-t*gMean).
// A dosage genotype needs nothing extra: any subject with a non-zero dosage is
// a carrier.
#pragma once

#include "util/math_helper.hpp" // M_PI fallback

#include <Eigen/Dense>
#include <cmath>
#include <cstdint>
#include <vector>

namespace ecgf {

// ══════════════════════════════════════════════════════════════════════
// Tabulated empirical CGF
// ══════════════════════════════════════════════════════════════════════

struct EmpiricalCgf {
    Eigen::VectorXd xGrid;   // L, strictly increasing (Cauchy-quantile spaced)
    Eigen::VectorXd yK0, slopeK0;
    Eigen::VectorXd yK1, slopeK1;
    Eigen::VectorXd yK2, slopeK2;
    int    nGrid    = 0;
    double invScale = 0.0;   // 1 / gridScale — inverts the Cauchy map in O(1)
    double Lp1      = 0.0;   // L + 1
    double rMin     = 0.0;   // exact min residual — K1(-inf)
    double rMax     = 0.0;   // exact max residual — K1(+inf)
    double varResid = 0.0;   // K0''(0) = population variance of the residuals

    bool empty() const {
        return nGrid == 0;
    }

    // Inverse of the Cauchy-quantile grid: no binary search needed.
    int binIdx(double v) const {
        const double f = (std::atan(v * invScale) * (1.0 / M_PI) + 0.5) * Lp1 - 1.0;
        int lo = static_cast<int>(f);
        if (lo < 0) lo = 0;
        if (lo >= nGrid - 1) lo = nGrid - 2;
        return lo;
    }

    double k0(double v) const {
        const double *x = xGrid.data();
        if (v <= x[0]) return yK0.data()[0];
        if (v >= x[nGrid - 1]) return yK0.data()[nGrid - 1];
        const int lo = binIdx(v);
        return yK0.data()[lo] + (v - x[lo]) * slopeK0.data()[lo];
    }

    // Fused K1 + K2: one index computation, both values from the same bin.
    void k1k2(double v, double &o1, double &o2) const {
        const double *x = xGrid.data();
        if (v <= x[0]) {
            o1 = yK1.data()[0];
            o2 = yK2.data()[0];
            return;
        }
        if (v >= x[nGrid - 1]) {
            o1 = yK1.data()[nGrid - 1];
            o2 = yK2.data()[nGrid - 1];
            return;
        }
        const int lo = binIdx(v);
        const double dx = v - x[lo];
        o1 = yK1.data()[lo] + dx * slopeK1.data()[lo];
        o2 = yK2.data()[lo] + dx * slopeK2.data()[lo];
    }
};

// Build the table.  `resid` is the full residual vector for one (chromosome,
// quantile); `nBins` controls the binning accuracy, `L` and `rangeMax` match
// SPACox's defaults (10000 points over [-100, 100]).
EmpiricalCgf buildEmpiricalCgf(
    const double *resid,
    int n,
    int nBins    = 4096,
    int L        = 10000,
    double rangeMax = 100.0
);

// ══════════════════════════════════════════════════════════════════════
// Per-marker carrier set
// ══════════════════════════════════════════════════════════════════════
//
// The score S = sum_i R_i * (G_i - gMean) splits into the carriers, kept one
// by one, and the non-carriers, who all sit at the same centred value -gMean
// and contribute n0 * K0(-t*gMean) in closed form.  Mean-imputed missing
// genotypes are non-zero and therefore land among the carriers, which keeps
// the sum exact.
struct CarrierSet {
    std::vector<double> gt;   // centred genotype of each carrier, G_i - gMean
    double g0 = 0.0;          // centred genotype of a non-carrier, -gMean
    double n0 = 0.0;          // number of non-carriers

    bool empty() const {
        return gt.empty() && n0 == 0.0;
    }

    // Gather from a union-space column, skipping subjects absent from this
    // phenotype (unionToLocal[i] == UINT32_MAX).
    void gather(
        const double *gUnion,
        const uint32_t *unionToLocal,
        uint32_t nUnion,
        double gMean
    ) {
        gt.clear();
        g0 = -gMean;
        n0 = 0.0;
        for (uint32_t i = 0; i < nUnion; ++i) {
            if (unionToLocal[i] == UINT32_MAX) continue;
            const double v = gUnion[i];
            if (v == 0.0) n0 += 1.0;
            else          gt.push_back(v - gMean);
        }
    }
};

// ══════════════════════════════════════════════════════════════════════
// CGF of the residual-randomised score for one marker
// ══════════════════════════════════════════════════════════════════════

// K'(t) - target and K''(t) in one pass.
inline void empK1K2(
    const EmpiricalCgf &tab,
    const CarrierSet &cs,
    double t,
    double target,
    double &outK1,
    double &outK2
) {
    double s1 = 0.0, s2 = 0.0, v1, v2;
    for (const double g : cs.gt) {
        tab.k1k2(t * g, v1, v2);
        s1 += g * v1;
        s2 += g * g * v2;
    }
    tab.k1k2(t * cs.g0, v1, v2);
    outK1 = s1 + cs.n0 * cs.g0 * v1 - target;
    outK2 = s2 + cs.n0 * cs.g0 * cs.g0 * v2;
}

inline double empK0(const EmpiricalCgf &tab, const CarrierSet &cs, double t) {
    double s = 0.0;
    for (const double g : cs.gt) s += tab.k0(t * g);
    return s + cs.n0 * tab.k0(t * cs.g0);
}

// K''(0) = Var(R) * sum_i (G_i - gMean)^2 — the variance the CGF encodes.
inline double empK2AtZero(const EmpiricalCgf &tab, const CarrierSet &cs) {
    double s = cs.n0 * cs.g0 * cs.g0;
    for (const double g : cs.gt) s += g * g;
    return tab.varResid * s;
}

// Support of the residual-randomised score: K'(t) is bounded, so the
// saddlepoint equation K'(zeta) = q has no solution outside [lower, upper].
inline void empSupport(
    const EmpiricalCgf &tab,
    const CarrierSet &cs,
    double &lower,
    double &upper
) {
    auto add = [&](double g, double w, double &up, double &lo) {
        up += w * g * (g > 0.0 ? tab.rMax : tab.rMin);
        lo += w * g * (g > 0.0 ? tab.rMin : tab.rMax);
    };
    double up = 0.0, lo = 0.0;
    for (const double g : cs.gt) add(g, 1.0, up, lo);
    add(cs.g0, cs.n0, up, lo);
    upper = up;
    lower = lo;
}

} // namespace ecgf
