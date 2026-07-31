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
// Two departures from src/spacox/spacox.cpp's buildCumulantTable, both driven
// by SPAsqr needing one table per chromosome AND per quantile (22 x ntaus)
// where SPACox needs one per phenotype:
//
//   1. The table is built from a *binned* residual histogram, turning the
//      O(L*N) double loop into O(N) binning + O(L*B).  Each bin contributes
//      through its own mean, which makes the approximation second order: the
//      relative error on M0 is about (1/2) * s_b^2 * v^2 with s_b the
//      within-bin sd (bin width / sqrt(12)).  With the default 4096 bins over
//      the residual range and the |v| <~ 10 region the saddlepoint actually
//      visits, that is below 1e-7.  The exact residual min/max are kept
//      separately so the K1 asymptotes (and hence the support bound) stay
//      exact.
//
//   2. The per-grid-point reduction ships the scalar / AVX2 / AVX-512 triple
//      required by CLAUDE.md, dispatched through simdLevel().  SPACox's
//      scalar-only loop is a one-off per phenotype; 22 x ntaus tables is not.
//
// Evaluation buckets subjects by quantised genotype rather than iterating
// over non-zero entries, so a Newton step costs O(#occupied buckets) instead
// of O(#non-zero).  Hard calls collapse to at most four values (0, 1, 2 and
// the mean-imputed 2*altFreq); dosages collapse onto the quantisation grid.
// There is deliberately no hard-call/dosage branch — quantisation subsumes
// both, and it is what makes the empirical CGF affordable for common
// variants, where SPACox's non-zero loop degenerates to O(N).
#pragma once

#include "util/math_helper.hpp" // M_PI fallback

#include <Eigen/Dense>
#include <cmath>
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
// Per-marker genotype bucketing
// ══════════════════════════════════════════════════════════════════════
//
// Collapses a post-impute genotype column into (centred value, multiplicity)
// pairs on a quantisation grid of step 1/qScale over the dosage range [0, 2].
// Reused across markers: the histogram is cleared through a touched-key list,
// so the per-marker cost is O(n) plus O(#occupied buckets), never O(qScale).
//
// Two details make the collapse accurate rather than merely cheap:
//
//   * a bucket is represented by the *mean* of the genotypes it holds, not by
//     its grid point, so the first-order quantisation error cancels;
//   * the leftover within-bucket spread is returned as residualSS() =
//     sum_i Gtilde_i^2 - sum_b n_b * Gtilde_b^2, which the evaluator adds back
//     as a closed-form quadratic.  That makes K''(0) exactly
//     Var(R) * sum_i Gtilde_i^2 no matter how coarse the grid is, and leaves
//     only the (negligible) within-bucket cubic term.
//
// Hard calls put every subject on an exact grid point, so both corrections
// vanish identically and the collapse is bit-exact with at most four buckets
// (0, 1, 2 and the mean-imputed value).
class GenoBucketizer {
  public:
    explicit GenoBucketizer(int qScale = 256)
        : m_q(qScale),
          m_cnt(static_cast<size_t>(2 * qScale) + 1, 0.0),
          m_sum(static_cast<size_t>(2 * qScale) + 1, 0.0) {
        m_touched.reserve(64);
        m_gt.reserve(64);
        m_ct.reserve(64);
    }

    // `g` is the post-impute genotype column (length n), `gMean` its mean.
    // Values are clamped into [0, 2] before quantisation; the readers already
    // clamp dosages, but a degenerate imputed value would otherwise index out
    // of bounds.
    int build(const double *g, int n, double gMean) {
        const double qd  = static_cast<double>(m_q);
        const int    hiK = 2 * m_q;
        double sumSq = 0.0;
        for (int i = 0; i < n; ++i) {
            double v = g[i];
            if (!(v > 0.0)) v = 0.0;       // also maps NaN to 0
            else if (v > 2.0) v = 2.0;
            sumSq += v * v;
            int k = static_cast<int>(v * qd + 0.5);
            if (k < 0) k = 0;
            else if (k > hiK) k = hiK;
            const size_t uk = static_cast<size_t>(k);
            if (m_cnt[uk] == 0.0) m_touched.push_back(k);
            m_cnt[uk] += 1.0;
            m_sum[uk] += v;
        }

        return finalise(n, sumSq, gMean);
    }

    // Same as build(), but reads straight out of a union-space column and
    // skips subjects absent from this phenotype.  Bucketing needs only the
    // multiset of genotype values, not the local ordering, so this avoids
    // materialising a scattered per-phenotype copy first — one pass over the
    // column instead of two.
    int buildMapped(
        const double *gUnion,
        const uint32_t *unionToLocal,
        uint32_t nUnion,
        double gMean
    ) {
        const double qd  = static_cast<double>(m_q);
        const int    hiK = 2 * m_q;
        double sumSq = 0.0;
        int    nSeen = 0;
        for (uint32_t i = 0; i < nUnion; ++i) {
            if (unionToLocal[i] == UINT32_MAX) continue;
            double v = gUnion[i];
            if (!(v > 0.0)) v = 0.0;
            else if (v > 2.0) v = 2.0;
            sumSq += v * v;
            ++nSeen;
            int k = static_cast<int>(v * qd + 0.5);
            if (k < 0) k = 0;
            else if (k > hiK) k = hiK;
            const size_t uk = static_cast<size_t>(k);
            if (m_cnt[uk] == 0.0) m_touched.push_back(k);
            m_cnt[uk] += 1.0;
            m_sum[uk] += v;
        }
        return finalise(nSeen, sumSq, gMean);
    }

  private:
    int finalise(int n, double sumSq, double gMean) {
        m_gt.clear();
        m_ct.clear();
        double bucketSS = 0.0;
        for (int k : m_touched) {
            const size_t uk = static_cast<size_t>(k);
            const double c  = m_cnt[uk];
            const double gt = m_sum[uk] / c - gMean;   // bucket mean, centred
            m_gt.push_back(gt);
            m_ct.push_back(c);
            bucketSS += c * gt * gt;
            m_cnt[uk] = 0.0;
            m_sum[uk] = 0.0;
        }
        m_touched.clear();

        // sum_i (g_i - gMean)^2 computed exactly from the raw column.
        const double exactSS = sumSq - static_cast<double>(n) * gMean * gMean;
        m_residSS = exactSS - bucketSS;
        if (m_residSS < 0.0) m_residSS = 0.0;   // guard rounding
        return static_cast<int>(m_gt.size());
    }

  public:
    const double *gtilde() const {
        return m_gt.data();
    }

    const double *count() const {
        return m_ct.data();
    }

    int size() const {
        return static_cast<int>(m_gt.size());
    }

    // Within-bucket sum of squares not captured by the bucket means.
    double residualSS() const {
        return m_residSS;
    }

  private:
    int m_q;
    std::vector<double> m_cnt, m_sum;
    std::vector<int>    m_touched;
    std::vector<double> m_gt, m_ct;
    double m_residSS = 0.0;
};

// ══════════════════════════════════════════════════════════════════════
// Bucketed CGF evaluation
// ══════════════════════════════════════════════════════════════════════

// The within-bucket spread left over by quantisation enters as a pure
// quadratic: for |v| small, K0(v) = (1/2) Var(R) v^2 + O(v^3) because the
// residuals are centred (K0'(0) = mean R = 0).  `resSS` is
// GenoBucketizer::residualSS(); it is identically zero for hard calls.

// K'(t) - target and K''(t) in one pass (one table lookup per bucket).
inline void bucketedK1K2(
    const EmpiricalCgf &tab,
    const double *gt,
    const double *ct,
    int nb,
    double resSS,
    double t,
    double target,
    double &outK1,
    double &outK2
) {
    double s1 = 0.0, s2 = 0.0;
    for (int b = 0; b < nb; ++b) {
        const double gtb = gt[b];
        double v1, v2;
        tab.k1k2(t * gtb, v1, v2);
        const double w = ct[b] * gtb;
        s1 += w * v1;
        s2 += w * gtb * v2;
    }
    const double q = tab.varResid * resSS;
    outK1 = s1 + q * t - target;
    outK2 = s2 + q;
}

inline double bucketedK0(
    const EmpiricalCgf &tab,
    const double *gt,
    const double *ct,
    int nb,
    double resSS,
    double t
) {
    double s = 0.0;
    for (int b = 0; b < nb; ++b)
        s += ct[b] * tab.k0(t * gt[b]);
    return s + 0.5 * tab.varResid * resSS * t * t;
}

// K''(0) = Var(R) * sum_i Gtilde_i^2 — the variance the CGF encodes.  With the
// residual term included this is exact regardless of the quantisation grid.
inline double bucketedK2AtZero(
    const EmpiricalCgf &tab,
    const double *gt,
    const double *ct,
    int nb,
    double resSS
) {
    double s = resSS;
    for (int b = 0; b < nb; ++b)
        s += ct[b] * gt[b] * gt[b];
    return tab.varResid * s;
}

// Support of the residual-randomised score: K'(t) is bounded, so the
// saddlepoint equation K'(zeta) = q has no solution outside
// [lower, upper].  O(#buckets), used to reject before Newton spins.
inline void bucketedSupport(
    const EmpiricalCgf &tab,
    const double *gt,
    const double *ct,
    int nb,
    double &lower,
    double &upper
) {
    double up = 0.0, lo = 0.0;
    for (int b = 0; b < nb; ++b) {
        const double w = ct[b] * gt[b];
        up += w * (gt[b] > 0.0 ? tab.rMax : tab.rMin);
        lo += w * (gt[b] > 0.0 ? tab.rMin : tab.rMax);
    }
    upper = up;
    lower = lo;
}

} // namespace ecgf
