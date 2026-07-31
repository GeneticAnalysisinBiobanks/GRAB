// empirical_cgf.cpp — table construction for the empirical residual CGF.
//
// See empirical_cgf.hpp for the model and for why this differs from
// src/spacox/spacox.cpp's buildCumulantTable.
#include "util/empirical_cgf.hpp"

#include "util/simd_dispatch.hpp"
#include "util/simd_math.hpp"

#include <algorithm>
#include <cmath>
#include <thread>

#if defined(__x86_64__) || defined(_M_X64)
#  include <immintrin.h>
#endif

namespace ecgf {
namespace {

// ── One grid point: reduce over bins ──────────────────────────────────
//
// Computes, with the exponent shifted by `m` for numerical stability,
//   S0 = sum_b w_b e^(r_b t - m)
//   S1 = sum_b w_b r_b e^(r_b t - m)
//   S2 = sum_b w_b r_b^2 e^(r_b t - m)
// The shift cancels out of K1 and K2 and contributes an additive `m` to K0,
// so any residual scale is safe (SPACox's unshifted loop overflows to inf
// once |r|*|t| exceeds ~709).
//
// `wb`, `wrb`, `wr2b` are the pre-multiplied weight arrays; `rb` the bin
// representatives.

struct BinSums {
    double s0, s1, s2;
};

BinSums reduceScalar(
    const double *rb, const double *wb, const double *wrb, const double *wr2b,
    int nb, double t, double m
) {
    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (int b = 0; b < nb; ++b) {
        const double e = std::exp(rb[b] * t - m);
        s0 += wb[b] * e;
        s1 += wrb[b] * e;
        s2 += wr2b[b] * e;
    }
    return {s0, s1, s2};
}

#if defined(__x86_64__) || defined(_M_X64)

__attribute__((target("avx2,fma")))
double hsum256(__m256d v) {
    __m128d lo = _mm256_castpd256_pd128(v);
    __m128d hi = _mm256_extractf128_pd(v, 1);
    lo = _mm_add_pd(lo, hi);
    return _mm_cvtsd_f64(_mm_add_sd(lo, _mm_unpackhi_pd(lo, lo)));
}

__attribute__((target("avx2,fma")))
BinSums reduceAvx2(
    const double *rb, const double *wb, const double *wrb, const double *wr2b,
    int nb, double t, double m
) {
    const __m256d vt = _mm256_set1_pd(t);
    const __m256d vm = _mm256_set1_pd(m);
    __m256d a0 = _mm256_setzero_pd(), a1 = _mm256_setzero_pd(), a2 = _mm256_setzero_pd();
    int b = 0;
    for (; b + 4 <= nb; b += 4) {
        const __m256d vr = _mm256_loadu_pd(rb + b);
        const __m256d ve = avx2_exp_pd(_mm256_sub_pd(_mm256_mul_pd(vr, vt), vm));
        a0 = _mm256_fmadd_pd(_mm256_loadu_pd(wb + b), ve, a0);
        a1 = _mm256_fmadd_pd(_mm256_loadu_pd(wrb + b), ve, a1);
        a2 = _mm256_fmadd_pd(_mm256_loadu_pd(wr2b + b), ve, a2);
    }
    BinSums out{hsum256(a0), hsum256(a1), hsum256(a2)};
    for (; b < nb; ++b) {
        const double e = std::exp(rb[b] * t - m);
        out.s0 += wb[b] * e;
        out.s1 += wrb[b] * e;
        out.s2 += wr2b[b] * e;
    }
    return out;
}

__attribute__((target("avx2,avx512f,avx512vl,fma")))
BinSums reduceAvx512(
    const double *rb, const double *wb, const double *wrb, const double *wr2b,
    int nb, double t, double m
) {
    const __m512d vt = _mm512_set1_pd(t);
    const __m512d vm = _mm512_set1_pd(m);
    __m512d a0 = _mm512_setzero_pd(), a1 = _mm512_setzero_pd(), a2 = _mm512_setzero_pd();
    int b = 0;
    for (; b + 8 <= nb; b += 8) {
        const __m512d vr = _mm512_loadu_pd(rb + b);
        const __m512d ve = avx512_exp_pd(_mm512_sub_pd(_mm512_mul_pd(vr, vt), vm));
        a0 = _mm512_fmadd_pd(_mm512_loadu_pd(wb + b), ve, a0);
        a1 = _mm512_fmadd_pd(_mm512_loadu_pd(wrb + b), ve, a1);
        a2 = _mm512_fmadd_pd(_mm512_loadu_pd(wr2b + b), ve, a2);
    }
    BinSums out{_mm512_reduce_add_pd(a0), _mm512_reduce_add_pd(a1), _mm512_reduce_add_pd(a2)};
    for (; b < nb; ++b) {
        const double e = std::exp(rb[b] * t - m);
        out.s0 += wb[b] * e;
        out.s1 += wrb[b] * e;
        out.s2 += wr2b[b] * e;
    }
    return out;
}

#endif // x86_64

using ReduceFn = BinSums (*)(const double *, const double *, const double *,
                             const double *, int, double, double);

ReduceFn pickReduce() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
        case SimdLevel::AVX512: return &reduceAvx512;
        case SimdLevel::AVX2:   return &reduceAvx2;
        default:                return &reduceScalar;
    }
#else
    return &reduceScalar;
#endif
}

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
    std::vector<double> rb, wb, wrb, wr2b;
    rb.reserve(static_cast<size_t>(nBins));
    wb.reserve(static_cast<size_t>(nBins));
    wrb.reserve(static_cast<size_t>(nBins));
    wr2b.reserve(static_cast<size_t>(nBins));
    const double invN = 1.0 / static_cast<double>(n);
    for (int k = 0; k < nBins; ++k) {
        const double c = binCnt[static_cast<size_t>(k)];
        if (c == 0.0) continue;
        const double r = binSum[static_cast<size_t>(k)] / c;
        const double w = c * invN;
        rb.push_back(r);
        wb.push_back(w);
        wrb.push_back(w * r);
        wr2b.push_back(w * r * r);
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
    Eigen::VectorXd yK0(L), yK1(L), yK2(L);
    const ReduceFn reduce = pickReduce();
    const double *rbp = rb.data(), *wbp = wb.data(), *wrbp = wrb.data(), *wr2bp = wr2b.data();

    auto fillRange = [&](int lo, int hi) {
        for (int i = lo; i < hi; ++i) {
            const double t = xGrid[i];
            // Exponent shift: the largest r*t over the bins.
            const double m = (t >= 0.0) ? ct.rMax * t : ct.rMin * t;
            const BinSums s = reduce(rbp, wbp, wrbp, wr2bp, nb, t, m);
            if (!(s.s0 > 0.0)) {   // unreachable after the shift; keep it total
                yK0[i] = m;
                yK1[i] = (t >= 0.0) ? ct.rMax : ct.rMin;
                yK2[i] = 0.0;
                continue;
            }
            const double k1 = s.s1 / s.s0;
            yK0[i] = m + std::log(s.s0);
            yK1[i] = k1;
            yK2[i] = s.s2 / s.s0 - k1 * k1;
        }
    };

    const unsigned hw = std::thread::hardware_concurrency();
    const int nThread = std::min<int>(static_cast<int>(hw ? hw : 1u), std::max(1, L / 512));
    if (nThread <= 1) {
        fillRange(0, L);
    } else {
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(nThread));
        const int chunk = (L + nThread - 1) / nThread;
        for (int w = 0; w < nThread; ++w) {
            const int lo = w * chunk;
            const int hi = std::min(L, lo + chunk);
            if (lo >= hi) break;
            pool.emplace_back(fillRange, lo, hi);
        }
        for (auto &th : pool) th.join();
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
