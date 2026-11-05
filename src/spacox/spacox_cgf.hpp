// spacox_cgf.hpp — SPACox's empirical cumulant-generating function.
//
// This is **tier 3** of the three-tier split described in
// dev-notes/methods/spa_unify/02_design.md: the saddlepoint solver and the
// Barndorff-Nielsen tail live in util/spa.hpp and are shared by all six SPA
// sites, but the CGF itself is SPACox's own and is deliberately not forced
// into util/spa_cgf.hpp's binomial variants.  SPACox does not model the
// genotype as Binomial(2, p); it tilts the *empirical distribution of the
// null-model residuals*,
//
//     K_R(t) = log( (1/N) * sum_j exp(t * r_j) )
//
// tabulated once per phenotype on a Cauchy-quantile grid and read back by
// piecewise-linear interpolation.  The score S = sum_i a_i R_i then has
//
//     K(t)   = sum_i K_R(t * a_i)
//     K'(t)  = sum_i a_i   * K_R'(t * a_i)
//     K''(t) = sum_i a_i^2 * K_R''(t * a_i)
//
// The header is separate from spacox.hpp (which pulls in engine/marker.hpp and
// geno_factory/geno_data.hpp) so that tests/spacox_cgf_test.cpp can exercise
// the table build and both evaluators without linking the engine.
//
// Two changes here relative to the original in-place implementation, both from
// the spa_unify audit:
//
//   P1  The stage-1 weights a_i = (g_i - 2*altFreq)/sqrt(Var S) are a function
//       of g_i alone, so for hard-called input they take at most four distinct
//       values.  buildGenoWeights collapses them into (weight, multiplicity)
//       buckets plus a leftover list for genuine dosages, which turns each
//       Newton iteration from O(nNonZero) into O(4 + nDosage).
//
//   P5  buildCumulantTable's 10 000 grid points are independent, so the outer
//       loop is split across std::thread workers (never OpenMP - see
//       CLAUDE.md), and the inner reduction over subjects is expressed as
//       Eigen array operations over an L1-resident block so that Eigen's own
//       vectorized exp is used.  The moment accumulation is additionally
//       reformulated to be overflow-free and cancellation-free; the previous
//       form produced a NaN at the grid endpoint for residual vectors with
//       max|r| above ~3.5 (see buildCumulantTable's comment).

#pragma once

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <thread>
#include <vector>

#include "util/spa.hpp"

// ======================================================================
// Empirical CGF interpolation table
// ======================================================================

struct CumulantTable {
    Eigen::VectorXd xGrid; // length L, strictly increasing
    int nGrid = 0;
    double invScale = 0.0; // 1/gridScale, for the O(1) Cauchy-inverse lookup
    double Lp1 = 0.0;      // L + 1
    Eigen::VectorXd yK0, slopeK0;
    Eigen::VectorXd yK1, slopeK1;
    Eigen::VectorXd yK2, slopeK2;
};

namespace spacox_cgf {

// ──────────────────────────────────────────────────────────────────────
// § 1  Table construction
// ──────────────────────────────────────────────────────────────────────

// Grid length and half-range are fixed by the method, not configurable.
constexpr int    kGridLen  = 10000;
constexpr double kGridHalf = 100.0;

// Number of subjects processed per inner block.  512 doubles is 4 kB, so the
// exponentiated-weight scratch stays in L1 while the three reductions read it
// back; the residual vector itself is streamed once per grid point.
constexpr int kInnerBlock = 512;

// Build the empirical-CGF table from a residual vector.
//
// PARALLELISM (P5).  The grid loop's iterations are independent: each reads
// xGrid[i] and the const residual pointer and writes only yK0[i], yK1[i],
// yK2[i].  Splitting it into contiguous ranges across `nthread` std::thread
// workers therefore leaves every table entry bit-identical to the serial
// build, which is exactly the property a parallelization of a numerical
// routine should have.  (Vectorizing the *inner* reduction does change the
// summation order and hence the last bits; that change is taken deliberately,
// see below.)
//
// NUMERICS.  The original accumulated raw moments and formed
//
//     M0 = mean(exp(r t)),  M1 = mean(r exp(r t)),  M2 = mean(r^2 exp(r t))
//     K = log M0,  K' = M1/M0,  K'' = (M0 M2 - M1 M1) / (M0 M0)
//
// which fails in two ways at the ends of the grid, where |t| reaches 100:
//
//   * OVERFLOW.  With max|r| = 4.14 (the Cox martingale residuals of the
//     bundled Time_Event fixture) the endpoint exponent is r t = 414, so
//     M0 = 1.63e176 and the products M0*M2 and M1*M1 both overflow to +inf.
//     Their difference is NaN, and because interpK2 clamps to yK2[0] for
//     v <= xGrid[0] while computeSlopes turns a NaN entry into a NaN slope,
//     the whole first bin - a span of 50 in t - returns NaN.
//   * CANCELLATION.  K'' is a variance, and at large |t| the exponentially
//     tilted distribution collapses onto a single residual, so the variance
//     tends to zero while both M2/M0 and (M1/M0)^2 tend to r_extreme^2.  The
//     difference loses every significant digit long before it overflows: at
//     t = -50 the old form returns 1.3135e-13 against a reference of
//     1.3529e-13, a relative error of 2.9e-2.
//
// Both are removed by (a) subtracting the maximum exponent before
// exponentiating, which is the standard log-sum-exp shift and makes every
// weight lie in (0, 1] with the largest exactly 1, and (b) accumulating the
// first two moments *about a centre* c chosen close to the tilted mean, so
// that K'' = A - B^2 is evaluated where B is small rather than where B^2
// nearly equals A.  The centre is the first-order tilted mean
// rMean + t * rVar, clamped to [rMin, rMax] so that it converges to the
// dominant residual in exactly the regime where the tilt degenerates.  Both
// c and the shift are exact functions of scalars computed once, so the inner
// loop pays two extra flops per grid point, not per subject.
//
// K itself is recovered as m + log(W/N), which is also better conditioned
// than log(M0) near t = 0: there M0 = 1 + O(t) and the logarithm of a double
// rounded to 1 + delta loses relative accuracy of order eps/delta, whereas
// here the O(t) part is carried in the exactly-representable shift.
inline CumulantTable buildCumulantTable(
    const Eigen::VectorXd &residuals,
    int nthread
) {
    constexpr int L = kGridLen;

    // Cauchy-quantile spacing: x0 = tan(pi*(k/(L+1) - 0.5)), rescaled so that
    // max|x| = kGridHalf.  Concentrates grid points near t = 0, where the
    // saddlepoint actually lands, while still covering the far tail.
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
        gridScale = kGridHalf / maxAbs;
        xGrid *= gridScale;
    }

    const int N = static_cast<int>(residuals.size());
    const double *rp = residuals.data();
    const double invN = 1.0 / static_cast<double>(N);

    // Scalars the inner loop needs but must not recompute: the extreme
    // residuals (which give the maximum exponent for either tilt direction)
    // and the mean and population variance (which give the first-order
    // tilted mean used as the centring constant).
    double rMin = rp[0], rMax = rp[0], rSum = 0.0;
    for (int j = 0; j < N; ++j) {
        const double r = rp[j];
        if (r < rMin) rMin = r;
        if (r > rMax) rMax = r;
        rSum += r;
    }
    const double rMean = rSum * invN;
    double rSs = 0.0;
    for (int j = 0; j < N; ++j) {
        const double d = rp[j] - rMean;
        rSs += d * d;
    }
    const double rVar = rSs * invN;

    Eigen::VectorXd yK0(L), yK1(L), yK2(L);
    const double *xp = xGrid.data();
    double *p0 = yK0.data(), *p1 = yK1.data(), *p2 = yK2.data();

    auto worker = [&](int iBegin, int iEnd) {
        alignas(64) double ubuf[kInnerBlock];
        for (int i = iBegin; i < iEnd; ++i) {
            const double t = xp[i];

            // Log-sum-exp shift: max_j (r_j * t) is attained at rMax when
            // t >= 0 and at rMin otherwise, so it costs one multiply.
            const double m = t * ((t >= 0.0) ? rMax : rMin);

            // Centre: first-order tilted mean, clamped to the residual range.
            double c = rMean + t * rVar;
            if (c < rMin) c = rMin;
            if (c > rMax) c = rMax;

            double W = 0.0, S1 = 0.0, S2 = 0.0;
            for (int j0 = 0; j0 < N; j0 += kInnerBlock) {
                const int nb = std::min(kInnerBlock, N - j0);
                Eigen::Map<const Eigen::ArrayXd> r(rp + j0, nb);
                Eigen::Map<Eigen::ArrayXd> u(ubuf, nb);
                u = (r * t - m).exp();
                W  += u.sum();
                S1 += (u * (r - c)).sum();
                S2 += (u * (r - c).square()).sum();
            }

            const double invW = 1.0 / W;
            const double B = S1 * invW;
            const double A = S2 * invW;
            p0[i] = m + std::log(W * invN);
            p1[i] = c + B;
            // A - B^2 is a variance and cannot be negative mathematically;
            // clamp so that a rounding-level negative can never reach the
            // GuardCurv test as a spurious "K'' <= 0".
            const double k2 = A - B * B;
            p2[i] = (k2 > 0.0) ? k2 : 0.0;
        }
    };

    int nthr = nthread;
    if (nthr < 1) nthr = 1;
    if (nthr > L) nthr = L;
    if (nthr == 1) {
        worker(0, L);
    } else {
        const int per = (L + nthr - 1) / nthr;
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(nthr) - 1);
        for (int k = 1; k < nthr; ++k) {
            const int a = k * per;
            if (a >= L) break;
            pool.emplace_back(worker, a, std::min(L, a + per));
        }
        worker(0, std::min(L, per));
        for (auto &th : pool) th.join();
    }

    // Piecewise-linear slopes.  L-1 divisions; left serial.
    auto computeSlopes = [](const Eigen::VectorXd &x,
                            const Eigen::VectorXd &y) -> Eigen::VectorXd {
        const int n = static_cast<int>(x.size());
        Eigen::VectorXd s(n - 1);
        for (int i = 0; i < n - 1; ++i)
            s[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        return s;
    };

    CumulantTable ct;
    ct.xGrid = std::move(xGrid);
    ct.nGrid = L;
    ct.invScale = 1.0 / gridScale;
    ct.Lp1 = static_cast<double>(L + 1);
    ct.yK0 = std::move(yK0);
    ct.slopeK0 = computeSlopes(ct.xGrid, ct.yK0);
    ct.yK1 = std::move(yK1);
    ct.slopeK1 = computeSlopes(ct.xGrid, ct.yK1);
    ct.yK2 = std::move(yK2);
    ct.slopeK2 = computeSlopes(ct.xGrid, ct.yK2);
    return ct;
}

// ──────────────────────────────────────────────────────────────────────
// § 2  Interpolation
// ──────────────────────────────────────────────────────────────────────

// Raw table triple at one abscissa: K_R, K_R', K_R''.
struct TableTriple {
    double k0, k1, k2;
};

// Flat pointer view, built once per evaluator call so that the inner loop
// does not reload seven Eigen data pointers per subject.
struct TableView {
    const double *xp, *y0, *s0, *y1, *s1, *y2, *s2;
    int nG;
    double invScale, Lp1;

    explicit TableView(const CumulantTable &ct)
        : xp(ct.xGrid.data()),
          y0(ct.yK0.data()), s0(ct.slopeK0.data()),
          y1(ct.yK1.data()), s1(ct.slopeK1.data()),
          y2(ct.yK2.data()), s2(ct.slopeK2.data()),
          nG(ct.nGrid), invScale(ct.invScale), Lp1(ct.Lp1) {}
};

// One bin lookup serving K, K' and K''.  The bin index is the analytic
// inverse of the Cauchy-quantile grid, so the lookup is O(1) rather than a
// binary search — one std::atan, which is why collapsing the per-subject walk
// into buckets (P1) is worth doing.
//
// The index clamp is written so that a non-finite `v` cannot reach the
// double-to-int conversion: both comparisons are false for NaN, which selects
// lo = 0 rather than executing undefined behaviour.  v is finite for every
// call the solver makes (t is checked finite before use and every weight is
// finite), so this is a belt rather than a live path.
inline TableTriple interpAll(const TableView &tv, double v) noexcept {
    if (v <= tv.xp[0])
        return TableTriple{tv.y0[0], tv.y1[0], tv.y2[0]};
    if (v >= tv.xp[tv.nG - 1]) {
        const int e = tv.nG - 1;
        return TableTriple{tv.y0[e], tv.y1[e], tv.y2[e]};
    }
    const double fIdx =
        (std::atan(v * tv.invScale) * (1.0 / M_PI) + 0.5) * tv.Lp1 - 1.0;
    int lo;
    if (fIdx > 0.0)
        lo = (fIdx < static_cast<double>(tv.nG - 1))
                 ? static_cast<int>(fIdx)
                 : tv.nG - 2;
    else
        lo = 0;
    const double dx = v - tv.xp[lo];
    return TableTriple{tv.y0[lo] + dx * tv.s0[lo],
                       tv.y1[lo] + dx * tv.s1[lo],
                       tv.y2[lo] + dx * tv.s2[lo]};
}

// ──────────────────────────────────────────────────────────────────────
// § 3  Genotype weights (P1)
// ──────────────────────────────────────────────────────────────────────

// The at-most-four distinct stage-1 weights, plus whatever did not fit.
//
// a_i = (g_i - twoMAF) / sqrtVarS depends on g_i alone.  Hard-called input
// gives g_i in {0, 1, 2}; a subject whose genotype was missing receives the
// mean imputation, which the engine writes as exactly 2.0 * altFreq
// (engine/marker.cpp), and twoMAF here is the same expression, so the
// equality test is exact and those subjects form a fourth bucket whose weight
// is exactly 0.  Genuine dosages — BGEN and VCF input, and any reader below
// the hard-call threshold — fall through to `extra`, one entry each, so the
// evaluator stays exactly correct for continuous genotypes; only the *speed*
// benefit is lost, in proportion to how many subjects carry a dosage.
//
// Bucket 0 is always the g == 0 bucket, reproducing the original code's `N0`
// collapse term in the same position and hence bit-identically.
struct GenoWeights {
    double bw[4];   // bucket weights
    double bc[4];   // bucket multiplicities, exact in double for N < 2^53
    int nBucket = 0;
    std::vector<double> extra;  // per-subject weights, multiplicity one
};

inline void buildGenoWeights(
    const double *g,
    int N,
    double twoMAF,
    double sqrtVarS,
    GenoWeights &out
) {
    int n0 = 0, n1 = 0, n2 = 0, nImp = 0;
    out.extra.clear();
    for (int i = 0; i < N; ++i) {
        const double v = g[i];
        if (v == 0.0)            ++n0;
        else if (v == 1.0)       ++n1;
        else if (v == 2.0)       ++n2;
        else if (v == twoMAF)    ++nImp;
        else                     out.extra.push_back((v - twoMAF) / sqrtVarS);
    }

    const double invSd = 1.0 / sqrtVarS;
    out.nBucket = 0;
    auto add = [&](int cnt, double gval) {
        if (cnt == 0) return;
        out.bw[out.nBucket] = (gval - twoMAF) * invSd;
        out.bc[out.nBucket] = static_cast<double>(cnt);
        ++out.nBucket;
    };
    // g == 0 first, so that the leading term matches the original N0 term.
    add(n0, 0.0);
    add(n1, 1.0);
    add(n2, 2.0);
    // The imputed bucket's weight is exactly zero; spelling it as
    // (twoMAF - twoMAF) * invSd keeps the construction uniform and yields
    // +0.0, which contributes nothing to K' or K'' and the fixed
    // K_R(0) offset to K — exactly what the per-subject loop did.
    add(nImp, twoMAF);
}

// ──────────────────────────────────────────────────────────────────────
// § 4  Cumulant evaluation
// ──────────────────────────────────────────────────────────────────────

namespace detail {

struct Acc {
    double k0 = 0.0, k1 = 0.0, k2 = 0.0;
};

// HasCounts == false is the dense stage-2 path: every weight has
// multiplicity one, so the multiply is elided rather than executed with 1.0.
// The association order of each term reproduces the original evaluator
// exactly — ((cnt * a) * k1) and (((cnt * a) * a) * k2) — so the g == 0
// bucket and the dense path are bit-identical to the pre-migration code at
// the same abscissa.
template <bool HasCounts>
inline void accumulate(
    const TableView &tv,
    double t,
    const double *w,
    const double *cnt,
    int n,
    Acc &acc
) noexcept {
    double k0 = acc.k0, k1 = acc.k1, k2 = acc.k2;
    for (int b = 0; b < n; ++b) {
        const double a = w[b];
        const TableTriple y = interpAll(tv, t * a);
        if (HasCounts) {
            const double c = cnt[b];
            k0 += c * y.k0;
            k1 += (c * a) * y.k1;
            k2 += ((c * a) * a) * y.k2;
        } else {
            k0 += y.k0;
            k1 += a * y.k1;
            k2 += (a * a) * y.k2;
        }
    }
    acc.k0 = k0;
    acc.k1 = k1;
    acc.k2 = k2;
}

} // namespace detail

// Stage 1: bucketed.  Cost is O(nBucket + extra.size()) per call, which for
// hard-called input is O(4) regardless of sample size.
inline spa::K012 evalBucketed(
    const TableView &tv,
    double t,
    const GenoWeights &gw,
    double s
) noexcept {
    detail::Acc acc;
    detail::accumulate<true>(tv, t, gw.bw, gw.bc, gw.nBucket, acc);
    if (!gw.extra.empty())
        detail::accumulate<false>(tv, t, gw.extra.data(), nullptr,
                                  static_cast<int>(gw.extra.size()), acc);
    return spa::K012{acc.k0, acc.k1 - s, acc.k2};
}

// Stage 2: dense.  The covariate-adjusted weights are a projection of the
// genotype onto the orthogonal complement of the design matrix, so they carry
// no bucket structure and every subject must be visited.  This runs only for
// markers whose stage-1 p is at or below --covar-p-threshold (5e-5 by
// default), so its cost is amortized over the whole scan.
inline spa::K012 evalDense(
    const TableView &tv,
    double t,
    const double *w,
    int n,
    double s
) noexcept {
    detail::Acc acc;
    detail::accumulate<false>(tv, t, w, nullptr, n, acc);
    return spa::K012{acc.k0, acc.k1 - s, acc.k2};
}

} // namespace spacox_cgf
