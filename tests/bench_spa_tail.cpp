// bench_spa_tail.cpp — What deleting the linear tail costs, and saves.
//
// log10p_unify Stage 3 removed `spa::bnTail` and `spa::combineTails`, the
// linear-scale twins of `spa::bnTailLog` and `spa::combineTailsLog`.  Until
// then every tier-3 adapter called BOTH per tail:
//
//     p[k]    = spa::bnTail   (zeta, s, K0, K2, lower, stLin);
//     logp[k] = spa::bnTailLog(zeta, s, K0, K2, lower, stLog);
//
// and both ran the whole of `spa::detail::tailAbscissa` — the finiteness
// checks, the K'' > 0 and temp >= 0 guards, a `sqrt`, the |w| singularity
// test, a second `sqrt`, a division, a `log` and a second division — diverging
// only at the last step, into `math::pnorm` on one side and `math::pnormLog`
// on the other.  The abscissa work was therefore performed exactly twice for
// every tail of every marker that reached the saddlepoint.
//
// This benchmark measures that.  The "before" arm is reconstructed here rather
// than imported, because the function it calls no longer exists: `linearTail`
// below is `bnTail` verbatim, kept in the benchmark where a second copy is a
// measurement device rather than a duplicate implementation.
//
// Three rows are reported, in nanoseconds per two-sided p-value (i.e. per
// marker that enters the SPA branch, two tails included):
//
//   before   two tails x (bnTail + bnTailLog), then combineTails and
//            combineTailsLog, then the hand splice — what Stage 2 shipped.
//   after    two tails x bnTailLog, then spa::assemble — what Stage 3 ships.
//
// A third arm measured the `spa::pFromNegLog10P` call that produced the P
// column on top of `after`.  Stage 8 deleted the column and the helper, so
// `after` IS the cost paid today and the arm was removed rather than retimed.
//
// WHAT THE NUMBER MEANS.  The saving lands only on markers that take the SPA
// branch — 2 % to 8 % of markers at the default --spa-z-threshold — and the
// tail is a small fraction of what such a marker costs, since the Newton loop
// evaluates the CGF over every outlier subject several times per tail.  So the
// end-to-end effect is small by construction, and this benchmark exists to say
// how small rather than to claim a speed-up.  The reason to delete the linear
// tail is decision D1, not throughput; the throughput change is a consequence.
//
// Build:  make bench   (or: make build/tests/bench_spa_tail && ./build/tests/bench_spa_tail)

#include "util/spa.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <random>
#include <vector>

namespace {

// `spa::bnTail` as it stood at log10p_unify Stage 2, reproduced so that the
// "before" arm measures the work that was actually done.  It is deliberately
// written against the same `spa::detail::tailAbscissa` the surviving tail uses,
// because that is the shared body whose double evaluation is the subject.
double linearTail(
    double zeta, double s, double K0, double K2, bool lowerTail,
    spa::Status &st
) {
    const double x = spa::detail::tailAbscissa(zeta, s, K0, K2, st);
    if (!std::isfinite(x)) return std::numeric_limits<double>::quiet_NaN();
    return math::pnorm(x, 0.0, 1.0, lowerTail);
}

// `spa::combineTails` as it stood at Stage 2, for the same reason.
struct LinearTwoSided {
    double p;
    double negLog10p;
    spa::Status status;
};

LinearTwoSided linearCombine(
    double pUpper, double pLower, spa::Status sUpper, spa::Status sLower
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const spa::Status st = spa::worseStatus(sUpper, sLower);
    if (!spa::statusIsUsable(st)) return LinearTwoSided{nan, nan, st};
    if (!std::isfinite(pUpper) || !std::isfinite(pLower))
        return LinearTwoSided{nan, nan, spa::Status::FallbackNonFinite};
    double p = pUpper + pLower;
    if (p > 1.0) p = 1.0;
    if (p < 0.0) p = 0.0;
    double neg = -std::log10(p);
    if (neg == 0.0) neg = 0.0;
    return LinearTwoSided{p, neg, st};
}

// A population of converged saddlepoints spanning the range the production
// adapters actually hand the tail: |w| from just outside the singularity band
// to well past the linear underflow, with v/w bounded away from 1 so the
// log(v/w) division is real work rather than a constant.
struct TailPoint {
    double zetaU, sU, K0U, K2U;
    double zetaL, sL, K0L, K2L;
    double zNorm;
};

std::vector<TailPoint> drawPopulation(int n) {
    std::mt19937_64 rng(20260801ULL);
    std::uniform_real_distribution<double> uw(0.05, 45.0);
    std::uniform_real_distribution<double> uk(0.2, 8.0);
    std::vector<TailPoint> out;
    out.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        // Gaussian construction: s = var*zeta, K0 = var*zeta^2/2 gives w
        // exactly, and an independently drawn K2 makes v/w != 1.
        const double w = uw(rng), var = uk(rng);
        const double zeta = w / std::sqrt(var);
        TailPoint t;
        t.zetaU = zeta;  t.sU = var * zeta;  t.K0U = 0.5 * var * zeta * zeta;
        t.K2U = uk(rng);
        t.zetaL = -zeta; t.sL = -var * zeta; t.K0L = t.K0U;
        t.K2L = uk(rng);
        t.zNorm = w;
        out.push_back(t);
    }
    return out;
}

using Clock = std::chrono::steady_clock;

double nsPer(Clock::time_point a, Clock::time_point b, long n) {
    return std::chrono::duration<double, std::nano>(b - a).count() /
           static_cast<double>(n);
}

} // namespace

int main() {
    const int N = 20000;
    const int REPS = 40;
    const std::vector<TailPoint> pop = drawPopulation(N);
    const long total = static_cast<long>(N) * REPS;

    double sink = 0.0;

    // ── before ────────────────────────────────────────────────────────
    Clock::time_point t0 = Clock::now();
    for (int r = 0; r < REPS; ++r)
        for (const TailPoint &t : pop) {
            spa::Status stLinU = spa::Status::SpaOk, stLogU = spa::Status::SpaOk;
            spa::Status stLinL = spa::Status::SpaOk, stLogL = spa::Status::SpaOk;
            const double pu = linearTail(t.zetaU, t.sU, t.K0U, t.K2U, false, stLinU);
            const double lu = spa::bnTailLog(t.zetaU, t.sU, t.K0U, t.K2U, false, stLogU);
            const double pl = linearTail(t.zetaL, t.sL, t.K0L, t.K2L, true, stLinL);
            const double ll = spa::bnTailLog(t.zetaL, t.sL, t.K0L, t.K2L, true, stLogL);
            const spa::Status su = spa::worseStatus(stLinU, stLogU);
            const spa::Status sl = spa::worseStatus(stLinL, stLogL);
            const LinearTwoSided lin = linearCombine(pu, pl, su, sl);
            const spa::Result lg = spa::combineTailsLog(lu, ll, su, sl);
            sink += lin.p + lg.negLog10p;
        }
    Clock::time_point t1 = Clock::now();

    // ── after ─────────────────────────────────────────────────────────
    for (int r = 0; r < REPS; ++r)
        for (const TailPoint &t : pop) {
            spa::Status su = spa::Status::SpaOk, sl = spa::Status::SpaOk;
            const double lu = spa::bnTailLog(t.zetaU, t.sU, t.K0U, t.K2U, false, su);
            const double ll = spa::bnTailLog(t.zetaL, t.sL, t.K0L, t.K2L, true, sl);
            const spa::Result lg = spa::assemble(lu, ll, su, sl, t.zNorm);
            sink += lg.negLog10p;
        }
    Clock::time_point t2 = Clock::now();

    const double before = nsPer(t0, t1, total);
    const double after = nsPer(t1, t2, total);

    std::printf("\nspa tail path, ns per two-sided p-value (%d points x %d reps)\n\n",
                N, REPS);
    std::printf("  %-42s %8s  %8s\n", "arm", "ns", "vs before");
    std::printf("  %-42s %8.2f  %8s\n",
                "before  2x(bnTail+bnTailLog)+2 assemblies", before, "1.00x");
    std::printf("  %-42s %8.2f  %7.2fx\n",
                "after   2x bnTailLog + assemble", after, before / after);
    std::printf("\n  tailAbscissa evaluations per two-sided p-value: 4 -> 2\n");
    std::printf("  Realized only on the %.0f-%.0f%% of markers that enter the SPA\n"
                "  branch, and only against the tail's own share of such a\n"
                "  marker's cost; see the header note.\n\n", 2.0, 8.0);

    if (!std::isfinite(sink)) std::printf("  (sink %g)\n", sink);
    return 0;
}
