// bench_spa_cgf.cpp — Microbenchmark for the diploid binomial CGF kernels.
//
// Stage 0 of the saddlepoint rework: records the "before" numbers against
// which Stages 2 onward are measured.  The kernel is the innermost loop of
// every saddlepoint method in GRAB; it executes once per outlier subject per
// Newton iteration per tail per marker per phenotype, so its cost dominates
// the SPA branch.
//
// Two axes are measured:
//
//   * form      — the production differencing algebra versus the stable
//                 algebra Stage 2 adopts.
//   * payload   — the full (K, K', K'') triple versus the (K', K'') pair the
//                 Newton loop actually consumes.  Every production kernel
//                 computes the logarithm on every iteration and discards it
//                 (src/spamix/common.cpp:44, src/spasqr/spasqr.cpp:88,
//                 src/wtcoxg/wtcoxg.cpp:79, src/localplus/spamixlocalp.cpp:890).
//
// Reported as nanoseconds per subject-evaluation.  The harness pins the
// working set below L2 so the measurement reflects arithmetic throughput
// rather than memory bandwidth, which is the regime the real kernels run in
// (outlier sets are small by construction).
//
// Build:  make bench       (or: make build/tests/bench_spa_cgf && ./build/tests/bench_spa_cgf)

#include "spa_reference.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

namespace {

// Prevent the optimizer from deleting a loop whose result is unused, without
// introducing a memory clobber that would also inhibit vectorization of the
// loop body itself.  Namespace scope rather than a function-local static so
// that -Wunused-but-set-variable does not fire; a volatile store is a side
// effect the compiler may not elide.
volatile double gSink = 0.0;

inline void sink(double v) { gSink = v; }

struct Result {
    double nsPerEval;
    double checksum;   // printed so the compiler cannot elide the work
};

template <class Fn>
Result timeKernel(
    const std::vector<double> &resid,
    double maf,
    int repeats,
    Fn &&fn
) {
    const int n = static_cast<int>(resid.size());

    // Warm-up: touch the data and let the branch predictors settle.
    double warm = 0.0;
    for (int i = 0; i < n; ++i) warm += fn(0.01, resid[i], maf);
    sink(warm);

    const auto t0 = std::chrono::steady_clock::now();
    double acc = 0.0;
    for (int rep = 0; rep < repeats; ++rep) {
        // Vary t across repeats so the saddlepoint argument sweeps the domain
        // the solver actually visits, rather than sitting at one value.
        const double t = 0.002 * static_cast<double>(rep % 64);
        for (int i = 0; i < n; ++i) acc += fn(t, resid[i], maf);
    }
    const auto t1 = std::chrono::steady_clock::now();
    sink(acc);

    const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
    return Result{ns / (static_cast<double>(repeats) * n), acc};
}

} // namespace

int main(int argc, char **argv) {
    int nOutlier = (argc > 1) ? std::atoi(argv[1]) : 512;
    int repeats  = (argc > 2) ? std::atoi(argv[2]) : 20000;
    if (nOutlier <= 0) nOutlier = 512;
    if (repeats  <= 0) repeats  = 20000;

    // Deterministic residuals: reproducibility across runs is required for the
    // benchmark to be comparable between stages.
    std::mt19937_64 rng(20260730u);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<double> resid(static_cast<size_t>(nOutlier));
    for (double &r : resid) r = nd(rng);

    std::printf("GRAB saddlepoint CGF microbenchmark\n");
    std::printf("  outliers = %d, repeats = %d, working set = %.1f KiB\n\n",
                nOutlier, repeats, nOutlier * sizeof(double) / 1024.0);

    std::printf("  %-34s %12s %10s\n", "kernel", "ns/subject", "rel");
    std::printf("  %-34s %12s %10s\n",
                "----------------------------------", "------------", "----------");

    double baseline = 0.0;

    struct Row { std::string name; Result res; };
    std::vector<Row> rows;

    for (double maf : {0.05, 0.30}) {
        char label[64];

        // Production algebra, full triple — what every kernel does today.
        std::snprintf(label, sizeof(label), "production K,K',K''  (MAF %.2f)", maf);
        rows.push_back({label, timeKernel(resid, maf, repeats,
            [](double t, double r, double p) {
                const sparef::Cgf c = sparef::cgfDiffForm(t, r, p);
                return c.K0 + c.K1 + c.K2;
            })});

        // Stable algebra, full triple.
        std::snprintf(label, sizeof(label), "stable     K,K',K''  (MAF %.2f)", maf);
        rows.push_back({label, timeKernel(resid, maf, repeats,
            [](double t, double r, double p) {
                const sparef::Cgf c = sparef::cgfStableForm(t, r, p, true);
                return c.K0 + c.K1 + c.K2;
            })});

        // Stable algebra, Newton payload only — no logarithm.
        std::snprintf(label, sizeof(label), "stable     K',K''    (MAF %.2f)", maf);
        rows.push_back({label, timeKernel(resid, maf, repeats,
            [](double t, double r, double p) {
                const sparef::Cgf c = sparef::cgfStableForm(t, r, p, false);
                return c.K1 + c.K2;
            })});
    }

    for (size_t i = 0; i < rows.size(); ++i) {
        if (i % 3 == 0) baseline = rows[i].res.nsPerEval;
        std::printf("  %-34s %12.3f %9.2fx\n",
                    rows[i].name.c_str(), rows[i].res.nsPerEval,
                    baseline / rows[i].res.nsPerEval);
        if (i % 3 == 2) std::printf("\n");
    }

    // Print the checksums so no kernel can be optimized away between builds.
    std::printf("  checksums:");
    for (const Row &r : rows) std::printf(" %.6e", r.res.checksum);
    std::printf("\n");

    return 0;
}
