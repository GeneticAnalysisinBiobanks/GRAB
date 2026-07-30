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
// ── Stage 2 addition ─────────────────────────────────────────────────
//
// PART I below is the Stage 0 harness, unchanged, timing the reference forms in
// tests/spa_reference.hpp through a per-subject inlined lambda.  Those rows are
// the recorded "before" numbers and must keep measuring the same thing.
//
// PART II times the REAL kernels in src/util/spa_cgf.{hpp,cpp}: three variants
// x two entry points x every SIMD tier the host supports.  These are whole-array
// reductions behind a non-inlinable function boundary, which is how the future
// consumers will call them, so the comparison against PART I is conservative
// rather than flattering: PART I's lambda is fully inlined into its loop.
//
// The stated Stage 2 target is >= 2x on the k12 payload relative to the Stage 0
// scalar production K,K',K'' figure, which is the row PART II's "vs prod"
// column divides by.
//
// Build:  make bench       (or: make build/tests/bench_spa_cgf && ./build/tests/bench_spa_cgf)

#include "spa_reference.hpp"

#include "util/spa_cgf.hpp"
#include "util/simd_dispatch.hpp"

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

// ──────────────────────────────────────────────────────────────────────
// Stage 2 harness: whole-array reductions
//
// `fn(t)` performs one complete reduction over all n subjects and returns a
// scalar derived from the result, so the timed quantity is the real kernel's
// per-subject cost including its call overhead and its tail handling.  The t
// sweep matches PART I exactly (0.002 * (rep % 64)) so the two parts drive the
// kernels over the same saddlepoint arguments.
// ──────────────────────────────────────────────────────────────────────
template <class Fn>
Result timeReduction(int n, int repeats, Fn &&fn) {
    double warm = 0.0;
    for (int rep = 0; rep < 16; ++rep) warm += fn(0.01);
    sink(warm);

    const auto t0 = std::chrono::steady_clock::now();
    double acc = 0.0;
    for (int rep = 0; rep < repeats; ++rep)
        acc += fn(0.002 * static_cast<double>(rep % 64));
    const auto t1 = std::chrono::steady_clock::now();
    sink(acc);

    const double ns = std::chrono::duration<double, std::nano>(t1 - t0).count();
    return Result{ns / (static_cast<double>(repeats) * n), acc};
}

const char *tierName(int level) {
    return (level >= 2) ? "avx512" : ((level >= 1) ? "avx2" : "scalar");
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
    std::printf("  outliers = %d, repeats = %d, working set = %.1f KiB\n",
                nOutlier, repeats, nOutlier * sizeof(double) / 1024.0);
    std::printf("  host SIMD tier: %s\n\n", tierName(static_cast<int>(simdLevel())));

    std::printf("PART I — reference forms (Stage 0 baseline, per-subject inlined lambda)\n\n");
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

    // Production K,K',K'' figures, kept per MAF: these are the denominators the
    // Stage 2 target is stated against.
    double prodAtMaf[2] = {0.0, 0.0};

    for (size_t i = 0; i < rows.size(); ++i) {
        if (i % 3 == 0) {
            baseline = rows[i].res.nsPerEval;
            prodAtMaf[i / 3] = baseline;
        }
        std::printf("  %-34s %12.3f %9.2fx\n",
                    rows[i].name.c_str(), rows[i].res.nsPerEval,
                    baseline / rows[i].res.nsPerEval);
        if (i % 3 == 2) std::printf("\n");
    }

    // Print the checksums so no kernel can be optimized away between builds.
    std::printf("  checksums:");
    for (const Row &r : rows) std::printf(" %.6e", r.res.checksum);
    std::printf("\n");

    // ══════════════════════════════════════════════════════════════════
    // PART II — the real kernels, every tier
    // ══════════════════════════════════════════════════════════════════

    const int level = static_cast<int>(simdLevel());

    // Per-individual allele frequencies and hapcounts for the two variants that
    // need them.  Deterministic, same seed discipline as the residuals.
    std::mt19937_64 rng2(20260730u);
    std::uniform_real_distribution<double> afd(1e-4, 0.5);
    std::uniform_int_distribution<int> hpd(0, 2);
    std::vector<double> af(static_cast<size_t>(nOutlier));
    std::vector<double> hap(static_cast<size_t>(nOutlier));
    for (int i = 0; i < nOutlier; ++i) {
        af[i]  = afd(rng2);
        hap[i] = static_cast<double>(hpd(rng2));
    }

    const double *R = resid.data();
    const double *A = af.data();
    const double *H = hap.data();
    const int N = nOutlier;

    struct Row2 { std::string variant, entry, tier; Result res; };
    std::vector<Row2> rows2;

    for (double maf : {0.05, 0.30}) {
        const double prod = prodAtMaf[(maf == 0.05) ? 0 : 1];
        (void)prod;

        // ── binomUniform ──
        rows2.push_back({"binomUniform", "k12", "scalar",
            timeReduction(N, repeats, [&](double t) {
                return spa_cgf::tier::binomUniformK12_scalar(t, R, N, maf).K2; })});
#if defined(__x86_64__) || defined(_M_X64)
        if (level >= 1)
            rows2.push_back({"binomUniform", "k12", "avx2",
                timeReduction(N, repeats, [&](double t) {
                    return spa_cgf::tier::binomUniformK12_avx2(t, R, N, maf).K2; })});
        if (level >= 2)
            rows2.push_back({"binomUniform", "k12", "avx512",
                timeReduction(N, repeats, [&](double t) {
                    return spa_cgf::tier::binomUniformK12_avx512(t, R, N, maf).K2; })});
#endif
        rows2.push_back({"binomUniform", "kFull", "dispatch",
            timeReduction(N, repeats, [&](double t) {
                return spa_cgf::binomUniformKFull(t, R, N, maf).K0; })});

        // ── binomIndiv (af per subject; maf unused) ──
        rows2.push_back({"binomIndiv", "k12", "scalar",
            timeReduction(N, repeats, [&](double t) {
                return spa_cgf::tier::binomIndivK12_scalar(t, R, A, N).K2; })});
#if defined(__x86_64__) || defined(_M_X64)
        if (level >= 1)
            rows2.push_back({"binomIndiv", "k12", "avx2",
                timeReduction(N, repeats, [&](double t) {
                    return spa_cgf::tier::binomIndivK12_avx2(t, R, A, N).K2; })});
        if (level >= 2)
            rows2.push_back({"binomIndiv", "k12", "avx512",
                timeReduction(N, repeats, [&](double t) {
                    return spa_cgf::tier::binomIndivK12_avx512(t, R, A, N).K2; })});
#endif
        rows2.push_back({"binomIndiv", "kFull", "dispatch",
            timeReduction(N, repeats, [&](double t) {
                return spa_cgf::binomIndivKFull(t, R, A, N).K0; })});

        // ── binomHapcount ──
        rows2.push_back({"binomHapcount", "k12", "scalar",
            timeReduction(N, repeats, [&](double t) {
                return spa_cgf::tier::binomHapcountK12_scalar(t, R, H, N, maf).K2; })});
#if defined(__x86_64__) || defined(_M_X64)
        if (level >= 1)
            rows2.push_back({"binomHapcount", "k12", "avx2",
                timeReduction(N, repeats, [&](double t) {
                    return spa_cgf::tier::binomHapcountK12_avx2(t, R, H, N, maf).K2; })});
        if (level >= 2)
            rows2.push_back({"binomHapcount", "k12", "avx512",
                timeReduction(N, repeats, [&](double t) {
                    return spa_cgf::tier::binomHapcountK12_avx512(t, R, H, N, maf).K2; })});
#endif
        rows2.push_back({"binomHapcount", "kFull", "dispatch",
            timeReduction(N, repeats, [&](double t) {
                return spa_cgf::binomHapcountKFull(t, R, H, N, maf).K0; })});

        // Emit this MAF's block, then move on.
        std::printf("\nPART II — real kernels (Stage 2), whole-array reductions, "
                    "MAF/q %.2f\n", maf);
        std::printf("  denominator: PART I production K,K',K'' = %.3f ns/subject\n\n",
                    prod);
        std::printf("  %-14s %-6s %-9s %12s %10s\n",
                    "variant", "entry", "tier", "ns/subject", "vs prod");
        std::printf("  %-14s %-6s %-9s %12s %10s\n",
                    "--------------", "------", "---------",
                    "------------", "----------");
        for (const Row2 &r : rows2)
            std::printf("  %-14s %-6s %-9s %12.3f %9.2fx\n",
                        r.variant.c_str(), r.entry.c_str(), r.tier.c_str(),
                        r.res.nsPerEval, prod / r.res.nsPerEval);

        std::printf("\n  checksums:");
        for (const Row2 &r : rows2) std::printf(" %.6e", r.res.checksum);
        std::printf("\n");

        rows2.clear();
    }

    return 0;
}
