// spacox_cgf_test.cpp — SPACox's empirical CGF: bucketing (P1) and the
// parallel, overflow-free table build (P5).
//
// src/spacox/spacox_cgf.hpp is header-only precisely so that this test needs
// no link against the marker engine; it includes Eigen and util/spa.hpp and
// nothing else from src/.
//
// What is asserted, and why each assertion exists:
//
//   1. BUCKETED == DENSE.  P1 replaces the per-subject walk over the non-zero
//      genotypes with a histogram over the at most four distinct weights.
//      The two must agree.  They are not required to agree bit-for-bit — the
//      bucketed form adds `count * (w * K')` once where the dense form adds
//      `w * K'` count times, which is a different association order and is in
//      fact the *more* accurate of the two — so the comparison is against a
//      long-double reference as well as against each other.
//
//   2. DOSAGES SURVIVE.  BGEN and VCF input, and any reader below the
//      hard-call threshold, produce genotypes outside {0, 1, 2, 2*altFreq}.
//      Those must land in the leftover list and still be evaluated exactly,
//      so the same equality is asserted on a genotype vector that is half
//      hard-called and half continuous, and on one that is entirely
//      continuous.
//
//   3. THE IMPUTED BUCKET.  A missing genotype is mean-imputed by the engine
//      as exactly 2.0 * altFreq, so its weight is exactly zero and it belongs
//      in a bucket rather than in the leftover list.  Asserted directly on
//      bucket occupancy, because getting this wrong is silent: the answer
//      stays correct and only the speed is lost.
//
//   4. THREADING IS BIT-EXACT.  P5 parallelizes the grid loop.  Each grid
//      point is computed from xGrid[i] and the residual vector alone, so the
//      table must be bit-identical for any thread count.  A parallelization
//      that changed the numbers would be a different routine, not a faster
//      one.
//
//   5. NO NON-FINITE TABLE ENTRY, AND THE ENDPOINT VALUE IS RIGHT.  The old
//      accumulation formed (M0*M2 - M1*M1)/(M0*M0); with max|r| = 4.14 the
//      grid endpoint has exponent 414, so M0 = 1.6e176 and both products
//      overflow, giving NaN — which computeSlopes then propagates across the
//      entire first bin, a span of 50 in t.  The replacement is checked
//      against a long-double reference over the whole grid, and the old form
//      is recomputed here so the test also demonstrates the defect rather
//      than merely asserting its absence.

#include <cmath>
#include <cstdint>
#include <random>
#include <vector>

#include <Eigen/Dense>

#include "spacox/spacox_cgf.hpp"
#include "tinytest.hpp"

namespace {

// ──────────────────────────────────────────────────────────────────────
// References
// ──────────────────────────────────────────────────────────────────────

// K, K', K'' of the empirical residual CGF at one abscissa, in long double
// with the log-sum-exp shift and a two-pass variance, i.e. as accurately as
// this can be done without arbitrary precision.
struct RefTriple {
    long double k0, k1, k2;
};

RefTriple refEmpiricalCgf(const Eigen::VectorXd &r, double t) {
    const int n = static_cast<int>(r.size());
    long double m = -std::numeric_limits<long double>::infinity();
    for (int j = 0; j < n; ++j) {
        const long double e = static_cast<long double>(r[j]) * static_cast<long double>(t);
        if (e > m) m = e;
    }
    long double W = 0.0L, M1 = 0.0L;
    std::vector<long double> u(static_cast<size_t>(n));
    for (int j = 0; j < n; ++j) {
        u[static_cast<size_t>(j)] =
            std::exp(static_cast<long double>(r[j]) * static_cast<long double>(t) - m);
        W += u[static_cast<size_t>(j)];
        M1 += u[static_cast<size_t>(j)] * static_cast<long double>(r[j]);
    }
    const long double mean = M1 / W;
    long double var = 0.0L;
    for (int j = 0; j < n; ++j) {
        const long double d = static_cast<long double>(r[j]) - mean;
        var += u[static_cast<size_t>(j)] * d * d;
    }
    var /= W;
    return RefTriple{m + std::log(W / static_cast<long double>(n)), mean, var};
}

// The pre-Stage-3 accumulation, kept verbatim so the defect is demonstrated
// and not merely asserted away.
struct OldTriple {
    double k0, k1, k2;
};

OldTriple oldMomentForm(const Eigen::VectorXd &r, double t) {
    const int n = static_cast<int>(r.size());
    double M0 = 0.0, M1 = 0.0, M2 = 0.0;
    for (int j = 0; j < n; ++j) {
        const double e = std::exp(r[j] * t);
        M0 += e;
        const double re = r[j] * e;
        M1 += re;
        M2 += r[j] * re;
    }
    const double invN = 1.0 / static_cast<double>(n);
    M0 *= invN;
    M1 *= invN;
    M2 *= invN;
    return OldTriple{std::log(M0), M1 / M0,
                     (M0 * M2 - M1 * M1) / (M0 * M0)};
}

// ──────────────────────────────────────────────────────────────────────
// Fixtures
// ──────────────────────────────────────────────────────────────────────

// Residual vectors with the shape SPACox actually sees.  `heavyTail`
// reproduces the property that triggered the overflow: max|r| above 4, which
// the Cox martingale residuals of the bundled Time_Event fixture reach
// (max|r| = 4.1357).
Eigen::VectorXd makeResiduals(int n, bool heavyTail, uint64_t seed) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> nd(0.0, 1.0);
    Eigen::VectorXd r(n);
    for (int i = 0; i < n; ++i) r[i] = nd(rng);
    if (heavyTail) {
        // Cox martingale residuals are bounded above by 1 and have a long
        // negative tail; mimic the *magnitude* rather than the exact law.
        r = (r.array() * 0.5).matrix();
        r[0] = -4.1357;
        r[1] = 0.9;
    }
    r.array() -= r.mean();  // null-model residuals sum to zero
    return r;
}

// Genotype vector with a controllable mix of hard calls, mean-imputed
// subjects and genuine dosages.
Eigen::VectorXd makeGenotypes(
    int n,
    double altFreq,
    int nImputed,
    int nDosage,
    uint64_t seed
) {
    std::mt19937_64 rng(seed);
    std::binomial_distribution<int> bd(2, altFreq);
    std::uniform_real_distribution<double> ud(0.0, 2.0);
    Eigen::VectorXd g(n);
    for (int i = 0; i < n; ++i) g[i] = static_cast<double>(bd(rng));
    for (int i = 0; i < nImputed && i < n; ++i) g[i] = 2.0 * altFreq;
    for (int i = 0; i < nDosage && (nImputed + i) < n; ++i)
        g[nImputed + i] = ud(rng);
    return g;
}

std::vector<double> weightVector(const Eigen::VectorXd &g, double twoMAF, double sd) {
    std::vector<double> w(static_cast<size_t>(g.size()));
    for (int i = 0; i < g.size(); ++i)
        w[static_cast<size_t>(i)] = (g[i] - twoMAF) / sd;
    return w;
}

}  // namespace

// ══════════════════════════════════════════════════════════════════════
// 1-3.  Bucketing
// ══════════════════════════════════════════════════════════════════════

TEST(bucketed_matches_dense_hardcalls) {
    const int N = 4000;
    const Eigen::VectorXd resid = makeResiduals(N, false, 20260730u);
    const CumulantTable ct = spacox_cgf::buildCumulantTable(resid, 4);
    const spacox_cgf::TableView tv(ct);

    for (double af : {0.01, 0.1, 0.3, 0.5, 0.85}) {
        const Eigen::VectorXd g = makeGenotypes(N, af, 0, 0, 99u);
        const double twoMAF = 2.0 * af;
        double ss = 0.0;
        for (int i = 0; i < N; ++i) ss += (g[i] - twoMAF) * (g[i] - twoMAF);
        const double sd = std::sqrt(ss);

        spacox_cgf::GenoWeights gw;
        spacox_cgf::buildGenoWeights(g.data(), N, twoMAF, sd, gw);
        // Pure hard calls with altFreq not in {0, 0.5, 1}: three buckets at
        // most, and never a leftover.
        CHECK(gw.extra.empty());
        CHECK(gw.nBucket >= 1 && gw.nBucket <= 4);

        const std::vector<double> w = weightVector(g, twoMAF, sd);
        for (double t : {-8.0, -2.5, -0.3, 0.0, 0.3, 2.5, 8.0}) {
            const spa::K012 b = spacox_cgf::evalBucketed(tv, t, gw, 1.5);
            const spa::K012 d =
                spacox_cgf::evalDense(tv, t, w.data(), N, 1.5);
            CHECK_REL(b.k0, d.k0, 1e-12);
            CHECK_REL(b.k1, d.k1, 1e-11);
            CHECK_REL(b.k2, d.k2, 1e-12);
        }
    }
}

TEST(bucketed_matches_dense_with_dosages) {
    const int N = 3000;
    const Eigen::VectorXd resid = makeResiduals(N, true, 20260731u);
    const CumulantTable ct = spacox_cgf::buildCumulantTable(resid, 8);
    const spacox_cgf::TableView tv(ct);
    const double af = 0.22;
    const double twoMAF = 2.0 * af;

    // Three regimes: no dosage, half dosage, all dosage.
    for (int nDos : {0, N / 2, N}) {
        const Eigen::VectorXd g = makeGenotypes(N, af, 17, nDos, 5150u);
        double ss = 0.0;
        for (int i = 0; i < N; ++i) ss += (g[i] - twoMAF) * (g[i] - twoMAF);
        const double sd = std::sqrt(ss);

        spacox_cgf::GenoWeights gw;
        spacox_cgf::buildGenoWeights(g.data(), N, twoMAF, sd, gw);
        const std::vector<double> w = weightVector(g, twoMAF, sd);

        // The classification must be exhaustive: every subject is either in a
        // bucket or in the leftover list, exactly once.
        double total = 0.0;
        for (int b = 0; b < gw.nBucket; ++b) total += gw.bc[b];
        CHECK_MSG(total + static_cast<double>(gw.extra.size()) ==
                      static_cast<double>(N),
                  "bucket multiplicities plus leftovers must total N");

        for (double t : {-6.0, -1.0, 0.7, 4.0}) {
            const spa::K012 b = spacox_cgf::evalBucketed(tv, t, gw, -2.0);
            const spa::K012 d = spacox_cgf::evalDense(tv, t, w.data(), N, -2.0);
            CHECK_REL(b.k0, d.k0, 1e-12);
            CHECK_REL(b.k1, d.k1, 1e-11);
            CHECK_REL(b.k2, d.k2, 1e-12);
        }
    }
}

TEST(imputed_subjects_form_a_zero_weight_bucket) {
    const int N = 500;
    const double af = 0.137;
    const double twoMAF = 2.0 * af;
    const int nImp = 23;
    const Eigen::VectorXd g = makeGenotypes(N, af, nImp, 0, 7u);
    double ss = 0.0;
    for (int i = 0; i < N; ++i) ss += (g[i] - twoMAF) * (g[i] - twoMAF);
    const double sd = std::sqrt(ss);

    spacox_cgf::GenoWeights gw;
    spacox_cgf::buildGenoWeights(g.data(), N, twoMAF, sd, gw);

    // The mean-imputed value is written by the engine as exactly 2.0*altFreq
    // and twoMAF is the same expression, so the equality test is exact and
    // these subjects must NOT fall through to the leftover list.
    CHECK_MSG(gw.extra.empty(), "mean-imputed subjects must be bucketed");
    bool sawZero = false;
    for (int b = 0; b < gw.nBucket; ++b) {
        if (gw.bw[b] == 0.0) {
            sawZero = true;
            CHECK(gw.bc[b] == static_cast<double>(nImp));
        }
    }
    CHECK_MSG(sawZero, "a zero-weight bucket must exist");
}

// ══════════════════════════════════════════════════════════════════════
// 4.  Threading
// ══════════════════════════════════════════════════════════════════════

TEST(table_is_bit_identical_across_thread_counts) {
    const Eigen::VectorXd resid = makeResiduals(2503, true, 20260732u);
    const CumulantTable a = spacox_cgf::buildCumulantTable(resid, 1);
    for (int nthr : {2, 7, 64, 4096}) {
        const CumulantTable b = spacox_cgf::buildCumulantTable(resid, nthr);
        REQUIRE(a.nGrid == b.nGrid);
        bool same = true;
        for (int i = 0; i < a.nGrid; ++i) {
            if (a.yK0[i] != b.yK0[i] || a.yK1[i] != b.yK1[i] ||
                a.yK2[i] != b.yK2[i] || a.xGrid[i] != b.xGrid[i])
                same = false;
        }
        CHECK_MSG(same, "parallel table build must be bit-identical to serial");
    }
}

// ══════════════════════════════════════════════════════════════════════
// 5.  Overflow, cancellation and accuracy of the table
// ══════════════════════════════════════════════════════════════════════

TEST(old_moment_form_overflows_at_the_grid_endpoint) {
    // Characterization, not a requirement: the point of this test is that the
    // defect was real.  With max|r| = 4.1357 the endpoint exponent is 413.5,
    // so M0 = 1.6e176 and M0*M2, M1*M1 both exceed DBL_MAX.
    const Eigen::VectorXd resid = makeResiduals(2503, true, 20260732u);
    const CumulantTable ct = spacox_cgf::buildCumulantTable(resid, 1);
    const double tEnd = ct.xGrid[0];
    const OldTriple old = oldMomentForm(resid, tEnd);
    CHECK_MSG(std::isnan(old.k2),
              "the pre-Stage-3 product form is expected to produce NaN here");
    CHECK_MSG(std::isfinite(ct.yK2[0]) && ct.yK2[0] >= 0.0,
              "the replacement must be finite and non-negative");

    // And the replacement must agree with the long-double reference, which is
    // where the old form's 2.9e-2 relative error at t = -50 also shows up.
    const RefTriple ref = refEmpiricalCgf(resid, tEnd);
    CHECK_REL(ct.yK2[0], static_cast<double>(ref.k2), 1e-9);
    CHECK_REL(ct.yK0[0], static_cast<double>(ref.k0), 1e-14);
    CHECK_REL(ct.yK1[0], static_cast<double>(ref.k1), 1e-14);
}

TEST(table_entries_are_all_finite_and_curvature_non_negative) {
    for (bool heavy : {false, true}) {
        const Eigen::VectorXd resid = makeResiduals(1777, heavy, 424242u);
        const CumulantTable ct = spacox_cgf::buildCumulantTable(resid, 16);
        int bad = 0;
        for (int i = 0; i < ct.nGrid; ++i) {
            if (!std::isfinite(ct.yK0[i]) || !std::isfinite(ct.yK1[i]) ||
                !std::isfinite(ct.yK2[i]) || ct.yK2[i] < 0.0)
                ++bad;
        }
        for (int i = 0; i < ct.nGrid - 1; ++i) {
            if (!std::isfinite(ct.slopeK0[i]) || !std::isfinite(ct.slopeK1[i]) ||
                !std::isfinite(ct.slopeK2[i]))
                ++bad;
        }
        CHECK_MSG(bad == 0, "no table entry or slope may be non-finite");
    }
}

TEST(table_accuracy_against_long_double_reference) {
    const Eigen::VectorXd resid = makeResiduals(1777, true, 424242u);
    const CumulantTable ct = spacox_cgf::buildCumulantTable(resid, 16);

    // Sample the grid, including both endpoints and the near-zero centre.
    //
    // K and K' are compared on a mixed relative/absolute criterion, and the
    // old form is measured alongside on the same criterion so the comparison
    // is like for like.  A purely relative test is the wrong instrument for
    // either quantity near t = 0: both pass through zero there (K(0) = 0 and
    // K'(t*) = 0 at the t* where the tilted mean crosses the residual mean, of
    // which there is one because the residuals sum to zero), and no algorithm
    // achieves a small *relative* error at a zero crossing.  What the
    // saddleprint needs is a small absolute error, since K enters
    // w = sqrt(2*(zeta*s - K)) additively against zeta*s = O(z^2) >= 4.
    //
    // K'' is genuinely multiplicative — it appears as v = zeta*sqrt(K'') — so
    // it is held to a relative bound, and it is exactly the quantity whose old
    // form lost every digit at large |t|.
    double worstK0 = 0.0, worstK1 = 0.0, worstK2 = 0.0;
    double oldWorstK0 = 0.0, oldWorstK1 = 0.0, oldWorstK2 = 0.0;
    int oldNonFinite = 0;
    auto mixed = [](double got, double ref, double scale) {
        if (std::isnan(got) || std::isnan(ref))
            return std::numeric_limits<double>::infinity();
        return std::fabs(got - ref) / std::fmax(std::fabs(ref), scale);
    };
    for (int i = 0; i < ct.nGrid; i += 37) {
        const double t = ct.xGrid[i];
        const RefTriple ref = refEmpiricalCgf(resid, t);
        const double r0 = static_cast<double>(ref.k0);
        const double r1 = static_cast<double>(ref.k1);
        const double r2 = static_cast<double>(ref.k2);
        worstK0 = std::max(worstK0, mixed(ct.yK0[i], r0, 1.0));
        worstK1 = std::max(worstK1, mixed(ct.yK1[i], r1, 1.0));
        worstK2 = std::max(worstK2, tinytest::relDiff(ct.yK2[i], r2));

        const OldTriple old = oldMomentForm(resid, t);
        if (!std::isfinite(old.k0) || !std::isfinite(old.k1) ||
            !std::isfinite(old.k2))
            ++oldNonFinite;
        oldWorstK0 = std::max(oldWorstK0, mixed(old.k0, r0, 1.0));
        oldWorstK1 = std::max(oldWorstK1, mixed(old.k1, r1, 1.0));
        oldWorstK2 = std::max(oldWorstK2, tinytest::relDiff(old.k2, r2));
    }
    std::printf("    worst error vs long double, new:  K %.3e  K' %.3e  K'' %.3e\n",
                worstK0, worstK1, worstK2);
    std::printf("    worst error vs long double, old:  K %.3e  K' %.3e  K'' %.3e"
                "   (%d non-finite samples)\n",
                oldWorstK0, oldWorstK1, oldWorstK2, oldNonFinite);
    CHECK(worstK0 <= 1e-14);
    CHECK(worstK1 <= 1e-14);
    CHECK(worstK2 <= 1e-11);
    // The replacement must not be worse than what it replaced, anywhere.
    CHECK(worstK0 <= oldWorstK0);
    CHECK(worstK2 <= oldWorstK2);
    CHECK_MSG(oldNonFinite > 0,
              "characterization: the old form is expected to go non-finite on "
              "this residual vector");
}

TINYTEST_MAIN
