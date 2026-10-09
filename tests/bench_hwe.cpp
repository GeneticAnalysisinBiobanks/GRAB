// bench_hwe.cpp — cost of the Hardy-Weinberg exact test per marker, for the
// implementation this repository used to carry and the one it calls now.
//
// Stage 9 of dev-notes/methods/log10p_unify (decision D7) replaced GRAB's own
// linear-domain SNPHWE2 transcription with plink2::HweLnP.  The two differ in
// asymptotic cost, not only in constant: the deleted version accumulated the
// whole tail term by term and is O(het_count), while HweLnP stops as soon as
// an increment falls below the 53-bit precision of the partial sum and jumps
// across the mode analytically, which is close to O(1) in the sample size.
// The benchmark exists so that claim is measured rather than asserted, and so
// that a future change to either side is caught by a number.
//
// The HWE test runs once per marker per QC group in every method GRAB ships,
// so at biobank scale the difference is not incidental: at n = 500 000 the old
// form costs about 0.1 ms per marker, which is roughly 20 CPU-minutes over ten
// million markers, entirely inside the QC path.
//
// Build:  make bench   (or: make build/tests/bench_hwe && ./build/tests/bench_hwe)

#include "geno_factory/hwe.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

namespace {

// The linear-domain exact test deleted from src/geno_factory/hwe.cpp by
// Stage 9, transcribed verbatim as the "before" side of the measurement.
double hweExactLinear(uint32_t obs_hets, uint32_t obs_hom1, uint32_t obs_hom2) {
    const int64_t obs_homc = std::max(obs_hom1, obs_hom2);
    const int64_t obs_homr = std::min(obs_hom1, obs_hom2);
    const int64_t rare = 2 * obs_homr + static_cast<int64_t>(obs_hets);
    const int64_t n = static_cast<int64_t>(obs_hets) + obs_homc + obs_homr;
    const int64_t obs = static_cast<int64_t>(obs_hets);
    if (n == 0) return 1.0;
    int64_t mid = (rare * (2 * n - rare)) / (2 * n);
    if ((rare & 1) ^ (mid & 1)) ++mid;
    {
        int64_t hr = (rare - mid) / 2;
        int64_t hc = n - mid - hr;
        if (mid + 2 <= rare && hr > 0 && 4.0 * hr * hc > (mid + 2.0) * (mid + 1.0)) {
            mid += 2;
        } else if (mid >= 2) {
            if (static_cast<double>(mid) * (mid - 1) > 4.0 * (hr + 1.0) * (hc + 1.0)) mid -= 2;
        }
    }
    const int64_t mid_homr = (rare - mid) / 2;
    const int64_t mid_homc = n - mid - mid_homr;
    double sum = 1.0, p = 0.0, thresh;
    if (obs <= mid) {
        {
            double prob = 1.0;
            int64_t cr = mid_homr, cc = mid_homc;
            for (int64_t h = mid; h > obs; h -= 2) {
                prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
                sum += prob; ++cr; ++cc;
            }
            thresh = prob; p = thresh;
            for (int64_t h = obs; h >= 2; h -= 2) {
                prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
                sum += prob; p += prob; ++cr; ++cc;
            }
        }
        {
            double prob = 1.0;
            int64_t cr = mid_homr, cc = mid_homc;
            for (int64_t h = mid; h <= rare - 2; h += 2) {
                prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
                sum += prob;
                if (prob <= thresh) p += prob;
                --cr; --cc;
            }
        }
    } else {
        {
            double prob = 1.0;
            int64_t cr = mid_homr, cc = mid_homc;
            for (int64_t h = mid; h < obs; h += 2) {
                prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
                sum += prob; --cr; --cc;
            }
            thresh = prob; p = thresh;
            for (int64_t h = obs; h <= rare - 2; h += 2) {
                prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
                sum += prob; p += prob; --cr; --cc;
            }
        }
        {
            double prob = 1.0;
            int64_t cr = mid_homr, cc = mid_homc;
            for (int64_t h = mid; h >= 2; h -= 2) {
                prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
                sum += prob;
                if (prob <= thresh) p += prob;
                ++cr; ++cc;
            }
        }
    }
    return std::min(p / sum, 1.0);
}

struct Counts { uint32_t het, homr, homc; };

// A per-n panel spanning the allele frequencies and the degrees of
// HWE departure a real cohort presents.  The het count is the driver of the
// deleted implementation's cost, so a panel dominated by rare variants would
// understate the difference and one dominated by common variants would
// overstate it; both are represented.
std::vector<Counts> panel(uint32_t n) {
    std::vector<Counts> out;
    const double mafs[] = {0.5, 0.3, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.001};
    const double fs[] = {1.0, 0.95, 0.9, 0.8, 0.6, 0.4, 0.2, 0.05,
                         1.05, 1.1, 1.2, 1.4, 1.6, 1.8};
    for (double q : mafs) {
        const double expHet = 2.0 * q * (1.0 - q) * n;
        const double rare = std::floor(2.0 * q * n + 0.5);
        for (double f : fs) {
            double het = std::min(std::floor(expHet * f + 0.5), rare);
            double homr = std::floor((rare - het) / 2.0);
            double homc = n - het - homr;
            if (homr < 0 || homc < 0) continue;
            out.push_back({static_cast<uint32_t>(het),
                           static_cast<uint32_t>(homr),
                           static_cast<uint32_t>(homc)});
        }
    }
    return out;
}

} // namespace

int main() {
    std::printf("Hardy-Weinberg exact test, nanoseconds per marker\n");
    std::printf("%10s %8s %16s %16s %10s\n",
                "n", "markers", "deleted linear", "plink2 HweLnP", "speedup");
    double sink = 0.0;
    for (uint32_t n : {500u, 5000u, 50000u, 500000u}) {
        const std::vector<Counts> p = panel(n);
        const int reps = (n >= 500000u) ? 3 : 20;

        // Warm the panel into cache; both sides then see the same state.
        for (const Counts &c : p) sink += hweNegLog10P(c.het, c.homr, c.homc);

        auto t0 = std::chrono::steady_clock::now();
        for (int r = 0; r < reps; ++r)
            for (const Counts &c : p) sink += hweExactLinear(c.het, c.homr, c.homc);
        auto t1 = std::chrono::steady_clock::now();
        for (int r = 0; r < reps; ++r)
            for (const Counts &c : p) sink += hweNegLog10P(c.het, c.homr, c.homc);
        auto t2 = std::chrono::steady_clock::now();

        const double denom = static_cast<double>(reps) * static_cast<double>(p.size());
        const double nsOld = std::chrono::duration<double, std::nano>(t1 - t0).count() / denom;
        const double nsNew = std::chrono::duration<double, std::nano>(t2 - t1).count() / denom;
        std::printf("%10u %8zu %16.1f %16.1f %9.1fx\n",
                    n, p.size(), nsOld, nsNew, nsOld / nsNew);
    }
    // Keep the loops from being optimized away without printing noise.
    if (!std::isfinite(sink)) std::printf("(sink %g)\n", sink);
    return 0;
}
