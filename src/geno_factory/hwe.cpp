// hwe.cpp — Hardy-Weinberg equilibrium exact test and QC stats
//
// Exact test (SNPHWE2): Wigginton JE, Cutler DJ, Abecasis GR (2005).
// Am J Hum Genet 76:887-893.
//
// This is the plink2 --hardy default method.  Always used for HWE
// computation — the chi-squared approximation has been removed because
// the exact test is fast enough for all practical sample sizes and
// avoids the information loss inherent in approximate methods.

#include "geno_factory/hwe.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

double HweExact(
    uint32_t obs_hets,
    uint32_t obs_hom1,
    uint32_t obs_hom2
) {
    const int64_t obs_homc = std::max(obs_hom1, obs_hom2);
    const int64_t obs_homr = std::min(obs_hom1, obs_hom2);
    const int64_t rare = 2 * obs_homr + static_cast<int64_t>(obs_hets);
    const int64_t n = static_cast<int64_t>(obs_hets) + obs_homc + obs_homr;
    const int64_t obs = static_cast<int64_t>(obs_hets);

    if (n == 0) return 1.0;

    // Mode of het-count distribution under HWE
    int64_t mid = (rare * (2 * n - rare)) / (2 * n);
    if ((rare & 1) ^ (mid & 1)) ++mid;

    // Adjust mid to sit at the actual peak
    {
        int64_t hr = (rare - mid) / 2;
        int64_t hc = n - mid - hr;
        if (mid + 2 <= rare && hr > 0 && 4.0 * hr * hc > (mid + 2.0) * (mid + 1.0)) {
            mid += 2;
        } else if (mid >= 2) {
            if (static_cast<double>(mid) * (mid - 1) > 4.0 * (hr + 1.0) * (hc + 1.0)) {
                mid -= 2;
            }
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
                sum += prob;
                ++cr;
                ++cc;
            }
            thresh = prob;
            p = thresh;
            for (int64_t h = obs; h >= 2; h -= 2) {
                prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
                sum += prob;
                p += prob;
                ++cr;
                ++cc;
            }
        }
        {
            double prob = 1.0;
            int64_t cr = mid_homr, cc = mid_homc;
            for (int64_t h = mid; h <= rare - 2; h += 2) {
                prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
                sum += prob;
                if (prob <= thresh) p += prob;
                --cr;
                --cc;
            }
        }
    } else {
        {
            double prob = 1.0;
            int64_t cr = mid_homr, cc = mid_homc;
            for (int64_t h = mid; h < obs; h += 2) {
                prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
                sum += prob;
                --cr;
                --cc;
            }
            thresh = prob;
            p = thresh;
            for (int64_t h = obs; h <= rare - 2; h += 2) {
                prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
                sum += prob;
                p += prob;
                --cr;
                --cc;
            }
        }
        {
            double prob = 1.0;
            int64_t cr = mid_homr, cc = mid_homc;
            for (int64_t h = mid; h >= 2; h -= 2) {
                prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
                sum += prob;
                if (prob <= thresh) p += prob;
                ++cr;
                ++cc;
            }
        }
    }

    return std::min(p / sum, 1.0);
}

GenoStats statsFromCounts(
    uint32_t nHom1,
    uint32_t nHet,
    uint32_t nHom2,
    uint32_t nMissing,
    uint32_t nSamples
) {
    const uint32_t nonMissing = nSamples - nMissing;
    GenoStats gs;
    gs.altCounts = 2 * nHom1 + nHet; // count A1 (ALT)
    gs.missingRate = static_cast<double>(nMissing) / nSamples;

    if (nonMissing == 0) {
        gs.altFreq = std::numeric_limits<double>::quiet_NaN();
        gs.hweP = std::numeric_limits<double>::quiet_NaN();
        gs.maf = std::numeric_limits<double>::quiet_NaN();
        gs.mac = 0;
        return gs;
    }

    gs.altFreq = static_cast<double>(gs.altCounts) / (2.0 * nonMissing);
    gs.maf = std::min(gs.altFreq, 1.0 - gs.altFreq);
    gs.mac = std::min(gs.altCounts, 2 * nonMissing - gs.altCounts);
    gs.hweP = HweExact(nHet, nHom1, nHom2);
    return gs;
}
