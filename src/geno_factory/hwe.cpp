// hwe.cpp — Hardy-Weinberg equilibrium exact test and QC stats
//
// Exact test (SNPHWE2): Wigginton JE, Cutler DJ, Abecasis GR (2005).
// Am J Hum Genet 76:887-893.
//
// This is the plink2 --hardy default method.  The evaluation is plink2's
// HweLnP, which returns ln(p) and is therefore the only form of the test that
// stays meaningful once p drops below the smallest normal double.  GRAB's own
// linear-domain transcription of the same recurrence (HweExact) has been
// deleted: it agreed with HweLnP to 1e-8 wherever it did not underflow, but it
// returned exactly 0 below p ~ 1e-300 and cost O(het_count) per marker where
// HweLnP costs ~O(1).

#include "geno_factory/hwe.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

namespace plink2 {
// Declared here rather than by including plink2_stats.h, which pulls in the
// whole pgenlib common header.  Definition: third_party/plink2-a.6.33/
// include/plink2_stats.cc, compiled into build/pgenlib/plink2_stats.o.
double HweLnP(
    int32_t obs_hets,
    int32_t obs_hom1,
    int32_t obs_hom2,
    uint32_t midp
);

} // namespace plink2

namespace {
// ln(10).  The change of base from HweLnP's natural log to the -log10 the
// LOG10P_HWE column carries.
constexpr double kLn10 = 2.30258509299404568402;

} // namespace

double hweNegLog10P(
    uint32_t obs_hets,
    uint32_t obs_hom1,
    uint32_t obs_hom2
) {
    // midp = 0 selects the non-mid-p exact test, which is the semantics the
    // deleted HweExact implemented and which plink2 --hardy uses by default.
    const double lnp = plink2::HweLnP(
        static_cast<int32_t>(obs_hets),
        static_cast<int32_t>(obs_hom1),
        static_cast<int32_t>(obs_hom2),
        /*midp=*/0
    );
    // HweLnP returns ln(p) <= 0 and returns exactly 0 for p = 1 (including the
    // rare_ct < 2 "no test" case).  Negating that would give -0.0, which the
    // output formatter prints as "-0"; return +0.0 so a marker in exact
    // Hardy-Weinberg proportions reads as 0.  This is a sign-of-zero
    // normalization, not a clamp: no non-zero value is altered, and a genuinely
    // positive ln p (which HweLnP does not produce) would still surface as a
    // negative LOG10P_HWE and fail invariant check C2.
    if (lnp == 0.0) return 0.0;
    return -lnp / kLn10;
}

GenoStats statsFromCounts(
    uint32_t nHomAlt,
    uint32_t nHet,
    uint32_t nHomRef,
    uint32_t nMissing,
    uint32_t nSamples
) {
    const uint32_t nonMissing = nSamples - nMissing;
    GenoStats gs;
    gs.altCounts = 2 * nHomAlt + nHet; // count of ALT allele
    gs.missingRate = static_cast<double>(nMissing) / nSamples;
    gs.fromDosage = false;             // this overload is the hard-call path

    if (nonMissing == 0) {
        gs.altFreq = std::numeric_limits<double>::quiet_NaN();
        gs.log10pHwe = std::numeric_limits<double>::quiet_NaN();
        gs.maf = std::numeric_limits<double>::quiet_NaN();
        gs.mac = 0;
        return gs;
    }

    const AlleleFreq af = alleleFreqFromTotal(gs.altCounts, nonMissing);
    gs.altFreq = af.altFreq;
    gs.maf = af.maf;
    gs.mac = af.mac;
    gs.log10pHwe = hweNegLog10P(nHet, nHomAlt, nHomRef);
    return gs;
}
