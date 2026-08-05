// hwe.hpp — Hardy-Weinberg equilibrium exact test and QC stats
//
// Single canonical implementation used by all genotype backends
// (PLINK, PGEN, VCF, BGEN) and the marker engine.
//
// Method: exact test (SNPHWE2), matching plink2 --hardy default.
// Wigginton JE, Cutler DJ, Abecasis GR (2005). Am J Hum Genet 76:887-893.
//
// The test is reported as -log10(p), not p (decision D7 of
// dev-notes/methods/log10p_unify/00_decisions.md).  The evaluation is
// delegated to plink2's HweLnP, which is already linked into the binary via
// build/pgenlib/plink2_stats.o; GRAB's own linear-domain copy of the same
// algorithm has been deleted.  The linear copy was correct but underflowed to
// exactly 0 for -log10 p >~ 300 (871 of 995 such markers on the panel measured
// in 01_numerics.md §5.2) and was O(het_count) rather than the ~O(1) of the
// plink2 version.
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

// Genotype class counts + derived QC stats.
struct GenoStats {
    double altFreq;
    uint32_t altCounts;
    double missingRate;
    double log10pHwe;
    double maf;
    uint32_t mac;
};

// plink2 --hard-call-threshold default.  A dosage is called to the nearest
// hard-call (0/1/2) when it lies within `thr` of that integer, and is
// otherwise treated as missing for HWE only.  This mirrors plink2, which
// derives every pgen's hard-call track from dosages once, at import time,
// via ApplyHardCallThresh, and then counts that hard-call track for --hardy.
inline constexpr double kHardCallThreshold = 0.1;

// Classify a dosage into a hard-call genotype for HWE.  Returns 0/1/2 when
// |d - round(d)| <= thr, or -1 (HWE-missing) otherwise.  NaN / out-of-range
// dosages return -1.  Exact integers 0.0/1.0/2.0 always pass (|d-r| == 0),
// so pure hard-call input reproduces the prior counting bit-for-bit.
inline int dosageHardcall(double d, double thr) {
    if (!(d >= 0.0 && d <= 2.0)) return -1;
    const double r = (d < 0.5) ? 0.0 : (d < 1.5 ? 1.0 : 2.0);
    return (std::fabs(d - r) <= thr) ? static_cast<int>(r) : -1;
}

// -log10 of the exact HWE p-value (plink2 --hardy, non-mid-p variant).
//
// The sole call into plink2::HweLnP in the tree: every producer of the
// LOG10P_HWE column goes through this function, so the ln -> -log10 change of
// base exists once rather than at each of the six count-to-stat sites.
//
// The returned value is >= 0 and is +0.0 (never -0.0) for a marker in exact
// agreement with Hardy-Weinberg proportions.  It is never +Inf: HweLnP works
// in the log domain throughout, so a p of 1e-30000 is returned as 30000, not
// as an underflowed 0.  Counts with fewer than two rare alleles have no test
// and return 0, which is the ln p = 0 that HweLnP itself reports for them and
// matches the deleted linear implementation's n == 0 -> p = 1.
double hweNegLog10P(
    uint32_t obs_hets,
    uint32_t obs_hom1,
    uint32_t obs_hom2
);

// Compute QC stats from genotype class counts.  Always uses the exact
// HWE test (plink2 default).  The first argument is the count of
// subjects homozygous for the ALT allele (the second allele in
// .pvar/.bim/.bgen/.vcf); altCounts is computed as 2*nHomAlt + nHet so
// that altFreq is the ALT allele frequency, not the REF allele
// frequency.  Callers must classify per-subject dosage 0/1/2 into
// nHomRef/nHet/nHomAlt before invoking this function.
GenoStats statsFromCounts(
    uint32_t nHomAlt,
    uint32_t nHet,
    uint32_t nHomRef,
    uint32_t nMissing,
    uint32_t nSamples
);
