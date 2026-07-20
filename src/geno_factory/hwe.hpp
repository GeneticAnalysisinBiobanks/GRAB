// hwe.hpp — Hardy-Weinberg equilibrium exact test and QC stats
//
// Single canonical implementation used by all genotype backends
// (PLINK, PGEN, VCF, BGEN) and the marker engine.
//
// Default method: exact test (SNPHWE2), matching plink2 --hardy default.
// Wigginton JE, Cutler DJ, Abecasis GR (2005). Am J Hum Genet 76:887-893.
#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

// Genotype class counts + derived QC stats.
struct GenoStats {
    double altFreq;
    uint32_t altCounts;
    double missingRate;
    double hweP;
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

// Exact HWE test (SNPHWE2).  O(het_count) time, O(1) auxiliary memory.
// This is the plink2 --hardy default method.
// Returns p-value in [0, 1].
double HweExact(
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
