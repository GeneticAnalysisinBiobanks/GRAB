// hwe.hpp — Hardy-Weinberg equilibrium exact test and QC stats
//
// Single canonical implementation used by all genotype backends
// (PLINK, PGEN, VCF, BGEN) and the marker engine.
//
// Default method: exact test (SNPHWE2), matching plink2 --hardy default.
// Wigginton JE, Cutler DJ, Abecasis GR (2005). Am J Hum Genet 76:887-893.
#pragma once

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
