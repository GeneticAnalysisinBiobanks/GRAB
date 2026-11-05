// lanc_convert_rfmix.hpp — RFMix .msp.tsv + phased query BCF/VCF -> .lanc
//
// Multi-chromosome converter: replaces the single-file .abed converter it
// supersedes with a plane-output writer driven by the merged .lanc format
// (see dev-notes/methods/recode_rfmix/03_merged_and_fixes.md).  All
// chromosomes are emitted as ordered segments of ONE {outPrefix}.lanc.
#pragma once

#include <string>

// Convert per-chromosome RFMix .msp.tsv files + matching phased query
// BCF/VCF files to ONE merged {outPrefix}.lanc + {outPrefix}.bim plus one
// shared {outPrefix}.fam.
//
//   genoPrefix:  --bcf or --vcf prefix; per-chromosome genotype files are
//                discovered by globbing {genoPrefix}*.bcf (expectBcf) or
//                {genoPrefix}*.vcf.gz / *.vcf.zst / *.vcf (otherwise).
//   expectBcf:   true  => genoPrefix files must be BCF2 (--bcf)
//                false => genoPrefix files must not be BCF2 (--vcf)
//                Mismatches are rejected with plink2-style wording, exactly
//                as the old (now-removed) .abed converter did.
//   mspPrefix:   --rfmix-msp prefix; per-chromosome RFMix output is
//                discovered by globbing {mspPrefix}*.msp.tsv.
//   outPrefix:   output prefix; writes ONE merged {outPrefix}.lanc (one
//                chromosome segment per matched chromosome, in numeric-aware
//                sorted order) and the companion merged {outPrefix}.bim (all
//                markers in that same chromosome order) plus one shared
//                {outPrefix}.fam.
//   compressionLevel: zstd level forwarded to LancWriter (default 3).
//   keepFile / removeFile: PLINK2-style --keep/--remove sample ID files;
//                applied once (the query cohort is shared across
//                chromosomes) via io/subject_filter.hpp::buildKeptIndices.
//   nthreads:    htslib decompression threads only (hts_set_threads).
//                LancWriter itself stays single-threaded for deterministic
//                frame ordering (see CLAUDE.md reproducibility rule).
//
// Requires that every chromosome's .msp.tsv agrees on K (number of
// ancestries, 1..15) and on the query sample-ID list and order (the query
// cohort is shared across chromosomes). Requires a 1:1 match between msp
// chromosomes and genotype-file chromosomes (by the same stripped token);
// throws a clear error if either side has an unmatched chromosome.
void convertRfmixToLanc(
    const std::string &genoPrefix,
    bool expectBcf,
    const std::string &mspPrefix,
    const std::string &outPrefix,
    int compressionLevel = 3,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    int nthreads = 1
);
