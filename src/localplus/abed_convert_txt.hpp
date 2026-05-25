// abed_convert_txt.hpp — Text-to-ABED conversion tool for admixed genotype data
//
// Reads gzipped text dosage/hapcount files (one pair per ancestry)
// and writes a unified .abed binary file + .bim + .fam.
#pragma once

#include <string>

// Convert extract_tracts output files to .abed binary format.
//
//   textPrefix: e.g. "data" — scans for {prefix}.anc0.dosage[.gz],
//               {prefix}.anc1.dosage[.gz], ... until the file is not found.
//   outPrefix:  output prefix for .abed/.bim/.fam
//   keepFile:   PLINK2-style sample ID file — only listed subjects are kept
//               (empty = keep all)
//   removeFile: PLINK2-style sample ID file — listed subjects are excluded
//               (empty = remove none)
//
// Text file format (from extract_tracts_fast_pgzip.py):
//   Header row:  CHROM  POS  ID  REF  ALT  SAMPLE1  SAMPLE2 ...
//   Data rows:   chr    pos  id  ref  alt  value1   value2   ...
// File naming: {prefix}.anc{k}.dosage[.gz] and {prefix}.anc{k}.hapcount[.gz]  (k=0,1,...)
void convertTextToAbed(
    const std::string &textPrefix,
    const std::string &outPrefix,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    int nthreads = 1
);
