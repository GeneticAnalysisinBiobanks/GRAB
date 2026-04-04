// admix_convert.hpp — Text-to-ABED conversion tool for admixed genotype data
//
// Reads gzipped text dosage/hapcount files (one pair per ancestry)
// and writes a unified .abed binary file + .bim + .fam.
#pragma once

#include <string>
#include <vector>

// Convert text-format admixed genotype files to .abed binary format.
//
//   textPrefix:  e.g. "simuAncestry" — expects files like
//                {textPrefix}{ancIdx}_Dosage.txt.gz, {textPrefix}{ancIdx}_HapCount.txt.gz
//   ancIdx:      ancestry indices e.g. {1, 2}
//   outPrefix:   output prefix for .abed/.bim/.fam
//
// Text file format:
//   Header row:  CHROM  POS  ID  REF  ALT  SAMPLE1  SAMPLE2 ...
//   Data rows:   chr    pos  id  ref  alt  value1   value2   ...
//   Values are integers: 0, 1, or 2 (dosage/hapcount)
//
// The first ancestry file's header defines the .bim and .fam.
// All subsequent files must have the same marker/sample order.
void convertTextToAbed(
    const std::string& textPrefix,
    const std::vector<int>& ancIdx,
    const std::string& outPrefix);

// Convert extract_tracts output files to .abed binary format.
//
//   textPrefix: e.g. "zz" — scans for {prefix}.anc0.dosage.txt[.gz],
//               {prefix}.anc1.dosage.txt[.gz], ... until the file is not found.
//   outPrefix:  output prefix for .abed/.bim/.fam (no extension appended here)
//
// Text file format (from extract_tracts_fast_pgzip.py):
//   Header row:  CHROM  POS  ID  REF  ALT  SAMPLE1  SAMPLE2 ...
//   Data rows:   chr    pos  id  ref  alt  value1   value2   ...
// Companion hapcount file: {prefix}.anc{k}.hapcount.txt[.gz]
void convertExtractTractsToAbed(
    const std::string& textPrefix,
    const std::string& outPrefix);
