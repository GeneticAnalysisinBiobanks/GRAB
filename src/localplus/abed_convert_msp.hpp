// abed_convert_msp.hpp — VCF/BCF + MSP (rfmix) to .abed converter
//
// Reads a phased VCF/BCF and a MSP local-ancestry file from rfmix,
// computes per-ancestry dosage and hapcount tracks, and writes an .abed
// binary file together with standard .bim and .fam files.
#pragma once

#include <string>

// Convert a phased VCF/BCF + MSP file to .abed/.bim/.fam.
//
//   vcfFile:    path to phased BCF/VCF/VCF.gz (used with htslib)
//   mspFile:    path to rfmix MSP file (.msp or .msp.tsv)
//   outPrefix:  output prefix — writes {outPrefix}.abed, .bim, .fam
//               (do NOT include the .abed extension in outPrefix)
//   keepFile:   PLINK2-style sample ID file — only listed subjects are kept
//               (empty = keep all)
//   removeFile: PLINK2-style sample ID file — listed subjects are excluded
//               (empty = remove none)
//
// MSP file format:
//   Line 1: #Subpopulation order/codes: 0=X\t1=Y\t...  (K inferred from this)
//   Line 2: #chm\tspos\tepos\tsgpos\tegpos\tn snps\tS0.0\tS0.1\tS1.0\tS1.1 ...
//   Data :  chrom\tspos\tepos\t...\tcall0\tcall1\t...  (2*N call columns)
//
// The number of ancestries K is auto-detected from line 1.
// Sample IDs are extracted from line 2 (even .0 suffix columns).
// Variants outside all MSP windows are silently skipped.
void convertVcfMspToAbed(
    const std::string &vcfFile,
    const std::string &mspFile,
    const std::string &outPrefix,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    int nthreads = 1
);
