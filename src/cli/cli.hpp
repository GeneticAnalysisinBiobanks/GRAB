// cli.hpp — Public interface for the GRAB command-line layer
#pragma once

#include <cstdint>
#include <string>

namespace cli {

struct Args {
    std::string method;
    std::string helpTopic;          // set when --help <topic> is used
    std::string residFile;
    std::string phenoFile;
    std::string covarFile;
    std::string covarName;          // comma-separated covariate column names
    std::string covarColNums;       // --covar-col-nums: 1-based column numbers/ranges
    std::string notCovar;           // --not-covar: exclude covariate(s) by name
    std::string binaryPheno;        // column name for binary phenotype
    std::string survPheno;          // "TIME:EVENT" survival phenotype
    std::string pcCols = "PC1,PC2,PC3,PC4"; // comma-separated PC column names (default: PC1,PC2,PC3,PC4)
    std::string quantPheno;         // column name for quantitative phenotype (SPAsqr)
    std::string ordinalPheno;       // column name for ordinal phenotype (POLMM)
    std::string spasqrTaus = "0.1,0.3,0.5,0.7,0.9";  // default tau levels (SPAsqr)
    std::string bfilePrefix;
    std::string pfilePrefix;        // --pfile (pgen/pvar/psam)
    std::string vcfFile;            // --vcf (vcf/bcf)
    std::string bgenFile;           // --bgen
    std::string refAfFile;
    std::string spGrmGrabFile;      // --sp-grm-grab
    std::string spGrmPlink2File;    // --sp-grm-plink2 (.grm.sp file)
    std::string indAfFile;          // --ind-af-coef
    std::string pairwiseIBDFile;
    std::string outputFile;
    std::string outPrefix;            // --out-prefix (GWAS multi-resid / --make-abed)
    std::string compression;          // --compression (gz|zst, default: empty = plain text)
    std::string extractFile;        // --extract (SNP include list)
    std::string excludeFile;        // --exclude (SNP exclude list)
    std::string admixBfilePrefix;   // --admix-bfile
    std::string admixPhiFile;       // --admix-phi
    std::string mspFile;            // --rfmix-msp
    std::string admixTextPrefix;    // --admix-text-prefix
    std::string keepFile;           // --keep (subject include list)
    std::string removeFile;         // --remove (subject exclude list)
    bool   calIndAfCoef      = false;
    bool   calPairwiseIBD    = false;
    bool   calAdmixPhi       = false;
    bool   makeAbed          = false;
    double minMafIBD         = 0.01;
    double refPrevalence     = -1.0;
    double cutoff            = 0.05;
    double spaCutoff         = 2.0;
    double pvalCovAdjCut     = 5e-5;
    double missingCutoff     = 0.1;
    double minMafCutoff      = 1e-4;
    double minMacCutoff      = 10.0;
    double outlierRatio      = 1.5;
    double outlierAbsBound   = 0.55;
    int    nthread           = 1;
    int    nSnpPerChunk      = 8192;
    int    compressionLevel  = 0;       // --compression-level (0 = library default)
    int    nClusters         = 0;       // 0 = auto (from --ref-af count)
    uint64_t seed             = 0;       // 0 = use std::random_device
};

// Entry point: parse argv, print help or dispatch the selected method.
int run(int argc, char* argv[]);

// ── Shared helpers (used across cli translation units) ─────────────

Args parseArgs(int argc, char* argv[]);
void printHelp(const std::string& topic);

} // namespace cli
