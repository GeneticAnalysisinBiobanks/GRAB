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
    std::string binaryPheno;        // column name for binary phenotype
    std::string survPheno;          // "TIME:EVENT" survival phenotype
    std::string pcCols;             // comma-separated PC column names
    std::string bfilePrefix;
    std::string refAfFile;
    std::string spGrmGrabFile;      // --sp-grm-grab
    std::string spGrmPlink2File;    // --sp-grm-plink2 (.grm.sp file)
    std::string indAfFile;          // --ind-af-coef
    std::string pairwiseIBDFile;
    std::string outputFile;
    bool   calIndAfCoef      = false;
    bool   calPairwiseIBD    = false;
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
    int    nClusters         = 0;       // 0 = auto (from --ref-af count)
    uint64_t seed             = 0;       // 0 = use std::random_device
};

// Entry point: parse argv, print help or dispatch the selected method.
int run(int argc, char* argv[]);

// ── Shared helpers (used across cli translation units) ─────────────

Args parseArgs(int argc, char* argv[]);
void printHelp(const std::string& topic);

} // namespace cli
