// spasqr.hpp — SPAsqr: SPA-squared (multi-tau quantile regression)
//
// Pure C++17 / Eigen / Boost port of ref_code/src/mtSPAsqr.h + SPAsqr.R.
//
// Workflow:
//   1. Load multi-column residual matrix (one column per tau)
//   2. Detect outliers per column (IQR-based, configurable)
//   3. Load sparse GRM and compute variance terms per column
//   4. Build one SPAGRMClass instance per tau
//   5. Per-marker: delegate to each SPAGRMClass, collect Z/P, compute CCT
#pragma once

#include <string>

void runSPAsqr(
    const std::string& residFile,
    const std::string& spgrmSaigeFile,
    const std::string& spgrmGctaPrefix,
    const std::string& bfilePrefix,
    const std::string& outputFile,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff);
