// spasqr.hpp — SPAsqr: SPA-squared (multi-tau quantile regression)
//
// Pure C++17 / Eigen / Boost port of ref_code/src/mtSPAsqr.h + SPAsqr.R.
//
// Workflow (--pheno + --pheno-quantitative + --spasqr-taus):
//   1. Load phenotype/covariates, run conquer quantile regression per tau
//   2. Build smoothed residual matrix: R(i,t) = tau - Phi(-resid(i)/h)
//   3. Detect outliers per column (IQR-based, configurable)
//   4. Load sparse GRM and compute variance terms per column
//   5. Build one SPAGRMClass instance per tau
//   6. Per-marker: delegate to each SPAGRMClass, collect Z/P, compute CCT
#pragma once

#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include <memory>
#include <string>
#include <vector>

// ── GRM entry for sparse GRM data ──────────────────────────────────
struct GRMEntry {
    uint32_t row, col;
    double value;
    double factor; // 1 for diagonal, 2 for off-diagonal
};

// Load GRM entries from disk and convert to flat GRMEntry vector.
std::vector<GRMEntry> loadGrmEntries(
    const std::vector<std::string> &subjOrder,
    const std::vector<std::string> &famIIDs,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile
);

// Build SPAsqrMethod from a pre-computed residual matrix and pre-loaded GRM entries.
std::unique_ptr<MethodBase> buildSPAsqrMethod(
    Eigen::MatrixXd &ResidMat,
    const std::vector<GRMEntry> &grmEntries,
    uint32_t nUsed,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    double minMafCutoff,
    double minMacCutoff,
    std::vector<std::string> tauLabels,
    std::vector<double> *outlierRatiosOut = nullptr
);

// Multi-phenotype entry point: loads data/GRM once, parallelizes
// ntraits × ntaus conquer fits with min(nthreads, ntraits × ntaus) workers.
void runSPAsqr(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<double> &taus,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    double spasqrTol = 1e-7,
    double spasqrH = -1.0,
    double spasqrHScale = -1.0,
    const std::string &keepFile = {},
    const std::string &removeFile = {}
);

// LOCO entry point: runs per-chromosome locoEngine with precomputed
// conquer fits (O1 optimization).
void runSPAsqrLoco(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<double> &taus,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &predListFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    double spasqrTol = 1e-7,
    double spasqrH = -1.0,
    double spasqrHScale = -1.0,
    const std::string &keepFile = {},
    const std::string &removeFile = {}
);
