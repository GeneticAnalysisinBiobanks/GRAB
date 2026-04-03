// spasqr.hpp — SPAsqr: SPA-squared (multi-tau quantile regression)
//
// Pure C++17 / Eigen / Boost port of ref_code/src/mtSPAsqr.h + SPAsqr.R.
//
// Legacy workflow (--null-resid):
//   1. Load multi-column residual matrix (one column per tau)
//   2. Detect outliers per column (IQR-based, configurable)
//   3. Load sparse GRM and compute variance terms per column
//   4. Build one SPAGRMClass instance per tau
//   5. Per-marker: delegate to each SPAGRMClass, collect Z/P, compute CCT
//
// Pheno workflow (--pheno + --pheno-quantitative + --spasqr-taus):
//   1. Load phenotype/covariates, run conquer quantile regression per tau
//   2. Build smoothed residual matrix: R(i,t) = tau - Phi(-resid(i)/h)
//   3-5. Same as above
#pragma once

#include "io/geno_data.hpp"
#include <vector>

void runSPAsqr(
    const std::string& residFile,
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const GenoSpec& geno,
    const std::string& outputFile,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff);

void runSPAsqrPheno(
    const std::string& phenoFile,
    const std::string& covarFile,
    const std::string& quantPhenoCol,
    const std::vector<std::string>& covarNames,
    const std::vector<double>& taus,
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const GenoSpec& geno,
    const std::string& outputFile,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff);
