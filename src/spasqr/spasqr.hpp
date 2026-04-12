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

#include "geno_factory/geno_data.hpp"
#include <vector>

void runSPAsqrPheno(const std::string &phenoFile,
                    const std::string &covarFile,
                    const std::string &quantPhenoCol,
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
                    const std::string &removeFile = {});
