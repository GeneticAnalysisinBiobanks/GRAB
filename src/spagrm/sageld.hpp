// sageld.hpp — SAGELD: G×E interaction analysis for longitudinal data
//              with sparse GRM relatedness correction
//
// Residual file column layout:
//   #IID  R_G  R_<E1>  R_Gx<E1>  [R_<E2>  R_Gx<E2>  ...]
//
// R_G is always the main-effect residual (sum of outcome residuals per IID).
// Each subsequent pair (R_<E>, R_Gx<E>) adds one G×E test.
// If no header row is present, columns are assigned by position and
// environments are named E1, E2, ...
//
// Output columns (single file, --out):
//   P_G  P_Gx<E1>  [P_Gx<E2>  ...]  Z_G  Z_Gx<E1>  [Z_Gx<E2>  ...]
#pragma once

#include "geno_factory/geno_data.hpp"

#include <string>

void runSAGELD(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &outputFile,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {}

);
