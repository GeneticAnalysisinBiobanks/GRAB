// sageld.hpp — SAGELD: G×E interaction analysis for longitudinal data
//              with sparse GRM relatedness correction
//
// Two input modes (the caller picks exactly one via residNames vs
// phenoNames):
//
//   Residual mode (residNames non-empty)
//     --pheno  : strict-format file with one row per IID and
//                columns #IID  R_G  R_<E1>  R_Gx<E1>  [R_<E2>  R_Gx<E2>  ...]
//     Output   : <outPrefix>.SAGELD              (single file)
//
//   Pheno mode (phenoNames non-empty)
//     --pheno     : long-format file (≥ 1 row per IID); the env column must
//                   appear under --covar-name as well so it enters X.
//     --pheno-name: list of outcome columns; one null model is fit per
//                   (outcome, env) pair via fitRandomSlopeML.
//     --sageld-x  : list of env column names; the random-effects design is
//                   Z_i = [1, E_i] within each subject.
//     Output      : <outPrefix>.<phenoName>.SAGELD   (per phenotype)
//
// In both modes the per-IID residuals (R_G, R_<E>, R_Gx<E>) feed into the
// same downstream GRM topology + SPAGRM null-model + marker-engine path.
#pragma once

#include "geno_factory/geno_data.hpp"

#include <string>
#include <vector>

void runSAGELD(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<std::string> &envNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    bool saveResid = false,
    const std::string &keepFile = {},
    const std::string &removeFile = {});
