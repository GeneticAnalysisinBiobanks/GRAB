// gt_prob.hpp — SPAGRM null model and full method runner
//
// Refactors SPAGRM.NullModel() + GRAB.mtMarker() from R/Rcpp to pure C++17.
//
// The null model computes:
//   - Outlier detection (IQR-based with adaptive ratio)
//   - Connected-component decomposition of GRM graph (union-find)
//   - Greedy family splitting (MaxNuminFam)
//   - Two-subject & three-subject outlier family structures
//   - Chow-Liu tree genotype probability arrays
//
// These feed into SPAGRMClass for per-marker SPA p-values.
#pragma once

#include "io/geno_data.hpp"

/// Run the full SPAGRM workflow: null model + marker-level SPA.
///
/// @param residFile        2-column residual file (SubjID  Resid), '#' header ok
/// @param spgrmGrabFile    3-column TSV (ID1 ID2 VALUE), '#' lines skipped
/// @param spgrmGctaFile   plink2 --make-grm-sparse .grm.sp file
/// @param pairwiseIBDFile  Tab-separated (ID1 ID2 pa pb pc) from runPairwiseIBD
/// @param bfilePrefix      PLINK binary prefix (.bed/.bim/.fam)
/// @param outputFile       Marker-level output file
/// @param spaCutoff        SPA z-score threshold (default 2.0)
/// @param nthreads         Number of worker threads
/// @param nSnpPerChunk     Markers per chunk for parallel engine
/// @param missingCutoff    Per-marker missing rate cutoff
/// @param minMafCutoff     Min minor allele frequency
/// @param minMacCutoff     Min minor allele count
void runSPAGRM(
    const std::string& residFile,
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const std::string& pairwiseIBDFile,
    const GenoSpec& geno,
    const std::string& outputFile,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff);
