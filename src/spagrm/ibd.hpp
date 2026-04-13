// ibd.hpp — Pairwise IBD estimation from PLINK genotypes + sparse GRM
//
// Replaces the R function getPairwiseIBD() with a pure C++17 / Eigen
// implementation.  Connected-component detection uses union-find instead
// of igraph.
//
// Output: tab-separated file with columns  ID1  ID2  pa  pb  pc
#pragma once

#include "geno_factory/geno_data.hpp"

/// Compute pairwise IBD (pa, pb, pc) for every off-diagonal pair in the
/// sparse GRM and write the result to `outputFile`.
///
/// @param spgrmGrabFile   3-column TSV (ID1 ID2 VALUE), '#' lines skipped.
/// @param spgrmGctaFile  plink2 --make-grm-sparse .grm.sp file.
/// @param bfilePrefix     PLINK binary prefix (.bed/.bim/.fam).
/// @param outputFile      Tab-separated output (ID1 ID2 pa pb pc).
/// @param minMafIBD       Minimum MAF for a marker to be used (default 0.01).
void runPairwiseIBD(
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outputFile,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    double minMafIBD = 0.01,
    int nthreads = 1
);
