// ibd.hpp — Pairwise IBD estimation from PLINK genotypes + sparse GRM
//
// Replaces the R function getPairwiseIBD() with a pure C++17 / Eigen
// implementation.  Connected-component detection uses union-find instead
// of igraph.
//
// Output: tab-separated file with columns  ID1  ID2  pa  pb  pc
#pragma once

#include <string>

/// Compute pairwise IBD (pa, pb, pc) for every off-diagonal pair in the
/// sparse GRM and write the result to `outputFile`.
///
/// @param sparseGrmFile  3-column TSV (ID1 ID2 VALUE), '#' lines skipped.
/// @param bfilePrefix    PLINK binary prefix (.bed/.bim/.fam).
/// @param outputFile     Tab-separated output (ID1 ID2 pa pb pc).
/// @param minMafIBD      Minimum MAF for a marker to be used (default 0.01).
void runPairwiseIBD(
    const std::string& sparseGrmFile,
    const std::string& bfilePrefix,
    const std::string& outputFile,
    double minMafIBD = 0.01);
