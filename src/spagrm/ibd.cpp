// ibd.cpp — Pairwise IBD estimation (C++17 / Eigen)
//
// Translates the R function getPairwiseIBD() into a single-pass streaming
// algorithm over the PLINK .bed file.
//
// Algorithm outline:
//   1. Load sparse GRM via SparseGRM; extract off-diagonal (related) pairs.
//   2. Create PlinkData for ALL subjects so that per-marker MAF is
//      computed from the full cohort (matching the R code's use of .frq).
//   3. Stream through every marker once.  For markers with MAF ≥ threshold,
//      mean-impute missing genotypes and accumulate weighted IBS0 statistics
//      for each related pair.
//   4. Derive (pa, pb, pc) from the accumulated statistics + GRM value.
//   5. Write output in the same pair order as the sparse GRM.

#include "spagrm/ibd.hpp"

#include "io/geno_data.hpp"
#include "io/subject_data.hpp"
#include "io/sparse_grm.hpp"
#include "util/logging.hpp"
#include "util/text_stream.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <numeric>
#include <string>
#include <vector>

#include <Eigen/Dense>





// ══════════════════════════════════════════════════════════════════════════════
// runPairwiseIBD
// ══════════════════════════════════════════════════════════════════════════════

void runPairwiseIBD(
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const GenoSpec& geno,
    const std::string& outputFile,
    double minMafIBD
) {
  const auto wallStart = std::chrono::steady_clock::now();
  const std::clock_t cpuStart = std::clock();

  // ── 1. Read sample IDs and load sparse GRM ──────────────────────
  std::vector<std::string> allIIDs = parseGenoIIDs(geno);
  const uint32_t nFam = static_cast<uint32_t>(allIIDs.size());
  infoMsg("Read %u subjects", nFam);

  SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile,
                                   allIIDs, allIIDs);
  infoMsg("Loaded sparse GRM: %zu entries", grm.nnz());

  // ── 2. Extract off-diagonal pairs (sorted by row, col — same as GRM) ─
  struct IndexedPair {
    uint32_t idx1, idx2;
    double   grmValue;
  };

  std::vector<IndexedPair> pairs;
  for (const auto& e : grm.entries()) {
    if (e.row != e.col)
      pairs.push_back({e.row, e.col, e.value});
  }
  infoMsg("Found %zu off-diagonal (related) pairs", pairs.size());

  if (pairs.empty()) {
    TextWriter writer(outputFile);
    writer.write("#ID1\tID2\tpa\tpb\tpc\n");
    infoMsg("No related pairs found — wrote empty output to %s", outputFile.c_str());
    return;
  }

  // ── 3. Build full bitmask and load genotype data ───────────────
  const size_t nMaskWords = (nFam + 63) / 64;
  std::vector<uint64_t> fullMask(nMaskWords, ~uint64_t(0));
  if (nFam % 64 != 0)
    fullMask.back() &= (uint64_t(1) << (nFam % 64)) - 1;

  // ── 4. Create genotype data for all subjects (full cohort MAF) ───
  auto genoData = makeGenoData(geno, fullMask, nFam, nFam, 1024);

  const uint32_t nSubj    = genoData->nSubjUsed();
  const uint32_t nMarkers = genoData->nMarkers();

  if (nMarkers > 200000) {
    infoMsg("WARNING: %u markers provided. "
            "IBD estimation is intended for LD-pruned markers (~100k). "
            "Consider pruning for speed and accuracy.", nMarkers);
  }

  infoMsg("IBD analysis: %u subjects, %u markers (MAF threshold %.4f)",
          nSubj, nMarkers, minMafIBD);

  // ── 5. Initialise per-pair accumulators ────────────────────────────
  //
  // For each marker m with MAF ≥ threshold:
  //   p_m         = altFreq       (allele frequency)
  //   pro_var_m   = 2 * (p*(1-p))^2
  //   w_m         = sqrt( pro_var / (1 - pro_var) )
  //
  // For each pair (i,j):
  //   d_m  = G_i[m] - G_j[m]          (genotype difference)
  //   x_m  = (|d-1| + |d+1| - 2) / pro_var   (IBS0 kernel)
  //
  //   sumXW += x_m * w_m
  //   sumW  += w_m
  //
  // Then  pc_raw = 0.5 * sumXW / sumW.

  const size_t nPairs = pairs.size();
  std::vector<double> sumXW(nPairs, 0.0);
  std::vector<double> sumW(nPairs, 0.0);

  // ── 6. Stream markers ──────────────────────────────────────────────
  auto cursor = genoData->makeCursor();

  Eigen::VectorXd genoVec(nSubj);
  std::vector<uint32_t> missingIdx;
  double altFreq, altCounts, missingRate, hweP, maf, mac;

  uint32_t nUsed = 0;  // markers actually used

  for (const auto& chunk : genoData->chunkIndices()) {
    cursor->beginSequentialBlock(chunk.front());
    for (uint64_t gi : chunk) {
      cursor->getGenotypes(gi, genoVec, altFreq, altCounts, missingRate,
                          hweP, maf, mac, missingIdx, /*exactHwe=*/false);

      if (maf < minMafIBD) continue;

      // Mean-impute missing genotypes
      const double meanGeno = 2.0 * altFreq;
      for (uint32_t mi : missingIdx)
        genoVec[mi] = meanGeno;

      // Pre-compute weights for this marker
      const double pq  = altFreq * (1.0 - altFreq);   // p*(1-p)
      const double pv  = 2.0 * pq * pq;               // pro_var
      const double opv = 1.0 - pv;                     // 1 - pro_var

      // Guard against degenerate markers (MAF ~0 or ~0.5)
      if (pv < 1e-30 || opv <= 0.0) continue;

      const double inv_pv = 1.0 / pv;
      const double w      = std::sqrt(pv / opv);

      // Accumulate for each pair
      for (size_t k = 0; k < nPairs; ++k) {
        const double d = genoVec[pairs[k].idx1] - genoVec[pairs[k].idx2];
        const double x = (std::abs(d - 1.0) + std::abs(d + 1.0) - 2.0) * inv_pv;
        sumXW[k] += x * w;
        sumW[k]  += w;
      }
      ++nUsed;
    }
  }

  infoMsg("IBD: used %u / %u markers (MAF >= %.4f)", nUsed, nMarkers, minMafIBD);

  // ── 7. Derive (pa, pb, pc) for each pair ───────────────────────────
  //
  // pc_raw = 0.5 * sumXW / sumW
  //
  // Constraints (from R code):
  //   upper_bound = (1 - grmValue)^2 - 1e-10
  //   lower_bound = 1 - 2 * grmValue
  //   pc = clamp(pc_raw, lower_bound, upper_bound)
  //   pb = 2 - 2*pc - 2*grmValue
  //   pa = 2*grmValue + pc - 1
  //   if (pb < 0) { pa += 0.5*pb; pb = 0; pc = 0; }

  struct IBDResult {
    uint32_t idx1, idx2;
    double pa, pb, pc;
  };
  std::vector<IBDResult> results(nPairs);

  for (size_t k = 0; k < nPairs; ++k) {
    double pc = (sumW[k] > 0.0) ? 0.5 * sumXW[k] / sumW[k] : 0.0;
    const double grm = pairs[k].grmValue;

    const double upper = (1.0 - grm) * (1.0 - grm) - 1e-10;
    const double lower = 1.0 - 2.0 * grm;
    pc = std::clamp(pc, lower, upper);

    double pb = 2.0 - 2.0 * pc - 2.0 * grm;
    double pa = 2.0 * grm + pc - 1.0;

    if (pb < 0.0) {
      pa += 0.5 * pb;
      pb = 0.0;
      pc = 0.0;
    }

    results[k] = {pairs[k].idx1, pairs[k].idx2, pa, pb, pc};
  }

  // ── 8. Write output (same pair order as the sparse GRM) ────────────
  TextWriter writer(outputFile);
  writer.write("#ID1\tID2\tpa\tpb\tpc\n");
  char buf[128];
  for (const auto& r : results) {
    std::string line;
    line.reserve(128);
    line += allIIDs[r.idx1]; line += '\t';
    line += allIIDs[r.idx2]; line += '\t';
    int n = std::snprintf(buf, sizeof(buf), "%.17g\t%.17g\t%.17g\n", r.pa, r.pb, r.pc);
    line.append(buf, n);
    writer.write(line);
  }

  infoMsg("Wrote %zu IBD records to %s", results.size(), outputFile.c_str());

  const double wallSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - wallStart).count();
  const double cpuSec = static_cast<double>(std::clock() - cpuStart) / CLOCKS_PER_SEC;
  infoMsg("Wall time: %.1f seconds, CPU time: %.1f seconds", wallSec, cpuSec);
}
