// ibd.cpp — Pairwise IBD estimation (C++17 / Eigen)
//
// Translates the R function getPairwiseIBD() into a single-pass streaming
// algorithm over the PLINK .bed file.
//
// Algorithm outline:
//   1. Parse sparse GRM to collect off-diagonal (related) pairs.
//   2. Read .fam to get all subject IDs; build an ID→index map.
//   3. Create PlinkData for ALL subjects so that per-marker MAF is
//      computed from the full cohort (matching the R code's use of .frq).
//   4. Stream through every marker once.  For markers with MAF ≥ threshold,
//      mean-impute missing genotypes and accumulate weighted IBS0 statistics
//      for each related pair.
//   5. Derive (pa, pb, pc) from the accumulated statistics + GRM value.
//   6. Write output.

#include "spagrm/ibd.hpp"

#include "io/plink.hpp"
#include "util/logging.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>


namespace {

// ── Sparse-GRM off-diagonal entry ────────────────────────────────────────────
struct RelatedPair {
  std::string id1, id2;
  double grmValue;
};

// Parse the 3-column sparse GRM file, keeping only off-diagonal entries.
std::vector<RelatedPair> parseGRMOffDiag(const std::string& filename) {
  std::ifstream ifs(filename);
  if (!ifs) throw std::runtime_error("Cannot open sparse GRM: " + filename);

  std::vector<RelatedPair> pairs;
  std::string line;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;

    const char*       p   = line.c_str();
    const char* const end = p + line.size();
    auto skipWS  = [&]() { while (p < end && (*p == ' ' || *p == '\t')) ++p; };
    auto nextTok = [&]() -> std::string {
      skipWS();
      const char* s = p;
      while (p < end && *p != ' ' && *p != '\t') ++p;
      return std::string(s, p);
    };

    skipWS();
    if (p >= end) continue;
    std::string id1 = nextTok();
    std::string id2 = nextTok();
    skipWS();
    if (id2.empty() || p >= end) continue;

    if (id1 == id2) continue;  // skip diagonal

    char* endPtr;
    double val = std::strtod(p, &endPtr);
    if (endPtr == p) continue;

    pairs.push_back({std::move(id1), std::move(id2), val});
  }
  return pairs;
}

// ── Read .fam IIDs (column 2) ────────────────────────────────────────────────
std::vector<std::string> readFamIIDs(const std::string& famFile) {
  std::ifstream ifs(famFile);
  if (!ifs) throw std::runtime_error("Cannot open " + famFile);
  std::vector<std::string> ids;
  std::string line;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    // .fam: FID IID ...
    const char*       p   = line.c_str();
    const char* const end = p + line.size();
    auto skipWS  = [&]() { while (p < end && (*p == ' ' || *p == '\t')) ++p; };
    auto nextTok = [&]() -> std::string {
      skipWS();
      const char* s = p;
      while (p < end && *p != ' ' && *p != '\t') ++p;
      return std::string(s, p);
    };
    nextTok();                // skip FID
    std::string iid = nextTok();
    if (!iid.empty()) ids.push_back(std::move(iid));
  }
  return ids;
}

} // anonymous namespace


// ══════════════════════════════════════════════════════════════════════════════
// runPairwiseIBD
// ══════════════════════════════════════════════════════════════════════════════

void runPairwiseIBD(
    const std::string& sparseGrmFile,
    const std::string& bfilePrefix,
    const std::string& outputFile,
    double minMafIBD)
{
  // ── 1. Parse sparse GRM → related pairs ────────────────────────────
  infoMsg("Parsing sparse GRM: %s", sparseGrmFile.c_str());
  std::vector<RelatedPair> pairs = parseGRMOffDiag(sparseGrmFile);
  infoMsg("Found %zu off-diagonal (related) pairs", pairs.size());

  if (pairs.empty()) {
    // Write header-only output
    std::ofstream ofs(outputFile);
    if (!ofs) throw std::runtime_error("Cannot write " + outputFile);
    ofs << "ID1\tID2\tpa\tpb\tpc\n";
    infoMsg("No related pairs found — wrote empty output to %s", outputFile.c_str());
    return;
  }

  // ── 2. Read .fam → all subject IIDs ────────────────────────────────
  const std::string famFile = bfilePrefix + ".fam";
  const std::string bimFile = bfilePrefix + ".bim";
  const std::string bedFile = bfilePrefix + ".bed";
  std::vector<std::string> allIIDs = readFamIIDs(famFile);
  infoMsg("Read %zu subjects from %s", allIIDs.size(), famFile.c_str());

  // Build IID → index map
  std::unordered_map<std::string, uint32_t> idMap;
  idMap.reserve(allIIDs.size());
  for (uint32_t i = 0; i < static_cast<uint32_t>(allIIDs.size()); ++i)
    idMap.emplace(allIIDs[i], i);

  // ── 3. Convert pair IDs to indices, drop pairs with unknown subjects ─
  struct IndexedPair {
    uint32_t idx1, idx2;
    double   grmValue;
    // string IDs kept for output
    const std::string* sid1;
    const std::string* sid2;
  };

  std::vector<IndexedPair> idxPairs;
  idxPairs.reserve(pairs.size());
  for (auto& rp : pairs) {
    auto it1 = idMap.find(rp.id1);
    auto it2 = idMap.find(rp.id2);
    if (it1 == idMap.end() || it2 == idMap.end()) continue;
    idxPairs.push_back({it1->second, it2->second, rp.grmValue,
                        &rp.id1, &rp.id2});
  }
  if (idxPairs.empty()) {
    std::ofstream ofs(outputFile);
    if (!ofs) throw std::runtime_error("Cannot write " + outputFile);
    ofs << "ID1\tID2\tpa\tpb\tpc\n";
    infoMsg("No related pairs matched .fam — wrote empty output");
    return;
  }
  infoMsg("IBD analysis: %zu pairs, %zu subjects",
          idxPairs.size(), allIIDs.size());

  // ── 4. Create PlinkData for all subjects (full cohort MAF) ─────────
  PlinkData pdata(
      bedFile, bimFile, famFile,
      allIIDs,
      "alt-first",     // allele coding doesn't matter — diffs cancel it
      {}, {}, {}, {},   // no marker filters
      1024);            // chunk size for sequential reading

  const uint32_t nSubj    = pdata.nSubjUsed();
  const uint32_t nMarkers = pdata.nMarkers();

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

  const size_t nPairs = idxPairs.size();
  std::vector<double> sumXW(nPairs, 0.0);
  std::vector<double> sumW(nPairs, 0.0);

  // ── 6. Stream markers ──────────────────────────────────────────────
  PlinkCursor cursor(
      pdata.bedFile(),
      nMarkers,
      pdata.nSubjInFile(),
      pdata.samplePosMap(),
      pdata.isAltFirst(),
      pdata.isIdentityMap());

  Eigen::VectorXd geno(nSubj);
  std::vector<uint32_t> missingIdx;
  double altFreq, altCounts, missingRate, hweP, maf, mac;

  uint32_t nUsed = 0;  // markers actually used

  for (const auto& chunk : pdata.chunkIndices()) {
    cursor.beginSequentialBlock(chunk.front());
    for (uint64_t gi : chunk) {
      cursor.getGenotypes(gi, geno, altFreq, altCounts, missingRate,
                          hweP, maf, mac, missingIdx, /*exactHwe=*/false);

      if (maf < minMafIBD) continue;

      // Mean-impute missing genotypes
      const double meanGeno = 2.0 * altFreq;
      for (uint32_t mi : missingIdx)
        geno[mi] = meanGeno;

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
        const double d = geno[idxPairs[k].idx1] - geno[idxPairs[k].idx2];
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
    const std::string* id1;
    const std::string* id2;
    double pa, pb, pc;
  };
  std::vector<IBDResult> results(nPairs);

  for (size_t k = 0; k < nPairs; ++k) {
    double pc = (sumW[k] > 0.0) ? 0.5 * sumXW[k] / sumW[k] : 0.0;
    const double grm = idxPairs[k].grmValue;

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

    results[k] = {idxPairs[k].sid1, idxPairs[k].sid2, pa, pb, pc};
  }

  // ── 8. Write output ────────────────────────────────────────────────
  std::ofstream ofs(outputFile);
  if (!ofs) throw std::runtime_error("Cannot write " + outputFile);
  ofs << "ID1\tID2\tpa\tpb\tpc\n";
  ofs.precision(12);
  for (const auto& r : results)
    ofs << *r.id1 << '\t' << *r.id2 << '\t'
        << r.pa   << '\t' << r.pb   << '\t' << r.pc << '\n';

  infoMsg("Wrote %zu IBD records to %s", results.size(), outputFile.c_str());
}
