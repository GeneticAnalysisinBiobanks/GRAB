// Marker4Engine.h -- Pure-C++ chunk-parallel marker execution engine interface
//
// Functions in this file:
//
//   [Data structures]
//     PlinkMarkerInfo    — one marker's metadata from a PLINK .bim file
//     ReaderConfig       — PLINK file paths and sample selection for the reader
//     MarkerFilterConfig — include/exclude lists and genomic ranges for marker filtering
//     WtRow              — per-marker external weight data for WtCoxG
//     LeafRow            — per-marker external weight data for LEAF
//     ExtraParams        — bundle of ReaderConfig + WtRow/LeafRow lookup maps
//
//   [Inline helpers]
//     numToStr           — format a double as a short string ("NA", "Inf", or %.6g)
//     cctPvalue          — Cauchy combination test aggregate p-value
//     appendCell         — append a tab-separated cell to an ostringstream
//     appendMeta         — append the 8 standard marker meta-columns
//
//   [Reader]
//     readPlinkMarkerInfo       — parse a PLINK .bim file into a PlinkMarkerInfo vector
//     readPlinkSampleIds        — parse a PLINK .fam file for sample IDs
//     validateRequestedSamples  — verify every requested sample exists in the .fam file
//     applyPlinkMarkerFilters   — filter markers by ID/range include-exclude lists
//     buildPlinkChunkIndexList  — partition markers into per-chromosome chunks
//
//   [Engine]
//     mainMarkerChunksCore — thread-pool dispatcher: reads genotypes, runs statistics,
//                            and writes results for all marker chunks

#ifndef MARKER4ENGINE_H
#define MARKER4ENGINE_H

#include <RcppArmadillo.h>
#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstdio>

namespace POLMM    { class POLMMClass; }
namespace SPACox   { class SPACoxClass; }
namespace SPAmix   { class SPAmixClass; }
namespace SPAGRM   { class SPAGRMClass; }
namespace SAGELD   { class SAGELDClass; }
namespace WtCoxG   { class WtCoxGClass; }
namespace SPAsqr   { class SPAsqrClass; }
namespace LEAF     { class LEAFClass; }
namespace SPAmixPlus { class SPAmixPlusClass; }

namespace Marker4 {

// ---- Data structures ----

struct PlinkMarkerInfo {
  std::string chrom;
  uint32_t    pos;
  std::string id;
  std::string ref;
  std::string alt;
  uint64_t    genoIndex;
};

struct ReaderConfig {
  std::string              genoType;
  std::string              bimFile;
  std::string              famFile;
  std::string              bedFile;
  std::vector<std::string> sampleInModel;
  std::string              alleleOrder;
};

struct MarkerFilterConfig {
  std::string idsToIncludeFile;
  std::string rangesToIncludeFile;
  std::string idsToExcludeFile;
  std::string rangesToExcludeFile;
};

struct WtRow {
  double AF_ref;
  double AN_ref;
  double TPR;
  double sigma2;
  double pvalue_bat;
  double w_ext;
  double var_ratio_w0;
  double var_ratio_int;
  double var_ratio_ext;
};

struct LeafRow {
  double AF_ref;
  double AN_ref;
  double TPR;
  double sigma2;
  double pvalue_bat;
  double w_ext;
  double var_ratio_w0;
  double var_ratio_int;
  double var_ratio_ext;
};

struct ExtraParams {
  ReaderConfig                                        reader;
  std::unordered_map<uint64_t, WtRow>                wtMap;
  std::vector<std::unordered_map<uint64_t, LeafRow>> leafMaps;
};

// ---- Inline helpers ----

// Format a double as a short string: "NA", "Inf", "-Inf", or %.6g.
inline std::string numToStr(double x) {
  if (std::isnan(x)) return "NA";
  if (std::isinf(x)) return (x > 0 ? "Inf" : "-Inf");
  char buf[32];
  std::snprintf(buf, sizeof(buf), "%.6g", x);
  return buf;
}

// Cauchy combination test: aggregate a vector of p-values into one p-value.
inline double cctPvalue(const std::vector<double>& pvals) {
  std::vector<double> p;
  p.reserve(pvals.size());
  for (double x : pvals)
    if (!std::isnan(x)) p.push_back(x);
  if (p.empty()) return std::numeric_limits<double>::quiet_NaN();
  double tStat = 0.0;
  for (double x : p) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) x = 0.999;
    tStat += std::tan((0.5 - x) * M_PI);
  }
  tStat /= static_cast<double>(p.size());
  if (tStat > 1e15) return (1.0 / tStat) / M_PI;
  return 0.5 - std::atan(tStat) / M_PI;
}

// Append one tab-separated value to an output row buffer.
inline void appendCell(std::ostringstream& oss, const std::string& v, bool first = false) {
  if (!first) oss << '\t';
  oss << v;
}
inline void appendCell(std::ostringstream& oss, double v, bool first = false) {
  appendCell(oss, numToStr(v), first);
}

// Append the 8 standard marker meta-columns (Marker, Chr, Pos, Ref, Alt,
// AltCount, AltFreq, MissRate) as the first cells of an output row.
inline void appendMeta(std::ostringstream& oss,
                       const std::string& marker,
                       const std::string& chr,
                       uint32_t pos,
                       const std::string& ref,
                       const std::string& alt,
                       double altCounts,
                       double altFreq,
                       double missingRate) {
  appendCell(oss, marker, true);
  appendCell(oss, chr);
  appendCell(oss, std::to_string(pos));
  appendCell(oss, ref);
  appendCell(oss, alt);
  appendCell(oss, altCounts);
  appendCell(oss, altFreq);
  appendCell(oss, missingRate);
}

// ---- Reader ----

// Parse a PLINK .bim file and return one PlinkMarkerInfo per marker line.
std::vector<PlinkMarkerInfo> readPlinkMarkerInfo(const std::string& bimFile,
                                                 const std::string& alleleOrder);

// Parse a PLINK .fam file and return the individual ID (column 2) of each sample.
std::vector<std::string> readPlinkSampleIds(const std::string& famFile);

// Throw if any sample in 'requested' is absent from 'available'.
void validateRequestedSamples(const std::vector<std::string>& requested,
                               const std::vector<std::string>& available);

// Apply ID/range include-exclude filters and return the surviving markers.
std::vector<PlinkMarkerInfo> applyPlinkMarkerFilters(
    const std::vector<PlinkMarkerInfo>& markers,
    const MarkerFilterConfig& filterConfig);

// Partition markers into chromosome-contiguous chunks of at most nMarkersEachChunk.
// Returns one inner vector of genoIndex values per chunk.
std::vector<std::vector<uint64_t>> buildPlinkChunkIndexList(
    const std::vector<PlinkMarkerInfo>& markers,
    int nMarkersEachChunk);

// ---- Engine ----

// Thread-pool dispatcher: for each chunk, reads genotypes from PLINK,
// computes per-marker statistics via the active method object, and writes
// results to outputFile (plain text or .gz).
void mainMarkerChunksCore(
  const std::string&                        method,
  const std::vector<std::vector<uint64_t>>& chunkMarkers,
  const std::string&                        outputFile,
  unsigned int                              nThreads,
  const std::string&                        imputeMethod,
  double                                    missingCutoff,
  double                                    minMafMarker,
  double                                    minMacMarker,
  const ExtraParams&                        extraParams
);

} // namespace Marker4

#endif // MARKER4ENGINE_H
