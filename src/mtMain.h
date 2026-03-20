// Marker4.h -- Pure-C++ chunk-parallel marker execution engine interface
//
// Functions in this file:
//
//   [Engine]
//     mainMarkerChunksCore — thread-pool dispatcher: reads genotypes, runs statistics,
//                            and writes results for all marker chunks

#ifndef MARKER4_H
#define MARKER4_H

#include <RcppArmadillo.h>
#include <string>
#include <vector>
#include <cstdint>
#include <sstream>
#include <cmath>
#include <limits>
#include <cstdio>

// Forward declarations — only the method class pointers used by the engine
namespace POLMM    { class POLMMClass; }
namespace SPACox   { class SPACoxClass; }
namespace SPAmix   { class SPAmixClass; }
namespace SPAGRM   { class SPAGRMClass; }
namespace SAGELD   { class SAGELDClass; }
namespace WtCoxG   { class WtCoxGClass; }
namespace SPAsqr   { class SPAsqrClass; }
namespace LEAF     { class LEAFClass; }
namespace SPAmixPlus { class SPAmixPlusClass; }
namespace PLINK4   { class PlinkData; }

namespace Marker4 {

// ---- Engine ----

// Thread-pool dispatcher: for each chunk, reads genotypes from PLINK,
// computes per-marker statistics via the active method object, and writes
// results to outputFile (plain text or .gz).
void mainMarkerChunksCore(
  const std::string&                        method,
  const std::vector<std::vector<uint64_t>>& chunkMarkers,
  const std::string&                        outputFile,
  unsigned int                              nthreads,
  const std::string&                        impute_method,
  double                                    missing_cutoff,
  double                                    min_maf_marker,
  double                                    min_mac_marker,
  const PLINK4::PlinkData&                  plinkData
);

} // namespace Marker4

#endif // MARKER4_H
