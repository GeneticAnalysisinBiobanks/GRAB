#ifndef MARKER4ENGINE_H
#define MARKER4ENGINE_H

// Marker4Engine.h -- Pure-C++ chunk-parallel marker execution engine
// No Rcpp types in this header.

#include <RcppArmadillo.h>

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <memory>


namespace POLMM { class POLMMClass; }
namespace SPACox { class SPACoxClass; }
namespace SPAmix { class SPAmixClass; }
namespace SPAGRM { class SPAGRMClass; }
namespace SAGELD { class SAGELDClass; }
namespace WtCoxG { class WtCoxGClass; }
namespace SPAsqr { class SPAsqrClass; }
namespace LEAF { class LEAFClass; }
namespace SPAmixPlus { class SPAmixPlusClass; }

namespace Marker4 {


// ---- Data structures ----

struct RangeFilter {
  std::string chrom;
  uint32_t start;
  uint32_t end;
};

struct PlinkMarkerInfo {
  std::string chrom;
  uint32_t pos;
  std::string id;
  std::string ref;
  std::string alt;
  uint64_t genoIndex;
};

struct ReaderConfig {
  std::string genoType;
  std::string bimFile;
  std::string famFile;
  std::string bedFile;
  std::vector<std::string> sampleInModel;
  std::string alleleOrder;
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
  ReaderConfig reader;
  std::unordered_map<uint64_t, WtRow> wtMap;
  std::vector<std::unordered_map<uint64_t, LeafRow>> leafMaps;
};


// Per-worker independent copies of method objects.
struct ThreadContext {
  std::unique_ptr<POLMM::POLMMClass> polmm;
  std::unique_ptr<SPACox::SPACoxClass> spacox;
  std::unique_ptr<SPAmix::SPAmixClass> spamix;
  std::unique_ptr<SPAmixPlus::SPAmixPlusClass> spamixPlus;
  std::unique_ptr<SPAGRM::SPAGRMClass> spagrm;
  std::unique_ptr<SAGELD::SAGELDClass> sageld;
  std::unique_ptr<WtCoxG::WtCoxGClass> wtcoxg;
  std::unique_ptr<SPAsqr::SPAsqrClass> spasqr;
  std::unique_ptr<LEAF::LEAFClass> leaf;
};


// ---- PLINK file helpers ----
std::vector<PlinkMarkerInfo> readPlinkMarkerInfo(const std::string& bimFile,
                                                 const std::string& alleleOrder);

std::vector<std::string> readPlinkSampleIds(const std::string& famFile);

void validateRequestedSamples(const std::vector<std::string>& requested,
                              const std::vector<std::string>& available);

std::vector<PlinkMarkerInfo> applyPlinkMarkerFilters(
    const std::vector<PlinkMarkerInfo>& markers,
    const MarkerFilterConfig& filterConfig);

std::vector<std::vector<uint64_t>> buildPlinkChunkIndexList(
    const std::vector<PlinkMarkerInfo>& markers,
    int nMarkersEachChunk);

// ---- Output helpers ----
std::string numToStr(double x);

std::string getHeader(const std::string& method,
                      const std::string& sageldMethod,
                      const std::vector<double>& taus,
                      int nPheno,
                      int nCluster);

double cctPvalue(const std::vector<double>& pvals);


// ---- Core engine ----
ThreadContext makeThreadContext(const std::string& method);

void mainMarkerChunksCore(
  const std::string& method,
  const std::vector<std::vector<uint64_t>>& chunkMarkers,
  const std::string& outputFile,
  unsigned int nThreads,
  const std::string& imputeMethod,
  double missingCutoff,
  double minMafMarker,
  double minMacMarker,
  const ExtraParams& extraParams
);

}

#endif
