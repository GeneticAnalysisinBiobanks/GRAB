// Marker4Engine.cpp -- Thread pool, PLINK reader, per-chunk worker, sequential writer,
//                      and per-method Rcpp entry points (via Marker4Bridge.h)
//
// Functions by role:
//
//   [File-local structs]
//     RangeFilter        — a genomic interval used for marker include/exclude filtering
//     ThreadContext      — per-worker copies of the active statistical method object
//
//   [IO helpers]          (static, file-local)
//     splitWhitespace          — split a line into whitespace-delimited tokens
//     readSingleColumnFile     — read a single-column text file (e.g. marker ID list)
//     readRangeFile            — parse a three-column chrom/start/end range file
//     markerInRanges           — test whether a marker falls inside any RangeFilter
//
//   [Reader]              (Marker4 namespace, external linkage)
//     readPlinkMarkerInfo       — parse .bim into a PlinkMarkerInfo vector
//     readPlinkSampleIds        — parse .fam into a sample ID vector
//     validateRequestedSamples  — verify all model samples exist in .fam
//     applyPlinkMarkerFilters   — apply ID/range include-exclude filters
//     buildPlinkChunkIndexList  — partition markers into per-chromosome chunks
//
//   [Output header]       (static)
//     getHeader           — produce the TSV column header for a given method
//
//   [Worker context]      (static)
//     makeThreadContext   — clone the global method object into a per-thread ThreadContext
//
//   [Engine]              (Marker4 namespace, external linkage)
//     mainMarkerChunksCore — coordinate the PLINK reader, worker pool, and writer thread
//
//   [Bridge functions]    (via #include "Marker4Bridge.h", compiled in this TU)
//     runEngineCore             — build filter/reader configs and invoke mainMarkerChunksCore
//     splitVec / splitUvec / splitMat / matRowsToVecs — reconstruct arma containers from R
//     runMarkerInCPP_POLMM      — create POLMMClass and dispatch to runEngineCore
//     runMarkerInCPP_SPACox     — create SPACoxClass and dispatch to runEngineCore
//     runMarkerInCPP_SPAmix     — create SPAmixClass and dispatch to runEngineCore
//     runMarkerInCPP_SPAmixPlus — create SPAmixPlusClass and dispatch to runEngineCore
//     runMarkerInCPP_SPAGRM     — create SPAGRMClass and dispatch to runEngineCore
//     runMarkerInCPP_SAGELD     — create SAGELDClass and dispatch to runEngineCore
//     runMarkerInCPP_WtCoxG     — create WtCoxGClass, populate wtMap, dispatch
//     runMarkerInCPP_SPAsqr     — create SPAsqrClass and dispatch to runEngineCore
//     runMarkerInCPP_LEAF       — create LEAFClass, populate leafMaps, dispatch

#include "Marker4Engine.h"
#include "PLINK4.h"
#include "POLMM.h"
#include "SPACox.h"
#include "SPAmix.h"
#include "SPAGRM.h"
#include "SAGELD.h"
#include "WtCoxG.h"
#include "SPAsqr.h"
#include "LEAF.h"
#include "SPAmixPlus.h"

#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <unordered_set>
#include <cmath>
#include <exception>
#include <cstdio>

#include <boost/math/distributions/chi_squared.hpp>
#include <zlib.h>


// ---- Global method pointers (one per active statistical method) ----
POLMM::POLMMClass*           ptr_gPOLMMobj      = nullptr;
SPACox::SPACoxClass*         ptr_gSPACoxobj     = nullptr;
SPAmix::SPAmixClass*         ptr_gSPAmixobj     = nullptr;
SPAGRM::SPAGRMClass*         ptr_gSPAGRMobj     = nullptr;
SAGELD::SAGELDClass*         ptr_gSAGELDobj     = nullptr;
WtCoxG::WtCoxGClass*         ptr_gWtCoxGobj     = nullptr;
SPAsqr::SPAsqrClass*         ptr_gSPAsqrobj     = nullptr;
LEAF::LEAFClass*             ptr_gLEAFobj       = nullptr;
SPAmixPlus::SPAmixPlusClass* ptr_gSPAmixPlusobj = nullptr;

namespace Marker4 {

// ---- File-local structs ----

// A genomic interval used for marker include/exclude range filtering.
struct RangeFilter {
  std::string chrom;
  uint32_t    start;
  uint32_t    end;
};

// Per-worker independent copies of the active method object.
struct ThreadContext {
  std::unique_ptr<POLMM::POLMMClass>           polmm;
  std::unique_ptr<SPACox::SPACoxClass>         spacox;
  std::unique_ptr<SPAmix::SPAmixClass>         spamix;
  std::unique_ptr<SPAmixPlus::SPAmixPlusClass> spamixPlus;
  std::unique_ptr<SPAGRM::SPAGRMClass>         spagrm;
  std::unique_ptr<SAGELD::SAGELDClass>         sageld;
  std::unique_ptr<WtCoxG::WtCoxGClass>         wtcoxg;
  std::unique_ptr<SPAsqr::SPAsqrClass>         spasqr;
  std::unique_ptr<LEAF::LEAFClass>             leaf;
};

// ---- IO helpers ----

static std::vector<std::string> splitWhitespace(const std::string& line) {
  std::istringstream iss(line);
  std::vector<std::string> tokens;
  std::string token;
  while (iss >> token) tokens.push_back(token);
  return tokens;
}

static std::vector<std::string> readSingleColumnFile(const std::string& path) {
  std::ifstream in(path.c_str());
  if (!in.is_open())
    throw std::runtime_error("Cannot open marker ID filter file: " + path);
  std::vector<std::string> values;
  std::string line;
  while (std::getline(in, line)) {
    auto tokens = splitWhitespace(line);
    if (!tokens.empty()) values.push_back(tokens[0]);
  }
  return values;
}

static std::vector<RangeFilter> readRangeFile(const std::string& path) {
  std::ifstream in(path.c_str());
  if (!in.is_open())
    throw std::runtime_error("Cannot open range filter file: " + path);
  std::vector<RangeFilter> ranges;
  std::string line;
  while (std::getline(in, line)) {
    auto tokens = splitWhitespace(line);
    if (tokens.empty()) continue;
    if (tokens.size() != 3)
      throw std::runtime_error("Range filter file should include exactly three whitespace-separated columns: " + path);
    ranges.push_back(RangeFilter{
      tokens[0],
      static_cast<uint32_t>(std::stoul(tokens[1])),
      static_cast<uint32_t>(std::stoul(tokens[2]))
    });
  }
  return ranges;
}

static bool markerInRanges(const PlinkMarkerInfo& marker,
                           const std::vector<RangeFilter>& ranges) {
  for (const auto& r : ranges) {
    if (marker.chrom == r.chrom && marker.pos >= r.start && marker.pos <= r.end)
      return true;
  }
  return false;
}

// ---- Reader implementations ----

std::vector<PlinkMarkerInfo> readPlinkMarkerInfo(const std::string& bimFile,
                                                 const std::string& alleleOrder) {
  std::ifstream in(bimFile.c_str());
  if (!in.is_open())
    throw std::runtime_error("Cannot open .bim file: " + bimFile);

  std::vector<PlinkMarkerInfo> markers;
  std::string line;
  uint64_t markerIndex = 0;
  while (std::getline(in, line)) {
    auto tokens = splitWhitespace(line);
    if (tokens.empty()) continue;
    if (tokens.size() != 6)
      throw std::runtime_error("PLINK .bim file should include 6 whitespace-separated columns: " + bimFile);

    std::string ref = tokens[5], alt = tokens[4];
    if (alleleOrder == "ref-first") { ref = tokens[4]; alt = tokens[5]; }
    std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
    std::transform(alt.begin(), alt.end(), alt.begin(), ::toupper);

    markers.push_back(PlinkMarkerInfo{
      tokens[0],
      static_cast<uint32_t>(std::stoul(tokens[3])),
      tokens[1], ref, alt, markerIndex
    });
    ++markerIndex;
  }
  return markers;
}

std::vector<std::string> readPlinkSampleIds(const std::string& famFile) {
  std::ifstream in(famFile.c_str());
  if (!in.is_open())
    throw std::runtime_error("Cannot open .fam file: " + famFile);
  std::vector<std::string> ids;
  std::string line;
  while (std::getline(in, line)) {
    auto tokens = splitWhitespace(line);
    if (tokens.empty()) continue;
    if (tokens.size() != 6)
      throw std::runtime_error("PLINK .fam file should include 6 whitespace-separated columns: " + famFile);
    ids.push_back(tokens[1]);
  }
  return ids;
}

void validateRequestedSamples(const std::vector<std::string>& requested,
                              const std::vector<std::string>& available) {
  std::unordered_set<std::string> avail(available.begin(), available.end());
  for (const auto& s : requested) {
    if (avail.find(s) == avail.end())
      throw std::runtime_error("At least one subject requested is not in PLINK file.");
  }
}

std::vector<PlinkMarkerInfo> applyPlinkMarkerFilters(
    const std::vector<PlinkMarkerInfo>& markers,
    const MarkerFilterConfig& filterConfig) {

  std::unordered_set<std::string> includeIds, excludeIds;
  std::vector<RangeFilter> includeRanges, excludeRanges;

  if (!filterConfig.idsToIncludeFile.empty()) {
    auto ids = readSingleColumnFile(filterConfig.idsToIncludeFile);
    includeIds.insert(ids.begin(), ids.end());
  }
  if (!filterConfig.rangesToIncludeFile.empty())
    includeRanges = readRangeFile(filterConfig.rangesToIncludeFile);
  if (!filterConfig.idsToExcludeFile.empty()) {
    auto ids = readSingleColumnFile(filterConfig.idsToExcludeFile);
    excludeIds.insert(ids.begin(), ids.end());
  }
  if (!filterConfig.rangesToExcludeFile.empty())
    excludeRanges = readRangeFile(filterConfig.rangesToExcludeFile);

  const bool anyInclude = !includeIds.empty() || !includeRanges.empty();
  const bool anyExclude = !excludeIds.empty() || !excludeRanges.empty();

  std::vector<PlinkMarkerInfo> filtered;
  filtered.reserve(markers.size());
  for (const auto& m : markers) {
    bool inc = !anyInclude;
    if (anyInclude)
      inc = (includeIds.count(m.id) > 0) || markerInRanges(m, includeRanges);
    if (!inc) continue;

    bool exc = false;
    if (anyExclude)
      exc = (excludeIds.count(m.id) > 0) || markerInRanges(m, excludeRanges);
    if (!exc) filtered.push_back(m);
  }
  return filtered;
}

std::vector<std::vector<uint64_t>> buildPlinkChunkIndexList(
    const std::vector<PlinkMarkerInfo>& markers,
    int nMarkersEachChunk) {

  std::vector<std::vector<uint64_t>> chunks;
  size_t start = 0;
  while (start < markers.size()) {
    const std::string chrom = markers[start].chrom;
    size_t chromEnd = start;
    while (chromEnd < markers.size() && markers[chromEnd].chrom == chrom) ++chromEnd;

    for (size_t cs = start; cs < chromEnd; cs += static_cast<size_t>(nMarkersEachChunk)) {
      size_t ce = std::min(cs + static_cast<size_t>(nMarkersEachChunk), chromEnd);
      std::vector<uint64_t> chunk;
      chunk.reserve(ce - cs);
      for (size_t i = cs; i < ce; ++i) chunk.push_back(markers[i].genoIndex);
      chunks.push_back(std::move(chunk));
    }
    start = chromEnd;
  }
  return chunks;
}


// ---- Output header ----

static const char* META_HEADER = "Marker\tChr\tPos\tRef\tAlt\tAltCount\tAltFreq\tMissRate";

// Build the TSV column header line for the given method.
static std::string getHeader(const std::string& method,
                      const std::string& sageldMethod,
                      const std::vector<double>& taus,
                      int nPheno,
                      int nCluster) {
  std::ostringstream oss;
  oss << META_HEADER;

  if (method == "POLMM") {
    oss << "\tPvalue\tbeta\tseBeta\tzScore";
  } else if (method == "SPACox") {
    oss << "\tPvalue\tzScore";
  } else if (method == "SPAmix" || method == "SPAmixPlus") {
    oss << "\tPheno\tPvalue\tzScore";
  } else if (method == "SPAGRM") {
    oss << "\tzScore\tPvalue\thwepval";
  } else if (method == "SAGELD") {
    if (sageldMethod == "GALLOP")
      oss << "\tMethod\tBeta_G\tBeta_GxE\tSE_G\tSE_GxE\tPvalue_G\tPvalue_GxE\thwepval";
    else
      oss << "\tMethod\tzScore_G\tzScore_GxE\tPvalue_G\tPvalue_GxE\thwepval";
  } else if (method == "WtCoxG") {
    oss << "\tWtCoxG.ext\tWtCoxG.noext\tzscore.ext\tzscore.noext"
        << "\tAF_ref\tAN_ref\tpvalue_bat";
  } else if (method == "SPAsqr") {
    oss << "\thwepval";
    for (double tau : taus) oss << "\tZ_tau" << numToStr(tau);
    for (double tau : taus) oss << "\tP_tau" << numToStr(tau);
    oss << "\tP_CCT";
  } else if (method == "LEAF") {
    oss << "\tmeta.p_ext\tmeta.p_noext";
    for (int i = 1; i <= nCluster; ++i)
      oss << "\tcl" << i << ".p_ext\tcl" << i << ".p_noext\tcl" << i << ".p_batch";
  } else {
    throw std::runtime_error("Unsupported method in Marker4: " + method);
  }
  return oss.str();
}


// ---- Worker context ----

// Clone the global method object into a per-thread ThreadContext.
static ThreadContext makeThreadContext(const std::string& method) {
  ThreadContext ctx;
  if (method == "POLMM") {
    if (!ptr_gPOLMMobj) throw std::runtime_error("POLMM object is not initialized.");
    ctx.polmm.reset(new POLMM::POLMMClass(*ptr_gPOLMMobj));
  } else if (method == "SPACox") {
    if (!ptr_gSPACoxobj) throw std::runtime_error("SPACox object is not initialized.");
    ctx.spacox.reset(new SPACox::SPACoxClass(*ptr_gSPACoxobj));
  } else if (method == "SPAmix") {
    if (!ptr_gSPAmixobj) throw std::runtime_error("SPAmix object is not initialized.");
    ctx.spamix.reset(new SPAmix::SPAmixClass(*ptr_gSPAmixobj));
  } else if (method == "SPAmixPlus") {
    if (!ptr_gSPAmixPlusobj) throw std::runtime_error("SPAmixPlus object is not initialized.");
    ctx.spamixPlus.reset(new SPAmixPlus::SPAmixPlusClass(*ptr_gSPAmixPlusobj));
  } else if (method == "SPAGRM") {
    if (!ptr_gSPAGRMobj) throw std::runtime_error("SPAGRM object is not initialized.");
    ctx.spagrm.reset(new SPAGRM::SPAGRMClass(*ptr_gSPAGRMobj));
  } else if (method == "SAGELD") {
    if (!ptr_gSAGELDobj) throw std::runtime_error("SAGELD object is not initialized.");
    ctx.sageld.reset(new SAGELD::SAGELDClass(*ptr_gSAGELDobj));
  } else if (method == "WtCoxG") {
    if (!ptr_gWtCoxGobj) throw std::runtime_error("WtCoxG object is not initialized.");
    ctx.wtcoxg.reset(new WtCoxG::WtCoxGClass(*ptr_gWtCoxGobj));
  } else if (method == "SPAsqr") {
    if (!ptr_gSPAsqrobj) throw std::runtime_error("SPAsqr object is not initialized.");
    ctx.spasqr.reset(new SPAsqr::SPAsqrClass(*ptr_gSPAsqrobj));
  } else if (method == "LEAF") {
    if (!ptr_gLEAFobj) throw std::runtime_error("LEAF object is not initialized.");
    ctx.leaf.reset(new LEAF::LEAFClass(*ptr_gLEAFobj));
  }
  return ctx;
}


// ---- Engine: mainMarkerChunksCore ----

void mainMarkerChunksCore(
  const std::string& method,
  const std::vector<std::vector<uint64_t>>& chunkMarkers,
  const std::string& outputFile,
  unsigned int nThreads,
  const std::string& impute_method,
  double missingCutoff,
  double minMafMarker,
  double minMacMarker,
  const ExtraParams& extra
) {
  if (chunkMarkers.empty())
    throw std::runtime_error("chunkMarkers is empty.");
  if (nThreads < 1)
    throw std::runtime_error("nThreads should be >= 1.");

  const ReaderConfig& reader = extra.reader;
  if (reader.genoType.empty())
    throw std::runtime_error("Reader configuration is required in Marker4 core.");

  int nPheno = 1;
  if (method == "SPAmix") nPheno = ptr_gSPAmixobj->getNpheno();
  if (method == "SPAmixPlus") nPheno = ptr_gSPAmixPlusobj->getNpheno();
  if (method == "SPAsqr") nPheno = ptr_gSPAsqrobj->get_ntaus();

  // Extract method-specific constants from class objects (not chunk-specific)
  const std::string sageldMethod = (method == "SAGELD" && ptr_gSAGELDobj) ? ptr_gSAGELDobj->getMethod() : "SAGELD";
  const std::vector<double> spasqrTaus = (method == "SPAsqr" && ptr_gSPAsqrobj) ? ptr_gSPAsqrobj->getTaus() : std::vector<double>{};
  const int leafNcluster = (method == "LEAF" && ptr_gLEAFobj) ? ptr_gLEAFobj->get_Ncluster() : 0;

  const std::string header = getHeader(method, sageldMethod, spasqrTaus, nPheno, leafNcluster);

  std::vector<std::string> chunkOutput(chunkMarkers.size());
  std::vector<char> chunkReady(chunkMarkers.size(), 0);
  std::atomic<size_t> nextChunk(0);

  std::exception_ptr workerError = nullptr;
  std::mutex errorMutex;
  std::mutex writeMutex;
  std::mutex logMutex;
  std::condition_variable writeCv;
  bool stopWriter = false;


  std::thread writerThread([&]() {
    const bool useGzip = outputFile.size() >= 3 &&
      outputFile.compare(outputFile.size() - 3, 3, ".gz") == 0;

    gzFile gz = nullptr;
    std::ofstream plainOut;

    if (useGzip) {
      gz = gzopen(outputFile.c_str(), "wb");
      if (!gz) {
        std::lock_guard<std::mutex> lock(errorMutex);
        workerError = std::make_exception_ptr(std::runtime_error("Cannot open gzip output file: " + outputFile));
        stopWriter = true;
        writeCv.notify_all();
        return;
      }
    } else {
      plainOut.open(outputFile.c_str());
      if (!plainOut.is_open()) {
        std::lock_guard<std::mutex> lock(errorMutex);
        workerError = std::make_exception_ptr(std::runtime_error("Cannot open output file: " + outputFile));
        stopWriter = true;
        writeCv.notify_all();
        return;
      }
    }

    auto writeStr = [&](const std::string& s) {
      if (useGzip) gzwrite(gz, s.data(), static_cast<unsigned>(s.size()));
      else plainOut << s;
    };

    writeStr(header + "\n");

    for (size_t i = 0; i < chunkOutput.size(); ++i) {
      std::unique_lock<std::mutex> lk(writeMutex);
      writeCv.wait(lk, [&]() { return chunkReady[i] || stopWriter; });
      if (!chunkReady[i]) break;
      writeStr(chunkOutput[i]);
      {
        std::lock_guard<std::mutex> lg(logMutex);
        std::fprintf(stderr, "[INFO] Writing finished: chunk %zu/%zu\n", i + 1, chunkOutput.size());
        std::fflush(stderr);
      }
    }

    if (useGzip) gzclose(gz);
    else plainOut.close();
  });


  auto workerFn = [&]() {
    try {
      ThreadContext ctx = makeThreadContext(method);

      if (reader.genoType != "PLINK")
        throw std::runtime_error("Unsupported reader genoType in Marker4 core: " + reader.genoType);
      std::unique_ptr<PLINK4::PlinkReader> plinkReader(new PLINK4::PlinkReader(
        reader.bimFile, reader.famFile, reader.bedFile,
        reader.sampleInModel, reader.alleleOrder
      ));

      while (true) {
        size_t cidx = nextChunk.fetch_add(1);
        if (cidx >= chunkMarkers.size()) break;

        const auto& gIdx = chunkMarkers[cidx];
        std::ostringstream out;

        if (!gIdx.empty()) plinkReader->beginSequentialBlock(gIdx.front());


        std::vector<WtRow> wtChunkRows;
        if (method == "WtCoxG") {
          wtChunkRows.reserve(gIdx.size());
          std::vector<double> AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext;
          AF_ref.reserve(gIdx.size()); AN_ref.reserve(gIdx.size());
          TPR.reserve(gIdx.size()); sigma2.reserve(gIdx.size());
          pvalue_bat.reserve(gIdx.size()); w_ext.reserve(gIdx.size());
          var_w0.reserve(gIdx.size()); var_int.reserve(gIdx.size()); var_ext.reserve(gIdx.size());

          for (auto gi : gIdx) {
            auto it = extra.wtMap.find(gi);
            wtChunkRows.push_back(it == extra.wtMap.end()
              ? WtRow{NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN} : it->second);
          }
          for (const auto& r : wtChunkRows) {
            AF_ref.push_back(r.AF_ref); AN_ref.push_back(r.AN_ref);
            TPR.push_back(r.TPR); sigma2.push_back(r.sigma2);
            pvalue_bat.push_back(r.pvalue_bat); w_ext.push_back(r.w_ext);
            var_w0.push_back(r.var_ratio_w0); var_int.push_back(r.var_ratio_int);
            var_ext.push_back(r.var_ratio_ext);
          }
          ctx.wtcoxg->updateMarkerInfo(AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext);
        }


        std::vector<std::vector<LeafRow>> leafChunkRows;
        if (method == "LEAF") {
          leafChunkRows.resize(leafNcluster);
          for (int c = 0; c < leafNcluster; ++c) {
            auto& rows = leafChunkRows[c];
            rows.reserve(gIdx.size());
            std::vector<double> AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext;
            AF_ref.reserve(gIdx.size()); AN_ref.reserve(gIdx.size());
            TPR.reserve(gIdx.size()); sigma2.reserve(gIdx.size());
            pvalue_bat.reserve(gIdx.size()); w_ext.reserve(gIdx.size());
            var_w0.reserve(gIdx.size()); var_int.reserve(gIdx.size()); var_ext.reserve(gIdx.size());

            for (auto gi : gIdx) {
              auto it = extra.leafMaps[c].find(gi);
              rows.push_back(it == extra.leafMaps[c].end()
                ? LeafRow{NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN} : it->second);
            }
            for (const auto& r : rows) {
              AF_ref.push_back(r.AF_ref); AN_ref.push_back(r.AN_ref);
              TPR.push_back(r.TPR); sigma2.push_back(r.sigma2);
              pvalue_bat.push_back(r.pvalue_bat); w_ext.push_back(r.w_ext);
              var_w0.push_back(r.var_ratio_w0); var_int.push_back(r.var_ratio_int);
              var_ext.push_back(r.var_ratio_ext);
            }
            ctx.leaf->updateMarkerInfo(c, AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext);
          }
        }


        for (size_t i = 0; i < gIdx.size(); ++i) {
          double altFreq = NAN, altCounts = NAN, missingRate = NAN, imputeInfo = NAN;
          std::vector<uint32_t> indexForMissing, indexForNonZero;
          std::string chr, ref, alt, marker;
          uint32_t pd = 0;

          arma::vec GVec = plinkReader->getOneMarker(gIdx[i], ref, alt, marker, pd, chr,
            altFreq, altCounts, missingRate, imputeInfo, true, indexForMissing, false, indexForNonZero, true);

          const double n = static_cast<double>(GVec.size());
          const double maf = std::min(altFreq, 1.0 - altFreq);
          const double mac = 2.0 * maf * n * (1.0 - missingRate);
          bool passQC = !((missingRate > missingCutoff) || (maf < minMafMarker) || (mac < minMacMarker));

          bool flip = false;
          if (passQC) {
            int nMissing = indexForMissing.size();
            double imputeG = 0;
            if (impute_method == "mean"){
              imputeG = 2 * altFreq;
            }

            if (impute_method == "minor"){
              if (altFreq > 0.5){
                imputeG = 2;
              }else{
                imputeG = 0;
              }
            }

            for (int i = 0; i < nMissing; i++){
              uint32_t index = indexForMissing.at(i);
              GVec.at(index) = imputeG;
            }

            if (method != "WtCoxG") {
              if (altFreq > 0.5){
                GVec = 2 - GVec;
                flip = true;
              }
            }
          }

          if (method == "SPAmix" || method == "SPAmixPlus") {
            std::vector<double> p(nPheno, NAN), z(nPheno, NAN);
            if (passQC) {
              if (method == "SPAmix") {
                ctx.spamix->getMarkerPval(GVec, altFreq);
                arma::vec pT = ctx.spamix->getpvalVec(), zT = ctx.spamix->getzScoreVec();
                for (int j = 0; j < nPheno; ++j) { p[j] = pT[j]; z[j] = zT[j]; }
              } else {
                ctx.spamixPlus->getMarkerPval(GVec, altFreq);
                arma::vec pT = ctx.spamixPlus->getpvalVec(), zT = ctx.spamixPlus->getzScoreVec();
                for (int j = 0; j < nPheno; ++j) { p[j] = pT[j]; z[j] = zT[j]; }
              }
            }
            for (int j = 0; j < nPheno; ++j) {
              std::ostringstream row;
              appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
              appendCell(row, "pheno_" + std::to_string(j + 1));
              appendCell(row, p[j]); appendCell(row, z[j]);
              out << row.str() << '\n';
            }
            continue;
          }


          if (method == "POLMM") {
            double Beta = NAN, seBeta = NAN, pval = NAN, zScore = NAN;
            if (passQC) {
              ctx.polmm->getMarkerPval(GVec, Beta, seBeta, pval, altFreq, zScore);
              Beta *= (1.0 - 2.0 * static_cast<double>(flip));
            }
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, pval);
            appendCell(row, Beta); appendCell(row, seBeta); appendCell(row, zScore);
            out << row.str() << '\n';
            continue;
          }


          if (method == "SPACox") {
            double pval = NAN, zScore = NAN;
            if (passQC) pval = ctx.spacox->getMarkerPval(GVec, altFreq, zScore);
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, pval); appendCell(row, zScore);
            out << row.str() << '\n';
            continue;
          }


          if (method == "SPAGRM") {
            double pval = NAN, zScore = NAN, hwepval = NAN;
            if (passQC) pval = ctx.spagrm->getMarkerPval(GVec, altFreq, zScore, hwepval);
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, zScore);
            appendCell(row, pval); appendCell(row, hwepval);
            out << row.str() << '\n';
            continue;
          }


          if (method == "SAGELD") {
            double hwepval = NAN;
            double zG = NAN, zGxE = NAN, pG = NAN, pGxE = NAN;
            double bG = NAN, bGxE = NAN, seG = NAN, seGxE = NAN;
            if (passQC) {
              ctx.sageld->getMarkerPval(GVec, altFreq, hwepval);
              arma::vec pT = ctx.sageld->getpvalVec(), zT = ctx.sageld->getzScoreVec();
              arma::vec bT = ctx.sageld->getBetaVec(), seT = ctx.sageld->getseBetaVec();
              zG = zT[0]; zGxE = zT[1]; pG = pT[0]; pGxE = pT[1];
              bG = bT[0] * (1.0 - 2.0 * static_cast<double>(flip));
              bGxE = bT[1] * (1.0 - 2.0 * static_cast<double>(flip));
              seG = seT[0]; seGxE = seT[1];
            }
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, sageldMethod);
            if (sageldMethod == "GALLOP") {
              appendCell(row, bG); appendCell(row, bGxE);
              appendCell(row, seG); appendCell(row, seGxE);
              appendCell(row, pG); appendCell(row, pGxE);
            } else {
              appendCell(row, zG); appendCell(row, zGxE);
              appendCell(row, pG); appendCell(row, pGxE);
            }
            appendCell(row, hwepval);
            out << row.str() << '\n';
            continue;
          }


          if (method == "WtCoxG") {
            double pExt = NAN, pNoext = NAN, zExt = NAN, zNoext = NAN;
            if (passQC) {
              arma::vec pT = ctx.wtcoxg->getpvalVec(GVec, static_cast<int>(i));
              arma::vec zT = ctx.wtcoxg->getZScoreVec();
              pExt = pT[0]; pNoext = pT[1]; zExt = zT[0]; zNoext = zT[1];
            }
            const WtRow& wr = wtChunkRows[i];
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, pExt); appendCell(row, pNoext);
            appendCell(row, zExt); appendCell(row, zNoext);
            appendCell(row, wr.AF_ref); appendCell(row, wr.AN_ref);
            appendCell(row, wr.pvalue_bat);
            out << row.str() << '\n';
            continue;
          }


          if (method == "SPAsqr") {
            double hwepval = NAN;
            std::vector<double> z(spasqrTaus.size(), NAN), p(spasqrTaus.size(), NAN);
            if (passQC) {
              arma::vec zT;
              arma::vec pT = ctx.spasqr->getMarkerPval(GVec, altFreq, zT, hwepval);
              for (size_t j = 0; j < p.size(); ++j) { p[j] = pT[j]; z[j] = zT[j]; }
            }
            double pCCT = cctPvalue(p);
            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, hwepval);
            for (double v : z) appendCell(row, v);
            for (double v : p) appendCell(row, v);
            appendCell(row, pCCT);
            out << row.str() << '\n';
            continue;
          }


          if (method == "LEAF") {
            std::vector<double> pExt(leafNcluster, NAN), pNoext(leafNcluster, NAN);
            std::vector<double> sExt(leafNcluster, NAN), sNoext(leafNcluster, NAN);

            if (passQC) {
              arma::vec all = ctx.leaf->getMarkerZSP(GVec, static_cast<int>(i));
              const int nOut = 2 * leafNcluster;
              for (int c = 0; c < leafNcluster; ++c) {
                sExt[c] = all[nOut + 2*c]; sNoext[c] = all[nOut + 2*c + 1];
                pExt[c] = all[2*nOut + 2*c]; pNoext[c] = all[2*nOut + 2*c + 1];
              }
            }

            auto metaP = [&](const std::vector<double>& scores, const std::vector<double>& pvals) {
              double sumScore = 0, sumVar = 0;
              boost::math::chi_squared dist(1.0);
              for (size_t c = 0; c < scores.size(); ++c) {
                if (std::isnan(scores[c]) || std::isnan(pvals[c]) || pvals[c] <= 0 || pvals[c] >= 1) continue;
                double chisq = boost::math::quantile(dist, 1.0 - pvals[c]);
                if (chisq < 1e-30) chisq = 1e-30;
                double var = (scores[c] * scores[c]) / chisq;
                if (std::isnan(var)) var = 0;
                sumScore += scores[c]; sumVar += var;
              }
              if (sumVar <= 0) return std::numeric_limits<double>::quiet_NaN();
              double z = sumScore / std::sqrt(sumVar);
              return std::erfc(std::fabs(z) / std::sqrt(2.0));
            };

            std::ostringstream row;
            appendMeta(row, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate);
            appendCell(row, metaP(sExt, pExt)); appendCell(row, metaP(sNoext, pNoext));
            for (int c = 0; c < leafNcluster; ++c) {
              appendCell(row, pExt[c]); appendCell(row, pNoext[c]);
              appendCell(row, leafChunkRows[c][i].pvalue_bat);
            }
            out << row.str() << '\n';
            continue;
          }
        }

        {
          std::lock_guard<std::mutex> lk(writeMutex);
          chunkOutput[cidx] = out.str();
          chunkReady[cidx] = 1;
        } {
          std::lock_guard<std::mutex> lg(logMutex);
          std::fprintf(stderr, "[INFO] Calculation finished: chunk %zu/%zu\n", cidx + 1, chunkMarkers.size());
          std::fflush(stderr);
        }
        writeCv.notify_all();
      }
    } catch (...) {
      { std::lock_guard<std::mutex> lk(errorMutex); if (!workerError) workerError = std::current_exception(); }
      { std::lock_guard<std::mutex> lk(writeMutex); stopWriter = true; }
      writeCv.notify_all();
    }
  };

  nThreads = std::min<unsigned int>(nThreads, static_cast<unsigned int>(chunkMarkers.size()));
  const bool runInline = (nThreads == 1 && (method == "WtCoxG" || method == "LEAF"));
  if (runInline) {
    workerFn();
  } else {
    std::vector<std::thread> workers;
    workers.reserve(nThreads);
    for (unsigned int t = 0; t < nThreads; ++t) workers.emplace_back(workerFn);
    for (auto& th : workers) th.join();
  }

  { std::lock_guard<std::mutex> lk(writeMutex); stopWriter = true; }
  writeCv.notify_all();
  writerThread.join();

  if (workerError) std::rethrow_exception(workerError);
}

} // namespace Marker4

// ---- Bridge functions (Rcpp entry points + support utilities) ----
// Included here so they are compiled as part of this translation unit.
// Marker4Bridge.h must be included from exactly one .cpp file.
#include "Marker4Bridge.h"
