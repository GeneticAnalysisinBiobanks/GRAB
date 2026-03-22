// mtMain.cpp -- Thread pool, per-chunk worker, sequential writer
//
// Functions by role:
//
//   [Inline helpers]
//     numToChars         — format a double directly into a char buffer
//     makeNaSuffix       — build precomputed NA suffix for fail-QC rows
//     formatLine         — format one output row (meta + numeric values)
//     formatLineNA       — format one fail-QC output row (meta + NA suffix)
//     formatLineWithPrefix — format one output row with a string prefix column
//
//   [File-local structs]
//     ThreadContext      — per-worker copies of the active statistical method object
//
//   [Output header]       (static)
//     getHeader           — produce the TSV column header for a given method
//
//   [Worker context]      (static)
//     makeThreadContext   — clone the global method object into a per-thread ThreadContext
//
//   [Engine]              (mtMarker namespace, external linkage)
//     mainMarkerChunksCore — coordinate the PLINK cursor, worker pool, and writer thread

#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <exception>
#include <cstdio>
#include <chrono>
#include <ctime>
#include <cstdarg>
#include <RcppArmadillo.h>
#include <zlib.h>

#include "mtPLINK.h"
#include "mtPOLMM.h"
#include "mtWtCoxG.h"
#include "mtLEAF.h"
#include "mtSPAGRM.h"
#include "mtSAGELD.h"
#include "mtSPAsqr.h"
#include "mtSPACox.h"
#include "mtSPAmix.h"
#include "mtSPAmixPlus.h"


// ---- Global method pointers (one per active statistical method) ----

mtPOLMMClass*           ptr_gPOLMMobj      = nullptr;
mtWtCoxGClass*         ptr_gWtCoxGobj     = nullptr;
mtLEAFClass*             ptr_gLEAFobj       = nullptr;
mtSPAGRMClass*         ptr_gSPAGRMobj     = nullptr;
mtSAGELDClass*         ptr_gSAGELDobj     = nullptr;
mtSPAsqrClass*         ptr_gSPAsqrobj     = nullptr;
mtSPACoxClass*         ptr_gSPACoxobj     = nullptr;
mtSPAmixClass*         ptr_gSPAmixobj     = nullptr;
mtSPAmixPlusClass* ptr_gSPAmixPlusobj = nullptr;


namespace {

// Timestamped info message (mirrors R's .message helper)
void infoMsg(const char* fmt, ...) {
  auto now = std::chrono::system_clock::now();
  auto t = std::chrono::system_clock::to_time_t(now);
  char ts[20];
  std::strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
  char buf[512];
  va_list args;
  va_start(args, fmt);
  std::vsnprintf(buf, sizeof(buf), fmt, args);
  va_end(args);
  Rprintf("[INFO] %s %s\n", ts, buf);
}

// Build a string of nCols tab-separated "NA" values (for fail-QC rows).
// Determines nResultCols based on method using global method pointers.
std::string makeNaSuffix(const std::string& method) {
  int nResultCols = 0;
  if      (method == "POLMM")      nResultCols = ptr_gPOLMMobj->resultSize();
  else if (method == "SPACox")     nResultCols = ptr_gSPACoxobj->resultSize();
  else if (method == "SPAmix")     nResultCols = ptr_gSPAmixobj->resultSize();
  else if (method == "SPAmixPlus") nResultCols = ptr_gSPAmixPlusobj->resultSize();
  else if (method == "SPAGRM")    nResultCols = ptr_gSPAGRMobj->resultSize();
  else if (method == "SAGELD")    nResultCols = ptr_gSAGELDobj->resultSize();
  else if (method == "WtCoxG")    nResultCols = ptr_gWtCoxGobj->resultSize();
  else if (method == "SPAsqr")    nResultCols = ptr_gSPAsqrobj->resultSize();
  else if (method == "LEAF")      nResultCols = ptr_gLEAFobj->resultSize();
  
  if (nResultCols <= 0) return {};
  // Each column: "\tNA" = 3 chars
  std::string s;
  s.reserve(3 * nResultCols);
  for (int i = 0; i < nResultCols; ++i) { s += '\t'; s += 'N'; s += 'A'; }
  return s;
}

// Format a double as a short string: "NA", "Inf", "-Inf", or %.6g.
// Writes directly into a char buffer and returns the length.
int numToChars(char* buf, double x) {
  if (std::isnan(x)) { buf[0]='N'; buf[1]='A'; return 2; }
  if (std::isinf(x)) {
    if (x > 0) { buf[0]='I'; buf[1]='n'; buf[2]='f'; return 3; }
    else { buf[0]='-'; buf[1]='I'; buf[2]='n'; buf[3]='f'; return 4; }
  }
  return std::snprintf(buf, 32, "%.6g", x);
}

// Format one output line: meta columns + numeric result values.
// Uses pre-reserved std::string for efficiency.
void formatLine(
  std::string& out,
  const std::string& marker,
  const std::string& chr,
  uint32_t pos,
  const std::string& ref,
  const std::string& alt,
  double altCounts,
  double altFreq,
  double missingRate,
  const std::vector<double>& vals
) {
  std::string s;
  s.reserve(512);
  char buf[32];
  int n;

  // Meta columns: Marker, Chr, Pos, Ref, Alt, AltCount, AltFreq, MissRate
  s += marker; s += '\t';
  s += chr;    s += '\t';
  n = std::snprintf(buf, sizeof(buf), "%u", pos);
  s.append(buf, n); s += '\t';
  s += ref;    s += '\t';
  s += alt;    s += '\t';
  n = numToChars(buf, altCounts); s.append(buf, n); s += '\t';
  n = numToChars(buf, altFreq);   s.append(buf, n); s += '\t';
  n = numToChars(buf, missingRate); s.append(buf, n);

  // Value columns
  for (double v : vals) {
    s += '\t';
    n = numToChars(buf, v);
    s.append(buf, n);
  }
  s += '\n';
  out += s;
}

// Format one output line for fail-QC markers using precomputed NA suffix.
void formatLineNA(
  std::string& out,
  const std::string& marker,
  const std::string& chr,
  uint32_t pos,
  const std::string& ref,
  const std::string& alt,
  double altCounts,
  double altFreq,
  double missingRate,
  const std::string& naSuffix
) {
  std::string s;
  s.reserve(512);
  char buf[32];
  int n;

  s += marker; s += '\t';
  s += chr;    s += '\t';
  n = std::snprintf(buf, sizeof(buf), "%u", pos);
  s.append(buf, n); s += '\t';
  s += ref;    s += '\t';
  s += alt;    s += '\t';
  n = numToChars(buf, altCounts); s.append(buf, n); s += '\t';
  n = numToChars(buf, altFreq);   s.append(buf, n); s += '\t';
  n = numToChars(buf, missingRate); s.append(buf, n);
  s += naSuffix;
  s += '\n';
  out += s;
}

// ---- File-local structs ----

// Per-worker independent copies of the active method object.
struct ThreadContext {
  std::unique_ptr<mtPOLMMClass>           polmm;
  std::unique_ptr<mtSPACoxClass>         spacox;
  std::unique_ptr<mtSPAmixClass>         spamix;
  std::unique_ptr<mtSPAmixPlusClass> spamixPlus;
  std::unique_ptr<mtSPAGRMClass>         spagrm;
  std::unique_ptr<mtSAGELDClass>         sageld;
  std::unique_ptr<mtWtCoxGClass>         wtcoxg;
  std::unique_ptr<mtSPAsqrClass>         spasqr;
  std::unique_ptr<mtLEAFClass>             leaf;
};


// ---- Output header ----

const char* META_HEADER = "Marker\tChr\tPos\tRef\tAlt\tAltCount\tAltFreq\tMissRate";

// Build the TSV column header line by dispatching to the active method object.
std::string getHeader(const std::string& method) {
  std::string cols;
  if (method == "POLMM")          cols = ptr_gPOLMMobj->getHeaderColumns();
  else if (method == "SPACox")    cols = ptr_gSPACoxobj->getHeaderColumns();
  else if (method == "SPAmix")    cols = ptr_gSPAmixobj->getHeaderColumns();
  else if (method == "SPAmixPlus") cols = ptr_gSPAmixPlusobj->getHeaderColumns();
  else if (method == "SPAGRM")   cols = ptr_gSPAGRMobj->getHeaderColumns();
  else if (method == "SAGELD")   cols = ptr_gSAGELDobj->getHeaderColumns();
  else if (method == "WtCoxG")   cols = ptr_gWtCoxGobj->getHeaderColumns();
  else if (method == "SPAsqr")   cols = ptr_gSPAsqrobj->getHeaderColumns();
  else if (method == "LEAF")     cols = ptr_gLEAFobj->getHeaderColumns();
  else throw std::runtime_error("Unsupported method in mtMarker: " + method);
  return std::string(META_HEADER) + cols;
}


// ---- Worker context ----

// Clone the global method object into a per-thread ThreadContext.
// Each method class stores per-marker mutable state (result vectors, scratch),
// so per-thread copies are required for thread safety.  The copy happens once
// per thread, not per chunk, so the cost is amortized.
ThreadContext makeThreadContext(const std::string& method) {
  ThreadContext ctx;
  if (method == "POLMM") {
    if (!ptr_gPOLMMobj) throw std::runtime_error("POLMM object is not initialized.");
    ctx.polmm.reset(new mtPOLMMClass(*ptr_gPOLMMobj));
  } else if (method == "SPACox") {
    if (!ptr_gSPACoxobj) throw std::runtime_error("SPACox object is not initialized.");
    ctx.spacox.reset(new mtSPACoxClass(*ptr_gSPACoxobj));
  } else if (method == "SPAmix") {
    if (!ptr_gSPAmixobj) throw std::runtime_error("SPAmix object is not initialized.");
    ctx.spamix.reset(new mtSPAmixClass(*ptr_gSPAmixobj));
  } else if (method == "SPAmixPlus") {
    if (!ptr_gSPAmixPlusobj) throw std::runtime_error("SPAmixPlus object is not initialized.");
    ctx.spamixPlus.reset(new mtSPAmixPlusClass(*ptr_gSPAmixPlusobj));
  } else if (method == "SPAGRM") {
    if (!ptr_gSPAGRMobj) throw std::runtime_error("SPAGRM object is not initialized.");
    ctx.spagrm.reset(new mtSPAGRMClass(*ptr_gSPAGRMobj));
  } else if (method == "SAGELD") {
    if (!ptr_gSAGELDobj) throw std::runtime_error("SAGELD object is not initialized.");
    ctx.sageld.reset(new mtSAGELDClass(*ptr_gSAGELDobj));
  } else if (method == "WtCoxG") {
    if (!ptr_gWtCoxGobj) throw std::runtime_error("WtCoxG object is not initialized.");
    ctx.wtcoxg.reset(new mtWtCoxGClass(*ptr_gWtCoxGobj));
  } else if (method == "SPAsqr") {
    if (!ptr_gSPAsqrobj) throw std::runtime_error("SPAsqr object is not initialized.");
    ctx.spasqr.reset(new mtSPAsqrClass(*ptr_gSPAsqrobj));
  } else if (method == "LEAF") {
    if (!ptr_gLEAFobj) throw std::runtime_error("LEAF object is not initialized.");
    ctx.leaf.reset(new mtLEAFClass(*ptr_gLEAFobj));
  }
  return ctx;
}

} // namespace


// ---- The Engine ----

void mtMarkerEngine(
  const std::string method,
  const std::string bedFile,
  const std::string bimFile,
  const std::string famFile,
  const std::string outputFile,
  const std::vector<std::string>& subjData,
  const std::string AlleleOrder,
  const int nMarkersEachChunk,
  const int nthreads,
  const std::string impute_method,
  const double missing_cutoff,
  const double min_maf_marker,
  const double min_mac_marker,
  const std::string IDsToIncludeFile,
  const std::string RangesToIncludeFile,
  const std::string IDsToExcludeFile,
  const std::string RangesToExcludeFile
) {

  mtPLINK::PlinkData plinkData(bedFile, bimFile, famFile, subjData, AlleleOrder);

  mtPLINK::MarkerFilterConfig filterConfig;
  filterConfig.IDsToIncludeFile    = IDsToIncludeFile;
  filterConfig.RangesToIncludeFile = RangesToIncludeFile;
  filterConfig.IDsToExcludeFile    = IDsToExcludeFile;
  filterConfig.RangesToExcludeFile = RangesToExcludeFile;

  std::vector<mtPLINK::MarkerInfo> markerInfo = plinkData.getFilteredMarkers(filterConfig);
  if (markerInfo.empty()) {
    throw std::runtime_error("No markers remain after PLINK marker filtering.");
  }

  std::vector<std::vector<uint64_t>> chunkIndices = mtPLINK::PlinkData::buildChunks(markerInfo, nMarkersEachChunk);
  const int effective_nthreads = std::min(nthreads, static_cast<int>(chunkIndices.size()));

  infoMsg("Number of subjects in the input file: %u", plinkData.nSubjInFile());
  infoMsg("Number of subjects to test: %u", plinkData.nSubjUsed());
  infoMsg("Number of markers in the input file: %u", plinkData.nMarkers());
  infoMsg("Number of markers to test: %zu", markerInfo.size());
  infoMsg("Number of markers in each chunk: %d", nMarkersEachChunk);
  infoMsg("Number of chunks for all markers: %zu", chunkIndices.size());
  infoMsg("Number of threads: %d", effective_nthreads);

  const std::string header = getHeader(method);
  std::vector<std::string> chunkOutput(chunkIndices.size());
  std::vector<char> chunkReady(chunkIndices.size(), 0);
  std::atomic<size_t> nextChunk(0);

  std::exception_ptr workerError = nullptr;
  std::mutex errorMutex;
  std::mutex writeMutex;
  std::condition_variable writeCv;
  std::atomic<bool> stopWriter{false};

  std::thread writerThread([&]() {
    const bool useGzip = outputFile.size() > 3 && outputFile.compare(outputFile.size() - 3, 3, ".gz") == 0;
    gzFile gz = nullptr;
    std::ofstream plainOut;

    if (useGzip) {
      gz = gzopen(outputFile.c_str(), "wb");
      if (!gz) {
        {
          std::lock_guard<std::mutex> lock(errorMutex);
          workerError = std::make_exception_ptr(std::runtime_error("Cannot open gzip output file: " + outputFile)); 
        }
        stopWriter.store(true);
        return;
      }
    } else {
      plainOut.open(outputFile.c_str());
      if (!plainOut.is_open()) {
        {
          std::lock_guard<std::mutex> lock(errorMutex);
          workerError = std::make_exception_ptr(std::runtime_error("Cannot open output file: " + outputFile));
        }
        stopWriter.store(true);
        return;
      }
    }

    auto writeStr = [&](const std::string& s) {
      if (useGzip) {
        gzwrite(gz, s.data(), static_cast<unsigned>(s.size()));
      } else {
        plainOut << s;
      }
    };

    writeStr(header + "\n");

    for (size_t i = 0; i < chunkOutput.size(); ++i) {
      std::unique_lock<std::mutex> lk(writeMutex);
      writeCv.wait(lk, [&]() { return chunkReady[i] || stopWriter; });
      if (!chunkReady[i]) {
        break;
      }
      writeStr(chunkOutput[i]);
      infoMsg("Writing finished: chunk %zu/%zu", i + 1, chunkOutput.size());
    }

    if (useGzip) {
      gzclose(gz);
    } else {
      plainOut.close();
    }
  });


  auto workerFn = [&]() {
    try {
      // Per-thread: copy method object (once) + open .bed file handle
      ThreadContext ctx = makeThreadContext(method);
      mtPLINK::PlinkCursor cursor(plinkData);

      // Precompute NA suffix for fail-QC lines (method-dependent column count)
      const std::string naSuffix = makeNaSuffix(method);

      while (true) {
        size_t cidx = nextChunk.fetch_add(1);
        if (cidx >= chunkIndices.size()) break;

        const auto& gIdx = chunkIndices[cidx];
        std::string out;
        out.reserve(gIdx.size() * 256);

        if (!gIdx.empty()) cursor.beginSequentialBlock(gIdx.front());

        if (method == "WtCoxG") ctx.wtcoxg->prepareChunk(gIdx);
        if (method == "LEAF")   ctx.leaf->prepareChunk(gIdx);

        for (size_t i = 0; i < gIdx.size(); ++i) {
          double altFreq = NAN, altCounts = NAN, missingRate = NAN;
          std::vector<uint32_t> indexForMissing;

          arma::vec GVec = cursor.getGenotypes(gIdx[i],
            altFreq, altCounts, missingRate, indexForMissing);

          const std::string& chr    = plinkData.chr(gIdx[i]);
          const std::string& ref    = plinkData.ref(gIdx[i]);
          const std::string& alt    = plinkData.alt(gIdx[i]);
          const std::string& marker = plinkData.markerId(gIdx[i]);
          uint32_t pd               = plinkData.pos(gIdx[i]);

          const double n = static_cast<double>(GVec.size());
          const double maf = std::min(altFreq, 1.0 - altFreq);
          const double mac = 2.0 * maf * n * (1.0 - missingRate);
          bool passQC = !((missingRate > missing_cutoff) || (maf < min_maf_marker) || (mac < min_mac_marker));

          bool flip = false;
          if (passQC) {
            int nMissing = indexForMissing.size();
            double imputeG = 2 * altFreq;
            for (int j = 0; j < nMissing; j++) {
              GVec.at(indexForMissing.at(j)) = imputeG;
            }
            if (method != "WtCoxG") {
              if (altFreq > 0.5) {
                GVec = 2 - GVec;
                flip = true;
              }
            }
          }

          // --- SPAmix ---
          if (method == "SPAmix") {
            if (!passQC) { formatLineNA(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, naSuffix); continue; }
            std::vector<double> rv = ctx.spamix->getResultVec(GVec, altFreq);
            formatLine(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, rv);
            continue;
          }

          // --- SPAmixPlus ---
          if (method == "SPAmixPlus") {
            if (!passQC) { formatLineNA(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, naSuffix); continue; }
            std::vector<double> rv = ctx.spamixPlus->getResultVec(GVec, altFreq);
            formatLine(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, rv);
            continue;
          }

          // --- POLMM ---
          if (method == "POLMM") {
            if (!passQC) { formatLineNA(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, naSuffix); continue; }
            std::vector<double> rv = ctx.polmm->getResultVec(GVec, altFreq);
            rv[1] *= (1.0 - 2.0 * static_cast<double>(flip)); // flip beta
            formatLine(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, rv);
            continue;
          }

          // --- SPACox ---
          if (method == "SPACox") {
            if (!passQC) { formatLineNA(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, naSuffix); continue; }
            std::vector<double> rv = ctx.spacox->getResultVec(GVec, altFreq);
            formatLine(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, rv);
            continue;
          }

          // --- SPAGRM ---
          if (method == "SPAGRM") {
            if (!passQC) { formatLineNA(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, naSuffix); continue; }
            std::vector<double> rv = ctx.spagrm->getResultVec(GVec, altFreq);
            formatLine(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, rv);
            continue;
          }

          // --- SAGELD ---
          if (method == "SAGELD") {
            if (!passQC) { formatLineNA(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, naSuffix); continue; }
            std::vector<double> rv = ctx.sageld->getResultVec(GVec, altFreq);
            formatLine(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, rv);
            continue;
          }

          // --- WtCoxG ---
          if (method == "WtCoxG") {
            if (!passQC) { formatLineNA(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, naSuffix); continue; }
            std::vector<double> rv = ctx.wtcoxg->getResultVec(GVec, static_cast<int>(i));
            formatLine(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, rv);
            continue;
          }

          // --- SPAsqr ---
          if (method == "SPAsqr") {
            if (!passQC) { formatLineNA(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, naSuffix); continue; }
            std::vector<double> rv = ctx.spasqr->getResultVec(GVec, altFreq);
            formatLine(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, rv);
            continue;
          }

          // --- LEAF ---
          if (method == "LEAF") {
            if (!passQC) { formatLineNA(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, naSuffix); continue; }
            std::vector<double> rv = ctx.leaf->getResultVec(GVec, static_cast<int>(i));
            formatLine(out, marker, chr, pd, ref, alt, altCounts, altFreq, missingRate, rv);
            continue;
          }
        }

        {
          std::lock_guard<std::mutex> lk(writeMutex);
          chunkOutput[cidx] = std::move(out);
          chunkReady[cidx] = 1;
        }
        infoMsg("Calculation finished: chunk %zu/%zu", cidx + 1, chunkIndices.size());
        writeCv.notify_all();
      }
    } catch (...) {
      { std::lock_guard<std::mutex> lk(errorMutex); if (!workerError) workerError = std::current_exception(); }
      { std::lock_guard<std::mutex> lk(writeMutex); stopWriter = true; }
      writeCv.notify_all();
    }
  };


  const bool runInline = (effective_nthreads == 1);
  if (runInline) {
    workerFn();
  } else {
    std::vector<std::thread> workers;
    workers.reserve(effective_nthreads);
    for (int t = 0; t < effective_nthreads; ++t) {
      workers.emplace_back(workerFn);
    }
    for (auto& th : workers) {
      th.join();
    }
  }

  { std::lock_guard<std::mutex> lk(writeMutex); stopWriter = true; }
  writeCv.notify_all();
  writerThread.join();

  if (workerError) std::rethrow_exception(workerError);
}
