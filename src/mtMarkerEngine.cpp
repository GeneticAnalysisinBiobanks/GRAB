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

mtPOLMMClass*           ptr_gPOLMMobj     = nullptr;
mtWtCoxGClass*         ptr_gWtCoxGobj     = nullptr;
mtLEAFClass*             ptr_gLEAFobj     = nullptr;
mtSPAGRMClass*         ptr_gSPAGRMobj     = nullptr;
mtSAGELDClass*         ptr_gSAGELDobj     = nullptr;
mtSPAsqrClass*         ptr_gSPAsqrobj     = nullptr;
mtSPACoxClass*         ptr_gSPACoxobj     = nullptr;
mtSPAmixClass*         ptr_gSPAmixobj     = nullptr;
mtSPAmixPlusClass*     ptr_gSPAmixPlusobj = nullptr;
PlinkData*              ptr_gPlinkDataObj = nullptr;
PlinkCursor*            ptr_gPlinkCursorObj = nullptr;


namespace {

// Thread-safe timestamped info message.
// Uses fprintf(stderr) — safe to call from any thread.
// Rprintf/REprintf are NOT thread-safe and must only be called from R's main thread.
static std::mutex g_logMutex;
void infoMsg(const char* fmt, ...) {
  char buf[512];
  va_list args;
  va_start(args, fmt);
  std::vsnprintf(buf, sizeof(buf), fmt, args);
  va_end(args);
  // localtime uses a static buffer — must be serialized.
  // Keep critical section minimal: snapshot time, format, and print under lock.
  std::lock_guard<std::mutex> lk(g_logMutex);
  auto t = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  char ts[20];
  std::strftime(ts, sizeof(ts), "%Y-%m-%d %H:%M:%S", std::localtime(&t));
  fprintf(stderr, "[INFO] %s %s\n", ts, buf);
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

// Append meta columns directly to out (no temp string allocation).
inline void appendMeta(
  std::string& out,
  char* buf,
  std::string_view marker,
  std::string_view chr,
  uint32_t pos,
  std::string_view ref,
  std::string_view alt,
  double altCounts,
  double altFreq,
  double missingRate,
  double hweP
) {
  int n;
  out += marker; out += '\t';
  out += chr;    out += '\t';
  n = std::snprintf(buf, 32, "%u", pos);
  out.append(buf, n); out += '\t';
  out += ref;    out += '\t';
  out += alt;    out += '\t';
  n = numToChars(buf, altCounts); out.append(buf, n); out += '\t';
  n = numToChars(buf, altFreq);   out.append(buf, n); out += '\t';
  n = numToChars(buf, missingRate); out.append(buf, n); out += '\t';
  n = numToChars(buf, hweP); out.append(buf, n);
}

// Format one output line: meta columns + numeric result values.
// Writes directly into out — no intermediate string allocation.
void formatLine(
  std::string& out,
  char* buf,
  std::string_view marker,
  std::string_view chr,
  uint32_t pos,
  std::string_view ref,
  std::string_view alt,
  double altCounts,
  double altFreq,
  double missingRate,
  double hweP,
  const std::vector<double>& vals
) {
  appendMeta(out, buf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP);
  int n;
  for (double v : vals) {
    out += '\t';
    n = numToChars(buf, v);
    out.append(buf, n);
  }
  out += '\n';
}

// Format one output line for fail-QC markers using precomputed NA suffix.
// Writes directly into out — no intermediate string allocation.
void formatLineNA(
  std::string& out,
  char* buf,
  std::string_view marker,
  std::string_view chr,
  uint32_t pos,
  std::string_view ref,
  std::string_view alt,
  double altCounts,
  double altFreq,
  double missingRate,
  double hweP,
  const std::string& naSuffix
) {
  appendMeta(out, buf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP);
  out += naSuffix;
  out += '\n';
}

// ---- File-local structs ----

// Per-worker independent copies of all shared state.
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
  std::unique_ptr<PlinkData>             plinkData;
  std::unique_ptr<PlinkCursor>           cursor;
  std::string                            naSuffix;
};


// ---- Output header ----

const char* META_HEADER = "Marker\tChr\tPos\tRef\tAlt\tAltCount\tAltFreq\tMissRate\tHWEpval";

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

  // Copy PlinkData and PlinkCursor for this thread
  if (!ptr_gPlinkDataObj) throw std::runtime_error("PlinkData is not initialized.");
  if (!ptr_gPlinkCursorObj) throw std::runtime_error("PlinkCursor is not initialized.");
  ctx.plinkData.reset(new PlinkData(*ptr_gPlinkDataObj));
  ctx.cursor.reset(new PlinkCursor(*ptr_gPlinkCursorObj));
  ctx.naSuffix = makeNaSuffix(method);
  return ctx;
}

} // namespace


// ---- The Engine ----

void mtMarkerEngine(
  const std::string method,
  const std::string outputFile,
  const int nthreads,
  const std::string impute_method,
  const double missing_cutoff,
  const double min_maf_marker,
  const double min_mac_marker
) {
  if (!ptr_gPlinkDataObj) throw std::runtime_error("PlinkData is not initialized.");
  const PlinkData& plinkRef = *ptr_gPlinkDataObj;

  const size_t nTotalChunks = plinkRef.chunkIndices().size();
  const int effective_nthreads = std::min(nthreads, static_cast<int>(nTotalChunks));

  infoMsg("Number of subjects in the input file: %u", plinkRef.nSubjInFile());
  infoMsg("Number of subjects to test: %u", plinkRef.nSubjUsed());
  infoMsg("Number of markers in the input file: %u", plinkRef.nMarkers());
  infoMsg("Number of markers to test: %zu", plinkRef.markerInfo().size());
  infoMsg("Number of chunks for all markers: %zu", nTotalChunks);
  infoMsg("Number of threads: %d", effective_nthreads);

  const std::string header = getHeader(method);
  std::vector<std::string> chunkOutput(nTotalChunks);
  std::vector<char> chunkReady(nTotalChunks, 0);
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
      // Move the chunk string out under the lock, then do I/O outside the lock.
      // This unblocks workers from storing subsequent completed chunks while we write.
      std::string tmp;
      {
        std::unique_lock<std::mutex> lk(writeMutex);
        writeCv.wait(lk, [&]() { return chunkReady[i] || stopWriter.load(); });
        if (!chunkReady[i]) break;  // writer-stop or error: abandon remaining chunks
        tmp = std::move(chunkOutput[i]);  // frees chunkOutput[i] memory immediately
      }
      writeStr(tmp);
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
      // Per-thread: copy all shared state (PlinkData + method + PlinkCursor)
      ThreadContext ctx = makeThreadContext(method);
      const PlinkData& pd = *ctx.plinkData;
      PlinkCursor& cursor = *ctx.cursor;
      const auto& localChunks = pd.chunkIndices();
      const std::string& naSuffix = ctx.naSuffix;

      while (!stopWriter.load(std::memory_order_relaxed)) {
        size_t cidx = nextChunk.fetch_add(1);
        if (cidx >= localChunks.size()) break;

        const auto& gIdx = localChunks[cidx];
        std::string out;
        out.reserve(gIdx.size() * 256);

        if (!gIdx.empty()) cursor.beginSequentialBlock(gIdx.front());

        if (method == "WtCoxG") ctx.wtcoxg->prepareChunk(gIdx);
        if (method == "LEAF")   ctx.leaf->prepareChunk(gIdx);

        // Per-thread reusable buffers (allocated once, reused across all markers)
        std::vector<uint32_t> indexForMissing;
        indexForMissing.reserve(pd.nSubjUsed() / 10);
        arma::vec GVec(pd.nSubjUsed());         // reused by getGenotypes
        std::vector<double> rv;                  // reused by getResultVec
        rv.reserve(16);
        char fmtBuf[32];  // scratch for numToChars / snprintf

        for (size_t i = 0; i < gIdx.size(); ++i) {
          double altFreq = NAN, altCounts = NAN, missingRate = NAN, hweP = NAN;
          double maf = NAN, mac = NAN;
          indexForMissing.clear();

          cursor.getGenotypes(gIdx[i], GVec,
            altFreq, altCounts, missingRate, hweP, maf, mac, indexForMissing);

          const std::string_view chr    = pd.chr(gIdx[i]);
          const std::string_view ref    = pd.ref(gIdx[i]);
          const std::string_view alt    = pd.alt(gIdx[i]);
          const std::string_view marker = pd.markerId(gIdx[i]);
          uint32_t pos                  = pd.pos(gIdx[i]);

          bool passQC = !((missingRate > missing_cutoff) || (maf < min_maf_marker) || (mac < min_mac_marker));

          bool flip = false;
          if (passQC) {
            const uint32_t nMissing = static_cast<uint32_t>(indexForMissing.size());
            const double imputeG = 2.0 * altFreq;
            double* gPtr = GVec.memptr();
            for (uint32_t j = 0; j < nMissing; ++j) {
              gPtr[indexForMissing[j]] = imputeG;
            }
            if (method != "WtCoxG" && altFreq > 0.5) {
              // In-place flip: avoids temporary arma::vec allocation
              const uint32_t n = GVec.n_elem;
              for (uint32_t k = 0; k < n; ++k) gPtr[k] = 2.0 - gPtr[k];
              flip = true;
            }
          }

          // --- SPAmix ---
          if (method == "SPAmix") {
            if (!passQC) { formatLineNA(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, naSuffix); continue; }
            ctx.spamix->getResultVec(GVec, altFreq, rv);
            formatLine(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, rv);
            continue;
          }

          // --- SPAmixPlus ---
          if (method == "SPAmixPlus") {
            if (!passQC) { formatLineNA(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, naSuffix); continue; }
            ctx.spamixPlus->getResultVec(GVec, altFreq, rv);
            formatLine(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, rv);
            continue;
          }

          // --- POLMM ---
          if (method == "POLMM") {
            if (!passQC) { formatLineNA(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, naSuffix); continue; }
            ctx.polmm->getResultVec(GVec, altFreq, rv);
            rv[1] *= (1.0 - 2.0 * static_cast<double>(flip)); // flip beta
            formatLine(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, rv);
            continue;
          }

          // --- SPACox ---
          if (method == "SPACox") {
            if (!passQC) { formatLineNA(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, naSuffix); continue; }
            ctx.spacox->getResultVec(GVec, altFreq, rv);
            formatLine(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, rv);
            continue;
          }

          // --- SPAGRM ---
          if (method == "SPAGRM") {
            if (!passQC) { formatLineNA(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, naSuffix); continue; }
            ctx.spagrm->getResultVec(GVec, altFreq, rv);
            formatLine(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, rv);
            continue;
          }

          // --- SAGELD ---
          if (method == "SAGELD") {
            if (!passQC) { formatLineNA(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, naSuffix); continue; }
            ctx.sageld->getResultVec(GVec, altFreq, rv);
            formatLine(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, rv);
            continue;
          }

          // --- WtCoxG ---
          if (method == "WtCoxG") {
            if (!passQC) { formatLineNA(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, naSuffix); continue; }
            ctx.wtcoxg->getResultVec(GVec, static_cast<int>(i), rv);
            formatLine(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, rv);
            continue;
          }

          // --- SPAsqr ---
          if (method == "SPAsqr") {
            if (!passQC) { formatLineNA(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, naSuffix); continue; }
            ctx.spasqr->getResultVec(GVec, altFreq, rv);
            formatLine(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, rv);
            continue;
          }

          // --- LEAF ---
          if (method == "LEAF") {
            if (!passQC) { formatLineNA(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, naSuffix); continue; }
            ctx.leaf->getResultVec(GVec, static_cast<int>(i), rv);
            formatLine(out, fmtBuf, marker, chr, pos, ref, alt, altCounts, altFreq, missingRate, hweP, rv);
            continue;
          }
        }

        {
          std::lock_guard<std::mutex> lk(writeMutex);
          chunkOutput[cidx] = std::move(out);
          chunkReady[cidx] = 1;
        }
        infoMsg("Calculation finished: chunk %zu/%zu", cidx + 1, nTotalChunks);
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
