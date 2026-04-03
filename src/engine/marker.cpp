// marker.cpp — Thread-pool marker engine (pure C++17 / Eigen)
//
// Architecture:
//   - N worker threads pull chunks via atomic counter (work-stealing).
//   - 1 writer thread drains completed chunks in order (plain or gzip).
//   - Each worker owns a MethodBase clone + an independent PlinkCursor.
//   - Per-thread buffers (GVec, indexForMissing, result, fmtBuf) are
//     allocated once and reused across all markers to minimize heap traffic,
//     which matters because this pipeline is memory-bound.

#include "engine/marker.hpp"
#include "io/geno_data.hpp"
#include "util/logging.hpp"

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
#include <cstdint>
#include <zlib.h>


namespace {

// ──────────────────────────────────────────────────────────────────────
// TSV formatting helpers — zero allocation in the hot path
// ──────────────────────────────────────────────────────────────────────

// Build a suffix of nCols tab-separated "NA" values for fail-QC rows.
std::string makeNaSuffix(int nResultCols) {
  if (nResultCols <= 0) return {};
  std::string s;
  s.reserve(3 * static_cast<size_t>(nResultCols));
  for (int i = 0; i < nResultCols; ++i) { s += '\t'; s += 'N'; s += 'A'; }
  return s;
}

// Format double: "NA" | "Inf" | "-Inf" | "%.6g".  Returns char count.
int numToChars(char* buf, double x) {
  if (std::isnan(x)) { buf[0] = 'N'; buf[1] = 'A'; return 2; }
  if (std::isinf(x)) {
    if (x > 0) { buf[0] = 'I'; buf[1] = 'n'; buf[2] = 'f'; return 3; }
    else { buf[0] = '-'; buf[1] = 'I'; buf[2] = 'n'; buf[3] = 'f'; return 4; }
  }
  return std::snprintf(buf, 32, "%.6g", x);
}

// Append 9 meta columns: CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P
inline void appendMeta(
    std::string& out, char* buf,
    std::string_view chrom, uint32_t pos, std::string_view id,
    std::string_view ref, std::string_view alt,
    double missRate, double altFreq, double mac, double hweP)
{
  int n;
  out += chrom;  out += '\t';
  n = std::snprintf(buf, 32, "%u", pos);
  out.append(buf, n); out += '\t';
  out += id;     out += '\t';
  out += ref;    out += '\t';
  out += alt;    out += '\t';
  n = numToChars(buf, missRate); out.append(buf, n); out += '\t';
  n = numToChars(buf, altFreq);  out.append(buf, n); out += '\t';
  n = numToChars(buf, mac);      out.append(buf, n); out += '\t';
  n = numToChars(buf, hweP); out.append(buf, n);
}

// Full result line: meta + tab-separated doubles.
void formatLine(
    std::string& out, char* buf,
    std::string_view chrom, uint32_t pos, std::string_view id,
    std::string_view ref, std::string_view alt,
    double missRate, double altFreq, double mac, double hweP,
    const std::vector<double>& vals)
{
  appendMeta(out, buf, chrom, pos, id, ref, alt,
             missRate, altFreq, mac, hweP);
  for (double v : vals) {
    out += '\t';
    int n = numToChars(buf, v);
    out.append(buf, n);
  }
  out += '\n';
}

// Fail-QC line: meta + precomputed NA suffix.
void formatLineNA(
    std::string& out, char* buf,
    std::string_view chrom, uint32_t pos, std::string_view id,
    std::string_view ref, std::string_view alt,
    double missRate, double altFreq, double mac, double hweP,
    const std::string& naSuffix)
{
  appendMeta(out, buf, chrom, pos, id, ref, alt,
             missRate, altFreq, mac, hweP);
  out += naSuffix;
  out += '\n';
}


// ──────────────────────────────────────────────────────────────────────
// Per-worker thread context
// ──────────────────────────────────────────────────────────────────────

struct ThreadContext {
  std::unique_ptr<MethodBase> method;   // cloned per thread
  std::unique_ptr<GenoCursor> cursor;   // per-thread genotype decoder
  std::string naSuffix;

  ThreadContext(const MethodBase& proto, const GenoMeta& gd)
    : method(proto.clone()),
      cursor(gd.makeCursor()),
      naSuffix(makeNaSuffix(proto.resultSize()))
  {}
};


// ──────────────────────────────────────────────────────────────────────
// Output header
// ──────────────────────────────────────────────────────────────────────

constexpr const char* META_HEADER =
    "CHROM\tPOS\tID\tREF\tALT\tMISS_RATE\tALT_FREQ\tMAC\tHWE_P";

std::string buildHeader(const MethodBase& method) {
  return std::string(META_HEADER) + method.getHeaderColumns();
}


// ──────────────────────────────────────────────────────────────────────
// Padded flag: one per chunk, prevents false sharing between workers
// ──────────────────────────────────────────────────────────────────────

struct alignas(64) PaddedFlag { char ready; };
static_assert(sizeof(PaddedFlag) == 64, "PaddedFlag must be 64 bytes");

} // anonymous namespace


// ══════════════════════════════════════════════════════════════════════
// Engine
// ══════════════════════════════════════════════════════════════════════

void markerEngine(
    const GenoMeta& genoData,
    const MethodBase& method,
    const std::string& outputFile,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    bool exactHwe)
{
  const auto wallStart = std::chrono::steady_clock::now();
  const std::clock_t cpuStart = std::clock();

  const size_t nTotalChunks = genoData.chunkIndices().size();
  const int effective_nthreads =
      std::min(nthreads, static_cast<int>(nTotalChunks));

  infoMsg("Number of subjects in the input file: %u", genoData.nSubjInFile());
  infoMsg("Number of subjects to test: %u", genoData.nSubjUsed());
  infoMsg("Number of markers in the input file: %u", genoData.nMarkers());
  infoMsg("Number of markers to test: %zu", genoData.markerInfo().size());
  infoMsg("Number of chunks for all markers: %zu", nTotalChunks);
  infoMsg("Start chunk-level parallel processing with %d worker threads.",
          effective_nthreads);

  const std::string header = buildHeader(method);

  // Per-chunk output buffer + cache-line-padded ready flag.
  std::vector<std::string> chunkOutput(nTotalChunks);
  std::vector<PaddedFlag> chunkReady(nTotalChunks, {0});
  std::atomic<size_t> nextChunk(0);

  std::exception_ptr workerError = nullptr;
  std::mutex errorMutex;
  std::mutex writeMutex;
  std::condition_variable writeCv;
  std::atomic<bool> stopWriter{false};

  // ── Writer thread ──────────────────────────────────────────────────
  std::thread writerThread([&]() {
    const bool useGzip =
        outputFile.size() > 3 &&
        outputFile.compare(outputFile.size() - 3, 3, ".gz") == 0;
    gzFile gz = nullptr;
    std::ofstream plainOut;

    if (useGzip) {
      gz = gzopen(outputFile.c_str(), "wb");
      if (!gz) {
        std::lock_guard<std::mutex> lock(errorMutex);
        workerError = std::make_exception_ptr(
            std::runtime_error("Cannot open gzip output: " + outputFile));
        stopWriter.store(true);
        return;
      }
      // 256 KB write buffer — reduces syscall frequency.
      gzbuffer(gz, 256u * 1024u);
    } else {
      plainOut.open(outputFile.c_str());
      if (!plainOut.is_open()) {
        std::lock_guard<std::mutex> lock(errorMutex);
        workerError = std::make_exception_ptr(
            std::runtime_error("Cannot open output: " + outputFile));
        stopWriter.store(true);
        return;
      }
    }

    auto writeStr = [&](const std::string& s) {
      if (useGzip)
        gzwrite(gz, s.data(), static_cast<unsigned>(s.size()));
      else
        plainOut << s;
    };

    writeStr(header + "\n");

    for (size_t i = 0; i < chunkOutput.size(); ++i) {
      std::string tmp;
      {
        std::unique_lock<std::mutex> lk(writeMutex);
        writeCv.wait(lk, [&]() {
          return chunkReady[i].ready || stopWriter.load();
        });
        if (!chunkReady[i].ready) break;
        tmp = std::move(chunkOutput[i]);
      }
      writeStr(tmp);
      infoMsg("Writing finished: chunk %zu/%zu", i + 1, chunkOutput.size());
    }

    if (useGzip)
      gzclose(gz);
    else
      plainOut.close();
  });

  // ── Worker function ────────────────────────────────────────────────
  auto workerFn = [&]() {
    try {
      ThreadContext ctx(method, genoData);
      const GenoMeta& gd = genoData;
      GenoCursor& cursor = *ctx.cursor;
      MethodBase& meth = *ctx.method;
      const auto& localChunks = gd.chunkIndices();
      const std::string& naSuffix = ctx.naSuffix;

      // Per-thread reusable buffers — allocated once, reused for every marker.
      Eigen::VectorXd GVec(gd.nSubjUsed());
      std::vector<uint32_t> indexForMissing;
      indexForMissing.reserve(gd.nSubjUsed() / 10);
      std::vector<double> rv;
      rv.reserve(16);
      char fmtBuf[32];

      while (!stopWriter.load(std::memory_order_relaxed)) {
        const size_t cidx = nextChunk.fetch_add(1);
        if (cidx >= localChunks.size()) break;

        const auto& gIdx = localChunks[cidx];
        std::string out;
        out.reserve(gIdx.size() * 256);

        if (!gIdx.empty()) cursor.beginSequentialBlock(gIdx.front());
        meth.prepareChunk(gIdx);

        for (size_t i = 0; i < gIdx.size(); ++i) {
          double altFreq = NAN, altCounts = NAN;
          double missingRate = NAN, hweP = NAN;
          double maf = NAN, mac = NAN;
          indexForMissing.clear();

          cursor.getGenotypes(gIdx[i], GVec,
              altFreq, altCounts, missingRate, hweP,
              maf, mac, indexForMissing, exactHwe);

          const std::string_view chr    = gd.chr(gIdx[i]);
          const std::string_view ref    = gd.ref(gIdx[i]);
          const std::string_view alt    = gd.alt(gIdx[i]);
          const std::string_view marker = gd.markerId(gIdx[i]);
          const uint32_t pos            = gd.pos(gIdx[i]);

          const bool passQC =
              !(missingRate > missingCutoff ||
                maf < minMafCutoff ||
                mac < minMacCutoff);

          if (!passQC) {
            formatLineNA(out, fmtBuf, chr, pos, marker,
                         alt/*REF=bim6*/, ref/*ALT=bim5*/,
                         missingRate, altFreq, mac, hweP, naSuffix);
            continue;
          }

          // Impute missing genotypes with mean (2 × altFreq).
          const double imputeG = 2.0 * altFreq;
          double* gPtr = GVec.data();
          const uint32_t nMissing =
              static_cast<uint32_t>(indexForMissing.size());
          for (uint32_t j = 0; j < nMissing; ++j)
            gPtr[indexForMissing[j]] = imputeG;

          rv.clear();
          meth.getResultVec(GVec, altFreq, static_cast<int>(i), rv);

          formatLine(out, fmtBuf, chr, pos, marker,
                     alt/*REF=bim6*/, ref/*ALT=bim5*/,
                     missingRate, altFreq, mac, hweP, rv);
        }

        {
          std::lock_guard<std::mutex> lk(writeMutex);
          chunkOutput[cidx] = std::move(out);
          chunkReady[cidx].ready = 1;
        }
        infoMsg("Calculation finished: chunk %zu/%zu",
                cidx + 1, nTotalChunks);
        writeCv.notify_all();
      }
    } catch (...) {
      {
        std::lock_guard<std::mutex> lk(errorMutex);
        if (!workerError) workerError = std::current_exception();
      }
      {
        std::lock_guard<std::mutex> lk(writeMutex);
        stopWriter = true;
      }
      writeCv.notify_all();
    }
  };

  // ── Launch ─────────────────────────────────────────────────────────
  if (effective_nthreads <= 1) {
    workerFn();
  } else {
    std::vector<std::thread> workers;
    workers.reserve(effective_nthreads);
    for (int t = 0; t < effective_nthreads; ++t)
      workers.emplace_back(workerFn);
    for (auto& th : workers)
      th.join();
  }

  { std::lock_guard<std::mutex> lk(writeMutex); stopWriter = true; }
  writeCv.notify_all();
  writerThread.join();

  if (workerError) std::rethrow_exception(workerError);

  const double wallSec = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - wallStart).count();
  const double cpuSec =
      static_cast<double>(std::clock() - cpuStart) / CLOCKS_PER_SEC;
  infoMsg("Output written to: %s", outputFile.c_str());
  infoMsg("Wall time: %.1f seconds, CPU time: %.1f seconds",
          wallSec, cpuSec);
}
