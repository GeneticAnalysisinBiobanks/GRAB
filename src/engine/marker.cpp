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
#include "util/text_stream.hpp"

#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <algorithm>
#include <cmath>
#include <exception>
#include <cstdio>
#include <cstdint>


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
    double missRate, double altFreq, double mac, double hweP
) {
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
    const std::vector<double>& vals
) {
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
    const std::string& naSuffix
) {
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
    bool exactHwe
) {
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
    try {
      TextWriter writer(outputFile);

      writer.write(header + "\n");

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
        writer.write(tmp);
        infoMsg("Writing finished: chunk %zu/%zu", i + 1, chunkOutput.size());
      }

      writer.close();
    } catch (...) {
      std::lock_guard<std::mutex> lock(errorMutex);
      workerError = std::current_exception();
      stopWriter.store(true);
    }
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

  infoMsg("Output written to: %s", outputFile.c_str());
}


// ──────────────────────────────────────────────────────────────────────
// MultiMethod implementation
// ──────────────────────────────────────────────────────────────────────

MultiMethod::MultiMethod(
    std::vector<std::unique_ptr<MethodBase>> methods,
    std::vector<std::string> residNames,
    std::vector<std::string> suffixes)
  : m_methods(std::move(methods)),
    m_residNames(std::move(residNames)),
    m_suffixes(std::move(suffixes))
{}

std::unique_ptr<MethodBase> MultiMethod::clone() const {
  std::vector<std::unique_ptr<MethodBase>> cloned;
  cloned.reserve(m_methods.size());
  for (const auto& m : m_methods) cloned.push_back(m->clone());
  return std::make_unique<MultiMethod>(
      std::move(cloned), m_residNames, m_suffixes);
}

int MultiMethod::resultSize() const {
  return static_cast<int>(m_methods.size() * m_suffixes.size());
}

std::string MultiMethod::getHeaderColumns() const {
  std::string h;
  for (const auto& name : m_residNames)
    for (const auto& suf : m_suffixes)
      h += "\t" + name + suf;
  return h;
}

void MultiMethod::prepareChunk(const std::vector<uint64_t>& gIndices) {
  for (auto& m : m_methods) m->prepareChunk(gIndices);
}

void MultiMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq, int markerInChunkIdx,
    std::vector<double>& result)
{
  std::vector<double> inner;
  for (auto& m : m_methods) {
    inner.clear();
    m->getResultVec(GVec, altFreq, markerInChunkIdx, inner);
    result.insert(result.end(), inner.begin(), inner.end());
  }
}


// ══════════════════════════════════════════════════════════════════════
// multiPhenoEngine — per-phenotype independent GWAS
// ══════════════════════════════════════════════════════════════════════

namespace {

// HWE exact test (SNPHWE2) — duplicated from plink.cpp to avoid linkage issues.
static double HweExact(uint32_t obs_hets, uint32_t obs_hom1, uint32_t obs_hom2) {
  const int64_t obs_homc = std::max(obs_hom1, obs_hom2);
  const int64_t obs_homr = std::min(obs_hom1, obs_hom2);
  const int64_t rare = 2 * obs_homr + static_cast<int64_t>(obs_hets);
  const int64_t n    = static_cast<int64_t>(obs_hets) + obs_homc + obs_homr;
  const int64_t obs  = static_cast<int64_t>(obs_hets);
  if (n == 0) return 1.0;
  int64_t mid = (rare * (2 * n - rare)) / (2 * n);
  if ((rare & 1) ^ (mid & 1)) ++mid;
  {
    int64_t hr = (rare - mid) / 2;
    int64_t hc = n - mid - hr;
    if (mid + 2 <= rare && hr > 0 &&
        4.0 * hr * hc > (mid + 2.0) * (mid + 1.0))
      mid += 2;
    else if (mid >= 2 &&
             static_cast<double>(mid) * (mid - 1) > 4.0 * (hr + 1.0) * (hc + 1.0))
      mid -= 2;
  }
  const int64_t mid_homr = (rare - mid) / 2;
  const int64_t mid_homc = n - mid - mid_homr;
  double sum = 1.0, p = 0.0, thresh;
  if (obs <= mid) {
    { double prob = 1.0; int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h > obs; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob; ++cr; ++cc;
      } thresh = prob; p = thresh;
      for (int64_t h = obs; h >= 2; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob; p += prob; ++cr; ++cc;
      }
    }
    { double prob = 1.0; int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h <= rare - 2; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob; if (prob <= thresh) p += prob; --cr; --cc;
      }
    }
  } else {
    { double prob = 1.0; int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h < obs; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob; --cr; --cc;
      } thresh = prob; p = thresh;
      for (int64_t h = obs; h <= rare - 2; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob; p += prob; --cr; --cc;
      }
    }
    { double prob = 1.0; int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h >= 2; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob; if (prob <= thresh) p += prob; ++cr; ++cc;
      }
    }
  }
  return std::min(p / sum, 1.0);
}

static double HweChiSq(uint32_t nHet, uint32_t nHom1, uint32_t nHom2) {
  const double n = static_cast<double>(nHet + nHom1 + nHom2);
  if (n == 0.0) return std::numeric_limits<double>::quiet_NaN();
  const double f = (2.0 * nHom1 + nHet) / (2.0 * n);
  const double g = 1.0 - f;
  const double Eh = 2.0 * f * g * n, E1 = f * f * n, E2 = g * g * n;
  if (E1 < 5.0 || Eh < 5.0 || E2 < 5.0)
    return std::numeric_limits<double>::quiet_NaN();
  const double d1 = nHom1 - E1, dh = nHet - Eh, d2 = nHom2 - E2;
  const double chi2 = d1*d1/E1 + dh*dh/Eh + d2*d2/E2;
  return std::erfc(std::sqrt(chi2 * 0.5));
}

// Compute allele stats from a per-phenotype genotype vector (doubles).
struct PhenoGenoStats {
  double altFreq, mac, missingRate, hweP;
};

static PhenoGenoStats statsFromGVec(
    const double* g, uint32_t n,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe)
{
  indexForMissing.clear();
  uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0;
  for (uint32_t i = 0; i < n; ++i) {
    const double v = g[i];
    if (std::isnan(v) || v < 0.0) { indexForMissing.push_back(i); continue; }
    if (v == 0.0)      ++nHomRef;
    else if (v == 1.0) ++nHet;
    else if (v == 2.0) ++nHomAlt;
    else               indexForMissing.push_back(i);
  }
  PhenoGenoStats s;
  uint32_t nonMissing = nHomRef + nHet + nHomAlt;
  if (nonMissing == 0) {
    s.altFreq = NAN; s.mac = NAN; s.missingRate = 1.0; s.hweP = NAN;
    return s;
  }
  uint32_t altCounts = 2 * nHomAlt + nHet;
  double n2 = 2.0 * nonMissing;
  s.altFreq = altCounts / n2;
  double maf = std::min(s.altFreq, 1.0 - s.altFreq);
  s.mac = maf * n2;
  s.missingRate = static_cast<double>(indexForMissing.size()) / n;
  if (exactHwe)
    s.hweP = HweExact(nHet, nHomAlt, nHomRef);
  else
    s.hweP = (nHet < 5 || nHomAlt < 5 || nHomRef < 5)
             ? HweExact(nHet, nHomAlt, nHomRef)
             : HweChiSq(nHet, nHomAlt, nHomRef);
  return s;
}

// Extract per-phenotype genotype vector from union vector.
static void extractPhenoGVec(
    const double* unionG, uint32_t nUnion,
    const uint32_t* unionToLocal, uint32_t nPheno,
    double* phenoG)
{
  for (uint32_t i = 0; i < nUnion; ++i) {
    uint32_t li = unionToLocal[i];
    if (li != UINT32_MAX)
      phenoG[li] = unionG[i];
  }
}

} // anonymous namespace


void multiPhenoEngine(
    const GenoMeta& genoData,
    std::vector<PhenoTask>& tasks,
    const std::string& outPrefix,
    const std::string& methodName,
    const std::string& compression,
    int compressionLevel,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    bool exactHwe
) {
  const size_t K = tasks.size();
  const size_t nTotalChunks = genoData.chunkIndices().size();
  const int effective_nthreads =
      std::min(nthreads, static_cast<int>(nTotalChunks));
  const uint32_t nUnion = genoData.nSubjUsed();

  infoMsg("Number of subjects in the input file: %u", genoData.nSubjInFile());
  infoMsg("Number of subjects in union mask: %u", nUnion);
  for (size_t p = 0; p < K; ++p)
    infoMsg("  Phenotype '%s': %u subjects",
            tasks[p].phenoName.c_str(), tasks[p].nUsed);
  infoMsg("Number of markers in the input file: %u", genoData.nMarkers());
  infoMsg("Number of markers to test: %zu", genoData.markerInfo().size());
  infoMsg("Number of chunks for all markers: %zu", nTotalChunks);
  infoMsg("Start chunk-level parallel processing with %d worker threads, %zu phenotypes.",
          effective_nthreads, K);

  // Build per-phenotype headers.
  std::vector<std::string> headers(K);
  for (size_t p = 0; p < K; ++p)
    headers[p] = std::string(META_HEADER) + tasks[p].method->getHeaderColumns() + "\n";

  // Build per-phenotype NA suffixes.
  std::vector<std::string> naSuffixes(K);
  for (size_t p = 0; p < K; ++p)
    naSuffixes[p] = makeNaSuffix(tasks[p].method->resultSize());

  // Build output file paths.
  auto writerMode = TextWriter::modeFromString(compression);
  std::vector<std::string> outPaths(K);
  for (size_t p = 0; p < K; ++p)
    outPaths[p] = TextWriter::buildOutputPath(outPrefix, tasks[p].phenoName,
                                              methodName, compression);

  // Max per-phenotype sample count (for buffer sizing).
  uint32_t maxPhenoN = 0;
  for (size_t p = 0; p < K; ++p)
    maxPhenoN = std::max(maxPhenoN, tasks[p].nUsed);

  // Per-chunk, per-phenotype output buffers + ready flags.
  // chunkOutput[chunk][pheno]
  std::vector<std::vector<std::string>> chunkOutput(nTotalChunks,
      std::vector<std::string>(K));
  std::vector<PaddedFlag> chunkReady(nTotalChunks, {0});
  std::atomic<size_t> nextChunk(0);

  std::exception_ptr workerError = nullptr;
  std::mutex errorMutex;
  std::mutex writeMutex;
  std::condition_variable writeCv;
  std::atomic<bool> stopWriter{false};

  // ── Writer thread ──────────────────────────────────────────────────
  std::thread writerThread([&]() {
    try {
      std::vector<TextWriter> writers;
      writers.reserve(K);
      for (size_t p = 0; p < K; ++p) {
        writers.emplace_back(outPaths[p], writerMode, compressionLevel);
        writers[p].write(headers[p]);
      }

      for (size_t i = 0; i < nTotalChunks; ++i) {
        std::vector<std::string> tmp(K);
        {
          std::unique_lock<std::mutex> lk(writeMutex);
          writeCv.wait(lk, [&]() {
            return chunkReady[i].ready || stopWriter.load();
          });
          if (!chunkReady[i].ready) break;
          for (size_t p = 0; p < K; ++p)
            tmp[p] = std::move(chunkOutput[i][p]);
        }
        for (size_t p = 0; p < K; ++p)
          writers[p].write(tmp[p]);
        infoMsg("Writing finished: chunk %zu/%zu", i + 1, nTotalChunks);
      }

      for (auto& w : writers) w.close();
    } catch (...) {
      std::lock_guard<std::mutex> lock(errorMutex);
      workerError = std::current_exception();
      stopWriter.store(true);
    }
  });

  // ── Worker function ────────────────────────────────────────────────
  auto workerFn = [&]() {
    try {
      // Per-thread state
      auto cursor = genoData.makeCursor();
      std::vector<std::unique_ptr<MethodBase>> methods(K);
      for (size_t p = 0; p < K; ++p)
        methods[p] = tasks[p].method->clone();

      Eigen::VectorXd GVec_union(nUnion);
      Eigen::VectorXd GVec_pheno(maxPhenoN);
      std::vector<uint32_t> unionMissing;
      unionMissing.reserve(nUnion / 10);
      std::vector<uint32_t> phenoMissing;
      phenoMissing.reserve(maxPhenoN / 10);
      std::vector<double> rv;
      rv.reserve(16);
      char fmtBuf[32];

      const auto& localChunks = genoData.chunkIndices();

      while (!stopWriter.load(std::memory_order_relaxed)) {
        const size_t cidx = nextChunk.fetch_add(1);
        if (cidx >= localChunks.size()) break;

        const auto& gIdx = localChunks[cidx];

        // Per-phenotype output buffers for this chunk.
        std::vector<std::string> phenoOut(K);
        for (size_t p = 0; p < K; ++p)
          phenoOut[p].reserve(gIdx.size() * 128);

        if (!gIdx.empty()) cursor->beginSequentialBlock(gIdx.front());
        for (size_t p = 0; p < K; ++p)
          methods[p]->prepareChunk(gIdx);

        for (size_t i = 0; i < gIdx.size(); ++i) {
          // 1. Decode genotypes for union mask
          double uAltFreq, uAltCounts, uMissRate, uHweP, uMaf, uMac;
          unionMissing.clear();
          cursor->getGenotypes(gIdx[i], GVec_union,
              uAltFreq, uAltCounts, uMissRate, uHweP,
              uMaf, uMac, unionMissing, exactHwe);

          const std::string_view chr    = genoData.chr(gIdx[i]);
          const std::string_view ref    = genoData.ref(gIdx[i]);
          const std::string_view alt    = genoData.alt(gIdx[i]);
          const std::string_view marker = genoData.markerId(gIdx[i]);
          const uint32_t pos            = genoData.pos(gIdx[i]);

          // 2. For each phenotype: extract, compute stats, QC, impute, run
          for (size_t p = 0; p < K; ++p) {
            const auto& task = tasks[p];
            const uint32_t nP = task.nUsed;

            // Extract per-phenotype genotypes
            extractPhenoGVec(GVec_union.data(), nUnion,
                             task.unionToLocal.data(), nP,
                             GVec_pheno.data());

            // Compute per-phenotype stats
            PhenoGenoStats gs = statsFromGVec(
                GVec_pheno.data(), nP, phenoMissing, exactHwe);

            double pMaf = std::min(gs.altFreq, 1.0 - gs.altFreq);
            const bool passQC =
                !(gs.missingRate > missingCutoff ||
                  pMaf < minMafCutoff ||
                  gs.mac < minMacCutoff);

            if (!passQC) {
              formatLineNA(phenoOut[p], fmtBuf, chr, pos, marker,
                           alt, ref,
                           gs.missingRate, gs.altFreq, gs.mac, gs.hweP,
                           naSuffixes[p]);
              continue;
            }

            // Impute missing with mean
            const double imputeG = 2.0 * gs.altFreq;
            double* gPtr = GVec_pheno.data();
            for (uint32_t j : phenoMissing)
              gPtr[j] = imputeG;

            rv.clear();
            methods[p]->getResultVec(
                GVec_pheno.head(nP), gs.altFreq,
                static_cast<int>(i), rv);

            formatLine(phenoOut[p], fmtBuf, chr, pos, marker,
                       alt, ref,
                       gs.missingRate, gs.altFreq, gs.mac, gs.hweP, rv);
          }
        }

        {
          std::lock_guard<std::mutex> lk(writeMutex);
          for (size_t p = 0; p < K; ++p)
            chunkOutput[cidx][p] = std::move(phenoOut[p]);
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

  for (size_t p = 0; p < K; ++p)
    infoMsg("Output written to: %s", outPaths[p].c_str());
}
