// loco.cpp — LOCO data loading + generic LOCO engine implementation

#include "engine/loco.hpp"
#include "engine/marker_impl.hpp"
#include "geno_factory/geno_data.hpp"
#include "util/logging.hpp"
#include "util/text_stream.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <mutex>
#include "util/text_scanner.hpp"
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace engine_impl;

// ══════════════════════════════════════════════════════════════════════
// LocoData — load Regenie step 1 LOCO predictions
// ══════════════════════════════════════════════════════════════════════

// Parse a .loco file directly into aligned vectors.
// Memory-maps the file and parses with strtod for zero-copy I/O —
// avoids ifstream/getline overhead on multi-MB lines (400K columns).
// Only autosomal chromosomes 1–22 are kept.
//
// Header: FID_IID  0_HG00096  0_HG00097  ...
// Data:   1  0.0306  0.271  ...
static std::unordered_map<std::string, Eigen::VectorXd>parseLocoFile(
    const std::string &path,
    const std::vector<std::string> &usedIIDs
) {
    // ── Memory-map the file for zero-copy I/O ───────────────────────
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0)
        throw std::runtime_error("Cannot open LOCO file: " + path);

    struct stat st;
    if (::fstat(fd, &st) != 0) {
        ::close(fd);
        throw std::runtime_error("Cannot stat LOCO file: " + path);
    }
    const size_t fileSize = static_cast<size_t>(st.st_size);
    if (fileSize == 0) {
        ::close(fd);
        throw std::runtime_error("Empty LOCO file: " + path);
    }

    void *mapped = ::mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
    ::close(fd);
    if (mapped == MAP_FAILED)
        throw std::runtime_error("Cannot mmap LOCO file: " + path);
    ::madvise(mapped, fileSize, MADV_SEQUENTIAL);

    // RAII guard: munmap on any exit path
    struct Guard { void *p; size_t n; ~Guard()
                   {
                       ::munmap(p, n);
                   }

    } guard{mapped, fileSize};

    const char *p   = static_cast<const char *>(mapped);
    const char *const fend = p + fileSize;
    const auto nUsed = static_cast<Eigen::Index>(usedIIDs.size());

    // ── Build IID → position index ──────────────────────────────────
    std::unordered_map<std::string, Eigen::Index> iidIdx;
    iidIdx.reserve(usedIIDs.size());
    for (size_t i = 0; i < usedIIDs.size(); ++i)
        iidIdx[usedIIDs[i]] = static_cast<Eigen::Index>(i);

    // ── Parse header → build colPos (zero-copy from mmap) ───────────
    // Skip first token ("FID_IID")
    while (p < fend && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') ++p;

    std::vector<Eigen::Index> colPos;
    colPos.reserve(500000);
    std::vector<bool> usedHit(usedIIDs.size(), false);
    std::string tok;
    tok.reserve(64); // reused buffer for hash lookups

    while (p < fend && *p != '\n' && *p != '\r') {
        while (p < fend && (*p == ' ' || *p == '\t')) ++p;
        if (p >= fend || *p == '\n' || *p == '\r') break;
        const char *s = p;
        while (p < fend && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') ++p;
        const size_t len = static_cast<size_t>(p - s);

        tok.assign(s, len);
        Eigen::Index pos = -1;
        auto it = iidIdx.find(tok);
        if (it != iidIdx.end()) {
            pos = it->second;
        } else {
            // Extract IID after first '_' (Regenie "FID_IID" format)
            const char *u = s;
            while (u < s + len && *u != '_') ++u;
            if (u < s + len) {
                tok.assign(u + 1, s + len - u - 1);
                it = iidIdx.find(tok);
                if (it != iidIdx.end())
                    pos = it->second;
            }
        }
        if (pos >= 0) usedHit[static_cast<size_t>(pos)] = true;
        colPos.push_back(pos);
    }
    if (p < fend && *p == '\r') ++p;
    if (p < fend && *p == '\n') ++p;

    if (colPos.empty())
        throw std::runtime_error("LOCO file has no subject columns: " + path);

    const size_t nCols = colPos.size();

    // Verify every usedIID was matched
    for (size_t i = 0; i < usedIIDs.size(); ++i) {
        if (!usedHit[i])
            throw std::runtime_error("Subject '" + usedIIDs[i] +
                                     "' not found in LOCO file: " + path);
    }

    // ── Parse data rows with strtod directly from mmap ──────────────
    std::unordered_map<std::string, Eigen::VectorXd> result;

    while (p < fend) {
        while (p < fend && (*p == '\n' || *p == '\r')) ++p;
        if (p >= fend) break;

        // Chromosome label
        const char *cs = p;
        while (p < fend && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') ++p;
        const size_t clen = static_cast<size_t>(p - cs);
        if (clen == 0) break;

        // Autosomal check: "1"–"22" only (skip chr 23, X, etc.)
        bool autosomal = false;
        if (clen <= 2) {
            int chrNum = 0;
            for (size_t i = 0; i < clen; ++i) {
                if (cs[i] < '0' || cs[i] > '9') { chrNum = 0; break; }
                chrNum = chrNum * 10 + (cs[i] - '0');
            }
            autosomal = (chrNum >= 1 && chrNum <= 22);
        }
        if (!autosomal) {
            while (p < fend && *p != '\n') ++p;
            if (p < fend) ++p;
            continue;
        }

        std::string chrStr(cs, clen);
        Eigen::VectorXd vec = Eigen::VectorXd::Zero(nUsed);
        for (size_t c = 0; c < nCols; ++c) {
            while (p < fend && (*p == ' ' || *p == '\t')) ++p;
            char *ep;
            double val = std::strtod(p, &ep);
            if (ep == p)
                throw std::runtime_error("LOCO file " + path + " chr " + chrStr +
                                         ": invalid value at column " + std::to_string(c + 1));
            p = ep;
            if (colPos[c] >= 0)
                vec[colPos[c]] = val;
        }

        result[chrStr] = std::move(vec);
        while (p < fend && *p != '\n') ++p;
        if (p < fend) ++p;
    }

    return result;
}

LocoData LocoData::load(
    const std::string &predListFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &usedIIDs,
    const std::vector<std::string> &famIIDs
) {
    // Parse pred.list: phenoName  locoFilePath
    std::ifstream ifs(predListFile);
    if (!ifs.is_open())
        throw std::runtime_error("Cannot open pred.list file: " + predListFile);

    std::unordered_map<std::string, std::string> predMap; // phenoName → locoFilePath
    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        text::TokenScanner lts(line);
        std::string pheno = lts.next();
        std::string locoPath = lts.next();
        if (pheno.empty() || locoPath.empty())
            throw std::runtime_error("Invalid line in pred.list: " + line);
        predMap[pheno] = locoPath;
    }

    // Validate all requested phenotypes are in pred.list
    for (const auto &pn : phenoNames) {
        if (predMap.find(pn) == predMap.end())
            throw std::runtime_error("Phenotype '" + pn +
                                     "' not found in pred.list file: " + predListFile);
    }

    LocoData result;

    for (const auto &pn : phenoNames) {
        const std::string &locoPath = predMap.at(pn);

        result.scores[pn] = parseLocoFile(locoPath, usedIIDs);

        if (result.scores[pn].size() < 22)
            throw std::runtime_error("LOCO file for '" + pn + "' contains " +
                                     std::to_string(result.scores[pn].size()) +
                                     " autosomal chromosomes (need at least 22): " + locoPath);
    }

    return result;
}

std::unordered_set<std::string> LocoData::availableChromosomes() const {
    if (scores.empty()) return {};

    // Start with chromosomes from first phenotype
    auto it = scores.begin();
    std::unordered_set<std::string> common;
    for (const auto &[chr, _] : it->second)
        common.insert(chr);

    // Intersect with remaining phenotypes
    for (++it; it != scores.end(); ++it) {
        std::unordered_set<std::string> next;
        for (const auto &chr : common) {
            if (it->second.count(chr))
                next.insert(chr);
        }
        common = std::move(next);
    }

    return common;
}

// ══════════════════════════════════════════════════════════════════════
// locoEngine — generic LOCO engine
// ══════════════════════════════════════════════════════════════════════

void locoEngine(
    const GenoMeta &genoData,
    const std::unordered_set<std::string> &locoChroms,
    const std::vector<std::string> &phenoNames,
    LocoTaskBuilder buildTasks,
    const std::string &outPrefix,
    const std::string &methodName,
    const std::string &compression,
    int compressionLevel,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff
) {
    const size_t K = phenoNames.size();
    const auto &allChunks = genoData.chunkIndices();
    const auto &markers = genoData.markerInfo();

    // ── 1. Extract unique chromosomes in encounter order ────────────
    std::vector<std::string> chrOrder;
    {
        std::unordered_set<std::string> seen;
        for (const auto &m : markers) {
            if (seen.insert(m.chrom).second)
                chrOrder.push_back(m.chrom);
        }
    }

    // ── 2. Intersect with LOCO data (O4) ────────────────────────────
    std::vector<std::string> activeChroms;
    for (const auto &chr : chrOrder) {
        if (locoChroms.count(chr)) {
            activeChroms.push_back(chr);
        } else {
            infoMsg("LOCO: skipping chromosome %s (not in pred.list)", chr.c_str());
        }
    }

    if (activeChroms.empty())
        throw std::runtime_error("LOCO: no chromosomes overlap between genotype data and pred.list");

    // ── 3. Partition chunks by chromosome ───────────────────────────
    // Chunks are already chromosome-aligned from buildChunks.
    // Build a map: chrStr → [startChunkIdx, endChunkIdx)
    struct ChrRange {
        size_t start, end;
    };

    std::unordered_map<std::string, ChrRange> chrChunks;
    {
        size_t ci = 0;
        while (ci < allChunks.size()) {
            const auto &firstIdx = allChunks[ci];
            if (firstIdx.empty()) { ++ci; continue; }
            std::string chr(genoData.chr(firstIdx[0]));
            size_t start = ci;
            // Find all consecutive chunks on the same chromosome
            while (ci < allChunks.size()) {
                const auto &chunk = allChunks[ci];
                if (chunk.empty()) { ++ci; continue; }
                if (std::string(genoData.chr(chunk[0])) != chr) break;
                ++ci;
            }
            chrChunks[chr] = {start, ci};
        }
    }

    // ── 4. Build tasks for first chromosome to get headers ──────────
    std::vector<PhenoTask> tasks;
    buildTasks(activeChroms[0], tasks);

    // Build per-phenotype headers and NA suffixes
    std::vector<std::string> headers(K);
    for (size_t p = 0; p < K; ++p)
        headers[p] = std::string(META_HEADER) + tasks[p].method->getHeaderColumns() + "\n";

    // Build output file paths
    auto writerMode = TextWriter::modeFromString(compression);
    std::vector<std::string> outPaths(K);
    for (size_t p = 0; p < K; ++p)
        outPaths[p] = TextWriter::buildOutputPath(outPrefix, phenoNames[p], methodName, compression);

    // ── 5. Open persistent writers and write headers ────────────────
    std::vector<TextWriter> writers;
    writers.reserve(K);
    for (size_t p = 0; p < K; ++p) {
        writers.emplace_back(outPaths[p], writerMode, compressionLevel);
        writers[p].write(headers[p]);
    }

    infoMsg("Number of markers to test: %zu", markers.size());
    infoMsg("Number of chunks for all markers: %zu", allChunks.size());
    infoMsg("LOCO: %zu chromosomes to process", activeChroms.size());

    // ── 6. Per-chromosome loop ──────────────────────────────────────
    for (size_t chrIdx = 0; chrIdx < activeChroms.size(); ++chrIdx) {
        const std::string &chr = activeChroms[chrIdx];
        auto rangeIt = chrChunks.find(chr);
        if (rangeIt == chrChunks.end()) {
            infoMsg("LOCO: chromosome %s has no chunks, skipping", chr.c_str());
            continue;
        }

        const size_t chunkStart = rangeIt->second.start;
        const size_t chunkEnd = rangeIt->second.end;
        const size_t nChrChunks = chunkEnd - chunkStart;

        // Build tasks for this chromosome (skip first chr — already built above)
        if (chrIdx > 0)
            buildTasks(chr, tasks);

        // Build per-phenotype NA suffixes for this chromosome's methods
        std::vector<std::string> naSuffixes(K);
        for (size_t p = 0; p < K; ++p)
            naSuffixes[p] = makeNaSuffix(tasks[p].method->resultSize());

        infoMsg("LOCO: chromosome %s — %zu chunks, %d workers (%zu/%zu)",
                chr.c_str(), nChrChunks,
                std::min(nthreads, static_cast<int>(nChrChunks)),
                chrIdx + 1, activeChroms.size());

        multiPhenoEngineRange(genoData, tasks, naSuffixes, chunkStart, chunkEnd, writers,
                              nthreads, missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);

        infoMsg("LOCO: chromosome %s finished (%zu/%zu)",
                chr.c_str(), chrIdx + 1, activeChroms.size());
    }

    // ── 7. Close writers ────────────────────────────────────────────
    for (auto &w : writers)
        w.close();

    for (size_t p = 0; p < K; ++p)
        infoMsg("Output written to: %s", outPaths[p].c_str());
}
