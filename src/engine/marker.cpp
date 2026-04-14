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
#include "engine/marker_impl.hpp"
#include "geno_factory/geno_data.hpp"
#include "util/logging.hpp"
#include "util/text_stream.hpp"

#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <mutex>
#include <thread>

using namespace engine_impl;

// ══════════════════════════════════════════════════════════════════════
// Engine
// ══════════════════════════════════════════════════════════════════════

void markerEngine(
    const GenoMeta &genoData,
    const MethodBase &method,
    const std::string &outputFile,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff
) {
    const size_t nTotalChunks = genoData.chunkIndices().size();
    const int effective_nthreads = std::min(nthreads, static_cast<int>(nTotalChunks));

    infoMsg("Number of subjects in the input file: %u", genoData.nSubjInFile());
    infoMsg("Number of subjects to test: %u", genoData.nSubjUsed());
    infoMsg("Number of markers in the input file: %u", genoData.nMarkers());
    infoMsg("Number of markers to test: %zu", genoData.markerInfo().size());
    infoMsg("Number of chunks for all markers: %zu", nTotalChunks);
    infoMsg("Start chunk-level parallel processing with %d worker threads.", effective_nthreads);

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
            const GenoMeta &gd = genoData;
            GenoCursor &cursor = *ctx.cursor;
            MethodBase &meth = *ctx.method;
            const auto &localChunks = gd.chunkIndices();
            const std::string &naSuffix = ctx.naSuffix;

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

                const auto &gIdx = localChunks[cidx];
                std::string out;
                out.reserve(gIdx.size() * 256);

                if (!gIdx.empty()) cursor.beginSequentialBlock(gIdx.front());
                meth.prepareChunk(gIdx);

                for (size_t i = 0; i < gIdx.size(); ++i) {
                    double altFreq = NAN, altCounts = NAN;
                    double missingRate = NAN, hweP = NAN;
                    double maf = NAN, mac = NAN;
                    indexForMissing.clear();

                    cursor.getGenotypes(gIdx[i], GVec, altFreq, altCounts, missingRate, hweP, maf, mac,
                                        indexForMissing);

                    const std::string_view chr = gd.chr(gIdx[i]);
                    const std::string_view ref = gd.ref(gIdx[i]);
                    const std::string_view alt = gd.alt(gIdx[i]);
                    const std::string_view marker = gd.markerId(gIdx[i]);
                    const uint32_t pos = gd.pos(gIdx[i]);

                    const bool passQC = !(missingRate > missingCutoff || maf < minMafCutoff || mac < minMacCutoff ||
                                          (hweCutoff > 0 && hweP < hweCutoff));

                    if (!passQC) {
                        formatLineNA(out, fmtBuf, chr, pos, marker, alt /*REF=bim6*/, ref /*ALT=bim5*/, missingRate,
                                     altFreq, mac, hweP, naSuffix);
                        continue;
                    }

                    // Impute missing genotypes with mean (2 × altFreq).
                    const double imputeG = 2.0 * altFreq;
                    double *gPtr = GVec.data();
                    const uint32_t nMissing = static_cast<uint32_t>(indexForMissing.size());
                    for (uint32_t j = 0; j < nMissing; ++j)
                        gPtr[indexForMissing[j]] = imputeG;

                    rv.clear();
                    meth.getResultVec(GVec, altFreq, static_cast<int>(i), rv);

                    formatLine(out, fmtBuf, chr, pos, marker, alt /*REF=bim6*/, ref /*ALT=bim5*/, missingRate, altFreq,
                               mac, hweP, rv);
                }

                {
                    std::lock_guard<std::mutex> lk(writeMutex);
                    chunkOutput[cidx] = std::move(out);
                    chunkReady[cidx].ready = 1;
                }
                infoMsg("Calculation finished: chunk %zu/%zu", cidx + 1, nTotalChunks);
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
        for (auto &th : workers)
            th.join();
    }

    {
        std::lock_guard<std::mutex> lk(writeMutex);
        stopWriter = true;
    }
    writeCv.notify_all();
    writerThread.join();

    if (workerError) std::rethrow_exception(workerError);

    infoMsg("Output written to: %s", outputFile.c_str());
}

// ──────────────────────────────────────────────────────────────────────
// MultiMethod implementation
// ──────────────────────────────────────────────────────────────────────

MultiMethod::MultiMethod(
    std::vector<std::unique_ptr<MethodBase> > methods,
    std::vector<std::string> residNames,
    std::vector<std::string> suffixes
)
    : m_methods(std::move(methods)),
      m_residNames(std::move(residNames)),
      m_suffixes(std::move(suffixes))
{
}

std::unique_ptr<MethodBase> MultiMethod::clone() const {
    std::vector<std::unique_ptr<MethodBase> > cloned;
    cloned.reserve(m_methods.size());
    for (const auto &m : m_methods)
        cloned.push_back(m->clone());
    return std::make_unique<MultiMethod>(std::move(cloned), m_residNames, m_suffixes);
}

int MultiMethod::resultSize() const {
    return static_cast<int>(m_methods.size() * m_suffixes.size());
}

std::string MultiMethod::getHeaderColumns() const {
    std::string h;
    for (const auto &name : m_residNames)
        for (const auto &suf : m_suffixes)
            h += "\t" + name + suf;
    return h;
}

void MultiMethod::prepareChunk(const std::vector<uint64_t> &gIndices) {
    for (auto &m : m_methods)
        m->prepareChunk(gIndices);
}

void MultiMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int markerInChunkIdx,
    std::vector<double> &result
) {
    std::vector<double> inner;
    for (auto &m : m_methods) {
        inner.clear();
        m->getResultVec(GVec, altFreq, markerInChunkIdx, inner);
        result.insert(result.end(), inner.begin(), inner.end());
    }
}

// ══════════════════════════════════════════════════════════════════════
// multiPhenoEngine — per-phenotype independent GWAS
// ══════════════════════════════════════════════════════════════════════

void multiPhenoEngine(
    const GenoMeta &genoData,
    std::vector<PhenoTask> &tasks,
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
    const size_t K = tasks.size();
    const size_t nTotalChunks = genoData.chunkIndices().size();
    const int effective_nthreads = std::min(nthreads, static_cast<int>(nTotalChunks));
    const uint32_t nUnion = genoData.nSubjUsed();

    infoMsg("Number of subjects in the input file: %u", genoData.nSubjInFile());
    infoMsg("Number of subjects in union mask: %u", nUnion);
    for (size_t p = 0; p < K; ++p)
        infoMsg("  Phenotype '%s': %u subjects", tasks[p].phenoName.c_str(), tasks[p].nUsed);
    infoMsg("Number of markers in the input file: %u", genoData.nMarkers());
    infoMsg("Number of markers to test: %zu", genoData.markerInfo().size());
    infoMsg("Number of chunks for all markers: %zu", nTotalChunks);
    infoMsg("Start chunk-level parallel processing with %d worker threads, %zu phenotypes.", effective_nthreads, K);

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
        outPaths[p] = TextWriter::buildOutputPath(outPrefix, tasks[p].phenoName, methodName, compression);

    // Per-chunk, per-phenotype output buffers + ready flags.
    // chunkOutput[chunk][pheno]
    std::vector<std::vector<std::string> > chunkOutput(nTotalChunks, std::vector<std::string>(K));
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

            for (auto &w : writers)
                w.close();
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
            std::vector<std::unique_ptr<MethodBase> > methods(K);
            for (size_t p = 0; p < K; ++p)
                methods[p] = tasks[p].method->clone();

            Eigen::VectorXd GVec_union(nUnion);
            std::vector<Eigen::VectorXd> GVec_phenos(K);
            for (size_t p = 0; p < K; ++p)
                GVec_phenos[p].resize(tasks[p].nUsed);
            std::vector<uint32_t> unionMissing;
            unionMissing.reserve(nUnion / 10);

            // Batched extract+stats buffers (allocated once, reused per marker)
            std::vector<const uint32_t *> utlPtrs(K);
            std::vector<uint32_t> nPhenoArr(K);
            std::vector<double *> phenoGPtrs(K);
            for (size_t p = 0; p < K; ++p) {
                utlPtrs[p] = tasks[p].unionToLocal.data();
                nPhenoArr[p] = tasks[p].nUsed;
                phenoGPtrs[p] = GVec_phenos[p].data();
            }
            std::vector<PhenoGenoStats> batchStats(K);
            std::vector<std::vector<uint32_t> > batchMissings(K);
            for (size_t p = 0; p < K; ++p)
                batchMissings[p].reserve(tasks[p].nUsed / 10);

            std::vector<double> rv;
            rv.reserve(16);
            char fmtBuf[32];

            const auto &localChunks = genoData.chunkIndices();

            while (!stopWriter.load(std::memory_order_relaxed)) {
                const size_t cidx = nextChunk.fetch_add(1);
                if (cidx >= localChunks.size()) break;

                const auto &gIdx = localChunks[cidx];

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
                    cursor->getGenotypes(gIdx[i], GVec_union, uAltFreq, uAltCounts, uMissRate, uHweP, uMaf, uMac,
                                         unionMissing);

                    const std::string_view chr = genoData.chr(gIdx[i]);
                    const std::string_view ref = genoData.ref(gIdx[i]);
                    const std::string_view alt = genoData.alt(gIdx[i]);
                    const std::string_view marker = genoData.markerId(gIdx[i]);
                    const uint32_t pos = genoData.pos(gIdx[i]);

                    // 2. Fused: extract all K phenotypes + compute stats in one pass
                    extractAndStatsBatched(
                        GVec_union.data(), nUnion, K,
                        utlPtrs.data(), nPhenoArr.data(),
                        phenoGPtrs.data(), batchStats.data(), batchMissings.data());

                    for (size_t p = 0; p < K; ++p) {
                        const PhenoGenoStats &gs = batchStats[p];

                        double pMaf = std::min(gs.altFreq, 1.0 - gs.altFreq);
                        const bool passQC = !(gs.missingRate > missingCutoff || pMaf < minMafCutoff ||
                                              gs.mac < minMacCutoff || (hweCutoff > 0 && gs.hweP < hweCutoff));

                        if (!passQC) {
                            formatLineNA(phenoOut[p], fmtBuf, chr, pos, marker, alt, ref, gs.missingRate, gs.altFreq,
                                         gs.mac, gs.hweP, naSuffixes[p]);
                            continue;
                        }

                        // Impute missing with mean
                        const double imputeG = 2.0 * gs.altFreq;
                        double *gPtr = GVec_phenos[p].data();
                        for (uint32_t j : batchMissings[p])
                            gPtr[j] = imputeG;

                        rv.clear();
                        methods[p]->getResultVec(GVec_phenos[p].head(tasks[p].nUsed), gs.altFreq, static_cast<int>(i),
                                                 rv);

                        formatLine(phenoOut[p], fmtBuf, chr, pos, marker, alt, ref, gs.missingRate, gs.altFreq, gs.mac,
                                   gs.hweP, rv);
                    }
                }

                {
                    std::lock_guard<std::mutex> lk(writeMutex);
                    for (size_t p = 0; p < K; ++p)
                        chunkOutput[cidx][p] = std::move(phenoOut[p]);
                    chunkReady[cidx].ready = 1;
                }
                infoMsg("Calculation finished: chunk %zu/%zu", cidx + 1, nTotalChunks);
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
        for (auto &th : workers)
            th.join();
    }

    {
        std::lock_guard<std::mutex> lk(writeMutex);
        stopWriter = true;
    }
    writeCv.notify_all();
    writerThread.join();

    if (workerError) std::rethrow_exception(workerError);

    for (size_t p = 0; p < K; ++p)
        infoMsg("Output written to: %s", outPaths[p].c_str());
}
