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
// multiPhenoEngineRange — shared worker+writer loop for a chunk range
// ══════════════════════════════════════════════════════════════════════

void multiPhenoEngineRange(
    const GenoMeta &genoData,
    std::vector<PhenoTask> &tasks,
    const std::vector<std::string> &naSuffixes,
    size_t chunkStart,
    size_t chunkEnd,
    std::vector<TextWriter> &writers,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff
) {
    const size_t K = tasks.size();
    const size_t nChunks = chunkEnd - chunkStart;
    const int effective_nthreads = std::min(nthreads, static_cast<int>(nChunks));
    const uint32_t nUnion = genoData.nSubjUsed();
    const auto &allChunks = genoData.chunkIndices();

    // ── Missingness-pattern batching ───────────────────────────────────
    // Group phenotypes with identical unionToLocal (= same non-missing
    // sample set).  All phenotypes in a batch share ONE genotype extraction
    // and ONE stats computation, eliminating (K - nBatches) scattered
    // extract passes and stats passes per marker.
    struct MissBatch {
        size_t rep;                  // representative phenotype index
        std::vector<size_t> members; // all phenotype indices with this pattern
    };

    std::vector<MissBatch> missBatches;
    {
        std::vector<bool> assigned(K, false);
        for (size_t p = 0; p < K; ++p) {
            if (assigned[p]) continue;
            MissBatch batch;
            batch.rep = p;
            batch.members.push_back(p);
            assigned[p] = true;
            for (size_t q = p + 1; q < K; ++q) {
                if (assigned[q]) continue;
                if (tasks[p].nUsed == tasks[q].nUsed &&
                    tasks[p].unionToLocal == tasks[q].unionToLocal) {
                    batch.members.push_back(q);
                    assigned[q] = true;
                }
            }
            missBatches.push_back(std::move(batch));
        }
    }
    const size_t nBatches = missBatches.size();
    if (nBatches < K)
        infoMsg("Missingness-pattern batching: %zu unique patterns for %zu phenotypes", nBatches, K);

    // Per-chunk, per-phenotype output buffers + ready flags.
    std::vector<std::vector<std::string> > chunkOutput(nChunks, std::vector<std::string>(K));
    std::vector<PaddedFlag> chunkReady(nChunks, {0});
    std::atomic<size_t> nextChunk(0);

    std::exception_ptr workerError = nullptr;
    std::mutex errorMutex;
    std::mutex writeMutex;
    std::condition_variable writeCv;
    std::atomic<bool> stopFlag{false};

    // ── Writer thread ──────────────────────────────────────────────────
    std::thread writerThread([&]() {
        try {
            for (size_t ci = 0; ci < nChunks; ++ci) {
                {
                    std::unique_lock<std::mutex> lk(writeMutex);
                    writeCv.wait(lk, [&]() {
                        return chunkReady[ci].ready || stopFlag.load();
                    });
                    if (!chunkReady[ci].ready) break;
                }
                for (size_t p = 0; p < K; ++p)
                    writers[p].write(chunkOutput[ci][p]);
                infoMsg("Writing finished: chunk %zu/%zu", ci + 1, nChunks);
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(errorMutex);
            if (!workerError) workerError = std::current_exception();
            stopFlag = true;
        }
    });

    // ── Worker function ────────────────────────────────────────────────
    auto workerFn = [&]() {
        try {
            auto cursor = genoData.makeCursor();
            std::vector<std::unique_ptr<MethodBase> > methods(K);
            for (size_t p = 0; p < K; ++p)
                methods[p] = tasks[p].method->clone();

            Eigen::VectorXd GVec_union(nUnion);
            std::vector<Eigen::VectorXd> GVec_phenos(K);
            for (size_t p = 0; p < K; ++p)
                GVec_phenos[p].resize(tasks[p].nUsed);

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

            // Difflist buffers for sparse genotype path
            const uint32_t maxDiffLen = nUnion / 8;
            std::vector<uint32_t> diffSampleIds(maxDiffLen + 2);
            std::vector<uint8_t> diffGenoCodes(maxDiffLen + 2);

            while (!stopFlag.load(std::memory_order_relaxed)) {
                const size_t localIdx = nextChunk.fetch_add(1);
                if (localIdx >= nChunks) break;

                const size_t globalIdx = chunkStart + localIdx;
                const auto &gIdx = allChunks[globalIdx];

                // Per-phenotype output buffers for this chunk.
                std::vector<std::string> phenoOut(K);
                for (size_t p = 0; p < K; ++p)
                    phenoOut[p].reserve(gIdx.size() * 128);

                if (!gIdx.empty()) cursor->beginSequentialBlock(gIdx.front());
                for (size_t p = 0; p < K; ++p)
                    methods[p]->prepareChunk(gIdx);

                for (size_t i = 0; i < gIdx.size(); ++i) {
                    const std::string_view chr = genoData.chr(gIdx[i]);
                    const std::string_view ref = genoData.ref(gIdx[i]);
                    const std::string_view alt = genoData.alt(gIdx[i]);
                    const std::string_view marker = genoData.markerId(gIdx[i]);
                    const uint32_t pos = genoData.pos(gIdx[i]);

                    // Try sparse read first; falls back to dense internally.
                    uint32_t diffLen = 0;
                    const uint32_t commonGeno = cursor->getGenotypesMaybeSparse(
                        gIdx[i], GVec_union, maxDiffLen,
                        diffSampleIds.data(), diffGenoCodes.data(), diffLen);

                    if (commonGeno != UINT32_MAX) {
                        // ── Sparse path ────────────────────────────────────
                        for (size_t b = 0; b < nBatches; ++b) {
                            const auto &batch = missBatches[b];
                            const size_t rep = batch.rep;
                            batchStats[rep] = sparseExtractAndStats(
                                phenoGPtrs[rep], nPhenoArr[rep], utlPtrs[rep],
                                commonGeno, diffSampleIds.data(),
                                diffGenoCodes.data(), diffLen, batchMissings[rep]);
                            for (size_t mi = 1; mi < batch.members.size(); ++mi) {
                                const size_t p = batch.members[mi];
                                std::copy(phenoGPtrs[rep],
                                          phenoGPtrs[rep] + nPhenoArr[rep],
                                          phenoGPtrs[p]);
                                batchStats[p] = batchStats[rep];
                                batchMissings[p] = batchMissings[rep];
                            }
                        }
                    } else {
                        // ── Dense path ─────────────────────────────────────
                        for (size_t b = 0; b < nBatches; ++b) {
                            const auto &batch = missBatches[b];
                            const size_t rep = batch.rep;
                            extractPhenoGVec(GVec_union.data(), nUnion,
                                             utlPtrs[rep], nPhenoArr[rep],
                                             phenoGPtrs[rep]);
                            batchStats[rep] = statsFromGVec(
                                phenoGPtrs[rep], nPhenoArr[rep],
                                batchMissings[rep]);
                            for (size_t mi = 1; mi < batch.members.size(); ++mi) {
                                const size_t p = batch.members[mi];
                                std::copy(phenoGPtrs[rep],
                                          phenoGPtrs[rep] + nPhenoArr[rep],
                                          phenoGPtrs[p]);
                                batchStats[p] = batchStats[rep];
                                batchMissings[p] = batchMissings[rep];
                            }
                        }
                    }

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
                        chunkOutput[localIdx][p] = std::move(phenoOut[p]);
                    chunkReady[localIdx].ready = 1;
                }
                infoMsg("Calculation finished: chunk %zu/%zu", localIdx + 1, nChunks);
                writeCv.notify_all();
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lk(errorMutex);
                if (!workerError) workerError = std::current_exception();
            }
            {
                std::lock_guard<std::mutex> lk(writeMutex);
                stopFlag = true;
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
        stopFlag = true;
    }
    writeCv.notify_all();
    writerThread.join();

    if (workerError) std::rethrow_exception(workerError);
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

    infoMsg("Number of subjects in the input file: %u", genoData.nSubjInFile());
    infoMsg("Number of subjects in union mask: %u", genoData.nSubjUsed());
    for (size_t p = 0; p < K; ++p)
        infoMsg("  Phenotype '%s': %u subjects", tasks[p].phenoName.c_str(), tasks[p].nUsed);
    infoMsg("Number of markers in the input file: %u", genoData.nMarkers());
    infoMsg("Number of markers to test: %zu", genoData.markerInfo().size());
    infoMsg("Number of chunks for all markers: %zu", nTotalChunks);
    infoMsg("Start chunk-level parallel processing with %d worker threads, %zu phenotypes.",
            std::min(nthreads, static_cast<int>(nTotalChunks)), K);

    // Build per-phenotype NA suffixes.
    std::vector<std::string> naSuffixes(K);
    for (size_t p = 0; p < K; ++p)
        naSuffixes[p] = makeNaSuffix(tasks[p].method->resultSize());

    // Build output file paths and open writers.
    auto writerMode = TextWriter::modeFromString(compression);
    std::vector<std::string> outPaths(K);
    std::vector<TextWriter> writers;
    writers.reserve(K);
    for (size_t p = 0; p < K; ++p) {
        outPaths[p] = TextWriter::buildOutputPath(outPrefix, tasks[p].phenoName, methodName, compression);
        writers.emplace_back(outPaths[p], writerMode, compressionLevel);
        std::string header = std::string(META_HEADER) + tasks[p].method->getHeaderColumns() + "\n";
        writers[p].write(header);
    }

    // Delegate to the range engine.
    multiPhenoEngineRange(genoData, tasks, naSuffixes, 0, nTotalChunks, writers,
                          nthreads, missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);

    for (auto &w : writers)
        w.close();

    for (size_t p = 0; p < K; ++p)
        infoMsg("Output written to: %s", outPaths[p].c_str());
}
