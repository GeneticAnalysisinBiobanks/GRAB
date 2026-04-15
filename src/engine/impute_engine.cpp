// impute_engine.cpp — Impute-mode multi-phenotype engine
//
// All phenotypes share one union GVec.  Missing phenotype residuals
// were 0-padded by padToUnionSpace(), so score = resid · GVec naturally
// skips absent subjects.  This eliminates K per-phenotype GVec buffers
// and extraction loops, dramatically reducing memory and cache pressure.

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
// imputeMultiPhenoEngine — wrapper for impute mode
//
// Pads every task's method to union space, then delegates to
// imputeMultiPhenoEngineRange.  Same interface as multiPhenoEngine.
// ══════════════════════════════════════════════════════════════════════

void imputeMultiPhenoEngine(
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
    const uint32_t nUnion = genoData.nSubjUsed();

    // Pad each task's method to union space.
    for (size_t p = 0; p < K; ++p) {
        tasks[p].method->padToUnionSpace(tasks[p].unionToLocal.data(), nUnion);
        tasks[p].nUsed = nUnion;
    }

    infoMsg("Number of subjects in the input file: %u", genoData.nSubjInFile());
    infoMsg("Number of subjects in union mask: %u", nUnion);
    infoMsg("Impute mode: all %zu phenotypes share union genotype vector", K);
    for (size_t p = 0; p < K; ++p)
        infoMsg("  Phenotype '%s': union-padded to %u subjects", tasks[p].phenoName.c_str(), nUnion);
    infoMsg("Number of markers in the input file: %u", genoData.nMarkers());
    infoMsg("Number of markers to test: %zu", genoData.markerInfo().size());
    infoMsg("Number of chunks for all markers: %zu", nTotalChunks);
    infoMsg("Start chunk-level parallel processing with %d worker threads, %zu phenotypes (impute).",
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

    // Delegate to the impute range engine.
    imputeMultiPhenoEngineRange(genoData, tasks, naSuffixes, 0, nTotalChunks, writers,
                                nthreads, missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);

    for (auto &w : writers)
        w.close();

    for (size_t p = 0; p < K; ++p)
        infoMsg("Output written to: %s", outPaths[p].c_str());
}

// ══════════════════════════════════════════════════════════════════════
// imputeMultiPhenoEngineRange — union GVec for impute mode
//
// All phenotypes share one union GVec.  Missing phenotype residuals
// were 0-padded by padToUnionSpace(), so score = resid · GVec naturally
// skips absent subjects.  This eliminates K per-phenotype GVec buffers
// and extraction loops, dramatically reducing memory and cache pressure.
// ══════════════════════════════════════════════════════════════════════

void imputeMultiPhenoEngineRange(
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

            // Single union GVec — shared across ALL phenotypes
            Eigen::VectorXd GVec_union(nUnion);

            // Single stats + missing buffer (union-level)
            std::vector<uint32_t> indexForMissing;
            indexForMissing.reserve(nUnion / 10);

            std::vector<double> rv;
            rv.reserve(16);
            char fmtBuf[32];

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

                    // Decode genotypes into union GVec
                    cursor->getGenotypesSimple(gIdx[i], GVec_union);

                    // Union-level stats (ONE call for all K phenotypes)
                    PhenoGenoStats gs = statsFromGVec(
                        GVec_union.data(), nUnion, indexForMissing);

                    double pMaf = std::min(gs.altFreq, 1.0 - gs.altFreq);
                    const bool passQC = !(gs.missingRate > missingCutoff || pMaf < minMafCutoff ||
                                          gs.mac < minMacCutoff || (hweCutoff > 0 && gs.hweP < hweCutoff));

                    if (!passQC) {
                        for (size_t p = 0; p < K; ++p)
                            formatLineNA(phenoOut[p], fmtBuf, chr, pos, marker, alt, ref,
                                         gs.missingRate, gs.altFreq, gs.mac, gs.hweP, naSuffixes[p]);
                        continue;
                    }

                    // Impute missing genotypes with mean (one pass)
                    const double imputeG = 2.0 * gs.altFreq;
                    double *gPtr = GVec_union.data();
                    for (uint32_t j : indexForMissing)
                        gPtr[j] = imputeG;

                    // Run each phenotype's method on the SAME union GVec
                    for (size_t p = 0; p < K; ++p) {
                        rv.clear();
                        methods[p]->getResultVec(GVec_union.head(nUnion), gs.altFreq,
                                                 static_cast<int>(i), rv);

                        formatLine(phenoOut[p], fmtBuf, chr, pos, marker, alt, ref,
                                   gs.missingRate, gs.altFreq, gs.mac, gs.hweP, rv);
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
