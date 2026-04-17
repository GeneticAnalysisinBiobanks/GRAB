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

            // ── Batch buffers (GEMV → GEMM optimisation) ───────────────
            const int batchB = meth.preferredBatchSize();

            struct MarkerMeta {
                double altFreq, missingRate, mac, hweP;
                std::string_view chr, ref, alt, marker;
                uint32_t pos;
                bool passQC;
            };

            Eigen::MatrixXd GBatch;
            std::vector<double> passAltFreqs;
            std::vector<std::vector<double> > batchResults;
            std::vector<MarkerMeta> batchMeta;

            if (batchB > 1) {
                GBatch.resize(gd.nSubjUsed(), batchB);
                passAltFreqs.reserve(batchB);
                batchMeta.resize(batchB);
            }

            while (!stopWriter.load(std::memory_order_relaxed)) {
                const size_t cidx = nextChunk.fetch_add(1);
                if (cidx >= localChunks.size()) break;

                const auto &gIdx = localChunks[cidx];
                std::string out;
                out.reserve(gIdx.size() * 256);

                if (!gIdx.empty()) cursor.beginSequentialBlock(gIdx.front());
                meth.prepareChunk(gIdx);

                if (batchB > 1) {
                    // ── Batched path: process batchB markers at a time ─────
                    const size_t B = static_cast<size_t>(batchB);
                    for (size_t bstart = 0; bstart < gIdx.size(); bstart += B) {
                        const size_t blen = std::min(B, gIdx.size() - bstart);

                        passAltFreqs.clear();
                        int passCount = 0;

                        for (size_t bi = 0; bi < blen; ++bi) {
                            const size_t i = bstart + bi;
                            double altFreq = NAN, altCounts = NAN;
                            double missingRate = NAN, hweP = NAN;
                            double maf = NAN, mac = NAN;
                            indexForMissing.clear();

                            cursor.getGenotypes(gIdx[i], GVec, altFreq, altCounts,
                                                missingRate, hweP, maf, mac, indexForMissing);

                            auto &mi = batchMeta[bi];
                            mi.altFreq = altFreq;
                            mi.missingRate = missingRate;
                            mi.mac = mac;
                            mi.hweP = hweP;
                            mi.chr = gd.chr(gIdx[i]);
                            mi.ref = gd.ref(gIdx[i]);
                            mi.alt = gd.alt(gIdx[i]);
                            mi.marker = gd.markerId(gIdx[i]);
                            mi.pos = gd.pos(gIdx[i]);
                            mi.passQC = !(missingRate > missingCutoff || maf < minMafCutoff ||
                                          mac < minMacCutoff || (hweCutoff > 0 && hweP < hweCutoff));

                            if (mi.passQC) {
                                const double imputeG = 2.0 * altFreq;
                                double *gPtr = GVec.data();
                                for (uint32_t j : indexForMissing)
                                    gPtr[j] = imputeG;
                                GBatch.col(passCount) = GVec;
                                passAltFreqs.push_back(altFreq);
                                ++passCount;
                            }
                        }

                        if (passCount > 0)
                            meth.getResultBatch(GBatch.leftCols(passCount), passAltFreqs, batchResults);

                        int ri = 0;
                        for (size_t bi = 0; bi < blen; ++bi) {
                            const auto &mi = batchMeta[bi];
                            if (!mi.passQC) {
                                formatLineNA(out, fmtBuf, mi.chr, mi.pos, mi.marker,
                                             mi.alt, mi.ref, mi.missingRate, mi.altFreq,
                                             mi.mac, mi.hweP, naSuffix);
                            } else {
                                formatLine(out, fmtBuf, mi.chr, mi.pos, mi.marker,
                                           mi.alt, mi.ref, mi.missingRate, mi.altFreq,
                                           mi.mac, mi.hweP, batchResults[ri++]);
                            }
                        }
                    }
                } else {
                    // ── Original single-marker path ────────────────────────
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

                        const bool passQC = !(missingRate > missingCutoff || maf < minMafCutoff ||
                                              mac < minMacCutoff || (hweCutoff > 0 && hweP < hweCutoff));

                        if (!passQC) {
                            formatLineNA(out, fmtBuf, chr, pos, marker, alt /*REF=bim6*/, ref /*ALT=bim5*/,
                                         missingRate, altFreq, mac, hweP, naSuffix);
                            continue;
                        }

                        const double imputeG = 2.0 * altFreq;
                        double *gPtr = GVec.data();
                        const uint32_t nMissing = static_cast<uint32_t>(indexForMissing.size());
                        for (uint32_t j = 0; j < nMissing; ++j)
                            gPtr[indexForMissing[j]] = imputeG;

                        rv.clear();
                        meth.getResultVec(GVec, altFreq, static_cast<int>(i), rv);

                        formatLine(out, fmtBuf, chr, pos, marker, alt /*REF=bim6*/, ref /*ALT=bim5*/,
                                   missingRate, altFreq, mac, hweP, rv);
                    }
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

void MultiMethod::getResultBatch(
    const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
    const std::vector<double> &altFreqs,
    std::vector<std::vector<double> > &results
) {
    const int B = static_cast<int>(GBatch.cols());
    results.resize(B);
    for (int b = 0; b < B; ++b)
        results[b].clear();

    std::vector<std::vector<double> > innerResults;
    for (auto &m : m_methods) {
        m->getResultBatch(GBatch, altFreqs, innerResults);
        for (int b = 0; b < B; ++b)
            results[b].insert(results[b].end(), innerResults[b].begin(), innerResults[b].end());
    }
}

int MultiMethod::preferredBatchSize() const {
    int best = 0;
    for (const auto &m : m_methods)
        best = std::max(best, m->preferredBatchSize());
    return best;
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

    // ── Fused union-level GEMM setup ───────────────────────────────────
    // Methods that support fused GEMM (Score = resid^T × G) contribute
    // their residual columns to a shared AugResid matrix.  The engine
    // does ONE GEMM per window for ALL fuseable phenotypes, eliminating
    // per-phenotype extraction entirely.
    //
    // AugResid layout (N_union × totalFusedCols + nFuseable):
    //   [residCols_0 | ... | residCols_{K-1} | mask_0 | ... | mask_{K-1}]
    //
    // The last nFuseable columns are 0/1 mask vectors for computing
    // per-phenotype genotype sums (gSum = mask^T × G).

    struct FusedPhenoInfo {
        size_t taskIdx;        // index into tasks[]
        int colOffset;         // start column in AugResid (residual columns)
        int nCols;             // fusedGemmColumns() for this phenotype
        int maskCol;           // column index for this phenotype's mask in AugResid
        uint32_t nUsed;        // this phenotype's sample count
    };

    std::vector<FusedPhenoInfo> fusedPhenos;
    std::vector<size_t> nonFusedPhenos;       // task indices of non-fuseable phenotypes
    int totalFusedCols = 0;

    for (size_t p = 0; p < K; ++p) {
        if (tasks[p].method->supportsFusedGemm()) {
            FusedPhenoInfo fi;
            fi.taskIdx = p;
            fi.colOffset = totalFusedCols;
            fi.nCols = tasks[p].method->fusedGemmColumns();
            fi.nUsed = tasks[p].nUsed;
            fi.maskCol = 0;  // set below
            totalFusedCols += fi.nCols;
            fusedPhenos.push_back(fi);
        } else {
            nonFusedPhenos.push_back(p);
        }
    }

    const size_t nFuseable = fusedPhenos.size();
    const int augCols = totalFusedCols + static_cast<int>(nFuseable);

    // Set mask column indices.
    for (size_t fi = 0; fi < nFuseable; ++fi)
        fusedPhenos[fi].maskCol = totalFusedCols + static_cast<int>(fi);

    // Build shared AugResid matrix (read-only, shared across all threads).
    // Also store per-fuseable residual sums.
    Eigen::MatrixXd AugResid;
    std::vector<double> allResidSums;    // length totalFusedCols

    const bool hasFused = nFuseable > 0;

    if (hasFused) {
        AugResid.setZero(nUnion, augCols);
        allResidSums.resize(totalFusedCols, 0.0);

        for (size_t fi = 0; fi < nFuseable; ++fi) {
            const auto &fp = fusedPhenos[fi];
            const size_t p = fp.taskIdx;

            // Fill residual columns (zero-padded for absent samples).
            tasks[p].method->fillUnionResiduals(
                AugResid.middleCols(fp.colOffset, fp.nCols),
                tasks[p].unionToLocal);

            // Fill residual sums.
            tasks[p].method->fillResidualSums(allResidSums.data() + fp.colOffset);

            // Fill mask column (1.0 for present, 0.0 for absent).
            for (uint32_t i = 0; i < nUnion; ++i)
                if (tasks[p].unionToLocal[i] != UINT32_MAX)
                    AugResid(i, fp.maskCol) = 1.0;
        }

        infoMsg("Fused GEMM: %zu fuseable phenotypes (%d residual cols + %zu mask cols = %d augmented cols), "
                "%zu non-fuseable phenotypes",
                nFuseable, totalFusedCols, nFuseable, augCols, nonFusedPhenos.size());
    }

    // ── Group fuseable phenotypes by identical subject sets (D1) ───────
    // Phenotypes sharing the same unionToLocal (same subjects) produce
    // identical statsFromUnionVec results, so we compute QC once per group.
    struct FusedStatsGroup {
        size_t repFi;                   // representative index in fusedPhenos[]
        std::vector<size_t> fiIndices;  // all fusedPhenos[] indices in this group
    };

    std::vector<FusedStatsGroup> fusedGroups;
    if (hasFused) {
        std::vector<bool> grouped(nFuseable, false);
        for (size_t fi = 0; fi < nFuseable; ++fi) {
            if (grouped[fi]) continue;
            FusedStatsGroup g;
            g.repFi = fi;
            g.fiIndices.push_back(fi);
            grouped[fi] = true;
            const size_t pi = fusedPhenos[fi].taskIdx;
            for (size_t fj = fi + 1; fj < nFuseable; ++fj) {
                if (grouped[fj]) continue;
                const size_t pj = fusedPhenos[fj].taskIdx;
                if (tasks[pi].nUsed == tasks[pj].nUsed &&
                    tasks[pi].unionToLocal == tasks[pj].unionToLocal) {
                    g.fiIndices.push_back(fj);
                    grouped[fj] = true;
                }
            }
            fusedGroups.push_back(std::move(g));
        }
        if (fusedGroups.size() < nFuseable)
            infoMsg("Fused QC groups: %zu groups for %zu fuseable phenotypes (%.0f%% QC savings)",
                    fusedGroups.size(), nFuseable,
                    100.0 * (1.0 - static_cast<double>(fusedGroups.size()) / static_cast<double>(nFuseable)));
    }

    // ── Missingness-pattern batching (non-fuseable phenotypes only) ────
    struct MissBatch {
        size_t rep;
        std::vector<size_t> members;
        std::vector<uint32_t> presentIndices;
        std::vector<uint64_t> presentMask;
        bool identity;
    };

    std::vector<MissBatch> missBatches;
    if (!nonFusedPhenos.empty()) {
        std::vector<bool> assigned(K, false);
        // Mark all fuseable phenotypes as assigned.
        for (const auto &fp : fusedPhenos)
            assigned[fp.taskIdx] = true;

        for (size_t ni = 0; ni < nonFusedPhenos.size(); ++ni) {
            const size_t p = nonFusedPhenos[ni];
            if (assigned[p]) continue;
            MissBatch batch;
            batch.rep = p;
            batch.members.push_back(p);
            assigned[p] = true;
            for (size_t nj = ni + 1; nj < nonFusedPhenos.size(); ++nj) {
                const size_t q = nonFusedPhenos[nj];
                if (assigned[q]) continue;
                if (tasks[p].nUsed == tasks[q].nUsed &&
                    tasks[p].unionToLocal == tasks[q].unionToLocal) {
                    batch.members.push_back(q);
                    assigned[q] = true;
                }
            }
            const uint32_t *utl = tasks[p].unionToLocal.data();
            batch.identity = isIdentityMapping(utl, nUnion, tasks[p].nUsed);
            if (!batch.identity) {
                batch.presentIndices = buildPresentIndices(utl, nUnion);
                batch.presentMask = buildPresentMask(utl, nUnion);
            }
            missBatches.push_back(std::move(batch));
        }
        infoMsg("Non-fuseable MissBatch: %zu unique patterns for %zu non-fuseable phenotypes",
                missBatches.size(), nonFusedPhenos.size());
    }
    const size_t nMissBatches = missBatches.size();

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
            char fmtBuf[32];

            // Difflist buffers for sparse genotype path
            const uint32_t maxDiffLen = nUnion / 8;
            std::vector<uint32_t> diffSampleIds(maxDiffLen + 2);
            std::vector<uint8_t> diffGenoCodes(maxDiffLen + 2);

            // ── Fused GEMM buffers ─────────────────────────────────────
            // Batch size: use preferredBatchSize from methods, default 16.
            int commonB = 16;
            for (size_t p = 0; p < K; ++p)
                commonB = std::max(commonB, methods[p]->preferredBatchSize());

            struct WindowMarkerMeta {
                std::string_view chr, ref, alt, marker;
                uint32_t pos;
            };

            const size_t B = static_cast<size_t>(commonB);

            // Union-level batch matrix: ONE per thread (N_union × B).
            Eigen::MatrixXd GBatch_union(nUnion, commonB);

            // Fused GEMM result: augCols × B (scores + gSums). Tiny.
            Eigen::MatrixXd allScoresAndSums;
            if (hasFused)
                allScoresAndSums.resize(augCols, commonB);

            std::vector<WindowMarkerMeta> winMarkers(B);

            // Per-fuseable-phenotype: QC stats per window marker.
            // winGSums[fi][bi] = genotype sum for fuseable pheno fi, marker bi
            // (from mask column dot product).
            std::vector<std::vector<double> > winGSums;
            std::vector<double> passAltFreqs;
            std::vector<double> passGSums;
            std::vector<std::vector<double> > fusedResultsBuf;

            if (hasFused) {
                winGSums.resize(nFuseable, std::vector<double>(B, 0.0));
                passAltFreqs.reserve(B);
                passGSums.reserve(B);
            }

            // Pre-allocated passScores buffer (A2): reused across windows.
            int maxFusedNCols = 0;
            Eigen::MatrixXd passScoresBuf;
            if (hasFused) {
                for (size_t fi = 0; fi < nFuseable; ++fi)
                    maxFusedNCols = std::max(maxFusedNCols, fusedPhenos[fi].nCols);
                passScoresBuf.resize(maxFusedNCols, commonB);
            }

            // ── Non-fuseable fallback buffers ──────────────────────────
            uint32_t maxNonFusedN = 0;
            std::vector<Eigen::VectorXd> GVec_phenos_nf;       // per non-fuseable phenotype
            std::vector<uint32_t> nfNPhenoArr;                  // nUsed for each
            std::vector<double *> nfPhenoGPtrs;
            std::vector<PhenoGenoStats> nfBatchStats;
            std::vector<std::vector<uint32_t> > nfBatchMissings;
            Eigen::MatrixXd GBatch_pheno_nf;
            std::vector<double> nfPassAltFreqs;
            std::vector<std::vector<double> > nfBatchResultsBuf;
            std::vector<PhenoGenoStats> nfWinStats;
            std::vector<char> nfWinPassQC;

            if (!nonFusedPhenos.empty()) {
                GVec_phenos_nf.resize(nonFusedPhenos.size());
                nfNPhenoArr.resize(nonFusedPhenos.size());
                nfPhenoGPtrs.resize(nonFusedPhenos.size());
                nfBatchStats.resize(nonFusedPhenos.size());
                nfBatchMissings.resize(nonFusedPhenos.size());

                for (size_t ni = 0; ni < nonFusedPhenos.size(); ++ni) {
                    const size_t p = nonFusedPhenos[ni];
                    nfNPhenoArr[ni] = tasks[p].nUsed;
                    GVec_phenos_nf[ni].resize(tasks[p].nUsed);
                    nfPhenoGPtrs[ni] = GVec_phenos_nf[ni].data();
                    nfBatchMissings[ni].reserve(tasks[p].nUsed / 10);
                    maxNonFusedN = std::max(maxNonFusedN, tasks[p].nUsed);
                }

                GBatch_pheno_nf.resize(maxNonFusedN, commonB);
                nfPassAltFreqs.reserve(B);
                nfWinStats.resize(B);
                nfWinPassQC.resize(B, 0);
            }

            // Shared result + format buffers.
            std::vector<double> rv;
            rv.reserve(16);

            while (!stopFlag.load(std::memory_order_relaxed)) {
                const size_t localIdx = nextChunk.fetch_add(1);
                if (localIdx >= nChunks) break;

                const size_t globalIdx = chunkStart + localIdx;
                const auto &gIdx = allChunks[globalIdx];

                std::vector<std::string> phenoOut(K);
                for (size_t p = 0; p < K; ++p)
                    phenoOut[p].reserve(gIdx.size() * 128);

                if (!gIdx.empty()) cursor->beginSequentialBlock(gIdx.front());
                for (size_t p = 0; p < K; ++p)
                    methods[p]->prepareChunk(gIdx);

                // ── Window loop ────────────────────────────────────────
                for (size_t wstart = 0; wstart < gIdx.size(); wstart += B) {
                    const size_t wlen = std::min(B, gIdx.size() - wstart);

                    // Phase 1: Decode wlen markers into GBatch_union.
                    for (size_t bi = 0; bi < wlen; ++bi) {
                        const size_t i = wstart + bi;
                        winMarkers[bi].chr = genoData.chr(gIdx[i]);
                        winMarkers[bi].ref = genoData.ref(gIdx[i]);
                        winMarkers[bi].alt = genoData.alt(gIdx[i]);
                        winMarkers[bi].marker = genoData.markerId(gIdx[i]);
                        winMarkers[bi].pos = genoData.pos(gIdx[i]);

                        uint32_t diffLen = 0;
                        const uint32_t commonGeno = cursor->getGenotypesMaybeSparse(
                            gIdx[i], GVec_union, maxDiffLen,
                            diffSampleIds.data(), diffGenoCodes.data(), diffLen);

                        if (commonGeno != UINT32_MAX) {
                            auto col = GBatch_union.col(static_cast<Eigen::Index>(bi));
                            col.setConstant(static_cast<double>(commonGeno));
                            for (uint32_t j = 0; j < diffLen; ++j) {
                                const uint32_t sid = diffSampleIds[j];
                                if (sid < nUnion) {
                                    const uint8_t gc = diffGenoCodes[j];
                                    col[sid] = (gc == 3)
                                        ? std::numeric_limits<double>::quiet_NaN()
                                        : static_cast<double>(gc);
                                }
                            }
                        } else {
                            GBatch_union.col(static_cast<Eigen::Index>(bi)) = GVec_union;
                        }
                    }

                    // Phase 1b: Impute NaN genotypes to 2*AF in GBatch_union.
                    // This prevents NaN from propagating through the fused GEMM.
                    // For imputed data (dosages) there are no NaNs → this is a no-op.
                    // For BED data with missingness, this matches the standard
                    // mean-imputation approach used by the non-fused path.
                    if (hasFused) {
                        for (size_t bi = 0; bi < wlen; ++bi) {
                            auto col = GBatch_union.col(static_cast<Eigen::Index>(bi));
                            double gSum = 0.0;
                            uint32_t nValid = 0;
                            for (Eigen::Index r = 0; r < static_cast<Eigen::Index>(nUnion); ++r) {
                                const double v = col[r];
                                if (!std::isnan(v)) {
                                    gSum += v;
                                    ++nValid;
                                }
                            }
                            if (nValid > 0 && nValid < nUnion) {
                                const double impVal = gSum / static_cast<double>(nValid);
                                for (Eigen::Index r = 0; r < static_cast<Eigen::Index>(nUnion); ++r) {
                                    if (std::isnan(col[r]))
                                        col[r] = impVal;
                                }
                            }
                        }
                    }

                    // Phase 2: Fused GEMM for all fuseable phenotypes.
                    if (hasFused) {
                        const auto wlenI = static_cast<Eigen::Index>(wlen);

                        // ONE GEMM: AugResid^T × GBatch_union → (augCols × wlen)
                        // Rows 0..totalFusedCols-1 = raw scores
                        // Rows totalFusedCols..augCols-1 = gSums (mask^T × G)
                        allScoresAndSums.leftCols(wlenI).noalias() =
                            AugResid.transpose() * GBatch_union.leftCols(wlenI);

                        // Phase 3: Per fuseable-phenotype group — shared QC + processScoreBatch.
                        // D1: Phenotypes with identical subjects share statsFromUnionVec results.
                        struct FusedMarkerQC {
                            double missingRate, altFreq, mac, hweP;
                            bool pass;
                        };

                        FusedMarkerQC wmQC[64];  // B ≤ 64

                        for (const auto &group : fusedGroups) {
                            // ── Compute QC once using the representative phenotype ──
                            const auto &repFp = fusedPhenos[group.repFi];
                            const double *maskCol = AugResid.col(repFp.maskCol).data();

                            passAltFreqs.clear();
                            passGSums.clear();
                            int passCount = 0;

                            for (size_t bi = 0; bi < wlen; ++bi) {
                                const double *unionCol = GBatch_union.col(
                                    static_cast<Eigen::Index>(bi)).data();

                                PhenoGenoStats gs = statsFromUnionVec(
                                    unionCol, maskCol, nUnion, repFp.nUsed);

                                const double gSum = allScoresAndSums(
                                    repFp.maskCol, static_cast<Eigen::Index>(bi));
                                const double twoN = 2.0 * static_cast<double>(repFp.nUsed);
                                if (gs.altFreq == 0.0 && gs.hweP == 1.0) {
                                    gs.altFreq = gSum / twoN;
                                    double maf = std::min(gs.altFreq, 1.0 - gs.altFreq);
                                    gs.mac = maf * twoN;
                                }

                                wmQC[bi].missingRate = gs.missingRate;
                                wmQC[bi].altFreq = gs.altFreq;
                                wmQC[bi].mac = gs.mac;
                                wmQC[bi].hweP = gs.hweP;

                                double maf = std::min(gs.altFreq, 1.0 - gs.altFreq);
                                wmQC[bi].pass = !(maf < minMafCutoff ||
                                                  gs.mac < minMacCutoff);

                                if (wmQC[bi].pass) {
                                    passAltFreqs.push_back(gs.altFreq);
                                    passGSums.push_back(gSum);
                                    ++passCount;
                                }
                            }

                            // ── Per member: compact scores, run SPA, format output ──
                            for (size_t fi : group.fiIndices) {
                                const auto &fp = fusedPhenos[fi];
                                const size_t p = fp.taskIdx;

                                if (passCount > 0) {
                                    // Compact pass-QC score columns into pre-allocated buffer.
                                    int pi = 0;
                                    for (size_t bi = 0; bi < wlen; ++bi) {
                                        if (wmQC[bi].pass) {
                                            passScoresBuf.col(pi).head(fp.nCols) =
                                                allScoresAndSums.block(
                                                    fp.colOffset, static_cast<Eigen::Index>(bi),
                                                    fp.nCols, 1);
                                            ++pi;
                                        }
                                    }

                                    methods[p]->processScoreBatch(
                                        passScoresBuf.topLeftCorner(fp.nCols, passCount),
                                        passGSums.data(), fp.nUsed,
                                        passAltFreqs, fusedResultsBuf);
                                }

                                // Format output lines.
                                int ri = 0;
                                for (size_t bi = 0; bi < wlen; ++bi) {
                                    const auto &wm = winMarkers[bi];
                                    if (wmQC[bi].pass) {
                                        formatLine(phenoOut[p], fmtBuf, wm.chr, wm.pos,
                                                   wm.marker, wm.alt, wm.ref,
                                                   wmQC[bi].missingRate,
                                                   wmQC[bi].altFreq,
                                                   wmQC[bi].mac,
                                                   wmQC[bi].hweP,
                                                   fusedResultsBuf[ri++]);
                                    } else {
                                        formatLineNA(phenoOut[p], fmtBuf, wm.chr, wm.pos,
                                                     wm.marker, wm.alt, wm.ref,
                                                     wmQC[bi].missingRate,
                                                     wmQC[bi].altFreq,
                                                     wmQC[bi].mac,
                                                     wmQC[bi].hweP,
                                                     naSuffixes[p]);
                                    }
                                }
                            }
                        }
                    }

                    // Phase 4: Non-fuseable phenotypes — MissBatch extraction path.
                    if (!nonFusedPhenos.empty()) {
                        for (size_t mb = 0; mb < nMissBatches; ++mb) {
                            const auto &batch = missBatches[mb];
                            // Find the non-fused index for the representative.
                            size_t repNfIdx = 0;
                            for (size_t ni = 0; ni < nonFusedPhenos.size(); ++ni)
                                if (nonFusedPhenos[ni] == batch.rep) { repNfIdx = ni; break; }

                            const uint32_t nP = nfNPhenoArr[repNfIdx];
                            const bool ident = batch.identity;
                            const uint32_t *pidx = ident ? nullptr : batch.presentIndices.data();
                            const uint64_t *pmask = ident ? nullptr : batch.presentMask.data();

                            nfPassAltFreqs.clear();
                            int passCount = 0;

                            for (size_t bi = 0; bi < wlen; ++bi) {
                                const double *unionCol = GBatch_union.col(
                                    static_cast<Eigen::Index>(bi)).data();
                                double *phenoG = nfPhenoGPtrs[repNfIdx];

                                extractPhenoFast(unionCol, nUnion, nP, ident,
                                                 pidx, pmask, phenoG);

                                nfBatchMissings[repNfIdx].clear();
                                PhenoGenoStats gs = statsFromGVec(phenoG, nP,
                                                                  nfBatchMissings[repNfIdx]);
                                nfWinStats[bi] = gs;

                                double pMaf = std::min(gs.altFreq, 1.0 - gs.altFreq);
                                const bool pass = !(gs.missingRate > missingCutoff ||
                                                    pMaf < minMafCutoff ||
                                                    gs.mac < minMacCutoff ||
                                                    (hweCutoff > 0 && gs.hweP < hweCutoff));
                                nfWinPassQC[bi] = pass ? 1 : 0;

                                if (pass) {
                                    const double imputeG = 2.0 * gs.altFreq;
                                    for (uint32_t j : nfBatchMissings[repNfIdx])
                                        phenoG[j] = imputeG;
                                    GBatch_pheno_nf.col(passCount).head(nP) =
                                        Eigen::Map<const Eigen::VectorXd>(phenoG, nP);
                                    nfPassAltFreqs.push_back(gs.altFreq);
                                    ++passCount;
                                }
                            }

                            for (size_t mi = 0; mi < batch.members.size(); ++mi) {
                                const size_t p = batch.members[mi];

                                if (passCount > 0)
                                    methods[p]->getResultBatch(
                                        GBatch_pheno_nf.topLeftCorner(nP, passCount),
                                        nfPassAltFreqs, nfBatchResultsBuf);

                                int ri = 0;
                                for (size_t bi = 0; bi < wlen; ++bi) {
                                    const auto &wm = winMarkers[bi];
                                    if (nfWinPassQC[bi]) {
                                        formatLine(phenoOut[p], fmtBuf, wm.chr, wm.pos,
                                                   wm.marker, wm.alt, wm.ref,
                                                   nfWinStats[bi].missingRate,
                                                   nfWinStats[bi].altFreq,
                                                   nfWinStats[bi].mac,
                                                   nfWinStats[bi].hweP,
                                                   nfBatchResultsBuf[ri++]);
                                    } else {
                                        formatLineNA(phenoOut[p], fmtBuf, wm.chr, wm.pos,
                                                     wm.marker, wm.alt, wm.ref,
                                                     nfWinStats[bi].missingRate,
                                                     nfWinStats[bi].altFreq,
                                                     nfWinStats[bi].mac,
                                                     nfWinStats[bi].hweP,
                                                     naSuffixes[p]);
                                    }
                                }
                            }
                        }
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
