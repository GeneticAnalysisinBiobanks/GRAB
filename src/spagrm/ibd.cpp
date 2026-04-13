// ibd.cpp — Pairwise IBD estimation (C++17 / Eigen)
//
// Translates the R function getPairwiseIBD() into a single-pass streaming
// algorithm over the PLINK .bed file.
//
// Algorithm outline:
//   1. Load sparse GRM via SparseGRM; extract off-diagonal (related) pairs.
//   2. Create PlinkData for ALL subjects so that per-marker MAF is
//      computed from the full cohort (matching the R code's use of .frq).
//   3. Stream through every marker once.  For markers with MAF ≥ threshold,
//      mean-impute missing genotypes and accumulate weighted IBS0 statistics
//      for each related pair.
//   4. Derive (pa, pb, pc) from the accumulated statistics + GRM value.
//   5. Write output in the same pair order as the sparse GRM.

#include "spagrm/ibd.hpp"

#include "geno_factory/geno_data.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "io/subject_set.hpp"
#include "util/logging.hpp"
#include "util/text_stream.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <mutex>
#include <numeric>
#include <string>
#include <thread>
#include <vector>

#include <Eigen/Dense>
#include <immintrin.h>
#include "util/simd_dispatch.hpp"

// ── IBD pair accumulation — multi-versioned SIMD kernels ─────────────

__attribute__((target("avx2,avx512f,avx512vl,fma")))
static void ibdAccPairs_avx512(
    const double *genoPtr,
    const int32_t *i1,
    const int32_t *i2,
    double *localXW,
    double *localW,
    double inv_pv,
    double w,
    size_t nPairs
) {
    const __m512d signmask = _mm512_set1_pd(-0.0);
    const __m512d vones = _mm512_set1_pd(1.0);
    const __m512d vtwos = _mm512_set1_pd(2.0);
    const __m512d vinv_pv = _mm512_set1_pd(inv_pv);
    const __m512d vw = _mm512_set1_pd(w);

    const size_t nPairs8 = nPairs & ~size_t(7);
    size_t k = 0;

    for (; k < nPairs8; k += 8) {
        __m256i vi1 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(i1 + k));
        __m256i vi2 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(i2 + k));
        __m512d g1 = _mm512_i32gather_pd(vi1, genoPtr, 8);
        __m512d g2 = _mm512_i32gather_pd(vi2, genoPtr, 8);

        __m512d d = _mm512_sub_pd(g1, g2);
        // |d - 1| + |d + 1| - 2
        __m512d dm1 = _mm512_sub_pd(d, vones);
        __m512d dp1 = _mm512_add_pd(d, vones);
        // abs via andnot(signmask, x) — clear sign bit
        __m512d abs_dm1 = (__m512d)_mm512_andnot_si512((__m512i)signmask, (__m512i)dm1);
        __m512d abs_dp1 = (__m512d)_mm512_andnot_si512((__m512i)signmask, (__m512i)dp1);
        __m512d x = _mm512_mul_pd(_mm512_sub_pd(_mm512_add_pd(abs_dm1, abs_dp1), vtwos), vinv_pv);

        _mm512_storeu_pd(localXW + k, _mm512_fmadd_pd(x, vw, _mm512_loadu_pd(localXW + k)));
        _mm512_storeu_pd(localW + k, _mm512_add_pd(_mm512_loadu_pd(localW + k), vw));
    }

    // Scalar tail
    for (; k < nPairs; ++k) {
        const double d = genoPtr[i1[k]] - genoPtr[i2[k]];
        const double x = (std::abs(d - 1.0) + std::abs(d + 1.0) - 2.0) * inv_pv;
        localXW[k] += x * w;
        localW[k] += w;
    }
}

__attribute__((target("avx2,fma")))
static void ibdAccPairs_avx2(
    const double *genoPtr,
    const int32_t *i1,
    const int32_t *i2,
    double *localXW,
    double *localW,
    double inv_pv,
    double w,
    size_t nPairs
) {
    const __m256d signmask = _mm256_set1_pd(-0.0);
    const __m256d vones = _mm256_set1_pd(1.0);
    const __m256d vtwos = _mm256_set1_pd(2.0);
    const __m256d vinv_pv = _mm256_set1_pd(inv_pv);
    const __m256d vw = _mm256_set1_pd(w);

    const size_t nPairs4 = nPairs & ~size_t(3);
    size_t k = 0;

    for (; k < nPairs4; k += 4) {
        __m128i vi1 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(i1 + k));
        __m128i vi2 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(i2 + k));
        __m256d g1 = _mm256_i32gather_pd(genoPtr, vi1, 8);
        __m256d g2 = _mm256_i32gather_pd(genoPtr, vi2, 8);

        __m256d d = _mm256_sub_pd(g1, g2);
        __m256d abs_dm1 = _mm256_andnot_pd(signmask, _mm256_sub_pd(d, vones));
        __m256d abs_dp1 = _mm256_andnot_pd(signmask, _mm256_add_pd(d, vones));
        __m256d x = _mm256_mul_pd(_mm256_sub_pd(_mm256_add_pd(abs_dm1, abs_dp1), vtwos), vinv_pv);

        _mm256_storeu_pd(localXW + k, _mm256_fmadd_pd(x, vw, _mm256_loadu_pd(localXW + k)));
        _mm256_storeu_pd(localW + k, _mm256_add_pd(_mm256_loadu_pd(localW + k), vw));
    }

    // Scalar tail
    for (; k < nPairs; ++k) {
        const double d = genoPtr[i1[k]] - genoPtr[i2[k]];
        const double x = (std::abs(d - 1.0) + std::abs(d + 1.0) - 2.0) * inv_pv;
        localXW[k] += x * w;
        localW[k] += w;
    }
}

static void ibdAccPairs_scalar(
    const double *genoPtr,
    const int32_t *i1,
    const int32_t *i2,
    double *localXW,
    double *localW,
    double inv_pv,
    double w,
    size_t nPairs
) {
    for (size_t k = 0; k < nPairs; ++k) {
        const double d = genoPtr[i1[k]] - genoPtr[i2[k]];
        const double x = (std::abs(d - 1.0) + std::abs(d + 1.0) - 2.0) * inv_pv;
        localXW[k] += x * w;
        localW[k] += w;
    }
}

using IbdAccFn = void (*)(const double *, const int32_t *, const int32_t *,
                          double *, double *, double, double, size_t);

static IbdAccFn pickIbdAccFn() {
    switch (simdLevel()) {
    case SimdLevel::AVX512: return ibdAccPairs_avx512;
    case SimdLevel::AVX2:   return ibdAccPairs_avx2;
    default:                return ibdAccPairs_scalar;
    }
}

static const IbdAccFn ibdAccPairs = pickIbdAccFn();

// ══════════════════════════════════════════════════════════════════════════════
// runPairwiseIBD
// ══════════════════════════════════════════════════════════════════════════════

void runPairwiseIBD(
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outputFile,
    const std::string &keepFile,
    const std::string &removeFile,
    double minMafIBD,
    int nthreads
) {
    // ── 1. Read sample IDs and apply keep/remove filters ────────────
    std::vector<std::string> allIIDs = parseGenoIIDs(geno);
    infoMsg("Read %u subjects from genotype file", static_cast<uint32_t>(allIIDs.size()));

    // Build filtered subject set via SubjectSet
    SubjectSet ss(std::move(allIIDs));
    ss.setKeepRemove(keepFile, removeFile);
    ss.setGenoLabel(geno.flagLabel());
    ss.finalize();

    const uint32_t nFam = ss.nFam();
    const uint32_t nFiltered = ss.nUsed();
    const auto &usedMask = ss.usedMask();
    auto filteredIIDs = ss.usedIIDs();

    // ── 2. Load sparse GRM against filtered subjects ─────────────────
    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, filteredIIDs, filteredIIDs);
    infoMsg("Loaded sparse GRM: %zu entries", grm.nnz());

    // ── 3. Extract off-diagonal pairs (sorted by row, col — same as GRM) ─
    struct IndexedPair {
        uint32_t idx1, idx2;
        double grmValue;
    };

    std::vector<IndexedPair> pairs;
    for (const auto &e : grm.entries()) {
        if (e.row != e.col) pairs.push_back({e.row, e.col, e.value});
    }
    infoMsg("Found %zu off-diagonal (related) pairs", pairs.size());

    if (pairs.empty()) {
        TextWriter writer(outputFile);
        writer.write("#ID1\tID2\tpa\tpb\tpc\n");
        infoMsg("No related pairs found — wrote empty output to %s", outputFile.c_str());
        return;
    }

    // ── 4. Create genotype data for filtered subjects ─────────────────
    auto genoData = makeGenoData(geno, usedMask, nFam, nFiltered, 1024);

    const uint32_t nSubj = genoData->nSubjUsed();
    const uint32_t nMarkers = genoData->nMarkers();

    if (nMarkers > 200000) {
        infoMsg("WARNING: %u markers provided. "
                "IBD estimation is intended for LD-pruned markers (~100k). "
                "Consider pruning for speed and accuracy.",
                nMarkers);
    }

    infoMsg("IBD analysis: %u subjects, %u markers (MAF threshold %.4f)", nSubj, nMarkers, minMafIBD);

    // ── 5. SoA pair indices for AVX2 gather ──────────────────────────
    const size_t nPairs = pairs.size();
    std::vector<int32_t> pairIdx1(nPairs), pairIdx2(nPairs);
    for (size_t k = 0; k < nPairs; ++k) {
        pairIdx1[k] = static_cast<int32_t>(pairs[k].idx1);
        pairIdx2[k] = static_cast<int32_t>(pairs[k].idx2);
    }

    // ── 6. Stream markers (multi-threaded, AVX2) ──────────────────────
    const auto &chunks = genoData->chunkIndices();
    const size_t nChunks = chunks.size();
    if (nthreads < 1) nthreads = 1;
    const int nt = std::min(nthreads, static_cast<int>(nChunks));

    struct ThreadAcc {
        std::vector<double> sumXW, sumW;
        uint32_t nUsed = 0;
    };

    std::vector<ThreadAcc> threadAccs(nt);
    for (auto &ta : threadAccs) {
        ta.sumXW.assign(nPairs, 0.0);
        ta.sumW.assign(nPairs, 0.0);
    }

    std::atomic<size_t> nextChunk(0);
    std::atomic<size_t> chunksDone(0);

    auto workerFn = [&](int tid) {
        auto cursor = genoData->makeCursor();
        Eigen::VectorXd genoVec(nSubj);
        std::vector<uint32_t> missingIdx;
        double altFreq, altCounts, missingRate, hweP, maf, mac;

        double *localXW = threadAccs[tid].sumXW.data();
        double *localW = threadAccs[tid].sumW.data();
        uint32_t localUsed = 0;

        const int32_t *i1 = pairIdx1.data();
        const int32_t *i2 = pairIdx2.data();

        while (true) {
            const size_t ci = nextChunk.fetch_add(1);
            if (ci >= nChunks) break;

            const auto &chunk = chunks[ci];
            if (!chunk.empty()) cursor->beginSequentialBlock(chunk.front());

            for (uint64_t gi : chunk) {
                cursor->getGenotypes(gi, genoVec, altFreq, altCounts, missingRate, hweP, maf, mac, missingIdx);
                if (maf < minMafIBD) continue;

                const double meanGeno = 2.0 * altFreq;
                for (uint32_t mi : missingIdx)
                    genoVec[mi] = meanGeno;

                const double pq = altFreq * (1.0 - altFreq);
                const double pv = 2.0 * pq * pq;
                const double opv = 1.0 - pv;
                if (pv < 1e-30 || opv <= 0.0) continue;

                const double inv_pv = 1.0 / pv;
                const double w = std::sqrt(pv / opv);
                const double *genoPtr = genoVec.data();

                ibdAccPairs(genoPtr, i1, i2, localXW, localW, inv_pv, w, nPairs);

                ++localUsed;
            }

            const size_t done = chunksDone.fetch_add(1) + 1;
            if (nChunks >= 10 && done % (nChunks / 10) == 0) infoMsg("  IBD progress: %zu/%zu chunks", done, nChunks);
        }
        threadAccs[tid].nUsed = localUsed;
    };

    infoMsg("IBD: streaming %u markers (%zu chunks) with %d threads", nMarkers, nChunks, nt);

    if (nt == 1) {
        workerFn(0);
    } else {
        std::vector<std::thread> workers;
        workers.reserve(nt - 1);
        for (int t = 1; t < nt; ++t)
            workers.emplace_back(workerFn, t);
        workerFn(0);
        for (auto &w : workers)
            w.join();
    }

    // Reduce thread-local accumulators
    uint32_t nUsed = 0;
    for (int t = 0; t < nt; ++t)
        nUsed += threadAccs[t].nUsed;
    for (int t = 1; t < nt; ++t) {
        const double *txw = threadAccs[t].sumXW.data();
        const double *tw = threadAccs[t].sumW.data();
        double *rxw = threadAccs[0].sumXW.data();
        double *rw = threadAccs[0].sumW.data();
        for (size_t k = 0; k < nPairs; ++k) {
            rxw[k] += txw[k];
            rw[k] += tw[k];
        }
    }

    infoMsg("IBD: used %u / %u markers (MAF >= %.4f)", nUsed, nMarkers, minMafIBD);

    // ── 7. Derive (pa, pb, pc) for each pair ───────────────────────────
    const double *sumXW = threadAccs[0].sumXW.data();
    const double *sumW = threadAccs[0].sumW.data();
    //
    // pc_raw = 0.5 * sumXW / sumW
    //
    // Constraints (from R code):
    //   upper_bound = (1 - grmValue)^2 - 1e-10
    //   lower_bound = 1 - 2 * grmValue
    //   pc = clamp(pc_raw, lower_bound, upper_bound)
    //   pb = 2 - 2*pc - 2*grmValue
    //   pa = 2*grmValue + pc - 1
    //   if (pb < 0) { pa += 0.5*pb; pb = 0; pc = 0; }

    struct IBDResult {
        uint32_t idx1, idx2;
        double pa, pb, pc;
    };

    std::vector<IBDResult> results(nPairs);

    for (size_t k = 0; k < nPairs; ++k) {
        double pc = (sumW[k] > 0.0) ? 0.5 * sumXW[k] / sumW[k] : 0.0;
        const double grm = pairs[k].grmValue;

        const double upper = (1.0 - grm) * (1.0 - grm) - 1e-10;
        const double lower = 1.0 - 2.0 * grm;
        pc = std::clamp(pc, lower, upper);

        double pb = 2.0 - 2.0 * pc - 2.0 * grm;
        double pa = 2.0 * grm + pc - 1.0;

        if (pb < 0.0) {
            pa += 0.5 * pb;
            pb = 0.0;
            pc = 0.0;
        }

        results[k] = {pairs[k].idx1, pairs[k].idx2, pa, pb, pc};
    }

    // ── 8. Write output (same pair order as the sparse GRM) ────────────
    TextWriter writer(outputFile);
    writer.write("#ID1\tID2\tpa\tpb\tpc\n");
    char buf[128];
    for (const auto &r : results) {
        std::string line;
        line.reserve(128);
        line += filteredIIDs[r.idx1];
        line += '\t';
        line += filteredIIDs[r.idx2];
        line += '\t';
        int n = std::snprintf(buf, sizeof(buf), "%.17g\t%.17g\t%.17g\n", r.pa, r.pb, r.pc);
        line.append(buf, n);
        writer.write(line);
    }

    infoMsg("Wrote %zu IBD records to %s", results.size(), outputFile.c_str());
}
