// spamixlocalp.cpp — SPAmixLocalPlus implementation
//
// Phi estimation + per-ancestry GWAS with streaming admix binary I/O.

#include "localplus/spamixlocalp.hpp"
#include "engine/marker.hpp"
#include "localplus/abed_io.hpp"
#include "spamix/common.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "io/subject_set.hpp"
#include "util/logging.hpp"
#include "util/null_model.hpp"
#include "util/text_scanner.hpp"
#include "util/text_stream.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/distributions/normal.hpp>

static constexpr double MIN_P_VALUE = std::numeric_limits<double>::min();

// ======================================================================
// Precomputed phi products — SoA layout for AVX2-friendly hot loop
// ======================================================================

RprodSoA buildRprodSoA(
    const PhiMatrices &phi,
    const Eigen::VectorXd &R
) {
    size_t total = phi.A.size() + phi.B.size() + phi.C.size() + phi.D.size();
    RprodSoA rp;
    rp.reserve(total);

    auto append = [&](const std::vector<PhiEntry> &src, double mult, uint8_t th_i, uint8_t th_j) {
        for (const auto &e : src)
            rp.push_back(e.i, e.j, mult * e.value * R[e.i] * R[e.j], th_i, th_j);
    };
    append(phi.A, 4.0, 2, 2);
    append(phi.B, 2.0, 2, 1);
    append(phi.C, 2.0, 1, 2);
    append(phi.D, 1.0, 1, 1);
    return rp;
}

// Build the multi-phenotype baked rprod table for one ancestry.
//
// rprod_packed[e * K_pheno + p] = mult · phi(i,j) · R[i,p] · R[j,p]
//
// Memory is allocated once and shared (const) across all worker threads.
// Layout is row-major over entries (K_pheno consecutive doubles per entry)
// so the hot scan loop reads K_pheno values sequentially — essential for
// SIMD throughput.
MultiPhenoRprodSoA buildMultiPhenoRprodSoA(
    const PhiMatrices &phi,
    const Eigen::Ref<const Eigen::MatrixXd> &R_mat
) {
    const size_t total = phi.A.size() + phi.B.size() + phi.C.size() + phi.D.size();
    const int K_pheno = static_cast<int>(R_mat.cols());

    MultiPhenoRprodSoA rp;
    rp.K_pheno = K_pheno;
    rp.idx_i.reserve(total);
    rp.idx_j.reserve(total);
    rp.target_hi.reserve(total);
    rp.target_hj.reserve(total);
    rp.rprod_packed.resize(total * static_cast<size_t>(K_pheno));

    size_t cursor = 0;
    auto append = [&](const std::vector<PhiEntry> &src,
                      double mult,
                      uint8_t th_i,
                      uint8_t th_j) {
        for (const auto &e : src) {
            rp.idx_i.push_back(e.i);
            rp.idx_j.push_back(e.j);
            rp.target_hi.push_back(th_i);
            rp.target_hj.push_back(th_j);
            double base = mult * e.value;
            for (int p = 0; p < K_pheno; ++p) {
                rp.rprod_packed[cursor * K_pheno + p] =
                    base * R_mat(e.i, p) * R_mat(e.j, p);
            }
            ++cursor;
        }
    };
    append(phi.A, 4.0, 2, 2);
    append(phi.B, 2.0, 2, 1);
    append(phi.C, 2.0, 1, 2);
    append(phi.D, 1.0, 1, 1);
    return rp;
}

// ── Runtime SIMD dispatch — phi variance functions ──────────────────
//
// Three tiers: AVX-512 (8 doubles / 512 bits), AVX2 (4 doubles / 256 bits),
// and scalar.  The resolver runs once at process startup via
// __builtin_cpu_supports() and caches the result; callers pay only a
// predictable branch per invocation.

#if defined(__x86_64__) || defined(_M_X64)
#  include <immintrin.h>
#endif
#include "util/simd_dispatch.hpp"

// ════════════════════════════════════════════════════════════════════════
// computeVarOffSoA — single-marker off-diagonal variance
// ════════════════════════════════════════════════════════════════════════

#if defined(__x86_64__) || defined(_M_X64)

// ── AVX-512: 8 phi entries per iteration ────────────────────────────

__attribute__((target("avx2,avx512f,avx512vl,fma")))
static double computeVarOffSoA_avx512(
    const RprodSoA &rp,
    const uint32_t *hInt
) {
    const size_t N = rp.size();
    if (N == 0) return 0.0;

    const uint32_t *__restrict ii = rp.idx_i.data();
    const uint32_t *__restrict jj = rp.idx_j.data();
    const double *__restrict rval = rp.rprod.data();
    const uint8_t *__restrict thi = rp.target_hi.data();
    const uint8_t *__restrict thj = rp.target_hj.data();

    __m512d acc = _mm512_setzero_pd();
    size_t p = 0;

    for (; p + 7 < N; p += 8) {
        // Load 8 index pairs (8×uint32 = 256 bits)
        __m256i vi = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(ii + p));
        __m256i vj = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(jj + p));

        // Gather hInt[i] and hInt[j] — 8×int32
        __m256i hi = _mm256_i32gather_epi32(reinterpret_cast<const int *>(hInt), vi, 4);
        __m256i hj = _mm256_i32gather_epi32(reinterpret_cast<const int *>(hInt), vj, 4);

        // Load 8 target bytes → 8×int32
        __m128i raw_thi = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(thi + p));
        __m128i raw_thj = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(thj + p));
        __m256i vthi = _mm256_cvtepu8_epi32(raw_thi);
        __m256i vthj = _mm256_cvtepu8_epi32(raw_thj);

        // Compare → __mmask8 (AVX-512VL)
        __mmask8 mi = _mm256_cmpeq_epi32_mask(hi, vthi);
        __mmask8 mj = _mm256_cmpeq_epi32_mask(hj, vthj);

        // Masked accumulate 8 rprod doubles
        __m512d rv = _mm512_loadu_pd(rval + p);
        acc = _mm512_mask_add_pd(acc, static_cast<__mmask8>(mi & mj), acc, rv);
    }

    double varOff = _mm512_reduce_add_pd(acc);

    // Scalar tail (0–7 remaining entries)
    for (; p < N; ++p)
        varOff += ((hInt[ii[p]] == thi[p]) & (hInt[jj[p]] == thj[p])) * rval[p];

    return varOff;
}

// ── AVX2: 4 phi entries per iteration ───────────────────────────────

__attribute__((target("avx2,fma")))
static double computeVarOffSoA_avx2(
    const RprodSoA &rp,
    const uint32_t *hInt
) {
    const size_t N = rp.size();
    if (N == 0) return 0.0;

    const uint32_t *__restrict ii = rp.idx_i.data();
    const uint32_t *__restrict jj = rp.idx_j.data();
    const double *__restrict rval = rp.rprod.data();
    const uint8_t *__restrict thi = rp.target_hi.data();
    const uint8_t *__restrict thj = rp.target_hj.data();

    __m256d acc = _mm256_setzero_pd();
    size_t p = 0;

    for (; p + 3 < N; p += 4) {
        __m128i vi = _mm_loadu_si128(reinterpret_cast<const __m128i *>(ii + p));
        __m128i vj = _mm_loadu_si128(reinterpret_cast<const __m128i *>(jj + p));

        __m128i hi = _mm_i32gather_epi32(reinterpret_cast<const int *>(hInt), vi, 4);
        __m128i hj = _mm_i32gather_epi32(reinterpret_cast<const int *>(hInt), vj, 4);

        __m128i raw_thi = _mm_cvtsi32_si128(*reinterpret_cast<const int *>(thi + p));
        __m128i raw_thj = _mm_cvtsi32_si128(*reinterpret_cast<const int *>(thj + p));
        __m128i vthi = _mm_cvtepu8_epi32(raw_thi);
        __m128i vthj = _mm_cvtepu8_epi32(raw_thj);

        __m128i cmpi = _mm_cmpeq_epi32(hi, vthi);
        __m128i cmpj = _mm_cmpeq_epi32(hj, vthj);
        __m128i mask32 = _mm_and_si128(cmpi, cmpj);

        __m256i mask64 = _mm256_cvtepi32_epi64(mask32);
        __m256d maskd = _mm256_castsi256_pd(mask64);

        __m256d rv = _mm256_loadu_pd(rval + p);
        acc = _mm256_add_pd(acc, _mm256_and_pd(maskd, rv));
    }

    __m128d lo = _mm256_castpd256_pd128(acc);
    __m128d hi128 = _mm256_extractf128_pd(acc, 1);
    __m128d sum2 = _mm_add_pd(lo, hi128);
    double tmp[2];
    _mm_storeu_pd(tmp, sum2);
    double varOff = tmp[0] + tmp[1];

    for (; p < N; ++p)
        varOff += ((hInt[ii[p]] == thi[p]) & (hInt[jj[p]] == thj[p])) * rval[p];

    return varOff;
}

#endif  // x86_64 SIMD variants — single-marker

// ── Scalar baseline ─────────────────────────────────────────────────

static double computeVarOffSoA_scalar(
    const RprodSoA &rp,
    const uint32_t *hInt
) {
    const size_t N = rp.size();
    double varOff = 0.0;
    for (size_t p = 0; p < N; ++p)
        varOff += ((hInt[rp.idx_i[p]] == rp.target_hi[p]) & (hInt[rp.idx_j[p]] == rp.target_hj[p])) * rp.rprod[p];
    return varOff;
}

// ── Public dispatch ─────────────────────────────────────────────────

double computeVarOffSoA(
    const RprodSoA &rp,
    const uint32_t *hInt,
    uint32_t                                                               /*nUsed*/
) {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return computeVarOffSoA_avx512(rp, hInt);
    case SimdLevel::AVX2:   return computeVarOffSoA_avx2(rp, hInt);
    default: break;
    }
#endif
    return computeVarOffSoA_scalar(rp, hInt);
}

// ════════════════════════════════════════════════════════════════════════
// computeVarOffSoABatch — scan phi entries ONCE for PHI_BATCH markers.
// Subject-major hIntSM layout: hIntSM[subj * PHI_BATCH + b].
// Each phi entry reads PHI_BATCH consecutive uint32s (1–2 cache lines),
// amortising DRAM bandwidth across all batch markers.
// ════════════════════════════════════════════════════════════════════════

#if defined(__x86_64__) || defined(_M_X64)

// ── AVX-512: one __m512d accumulator covers all 8 batch items ───────

__attribute__((target("avx2,avx512f,avx512vl,fma")))
static void computeVarOffSoABatch_avx512(
    const RprodSoA &rp,
    const uint32_t *hIntSM,
    int batchLen,
    double *varOff
) {
    const size_t N = rp.size();
    std::fill(varOff, varOff + batchLen, 0.0);
    if (N == 0) return;

    const uint32_t *__restrict ii = rp.idx_i.data();
    const uint32_t *__restrict jj = rp.idx_j.data();
    const double *__restrict rval = rp.rprod.data();
    const uint8_t *__restrict thi = rp.target_hi.data();
    const uint8_t *__restrict thj = rp.target_hj.data();

    // Single 512-bit accumulator for all 8 batch items
    __m512d acc = _mm512_setzero_pd();

    for (size_t p = 0; p < N; ++p) {
        const uint32_t *hi_ptr = hIntSM + ii[p] * PHI_BATCH;
        const uint32_t *hj_ptr = hIntSM + jj[p] * PHI_BATCH;

        // Load 8 hInt values per subject (256 bits)
        __m256i hi8 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(hi_ptr));
        __m256i hj8 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(hj_ptr));

        // Broadcast targets
        __m256i vthi = _mm256_set1_epi32(static_cast<int>(thi[p]));
        __m256i vthj = _mm256_set1_epi32(static_cast<int>(thj[p]));

        // Compare → __mmask8 (AVX-512VL)
        __mmask8 mi = _mm256_cmpeq_epi32_mask(hi8, vthi);
        __mmask8 mj = _mm256_cmpeq_epi32_mask(hj8, vthj);

        // Broadcast rprod and masked-accumulate
        __m512d rv = _mm512_set1_pd(rval[p]);
        acc = _mm512_mask_add_pd(acc, static_cast<__mmask8>(mi & mj), acc, rv);
    }

    double tmp[8];
    _mm512_storeu_pd(tmp, acc);
    for (int b = 0; b < batchLen; ++b)
        varOff[b] = tmp[b];
}

// ── AVX2: two __m256d accumulators for batch[0..3] and [4..7] ───────

__attribute__((target("avx2,fma")))
static void computeVarOffSoABatch_avx2(
    const RprodSoA &rp,
    const uint32_t *hIntSM,
    int batchLen,
    double *varOff
) {
    const size_t N = rp.size();
    std::fill(varOff, varOff + batchLen, 0.0);
    if (N == 0) return;

    const uint32_t *__restrict ii = rp.idx_i.data();
    const uint32_t *__restrict jj = rp.idx_j.data();
    const double *__restrict rval = rp.rprod.data();
    const uint8_t *__restrict thi = rp.target_hi.data();
    const uint8_t *__restrict thj = rp.target_hj.data();

    __m256d acc0 = _mm256_setzero_pd();
    __m256d acc1 = _mm256_setzero_pd();

    for (size_t p = 0; p < N; ++p) {
        const uint32_t *hi_ptr = hIntSM + ii[p] * PHI_BATCH;
        const uint32_t *hj_ptr = hIntSM + jj[p] * PHI_BATCH;
        __m128i vthi = _mm_set1_epi32(static_cast<int>(thi[p]));
        __m128i vthj = _mm_set1_epi32(static_cast<int>(thj[p]));
        __m256d rv = _mm256_set1_pd(rval[p]);

        // Batch items 0–3
        {
            __m128i hi4 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(hi_ptr));
            __m128i hj4 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(hj_ptr));
            __m128i mask32 = _mm_and_si128(_mm_cmpeq_epi32(hi4, vthi), _mm_cmpeq_epi32(hj4, vthj));
            __m256d maskd = _mm256_castsi256_pd(_mm256_cvtepi32_epi64(mask32));
            acc0 = _mm256_add_pd(acc0, _mm256_and_pd(maskd, rv));
        }
        // Batch items 4–7
        {
            __m128i hi4 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(hi_ptr + 4));
            __m128i hj4 = _mm_loadu_si128(reinterpret_cast<const __m128i *>(hj_ptr + 4));
            __m128i mask32 = _mm_and_si128(_mm_cmpeq_epi32(hi4, vthi), _mm_cmpeq_epi32(hj4, vthj));
            __m256d maskd = _mm256_castsi256_pd(_mm256_cvtepi32_epi64(mask32));
            acc1 = _mm256_add_pd(acc1, _mm256_and_pd(maskd, rv));
        }
    }

    double tmp[8];
    _mm256_storeu_pd(tmp, acc0);
    _mm256_storeu_pd(tmp + 4, acc1);
    for (int b = 0; b < batchLen; ++b)
        varOff[b] = tmp[b];
}

#endif  // x86_64 SIMD variants — batch

// ── Scalar baseline ─────────────────────────────────────────────────

static void computeVarOffSoABatch_scalar(
    const RprodSoA &rp,
    const uint32_t *hIntSM,
    int batchLen,
    double *varOff
) {
    const size_t N = rp.size();
    std::fill(varOff, varOff + batchLen, 0.0);
    for (size_t p = 0; p < N; ++p) {
        uint32_t i = rp.idx_i[p], j = rp.idx_j[p];
        double r = rp.rprod[p];
        uint32_t th_i = rp.target_hi[p], th_j = rp.target_hj[p];
        const uint32_t *hi_ptr = hIntSM + i * PHI_BATCH;
        const uint32_t *hj_ptr = hIntSM + j * PHI_BATCH;
        for (int b = 0; b < batchLen; ++b)
            varOff[b] += ((hi_ptr[b] == th_i) & (hj_ptr[b] == th_j)) * r;
    }
}

// ── Public dispatch ─────────────────────────────────────────────────

void computeVarOffSoABatch(
    const RprodSoA &rp,
    const uint32_t *hIntSM,
    uint32_t /*nUsed*/,
    int batchLen,
    double *varOff
) {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: computeVarOffSoABatch_avx512(rp, hIntSM, batchLen, varOff); return;
    case SimdLevel::AVX2:   computeVarOffSoABatch_avx2(rp, hIntSM, batchLen, varOff); return;
    default: break;
    }
#endif
    computeVarOffSoABatch_scalar(rp, hIntSM, batchLen, varOff);
}

// ════════════════════════════════════════════════════════════════════════
// computeVarOffMultiPhenoBatch — multi-phenotype fused phi scan.
//
// One pass over phi entries serves all K_pheno phenotypes.  The match mask
// (hi == target_hi & hj == target_hj) depends only on (entry, batch lane);
// it is computed once and reused.  rprod_packed[e * K_pheno + p] is read
// sequentially (K_pheno consecutive doubles per entry), so the per-
// phenotype cost is one cache-resident scalar load + broadcast + masked-
// add.  This avoids the random R-matrix gather pattern that would result
// from computing R[i,p] * R[j,p] inline.
//
// Output layout: varOff[b * K_pheno + p].
// ════════════════════════════════════════════════════════════════════════

// Scalar fallback (also used as the K_pheno > MAX_KP slow path).
static void computeVarOffMultiPhenoBatch_scalar(
    const MultiPhenoRprodSoA &rp,
    const uint32_t *hIntSM,
    int batchLen,
    double *varOff
) {
    const size_t E = rp.nEntries();
    const int K_pheno = rp.K_pheno;
    std::fill(varOff, varOff + static_cast<size_t>(batchLen) * K_pheno, 0.0);
    if (E == 0) return;

    const double *__restrict pp = rp.rprod_packed.data();

    for (size_t e = 0; e < E; ++e) {
        uint32_t i = rp.idx_i[e];
        uint32_t j = rp.idx_j[e];
        uint32_t ti = rp.target_hi[e];
        uint32_t tj = rp.target_hj[e];
        const uint32_t *hi_ptr = hIntSM + static_cast<size_t>(i) * PHI_BATCH;
        const uint32_t *hj_ptr = hIntSM + static_cast<size_t>(j) * PHI_BATCH;
        const double *rp_e = pp + e * K_pheno;

        for (int p = 0; p < K_pheno; ++p) {
            double term = rp_e[p];
            for (int b = 0; b < batchLen; ++b) {
                int matched = (hi_ptr[b] == ti) & (hj_ptr[b] == tj);
                varOff[static_cast<size_t>(b) * K_pheno + p] += matched * term;
            }
        }
    }
}

#if defined(__x86_64__) || defined(_M_X64)

// AVX-512 path — one 512-bit accumulator per phenotype covers all 8 batch
// lanes.  Cap at 16 phenotypes for register-pressure safety.
static constexpr int MULTI_PHENO_AVX512_MAX_KP = 16;

__attribute__((target("avx2,avx512f,avx512vl,fma")))
static void computeVarOffMultiPhenoBatch_avx512(
    const MultiPhenoRprodSoA &rp,
    const uint32_t *hIntSM,
    int batchLen,
    double *varOff
) {
    const int K_pheno = rp.K_pheno;
    if (K_pheno > MULTI_PHENO_AVX512_MAX_KP) {
        computeVarOffMultiPhenoBatch_scalar(rp, hIntSM, batchLen, varOff);
        return;
    }
    const size_t E = rp.nEntries();
    std::fill(varOff, varOff + static_cast<size_t>(batchLen) * K_pheno, 0.0);
    if (E == 0) return;

    const uint32_t *__restrict ii  = rp.idx_i.data();
    const uint32_t *__restrict jj  = rp.idx_j.data();
    const uint8_t  *__restrict thi = rp.target_hi.data();
    const uint8_t  *__restrict thj = rp.target_hj.data();
    const double   *__restrict pp  = rp.rprod_packed.data();

    __m512d acc[MULTI_PHENO_AVX512_MAX_KP];
    for (int p = 0; p < K_pheno; ++p) acc[p] = _mm512_setzero_pd();

    for (size_t e = 0; e < E; ++e) {
        uint32_t i = ii[e];
        uint32_t j = jj[e];
        const uint32_t *hi_ptr = hIntSM + static_cast<size_t>(i) * PHI_BATCH;
        const uint32_t *hj_ptr = hIntSM + static_cast<size_t>(j) * PHI_BATCH;

        __m256i hi8 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(hi_ptr));
        __m256i hj8 = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(hj_ptr));
        __m256i vthi = _mm256_set1_epi32(static_cast<int>(thi[e]));
        __m256i vthj = _mm256_set1_epi32(static_cast<int>(thj[e]));
        __mmask8 mi = _mm256_cmpeq_epi32_mask(hi8, vthi);
        __mmask8 mj = _mm256_cmpeq_epi32_mask(hj8, vthj);
        __mmask8 mask = static_cast<__mmask8>(mi & mj);

        const double *rp_e = pp + e * K_pheno;
        for (int p = 0; p < K_pheno; ++p) {
            __m512d sv = _mm512_set1_pd(rp_e[p]);
            acc[p] = _mm512_mask_add_pd(acc[p], mask, acc[p], sv);
        }
    }

    double tmp[8];
    for (int p = 0; p < K_pheno; ++p) {
        _mm512_storeu_pd(tmp, acc[p]);
        for (int b = 0; b < batchLen; ++b)
            varOff[static_cast<size_t>(b) * K_pheno + p] = tmp[b];
    }
}

// AVX2 path: K_pheno × 2 accumulators of __m256d cover batch[0..3] and [4..7].
// Cap at 8 phenotypes for YMM register pressure (16 YMMs available).
static constexpr int MULTI_PHENO_AVX2_MAX_KP = 8;

__attribute__((target("avx2,fma")))
static void computeVarOffMultiPhenoBatch_avx2(
    const MultiPhenoRprodSoA &rp,
    const uint32_t *hIntSM,
    int batchLen,
    double *varOff
) {
    const int K_pheno = rp.K_pheno;
    if (K_pheno > MULTI_PHENO_AVX2_MAX_KP) {
        computeVarOffMultiPhenoBatch_scalar(rp, hIntSM, batchLen, varOff);
        return;
    }
    const size_t E = rp.nEntries();
    std::fill(varOff, varOff + static_cast<size_t>(batchLen) * K_pheno, 0.0);
    if (E == 0) return;

    const uint32_t *__restrict ii  = rp.idx_i.data();
    const uint32_t *__restrict jj  = rp.idx_j.data();
    const uint8_t  *__restrict thi = rp.target_hi.data();
    const uint8_t  *__restrict thj = rp.target_hj.data();
    const double   *__restrict pp  = rp.rprod_packed.data();

    __m256d acc0[MULTI_PHENO_AVX2_MAX_KP];
    __m256d acc1[MULTI_PHENO_AVX2_MAX_KP];
    for (int p = 0; p < K_pheno; ++p) {
        acc0[p] = _mm256_setzero_pd();
        acc1[p] = _mm256_setzero_pd();
    }

    for (size_t e = 0; e < E; ++e) {
        uint32_t i = ii[e];
        uint32_t j = jj[e];
        const uint32_t *hi_ptr = hIntSM + static_cast<size_t>(i) * PHI_BATCH;
        const uint32_t *hj_ptr = hIntSM + static_cast<size_t>(j) * PHI_BATCH;

        __m128i vthi = _mm_set1_epi32(static_cast<int>(thi[e]));
        __m128i vthj = _mm_set1_epi32(static_cast<int>(thj[e]));

        __m128i hi4a = _mm_loadu_si128(reinterpret_cast<const __m128i *>(hi_ptr));
        __m128i hj4a = _mm_loadu_si128(reinterpret_cast<const __m128i *>(hj_ptr));
        __m128i mask32a = _mm_and_si128(_mm_cmpeq_epi32(hi4a, vthi),
                                        _mm_cmpeq_epi32(hj4a, vthj));
        __m256d maskd_a = _mm256_castsi256_pd(_mm256_cvtepi32_epi64(mask32a));

        __m128i hi4b = _mm_loadu_si128(reinterpret_cast<const __m128i *>(hi_ptr + 4));
        __m128i hj4b = _mm_loadu_si128(reinterpret_cast<const __m128i *>(hj_ptr + 4));
        __m128i mask32b = _mm_and_si128(_mm_cmpeq_epi32(hi4b, vthi),
                                        _mm_cmpeq_epi32(hj4b, vthj));
        __m256d maskd_b = _mm256_castsi256_pd(_mm256_cvtepi32_epi64(mask32b));

        const double *rp_e = pp + e * K_pheno;
        for (int p = 0; p < K_pheno; ++p) {
            __m256d sv = _mm256_set1_pd(rp_e[p]);
            acc0[p] = _mm256_add_pd(acc0[p], _mm256_and_pd(maskd_a, sv));
            acc1[p] = _mm256_add_pd(acc1[p], _mm256_and_pd(maskd_b, sv));
        }
    }

    double tmp[8];
    for (int p = 0; p < K_pheno; ++p) {
        _mm256_storeu_pd(tmp, acc0[p]);
        _mm256_storeu_pd(tmp + 4, acc1[p]);
        for (int b = 0; b < batchLen; ++b)
            varOff[static_cast<size_t>(b) * K_pheno + p] = tmp[b];
    }
}

#endif // x86_64 SIMD variants — multi-pheno batch

void computeVarOffMultiPhenoBatch(
    const MultiPhenoRprodSoA &rp,
    const uint32_t *hIntSM,
    int batchLen,
    double *varOff
) {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512:
        computeVarOffMultiPhenoBatch_avx512(rp, hIntSM, batchLen, varOff);
        return;
    case SimdLevel::AVX2:
        computeVarOffMultiPhenoBatch_avx2(rp, hIntSM, batchLen, varOff);
        return;
    default: break;
    }
#endif
    computeVarOffMultiPhenoBatch_scalar(rp, hIntSM, batchLen, varOff);
}

// ======================================================================
// Phi estimation — streaming through .abed, no full-matrix materialization
// ======================================================================

static constexpr double PHI_MAF_CUTOFF = 0.01;

PhiMatrices estimatePhiOneAncestry(
    const AdmixData &admixData,
    const SparseGRM &grm,
    int ancIdx,
    int nthreads
) {
    const double mafCutoff = PHI_MAF_CUTOFF;
    uint32_t nUsed = admixData.nSubjUsed();
    uint32_t nMarkers = admixData.nMarkers();
    if (nthreads < 1) nthreads = 1;

    // Build directed pair list from GRM entries (both directions for off-diagonal)
    struct DPair {
        uint32_t i, j;
    };

    std::vector<DPair> pairs;
    pairs.reserve(grm.nnz() * 2);
    for (const auto &e : grm.entries()) {
        if (e.row != e.col) {
            pairs.push_back({e.row, e.col});
            pairs.push_back({e.col, e.row});
        }
    }
    const size_t nPairs = pairs.size();

    // Per-thread accumulator structure
    struct PairAcc {
        double ratioSum;
        uint32_t validCount;
    };

    // Worker function: processes markers [mBegin, mEnd) using its own cursor
    // and accumulates into local vectors.
    struct ThreadResult {
        std::vector<PairAcc> accA, accB, accC, accD;
        uint32_t snpsPassed = 0;
    };

    auto workerFn = [&](int workerId, uint32_t mBegin, uint32_t mEnd) -> ThreadResult {
        ThreadResult res;
        res.accA.assign(nPairs, {0.0, 0});
        res.accB.assign(nPairs, {0.0, 0});
        res.accC.assign(nPairs, {0.0, 0});
        res.accD.assign(nPairs, {0.0, 0});

        auto cursor = admixData.makeCursor();
        cursor->beginSequentialBlock(mBegin);
        Eigen::VectorXd dosage(nUsed), hapcount(nUsed);
        double mafHigh = 1.0 - mafCutoff;
        bool checkMissing = !admixData.hasNoMissing();
        uint32_t total = mEnd > mBegin ? mEnd - mBegin : 0;

        for (uint32_t m = mBegin; m < mEnd; ++m) {
            double q = cursor->getAdmixGenotypes(m, ancIdx, dosage, hapcount);

            if (q <= mafCutoff || q >= mafHigh) continue;
            ++res.snpsPassed;

            double qTerm = q * (1.0 - q);

            for (size_t p = 0; p < nPairs; ++p) {
                uint32_t i = pairs[p].i;
                uint32_t j = pairs[p].j;
                if (i >= nUsed || j >= nUsed) continue;

                double hi = hapcount[i];
                double hj = hapcount[j];
                double gi = dosage[i];
                double gj = dosage[j];

                if (checkMissing) {
                    if (!std::isfinite(hi) || !std::isfinite(hj) || !std::isfinite(gi) ||
                        !std::isfinite(gj)) continue;
                }

                int hiInt = static_cast<int>(std::round(hi));
                int hjInt = static_cast<int>(std::round(hj));

                double denom = hi * hj * qTerm;
                if (std::abs(denom) < 1e-15) continue;

                double numer = (gi - hi * q) * (gj - hj * q);
                double ratio = numer / denom;

                if (hiInt == 2 && hjInt == 2) {
                    res.accA[p].ratioSum += ratio;
                    res.accA[p].validCount++;
                } else if (hiInt == 2 && hjInt == 1) {
                    res.accB[p].ratioSum += ratio;
                    res.accB[p].validCount++;
                } else if (hiInt == 1 && hjInt == 2) {
                    res.accC[p].ratioSum += ratio;
                    res.accC[p].validCount++;
                } else if (hiInt == 1 && hjInt == 1) {
                    res.accD[p].ratioSum += ratio;
                    res.accD[p].validCount++;
                }
            }

            uint32_t finished = m - mBegin + 1;
            if (finished % 10000 == 0 || finished == total) infoMsg(
                    "  worker %d finished %u/%u markers",
                    workerId,
                    finished,
                    total
            );
        }
        return res;
    };

    // Split markers across threads
    std::vector<ThreadResult> results(nthreads);
    std::vector<uint32_t> workerMarkerCount(nthreads);
    uint32_t chunkSize = (nMarkers + nthreads - 1) / nthreads;
    for (int t = 0; t < nthreads; ++t) {
        uint32_t mBegin = t * chunkSize;
        uint32_t mEnd = std::min(mBegin + chunkSize, nMarkers);
        workerMarkerCount[t] = (mEnd > mBegin ? mEnd - mBegin : 0);
    }
    infoMsg("  workers 0-%d will each process ~%u markers", nthreads - 1, chunkSize);
    if (nthreads == 1) {
        results[0] = workerFn(0, 0, nMarkers);
    } else {
        std::vector<std::thread> workers;
        workers.reserve(nthreads - 1);
        for (int t = 0; t < nthreads - 1; ++t) {
            uint32_t mBegin = t * chunkSize;
            uint32_t mEnd = std::min(mBegin + chunkSize, nMarkers);
            workers.emplace_back(
                [&results, &workerFn, t, mBegin, mEnd]() {
                results[t] = workerFn(t, mBegin, mEnd);
            }
            );
        }
        // Main thread handles last chunk
        {
            uint32_t mBegin = (nthreads - 1) * chunkSize;
            uint32_t mEnd = nMarkers;
            results[nthreads - 1] = workerFn(nthreads - 1, mBegin, mEnd);
        }
        for (auto &w : workers)
            w.join();
    }

    // Merge accumulators
    std::vector<PairAcc> accA(nPairs, {0.0, 0});
    std::vector<PairAcc> accB(nPairs, {0.0, 0});
    std::vector<PairAcc> accC(nPairs, {0.0, 0});
    std::vector<PairAcc> accD(nPairs, {0.0, 0});
    uint32_t snpsPassed = 0;

    for (int t = 0; t < nthreads; ++t) {
        snpsPassed += results[t].snpsPassed;
        for (size_t p = 0; p < nPairs; ++p) {
            accA[p].ratioSum += results[t].accA[p].ratioSum;
            accA[p].validCount += results[t].accA[p].validCount;
            accB[p].ratioSum += results[t].accB[p].ratioSum;
            accB[p].validCount += results[t].accB[p].validCount;
            accC[p].ratioSum += results[t].accC[p].ratioSum;
            accC[p].validCount += results[t].accC[p].validCount;
            accD[p].ratioSum += results[t].accD[p].ratioSum;
            accD[p].validCount += results[t].accD[p].validCount;
        }
    }

    infoMsg("  Phi estimation: %u markers passed MAF filter (anc%d)", snpsPassed, ancIdx);

    // Build PhiMatrices from accumulators
    PhiMatrices result;
    auto buildEntries = [&](const std::vector<PairAcc> &acc, std::vector<PhiEntry> &out) {
        for (size_t p = 0; p < nPairs; ++p) {
            if (acc[p].validCount > 0) {
                double phi = acc[p].ratioSum / acc[p].validCount;
                out.push_back({pairs[p].i, pairs[p].j, phi});
            }
        }
    };
    buildEntries(accA, result.A);
    buildEntries(accB, result.B);
    buildEntries(accC, result.C);
    buildEntries(accD, result.D);

    infoMsg(
        "  Phi entries: A=%zu, B=%zu, C=%zu, D=%zu",
        result.A.size(),
        result.B.size(),
        result.C.size(),
        result.D.size()
    );

    return result;
}

// ======================================================================
// Variance computation with phi matrices
// ======================================================================

double computePhiVariance(
    const Eigen::VectorXd &R,
    const Eigen::VectorXd &hapcount,
    double q,
    const PhiMatrices &phi
) {
    double qTerm = q * (1.0 - q);
    double var = 0.0;

    // Off-diagonal terms (phi already has bidirectional pairs)
    // Scenario A: h_i=2, h_j=2 → multiplier = 4
    for (const auto &e : phi.A) {
        if (std::abs(hapcount[e.i] - 2.0) < 0.5 && std::abs(hapcount[e.j] - 2.0) < 0.5)var += 4.0 * qTerm * e.value * R[
                e.i] * R[e.j];
    }
    // Scenario B: h_i=2, h_j=1 → multiplier = 2
    for (const auto &e : phi.B) {
        if (std::abs(hapcount[e.i] - 2.0) < 0.5 && std::abs(hapcount[e.j] - 1.0) < 0.5)var += 2.0 * qTerm * e.value * R[
                e.i] * R[e.j];
    }
    // Scenario C: h_i=1, h_j=2 → multiplier = 2
    for (const auto &e : phi.C) {
        if (std::abs(hapcount[e.i] - 1.0) < 0.5 && std::abs(hapcount[e.j] - 2.0) < 0.5)var += 2.0 * qTerm * e.value * R[
                e.i] * R[e.j];
    }
    // Scenario D: h_i=1, h_j=1 → multiplier = 1
    for (const auto &e : phi.D) {
        if (std::abs(hapcount[e.i] - 1.0) < 0.5 && std::abs(hapcount[e.j] - 1.0) < 0.5)var += 1.0 * qTerm * e.value * R[
                e.i] * R[e.j];
    }

    // Diagonal: sum_i R_i^2 * h_i * q(1-q)
    for (int i = 0; i < static_cast<int>(R.size()); ++i) {
        var += R[i] * R[i] * hapcount[i] * qTerm;
    }

    return var;
}

// ======================================================================
// SPA p-value with outlier split (adapted from reference implementation)
// ======================================================================

namespace {

// Fused CGF derivatives — single exp() per call instead of 3 separate calls.
// K0 = h * log((1-p) + p*e^t),  K1 = h * p*e^t / base,  K2 = h * p*e^t*(1-p) / base^2
inline void kG012Local(
    double t,
    double maf,
    double h,
    double &K0out,
    double &K1out,
    double &K2out
) {
    double et = std::exp(std::clamp(t, -700.0, 700.0));
    double base = (1.0 - maf) + maf * et;
    if (base > 1e-15) {
        double pe = maf * et;     // p * e^t
        double q1p = (1.0 - maf); // (1-p)
        double bsq = base * base;
        K0out = h * std::log(base);
        K1out = h * pe / base;
        K2out = h * pe * q1p / bsq;
    } else {
        K0out = -std::numeric_limits<double>::infinity();
        K1out = 0.0;
        K2out = 0.0;
    }
}

// Newton-Raphson root finding for K'(t) = s
struct RootResult {
    double root;
    bool converged;
};

RootResult findRoot(
    double s,
    const double *rOut,
    const double *hOut,
    int nOut,
    double q,
    double meanNorm,
    double
    varNorm
) {
    static constexpr double tol = 0.001;
    static constexpr int maxIter = 100;
    double initVals[] = {0.0, -1.0, 1.0, -2.0, 2.0};

    for (double init : initVals) {
        double t = init;
        bool conv = false;
        for (int iter = 0; iter < maxIter; ++iter) {
            double K1 = 0.0, K2 = 0.0;
            for (int i = 0; i < nOut; ++i) {
                double tR = std::clamp(t * rOut[i], -700.0, 700.0);
                double k0i, k1i, k2i;
                kG012Local(tR, q, hOut[i], k0i, k1i, k2i);
                K1 += rOut[i] * k1i;
                K2 += rOut[i] * rOut[i] * k2i;
            }
            double K1total = K1 + meanNorm + varNorm * t - s;
            double K2total = K2 + varNorm;

            if (std::abs(K1total) < tol) {
                conv = true;
                break;
            }
            if (K2total <= 1e-10) break;

            double dt = std::clamp(-K1total / K2total, -2.0, 2.0);
            t = std::clamp(t + dt, -20.0, 20.0);
            if (std::abs(dt) < tol) {
                conv = true;
                break;
            }
        }
        if (conv) return {t, true};
    }
    return {0.0, false};
}

// Lugannani-Rice tail probability
double lugannamiRicePval(
    double zeta,
    double s,
    double q,
    const double *rOut,
    const double *hOut,
    int nOut,
    double meanNorm,
    double varNorm,
    bool upperTail
) {
    double K0 = 0.0, K2 = 0.0;
    for (int i = 0; i < nOut; ++i) {
        double tR = std::clamp(zeta * rOut[i], -700.0, 700.0);
        double k0i, k1i, k2i;
        kG012Local(tR, q, hOut[i], k0i, k1i, k2i);
        K0 += k0i;
        K2 += rOut[i] * rOut[i] * k2i;
    }
    K0 += meanNorm * zeta + 0.5 * varNorm * zeta * zeta;
    K2 += varNorm;

    double temp = zeta * s - K0;
    if (temp <= 0 || K2 <= 0) {
        // Fallback to normal
        double z = std::abs(s) / std::sqrt(K2 > 0 ? K2 : 1.0);
        boost::math::normal_distribution<> norm;
        return boost::math::cdf(boost::math::complement(norm, z));
    }

    double w = std::copysign(std::sqrt(2.0 * temp), zeta);
    double v = zeta * std::sqrt(K2);

    if (std::abs(w) < 1e-12 || std::abs(v) < 1e-12 || !std::isfinite(w) || !std::isfinite(v)) return MIN_P_VALUE;

    double lr_arg = w + (1.0 / w) * std::log(v / w);
    if (!std::isfinite(lr_arg)) return MIN_P_VALUE;

    boost::math::normal_distribution<> norm;
    if (upperTail) {
        return boost::math::cdf(boost::math::complement(norm, lr_arg));
    } else {
        return boost::math::cdf(norm, lr_arg);
    }
}

} // anonymous namespace

std::pair<double, double> spaLocalPval(
    double S,
    double sMean,
    double varDiag,
    const Eigen::Ref<const Eigen::VectorXd> &R,
    const Eigen::Ref<const Eigen::VectorXd> &hapcount,
    double q,
    double varS,
    const OutlierData &outlier,
    double spaCutoff
) {
    double z = (varS > 0.0) ? (S - sMean) / std::sqrt(varS) : 0.0;
    boost::math::normal_distribution<> norm;
    double pNorm = 2.0 * boost::math::cdf(boost::math::complement(norm, std::abs(z)));
    if (pNorm <= 0.0) pNorm = MIN_P_VALUE;

    if (std::abs(z) < spaCutoff || outlier.posOutlier.empty()) {
        return {pNorm, pNorm};
    }

    // varDiag already computed by caller
    double varRatio = (varS > 0.0) ? varDiag / varS : 1.0;
    double sNew = S * std::sqrt(varRatio);
    double sMeanNew = sMean * std::sqrt(varRatio);

    // Extract outlier data
    int nOut = static_cast<int>(outlier.posOutlier.size());
    std::vector<double> rOut(nOut), hOut(nOut);
    double meanNorm = 0.0, varNorm = 0.0;

    for (int i = 0; i < nOut; ++i) {
        uint32_t idx = outlier.posOutlier[i];
        rOut[i] = R[idx];
        hOut[i] = hapcount[idx];
    }
    for (uint32_t idx : outlier.posNonOutlier) {
        meanNorm += R[idx] * q * hapcount[idx];
        varNorm += R[idx] * R[idx] * q * (1.0 - q) * hapcount[idx];
    }

    double sUpper = std::max(sNew, 2.0 * sMeanNew - sNew);
    double sLower = std::min(sNew, 2.0 * sMeanNew - sNew);

    auto rootUpper = findRoot(sUpper, rOut.data(), hOut.data(), nOut, q, meanNorm, varNorm);
    auto rootLower = findRoot(sLower, rOut.data(), hOut.data(), nOut, q, meanNorm, varNorm);

    double pval = 0.0;
    if (rootUpper.converged) {
        pval += lugannamiRicePval(rootUpper.root, sUpper, q, rOut.data(), hOut.data(), nOut, meanNorm, varNorm, true);
    }
    if (rootLower.converged) {
        pval += lugannamiRicePval(rootLower.root, sLower, q, rOut.data(), hOut.data(), nOut, meanNorm, varNorm, false);
    }

    if (!rootUpper.converged && !rootLower.converged) pval = std::numeric_limits<double>::quiet_NaN();

    if (std::isfinite(pval)) {
        if (pval <= 0.0) pval = MIN_P_VALUE;
        if (pval > 1.0) pval = 1.0;
    }

    return {pval, pNorm};
}

// ======================================================================
// Phi I/O — wide format (single file, all ancestries)
//
// Header: idx1\tidx2\tanc0_A\tanc0_B\tanc0_C\tanc0_D\tanc1_A\t...
// Indices are 0-based into .fam order.  No .id sidecar needed.
// ======================================================================

static void writePhiWide(
    const std::string &path,
    const std::vector<PhiMatrices> &allPhi,
    int K
) {
    // Collect union of all (i,j) pairs across ancestries and scenarios
    struct PairKey {
        uint32_t i, j;
    };

    auto pairHash = [](const PairKey &p) {
        return std::hash<uint64_t>{}(uint64_t(p.i) << 32 | p.j);
    };
    auto pairEq = [](const PairKey &a, const PairKey &b) {
        return a.i == b.i && a.j == b.j;
    };
    std::unordered_map<PairKey, std::vector<double>, decltype(pairHash), decltype(pairEq)> rows(0, pairHash, pairEq);

    int nCols = K * 4; // A,B,C,D per ancestry

    auto insertEntries = [&](const std::vector<PhiEntry> &entries, int colIdx) {
        for (const auto &e : entries) {
            auto [it, inserted] =
                rows.try_emplace(
                    {e.i, e.j},
                    std::vector<double>(
                        nCols,
                        std::numeric_limits<double>::quiet_NaN()
                    )
                );
            it->second[colIdx] = e.value;
        }
    };

    for (int k = 0; k < K; ++k) {
        insertEntries(allPhi[k].A, k * 4 + 0);
        insertEntries(allPhi[k].B, k * 4 + 1);
        insertEntries(allPhi[k].C, k * 4 + 2);
        insertEntries(allPhi[k].D, k * 4 + 3);
    }

    // Write
    TextWriter out(path);

    // Header
    std::string hdr = "idx1\tidx2";
    const char *scenarioTag[] = {"_A", "_B", "_C", "_D"};
    for (int k = 0; k < K; ++k)
        for (int s = 0; s < 4; ++s)
            hdr += std::string("\tanc") + std::to_string(k) + scenarioTag[s];
    hdr += '\n';
    out.write(hdr);

    // Data rows
    char buf[64];
    for (const auto &[pair, vals] : rows) {
        std::string line;
        line.reserve(64 + 16 * nCols);
        std::snprintf(buf, sizeof(buf), "%u\t%u", pair.i, pair.j);
        line += buf;
        for (int c = 0; c < nCols; ++c) {
            if (std::isnan(vals[c]))line += "\tNA";
            else {
                std::snprintf(buf, sizeof(buf), "\t%.17g", vals[c]);
                line += buf;
            }
        }
        line += '\n';
        out.write(line);
    }

    infoMsg("Phi written: %zu pairs x %d columns -> %s", rows.size(), nCols, path.c_str());
}

static std::vector<PhiMatrices> readPhiWide(
    const std::string &path,
    int K
) {
    TextReader reader(path);

    // Parse header to determine column mapping
    std::string header;
    if (!reader.getline(header)) throw std::runtime_error("Cannot read phi file header: " + path);
    // Validate header: first token must be idx1 or #idx1
    {
        text::TokenScanner hts(header);
        auto firstCol = hts.nextView();
        if (firstCol != "idx1" && firstCol != "#idx1")throw std::runtime_error(path +
                                                                               ": header must start with idx1 or #idx1, got '"
                                                                               + std::string(firstCol) + "'");
    }

    int nCols = K * 4;
    std::vector<PhiMatrices> result(K);

    std::string line;
    while (reader.getline(line)) {
        if (line.empty()) continue;
        text::TokenScanner ts(line);
        ts.skipWS();
        char *ep;
        uint32_t i = static_cast<uint32_t>(std::strtoul(ts.pos(), &ep, 10));
        ts.p = ep;
        ts.skipWS();
        uint32_t j = static_cast<uint32_t>(std::strtoul(ts.pos(), &ep, 10));
        ts.p = ep;

        for (int c = 0; c < nCols; ++c) {
            auto sv = ts.nextView();
            if (sv.empty() || sv == "NA") continue;
            double val = std::strtod(sv.data(), &ep);
            int k = c / 4;
            int s = c % 4;
            PhiEntry entry{i, j, val};
            switch (s) {
            case 0:
                result[k].A.push_back(entry);
                break;
            case 1:
                result[k].B.push_back(entry);
                break;
            case 2:
                result[k].C.push_back(entry);
                break;
            case 3:
                result[k].D.push_back(entry);
                break;
            }
        }
    }

    for (int k = 0; k < K; ++k)
        infoMsg(
            "  anc%d phi: A=%zu B=%zu C=%zu D=%zu",
            k,
            result[k].A.size(),
            result[k].B.size(),
            result[k].C.size(),
            result[k].D.size()
        );

    return result;
}

// ======================================================================
// runPhiEstimation — full pipeline (writes single wide phi file)
// ======================================================================

void runPhiEstimation(
    const std::string &admixPrefix,
    const std::string &grmGrabFile,
    const std::string &grmGctaFile,
    const std::string &phiOutputFile,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &extractFile,
    const std::string &excludeFile,
    int nthreads
) {
    infoMsg("=== SPAmixLocalPlus Phi Estimation ===");

    // Load .fam IIDs (genotype subject list)
    auto famIIDs = parseFamIIDs(admixPrefix + ".fam");
    uint32_t nFam = static_cast<uint32_t>(famIIDs.size());
    infoMsg("Genotype (.fam): %u subjects", nFam);

    // Build filtered subject set via SubjectSet:
    // genotype → GRM (via parseSubjectIDs) → keep → remove
    SubjectSet ss(famIIDs);

    auto grmIDs = SparseGRM::parseSubjectIDs(grmGrabFile, grmGctaFile, famIIDs);
    if (!grmIDs.empty()) ss.setGrmSubjects(std::move(grmIDs));

    ss.setGrmLabel(grmFlagLabel(grmGrabFile, grmGctaFile));
    ss.setKeepRemove(keepFile, removeFile);
    ss.finalize();

    const uint32_t nUsed = ss.nUsed();
    const auto &usedMask = ss.usedMask();
    auto usedIIDs = ss.usedIIDs();

    // Load GRM against filtered subjects only
    auto grm = SparseGRM::load(grmGrabFile, grmGctaFile, usedIIDs, famIIDs);
    infoMsg("GRM loaded: %u subjects (dimension), %zu entries", grm.nSubjects(), grm.nnz());

    // Load admix data using intersection mask
    AdmixData admixData(admixPrefix, usedMask, nFam, nUsed, extractFile, excludeFile);
    int K = admixData.nAncestries();
    infoMsg("Ancestries: %d, Markers: %u, Subjects: %u", K, admixData.nMarkers(), admixData.nSubjUsed());

    // Estimate phi for all ancestries
    std::vector<PhiMatrices> allPhi(K);
    for (int k = 0; k < K; ++k) {
        infoMsg("Estimating phi for anc%d (%d thread%s)...", k, nthreads, nthreads > 1 ? "s" : "");
        allPhi[k] = estimatePhiOneAncestry(admixData, grm, k, nthreads);
    }

    // Write single wide file
    writePhiWide(phiOutputFile, allPhi, K);

    infoMsg("Output written to: %s", phiOutputFile.c_str());
}

// ======================================================================
// Unified GWAS — chunk-parallel × intra-chunk multi-phenotype loop.
//
// One pass over all admix chunks serves K_pheno phenotypes simultaneously.
// Each chunk:
//   - decoded ONCE per worker, shared across phenotypes
//   - Phase 1: per-(b, k) anc-shared scalars (maf / mac / missRate, hInt fill),
//              then three fused GEMMs (S, H·R, H·R^2) over (K_anc × K_pheno)
//   - Phase 2: one K_pheno-fused phi scan per ancestry
//              (kernel: computeVarOffMultiPhenoBatch)
//   - Phase 3: per-(b, p, k) SPA branch + output formatting into K_pheno
//              per-phenotype buffers
//
// Output: K_pheno files written in marker (chunk) order, byte-identical at
// the "%.6g" level to running each phenotype through the single-phenotype
// path.  Internal floating-point summation order in S / H·R / H·R^2 follows
// Eigen GEMM (small ULP-level difference from a hand-rolled per-(b,k,p)
// dot product); the phi base table is shared rather than R-baked, so the
// off-diagonal variance hot path differs only in inner multiplication
// order (mult·phi·R_p[i]·R_p[j] vs. (mult·phi·R[i]·R[j]) precomputed).
// ======================================================================

static void runUnifiedGWAS(
    const AdmixData &admixData,
    const std::vector<MultiPhenoRprodSoA> &rphi,
    const Eigen::MatrixXd &R_mat,
    const Eigen::MatrixXd &R2_mat,
    const std::vector<OutlierData> &outliers,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &outFiles,
    double spaCutoff,
    double missingCutoff,
    double mafCutoff,
    double macCutoff,
    const std::string &compression,
    int compressionLevel,
    int nthreads
) {
    const int K = admixData.nAncestries();
    const uint32_t nUsed = admixData.nSubjUsed();
    const int K_pheno = static_cast<int>(outliers.size());
    const auto &markerInfo = admixData.markerInfo();

    if (R_mat.rows() != static_cast<Eigen::Index>(nUsed) ||
        R2_mat.rows() != static_cast<Eigen::Index>(nUsed) ||
        R_mat.cols() != K_pheno || R2_mat.cols() != K_pheno) {
        throw std::runtime_error("runUnifiedGWAS: R_mat/R2_mat dimensions inconsistent with nUsed / K_pheno");
    }
    if (static_cast<int>(phenoNames.size()) != K_pheno ||
        static_cast<int>(outFiles.size()) != K_pheno) {
        throw std::runtime_error("runUnifiedGWAS: phenoNames / outFiles length must match K_pheno");
    }
    if (static_cast<int>(rphi.size()) != K) {
        throw std::runtime_error("runUnifiedGWAS: rphi size must equal nAncestries()");
    }
    for (int k = 0; k < K; ++k) {
        if (rphi[k].K_pheno != K_pheno) {
            throw std::runtime_error("runUnifiedGWAS: rphi[k].K_pheno must match phenotype count");
        }
    }

    // ── Output writers, one per phenotype ─────────────────────────────
    std::vector<TextWriter> writers;
    writers.reserve(K_pheno);
    for (int p = 0; p < K_pheno; ++p)
        writers.emplace_back(outFiles[p], TextWriter::modeFromString(compression), compressionLevel);

    // Header shared across phenotypes
    std::string header = "CHROM\tPOS\tID\tREF\tALT";
    for (int k = 0; k < K; ++k) {
        std::string pfx = "\tanc" + std::to_string(k) + "_";
        header += pfx + "MISS_RATE";
        header += pfx + "ALT_FREQ";
        header += pfx + "MAC";
        header += pfx + "P";
        header += pfx + "BETA";
        header += pfx + "SE";
    }
    header += '\n';
    for (int p = 0; p < K_pheno; ++p) writers[p].write(header);

    // ── Chunk parallel state ──────────────────────────────────────────
    const auto &chunks = admixData.chunkIndices();
    std::atomic<size_t> nextChunk{0};
    std::atomic<size_t> chunksCompleted{0};
    const size_t nChunks = chunks.size();
    const uint32_t nMarkers = admixData.nMarkers();
    const bool noMissing = admixData.hasNoMissing();

    struct PaddedFlag {
        alignas(64) int ready = 0;
    };

    // chunkOutputs[p][ci] = chunk ci's output buffer for phenotype p
    std::vector<std::vector<std::string> > chunkOutputs(K_pheno);
    for (int p = 0; p < K_pheno; ++p) chunkOutputs[p].resize(nChunks);
    std::vector<PaddedFlag> chunkReady(nChunks);

    std::mutex writeMutex;
    std::condition_variable writeCv;
    bool stopWriter = false;

    // Single writer thread: serializes emission per chunk in chunk order,
    // and across phenotypes in stable phenotype-index order within each chunk.
    std::thread writer([&]() {
        std::vector<std::string> mine(K_pheno);
        for (size_t ci = 0; ci < nChunks; ++ci) {
            {
                std::unique_lock<std::mutex> lk(writeMutex);
                writeCv.wait(lk, [&] { return chunkReady[ci].ready || stopWriter; });
                if (!chunkReady[ci].ready) continue;
                for (int p = 0; p < K_pheno; ++p)
                    mine[p] = std::move(chunkOutputs[p][ci]);
            }
            for (int p = 0; p < K_pheno; ++p) {
                writers[p].write(mine[p]);
                mine[p].clear();
                mine[p].shrink_to_fit();
            }
        }
    });


    // ── Worker function ───────────────────────────────────────────────
    auto workerFn = [&]() {
        auto cursor = admixData.makeCursor();
        char fmtBuf[64];

        // Per-batch decoded genotype storage (shared across phenotypes).
        // Allocated as ONE contiguous matrix per kind, with each batch slot
        // occupying K consecutive columns.  This lets Phase 1B issue a
        // single large GEMM per matrix (replacing PHI_BATCH×3 small GEMMs
        // whose per-call dispatch overhead otherwise dominates wall time
        // when K, K_pheno are small).
        Eigen::MatrixXd bDosBig(nUsed, static_cast<Eigen::Index>(K) * PHI_BATCH);
        Eigen::MatrixXd bHapBig(nUsed, static_cast<Eigen::Index>(K) * PHI_BATCH);

        // Subject-major hInt for batched phi scan: hIntSM[s * PHI_BATCH + b]
        std::vector<std::vector<uint32_t> > hIntPerAnc(
            K, std::vector<uint32_t>(static_cast<size_t>(nUsed) * PHI_BATCH, 0));

        // Per-batch per-ancestry shared scalar results (no S/sMean/diagVar —
        // those are now phenotype-dependent and live in the GEMM result blocks)
        struct AncScalar {
            double maf, missRate, dosSum, mac, q;
            bool pass;
        };
        std::vector<std::vector<AncScalar> > bAncS(PHI_BATCH, std::vector<AncScalar>(K));

        // GEMM result blocks: one big result per kind, of shape
        // (K * PHI_BATCH) × K_pheno.  Row index = b * K + k.
        Eigen::MatrixXd S_all(static_cast<Eigen::Index>(K) * PHI_BATCH, K_pheno);
        Eigen::MatrixXd HR_all(static_cast<Eigen::Index>(K) * PHI_BATCH, K_pheno);
        Eigen::MatrixXd HR2_all(static_cast<Eigen::Index>(K) * PHI_BATCH, K_pheno);

        // Phase 2 outputs: per-ancestry kernel buffer + 3-D accumulator
        std::vector<double> varOffKbuf(static_cast<size_t>(PHI_BATCH) * K_pheno);
        std::vector<double> varOffAll(static_cast<size_t>(PHI_BATCH) * K * K_pheno);

        // Per-phenotype chunk-output buffers (worker-local, transferred to
        // chunkOutputs[p][ci] under writeMutex at chunk completion).
        std::vector<std::string> bufPerPheno(K_pheno);

        for (size_t ci = nextChunk.fetch_add(1); ci < nChunks; ci = nextChunk.fetch_add(1)) {
            const auto &gIndices = chunks[ci];
            for (int p = 0; p < K_pheno; ++p) {
                bufPerPheno[p].clear();
                bufPerPheno[p].reserve(gIndices.size() * (80 + 90 * K));
            }

            cursor->beginSequentialBlock(gIndices.front());

            // Process markers in mini-batches of PHI_BATCH
            for (size_t mi = 0; mi < gIndices.size(); mi += PHI_BATCH) {
                int batchLen = static_cast<int>(std::min(
                    static_cast<size_t>(PHI_BATCH),
                    gIndices.size() - mi));

                // ── Phase 1A: decode + anc-shared scalars + hInt fill ──
                for (int b = 0; b < batchLen; ++b) {
                    uint64_t localIdx = gIndices[mi + b];
                    auto bDosView = bDosBig.middleCols(b * K, K);
                    auto bHapView = bHapBig.middleCols(b * K, K);
                    cursor->getAllAncestries(localIdx, bDosView, bHapView);

                    for (int k = 0; k < K; ++k) {
                        auto dosCol = bDosView.col(k);
                        auto hapCol = bHapView.col(k);

                        double hapSum = 0.0, dosSum = 0.0;
                        uint32_t nMissing = 0;
                        if (noMissing) {
                            dosSum = dosCol.sum();
                            hapSum = hapCol.sum();
                        } else {
                            for (uint32_t s = 0; s < nUsed; ++s) {
                                if (!std::isfinite(dosCol[s]) || !std::isfinite(hapCol[s])) {
                                    ++nMissing;
                                    dosCol[s] = 0.0;
                                    hapCol[s] = 0.0;
                                } else {
                                    dosSum += dosCol[s];
                                    hapSum += hapCol[s];
                                }
                            }
                        }

                        AncScalar &as = bAncS[b][k];
                        as.missRate = static_cast<double>(nMissing) / nUsed;
                        as.maf = (hapSum > 0) ? dosSum / hapSum : 0.0;
                        as.mac = std::min(dosSum, hapSum - dosSum);
                        as.dosSum = dosSum;

                        if (as.missRate > missingCutoff || as.maf < mafCutoff ||
                            as.maf > (1.0 - mafCutoff) || as.mac < macCutoff) {
                            as.pass = false;
                            as.q = 0.0;
                            // Zero out hInt slot so it won't match any scenario
                            uint32_t *dst = hIntPerAnc[k].data() + static_cast<size_t>(b);
                            for (uint32_t s = 0; s < nUsed; ++s)
                                dst[s * PHI_BATCH] = 0;
                            // GEMM still produces S_blk[b](k, *) etc. for the failed
                            // (b, k) column, but Phase 3 skips it via as.pass == false,
                            // so no further sanitization is needed here.
                            continue;
                        }
                        as.pass = true;
                        as.q = as.maf;

                        // Fill subject-major hInt for this batch slot
                        uint32_t *dst = hIntPerAnc[k].data() + static_cast<size_t>(b);
                        if (noMissing) {
                            for (uint32_t s = 0; s < nUsed; ++s)
                                dst[s * PHI_BATCH] = static_cast<uint32_t>(hapCol[s]);
                        } else {
                            for (uint32_t s = 0; s < nUsed; ++s) {
                                double h = hapCol[s];
                                dst[s * PHI_BATCH] = std::isfinite(h) ? static_cast<uint32_t>(h) : 0;
                            }
                        }
                    }
                }
                // Zero out unused trailing batch slots (no false matches)
                for (int b = batchLen; b < PHI_BATCH; ++b) {
                    for (int k = 0; k < K; ++k) {
                        uint32_t *dst = hIntPerAnc[k].data() + static_cast<size_t>(b);
                        for (uint32_t s = 0; s < nUsed; ++s)
                            dst[s * PHI_BATCH] = 0;
                    }
                }

                // ── Phase 1B: three big fused GEMMs across (PHI_BATCH × K_anc × K_pheno) ──
                // S_all   = bDosBig^T · R_mat    (K·PHI_BATCH × K_pheno)
                // HR_all  = bHapBig^T · R_mat    (K·PHI_BATCH × K_pheno)
                // HR2_all = bHapBig^T · R2_mat   (K·PHI_BATCH × K_pheno)
                //
                // The single large product replaces PHI_BATCH × 3 = 24 small
                // 3×3 GEMM dispatches; with N=10k, K=3, K_pheno=3 the per-call
                // dispatch overhead of the small variant dominates wall time,
                // whereas the large GEMM benefits from Eigen's blocked kernel.
                // Trailing batch slots (b ≥ batchLen) hold stale data; their
                // GEMM output rows are never read in Phase 3.
                S_all.noalias()   = bDosBig.transpose() * R_mat;
                HR_all.noalias()  = bHapBig.transpose() * R_mat;
                HR2_all.noalias() = bHapBig.transpose() * R2_mat;

                // ── Phase 2: K_pheno-fused phi scan, one pass per ancestry ──
                for (int k = 0; k < K; ++k) {
                    computeVarOffMultiPhenoBatch(
                        rphi[k], hIntPerAnc[k].data(),
                        batchLen, varOffKbuf.data());
                    for (int b = 0; b < batchLen; ++b) {
                        for (int p = 0; p < K_pheno; ++p) {
                            varOffAll[(static_cast<size_t>(b) * K + k) * K_pheno + p] =
                                varOffKbuf[static_cast<size_t>(b) * K_pheno + p];
                        }
                    }
                }

                // ── Phase 3: SPA branch + per-phenotype output formatting ──
                for (int b = 0; b < batchLen; ++b) {
                    uint64_t localIdx = gIndices[mi + b];
                    const auto &mInfo = markerInfo[localIdx];

                    // Compose marker meta once per (b)
                    char metaBuf[128];
                    int metaLen = std::snprintf(
                        metaBuf, sizeof(metaBuf), "%s\t%u\t%s\t%s\t%s",
                        mInfo.chrom.c_str(), mInfo.pos, mInfo.id.c_str(),
                        mInfo.ref.c_str(), mInfo.alt.c_str());

                    for (int p = 0; p < K_pheno; ++p) {
                        std::string &buf = bufPerPheno[p];
                        buf.append(metaBuf, metaLen);

                        for (int k = 0; k < K; ++k) {
                            const AncScalar &as = bAncS[b][k];

                            // MISS_RATE / ALT_FREQ / MAC (shared across phenotypes)
                            buf += '\t';
                            std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", as.missRate);
                            buf += fmtBuf;
                            buf += '\t';
                            std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", as.maf);
                            buf += fmtBuf;
                            buf += '\t';
                            std::snprintf(fmtBuf, sizeof(fmtBuf), "%.0f", as.mac);
                            buf += fmtBuf;

                            if (!as.pass) {
                                buf += "\tNA\tNA\tNA";
                                continue;
                            }

                            double q     = as.q;
                            double qTerm = q * (1.0 - q);
                            Eigen::Index gemmRow = static_cast<Eigen::Index>(b) * K + k;
                            double S       = S_all(gemmRow, p);
                            double sMean   = q * HR_all(gemmRow, p);
                            double diagVar = qTerm * HR2_all(gemmRow, p);
                            // varOff kernel returns sum(mult · phi · R_p[i] · R_p[j]) for matching
                            // entries; the q(1-q) factor is applied here, matching the single-
                            // phenotype path (see git history of runUnifiedGWAS).
                            double varOffRaw = varOffAll[(static_cast<size_t>(b) * K + k) * K_pheno + p];
                            double varS      = varOffRaw * qTerm + diagVar;

                            auto hapCol = bHapBig.col(b * K + k);
                            auto [pSpa, pNorm] = spaLocalPval(
                                S, sMean, diagVar,
                                R_mat.col(p), hapCol,
                                q, varS, outliers[p], spaCutoff);
                            (void)pNorm;
                            double betaG = (varS > 0.0) ? (S - sMean) / varS
                                                         : std::numeric_limits<double>::quiet_NaN();
                            double seG   = (varS > 0.0) ? 1.0 / std::sqrt(varS)
                                                         : std::numeric_limits<double>::quiet_NaN();

                            buf += '\t';
                            std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", pSpa);
                            buf += fmtBuf;
                            buf += '\t';
                            if (std::isfinite(betaG)) {
                                std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", betaG);
                                buf += fmtBuf;
                            } else {
                                buf += "NA";
                            }
                            buf += '\t';
                            if (std::isfinite(seG)) {
                                std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", seG);
                                buf += fmtBuf;
                            } else {
                                buf += "NA";
                            }
                        }
                        buf += '\n';
                    }
                }
            } // end mini-batch loop

            // Commit per-phenotype chunk outputs to the shared queue
            {
                std::lock_guard<std::mutex> lk(writeMutex);
                for (int p = 0; p < K_pheno; ++p)
                    chunkOutputs[p][ci] = std::move(bufPerPheno[p]);
                chunkReady[ci].ready = 1;
            }
            writeCv.notify_all();

            // Progress logging at ~25 / 50 / 75 %
            size_t done = chunksCompleted.fetch_add(1) + 1;
            if (nChunks >= 20) {
                size_t q1 = nChunks / 4, q2 = nChunks / 2, q3 = nChunks * 3 / 4;
                if (done == q1 || done == q2 || done == q3) {
                    uint32_t markersDone = static_cast<uint32_t>(
                        static_cast<uint64_t>(done) * nMarkers / nChunks);
                    infoMsg(
                        "    %u / %u markers (~%u%%)",
                        markersDone, nMarkers,
                        static_cast<unsigned>(done * 100 / nChunks));
                }
            }
        }
    };

    // Launch workers
    int nWorkers = std::max(1, nthreads);
    std::vector<std::thread> workers;
    workers.reserve(nWorkers);
    for (int t = 0; t < nWorkers; ++t)
        workers.emplace_back(workerFn);

    for (auto &w : workers) w.join();
    {
        std::lock_guard<std::mutex> lk(writeMutex);
        stopWriter = true;
    }
    writeCv.notify_all();
    writer.join();

    for (int p = 0; p < K_pheno; ++p)
        infoMsg("  Phenotype '%s': %u markers processed -> %s",
                phenoNames[p].c_str(), admixData.nMarkers(), outFiles[p].c_str());
}

// ======================================================================
// runSPAmixLocalPlus — main entry point
// ======================================================================

void runSPAmixLocalPlus(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &admixPrefix,
    const std::string &admixPhiFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierRatio,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &extractFile,
    const std::string &excludeFile,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    bool saveResid
) {
    infoMsg("=== SPAmixLocalPlus GWAS ===");

    // Fit-path detection: when --pheno-name is supplied, the null-model
    // engine fits residuals from the phenotype columns + covariates;
    // otherwise loadResidOne consumes pre-computed residuals via --resid-name.
    const bool fitPath = !phenoNameSpec.empty();
    nullmodel::RegressionModel regModel{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        regModel = nullmodel::parseRegressionModel(regressionModelStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
        infoMsg("SPAmixLocalPlus: fitting %s null model for %zu phenotype(s)",
                nullmodel::regressionModelName(regModel), phenoSpecs.size());
    }

    // Load subjects
    auto famIIDs = parseFamIIDs(admixPrefix + ".fam");
    SubjectData sd(famIIDs);
    if (fitPath) {
        // Single-pass load: gather every column the null-model engine needs
        // from the pheno file.  Covariates are read from the pheno file
        // when --covar is absent.
        std::vector<std::string> wanted;
        auto add = [&](const std::string &name) {
            if (name.empty()) return;
            if (std::find(wanted.begin(), wanted.end(), name) == wanted.end())
                wanted.push_back(name);
        };
        for (const auto &name : nullmodel::columnsNeeded(phenoSpecs)) add(name);
        if (covarFile.empty())
            for (const auto &name : covarNames) add(name);
        sd.loadPhenoFile(phenoFile, wanted);
    } else {
        sd.loadResidOne(phenoFile, residNames);
    }
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.finalize();

    uint32_t nUsed = sd.nUsed();
    infoMsg("Subjects: %u in .fam, %u used", sd.nFam(), nUsed);

    if (fitPath) {
        Eigen::MatrixXd covarUnion;
        if (!covarNames.empty()) {
            covarUnion = sd.getColumns(covarNames);
        } else {
            covarUnion.resize(sd.nUsed(), 0);
        }
        nullmodel::EngineOptions eo;
        eo.nthreads = nthread;
        auto fits = nullmodel::fitAll(sd, phenoSpecs, regModel, covarUnion, eo);
        std::vector<Eigen::VectorXd> rs;
        std::vector<std::string> ns;
        rs.reserve(fits.size());
        ns.reserve(fits.size());
        for (auto &f : fits) {
            infoMsg("  Fitted '%s': %d subjects after NaN removal",
                    f.name.c_str(), f.nUsedRows);
            rs.push_back(std::move(f.residuals));
            ns.push_back(f.name);
        }
        if (saveResid) {
            std::vector<nullmodel::NullModelFit> dumpFits(rs.size());
            for (size_t i = 0; i < rs.size(); ++i) {
                dumpFits[i].name = ns[i];
                dumpFits[i].residuals = rs[i];
                dumpFits[i].nUsedRows = static_cast<int>(rs[i].size());
            }
            nullmodel::writeResidualsFile(outPrefix + ".null.resid", sd, dumpFits);
        }
        sd.setResidualsFromFit(std::move(rs), std::move(ns));
    }

    // Load admix data
    AdmixData admixData(admixPrefix, sd.usedMask(), sd.nFam(), nUsed, extractFile, excludeFile, nSnpPerChunk);
    int K = admixData.nAncestries();
    infoMsg("Ancestries: %d, Markers: %u", K, admixData.nMarkers());

    // Load phi from wide file (all ancestries at once)
    auto allPhi = readPhiWide(admixPhiFile, K);

    // Phenotype set
    const int nRC = sd.residOneCols();
    if (nRC > 1) infoMsg("Multi-column residual file: %d phenotypes", nRC);
    auto phenoInfos = sd.buildPerColumnMasks();

    // R_mat (N × K_pheno) at union dimension, NaN→0 filled.
    Eigen::MatrixXd R_mat;
    if (nRC > 1) {
        R_mat = sd.residMatrix();
        for (Eigen::Index c = 0; c < R_mat.cols(); ++c)
            for (Eigen::Index s = 0; s < R_mat.rows(); ++s)
                if (std::isnan(R_mat(s, c))) R_mat(s, c) = 0.0;
    } else {
        R_mat.resize(sd.residuals().size(), 1);
        R_mat.col(0) = sd.residuals();
    }
    Eigen::MatrixXd R2_mat = R_mat.array().square().matrix();

    // Build multi-phenotype rprod tables once per ancestry; shared (const)
    // across all phenotypes and worker threads.  rprod_packed bakes
    // `mult · phi · R[i,p] · R[j,p]` so the kernel hot loop reads K_pheno
    // sequential doubles per entry — no random R gather.
    std::vector<MultiPhenoRprodSoA> rphi(K);
    for (int k = 0; k < K; ++k) rphi[k] = buildMultiPhenoRprodSoA(allPhi[k], R_mat);

    // Per-phenotype outliers, names, output paths
    std::vector<OutlierData>  outliers(nRC);
    std::vector<std::string>  phenoNames(nRC);
    std::vector<std::string>  outFiles(nRC);
    for (int p = 0; p < nRC; ++p) {
        const auto &pi = phenoInfos[p];
        Eigen::VectorXd col = R_mat.col(p);
        outliers[p]   = detectOutliers(col, outlierRatio);
        phenoNames[p] = pi.name;
        outFiles[p]   = TextWriter::buildOutputPath(outPrefix, pi.name, "LocalP", compression);
        infoMsg(
            "  Phenotype '%s': %u subjects, %u markers, %d ancestries -> %s",
            pi.name.c_str(), pi.nUsed, admixData.nMarkers(), K, outFiles[p].c_str());
    }

    runUnifiedGWAS(
        admixData,
        rphi,
        R_mat,
        R2_mat,
        outliers,
        phenoNames,
        outFiles,
        spaCutoff,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        compression,
        compressionLevel,
        nthread
    );
}
