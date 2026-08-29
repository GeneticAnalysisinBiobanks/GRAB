// bgen.cpp — BgenData and BgenCursor implementations
//
// Wraps the bgen reference implementation behind GenoMeta / GenoCursor.
// Each BgenCursor owns its own ifstream for thread safety.

#include "geno_factory/bgen.hpp"
#include "geno_factory/hwe.hpp"
#include "geno_factory/variant_filter.hpp"
#include "util/logging.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "genfile/bgen/bgen.hpp"
#include "libdeflate.h"
#include "sqlite3.h"
#include "util/simd_dispatch.hpp"
#include "zstd.h"

#include <sys/stat.h>

#if defined(__x86_64__) || defined(_M_X64)
#include <immintrin.h>
#endif

namespace {

// True if a regular file exists and its size (bytes) is returned via *size.
inline bool statFileSize(const std::string &path, uint64_t *size) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return false;
    if (size) *size = static_cast<uint64_t>(st.st_size);
    return true;
}

// ──────────────────────────────────────────────────────────────────────────────
// DosageSetter — callback that bgen's read_and_parse_genotype_data_block uses
// to deliver per-sample genotype probabilities.  We convert on the fly to a
// per-subject dosage, applying the usedMask.
//
// Supports both BGEN-1.2 unphased and phased layouts:
//
//   Unphased biallelic diploid (Z=3, ePerUnorderedGenotype):
//     idx 0 = P(AA)  [hom first allele]
//     idx 1 = P(AB)  [het]
//     idx 2 = P(BB)  [hom second allele]
//     dosage of second allele = idx_1 + 2*idx_2
//
//   Phased biallelic diploid (Z=4, ePerPhasedHaplotypePerAllele):
//     idx 0 = P(hap0 = first allele)
//     idx 1 = P(hap0 = second allele)
//     idx 2 = P(hap1 = first allele)
//     idx 3 = P(hap1 = second allele)
//     dosage of second allele = idx_1 + idx_3
//
// `countFirstAllele = false` (default) → dosage of second allele (alleles[1]).
// `countFirstAllele = true`            → dosage of first allele  (alleles[0]),
//   used by --bgen FILE ref-last (or ref-unknown) to compensate for plink2's
//   default BGEN export which writes ALT as the first allele.
// ──────────────────────────────────────────────────────────────────────────────

struct DosageSetter {
    // Outputs
    double *out = nullptr;
    std::vector<uint32_t> *missingIdx = nullptr; // may be nullptr for simple path

    // Config
    const std::vector<uint64_t> *usedMask = nullptr;
    uint32_t nSamplesInFile = 0;
    bool allUsed = true;
    bool countFirstAllele = false; // true ⇒ count alleles[0], else count alleles[1]

    // ⑤ fast-path inputs (used only by fastParseV12DiploidBiallelic — the
    // genfile setter-callback fallback below ignores them):
    const std::vector<uint32_t> *usedIndices = nullptr; // ascending used indices
    double *scratchAll = nullptr; // scratch sized nSamplesInFile, used when !allUsed
    uint32_t nUsed = 0;           // == usedIndices->size() when !allUsed
    bool needStats = true;        // false ⇒ getGenotypesSimple (fill out only)

    // Per-sample state
    uint32_t outIdx = 0;
    uint32_t currentSample = 0;
    bool currentUsed = true;
    bool phasedCurrent = false;
    uint32_t expectedEntries = 0;
    uint32_t probIdx = 0;
    double probSum = 0.0;       // running dosage of the "counted" allele
    bool sampleMissing = false; // any MissingValue seen for this sample
    bool sampleFinalized = false;

    // Counts
    uint32_t nHomRef = 0;
    uint32_t nHet = 0;
    uint32_t nHomAlt = 0;
    uint32_t nMissing = 0;
    double dosageSum = 0.0;  // Σ dosage over non-missing used samples (for dosage AF)
    // Set when any decoded dosage is not exactly 0/1/2.  BGEN quantises
    // probabilities to 8 bits, so a hard-call variant decodes to exact
    // integers and this stays false; recording it here costs one test on a
    // value already in a register, where the alternative is a second pass
    // over the column further downstream.
    bool anyFractional = false;
    // plink2 hard-call threshold: a non-missing dosage counts toward the HWE
    // hom/het buckets only when |dosage - round| <= hardCallThr; otherwise it
    // is HWE-uncertain and excluded (still summed into dosageSum for AF).
    double hardCallThr = kHardCallThreshold;

    void initialise(
        std::size_t /*nsamples*/,
        std::size_t                                       /*nalleles*/
    ) {
        outIdx = 0;
        nHomRef = nHet = nHomAlt = nMissing = 0;
        dosageSum = 0.0;
        anyFractional = false;
    }

    bool set_sample(std::size_t i) {
        currentSample = static_cast<uint32_t>(i);
        currentUsed = allUsed || (((*usedMask)[currentSample / 64] >> (currentSample % 64)) & 1);
        probIdx = 0;
        probSum = 0.0;
        sampleMissing = false;
        sampleFinalized = false;
        return true; // always receive data to advance correctly
    }

    void set_number_of_entries(
        std::size_t /*ploidy*/,
        std::size_t Z,
        genfile::OrderType order,
        genfile::ValueType                        /*valueType*/
    ) {
        phasedCurrent = (order == genfile::ePerPhasedHaplotypePerAllele);
        expectedEntries = static_cast<uint32_t>(Z);
    }

    void finalizeSample() {
        if (sampleFinalized) return;
        sampleFinalized = true;
        if (!currentUsed) return;
        if (sampleMissing) {
            out[outIdx] = std::numeric_limits<double>::quiet_NaN();
            if (missingIdx) missingIdx->push_back(outIdx);
            ++nMissing;
        } else {
            double dosage = countFirstAllele ? (2.0 - probSum) : probSum;
            // Clamp to [0, 2] in case of small numerical drift in phased branch.
            if (dosage < 0.0) dosage = 0.0;
            else if (dosage > 2.0) dosage = 2.0;
            out[outIdx] = dosage;
            dosageSum += dosage;
            if (dosage != 0.0 && dosage != 1.0 && dosage != 2.0) anyFractional = true;
            const int hc = dosageHardcall(dosage, hardCallThr);
            if (hc == 0) ++nHomRef;
            else if (hc == 1) ++nHet;
            else if (hc == 2) ++nHomAlt;
            // else: HWE-uncertain — excluded from hom/het (still in dosageSum).
        }
        ++outIdx;
    }

    void set_value(
        uint32_t idx,
        double value
    ) {
        if (currentUsed && !sampleMissing) {
            if (phasedCurrent) {
                // Z=4 phased: idx 1 = P(hap0=allele1), idx 3 = P(hap1=allele1).
                if (idx == 1 || idx == 3) probSum += value;
            } else {
                // Z=3 unphased: idx 1 = P(het), idx 2 = P(hom allele1).
                if (idx == 1)
                    probSum += value;
                else if (idx == 2)
                    probSum += 2.0 * value;
            }
        }
        ++probIdx;
        if (probIdx >= expectedEntries) finalizeSample();
    }

    void set_value(
        uint32_t /*idx*/,
        genfile::MissingValue
    ) {
        sampleMissing = true;
        ++probIdx;
        if (probIdx >= expectedEntries) finalizeSample();
    }

    void finalise() {
    }

};

// Lookup table for the 8-bit probability decode: t[b] == double(b)/255.0,
// computed once.  Replacing the per-sample division with a table load is
// bit-for-bit identical (same operation, same rounding) and noticeably faster
// in the inner loop over hundreds of thousands of samples.
inline const double *prob8Table() {
    static const std::array<double, 256> table = [] {
        std::array<double, 256> a{};
        for (int i = 0; i < 256; ++i) a[i] = static_cast<double>(i) / 255.0;
        return a;
    }();
    return table.data();
}

// ── Per-sample decode kernels ───────────────────────────────────────────────
// Decode ALL N file samples of a v12 diploid/biallelic/unphased block into
// dst[0..N): dst[i] = dosage of the counted allele, or NaN when sample i is
// missing (ploidy byte high bit set).  No used-sample filtering and no
// statistics — the caller gathers used samples and computes QC in a second pass.
//
// The per-sample formula is: v = byte/255 (8-bit, via prob8Table) or u16/65535
// (16-bit); p2 = max((1 - v1) - v2, 0); dosage = v2 + 2*p2 (count second
// allele), or 2 - that under ref-last (countFirst); clamp to [0,2].  The
// _avx2/_avx512 variants below reproduce this bit-for-bit: the divisions are the
// same correctly-rounded IEEE operations, and the only multiply is *2.0 which is
// exact (a pure exponent increment), so FMA-vs-mul+add cannot change the result.

// Scalar per-sample helpers, reused verbatim by the SIMD remainder loops so the
// tail is bit-identical to the SIMD body.
inline double decodeDosage8(
    const genfile::byte_t *probs,
    const genfile::byte_t *ploidy,
    uint32_t i,
    bool countFirst,
    const double *t8
) {
    if (ploidy[i] & 0x80) return std::numeric_limits<double>::quiet_NaN();
    const genfile::byte_t *p = probs + static_cast<size_t>(i) * 2;
    const double v1 = t8[p[0]];
    const double v2 = t8[p[1]];
    const double p2 = std::max(1.0 - v1 - v2, 0.0);
    double dosage = v2 + 2.0 * p2;
    if (countFirst) dosage = 2.0 - dosage;
    if (dosage < 0.0) dosage = 0.0;
    else if (dosage > 2.0) dosage = 2.0;
    return dosage;
}

inline double decodeDosage16(
    const genfile::byte_t *probs,
    const genfile::byte_t *ploidy,
    uint32_t i,
    bool countFirst
) {
    if (ploidy[i] & 0x80) return std::numeric_limits<double>::quiet_NaN();
    const genfile::byte_t *p = probs + static_cast<size_t>(i) * 4;
    const uint32_t u1 = static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8);
    const uint32_t u2 = static_cast<uint32_t>(p[2]) | (static_cast<uint32_t>(p[3]) << 8);
    const double v1 = static_cast<double>(u1) / 65535.0;
    const double v2 = static_cast<double>(u2) / 65535.0;
    const double p2 = std::max(1.0 - v1 - v2, 0.0);
    double dosage = v2 + 2.0 * p2;
    if (countFirst) dosage = 2.0 - dosage;
    if (dosage < 0.0) dosage = 0.0;
    else if (dosage > 2.0) dosage = 2.0;
    return dosage;
}

void parseKernel8_scalar(
    const genfile::byte_t *probs,
    const genfile::byte_t *ploidy,
    uint32_t N,
    bool countFirst,
    double *dst
) {
    const double *const t8 = prob8Table();
    for (uint32_t i = 0; i < N; ++i) dst[i] = decodeDosage8(probs, ploidy, i, countFirst, t8);
}

void parseKernel16_scalar(
    const genfile::byte_t *probs,
    const genfile::byte_t *ploidy,
    uint32_t N,
    bool countFirst,
    double *dst
) {
    for (uint32_t i = 0; i < N; ++i) dst[i] = decodeDosage16(probs, ploidy, i, countFirst);
}

#if defined(__x86_64__) || defined(_M_X64)

// ② AVX2 / AVX-512 variants.  Double precision throughout; byte→double via
// cvtepu8/epu16 → cvtepi32_pd → div, matching the scalar table/division exactly.
// The *2.0 is exact, so explicit mul+add reproduces the scalar result.  Missing
// samples (ploidy high bit) are blended to NaN.  Remainder handled by the scalar
// per-sample helper, identical to the SIMD body.

__attribute__((target("avx2,fma")))
void parseKernel8_avx2(
    const genfile::byte_t *probs,
    const genfile::byte_t *ploidy,
    uint32_t N,
    bool countFirst,
    double *dst
) {
    const __m256d v255 = _mm256_set1_pd(255.0);
    const __m256d vone = _mm256_set1_pd(1.0);
    const __m256d vtwo = _mm256_set1_pd(2.0);
    const __m256d vzero = _mm256_setzero_pd();
    const __m256d vnan = _mm256_set1_pd(std::numeric_limits<double>::quiet_NaN());
    const __m128i evenMask =
        _mm_setr_epi8(0, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i oddMask =
        _mm_setr_epi8(1, 3, 5, 7, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i v0x80 = _mm_set1_epi32(0x80);
    uint32_t i = 0;
    for (; i + 4 <= N; i += 4) {
        const __m128i raw =
            _mm_loadl_epi64(reinterpret_cast<const __m128i *>(probs + static_cast<size_t>(i) * 2));
        const __m256d v1 = _mm256_div_pd(
            _mm256_cvtepi32_pd(_mm_cvtepu8_epi32(_mm_shuffle_epi8(raw, evenMask))), v255);
        const __m256d v2 = _mm256_div_pd(
            _mm256_cvtepi32_pd(_mm_cvtepu8_epi32(_mm_shuffle_epi8(raw, oddMask))), v255);
        const __m256d p2 =
            _mm256_max_pd(_mm256_sub_pd(_mm256_sub_pd(vone, v1), v2), vzero);
        __m256d dosage = _mm256_add_pd(v2, _mm256_mul_pd(vtwo, p2));
        if (countFirst) dosage = _mm256_sub_pd(vtwo, dosage);
        dosage = _mm256_min_pd(_mm256_max_pd(dosage, vzero), vtwo);
        const __m128i pl =
            _mm_loadl_epi64(reinterpret_cast<const __m128i *>(ploidy + i)); // low 4 used
        const __m128i miss =
            _mm_cmpgt_epi32(_mm_and_si128(_mm_cvtepu8_epi32(pl), v0x80), _mm_setzero_si128());
        dosage = _mm256_blendv_pd(dosage, vnan,
                                  _mm256_castsi256_pd(_mm256_cvtepi32_epi64(miss)));
        _mm256_storeu_pd(dst + i, dosage);
    }
    const double *const t8 = prob8Table();
    for (; i < N; ++i) dst[i] = decodeDosage8(probs, ploidy, i, countFirst, t8);
}

__attribute__((target("avx2,avx512f,avx512vl,fma")))
void parseKernel8_avx512(
    const genfile::byte_t *probs,
    const genfile::byte_t *ploidy,
    uint32_t N,
    bool countFirst,
    double *dst
) {
    const __m512d v255 = _mm512_set1_pd(255.0);
    const __m512d vone = _mm512_set1_pd(1.0);
    const __m512d vtwo = _mm512_set1_pd(2.0);
    const __m512d vzero = _mm512_setzero_pd();
    const __m512d vnan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
    const __m128i evenMask =
        _mm_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m128i oddMask =
        _mm_setr_epi8(1, 3, 5, 7, 9, 11, 13, 15, -1, -1, -1, -1, -1, -1, -1, -1);
    const __m256i v0x80 = _mm256_set1_epi32(0x80);
    uint32_t i = 0;
    for (; i + 8 <= N; i += 8) {
        const __m128i raw =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(probs + static_cast<size_t>(i) * 2));
        const __m512d v1 = _mm512_div_pd(
            _mm512_cvtepi32_pd(_mm256_cvtepu8_epi32(_mm_shuffle_epi8(raw, evenMask))), v255);
        const __m512d v2 = _mm512_div_pd(
            _mm512_cvtepi32_pd(_mm256_cvtepu8_epi32(_mm_shuffle_epi8(raw, oddMask))), v255);
        const __m512d p2 =
            _mm512_max_pd(_mm512_sub_pd(_mm512_sub_pd(vone, v1), v2), vzero);
        __m512d dosage = _mm512_add_pd(v2, _mm512_mul_pd(vtwo, p2));
        if (countFirst) dosage = _mm512_sub_pd(vtwo, dosage);
        dosage = _mm512_min_pd(_mm512_max_pd(dosage, vzero), vtwo);
        const __m128i pl =
            _mm_loadl_epi64(reinterpret_cast<const __m128i *>(ploidy + i)); // 8 bytes
        const __mmask8 miss =
            _mm256_test_epi32_mask(_mm256_cvtepu8_epi32(pl), v0x80);
        dosage = _mm512_mask_blend_pd(miss, dosage, vnan);
        _mm512_storeu_pd(dst + i, dosage);
    }
    const double *const t8 = prob8Table();
    for (; i < N; ++i) dst[i] = decodeDosage8(probs, ploidy, i, countFirst, t8);
}

__attribute__((target("avx2,fma")))
void parseKernel16_avx2(
    const genfile::byte_t *probs,
    const genfile::byte_t *ploidy,
    uint32_t N,
    bool countFirst,
    double *dst
) {
    const __m256d v65535 = _mm256_set1_pd(65535.0);
    const __m256d vone = _mm256_set1_pd(1.0);
    const __m256d vtwo = _mm256_set1_pd(2.0);
    const __m256d vzero = _mm256_setzero_pd();
    const __m256d vnan = _mm256_set1_pd(std::numeric_limits<double>::quiet_NaN());
    const __m128i v0xffff = _mm_set1_epi32(0xFFFF);
    const __m128i v0x80 = _mm_set1_epi32(0x80);
    const int *const base = reinterpret_cast<const int *>(probs);
    uint32_t i = 0;
    for (; i + 4 <= N; i += 4) {
        // One 32-bit dword per sample (stride 4) packs both u16 probabilities.
        const __m128i idx = _mm_setr_epi32(i, i + 1, i + 2, i + 3);
        const __m128i dw = _mm_i32gather_epi32(base, idx, 4);
        const __m256d v1 =
            _mm256_div_pd(_mm256_cvtepi32_pd(_mm_and_si128(dw, v0xffff)), v65535);
        const __m256d v2 =
            _mm256_div_pd(_mm256_cvtepi32_pd(_mm_srli_epi32(dw, 16)), v65535);
        const __m256d p2 =
            _mm256_max_pd(_mm256_sub_pd(_mm256_sub_pd(vone, v1), v2), vzero);
        __m256d dosage = _mm256_add_pd(v2, _mm256_mul_pd(vtwo, p2));
        if (countFirst) dosage = _mm256_sub_pd(vtwo, dosage);
        dosage = _mm256_min_pd(_mm256_max_pd(dosage, vzero), vtwo);
        const __m128i pl =
            _mm_loadl_epi64(reinterpret_cast<const __m128i *>(ploidy + i));
        const __m128i miss =
            _mm_cmpgt_epi32(_mm_and_si128(_mm_cvtepu8_epi32(pl), v0x80), _mm_setzero_si128());
        dosage = _mm256_blendv_pd(dosage, vnan,
                                  _mm256_castsi256_pd(_mm256_cvtepi32_epi64(miss)));
        _mm256_storeu_pd(dst + i, dosage);
    }
    for (; i < N; ++i) dst[i] = decodeDosage16(probs, ploidy, i, countFirst);
}

__attribute__((target("avx2,avx512f,avx512vl,fma")))
void parseKernel16_avx512(
    const genfile::byte_t *probs,
    const genfile::byte_t *ploidy,
    uint32_t N,
    bool countFirst,
    double *dst
) {
    const __m512d v65535 = _mm512_set1_pd(65535.0);
    const __m512d vone = _mm512_set1_pd(1.0);
    const __m512d vtwo = _mm512_set1_pd(2.0);
    const __m512d vzero = _mm512_setzero_pd();
    const __m512d vnan = _mm512_set1_pd(std::numeric_limits<double>::quiet_NaN());
    const __m256i v0xffff = _mm256_set1_epi32(0xFFFF);
    const __m256i v0x80 = _mm256_set1_epi32(0x80);
    const int *const base = reinterpret_cast<const int *>(probs);
    uint32_t i = 0;
    for (; i + 8 <= N; i += 8) {
        const __m256i idx = _mm256_setr_epi32(i, i + 1, i + 2, i + 3, i + 4, i + 5, i + 6, i + 7);
        const __m256i dw = _mm256_i32gather_epi32(base, idx, 4);
        const __m512d v1 =
            _mm512_div_pd(_mm512_cvtepi32_pd(_mm256_and_si256(dw, v0xffff)), v65535);
        const __m512d v2 =
            _mm512_div_pd(_mm512_cvtepi32_pd(_mm256_srli_epi32(dw, 16)), v65535);
        const __m512d p2 =
            _mm512_max_pd(_mm512_sub_pd(_mm512_sub_pd(vone, v1), v2), vzero);
        __m512d dosage = _mm512_add_pd(v2, _mm512_mul_pd(vtwo, p2));
        if (countFirst) dosage = _mm512_sub_pd(vtwo, dosage);
        dosage = _mm512_min_pd(_mm512_max_pd(dosage, vzero), vtwo);
        const __m128i pl =
            _mm_loadl_epi64(reinterpret_cast<const __m128i *>(ploidy + i));
        const __mmask8 miss =
            _mm256_test_epi32_mask(_mm256_cvtepu8_epi32(pl), v0x80);
        dosage = _mm512_mask_blend_pd(miss, dosage, vnan);
        _mm512_storeu_pd(dst + i, dosage);
    }
    for (; i < N; ++i) dst[i] = decodeDosage16(probs, ploidy, i, countFirst);
}

#endif // x86_64 SIMD variants

// Runtime dispatch: resolve each kernel once at startup to the widest available
// SIMD tier, with a mandatory AVX2 tier and a scalar fallback (ARM / pre-AVX2).
using ParseKernelFn = void (*)(const genfile::byte_t *, const genfile::byte_t *,
                               uint32_t, bool, double *);

inline ParseKernelFn pickParseKernel8() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return parseKernel8_avx512;
    case SimdLevel::AVX2: return parseKernel8_avx2;
    default: break;
    }
#endif
    return parseKernel8_scalar;
}

inline ParseKernelFn pickParseKernel16() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return parseKernel16_avx512;
    case SimdLevel::AVX2: return parseKernel16_avx2;
    default: break;
    }
#endif
    return parseKernel16_scalar;
}

const ParseKernelFn g_parseKernel8 = pickParseKernel8();
const ParseKernelFn g_parseKernel16 = pickParseKernel16();

// P3: specialised parser for the most common BGEN v12 case — diploid,
// biallelic, unphased, 8- or 16-bit probabilities — bypassing genfile's
// per-sample / per-entry setter-callback machinery.  Reproduces genfile +
// DosageSetter bit-for-bit (the kernels above hold the per-sample formula).
//
// ⑤ decode-all-then-compact: step 1 decodes every sample (into out itself when
// allUsed, else into setter.scratchAll); step 2 gathers the used samples in
// ascending file order into out[0..nUsed) and, when setter.needStats, computes
// QC.  dosageSum is summed in that same ascending order, so it is bit-identical
// to the previous single-loop accumulation; getGenotypesSimple sets
// needStats=false and skips the QC pass entirely.  Returns false (caller falls
// back to genfile) if the probability region is shorter than expected.
inline bool fastParseV12DiploidBiallelic(
    const genfile::bgen::v12::GenotypeDataBlock &pack,
    DosageSetter &setter
) {
    const uint32_t N = pack.numberOfSamples;
    const uint32_t stride = (pack.bits == 8) ? 2u : 4u; // bytes/sample (2 probs)
    const genfile::byte_t *probs = pack.buffer;
    if (probs + static_cast<size_t>(N) * stride > pack.end)
        return false; // truncated / unexpected → let genfile handle it

    const genfile::byte_t *ploidy = pack.ploidy;
    const bool countFirst = setter.countFirstAllele;
    double *const out = setter.out;

    // Step 1 — decode every sample into dst (NaN for missing).  allUsed writes
    // straight into out (nUsed == N); the subset case decodes into scratchAll.
    double *const dst = setter.allUsed ? out : setter.scratchAll;
    if (pack.bits == 8) g_parseKernel8(probs, ploidy, N, countFirst, dst);
    else g_parseKernel16(probs, ploidy, N, countFirst, dst);

    // Step 2 — gather used samples (ascending file order) + optional QC.
    std::vector<uint32_t> *const missingIdx = setter.missingIdx;
    uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0, nMissing = 0;
    double dosageSum = 0.0;
    bool anyFractional = false;

    if (setter.allUsed) {
        if (setter.needStats) {
            for (uint32_t i = 0; i < N; ++i) {
                const double d = dst[i]; // == out[i]
                if (std::isnan(d)) {
                    if (missingIdx) missingIdx->push_back(i);
                    ++nMissing;
                    continue;
                }
                dosageSum += d;
                if (d != 0.0 && d != 1.0 && d != 2.0) anyFractional = true;
                const int hc = dosageHardcall(d, setter.hardCallThr);
                if (hc == 0) ++nHomRef;
                else if (hc == 1) ++nHet;
                else if (hc == 2) ++nHomAlt;
                // else: HWE-uncertain — excluded (still summed for AF).
            }
        }
    } else {
        const uint32_t *const ui = setter.usedIndices->data();
        const uint32_t nU = setter.nUsed;
        if (setter.needStats) {
            for (uint32_t k = 0; k < nU; ++k) {
                const double d = dst[ui[k]];
                out[k] = d;
                if (std::isnan(d)) {
                    if (missingIdx) missingIdx->push_back(k);
                    ++nMissing;
                    continue;
                }
                dosageSum += d;
                if (d != 0.0 && d != 1.0 && d != 2.0) anyFractional = true;
                const int hc = dosageHardcall(d, setter.hardCallThr);
                if (hc == 0) ++nHomRef;
                else if (hc == 1) ++nHet;
                else if (hc == 2) ++nHomAlt;
                // else: HWE-uncertain — excluded (still summed for AF).
            }
        } else {
            for (uint32_t k = 0; k < nU; ++k) out[k] = dst[ui[k]];
        }
    }

    setter.nHomRef = nHomRef;
    setter.nHet = nHet;
    setter.nHomAlt = nHomAlt;
    setter.nMissing = nMissing;
    setter.dosageSum = dosageSum;
    setter.anyFractional = anyFractional;
    return true;
}

// Skip the current variant's genotype data block by SEEKING past it, instead of
// streaming the bytes through the buffer as genfile's ignore_genotype_data_block
// does (std::istream::ignore reads/discards every byte).  Mirrors that function's
// framing exactly — for a compressed block the leading 4-byte payload size is
// consumed, then the payload is skipped; for the uncompressed layout-1 case the
// fixed 6*N bytes are skipped — so the stream lands at the identical position and
// the recorded marker metadata is byte-for-byte unchanged.  Used on the metadata
// first pass and the cursor's forward advance, where the genotype bytes are not
// needed; turns an O(file-size) read into O(num-variants) seeks.
inline void skipGenotypeDataBlock(std::istream &s, const genfile::bgen::Context &ctx) {
    if ((ctx.flags & genfile::bgen::e_CompressedSNPBlocks) != genfile::bgen::e_NoCompression) {
        unsigned char b[4];
        s.read(reinterpret_cast<char *>(b), 4);
        const uint32_t sz = static_cast<uint32_t>(b[0]) |
                            (static_cast<uint32_t>(b[1]) << 8) |
                            (static_cast<uint32_t>(b[2]) << 16) |
                            (static_cast<uint32_t>(b[3]) << 24);
        if (sz > 0) s.seekg(static_cast<std::streamoff>(sz), std::ios::cur);
    } else {
        s.seekg(static_cast<std::streamoff>(6) * ctx.number_of_samples, std::ios::cur);
    }
}

// Little-endian 4-byte read from an in-memory buffer (BGEN is little-endian).
inline uint32_t readLE32(const genfile::byte_t *p) {
    return static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8) |
           (static_cast<uint32_t>(p[2]) << 16) | (static_cast<uint32_t>(p[3]) << 24);
}

// ③: in-memory analogue of genfile::bgen::read_snp_identifying_data's framing.
// Advances past the variant's identifying-data block (whose values are already
// captured in the marker metadata, so they are discarded here) and returns a
// pointer to the first byte of the genotype data block.  Mirrors
// read_snp_identifying_data exactly (genfile bgen.hpp:466-519): layout-2 reads
// uint16-length SNPID/RSID/chromosome, uint32 position, uint16 number-of-alleles,
// then per allele a uint32-length string; layout-0/1 prepend a uint32 sample
// count and fix the allele count at 2.  Returns nullptr if the framing would
// overrun [p, end) (treated as a malformed record by the caller).
inline const genfile::byte_t *skipIdentifyingDataInMemory(
    const genfile::byte_t *p,
    const genfile::byte_t *end,
    const genfile::bgen::Context &ctx
) {
    const uint32_t layout = ctx.flags & genfile::bgen::e_Layout;
    auto have = [&](size_t n) { return static_cast<size_t>(end - p) >= n; };
    auto u16 = [&](uint16_t &v) -> bool {
        if (!have(2)) return false;
        v = static_cast<uint16_t>(static_cast<uint16_t>(p[0]) |
                                  (static_cast<uint16_t>(p[1]) << 8));
        p += 2;
        return true;
    };
    auto u32 = [&](uint32_t &v) -> bool {
        if (!have(4)) return false;
        v = readLE32(p);
        p += 4;
        return true;
    };

    uint16_t len16 = 0;
    uint32_t len32 = 0;
    uint16_t nAlleles = 2;
    if (layout == genfile::bgen::e_Layout1 || layout == genfile::bgen::e_Layout0) {
        uint32_t nSamples;
        if (!u32(nSamples)) return nullptr; // number_of_samples (discarded)
    }
    if (!u16(len16) || !have(len16)) return nullptr;
    p += len16; // SNPID
    if (!u16(len16) || !have(len16)) return nullptr;
    p += len16; // RSID
    if (!u16(len16) || !have(len16)) return nullptr;
    p += len16; // chromosome
    if (!have(4)) return nullptr;
    p += 4; // SNP_position
    if (layout == genfile::bgen::e_Layout2) {
        if (!u16(nAlleles)) return nullptr;
    }
    for (uint16_t a = 0; a < nAlleles; ++a) {
        if (!u32(len32) || !have(len32)) return nullptr;
        p += len32; // allele
    }
    return p;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════════════
// BgenData
// ══════════════════════════════════════════════════════════════════════════════

bool BgenData::loadFromIndex(
    const std::string &bgiPath,
    const std::unordered_set<std::string> &chrFilter,
    uint32_t &nSkippedMultiallelic
) {
    sqlite3 *db = nullptr;
    if (sqlite3_open_v2(bgiPath.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) {
        if (db) sqlite3_close(db);
        return false;
    }

    // Staleness guard: bgenix's Metadata table records the indexed .bgen's
    // size; if it disagrees with the actual file the byte offsets cannot be
    // trusted, so fall back to a fresh scan.  Best-effort — a missing column or
    // table simply skips the check.
    {
        sqlite3_stmt *meta = nullptr;
        if (sqlite3_prepare_v2(db, "SELECT file_size FROM Metadata LIMIT 1", -1,
                               &meta, nullptr) == SQLITE_OK &&
            sqlite3_step(meta) == SQLITE_ROW) {
            const uint64_t indexed = static_cast<uint64_t>(sqlite3_column_int64(meta, 0));
            uint64_t actual = 0;
            if (indexed != 0 && statFileSize(m_bgenFile, &actual) && indexed != actual) {
                sqlite3_finalize(meta);
                sqlite3_close(db);
                warnMsg("BGEN: ignoring stale index %s (indexed file_size %llu != "
                        "actual %llu); falling back to scan",
                        bgiPath.c_str(),
                        static_cast<unsigned long long>(indexed),
                        static_cast<unsigned long long>(actual));
                return false;
            }
        }
        sqlite3_finalize(meta);
    }

    sqlite3_stmt *st = nullptr;
    const char *sql =
        "SELECT chromosome, position, rsid, number_of_alleles, allele1, allele2, "
        "file_start_position, size_in_bytes FROM Variant ORDER BY file_start_position";
    if (sqlite3_prepare_v2(db, sql, -1, &st, nullptr) != SQLITE_OK) {
        sqlite3_close(db);
        return false;
    }

    uint64_t ordinal = 0;
    int rc;
    while ((rc = sqlite3_step(st)) == SQLITE_ROW) {
        // Skip non-biallelic to match the scan path (the cursor decodes only
        // the first two alleles).
        if (sqlite3_column_int(st, 3) != 2) {
            ++nSkippedMultiallelic;
            continue;
        }
        const unsigned char *chrTxt = sqlite3_column_text(st, 0);
        const uint32_t position = static_cast<uint32_t>(sqlite3_column_int64(st, 1));
        const unsigned char *rsidTxt = sqlite3_column_text(st, 2);
        const unsigned char *a1 = sqlite3_column_text(st, 4);
        const unsigned char *a2 = sqlite3_column_text(st, 5);
        const uint64_t fpos = static_cast<uint64_t>(sqlite3_column_int64(st, 6));
        const uint64_t rsize = static_cast<uint64_t>(sqlite3_column_int64(st, 7));

        std::string chromosome = chrTxt ? reinterpret_cast<const char *>(chrTxt) : "";
        std::string rsid = rsidTxt ? reinterpret_cast<const char *>(rsidTxt) : "";
        std::string allele1 = a1 ? reinterpret_cast<const char *>(a1) : "";
        std::string allele2 = a2 ? reinterpret_cast<const char *>(a2) : "";
        // .bgi stores allele1 = first allele, allele2 = second allele; map to
        // REF/ALT per the user's --bgen orientation, identical to the scan path.
        const std::string &refAllele = m_altFirst ? allele2 : allele1;
        const std::string &altAllele = m_altFirst ? allele1 : allele2;

        m_chr.push_back(chromosome);
        m_pos.push_back(position);
        m_markerId.push_back(rsid);
        m_ref.push_back(refAllele);
        m_alt.push_back(altAllele);
        m_fileOffset.push_back(fpos);
        m_recordSize.push_back(rsize);
        if (chrFilter.empty() || chrFilter.count(chromosome))
            m_markerInfo.push_back({chromosome, position, rsid, refAllele,
                                    altAllele, ordinal});
        ++ordinal;
    }
    const bool ok = (rc == SQLITE_DONE);
    sqlite3_finalize(st);
    sqlite3_close(db);

    if (!ok) {
        // Mid-iteration error: discard any partial metadata and fall back.
        m_chr.clear();
        m_pos.clear();
        m_markerId.clear();
        m_ref.clear();
        m_alt.clear();
        m_fileOffset.clear();
        m_recordSize.clear();
        m_markerInfo.clear();
        nSkippedMultiallelic = 0;
        return false;
    }
    return true;
}

BgenData::BgenData(
    std::string bgenFile,
    const std::vector<uint64_t> &usedMask,
    uint32_t nSamplesInFile,
    uint32_t nUsed,
    std::unordered_set<std::string> chrFilter,
    std::string extractFile,
    std::string excludeFile,
    int nMarkersEachChunk,
    bool altFirst
)
    : m_bgenFile(std::move(bgenFile)),
      m_nSubjInFile(nSamplesInFile),
      m_nSubjUsed(nUsed),
      m_usedMask(usedMask),
      m_altFirst(altFirst)
{
    m_allUsed = (nUsed == nSamplesInFile);

    // ⑤: precompute the ascending list of used file-sample indices once, so the
    // cursor's parse loop can decode all samples and then gather the used ones,
    // instead of testing the used-mask bit per sample.  Empty when all used.
    if (!m_allUsed) {
        m_usedIndices.reserve(nUsed);
        for (uint32_t i = 0; i < nSamplesInFile; ++i)
            if ((m_usedMask[i >> 6] >> (i & 63)) & 1ULL)
                m_usedIndices.push_back(i);
    }

    // ---- Open and read header ----
    std::ifstream stream(m_bgenFile, std::ios::binary);
    if (!stream.is_open()) throw std::runtime_error("Cannot open BGEN file: " + m_bgenFile);

    genfile::bgen::uint32_t offset;
    genfile::bgen::read_offset(stream, &offset);

    genfile::bgen::Context context;
    genfile::bgen::read_header_block(stream, &context);

    if (context.number_of_samples != m_nSubjInFile)
        throw std::runtime_error("BGEN sample count (" + std::to_string(context.number_of_samples) +
                                 ") does not match expected (" + std::to_string(m_nSubjInFile) + ")");

    m_bgenFlags = context.flags;

    // Skip sample identifiers if present
    if (context.flags & genfile::bgen::e_SampleIdentifiers) {
        genfile::bgen::read_sample_identifier_block(stream, context, [](std::string const &) {
        });
    }

    // Seek to first variant
    stream.seekg(offset + 4);
    m_dataOffset = stream.tellg();

    // ---- Marker metadata: prefer the .bgi (bgenix) index, else linear scan ----
    uint32_t nSkippedMultiallelic = 0;
    const std::string bgiPath = m_bgenFile + ".bgi";
    bool usedIndex = false;
    {
        uint64_t bgiSize = 0;
        if (statFileSize(bgiPath, &bgiSize) && bgiSize > 0)
            usedIndex = loadFromIndex(bgiPath, chrFilter, nSkippedMultiallelic);
    }

    if (usedIndex) {
        infoMsg("  BGEN: read marker metadata from index %s", bgiPath.c_str());
    } else {
        // Linear first pass: read each variant's identifying data and SEEK past
        // its genotype block, recording the byte offset of every biallelic
        // record so the cursor can later seek to it directly.
        std::string SNPID, RSID, chromosome;
        genfile::bgen::uint32_t position;
        std::vector<std::string> alleles;
        uint64_t variantIdx = 0;

        for (;;) {
            const std::streampos varStart = stream.tellg();
            if (!genfile::bgen::read_snp_identifying_data(
                    stream, context, &SNPID, &RSID, &chromosome, &position,
                    [&](std::size_t n) { alleles.resize(n); },
                    [&](std::size_t i, std::string const &a) { alleles[i] = a; }))
                break;

            // Skip non-biallelic: the BGEN cursor decodes only the first two
            // alleles, so retaining a multi-allelic record would silently
            // mis-report carriers of ALT2/ALT3/... .  Count the skips and
            // surface them in a warning after the first pass completes.
            if (alleles.size() != 2) {
                ++nSkippedMultiallelic;
                skipGenotypeDataBlock(stream, context);
                continue;
            }

            m_chr.push_back(chromosome);
            m_pos.push_back(position);
            // Use RSID as marker ID; fall back to SNPID
            m_markerId.push_back(RSID.empty() ? SNPID : RSID);
            // BGEN itself only encodes "first allele" / "second allele".  The
            // user declares the REF/ALT convention at the CLI via
            // --bgen FILE {ref-first|ref-last|ref-unknown}:
            //   ref-first   → alleles[0] = REF (m_altFirst = false; default of
            //                 IMPUTE / qctool / UK Biobank)
            //   ref-last    → alleles[0] = ALT (m_altFirst = true; default of
            //                 plink2 --export bgen-1.x)
            //   ref-unknown → treated as ref-last with a CLI-level warning.
            const std::string &refAllele = m_altFirst ? alleles[1] : alleles[0];
            const std::string &altAllele = m_altFirst ? alleles[0] : alleles[1];
            m_ref.push_back(refAllele);
            m_alt.push_back(altAllele);
            m_fileOffset.push_back(static_cast<uint64_t>(varStart));
            if (chrFilter.empty() || chrFilter.count(chromosome))
                m_markerInfo.push_back({chromosome, position,
                                        RSID.empty() ? SNPID : RSID, refAllele,
                                        altAllele, variantIdx});
            ++variantIdx;

            skipGenotypeDataBlock(stream, context);
            // ③: whole-record size = (post-genotype-block position) - varStart;
            // the seek lands exactly at the next record start (or EOF), so this
            // is id-block + genotype-block bytes, parallel to m_fileOffset.
            m_recordSize.push_back(static_cast<uint64_t>(stream.tellg()) -
                                   static_cast<uint64_t>(varStart));
        }
    }

    // ---- Apply --extract / --exclude filter ----
    geno_factory::filterMarkersByIds(m_markerInfo, extractFile, excludeFile);

    m_nMarkers = static_cast<uint32_t>(m_markerInfo.size());
    m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);

    infoMsg("  BGEN: read %u biallelic variants from %s", m_nMarkers, m_bgenFile.c_str());

    if (nSkippedMultiallelic > 0) {
        warnMsg("BGEN: skipped %u multi-allelic variant(s) in %s; "
                "the BGEN reader handles only biallelic records. "
                "Split them with plink2 --export bgen-1.x 'multiallelics-already-joined=split' "
                "(or an equivalent qctool pipeline) to retain them as separate biallelic records.",
                nSkippedMultiallelic, m_bgenFile.c_str());
    }
}

BgenData::~BgenData() = default;

std::vector<std::vector<uint64_t> > BgenData::buildChunks(
    const std::vector<MarkerInfo> &markers,
    int chunkSize
) {
    std::vector<std::vector<uint64_t> > chunks;
    if (markers.empty()) return chunks;
    std::vector<uint64_t> cur;
    cur.reserve(chunkSize);
    std::string prevChr = markers[0].chrom;
    for (uint32_t i = 0; i < markers.size(); ++i) {
        if (markers[i].chrom != prevChr || static_cast<int>(cur.size()) >= chunkSize) {
            if (!cur.empty()) chunks.push_back(std::move(cur));
            cur.clear();
            cur.reserve(chunkSize);
            prevChr = markers[i].chrom;
        }
        cur.push_back(markers[i].genoIndex);
    }
    if (!cur.empty()) chunks.push_back(std::move(cur));
    return chunks;
}

std::unique_ptr<GenoCursor> BgenData::makeCursor() const {
    return std::make_unique<BgenCursor>(*this);
}

// ══════════════════════════════════════════════════════════════════════════════
// BgenCursor::Impl
// ══════════════════════════════════════════════════════════════════════════════

struct BgenCursor::Impl {
    std::ifstream stream;
    genfile::bgen::Context context;

    // Scratch buffers for decompression
    std::vector<genfile::byte_t> buffer1;
    std::vector<genfile::byte_t> buffer2;

    // Sample subsetting
    const std::vector<uint64_t> *usedMask = nullptr;
    uint32_t nSamplesInFile = 0;
    uint32_t nUsed = 0;
    bool allUsed = true;
    bool altFirst = false;

    // ⑤ — ascending used-index list (borrowed from parent) + scratch holding all
    // N decoded dosages before compaction (sized to nSamplesInFile when subset).
    const std::vector<uint32_t> *usedIndices = nullptr;
    std::vector<double> scratchAll;

    // ④ — capped per-cursor read-ahead window holding a contiguous file span.
    // Records within a chunk are decoded in ascending, contiguous file order, so
    // one window fill serves many subsequent records, coalescing O(records)
    // small reads into O(span / WINDOW_CAP) large sequential reads.
    static constexpr size_t WINDOW_CAP = size_t(8) << 20; // 8 MB
    std::vector<genfile::byte_t> window;
    uint64_t winFileStart = 0; // file offset of window[0]
    size_t winValid = 0;       // bytes currently valid in window
    uint64_t fileSize = 0;     // cached .bgen size, to clamp reads at EOF

    // Return a pointer to the `recSize` whole-record bytes for the variant whose
    // record starts at byte `recStart`.  Served from the window when fully
    // covered; otherwise refill with a capped span starting at recStart (grown
    // to recSize for a rare oversized record), replacing the old per-variant
    // seek + read_snp_identifying_data re-parse with one read per ~WINDOW_CAP.
    // Correct for arbitrary recStart (absolute offsets), so out-of-order or
    // non-primed callers still get the right bytes — only coalescing is lost.
    const genfile::byte_t *getRecordBytes(uint64_t recStart, uint64_t recSize) {
        if (recStart >= winFileStart &&
            recStart + recSize <= winFileStart + winValid)
            return window.data() + (recStart - winFileStart);

        size_t want = WINDOW_CAP;
        if (recSize > want) want = recSize; // ensure this record fits
        if (fileSize && recStart + want > fileSize)
            want = static_cast<size_t>(fileSize - recStart); // clamp at EOF
        if (want < recSize)
            throw std::runtime_error("BGEN: variant record at offset " +
                                     std::to_string(recStart) + " exceeds file size");
        if (window.size() < want) window.resize(want);

        stream.clear(); // drop any eof/fail set by a previous block
        stream.seekg(static_cast<std::streamoff>(recStart));
        stream.read(reinterpret_cast<char *>(window.data()),
                    static_cast<std::streamsize>(want));
        winFileStart = recStart;
        winValid = static_cast<size_t>(stream.gcount());
        if (winValid < recSize)
            throw std::runtime_error("BGEN: short read for variant record at offset " +
                                     std::to_string(recStart));
        return window.data();
    }

    // Per-cursor (per-thread) libdeflate zlib decompressor; allocated lazily.
    libdeflate_decompressor *deflate = nullptr;

    // P2: decompress the layout-2 genotype block [compData, compData+compLen)
    // into buffer2 using libdeflate (zlib) or the vendored zstd, or copy when
    // uncompressed.  The decompressed bytes are uniquely defined by the input,
    // so the result is identical to genfile's zlib/zstd uncompress.  Returns
    // false (leaving the caller to use genfile's uncompress) for any codec we do
    // not specialise.  compData points into the in-memory record buffer (③).
    bool decompressLayout2(const genfile::byte_t *compData, size_t compLen) {
        const uint32_t comp = context.flags & genfile::bgen::e_CompressedSNPBlocks;
        if (comp == genfile::bgen::e_NoCompression) {
            buffer2.assign(compData, compData + compLen);
            return true;
        }
        if (compLen < 4) return false;
        const unsigned char *b = reinterpret_cast<const unsigned char *>(compData);
        const uint32_t usize = static_cast<uint32_t>(b[0]) |
                               (static_cast<uint32_t>(b[1]) << 8) |
                               (static_cast<uint32_t>(b[2]) << 16) |
                               (static_cast<uint32_t>(b[3]) << 24);
        buffer2.resize(usize);
        const void *in = compData + 4;
        const size_t inLen = compLen - 4;
        if (comp == genfile::bgen::e_ZlibCompression) {
            if (!deflate) deflate = libdeflate_alloc_decompressor();
            if (!deflate || usize == 0) return false;
            size_t got = 0;
            const libdeflate_result r = libdeflate_zlib_decompress(
                deflate, in, inLen, buffer2.data(), usize, &got);
            return (r == LIBDEFLATE_SUCCESS && got == usize);
        }
        if (comp == genfile::bgen::e_ZstdCompression) {
            if (usize == 0) return false;
            const size_t got = ZSTD_decompress(buffer2.data(), usize, in, inLen);
            return (!ZSTD_isError(got) && got == usize);
        }
        return false;
    }

    // Read + decompress + parse the variant's genotype block into the setter.
    // ③: the whole record is read once into recordBuf (or served from the ④
    // window); the identifying-data block is skipped in memory and the genotype
    // block located there, then routed through libdeflate (P2) + the specialised
    // parser (P3), falling back to genfile's uncompress/parse for any other
    // layout, codec, ploidy, phasing, or bit depth — byte-identical in every
    // case.  recStart / recSize come from BgenData::fileOffset / recordSize.
    void readDecompressParse(DosageSetter &setter, uint64_t recStart, uint64_t recSize) {
        const genfile::byte_t *rec = getRecordBytes(recStart, recSize);
        const genfile::byte_t *recEnd = rec + recSize;

        // ③: skip the identifying-data block in memory (its values are already
        // in the marker metadata) to land on the genotype data block.
        const genfile::byte_t *gp = skipIdentifyingDataInMemory(rec, recEnd, context);
        if (!gp)
            throw std::runtime_error("BGEN: malformed variant record at offset " +
                                     std::to_string(recStart));

        // Genotype-block framing, mirroring read_genotype_data_block: a layout-2
        // block or any compressed block is prefixed by a 4-byte payload size; an
        // uncompressed layout-0/1 block is 6*N bytes with no prefix.
        const uint32_t layout = context.flags & genfile::bgen::e_Layout;
        const bool compressed = (context.flags & genfile::bgen::e_CompressedSNPBlocks) !=
                                genfile::bgen::e_NoCompression;
        uint32_t payloadSize;
        const genfile::byte_t *dataStart;
        if (layout == genfile::bgen::e_Layout2 || compressed) {
            if (gp + 4 > recEnd)
                throw std::runtime_error("BGEN: truncated genotype block at offset " +
                                         std::to_string(recStart));
            payloadSize = readLE32(gp);
            dataStart = gp + 4;
        } else {
            payloadSize = 6u * context.number_of_samples;
            dataStart = gp;
        }
        if (dataStart + payloadSize > recEnd)
            throw std::runtime_error("BGEN: genotype block overruns record at offset " +
                                     std::to_string(recStart));

        bool haveBuffer2 = false;
        if (layout == genfile::bgen::e_Layout2)
            haveBuffer2 = decompressLayout2(dataStart, payloadSize); // P2
        if (!haveBuffer2) {
            buffer1.assign(dataStart, dataStart + payloadSize);
            genfile::bgen::uncompress_probability_data(context, buffer1, &buffer2);
        }

        bool parsed = false;
        if (layout == genfile::bgen::e_Layout2) {
            genfile::bgen::v12::GenotypeDataBlock pack(
                context, buffer2.data(), buffer2.data() + buffer2.size());
            if (pack.numberOfAlleles == 2 && pack.ploidyExtent[0] == 2 &&
                pack.ploidyExtent[1] == 2 && !pack.phased &&
                (pack.bits == 8 || pack.bits == 16))
                parsed = fastParseV12DiploidBiallelic(pack, setter); // P3
        }
        if (!parsed)
            genfile::bgen::parse_probability_data(
                buffer2.data(), buffer2.data() + buffer2.size(), context, setter);
    }

    ~Impl() {
        if (deflate) libdeflate_free_decompressor(deflate);
    }
};

BgenCursor::BgenCursor(const BgenData &parent)
    : m_parent(parent),
      m_impl(std::make_unique<Impl>())
{
    auto &impl = *m_impl;
    impl.nSamplesInFile = parent.nSubjInFile();
    impl.nUsed = parent.nSubjUsed();
    impl.allUsed = parent.allUsed();
    impl.usedMask = &parent.usedMask();
    impl.altFirst = parent.altFirst();
    impl.usedIndices = &parent.usedIndices();
    if (!impl.allUsed) impl.scratchAll.resize(impl.nSamplesInFile);

    // Open file and read header to build context
    impl.stream.open(parent.bgenFile(), std::ios::binary);
    if (!impl.stream.is_open()) throw std::runtime_error("BGEN cursor: cannot open " + parent.bgenFile());

    // ④: cache the file size so the read-ahead window can clamp its span at EOF.
    { uint64_t fsz = 0; if (statFileSize(parent.bgenFile(), &fsz)) impl.fileSize = fsz; }

    genfile::bgen::uint32_t offset;
    genfile::bgen::read_offset(impl.stream, &offset);
    (void)offset; // consumed only to position the stream for the header block
    genfile::bgen::read_header_block(impl.stream, &impl.context);

    // Skip sample identifiers
    if (impl.context.flags & genfile::bgen::e_SampleIdentifiers) {
        genfile::bgen::read_sample_identifier_block(impl.stream, impl.context, [](std::string const &) {
        });
    }
    // No seek to the first variant: the cursor seeks directly to each marker's
    // byte offset (BgenData::fileOffset) on every getGenotypes call.
}

BgenCursor::~BgenCursor() = default;

void BgenCursor::beginSequentialBlock(uint64_t /*firstMarker*/) {
    // ④: invalidate the read-ahead window at each chunk boundary (a worker may
    // be handed a non-adjacent chunk) and clear any residual eof/fail state so
    // the next refill's seekg succeeds.  getRecordBytes then fills the window
    // from the chunk's first record.
    m_impl->stream.clear();
    m_impl->winValid = 0;
}

GenoStats BgenCursor::getGenotypes(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out,
    std::vector<uint32_t> &indexForMissing
) {
    auto &impl = *m_impl;
    indexForMissing.clear();

    DosageSetter setter;
    setter.out = out.data();
    setter.missingIdx = &indexForMissing;
    setter.usedMask = impl.usedMask;
    setter.nSamplesInFile = impl.nSamplesInFile;
    setter.allUsed = impl.allUsed;
    setter.countFirstAllele = impl.altFirst;
    setter.usedIndices = impl.usedIndices;
    setter.scratchAll = impl.scratchAll.empty() ? nullptr : impl.scratchAll.data();
    setter.nUsed = impl.nUsed;
    setter.needStats = true; // QC stats required
    setter.hardCallThr = m_parent.hardCallThreshold;

    impl.readDecompressParse(setter, m_parent.fileOffset(gIndex),
                             m_parent.recordSize(gIndex));

    const uint32_t nUsed = impl.nUsed;
    // statsFromCounts treats its first count argument as "hom of the ALT
    // allele" and returns altCounts = 2*nHomAlt + nHet.  Our DosageSetter
    // already classified the per-sample dosage of the ALT allele
    // (alleles[1] under --bgen FILE ref-first, alleles[0] under ref-last
    // or ref-unknown), so setter.nHomAlt is the count of subjects with
    // dosage ≈ 2 (hom ALT).
    GenoStats gs = statsFromCounts(setter.nHomAlt, setter.nHet, setter.nHomRef, setter.nMissing, nUsed);
    gs.fromDosage = setter.anyFractional;

    // Dosage-aware AF: compute the allele frequency from the dosage sum rather
    // than the thresholded hard-call counts.  For hard calls the dosage sum
    // equals 2*nHomAlt + nHet exactly, so altFreq/altCounts/mac are unchanged;
    // for true dosages this is the correct (un-binned) frequency.  log10pHwe retains
    // the count-based value from statsFromCounts, where the counts are the
    // plink2 hard-call-threshold classification (dosageHardcall) — matching
    // plink2 --hardy; missingRate stays NaN-only.
    const uint32_t nNonMissing = nUsed - setter.nMissing;
    if (nNonMissing > 0) {
        const AlleleFreq af = alleleFreqFromTotal(setter.dosageSum, nNonMissing);
        gs.altCounts = setter.dosageSum;
        gs.altFreq = af.altFreq;
        gs.maf = af.maf;
        gs.mac = af.mac;
    }
    return gs;
}

void BgenCursor::getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out
) {
    auto &impl = *m_impl;

    DosageSetter setter;
    setter.out = out.data();
    setter.missingIdx = nullptr;
    setter.usedMask = impl.usedMask;
    setter.nSamplesInFile = impl.nSamplesInFile;
    setter.allUsed = impl.allUsed;
    setter.countFirstAllele = impl.altFirst;
    setter.usedIndices = impl.usedIndices;
    setter.scratchAll = impl.scratchAll.empty() ? nullptr : impl.scratchAll.data();
    setter.nUsed = impl.nUsed;
    setter.needStats = false; // fill out only; no QC

    impl.readDecompressParse(setter, m_parent.fileOffset(gIndex),
                             m_parent.recordSize(gIndex));
}
