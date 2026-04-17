// marker_impl.hpp — Shared helpers for marker engines (formatting, stats)
//
// Used by marker.cpp (markerEngine, multiPhenoEngine) and loco.cpp (locoEngine).
// Not part of the public API — include only from engine translation units.
#pragma once

#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "geno_factory/hwe.hpp"
#include "util/simd_dispatch.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <immintrin.h>
#include <string>
#include <string_view>
#include <vector>

// Forward-declare plink2's fast double-to-string (6 sig figs, ~5-10× faster
// than snprintf("%.6g")).  Already compiled and linked via pgenlib.
namespace plink2 {
char *dtoa_g(
    double dxx,
    char *start
);

} // namespace plink2

namespace engine_impl {

// ──────────────────────────────────────────────────────────────────────
// TSV formatting helpers — zero allocation in the hot path
// ──────────────────────────────────────────────────────────────────────

// Build a suffix of nCols tab-separated "NA" values for fail-QC rows.
inline std::string makeNaSuffix(int nResultCols) {
    if (nResultCols <= 0) return {};
    std::string s;
    s.reserve(3 * static_cast<size_t>(nResultCols));
    for (int i = 0; i < nResultCols; ++i) {
        s += '\t';
        s += 'N';
        s += 'A';
    }
    return s;
}

// Format double: "NA" | "Inf" | "-Inf" | 6-sig-fig decimal/scientific.
// Returns char count.  Uses plink2's dtoa_g for ~5-10× speedup over snprintf.
inline int numToChars(
    char *buf,
    double x
) {
    if (std::isnan(x)) {
        buf[0] = 'N';
        buf[1] = 'A';
        return 2;
    }
    if (std::isinf(x)) {
        if (x > 0) {
            buf[0] = 'I';
            buf[1] = 'n';
            buf[2] = 'f';
            return 3;
        } else {
            buf[0] = '-';
            buf[1] = 'I';
            buf[2] = 'n';
            buf[3] = 'f';
            return 4;
        }
    }
    char *end = plink2::dtoa_g(x, buf);
    return static_cast<int>(end - buf);
}

// Append 9 meta columns: CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P
inline void appendMeta(
    std::string &out,
    char *buf,
    std::string_view chrom,
    uint32_t pos,
    std::string_view id,
    std::string_view ref,
    std::string_view alt,
    double missRate,
    double altFreq,
    double mac,
    double hweP
) {
    int n;
    out += chrom;
    out += '\t';
    n = std::snprintf(buf, 32, "%u", pos);
    out.append(buf, n);
    out += '\t';
    out += id;
    out += '\t';
    out += ref;
    out += '\t';
    out += alt;
    out += '\t';
    n = numToChars(buf, missRate);
    out.append(buf, n);
    out += '\t';
    n = numToChars(buf, altFreq);
    out.append(buf, n);
    out += '\t';
    n = numToChars(buf, mac);
    out.append(buf, n);
    out += '\t';
    n = numToChars(buf, hweP);
    out.append(buf, n);
}

// Full result line: meta + tab-separated doubles.
inline void formatLine(
    std::string &out,
    char *buf,
    std::string_view chrom,
    uint32_t pos,
    std::string_view id,
    std::string_view ref,
    std::string_view alt,
    double missRate,
    double altFreq,
    double mac,
    double hweP,
    const std::vector<double> &vals
) {
    appendMeta(out, buf, chrom, pos, id, ref, alt, missRate, altFreq, mac, hweP);
    for (double v : vals) {
        out += '\t';
        int n = numToChars(buf, v);
        out.append(buf, n);
    }
    out += '\n';
}

// Fail-QC line: meta + precomputed NA suffix.
inline void formatLineNA(
    std::string &out,
    char *buf,
    std::string_view chrom,
    uint32_t pos,
    std::string_view id,
    std::string_view ref,
    std::string_view alt,
    double missRate,
    double altFreq,
    double mac,
    double hweP,
    const std::string &naSuffix
) {
    appendMeta(out, buf, chrom, pos, id, ref, alt, missRate, altFreq, mac, hweP);
    out += naSuffix;
    out += '\n';
}

// ──────────────────────────────────────────────────────────────────────
// Output header
// ──────────────────────────────────────────────────────────────────────

constexpr const char *META_HEADER = "CHROM\tPOS\tID\tREF\tALT\tMISS_RATE\tALT_FREQ\tMAC\tHWE_P";

inline std::string buildHeader(const MethodBase &method) {
    return std::string(META_HEADER) + method.getHeaderColumns();
}

// ──────────────────────────────────────────────────────────────────────
// Padded flag: one per chunk, prevents false sharing between workers
// ──────────────────────────────────────────────────────────────────────

struct alignas(64) PaddedFlag {
    char ready;
};

static_assert(sizeof(PaddedFlag) == 64, "PaddedFlag must be 64 bytes");

// ──────────────────────────────────────────────────────────────────────
// Per-phenotype genotype stats and extraction
// ──────────────────────────────────────────────────────────────────────

struct PhenoGenoStats {
    double altFreq, mac, missingRate, hweP;
};

// Compute per-phenotype genotype stats from union-level genotype vector
// and a mask column (1.0 = present, 0.0 = absent).
// Used by the fused GEMM path where per-phenotype extraction is skipped.
inline PhenoGenoStats statsFromUnionVec(
    const double *unionG,
    const double *mask,
    uint32_t nUnion,
    uint32_t nUsed
) {
    uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0, nMissing = 0;
    for (uint32_t i = 0; i < nUnion; ++i) {
        if (mask[i] == 0.0) continue;  // absent from this phenotype
        const double v = unionG[i];
        if (v == 0.0)
            ++nHomRef;
        else if (v == 1.0)
            ++nHet;
        else if (v == 2.0)
            ++nHomAlt;
        else
            ++nMissing;  // NaN or dosage
    }
    PhenoGenoStats s;
    uint32_t nonMissing = nHomRef + nHet + nHomAlt;
    if (nonMissing == 0) {
        // Dosage data: compute from gSum (not discrete genotypes).
        // No HWE for dosages → set to 1.0, compute AF from gSum.
        s.missingRate = 0.0;
        s.hweP = 1.0;
        // altFreq and mac will be overwritten by caller from gSum.
        s.altFreq = 0.0;
        s.mac = 0.0;
        return s;
    }
    uint32_t altCounts = 2 * nHomAlt + nHet;
    double n2 = 2.0 * nonMissing;
    s.altFreq = altCounts / n2;
    double maf = std::min(s.altFreq, 1.0 - s.altFreq);
    s.mac = maf * n2;
    s.missingRate = static_cast<double>(nMissing) / static_cast<double>(nUsed);
    s.hweP = HweExact(nHet, nHomAlt, nHomRef);
    return s;
}

inline PhenoGenoStats statsFromGVec(
    const double *g,
    uint32_t n,
    std::vector<uint32_t> &indexForMissing
) {
    indexForMissing.clear();
    uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0;
    for (uint32_t i = 0; i < n; ++i) {
        const double v = g[i];
        if (std::isnan(v) || v < 0.0) {
            indexForMissing.push_back(i);
            continue;
        }
        if (v == 0.0)
            ++nHomRef;
        else if (v == 1.0)
            ++nHet;
        else if (v == 2.0)
            ++nHomAlt;
        else
            indexForMissing.push_back(i);
    }
    PhenoGenoStats s;
    uint32_t nonMissing = nHomRef + nHet + nHomAlt;
    if (nonMissing == 0) {
        s.altFreq = NAN;
        s.mac = NAN;
        s.missingRate = 1.0;
        s.hweP = NAN;
        return s;
    }
    uint32_t altCounts = 2 * nHomAlt + nHet;
    double n2 = 2.0 * nonMissing;
    s.altFreq = altCounts / n2;
    double maf = std::min(s.altFreq, 1.0 - s.altFreq);
    s.mac = maf * n2;
    s.missingRate = static_cast<double>(indexForMissing.size()) / n;
    s.hweP = HweExact(nHet, nHomAlt, nHomRef);
    return s;
}

// Extract per-phenotype genotype vector from union vector.
inline void extractPhenoGVec(
    const double *unionG,
    uint32_t nUnion,
    const uint32_t *unionToLocal,
    uint32_t /*nPheno*/,
    double *phenoG
) {
    for (uint32_t i = 0; i < nUnion; ++i) {
        uint32_t li = unionToLocal[i];
        if (li != UINT32_MAX) phenoG[li] = unionG[i];
    }
}

// ──────────────────────────────────────────────────────────────────────
// Precomputed gather-index extraction: O(nUsed) instead of O(nUnion).
//
// presentIndices[j] = union position of the j-th present subject.
// Because unionToLocal is monotonically increasing where != UINT32_MAX,
// phenoG[j] = unionG[presentIndices[j]] — sequential write, gathered read.
// ──────────────────────────────────────────────────────────────────────

// Build the gather-index list from unionToLocal. Called once per MissBatch.
inline std::vector<uint32_t> buildPresentIndices(
    const uint32_t *unionToLocal,
    uint32_t nUnion
) {
    std::vector<uint32_t> idx;
    idx.reserve(nUnion / 2);
    for (uint32_t i = 0; i < nUnion; ++i)
        if (unionToLocal[i] != UINT32_MAX)
            idx.push_back(i);
    return idx;
}

// Check if unionToLocal is identity (nUsed == nUnion, no gaps).
inline bool isIdentityMapping(
    const uint32_t *unionToLocal,
    uint32_t nUnion,
    uint32_t nUsed
) {
    return nUsed == nUnion;
}

// Gather-index extraction: branchless, sequential output.
inline void extractGather(
    const double *unionG,
    const uint32_t *presentIndices,
    uint32_t nUsed,
    double *phenoG
) {
    for (uint32_t j = 0; j < nUsed; ++j)
        phenoG[j] = unionG[presentIndices[j]];
}

// AVX-512 compress-store extraction using precomputed bitmask.
// presentMask is a packed bitmask: bit i set ↔ union subject i is present.
// Processes 8 doubles per iteration using VCOMPRESSSTOREPD.
__attribute__((target("avx512f,avx512vl")))
inline void extractCompress_avx512(
    const double *unionG,
    const uint64_t *presentMask,
    uint32_t nUnion,
    double *phenoG
) {
    uint32_t outIdx = 0;
    const uint32_t nUnion8 = nUnion & ~uint32_t(7);
    for (uint32_t i = 0; i < nUnion8; i += 8) {
        __m512d vals = _mm512_loadu_pd(&unionG[i]);
        __mmask8 mask = static_cast<__mmask8>(
            (presentMask[i >> 6] >> (i & 63)) & 0xFF);
        _mm512_mask_compressstoreu_pd(&phenoG[outIdx], mask, vals);
        outIdx += static_cast<uint32_t>(_mm_popcnt_u32(mask));
    }
    // Scalar tail
    for (uint32_t i = nUnion8; i < nUnion; ++i) {
        if ((presentMask[i >> 6] >> (i & 63)) & 1)
            phenoG[outIdx++] = unionG[i];
    }
}

// Build bitmask from unionToLocal. Called once per MissBatch.
inline std::vector<uint64_t> buildPresentMask(
    const uint32_t *unionToLocal,
    uint32_t nUnion
) {
    const size_t nWords = (static_cast<size_t>(nUnion) + 63) / 64;
    std::vector<uint64_t> mask(nWords, 0);
    for (uint32_t i = 0; i < nUnion; ++i)
        if (unionToLocal[i] != UINT32_MAX)
            mask[i >> 6] |= uint64_t(1) << (i & 63);
    return mask;
}

// Runtime-dispatched fast extraction.
inline void extractPhenoFast(
    const double *unionG,
    uint32_t nUnion,
    uint32_t nUsed,
    bool identity,
    const uint32_t *presentIndices,
    const uint64_t *presentMask,
    double *phenoG
) {
    if (identity) {
        std::memcpy(phenoG, unionG, static_cast<size_t>(nUsed) * sizeof(double));
        return;
    }
    if (simdLevel() >= SimdLevel::AVX512) {
        extractCompress_avx512(unionG, presentMask, nUnion, phenoG);
    } else {
        extractGather(unionG, presentIndices, nUsed, phenoG);
    }
}

// ──────────────────────────────────────────────────────────────────────
// Fused extract + stats for K phenotypes in a single pass over unionG.
//
// Reads GVec_union once instead of K times, accumulating per-phenotype
// genotype counts and missing indices inline.  Eliminates K separate
// statsFromGVec passes over the per-phenotype buffers.
//
// Memory-bandwidth savings per marker: ~(2K-2) × nUnion × 8 bytes
// (K-1 extra reads of GVec_union + K reads of GVec_pheno for stats).
// ──────────────────────────────────────────────────────────────────────

inline void extractAndStatsBatched(
    const double *unionG,
    uint32_t nUnion,
    size_t K,
    const uint32_t *const *unionToLocalPtrs,
    const uint32_t *nPhenoArr,
    double *const *phenoGPtrs,
    PhenoGenoStats *outStats,
    std::vector<uint32_t> *outMissings
) {
    // Per-phenotype accumulator
    struct Acc {
        uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0;
    };

    // Use small stack buffer for typical K; fall back to heap for large K
    constexpr size_t kStackMax = 64;
    Acc stackAcc[kStackMax];
    Acc *acc = (K <= kStackMax) ? stackAcc : new Acc[K]{};
    if (K <= kStackMax) {
        for (size_t p = 0; p < K; ++p)
            acc[p] = {};
    }
    for (size_t p = 0; p < K; ++p)
        outMissings[p].clear();

    // ── Single pass over union genotypes ────────────────────────────
    for (uint32_t i = 0; i < nUnion; ++i) {
        const double val = unionG[i];
        for (size_t p = 0; p < K; ++p) {
            const uint32_t li = unionToLocalPtrs[p][i];
            if (li == UINT32_MAX) continue;
            phenoGPtrs[p][li] = val;
            if (val == 0.0) {
                ++acc[p].nHomRef;
            } else if (val == 1.0) {
                ++acc[p].nHet;
            } else if (val == 2.0) {
                ++acc[p].nHomAlt;
            } else {
                outMissings[p].push_back(li);
            }
        }
    }

    // ── Compute final stats ─────────────────────────────────────────
    for (size_t p = 0; p < K; ++p) {
        const uint32_t nonMissing = acc[p].nHomRef + acc[p].nHet + acc[p].nHomAlt;
        PhenoGenoStats &s = outStats[p];
        if (nonMissing == 0) {
            s.altFreq = NAN;
            s.mac = NAN;
            s.missingRate = 1.0;
            s.hweP = NAN;
        } else {
            const uint32_t nP = nPhenoArr[p];
            const uint32_t altCounts = 2 * acc[p].nHomAlt + acc[p].nHet;
            const double n2 = 2.0 * nonMissing;
            s.altFreq = altCounts / n2;
            const double maf = std::min(s.altFreq, 1.0 - s.altFreq);
            s.mac = maf * n2;
            s.missingRate = static_cast<double>(outMissings[p].size()) / nP;
            s.hweP = HweExact(acc[p].nHet, acc[p].nHomAlt, acc[p].nHomRef);
        }
    }

    if (K > kStackMax) delete[] acc;
}

// ──────────────────────────────────────────────────────────────────────
// Sparse genotype: fill + scatter + stats for one phenotype batch
//
// Given a difflist (sparse representation from pgenlib), fills the phenotype
// vector with the common genotype, scatters difflist entries, and computes
// per-phenotype stats — all in O(diffLen) rather than O(nPheno).
// ──────────────────────────────────────────────────────────────────────

inline PhenoGenoStats sparseExtractAndStats(
    double *phenoG,
    uint32_t nPheno,
    const uint32_t *unionToLocal,
    uint32_t commonGeno,
    const uint32_t *diffSampleIds,
    const uint8_t *diffGenoCodes,
    uint32_t diffLen,
    std::vector<uint32_t> &indexForMissing
) {
    indexForMissing.clear();

    // Fill with common genotype (sequential write, very cache-friendly)
    const double commonD = static_cast<double>(commonGeno);
    std::fill(phenoG, phenoG + nPheno, commonD);

    // Count genotypes from difflist entries within this phenotype
    uint32_t counts[4] = {0, 0, 0, 0}; // [0]=hom_ref, [1]=het, [2]=hom_alt, [3]=missing
    for (uint32_t j = 0; j < diffLen; ++j) {
        const uint32_t li = unionToLocal[diffSampleIds[j]];
        if (li == UINT32_MAX) continue;
        const uint8_t gc = diffGenoCodes[j];
        if (gc == 3) {
            phenoG[li] = std::numeric_limits<double>::quiet_NaN();
            indexForMissing.push_back(li);
        } else {
            phenoG[li] = static_cast<double>(gc);
        }
        counts[gc]++;
    }

    // All remaining samples have the common genotype
    const uint32_t nDiffInPheno = counts[0] + counts[1] + counts[2] + counts[3];
    counts[commonGeno] += nPheno - nDiffInPheno;

    // Compute stats from counts
    PhenoGenoStats s;
    const uint32_t nonMissing = counts[0] + counts[1] + counts[2];
    if (nonMissing == 0) {
        s.altFreq = NAN;
        s.mac = NAN;
        s.missingRate = 1.0;
        s.hweP = NAN;
        return s;
    }
    const uint32_t altCounts = 2 * counts[2] + counts[1];
    const double n2 = 2.0 * nonMissing;
    s.altFreq = altCounts / n2;
    const double maf = std::min(s.altFreq, 1.0 - s.altFreq);
    s.mac = maf * n2;
    s.missingRate = static_cast<double>(counts[3]) / nPheno;
    s.hweP = HweExact(counts[1], counts[2], counts[0]);
    return s;
}

// ──────────────────────────────────────────────────────────────────────
// Per-worker thread context (used by markerEngine)
// ──────────────────────────────────────────────────────────────────────

struct ThreadContext {
    std::unique_ptr<MethodBase> method; // cloned per thread
    std::unique_ptr<GenoCursor> cursor; // per-thread genotype decoder
    std::string naSuffix;

    ThreadContext(
        const MethodBase &proto,
        const GenoMeta &gd
    )
        : method(proto.clone()),
          cursor(gd.makeCursor()),
          naSuffix(makeNaSuffix(proto.resultSize()))
    {
    }

};

} // namespace engine_impl
