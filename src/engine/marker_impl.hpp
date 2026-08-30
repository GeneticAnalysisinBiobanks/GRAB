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
#if defined(__x86_64__) || defined(_M_X64)
#  include <immintrin.h>
#endif
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

// Append 9 meta columns: CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC LOG10P_HWE
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
    double log10pHwe
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
    n = numToChars(buf, log10pHwe);
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
    double log10pHwe,
    const std::vector<double> &vals
) {
    appendMeta(out, buf, chrom, pos, id, ref, alt, missRate, altFreq, mac, log10pHwe);
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
    double log10pHwe,
    const std::string &naSuffix
) {
    appendMeta(out, buf, chrom, pos, id, ref, alt, missRate, altFreq, mac, log10pHwe);
    out += naSuffix;
    out += '\n';
}

// ──────────────────────────────────────────────────────────────────────
// Output header
// ──────────────────────────────────────────────────────────────────────

// LOG10P_HWE is -log10 of the exact HWE p-value (decision D7 of
// dev-notes/methods/log10p_unify/00_decisions.md): 0 means p = 1, larger means
// stronger departure from Hardy-Weinberg proportions, and the column stays a
// magnitude rather than saturating at 0 the way the former linear HWE_P did
// below p ~ 1e-300.  NA when no subject has a hard call.
constexpr const char *META_HEADER = "CHROM\tPOS\tID\tREF\tALT\tMISS_RATE\tALT_FREQ\tMAC\tLOG10P_HWE";

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

// --hwe is an *input* threshold and stays on the linear p scale (decision D8
// of dev-notes/methods/log10p_unify/00_decisions.md): 1e-6 is the habitual
// spelling and 6 is not.  The QC statistic is now -log10(p), so the comparison
// is carried out against -log10(cutoff) and its direction is reversed:
// "p below the cutoff" becomes "-log10(p) above -log10(cutoff)".  A marker
// with no hard-called subject has log10pHwe = NaN, and NaN fails both the old
// and the new comparison, so such a marker is still not rejected on HWE.
// Returns 0 when the filter is disabled (cutoff <= 0), which no caller reads.
inline double hweLog10CutoffOf(double hweCutoff) {
    return (hweCutoff > 0.0) ? -std::log10(hweCutoff) : 0.0;
}

struct PhenoGenoStats {
    double altFreq, mac, missingRate, log10pHwe;
    // Sum of squared post-impute genotypes over mask-included subjects.
    // Used by SPAsqr to compute empirical Var(G) = (sumSq − sum²/n) / (n−1).
    double sumSq;
    // Sum of post-impute genotypes over mask-included subjects (= mask^T·G).
    // Replaces the GEMM mask column on the fused path; consumed as the per-
    // phenotype genotype sum and (for dosage markers) to recompute altFreq/mac.
    double gSum;
    // True ⇒ this is a dosage marker: altFreq/mac were left as sentinels and
    // the caller must recompute them from the GEMM gSum (see
    // statsFromUnionVecMultiGroup).
    bool fromDosage;
};

// Single-pass multi-group QC.  The fused path previously called
// statsFromUnionVec once per (QC group × marker), so the same N_union-length
// genotype column was re-streamed and the per-subject 0/1/2 classification was
// re-run for every group (typically dozens at biobank scale).  That re-scan
// dominated the non-GEMM runtime: it scaled genotype DRAM / dTLB traffic and
// the mispredicted classification branch by a factor of nGroups.
//
// This routine reads each column ONCE and classifies each subject ONCE, then
// fans the contribution into every group's accumulators through that group's
// 1-bit subject mask.  Genotype memory traffic and classification branches
// drop by a factor of nGroups; the per-subject inner loop touches only the
// small L1-resident accumulator arrays and the compact per-group mask words.
//
// For hard-call data with no missingness this is byte-identical to the prior
// per-group scan (no tolerance); dosage and genuine-missing subjects now follow
// the plink2 hard-call threshold (see the classification below).  Ordering:
//   - gSum and sumSq are accumulated per group in increasing subject index i
//     (i is the outer loop), matching the original per-group scan order.  For
//     excluded subjects the branchless mask multiply adds v·0 = +0.0; post-
//     impute genotypes are non-negative, so each group's running sum is ≥ 0
//     and "x + 0.0 == x" holds bit-exactly — the masked partial sums are thus
//     identical to skipping the excluded subjects outright.
//   - the integer class counts (and hweNegLog10P, which depends only on them) are
//     order-independent.
//
// out[g * outStride] receives group g's stats; isDosage is the per-marker flag
// (identical for all groups of one marker).  accGSum/accSumSq (length nGroups)
// and accCnt (length nGroups*4: class counts 0..3) are caller-owned scratch,
// reused across markers; this routine zeroes them on entry.
inline void statsFromUnionVecMultiGroup(
    const double *unionG,
    const uint64_t *const *groupMaskBits,
    const uint32_t *groupNUsed,
    uint32_t nGroups,
    uint32_t nUnion,
    bool isDosage,
    double thr,
    const uint64_t *origMissBits,
    double *accGSum,
    double *accSumSq,
    uint32_t *accCnt,
    PhenoGenoStats *out,
    size_t outStride
) {
    for (uint32_t g = 0; g < nGroups; ++g) {
        accGSum[g] = 0.0;
        accSumSq[g] = 0.0;
        accCnt[g * 4 + 0] = 0;
        accCnt[g * 4 + 1] = 0;
        accCnt[g * 4 + 2] = 0;
        accCnt[g * 4 + 3] = 0;
    }

    for (uint32_t i = 0; i < nUnion; ++i) {
        const double v = unionG[i];
        const double vv = v * v;
        // Classify once per subject (mask-independent) for the HWE counts.
        // origMissBits records subjects that were NaN before the caller's mean-
        // imputation: these are genuine-missing (cls 3 — excluded from hom/het,
        // counted for missingRate).  Non-missing subjects use the plink2 hard-
        // call threshold: within thr of an integer → that hard-call (cls 0/1/2),
        // otherwise HWE-uncertain (cls 4 — kept in gSum/sumSq for AF but excluded
        // from hom/het and from missingRate; not stored in accCnt).  For pure
        // hard-calls with no missing, cls == v exactly, reproducing the previous
        // 0/1/2 counting bit-for-bit.  cls indexes accCnt: 0 hom-ref, 1 het,
        // 2 hom-alt, 3 genuine-missing.
        uint32_t cls;
        const bool origMiss = origMissBits &&
            ((origMissBits[i >> 6] >> (i & 63)) & 1ULL);
        if (origMiss) {
            cls = 3;
        } else {
            const int hc = dosageHardcall(v, thr);
            cls = (hc < 0) ? 4u : static_cast<uint32_t>(hc);
        }
        const uint32_t w = i >> 6;
        const uint32_t b = i & 63;
        for (uint32_t g = 0; g < nGroups; ++g) {
            const uint64_t m = (groupMaskBits[g][w] >> b) & 1ULL;
            const double md = static_cast<double>(m);
            accGSum[g] += v * md;       // += v when masked-in, += +0.0 otherwise
            accSumSq[g] += vv * md;
            if (cls < 4) accCnt[g * 4 + cls] += static_cast<uint32_t>(m);
        }
    }

    for (uint32_t g = 0; g < nGroups; ++g) {
        const uint32_t nHomRef = accCnt[g * 4 + 0];
        const uint32_t nHet = accCnt[g * 4 + 1];
        const uint32_t nHomAlt = accCnt[g * 4 + 2];
        const uint32_t nMissing = accCnt[g * 4 + 3];  // genuine-missing only
        PhenoGenoStats &s = out[g * outStride];
        s.sumSq = accSumSq[g];
        s.gSum = accGSum[g];
        const uint32_t nonMissing = nHomRef + nHet + nHomAlt;
        if (nonMissing == 0) {
            // No hard-called subjects (all missing/uncertain): keep the prior
            // all-missing sentinel; the caller recomputes AF from gSum (= 0).
            s.missingRate = 0.0;
            s.log10pHwe = 0.0;   // -log10(1): the prior all-missing sentinel
            s.altFreq = 0.0;
            s.mac = 0.0;
            s.fromDosage = true;
            continue;
        }
        // HWE from the thresholded hard-call counts (plink2 --hardy); this is
        // computed for dosage markers too.  missingRate counts genuine-missing
        // only — HWE-uncertain dosages are excluded from it (they remain in the
        // dosage sum used for AF).
        s.log10pHwe = hweNegLog10P(nHet, nHomAlt, nHomRef);
        s.missingRate =
            static_cast<double>(nMissing) / static_cast<double>(groupNUsed[g]);
        if (isDosage) {
            // AF/MAC recomputed from the GEMM gSum by the caller (the un-binned
            // dosage frequency, not the hard-call count ratio).
            s.altFreq = 0.0;
            s.mac = 0.0;
            s.fromDosage = true;
        } else {
            const uint32_t altCounts = 2 * nHomAlt + nHet;
            const double n2 = 2.0 * nonMissing;
            s.altFreq = altCounts / n2;
            const double maf = std::min(s.altFreq, 1.0 - s.altFreq);
            s.mac = maf * n2;
            s.fromDosage = false;
        }
    }
}

inline PhenoGenoStats statsFromGVec(
    const double *g,
    uint32_t n,
    double thr,
    std::vector<uint32_t> &indexForMissing
) {
    indexForMissing.clear();
    uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0, nNonMissing = 0;
    double dosageSum = 0.0;
    bool isDosage = false;
    for (uint32_t i = 0; i < n; ++i) {
        const double v = g[i];
        if (std::isnan(v) || v < 0.0) {
            indexForMissing.push_back(i);
            continue;
        }
        dosageSum += v;
        ++nNonMissing;
        // AF track: any fractional value marks a dosage marker (AF then comes
        // from the dosage sum, not the hard-call counts).  Unchanged detection.
        if (v != 0.0 && v != 1.0 && v != 2.0)
            isDosage = true;
        // HWE track: classify via the plink2 hard-call threshold.  Exact 0/1/2
        // always pass (|v-round| == 0); dosages farther than thr from an integer
        // are HWE-uncertain and excluded from hom/het (still in the dosage sum).
        const int hc = dosageHardcall(v, thr);
        if (hc == 0) ++nHomRef;
        else if (hc == 1) ++nHet;
        else if (hc == 2) ++nHomAlt;
    }
    PhenoGenoStats s;
    s.sumSq = 0.0;  // unused on the non-fused path
    s.gSum = 0.0;   // unused on the non-fused path (kept initialized)
    s.fromDosage = isDosage;
    if (nNonMissing == 0) {
        s.altFreq = NAN;
        s.mac = NAN;
        s.missingRate = 1.0;
        s.log10pHwe = NAN;
        return s;
    }
    double n2 = 2.0 * nNonMissing;
    // HWE from the thresholded hard-call counts (plink2 --hardy); genuine-missing
    // (NaN) is already excluded, and HWE-uncertain dosages were excluded above.
    // Computed for dosage markers too (previously skipped as NaN).
    s.log10pHwe = hweNegLog10P(nHet, nHomAlt, nHomRef);
    s.missingRate = static_cast<double>(indexForMissing.size()) / n;
    if (isDosage) {
        // Dosage marker: AF/MAC from the un-binned dosage sum.
        s.altFreq = dosageSum / n2;
        double maf = std::min(s.altFreq, 1.0 - s.altFreq);
        s.mac = maf * n2;
    } else {
        uint32_t altCounts = 2 * nHomAlt + nHet;
        s.altFreq = altCounts / n2;
        double maf = std::min(s.altFreq, 1.0 - s.altFreq);
        s.mac = maf * n2;
    }
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

#if defined(__x86_64__) || defined(_M_X64)
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
#endif

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
#if defined(__x86_64__) || defined(_M_X64)
    if (simdLevel() >= SimdLevel::AVX512) {
        extractCompress_avx512(unionG, presentMask, nUnion, phenoG);
    } else {
        extractGather(unionG, presentIndices, nUsed, phenoG);
    }
#else
    (void)presentMask; (void)nUnion;
    extractGather(unionG, presentIndices, nUsed, phenoG);
#endif
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
            s.log10pHwe = NAN;
        } else {
            const uint32_t nP = nPhenoArr[p];
            const uint32_t altCounts = 2 * acc[p].nHomAlt + acc[p].nHet;
            const double n2 = 2.0 * nonMissing;
            s.altFreq = altCounts / n2;
            const double maf = std::min(s.altFreq, 1.0 - s.altFreq);
            s.mac = maf * n2;
            s.missingRate = static_cast<double>(outMissings[p].size()) / nP;
            s.log10pHwe = hweNegLog10P(acc[p].nHet, acc[p].nHomAlt, acc[p].nHomRef);
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
        s.log10pHwe = NAN;
        return s;
    }
    const uint32_t altCounts = 2 * counts[2] + counts[1];
    const double n2 = 2.0 * nonMissing;
    s.altFreq = altCounts / n2;
    const double maf = std::min(s.altFreq, 1.0 - s.altFreq);
    s.mac = maf * n2;
    s.missingRate = static_cast<double>(counts[3]) / nPheno;
    s.log10pHwe = hweNegLog10P(counts[1], counts[2], counts[0]);
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
