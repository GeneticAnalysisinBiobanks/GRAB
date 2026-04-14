// marker_impl.hpp — Shared helpers for marker engines (formatting, stats)
//
// Used by marker.cpp (markerEngine, multiPhenoEngine) and loco.cpp (locoEngine).
// Not part of the public API — include only from engine translation units.
#pragma once

#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "geno_factory/hwe.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>
#include <string_view>
#include <vector>

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

// Format double: "NA" | "Inf" | "-Inf" | "%.6g".  Returns char count.
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
    return std::snprintf(buf, 32, "%.6g", x);
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
