// plink.hpp — PLINK binary genotype reader (pure C++17 / Eigen)
//
// plink2-style bitmask design: subjects are never reordered.  A bitmask
// marks which .fam subjects are "used"; decoding writes only those
// subjects into a dense output vector in .fam order.
#pragma once

#include "geno_factory/geno_data.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

// ======== PlinkData: shared, constructed once on main thread ========
// Parses .bim, counts .fam lines, validates .bed, applies marker
// filters, builds chunks.  Stores a bitmask of used subjects.
// All accessors are const — safe to share across worker threads.

class PlinkData : public GenoMeta {
  public:
    // `usedMask` has ceil(nFam/64) words; bit i set ↔ .fam subject i is used.
    // `nFam` = total .fam lines.  `nUsed` = popcount(usedMask).
    PlinkData(
        std::string bedFile,
        std::string bimFile,
        std::string famFile,
        const std::vector<uint64_t> &usedMask,
        uint32_t nFam,
        uint32_t nUsed,
        std::string IDsToIncludeFile = {},
        std::string RangesToIncludeFile = {},
        std::string IDsToExcludeFile = {},
        std::string RangesToExcludeFile = {},
        std::unordered_set<std::string> chrFilter = {},
        int nMarkersEachChunk = 1024
    );

    // MarkerInfo is inherited from GenoMeta.

    // ---- GenoMeta overrides ----
    uint32_t nMarkers() const override {
        return m_nMarkers;
    }

    uint32_t nSubjUsed() const override {
        return m_nSubjUsed;
    }

    uint32_t nSubjInFile() const override {
        return m_nSubjInFile;
    }

    std::string_view chr(uint64_t i) const override {
        return m_chr[i];
    }

    std::string_view markerId(uint64_t i) const override {
        return m_markerId[i];
    }

    uint32_t pos(uint64_t i) const override {
        return m_pos[i];
    }

    std::string_view ref(uint64_t i) const override {
        return m_ref[i];
    }

    std::string_view alt(uint64_t i) const override {
        return m_alt[i];
    }

    const std::vector<MarkerInfo> &markerInfo() const override {
        return m_markerInfo;
    }

    const std::vector<std::vector<uint64_t> > &chunkIndices() const override {
        return m_chunkIndices;
    }

    std::unique_ptr<GenoCursor> makeCursor() const override;

    // ---- PLINK-specific accessors ----
    uint64_t bytesPerMarker() const {
        return m_bytesPerMarker;
    }

    const std::string &bedFile() const {
        return m_bedFile;
    }

    bool allUsed() const {
        return m_allUsed;
    }

    const std::vector<uint64_t> &usedMask() const {
        return m_usedMask;
    }

    uint32_t nMaskWords() const {
        return static_cast<uint32_t>(m_usedMask.size());
    }

  private:
    static std::vector<MarkerInfo> getFilteredMarkers(
        const std::vector<std::string> &chr,
        const std::vector<uint32_t> &pos,
        const std::vector<std::string> &markerId,
        const std::vector<std::string> &ref,
        const std::vector<std::string> &alt,
        const std::string &IDsToIncludeFile,
        const std::string &RangesToIncludeFile,
        const std::string &IDsToExcludeFile,
        const std::string &RangesToExcludeFile,
        const std::unordered_set<std::string> &chrFilter
    );

    static std::vector<std::vector<uint64_t> > buildChunks(
        const std::vector<MarkerInfo> &markers,
        int chunkSize
    );

    std::string m_bedFile;
    bool m_allUsed;
    uint32_t m_nSubjInFile;
    uint32_t m_nSubjUsed;
    uint32_t m_nMarkers;
    uint64_t m_bytesPerMarker;
    std::vector<uint64_t> m_usedMask;
    std::vector<std::string> m_chr;
    std::vector<std::string> m_markerId;
    std::vector<uint32_t> m_pos;
    std::vector<std::string> m_ref;
    std::vector<std::string> m_alt;
    std::vector<MarkerInfo> m_markerInfo;
    std::vector<std::vector<uint64_t> > m_chunkIndices;
};

// ======== PlinkCursor: per-thread, lightweight ========
// Each worker thread creates its own PlinkCursor with an independent
// file handle and scratch buffers.  Uses bitmask to decode only used
// subjects into a dense output vector.

class PlinkCursor : public GenoCursor {
  public:
    PlinkCursor(
        const std::string &bedFile,
        uint32_t nBimLines,
        uint32_t nFamLines,
        const std::vector<uint64_t> &usedMask,
        uint32_t nUsed,
        bool allUsed
    );

    PlinkCursor(const PlinkCursor &other);

    // Prepare sequential reading starting from a given marker index.
    void beginSequentialBlock(uint64_t firstMarker) override;

    // Decode genotype for marker gIndex into caller-owned Eigen vector.
    // Returns QC statistics through the output parameters.
    void getGenotypes(
        uint64_t gIndex,
        Eigen::Ref<Eigen::VectorXd> out,
        double &altFreq,
        double &altCounts,
        double &missingRate,
        double &hweP,
        double &maf,
        double &mac,
        std::vector<uint32_t> &indexForMissing
    ) override;

    // Lightweight variant: genotype vector only, missing → NaN, no QC stats.
    void getGenotypesSimple(
        uint64_t gIndex,
        Eigen::Ref<Eigen::VectorXd> out
    ) override;

  private:
    void loadBlock(uint64_t startMarker);

    const uint8_t *readMarkerPtr(uint64_t gIndex);

    // ---- All-used fast path (identity decode, no scatter) ----
    void getGenotypesAllUsed(
        const uint8_t *raw,
        uint32_t n,
        Eigen::Ref<Eigen::VectorXd> out,
        double &altFreq,
        double &altCounts,
        double &missingRate,
        double &hweP,
        double &maf,
        double &mac,
        std::vector<uint32_t> &indexForMissing
    );

    // ---- Bitmask path (decode only set-bit subjects) ----
    void getGenotypesMasked(
        const uint8_t *raw,
        Eigen::Ref<Eigen::VectorXd> out,
        double &altFreq,
        double &altCounts,
        double &missingRate,
        double &hweP,
        double &maf,
        double &mac,
        std::vector<uint32_t> &indexForMissing
    );

    std::string m_bedFile;
    std::ifstream m_bedStream;
    uint32_t m_nSubjInFile;
    uint32_t m_bytesPerMarker;
    uint32_t m_nMarkers;
    uint32_t m_nUsed;
    std::vector<uint64_t> m_usedMask;
    bool m_allUsed;

    // Pre-allocated scratch buffer (one per cursor)
    std::vector<uint8_t> m_rawBytes;

    // Word-aligned scratch for plink2 SIMD primitives
    std::vector<uintptr_t> m_alignedBed;

    // Block buffer for sequential reads
    std::vector<uint8_t> m_blockBytes;
    uint64_t m_blockStart = 0;
    uint64_t m_blockEnd = 0;
    bool m_hasBlock = false;
    bool m_hasSeqCursor = false;
    uint64_t m_nextMarker = 0;

    static constexpr uint64_t BLOCK_CAPACITY = 512;
};
