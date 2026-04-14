// bgen.hpp — BGEN genotype reader
//
// BgenData:   GenoMeta implementation for .bgen files
// BgenCursor: per-thread GenoCursor using the low-level bgen API
//
// Uses the bgen reference library from third_party/bgen-1.2.0.
// Supports Layout 1 and Layout 2 (with zlib or zstd compression).
// Only biallelic diploid variants are processed; multiallelic are skipped.
#pragma once

#include "geno_factory/geno_data.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

class BgenData : public GenoMeta {
  public:
    // bgenFile: path to .bgen
    // sampleFile: optional path to .sample file (for sample IDs); empty = use embedded IDs
    // usedMask, nSamplesInFile, nUsed: from SubjectData
    BgenData(
        std::string bgenFile,
        const std::vector<uint64_t> &usedMask,
        uint32_t nSamplesInFile,
        uint32_t nUsed,
        std::unordered_set<std::string> chrFilter = {},
        int nMarkersEachChunk = 1024
    );

    ~BgenData() override;

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

    // ---- BGEN-specific accessors ----
    const std::string &bgenFile() const {
        return m_bgenFile;
    }

    bool allUsed() const {
        return m_allUsed;
    }

    const std::vector<uint64_t> &usedMask() const {
        return m_usedMask;
    }

    // File byte offset of the first variant data (after header + sample block).
    std::streampos dataOffset() const {
        return m_dataOffset;
    }

    uint32_t bgenFlags() const {
        return m_bgenFlags;
    }

  private:
    static std::vector<std::vector<uint64_t> > buildChunks(
        const std::vector<MarkerInfo> &markers,
        int chunkSize
    );

    std::string m_bgenFile;
    bool m_allUsed;
    uint32_t m_nSubjInFile;
    uint32_t m_nSubjUsed;
    uint32_t m_nMarkers;

    std::vector<uint64_t> m_usedMask;

    std::vector<std::string> m_chr;
    std::vector<std::string> m_markerId;
    std::vector<uint32_t> m_pos;
    std::vector<std::string> m_ref;
    std::vector<std::string> m_alt;
    std::vector<MarkerInfo> m_markerInfo;
    std::vector<std::vector<uint64_t> > m_chunkIndices;

    std::streampos m_dataOffset = 0;
    uint32_t m_bgenFlags = 0;
};

class BgenCursor : public GenoCursor {
  public:
    BgenCursor(const BgenData &parent);
    ~BgenCursor() override;

    BgenCursor(const BgenCursor &) = delete;
    BgenCursor &operator=(const BgenCursor &) = delete;

    void beginSequentialBlock(uint64_t firstMarker) override;

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

    void getGenotypesSimple(
        uint64_t gIndex,
        Eigen::Ref<Eigen::VectorXd> out
    ) override;

  private:
    const BgenData &m_parent;

    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
