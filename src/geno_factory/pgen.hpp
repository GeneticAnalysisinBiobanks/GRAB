// pgen.hpp — PLINK 2 pgen/pvar genotype reader
//
// PgenData:   GenoMeta implementation for .pgen + .pvar [+ .psam]
// PgenCursor: per-thread GenoCursor backed by a pgenlib PgenReader
//
// Uses pgenlib from third_party/plink2-a.6.33/include.
// Reads variant metadata from .pvar (tab/space-delimited, with header).
// Sample count comes from the .pgen header; the caller's usedMask
// (built by SubjectData) defines which samples to include.
#pragma once

#include "geno_factory/geno_data.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

// Forward-declare pgenlib types to avoid leaking the heavy header.
// The .cpp includes the actual pgenlib headers.
namespace plink2 {
struct PgenFileInfoStruct;
}
using PgenFileInfo = plink2::PgenFileInfoStruct;

class PgenData : public GenoMeta {
  public:
    // pgenFile: path to .pgen
    // pvarFile: path to .pvar (tab-delimited with #CHROM POS ID REF ALT header)
    // usedMask, nSamplesInFile, nUsed: from SubjectData (same semantics as PlinkData)
    PgenData(
        std::string pgenFile,
        std::string pvarFile,
        const std::vector<uint64_t> &usedMask,
        uint32_t nSamplesInFile,
        uint32_t nUsed,
        std::unordered_set<std::string> chrFilter = {},
        int nMarkersEachChunk = 1024
    );

    ~PgenData() override;

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

    // ---- Pgen-specific accessors ----
    const std::string &pgenFile() const {
        return m_pgenFile;
    }

    bool allUsed() const {
        return m_allUsed;
    }

    const std::vector<uint64_t> &usedMask() const {
        return m_usedMask;
    }

    // Shared pgfi for all cursors (read-only after init).
    PgenFileInfo *pgfi() const {
        return m_pgfi.get();
    }

    uint32_t maxVrecWidth() const {
        return m_maxVrecWidth;
    }

    std::size_t pgrAllocCachelineCt() const {
        return m_pgrAllocCachelineCt;
    }

  private:
    void parsePvar(const std::string &pvarFile);

    static std::vector<std::vector<uint64_t> > buildChunks(
        const std::vector<MarkerInfo> &markers,
        int chunkSize
    );

    std::string m_pgenFile;
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

    // pgenlib state
    struct PgfiDeleter {
        void operator()(PgenFileInfo *p) const;

    };

    std::unique_ptr<PgenFileInfo, PgfiDeleter> m_pgfi;
    std::unique_ptr<unsigned char, void (*)(void *)> m_pgfiAlloc;
    uint32_t m_maxVrecWidth = 0;
    std::size_t m_pgrAllocCachelineCt = 0;
};

class PgenCursor : public GenoCursor {
  public:
    PgenCursor(const PgenData &parent);
    ~PgenCursor() override;

    PgenCursor(const PgenCursor &) = delete;
    PgenCursor &operator=(const PgenCursor &) = delete;

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

    uint32_t getGenotypesMaybeSparse(
        uint64_t gIndex,
        Eigen::Ref<Eigen::VectorXd> out,
        uint32_t maxLen,
        uint32_t *diffSampleIds,
        uint8_t *diffGenoCodes,
        uint32_t &diffLen
    ) override;

  private:
    const PgenData &m_parent;

    // Per-thread pgenlib state (opaque, allocated in .cpp)
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
