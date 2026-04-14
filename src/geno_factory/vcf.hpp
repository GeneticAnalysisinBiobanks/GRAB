// vcf.hpp — VCF/BCF genotype reader via htslib
//
// VcfData:   GenoMeta implementation for .vcf/.vcf.gz/.bcf
// VcfCursor: per-thread GenoCursor backed by independent htsFile handles
//
// Reads all biallelic variants into memory at construction time (variant
// metadata only; genotypes are re-read per cursor).  This two-pass design
// keeps the metadata read fast while allowing each worker thread to seek
// independently.
#pragma once

#include "geno_factory/geno_data.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

class VcfData : public GenoMeta {
  public:
    // vcfFile: path to .vcf / .vcf.gz / .bcf
    // usedMask, nSamplesInFile, nUsed: from SubjectData
    VcfData(
        std::string vcfFile,
        const std::vector<uint64_t> &usedMask,
        uint32_t nSamplesInFile,
        uint32_t nUsed,
        std::unordered_set<std::string> chrFilter = {},
        int nMarkersEachChunk = 1024
    );

    ~VcfData() override;

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

    // ---- VCF-specific accessors ----
    const std::string &vcfFile() const {
        return m_vcfFile;
    }

    bool allUsed() const {
        return m_allUsed;
    }

    const std::vector<uint64_t> &usedMask() const {
        return m_usedMask;
    }

  private:
    static std::vector<std::vector<uint64_t> > buildChunks(
        const std::vector<MarkerInfo> &markers,
        int chunkSize
    );

    std::string m_vcfFile;
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
};

class VcfCursor : public GenoCursor {
  public:
    VcfCursor(const VcfData &parent);
    ~VcfCursor() override;

    VcfCursor(const VcfCursor &) = delete;
    VcfCursor &operator=(const VcfCursor &) = delete;

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
    const VcfData &m_parent;

    struct Impl;
    std::unique_ptr<Impl> m_impl;
};
