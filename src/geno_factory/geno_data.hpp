// geno_data.hpp — Abstract interfaces for genotype data backends
//
// GenoMeta:   shared metadata + cursor factory (constructed once, main thread)
// GenoCursor: per-thread lightweight genotype decoder
//
// Concrete backends (PlinkData, PgenData, VcfData, BgenData)
// derive from these interfaces.  The marker engine and helper functions
// operate through the abstract interface so they work with any backend.
#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

class GenoCursor;

// ======================================================================
// GenoSpec — lightweight descriptor for the genotype data source
// ======================================================================

enum class GenoFormat { Plink, Pgen, Vcf, Bgen };

struct GenoSpec {
    GenoFormat format;
    std::string path;        // bfilePrefix for Plink/Pgen, full file for Vcf/Bgen
    std::string extractFile; // --extract: SNP include list
    std::string excludeFile; // --exclude: SNP exclude list
    std::unordered_set<std::string> chrFilter; // --chr: chromosome include set

    // Return a descriptive flag label, e.g. "--bfile d_bed".
    std::string flagLabel() const {
        const char *flag = nullptr;
        switch (format) {
        case GenoFormat::Plink: flag = "--bfile"; break;
        case GenoFormat::Pgen:  flag = "--pfile"; break;
        case GenoFormat::Vcf:   flag = "--vcf";   break;
        case GenoFormat::Bgen:  flag = "--bgen";  break;
        }
        return std::string(flag) + " " + path;
    }

};

// Return a descriptive label for the GRM source, e.g. "--sp-grm-grab e_grm.grab".
inline std::string grmFlagLabel(
    const std::string &grabFile,
    const std::string &gctaFile
) {
    if (!grabFile.empty()) return "--sp-grm-grab " + grabFile;
    if (!gctaFile.empty()) return "--sp-grm-plink2 " + gctaFile;
    return "";
}

// Extract all sample IIDs from the genotype file metadata
// (plink .fam, pgen .psam, VCF header, BGEN sample block).
std::vector<std::string> parseGenoIIDs(const GenoSpec &spec);

// Construct the appropriate GenoMeta backend from a spec + subject mask.
std::unique_ptr<class GenoMeta> makeGenoData(
    const GenoSpec &spec,
    const std::vector<uint64_t> &usedMask,
    uint32_t nSamplesInFile,
    uint32_t nUsed,
    int nMarkersEachChunk = 1024
);

// ======================================================================
// GenoMeta — shared, read-only genotype metadata + cursor factory
// ======================================================================

class GenoMeta {
  public:
    struct MarkerInfo {
        std::string chrom;
        uint32_t pos;
        std::string id;
        std::string ref;
        std::string alt;
        uint64_t genoIndex; // backend-specific row index
    };

    virtual ~GenoMeta() = default;

    virtual uint32_t nMarkers() const = 0;

    virtual uint32_t nSubjUsed() const = 0;

    virtual uint32_t nSubjInFile() const = 0;

    virtual std::string_view chr(uint64_t i) const = 0;

    virtual std::string_view markerId(uint64_t i) const = 0;

    virtual uint32_t pos(uint64_t i) const = 0;

    virtual std::string_view ref(uint64_t i) const = 0;

    virtual std::string_view alt(uint64_t i) const = 0;

    virtual const std::vector<MarkerInfo> &markerInfo() const = 0;

    virtual const std::vector<std::vector<uint64_t> > &chunkIndices() const = 0;

    // Factory: create a per-thread cursor for genotype decoding.
    virtual std::unique_ptr<GenoCursor> makeCursor() const = 0;

};

// ======================================================================
// GenoCursor — per-thread, lightweight genotype decoder
// ======================================================================

class GenoCursor {
  public:
    virtual ~GenoCursor() = default;

    // Prepare sequential reading starting from a given marker index.
    virtual void beginSequentialBlock(uint64_t firstMarker) = 0;

    // Decode genotype for marker gIndex into caller-owned Eigen vector.
    // Returns QC statistics through the output parameters.
    virtual void getGenotypes(
        uint64_t gIndex,
        Eigen::Ref<Eigen::VectorXd> out,
        double &altFreq,
        double &altCounts,
        double &missingRate,
        double &hweP,
        double &maf,
        double &mac,
        std::vector<uint32_t> &indexForMissing
    ) = 0;

    // Lightweight variant: genotype vector only, missing → NaN, no QC stats.
    virtual void getGenotypesSimple(
        uint64_t gIndex,
        Eigen::Ref<Eigen::VectorXd> out
    ) = 0;

    // Optional sparse genotype read via pgenlib difflist.
    //
    // If the variant qualifies as sparse (difflist_len <= maxLen):
    //   Returns commonGeno (0, 1, or 2).  diffSampleIds[0..diffLen-1]
    //   contains subsetted sample indices; diffGenoCodes[0..diffLen-1]
    //   contains 2-bit genotype codes (0/1/2/3).  `out` is NOT populated.
    //
    // If dense (commonGeno == UINT32_MAX):
    //   `out` is populated as in getGenotypesSimple.  Diff arrays ignored.
    //
    // Buffers must be sized for at least maxLen entries.
    // Default: always returns dense.
    virtual uint32_t getGenotypesMaybeSparse(
        uint64_t gIndex,
        Eigen::Ref<Eigen::VectorXd> out,
        uint32_t maxLen,
        uint32_t *diffSampleIds,
        uint8_t *diffGenoCodes,
        uint32_t &diffLen
    ) {
        getGenotypesSimple(gIndex, out);
        return UINT32_MAX;
    }

};
