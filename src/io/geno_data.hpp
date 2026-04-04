// geno_data.hpp — Abstract interfaces for genotype data backends
//
// GenoMeta:   shared metadata + cursor factory (constructed once, main thread)
// GenoCursor: per-thread lightweight genotype decoder
//
// Concrete backends (PlinkData, PgenData, VcfData, BgenData)
// derive from these interfaces.  The marker engine and helper functions
// operate through the abstract interface so they work with any backend.
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>
#include <Eigen/Dense>

class GenoCursor;


// ======================================================================
// GenoSpec — lightweight descriptor for the genotype data source
// ======================================================================

enum class GenoFormat { Plink, Pgen, Vcf, Bgen };

struct GenoSpec {
    GenoFormat  format;
    std::string path;  // bfilePrefix for Plink/Pgen, full file for Vcf/Bgen
    std::string extractFile;  // --extract: SNP include list
    std::string excludeFile;  // --exclude: SNP exclude list
};

// Extract all sample IIDs from the genotype file metadata
// (plink .fam, pgen .psam, VCF header, BGEN sample block).
std::vector<std::string> parseGenoIIDs(const GenoSpec& spec);

// Construct the appropriate GenoMeta backend from a spec + subject mask.
std::unique_ptr<class GenoMeta> makeGenoData(
    const GenoSpec& spec,
    const std::vector<uint64_t>& usedMask,
    uint32_t nSamplesInFile,
    uint32_t nUsed,
    int nMarkersEachChunk = 1024);


// ======================================================================
// GenoMeta — shared, read-only genotype metadata + cursor factory
// ======================================================================

class GenoMeta {
public:
  struct MarkerInfo {
    std::string chrom;
    uint32_t    pos;
    std::string id;
    std::string ref;
    std::string alt;
    uint64_t    genoIndex;  // backend-specific row index
  };

  virtual ~GenoMeta() = default;

  virtual uint32_t nMarkers()    const = 0;
  virtual uint32_t nSubjUsed()   const = 0;
  virtual uint32_t nSubjInFile() const = 0;

  virtual std::string_view chr(uint64_t i)      const = 0;
  virtual std::string_view markerId(uint64_t i) const = 0;
  virtual uint32_t         pos(uint64_t i)      const = 0;
  virtual std::string_view ref(uint64_t i)      const = 0;
  virtual std::string_view alt(uint64_t i)      const = 0;

  virtual const std::vector<MarkerInfo>&            markerInfo()   const = 0;
  virtual const std::vector<std::vector<uint64_t>>& chunkIndices() const = 0;

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
    double& altFreq,
    double& altCounts,
    double& missingRate,
    double& hweP,
    double& maf,
    double& mac,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe = true) = 0;

  // Lightweight variant: genotype vector only, missing → NaN, no QC stats.
  virtual void getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out) = 0;
};
