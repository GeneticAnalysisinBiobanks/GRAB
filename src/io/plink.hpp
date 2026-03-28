// plink.hpp — PLINK binary genotype reader (pure C++17 / Eigen)
#pragma once

#include <cstdint>
#include <cmath>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <Eigen/Dense>


// ======== PlinkData: shared, constructed once on main thread ========
// Parses .bim/.fam, validates .bed, applies marker filters, builds chunks.
// All accessors are const — safe to share across worker threads without locks.

class PlinkData {
public:

  PlinkData(
    std::string bedFile,
    std::string bimFile,
    std::string famFile,
    std::vector<std::string> subjData,
    std::string AlleleOrder,
    std::string IDsToIncludeFile   = {},
    std::string RangesToIncludeFile = {},
    std::string IDsToExcludeFile   = {},
    std::string RangesToExcludeFile = {},
    int nMarkersEachChunk = 1024
  );

  struct MarkerInfo {
    std::string chrom;
    uint32_t    pos;
    std::string id;
    std::string ref;
    std::string alt;
    uint64_t    genoIndex;  // 0-based .bim line = row in .bed
  };

  // ---- Accessors ----
  uint32_t nMarkers()       const { return m_nMarkers; }
  uint32_t nSubjUsed()      const { return m_nSubjUsed; }
  uint32_t nSubjInFile()    const { return m_nSubjInFile; }
  uint64_t bytesPerMarker() const { return m_bytesPerMarker; }
  const std::string& bedFile()  const { return m_bedFile; }
  bool isAltFirst()             const { return m_altFirst; }

  // True when posMap is an identity permutation (all subjects used, in order).
  // Enables a SIMD fast-path that skips the scatter/gather step entirely.
  bool isIdentityMap()          const { return m_identityMap; }

  const std::vector<uint32_t>& samplePosMap() const { return m_samplePosMap; }

  // Per-marker metadata by genoIndex (0-based .bim line number)
  std::string_view chr(uint64_t i)      const { return m_chr[i]; }
  std::string_view markerId(uint64_t i) const { return m_markerId[i]; }
  uint32_t         pos(uint64_t i)      const { return m_pos[i]; }
  std::string_view ref(uint64_t i)      const { return m_ref[i]; }
  std::string_view alt(uint64_t i)      const { return m_alt[i]; }

  const std::vector<MarkerInfo>&            markerInfo()   const { return m_markerInfo; }
  const std::vector<std::vector<uint64_t>>& chunkIndices() const { return m_chunkIndices; }

private:

  static std::vector<MarkerInfo> getFilteredMarkers(
    const std::vector<std::string>& chr,
    const std::vector<uint32_t>& pos,
    const std::vector<std::string>& markerId,
    const std::vector<std::string>& ref,
    const std::vector<std::string>& alt,
    const std::string& IDsToIncludeFile,
    const std::string& RangesToIncludeFile,
    const std::string& IDsToExcludeFile,
    const std::string& RangesToExcludeFile
  );

  static std::vector<std::vector<uint64_t>> buildChunks(
    const std::vector<MarkerInfo>& markers, int chunkSize
  );

  std::string              m_bedFile;
  bool                     m_altFirst;
  bool                     m_identityMap;
  uint32_t                 m_nSubjInFile;
  uint32_t                 m_nSubjUsed;
  uint32_t                 m_nMarkers;
  uint64_t                 m_bytesPerMarker;
  std::vector<uint32_t>    m_samplePosMap;
  std::vector<std::string> m_chr;
  std::vector<std::string> m_markerId;
  std::vector<uint32_t>    m_pos;
  std::vector<std::string> m_ref;
  std::vector<std::string> m_alt;
  std::vector<MarkerInfo>  m_markerInfo;
  std::vector<std::vector<uint64_t>> m_chunkIndices;
};


// ======== PlinkCursor: per-thread, lightweight ========
// Each worker thread creates its own PlinkCursor with an independent
// file handle and scratch buffers.  The copy constructor re-opens the
// .bed file so copies are fully independent.

class PlinkCursor {
public:
  PlinkCursor(const std::string& bedFile,
              uint32_t nBimLines,
              uint32_t nFamLines,
              const std::vector<uint32_t>& posMap,
              bool altFirst,
              bool identityMap);

  PlinkCursor(const PlinkCursor& other);

  // Prepare sequential reading starting from a given marker index.
  void beginSequentialBlock(uint64_t firstMarker);

  // Decode genotype for marker gIndex into caller-owned Eigen vector.
  // Returns QC statistics through the output parameters.
  void getGenotypes(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    double& hweP,
    double& maf,
    double& mac,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe = true
  );

  // Lightweight variant: genotype vector only, missing → NaN, no QC stats.
  // Use this when you only need the genotype values (e.g., AF scans).
  void getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out
  );

private:
  void loadBlock(uint64_t startMarker);
  const uint8_t* readMarkerPtr(uint64_t gIndex);

  // ---- Identity-map fast path (no repack, no scatter) ----
  void getGenotypesIdentity(
    const uint8_t* raw,
    uint32_t n,
    Eigen::Ref<Eigen::VectorXd> out,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    double& hweP,
    double& maf,
    double& mac,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe
  );

  // ---- Subset path (repack + AVX2 stats + scalar decode) ----
  void getGenotypesSubset(
    const uint8_t* raw,
    uint32_t n,
    Eigen::Ref<Eigen::VectorXd> out,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    double& hweP,
    double& maf,
    double& mac,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe
  );

  std::string   m_bedFile;
  std::ifstream m_bedStream;
  uint32_t      m_nSubjInFile;
  uint32_t      m_bytesPerMarker;
  uint32_t      m_nMarkers;
  std::vector<uint32_t> m_posMap;
  bool          m_altFirst;
  bool          m_identityMap;

  // Pre-allocated scratch buffers (one per cursor, not per call)
  std::vector<uint8_t> m_rawBytes;
  std::vector<uint8_t> m_compactGeno;

  // Block buffer for sequential reads
  std::vector<uint8_t> m_blockBytes;
  uint64_t m_blockStart = 0;
  uint64_t m_blockEnd   = 0;
  bool     m_hasBlock   = false;
  bool     m_hasSeqCursor = false;
  uint64_t m_nextMarker   = 0;

  static constexpr uint64_t BLOCK_CAPACITY = 512;
};
