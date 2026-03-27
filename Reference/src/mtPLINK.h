#ifndef PLINK_H
#define PLINK_H

// mtPLINK.h -- PLINK binary genotype reader

#include <cstdint>
#include <cmath>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include <memory>
#include <RcppArmadillo.h>


// ======== PlinkData: shared, constructed once on main thread ========

class PlinkData {

public:
  // Constructor — parses .bim/.fam, validates .bed, applies filters, builds chunks.
  // Throws if no markers remain after filtering.
  PlinkData(
    std::string bedFile,
    std::string bimFile,
    std::string famFile,
    std::vector<std::string> subjData,
    std::string AlleleOrder,
    std::string IDsToIncludeFile,
    std::string RangesToIncludeFile,
    std::string IDsToExcludeFile,
    std::string RangesToExcludeFile,
    int nMarkersEachChunk = 1024
  );

  // Marker metadata from a .bim file (alleles normalized by AlleleOrder).
  struct MarkerInfo {
    std::string chrom;
    uint32_t    pos;
    std::string id;
    std::string ref;
    std::string alt;
    uint64_t    genoIndex;   // 0-based line number in .bim = row in .bed
  };

  // ---- Accessors ----
  uint32_t nMarkers()       const { return m_nMarkers; }
  uint32_t nSubjUsed()      const { return m_nSubjUsed; }
  uint32_t nSubjInFile()    const { return m_nSubjInFile; }
  uint64_t bytesPerMarker() const { return m_bytesPerMarker; }
  const std::string& bedFile()  const { return m_bedFile; }
  bool isAltFirst()             const { return m_altFirst; }
  const std::vector<uint32_t>& samplePosMap() const { return m_samplePosMap; }

  // Per-marker metadata by genoIndex
  std::string_view chr(uint64_t i)      const { return m_chr[i]; }
  std::string_view markerId(uint64_t i) const { return m_markerId[i]; }
  uint32_t         pos(uint64_t i)      const { return m_pos[i]; }
  std::string_view ref(uint64_t i)      const { return m_ref[i]; }
  std::string_view alt(uint64_t i)      const { return m_alt[i]; }

  // Filtered marker list and chunk indices built during construction
  const std::vector<MarkerInfo>&              markerInfo()   const { return m_markerInfo; }
  const std::vector<std::vector<uint64_t>>&   chunkIndices() const { return m_chunkIndices; }

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

  // ---- Members ----
  std::string              m_bedFile;
  bool                     m_altFirst;
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

class PlinkCursor {
public:
    PlinkCursor(const std::string& bedFile,
                uint32_t nBimLines,
                uint32_t nFamLines,
                std::vector<uint32_t> posMap,
                bool altFirst);

    // Copy constructor: copies all state but re-opens an independent file handle.
    PlinkCursor(const PlinkCursor& other);

    const uint8_t* readMarkerPtr(uint64_t gIndex);

    void getGenotypes(
        uint64_t gIndex,
        arma::vec& out,
        double& altFreq,
        double& altCounts,
        double& missingRate,
        double& hweP,
        double& maf,
        double& mac,
        std::vector<uint32_t>& indexForMissing,
        bool exactHwe = true
    );

    void beginSequentialBlock(uint64_t firstMarker);

private:
    void loadBlock(uint64_t startMarker);
    // ===== file =====
    std::string m_bedFile;
    std::ifstream m_bedStream;

    // ===== dimensions =====
    uint32_t m_nSubjInFile;
    uint32_t m_bytesPerMarker;
    uint32_t m_nMarkers;

    // ===== sample mapping =====
    std::vector<uint32_t> m_posMap;
    bool m_altFirst;

    // ===== buffer =====
    std::vector<uint8_t> m_rawBytes;
    std::vector<uint8_t> m_compactGeno;

    // ===== block buffer =====
    std::vector<uint8_t> m_blockBytes;
    uint64_t m_blockStart = 0;
    uint64_t m_blockEnd = 0;
    bool m_hasBlock = false;

    // ===== sequential read =====
    bool m_hasSeqCursor = false;
    uint64_t m_nextMarker = 0;

    // ===== constants =====
    static constexpr uint64_t BLOCK_CAPACITY = 256;
};

#endif
