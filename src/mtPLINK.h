#ifndef PLINK_H
#define PLINK_H

// mtPLINK.h -- PLINK binary genotype reader
//
//   PlinkData   — shared across all threads (read-only after construction).
//                 Parses .bim/.fam once, builds sample-position map,
//                 provides marker filtering and chunk building.
//
//   PlinkCursor — per-thread lightweight .bed reader.
//                 Opens its own file handle and read buffer.
//                 References shared PlinkData for metadata.

#include <RcppArmadillo.h>

#include <cstdint>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

namespace mtPLINK {

// Marker metadata from a .bim file (alleles normalized by AlleleOrder).
struct MarkerInfo {
  std::string chrom;
  uint32_t    pos;
  std::string id;
  std::string ref;
  std::string alt;
  uint64_t    genoIndex;   // 0-based line number in .bim = row in .bed
};

// Marker include/exclude filter specification.
struct MarkerFilterConfig {
  std::string IDsToIncludeFile;
  std::string RangesToIncludeFile;
  std::string IDsToExcludeFile;
  std::string RangesToExcludeFile;
};

// ======== PlinkData: shared, constructed once on main thread ========

class PlinkData {
private:
  std::string m_bedFile;
  bool        m_altFirst;            // true when AlleleOrder == "alt-first"

  uint32_t m_nSubjInFile = 0;     // total samples in .fam
  uint32_t m_nSubjUsed   = 0;     // selected samples for analysis
  uint32_t m_nMarkers       = 0;     // total markers in .bim
  uint64_t m_bytesPerMarker = 0;     // ceil(nSubjInFile / 4)

  std::vector<uint32_t> m_samplePosMap;  // model-sample i → .fam row

  // Per-marker metadata (indexed by genoIndex = BIM line number)
  std::vector<std::string> m_chr;
  std::vector<std::string> m_markerId;
  std::vector<uint32_t>    m_pos;
  std::vector<std::string> m_ref;    // ref allele (normalized by AlleleOrder)
  std::vector<std::string> m_alt;    // alt allele (normalized by AlleleOrder)

public:
  PlinkData(
    const std::string& bedFile,
    const std::string& bimFile,
    const std::string& famFile,
    const std::vector<std::string>& subjData,
    const std::string& AlleleOrder
  );

  // ---- Accessors ----
  uint32_t nMarkers()       const { return m_nMarkers; }
  uint32_t nSubjUsed()   const { return m_nSubjUsed; }
  uint32_t nSubjInFile() const { return m_nSubjInFile; }
  uint64_t bytesPerMarker() const { return m_bytesPerMarker; }
  const std::string& bedFile() const { return m_bedFile; }
  bool isAltFirst()         const { return m_altFirst; }
  const std::vector<uint32_t>& samplePosMap() const { return m_samplePosMap; }

  // Per-marker metadata by genoIndex
  const std::string& chr(uint64_t i)      const { return m_chr[i]; }
  const std::string& markerId(uint64_t i) const { return m_markerId[i]; }
  uint32_t           pos(uint64_t i)      const { return m_pos[i]; }
  const std::string& ref(uint64_t i)      const { return m_ref[i]; }
  const std::string& alt(uint64_t i)      const { return m_alt[i]; }

  // ---- Filtering & chunking (called once on main thread) ----
  std::vector<MarkerInfo> getFilteredMarkers(
    const MarkerFilterConfig& filter
  ) const;

  static std::vector<std::vector<uint64_t>> buildChunks(
    const std::vector<MarkerInfo>& markers, int chunkSize
  );
};

// ======== PlinkCursor: per-thread, lightweight ========

class PlinkCursor {
private:
  const PlinkData& m_data;
  std::ifstream    m_bedStream;

  // Block buffering
  static constexpr uint64_t BLOCK_CAPACITY = 64;
  std::vector<unsigned char> m_blockBytes;
  uint64_t m_blockStart = 0;
  uint64_t m_blockEnd   = 0;
  bool     m_hasBlock   = false;

  // Sequential cursor
  uint64_t m_nextMarker   = 0;
  bool     m_hasSeqCursor = false;

  // Single-marker raw bytes buffer
  std::vector<unsigned char> m_rawBytes;

  void loadBlock(uint64_t startMarker);
  void readMarkerBytes(uint64_t gIndex);

  static int decodeGenotype(unsigned char byte, int pos) {
    return (byte >> (pos << 1)) & 0x3;
  }

public:
  explicit PlinkCursor(const PlinkData& data);

  // Prepare for sequential reading starting at firstMarker.
  void beginSequentialBlock(uint64_t firstMarker);

  // Read genotypes for one marker.  Returns alt-allele-count vector
  // (length = nSubjUsed).  Missing positions are coded -1 and
  // listed in indexForMissing.
  arma::vec getGenotypes(
    uint64_t gIndex,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    std::vector<uint32_t>& indexForMissing
  );
};

} // namespace mtPLINK

#endif
