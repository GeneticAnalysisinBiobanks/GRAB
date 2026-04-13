// abed_io.hpp — Admixed ancestry binary genotype reader (.abed format)
//
// .abed format: 8-byte header + SNP-major 2K tracks per marker
//   Header: magic[2]{0xAD,0x4D} + version(0x02) + nAnc(1) + nSamples(4)
//   nAnc bit 7 = NO_MISSING flag; K = nAnc & 0x7F.
//   nMarkers is derived from the companion .bim file.
//
//   Each marker stores 2*K tracks: [dosage_anc0][hapcount_anc0]...[dosage_anc(K-1)][hapcount_anc(K-1)]
//   Each track is ceil(N/4) bytes using PLINK-compatible 2-bit encoding.
//   Values: 00→0, 10→1, 11→2, 01→missing  (same as PLINK BED)
//
// Shares standard .bim and .fam files with PLINK.
// Uses the same bitmask/cursor pattern as PlinkData/PlinkCursor.
#pragma once

#include "geno_factory/geno_data.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <vector>

struct BGZF; // forward declare (htslib/bgzf.h)

// ── .abed file header (8 bytes) ───────────────────────────────────────

struct AbedHeader {
    uint8_t magic[2];  // {0xAD, 0x4D}
    uint8_t version;   // 0x02
    uint8_t nAnc;      // K in bits 0-6; bit 7 = NO_MISSING flag
    uint32_t nSamples; // N (little-endian)
};

static_assert(sizeof(AbedHeader) == 8, "AbedHeader must be 8 bytes");

static constexpr uint8_t ABED_MAGIC_0 = 0xAD;
static constexpr uint8_t ABED_MAGIC_1 = 0x4D;
static constexpr uint8_t ABED_VERSION = 0x02;
static constexpr uint64_t ABED_HEADER_SIZE = 8;

// Flag encoded in nAnc high bit
static constexpr uint8_t ABED_NANC_FLAG_NO_MISSING = 0x80; // bit 7

// ======================================================================
// AdmixData: shared metadata + cursor factory
// ======================================================================

class AdmixData {
  public:
    using MarkerInfo = GenoMeta::MarkerInfo;

    // Construct from prefix: reads {prefix}.abed, {prefix}.bim, {prefix}.fam
    // usedMask: ceil(nFam/64) words, bit i set ↔ subject i included
    AdmixData(
        const std::string &prefix,
        const std::vector<uint64_t> &usedMask,
        uint32_t nFam,
        uint32_t nUsed,
        const std::string &extractFile = {},
        const std::string &excludeFile = {},
        int nMarkersEachChunk = 1024
    );

    uint16_t nAncestries() const {
        return m_nAnc;
    }

    uint32_t nMarkers() const {
        return m_nMarkers;
    }

    uint32_t nSubjUsed() const {
        return m_nSubjUsed;
    }

    uint32_t nSubjInFile() const {
        return m_nSubjInFile;
    }

    uint64_t bytesPerTrack() const {
        return m_bytesPerTrack;
    }

    uint64_t bytesPerMarker() const {
        return m_bytesPerMarker;
    }

    bool allUsed() const {
        return m_allUsed;
    }

    bool hasNoMissing() const {
        return m_noMissing;
    }

    const std::string &abedFile() const {
        return m_abedFile;
    }

    const std::vector<uint64_t> &usedMask() const {
        return m_usedMask;
    }

    const std::vector<MarkerInfo> &markerInfo() const {
        return m_markerInfo;
    }

    const std::vector<std::vector<uint64_t> > &chunkIndices() const {
        return m_chunkIndices;
    }

    std::string_view chr(uint64_t i) const {
        return m_markerInfo[i].chrom;
    }

    std::string_view markerId(uint64_t i) const {
        return m_markerInfo[i].id;
    }

    uint32_t pos(uint64_t i) const {
        return m_markerInfo[i].pos;
    }

    std::string_view ref(uint64_t i) const {
        return m_markerInfo[i].ref;
    }

    std::string_view alt(uint64_t i) const {
        return m_markerInfo[i].alt;
    }

    // Create per-thread cursor
    std::unique_ptr<class AdmixCursor> makeCursor() const;

  private:
    static std::vector<MarkerInfo>filterMarkers(
        const std::vector<MarkerInfo> &all,
        const std::string &extractFile,
        const std::string &excludeFile
    );

    static std::vector<std::vector<uint64_t> > buildChunks(
        const std::vector<MarkerInfo> &markers,
        int chunkSize
    );

    std::string m_abedFile;
    uint16_t m_nAnc;
    uint32_t m_nSubjInFile;
    uint32_t m_nSubjUsed;
    uint32_t m_nMarkers;
    uint64_t m_bytesPerTrack;  // ceil(nFam / 4)
    uint64_t m_bytesPerMarker; // 2 * nAnc * bytesPerTrack
    bool m_allUsed;
    bool m_noMissing;
    std::vector<uint64_t> m_usedMask;
    std::vector<MarkerInfo> m_markerInfo;
    std::vector<std::vector<uint64_t> > m_chunkIndices;
};

// ======================================================================
// AdmixCursor: per-thread, lightweight decoder
// ======================================================================

class AdmixCursor {
  public:
    AdmixCursor(const AdmixData &data);
    ~AdmixCursor();

    AdmixCursor(const AdmixCursor &) = delete;
    AdmixCursor &operator=(const AdmixCursor &) = delete;

    // Prepare sequential reading starting from a given marker index.
    void beginSequentialBlock(uint64_t firstMarker);

    // Decode one ancestry track for marker at genoIndex.
    //   ancIdx: 0..K-1
    //   dosageOut: receives dosage values (N_used elements)
    //   hapcountOut: receives hapcount values (N_used elements)
    //   Returns: altFreq (dosage sum / hapcount sum)
    double getAdmixGenotypes(
        uint64_t genoIndex,
        int ancIdx,
        Eigen::Ref<Eigen::VectorXd> dosageOut,
        Eigen::Ref<Eigen::VectorXd> hapcountOut
    );

    // Decode ALL ancestries at once for a given marker.
    //   dosageMatrix: N_used × K matrix, column k = dosage for ancestry k
    //   hapcountMatrix: N_used × K matrix, column k = hapcount for ancestry k
    void getAllAncestries(
        uint64_t genoIndex,
        Eigen::Ref<Eigen::MatrixXd> dosageMatrix,
        Eigen::Ref<Eigen::MatrixXd> hapcountMatrix
    );

  private:
    void loadBlock(uint64_t startMarker);

    const uint8_t *readTrackPtr(
        uint64_t genoIndex,
        int trackIdx
    );

    // Decode a single track from raw bytes into a dense vector.
    // Handles bitmask (subset) or all-used (identity) path.
    void decodeTrack(
        const uint8_t *raw,
        Eigen::Ref<Eigen::VectorXd> out
    );

    const AdmixData &m_data;
    BGZF *m_bgzf;
    uint32_t m_nUsed;
    uint64_t m_bytesPerTrack;
    uint64_t m_bytesPerMarker;
    int m_nTracks; // 2 * K

    // Per-cursor scratch buffer
    std::vector<uint8_t> m_rawBytes;

    // Block buffer for sequential reads
    std::vector<uint8_t> m_blockBytes;
    uint64_t m_blockStart = 0;
    uint64_t m_blockEnd = 0;
    bool m_hasBlock = false;
    bool m_hasSeqCursor = false;
    uint64_t m_nextMarker = 0;

    static constexpr uint64_t BLOCK_CAPACITY = 512; // markers per block
};

// ======================================================================
// .abed writer — for conversion tools
// ======================================================================

class AbedWriter {
  public:
    // Create a new .abed v2 file (8-byte header, no nMarkers).
    //   noMissing: if true, sets NO_MISSING flag in nAnc high bit
    AbedWriter(
        const std::string &filename,
        uint8_t nAnc,
        uint32_t nSamples,
        bool noMissing = false,
        int nthreads = 1
    );

    // Write one track (2-bit packed, ceil(nSamples/4) bytes).
    void writeTrack(
        const uint8_t *data,
        uint64_t nBytes
    );

    // Encode integer values (0, 1, 2 or -1 for missing) into 2-bit packed format.
    // Returns packed bytes. `values` has nSamples elements.
    static std::vector<uint8_t> encode(const std::vector<int8_t> &values);

    void close();

  private:
    BGZF *m_bgzf;
    std::string m_filename;
    uint32_t m_nSamples;
    uint8_t m_nAnc;
    bool m_noMissing;
    bool m_headerWritten;

    void ensureHeader();

};
