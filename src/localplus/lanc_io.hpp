// lanc_io.hpp — Local-ancestry plane-separated binary format (.lanc) reader/writer
//
// The .lanc v1 format is the plane-separated successor to .abed.  Instead of
// storing 2*K dense PLINK-2-bit tracks per marker, it stores two or three
// SNP-major bit-planes under block-framed zstd, from which the consumer-visible
// per-ancestry dosage_k / hapcount_k matrices are DERIVED at decode time:
//
//   REGION A — ANCESTRY : nSamples bytes/marker, byte(i) = (hap1<<4)|(hap0&0xF)
//                         nibble in 0..K-1 => ancestry; 0xF => unassigned
//   REGION B — ALLELE   : ceil(nSamples/4) bytes/marker, 2 bits/individual
//                         bit0 = hap0 ALT, bit1 = hap1 ALT (phased)
//   REGION M — MISSING  : present iff flags.NO_MISSING==0; same layout as B
//                         bit0 = hap0 missing, bit1 = hap1 missing
//
// Each region is a run of ceil(nMarkers/blockLen) INDEPENDENT zstd frames, one
// per block of blockLen consecutive markers.  A footer index maps a marker to
// its frame for O(1) random access without a sidecar.
//
// Per-chromosome file {prefix}.chr{N}.lanc, companion {prefix}.chr{N}.bim
// (standard PLINK BIM), single shared {prefix}.fam.
//
// The public decode contract (getAllAncestries / getAdmixGenotypes) mirrors the
// old (now-removed) .abed reader's cursor contract, so the spamixlocalp
// rewiring onto LancData/LancCursor was mechanical.
#pragma once

#include "geno_factory/geno_data.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

// ── .lanc format constants ────────────────────────────────────────────

static constexpr uint8_t  LANC_MAGIC_0 = 0x4C; // 'L'
static constexpr uint8_t  LANC_MAGIC_1 = 0x41; // 'A'
static constexpr uint8_t  LANC_MAGIC_2 = 0x4E; // 'N'
static constexpr uint8_t  LANC_MAGIC_3 = 0x43; // 'C'
static constexpr uint8_t  LANC_VERSION = 1;
static constexpr uint64_t LANC_HEADER_SIZE = 24;
static constexpr uint8_t  LANC_FLAG_NO_MISSING = 0x01; // flags bit0
static constexpr uint8_t  LANC_ANC_UNASSIGNED = 0x0F;  // reserved nibble
static constexpr uint16_t LANC_DEFAULT_BLOCKLEN = 512;
static constexpr int      LANC_DEFAULT_ZSTD_LEVEL = 3;

// ── SIMD plane-unpack kernels (runtime-dispatched) ─────────────────────
//
// The decode path (LancCursor) expands each packed plane into flat per-
// individual uint8 lanes before the compare-scatter over ancestries:
//
//   unpackNibbles(in, n, lo, hi) : lo[i] = in[i] & 0x0F  (hap0 ancestry),
//                                  hi[i] = in[i] >> 4     (hap1 ancestry).
//   unpack2bit(in, n, bit0, bit1): 2-bits/individual stream ->
//                                  bit0[i] = low bit (hap0), bit1[i] = high
//                                  bit (hap1); output bytes are 0/1.  Serves
//                                  BOTH the allele plane and the missing plane.
//
// The dispatched entry points resolve the widest supported variant once via
// simdLevel() (the SPAsqr pattern).  The scalar variant compiles on every
// architecture and is the ARM/NEON path.  The per-level variants are exposed
// so the equality test (tests/lanc_simd_test.cpp) can assert byte-identity of
// each SIMD tier against the scalar reference.
namespace lanc_simd {

// Dispatched entry points (resolve to the widest supported variant).
void unpackNibbles(const uint8_t *in, int n, uint8_t *lo, uint8_t *hi);
void unpack2bit(const uint8_t *in, int n, uint8_t *bit0, uint8_t *bit1);

// Scalar variants (always available; also the ARM path).
void unpackNibbles_scalar(const uint8_t *in, int n, uint8_t *lo, uint8_t *hi);
void unpack2bit_scalar(const uint8_t *in, int n, uint8_t *bit0, uint8_t *bit1);

#if defined(__x86_64__) || defined(_M_X64)
void unpackNibbles_avx2(const uint8_t *in, int n, uint8_t *lo, uint8_t *hi);
void unpackNibbles_avx512(const uint8_t *in, int n, uint8_t *lo, uint8_t *hi);
void unpack2bit_avx2(const uint8_t *in, int n, uint8_t *bit0, uint8_t *bit1);
void unpack2bit_avx512(const uint8_t *in, int n, uint8_t *bit0, uint8_t *bit1);
#endif

// Testing helper: current runtime SIMD tier (0 = scalar, 1 = AVX2, 2 = AVX512).
int simdLevelValue();

} // namespace lanc_simd

// ======================================================================
// LancWriter — framed-zstd .lanc writer (single-threaded, deterministic)
// ======================================================================
//
// Usage:
//   LancWriter w(path, K, nSamples[, blockLen][, zstdLevel]);
//   for each marker: w.addMarker(anc0, anc1, alt0, alt1[, miss0, miss1]);
//   w.close();
//
// Frames are accumulated in memory (per-chromosome bounded) and flushed at
// close() as: 24-byte header | all A frames | all B frames | [all M frames] |
// footer index.  NO_MISSING is decided at close(): if no missing bit was set
// across all markers, the M region is omitted entirely.
class LancWriter {
  public:
    LancWriter(
        const std::string &filename,
        uint8_t K,
        uint32_t nSamples,
        uint16_t blockLen = LANC_DEFAULT_BLOCKLEN,
        int zstdLevel = LANC_DEFAULT_ZSTD_LEVEL
    );

    ~LancWriter();

    LancWriter(const LancWriter &) = delete;
    LancWriter &operator=(const LancWriter &) = delete;

    // Append one marker.  Each array has nSamples entries (one per individual):
    //   anc0[i], anc1[i] : ancestry code of hap0/hap1 in 0..K-1, or 0xF (unassigned)
    //   alt0[i], alt1[i] : ALT allele bit (0=REF, 1=ALT), phased
    //   miss0, miss1     : allele-missing bit (0/1); nullptr => no missing.
    // When a missing bit is set, the matching REGION B allele bit is written 0.
    void addMarker(
        const uint8_t *anc0,
        const uint8_t *anc1,
        const uint8_t *alt0,
        const uint8_t *alt1,
        const uint8_t *miss0 = nullptr,
        const uint8_t *miss1 = nullptr
    );

    // Finalize: flush the last partial block, write header + regions + footer.
    void close();

    uint32_t nMarkersWritten() const {
        return m_nMarkers;
    }

  private:
    void flushBlock();

    std::string m_filename;
    uint8_t     m_K;
    uint32_t    m_nSamples;
    uint16_t    m_blockLen;
    int         m_zstdLevel;

    uint64_t m_bytesPerMarkerB; // ceil(nSamples/4)
    uint32_t m_markersInBlock = 0;
    uint32_t m_nMarkers = 0;
    bool     m_anyMissing = false;
    bool     m_closed = false;

    // Raw per-region block buffers (cleared at each flush).
    std::vector<uint8_t> m_rawA;
    std::vector<uint8_t> m_rawB;
    std::vector<uint8_t> m_rawM;

    // Accumulated compressed frame blobs, one entry per block.
    std::vector<std::vector<uint8_t> > m_framesA;
    std::vector<std::vector<uint8_t> > m_framesB;
    std::vector<std::vector<uint8_t> > m_framesM;
};

// ======================================================================
// LancData — shared metadata + multi-chr marker list + cursor factory
// ======================================================================

class LancData {
  public:
    using MarkerInfo = GenoMeta::MarkerInfo;

    // Per-chromosome file record (header + footer index + provenance).
    struct ChrFileInfo {
        std::string path;
        uint8_t  K;
        bool     noMissing;
        uint32_t nSamples;
        uint32_t nMarkers;
        uint16_t blockLen;
        uint16_t zstdLevel;
        uint32_t nFrames;
        uint64_t baseA, baseB, baseM;
        uint64_t footerStart;
        std::vector<uint64_t> offA, offB, offM;
        uint64_t firstGlobalIdx; // global genoIndex of this file's first marker
    };

    // Construct from prefix: globs {prefix}.chr*.lanc, sorts by chromosome,
    // opens each header + companion .bim, assigns global genoIndex across files.
    //   usedMask: ceil(nFam/64) words, bit i set <=> subject i included.
    LancData(
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

    bool allUsed() const {
        return m_allUsed;
    }

    // True only if every chromosome file is NO_MISSING.  The cursor still
    // consults the per-file flag when decoding.
    bool hasNoMissing() const {
        return m_noMissingAll;
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

    const std::vector<ChrFileInfo> &chrFiles() const {
        return m_chrFiles;
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

    // Map a global genoIndex to its chromosome file and local marker index.
    // Returns the file index; sets localIdx = g - firstGlobalIdx.
    uint32_t fileIndexForGeno(
        uint64_t g,
        uint64_t &localIdx
    ) const;

    // Create a per-thread cursor.
    std::unique_ptr<class LancCursor> makeCursor() const;

  private:
    static std::vector<std::vector<uint64_t> > buildChunks(
        const std::vector<MarkerInfo> &markers,
        int chunkSize
    );

    uint16_t m_nAnc = 0;
    uint32_t m_nSubjInFile = 0;
    uint32_t m_nSubjUsed = 0;
    uint32_t m_nMarkers = 0;
    bool     m_allUsed = false;
    bool     m_noMissingAll = true;
    std::vector<uint64_t> m_usedMask;
    std::vector<MarkerInfo> m_markerInfo;
    std::vector<std::vector<uint64_t> > m_chunkIndices;
    std::vector<ChrFileInfo> m_chrFiles;
};

// ======================================================================
// LancCursor — per-thread, per-region framed scalar decoder
// ======================================================================

class LancCursor {
  public:
    explicit LancCursor(const LancData &data);
    ~LancCursor();

    LancCursor(const LancCursor &) = delete;
    LancCursor &operator=(const LancCursor &) = delete;

    // Hint: prepare a forward sequential run starting at a marker index.  The
    // per-region block cache makes correctness independent of this call.
    void beginSequentialBlock(uint64_t firstMarker);

    // Decode one ancestry column for the marker at the given (filtered) local
    // index.  dosageOut/hapcountOut receive N_used values.
    // Returns altFreq = sum(dosageOut)/sum(hapcountOut) (0 if hapcount sum <= 0).
    double getAdmixGenotypes(
        uint64_t markerLocalIdx,
        int ancIdx,
        Eigen::Ref<Eigen::VectorXd> dosageOut,
        Eigen::Ref<Eigen::VectorXd> hapcountOut
    );

    // Decode all K ancestries for the marker at the given (filtered) local
    // index into N_used x K matrices (column k = ancestry k).
    void getAllAncestries(
        uint64_t markerLocalIdx,
        Eigen::Ref<Eigen::MatrixXd> dosageMatrix,
        Eigen::Ref<Eigen::MatrixXd> hapcountMatrix
    );

  private:
    // Region selector for the frame cache.
    enum Region { REGION_A = 0, REGION_B = 1, REGION_M = 2 };

    struct RegionCache {
        int      fileIdx = -1;
        uint64_t block = UINT64_MAX;
        uint32_t markersInBuf = 0;
        std::vector<uint8_t> buf;
    };

    void ensureFileOpen(uint32_t fileIdx);

    // Ensure the region's frame for (fileIdx, block) is decompressed & cached.
    void ensureRegionBlock(
        Region region,
        uint32_t fileIdx,
        uint64_t block
    );

    // Resolve the marker to its plane byte pointers (ensures blocks are cached).
    // Returns the ancestry byte pointer; alleleBase/missBase point at the
    // marker's REGION B / REGION M raw bytes (missBase null if noMissing).
    const uint8_t *loadMarkerPlanes(
        uint64_t markerLocalIdx,
        const uint8_t *&alleleBase,
        const uint8_t *&missBase,
        bool &noMissing
    );

    const LancData &m_data;
    uint16_t m_K;
    uint32_t m_nSamples; // subjects in file
    uint32_t m_nUsed;
    bool m_allUsed;

    std::ifstream m_ifs;
    int m_openFileIdx = -1;

    RegionCache m_cA, m_cB, m_cM;
    std::vector<uint8_t> m_compBuf; // scratch for compressed frame reads

    // Preallocated per-cursor unpack scratch (length nSamplesInFile) so the hot
    // decode loop performs no allocation.  Filled by the dispatched SIMD kernels
    // once per marker, then read by the compare-scatter over ancestries.
    // m_mbit0/m_mbit1 are used only for files that carry a missing plane.
    std::vector<uint8_t> m_lo, m_hi;       // ancestry nibbles: hap0 (c0), hap1 (c1)
    std::vector<uint8_t> m_bit0, m_bit1;   // allele plane:  hap0 ALT, hap1 ALT
    std::vector<uint8_t> m_mbit0, m_mbit1; // missing plane: hap0 miss, hap1 miss
};
