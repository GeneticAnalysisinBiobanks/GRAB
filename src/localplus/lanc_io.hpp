// lanc_io.hpp — Local-ancestry plane-separated binary format (.lanc) reader/writer
//
// The .lanc format is the plane-separated successor to .abed.  Instead of
// storing 2*K dense PLINK-2-bit tracks per marker, it stores two or three
// SNP-major bit-planes under block-framed zstd, from which the consumer-visible
// per-ancestry dosage_k / hapcount_k matrices are DERIVED at decode time:
//
//   REGION A — ANCESTRY : nSamples bytes/marker, byte(i) = (hap1<<4)|(hap0&0xF)
//                         nibble in 0..K-1 => ancestry; 0xF => unassigned
//   REGION B — ALLELE   : ceil(nSamples/4) bytes/marker, 2 bits/individual
//                         bit0 = hap0 ALT, bit1 = hap1 ALT (phased)
//   REGION M — MISSING  : present iff a segment carries a missing bit; same
//                         layout as B (bit0 = hap0 missing, bit1 = hap1 missing)
//
// ── Merged single-file layout ──────────────────────────────────────────
//
// One {prefix}.lanc holds every chromosome, a companion merged {prefix}.bim
// lists all markers in chromosome order, and one {prefix}.fam is shared.  A
// "segment" is one chromosome's contiguous marker run; the writer flushes the
// pending zstd block at each chromosome boundary so a frame never straddles two
// chromosomes.  Little-endian throughout.
//
//   HEADER (32 bytes)
//      0  4  magic     "LANC"
//      4  1  version   1  (grab2 is alpha; the byte is bumped only at release)
//      5  1  K         1..15
//      6  1  flags     bit0 NO_MISSING_ALL (1 => no segment carries a missing plane)
//      7  1  reserved  0
//      8  4  nSamples  u32
//     12  4  nMarkers  u32  total across all segments
//     16  2  blockLen  u16  (512)
//     18  2  zstdLevel u16
//     20  4  nSeg      u32  number of chromosome segments
//     24  8  reserved  0
//
//   BODY (segments in .bim chromosome order)
//     segment 0: [A frames][B frames][M frames iff the segment has missing]
//     segment 1: ...
//     Each region = ceil(nMarkers_s / blockLen) INDEPENDENT zstd frames; the
//     last frame of each region in each segment may be partial.
//
//   FOOTER (per-segment directory, then the trailer)
//     for s in 0..nSeg-1:
//        u32 nMarkers_s
//        u8  noMissing_s
//        u8  reserved
//        u32 nFrames_s                    (= ceil(nMarkers_s / blockLen))
//        u64 baseA_s, baseB_s, baseM_s    (file-absolute; baseM_s = 0 if noMissing_s)
//        u64 offA_s[nFrames_s]
//        u64 offB_s[nFrames_s]
//        u64 offM_s[nFrames_s]            (only if !noMissing_s)
//     u64 footerStart                     (file-absolute; the LAST 8 bytes of the file)
//
// Open sequence: read the 32-byte header; seek end-8; read footerStart; seek
// there; load the nSeg segment descriptors.  A global genoIndex runs
// 0..nMarkers-1 across segments in order (and survives filterMarkersByIds).
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
static constexpr uint64_t LANC_HEADER_SIZE = 32;
static constexpr uint8_t  LANC_FLAG_NO_MISSING = 0x01; // flags bit0 (NO_MISSING_ALL)
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
// LancWriter — chromosome-aware framed-zstd .lanc writer (deterministic)
// ======================================================================
//
// Usage:
//   LancWriter w(path, K, nSamples[, blockLen][, zstdLevel]);
//   for each chromosome:
//     w.beginChromosome();
//     for each marker: w.addMarker(anc0, anc1, alt0, alt1[, miss0, miss1]);
//   w.close();
//
// One writer emits the whole merged {prefix}.lanc.  The current segment's
// compressed frames are accumulated in memory and streamed to the file (A, then
// B, then optionally M) whenever a segment is finalized — either at the next
// beginChromosome() or at close() — so peak memory stays bounded to one
// chromosome.  A 32-byte placeholder header is written up front and overwritten
// (seekp 0) at close() with the final flags/nMarkers/nSeg.  A segment is
// NO_MISSING iff no missing bit was set within it; the header NO_MISSING_ALL
// flag is set iff every segment is NO_MISSING.
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

    // Start a new chromosome segment.  If the current segment already holds >=1
    // marker it is finalized first (partial block flushed; A|B|[M] frames
    // written; footer entry recorded); then per-segment state is reset.  The
    // first call (empty current segment) simply opens segment 0.
    void beginChromosome();

    // Append one marker to the CURRENT segment.  Each array has nSamples entries
    // (one per individual):
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

    // Finalize: flush the last segment, write the footer directory, then patch
    // the 32-byte header.
    void close();

    uint32_t nMarkersWritten() const {
        return m_nMarkers;
    }

  private:
    // Per-segment footer descriptor, populated as each segment is finalized.
    struct SegDesc {
        uint32_t nMarkers = 0;
        bool     noMissing = true;
        uint32_t nFrames = 0;
        uint64_t baseA = 0, baseB = 0, baseM = 0;
        std::vector<uint64_t> offA, offB, offM;
    };

    void flushBlock();      // compress the current partial/full block
    void finalizeSegment(); // flush + write the current segment's frames
    void emit(const uint8_t *data, size_t n);

    std::string   m_filename;
    std::ofstream m_out;
    uint8_t       m_K;
    uint32_t      m_nSamples;
    uint16_t      m_blockLen;
    int           m_zstdLevel;

    uint64_t m_bytesPerMarkerB; // ceil(nSamples/4)
    uint64_t m_pos = 0;         // current file write position
    uint32_t m_nMarkers = 0;    // total markers across all segments
    bool     m_noMissingAll = true;
    bool     m_closed = false;

    // Current-segment state (reset at each finalize).
    uint32_t m_markersInBlock = 0;
    uint32_t m_segMarkers = 0;
    bool     m_segAnyMissing = false;

    // Raw per-region block buffers (cleared at each flush).
    std::vector<uint8_t> m_rawA;
    std::vector<uint8_t> m_rawB;
    std::vector<uint8_t> m_rawM;

    // Current segment's compressed frame blobs, one entry per block.
    std::vector<std::vector<uint8_t> > m_framesA;
    std::vector<std::vector<uint8_t> > m_framesB;
    std::vector<std::vector<uint8_t> > m_framesM;

    // Finalized segment descriptors (written to the footer at close()).
    std::vector<SegDesc> m_segs;
};

// ======================================================================
// LancData — shared metadata + multi-chr marker list + cursor factory
// ======================================================================

class LancData {
  public:
    using MarkerInfo = GenoMeta::MarkerInfo;

    // Per-chromosome segment descriptor (footer directory entry, file-absolute
    // offsets in the single merged {prefix}.lanc).  K / nSamples / blockLen /
    // zstdLevel are global (header) and live on LancData, not here.
    struct SegmentInfo {
        bool     noMissing;
        uint32_t nMarkers;
        uint32_t nFrames;
        uint64_t baseA, baseB, baseM;
        std::vector<uint64_t> offA, offB, offM;
        uint64_t firstGlobalIdx; // global genoIndex of this segment's first marker
    };

    // Construct from prefix: opens the single {prefix}.lanc (header + footer
    // segment directory) and the merged {prefix}.bim, assigns a global
    // genoIndex 0..nMarkers-1 across segments in chromosome order, and validates
    // that segment boundaries coincide with .bim chromosome changes.
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

    // True only if every segment is NO_MISSING (header NO_MISSING_ALL flag).
    // The cursor still consults the per-segment flag when decoding.
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

    const std::vector<SegmentInfo> &segments() const {
        return m_segs;
    }

    const std::string &lancPath() const {
        return m_lancPath;
    }

    uint16_t blockLen() const {
        return m_blockLen;
    }

    uint64_t footerStart() const {
        return m_footerStart;
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

    // Map a global genoIndex to its segment and local marker index.
    // Returns the segment index; sets localIdx = g - firstGlobalIdx.
    uint32_t segmentIndexForGeno(
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

    std::string m_lancPath;
    uint16_t m_nAnc = 0;
    uint32_t m_nSubjInFile = 0;
    uint32_t m_nSubjUsed = 0;
    uint32_t m_nMarkers = 0;
    uint16_t m_blockLen = 0;
    uint16_t m_zstdLevel = 0;
    uint64_t m_footerStart = 0;
    bool     m_allUsed = false;
    bool     m_noMissingAll = true;
    std::vector<uint64_t> m_usedMask;
    std::vector<MarkerInfo> m_markerInfo;
    std::vector<std::vector<uint64_t> > m_chunkIndices;
    std::vector<SegmentInfo> m_segs;
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
        int      segIdx = -1;
        uint64_t block = UINT64_MAX;
        uint32_t markersInBuf = 0;
        std::vector<uint8_t> buf;
    };

    void ensureOpen();

    // Ensure the region's frame for (segIdx, block) is decompressed & cached.
    void ensureRegionBlock(
        Region region,
        uint32_t segIdx,
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

    std::ifstream m_ifs; // one open handle on the single {prefix}.lanc

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
