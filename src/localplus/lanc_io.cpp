// lanc_io.cpp — .lanc plane-separated local-ancestry reader/writer (scalar decode)

#include "localplus/lanc_io.hpp"
#include "geno_factory/variant_filter.hpp"
#include "util/logging.hpp"
#include "util/text_scanner.hpp"

#if defined(__x86_64__) || defined(_M_X64)
#  include <immintrin.h>
#endif
#include "util/simd_dispatch.hpp"

#include <zstd.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace {

// ── Little-endian scalar (de)serialization helpers ──────────────────────
// Written byte-wise so the on-disk format is host-endianness independent.

inline void putU16(std::vector<uint8_t> &v, uint16_t x) {
    v.push_back(static_cast<uint8_t>(x & 0xFF));
    v.push_back(static_cast<uint8_t>((x >> 8) & 0xFF));
}

inline void putU32(std::vector<uint8_t> &v, uint32_t x) {
    for (int i = 0; i < 4; ++i)
        v.push_back(static_cast<uint8_t>((x >> (8 * i)) & 0xFF));
}

inline void putU64(std::vector<uint8_t> &v, uint64_t x) {
    for (int i = 0; i < 8; ++i)
        v.push_back(static_cast<uint8_t>((x >> (8 * i)) & 0xFF));
}

inline uint16_t getU16(const uint8_t *p) {
    return static_cast<uint16_t>(p[0]) | (static_cast<uint16_t>(p[1]) << 8);
}

inline uint32_t getU32(const uint8_t *p) {
    uint32_t x = 0;
    for (int i = 0; i < 4; ++i)
        x |= static_cast<uint32_t>(p[i]) << (8 * i);
    return x;
}

inline uint64_t getU64(const uint8_t *p) {
    uint64_t x = 0;
    for (int i = 0; i < 8; ++i)
        x |= static_cast<uint64_t>(p[i]) << (8 * i);
    return x;
}

// ── BIM parsing (mirrors the .abed reader; standard PLINK BIM) ──────────
std::vector<GenoMeta::MarkerInfo> parseBimFile(const std::string &bimFile) {
    std::ifstream ifs(bimFile);
    if (!ifs) throw std::runtime_error("Cannot open BIM file: " + bimFile);

    std::vector<GenoMeta::MarkerInfo> markers;
    std::string line;
    uint64_t idx = 0;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        text::TokenScanner ts(line);
        std::string chrom = ts.next();
        std::string id = ts.next();
        ts.nextView(); // skip cM
        ts.skipWS();
        char *ep;
        uint32_t pos = static_cast<uint32_t>(std::strtoul(ts.pos(), &ep, 10));
        ts.p = ep;
        std::string ref = ts.next();
        std::string alt = ts.next();
        if (chrom.empty() || id.empty())
            throw std::runtime_error("Malformed BIM line: " + line);
        markers.push_back({std::move(chrom), pos, std::move(id), std::move(ref), std::move(alt), idx++});
    }
    return markers;
}

} // namespace

// ======================================================================
// SIMD plane-unpack kernels — scalar + AVX2 + AVX-512, runtime-dispatched.
//
// Each variant is byte-identical to the scalar reference for arbitrary n:
// a SIMD body over whole lanes followed by a scalar remainder tail.  The
// AVX-512 tier includes an AVX2 cleanup so it never falls back to more
// scalar work than the AVX2 tier.  Dispatch is resolved once at first use
// via a static function pointer selected by simdLevel() (SPAsqr pattern).
// ======================================================================

namespace lanc_simd {

// ── Kernel 1: nibble unpack — scalar reference ─────────────────────────
void unpackNibbles_scalar(const uint8_t *in, int n, uint8_t *lo, uint8_t *hi) {
    for (int i = 0; i < n; ++i) {
        lo[i] = static_cast<uint8_t>(in[i] & 0x0F);
        hi[i] = static_cast<uint8_t>(in[i] >> 4);
    }
}

// ── Kernel 2: 2-bit unpack — scalar reference ──────────────────────────
// Individual i lives in byte i/4 at bit positions [2*(i%4), 2*(i%4)+1].
// The final byte may hold fewer than 4 individuals; the loop bound (i < n)
// covers partial trailing bytes.
void unpack2bit_scalar(const uint8_t *in, int n, uint8_t *bit0, uint8_t *bit1) {
    for (int i = 0; i < n; ++i) {
        const uint8_t byte = in[i >> 2];
        const int shift = (i & 3) * 2;
        bit0[i] = static_cast<uint8_t>((byte >> shift) & 1);
        bit1[i] = static_cast<uint8_t>((byte >> (shift + 1)) & 1);
    }
}

#if defined(__x86_64__) || defined(_M_X64)

// ── Kernel 1: nibble unpack — AVX2 (32 bytes/iter) ─────────────────────
// Per-byte low nibble = v & 0x0F; per-byte high nibble = (v >> 4) & 0x0F.
// A 16-bit/32-bit-granular shift borrows bits across the byte boundary, but
// the subsequent per-byte AND 0x0F discards them, so an epi32 shift is exact.
__attribute__((target("avx2")))
void unpackNibbles_avx2(const uint8_t *in, int n, uint8_t *lo, uint8_t *hi) {
    const __m256i mask = _mm256_set1_epi8(0x0F);
    int i = 0;
    for (; i + 32 <= n; i += 32) {
        const __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(in + i));
        const __m256i vlo = _mm256_and_si256(v, mask);
        const __m256i vhi = _mm256_and_si256(_mm256_srli_epi32(v, 4), mask);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(lo + i), vlo);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(hi + i), vhi);
    }
    for (; i < n; ++i) {
        lo[i] = static_cast<uint8_t>(in[i] & 0x0F);
        hi[i] = static_cast<uint8_t>(in[i] >> 4);
    }
}

// ── Kernel 1: nibble unpack — AVX-512 (64 bytes/iter) ──────────────────
__attribute__((target("avx2,avx512f,avx512vl")))
void unpackNibbles_avx512(const uint8_t *in, int n, uint8_t *lo, uint8_t *hi) {
    const __m512i mask512 = _mm512_set1_epi32(0x0F0F0F0F);
    int i = 0;
    for (; i + 64 <= n; i += 64) {
        const __m512i v = _mm512_loadu_si512(reinterpret_cast<const void *>(in + i));
        const __m512i vlo = _mm512_and_si512(v, mask512);
        const __m512i vhi = _mm512_and_si512(_mm512_srli_epi32(v, 4), mask512);
        _mm512_storeu_si512(reinterpret_cast<void *>(lo + i), vlo);
        _mm512_storeu_si512(reinterpret_cast<void *>(hi + i), vhi);
    }
    // AVX2 cleanup for the 32..63 remainder.
    const __m256i mask256 = _mm256_set1_epi8(0x0F);
    for (; i + 32 <= n; i += 32) {
        const __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(in + i));
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(lo + i), _mm256_and_si256(v, mask256));
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(hi + i),
                            _mm256_and_si256(_mm256_srli_epi32(v, 4), mask256));
    }
    for (; i < n; ++i) {
        lo[i] = static_cast<uint8_t>(in[i] & 0x0F);
        hi[i] = static_cast<uint8_t>(in[i] >> 4);
    }
}

// ── Kernel 2: 2-bit unpack — AVX2 (32 individuals/iter) ────────────────
// For a 32-individual block (i%4==0) the 8 covering input bytes are broadcast
// and byte-replicated 4x via pshufb, so lane t holds in[base + t/4].  Each
// lane is then AND-tested against a per-lane single-bit mask:
//   bit0 tests bit 2*(t%4)  -> masks {0x01,0x04,0x10,0x40}
//   bit1 tests bit 2*(t%4)+1-> masks {0x02,0x08,0x20,0x80}
// cmpeq(masked, mask) yields 0xFF when the bit is set, then AND 1 -> {0,1}.
__attribute__((target("avx2")))
void unpack2bit_avx2(const uint8_t *in, int n, uint8_t *bit0, uint8_t *bit1) {
    const __m256i mask0 = _mm256_set1_epi32(0x40100401);              // bytes {1,4,16,64}
    const __m256i mask1 = _mm256_set1_epi32(static_cast<int>(0x80200802u)); // bytes {2,8,32,128}
    const __m256i one   = _mm256_set1_epi8(1);
    const __m256i ctrl  = _mm256_set_epi8(
        7, 7, 7, 7, 6, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4,
        3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0);
    int i = 0;
    for (; i + 32 <= n; i += 32) {
        const __m128i x     = _mm_loadl_epi64(reinterpret_cast<const __m128i *>(in + (i >> 2)));
        const __m256i bcast = _mm256_broadcastsi128_si256(x);
        const __m256i rep   = _mm256_shuffle_epi8(bcast, ctrl);
        const __m256i b0 =
            _mm256_and_si256(_mm256_cmpeq_epi8(_mm256_and_si256(rep, mask0), mask0), one);
        const __m256i b1 =
            _mm256_and_si256(_mm256_cmpeq_epi8(_mm256_and_si256(rep, mask1), mask1), one);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(bit0 + i), b0);
        _mm256_storeu_si256(reinterpret_cast<__m256i *>(bit1 + i), b1);
    }
    for (; i < n; ++i) {
        const uint8_t byte = in[i >> 2];
        const int shift = (i & 3) * 2;
        bit0[i] = static_cast<uint8_t>((byte >> shift) & 1);
        bit1[i] = static_cast<uint8_t>((byte >> (shift + 1)) & 1);
    }
}

// ── Kernel 2: 2-bit unpack — AVX-512 (16 individuals/lane group) ───────
// For a 16-individual group (i%4==0) the 4 covering input bytes form a 32-bit
// word W broadcast to 16 lanes.  For lane t the wanted bit sits at W-bit 2*t
// (bit0) or 2*t+1 (bit1), because 8*(t/4)+2*(t%4) == 2*t.  A per-lane variable
// shift (VPSRLVD) followed by AND 1 extracts it; VPMOVDB narrows the 16 dword
// results to 16 bytes.  The 64/iter loop unrolls this over four groups.
__attribute__((target("avx2,avx512f,avx512vl")))
void unpack2bit_avx512(const uint8_t *in, int n, uint8_t *bit0, uint8_t *bit1) {
    const __m512i sh0 = _mm512_set_epi32(30, 28, 26, 24, 22, 20, 18, 16,
                                         14, 12, 10, 8, 6, 4, 2, 0);
    const __m512i sh1 = _mm512_set_epi32(31, 29, 27, 25, 23, 21, 19, 17,
                                         15, 13, 11, 9, 7, 5, 3, 1);
    const __m512i one = _mm512_set1_epi32(1);
    int i = 0;
    for (; i + 64 <= n; i += 64) {
        for (int g = 0; g < 4; ++g) {
            uint32_t W;
            std::memcpy(&W, in + (i >> 2) + g * 4, 4);
            const __m512i vW = _mm512_set1_epi32(static_cast<int>(W));
            const __m512i v0 = _mm512_and_si512(_mm512_srlv_epi32(vW, sh0), one);
            const __m512i v1 = _mm512_and_si512(_mm512_srlv_epi32(vW, sh1), one);
            _mm_storeu_si128(reinterpret_cast<__m128i *>(bit0 + i + g * 16), _mm512_cvtepi32_epi8(v0));
            _mm_storeu_si128(reinterpret_cast<__m128i *>(bit1 + i + g * 16), _mm512_cvtepi32_epi8(v1));
        }
    }
    // 16-individual cleanup for the 16..63 remainder.
    for (; i + 16 <= n; i += 16) {
        uint32_t W;
        std::memcpy(&W, in + (i >> 2), 4);
        const __m512i vW = _mm512_set1_epi32(static_cast<int>(W));
        const __m512i v0 = _mm512_and_si512(_mm512_srlv_epi32(vW, sh0), one);
        const __m512i v1 = _mm512_and_si512(_mm512_srlv_epi32(vW, sh1), one);
        _mm_storeu_si128(reinterpret_cast<__m128i *>(bit0 + i), _mm512_cvtepi32_epi8(v0));
        _mm_storeu_si128(reinterpret_cast<__m128i *>(bit1 + i), _mm512_cvtepi32_epi8(v1));
    }
    for (; i < n; ++i) {
        const uint8_t byte = in[i >> 2];
        const int shift = (i & 3) * 2;
        bit0[i] = static_cast<uint8_t>((byte >> shift) & 1);
        bit1[i] = static_cast<uint8_t>((byte >> (shift + 1)) & 1);
    }
}

#endif // x86_64 SIMD variants

// ── Runtime dispatch: resolve the widest supported variant once ────────
namespace {

using UnpackFn = void (*)(const uint8_t *, int, uint8_t *, uint8_t *);

UnpackFn pickUnpackNibbles() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return unpackNibbles_avx512;
    case SimdLevel::AVX2:   return unpackNibbles_avx2;
    default: break;
    }
#endif
    return unpackNibbles_scalar;
}

UnpackFn pickUnpack2bit() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return unpack2bit_avx512;
    case SimdLevel::AVX2:   return unpack2bit_avx2;
    default: break;
    }
#endif
    return unpack2bit_scalar;
}

const UnpackFn g_unpackNibbles = pickUnpackNibbles();
const UnpackFn g_unpack2bit    = pickUnpack2bit();

} // namespace

void unpackNibbles(const uint8_t *in, int n, uint8_t *lo, uint8_t *hi) {
    g_unpackNibbles(in, n, lo, hi);
}

void unpack2bit(const uint8_t *in, int n, uint8_t *bit0, uint8_t *bit1) {
    g_unpack2bit(in, n, bit0, bit1);
}

int simdLevelValue() {
    return static_cast<int>(simdLevel());
}

} // namespace lanc_simd

// ======================================================================
// LancWriter
// ======================================================================

LancWriter::LancWriter(
    const std::string &filename,
    uint8_t K,
    uint32_t nSamples,
    uint16_t blockLen,
    int zstdLevel
)
    : m_filename(filename),
      m_K(K),
      m_nSamples(nSamples),
      m_blockLen(blockLen),
      m_zstdLevel(zstdLevel)
{
    if (K < 1 || K > 15)
        throw std::runtime_error("LancWriter: K must be in 1..15 (got " + std::to_string(K) + ")");
    if (blockLen == 0)
        throw std::runtime_error("LancWriter: blockLen must be > 0");
    m_bytesPerMarkerB = (static_cast<uint64_t>(nSamples) + 3) / 4;

    m_out.open(m_filename, std::ios::binary | std::ios::out | std::ios::trunc);
    if (!m_out) throw std::runtime_error("LancWriter: cannot create " + m_filename);

    // Write a 32-byte placeholder header; it is overwritten (seekp 0) at close()
    // with the final flags / nMarkers / nSeg.  Overwriting exactly these 32 bytes
    // with deterministic content keeps the output byte-reproducible.
    std::vector<uint8_t> placeholder(LANC_HEADER_SIZE, 0);
    emit(placeholder.data(), placeholder.size());
}

LancWriter::~LancWriter()
{
    if (!m_closed) {
        try {
            close();
        } catch (...) {
            // best-effort flush on destruction
        }
    }
}

void LancWriter::emit(const uint8_t *data, size_t n) {
    m_out.write(reinterpret_cast<const char *>(data), static_cast<std::streamsize>(n));
    m_pos += n;
}

void LancWriter::beginChromosome() {
    // Finalize the pending segment (if any); then start a fresh one.  The first
    // call has an empty current segment and only resets state.
    if (m_segMarkers > 0) finalizeSegment();
    m_rawA.clear();
    m_rawB.clear();
    m_rawM.clear();
    m_framesA.clear();
    m_framesB.clear();
    m_framesM.clear();
    m_markersInBlock = 0;
    m_segMarkers = 0;
    m_segAnyMissing = false;
}

void LancWriter::addMarker(
    const uint8_t *anc0,
    const uint8_t *anc1,
    const uint8_t *alt0,
    const uint8_t *alt1,
    const uint8_t *miss0,
    const uint8_t *miss1
) {
    if (m_closed)
        throw std::runtime_error("LancWriter::addMarker after close()");

    const uint32_t N = m_nSamples;

    // REGION A — one byte per individual: (hap1<<4)|(hap0&0xF)
    size_t baseA = m_rawA.size();
    m_rawA.resize(baseA + N);
    for (uint32_t i = 0; i < N; ++i)
        m_rawA[baseA + i] = static_cast<uint8_t>(((anc1[i] & 0x0F) << 4) | (anc0[i] & 0x0F));

    // REGION B (allele) and REGION M (missing) — 2 bits per individual.
    const size_t bpm = static_cast<size_t>(m_bytesPerMarkerB);
    size_t baseB = m_rawB.size();
    size_t baseM = m_rawM.size();
    m_rawB.resize(baseB + bpm, 0);
    m_rawM.resize(baseM + bpm, 0);
    for (uint32_t i = 0; i < N; ++i) {
        uint8_t m0 = miss0 ? static_cast<uint8_t>(miss0[i] & 1) : 0;
        uint8_t m1 = miss1 ? static_cast<uint8_t>(miss1[i] & 1) : 0;
        uint8_t a0 = m0 ? 0 : static_cast<uint8_t>(alt0[i] & 1); // 0 placeholder when missing
        uint8_t a1 = m1 ? 0 : static_cast<uint8_t>(alt1[i] & 1);
        uint8_t bval = static_cast<uint8_t>(a0 | (a1 << 1));
        uint8_t mval = static_cast<uint8_t>(m0 | (m1 << 1));
        if (mval) m_segAnyMissing = true;
        uint32_t shift = (i & 3) * 2;
        m_rawB[baseB + (i >> 2)] |= static_cast<uint8_t>(bval << shift);
        m_rawM[baseM + (i >> 2)] |= static_cast<uint8_t>(mval << shift);
    }

    ++m_markersInBlock;
    ++m_segMarkers;
    ++m_nMarkers;
    if (m_markersInBlock == m_blockLen)
        flushBlock();
}

void LancWriter::flushBlock() {
    if (m_markersInBlock == 0) return;

    auto compress = [this](const std::vector<uint8_t> &src, std::vector<std::vector<uint8_t> > &frames) {
        size_t bound = ZSTD_compressBound(src.size());
        std::vector<uint8_t> dst(bound);
        size_t csz = ZSTD_compress(dst.data(), bound, src.data(), src.size(), m_zstdLevel);
        if (ZSTD_isError(csz))
            throw std::runtime_error(std::string("LancWriter: ZSTD_compress failed: ") + ZSTD_getErrorName(csz));
        dst.resize(csz);
        frames.push_back(std::move(dst));
    };

    compress(m_rawA, m_framesA);
    compress(m_rawB, m_framesB);
    compress(m_rawM, m_framesM);

    m_rawA.clear();
    m_rawB.clear();
    m_rawM.clear();
    m_markersInBlock = 0;
}

void LancWriter::finalizeSegment() {
    flushBlock(); // compress the segment's trailing partial block

    const bool noMissing = !m_segAnyMissing;
    const uint32_t nFrames = static_cast<uint32_t>(m_framesA.size());
    if (m_framesB.size() != nFrames || m_framesM.size() != nFrames)
        throw std::runtime_error("LancWriter: internal frame-count mismatch");

    SegDesc seg;
    seg.nMarkers = m_segMarkers;
    seg.noMissing = noMissing;
    seg.nFrames = nFrames;

    // ── Region A frames (file-absolute offsets) ──
    seg.baseA = m_pos;
    seg.offA.resize(nFrames);
    for (uint32_t f = 0; f < nFrames; ++f) {
        seg.offA[f] = m_pos;
        emit(m_framesA[f].data(), m_framesA[f].size());
    }

    // ── Region B frames ──
    seg.baseB = m_pos;
    seg.offB.resize(nFrames);
    for (uint32_t f = 0; f < nFrames; ++f) {
        seg.offB[f] = m_pos;
        emit(m_framesB[f].data(), m_framesB[f].size());
    }

    // ── Region M frames (only if the segment saw any missing bit) ──
    seg.baseM = 0;
    if (!noMissing) {
        seg.baseM = m_pos;
        seg.offM.resize(nFrames);
        for (uint32_t f = 0; f < nFrames; ++f) {
            seg.offM[f] = m_pos;
            emit(m_framesM[f].data(), m_framesM[f].size());
        }
        m_noMissingAll = false;
    }

    m_segs.push_back(std::move(seg));

    // Release this segment's frame memory (peak stays bounded to one segment).
    m_framesA.clear();
    m_framesB.clear();
    m_framesM.clear();
}

void LancWriter::close() {
    if (m_closed) return;
    m_closed = true;

    if (m_segMarkers > 0) finalizeSegment(); // finalize the last segment

    // ── Footer: per-segment directory, then the u64 footerStart trailer ──
    const uint64_t footerStart = m_pos;
    std::vector<uint8_t> footer;
    for (const SegDesc &seg : m_segs) {
        putU32(footer, seg.nMarkers);
        footer.push_back(seg.noMissing ? 1 : 0);
        footer.push_back(0); // reserved
        putU32(footer, seg.nFrames);
        putU64(footer, seg.baseA);
        putU64(footer, seg.baseB);
        putU64(footer, seg.baseM);
        for (uint32_t f = 0; f < seg.nFrames; ++f) putU64(footer, seg.offA[f]);
        for (uint32_t f = 0; f < seg.nFrames; ++f) putU64(footer, seg.offB[f]);
        if (!seg.noMissing)
            for (uint32_t f = 0; f < seg.nFrames; ++f) putU64(footer, seg.offM[f]);
    }
    putU64(footer, footerStart); // final 8 bytes of the file
    emit(footer.data(), footer.size());

    // ── Patch the 32-byte header in place ──
    std::vector<uint8_t> header;
    header.reserve(LANC_HEADER_SIZE);
    header.push_back(LANC_MAGIC_0);
    header.push_back(LANC_MAGIC_1);
    header.push_back(LANC_MAGIC_2);
    header.push_back(LANC_MAGIC_3);
    header.push_back(LANC_VERSION);
    header.push_back(m_K);
    header.push_back(m_noMissingAll ? LANC_FLAG_NO_MISSING : 0);
    header.push_back(0); // reserved
    putU32(header, m_nSamples);
    putU32(header, m_nMarkers);
    putU16(header, m_blockLen);
    putU16(header, static_cast<uint16_t>(m_zstdLevel));
    putU32(header, static_cast<uint32_t>(m_segs.size()));
    putU32(header, 0); // reserved (upper half of the 8-byte reserved field)
    putU32(header, 0); // reserved (lower half)
    if (header.size() != LANC_HEADER_SIZE)
        throw std::runtime_error("LancWriter: header size assertion failed");

    m_out.flush();
    m_out.seekp(0, std::ios::beg);
    m_out.write(reinterpret_cast<const char *>(header.data()), static_cast<std::streamsize>(header.size()));
    m_out.flush();
    if (!m_out) throw std::runtime_error("LancWriter: write error on " + m_filename);
    m_out.close();

    m_segs.clear();
}

// ======================================================================
// LancData
// ======================================================================

LancData::LancData(
    const std::string &prefix,
    const std::vector<uint64_t> &usedMask,
    uint32_t nFam,
    uint32_t nUsed,
    const std::string &extractFile,
    const std::string &excludeFile,
    int nMarkersEachChunk
) {
    m_lancPath = prefix + ".lanc";
    const std::string bimFile = prefix + ".bim";

    m_nSubjInFile = nFam;
    m_nSubjUsed = nUsed;
    m_usedMask = usedMask;

    // ── Header (32 bytes) ──
    std::ifstream ifs(m_lancPath, std::ios::binary);
    if (!ifs) throw std::runtime_error("LancData: cannot open " + m_lancPath);

    uint8_t hbuf[LANC_HEADER_SIZE];
    ifs.read(reinterpret_cast<char *>(hbuf), LANC_HEADER_SIZE);
    if (ifs.gcount() != static_cast<std::streamsize>(LANC_HEADER_SIZE))
        throw std::runtime_error("LancData: file too short for header: " + m_lancPath);
    if (hbuf[0] != LANC_MAGIC_0 || hbuf[1] != LANC_MAGIC_1 || hbuf[2] != LANC_MAGIC_2 ||
        hbuf[3] != LANC_MAGIC_3)
        throw std::runtime_error("LancData: invalid magic in " + m_lancPath);
    if (hbuf[4] != LANC_VERSION)
        throw std::runtime_error("LancData: unsupported version " + std::to_string(hbuf[4]) + " in " + m_lancPath);

    m_nAnc = hbuf[5];
    m_noMissingAll = (hbuf[6] & LANC_FLAG_NO_MISSING) != 0;
    const uint32_t nSamples = getU32(hbuf + 8);
    const uint32_t headerMarkers = getU32(hbuf + 12);
    m_blockLen = getU16(hbuf + 16);
    m_zstdLevel = getU16(hbuf + 18);
    const uint32_t nSeg = getU32(hbuf + 20);

    if (m_nAnc < 1 || m_nAnc > 15)
        throw std::runtime_error("LancData: invalid K in " + m_lancPath);
    if (nSamples != nFam)
        throw std::runtime_error("LancData: " + m_lancPath + " has " + std::to_string(nSamples) +
                                 " samples but .fam has " + std::to_string(nFam));
    if (m_blockLen == 0)
        throw std::runtime_error("LancData: invalid blockLen in " + m_lancPath);
    if (nSeg == 0)
        throw std::runtime_error("LancData: file declares zero segments: " + m_lancPath);

    // ── Footer trailer: footerStart lives in the last 8 bytes ──
    ifs.seekg(-8, std::ios::end);
    uint8_t fbuf8[8];
    ifs.read(reinterpret_cast<char *>(fbuf8), 8);
    if (ifs.gcount() != 8)
        throw std::runtime_error("LancData: cannot read footer trailer in " + m_lancPath);
    m_footerStart = getU64(fbuf8);

    // ── Segment directory ──
    ifs.seekg(static_cast<std::streamoff>(m_footerStart), std::ios::beg);
    m_segs.reserve(nSeg);
    uint64_t globalCounter = 0;
    for (uint32_t s = 0; s < nSeg; ++s) {
        uint8_t sbuf[34]; // nMarkers(4)+noMissing(1)+reserved(1)+nFrames(4)+baseA/B/M(24)
        ifs.read(reinterpret_cast<char *>(sbuf), 34);
        if (ifs.gcount() != 34)
            throw std::runtime_error("LancData: cannot read segment descriptor " + std::to_string(s) +
                                     " in " + m_lancPath);
        SegmentInfo seg{};
        seg.nMarkers = getU32(sbuf + 0);
        seg.noMissing = sbuf[4] != 0;
        // sbuf[5] reserved
        seg.nFrames = getU32(sbuf + 6);
        seg.baseA = getU64(sbuf + 10);
        seg.baseB = getU64(sbuf + 18);
        seg.baseM = getU64(sbuf + 26);

        const uint32_t expFrames =
            static_cast<uint32_t>((static_cast<uint64_t>(seg.nMarkers) + m_blockLen - 1) / m_blockLen);
        if (seg.nFrames != expFrames)
            throw std::runtime_error("LancData: frame-count mismatch in segment " + std::to_string(s) +
                                     " of " + m_lancPath);

        const int nTables = seg.noMissing ? 2 : 3;
        std::vector<uint8_t> tbuf(static_cast<size_t>(nTables) * seg.nFrames * 8);
        ifs.read(reinterpret_cast<char *>(tbuf.data()), static_cast<std::streamsize>(tbuf.size()));
        if (ifs.gcount() != static_cast<std::streamsize>(tbuf.size()))
            throw std::runtime_error("LancData: cannot read frame tables for segment " + std::to_string(s) +
                                     " in " + m_lancPath);
        seg.offA.resize(seg.nFrames);
        seg.offB.resize(seg.nFrames);
        size_t p = 0;
        for (uint32_t f = 0; f < seg.nFrames; ++f, p += 8) seg.offA[f] = getU64(tbuf.data() + p);
        for (uint32_t f = 0; f < seg.nFrames; ++f, p += 8) seg.offB[f] = getU64(tbuf.data() + p);
        if (!seg.noMissing) {
            seg.offM.resize(seg.nFrames);
            for (uint32_t f = 0; f < seg.nFrames; ++f, p += 8) seg.offM[f] = getU64(tbuf.data() + p);
        }

        seg.firstGlobalIdx = globalCounter;
        globalCounter += seg.nMarkers;
        m_segs.push_back(std::move(seg));
    }
    ifs.close();

    // ── Σ nMarkers_s == header nMarkers ──
    if (globalCounter != headerMarkers)
        throw std::runtime_error("LancData: segment marker counts (" + std::to_string(globalCounter) +
                                 ") != header nMarkers (" + std::to_string(headerMarkers) + ") in " + m_lancPath);

    // ── Merged .bim (all markers, chromosome order) ──
    m_markerInfo = parseBimFile(bimFile);
    if (m_markerInfo.size() != headerMarkers)
        throw std::runtime_error("LancData: .bim rows (" + std::to_string(m_markerInfo.size()) +
                                 ") != header nMarkers (" + std::to_string(headerMarkers) + ") for " + bimFile);
    for (uint64_t i = 0; i < m_markerInfo.size(); ++i) m_markerInfo[i].genoIndex = i;

    // ── Validate segment boundaries coincide with .bim chromosome changes ──
    //   Each segment's nMarkers_s must equal a run of consecutive .bim rows that
    //   all share one chromosome, and adjacent segments must differ in chrom.
    {
        size_t rowStart = 0;
        for (uint32_t s = 0; s < nSeg; ++s) {
            const uint32_t nm = m_segs[s].nMarkers;
            if (rowStart + nm > m_markerInfo.size())
                throw std::runtime_error("LancData: segment " + std::to_string(s) +
                                         " overruns the .bim marker list");
            const std::string &segChrom = m_markerInfo[rowStart].chrom;
            for (size_t i = rowStart; i < rowStart + nm; ++i)
                if (m_markerInfo[i].chrom != segChrom)
                    throw std::runtime_error("LancData: segment " + std::to_string(s) +
                                             " spans multiple chromosomes in " + bimFile + " ('" + segChrom +
                                             "' vs '" + m_markerInfo[i].chrom + "')");
            if (rowStart > 0 && m_markerInfo[rowStart - 1].chrom == segChrom)
                throw std::runtime_error("LancData: segment boundary at .bim row " + std::to_string(rowStart) +
                                         " does not align with a chromosome change (chrom '" + segChrom +
                                         "' continues from the previous segment) in " + bimFile);
            rowStart += nm;
        }
        if (rowStart != m_markerInfo.size())
            throw std::runtime_error("LancData: segment marker counts do not cover the .bim rows in " + bimFile);
    }

    // ── Apply --extract / --exclude to the global marker list ──
    geno_factory::filterMarkersByIds(m_markerInfo, extractFile, excludeFile);
    m_nMarkers = static_cast<uint32_t>(m_markerInfo.size());

    // ── Chunks (break at chromosome boundaries) ──
    m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);

    // ── all-used fast path detection ──
    uint64_t maskBits = 0;
    for (auto w : m_usedMask) maskBits += static_cast<uint64_t>(__builtin_popcountll(w));
    m_allUsed = (maskBits == m_nSubjInFile);

    infoMsg("Loaded LANC: %u ancestries x %u markers x %u samples (%u used), %zu segments",
            m_nAnc, m_nMarkers, m_nSubjInFile, m_nSubjUsed, m_segs.size());
}

std::vector<std::vector<uint64_t> > LancData::buildChunks(
    const std::vector<MarkerInfo> &markers,
    int chunkSize
) {
    std::vector<std::vector<uint64_t> > chunks;
    if (markers.empty()) return chunks;

    std::vector<uint64_t> cur;
    cur.reserve(chunkSize);
    std::string prevChr = markers[0].chrom;

    for (size_t i = 0; i < markers.size(); ++i) {
        if (markers[i].chrom != prevChr) {
            if (!cur.empty()) {
                chunks.push_back(std::move(cur));
                cur.clear();
                cur.reserve(chunkSize);
            }
            prevChr = markers[i].chrom;
        }
        cur.push_back(i);
        if (static_cast<int>(cur.size()) >= chunkSize) {
            chunks.push_back(std::move(cur));
            cur.clear();
            cur.reserve(chunkSize);
        }
    }
    if (!cur.empty()) chunks.push_back(std::move(cur));
    return chunks;
}

uint32_t LancData::segmentIndexForGeno(
    uint64_t g,
    uint64_t &localIdx
) const {
    // m_segs is ordered by firstGlobalIdx (contiguous, ascending).
    size_t lo = 0, hi = m_segs.size();
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (m_segs[mid].firstGlobalIdx <= g)
            lo = mid + 1;
        else
            hi = mid;
    }
    if (lo == 0)
        throw std::runtime_error("LancData: genoIndex out of range");
    uint32_t idx = static_cast<uint32_t>(lo - 1);
    localIdx = g - m_segs[idx].firstGlobalIdx;
    return idx;
}

std::unique_ptr<LancCursor> LancData::makeCursor() const {
    return std::make_unique<LancCursor>(*this);
}

// ======================================================================
// LancCursor
// ======================================================================

LancCursor::LancCursor(const LancData &data)
    : m_data(data),
      m_K(data.nAncestries()),
      m_nSamples(data.nSubjInFile()),
      m_nUsed(data.nSubjUsed()),
      m_allUsed(data.allUsed())
{
    // Preallocate the SIMD unpack scratch (one byte per in-file individual per
    // plane lane) so the per-marker decode never allocates.
    m_lo.resize(m_nSamples);
    m_hi.resize(m_nSamples);
    m_bit0.resize(m_nSamples);
    m_bit1.resize(m_nSamples);
    if (!m_data.hasNoMissing()) {
        m_mbit0.resize(m_nSamples);
        m_mbit1.resize(m_nSamples);
    }
}

LancCursor::~LancCursor()
{
    if (m_ifs.is_open()) m_ifs.close();
}

void LancCursor::beginSequentialBlock(uint64_t /*firstMarker*/) {
    // The per-region (segIdx, block) cache makes decode correctness independent
    // of this hint; nothing further is required in the scalar path.
}

void LancCursor::ensureOpen() {
    if (m_ifs.is_open()) return;
    m_ifs.clear();
    m_ifs.open(m_data.lancPath(), std::ios::binary);
    if (!m_ifs)
        throw std::runtime_error("LancCursor: cannot open " + m_data.lancPath());
}

void LancCursor::ensureRegionBlock(
    Region region,
    uint32_t segIdx,
    uint64_t block
) {
    RegionCache &c = (region == REGION_A) ? m_cA : (region == REGION_B) ? m_cB : m_cM;
    if (c.segIdx == static_cast<int>(segIdx) && c.block == block) return; // cache hit

    const auto &segs = m_data.segments();
    const LancData::SegmentInfo &seg = segs[segIdx];
    const uint16_t blockLen = m_data.blockLen();

    // segEnd_s: the file-absolute end of this segment's compressed body — the
    // next segment's baseA, or footerStart for the last segment.  It bounds the
    // last frame of the final region present in the segment.
    const uint64_t segEnd =
        (segIdx + 1 < segs.size()) ? segs[segIdx + 1].baseA : m_data.footerStart();

    const std::vector<uint64_t> *offs = nullptr;
    uint64_t fallbackEnd = 0;
    uint64_t rawPerMarker = 0;
    switch (region) {
    case REGION_A:
        offs = &seg.offA;
        fallbackEnd = seg.baseB;
        rawPerMarker = m_nSamples;
        break;
    case REGION_B:
        offs = &seg.offB;
        fallbackEnd = seg.noMissing ? segEnd : seg.baseM;
        rawPerMarker = (static_cast<uint64_t>(m_nSamples) + 3) / 4;
        break;
    case REGION_M:
        offs = &seg.offM;
        fallbackEnd = segEnd;
        rawPerMarker = (static_cast<uint64_t>(m_nSamples) + 3) / 4;
        break;
    }

    const uint64_t start = (*offs)[block];
    const uint64_t end = (block + 1 < seg.nFrames) ? (*offs)[block + 1] : fallbackEnd;
    if (end < start)
        throw std::runtime_error("LancCursor: corrupt frame offsets");
    const size_t compressedSize = static_cast<size_t>(end - start);

    const uint64_t markersInBuf =
        std::min<uint64_t>(blockLen, seg.nMarkers - block * blockLen);
    const size_t rawSize = static_cast<size_t>(markersInBuf * rawPerMarker);

    ensureOpen();
    m_compBuf.resize(compressedSize);
    m_ifs.clear();
    m_ifs.seekg(static_cast<std::streamoff>(start), std::ios::beg);
    m_ifs.read(reinterpret_cast<char *>(m_compBuf.data()), static_cast<std::streamsize>(compressedSize));
    if (m_ifs.gcount() != static_cast<std::streamsize>(compressedSize))
        throw std::runtime_error("LancCursor: short read on frame in " + m_data.lancPath());

    c.buf.resize(rawSize);
    size_t dsz = ZSTD_decompress(c.buf.data(), rawSize, m_compBuf.data(), compressedSize);
    if (ZSTD_isError(dsz))
        throw std::runtime_error(std::string("LancCursor: ZSTD_decompress failed: ") + ZSTD_getErrorName(dsz));
    if (dsz != rawSize)
        throw std::runtime_error("LancCursor: decompressed size mismatch in " + m_data.lancPath());

    c.segIdx = static_cast<int>(segIdx);
    c.block = block;
    c.markersInBuf = static_cast<uint32_t>(markersInBuf);
}

const uint8_t *LancCursor::loadMarkerPlanes(
    uint64_t markerLocalIdx,
    const uint8_t *&alleleBase,
    const uint8_t *&missBase,
    bool &noMissing
) {
    const uint64_t g = m_data.markerInfo()[markerLocalIdx].genoIndex;
    uint64_t localIdx = 0;
    const uint32_t segIdx = m_data.segmentIndexForGeno(g, localIdx);
    const LancData::SegmentInfo &seg = m_data.segments()[segIdx];
    const uint16_t blockLen = m_data.blockLen();
    const uint64_t block = localIdx / blockLen;
    const uint64_t within = localIdx % blockLen;
    noMissing = seg.noMissing;

    ensureRegionBlock(REGION_A, segIdx, block);
    ensureRegionBlock(REGION_B, segIdx, block);
    if (!noMissing) ensureRegionBlock(REGION_M, segIdx, block);

    const uint64_t bpmB = (static_cast<uint64_t>(m_nSamples) + 3) / 4;
    const uint8_t *ancByte = m_cA.buf.data() + within * m_nSamples;
    alleleBase = m_cB.buf.data() + within * bpmB;
    missBase = noMissing ? nullptr : (m_cM.buf.data() + within * bpmB);
    return ancByte;
}

void LancCursor::getAllAncestries(
    uint64_t markerLocalIdx,
    Eigen::Ref<Eigen::MatrixXd> dosageMatrix,
    Eigen::Ref<Eigen::MatrixXd> hapcountMatrix
) {
    const uint8_t *alleleBase = nullptr;
    const uint8_t *missBase = nullptr;
    bool noMissing = true;
    const uint8_t *ancByte = loadMarkerPlanes(markerLocalIdx, alleleBase, missBase, noMissing);

    const int N = static_cast<int>(m_nSamples);
    const int K = static_cast<int>(m_K);
    const double NaN = std::numeric_limits<double>::quiet_NaN();

    // Unpack each packed plane into flat per-individual lanes once; the
    // compare-scatter below then reads the plain uint8 arrays (Plan A §2.1).
    lanc_simd::unpackNibbles(ancByte, N, m_lo.data(), m_hi.data());
    lanc_simd::unpack2bit(alleleBase, N, m_bit0.data(), m_bit1.data());
    if (!noMissing)
        lanc_simd::unpack2bit(missBase, N, m_mbit0.data(), m_mbit1.data());

    auto decode = [&](uint32_t s, uint32_t row) {
        const int c0 = m_lo[s];
        const int c1 = m_hi[s];
        const uint8_t a0 = m_bit0[s];
        const uint8_t a1 = m_bit1[s];
        uint8_t m0 = 0, m1 = 0;
        if (!noMissing) {
            m0 = m_mbit0[s];
            m1 = m_mbit1[s];
        }
        for (int k = 0; k < K; ++k) {
            dosageMatrix(row, k) = 0.0;
            hapcountMatrix(row, k) = 0.0;
        }
        if (c0 < K) {
            hapcountMatrix(row, c0) += 1.0;
            dosageMatrix(row, c0) += (m0 ? NaN : static_cast<double>(a0));
        }
        if (c1 < K) {
            hapcountMatrix(row, c1) += 1.0;
            dosageMatrix(row, c1) += (m1 ? NaN : static_cast<double>(a1));
        }
    };

    if (m_allUsed) {
        for (uint32_t s = 0; s < m_nSamples; ++s) decode(s, s);
    } else {
        const auto &mask = m_data.usedMask();
        uint32_t di = 0;
        for (uint32_t w = 0; w < static_cast<uint32_t>(mask.size()); ++w) {
            uint64_t bits = mask[w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                decode(w * 64 + static_cast<uint32_t>(bit), di++);
                bits &= bits - 1;
            }
        }
    }
}

double LancCursor::getAdmixGenotypes(
    uint64_t markerLocalIdx,
    int ancIdx,
    Eigen::Ref<Eigen::VectorXd> dosageOut,
    Eigen::Ref<Eigen::VectorXd> hapcountOut
) {
    const uint8_t *alleleBase = nullptr;
    const uint8_t *missBase = nullptr;
    bool noMissing = true;
    const uint8_t *ancByte = loadMarkerPlanes(markerLocalIdx, alleleBase, missBase, noMissing);

    const int N = static_cast<int>(m_nSamples);
    const double NaN = std::numeric_limits<double>::quiet_NaN();

    lanc_simd::unpackNibbles(ancByte, N, m_lo.data(), m_hi.data());
    lanc_simd::unpack2bit(alleleBase, N, m_bit0.data(), m_bit1.data());
    if (!noMissing)
        lanc_simd::unpack2bit(missBase, N, m_mbit0.data(), m_mbit1.data());

    auto decode = [&](uint32_t s, uint32_t row) {
        const int c0 = m_lo[s];
        const int c1 = m_hi[s];
        const uint8_t a0 = m_bit0[s];
        const uint8_t a1 = m_bit1[s];
        uint8_t m0 = 0, m1 = 0;
        if (!noMissing) {
            m0 = m_mbit0[s];
            m1 = m_mbit1[s];
        }
        double d = 0.0, h = 0.0;
        if (c0 == ancIdx) {
            h += 1.0;
            d += (m0 ? NaN : static_cast<double>(a0));
        }
        if (c1 == ancIdx) {
            h += 1.0;
            d += (m1 ? NaN : static_cast<double>(a1));
        }
        dosageOut(row) = d;
        hapcountOut(row) = h;
    };

    if (m_allUsed) {
        for (uint32_t s = 0; s < m_nSamples; ++s) decode(s, s);
    } else {
        const auto &mask = m_data.usedMask();
        uint32_t di = 0;
        for (uint32_t w = 0; w < static_cast<uint32_t>(mask.size()); ++w) {
            uint64_t bits = mask[w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                decode(w * 64 + static_cast<uint32_t>(bit), di++);
                bits &= bits - 1;
            }
        }
    }

    const double dosSum = dosageOut.sum();
    const double hapSum = hapcountOut.sum();
    return (hapSum > 0.0) ? dosSum / hapSum : 0.0;
}
