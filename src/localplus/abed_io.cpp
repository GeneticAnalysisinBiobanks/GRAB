// abed_io.cpp — Admixed ancestry binary genotype reader/writer (.abed format)

#include "localplus/abed_io.hpp"
#include "io/subject_data.hpp" // parseFamIIDs, parseBimLines
#include "util/logging.hpp"
#include "util/text_scanner.hpp"

extern "C" {
#include <htslib/bgzf.h>
}

#include <immintrin.h>
#include "util/simd_dispatch.hpp"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <unordered_set>

// ── BIM parsing (shared with plink.hpp) ─────────────────────────────

static std::vector<GenoMeta::MarkerInfo> parseBimFile(const std::string &bimFile) {
    std::ifstream ifs(bimFile);
    if (!ifs) throw std::runtime_error("Cannot open BIM file: " + bimFile);

    std::vector<GenoMeta::MarkerInfo> markers;
    std::string line;
    uint64_t idx = 0;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        text::TokenScanner ts(line);
        std::string chrom = ts.next();
        std::string id    = ts.next();
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

// ======================================================================
// AdmixData
// ======================================================================

AdmixData::AdmixData(
    const std::string &prefix,
    const std::vector<uint64_t> &usedMask,
    uint32_t nFam,
    uint32_t nUsed,
    const std::string &extractFile,
    const std::string &excludeFile,
    int nMarkersEachChunk
)
{
    m_abedFile = prefix + ".abed";
    std::string bimFile = prefix + ".bim";
    std::string famFile = prefix + ".fam";

    // Read and validate 8-byte header
    BGZF *bgz = bgzf_open(m_abedFile.c_str(), "r");
    if (!bgz) throw std::runtime_error("Cannot open ABED file: " + m_abedFile);

    AbedHeader hdr{};
    if (bgzf_read(bgz, &hdr, sizeof(hdr)) != sizeof(hdr)) {
        bgzf_close(bgz);
        throw std::runtime_error("ABED file too short for header: " + m_abedFile);
    }
    bgzf_close(bgz);

    if (hdr.magic[0] != ABED_MAGIC_0 || hdr.magic[1] != ABED_MAGIC_1)
        throw std::runtime_error("Invalid ABED magic bytes in " + m_abedFile);
    if (hdr.version != ABED_VERSION)
        throw std::runtime_error("Unsupported ABED version " + std::to_string(hdr.version) + " in " + m_abedFile);

    m_nAnc = hdr.nAnc & 0x7F;
    m_noMissing = (hdr.nAnc & ABED_NANC_FLAG_NO_MISSING) != 0;
    m_nSubjInFile = hdr.nSamples;

    if (m_nSubjInFile != nFam)
        throw std::runtime_error("ABED header says " + std::to_string(m_nSubjInFile) + " samples but .fam has " +
                                 std::to_string(nFam));

    m_nSubjUsed = nUsed;
    m_usedMask = usedMask;
    m_bytesPerTrack = (static_cast<uint64_t>(m_nSubjInFile) + 3) / 4;
    m_bytesPerMarker = 2ULL * m_nAnc * m_bytesPerTrack;

    // Check all used
    uint64_t maskBits = 0;
    for (auto w : m_usedMask)
        maskBits += __builtin_popcountll(w);
    m_allUsed = (maskBits == m_nSubjInFile);

    // Parse .bim — marker count derived from .bim, not header
    auto allMarkers = parseBimFile(bimFile);

    // Apply extract/exclude filters
    m_markerInfo = filterMarkers(allMarkers, extractFile, excludeFile);
    m_nMarkers = static_cast<uint32_t>(m_markerInfo.size());

    // Build chunks
    m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);

    infoMsg("Loaded ABED: %u ancestries x %u markers x %u samples (%u used)", m_nAnc, m_nMarkers, m_nSubjInFile,
            m_nSubjUsed);
}

std::vector<GenoMeta::MarkerInfo> AdmixData::filterMarkers(
    const std::vector<MarkerInfo> &all,
    const std::string &extractFile,
    const std::string &excludeFile
) {
    std::unordered_set<std::string> includeSet, excludeSet;

    auto readIDFile = [](const std::string &path, std::unordered_set<std::string> &out) {
        std::ifstream ifs(path);
        if (!ifs) throw std::runtime_error("Cannot open SNP list file: " + path);
        std::string id;
        while (ifs >> id) {
            if (!id.empty() && id[0] != '#') out.insert(std::move(id));
        }
    };

    if (!extractFile.empty()) readIDFile(extractFile, includeSet);
    if (!excludeFile.empty()) readIDFile(excludeFile, excludeSet);

    if (includeSet.empty() && excludeSet.empty()) return all;

    std::vector<MarkerInfo> filtered;
    filtered.reserve(all.size());
    for (const auto &m : all) {
        if (!includeSet.empty() && includeSet.find(m.id) == includeSet.end()) continue;
        if (!excludeSet.empty() && excludeSet.find(m.id) != excludeSet.end()) continue;
        filtered.push_back(m);
    }
    return filtered;
}

std::vector<std::vector<uint64_t> > AdmixData::buildChunks(
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

std::unique_ptr<AdmixCursor> AdmixData::makeCursor() const {
    return std::make_unique<AdmixCursor>(*this);
}

// ======================================================================
// AdmixCursor
// ======================================================================

AdmixCursor::AdmixCursor(const AdmixData &data)
    : m_data(data),
      m_bgzf(nullptr),
      m_nUsed(data.nSubjUsed()),
      m_bytesPerTrack(data.bytesPerTrack()),
      m_bytesPerMarker(data.bytesPerMarker()),
      m_nTracks(2 * data.nAncestries())
{
    m_bgzf = bgzf_open(data.abedFile().c_str(), "r");
    if (!m_bgzf) throw std::runtime_error("AdmixCursor: cannot open " + data.abedFile());

    // Load .gzi index for random access via bgzf_useek (required for BGZF files;
    // silently ignored for uncompressed legacy .abed files)
    int idxRc = bgzf_index_load(m_bgzf, data.abedFile().c_str(), ".gzi");
    (void)idxRc;

    m_rawBytes.resize(m_bytesPerTrack);
}

AdmixCursor::~AdmixCursor()
{
    if (m_bgzf) bgzf_close(m_bgzf);
}

void AdmixCursor::beginSequentialBlock(uint64_t firstMarker) {
    m_hasSeqCursor = true;
    m_nextMarker = firstMarker;
    m_hasBlock = false;
}

void AdmixCursor::loadBlock(uint64_t startMarker) {
    uint64_t nTotal = m_data.markerInfo().size();
    uint64_t end = std::min(startMarker + BLOCK_CAPACITY, nTotal);
    uint64_t count = end - startMarker;
    if (count == 0) {
        m_hasBlock = false;
        return;
    }

    // Resolve genoIndex of first marker for file offset
    uint64_t firstGenoIdx = m_data.markerInfo()[startMarker].genoIndex;
    uint64_t offset = ABED_HEADER_SIZE + firstGenoIdx * m_bytesPerMarker;

    m_blockBytes.resize(count * m_bytesPerMarker);
    if (bgzf_useek(m_bgzf, static_cast<off_t>(offset), SEEK_SET) < 0)
        throw std::runtime_error("AdmixCursor: bgzf_useek failed in loadBlock");
    if (bgzf_read(m_bgzf, m_blockBytes.data(), count * m_bytesPerMarker) !=
        static_cast<ssize_t>(count * m_bytesPerMarker))
        throw std::runtime_error("AdmixCursor: bgzf_read failed in loadBlock");

    m_blockStart = startMarker;
    m_blockEnd = end;
    m_hasBlock = true;
}

const uint8_t *AdmixCursor::readTrackPtr(
    uint64_t markerLocalIdx,
    int trackIdx
) {
    // Check block cache
    if (m_hasBlock && markerLocalIdx >= m_blockStart && markerLocalIdx < m_blockEnd) {
        // Check if genoIndices are contiguous in the block
        uint64_t genoIdx = m_data.markerInfo()[markerLocalIdx].genoIndex;
        uint64_t baseGenoIdx = m_data.markerInfo()[m_blockStart].genoIndex;
        uint64_t relMarker = genoIdx - baseGenoIdx;
        return m_blockBytes.data() + relMarker * m_bytesPerMarker + trackIdx * m_bytesPerTrack;
    }

    // Try loading a block if sequential
    if (m_hasSeqCursor && markerLocalIdx == m_nextMarker) {
        loadBlock(markerLocalIdx);
        m_nextMarker = m_blockEnd;
        if (m_hasBlock) {
            return m_blockBytes.data() + trackIdx * m_bytesPerTrack;
        }
    }

    // Fallback: single-marker read
    uint64_t genoIdx = m_data.markerInfo()[markerLocalIdx].genoIndex;
    uint64_t offset = ABED_HEADER_SIZE + genoIdx * m_bytesPerMarker + trackIdx * m_bytesPerTrack;
    if (bgzf_useek(m_bgzf, static_cast<off_t>(offset), SEEK_SET) < 0)
        throw std::runtime_error("AdmixCursor: bgzf_useek failed");
    if (bgzf_read(m_bgzf, m_rawBytes.data(), m_bytesPerTrack) != static_cast<ssize_t>(m_bytesPerTrack))
        throw std::runtime_error("AdmixCursor: bgzf_read failed");
    return m_rawBytes.data();
}

// ── AVX2 helper for no-missing 2-bit decode ─────────────────────────
// Returns number of full bytes decoded (caller handles the scalar tail).
__attribute__((target("avx2")))
static uint32_t decodeNoMissAVX2(
    const uint8_t *raw,
    double *p,
    uint32_t fullBytes
) {
    const __m128i vone = _mm_set1_epi32(1);
    const __m128i vzero = _mm_setzero_si128();
    uint32_t di = 0, b = 0;

    for (; b + 7 < fullBytes; b += 8) {
        for (int bi = 0; bi < 8; ++bi) {
            uint8_t byte = raw[b + bi];
            __m128i vc = _mm_setr_epi32((byte) & 3, (byte >> 2) & 3, (byte >> 4) & 3, (byte >> 6) & 3);
            __m128i vi = _mm_max_epi32(_mm_sub_epi32(vc, vone), vzero);
            _mm256_storeu_pd(p + di, _mm256_cvtepi32_pd(vi));
            di += 4;
        }
    }
    return b;
}

void AdmixCursor::decodeTrack(
    const uint8_t *raw,
    Eigen::Ref<Eigen::VectorXd> out
) {
    // Lookup table: 00→0, 10→1, 11→2, 01→missing(NaN)
    static const double kDecode[4] = {0.0, std::numeric_limits<double>::quiet_NaN(), 1.0, 2.0};
    // No-missing variant: 01→0 (treat as ref) — avoids NaN in the output
    static const double kDecodeNoMiss[4] = {0.0, 0.0, 1.0, 2.0};

    const bool noMissing = m_data.hasNoMissing();

    if (m_data.allUsed()) {
        uint32_t n = m_data.nSubjInFile();
        double *__restrict p = out.data();

        if (noMissing) {
            // 2-bit code mapping: 0→0, 1→0, 2→1, 3→2 = max(0, code - 1)
            uint32_t fullBytes = n / 4;
            uint32_t di = 0;
            uint32_t b = 0;

            if (simdLevel() >= SimdLevel::AVX2) {
                b = decodeNoMissAVX2(raw, p, fullBytes);
                di = b * 4;
            }

            // Scalar tail
            for (; b < fullBytes; ++b) {
                uint8_t byte = raw[b];
                p[di++] = kDecodeNoMiss[(byte) & 3];
                p[di++] = kDecodeNoMiss[(byte >> 2) & 3];
                p[di++] = kDecodeNoMiss[(byte >> 4) & 3];
                p[di++] = kDecodeNoMiss[(byte >> 6) & 3];
            }
            uint32_t rem = n % 4;
            if (rem > 0) {
                uint8_t byte = raw[fullBytes];
                for (uint32_t r = 0; r < rem; ++r)
                    p[di++] = kDecodeNoMiss[(byte >> (2 * r)) & 3];
            }
            return;
        }

        // ── Scalar path (with missing values) ──
        uint32_t fullBytes = n / 4;
        uint32_t rem = n % 4;
        uint32_t di = 0;
        for (uint32_t b = 0; b < fullBytes; ++b) {
            uint8_t byte = raw[b];
            p[di++] = kDecode[(byte) & 3];
            p[di++] = kDecode[(byte >> 2) & 3];
            p[di++] = kDecode[(byte >> 4) & 3];
            p[di++] = kDecode[(byte >> 6) & 3];
        }
        if (rem > 0) {
            uint8_t byte = raw[fullBytes];
            for (uint32_t r = 0; r < rem; ++r)
                p[di++] = kDecode[(byte >> (2 * r)) & 3];
        }
    } else {
        // Bitmask path: decode only used subjects
        const double *lut = noMissing ? kDecodeNoMiss : kDecode;
        const auto &mask = m_data.usedMask();
        uint32_t nWords = static_cast<uint32_t>(mask.size());
        uint32_t di = 0;
        for (uint32_t w = 0; w < nWords; ++w) {
            uint64_t bits = mask[w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                uint32_t subjIdx = w * 64 + bit;
                uint32_t byteIdx = subjIdx / 4;
                uint32_t shift = (subjIdx % 4) * 2;
                uint8_t code = (raw[byteIdx] >> shift) & 3;
                out[di++] = lut[code];
                bits &= bits - 1; // clear lowest set bit
            }
        }
    }
}

double AdmixCursor::getAdmixGenotypes(
    uint64_t markerLocalIdx,
    int ancIdx,
    Eigen::Ref<Eigen::VectorXd> dosageOut,
    Eigen::Ref<Eigen::VectorXd> hapcountOut
) {
    int dosageTrack = 2 * ancIdx;
    int hapTrack = 2 * ancIdx + 1;

    const uint8_t *dosageRaw = readTrackPtr(markerLocalIdx, dosageTrack);
    decodeTrack(dosageRaw, dosageOut);

    const uint8_t *hapRaw = readTrackPtr(markerLocalIdx, hapTrack);
    decodeTrack(hapRaw, hapcountOut);

    // Compute allele frequency: sum(dosage) / sum(hapcount)
    double dosSum = 0.0, hapSum = 0.0;
    if (m_data.hasNoMissing()) {
        dosSum = dosageOut.sum();
        hapSum = hapcountOut.sum();
    } else {
        for (int i = 0; i < static_cast<int>(m_nUsed); ++i) {
            if (std::isfinite(dosageOut[i]) && std::isfinite(hapcountOut[i])) {
                dosSum += dosageOut[i];
                hapSum += hapcountOut[i];
            }
        }
    }
    return (hapSum > 0.0) ? dosSum / hapSum : 0.0;
}

void AdmixCursor::getAllAncestries(
    uint64_t markerLocalIdx,
    Eigen::Ref<Eigen::MatrixXd> dosageMatrix,
    Eigen::Ref<Eigen::MatrixXd> hapcountMatrix
) {
    int K = m_data.nAncestries();
    for (int k = 0; k < K; ++k) {
        getAdmixGenotypes(markerLocalIdx, k, dosageMatrix.col(k), hapcountMatrix.col(k));
    }
}

// ======================================================================
// AbedWriter
// ======================================================================

AbedWriter::AbedWriter(
    const std::string &filename,
    uint8_t nAnc,
    uint32_t nSamples,
    bool noMissing,
    int nthreads
)
    : m_bgzf(nullptr),
      m_filename(filename),
      m_nSamples(nSamples),
      m_nAnc(nAnc),
      m_noMissing(noMissing),
      m_headerWritten(false)
{
    m_bgzf = bgzf_open(filename.c_str(), "w");
    if (!m_bgzf) throw std::runtime_error("Cannot create ABED file: " + filename);

    if (nthreads > 1) bgzf_mt(m_bgzf, nthreads, 0);

    if (bgzf_index_build_init(m_bgzf) < 0) {
        bgzf_close(m_bgzf);
        m_bgzf = nullptr;
        throw std::runtime_error("Cannot init BGZF index for: " + filename);
    }
}

void AbedWriter::ensureHeader() {
    if (m_headerWritten) return;

    AbedHeader hdr{};
    hdr.magic[0] = ABED_MAGIC_0;
    hdr.magic[1] = ABED_MAGIC_1;
    hdr.version = ABED_VERSION;
    hdr.nAnc = m_nAnc | (m_noMissing ? ABED_NANC_FLAG_NO_MISSING : 0);
    hdr.nSamples = m_nSamples;

    if (bgzf_write(m_bgzf, &hdr, sizeof(hdr)) < 0)
        throw std::runtime_error("Failed to write ABED header: " + m_filename);
    m_headerWritten = true;
}

void AbedWriter::writeTrack(
    const uint8_t *data,
    uint64_t nBytes
) {
    ensureHeader();
    if (bgzf_write(m_bgzf, data, nBytes) < 0) throw std::runtime_error("Failed to write ABED track");
}

// ── AVX2 helper for PLINK encode (32 int8 samples → 8 packed bytes) ─────────
__attribute__((target("avx2")))
static uint64_t encodeSamplesAVX2(
    const int8_t *src,
    uint8_t *packed,
    uint64_t n
) {
    const __m256i vlut = _mm256_setr_epi8(0, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 3, 1, 1, 1, 1, 1, 1, 1,
                                          1, 1, 1, 1, 1, 1);
    const __m256i vmask_lo = _mm256_set1_epi8(0x0F);

    uint64_t i = 0;
    uint64_t ob = 0;

    for (; i + 31 < n; i += 32, ob += 8) {
        __m256i v = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(src + i));
        __m256i codes = _mm256_shuffle_epi8(vlut, _mm256_and_si256(v, vmask_lo));

        const __m256i vmul = _mm256_set1_epi16(0x0401);
        __m256i pairs = _mm256_maddubs_epi16(codes, vmul);

        const __m256i vpack16 = _mm256_setr_epi8(0, 2, 4, 6, 8, 10, 12, 14, -1, -1, -1, -1, -1, -1, -1, -1, 0, 2, 4, 6,
                                                 8, 10, 12, 14, -1, -1, -1, -1, -1, -1, -1, -1);
        __m256i narrowed = _mm256_shuffle_epi8(pairs, vpack16);

        __m128i lo128 = _mm256_castsi256_si128(narrowed);
        __m128i hi128 = _mm256_extracti128_si256(narrowed, 1);
        const __m128i vmul2 = _mm_set1_epi16(0x1001);
        __m128i q_lo = _mm_maddubs_epi16(lo128, vmul2);
        __m128i q_hi = _mm_maddubs_epi16(hi128, vmul2);
        const __m128i vpack8 = _mm_setr_epi8(0, 2, 4, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1);
        __m128i r_lo = _mm_shuffle_epi8(q_lo, vpack8);
        __m128i r_hi = _mm_shuffle_epi8(q_hi, vpack8);
        __m128i result = _mm_unpacklo_epi32(r_lo, r_hi);
        _mm_storel_epi64(reinterpret_cast<__m128i *>(packed + ob), result);
    }
    return i;
}

std::vector<uint8_t> AbedWriter::encode(const std::vector<int8_t> &values) {
    // PLINK encoding: 0→0b00, 1→0b10, 2→0b11, -1(missing)→0b01
    uint64_t n = values.size();
    uint64_t nBytes = (n + 3) / 4;
    std::vector<uint8_t> packed(nBytes, 0);

    uint64_t i = 0;

    if (simdLevel() >= SimdLevel::AVX2)
        i = encodeSamplesAVX2(values.data(), packed.data(), n);

    // ── Scalar tail ───────────────────────────────────────────────
    for (; i < n; ++i) {
        uint8_t code;
        switch (values[i]) {
        case 0:
            code = 0b00;
            break;
        case 1:
            code = 0b10;
            break;
        case 2:
            code = 0b11;
            break;
        default:
            code = 0b01;
            break; // missing
        }
        uint64_t byteIdx = i / 4;
        uint32_t shift = (i % 4) * 2;
        packed[byteIdx] |= (code << shift);
    }
    return packed;
}

void AbedWriter::close() {
    if (m_bgzf) {
        ensureHeader();
        if (bgzf_index_dump(m_bgzf, m_filename.c_str(), ".gzi") < 0)
            infoMsg("Warning: failed to write BGZF index: %s.gzi", m_filename.c_str());
        bgzf_close(m_bgzf);
        m_bgzf = nullptr;
    }
}
