// admix.cpp — Admixed ancestry binary genotype reader/writer (.abed format)

#include "io/admix.hpp"
#include "io/subject_data.hpp"   // parseFamIIDs, parseBimLines
#include "util/logging.hpp"

extern "C" {
#include <htslib/bgzf.h>
}

// x86 SIMD header
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
#include <immintrin.h>
#endif

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>


// ── BIM parsing (shared with plink.hpp) ─────────────────────────────

static std::vector<GenoMeta::MarkerInfo> parseBimFile(const std::string& bimFile) {
    std::ifstream ifs(bimFile);
    if (!ifs) throw std::runtime_error("Cannot open BIM file: " + bimFile);

    std::vector<GenoMeta::MarkerInfo> markers;
    std::string line;
    uint64_t idx = 0;
    while (std::getline(ifs, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string chrom, id, cm, ref, alt;
        uint32_t pos = 0;
        ss >> chrom >> id >> cm >> pos >> ref >> alt;
        if (ss.fail())
            throw std::runtime_error("Malformed BIM line: " + line);
        markers.push_back({std::move(chrom), pos, std::move(id),
                           std::move(ref), std::move(alt), idx++});
    }
    return markers;
}


// ======================================================================
// AdmixData
// ======================================================================

AdmixData::AdmixData(
    const std::string& prefix,
    const std::vector<uint64_t>& usedMask,
    uint32_t nFam,
    uint32_t nUsed,
    const std::string& extractFile,
    const std::string& excludeFile,
    int nMarkersEachChunk
) {
    m_abedFile = prefix + ".abed";
    std::string bimFile = prefix + ".bim";
    std::string famFile = prefix + ".fam";

    // Read and validate header
    BGZF* bgz = bgzf_open(m_abedFile.c_str(), "r");
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
        throw std::runtime_error("Unsupported ABED version in " + m_abedFile);
    if (hdr.mode != ABED_MODE_SNP_MAJOR)
        throw std::runtime_error("Only SNP-major mode supported in " + m_abedFile);

    m_nAnc        = hdr.nAncestries;
    m_nSubjInFile = hdr.nSamples;
    m_noMissing   = (hdr.reserved & ABED_FLAG_NO_MISSING) != 0;
    uint32_t fileMark = hdr.nMarkers;

    if (m_nSubjInFile != nFam)
        throw std::runtime_error("ABED header says " + std::to_string(m_nSubjInFile) +
                                 " samples but .fam has " + std::to_string(nFam));

    m_nSubjUsed      = nUsed;
    m_usedMask       = usedMask;
    m_bytesPerTrack  = (static_cast<uint64_t>(m_nSubjInFile) + 3) / 4;
    m_bytesPerMarker = 2ULL * m_nAnc * m_bytesPerTrack;

    // Check all used
    uint64_t maskBits = 0;
    for (auto w : m_usedMask) maskBits += __builtin_popcountll(w);
    m_allUsed = (maskBits == m_nSubjInFile);

    // Parse .bim
    auto allMarkers = parseBimFile(bimFile);
    if (allMarkers.size() != fileMark)
        throw std::runtime_error("ABED header says " + std::to_string(fileMark) +
                                 " markers but .bim has " + std::to_string(allMarkers.size()));

    // Apply extract/exclude filters
    m_markerInfo = filterMarkers(allMarkers, extractFile, excludeFile);
    m_nMarkers   = static_cast<uint32_t>(m_markerInfo.size());

    // Build chunks
    m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);

    infoMsg("Loaded ABED: %u ancestries x %u markers x %u samples (%u used)",
            m_nAnc, m_nMarkers, m_nSubjInFile, m_nSubjUsed);
}

std::vector<GenoMeta::MarkerInfo> AdmixData::filterMarkers(
    const std::vector<MarkerInfo>& all,
    const std::string& extractFile,
    const std::string& excludeFile
) {
    std::unordered_set<std::string> includeSet, excludeSet;

    auto readIDFile = [](const std::string& path, std::unordered_set<std::string>& out) {
        std::ifstream ifs(path);
        if (!ifs) throw std::runtime_error("Cannot open SNP list file: " + path);
        std::string id;
        while (ifs >> id) {
            if (!id.empty() && id[0] != '#')
                out.insert(std::move(id));
        }
    };

    if (!extractFile.empty()) readIDFile(extractFile, includeSet);
    if (!excludeFile.empty()) readIDFile(excludeFile, excludeSet);

    if (includeSet.empty() && excludeSet.empty()) return all;

    std::vector<MarkerInfo> filtered;
    filtered.reserve(all.size());
    for (const auto& m : all) {
        if (!includeSet.empty() && includeSet.find(m.id) == includeSet.end()) continue;
        if (!excludeSet.empty() && excludeSet.find(m.id) != excludeSet.end()) continue;
        filtered.push_back(m);
    }
    return filtered;
}

std::vector<std::vector<uint64_t>> AdmixData::buildChunks(
    const std::vector<MarkerInfo>& markers, int chunkSize
) {
    std::vector<std::vector<uint64_t>> chunks;
    if (markers.empty()) return chunks;

    std::vector<uint64_t> cur;
    cur.reserve(chunkSize);
    std::string prevChr = markers[0].chrom;

    for (size_t i = 0; i < markers.size(); ++i) {
        if (markers[i].chrom != prevChr) {
            if (!cur.empty()) { chunks.push_back(std::move(cur)); cur.clear(); cur.reserve(chunkSize); }
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

AdmixCursor::AdmixCursor(const AdmixData& data)
    : m_data(data)
    , m_bgzf(nullptr)
    , m_nUsed(data.nSubjUsed())
    , m_bytesPerTrack(data.bytesPerTrack())
    , m_bytesPerMarker(data.bytesPerMarker())
    , m_nTracks(2 * data.nAncestries()
) {
    m_bgzf = bgzf_open(data.abedFile().c_str(), "r");
    if (!m_bgzf)
        throw std::runtime_error("AdmixCursor: cannot open " + data.abedFile());

    // Load .gzi index for random access via bgzf_useek (required for BGZF files;
    // silently ignored for uncompressed legacy .abed files)
    int idxRc = bgzf_index_load(m_bgzf, data.abedFile().c_str(), ".gzi");
    (void)idxRc;

    m_rawBytes.resize(m_bytesPerTrack);
}

AdmixCursor::~AdmixCursor() {
    if (m_bgzf) bgzf_close(m_bgzf);
}

void AdmixCursor::beginSequentialBlock(uint64_t firstMarker) {
    m_hasSeqCursor = true;
    m_nextMarker   = firstMarker;
    m_hasBlock     = false;
}

void AdmixCursor::loadBlock(uint64_t startMarker) {
    uint64_t nTotal = m_data.markerInfo().size();
    uint64_t end = std::min(startMarker + BLOCK_CAPACITY, nTotal);
    uint64_t count = end - startMarker;
    if (count == 0) { m_hasBlock = false; return; }

    // Resolve genoIndex of first marker for file offset
    uint64_t firstGenoIdx = m_data.markerInfo()[startMarker].genoIndex;
    uint64_t offset = ABED_HEADER_SIZE + firstGenoIdx * m_bytesPerMarker;

    m_blockBytes.resize(count * m_bytesPerMarker);
    if (bgzf_useek(m_bgzf, static_cast<off_t>(offset), SEEK_SET) < 0)
        throw std::runtime_error("AdmixCursor: bgzf_useek failed in loadBlock");
    if (bgzf_read(m_bgzf, m_blockBytes.data(),
                  count * m_bytesPerMarker)
        != static_cast<ssize_t>(count * m_bytesPerMarker))
        throw std::runtime_error("AdmixCursor: bgzf_read failed in loadBlock");

    m_blockStart = startMarker;
    m_blockEnd   = end;
    m_hasBlock   = true;
}

const uint8_t* AdmixCursor::readTrackPtr(uint64_t markerLocalIdx, int trackIdx) {
    // Check block cache
    if (m_hasBlock && markerLocalIdx >= m_blockStart && markerLocalIdx < m_blockEnd) {
        // Check if genoIndices are contiguous in the block
        uint64_t genoIdx = m_data.markerInfo()[markerLocalIdx].genoIndex;
        uint64_t baseGenoIdx = m_data.markerInfo()[m_blockStart].genoIndex;
        uint64_t relMarker = genoIdx - baseGenoIdx;
        return m_blockBytes.data() + relMarker * m_bytesPerMarker
                                   + trackIdx * m_bytesPerTrack;
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
    uint64_t offset = ABED_HEADER_SIZE + genoIdx * m_bytesPerMarker
                    + trackIdx * m_bytesPerTrack;
    if (bgzf_useek(m_bgzf, static_cast<off_t>(offset), SEEK_SET) < 0)
        throw std::runtime_error("AdmixCursor: bgzf_useek failed");
    if (bgzf_read(m_bgzf, m_rawBytes.data(), m_bytesPerTrack)
        != static_cast<ssize_t>(m_bytesPerTrack))
        throw std::runtime_error("AdmixCursor: bgzf_read failed");
    return m_rawBytes.data();
}

void AdmixCursor::decodeTrack(const uint8_t* raw,
                              Eigen::Ref<Eigen::VectorXd> out
) {
    // Lookup table: 00→0, 10→1, 11→2, 01→missing(NaN)
    static const double kDecode[4] = {0.0, std::numeric_limits<double>::quiet_NaN(), 1.0, 2.0};
    // No-missing variant: 01→0 (treat as ref) — avoids NaN in the output
    static const double kDecodeNoMiss[4] = {0.0, 0.0, 1.0, 2.0};

    const bool noMissing = m_data.hasNoMissing();

    if (m_data.allUsed()) {
        uint32_t n = m_data.nSubjInFile();
        double* __restrict p = out.data();

#if defined(__AVX2__)
        if (noMissing) {
            // ── AVX2 fast path: 8 bytes → 32 samples per iteration ──
            // Map 2-bit codes {0,1,2,3} → doubles {0.0, 0.0, 1.0, 2.0}
            // Strategy: extract 2-bit codes → 32-bit ints, apply LUT, cvt to double
            //
            // For each byte, 4 samples are extracted, then widened to int32 and
            // converted to double via _mm256_cvtepi32_pd (4 doubles per 128-bit int).
            // Process 8 bytes (32 samples) at a time via 8 groups of 4 cvt ops.

            // LUT: geno code → integer value (0→0, 1→0, 2→1, 3→2)
            alignas(16) static const int32_t lut[4] = {0, 0, 1, 2};

            uint32_t fullBytes = n / 4;
            uint32_t di = 0;

            uint32_t b = 0;
            for (; b + 7 < fullBytes; b += 8) {
                for (int bi = 0; bi < 8; ++bi) {
                    uint8_t byte = raw[b + bi];
                    // Extract 4 codes, look up in LUT, convert to double
                    alignas(16) int32_t codes[4] = {
                        (byte     ) & 3,
                        (byte >> 2) & 3,
                        (byte >> 4) & 3,
                        (byte >> 6) & 3
                    };
                    __m128i vcodes = _mm_load_si128(reinterpret_cast<const __m128i*>(codes));
                    // Gather from LUT: result[i] = lut[codes[i]]
                    __m128i vvals = _mm_i32gather_epi32(lut, vcodes, 4);
                    // Convert int32 → double (4 values)
                    __m256d vd = _mm256_cvtepi32_pd(vvals);
                    _mm256_storeu_pd(p + di, vd);
                    di += 4;
                }
            }
            // Scalar tail
            for (; b < fullBytes; ++b) {
                uint8_t byte = raw[b];
                p[di++] = kDecodeNoMiss[(byte     ) & 3];
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
#endif // __AVX2__

        // ── Scalar path ──
        const double* lut = noMissing ? kDecodeNoMiss : kDecode;
        uint32_t fullBytes = n / 4;
        uint32_t rem = n % 4;
        uint32_t di = 0;
        for (uint32_t b = 0; b < fullBytes; ++b) {
            uint8_t byte = raw[b];
            p[di++] = lut[(byte     ) & 3];
            p[di++] = lut[(byte >> 2) & 3];
            p[di++] = lut[(byte >> 4) & 3];
            p[di++] = lut[(byte >> 6) & 3];
        }
        if (rem > 0) {
            uint8_t byte = raw[fullBytes];
            for (uint32_t r = 0; r < rem; ++r)
                p[di++] = lut[(byte >> (2 * r)) & 3];
        }
    } else {
        // Bitmask path: decode only used subjects
        const double* lut = noMissing ? kDecodeNoMiss : kDecode;
        const auto& mask = m_data.usedMask();
        uint32_t nWords = static_cast<uint32_t>(mask.size());
        uint32_t di = 0;
        for (uint32_t w = 0; w < nWords; ++w) {
            uint64_t bits = mask[w];
            while (bits) {
                int bit = __builtin_ctzll(bits);
                uint32_t subjIdx = w * 64 + bit;
                uint32_t byteIdx = subjIdx / 4;
                uint32_t shift   = (subjIdx % 4) * 2;
                uint8_t code = (raw[byteIdx] >> shift) & 3;
                out[di++] = lut[code];
                bits &= bits - 1;  // clear lowest set bit
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
    int dosageTrack  = 2 * ancIdx;
    int hapTrack     = 2 * ancIdx + 1;

    const uint8_t* dosageRaw = readTrackPtr(markerLocalIdx, dosageTrack);
    decodeTrack(dosageRaw, dosageOut);

    const uint8_t* hapRaw = readTrackPtr(markerLocalIdx, hapTrack);
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
        getAdmixGenotypes(markerLocalIdx, k,
                          dosageMatrix.col(k), hapcountMatrix.col(k));
    }
}


// ======================================================================
// AbedWriter
// ======================================================================

AbedWriter::AbedWriter(const std::string& filename, uint16_t nAnc,
                       uint32_t nSamples, uint32_t nMarkers,
                       uint16_t flags)
    : m_bgzf(nullptr)
    , m_filename(filename)
    , m_nSamples(nSamples)
    , m_flags(flags)
    , m_nAnc(nAnc)
    , m_nMarkers(nMarkers)
    , m_headerWritten(false)
{
    m_bgzf = bgzf_open(filename.c_str(), "w");
    if (!m_bgzf)
        throw std::runtime_error("Cannot create ABED file: " + filename);

    if (bgzf_index_build_init(m_bgzf) < 0) {
        bgzf_close(m_bgzf);
        m_bgzf = nullptr;
        throw std::runtime_error("Cannot init BGZF index for: " + filename);
    }
}

void AbedWriter::ensureHeader() {
    if (m_headerWritten) return;

    AbedHeader hdr{};
    hdr.magic[0]    = ABED_MAGIC_0;
    hdr.magic[1]    = ABED_MAGIC_1;
    hdr.version     = ABED_VERSION;
    hdr.mode        = ABED_MODE_SNP_MAJOR;
    hdr.nAncestries = m_nAnc;
    hdr.nSamples    = m_nSamples;
    hdr.nMarkers    = m_nMarkers;
    hdr.reserved    = m_flags;

    if (bgzf_write(m_bgzf, &hdr, sizeof(hdr)) < 0)
        throw std::runtime_error("Failed to write ABED header: " + m_filename);
    m_headerWritten = true;
}

void AbedWriter::writeTrack(const uint8_t* data, uint64_t nBytes) {
    ensureHeader();
    if (bgzf_write(m_bgzf, data, nBytes) < 0)
        throw std::runtime_error("Failed to write ABED track");
}

std::vector<uint8_t> AbedWriter::encode(const std::vector<int8_t>& values) {
    // PLINK encoding: 0→0b00, 1→0b10, 2→0b11, -1(missing)→0b01
    uint64_t n = values.size();
    uint64_t nBytes = (n + 3) / 4;
    std::vector<uint8_t> packed(nBytes, 0);

    for (uint64_t i = 0; i < n; ++i) {
        uint8_t code;
        switch (values[i]) {
            case 0:  code = 0b00; break;
            case 1:  code = 0b10; break;
            case 2:  code = 0b11; break;
            default: code = 0b01; break;  // missing
        }
        uint64_t byteIdx = i / 4;
        uint32_t shift   = (i % 4) * 2;
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
