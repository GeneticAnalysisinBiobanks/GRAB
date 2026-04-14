// plink.cpp — PlinkData and PlinkCursor implementations (pure C++17 / Eigen)

#include "geno_factory/plink.hpp"
#include "geno_factory/hwe.hpp"

// plink2 SIMD primitives for BED decode + counting
#include "pgenlib_misc.h"

#include "util/text_scanner.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace {

// ──────────────────────────────────────────────────────────────────────
// File-local helpers
// ──────────────────────────────────────────────────────────────────────

std::vector<std::string> splitWhitespace(const std::string &line) {
    text::TokenScanner ts(line);
    std::vector<std::string> tokens;
    while (!ts.atEnd()) {
        auto sv = ts.nextView();
        if (sv.empty()) break;
        tokens.emplace_back(sv);
    }
    return tokens;
}

std::vector<std::string> readSingleColumnFile(const std::string &path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open filter file: " + path);
    std::vector<std::string> values;
    std::string line;
    while (std::getline(in, line)) {
        auto tokens = splitWhitespace(line);
        if (!tokens.empty()) values.push_back(std::move(tokens[0]));
    }
    return values;
}

struct RangeFilter {
    std::string chrom;
    uint32_t start;
    uint32_t end;
};

std::vector<RangeFilter> readRangeFile(const std::string &path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open range filter file: " + path);
    std::vector<RangeFilter> ranges;
    std::string line;
    while (std::getline(in, line)) {
        auto tokens = splitWhitespace(line);
        if (tokens.empty()) continue;
        if (tokens.size() != 3) throw std::runtime_error("Range filter file needs exactly 3 columns: " + path);
        ranges.push_back(
            {tokens[0], static_cast<uint32_t>(std::stoul(tokens[1])), static_cast<uint32_t>(std::stoul(tokens[2]))});
    }
    return ranges;
}

bool markerInRanges(
    const PlinkData::MarkerInfo &m,
    const std::vector<RangeFilter> &ranges
) {
    for (const auto &r : ranges)
        if (m.chrom == r.chrom && m.pos >= r.start && m.pos <= r.end) return true;
    return false;
}

// BED 2-bit genotype → count of bim col5 (A1 = ALT in output)
// BED codes: 0 = hom A1, 1 = missing, 2 = het, 3 = hom A2
static constexpr int GENO_BIM5[4] = {2, -1, 1, 0}; // count A1 (bim col5 = ALT)

// ──────────────────────────────────────────────────────────────────────
// BED-code-to-double lookup table for plink2's GenoarrLookup16x8bx2.
// BED codes: 0=homA1→2.0, 1=miss→NaN, 2=het→1.0, 3=homA2→0.0.
// The table has 16 entries (one per nibble = two consecutive 2-bit genotypes),
// each entry is a pair of doubles (16 bytes).
static const double kBedDoublePairs[32] __attribute__((aligned(16))) =
    PAIR_TABLE16(2.0, std::numeric_limits<double>::quiet_NaN(), 1.0, 0.0);

// ──────────────────────────────────────────────────────────────────────
// Simple masked decode: bitmask scatter, no QC stats, missing → NaN.
// ──────────────────────────────────────────────────────────────────────

static const double kNaN = std::numeric_limits<double>::quiet_NaN();

static void decodeMaskedSimple(
    const uint8_t *__restrict raw,
    const uint64_t *__restrict usedMask,
    uint32_t nFam,
    const int *genoMap,
    double *__restrict out
) {
    uint32_t denseIdx = 0;
    const uint32_t nWords = (nFam + 63) / 64;

    for (uint32_t w = 0; w < nWords; ++w) {
        uint64_t mask = usedMask[w];
        while (mask) {
            const uint32_t bit = static_cast<uint32_t>(__builtin_ctzll(mask));
            const uint32_t famIdx = w * 64 + bit;
            const int code = (raw[famIdx >> 2] >> ((famIdx & 3) << 1)) & 3;
            const int g = genoMap[code];
            out[denseIdx++] = g < 0 ? kNaN : static_cast<double>(g);
            mask &= mask - 1;
        }
    }
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// PlinkData
// ══════════════════════════════════════════════════════════════════════

PlinkData::PlinkData(
    std::string bedFile,
    std::string bimFile,
    std::string famFile,
    const std::vector<uint64_t> &usedMask,
    uint32_t nFam,
    uint32_t nUsed,
    std::string IDsToIncludeFile,
    std::string RangesToIncludeFile,
    std::string IDsToExcludeFile,
    std::string RangesToExcludeFile,
    std::unordered_set<std::string> chrFilter,
    int nMarkersEachChunk
)
    : m_bedFile(std::move(bedFile)),
      m_allUsed(nUsed == nFam),
      m_nSubjInFile(nFam),
      m_nSubjUsed(nUsed),
      m_usedMask(usedMask)
{
    // ---- Parse .bim ----
    {
        std::ifstream in(bimFile);
        if (!in.is_open()) throw std::runtime_error("Cannot open PLINK bim file: " + bimFile);
        std::string line;
        while (std::getline(in, line)) {
            auto tokens = splitWhitespace(line);
            if (tokens.empty()) continue;
            if (tokens.size() != 6) throw std::runtime_error("PLINK .bim file needs 6 columns: " + bimFile);
            m_chr.push_back(tokens[0]);
            m_markerId.push_back(tokens[1]);
            m_pos.push_back(static_cast<uint32_t>(std::stoul(tokens[3])));
            std::string a1 = tokens[4], a2 = tokens[5];
            std::transform(a1.begin(), a1.end(), a1.begin(), ::toupper);
            std::transform(a2.begin(), a2.end(), a2.begin(), ::toupper);
            m_ref.push_back(std::move(a1)); // bim col5 = A1 = REF in plink
            m_alt.push_back(std::move(a2)); // bim col6 = A2
        }
        m_nMarkers = static_cast<uint32_t>(m_chr.size());
    }

    // ---- Count .fam lines for validation ----
    {
        std::ifstream in(famFile);
        if (!in.is_open()) throw std::runtime_error("Cannot open PLINK fam file: " + famFile);
        uint32_t famCount = 0;
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty() && line.find_first_not_of(" \t\r\n") != std::string::npos) ++famCount;
        }
        if (famCount != nFam)
            throw std::runtime_error("PLINK .fam line count (" + std::to_string(famCount) + ") does not match nFam (" +
                                     std::to_string(nFam) + ")");
    }
    m_bytesPerMarker = (m_nSubjInFile + 3) / 4;

    // ---- Validate .bed header ----
    {
        std::ifstream bed(m_bedFile, std::ios::binary);
        if (!bed.is_open()) throw std::runtime_error("Cannot open PLINK bed file: " + m_bedFile);
        char magic[3] = {0};
        bed.read(magic, 3);
        if (!bed.good() || magic[0] != 0x6C || magic[1] != 0x1B || magic[2] != 0x01)
            throw std::runtime_error("Invalid or unsupported PLINK bed file format");
    }

    // ---- Filter markers and build chunks ----
    m_markerInfo = getFilteredMarkers(m_chr, m_pos, m_markerId, m_ref, m_alt, IDsToIncludeFile, RangesToIncludeFile,
                                      IDsToExcludeFile, RangesToExcludeFile, chrFilter);
    if (m_markerInfo.empty()) throw std::runtime_error("No markers remain after PLINK marker filtering.");
    m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);
}

std::vector<PlinkData::MarkerInfo> PlinkData::getFilteredMarkers(
    const std::vector<std::string> &chr,
    const std::vector<uint32_t> &pos,
    const std::vector<std::string> &markerId,
    const std::vector<std::string> &ref,
    const std::vector<std::string> &alt,
    const std::string &IDsToIncludeFile,
    const std::string &RangesToIncludeFile,
    const std::string &IDsToExcludeFile,
    const std::string &RangesToExcludeFile,
    const std::unordered_set<std::string> &chrFilter
) {
    const uint32_t nMarkers = static_cast<uint32_t>(chr.size());
    std::vector<PlinkData::MarkerInfo> all;
    all.reserve(nMarkers);
    for (uint32_t i = 0; i < nMarkers; ++i) {
        all.push_back({chr[i], pos[i], markerId[i], ref[i], alt[i], static_cast<uint64_t>(i)});
    }

    std::unordered_set<std::string> includeIds, excludeIds;
    std::vector<RangeFilter> includeRanges, excludeRanges;

    if (!IDsToIncludeFile.empty()) {
        auto ids = readSingleColumnFile(IDsToIncludeFile);
        includeIds.insert(ids.begin(), ids.end());
    }
    if (!RangesToIncludeFile.empty()) includeRanges = readRangeFile(RangesToIncludeFile);
    if (!IDsToExcludeFile.empty()) {
        auto ids = readSingleColumnFile(IDsToExcludeFile);
        excludeIds.insert(ids.begin(), ids.end());
    }
    if (!RangesToExcludeFile.empty()) excludeRanges = readRangeFile(RangesToExcludeFile);

    const bool anyChr     = !chrFilter.empty();
    const bool anyInclude = !includeIds.empty() || !includeRanges.empty();
    const bool anyExclude = !excludeIds.empty() || !excludeRanges.empty();
    if (!anyChr && !anyInclude && !anyExclude) return all;

    std::vector<PlinkData::MarkerInfo> filtered;
    filtered.reserve(all.size());
    for (const auto &m : all) {
        if (anyChr && chrFilter.count(m.chrom) == 0) continue;
        bool inc = !anyInclude;
        if (anyInclude) inc = (includeIds.count(m.id) > 0) || markerInRanges(m, includeRanges);
        if (!inc) continue;
        bool exc = false;
        if (anyExclude) exc = (excludeIds.count(m.id) > 0) || markerInRanges(m, excludeRanges);
        if (!exc) filtered.push_back(m);
    }
    return filtered;
}

std::vector<std::vector<uint64_t> > PlinkData::buildChunks(
    const std::vector<MarkerInfo> &markers,
    int chunkSize
) {
    std::vector<std::vector<uint64_t> > chunks;
    size_t start = 0;
    while (start < markers.size()) {
        const std::string &chrom = markers[start].chrom;
        size_t chromEnd = start;
        while (chromEnd < markers.size() && markers[chromEnd].chrom == chrom)
            ++chromEnd;

        for (size_t cs = start; cs < chromEnd; cs += static_cast<size_t>(chunkSize)) {
            size_t ce = std::min(cs + static_cast<size_t>(chunkSize), chromEnd);
            std::vector<uint64_t> chunk;
            chunk.reserve(ce - cs);
            for (size_t i = cs; i < ce; ++i)
                chunk.push_back(markers[i].genoIndex);
            chunks.push_back(std::move(chunk));
        }
        start = chromEnd;
    }
    return chunks;
}

// ══════════════════════════════════════════════════════════════════════
// PlinkCursor
// ══════════════════════════════════════════════════════════════════════

PlinkCursor::PlinkCursor(
    const std::string &bedFile,
    uint32_t nBimLines,
    uint32_t nFamLines,
    const std::vector<uint64_t> &usedMask,
    uint32_t nUsed,
    bool allUsed
)
    : m_bedFile(bedFile),
      m_nSubjInFile(nFamLines),
      m_bytesPerMarker((nFamLines + 3) / 4),
      m_nMarkers(nBimLines),
      m_nUsed(nUsed),
      m_usedMask(usedMask),
      m_allUsed(allUsed),
      m_rawBytes(m_bytesPerMarker),
      m_alignedBed((m_bytesPerMarker + sizeof(uintptr_t) - 1) / sizeof(uintptr_t))
{
    m_bedStream.open(bedFile, std::ios::binary);
    if (!m_bedStream.is_open()) throw std::runtime_error("Cannot open PLINK bed file: " + bedFile);
    char magic[3] = {0};
    m_bedStream.read(magic, 3);
    if (!m_bedStream.good() || magic[0] != 0x6C || magic[1] != 0x1B || magic[2] != 0x01)
        throw std::runtime_error("Invalid or unsupported PLINK bed file format: " + bedFile);
}

PlinkCursor::PlinkCursor(const PlinkCursor &other)
    : m_bedFile(other.m_bedFile),
      m_nSubjInFile(other.m_nSubjInFile),
      m_bytesPerMarker(other.m_bytesPerMarker),
      m_nMarkers(other.m_nMarkers),
      m_nUsed(other.m_nUsed),
      m_usedMask(other.m_usedMask),
      m_allUsed(other.m_allUsed),
      m_rawBytes(other.m_rawBytes.size()),
      m_alignedBed(other.m_alignedBed.size())
{
    m_bedStream.open(m_bedFile, std::ios::binary);
    if (!m_bedStream.is_open()) throw std::runtime_error("Cannot open PLINK bed file: " + m_bedFile);
    char magic[3] = {0};
    m_bedStream.read(magic, 3);
    if (!m_bedStream.good() || magic[0] != 0x6C || magic[1] != 0x1B || magic[2] != 0x01)
        throw std::runtime_error("Invalid or unsupported PLINK bed file format: " + m_bedFile);
}

void PlinkCursor::beginSequentialBlock(uint64_t firstMarker) {
    if (firstMarker >= m_nMarkers) throw std::runtime_error("PLINK marker index out of range in beginSequentialBlock.");
    m_bedStream.clear();
    m_bedStream.seekg(3 + static_cast<std::streamoff>(m_bytesPerMarker) * firstMarker);
    if (!m_bedStream.good()) throw std::runtime_error("Failed to seek PLINK bed file.");
    m_hasSeqCursor = true;
    m_nextMarker = firstMarker;
    m_hasBlock = false;
}

void PlinkCursor::loadBlock(uint64_t startMarker) {
    if (startMarker >= m_nMarkers) throw std::runtime_error("PLINK marker index out of range when loading block.");

    const uint64_t nRemain = m_nMarkers - startMarker;
    const uint64_t nRead = std::min(BLOCK_CAPACITY, nRemain);
    const uint64_t nBytes = nRead * m_bytesPerMarker;

    if (m_blockBytes.size() < static_cast<size_t>(nBytes)) m_blockBytes.resize(static_cast<size_t>(nBytes));

    if (!m_hasSeqCursor || m_nextMarker != startMarker) {
        m_bedStream.clear();
        m_bedStream.seekg(3 + static_cast<std::streamoff>(m_bytesPerMarker) * startMarker);
        if (!m_bedStream.good()) throw std::runtime_error("Failed to seek PLINK bed file for block read.");
    }

    m_bedStream.read(reinterpret_cast<char *>(m_blockBytes.data()), static_cast<std::streamsize>(nBytes));
    if (!m_bedStream.good()) throw std::runtime_error("Failed to read PLINK marker block.");

    m_blockStart = startMarker;
    m_blockEnd = startMarker + nRead;
    m_hasBlock = true;
    m_hasSeqCursor = true;
    m_nextMarker = m_blockEnd;
}

const uint8_t *PlinkCursor::readMarkerPtr(uint64_t gIndex) {
    const uint64_t bpm = m_bytesPerMarker;

    if (m_hasBlock && gIndex >= m_blockStart && gIndex < m_blockEnd) {
        return m_blockBytes.data() + (gIndex - m_blockStart) * bpm;
    }

    if (m_hasSeqCursor && gIndex == m_nextMarker) {
        loadBlock(gIndex);
        return m_blockBytes.data();
    }

    const auto filePos = static_cast<std::streamoff>(3 + bpm * gIndex);
    m_bedStream.clear();
    m_bedStream.seekg(filePos);
    if (!m_bedStream.good()) throw std::runtime_error("seek failed");

    m_bedStream.read(reinterpret_cast<char *>(m_rawBytes.data()), bpm);
    if (!m_bedStream.good()) throw std::runtime_error("read failed");

    m_hasBlock = false;
    m_hasSeqCursor = true;
    m_nextMarker = gIndex + 1;

    return m_rawBytes.data();
}

// ── All-used fast path ───────────────────────────────────────────────
// All subjects in file order → plink2 SIMD decode + count.

void PlinkCursor::getGenotypesAllUsed(
    const uint8_t *raw,
    uint32_t n,
    Eigen::Ref<Eigen::VectorXd> out,
    double &altFreq,
    double &altCounts,
    double &missingRate,
    double &hweP,
    double &maf,
    double &mac,
    std::vector<uint32_t> &indexForMissing
) {
    indexForMissing.clear();

    // Copy BED bytes to word-aligned buffer; trailing bytes stay zero.
    m_alignedBed.back() = 0;
    std::memcpy(m_alignedBed.data(), raw, m_bytesPerMarker);
    const uintptr_t *bedWords = m_alignedBed.data();

    // Vectorized counting via plink2 carry-save-adder.
    // Counts are in BED code order: [0]=homA1, [1]=miss, [2]=het, [3]=homA2.
    std::array<uint32_t, 4> bedCounts;
    plink2::GenoarrCountFreqsUnsafe(bedWords, n, bedCounts);

    // Vectorized decode via plink2 SIMD lookup table.
    plink2::GenoarrLookup16x8bx2(bedWords, kBedDoublePairs, n, out.data());

    // Collect missing indices (BED code 1 → NaN) — skip scan when none missing.
    if (bedCounts[1] > 0) {
        for (uint32_t i = 0; i < n; ++i)
            if (std::isnan(out[i])) indexForMissing.push_back(i);
    }

    // BED code mapping: [0]=homA1, [1]=miss, [2]=het, [3]=homA2
    GenoStats gs = statsFromCounts(bedCounts[0], bedCounts[2], bedCounts[3], bedCounts[1], n);
    altFreq = gs.altFreq;
    altCounts = gs.altCounts;
    missingRate = gs.missingRate;
    hweP = gs.hweP;
    maf = gs.maf;
    mac = gs.mac;
}

// ── Bitmask path ─────────────────────────────────────────────────────
// Subset selection via bitmask → plink2 vectorized counting +
// scalar scatter for decode.

void PlinkCursor::getGenotypesMasked(
    const uint8_t *raw,
    Eigen::Ref<Eigen::VectorXd> out,
    double &altFreq,
    double &altCounts,
    double &missingRate,
    double &hweP,
    double &maf,
    double &mac,
    std::vector<uint32_t> &indexForMissing
) {
    indexForMissing.clear();

    // Copy BED bytes to word-aligned buffer for plink2.
    m_alignedBed.back() = 0;
    std::memcpy(m_alignedBed.data(), raw, m_bytesPerMarker);

    // Vectorized subset counting via plink2 popcount-based loop.
    // Counts are in BED code order: [0]=homA1, [1]=miss, [2]=het, [3]=homA2.
    std::array<uint32_t, 4> bedCounts;
    plink2::GenoarrCountSubsetFreqs2(m_alignedBed.data(), reinterpret_cast<const uintptr_t *>(m_usedMask.data()),
                                     m_nSubjInFile, m_nUsed, bedCounts);

    // Scatter decode: write only used subjects into dense output.
    const int *genoMap = GENO_BIM5;
    uint32_t denseIdx = 0;
    const uint32_t nWords = (m_nSubjInFile + 63) / 64;
    for (uint32_t w = 0; w < nWords; ++w) {
        uint64_t mask = m_usedMask[w];
        while (mask) {
            const uint32_t bit = static_cast<uint32_t>(__builtin_ctzll(mask));
            const uint32_t famIdx = w * 64 + bit;
            const int code = (raw[famIdx >> 2] >> ((famIdx & 3) << 1)) & 3;
            const int g = genoMap[code];
            if (g < 0) {
                out[denseIdx] = kNaN;
                indexForMissing.push_back(denseIdx);
            } else {
                out[denseIdx] = static_cast<double>(g);
            }
            ++denseIdx;
            mask &= mask - 1;
        }
    }

    // BED code mapping: [0]=homA1, [1]=miss, [2]=het, [3]=homA2
    GenoStats gs = statsFromCounts(bedCounts[0], bedCounts[2], bedCounts[3], bedCounts[1], m_nUsed);
    altFreq = gs.altFreq;
    altCounts = gs.altCounts;
    missingRate = gs.missingRate;
    hweP = gs.hweP;
    maf = gs.maf;
    mac = gs.mac;
}

// ── Public entry point ───────────────────────────────────────────────

void PlinkCursor::getGenotypes(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out,
    double &altFreq,
    double &altCounts,
    double &missingRate,
    double &hweP,
    double &maf,
    double &mac,
    std::vector<uint32_t> &indexForMissing
) {
    const uint8_t *raw = readMarkerPtr(gIndex);

    if (m_allUsed) {
        getGenotypesAllUsed(raw, m_nSubjInFile, out, altFreq, altCounts, missingRate, hweP, maf, mac, indexForMissing);
    } else {
        getGenotypesMasked(raw, out, altFreq, altCounts, missingRate, hweP, maf, mac, indexForMissing);
    }
}

// ── Simple entry point: genotype vector only, missing → NaN ──────────

void PlinkCursor::getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out
) {
    const uint8_t *raw = readMarkerPtr(gIndex);

    if (m_allUsed) {
        // Vectorized decode via plink2 SIMD lookup table.
        m_alignedBed.back() = 0;
        std::memcpy(m_alignedBed.data(), raw, m_bytesPerMarker);
        plink2::GenoarrLookup16x8bx2(m_alignedBed.data(), kBedDoublePairs, m_nSubjInFile, out.data());
    } else {
        decodeMaskedSimple(raw, m_usedMask.data(), m_nSubjInFile, GENO_BIM5, out.data());
    }
}

// ══════════════════════════════════════════════════════════════════════
// GenoMeta factory
// ══════════════════════════════════════════════════════════════════════

std::unique_ptr<GenoCursor> PlinkData::makeCursor() const {
    return std::make_unique<PlinkCursor>(bedFile(), nMarkers(), nSubjInFile(), usedMask(), nSubjUsed(), allUsed());
}
