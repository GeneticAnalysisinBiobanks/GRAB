// PLINK4.cpp -- PlinkData and PlinkCursor implementations

#include "PLINK4.h"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <sstream>

namespace {

// ---- File-local helpers ----

std::vector<std::string> splitWhitespace(const std::string& line) {
  std::istringstream iss(line);
  std::vector<std::string> tokens;
  std::string tok;
  while (iss >> tok) tokens.push_back(tok);
  return tokens;
}

std::vector<std::string> readSingleColumnFile(const std::string& path) {
  std::ifstream in(path);
  if (!in.is_open())
    throw std::runtime_error("Cannot open filter file: " + path);
  std::vector<std::string> values;
  std::string line;
  while (std::getline(in, line)) {
    auto tokens = splitWhitespace(line);
    if (!tokens.empty()) values.push_back(tokens[0]);
  }
  return values;
}

struct RangeFilter {
  std::string chrom;
  uint32_t start;
  uint32_t end;
};

std::vector<RangeFilter> readRangeFile(const std::string& path) {
  std::ifstream in(path);
  if (!in.is_open())
    throw std::runtime_error("Cannot open range filter file: " + path);
  std::vector<RangeFilter> ranges;
  std::string line;
  while (std::getline(in, line)) {
    auto tokens = splitWhitespace(line);
    if (tokens.empty()) continue;
    if (tokens.size() != 3)
      throw std::runtime_error("Range filter file needs exactly 3 columns: " + path);
    ranges.push_back({
      tokens[0],
      static_cast<uint32_t>(std::stoul(tokens[1])),
      static_cast<uint32_t>(std::stoul(tokens[2]))
    });
  }
  return ranges;
}

bool markerInRanges(const PLINK4::MarkerInfo& m, const std::vector<RangeFilter>& ranges) {
  for (const auto& r : ranges)
    if (m.chrom == r.chrom && m.pos >= r.start && m.pos <= r.end) return true;
  return false;
}

// BED 2-bit genotype → alt allele count (indexed by 2-bit BED code 0..3)
// BED codes: 0 = hom A1, 1 = missing, 2 = het, 3 = hom A2
static const int GENO_ALT_FIRST[4] = { 2, -1,  1,  0};  // alt=A1
static const int GENO_REF_FIRST[4] = { 0, -1,  1,  2};  // alt=A2

} // anonymous namespace


namespace PLINK4 {

// ==================== PlinkData ====================

PlinkData::PlinkData(
    const std::string& bedFile,
    const std::string& bimFile,
    const std::string& famFile,
    const std::vector<std::string>& requestedSamples,
    const std::string& AlleleOrder)
  : m_bedFile(bedFile),
    m_altFirst(AlleleOrder == "alt-first")
{
  // ---- Parse .bim ----
  {
    std::ifstream in(bimFile);
    if (!in.is_open())
      throw std::runtime_error("Cannot open PLINK bim file: " + bimFile);
    std::string line;
    while (std::getline(in, line)) {
      auto tokens = splitWhitespace(line);
      if (tokens.empty()) continue;
      if (tokens.size() != 6)
        throw std::runtime_error("PLINK .bim file needs 6 columns: " + bimFile);
      m_chr.push_back(tokens[0]);
      m_markerId.push_back(tokens[1]);
      m_pos.push_back(static_cast<uint32_t>(std::stoul(tokens[3])));
      // A1 = tokens[4], A2 = tokens[5]
      std::string a1 = tokens[4], a2 = tokens[5];
      std::transform(a1.begin(), a1.end(), a1.begin(), ::toupper);
      std::transform(a2.begin(), a2.end(), a2.begin(), ::toupper);
      if (m_altFirst) {          // alt-first: alt=A1, ref=A2
        m_alt.push_back(a1);
        m_ref.push_back(a2);
      } else {                   // ref-first: ref=A1, alt=A2
        m_ref.push_back(a1);
        m_alt.push_back(a2);
      }
    }
    m_nMarkers = static_cast<uint32_t>(m_chr.size());
  }

  // ---- Parse .fam ----
  std::vector<std::string> famSamples;
  {
    std::ifstream in(famFile);
    if (!in.is_open())
      throw std::runtime_error("Cannot open PLINK fam file: " + famFile);
    std::string line;
    while (std::getline(in, line)) {
      auto tokens = splitWhitespace(line);
      if (tokens.empty()) continue;
      if (tokens.size() != 6)
        throw std::runtime_error("PLINK .fam file needs 6 columns: " + famFile);
      famSamples.push_back(tokens[1]);
    }
  }
  m_nSamplesInFile = static_cast<uint32_t>(famSamples.size());
  m_bytesPerMarker = (m_nSamplesInFile + 3) / 4;

  // ---- Build sample position map + validate ----
  m_nSamplesUsed = static_cast<uint32_t>(requestedSamples.size());
  std::unordered_map<std::string, uint32_t> famPosLookup;
  famPosLookup.reserve(famSamples.size());
  for (uint32_t i = 0; i < famSamples.size(); ++i)
    famPosLookup.emplace(famSamples[i], i);

  m_samplePosMap.resize(m_nSamplesUsed);
  for (uint32_t i = 0; i < m_nSamplesUsed; ++i) {
    auto it = famPosLookup.find(requestedSamples[i]);
    if (it == famPosLookup.end())
      throw std::runtime_error("Subject not found in PLINK file: " + requestedSamples[i]);
    m_samplePosMap[i] = it->second;
  }

  // ---- Validate .bed header ----
  {
    std::ifstream bed(m_bedFile, std::ios::binary);
    if (!bed.is_open())
      throw std::runtime_error("Cannot open PLINK bed file: " + m_bedFile);
      
    char magic[3] = {0};
    bed.read(magic, 3);
    if (!bed.good() || magic[0] != 0x6C || magic[1] != 0x1B || magic[2] != 0x01)
        throw std::runtime_error("Invalid or unsupported PLINK bed file format");
  }
}

std::vector<MarkerInfo> PlinkData::getFilteredMarkers(
  const MarkerFilterConfig& filter
) const {

  // Build full marker list
  std::vector<MarkerInfo> all;
  all.reserve(m_nMarkers);
  for (uint32_t i = 0; i < m_nMarkers; ++i) {
    all.push_back({m_chr[i], m_pos[i], m_markerId[i],
                   m_ref[i], m_alt[i], static_cast<uint64_t>(i)});
  }

  // Load filter sets
  std::unordered_set<std::string> includeIds, excludeIds;
  std::vector<RangeFilter> includeRanges, excludeRanges;

  if (!filter.IDsToIncludeFile.empty()) {
    auto ids = readSingleColumnFile(filter.IDsToIncludeFile);
    includeIds.insert(ids.begin(), ids.end());
  }
  if (!filter.RangesToIncludeFile.empty())
    includeRanges = readRangeFile(filter.RangesToIncludeFile);
  if (!filter.IDsToExcludeFile.empty()) {
    auto ids = readSingleColumnFile(filter.IDsToExcludeFile);
    excludeIds.insert(ids.begin(), ids.end());
  }
  if (!filter.RangesToExcludeFile.empty())
    excludeRanges = readRangeFile(filter.RangesToExcludeFile);

  const bool anyInclude = !includeIds.empty() || !includeRanges.empty();
  const bool anyExclude = !excludeIds.empty() || !excludeRanges.empty();
  if (!anyInclude && !anyExclude) return all;

  std::vector<MarkerInfo> filtered;
  filtered.reserve(all.size());
  for (const auto& m : all) {
    bool inc = !anyInclude;
    if (anyInclude)
      inc = (includeIds.count(m.id) > 0) || markerInRanges(m, includeRanges);
    if (!inc) continue;
    bool exc = false;
    if (anyExclude)
      exc = (excludeIds.count(m.id) > 0) || markerInRanges(m, excludeRanges);
    if (!exc) filtered.push_back(m);
  }
  return filtered;
}

std::vector<std::vector<uint64_t>> PlinkData::buildChunks(
    const std::vector<MarkerInfo>& markers, int chunkSize) {
  std::vector<std::vector<uint64_t>> chunks;
  size_t start = 0;
  while (start < markers.size()) {
    const std::string& chrom = markers[start].chrom;
    size_t chromEnd = start;
    while (chromEnd < markers.size() && markers[chromEnd].chrom == chrom) ++chromEnd;

    for (size_t cs = start; cs < chromEnd;
         cs += static_cast<size_t>(chunkSize)) {
      size_t ce = std::min(cs + static_cast<size_t>(chunkSize), chromEnd);
      std::vector<uint64_t> chunk;
      chunk.reserve(ce - cs);
      for (size_t i = cs; i < ce; ++i) chunk.push_back(markers[i].genoIndex);
      chunks.push_back(std::move(chunk));
    }
    start = chromEnd;
  }
  return chunks;
}


// ==================== PlinkCursor ====================

PlinkCursor::PlinkCursor(const PlinkData& data)
  : m_data(data),
    m_rawBytes(data.bytesPerMarker())
{
  m_bedStream.open(data.bedFile(), std::ios::binary);
  if (!m_bedStream.is_open())
    throw std::runtime_error("Cannot open PLINK bed file: " + data.bedFile());
  // Skip 3-byte header
  m_bedStream.seekg(3);
  if (!m_bedStream.good())
    throw std::runtime_error("Failed to read PLINK bed header: " + data.bedFile());
}

void PlinkCursor::beginSequentialBlock(uint64_t firstMarker) {
  if (firstMarker >= m_data.nMarkers())
    throw std::runtime_error("PLINK marker index out of range in beginSequentialBlock.");
  const uint64_t pos = 3 + m_data.bytesPerMarker() * firstMarker;
  m_bedStream.clear();
  m_bedStream.seekg(pos);
  if (!m_bedStream.good())
    throw std::runtime_error("Failed to seek PLINK bed file.");
  m_hasSeqCursor = true;
  m_nextMarker = firstMarker;
  m_hasBlock = false;
}

void PlinkCursor::loadBlock(uint64_t startMarker) {
  if (startMarker >= m_data.nMarkers())
    throw std::runtime_error("PLINK marker index out of range when loading block.");

  const uint64_t nRemain = static_cast<uint64_t>(m_data.nMarkers()) - startMarker;
  const uint64_t nRead = std::min(BLOCK_CAPACITY, nRemain);
  const uint64_t nBytes = nRead * m_data.bytesPerMarker();

  if (m_blockBytes.size() != nBytes)
    m_blockBytes.resize(static_cast<size_t>(nBytes));

  if (!m_hasSeqCursor || m_nextMarker != startMarker) {
    const uint64_t pos = 3 + m_data.bytesPerMarker() * startMarker;
    m_bedStream.clear();
    m_bedStream.seekg(pos);
    if (!m_bedStream.good())
      throw std::runtime_error("Failed to seek PLINK bed file for block read.");
  }

  m_bedStream.read(reinterpret_cast<char*>(m_blockBytes.data()),
                   static_cast<std::streamsize>(nBytes));
  if (!m_bedStream.good())
    throw std::runtime_error("Failed to read PLINK marker block.");

  m_blockStart = startMarker;
  m_blockEnd = startMarker + nRead;
  m_hasBlock = true;
  m_hasSeqCursor = true;
  m_nextMarker = m_blockEnd;
}

void PlinkCursor::readMarkerBytes(uint64_t gIndex) {
  const uint64_t bpm = m_data.bytesPerMarker();

  // Try buffered block first
  if (m_hasBlock && gIndex >= m_blockStart && gIndex < m_blockEnd) {
    const uint64_t off = (gIndex - m_blockStart) * bpm;
    std::copy(m_blockBytes.begin() + static_cast<std::ptrdiff_t>(off),
              m_blockBytes.begin() + static_cast<std::ptrdiff_t>(off + bpm),
              m_rawBytes.begin());
    return;
  }

  // Try sequential block load
  if (m_hasSeqCursor && gIndex == m_nextMarker) {
    loadBlock(gIndex);
    const uint64_t off = (gIndex - m_blockStart) * bpm;
    std::copy(m_blockBytes.begin() + static_cast<std::ptrdiff_t>(off),
              m_blockBytes.begin() + static_cast<std::ptrdiff_t>(off + bpm),
              m_rawBytes.begin());
    return;
  }

  // Fallback: single-marker read
  const uint64_t filePos = 3 + bpm * gIndex;
  m_bedStream.clear();
  m_bedStream.seekg(filePos);
  if (!m_bedStream.good())
    throw std::runtime_error("Failed to seek PLINK bed file.");
  m_bedStream.read(reinterpret_cast<char*>(m_rawBytes.data()), bpm);
  if (!m_bedStream.good())
    throw std::runtime_error("Failed to read PLINK marker bytes.");
  m_hasBlock = false;
  m_hasSeqCursor = true;
  m_nextMarker = gIndex + 1;
}

arma::vec PlinkCursor::getGenotypes(
    uint64_t gIndex,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    std::vector<uint32_t>& indexForMissing
) {
  readMarkerBytes(gIndex);

  const uint32_t n = m_data.nSamplesUsed();
  const auto& posMap = m_data.samplePosMap();
  const int* genoMap = m_data.isAltFirst() ? GENO_ALT_FIRST : GENO_REF_FIRST;

  arma::vec geno(n);
  indexForMissing.clear();
  int sumAlt = 0;
  int numMissing = 0;

  for (uint32_t i = 0; i < n; ++i) {
    const uint32_t fam_pos = posMap[i];
    const int raw = decodeGenotype(m_rawBytes[fam_pos / 4], fam_pos % 4);
    const int g = genoMap[raw];

    if (g < 0) {   // missing
      geno[i] = -1;
      ++numMissing;
      indexForMissing.push_back(i);
    } else {
      geno[i] = g;
      sumAlt += g;
    }
  }

  const int count = static_cast<int>(n) - numMissing;
  missingRate = static_cast<double>(numMissing) / static_cast<double>(n);
  altCounts = static_cast<double>(sumAlt);
  altFreq = (count > 0) ? altCounts / (2.0 * count) : 0.0;

  return geno;
}

} // namespace PLINK4
