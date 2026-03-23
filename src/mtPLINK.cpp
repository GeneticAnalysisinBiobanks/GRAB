// mtPLINK.cpp -- PlinkData and PlinkCursor implementations

#include "mtPLINK.h"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <sstream>
#include <mutex>
#include <immintrin.h>
#include <cstdint>
#include <limits>


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

bool markerInRanges(const PlinkData::MarkerInfo& m, const std::vector<RangeFilter>& ranges) {
  for (const auto& r : ranges)
    if (m.chrom == r.chrom && m.pos >= r.start && m.pos <= r.end) return true;
  return false;
}

// BED 2-bit genotype → alt allele count (indexed by 2-bit BED code 0..3)
// BED codes: 0 = hom A1, 1 = missing, 2 = het, 3 = hom A2
static const int GENO_ALT_FIRST[4] = { 2, -1,  1,  0};  // alt=A1
static const int GENO_REF_FIRST[4] = { 0, -1,  1,  2};  // alt=A2

// Repack subset samples from full-cohort BED into compact 2-bit buffer.
static void repackSubset(const uint8_t* __restrict raw,
                         const uint32_t* __restrict posMap,
                         uint8_t* __restrict compact,
                         uint32_t nSubset) {
  const uint32_t compactBytes = (nSubset + 3) / 4;
  std::fill(compact, compact + compactBytes, uint8_t(0));
  for (uint32_t i = 0; i < nSubset; ++i) {
    const uint32_t src = posMap[i];
    const int code = (raw[src >> 2] >> ((src & 3) << 1)) & 3;
    compact[i >> 2] |= static_cast<uint8_t>(code << ((i & 3) << 1));
  }
}


// ==================== AVX2 nibble-LUT ====================

alignas(32) static uint8_t LUT_HOMREF[32];
alignas(32) static uint8_t LUT_HET[32];
alignas(32) static uint8_t LUT_HOMALT[32];
alignas(32) static uint8_t LUT_MISS[32];

static void init_lut() {
  static std::once_flag flag;
  std::call_once(flag, []() {
    for (int n = 0; n < 16; ++n) {
      uint8_t homRef = 0, het = 0, homAlt = 0, miss = 0;
      for (int s = 0; s < 2; ++s) {
        int code = (n >> (2*s)) & 3;
        if (code == 0) homRef++;
        else if (code == 1) miss++;
        else if (code == 2) het++;
        else homAlt++;
      }
      LUT_HOMREF[n] = homRef;
      LUT_HET[n]    = het;
      LUT_HOMALT[n] = homAlt;
      LUT_MISS[n]   = miss;
    }
    for (int i = 16; i < 32; ++i) {
      LUT_HOMREF[i] = LUT_HOMREF[i - 16];
      LUT_HET[i]    = LUT_HET[i - 16];
      LUT_HOMALT[i] = LUT_HOMALT[i - 16];
      LUT_MISS[i]   = LUT_MISS[i - 16];
    }
  });
}

// ==================== COUNT ====================

static inline void count_geno_classes(
  const uint8_t* __restrict geno_bytes,
  uint32_t n_samples,
  uint32_t out[4]
) {
  init_lut();

  const uint32_t n_bytes = (n_samples + 3) / 4;
  uint64_t nHomRef = 0, nHet = 0, nHomAlt = 0, nMissing = 0;
  const __m256i mask_lo = _mm256_set1_epi8(0x0F);
  const __m256i zero = _mm256_setzero_si256();

  const __m256i lut_homref = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(LUT_HOMREF));
  const __m256i lut_het = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(LUT_HET));
  const __m256i lut_homalt = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(LUT_HOMALT));
  const __m256i lut_miss = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(LUT_MISS));

  uint32_t i = 0;
  uint32_t limit = n_bytes & ~63u;

  __m256i acc_h = _mm256_setzero_si256();
  __m256i acc_e = _mm256_setzero_si256();
  __m256i acc_a = _mm256_setzero_si256();
  __m256i acc_m = _mm256_setzero_si256();

  int iter = 0;

  for (; i < limit; i += 64) {
    __m256i v0 = _mm256_loadu_si256((const __m256i*)(geno_bytes + i));
    __m256i v1 = _mm256_loadu_si256((const __m256i*)(geno_bytes + i + 32));

#define PROCESS_VEC(v) \
    do { \
      __m256i lo = _mm256_and_si256(v, mask_lo); \
      __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), mask_lo); \
      acc_h = _mm256_add_epi8(acc_h, \
          _mm256_add_epi8( \
            _mm256_shuffle_epi8(lut_homref, lo), \
            _mm256_shuffle_epi8(lut_homref, hi))); \
      acc_e = _mm256_add_epi8(acc_e, \
          _mm256_add_epi8( \
            _mm256_shuffle_epi8(lut_het, lo), \
            _mm256_shuffle_epi8(lut_het, hi))); \
      acc_a = _mm256_add_epi8(acc_a, \
          _mm256_add_epi8( \
            _mm256_shuffle_epi8(lut_homalt, lo), \
            _mm256_shuffle_epi8(lut_homalt, hi))); \
      acc_m = _mm256_add_epi8(acc_m, \
          _mm256_add_epi8( \
            _mm256_shuffle_epi8(lut_miss, lo), \
            _mm256_shuffle_epi8(lut_miss, hi))); \
    } while (0)

    PROCESS_VEC(v0);
    PROCESS_VEC(v1);

#undef PROCESS_VEC

    ++iter;

    // Flush accumulators before uint8 overflow.
    // Max increment per lane per iteration = 8 (two PROCESS_VECs × max 4 each).
    // Safe limit: floor(255 / 8) = 31 iterations.
    if (iter >= 31) {
      __m256i s;

      s = _mm256_sad_epu8(acc_h, zero);
      uint64_t tmp[4];
      _mm256_storeu_si256((__m256i*)tmp, s);
      nHomRef += tmp[0] + tmp[1] + tmp[2] + tmp[3];

      s = _mm256_sad_epu8(acc_e, zero);
      _mm256_storeu_si256((__m256i*)tmp, s);
      nHet += tmp[0] + tmp[1] + tmp[2] + tmp[3];

      s = _mm256_sad_epu8(acc_a, zero);
      _mm256_storeu_si256((__m256i*)tmp, s);
      nHomAlt += tmp[0] + tmp[1] + tmp[2] + tmp[3];

      s = _mm256_sad_epu8(acc_m, zero);
      _mm256_storeu_si256((__m256i*)tmp, s);
      nMissing += tmp[0] + tmp[1] + tmp[2] + tmp[3];

      acc_h = acc_e = acc_a = acc_m = zero;
      iter = 0;
    }
  }

  // final flush
  if (iter > 0) {
    __m256i s;
    uint64_t tmp[4];

    s = _mm256_sad_epu8(acc_h, zero);
    _mm256_storeu_si256((__m256i*)tmp, s);
    nHomRef += tmp[0] + tmp[1] + tmp[2] + tmp[3];

    s = _mm256_sad_epu8(acc_e, zero);
    _mm256_storeu_si256((__m256i*)tmp, s);
    nHet += tmp[0] + tmp[1] + tmp[2] + tmp[3];

    s = _mm256_sad_epu8(acc_a, zero);
    _mm256_storeu_si256((__m256i*)tmp, s);
    nHomAlt += tmp[0] + tmp[1] + tmp[2] + tmp[3];

    s = _mm256_sad_epu8(acc_m, zero);
    _mm256_storeu_si256((__m256i*)tmp, s);
    nMissing += tmp[0] + tmp[1] + tmp[2] + tmp[3];
  }

  // tail
  for (; i < n_bytes; ++i) {
    uint8_t byte = geno_bytes[i];
    for (int s = 0; s < 4; ++s) {
      int code = (byte >> (2*s)) & 3;
      if (code == 0) nHomRef++;
      else if (code == 1) nMissing++;
      else if (code == 2) nHet++;
      else nHomAlt++;
    }
  }

  // padding fix
  uint32_t valid_last = n_samples % 4;
  if (valid_last != 0) {
    uint8_t last = geno_bytes[n_bytes - 1];
    for (uint32_t s = valid_last; s < 4; ++s) {
      int code = (last >> (2*s)) & 3;
      if (code == 0) --nHomRef;
      else if (code == 1) --nMissing;
      else if (code == 2) --nHet;
      else --nHomAlt;
    }
  }

  out[0] = (uint32_t)nHomRef;
  out[1] = (uint32_t)nHet;
  out[2] = (uint32_t)nHomAlt;
  out[3] = (uint32_t)nMissing;
}

// ==================== HWE tests ====================

// Chi-squared HWE approximation — O(1), chi²(1 df) upper-tail p = erfc(sqrt(chi2/2)).
// Returns NaN when any expected genotype count < 5 (chi² not applicable).
static double HweChiSq(uint32_t nHet, uint32_t nHom1, uint32_t nHom2) {
  const double n  = static_cast<double>(nHet + nHom1 + nHom2);
  if (n == 0.0) return std::numeric_limits<double>::quiet_NaN();
  const double f  = (2.0 * nHom1 + nHet) / (2.0 * n);  // allele A1 freq
  const double g  = 1.0 - f;
  const double Eh = 2.0 * f * g * n;
  const double E1 = f * f * n;
  const double E2 = g * g * n;
  // Chi-sq approximation requires all expected counts >= 5
  if (E1 < 5.0 || Eh < 5.0 || E2 < 5.0)
    return std::numeric_limits<double>::quiet_NaN();
  const double d1 = static_cast<double>(nHom1) - E1;
  const double dh = static_cast<double>(nHet)  - Eh;
  const double d2 = static_cast<double>(nHom2) - E2;
  const double chi2 = d1*d1/E1 + dh*dh/Eh + d2*d2/E2;
  return std::erfc(std::sqrt(chi2 * 0.5));
}

// Exact HWE test — matches plink2 --hardy default (SNPHWE2).
// Wigginton JE, Cutler DJ, Abecasis GR (2005). A Note on Exact Tests of
// Hardy-Weinberg Equilibrium. Am J Hum Genet 76:887-893.
//
// Optimized for repeated calls (e.g. 30M markers): no heap allocation, O(1)
// auxiliary memory. Walks the "obs side" first to derive thresh inline, then
// walks the opposite side and accumulates p on the fly.
static double HweExact(uint32_t obs_hets, uint32_t obs_hom1, uint32_t obs_hom2) {
  const int64_t obs_homc = std::max(obs_hom1, obs_hom2);
  const int64_t obs_homr = std::min(obs_hom1, obs_hom2);
  const int64_t rare = 2 * obs_homr + (int64_t)obs_hets;
  const int64_t n    = (int64_t)obs_hets + obs_homc + obs_homr;
  const int64_t obs  = (int64_t)obs_hets;

  if (n == 0) return 1.0;

  // Mode of het-count distribution under HWE
  int64_t mid = (rare * (2 * n - rare)) / (2 * n);
  if ((rare & 1) ^ (mid & 1)) ++mid;

  // The floor+parity estimate can land one grid step from the true mode.
  // Adjust mid to sit at the actual peak of the distribution.
  {
    int64_t hr = (rare - mid) / 2;
    int64_t hc = n - mid - hr;
    if (mid + 2 <= rare && hr > 0 &&
        4.0 * hr * hc > (mid + 2.0) * (mid + 1.0)) {
      mid += 2;
    } else if (mid >= 2) {
      if (static_cast<double>(mid) * (mid - 1) >
          4.0 * (hr + 1.0) * (hc + 1.0)) {
        mid -= 2;
      }
    }
  }

  const int64_t mid_homr = (rare - mid) / 2;
  const int64_t mid_homc = n - mid - mid_homr;

  // All probabilities are relative to probs[mid] = 1 (unnormalized).
  // thresh = probs[obs].  p = sum of probs[h] <= thresh.  sum = normalizer.
  double sum = 1.0, p = 0.0, thresh;

  if (obs <= mid) {
    // --- Downward walk (obs side): derive thresh, then count obs's outer tail ---
    {
      double prob = 1.0;
      int64_t cr = mid_homr, cc = mid_homc;
      // Phase 1: mid → obs+2, probs > thresh; contribute to sum only
      for (int64_t h = mid; h > obs; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob;
        ++cr; ++cc;
      }
      thresh = prob;   // prob[obs]
      p      = thresh; // obs itself is always counted
      // Phase 2: obs → 2, probs <= thresh; contribute to both sum and p
      for (int64_t h = obs; h >= 2; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob;
        p   += prob;
        ++cr; ++cc;
      }
    }
    // --- Upward walk (other side): monotone decrease; count tail <= thresh ---
    {
      double prob = 1.0;
      int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h <= rare - 2; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob;
        if (prob <= thresh) p += prob;
        --cr; --cc;
      }
    }
  } else {
    // --- Upward walk (obs side): derive thresh, then count obs's outer tail ---
    {
      double prob = 1.0;
      int64_t cr = mid_homr, cc = mid_homc;
      // Phase 1: mid → obs-2, probs > thresh; contribute to sum only
      for (int64_t h = mid; h < obs; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob;
        --cr; --cc;
      }
      thresh = prob;   // prob[obs]
      p      = thresh; // obs itself is always counted
      // Phase 2: obs → rare, probs <= thresh; contribute to both sum and p
      for (int64_t h = obs; h <= rare - 2; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob;
        p   += prob;
        --cr; --cc;
      }
    }
    // --- Downward walk (other side): monotone decrease; count tail <= thresh ---
    {
      double prob = 1.0;
      int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h >= 2; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob;
        if (prob <= thresh) p += prob;
        ++cr; ++cc;
      }
    }
  }

  return std::min(p / sum, 1.0);
}


// ==================== compute_marker_stats_avx2 ====================

struct GenoStats {
  double   altFreq;
  uint32_t altCounts;
  double   missingRate;
  double   hweP;
  double   maf;
  uint32_t mac;
};

// // Compute genotype counts and derived statistics from packed BED bytes.
// // Uses AVX2 intrinsics for the main loop with scalar tail fallback.
GenoStats compute_marker_stats_avx2(const uint8_t* geno_bytes, uint32_t n_samples, bool altFirst, bool exactHwe) {
  uint32_t counts[4];  // homA1(code0), het(code2), homA2(code3), missing(code1)
  count_geno_classes(geno_bytes, n_samples, counts);

  const uint32_t nHomA1   = counts[0];
  const uint32_t nHet     = counts[1];
  const uint32_t nHomA2   = counts[2];
  const uint32_t nMissing = counts[3];
  const uint32_t nonMissing = n_samples - nMissing;

  GenoStats gs;
  // altCounts is derived from AVX2-computed genotype class counts
  gs.altCounts = altFirst ? (2 * nHomA1 + nHet) : (nHet + 2 * nHomA2);
  gs.missingRate = static_cast<double>(nMissing) / n_samples;

  if (nonMissing == 0) {
    gs.altFreq = std::numeric_limits<double>::quiet_NaN();
    gs.hweP    = std::numeric_limits<double>::quiet_NaN();
    gs.maf     = std::numeric_limits<double>::quiet_NaN();
    gs.mac     = 0;
    return gs;
  }

  gs.altFreq = static_cast<double>(gs.altCounts) / (2.0 * nonMissing);
  gs.maf = std::min(gs.altFreq, 1.0 - gs.altFreq);
  gs.mac = std::min(gs.altCounts, 2 * nonMissing - gs.altCounts);
  if (exactHwe) {
    gs.hweP = HweExact(nHet, nHomA1, nHomA2);
  } else {
    if (nHet < 5.0 || nHomA1 < 5.0 || nHomA2 < 5.0) {
      gs.hweP = HweExact(nHet, nHomA1, nHomA2);
    } else {
      gs.hweP = HweChiSq(nHet, nHomA1, nHomA2);
    }
  }
  return gs;
}

} // anonymous namespace


// ==================== PlinkData ====================

PlinkData::PlinkData(
    std::string bedFile,
    std::string bimFile,
    std::string famFile,
    std::vector<std::string> subjData,
    std::string AlleleOrder,
    std::string IDsToIncludeFile,
    std::string RangesToIncludeFile,
    std::string IDsToExcludeFile,
    std::string RangesToExcludeFile,
    int nMarkersEachChunk)
  : m_bedFile(std::move(bedFile)),
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
  std::vector<std::string> famIIDs;
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
      famIIDs.push_back(tokens[1]);
    }
  }
  m_nSubjInFile    = static_cast<uint32_t>(famIIDs.size());
  m_bytesPerMarker = (m_nSubjInFile + 3) / 4;

  // ---- Build sample position map + validate ----
  m_nSubjUsed = static_cast<uint32_t>(subjData.size());
  std::unordered_map<std::string, uint32_t> famPosLookup;
  famPosLookup.reserve(famIIDs.size());
  for (uint32_t i = 0; i < famIIDs.size(); ++i)
    famPosLookup.emplace(famIIDs[i], i);

  m_samplePosMap.resize(m_nSubjUsed);
  for (uint32_t i = 0; i < m_nSubjUsed; ++i) {
    auto it = famPosLookup.find(subjData[i]);
    if (it == famPosLookup.end())
      throw std::runtime_error("Subject not found in PLINK file: " + subjData[i]);
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

  // ---- Filter markers and build chunks ----
  m_markerInfo = getFilteredMarkers(m_chr, m_pos, m_markerId, m_ref, m_alt,
                                    IDsToIncludeFile, RangesToIncludeFile,
                                    IDsToExcludeFile, RangesToExcludeFile);
  if (m_markerInfo.empty())
    throw std::runtime_error("No markers remain after PLINK marker filtering.");
  m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);
}

std::vector<PlinkData::MarkerInfo> PlinkData::getFilteredMarkers(
    const std::vector<std::string>& chr,
    const std::vector<uint32_t>& pos,
    const std::vector<std::string>& markerId,
    const std::vector<std::string>& ref,
    const std::vector<std::string>& alt,
    const std::string& IDsToIncludeFile,
    const std::string& RangesToIncludeFile,
    const std::string& IDsToExcludeFile,
    const std::string& RangesToExcludeFile) {

  // Build full marker list
  const uint32_t nMarkers = static_cast<uint32_t>(chr.size());
  std::vector<PlinkData::MarkerInfo> all;
  all.reserve(nMarkers);
  for (uint32_t i = 0; i < nMarkers; ++i) {
    all.push_back({chr[i], pos[i], markerId[i],
                   ref[i], alt[i], static_cast<uint64_t>(i)});
  }

  // Load filter sets
  std::unordered_set<std::string> includeIds, excludeIds;
  std::vector<RangeFilter> includeRanges, excludeRanges;

  if (!IDsToIncludeFile.empty()) {
    auto ids = readSingleColumnFile(IDsToIncludeFile);
    includeIds.insert(ids.begin(), ids.end());
  }
  if (!RangesToIncludeFile.empty())
    includeRanges = readRangeFile(RangesToIncludeFile);
  if (!IDsToExcludeFile.empty()) {
    auto ids = readSingleColumnFile(IDsToExcludeFile);
    excludeIds.insert(ids.begin(), ids.end());
  }
  if (!RangesToExcludeFile.empty())
    excludeRanges = readRangeFile(RangesToExcludeFile);

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

PlinkCursor::PlinkCursor(const std::string& bedFile,
                         uint32_t nBimLines, uint32_t nFamLines,
                         std::vector<uint32_t> posMap, bool altFirst)
  : m_bedFile(bedFile),
    m_nSubjInFile(nFamLines),
    m_bytesPerMarker((nFamLines + 3) / 4),
    m_nMarkers(nBimLines),
    m_posMap(std::move(posMap)),
    m_altFirst(altFirst),
    m_rawBytes(m_bytesPerMarker),
    m_compactGeno((m_posMap.size() + 3) / 4)
{
  m_bedStream.open(bedFile, std::ios::binary);
  if (!m_bedStream.is_open())
    throw std::runtime_error("Cannot open PLINK bed file: " + bedFile);
  // Validate 3-byte header
  char magic[3] = {0};
  m_bedStream.read(magic, 3);
  if (!m_bedStream.good() || magic[0] != 0x6C || magic[1] != 0x1B || magic[2] != 0x01)
    throw std::runtime_error("Invalid or unsupported PLINK bed file format: " + bedFile);
}

PlinkCursor::PlinkCursor(const PlinkCursor& other)
  : m_bedFile(other.m_bedFile),
    m_nSubjInFile(other.m_nSubjInFile),
    m_bytesPerMarker(other.m_bytesPerMarker),
    m_nMarkers(other.m_nMarkers),
    m_posMap(other.m_posMap),
    m_altFirst(other.m_altFirst),
    m_rawBytes(other.m_rawBytes.size()),
    m_compactGeno(other.m_compactGeno.size())
{
  m_bedStream.open(m_bedFile, std::ios::binary);
  if (!m_bedStream.is_open())
    throw std::runtime_error("Cannot open PLINK bed file: " + m_bedFile);
  char magic[3] = {0};
  m_bedStream.read(magic, 3);
  if (!m_bedStream.good() || magic[0] != 0x6C || magic[1] != 0x1B || magic[2] != 0x01)
    throw std::runtime_error("Invalid or unsupported PLINK bed file format: " + m_bedFile);
}

void PlinkCursor::beginSequentialBlock(uint64_t firstMarker) {
  if (firstMarker >= m_nMarkers)
    throw std::runtime_error("PLINK marker index out of range in beginSequentialBlock.");
  m_bedStream.clear();
  m_bedStream.seekg(3 + static_cast<uint64_t>(m_bytesPerMarker) * firstMarker);
  if (!m_bedStream.good())
    throw std::runtime_error("Failed to seek PLINK bed file.");
  m_hasSeqCursor = true;
  m_nextMarker = firstMarker;
  m_hasBlock = false;
}

void PlinkCursor::loadBlock(uint64_t startMarker) {
  if (startMarker >= m_nMarkers)
    throw std::runtime_error("PLINK marker index out of range when loading block.");

  const uint64_t nRemain = m_nMarkers - startMarker;
  const uint64_t nRead = std::min(BLOCK_CAPACITY, nRemain);
  const uint64_t nBytes = nRead * m_bytesPerMarker;

  if (m_blockBytes.size() != nBytes)
    m_blockBytes.resize(static_cast<size_t>(nBytes));

  if (!m_hasSeqCursor || m_nextMarker != startMarker) {
    m_bedStream.clear();
    m_bedStream.seekg(3 + static_cast<uint64_t>(m_bytesPerMarker) * startMarker);
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

const uint8_t* PlinkCursor::readMarkerPtr(uint64_t gIndex) {
  const uint64_t bpm = m_bytesPerMarker;

  if (m_hasBlock && gIndex >= m_blockStart && gIndex < m_blockEnd) {
    return m_blockBytes.data() + (gIndex - m_blockStart) * bpm;
  }

  if (m_hasSeqCursor && gIndex == m_nextMarker) {
    loadBlock(gIndex);
    return m_blockBytes.data();
  }

  const uint64_t filePos = 3 + bpm * gIndex;
  m_bedStream.clear();
  m_bedStream.seekg(filePos);
  if (!m_bedStream.good())
    throw std::runtime_error("seek failed");

  m_bedStream.read(reinterpret_cast<char*>(m_rawBytes.data()), bpm);
  if (!m_bedStream.good())
    throw std::runtime_error("read failed");

  m_hasBlock = false;
  m_hasSeqCursor = true;
  m_nextMarker = gIndex + 1;

  return m_rawBytes.data();
}

void PlinkCursor::getGenotypes(
    uint64_t gIndex,
    arma::vec& out,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    double& hweP,
    double& maf,
    double& mac,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe
) {
  const uint8_t* raw = readMarkerPtr(gIndex);
  const uint32_t n = static_cast<uint32_t>(m_posMap.size());

  // Step 1: Repack subset into compact 2-bit buffer
  repackSubset(raw, m_posMap.data(), m_compactGeno.data(), n);

  // Step 2: AVX2 stats on subset
  GenoStats gs = compute_marker_stats_avx2(m_compactGeno.data(), n, m_altFirst, exactHwe);
  altFreq     = gs.altFreq;
  altCounts   = gs.altCounts;
  missingRate = gs.missingRate;
  hweP        = gs.hweP;
  maf         = gs.maf;
  mac         = gs.mac;

  // Step 3: Decode compact 2-bit to float (reuse caller's buffer)
  const int* genoMap = m_altFirst ? GENO_ALT_FIRST : GENO_REF_FIRST;
  out.set_size(n);
  indexForMissing.clear();

  for (uint32_t i = 0; i < n; ++i) {
    const int code = (m_compactGeno[i >> 2] >> ((i & 3) << 1)) & 3;
    const int g = genoMap[code];
    if (g < 0) {
      out[i] = -1;
      indexForMissing.push_back(i);
    } else {
      out[i] = g;
    }
  }
}
