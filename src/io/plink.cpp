// plink.cpp — PlinkData and PlinkCursor implementations (pure C++17 / Eigen)

#include "io/plink.hpp"

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

// ──────────────────────────────────────────────────────────────────────
// File-local helpers
// ──────────────────────────────────────────────────────────────────────

std::vector<std::string> splitWhitespace(const std::string& line) {
  std::istringstream iss(line);
  std::vector<std::string> tokens;
  std::string tok;
  while (iss >> tok) tokens.push_back(std::move(tok));
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
    if (!tokens.empty()) values.push_back(std::move(tokens[0]));
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

bool markerInRanges(const PlinkData::MarkerInfo& m,
                    const std::vector<RangeFilter>& ranges) {
  for (const auto& r : ranges)
    if (m.chrom == r.chrom && m.pos >= r.start && m.pos <= r.end) return true;
  return false;
}

// BED 2-bit genotype → alt allele count (indexed by 2-bit BED code 0..3)
// BED codes: 0 = hom A1, 1 = missing, 2 = het, 3 = hom A2
static constexpr int GENO_ALT_FIRST[4] = { 2, -1,  1,  0};  // alt=A1
static constexpr int GENO_REF_FIRST[4] = { 0, -1,  1,  2};  // alt=A2

// Repack subset samples from full-cohort BED into compact 2-bit buffer.
static void repackSubset(const uint8_t* __restrict raw,
                         const uint32_t* __restrict posMap,
                         uint8_t* __restrict compact,
                         uint32_t nSubset) {
  const uint32_t compactBytes = (nSubset + 3) / 4;
  std::memset(compact, 0, compactBytes);
  for (uint32_t i = 0; i < nSubset; ++i) {
    const uint32_t src = posMap[i];
    const int code = (raw[src >> 2] >> ((src & 3) << 1)) & 3;
    compact[i >> 2] |= static_cast<uint8_t>(code << ((i & 3) << 1));
  }
}


// ──────────────────────────────────────────────────────────────────────
// AVX2 nibble-LUT genotype counting
// ──────────────────────────────────────────────────────────────────────

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
        if      (code == 0) homRef++;
        else if (code == 1) miss++;
        else if (code == 2) het++;
        else                homAlt++;
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

// Flush an AVX2 uint8 accumulator via SAD → 64-bit horizontal sum.
static inline uint64_t hsum_sad(const __m256i& acc) {
  const __m256i s = _mm256_sad_epu8(acc, _mm256_setzero_si256());
  alignas(32) uint64_t tmp[4];
  _mm256_store_si256(reinterpret_cast<__m256i*>(tmp), s);
  return tmp[0] + tmp[1] + tmp[2] + tmp[3];
}

static void count_geno_classes(
  const uint8_t* __restrict geno_bytes,
  uint32_t n_samples,
  uint32_t out[4])
{
  init_lut();

  const uint32_t n_bytes = (n_samples + 3) / 4;
  uint64_t nHomRef = 0, nHet = 0, nHomAlt = 0, nMissing = 0;
  const __m256i mask_lo = _mm256_set1_epi8(0x0F);

  const __m256i lut_homref = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(LUT_HOMREF));
  const __m256i lut_het    = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(LUT_HET));
  const __m256i lut_homalt = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(LUT_HOMALT));
  const __m256i lut_miss   = _mm256_load_si256(
      reinterpret_cast<const __m256i*>(LUT_MISS));

  __m256i acc_h = _mm256_setzero_si256();
  __m256i acc_e = _mm256_setzero_si256();
  __m256i acc_a = _mm256_setzero_si256();
  __m256i acc_m = _mm256_setzero_si256();

  uint32_t i = 0;
  const uint32_t limit = n_bytes & ~63u;
  int iter = 0;

  for (; i < limit; i += 64) {
    __m256i v0 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(geno_bytes + i));
    __m256i v1 = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(geno_bytes + i + 32));

#define PROCESS_VEC(v) do { \
    __m256i lo = _mm256_and_si256(v, mask_lo); \
    __m256i hi = _mm256_and_si256(_mm256_srli_epi16(v, 4), mask_lo); \
    acc_h = _mm256_add_epi8(acc_h, _mm256_add_epi8( \
        _mm256_shuffle_epi8(lut_homref, lo), _mm256_shuffle_epi8(lut_homref, hi))); \
    acc_e = _mm256_add_epi8(acc_e, _mm256_add_epi8( \
        _mm256_shuffle_epi8(lut_het, lo), _mm256_shuffle_epi8(lut_het, hi))); \
    acc_a = _mm256_add_epi8(acc_a, _mm256_add_epi8( \
        _mm256_shuffle_epi8(lut_homalt, lo), _mm256_shuffle_epi8(lut_homalt, hi))); \
    acc_m = _mm256_add_epi8(acc_m, _mm256_add_epi8( \
        _mm256_shuffle_epi8(lut_miss, lo), _mm256_shuffle_epi8(lut_miss, hi))); \
  } while (0)

    PROCESS_VEC(v0);
    PROCESS_VEC(v1);
#undef PROCESS_VEC

    if (++iter >= 31) {
      nHomRef  += hsum_sad(acc_h);
      nHet     += hsum_sad(acc_e);
      nHomAlt  += hsum_sad(acc_a);
      nMissing += hsum_sad(acc_m);
      acc_h = acc_e = acc_a = acc_m = _mm256_setzero_si256();
      iter = 0;
    }
  }

  // flush remaining accumulator
  if (iter > 0) {
    nHomRef  += hsum_sad(acc_h);
    nHet     += hsum_sad(acc_e);
    nHomAlt  += hsum_sad(acc_a);
    nMissing += hsum_sad(acc_m);
  }

  // scalar tail
  for (; i < n_bytes; ++i) {
    uint8_t byte = geno_bytes[i];
    for (int s = 0; s < 4; ++s) {
      int code = (byte >> (2*s)) & 3;
      if      (code == 0) nHomRef++;
      else if (code == 1) nMissing++;
      else if (code == 2) nHet++;
      else                nHomAlt++;
    }
  }

  // undo padding genotypes in the last byte
  uint32_t valid_last = n_samples & 3;
  if (valid_last != 0) {
    uint8_t last = geno_bytes[n_bytes - 1];
    for (uint32_t s = valid_last; s < 4; ++s) {
      int code = (last >> (2*s)) & 3;
      if      (code == 0) --nHomRef;
      else if (code == 1) --nMissing;
      else if (code == 2) --nHet;
      else                --nHomAlt;
    }
  }

  out[0] = static_cast<uint32_t>(nHomRef);
  out[1] = static_cast<uint32_t>(nHet);
  out[2] = static_cast<uint32_t>(nHomAlt);
  out[3] = static_cast<uint32_t>(nMissing);
}


// ──────────────────────────────────────────────────────────────────────
// HWE tests
// ──────────────────────────────────────────────────────────────────────

// Chi-squared approximation: O(1), valid when expected counts >= 5.
static double HweChiSq(uint32_t nHet, uint32_t nHom1, uint32_t nHom2) {
  const double n  = static_cast<double>(nHet + nHom1 + nHom2);
  if (n == 0.0) return std::numeric_limits<double>::quiet_NaN();
  const double f  = (2.0 * nHom1 + nHet) / (2.0 * n);
  const double g  = 1.0 - f;
  const double Eh = 2.0 * f * g * n;
  const double E1 = f * f * n;
  const double E2 = g * g * n;
  if (E1 < 5.0 || Eh < 5.0 || E2 < 5.0)
    return std::numeric_limits<double>::quiet_NaN();
  const double d1 = static_cast<double>(nHom1) - E1;
  const double dh = static_cast<double>(nHet)  - Eh;
  const double d2 = static_cast<double>(nHom2) - E2;
  const double chi2 = d1*d1/E1 + dh*dh/Eh + d2*d2/E2;
  return std::erfc(std::sqrt(chi2 * 0.5));
}

// Exact HWE test (SNPHWE2). O(het_count) time, O(1) auxiliary memory.
// Wigginton JE, Cutler DJ, Abecasis GR (2005). Am J Hum Genet 76:887-893.
static double HweExact(uint32_t obs_hets, uint32_t obs_hom1, uint32_t obs_hom2) {
  const int64_t obs_homc = std::max(obs_hom1, obs_hom2);
  const int64_t obs_homr = std::min(obs_hom1, obs_hom2);
  const int64_t rare = 2 * obs_homr + static_cast<int64_t>(obs_hets);
  const int64_t n    = static_cast<int64_t>(obs_hets) + obs_homc + obs_homr;
  const int64_t obs  = static_cast<int64_t>(obs_hets);

  if (n == 0) return 1.0;

  // Mode of het-count distribution under HWE
  int64_t mid = (rare * (2 * n - rare)) / (2 * n);
  if ((rare & 1) ^ (mid & 1)) ++mid;

  // Adjust mid to sit at the actual peak
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

  double sum = 1.0, p = 0.0, thresh;

  if (obs <= mid) {
    {
      double prob = 1.0;
      int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h > obs; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob;
        ++cr; ++cc;
      }
      thresh = prob;
      p      = thresh;
      for (int64_t h = obs; h >= 2; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob;
        p   += prob;
        ++cr; ++cc;
      }
    }
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
    {
      double prob = 1.0;
      int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h < obs; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob;
        --cr; --cc;
      }
      thresh = prob;
      p      = thresh;
      for (int64_t h = obs; h <= rare - 2; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob;
        p   += prob;
        --cr; --cc;
      }
    }
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

// Compute genotype class counts + derived QC stats from packed BED bytes.
struct GenoStats {
  double   altFreq;
  uint32_t altCounts;
  double   missingRate;
  double   hweP;
  double   maf;
  uint32_t mac;
};

static GenoStats computeStats(const uint8_t* geno_bytes, uint32_t n_samples,
                               bool altFirst, bool exactHwe) {
  uint32_t counts[4];
  count_geno_classes(geno_bytes, n_samples, counts);

  const uint32_t nHomA1   = counts[0];
  const uint32_t nHet     = counts[1];
  const uint32_t nHomA2   = counts[2];
  const uint32_t nMissing = counts[3];
  const uint32_t nonMissing = n_samples - nMissing;

  GenoStats gs;
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
    gs.hweP = (nHet < 5 || nHomA1 < 5 || nHomA2 < 5)
                ? HweExact(nHet, nHomA1, nHomA2)
                : HweChiSq(nHet, nHomA1, nHomA2);
  }
  return gs;
}


// ──────────────────────────────────────────────────────────────────────
// Identity-map fast decode: 2-bit BED → double[n] using 4-wide LUT
//
// When all subjects are used in file order, skip repackSubset entirely
// and decode directly from the raw BED row.  Each byte holds 4 genotypes;
// we decode 4 doubles at a time via a scalar LUT.
// ──────────────────────────────────────────────────────────────────────

static void decodeBedsIdentity(const uint8_t* __restrict bed,
                               uint32_t n,
                               const int* genoMap,
                               double* __restrict out,
                               std::vector<uint32_t>& indexForMissing)
{
  const uint32_t nFullBytes = n / 4;
  uint32_t idx = 0;

  for (uint32_t b = 0; b < nFullBytes; ++b) {
    const uint8_t byte = bed[b];
    // Unrolled 4-sample decode from one byte
    const int g0 = genoMap[(byte      ) & 3];
    const int g1 = genoMap[(byte >>  2) & 3];
    const int g2 = genoMap[(byte >>  4) & 3];
    const int g3 = genoMap[(byte >>  6) & 3];

    // Branchless: write genotype, conditionally push missing index.
    // Predictable branch: missingness is typically < 5%.
    out[idx] = g0; if (g0 < 0) { out[idx] = -1; indexForMissing.push_back(idx); } ++idx;
    out[idx] = g1; if (g1 < 0) { out[idx] = -1; indexForMissing.push_back(idx); } ++idx;
    out[idx] = g2; if (g2 < 0) { out[idx] = -1; indexForMissing.push_back(idx); } ++idx;
    out[idx] = g3; if (g3 < 0) { out[idx] = -1; indexForMissing.push_back(idx); } ++idx;
  }

  // tail samples (last partial byte)
  if (idx < n) {
    const uint8_t byte = bed[nFullBytes];
    for (uint32_t s = 0; idx < n; ++s, ++idx) {
      const int g = genoMap[(byte >> (2*s)) & 3];
      if (g < 0) {
        out[idx] = -1;
        indexForMissing.push_back(idx);
      } else {
        out[idx] = g;
      }
    }
  }
}

// Decode packed 2-bit compact subset buffer → double (scalar, for subset path).
static void decodeCompact(const uint8_t* __restrict compact,
                          uint32_t n,
                          const int* genoMap,
                          double* __restrict out,
                          std::vector<uint32_t>& indexForMissing)
{
  for (uint32_t i = 0; i < n; ++i) {
    const int code = (compact[i >> 2] >> ((i & 3) << 1)) & 3;
    const int g = genoMap[code];
    if (g < 0) {
      out[i] = -1;
      indexForMissing.push_back(i);
    } else {
      out[i] = g;
    }
  }
}

// ──────────────────────────────────────────────────────────────────────
// Simple decode helpers: no QC stats, missing → NaN, no indexForMissing.
// Used by getGenotypesSimple() to skip the count_geno_classes() AVX2 pass.
// ──────────────────────────────────────────────────────────────────────

static const double kNaN = std::numeric_limits<double>::quiet_NaN();

static void decodeBedsIdentitySimple(const uint8_t* __restrict bed,
                                     uint32_t n,
                                     const int* genoMap,
                                     double* __restrict out)
{
  const uint32_t nFullBytes = n / 4;
  uint32_t idx = 0;
  for (uint32_t b = 0; b < nFullBytes; ++b) {
    const uint8_t byte = bed[b];
    const int g0 = genoMap[(byte      ) & 3];
    const int g1 = genoMap[(byte >>  2) & 3];
    const int g2 = genoMap[(byte >>  4) & 3];
    const int g3 = genoMap[(byte >>  6) & 3];
    out[idx++] = g0 < 0 ? kNaN : static_cast<double>(g0);
    out[idx++] = g1 < 0 ? kNaN : static_cast<double>(g1);
    out[idx++] = g2 < 0 ? kNaN : static_cast<double>(g2);
    out[idx++] = g3 < 0 ? kNaN : static_cast<double>(g3);
  }
  if (idx < n) {
    const uint8_t byte = bed[nFullBytes];
    for (uint32_t s = 0; idx < n; ++s, ++idx) {
      const int g = genoMap[(byte >> (2*s)) & 3];
      out[idx] = g < 0 ? kNaN : static_cast<double>(g);
    }
  }
}

static void decodeCompactSimple(const uint8_t* __restrict compact,
                                uint32_t n,
                                const int* genoMap,
                                double* __restrict out)
{
  for (uint32_t i = 0; i < n; ++i) {
    const int code = (compact[i >> 2] >> ((i & 3) << 1)) & 3;
    const int g = genoMap[code];
    out[i] = g < 0 ? kNaN : static_cast<double>(g);
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
    std::vector<std::string> subjData,
    std::string AlleleOrder,
    std::string IDsToIncludeFile,
    std::string RangesToIncludeFile,
    std::string IDsToExcludeFile,
    std::string RangesToExcludeFile,
    int nMarkersEachChunk)
  : m_bedFile(std::move(bedFile)),
    m_altFirst(AlleleOrder == "alt-first"),
    m_identityMap(false)
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
      std::string a1 = tokens[4], a2 = tokens[5];
      std::transform(a1.begin(), a1.end(), a1.begin(), ::toupper);
      std::transform(a2.begin(), a2.end(), a2.begin(), ::toupper);
      if (m_altFirst) {
        m_alt.push_back(std::move(a1));
        m_ref.push_back(std::move(a2));
      } else {
        m_ref.push_back(std::move(a1));
        m_alt.push_back(std::move(a2));
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

  // ---- Build sample position map ----
  m_nSubjUsed = static_cast<uint32_t>(subjData.size());
  std::unordered_map<std::string, uint32_t> famPosLookup;
  famPosLookup.reserve(famIIDs.size());
  for (uint32_t i = 0; i < static_cast<uint32_t>(famIIDs.size()); ++i)
    famPosLookup.emplace(famIIDs[i], i);

  m_samplePosMap.resize(m_nSubjUsed);
  for (uint32_t i = 0; i < m_nSubjUsed; ++i) {
    auto it = famPosLookup.find(subjData[i]);
    if (it == famPosLookup.end())
      throw std::runtime_error("Subject not found in PLINK file: " + subjData[i]);
    m_samplePosMap[i] = it->second;
  }

  // Detect identity mapping: subjects are all subjects in file order.
  // This enables a SIMD fast-path that skips repackSubset entirely.
  if (m_nSubjUsed == m_nSubjInFile) {
    m_identityMap = true;
    for (uint32_t i = 0; i < m_nSubjUsed; ++i) {
      if (m_samplePosMap[i] != i) { m_identityMap = false; break; }
    }
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
    const std::string& RangesToExcludeFile)
{
  const uint32_t nMarkers = static_cast<uint32_t>(chr.size());
  std::vector<PlinkData::MarkerInfo> all;
  all.reserve(nMarkers);
  for (uint32_t i = 0; i < nMarkers; ++i) {
    all.push_back({chr[i], pos[i], markerId[i],
                   ref[i], alt[i], static_cast<uint64_t>(i)});
  }

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

  std::vector<PlinkData::MarkerInfo> filtered;
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


// ══════════════════════════════════════════════════════════════════════
// PlinkCursor
// ══════════════════════════════════════════════════════════════════════

PlinkCursor::PlinkCursor(const std::string& bedFile,
                         uint32_t nBimLines, uint32_t nFamLines,
                         const std::vector<uint32_t>& posMap, bool altFirst,
                         bool identityMap)
  : m_bedFile(bedFile),
    m_nSubjInFile(nFamLines),
    m_bytesPerMarker((nFamLines + 3) / 4),
    m_nMarkers(nBimLines),
    m_posMap(posMap),
    m_altFirst(altFirst),
    m_identityMap(identityMap),
    m_rawBytes(m_bytesPerMarker),
    m_compactGeno(identityMap ? 0 : (posMap.size() + 3) / 4)
{
  m_bedStream.open(bedFile, std::ios::binary);
  if (!m_bedStream.is_open())
    throw std::runtime_error("Cannot open PLINK bed file: " + bedFile);
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
    m_identityMap(other.m_identityMap),
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
  m_bedStream.seekg(3 + static_cast<std::streamoff>(m_bytesPerMarker) * firstMarker);
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

  if (m_blockBytes.size() < static_cast<size_t>(nBytes))
    m_blockBytes.resize(static_cast<size_t>(nBytes));

  if (!m_hasSeqCursor || m_nextMarker != startMarker) {
    m_bedStream.clear();
    m_bedStream.seekg(3 + static_cast<std::streamoff>(m_bytesPerMarker) * startMarker);
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

  const auto filePos = static_cast<std::streamoff>(3 + bpm * gIndex);
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


// ── Identity-map fast path ───────────────────────────────────────────
// All subjects in file order → compute stats and decode from raw BED
// directly.  Saves one full memcpy (repackSubset) per marker, which at
// 500K subjects is ~125 KB saved per marker × millions of markers.

void PlinkCursor::getGenotypesIdentity(
    const uint8_t* raw,
    uint32_t n,
    Eigen::Ref<Eigen::VectorXd> out,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    double& hweP,
    double& maf,
    double& mac,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe)
{
  // Stats directly on the raw BED row (no repack needed)
  GenoStats gs = computeStats(raw, n, m_altFirst, exactHwe);
  altFreq     = gs.altFreq;
  altCounts   = gs.altCounts;
  missingRate = gs.missingRate;
  hweP        = gs.hweP;
  maf         = gs.maf;
  mac         = gs.mac;

  // Decode 2-bit → double directly from raw BED
  const int* genoMap = m_altFirst ? GENO_ALT_FIRST : GENO_REF_FIRST;
  indexForMissing.clear();
  decodeBedsIdentity(raw, n, genoMap, out.data(), indexForMissing);
}


// ── Subset path ──────────────────────────────────────────────────────
// Sparse subject selection → repack into compact buffer, then stats + decode.

void PlinkCursor::getGenotypesSubset(
    const uint8_t* raw,
    uint32_t n,
    Eigen::Ref<Eigen::VectorXd> out,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    double& hweP,
    double& maf,
    double& mac,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe)
{
  repackSubset(raw, m_posMap.data(), m_compactGeno.data(), n);

  GenoStats gs = computeStats(m_compactGeno.data(), n, m_altFirst, exactHwe);
  altFreq     = gs.altFreq;
  altCounts   = gs.altCounts;
  missingRate = gs.missingRate;
  hweP        = gs.hweP;
  maf         = gs.maf;
  mac         = gs.mac;

  const int* genoMap = m_altFirst ? GENO_ALT_FIRST : GENO_REF_FIRST;
  indexForMissing.clear();
  decodeCompact(m_compactGeno.data(), n, genoMap, out.data(), indexForMissing);
}


// ── Public entry point ───────────────────────────────────────────────

void PlinkCursor::getGenotypes(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    double& hweP,
    double& maf,
    double& mac,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe)
{
  const uint8_t* raw = readMarkerPtr(gIndex);
  const uint32_t n = static_cast<uint32_t>(m_posMap.size());

  if (m_identityMap) {
    getGenotypesIdentity(raw, n, out, altFreq, altCounts,
                         missingRate, hweP, maf, mac,
                         indexForMissing, exactHwe);
  } else {
    getGenotypesSubset(raw, n, out, altFreq, altCounts,
                       missingRate, hweP, maf, mac,
                       indexForMissing, exactHwe);
  }
}

// ── Simple entry point: genotype vector only, missing → NaN ──────────
// Skips count_geno_classes (AVX2 counting pass) and indexForMissing
// book-keeping entirely. Use for AF scans where only the genotype values
// are needed.

void PlinkCursor::getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out)
{
  const uint8_t* raw = readMarkerPtr(gIndex);
  const uint32_t n = static_cast<uint32_t>(m_posMap.size());
  const int* genoMap = m_altFirst ? GENO_ALT_FIRST : GENO_REF_FIRST;

  if (m_identityMap) {
    decodeBedsIdentitySimple(raw, n, genoMap, out.data());
  } else {
    repackSubset(raw, m_posMap.data(), m_compactGeno.data(), n);
    decodeCompactSimple(m_compactGeno.data(), n, genoMap, out.data());
  }
}
