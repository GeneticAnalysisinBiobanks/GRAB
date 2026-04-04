// bgen.cpp — BgenData and BgenCursor implementations
//
// Wraps the bgen reference implementation behind GenoMeta / GenoCursor.
// Each BgenCursor owns its own ifstream for thread safety.

#include "io/bgen.hpp"
#include "util/logging.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "genfile/bgen/bgen.hpp"

namespace {

// ──────────────────────────────────────────────────────────────────────────────
// HWE helpers — same algorithms as plink.cpp
// ──────────────────────────────────────────────────────────────────────────────

static double HweChiSq(uint32_t nHet, uint32_t nHom1, uint32_t nHom2) {
  const double n = static_cast<double>(nHet + nHom1 + nHom2);
  const double p = (2.0 * nHom1 + nHet) / (2.0 * n);
  const double q = 1.0 - p;
  const double E1 = n * p * p;
  const double Eh = 2.0 * n * p * q;
  const double E2 = n * q * q;
  if (E1 < 1e-15 || Eh < 1e-15 || E2 < 1e-15) return 1.0;
  const double d1 = static_cast<double>(nHom1) - E1;
  const double dh = static_cast<double>(nHet)  - Eh;
  const double d2 = static_cast<double>(nHom2) - E2;
  const double chi2 = d1*d1/E1 + dh*dh/Eh + d2*d2/E2;
  return std::erfc(std::sqrt(chi2 * 0.5));
}

static double HweExact(uint32_t obs_hets, uint32_t obs_hom1, uint32_t obs_hom2) {
  const int64_t obs_homc = std::max(obs_hom1, obs_hom2);
  const int64_t obs_homr = std::min(obs_hom1, obs_hom2);
  const int64_t rare = 2 * obs_homr + static_cast<int64_t>(obs_hets);
  const int64_t n    = static_cast<int64_t>(obs_hets) + obs_homc + obs_homr;
  const int64_t obs  = static_cast<int64_t>(obs_hets);
  if (n == 0) return 1.0;
  int64_t mid = (rare * (2 * n - rare)) / (2 * n);
  if ((rare & 1) ^ (mid & 1)) ++mid;
  {
    int64_t hr = (rare - mid) / 2;
    int64_t hc = n - mid - hr;
    if (mid + 2 <= rare && hr > 0 &&
        4.0 * hr * hc > (mid + 2.0) * (mid + 1.0))
      mid += 2;
    else if (mid >= 2 &&
             static_cast<double>(mid) * (mid - 1) > 4.0 * (hr + 1.0) * (hc + 1.0))
      mid -= 2;
  }
  const int64_t mid_homr = (rare - mid) / 2;
  const int64_t mid_homc = n - mid - mid_homr;
  double sum = 1.0, p = 0.0, thresh;
  if (obs <= mid) {
    {
      double prob = 1.0; int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h > obs; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob; ++cr; ++cc;
      }
      thresh = prob; p = thresh;
      for (int64_t h = obs; h >= 2; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob; p += prob; ++cr; ++cc;
      }
    }
    {
      double prob = 1.0; int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h <= rare - 2; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob; if (prob <= thresh) p += prob;
        --cr; --cc;
      }
    }
  } else {
    {
      double prob = 1.0; int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h < obs; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob; --cr; --cc;
      }
      thresh = prob; p = thresh;
      for (int64_t h = obs; h <= rare - 2; h += 2) {
        prob *= 4.0 * cr * cc / ((h + 2.0) * (h + 1.0));
        sum += prob; p += prob; --cr; --cc;
      }
    }
    {
      double prob = 1.0; int64_t cr = mid_homr, cc = mid_homc;
      for (int64_t h = mid; h >= 2; h -= 2) {
        prob *= h * (h - 1.0) / (4.0 * (cr + 1.0) * (cc + 1.0));
        sum += prob; if (prob <= thresh) p += prob;
        ++cr; ++cc;
      }
    }
  }
  return std::min(p / sum, 1.0);
}

struct GenoStats {
  double   altFreq;
  uint32_t altCounts;
  double   missingRate;
  double   hweP;
  double   maf;
  uint32_t mac;
};

static GenoStats statsFromCounts(uint32_t nHomRef, uint32_t nHet, uint32_t nHomAlt,
                                  uint32_t nMissing, uint32_t nSamples,
                                  bool exactHwe) {
  const uint32_t nonMissing = nSamples - nMissing;
  GenoStats gs;
  gs.altCounts = 2 * nHomAlt + nHet;
  gs.missingRate = static_cast<double>(nMissing) / nSamples;
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
  if (exactHwe)
    gs.hweP = HweExact(nHet, nHomAlt, nHomRef);
  else
    gs.hweP = (nHet < 5 || nHomAlt < 5 || nHomRef < 5)
              ? HweExact(nHet, nHomAlt, nHomRef)
              : HweChiSq(nHet, nHomAlt, nHomRef);
  return gs;
}


// ──────────────────────────────────────────────────────────────────────────────
// DosageSetter — callback that bgen's read_and_parse_genotype_data_block uses
// to deliver per-sample genotype probabilities.  We convert on the fly to
// ALT dosage, applying the usedMask.
// ──────────────────────────────────────────────────────────────────────────────

struct DosageSetter {
  // Outputs
  double*                out         = nullptr;
  std::vector<uint32_t>* missingIdx  = nullptr;  // may be nullptr for simple path

  // Config
  const std::vector<uint64_t>* usedMask = nullptr;
  uint32_t nSamplesInFile = 0;
  bool     allUsed        = true;

  // Per-sample state
  uint32_t outIdx         = 0;
  uint32_t currentSample  = 0;
  bool     currentUsed    = true;
  uint32_t probIdx        = 0;
  double   probSum        = 0.0;      // P(AB) + 2*P(BB)

  // Counts
  uint32_t nHomRef  = 0;
  uint32_t nHet     = 0;
  uint32_t nHomAlt  = 0;
  uint32_t nMissing = 0;

  void initialise(std::size_t /*nsamples*/, std::size_t /*nalleles*/) {
    outIdx = 0;
    nHomRef = nHet = nHomAlt = nMissing = 0;
  }

  bool set_sample(std::size_t i) {
    currentSample = static_cast<uint32_t>(i);
    currentUsed = allUsed ||
      (((*usedMask)[currentSample / 64] >> (currentSample % 64)) & 1);
    probIdx = 0;
    probSum = 0.0;
    return true;  // always receive data to advance correctly
  }

  void set_number_of_entries(
      std::size_t /*ploidy*/, std::size_t /*Z*/,
      genfile::OrderType /*order*/, genfile::ValueType /*valueType*/) {}

  void set_value(uint32_t idx, double value) {
    if (!currentUsed) return;
    // For biallelic diploid: idx 0 = P(RR), idx 1 = P(RA), idx 2 = P(AA)
    if (idx == 1) probSum += value;
    else if (idx == 2) probSum += 2.0 * value;
    ++probIdx;

    if (idx == 2) {
      // Final probability for this sample
      out[outIdx] = probSum;
      if (probSum < 0.5) ++nHomRef;
      else if (probSum > 1.5) ++nHomAlt;
      else ++nHet;
      ++outIdx;
    }
  }

  void set_value(uint32_t /*idx*/, genfile::MissingValue) {
    if (!currentUsed) return;
    ++probIdx;
    // Mark missing at last index
    if (probIdx >= 3) {
      out[outIdx] = std::numeric_limits<double>::quiet_NaN();
      if (missingIdx) missingIdx->push_back(outIdx);
      ++nMissing;
      ++outIdx;
    }
  }

  void finalise() {}
};

} // anonymous namespace


// ══════════════════════════════════════════════════════════════════════════════
// BgenData
// ══════════════════════════════════════════════════════════════════════════════

BgenData::BgenData(
    std::string bgenFile,
    const std::vector<uint64_t>& usedMask,
    uint32_t nSamplesInFile,
    uint32_t nUsed,
    int nMarkersEachChunk)
  : m_bgenFile(std::move(bgenFile)),
    m_nSubjInFile(nSamplesInFile),
    m_nSubjUsed(nUsed),
    m_usedMask(usedMask
) {
  m_allUsed = (nUsed == nSamplesInFile);

  // ---- Open and read header ----
  std::ifstream stream(m_bgenFile, std::ios::binary);
  if (!stream.is_open())
    throw std::runtime_error("Cannot open BGEN file: " + m_bgenFile);

  genfile::bgen::uint32_t offset;
  genfile::bgen::read_offset(stream, &offset);

  genfile::bgen::Context context;
  genfile::bgen::read_header_block(stream, &context);

  if (context.number_of_samples != m_nSubjInFile)
    throw std::runtime_error("BGEN sample count (" +
      std::to_string(context.number_of_samples) +
      ") does not match expected (" +
      std::to_string(m_nSubjInFile) + ")");

  m_bgenFlags = context.flags;

  // Skip sample identifiers if present
  if (context.flags & genfile::bgen::e_SampleIdentifiers) {
    genfile::bgen::read_sample_identifier_block(
      stream, context, [](std::string const&) {});
  }

  // Seek to first variant
  stream.seekg(offset + 4);
  m_dataOffset = stream.tellg();

  // ---- First pass: read variant metadata (skip genotype data) ----
  std::string SNPID, RSID, chromosome;
  genfile::bgen::uint32_t position;
  std::vector<std::string> alleles;
  uint64_t variantIdx = 0;

  while (genfile::bgen::read_snp_identifying_data(
    stream, context, &SNPID, &RSID, &chromosome, &position,
    [&](std::size_t n) { alleles.resize(n); },
    [&](std::size_t i, std::string const& a) { alleles[i] = a; })) {

    // Skip non-biallelic
    if (alleles.size() != 2) {
      genfile::bgen::ignore_genotype_data_block(stream, context);
      continue;
    }

    m_chr.push_back(chromosome);
    m_pos.push_back(position);
    // Use RSID as marker ID; fall back to SNPID
    m_markerId.push_back(RSID.empty() ? SNPID : RSID);
    // BGEN convention: allele[0] = first allele (typically REF), allele[1] = ALT
    m_ref.push_back(alleles[0]);
    m_alt.push_back(alleles[1]);
    m_markerInfo.push_back({chromosome, position,
                            RSID.empty() ? SNPID : RSID,
                            alleles[0], alleles[1], variantIdx});
    ++variantIdx;

    genfile::bgen::ignore_genotype_data_block(stream, context);
  }

  m_nMarkers = static_cast<uint32_t>(m_markerInfo.size());
  m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);

  infoMsg("  BGEN: read %u biallelic variants from %s",
          m_nMarkers, m_bgenFile.c_str());
}

BgenData::~BgenData() = default;

std::vector<std::vector<uint64_t>> BgenData::buildChunks(
    const std::vector<MarkerInfo>& markers, int chunkSize) {
  std::vector<std::vector<uint64_t>> chunks;
  if (markers.empty()) return chunks;
  std::vector<uint64_t> cur;
  cur.reserve(chunkSize);
  std::string prevChr = markers[0].chrom;
  for (uint32_t i = 0; i < markers.size(); ++i) {
    if (markers[i].chrom != prevChr ||
        static_cast<int>(cur.size()) >= chunkSize) {
      if (!cur.empty()) chunks.push_back(std::move(cur));
      cur.clear();
      cur.reserve(chunkSize);
      prevChr = markers[i].chrom;
    }
    cur.push_back(markers[i].genoIndex);
  }
  if (!cur.empty()) chunks.push_back(std::move(cur));
  return chunks;
}

std::unique_ptr<GenoCursor> BgenData::makeCursor() const {
  return std::make_unique<BgenCursor>(*this);
}


// ══════════════════════════════════════════════════════════════════════════════
// BgenCursor::Impl
// ══════════════════════════════════════════════════════════════════════════════

struct BgenCursor::Impl {
  std::ifstream          stream;
  genfile::bgen::Context context;

  // Scratch buffers for decompression
  std::vector<genfile::byte_t> buffer1;
  std::vector<genfile::byte_t> buffer2;

  // Sample subsetting
  const std::vector<uint64_t>* usedMask = nullptr;
  uint32_t nSamplesInFile = 0;
  uint32_t nUsed          = 0;
  bool     allUsed        = true;

  // Sequential reading state
  uint64_t currentBiallelicIdx = 0;  // index among biallelic variants seen so far

  // Advance to the biallelic variant at gIndex
  bool advanceTo(uint64_t gIndex) {
    std::string SNPID, RSID, chromosome;
    genfile::bgen::uint32_t position;
    std::vector<std::string> alleles;

    while (currentBiallelicIdx <= gIndex) {
      if (!genfile::bgen::read_snp_identifying_data(
          stream, context, &SNPID, &RSID, &chromosome, &position,
          [&](std::size_t n) { alleles.resize(n); },
          [&](std::size_t i, std::string const& a) { alleles[i] = a; }))
        return false;

      if (alleles.size() != 2) {
        genfile::bgen::ignore_genotype_data_block(stream, context);
        continue;
      }

      if (currentBiallelicIdx == gIndex)
        return true;

      genfile::bgen::ignore_genotype_data_block(stream, context);
      ++currentBiallelicIdx;
    }
    return false;
  }
};


BgenCursor::BgenCursor(const BgenData& parent)
  : m_parent(parent), m_impl(std::make_unique<Impl>()
) {
  auto& impl = *m_impl;
  impl.nSamplesInFile = parent.nSubjInFile();
  impl.nUsed          = parent.nSubjUsed();
  impl.allUsed        = parent.allUsed();
  impl.usedMask       = &parent.usedMask();

  // Open file and read header to build context
  impl.stream.open(parent.bgenFile(), std::ios::binary);
  if (!impl.stream.is_open())
    throw std::runtime_error("BGEN cursor: cannot open " + parent.bgenFile());

  genfile::bgen::uint32_t offset;
  genfile::bgen::read_offset(impl.stream, &offset);
  genfile::bgen::read_header_block(impl.stream, &impl.context);

  // Skip sample identifiers
  if (impl.context.flags & genfile::bgen::e_SampleIdentifiers) {
    genfile::bgen::read_sample_identifier_block(
      impl.stream, impl.context, [](std::string const&) {});
  }

  impl.stream.seekg(offset + 4);
  impl.currentBiallelicIdx = 0;
}

BgenCursor::~BgenCursor() = default;

void BgenCursor::beginSequentialBlock(uint64_t /*firstMarker*/) {
  // Reopen and seek to data start
  auto& impl = *m_impl;
  impl.stream.close();
  impl.stream.open(m_parent.bgenFile(), std::ios::binary);
  if (!impl.stream.is_open())
    throw std::runtime_error("BGEN cursor: cannot reopen file");

  genfile::bgen::uint32_t offset;
  genfile::bgen::read_offset(impl.stream, &offset);
  genfile::bgen::read_header_block(impl.stream, &impl.context);

  if (impl.context.flags & genfile::bgen::e_SampleIdentifiers) {
    genfile::bgen::read_sample_identifier_block(
      impl.stream, impl.context, [](std::string const&) {});
  }

  impl.stream.seekg(offset + 4);
  impl.currentBiallelicIdx = 0;
}

void BgenCursor::getGenotypes(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out,
    double& altFreq,
    double& altCounts,
    double& missingRate,
    double& hweP,
    double& maf,
    double& mac,
    std::vector<uint32_t>& indexForMissing,
    bool exactHwe
) {
  auto& impl = *m_impl;
  indexForMissing.clear();

  if (!impl.advanceTo(gIndex))
    throw std::runtime_error("BGEN: failed to advance to variant " +
                             std::to_string(gIndex));

  DosageSetter setter;
  setter.out            = out.data();
  setter.missingIdx     = &indexForMissing;
  setter.usedMask       = impl.usedMask;
  setter.nSamplesInFile = impl.nSamplesInFile;
  setter.allUsed        = impl.allUsed;

  genfile::bgen::read_and_parse_genotype_data_block(
    impl.stream, impl.context, setter, &impl.buffer1, &impl.buffer2);

  const uint32_t nUsed = impl.nUsed;
  GenoStats gs = statsFromCounts(setter.nHomRef, setter.nHet, setter.nHomAlt,
                                  setter.nMissing, nUsed, exactHwe);
  altFreq     = gs.altFreq;
  altCounts   = gs.altCounts;
  missingRate = gs.missingRate;
  hweP        = gs.hweP;
  maf         = gs.maf;
  mac         = gs.mac;

  ++impl.currentBiallelicIdx;
}

void BgenCursor::getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out
) {
  auto& impl = *m_impl;

  if (!impl.advanceTo(gIndex))
    throw std::runtime_error("BGEN: failed to advance to variant " +
                             std::to_string(gIndex));

  DosageSetter setter;
  setter.out            = out.data();
  setter.missingIdx     = nullptr;
  setter.usedMask       = impl.usedMask;
  setter.nSamplesInFile = impl.nSamplesInFile;
  setter.allUsed        = impl.allUsed;

  genfile::bgen::read_and_parse_genotype_data_block(
    impl.stream, impl.context, setter, &impl.buffer1, &impl.buffer2);

  ++impl.currentBiallelicIdx;
}
