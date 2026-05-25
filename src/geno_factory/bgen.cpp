// bgen.cpp — BgenData and BgenCursor implementations
//
// Wraps the bgen reference implementation behind GenoMeta / GenoCursor.
// Each BgenCursor owns its own ifstream for thread safety.

#include "geno_factory/bgen.hpp"
#include "geno_factory/hwe.hpp"
#include "geno_factory/variant_filter.hpp"
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
// DosageSetter — callback that bgen's read_and_parse_genotype_data_block uses
// to deliver per-sample genotype probabilities.  We convert on the fly to a
// per-subject dosage, applying the usedMask.
//
// Supports both BGEN-1.2 unphased and phased layouts:
//
//   Unphased biallelic diploid (Z=3, ePerUnorderedGenotype):
//     idx 0 = P(AA)  [hom first allele]
//     idx 1 = P(AB)  [het]
//     idx 2 = P(BB)  [hom second allele]
//     dosage of second allele = idx_1 + 2*idx_2
//
//   Phased biallelic diploid (Z=4, ePerPhasedHaplotypePerAllele):
//     idx 0 = P(hap0 = first allele)
//     idx 1 = P(hap0 = second allele)
//     idx 2 = P(hap1 = first allele)
//     idx 3 = P(hap1 = second allele)
//     dosage of second allele = idx_1 + idx_3
//
// `countFirstAllele = false` (default) → dosage of second allele (alleles[1]).
// `countFirstAllele = true`            → dosage of first allele  (alleles[0]),
//   used by --bgen FILE ref-last (or ref-unknown) to compensate for plink2's
//   default BGEN export which writes ALT as the first allele.
// ──────────────────────────────────────────────────────────────────────────────

struct DosageSetter {
    // Outputs
    double *out = nullptr;
    std::vector<uint32_t> *missingIdx = nullptr; // may be nullptr for simple path

    // Config
    const std::vector<uint64_t> *usedMask = nullptr;
    uint32_t nSamplesInFile = 0;
    bool allUsed = true;
    bool countFirstAllele = false; // true ⇒ count alleles[0], else count alleles[1]

    // Per-sample state
    uint32_t outIdx = 0;
    uint32_t currentSample = 0;
    bool currentUsed = true;
    bool phasedCurrent = false;
    uint32_t expectedEntries = 0;
    uint32_t probIdx = 0;
    double probSum = 0.0;       // running dosage of the "counted" allele
    bool sampleMissing = false; // any MissingValue seen for this sample
    bool sampleFinalized = false;

    // Counts
    uint32_t nHomRef = 0;
    uint32_t nHet = 0;
    uint32_t nHomAlt = 0;
    uint32_t nMissing = 0;

    void initialise(
        std::size_t /*nsamples*/,
        std::size_t                                       /*nalleles*/
    ) {
        outIdx = 0;
        nHomRef = nHet = nHomAlt = nMissing = 0;
    }

    bool set_sample(std::size_t i) {
        currentSample = static_cast<uint32_t>(i);
        currentUsed = allUsed || (((*usedMask)[currentSample / 64] >> (currentSample % 64)) & 1);
        probIdx = 0;
        probSum = 0.0;
        sampleMissing = false;
        sampleFinalized = false;
        return true; // always receive data to advance correctly
    }

    void set_number_of_entries(
        std::size_t /*ploidy*/,
        std::size_t Z,
        genfile::OrderType order,
        genfile::ValueType                        /*valueType*/
    ) {
        phasedCurrent = (order == genfile::ePerPhasedHaplotypePerAllele);
        expectedEntries = static_cast<uint32_t>(Z);
    }

    void finalizeSample() {
        if (sampleFinalized) return;
        sampleFinalized = true;
        if (!currentUsed) return;
        if (sampleMissing) {
            out[outIdx] = std::numeric_limits<double>::quiet_NaN();
            if (missingIdx) missingIdx->push_back(outIdx);
            ++nMissing;
        } else {
            double dosage = countFirstAllele ? (2.0 - probSum) : probSum;
            // Clamp to [0, 2] in case of small numerical drift in phased branch.
            if (dosage < 0.0) dosage = 0.0;
            else if (dosage > 2.0) dosage = 2.0;
            out[outIdx] = dosage;
            if (dosage < 0.5)
                ++nHomRef;
            else if (dosage > 1.5)
                ++nHomAlt;
            else
                ++nHet;
        }
        ++outIdx;
    }

    void set_value(
        uint32_t idx,
        double value
    ) {
        if (currentUsed && !sampleMissing) {
            if (phasedCurrent) {
                // Z=4 phased: idx 1 = P(hap0=allele1), idx 3 = P(hap1=allele1).
                if (idx == 1 || idx == 3) probSum += value;
            } else {
                // Z=3 unphased: idx 1 = P(het), idx 2 = P(hom allele1).
                if (idx == 1)
                    probSum += value;
                else if (idx == 2)
                    probSum += 2.0 * value;
            }
        }
        ++probIdx;
        if (probIdx >= expectedEntries) finalizeSample();
    }

    void set_value(
        uint32_t /*idx*/,
        genfile::MissingValue
    ) {
        sampleMissing = true;
        ++probIdx;
        if (probIdx >= expectedEntries) finalizeSample();
    }

    void finalise() {
    }

};

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════════════
// BgenData
// ══════════════════════════════════════════════════════════════════════════════

BgenData::BgenData(
    std::string bgenFile,
    const std::vector<uint64_t> &usedMask,
    uint32_t nSamplesInFile,
    uint32_t nUsed,
    std::unordered_set<std::string> chrFilter,
    std::string extractFile,
    std::string excludeFile,
    int nMarkersEachChunk,
    bool altFirst
)
    : m_bgenFile(std::move(bgenFile)),
      m_nSubjInFile(nSamplesInFile),
      m_nSubjUsed(nUsed),
      m_usedMask(usedMask),
      m_altFirst(altFirst)
{
    m_allUsed = (nUsed == nSamplesInFile);

    // ---- Open and read header ----
    std::ifstream stream(m_bgenFile, std::ios::binary);
    if (!stream.is_open()) throw std::runtime_error("Cannot open BGEN file: " + m_bgenFile);

    genfile::bgen::uint32_t offset;
    genfile::bgen::read_offset(stream, &offset);

    genfile::bgen::Context context;
    genfile::bgen::read_header_block(stream, &context);

    if (context.number_of_samples != m_nSubjInFile)
        throw std::runtime_error("BGEN sample count (" + std::to_string(context.number_of_samples) +
                                 ") does not match expected (" + std::to_string(m_nSubjInFile) + ")");

    m_bgenFlags = context.flags;

    // Skip sample identifiers if present
    if (context.flags & genfile::bgen::e_SampleIdentifiers) {
        genfile::bgen::read_sample_identifier_block(stream, context, [](std::string const &) {
        });
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
               stream, context, &SNPID, &RSID, &chromosome, &position, [&](std::size_t n) {
        alleles.resize(n);
    },
               [&](std::size_t i, std::string const &a) {
        alleles[i] = a;
    })) {

        // Skip non-biallelic
        if (alleles.size() != 2) {
            genfile::bgen::ignore_genotype_data_block(stream, context);
            continue;
        }

        m_chr.push_back(chromosome);
        m_pos.push_back(position);
        // Use RSID as marker ID; fall back to SNPID
        m_markerId.push_back(RSID.empty() ? SNPID : RSID);
        // BGEN itself only encodes "first allele" / "second allele".  The
        // user declares the REF/ALT convention at the CLI via
        // --bgen FILE {ref-first|ref-last|ref-unknown}:
        //   ref-first   → alleles[0] = REF (m_altFirst = false; default of
        //                 IMPUTE / qctool / UK Biobank)
        //   ref-last    → alleles[0] = ALT (m_altFirst = true; default of
        //                 plink2 --export bgen-1.x)
        //   ref-unknown → treated as ref-last with a CLI-level warning.
        const std::string &refAllele = m_altFirst ? alleles[1] : alleles[0];
        const std::string &altAllele = m_altFirst ? alleles[0] : alleles[1];
        m_ref.push_back(refAllele);
        m_alt.push_back(altAllele);
        if (chrFilter.empty() || chrFilter.count(chromosome))
            m_markerInfo.push_back({chromosome, position, RSID.empty() ? SNPID : RSID,
                                    refAllele, altAllele,
                                    variantIdx});
        ++variantIdx;

        genfile::bgen::ignore_genotype_data_block(stream, context);
    }

    // ---- Apply --extract / --exclude filter ----
    geno_factory::filterMarkersByIds(m_markerInfo, extractFile, excludeFile);

    m_nMarkers = static_cast<uint32_t>(m_markerInfo.size());
    m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);

    infoMsg("  BGEN: read %u biallelic variants from %s", m_nMarkers, m_bgenFile.c_str());
}

BgenData::~BgenData() = default;

std::vector<std::vector<uint64_t> > BgenData::buildChunks(
    const std::vector<MarkerInfo> &markers,
    int chunkSize
) {
    std::vector<std::vector<uint64_t> > chunks;
    if (markers.empty()) return chunks;
    std::vector<uint64_t> cur;
    cur.reserve(chunkSize);
    std::string prevChr = markers[0].chrom;
    for (uint32_t i = 0; i < markers.size(); ++i) {
        if (markers[i].chrom != prevChr || static_cast<int>(cur.size()) >= chunkSize) {
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
    std::ifstream stream;
    genfile::bgen::Context context;

    // Scratch buffers for decompression
    std::vector<genfile::byte_t> buffer1;
    std::vector<genfile::byte_t> buffer2;

    // Sample subsetting
    const std::vector<uint64_t> *usedMask = nullptr;
    uint32_t nSamplesInFile = 0;
    uint32_t nUsed = 0;
    bool allUsed = true;
    bool altFirst = false;

    // Sequential reading state
    uint64_t currentBiallelicIdx = 0; // index among biallelic variants seen so far

    // Advance to the biallelic variant at gIndex
    bool advanceTo(uint64_t gIndex) {
        std::string SNPID, RSID, chromosome;
        genfile::bgen::uint32_t position;
        std::vector<std::string> alleles;

        while (currentBiallelicIdx <= gIndex) {
            if (!genfile::bgen::read_snp_identifying_data(
                    stream, context, &SNPID, &RSID, &chromosome, &position, [&](std::size_t n) {
                alleles.resize(n);
            },
                    [&](std::size_t i, std::string const &a) {
                alleles[i] = a;
            }))
                return false;

            if (alleles.size() != 2) {
                genfile::bgen::ignore_genotype_data_block(stream, context);
                continue;
            }

            if (currentBiallelicIdx == gIndex) return true;

            genfile::bgen::ignore_genotype_data_block(stream, context);
            ++currentBiallelicIdx;
        }
        return false;
    }

};

BgenCursor::BgenCursor(const BgenData &parent)
    : m_parent(parent),
      m_impl(std::make_unique<Impl>())
{
    auto &impl = *m_impl;
    impl.nSamplesInFile = parent.nSubjInFile();
    impl.nUsed = parent.nSubjUsed();
    impl.allUsed = parent.allUsed();
    impl.usedMask = &parent.usedMask();
    impl.altFirst = parent.altFirst();

    // Open file and read header to build context
    impl.stream.open(parent.bgenFile(), std::ios::binary);
    if (!impl.stream.is_open()) throw std::runtime_error("BGEN cursor: cannot open " + parent.bgenFile());

    genfile::bgen::uint32_t offset;
    genfile::bgen::read_offset(impl.stream, &offset);
    genfile::bgen::read_header_block(impl.stream, &impl.context);

    // Skip sample identifiers
    if (impl.context.flags & genfile::bgen::e_SampleIdentifiers) {
        genfile::bgen::read_sample_identifier_block(impl.stream, impl.context, [](std::string const &) {
        });
    }

    impl.stream.seekg(offset + 4);
    impl.currentBiallelicIdx = 0;
}

BgenCursor::~BgenCursor() = default;

void BgenCursor::beginSequentialBlock(uint64_t /*firstMarker*/) {
    // Reopen and seek to data start
    auto &impl = *m_impl;
    impl.stream.close();
    impl.stream.open(m_parent.bgenFile(), std::ios::binary);
    if (!impl.stream.is_open()) throw std::runtime_error("BGEN cursor: cannot reopen file");

    genfile::bgen::uint32_t offset;
    genfile::bgen::read_offset(impl.stream, &offset);
    genfile::bgen::read_header_block(impl.stream, &impl.context);

    if (impl.context.flags & genfile::bgen::e_SampleIdentifiers) {
        genfile::bgen::read_sample_identifier_block(impl.stream, impl.context, [](std::string const &) {
        });
    }

    impl.stream.seekg(offset + 4);
    impl.currentBiallelicIdx = 0;
}

void BgenCursor::getGenotypes(
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
    auto &impl = *m_impl;
    indexForMissing.clear();

    if (!impl.advanceTo(gIndex))
        throw std::runtime_error("BGEN: failed to advance to variant " + std::to_string(gIndex));

    DosageSetter setter;
    setter.out = out.data();
    setter.missingIdx = &indexForMissing;
    setter.usedMask = impl.usedMask;
    setter.nSamplesInFile = impl.nSamplesInFile;
    setter.allUsed = impl.allUsed;
    setter.countFirstAllele = impl.altFirst;

    genfile::bgen::read_and_parse_genotype_data_block(impl.stream, impl.context, setter, &impl.buffer1, &impl.buffer2);

    const uint32_t nUsed = impl.nUsed;
    // statsFromCounts treats its first count argument as "hom of the ALT
    // allele" and returns altCounts = 2*nHomAlt + nHet.  Our DosageSetter
    // already classified the per-sample dosage of the ALT allele
    // (alleles[1] under --bgen FILE ref-first, alleles[0] under ref-last
    // or ref-unknown), so setter.nHomAlt is the count of subjects with
    // dosage ≈ 2 (hom ALT).
    GenoStats gs = statsFromCounts(setter.nHomAlt, setter.nHet, setter.nHomRef, setter.nMissing, nUsed);
    altFreq = gs.altFreq;
    altCounts = gs.altCounts;
    missingRate = gs.missingRate;
    hweP = gs.hweP;
    maf = gs.maf;
    mac = gs.mac;

    ++impl.currentBiallelicIdx;
}

void BgenCursor::getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out
) {
    auto &impl = *m_impl;

    if (!impl.advanceTo(gIndex))
        throw std::runtime_error("BGEN: failed to advance to variant " + std::to_string(gIndex));

    DosageSetter setter;
    setter.out = out.data();
    setter.missingIdx = nullptr;
    setter.usedMask = impl.usedMask;
    setter.nSamplesInFile = impl.nSamplesInFile;
    setter.allUsed = impl.allUsed;
    setter.countFirstAllele = impl.altFirst;

    genfile::bgen::read_and_parse_genotype_data_block(impl.stream, impl.context, setter, &impl.buffer1, &impl.buffer2);

    ++impl.currentBiallelicIdx;
}
