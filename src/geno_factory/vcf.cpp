// vcf.cpp — VcfData and VcfCursor implementations
//
// Wraps htslib VCF/BCF reader behind the GenoMeta / GenoCursor interfaces.
// Each VcfCursor owns its own htsFile + header + record for thread safety.

#include "geno_factory/vcf.hpp"
#include "geno_factory/hwe.hpp"
#include "util/logging.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>

extern "C" {
#include <htslib/hts.h>
#include <htslib/vcf.h>
}

// ══════════════════════════════════════════════════════════════════════════════
// VcfData
// ══════════════════════════════════════════════════════════════════════════════

VcfData::VcfData(
    std::string vcfFile,
    const std::vector<uint64_t> &usedMask,
    uint32_t nSamplesInFile,
    uint32_t nUsed,
    std::unordered_set<std::string> chrFilter,
    int nMarkersEachChunk
)
    : m_vcfFile(std::move(vcfFile)),
      m_nSubjInFile(nSamplesInFile),
      m_nSubjUsed(nUsed),
      m_usedMask(usedMask)
{
    m_allUsed = (nUsed == nSamplesInFile);

    // ---- First pass: read all variant metadata ----
    htsFile *fp = hts_open(m_vcfFile.c_str(), "r");
    if (!fp) throw std::runtime_error("Cannot open VCF/BCF file: " + m_vcfFile);

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        hts_close(fp);
        throw std::runtime_error("Cannot read VCF/BCF header: " + m_vcfFile);
    }

    // Verify sample count
    const int nSamplesHeader = bcf_hdr_nsamples(hdr);
    if (static_cast<uint32_t>(nSamplesHeader) != m_nSubjInFile) {
        bcf_hdr_destroy(hdr);
        hts_close(fp);
        throw std::runtime_error("VCF sample count (" + std::to_string(nSamplesHeader) + ") does not match expected (" +
                                 std::to_string(m_nSubjInFile) + ")");
    }

    bcf1_t *rec = bcf_init();
    uint64_t variantIdx = 0;

    while (bcf_read(fp, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_STR); // unpack chrom/pos/id/alleles only

        // Skip non-biallelic
        if (rec->n_allele != 2) continue;

        const char *chromName = bcf_hdr_id2name(hdr, rec->rid);
        std::string chrom = chromName ? chromName : ".";
        uint32_t pos = static_cast<uint32_t>(rec->pos + 1); // htslib is 0-based
        std::string id = rec->d.id ? rec->d.id : ".";
        std::string refAllele = rec->d.allele[0];
        std::string altAllele = rec->d.allele[1];

        m_chr.push_back(chrom);
        m_pos.push_back(pos);
        m_markerId.push_back(id);
        m_ref.push_back(refAllele);
        m_alt.push_back(altAllele);
        if (chrFilter.empty() || chrFilter.count(chrom))
            m_markerInfo.push_back({chrom, pos, id, refAllele, altAllele, variantIdx});
        ++variantIdx;
    }

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    hts_close(fp);

    m_nMarkers = static_cast<uint32_t>(m_markerInfo.size());
    m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);

    infoMsg("  VCF: read %u biallelic variants from %s", m_nMarkers, m_vcfFile.c_str());
}

VcfData::~VcfData() = default;

std::vector<std::vector<uint64_t> > VcfData::buildChunks(
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

std::unique_ptr<GenoCursor> VcfData::makeCursor() const {
    return std::make_unique<VcfCursor>(*this);
}

// ══════════════════════════════════════════════════════════════════════════════
// VcfCursor::Impl
// ══════════════════════════════════════════════════════════════════════════════

struct VcfCursor::Impl {
    htsFile *fp = nullptr;
    bcf_hdr_t *hdr = nullptr;
    bcf1_t *rec = nullptr;

    // bcf_get_genotypes buffer
    int32_t *gtArr = nullptr;
    int ngtArr = 0;

    // bcf_get_format_float buffer for DS
    float *dsArr = nullptr;
    int ndsArr = 0;

    // Sample subsetting
    const std::vector<uint64_t> *usedMask = nullptr;
    uint32_t nSamplesInFile = 0;
    uint32_t nUsed = 0;
    bool allUsed = true;

    // Sequential reading state
    uint64_t currentRecordIdx = 0; // 0-based record index among biallelic variants

    // Advance to the record matching gIndex (index among biallelic variants).
    bool advanceTo(uint64_t gIndex) {
        while (currentRecordIdx <= gIndex) {
            if (bcf_read(fp, hdr, rec) != 0) return false;
            bcf_unpack(rec, BCF_UN_ALL);
            if (rec->n_allele != 2) continue; // skip non-biallelic
            if (currentRecordIdx == gIndex) return true;
            ++currentRecordIdx;
        }
        return false; // should not happen if gIndex is valid
    }

    ~Impl()
    {
        std::free(gtArr);
        std::free(dsArr);
        if (rec) bcf_destroy(rec);
        if (hdr) bcf_hdr_destroy(hdr);
        if (fp) hts_close(fp);
    }

};

VcfCursor::VcfCursor(const VcfData &parent)
    : m_parent(parent),
      m_impl(std::make_unique<Impl>())
{
    auto &impl = *m_impl;
    impl.nSamplesInFile = parent.nSubjInFile();
    impl.nUsed = parent.nSubjUsed();
    impl.allUsed = parent.allUsed();
    impl.usedMask = &parent.usedMask();

    impl.fp = hts_open(parent.vcfFile().c_str(), "r");
    if (!impl.fp) throw std::runtime_error("VCF cursor: cannot open " + parent.vcfFile());

    impl.hdr = bcf_hdr_read(impl.fp);
    if (!impl.hdr) throw std::runtime_error("VCF cursor: cannot read header from " + parent.vcfFile());

    impl.rec = bcf_init();
    impl.currentRecordIdx = 0;
}

VcfCursor::~VcfCursor() = default;

void VcfCursor::beginSequentialBlock(uint64_t /*firstMarker*/) {
    // Re-open or rewind for sequential reading.
    // For simplicity, we rewind to the beginning.
    auto &impl = *m_impl;
    // Close and reopen
    if (impl.rec) {
        bcf_destroy(impl.rec);
        impl.rec = nullptr;
    }
    if (impl.hdr) {
        bcf_hdr_destroy(impl.hdr);
        impl.hdr = nullptr;
    }
    if (impl.fp) {
        hts_close(impl.fp);
        impl.fp = nullptr;
    }

    impl.fp = hts_open(m_parent.vcfFile().c_str(), "r");
    if (!impl.fp) throw std::runtime_error("VCF cursor: cannot reopen file");
    impl.hdr = bcf_hdr_read(impl.fp);
    if (!impl.hdr) throw std::runtime_error("VCF cursor: cannot re-read header");
    impl.rec = bcf_init();
    impl.currentRecordIdx = 0;
}

void VcfCursor::getGenotypes(
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
    const uint32_t nUsed = impl.nUsed;
    indexForMissing.clear();

    if (!impl.advanceTo(gIndex))
        throw std::runtime_error("VCF: failed to advance to variant " + std::to_string(gIndex));

    // Try to extract dosages: prefer DS > GP > GT

    // Try DS field first
    int nds = bcf_get_format_float(impl.hdr, impl.rec, "DS", &impl.dsArr, &impl.ndsArr);
    if (nds > 0) {
        uint32_t outIdx = 0;
        uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0, nMissing = 0;
        for (uint32_t i = 0; i < impl.nSamplesInFile; ++i) {
            bool used = impl.allUsed || (((*impl.usedMask)[i / 64] >> (i % 64)) & 1);
            if (!used) continue;
            float ds = impl.dsArr[i];
            if (bcf_float_is_missing(ds) || bcf_float_is_vector_end(ds)) {
                out[outIdx] = std::numeric_limits<double>::quiet_NaN();
                indexForMissing.push_back(outIdx);
                ++nMissing;
            } else {
                out[outIdx] = static_cast<double>(ds);
                // Classify for HWE
                if (ds < 0.5)
                    ++nHomRef;
                else if (ds > 1.5)
                    ++nHomAlt;
                else
                    ++nHet;
            }
            ++outIdx;
        }
        GenoStats gs = statsFromCounts(nHomRef, nHet, nHomAlt, nMissing, nUsed);
        altFreq = gs.altFreq;
        altCounts = gs.altCounts;
        missingRate = gs.missingRate;
        hweP = gs.hweP;
        maf = gs.maf;
        mac = gs.mac;
        ++impl.currentRecordIdx;
        return;
    }

    // Fall back to GT field
    int ngt = bcf_get_genotypes(impl.hdr, impl.rec, &impl.gtArr, &impl.ngtArr);
    if (ngt <= 0) throw std::runtime_error("VCF: no GT or DS field at variant " + std::to_string(gIndex));

    const int maxPloidy = ngt / static_cast<int>(impl.nSamplesInFile);
    uint32_t outIdx = 0;
    uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0, nMissing = 0;

    for (uint32_t i = 0; i < impl.nSamplesInFile; ++i) {
        bool used = impl.allUsed || (((*impl.usedMask)[i / 64] >> (i % 64)) & 1);
        if (!used) continue;

        const int32_t *ptr = impl.gtArr + i * maxPloidy;
        int dosage = 0;
        bool missing = false;
        int nAlleles = 0;

        for (int j = 0; j < maxPloidy; ++j) {
            if (ptr[j] == bcf_int32_vector_end) break;
            if (bcf_gt_is_missing(ptr[j])) {
                missing = true;
                break;
            }
            dosage += bcf_gt_allele(ptr[j]); // 0 for ref, 1 for alt1
            ++nAlleles;
        }

        if (missing || nAlleles == 0) {
            out[outIdx] = std::numeric_limits<double>::quiet_NaN();
            indexForMissing.push_back(outIdx);
            ++nMissing;
        } else {
            out[outIdx] = static_cast<double>(dosage);
            if (dosage == 0)
                ++nHomRef;
            else if (dosage == 2)
                ++nHomAlt;
            else
                ++nHet;
        }
        ++outIdx;
    }

    GenoStats gs = statsFromCounts(nHomRef, nHet, nHomAlt, nMissing, nUsed);
    altFreq = gs.altFreq;
    altCounts = gs.altCounts;
    missingRate = gs.missingRate;
    hweP = gs.hweP;
    maf = gs.maf;
    mac = gs.mac;
    ++impl.currentRecordIdx;
}

void VcfCursor::getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out
) {
    auto &impl = *m_impl;

    if (!impl.advanceTo(gIndex))
        throw std::runtime_error("VCF: failed to advance to variant " + std::to_string(gIndex));

    // Try DS first
    int nds = bcf_get_format_float(impl.hdr, impl.rec, "DS", &impl.dsArr, &impl.ndsArr);
    if (nds > 0) {
        uint32_t outIdx = 0;
        for (uint32_t i = 0; i < impl.nSamplesInFile; ++i) {
            bool used = impl.allUsed || (((*impl.usedMask)[i / 64] >> (i % 64)) & 1);
            if (!used) continue;
            float ds = impl.dsArr[i];
            out[outIdx] = (bcf_float_is_missing(ds) || bcf_float_is_vector_end(ds))
                              ? std::numeric_limits<double>::quiet_NaN()
                              : static_cast<double>(ds);
            ++outIdx;
        }
        ++impl.currentRecordIdx;
        return;
    }

    // Fall back to GT
    int ngt = bcf_get_genotypes(impl.hdr, impl.rec, &impl.gtArr, &impl.ngtArr);
    if (ngt <= 0) throw std::runtime_error("VCF: no GT or DS field at variant " + std::to_string(gIndex));

    const int maxPloidy = ngt / static_cast<int>(impl.nSamplesInFile);
    uint32_t outIdx = 0;

    for (uint32_t i = 0; i < impl.nSamplesInFile; ++i) {
        bool used = impl.allUsed || (((*impl.usedMask)[i / 64] >> (i % 64)) & 1);
        if (!used) continue;

        const int32_t *ptr = impl.gtArr + i * maxPloidy;
        int dosage = 0;
        bool missing = false;

        for (int j = 0; j < maxPloidy; ++j) {
            if (ptr[j] == bcf_int32_vector_end) break;
            if (bcf_gt_is_missing(ptr[j])) {
                missing = true;
                break;
            }
            dosage += bcf_gt_allele(ptr[j]);
        }

        out[outIdx] = missing ? std::numeric_limits<double>::quiet_NaN() : static_cast<double>(dosage);
        ++outIdx;
    }
    ++impl.currentRecordIdx;
}
