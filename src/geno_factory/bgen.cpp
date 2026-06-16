// bgen.cpp — BgenData and BgenCursor implementations
//
// Wraps the bgen reference implementation behind GenoMeta / GenoCursor.
// Each BgenCursor owns its own ifstream for thread safety.

#include "geno_factory/bgen.hpp"
#include "geno_factory/hwe.hpp"
#include "geno_factory/variant_filter.hpp"
#include "util/logging.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "genfile/bgen/bgen.hpp"
#include "libdeflate.h"
#include "sqlite3.h"
#include "zstd.h"

#include <sys/stat.h>

namespace {

// True if a regular file exists and its size (bytes) is returned via *size.
inline bool statFileSize(const std::string &path, uint64_t *size) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return false;
    if (size) *size = static_cast<uint64_t>(st.st_size);
    return true;
}

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
    double dosageSum = 0.0;  // Σ dosage over non-missing used samples (for dosage AF)

    void initialise(
        std::size_t /*nsamples*/,
        std::size_t                                       /*nalleles*/
    ) {
        outIdx = 0;
        nHomRef = nHet = nHomAlt = nMissing = 0;
        dosageSum = 0.0;
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
            dosageSum += dosage;
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

// Lookup table for the 8-bit probability decode: t[b] == double(b)/255.0,
// computed once.  Replacing the per-sample division with a table load is
// bit-for-bit identical (same operation, same rounding) and noticeably faster
// in the inner loop over hundreds of thousands of samples.
inline const double *prob8Table() {
    static const std::array<double, 256> table = [] {
        std::array<double, 256> a{};
        for (int i = 0; i < 256; ++i) a[i] = static_cast<double>(i) / 255.0;
        return a;
    }();
    return table.data();
}

// P3: specialised in-place parser for the most common BGEN v12 case — diploid,
// biallelic, unphased, 8- or 16-bit probabilities — bypassing genfile's
// per-sample / per-entry setter-callback machinery.  Reproduces genfile +
// DosageSetter bit-for-bit: v1 = P(g=0), v2 = P(g=1) read as byte/255 or
// u16/65535 (matching genfile's SpecialisedBitParser), dosage(ALT) =
// v2 + 2*max(1-v1-v2,0), or 2 - that under ref-last (countFirstAllele), then
// clamped to [0,2]; per-sample missingness from the ploidy byte's high bit; only
// used samples advance the output index (as DosageSetter::finalizeSample does).
// Returns false (caller falls back to genfile) if the probability region is
// shorter than expected.
inline bool fastParseV12DiploidBiallelic(
    const genfile::bgen::v12::GenotypeDataBlock &pack,
    DosageSetter &setter
) {
    const uint32_t N = pack.numberOfSamples;
    const uint32_t stride = (pack.bits == 8) ? 2u : 4u; // bytes/sample (2 probs)
    const genfile::byte_t *probs = pack.buffer;
    if (probs + static_cast<size_t>(N) * stride > pack.end)
        return false; // truncated / unexpected → let genfile handle it

    const genfile::byte_t *ploidy = pack.ploidy;
    const uint64_t *mask = setter.usedMask ? setter.usedMask->data() : nullptr;
    const bool allUsed = setter.allUsed;
    const bool countFirst = setter.countFirstAllele;
    double *const out = setter.out;
    std::vector<uint32_t> *const missingIdx = setter.missingIdx;
    const bool eightBit = (pack.bits == 8);
    const double *const t8 = eightBit ? prob8Table() : nullptr;

    uint32_t outIdx = 0;
    uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0, nMissing = 0;
    double dosageSum = 0.0;

    for (uint32_t i = 0; i < N; ++i) {
        const bool used = allUsed || ((mask[i >> 6] >> (i & 63)) & 1ULL);
        if (!used) continue; // non-used samples do not advance the output index
        if (ploidy[i] & 0x80) {
            out[outIdx] = std::numeric_limits<double>::quiet_NaN();
            if (missingIdx) missingIdx->push_back(outIdx);
            ++nMissing;
            ++outIdx;
            continue;
        }
        const genfile::byte_t *p = probs + static_cast<size_t>(i) * stride;
        double v1, v2;
        if (eightBit) {
            v1 = t8[p[0]];
            v2 = t8[p[1]];
        } else {
            const uint32_t u1 = static_cast<uint32_t>(p[0]) | (static_cast<uint32_t>(p[1]) << 8);
            const uint32_t u2 = static_cast<uint32_t>(p[2]) | (static_cast<uint32_t>(p[3]) << 8);
            v1 = static_cast<double>(u1) / 65535.0;
            v2 = static_cast<double>(u2) / 65535.0;
        }
        const double p2 = std::max(1.0 - v1 - v2, 0.0);
        double dosage = v2 + 2.0 * p2;          // dosage of the second allele
        if (countFirst) dosage = 2.0 - dosage;  // ref-last: count the first allele
        if (dosage < 0.0) dosage = 0.0;
        else if (dosage > 2.0) dosage = 2.0;
        out[outIdx] = dosage;
        dosageSum += dosage;
        if (dosage < 0.5) ++nHomRef;
        else if (dosage > 1.5) ++nHomAlt;
        else ++nHet;
        ++outIdx;
    }

    setter.nHomRef = nHomRef;
    setter.nHet = nHet;
    setter.nHomAlt = nHomAlt;
    setter.nMissing = nMissing;
    setter.dosageSum = dosageSum;
    return true;
}

// Skip the current variant's genotype data block by SEEKING past it, instead of
// streaming the bytes through the buffer as genfile's ignore_genotype_data_block
// does (std::istream::ignore reads/discards every byte).  Mirrors that function's
// framing exactly — for a compressed block the leading 4-byte payload size is
// consumed, then the payload is skipped; for the uncompressed layout-1 case the
// fixed 6*N bytes are skipped — so the stream lands at the identical position and
// the recorded marker metadata is byte-for-byte unchanged.  Used on the metadata
// first pass and the cursor's forward advance, where the genotype bytes are not
// needed; turns an O(file-size) read into O(num-variants) seeks.
inline void skipGenotypeDataBlock(std::istream &s, const genfile::bgen::Context &ctx) {
    if ((ctx.flags & genfile::bgen::e_CompressedSNPBlocks) != genfile::bgen::e_NoCompression) {
        unsigned char b[4];
        s.read(reinterpret_cast<char *>(b), 4);
        const uint32_t sz = static_cast<uint32_t>(b[0]) |
                            (static_cast<uint32_t>(b[1]) << 8) |
                            (static_cast<uint32_t>(b[2]) << 16) |
                            (static_cast<uint32_t>(b[3]) << 24);
        if (sz > 0) s.seekg(static_cast<std::streamoff>(sz), std::ios::cur);
    } else {
        s.seekg(static_cast<std::streamoff>(6) * ctx.number_of_samples, std::ios::cur);
    }
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════════════
// BgenData
// ══════════════════════════════════════════════════════════════════════════════

bool BgenData::loadFromIndex(
    const std::string &bgiPath,
    const std::unordered_set<std::string> &chrFilter,
    uint32_t &nSkippedMultiallelic
) {
    sqlite3 *db = nullptr;
    if (sqlite3_open_v2(bgiPath.c_str(), &db, SQLITE_OPEN_READONLY, nullptr) != SQLITE_OK) {
        if (db) sqlite3_close(db);
        return false;
    }

    // Staleness guard: bgenix's Metadata table records the indexed .bgen's
    // size; if it disagrees with the actual file the byte offsets cannot be
    // trusted, so fall back to a fresh scan.  Best-effort — a missing column or
    // table simply skips the check.
    {
        sqlite3_stmt *meta = nullptr;
        if (sqlite3_prepare_v2(db, "SELECT file_size FROM Metadata LIMIT 1", -1,
                               &meta, nullptr) == SQLITE_OK &&
            sqlite3_step(meta) == SQLITE_ROW) {
            const uint64_t indexed = static_cast<uint64_t>(sqlite3_column_int64(meta, 0));
            uint64_t actual = 0;
            if (indexed != 0 && statFileSize(m_bgenFile, &actual) && indexed != actual) {
                sqlite3_finalize(meta);
                sqlite3_close(db);
                warnMsg("BGEN: ignoring stale index %s (indexed file_size %llu != "
                        "actual %llu); falling back to scan",
                        bgiPath.c_str(),
                        static_cast<unsigned long long>(indexed),
                        static_cast<unsigned long long>(actual));
                return false;
            }
        }
        sqlite3_finalize(meta);
    }

    sqlite3_stmt *st = nullptr;
    const char *sql =
        "SELECT chromosome, position, rsid, number_of_alleles, allele1, allele2, "
        "file_start_position FROM Variant ORDER BY file_start_position";
    if (sqlite3_prepare_v2(db, sql, -1, &st, nullptr) != SQLITE_OK) {
        sqlite3_close(db);
        return false;
    }

    uint64_t ordinal = 0;
    int rc;
    while ((rc = sqlite3_step(st)) == SQLITE_ROW) {
        // Skip non-biallelic to match the scan path (the cursor decodes only
        // the first two alleles).
        if (sqlite3_column_int(st, 3) != 2) {
            ++nSkippedMultiallelic;
            continue;
        }
        const unsigned char *chrTxt = sqlite3_column_text(st, 0);
        const uint32_t position = static_cast<uint32_t>(sqlite3_column_int64(st, 1));
        const unsigned char *rsidTxt = sqlite3_column_text(st, 2);
        const unsigned char *a1 = sqlite3_column_text(st, 4);
        const unsigned char *a2 = sqlite3_column_text(st, 5);
        const uint64_t fpos = static_cast<uint64_t>(sqlite3_column_int64(st, 6));

        std::string chromosome = chrTxt ? reinterpret_cast<const char *>(chrTxt) : "";
        std::string rsid = rsidTxt ? reinterpret_cast<const char *>(rsidTxt) : "";
        std::string allele1 = a1 ? reinterpret_cast<const char *>(a1) : "";
        std::string allele2 = a2 ? reinterpret_cast<const char *>(a2) : "";
        // .bgi stores allele1 = first allele, allele2 = second allele; map to
        // REF/ALT per the user's --bgen orientation, identical to the scan path.
        const std::string &refAllele = m_altFirst ? allele2 : allele1;
        const std::string &altAllele = m_altFirst ? allele1 : allele2;

        m_chr.push_back(chromosome);
        m_pos.push_back(position);
        m_markerId.push_back(rsid);
        m_ref.push_back(refAllele);
        m_alt.push_back(altAllele);
        m_fileOffset.push_back(fpos);
        if (chrFilter.empty() || chrFilter.count(chromosome))
            m_markerInfo.push_back({chromosome, position, rsid, refAllele,
                                    altAllele, ordinal});
        ++ordinal;
    }
    const bool ok = (rc == SQLITE_DONE);
    sqlite3_finalize(st);
    sqlite3_close(db);

    if (!ok) {
        // Mid-iteration error: discard any partial metadata and fall back.
        m_chr.clear();
        m_pos.clear();
        m_markerId.clear();
        m_ref.clear();
        m_alt.clear();
        m_fileOffset.clear();
        m_markerInfo.clear();
        nSkippedMultiallelic = 0;
        return false;
    }
    return true;
}

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

    // ---- Marker metadata: prefer the .bgi (bgenix) index, else linear scan ----
    uint32_t nSkippedMultiallelic = 0;
    const std::string bgiPath = m_bgenFile + ".bgi";
    bool usedIndex = false;
    {
        uint64_t bgiSize = 0;
        if (statFileSize(bgiPath, &bgiSize) && bgiSize > 0)
            usedIndex = loadFromIndex(bgiPath, chrFilter, nSkippedMultiallelic);
    }

    if (usedIndex) {
        infoMsg("  BGEN: read marker metadata from index %s", bgiPath.c_str());
    } else {
        // Linear first pass: read each variant's identifying data and SEEK past
        // its genotype block, recording the byte offset of every biallelic
        // record so the cursor can later seek to it directly.
        std::string SNPID, RSID, chromosome;
        genfile::bgen::uint32_t position;
        std::vector<std::string> alleles;
        uint64_t variantIdx = 0;

        for (;;) {
            const std::streampos varStart = stream.tellg();
            if (!genfile::bgen::read_snp_identifying_data(
                    stream, context, &SNPID, &RSID, &chromosome, &position,
                    [&](std::size_t n) { alleles.resize(n); },
                    [&](std::size_t i, std::string const &a) { alleles[i] = a; }))
                break;

            // Skip non-biallelic: the BGEN cursor decodes only the first two
            // alleles, so retaining a multi-allelic record would silently
            // mis-report carriers of ALT2/ALT3/... .  Count the skips and
            // surface them in a warning after the first pass completes.
            if (alleles.size() != 2) {
                ++nSkippedMultiallelic;
                skipGenotypeDataBlock(stream, context);
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
            m_fileOffset.push_back(static_cast<uint64_t>(varStart));
            if (chrFilter.empty() || chrFilter.count(chromosome))
                m_markerInfo.push_back({chromosome, position,
                                        RSID.empty() ? SNPID : RSID, refAllele,
                                        altAllele, variantIdx});
            ++variantIdx;

            skipGenotypeDataBlock(stream, context);
        }
    }

    // ---- Apply --extract / --exclude filter ----
    geno_factory::filterMarkersByIds(m_markerInfo, extractFile, excludeFile);

    m_nMarkers = static_cast<uint32_t>(m_markerInfo.size());
    m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);

    infoMsg("  BGEN: read %u biallelic variants from %s", m_nMarkers, m_bgenFile.c_str());

    if (nSkippedMultiallelic > 0) {
        warnMsg("BGEN: skipped %u multi-allelic variant(s) in %s; "
                "the BGEN reader handles only biallelic records. "
                "Split them with plink2 --export bgen-1.x 'multiallelics-already-joined=split' "
                "(or an equivalent qctool pipeline) to retain them as separate biallelic records.",
                nSkippedMultiallelic, m_bgenFile.c_str());
    }
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

    // Scratch for the per-variant identifying-data read when seeking (reused
    // across calls to avoid per-marker allocations).
    std::string idSNPID, idRSID, idChr;
    genfile::bgen::uint32_t idPos = 0;
    std::vector<std::string> idAlleles;

    // Seek to the variant whose identifying-data block begins at byte `off`,
    // then read (and discard) the identifying data so the stream is left
    // positioned at the genotype block, ready for readDecompressParse().  The
    // byte offset comes from BgenData::fileOffset (.bgi index or recorded scan
    // position), so no forward scan is needed.
    void seekToVariant(uint64_t off) {
        stream.clear(); // drop any eof/fail set by a previous block
        stream.seekg(static_cast<std::streamoff>(off));
        if (!genfile::bgen::read_snp_identifying_data(
                stream, context, &idSNPID, &idRSID, &idChr, &idPos,
                [&](std::size_t n) { idAlleles.resize(n); },
                [&](std::size_t i, std::string const &a) { idAlleles[i] = a; }))
            throw std::runtime_error("BGEN: failed to read variant at offset " +
                                     std::to_string(off));
    }

    // Per-cursor (per-thread) libdeflate zlib decompressor; allocated lazily.
    libdeflate_decompressor *deflate = nullptr;

    // P2: decompress the layout-2 genotype block in buffer1 into buffer2 using
    // libdeflate (zlib) or the vendored zstd, or copy when uncompressed.  The
    // decompressed bytes are uniquely defined by the stream, so the result is
    // identical to genfile's zlib/zstd uncompress.  Returns false (leaving the
    // caller to use genfile's uncompress) for any codec we do not specialise.
    bool decompressLayout2() {
        const uint32_t comp = context.flags & genfile::bgen::e_CompressedSNPBlocks;
        if (comp == genfile::bgen::e_NoCompression) {
            buffer2.assign(buffer1.begin(), buffer1.end());
            return true;
        }
        if (buffer1.size() < 4) return false;
        const unsigned char *b = reinterpret_cast<const unsigned char *>(buffer1.data());
        const uint32_t usize = static_cast<uint32_t>(b[0]) |
                               (static_cast<uint32_t>(b[1]) << 8) |
                               (static_cast<uint32_t>(b[2]) << 16) |
                               (static_cast<uint32_t>(b[3]) << 24);
        buffer2.resize(usize);
        const void *in = buffer1.data() + 4;
        const size_t inLen = buffer1.size() - 4;
        if (comp == genfile::bgen::e_ZlibCompression) {
            if (!deflate) deflate = libdeflate_alloc_decompressor();
            if (!deflate || usize == 0) return false;
            size_t got = 0;
            const libdeflate_result r = libdeflate_zlib_decompress(
                deflate, in, inLen, buffer2.data(), usize, &got);
            return (r == LIBDEFLATE_SUCCESS && got == usize);
        }
        if (comp == genfile::bgen::e_ZstdCompression) {
            if (usize == 0) return false;
            const size_t got = ZSTD_decompress(buffer2.data(), usize, in, inLen);
            return (!ZSTD_isError(got) && got == usize);
        }
        return false;
    }

    // Read + decompress + parse the current variant's genotype block into the
    // setter.  Mirrors genfile's read_and_parse_genotype_data_block but routes
    // the common v12 case through libdeflate (P2) + the specialised parser (P3),
    // falling back to genfile's uncompress/parse for any other layout, codec,
    // ploidy, phasing, or bit depth — byte-identical in every case.
    void readDecompressParse(DosageSetter &setter) {
        genfile::bgen::read_genotype_data_block(stream, context, &buffer1);

        const uint32_t layout = context.flags & genfile::bgen::e_Layout;
        bool haveBuffer2 = false;
        if (layout == genfile::bgen::e_Layout2)
            haveBuffer2 = decompressLayout2(); // P2: libdeflate / zstd / copy
        if (!haveBuffer2)
            genfile::bgen::uncompress_probability_data(context, buffer1, &buffer2);

        bool parsed = false;
        if (layout == genfile::bgen::e_Layout2) {
            genfile::bgen::v12::GenotypeDataBlock pack(
                context, buffer2.data(), buffer2.data() + buffer2.size());
            if (pack.numberOfAlleles == 2 && pack.ploidyExtent[0] == 2 &&
                pack.ploidyExtent[1] == 2 && !pack.phased &&
                (pack.bits == 8 || pack.bits == 16))
                parsed = fastParseV12DiploidBiallelic(pack, setter); // P3
        }
        if (!parsed)
            genfile::bgen::parse_probability_data(
                buffer2.data(), buffer2.data() + buffer2.size(), context, setter);
    }

    ~Impl() {
        if (deflate) libdeflate_free_decompressor(deflate);
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
    (void)offset; // consumed only to position the stream for the header block
    genfile::bgen::read_header_block(impl.stream, &impl.context);

    // Skip sample identifiers
    if (impl.context.flags & genfile::bgen::e_SampleIdentifiers) {
        genfile::bgen::read_sample_identifier_block(impl.stream, impl.context, [](std::string const &) {
        });
    }
    // No seek to the first variant: the cursor seeks directly to each marker's
    // byte offset (BgenData::fileOffset) on every getGenotypes call.
}

BgenCursor::~BgenCursor() = default;

void BgenCursor::beginSequentialBlock(uint64_t /*firstMarker*/) {
    // Direct-seek cursor: nothing to prime per chunk.  Just clear any residual
    // eof/fail state so the next seekg succeeds; the genotype reads themselves
    // seek to the exact byte offset of each marker.
    m_impl->stream.clear();
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

    impl.seekToVariant(m_parent.fileOffset(gIndex));

    DosageSetter setter;
    setter.out = out.data();
    setter.missingIdx = &indexForMissing;
    setter.usedMask = impl.usedMask;
    setter.nSamplesInFile = impl.nSamplesInFile;
    setter.allUsed = impl.allUsed;
    setter.countFirstAllele = impl.altFirst;

    impl.readDecompressParse(setter);

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

    // Dosage-aware AF: compute the allele frequency from the dosage sum rather
    // than the 0.5/1.5-binned hard-call counts.  For hard calls the dosage sum
    // equals 2*nHomAlt + nHet exactly, so altFreq/altCounts/mac are unchanged;
    // for true dosages this is the correct (un-binned) frequency.  hweP and
    // missingRate retain the count-based values (HWE is not defined on dosages,
    // but hard calls reproduce the same hweP).
    const uint32_t nNonMissing = nUsed - setter.nMissing;
    if (nNonMissing > 0) {
        const double n2 = 2.0 * static_cast<double>(nNonMissing);
        altCounts = setter.dosageSum;
        altFreq = setter.dosageSum / n2;
        maf = std::min(altFreq, 1.0 - altFreq);
        mac = maf * n2;
    }
}

void BgenCursor::getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out
) {
    auto &impl = *m_impl;

    impl.seekToVariant(m_parent.fileOffset(gIndex));

    DosageSetter setter;
    setter.out = out.data();
    setter.missingIdx = nullptr;
    setter.usedMask = impl.usedMask;
    setter.nSamplesInFile = impl.nSamplesInFile;
    setter.allUsed = impl.allUsed;
    setter.countFirstAllele = impl.altFirst;

    impl.readDecompressParse(setter);
}
