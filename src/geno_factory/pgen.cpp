// pgen.cpp — PgenData and PgenCursor implementations
//
// Wraps pgenlib (plink2 pgen reader) behind the GenoMeta / GenoCursor
// interfaces.  Each PgenCursor owns an independent PgenReader (thread-safe).

#include "geno_factory/pgen.hpp"
#include "geno_factory/hwe.hpp"
#include "util/logging.hpp"

#include "util/text_scanner.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>

// ── pgenlib headers (C++ mode, plink2 namespace) ─────────────────────────────
#include "pgenlib_ffi_support.h"
#include "pgenlib_read.h"

namespace {

// ──────────────────────────────────────────────────────────────────────────────
// .pvar parser — tab-delimited with header line starting with #CHROM
// Columns: #CHROM  POS  ID  REF  ALT  [...]
// ──────────────────────────────────────────────────────────────────────────────

struct PvarRecord {
    std::string chrom;
    uint32_t pos;
    std::string id;
    std::string ref;
    std::string alt;
};

std::vector<PvarRecord> parsePvarFile(const std::string &path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open .pvar file: " + path);

    std::vector<PvarRecord> records;
    std::string line;
    bool headerSeen = false;
    int colChrom = -1, colPos = -1, colId = -1, colRef = -1, colAlt = -1;

    while (std::getline(in, line)) {
        if (line.empty()) continue;
        // Skip comment lines
        if (line.size() >= 2 && line[0] == '#' && line[1] == '#') continue;

        // Parse header
        if (!headerSeen && line[0] == '#') {
            headerSeen = true;
            text::TokenScanner ts(line);
            int col = 0;
            while (!ts.atEnd()) {
                auto sv = ts.nextView();
                if (sv.empty()) break;
                if (sv == "#CHROM" || sv == "CHROM")
                    colChrom = col;
                else if (sv == "POS")
                    colPos = col;
                else if (sv == "ID")
                    colId = col;
                else if (sv == "REF")
                    colRef = col;
                else if (sv == "ALT")
                    colAlt = col;
                ++col;
            }
            if (colChrom < 0 || colPos < 0 || colId < 0 || colRef < 0 || colAlt < 0)
                throw std::runtime_error(".pvar header missing required columns (#CHROM POS ID REF ALT)");
            continue;
        }

        // If no header was found, assume default column order
        if (!headerSeen) {
            headerSeen = true;
            colChrom = 0;
            colPos = 1;
            colId = 2;
            colRef = 3;
            colAlt = 4;
            // Fall through to parse this line as data
        }

        // Parse data line with zero-copy token scanning
        int maxCol = std::max({colChrom, colPos, colId, colRef, colAlt});
        text::TokenScanner ts(line);
        std::vector<std::string_view> tokViews;
        tokViews.reserve(maxCol + 2);
        while (!ts.atEnd()) {
            auto sv = ts.nextView();
            if (sv.empty()) break;
            tokViews.push_back(sv);
            if (static_cast<int>(tokViews.size()) > maxCol) break; // enough columns
        }
        if (static_cast<int>(tokViews.size()) <= maxCol)
            throw std::runtime_error(".pvar line has too few columns: " + line);

        PvarRecord r;
        r.chrom = std::string(tokViews[colChrom]);
        r.pos = static_cast<uint32_t>(std::stoul(std::string(tokViews[colPos])));
        r.id = std::string(tokViews[colId]);
        r.ref = std::string(tokViews[colRef]);
        r.alt = std::string(tokViews[colAlt]);
        // For multiallelic with comma-separated ALT, take only the first
        auto comma = r.alt.find(',');
        if (comma != std::string::npos) r.alt.resize(comma);
        records.push_back(std::move(r));
    }

    if (records.empty()) throw std::runtime_error(".pvar file has no data records: " + path);

    return records;
}

} // anonymous namespace

// Allocate zero-initialized memory with 32-byte alignment (required by pgenlib AVX2).
static uintptr_t* allocAligned(size_t nWords) {
    if (nWords == 0) return nullptr;
    const size_t nBytes = nWords * sizeof(uintptr_t);
    void* ptr = nullptr;
    if (posix_memalign(&ptr, 32, nBytes) != 0)
        throw std::bad_alloc();
    std::memset(ptr, 0, nBytes);
    return static_cast<uintptr_t*>(ptr);
}

// Free a buffer allocated by plink2::cachealigned_malloc / aligned_malloc.
// These store the original malloc pointer at aligned_ptr[-1].
static void freeCacheAligned(void* p) {
    if (p) plink2::aligned_free(p);
}

// ══════════════════════════════════════════════════════════════════════════════
// PgenData
// ══════════════════════════════════════════════════════════════════════════════

void PgenData::PgfiDeleter::operator()(PgenFileInfo *p) const {
    if (p) {
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgfi(p, &reterr);
        delete p;
    }
}

PgenData::PgenData(
    std::string pgenFile,
    std::string pvarFile,
    const std::vector<uint64_t> &usedMask,
    uint32_t nSamplesInFile,
    uint32_t nUsed,
    std::unordered_set<std::string> chrFilter,
    int nMarkersEachChunk
)
    : m_pgenFile(std::move(pgenFile)),
      m_nSubjInFile(nSamplesInFile),
      m_nSubjUsed(nUsed),
      m_usedMask(usedMask),
      m_pgfiAlloc(nullptr, freeCacheAligned)
{
    m_allUsed = (nUsed == nSamplesInFile);

    // ---- Parse .pvar ----
    auto pvarRecords = parsePvarFile(pvarFile);
    m_nMarkers = static_cast<uint32_t>(pvarRecords.size());
    m_chr.reserve(m_nMarkers);
    m_pos.reserve(m_nMarkers);
    m_markerId.reserve(m_nMarkers);
    m_ref.reserve(m_nMarkers);
    m_alt.reserve(m_nMarkers);
    m_markerInfo.reserve(m_nMarkers);

    for (uint32_t i = 0; i < m_nMarkers; ++i) {
        auto &r = pvarRecords[i];
        m_chr.push_back(r.chrom);
        m_pos.push_back(r.pos);
        m_markerId.push_back(r.id);
        m_ref.push_back(r.ref);
        m_alt.push_back(r.alt);
        // For pgen: REF = pvar REF, ALT = pvar ALT.
        // Genotype coding: 0 = hom_ref (0 alt), 1 = het (1 alt), 2 = hom_alt (2 alt).
        // This matches our convention: always count ALT.
        m_markerInfo.push_back({r.chrom, r.pos, r.id, r.ref, r.alt, i});
    }

    // ---- Apply --chr filter ----
    if (!chrFilter.empty()) {
        std::vector<MarkerInfo> filtered;
        filtered.reserve(m_markerInfo.size());
        for (const auto &m : m_markerInfo)
            if (chrFilter.count(m.chrom)) filtered.push_back(m);
        m_markerInfo = std::move(filtered);
    }

    m_chunkIndices = buildChunks(m_markerInfo, nMarkersEachChunk);

    // ---- Initialize pgenlib (two-phase) ----
    using namespace plink2;

    auto *pgfi_raw = new PgenFileInfo{};
    PreinitPgfi(pgfi_raw);
    m_pgfi.reset(pgfi_raw);

    PgenHeaderCtrl header_ctrl;
    uintptr_t pgfi_alloc_cacheline_ct;
    char errstr_buf[kPglErrstrBufBlen];

    PglErr err = PgfiInitPhase1(m_pgenFile.c_str(),
                                nullptr,       // pgi_fname auto-detect
                                m_nMarkers,    // raw_variant_ct
                                m_nSubjInFile, // raw_sample_ct
                                &header_ctrl, pgfi_raw, &pgfi_alloc_cacheline_ct, errstr_buf);
    if (err != kPglRetSuccess) throw std::runtime_error(std::string("pgen init phase 1 failed: ") + errstr_buf);

    // Verify counts match
    if (pgfi_raw->raw_sample_ct != m_nSubjInFile)
        throw std::runtime_error("pgen sample count (" + std::to_string(pgfi_raw->raw_sample_ct) +
                                 ") does not match expected (" + std::to_string(m_nSubjInFile) + ")");
    if (pgfi_raw->raw_variant_ct != m_nMarkers)
        throw std::runtime_error("pgen variant count (" + std::to_string(pgfi_raw->raw_variant_ct) +
                                 ") does not match .pvar (" + std::to_string(m_nMarkers) + ")");

    // Allocate pgfi memory
    unsigned char *pgfi_alloc_raw = nullptr;
    if (pgfi_alloc_cacheline_ct > 0) {
        if (cachealigned_malloc(pgfi_alloc_cacheline_ct * kCacheline, &pgfi_alloc_raw))
            throw std::runtime_error("pgen: pgfi memory allocation failed");
    }
    m_pgfiAlloc.reset(pgfi_alloc_raw);

    // Phase 2
    uintptr_t pgr_alloc_cacheline_ct;
    err = PgfiInitPhase2(header_ctrl,
                         0, // allele_cts_already_loaded
                         0, // nonref_flags_already_loaded
                         0, // use_blockload (mode 2 = per-variant fread)
                         0, // vblock_idx_start
                         m_nMarkers, &m_maxVrecWidth, pgfi_raw, pgfi_alloc_raw, &pgr_alloc_cacheline_ct, errstr_buf);
    if (err != kPglRetSuccess) throw std::runtime_error(std::string("pgen init phase 2 failed: ") + errstr_buf);

    m_pgrAllocCachelineCt = pgr_alloc_cacheline_ct;

    // Close the shared file handle so that each PgenCursor (PgrInit) opens
    // its own independent FILE*.  Without this, the very first PgrInit call
    // steals shared_ff, and concurrent PgrInit calls from worker threads
    // can race on the same pointer (data race → double-fclose / corrupt reads).
    if (pgfi_raw->shared_ff) {
        fclose(pgfi_raw->shared_ff);
        pgfi_raw->shared_ff = nullptr;
    }
}

PgenData::~PgenData() = default;

std::vector<std::vector<uint64_t> > PgenData::buildChunks(
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

std::unique_ptr<GenoCursor> PgenData::makeCursor() const {
    return std::make_unique<PgenCursor>(*this);
}

// ══════════════════════════════════════════════════════════════════════════════
// PgenCursor::Impl — per-thread pgenlib reader state
// ══════════════════════════════════════════════════════════════════════════════

struct PgenCursor::Impl {
    plink2::PgenReader pgr;
    std::unique_ptr<unsigned char, void (*)(void *)> pgrAlloc;
    plink2::PgrSampleSubsetIndex pssi;

    // Genotype scratch buffers — must be 32-byte aligned for pgenlib AVX2 SIMD.
    // Sizes are rounded up to VecW (kWordsPerVec) boundary so SIMD loops
    // that process in 32-byte chunks never overrun.
    uintptr_t* genovec = nullptr;         // packed 2-bit genotypes
    uintptr_t* dosagePresent = nullptr;   // dosage presence bitarray
    std::vector<uint16_t> dosageMain;     // raw dosage values (no SIMD alignment needed)

    // Difflist scratch buffers for sparse read
    uintptr_t* diffRaregeno = nullptr;    // 2-bit packed rare genotype codes
    std::vector<uint32_t> diffSampleIdsRaw;

    // Sample subsetting support
    uintptr_t* sampleInclude = nullptr;   // pgenlib reads with SIMD
    std::vector<uint32_t> cumulativePopcounts;

    uint32_t sampleCt;
    uint32_t rawSampleCt;
    bool allUsed;

    Impl()
        : pgrAlloc(nullptr, freeCacheAligned)
    {
    }

    ~Impl()
    {
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgr(&pgr, &reterr);
        free(genovec);
        free(dosagePresent);
        free(diffRaregeno);
        free(sampleInclude);
    }

};

PgenCursor::PgenCursor(const PgenData &parent)
    : m_parent(parent),
      m_impl(std::make_unique<Impl>())
{
    using namespace plink2;

    auto &impl = *m_impl;
    impl.rawSampleCt = parent.nSubjInFile();
    impl.sampleCt = parent.nSubjUsed();
    impl.allUsed = parent.allUsed();

    // Allocate genovec buffer: rounded up to VecW boundary for AVX2 safety
    const uint32_t rawSampleCtl2 = plink2::NypCtToVecCt(impl.rawSampleCt) * plink2::kWordsPerVec;
    impl.genovec = allocAligned(rawSampleCtl2);

    // Dosage buffers (sized for raw_sample_ct)
    const uint32_t rawSampleCtaw = plink2::BitCtToAlignedWordCt(impl.rawSampleCt);
    impl.dosagePresent = allocAligned(rawSampleCtaw);
    impl.dosageMain.resize(impl.rawSampleCt, 0);

    // Difflist buffers for sparse read — round up to VecW boundary
    const uint32_t maxDifflistLen = impl.sampleCt / 8;
    const uint32_t diffRaregenoWords = plink2::NypCtToVecCt(maxDifflistLen + 1) * plink2::kWordsPerVec;
    impl.diffRaregeno = allocAligned(diffRaregenoWords);
    impl.diffSampleIdsRaw.resize(maxDifflistLen + 2, 0); // +2 for pgenlib sentinel

    // Sample include bitmask + cumulative popcounts
    const uint32_t rawSampleCtl = BitCtToWordCt(impl.rawSampleCt);
    if (!impl.allUsed) {
        const uint32_t sampleIncludeWords = plink2::BitCtToVecCt(impl.rawSampleCt) * plink2::kWordsPerVec;
        impl.sampleInclude = allocAligned(sampleIncludeWords);
        // Convert our uint64_t usedMask to pgenlib's uintptr_t bitarray
        const auto &mask = parent.usedMask();
        static_assert(sizeof(uintptr_t) == sizeof(uint64_t), "PgenCursor assumes LP64 (64-bit uintptr_t)");
        std::memcpy(impl.sampleInclude, mask.data(),
                    std::min(static_cast<size_t>(rawSampleCtl), mask.size()) * sizeof(uintptr_t));

        impl.cumulativePopcounts.resize(rawSampleCtl, 0);
        FillCumulativePopcounts(impl.sampleInclude, rawSampleCtl, impl.cumulativePopcounts.data());
    }

    // Initialize PgenReader
    PreinitPgr(&impl.pgr);
    unsigned char *pgr_alloc_raw = nullptr;
    if (parent.pgrAllocCachelineCt() > 0) {
        if (cachealigned_malloc(parent.pgrAllocCachelineCt() * kCacheline, &pgr_alloc_raw))
            throw std::runtime_error("pgen: per-thread pgr allocation failed");
    }
    impl.pgrAlloc.reset(pgr_alloc_raw);

    PglErr err = PgrInit(parent.pgenFile().c_str(), parent.maxVrecWidth(), parent.pgfi(), &impl.pgr, pgr_alloc_raw);
    if (err != kPglRetSuccess)
        throw std::runtime_error("pgen: PgrInit failed (error " + std::to_string(static_cast<int>(err)) + ")");

    // Set sample subset
    if (impl.allUsed)
        PgrClearSampleSubsetIndex(&impl.pgr, &impl.pssi);
    else
        PgrSetSampleSubsetIndex(impl.cumulativePopcounts.data(), &impl.pgr, &impl.pssi);
}

PgenCursor::~PgenCursor() = default;

void PgenCursor::beginSequentialBlock(uint64_t /*firstMarker*/) {
    // pgenlib handles seeking internally; nothing to do.
}

void PgenCursor::getGenotypes(
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
    using namespace plink2;
    auto &impl = *m_impl;
    const uint32_t sampleCt = impl.sampleCt;
    indexForMissing.clear();

    // Read genotypes with dosage
    uint32_t dosage_ct = 0;
    PglErr err =
        PgrGetD(impl.allUsed ? nullptr : impl.sampleInclude, impl.pssi, sampleCt, static_cast<uint32_t>(gIndex),
                &impl.pgr, impl.genovec, impl.dosagePresent, impl.dosageMain.data(), &dosage_ct);

    if (err != kPglRetSuccess) throw std::runtime_error("pgen: PgrGetD failed at variant " + std::to_string(gIndex));

    // Decode genotypes to doubles and collect counts
    // pgen 2-bit codes: 0=hom_ref, 1=het, 2=hom_alt, 3=missing
    // We want ALT allele dosage: 0, 1, 2, NaN
    uint32_t nHomRef = 0, nHet = 0, nHomAlt = 0, nMissing = 0;
    const uintptr_t *gv = impl.genovec;

    for (uint32_t i = 0; i < sampleCt; ++i) {
        const uint32_t code = (gv[i / kBitsPerWordD2] >> (2 * (i % kBitsPerWordD2))) & 3;
        switch (code) {
        case 0:
            out[i] = 0.0;
            ++nHomRef;
            break;
        case 1:
            out[i] = 1.0;
            ++nHet;
            break;
        case 2:
            out[i] = 2.0;
            ++nHomAlt;
            break;
        default:
            out[i] = std::numeric_limits<double>::quiet_NaN();
            indexForMissing.push_back(i);
            ++nMissing;
            break;
        }
    }

    // Overlay dosage values where present
    if (dosage_ct > 0) {
        const uintptr_t *dp = impl.dosagePresent;
        const uint16_t *dm = impl.dosageMain.data();
        uint32_t dIdx = 0;
        for (uint32_t i = 0; i < sampleCt && dIdx < dosage_ct; ++i) {
            const uint32_t word_idx = i / kBitsPerWord;
            const uint32_t bit_idx = i % kBitsPerWord;
            if ((dp[word_idx] >> bit_idx) & 1) {
                // dosage_main values: 0..32768 maps to 0.0..2.0
                out[i] = static_cast<double>(dm[dIdx]) / 16384.0;
                ++dIdx;
            }
        }
    }

    // Compute QC stats
    GenoStats gs = statsFromCounts(nHomRef, nHet, nHomAlt, nMissing, sampleCt);
    altFreq = gs.altFreq;
    altCounts = gs.altCounts;
    missingRate = gs.missingRate;
    hweP = gs.hweP;
    maf = gs.maf;
    mac = gs.mac;
}

void PgenCursor::getGenotypesSimple(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out
) {
    using namespace plink2;
    auto &impl = *m_impl;
    const uint32_t sampleCt = impl.sampleCt;

    PglErr err = PgrGet(impl.allUsed ? nullptr : impl.sampleInclude, impl.pssi, sampleCt,
                        static_cast<uint32_t>(gIndex), &impl.pgr, impl.genovec);

    if (err != kPglRetSuccess) throw std::runtime_error("pgen: PgrGet failed at variant " + std::to_string(gIndex));

    const uintptr_t *gv = impl.genovec;
    for (uint32_t i = 0; i < sampleCt; ++i) {
        const uint32_t code = (gv[i / kBitsPerWordD2] >> (2 * (i % kBitsPerWordD2))) & 3;
        switch (code) {
        case 0:
            out[i] = 0.0;
            break;
        case 1:
            out[i] = 1.0;
            break;
        case 2:
            out[i] = 2.0;
            break;
        default:
            out[i] = std::numeric_limits<double>::quiet_NaN();
            break;
        }
    }
}

uint32_t PgenCursor::getGenotypesMaybeSparse(
    uint64_t gIndex,
    Eigen::Ref<Eigen::VectorXd> out,
    uint32_t maxLen,
    uint32_t *diffSampleIds,
    uint8_t *diffGenoCodes,
    uint32_t &diffLen
) {
    using namespace plink2;
    auto &impl = *m_impl;
    const uint32_t sampleCt = impl.sampleCt;

    // Cap maxLen at our pre-allocated buffer size (sampleCt / 8)
    const uint32_t bufMaxLen = sampleCt / 8;
    if (maxLen > bufMaxLen) maxLen = bufMaxLen;

    uint32_t commonGeno = UINT32_MAX;
    diffLen = 0;

    PglErr err = PgrGetDifflistOrGenovec(
        impl.allUsed ? nullptr : impl.sampleInclude,
        impl.pssi, sampleCt, maxLen,
        static_cast<uint32_t>(gIndex), &impl.pgr,
        impl.genovec, &commonGeno,
        impl.diffRaregeno, impl.diffSampleIdsRaw.data(), &diffLen);

    if (err != kPglRetSuccess)
        throw std::runtime_error("pgen: PgrGetDifflistOrGenovec failed at variant " +
                                 std::to_string(gIndex));

    if (commonGeno == UINT32_MAX) {
        // Dense: expand genovec to doubles
        const uintptr_t *gv = impl.genovec;
        for (uint32_t i = 0; i < sampleCt; ++i) {
            const uint32_t code = (gv[i / kBitsPerWordD2] >> (2 * (i % kBitsPerWordD2))) & 3;
            out[i] = (code < 3) ? static_cast<double>(code)
                                : std::numeric_limits<double>::quiet_NaN();
        }
        return UINT32_MAX;
    }

    // Sparse: unpack 2-bit raregeno to uint8 + copy sample IDs
    const uintptr_t *rg = impl.diffRaregeno;
    for (uint32_t j = 0; j < diffLen; ++j) {
        diffGenoCodes[j] = static_cast<uint8_t>(
            (rg[j / kBitsPerWordD2] >> (2 * (j % kBitsPerWordD2))) & 3);
        diffSampleIds[j] = impl.diffSampleIdsRaw[j];
    }

    return commonGeno;
}
