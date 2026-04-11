// subject_set.cpp — Lightweight genotype→GRM→keep→remove subject filter

#include "subj_reader/subject_set.hpp"
#include "subj_reader/subject_filter.hpp"
#include "util/logging.hpp"

#include <cstdio>
#include <stdexcept>

SubjectSet::SubjectSet(std::vector<std::string> famIIDs)
    : m_nFam(static_cast<uint32_t>(famIIDs.size())), m_famIIDs(std::move(famIIDs)) {}

void SubjectSet::setGrmSubjects(std::unordered_set<std::string> grmIDs) {
    if (m_finalized) throw std::runtime_error("SubjectSet::setGrmSubjects called after finalize");
    m_grmSubjects = std::move(grmIDs);
}

void SubjectSet::setKeepRemove(const std::string &keepFile, const std::string &removeFile) {
    if (m_finalized) throw std::runtime_error("SubjectSet::setKeepRemove called after finalize");
    m_keepFile = keepFile;
    m_removeFile = removeFile;
}

void SubjectSet::finalize() {
    if (m_finalized) throw std::runtime_error("SubjectSet::finalize already called");

    // ── Load --keep / --remove subject ID sets ──────────────────────────
    const bool hasKeep = !m_keepFile.empty();
    const bool hasRemove = !m_removeFile.empty();
    std::unordered_set<std::string> keepSet, removeSet;
    if (hasKeep) {
        keepSet = parseSubjectIDFile(m_keepFile);
        infoMsg("--keep: %zu subjects listed in %s", keepSet.size(), m_keepFile.c_str());
    }
    if (hasRemove) {
        removeSet = parseSubjectIDFile(m_removeFile);
        infoMsg("--remove: %zu subjects listed in %s", removeSet.size(), m_removeFile.c_str());
    }

    // ── Intersect with .fam, build bitmask ──────────────────────────────
    const uint32_t nWords = (m_nFam + 63) / 64;
    m_usedMask.assign(nWords, 0ULL);

    // Pipeline counters
    const uint32_t nGeno = m_nFam;
    uint32_t nGrm = 0;
    uint32_t nKeep = 0;
    uint32_t nRemove = 0;
    const bool hasGrm = !m_grmSubjects.empty();

    m_usedFamIndices.clear();
    m_usedFamIndices.reserve(m_nFam);

    for (uint32_t f = 0; f < m_nFam; ++f) {
        const auto &iid = m_famIIDs[f];

        // ── GRM intersection ───────────────────────────────────────────
        if (hasGrm && !m_grmSubjects.count(iid)) continue;
        ++nGrm;

        // ── --keep filter ──────────────────────────────────────────────
        if (hasKeep && !keepSet.count(iid)) continue;
        ++nKeep;

        // ── --remove filter ────────────────────────────────────────────
        if (hasRemove && removeSet.count(iid)) continue;
        ++nRemove;

        // Subject passes all checks
        m_usedMask[f / 64] |= 1ULL << (f % 64);
        m_usedFamIndices.push_back(f);
    }
    m_nUsed = static_cast<uint32_t>(m_usedFamIndices.size());

    if (m_nUsed == 0) throw std::runtime_error("SubjectSet: no subjects remain after filtering");

    // ── Log subject pipeline table ──────────────────────────────────────
    const uint32_t nAfterGrm = hasGrm ? nGrm : nGeno;
    const uint32_t nAfterKeep = hasKeep ? nKeep : nAfterGrm;
    const uint32_t nAfterRemove = hasRemove ? nRemove : nAfterKeep;
    (void)nAfterRemove; // last step count is m_nUsed

    std::fprintf(
        stderr,
        "\n\xe2\x94\x80\xe2\x94\x80 Subject pipeline "
        "\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80"
        "\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80"
        "\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80"
        "\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\n"
        "  Step              Count  Description\n"
        "  genotype         %6u  .fam / .psam / VCF / BGEN subjects\n",
        nGeno);
    if (hasGrm)
        std::fprintf(stderr, "  \xe2\x88\xa9 GRM            %6u  sparse GRM subjects\n", nGrm);
    else
        std::fprintf(stderr, "  \xe2\x88\xa9 GRM                 -  (no GRM provided)\n");
    std::fprintf(stderr,
                 "  \xe2\x88\xa9 keep           %6u  %s\n"
                 "  \\ remove         %6u  %s\n",
                 nAfterKeep, hasKeep ? "--keep filter" : "(no --keep provided)", nAfterRemove,
                 hasRemove ? "--remove filter" : "(no --remove provided)");
    std::fprintf(stderr,
                 "\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2"
                 "\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94"
                 "\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80"
                 "\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2"
                 "\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94"
                 "\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80"
                 "\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\xe2\x94\x80\n\n");

    m_finalized = true;
}

std::vector<std::string> SubjectSet::usedIIDs() const {
    std::vector<std::string> result;
    result.reserve(m_nUsed);
    for (uint32_t idx : m_usedFamIndices)
        result.push_back(m_famIIDs[idx]);
    return result;
}

void SubjectSet::narrowMask(const std::vector<uint64_t> &newMask, uint32_t newNUsed) {
    m_usedMask = newMask;
    m_nUsed = newNUsed;
    m_usedFamIndices.clear();
    m_usedFamIndices.reserve(m_nUsed);
    for (uint32_t w = 0; w < static_cast<uint32_t>(m_usedMask.size()); ++w) {
        uint64_t bits = m_usedMask[w];
        while (bits) {
            int b = __builtin_ctzll(bits);
            m_usedFamIndices.push_back(w * 64 + static_cast<uint32_t>(b));
            bits &= bits - 1;
        }
    }
}
