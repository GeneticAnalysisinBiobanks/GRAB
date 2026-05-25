// subject_set.cpp — Lightweight genotype→GRM→keep→remove subject filter

#include "io/subject_set.hpp"
#include "io/subject_filter.hpp"
#include "util/logging.hpp"

#include <cstdio>
#include <stdexcept>

SubjectSet::SubjectSet(std::vector<std::string> famIIDs)
    : m_nFam(static_cast<uint32_t>(famIIDs.size())),
      m_famIIDs(std::move(famIIDs))
{
}

void SubjectSet::setGrmSubjects(std::unordered_set<std::string> grmIDs) {
    if (m_finalized) throw std::runtime_error("SubjectSet::setGrmSubjects called after finalize");
    m_grmSubjects = std::move(grmIDs);
}

void SubjectSet::setKeepRemove(
    const std::string &keepFile,
    const std::string &removeFile
) {
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

    // ── Log subject pipeline ──────────────────────────────────────────
    const uint32_t nAfterGrm = hasGrm ? nGrm : nGeno;
    const uint32_t nAfterKeep = hasKeep ? nKeep : nAfterGrm;
    const uint32_t nAfterRemove = hasRemove ? nRemove : nAfterKeep;
    (void)nAfterRemove;

    const std::string genoLabel = m_genoLabel.empty() ? "genotype file" : m_genoLabel;
    const std::string grmLabel = m_grmLabel.empty() ? "sparse GRM" : m_grmLabel;

    std::fprintf(stderr, "\n");
    std::fprintf(stderr, "  %u subjects in %s\n", nGeno, genoLabel.c_str());
    if (hasGrm)
        std::fprintf(stderr, "  %u subjects in %s\n", nGrm, grmLabel.c_str());
    else
        std::fprintf(stderr, "  (no GRM provided, skipping GRM intersection)\n");
    if (hasKeep)
        std::fprintf(stderr, "  %u subjects in --keep %s\n", nAfterKeep, m_keepFile.c_str());
    else
        std::fprintf(stderr, "  (no --keep provided)\n");
    if (hasRemove)
        std::fprintf(stderr, "  %u subjects after --remove %s\n", nAfterRemove, m_removeFile.c_str());
    else
        std::fprintf(stderr, "  (no --remove provided)\n");
    std::fprintf(stderr, "  %u subjects used in analysis\n\n", m_nUsed);

    m_finalized = true;
}

std::vector<std::string> SubjectSet::usedIIDs() const {
    std::vector<std::string> result;
    result.reserve(m_nUsed);
    for (uint32_t idx : m_usedFamIndices)
        result.push_back(m_famIIDs[idx]);
    return result;
}

void SubjectSet::narrowMask(
    const std::vector<uint64_t> &newMask,
    uint32_t newNUsed
) {
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
