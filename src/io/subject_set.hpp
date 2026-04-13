// subject_set.hpp — Lightweight genotype→GRM→keep→remove subject filter
//
// Handles the subject intersection pipeline that is common to all analysis
// modes: intersect .fam IIDs with GRM subjects, apply --keep/--remove
// filters, and build a 64-bit bitmask of used subjects.
//
// SubjectData owns a SubjectSet (composition) and adds the phenotype/
// residual intersection on top.  Utility modes like runPairwiseIBD and
// runPhiEstimation use SubjectSet directly (they have no pheno/resid).
#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

class SubjectSet {
  public:
    // Construct from .fam IIDs (genotype subject list).
    explicit SubjectSet(std::vector<std::string> famIIDs);

    // ── Pre-finalize configuration ─────────────────────────────────────

    // Set GRM subject IDs for intersection.  Empty set → skip GRM step.
    void setGrmSubjects(std::unordered_set<std::string> grmIDs);

    // Set --keep / --remove filter files.  Empty path → skip that filter.
    void setKeepRemove(
        const std::string &keepFile,
        const std::string &removeFile
    );

    // Set descriptive labels for pipeline logging (e.g., "--bfile d_bed").
    // If not set, generic descriptions are used.
    void setGenoLabel(const std::string &label) {
        m_genoLabel = label;
    }

    void setGrmLabel(const std::string &label) {
        m_grmLabel = label;
    }

    // ── Build bitmask ──────────────────────────────────────────────────
    // Applies the pipeline: genotype → GRM → keep → remove.
    // Prints the "Subject pipeline" log table.
    // Must be called exactly once.
    void finalize();

    // ── Accessors (valid after finalize) ───────────────────────────────

    uint32_t nFam() const {
        return m_nFam;
    }

    uint32_t nUsed() const {
        return m_nUsed;
    }

    const std::vector<std::string> &famIIDs() const {
        return m_famIIDs;
    }

    const std::vector<uint64_t> &usedMask() const {
        return m_usedMask;
    }

    uint32_t nMaskWords() const {
        return static_cast<uint32_t>(m_usedMask.size());
    }

    // Dense index → .fam index mapping.
    const std::vector<uint32_t> &usedFamIndices() const {
        return m_usedFamIndices;
    }

    // IIDs of used subjects in .fam order.
    std::vector<std::string> usedIIDs() const;

    // ── Post-finalize mask manipulation ────────────────────────────────
    // Remove a subject from the used set (clear its bit, rebuild indices).
    // Used by SubjectData when phenotype/residual filtering narrows the set.
    void narrowMask(
        const std::vector<uint64_t> &newMask,
        uint32_t newNUsed
    );

  private:
    uint32_t m_nFam = 0;
    uint32_t m_nUsed = 0;
    bool m_finalized = false;

    std::vector<std::string> m_famIIDs;
    std::vector<uint64_t> m_usedMask;
    std::vector<uint32_t> m_usedFamIndices;

    // Pre-finalize state
    std::string m_keepFile;
    std::string m_removeFile;
    std::unordered_set<std::string> m_grmSubjects;
    std::string m_genoLabel; // e.g., "--bfile d_bed"
    std::string m_grmLabel;  // e.g., "--sp-grm-grab e_grm.grab"
};
