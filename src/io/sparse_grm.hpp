// sparse_grm.hpp — Sparse GRM reader and quadratic-form evaluator
//
// Supports two formats:
//   GRAB: 3-column TSV  (ID1 \t ID2 \t value), '#'-header skipped
//   GCTA:  plink2 --make-grm-sparse .grm.sp file + optional .grm.id
//          .grm.id has FID IID per line; .grm.sp has 0-based-idx1 idx2 value
//
// Dedicated method for the quadratic form  x^T G x  where G is the sparse GRM
// restricted to a given set of subjects.  This is needed by WtCoxG (variance
// ratios) and will be reused by SPAGRM, LEAF, etc.
#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class SparseGRM {
  public:
    // Storage: COO (coordinate) format — lower triangle + diagonal only,
    // kept sorted by (row, col) for cache locality when iterating.
    // quadForm / bilinearForm handle the 2× symmetry factor internally.
    struct Entry {
        uint32_t row; // index into subjectOrder
        uint32_t col; // index into subjectOrder  (col <= row)
        double value;
    };

    // ── GRAB format ────────────────────────────────────────────────────
    // Parse the 3-column file.  Subjects not in `subjectOrder` are silently
    // dropped.  `subjectOrder` defines the canonical index mapping so that
    // quadratic-form vectors can use the same ordering.
    // Stores lower triangle + diagonal only (half-storage).
    SparseGRM(
        const std::string &filename,
        const std::vector<std::string> &subjectOrder
    );

    // ── GCTA format (plink2 --make-grm-sparse) ────────────────────────
    // Reads `spFile` (the .grm.sp file) and derives the .grm.id path from it.
    // Subjects not in `subjectOrder` are silently dropped.
    // If .grm.id does not exist and `famIIDs` is provided, the 0-based
    // indices in .grm.sp are assumed to correspond to `famIIDs` order.
    // Stores lower triangle + diagonal only (half-storage).
    static SparseGRM fromGCTA(
        const std::string &spFile,
        const std::vector<std::string> &subjectOrder,
        const std::vector<std::string> &famIIDs = {}

    );

    // Number of subjects that appear in both the file and subjectOrder.
    uint32_t nSubjects() const {
        return m_nSubj;
    }

    // Number of stored entries (lower triangle + diagonal, no double-counting).
    size_t nnz() const {
        return m_entries.size();
    }

    // Compute  x^T G x  using half-storage:
    //   sum_diag(value * x[i]²) + 2 * sum_offdiag(value * x[i] * x[j]).
    // `x` must have size == nSubjects() and use the same order as subjectOrder.
    double quadForm(
        const double *x,
        uint32_t n
    ) const;

    // Compute  x^T G y  (full bilinear form, symmetric G):
    //   sum_diag(value * x[i]*y[i]) + sum_offdiag(value * (x[i]*y[j] + x[j]*y[i])).
    double bilinearForm(
        const double *x,
        const double *y,
        uint32_t n
    ) const;

    // Matrix-vector multiply:  result = G * x
    // Handles lower-tri symmetry internally (result is size n).
    void multiply(
        const double *x,
        double *result,
        uint32_t n
    ) const;

    // SPAmixPlus variance:  2 * sum(grm * R[i]*R[j]) - R·R
    //   (summing over lower-tri + diagonal entries).
    double spaVariance(
        const double *R,
        uint32_t n
    ) const;

    const std::vector<Entry> &entries() const {
        return m_entries;
    }

    // Cached diagonal of G (self-relatedness per subject).
    // diagonal()[i] == G(i,i).  Default 0.0 if subject has no diagonal entry.
    const std::vector<double> &diagonal() const {
        return m_diagonal;
    }

    // ── Shared GCTA .grm.id reader ────────────────────────────────────
    // Given a .grm.sp file path, derive the .grm.id path and read IIDs.
    // If the .grm.id file does not exist, returns an empty vector.
    static std::vector<std::string> readGctaIIDs(const std::string &spFile);

    // ── Convenience loader ─────────────────────────────────────────────
    // Dispatches to fromGCTA() or the text constructor based on which
    // path is non-empty.  Exactly one of grabFile / gctaFile must be set.
    static SparseGRM load(
        const std::string &grabFile,
        const std::string &gctaFile,
        const std::vector<std::string> &subjectOrder,
        const std::vector<std::string> &famIIDs = {}

    );

    // ── From pre-built entries (re-indexing) ──────────────────────────
    static SparseGRM fromEntries(
        uint32_t nSubj,
        std::vector<Entry> entries
    );

    // ── Lightweight subject-ID scanner ─────────────────────────────────
    // Scans GRM file(s) and returns the set of unique subject IDs that
    // appear in the file, WITHOUT loading any numeric entries.
    // Used for the subject intersection pipeline (GRM step).
    //   grabFile  — GRAB format path (empty to skip)
    //   gctaFile  — GCTA/plink2 .grm.sp path (empty to skip)
    //   famIIDs   — .fam IIDs (fallback when GCTA .grm.id is absent)
    static std::unordered_set<std::string> parseSubjectIDs(
        const std::string &grabFile,
        const std::string &gctaFile,
        const std::vector<std::string> &famIIDs = {}

    );

  private:
    SparseGRM() = default; // used by fromGCTA factory
    void buildDiagonal();  // populate m_diagonal from m_entries

    uint32_t m_nSubj = 0;
    std::vector<Entry> m_entries;   // lower-tri + diagonal COO
    std::vector<double> m_diagonal; // cached G(i,i) per subject
};
