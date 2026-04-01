// sparse_grm.hpp — Sparse GRM reader and quadratic-form evaluator
//
// Supports two formats:
//   SAIGE: 3-column TSV  (ID1 \t ID2 \t value), '#'-header skipped
//   GCTA:  plink2 --make-grm-sparse output (PREFIX.grm.sp + PREFIX.grm.id)
//          .grm.id has FID IID per line; .grm.sp has 0-based-idx1 idx2 value
//
// Dedicated method for the quadratic form  x^T G x  where G is the sparse GRM
// restricted to a given set of subjects.  This is needed by WtCoxG (variance
// ratios) and will be reused by SPAGRM, LEAF, etc.
#pragma once

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

class SparseGRM {
public:
  // Storage: COO (coordinate) format kept sorted by (row, col) for cache
  // locality when iterating.
  struct Entry {
    uint32_t row;   // index into m_ids
    uint32_t col;   // index into m_ids
    double   value;
  };

  // ── SAIGE format ───────────────────────────────────────────────────
  // Parse the 3-column file.  Subjects not in `subjectOrder` are silently
  // dropped.  `subjectOrder` defines the canonical index mapping so that
  // quadratic-form vectors can use the same ordering.
  //
  // If symmetrize=true (default), stores both (i,j) and (j,i) for every
  // off-diagonal entry so that quadForm gives the full x^T G x.
  // If symmetrize=false, stores entries exactly as they appear in the file
  // (lower triangle + diagonal).  Use spaVariance() with raw entries.
  SparseGRM(const std::string& filename,
            const std::vector<std::string>& subjectOrder,
            bool symmetrize = true);

  // ── GCTA format (plink2 --make-grm-sparse) ────────────────────────
  // Reads PREFIX.grm.id (FID IID) and PREFIX.grm.sp (idx1 idx2 value).
  // Subjects not in `subjectOrder` are silently dropped.
  static SparseGRM fromGCTA(const std::string& prefix,
                            const std::vector<std::string>& subjectOrder,
                            bool symmetrize = true);

  // Number of subjects that appear in both the file and subjectOrder.
  uint32_t nSubjects() const { return m_nSubj; }

  // Number of non-zero entries (counting both halves after symmetrisation).
  size_t nnz() const { return m_entries.size(); }

  // Compute  sum_{(i,j) in GRM} value * x[i] * x[j].
  // `x` must have size == nSubjects() and use the same order as subjectOrder.
  // Off-diagonal entries are counted twice (i,j) and (j,i) automatically.
  double quadForm(const double* x, uint32_t n) const;

  // Compute  sum_{(i,j) in GRM} value * x[i] * y[j].
  // Symmetric: returns the full bilinear form x^T G y.
  double bilinearForm(const double* x, const double* y, uint32_t n) const;

  // SPAmixPlus variance:  2 * sum_raw(grm * R[i]*R[j]) - R·R.
  // Requires raw (non-symmetrised) entries (symmetrize=false at construction).
  double spaVariance(const double* R, uint32_t n) const;

  const std::vector<Entry>& entries() const { return m_entries; }

  // ── Shared GCTA .grm.id reader ────────────────────────────────────
  // Read PREFIX.grm.id, return vector of IIDs (one per file line, IID = col 2).
  static std::vector<std::string> readGctaIIDs(const std::string& prefix);

private:
  SparseGRM() = default;  // used by fromGCTA factory
  uint32_t             m_nSubj = 0;
  std::vector<Entry>   m_entries;   // symmetrised COO
};
