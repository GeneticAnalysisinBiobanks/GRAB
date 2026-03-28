// sparse_grm.hpp — Sparse GRM reader and quadratic-form evaluator
//
// File format: 3-column TSV  (ID1 \t ID2 \t value)
// - Symmetric: only lower triangle (or arbitrary half) stored.
// - Diagonal entries have ID1 == ID2.
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

  // Parse the 3-column file.  Subjects not in `subjectOrder` are silently
  // dropped.  `subjectOrder` defines the canonical index mapping so that
  // quadratic-form vectors can use the same ordering.
  SparseGRM(const std::string& filename,
            const std::vector<std::string>& subjectOrder);

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

  const std::vector<Entry>& entries() const { return m_entries; }

private:
  uint32_t             m_nSubj = 0;
  std::vector<Entry>   m_entries;   // symmetrised COO
};
