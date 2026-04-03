// marker.hpp — Thread-pool marker-level association engine (pure C++17 / Eigen)
#pragma once

#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>

class GenoMeta;

// ======================================================================
// MethodBase — abstract interface for statistical methods.
//
// Each concrete method (SPAGRM, SPACox, WtCoxG, LEAF, …) derives from
// this and implements clone(), resultSize(), getHeaderColumns(), and
// getResultVec().  The engine clones one copy per worker thread.
// ======================================================================

class MethodBase {
public:
  virtual ~MethodBase() = default;

  // Deep-copy for per-thread isolation (mutable scratch state, etc.).
  virtual std::unique_ptr<MethodBase> clone() const = 0;

  // Number of result columns after the 9 meta columns.
  virtual int resultSize() const = 0;

  // Tab-separated column headers for the result columns (leading tab included).
  virtual std::string getHeaderColumns() const = 0;

  // Called once per chunk before processing markers.  Default: no-op.
  virtual void prepareChunk(const std::vector<uint64_t>& /*gIndices*/) {}

  // Per-marker analysis.  The engine handles QC and imputation before
  // calling this.  `result` is cleared before each call.
  //
  //   GVec             — imputed genotype vector (counts bim col5 = ALT allele)
  //   altFreq          — ALT allele frequency (= freq of bim col5)
  //   markerInChunkIdx — 0-based index within current chunk
  //   result           — output: method-specific result values
  virtual void getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int markerInChunkIdx,
    std::vector<double>& result) = 0;
};


// ======================================================================
// Engine entry point
// ======================================================================

void markerEngine(
  const GenoMeta& genoData,
  const MethodBase& method,
  const std::string& outputFile,
  int nthreads,
  double missingCutoff,
  double minMafCutoff,
  double minMacCutoff,
  bool exactHwe);
