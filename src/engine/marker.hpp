// marker.hpp — Thread-pool marker-level association engine (pure C++17 / Eigen)
#pragma once

#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>

class PlinkData;

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

  // If true, the engine will NOT flip GVec when altFreq > 0.5.
  virtual bool skipFlip() const { return false; }

  // Per-marker analysis.  The engine handles QC, imputation, and flip
  // before calling this.  `result` is cleared before each call.
  //
  //   GVec             — imputed genotype vector (may be flipped)
  //   altFreq          — alternate allele frequency (pre-flip)
  //   markerInChunkIdx — 0-based index within current chunk
  //   flipped          — true if we flipped the genotype vector
  //   result           — output: method-specific result values
  virtual void getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int markerInChunkIdx,
    bool flipped,
    std::vector<double>& result) = 0;
};


// ======================================================================
// Engine entry point
// ======================================================================

void markerEngine(
  const PlinkData& plinkData,
  const MethodBase& method,
  const std::string& outputFile,
  int nthreads,
  double missingCutoff,
  double minMafMarker,
  double minMacMarker,
  bool exactHwe);
