// leaf.hpp — LEAF: cluster-stratified WtCoxG with ancestry-matched reference AFs
//
// Workflow:
//   1. Load per-cluster resid files → combined subject list + cluster indices
//   2. Load PLINK data with combined subjects
//   3. Load multi-population ref-af file (#CHROM POS A1 A2 A1F_POP1 N_POP1 ...)
//   4. Exact-match markers by (chr, pos, a1, a2) against .bim
//   5. Per cluster: compute internal AF → summix → ancestry-matched AF_ref/N_ref
//   6. Per cluster: batch-effect testing → refInfoMap
//   7. Create LEAFMethod (N WtCoxGMethod objects) → markerEngine → meta-analysis
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "engine/marker.hpp"

class WtCoxGMethod;
struct WtCoxGRefInfo;

// ======================================================================
// Multi-population reference AF record  (4+2K columns, '#' header line)
// Format: #CHROM  POS  A1  A2  A1F_POP1  N_POP1  [A1F_POP2  N_POP2  ...]
// ======================================================================

struct MultiPopRefAf {
  struct Record {
    std::string chrom;
    uint32_t    pos;
    std::string a1;     // same meaning as .bim A1
    std::string a2;     // same meaning as .bim A2
    std::vector<double> popAF;  // A1 frequency per reference population
    std::vector<double> popN;   // sample size N per reference population
  };
  int                 nPop = 0;    // number of reference populations
  std::vector<Record> records;
};

// Parse 4+2K column ref-af file; header line (starting with '#') is skipped.
// N_POPk is sample size (allele count AN = 2*N).
MultiPopRefAf loadMultiPopRefAfFile(const std::string& filename);


// ======================================================================
// Matched marker with multi-population allele frequencies
// ======================================================================

struct MultiPopMatchedMarker {
  uint64_t genoIndex;
  std::vector<double> popAF;  // per-pop A1 frequency (no flip)
  std::vector<double> popN;   // per-pop sample size N
};

// Exact match by (chr, pos, a1, a2) between .bim and ref-af.
// No allele flipping. Unmatched markers are dropped.
std::vector<MultiPopMatchedMarker> matchMultiPopMarkers(
    const class PlinkData& plinkData,
    const MultiPopRefAf&   refAf);


// ======================================================================
// Summix — ancestry proportion estimation
// ======================================================================

// Estimate ancestry proportions by minimising ||D*p − o||²
// subject to p ≥ 0, Σp = 1.
// Uses exact active-set enumeration (K ≤ 6 populations, up to 2^K subsets)
// with KKT / Lagrange-multiplier equality-constrained least squares.
//   observedAF: per-marker internal AF   (nSNP)
//   refAF:      reference pop AFs        (nSNP × nPop, nPop ≤ 6)
// Returns: proportion vector (nPop), sums to 1.
Eigen::VectorXd summixEstimate(
    const Eigen::VectorXd& observedAF,
    const Eigen::MatrixXd& refAF);


// ======================================================================
// LEAFMethod — MethodBase implementation for cluster-stratified WtCoxG
// ======================================================================

class LEAFMethod : public MethodBase {
public:
  // Constructs from per-cluster WtCoxGMethod objects and index arrays.
  //   clusterMethods:  one WtCoxGMethod per cluster (ownership transferred)
  //   clusterIndices:  per-cluster indices into the full subject vector
  LEAFMethod(
    std::vector<std::unique_ptr<WtCoxGMethod>> clusterMethods,
    std::vector<std::vector<uint32_t>>         clusterIndices);

  // ---- MethodBase interface ----
  std::unique_ptr<MethodBase> clone() const override;
  int  resultSize()       const override;
  std::string getHeaderColumns() const override;
  bool skipFlip()         const override { return true; }
  void prepareChunk(const std::vector<uint64_t>& gIndices) override;
  void getResultVec(Eigen::Ref<Eigen::VectorXd> GVec,
                    double altFreq, int markerInChunkIdx,
                    bool flipped, std::vector<double>& result) override;

private:
  int m_nCluster;
  std::vector<std::unique_ptr<WtCoxGMethod>> m_clusterMethods;
  std::vector<std::vector<uint32_t>>         m_clusterIndices;
  std::vector<Eigen::VectorXd>               m_clusterGVec;  // pre-allocated buffers
};


// ======================================================================
// Top-level orchestration — called from main()
// ======================================================================

void runLEAF(
    const std::vector<std::string>& residFiles,     // one per cluster
    const std::string& bfilePrefix,
    const std::string& refAfFile,
    const std::string& sparseGrmFile,               // empty = no GRM
    const std::string& outputFile,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff);
