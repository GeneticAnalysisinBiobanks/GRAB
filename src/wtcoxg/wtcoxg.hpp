// wtcoxg.hpp — WtCoxG: Weighted Cox-type G-test (Li et al., 2025)
//
// Phases:
//   1. Data loading & marker matching  (resid-file, bfile, ref-af-file)
//   2. Batch-effect testing & parameter estimation (TPR, sigma2, w.ext, var ratios)
//   3. Marker-level SPA tests via markerEngine (MethodBase interface)
#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>
#include <memory>
#include <array>
#include <cmath>
#include <limits>
#include <Eigen/Dense>

#include "engine/marker.hpp"
#include "io/resid_file.hpp"

class PlinkData;
class SparseGRM;

// ======================================================================
// Per-marker reference info  (populated during Phase 2)
// ======================================================================

struct WtCoxGRefInfo {
  double AF_ref       = std::numeric_limits<double>::quiet_NaN();
  double N_ref        = std::numeric_limits<double>::quiet_NaN();  // sample size
  double TPR          = std::numeric_limits<double>::quiet_NaN();
  double sigma2       = std::numeric_limits<double>::quiet_NaN();
  double pvalue_bat   = std::numeric_limits<double>::quiet_NaN();
  double w_ext        = std::numeric_limits<double>::quiet_NaN();
  double var_ratio_w0 = 1.0;
  double var_ratio_int = 1.0;
  double var_ratio_ext = 1.0;
};

// ======================================================================
// WtCoxGMethod — implements MethodBase for marker engine
// ======================================================================

class WtCoxGMethod : public MethodBase {
public:
  // Construct from pre-computed null-model quantities.
  //   R:           martingale residuals (nSubj)
  //   w:           sampling weights     (nSubj)
  //   cutoff:      batch-effect p-value threshold (e.g. 0.05)
  //   SPA_Cutoff:  z-score threshold to switch from normal to SPA
  //   refMap:      genoIndex → WtCoxGRefInfo (from Phase 2)
  WtCoxGMethod(Eigen::VectorXd R,
               Eigen::VectorXd w,
               double cutoff,
               double SPA_Cutoff,
               std::shared_ptr<const std::unordered_map<uint64_t, WtCoxGRefInfo>> refMap);

  // ---- MethodBase interface ----
  std::unique_ptr<MethodBase> clone() const override;
  int resultSize() const override { return 7; }
  std::string getHeaderColumns() const override;
  bool skipFlip() const override { return true; }  // WtCoxG handles alleles itself
  void prepareChunk(const std::vector<uint64_t>& gIndices) override;
  void getResultVec(Eigen::Ref<Eigen::VectorXd> GVec,
                    double altFreq, int markerInChunkIdx,
                    bool flipped, std::vector<double>& result) override;

  // For LEAF: compute ext/noext results with raw scores for meta-analysis.
  struct DualResult {
    double p_ext, p_noext;
    double score_ext, score_noext;
  };
  DualResult computeDual(Eigen::Ref<Eigen::VectorXd> GVec, int markerInChunkIdx);

  // Access per-chunk ref info (for LEAF meta-analysis).
  const WtCoxGRefInfo& chunkRefInfoAt(int idx) const { return m_chunkRefInfo[idx]; }

private:
  struct WtResult { double pval; double score; double zscore; };

  // Core SPA test for one marker with external adjustment.
  WtResult wtCoxGTest(
    const Eigen::Ref<const Eigen::VectorXd>& g,
    double p_bat, double TPR, double sigma2, double b,
    double var_ratio_int, double var_ratio_w0, double var_ratio_w1,
    double var_ratio0, double var_ratio1,
    double mu_ext, double n_ext, double p_cut) const;

  // Members (const after construction, except per-chunk scratch)
  Eigen::VectorXd m_R;
  Eigen::VectorXd m_w;
  Eigen::VectorXd m_w1;      // w / (2 * sum(w))
  double          m_meanR;
  double          m_sumR;
  double          m_cutoff;
  double          m_SPA_Cutoff;
  std::shared_ptr<const std::unordered_map<uint64_t, WtCoxGRefInfo>> m_refMap;

  // Per-chunk scratch (rebuilt in prepareChunk)
  std::vector<WtCoxGRefInfo> m_chunkRefInfo;
  std::array<double, 2>      m_scoreArr;
  std::array<double, 2>      m_zScoreArr;
};


// ======================================================================
// Phase 1 — Data loading & marker matching
// ======================================================================

struct RefAfRecord {
  std::string chrom;
  uint32_t    pos;
  std::string a1;
  std::string a2;
  double      AF_ref;  // A1 allele frequency
  double      N_ref;   // sample size (file column N)
};

// Parse a 6-column ref-af file: #CHROM  POS  A1  A2  A1F  N
// Header line (starting with '#') is skipped. N is sample size (AN = 2*N).
std::vector<RefAfRecord> loadRefAfFile(const std::string& filename);

// Match bim markers to ref-af records by (chr, bp, a1, a2) — exact only.
// Returns a map: genoIndex → {AF_ref, N_ref, ...} for matched markers.
// mu0, mu1, n0, n1 per marker are computed from the genotype data.
struct MatchedMarkerInfo {
  uint64_t genoIndex;
  double   AF_ref;
  double   N_ref;   // sample size
  double   mu0;        // control allele freq
  double   mu1;        // case allele freq
  double   n0;         // effective control count
  double   n1;         // effective case count
  double   mu_int;     // internal MAF = 0.5*mu0 + 0.5*mu1 (folded ≤0.5)
};

std::vector<MatchedMarkerInfo> matchMarkers(
    const PlinkData& plinkData,
    const std::vector<RefAfRecord>& refAf);

// Scan genotypes to compute per-marker case/control allele frequencies.
// Populates mu0, mu1, n0, n1, mu_int in each MatchedMarkerInfo.
void computeMarkerStats(
    std::vector<MatchedMarkerInfo>& matched,
    const PlinkData& plinkData,
    const ResidData& resid);


// ======================================================================
// Phase 2 — Batch-effect testing and parameter estimation
// ======================================================================

// Full batch-effect QC pipeline:
//   1. Compute batch p-values per marker
//   2. Compute variance ratios from sparse GRM (if provided)
//   3. Estimate TPR, sigma2 per MAF group (Nelder-Mead)
//   4. Compute optimal w.ext per MAF group (Brent 1D)
//   5. Compute var_ratio_ext per MAF group from sparse GRM (if provided)
//
// Populates refInfo for each matched marker and returns a shared map.
std::shared_ptr<std::unordered_map<uint64_t, WtCoxGRefInfo>>
testBatchEffects(
    const std::vector<MatchedMarkerInfo>& matched,
    const ResidData& resid,
    const SparseGRM* grm,       // nullptr if no GRM
    double refPrevalence,
    double cutoff);


// ======================================================================
// Top-level orchestration — called from main()
// ======================================================================

void runWtCoxG(
    const std::string& residFile,
    const std::string& bfilePrefix,
    const std::string& refAfFile,
    const std::string& sparseGrmFile,  // empty = no GRM
    const std::string& outputFile,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff);
