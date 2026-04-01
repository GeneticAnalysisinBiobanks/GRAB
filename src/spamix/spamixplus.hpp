// spamixplus.hpp — SPAmixPlus: SPA for admixed populations with sparse GRM
//
// Port of ref_code/src/mtSPAmixPlus into the pure C++17 / Eigen framework.
//
// Differences from SPAmix:
//   - Variance uses sparse GRM:  2·Σ_raw grm·R[i]·R[j] − R·R
//     instead of the diagonal   Σ resid²·2·AF·(1−AF)
//   - SPA tail probability is applied to variance-ratio-corrected S
//   - AF models are always pre-computed (from file or on-the-fly at startup)
//
// Output columns: [Pvalue, zScore]  (same as SPAmix)
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "engine/marker.hpp"
#include "io/sparse_grm.hpp"
#include "spamix/common.hpp"
#include "spamix/indiv_af.hpp"

// ======================================================================
// SPAmixPlusMethod — MethodBase implementation
// ======================================================================

class SPAmixPlusMethod : public MethodBase {
public:
  // Pre-computed AF mode
  SPAmixPlusMethod(
      const Eigen::VectorXd& residuals,
      const Eigen::VectorXd& resid2,
      const Eigen::MatrixXd& onePlusPCs,
      const OutlierData& outlier,
      double spaCutoff,
      const SparseGRM& grm,
      const std::vector<AFModel>& afModels,
      const std::vector<uint32_t>& genoToFlat);

  // On-the-fly AF mode (single PLINK pass — AF computed during marker testing)
  SPAmixPlusMethod(
      const Eigen::VectorXd& residuals,
      const Eigen::VectorXd& resid2,
      const Eigen::MatrixXd& onePlusPCs,
      const OutlierData& outlier,
      double spaCutoff,
      const SparseGRM& grm,
      const Eigen::MatrixXd& XtX_inv_Xt,
      const Eigen::VectorXd& sqrt_XtX_inv_diag,
      int nPC);

  std::unique_ptr<MethodBase> clone() const override;
  int resultSize() const override { return 2; }
  std::string getHeaderColumns() const override { return "\tPvalue\tzScore"; }
  bool skipFlip() const override { return true; }

  void prepareChunk(const std::vector<uint64_t>& gIndices) override;
  void getResultVec(Eigen::Ref<Eigen::VectorXd> GVec,
                    double altFreq, int markerInChunkIdx,
                    bool flipped, std::vector<double>& result) override;

private:
  double getMarkerPval(const Eigen::Ref<const Eigen::VectorXd>& GVec,
                       double altFreq, double& zScore);

  // Read-only shared data (stable references — owner outlives all clones)
  const Eigen::VectorXd&        m_resid;
  const Eigen::VectorXd&        m_resid2;        // shared across clones
  const Eigen::MatrixXd&        m_onePlusPCs;
  const OutlierData&            m_outlier;
  double                        m_spaCutoff;
  const SparseGRM&              m_grm;
  int                           m_N;
  int                           m_nPC;

  // Pre-computed AF (non-null in pre-computed mode)
  const std::vector<AFModel>*   m_afModels;
  const std::vector<uint32_t>*  m_genoToFlat;

  // On-the-fly AF (non-null in on-the-fly mode)
  const Eigen::MatrixXd*        m_XtX_inv_Xt;
  const Eigen::VectorXd*        m_sqrt_XtX_inv_diag;

  // Per-thread scratch (mutable, freshly allocated in clone)
  Eigen::VectorXd m_AFVec;
  Eigen::VectorXd m_R_new;
  Eigen::VectorXd m_mafOutlier;
  Eigen::VectorXd m_mafNonOutlier;

  // Chunk state
  std::vector<uint64_t> m_chunkGenoIndices;
};

// ======================================================================
// Orchestration
// ======================================================================

void runSPAmixPlus(
    const std::string& residFile,
    const std::string& eigenVecsFile,
    const std::string& bfilePrefix,
    const std::string& spgrmSaigeFile,
    const std::string& spgrmGctaPrefix,
    const std::string& afFile,          // empty → compute on-the-fly
    const std::string& outputFile,
    double spaCutoff,
    double outlierRatio,
    int    nthread,
    int    nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff);
