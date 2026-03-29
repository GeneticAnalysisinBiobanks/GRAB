// spamix.hpp — SPAmix: Saddlepoint Approximation for Admixed Populations
//
// Port of ref_code/src/mtSPAmix into the pure C++17 / Eigen framework.
//
// Workflow:
//   1. Detect outlier residuals (IQR-based)
//   2. Per-marker: estimate per-individual MAF using PCs (linear model / logistic)
//   3. Score test → normal approx or SPA with outlier/non-outlier split
//   4. Output: [Pvalue, zScore]
#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>
#include <Eigen/Dense>

#include "engine/marker.hpp"
#include "spamix/common.hpp"
#include "spamix/indiv_af.hpp"

// ======================================================================
// SPAmixMethod — MethodBase implementation
// ======================================================================

class SPAmixMethod : public MethodBase {
public:
  SPAmixMethod(const Eigen::VectorXd& residuals,
               const Eigen::VectorXd& resid2,
               const Eigen::MatrixXd& onePlusPCs,
               const Eigen::MatrixXd& XtX_inv_Xt,
               const Eigen::VectorXd& sqrt_XtX_inv_diag,
               const OutlierData& outlier,
               double spaCutoff);

  std::unique_ptr<MethodBase> clone() const override;
  int resultSize() const override { return 2; }
  std::string getHeaderColumns() const override { return "\tPvalue\tzScore"; }
  bool skipFlip() const override { return true; }
  void getResultVec(Eigen::Ref<Eigen::VectorXd> GVec,
                    double altFreq, int markerInChunkIdx,
                    bool flipped, std::vector<double>& result) override;

private:
  // Score test + SPA.  Returns p-value; writes zScore.
  double getMarkerPval(const Eigen::Ref<const Eigen::VectorXd>& GVec,
                       double altFreq, double& zScore);

  // Read-only shared data (stable references — owner outlives all clones)
  const Eigen::VectorXd& m_resid;
  const Eigen::VectorXd& m_resid2;
  const Eigen::MatrixXd& m_onePlusPCs;
  const Eigen::MatrixXd& m_XtX_inv_Xt;
  const Eigen::VectorXd& m_sqrt_XtX_inv_diag;
  const OutlierData&     m_outlier;
  int    m_N;
  int    m_nPC;
  double m_spaCutoff;

  // Per-thread scratch (mutable, freshly allocated in clone)
  Eigen::VectorXd m_AFVec;
  Eigen::VectorXd m_mafOutlier;
  Eigen::VectorXd m_mafNonOutlier;
};

// ======================================================================
// Orchestration
// ======================================================================

void runSPAmix(
    const std::string& residFile,
    const std::string& eigenVecsFile,
    const std::string& bfilePrefix,
    const std::string& outputFile,
    double spaCutoff,
    double outlierRatio,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff);
