// spagrm.hpp — Saddlepoint approximation with GRM (pure C++17 / Eigen / Boost)
//
// Translated from mtSPAGRM.h / mtSPAGRM.cpp (RcppArmadillo) to Eigen.
// Provides the nsSPAGRM namespace (mgf, fastGetRoot, getProbSpa) and
// the SPAGRMClass marker-level evaluator used by SPAsqr and SPAGRM.
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>
#include <memory>
#include <Eigen/Dense>

#include "engine/marker.hpp"   // MethodBase
#include "io/geno_data.hpp"    // GenoSpec

namespace nsSPAGRM {

// Per-marker updated data for three-or-more-subject families.
struct UpdatedThreeSubj {
  std::vector<double> stand_S;
  std::vector<double> arr_prob;
};

// Family data for one tau (or one residual column).
struct FamilyData {
  Eigen::VectorXd                        resid_unrelated_outliers;
  std::vector<std::array<double, 2>>     twoSubj_resid;
  std::vector<std::vector<double>>       twoSubj_rho;
  std::vector<std::vector<double>>       threeSubj_standS;
  std::vector<Eigen::MatrixXd>           threeSubj_CLT;
};

// Total number of elements in the mgf output vectors.
inline size_t mgfOutputSize(
    size_t n_unrelated,
    const std::vector<std::vector<double>>& TwoSubj_rho,
    size_t nThreeSubj) {
  size_t sz = n_unrelated;
  for (const auto& r : TwoSubj_rho) sz += r.size();
  sz += nThreeSubj;
  return sz;
}

// Pre-allocated scratch workspace (one per SPAGRMClass instance).
struct MgfWorkspace {
  Eigen::VectorXd MGF0, MGF1, MGF2;
  Eigen::VectorXd temp;
  // Unrelated-outlier intermediates
  Eigen::VectorXd ul_lambda, ul_alpha, ul_alpha_1, ul_alpha_2;

  MgfWorkspace() = default;
  MgfWorkspace(Eigen::Index mgfSz, Eigen::Index n_unrelated)
    : MGF0(mgfSz), MGF1(mgfSz), MGF2(mgfSz), temp(mgfSz),
      ul_lambda(n_unrelated), ul_alpha(n_unrelated),
      ul_alpha_1(n_unrelated), ul_alpha_2(n_unrelated)
  {}
};

// MGF and its first two derivatives → ws.MGF0/MGF1/MGF2.
void mgf(
    double t,
    const Eigen::VectorXd& resid_unrelated_outliers,
    const std::vector<std::array<double, 2>>& TwoSubj_resid,
    const std::vector<std::vector<double>>& TwoSubj_rho,
    const std::vector<UpdatedThreeSubj>& threeSubj,
    double MAF,
    MgfWorkspace& ws);

// Newton-Raphson root finder for the CGF equation.
double fastGetRoot(
    const Eigen::VectorXd& resid_unrelated_outliers,
    const std::vector<std::array<double, 2>>& TwoSubj_resid,
    const std::vector<std::vector<double>>& TwoSubj_rho,
    const std::vector<UpdatedThreeSubj>& threeSubj,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score, double MAF,
    double init_t, double tol,
    MgfWorkspace& ws,
    int maxiter = 50);

// Saddlepoint approximation tail probability.
double getProbSpa(
    const Eigen::VectorXd& resid_unrelated_outliers,
    const std::vector<std::array<double, 2>>& TwoSubj_resid,
    const std::vector<std::vector<double>>& TwoSubj_rho,
    const std::vector<UpdatedThreeSubj>& threeSubj,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score, double MAF,
    bool lower_tail, double zeta, double tol,
    MgfWorkspace& ws);

} // namespace nsSPAGRM


// ══════════════════════════════════════════════════════════════════════
// SPAGRMClass — per-column (per-tau) marker evaluator
// ══════════════════════════════════════════════════════════════════════

class SPAGRMClass {
public:
  SPAGRMClass(
      Eigen::VectorXd resid,
      double sum_R_nonOutlier,
      double R_GRM_R_nonOutlier,
      double R_GRM_R_TwoSubjOutlier,
      double R_GRM_R,
      std::vector<double> MAF_interval,
      nsSPAGRM::FamilyData fam,
      double SPA_Cutoff,
      double zeta,
      double tol);

  // Deep copy (including scratch state) for per-thread isolation.
  SPAGRMClass(const SPAGRMClass& o);
  SPAGRMClass& operator=(const SPAGRMClass&) = delete;

  double getMarkerPval(
      const Eigen::VectorXd& GVec,
      double altFreq,
      double& zScore);

private:
  Eigen::VectorXd m_resid;
  Eigen::VectorXd m_resid_unrelated_outliers;
  double m_sum_unrelated_outliers2;
  double m_sum_R_nonOutlier;
  double m_R_GRM_R_nonOutlier;
  double m_R_GRM_R_TwoSubjOutlier;
  double m_R_GRM_R;
  std::vector<double> m_MAF_interval;
  std::vector<std::array<double, 2>> m_TwoSubj_resid_list;
  std::vector<std::vector<double>> m_TwoSubj_rho_list;
  std::vector<std::vector<double>> m_ThreeSubj_standS_list;
  std::vector<Eigen::MatrixXd> m_ThreeSubj_CLT_list;

  double m_SPA_Cutoff;
  double m_zeta;
  double m_tol;

  nsSPAGRM::MgfWorkspace m_workspace;
  std::vector<nsSPAGRM::UpdatedThreeSubj> m_threeSubj_scratch;
};


// ══════════════════════════════════════════════════════════════════════
// SPAGRMMethod — MethodBase adapter wrapping a single SPAGRMClass
// ══════════════════════════════════════════════════════════════════════

class SPAGRMMethod : public MethodBase {
public:
  explicit SPAGRMMethod(SPAGRMClass spagrm) : m_spagrm(std::move(spagrm)) {}

  std::unique_ptr<MethodBase> clone() const override {
    return std::make_unique<SPAGRMMethod>(*this);
  }

  int resultSize() const override { return 2; }

  std::string getHeaderColumns() const override {
    return "\tP\tZ";
  }

  void getResultVec(
      Eigen::Ref<Eigen::VectorXd> GVec,
      double altFreq,
      int /*markerInChunkIdx*/,
      std::vector<double>& result) override
  {
    result.clear();
    double z;
    double p = m_spagrm.getMarkerPval(GVec, altFreq, z);
    result.push_back(p);
    result.push_back(z);
  }

private:
  SPAGRMClass m_spagrm;
};


// ══════════════════════════════════════════════════════════════════════
// runSPAGRM — full workflow entry point
// ══════════════════════════════════════════════════════════════════════

void runSPAGRM(
    const std::string& residFile,
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const std::string& pairwiseIBDFile,
    const GenoSpec& geno,
    const std::string& outPrefix,
    const std::string& compression,
    int compressionLevel,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    const std::string& keepFile = {},
    const std::string& removeFile = {});
