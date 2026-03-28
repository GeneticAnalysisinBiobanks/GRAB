#ifndef SAGELD_H
#define SAGELD_H

// SAGELD.h -- SAGE linkage-disequilibrium / GALLOP family-based association testing

#include <RcppArmadillo.h>
#include <limits>
#include <boost/math/distributions/normal.hpp>
#include <vector>
#include "mtSPAGRM.h"

class mtSAGELDClass {
public:

  // ---- Public types ----
  struct TwoSubjFamily {
    std::array<double, 2> Resid;
    std::vector<double> Rho;
    std::array<double, 2> Resid_G;
    std::array<double, 2> Resid_GxE;
  };

  struct ThreeSubjFamily {
    arma::mat CLT;
    std::vector<double> stand_S;
    std::vector<double> stand_S_G;
    std::vector<double> stand_S_GxE;
  };

  // ---- Public interface ----
  mtSAGELDClass(
    std::string Method,
    arma::mat XTs,
    arma::mat SS,
    arma::mat AtS,
    arma::mat Q,
    arma::mat A21,
    arma::mat TTs,
    arma::mat Tys,
    arma::vec sol,
    arma::vec blups,
    double sig,
    arma::vec resid,
    arma::vec resid_G,
    arma::vec resid_GxE,
    arma::vec resid_E,
    arma::vec resid_unrelated_outliers,
    arma::vec resid_unrelated_outliers_G,
    arma::vec resid_unrelated_outliers_GxE,
    double sum_R_nonOutlier,
    double sum_R_nonOutlier_G,
    double sum_R_nonOutlier_GxE,
    double R_GRM_R,
    double R_GRM_R_G,
    double R_GRM_R_GxE,
    double R_GRM_R_G_GxE,
    double R_GRM_R_E,
    double R_GRM_R_nonOutlier,
    double R_GRM_R_nonOutlier_G,
    double R_GRM_R_nonOutlier_GxE,
    double R_GRM_R_nonOutlier_G_GxE,
    double R_GRM_R_TwoSubjOutlier,
    double R_GRM_R_TwoSubjOutlier_G,
    double R_GRM_R_TwoSubjOutlier_GxE,
    double R_GRM_R_TwoSubjOutlier_G_GxE,
    std::vector<TwoSubjFamily> TwoSubj_list,
    std::vector<ThreeSubjFamily> ThreeSubj_list,
    std::vector<double> MAF_interval,
    double zScoreE_cutoff,
    double SPA_Cutoff,
    double zeta,
    double tol
  );

  const std::string& getMethod() const { return m_Method; }

  std::string getHeaderColumns() const {
    if (m_Method == "GALLOP")
      return "\tisGALLOP\tBeta_G\tBeta_GxE\tSE_G\tSE_GxE\tPvalue_G\tPvalue_GxE";
    else
      return "\tisGALLOP\tzScore_G\tzScore_GxE\tPvalue_G\tPvalue_GxE";
  }

  int resultSize() const { return (m_Method == "GALLOP") ? 7 : 5; }

  // GALLOP: fills rv with [isGALLOP, bG, bGxE, seG, seGxE, pG, pGxE] (7)
  // else:   fills rv with [isGALLOP, zG, zGxE, pG, pGxE]              (5)
  void getResultVec(const arma::vec& GVec, double altFreq, std::vector<double>& rv) {
    getMarkerPval(GVec, altFreq);
    rv.clear();
    if (m_Method == "GALLOP") {
      double flipSign = (altFreq > 0.5) ? -1.0 : 1.0;
      rv.push_back(1.0);
      rv.push_back(m_BetaVec[0] * flipSign);
      rv.push_back(m_BetaVec[1] * flipSign);
      rv.push_back(m_seBetaVec[0]);
      rv.push_back(m_seBetaVec[1]);
      rv.push_back(m_pvalVec[0]);
      rv.push_back(m_pvalVec[1]);
    } else {
      rv.push_back(0.0);
      rv.push_back(m_zScoreVec[0]);
      rv.push_back(m_zScoreVec[1]);
      rv.push_back(m_pvalVec[0]);
      rv.push_back(m_pvalVec[1]);
    }
  }

private:

  // ---- Members ----
  const std::string m_Method;

  const arma::mat m_XTs;
  const arma::mat m_SS;
  const arma::mat m_AtS;
  const arma::mat m_Q;
  const arma::mat m_A21;
  const arma::mat m_TTs;
  const arma::mat m_Tys;
  const arma::vec m_sol;
  const arma::vec m_blups;
  const double m_sig;
  const double m_ncov;

  const arma::vec m_resid;
  const arma::vec m_resid_G;
  const arma::vec m_resid_GxE;
  const arma::vec m_resid_E;

  const arma::vec m_resid_unrelated_outliers;
  const arma::vec m_resid_unrelated_outliers_G;
  const arma::vec m_resid_unrelated_outliers_GxE;

  const double m_sum_unrelated_outliers2;
  const double m_sum_unrelated_outliers_G2;
  const double m_sum_unrelated_outliers_GxE2;
  const double m_sum_unrelated_outliers_G_GxE2;

  const double m_sum_R_nonOutlier;
  const double m_sum_R_nonOutlier_G;
  const double m_sum_R_nonOutlier_GxE;

  const double m_R_GRM_R;
  const double m_R_GRM_R_G;
  const double m_R_GRM_R_GxE;
  const double m_R_GRM_R_G_GxE;
  const double m_R_GRM_R_E;
  const double m_R_GRM_R_nonOutlier;
  const double m_R_GRM_R_nonOutlier_G;
  const double m_R_GRM_R_nonOutlier_GxE;
  const double m_R_GRM_R_nonOutlier_G_GxE;
  const double m_R_GRM_R_TwoSubjOutlier;
  const double m_R_GRM_R_TwoSubjOutlier_G;
  const double m_R_GRM_R_TwoSubjOutlier_GxE;
  const double m_R_GRM_R_TwoSubjOutlier_G_GxE;

  const std::vector<TwoSubjFamily> m_TwoSubj_list;
  const std::vector<ThreeSubjFamily> m_ThreeSubj_list;

  // Extracted from m_TwoSubj_list for use with nsSPAGRM free functions
  std::vector<std::array<double, 2>> m_TwoSubj_resid_list;
  std::vector<std::vector<double>> m_TwoSubj_rho_list;

  const std::vector<double> m_MAF_interval;
  const double m_zScoreE_cutoff;
  const double m_SPA_Cutoff;
  const double m_zeta;
  const double m_tol;

  // Per-instance SPA scratch — pre-allocated at construction, reused across all marker calls.
  // Thread-safe: each worker thread owns its own mtSAGELDClass copy.
  nsSPAGRM::MgfWorkspace m_workspace;
  // Three-subject scratch: stand_S pre-set for path-1; overwritten per marker in path-2.
  std::vector<nsSPAGRM::UpdatedThreeSubj> m_threeSubj_scratch;
  // Path-2 (lambda_i) scratch: avoids per-marker heap allocation for adjusted residuals.
  arma::vec m_resid_outliers_i;
  std::vector<std::array<double, 2>> m_twoResid_i;

  arma::vec m_pvalVec;
  arma::vec m_zScoreVec;
  arma::vec m_BetaVec;
  arma::vec m_seBetaVec;

  double getMarkerPval(
    const arma::vec& GVec,
    double altFreq
  );

};

#endif
