#ifndef SAGELD_H
#define SAGELD_H

// SAGELD.h -- SAGE linkage-disequilibrium / GALLOP family-based association testing

#include <RcppArmadillo.h>
#include <limits>
#include <boost/math/distributions/normal.hpp>
#include <vector>

class mtSAGELDClass {
public:

  // ---- Public types ----
  struct TwoSubjFamily {
    arma::vec Resid;
    arma::vec Rho;
    arma::vec Resid_G;
    arma::vec Resid_GxE;
  };

  struct ThreeSubjFamily {
    arma::mat CLT;
    arma::vec stand_S;
    arma::vec stand_S_G;
    arma::vec stand_S_GxE;
  };

  const std::string& getMethod() const { return m_Method; }

  std::string getHeaderColumns() const {
    if (m_Method == "GALLOP")
      return "\tisGALLOP\tBeta_G\tBeta_GxE\tSE_G\tSE_GxE\tPvalue_G\tPvalue_GxE";
    else
      return "\tisGALLOP\tzScore_G\tzScore_GxE\tPvalue_G\tPvalue_GxE";
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

  const arma::vec m_R_GRM_R;               // [0]=main, [1]=G, [2]=GxE, [3]=G_GxE, [4]=E
  const arma::vec m_R_GRM_R_nonOutlier;    // [0]=main, [1]=G, [2]=GxE, [3]=G_GxE
  const arma::vec m_R_GRM_R_TwoSubjOutlier; // [0]=main, [1]=G, [2]=GxE, [3]=G_GxE

  const std::vector<TwoSubjFamily> m_TwoSubj_list;
  const std::vector<ThreeSubjFamily> m_ThreeSubj_list;

  // Extracted from m_TwoSubj_list for use with nsSPAGRM free functions
  std::vector<arma::vec> m_TwoSubj_resid_list;
  std::vector<arma::vec> m_TwoSubj_rho_list;

  const arma::vec m_MAF_interval;
  const double m_zScoreE_cutoff;
  const double m_SPA_Cutoff;
  const double m_zeta;
  const double m_tol;

  arma::vec m_pvalVec;
  arma::vec m_zScoreVec;
  arma::vec m_BetaVec;
  arma::vec m_seBetaVec;

public:

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
    arma::vec R_GRM_R,
    arma::vec R_GRM_R_nonOutlier,
    arma::vec R_GRM_R_TwoSubjOutlier,
    std::vector<TwoSubjFamily> TwoSubj_list,
    std::vector<ThreeSubjFamily> ThreeSubj_list,
    arma::vec MAF_interval,
    double zScoreE_cutoff,
    double SPA_Cutoff,
    double zeta,
    double tol
  );

  arma::vec getpvalVec(){return m_pvalVec;}
  arma::vec getzScoreVec(){return m_zScoreVec;}
  arma::vec getBetaVec(){return m_BetaVec;}
  arma::vec getseBetaVec(){return m_seBetaVec;}

  // GALLOP: returns [isGALLOP, bG, bGxE, seG, seGxE, pG, pGxE] (7)
  // else:   returns [isGALLOP, zG, zGxE, pG, pGxE]              (5)
  std::vector<double> getResultVec(arma::vec GVec, double altFreq) {
    getMarkerPval(std::move(GVec), altFreq);
    if (m_Method == "GALLOP") {
      double flipSign = (altFreq > 0.5) ? -1.0 : 1.0;
      return {1.0, m_BetaVec[0] * flipSign, m_BetaVec[1] * flipSign,
              m_seBetaVec[0], m_seBetaVec[1],
              m_pvalVec[0], m_pvalVec[1]};
    }
    return {0.0, m_zScoreVec[0], m_zScoreVec[1], m_pvalVec[0], m_pvalVec[1]};
  }

  int resultSize() const { return (m_Method == "GALLOP") ? 7 : 5; }

  double getMarkerPval(
    arma::vec GVec,
    double altFreq
  );

};

#endif
