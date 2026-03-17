#ifndef SAGELD_H
#define SAGELD_H

// SAGELD.h -- SAGE linkage-disequilibrium / GALLOP family-based association testing

#include <RcppArmadillo.h>
#include <limits>
#include <boost/math/distributions/normal.hpp>

#include "UTIL.h"
#include <vector>

namespace SAGELD{

class SAGELDClass {
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


  // Per-marker updated data for three-or-more-subject families.
  struct UpdatedThreeSubj {
    arma::vec stand_S;
    arma::vec arr_prob;
  };

  const std::string& getMethod() const { return m_Method; }

private:


  // ---- Members ----
  std::string m_Method;
  arma::mat m_XTs;
  arma::mat m_SS;
  arma::mat m_AtS;
  arma::mat m_Q;
  arma::mat m_A21;
  arma::mat m_TTs;
  arma::mat m_Tys;
  arma::vec m_sol;
  arma::vec m_blups;
  double m_sig;
  double m_ncov;

  arma::vec m_resid;
  arma::vec m_resid_G;
  arma::vec m_resid_GxE;
  arma::vec m_resid_E;

  arma::vec m_resid_unrelated_outliers;
  arma::vec m_resid_unrelated_outliers_G;
  arma::vec m_resid_unrelated_outliers_GxE;

  double m_sum_unrelated_outliers2;
  double m_sum_unrelated_outliers_G2;
  double m_sum_unrelated_outliers_GxE2;
  double m_sum_unrelated_outliers_G_GxE2;

  double m_sum_R_nonOutlier;
  double m_sum_R_nonOutlier_G;
  double m_sum_R_nonOutlier_GxE;

  arma::vec m_R_GRM_R;               // [0]=main, [1]=G, [2]=GxE, [3]=G_GxE, [4]=E
  arma::vec m_R_GRM_R_nonOutlier;    // [0]=main, [1]=G, [2]=GxE, [3]=G_GxE
  arma::vec m_R_GRM_R_TwoSubjOutlier; // [0]=main, [1]=G, [2]=GxE, [3]=G_GxE

  std::vector<TwoSubjFamily> m_TwoSubj_list;
  std::vector<ThreeSubjFamily> m_ThreeSubj_list;

  arma::vec m_MAF_interval;
  double m_zScoreE_cutoff;
  double m_SPA_Cutoff;
  double m_zeta;
  double m_tol;

  arma::vec m_pvalVec;
  arma::vec m_zScoreVec;
  arma::vec m_BetaVec;
  arma::vec m_seBetaVec;

public:


  // ---- Public interface ----
  SAGELDClass(std::string Method,
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
              double tol);

  arma::vec getpvalVec(){return m_pvalVec;}
  arma::vec getzScoreVec(){return m_zScoreVec;}
  arma::vec getBetaVec(){return m_BetaVec;}
  arma::vec getseBetaVec(){return m_seBetaVec;}

  arma::mat mgf(double t,
                    const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
                    double MAF);

  arma::mat mgf(double t,
                    const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
                    double MAF,
                    double lambda_i);

  double fastGetRoot(const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
                         double Score,
                         double MAF,
                         double init_t,
                         double tol,
                         int maxiter = 50);

  double fastGetRoot(const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
                         double m_R_GRM_R_nonOutlier_i,
                         double lambda_i,
                         double Score,
                         double MAF,
                         double init_t,
                         double tol,
                         int maxiter = 50);

  double getProbSpa(const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
                     double Score,
                     double MAF,
                     bool lower_tail,
                     double zeta,
                     double tol);

  double getProbSpa(const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
                     double m_R_GRM_R_nonOutlier_i,
                     double lambda_i,
                     double Score,
                     double MAF,
                     bool lower_tail,
                     double zeta,
                     double tol);

  double getMarkerPval(arma::vec GVec,
                       double altFreq,
                       double& hwepval,
                       double hwepvalCutoff = 0.1);

};

}

#endif
