
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "SAGELD.hpp"

namespace SAGELD {

SAGELDClass::SAGELDClass(std::string t_Method,
                         arma::mat t_XTs,
                         arma::mat t_SS,
                         arma::mat t_AtS,
                         arma::mat t_Q,
                         arma::mat t_A21,
                         arma::mat t_TTs,
                         arma::mat t_Tys,
                         arma::vec t_sol,
                         arma::vec t_blups,
                         double t_sig,
                         arma::vec t_resid,
                         arma::vec t_resid_G,
                         arma::vec t_resid_GxE,
                         arma::vec t_resid_E,
                         arma::vec t_resid_unrelated_outliers,
                         arma::vec t_resid_unrelated_outliers_G,
                         arma::vec t_resid_unrelated_outliers_GxE,
                         double t_sum_R_nonOutlier,
                         double t_sum_R_nonOutlier_G,
                         double t_sum_R_nonOutlier_GxE,
                         double t_R_GRM_R,
                         double t_R_GRM_R_G,
                         double t_R_GRM_R_GxE,
                         double t_R_GRM_R_G_GxE,
                         double t_R_GRM_R_E,
                         double t_R_GRM_R_nonOutlier,
                         double t_R_GRM_R_nonOutlier_G,
                         double t_R_GRM_R_nonOutlier_GxE,
                         double t_R_GRM_R_nonOutlier_G_GxE,
                         double t_R_GRM_R_TwoSubjOutlier,
                         double t_R_GRM_R_TwoSubjOutlier_G,
                         double t_R_GRM_R_TwoSubjOutlier_GxE,
                         double t_R_GRM_R_TwoSubjOutlier_G_GxE,
                         Rcpp::List t_TwoSubj_list,
                         Rcpp::List t_ThreeSubj_list,
                         arma::vec t_MAF_interval,
                         double t_zScoreE_cutoff,
                         double t_SPA_Cutoff,
                         double t_zeta,
                         double t_tol)
{
  m_Method = t_Method;
  
  m_XTs = t_XTs;
  m_SS = t_SS;
  m_AtS = t_AtS;
  m_Q = t_Q;
  m_A21 = t_A21;
  m_TTs = t_TTs;
  m_Tys = t_Tys;
  m_sol = t_sol;
  m_blups = t_blups;
  m_sig = t_sig;
  m_ncov = t_sol.n_elem;
  
  m_resid = t_resid;
  m_resid_G = t_resid_G;
  m_resid_GxE = t_resid_GxE;
  m_resid_E = t_resid_E;
  
  m_resid_unrelated_outliers = t_resid_unrelated_outliers;
  m_resid_unrelated_outliers_G = t_resid_unrelated_outliers_G;
  m_resid_unrelated_outliers_GxE = t_resid_unrelated_outliers_GxE;
  
  m_sum_unrelated_outliers2 = sum(t_resid_unrelated_outliers % t_resid_unrelated_outliers);
  m_sum_unrelated_outliers_G2 = sum(t_resid_unrelated_outliers_G % t_resid_unrelated_outliers_G);
  m_sum_unrelated_outliers_GxE2 = sum(t_resid_unrelated_outliers_GxE % t_resid_unrelated_outliers_GxE);
  m_sum_unrelated_outliers_G_GxE2 = 2 * sum(t_resid_unrelated_outliers_G % t_resid_unrelated_outliers_GxE);
  
  m_sum_R_nonOutlier = t_sum_R_nonOutlier;
  m_sum_R_nonOutlier_G = t_sum_R_nonOutlier_G;
  m_sum_R_nonOutlier_GxE = t_sum_R_nonOutlier_GxE;
  
  m_R_GRM_R = t_R_GRM_R;
  m_R_GRM_R_G = t_R_GRM_R_G;
  m_R_GRM_R_GxE = t_R_GRM_R_GxE;
  m_R_GRM_R_G_GxE = t_R_GRM_R_G_GxE;
  m_R_GRM_R_E = t_R_GRM_R_E;
  
  m_R_GRM_R_nonOutlier = t_R_GRM_R_nonOutlier;
  m_R_GRM_R_nonOutlier_G = t_R_GRM_R_nonOutlier_G;
  m_R_GRM_R_nonOutlier_GxE = t_R_GRM_R_nonOutlier_GxE;
  m_R_GRM_R_nonOutlier_G_GxE = t_R_GRM_R_nonOutlier_G_GxE;
  
  m_R_GRM_R_TwoSubjOutlier = t_R_GRM_R_TwoSubjOutlier;
  m_R_GRM_R_TwoSubjOutlier_G = t_R_GRM_R_TwoSubjOutlier_G;
  m_R_GRM_R_TwoSubjOutlier_GxE = t_R_GRM_R_TwoSubjOutlier_GxE;
  m_R_GRM_R_TwoSubjOutlier_G_GxE = t_R_GRM_R_TwoSubjOutlier_G_GxE;
  
  m_TwoSubj_list = t_TwoSubj_list;
  m_ThreeSubj_list = t_ThreeSubj_list;
  m_MAF_interval = t_MAF_interval;
  m_zScoreE_cutoff = t_zScoreE_cutoff;
  
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_zeta = t_zeta;
  m_tol = t_tol;
  
  m_pvalVec.resize(2);
  m_zScoreVec.resize(2);
  m_BetaVec.resize(2);
  m_seBetaVec.resize(2);
  
}
}