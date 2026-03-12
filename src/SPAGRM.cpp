
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "SPAGRM.h"

namespace SPAGRM {

SPAGRMClass::SPAGRMClass(arma::vec t_resid,
                         arma::vec t_resid_unrelated_outliers,
                         double t_sum_R_nonOutlier,
                         double t_R_GRM_R_nonOutlier,
                         double t_R_GRM_R_TwoSubjOutlier,
                         double t_R_GRM_R,
                         arma::vec t_MAF_interval,
                         Rcpp::List t_TwoSubj_list,
                         Rcpp::List t_ThreeSubj_list,
                         double t_SPA_Cutoff,
                         double t_zeta,
                         double t_tol)
{
  m_resid = t_resid;
  m_resid_unrelated_outliers = t_resid_unrelated_outliers;
  m_sum_unrelated_outliers2 = sum(t_resid_unrelated_outliers % t_resid_unrelated_outliers);
  m_sum_R_nonOutlier = t_sum_R_nonOutlier;
  m_R_GRM_R_nonOutlier = t_R_GRM_R_nonOutlier;
  m_R_GRM_R_TwoSubjOutlier = t_R_GRM_R_TwoSubjOutlier;
  m_R_GRM_R = t_R_GRM_R;
  m_MAF_interval = t_MAF_interval;

  m_TwoSubj_resid_list.clear();
  m_TwoSubj_rho_list.clear();
  m_TwoSubj_resid_list.reserve(t_TwoSubj_list.size());
  m_TwoSubj_rho_list.reserve(t_TwoSubj_list.size());
  for (int i = 0; i < t_TwoSubj_list.size(); i++) {
    Rcpp::List two_i = Rcpp::as<Rcpp::List>(t_TwoSubj_list[i]);
    m_TwoSubj_resid_list.push_back(Rcpp::as<arma::vec>(two_i["Resid"]));
    m_TwoSubj_rho_list.push_back(Rcpp::as<arma::vec>(two_i["Rho"]));
  }

  m_ThreeSubj_standS_list.clear();
  m_ThreeSubj_CLT_list.clear();
  m_ThreeSubj_standS_list.reserve(t_ThreeSubj_list.size());
  m_ThreeSubj_CLT_list.reserve(t_ThreeSubj_list.size());
  for (int i = 0; i < t_ThreeSubj_list.size(); i++) {
    Rcpp::List three_i = Rcpp::as<Rcpp::List>(t_ThreeSubj_list[i]);
    m_ThreeSubj_standS_list.push_back(Rcpp::as<arma::vec>(three_i["stand.S"]));
    m_ThreeSubj_CLT_list.push_back(Rcpp::as<arma::mat>(three_i["CLT"]));
  }
  
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_zeta = t_zeta;
  m_tol = t_tol;
}
}