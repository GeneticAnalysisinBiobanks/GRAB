
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "SPAGRM.hpp"

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
  m_TwoSubj_list = t_TwoSubj_list;
  m_ThreeSubj_list = t_ThreeSubj_list;
  
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_zeta = t_zeta;
  m_tol = t_tol;
}
}