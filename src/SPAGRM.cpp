
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
                         std::vector<arma::vec> t_TwoSubj_resid,
                         std::vector<arma::vec> t_TwoSubj_rho,
                         std::vector<arma::vec> t_ThreeSubj_standS,
                         std::vector<arma::mat> t_ThreeSubj_CLT,
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

  m_TwoSubj_resid_list = std::move(t_TwoSubj_resid);
  m_TwoSubj_rho_list   = std::move(t_TwoSubj_rho);
  m_ThreeSubj_standS_list = std::move(t_ThreeSubj_standS);
  m_ThreeSubj_CLT_list    = std::move(t_ThreeSubj_CLT);

  m_SPA_Cutoff = t_SPA_Cutoff;
  m_zeta = t_zeta;
  m_tol = t_tol;
}
}