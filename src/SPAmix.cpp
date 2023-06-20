
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "SPAmix.hpp"

namespace SPAmix {

SPAmixClass::SPAmixClass(arma::mat t_resid,
                         arma::mat t_PCs,
                         int t_N,
                         double t_SPA_Cutoff,
                         Rcpp::List t_outlierList)
                         // arma::uvec t_posOutlier,
                         // arma::uvec t_posNonOutlier)
{
  m_resid = t_resid;
  // m_resid2 = pow(m_resid, 2);
  m_Npheno = t_resid.n_cols;
  
  m_pvalVec.resize(m_Npheno);
  m_zScoreVec.resize(m_Npheno);
  
  m_N = t_N;
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_PCs = t_PCs;
  m_outlierList = t_outlierList;
  
  // m_posOutlier = t_posOutlier;
  // m_posNonOutlier = t_posNonOutlier;
  
  // m_R_outlier = m_resid(m_posOutlier);  // what if no residuals are outlier?? check later (2023-04-23)
  // m_R_nonOutlier = m_resid(m_posNonOutlier);  
  
  arma::mat X = arma::join_horiz(arma::ones(t_N), t_PCs);
  arma::mat X_t = X.t();
  arma::mat XTX = X_t * X;
  arma::mat XTX_inv = arma::inv(XTX);
  arma::vec XTX_inv_diag = XTX_inv.diag();
  
  m_onePlusPCs = X;
  m_sqrt_XTX_inv_diag = arma::sqrt(XTX_inv_diag);
  
  m_diffTime1.zeros(2);
  m_diffTime2.zeros(2);
}
}