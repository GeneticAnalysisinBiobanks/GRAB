
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "SPAmix.h"

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
  
  // Convert Rcpp::List to std::vector<OutlierData> (Thread-safe storage)
  m_outlierVec.resize(m_Npheno);
  for(int i = 0; i < m_Npheno; ++i) {
      Rcpp::List l = t_outlierList[i];
      m_outlierVec[i].posValue = Rcpp::as<arma::uvec>(l["posValue"]);
      m_outlierVec[i].posOutlier = Rcpp::as<arma::uvec>(l["posOutlier"]);
      m_outlierVec[i].posNonOutlier = Rcpp::as<arma::uvec>(l["posNonOutlier"]);
      m_outlierVec[i].resid = Rcpp::as<arma::vec>(l["resid"]);
      m_outlierVec[i].resid2 = Rcpp::as<arma::vec>(l["resid2"]);
      m_outlierVec[i].residOutlier = Rcpp::as<arma::vec>(l["residOutlier"]);
      m_outlierVec[i].residNonOutlier = Rcpp::as<arma::vec>(l["residNonOutlier"]);
      m_outlierVec[i].resid2NonOutlier = Rcpp::as<arma::vec>(l["resid2NonOutlier"]);
  }
  
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