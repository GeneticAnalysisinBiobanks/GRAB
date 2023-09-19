
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "WtSPAG.hpp"

namespace WtSPAG {

WtSPAGClass::WtSPAGClass(arma::mat t_mresid,
                         // arma::vec t_weight,
                         int t_N,
                         double t_SPA_Cutoff,
                         Rcpp::List t_outlierList)
{
  m_mresid = t_mresid;
  // m_weight = t_weight;
  m_N = t_N;
  m_SPA_Cutoff = t_SPA_Cutoff;
  
  arma::uvec posOutlier = t_outlierList["posOutlier"];
  arma::uvec posNonOutlier = t_outlierList["posNonOutlier"];
 // m_N_NonOutlier = m_posNonOutlier.size();
  
  // m_resid is the same as m_mresid, might be updated later 2023-08-31
  arma::vec resid = t_outlierList["resid"];
  arma::vec resid2 = t_outlierList["resid2"];
  arma::vec residOutlier = t_outlierList["residOutlier"];
  arma::vec residNonOutlier = t_outlierList["residNonOutlier"];
  arma::vec resid2NonOutlier = t_outlierList["resid2NonOutlier"];
  
  m_posOutlier = posOutlier;
  m_posNonOutlier = posNonOutlier;
  m_N_NonOutlier = m_posNonOutlier.size();
  
  m_resid = resid;
  m_resid2 = resid2;
  m_residOutlier = residOutlier;
  m_residNonOutlier = residNonOutlier;
  m_resid2NonOutlier = resid2NonOutlier;
  
  m_sum_resid2 = sum(m_resid2);
  m_sum_resid = sum(m_resid);
  m_sum_resid2NonOutlier = sum(m_resid2NonOutlier);
  m_sum_residNonOutlier = sum(m_residNonOutlier);
  
}
}
