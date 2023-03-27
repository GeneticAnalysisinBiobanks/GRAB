
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "SPAmix.hpp"

namespace SPAmix {

SPAmixClass::SPAmixClass(arma::vec t_resid,
                         arma::mat t_XinvXX,
                         arma::mat t_tX,
                         arma::mat t_PCs,
                         int t_N,
                         double t_SPA_Cutoff)
{
  m_resid = t_resid;
  m_resid2 = pow(m_resid, 2);
  m_XinvXX = t_XinvXX;
  m_tX = t_tX;
  m_N = t_N;
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_PCs = t_PCs;
}
}