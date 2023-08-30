
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "WtSPAG.hpp"

namespace WtSPAG {

WtSPAGClass::WtSPAGClass(arma::mat t_mresid,
                         arma::vec t_weight,
                         int t_N,
                         double t_SPA_Cutoff)
{
  m_mresid = t_mresid;
  m_weight = t_weight;
  m_N = t_N;
  m_SPA_Cutoff = t_SPA_Cutoff;
}
}
