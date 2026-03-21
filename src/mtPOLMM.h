#ifndef POLMM_H
#define POLMM_H

// POLMM.h -- Proportional-odds logistic mixed model for ordinal categorical traits

#include <RcppArmadillo.h>

class mtPOLMMClass {
private:
  const arma::mat m_muMat;
  const arma::mat m_iRMat;
  const int       m_n, m_J, m_p;
  const double    m_varRatio;
  const double    m_SPA_Cutoff;

  arma::mat m_XXR_Psi_RX;
  arma::mat m_XR_Psi_R;
  arma::vec m_RymuVec;
  arma::vec m_RPsiR;

  arma::mat getPsixMat(const arma::mat& xMat) const;
  arma::vec getadjGFast(const arma::vec& GVec) const;
  double    getStatFast(const arma::vec& adjGVec) const;
  arma::vec getVarWVec(const arma::vec& adjGVec) const;

public:
  mtPOLMMClass(
    arma::mat muMat,
    arma::mat iRMat,
    arma::mat Cova,
    arma::uvec yVec,
    double varRatio,
    double SPA_Cutoff
  );

  void getMarkerPval(
    arma::vec GVec,
    double& Beta,
    double& seBeta,
    double& pval,
    double altFreq,
    double& zScore
  ) const;
};

#endif
