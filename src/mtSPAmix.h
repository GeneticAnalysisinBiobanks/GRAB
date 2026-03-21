#ifndef SPAMIX_H
#define SPAMIX_H

// SPAmix.h -- Saddlepoint approximation for admixed populations

#include <RcppArmadillo.h>
#include <boost/math/distributions/students_t.hpp>

namespace SPAmixSpace {

double getProbSpaG(
  const arma::vec MAF_outlier,
  const arma::vec residOutlier,
  double s,
  bool lower_tail,
  double mean_nonOutlier,
  double var_nonOutlier
);

arma::vec logistic_regression_beta(const arma::mat& X, const arma::vec& y);

} // namespace SPAmixSpace

class mtSPAmixClass {
public:
  // Per-phenotype outlier partition for SPA.
  struct OutlierData {
    arma::uvec posValue;
    arma::uvec posOutlier;
    arma::uvec posNonOutlier;
    arma::vec resid;
    arma::vec resid2;
    arma::vec residOutlier;
    arma::vec residNonOutlier;
    arma::vec resid2NonOutlier;
  };

private:

  const arma::mat m_resid;
  const arma::mat m_onePlusPCs;

  const int m_N;
  const int m_Npheno;
  const double m_SPA_Cutoff;

  const arma::mat m_PCs;
  const arma::vec m_sqrt_XTX_inv_diag;

  const std::vector<OutlierData> m_outlierVec;
  arma::vec m_pvalVec;
  arma::vec m_zScoreVec;

public:

  mtSPAmixClass(
    arma::mat resid,
    arma::mat PCs,
    int N,
    double SPA_Cutoff,
    std::vector<mtSPAmixClass::OutlierData> outlierVec
  );

  int getNpheno(){return m_Npheno;}
  arma::vec getpvalVec(){return m_pvalVec;}
  arma::vec getzScoreVec(){return m_zScoreVec;}

  double getMarkerPval(arma::vec GVec, double altFreq);

private:

  arma::vec fit_lm(const arma::vec& g, arma::vec& pvalues);

  arma::vec getMafEst(
    arma::vec g,
    double altFreq,
    double MAC_cutoff = 20,
    double PCs_pvalue_cutoff = 0.05,
    double MAF_est_negative_ratio_cutoff = 0.1
  );

};

#endif
