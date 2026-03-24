#ifndef SPAMIX_H
#define SPAMIX_H

// SPAmix.h -- Saddlepoint approximation for admixed populations

#include <RcppArmadillo.h>
#include <boost/math/distributions/students_t.hpp>

namespace nsSPAmix {

double getProbSpaG(
  const arma::vec& MAF_outlier,
  const arma::vec& residOutlier,
  double s,
  bool lower_tail,
  double mean_nonOutlier,
  double var_nonOutlier
);

arma::vec logistic_regression_beta(const arma::mat& X, const arma::vec& y);

} // namespace nsSPAmix

class mtSPAmixClass {
public:
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

  const arma::vec m_resid;
  const arma::mat m_onePlusPCs;

  const int m_N;
  const double m_SPA_Cutoff;

  const arma::mat m_PCs;
  const arma::vec m_sqrt_XTX_inv_diag;

  const OutlierData m_outlier;
  double m_pval;
  double m_zScore;

public:

  mtSPAmixClass(
    arma::vec resid,
    arma::mat PCs,
    int N,
    double SPA_Cutoff,
    OutlierData outlier
  );

  void getResultVec(const arma::vec& GVec, double altFreq, std::vector<double>& rv) {
    getMarkerPval(GVec, altFreq);
    rv.clear();
    rv.push_back(m_pval);
    rv.push_back(m_zScore);
  }

  int resultSize() const { return 2; }

  std::string getHeaderColumns() const {
    return "\tPvalue\tzScore";
  }

  double getMarkerPval(const arma::vec& GVec, double altFreq);

private:

  arma::vec fit_lm(const arma::vec& g, arma::vec& pvalues);

  arma::vec getMafEst(
    const arma::vec& g,
    double altFreq,
    double MAC_cutoff = 20,
    double PCs_pvalue_cutoff = 0.05,
    double MAF_est_negative_ratio_cutoff = 0.1
  );

};

#endif
