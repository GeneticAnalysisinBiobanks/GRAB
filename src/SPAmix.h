#ifndef SPAMIX_H
#define SPAMIX_H

// SPAmix.h -- Saddlepoint approximation for admixed populations

#include <RcppArmadillo.h>
#include <boost/math/distributions/students_t.hpp>

namespace SPAmix{

class SPAmixClass {
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


  // Newton root-finder convergence result.
  struct RootResult {
      double root;
      int iter;
      bool converge;
      double K2;
  };

private:


  arma::mat m_resid;


  arma::mat m_onePlusPCs;

  int m_N;
  int m_Npheno;

  double m_SPA_Cutoff;
  arma::mat m_PCs;
  arma::vec m_sqrt_XTX_inv_diag;
  arma::vec m_diffTime1, m_diffTime2;

  std::vector<OutlierData> m_outlierVec;

  arma::vec m_pvalVec;
  arma::vec m_zScoreVec;


public:

  SPAmixClass(arma::mat resid,
              arma::mat PCs,
              int N,
              double SPA_Cutoff,
              std::vector<SPAmixClass::OutlierData> outlierVec);

  arma::vec getTestTime1(){return m_diffTime1;}
  arma::vec getTestTime2(){return m_diffTime2;}
  int getNpheno(){return m_Npheno;}
  arma::vec getpvalVec(){return m_pvalVec;}
  arma::vec getzScoreVec(){return m_zScoreVec;}

  arma::vec M_G0(arma::vec t, arma::vec MAF);
  arma::vec M_G1(arma::vec t, arma::vec MAF);
  arma::vec M_G2(arma::vec t, arma::vec MAF);
  arma::vec K_G0(arma::vec t, arma::vec MAF);
  arma::vec K_G1(arma::vec t, arma::vec MAF);
  arma::vec K_G2(arma::vec t, arma::vec MAF);
  double H_org(double t, arma::vec R, const arma::vec& MAFVec);
  double H1_adj(double t, arma::vec R, const double& s, const arma::vec& MAFVec);
  double H2(double t, arma::vec R, const arma::vec& MAFVec);
  arma::vec Horg_H2(double t, arma::vec R, const arma::vec MAFVec);
  arma::vec H1_adj_H2(double t, arma::vec R, double s, const arma::vec MAFVec);
  RootResult fastGetRootK1(double initX,
                            const double& s,
                            const arma::vec MAF_outlier,
                            double mean_nonOutlier,
                            double var_nonOutlier,
                            const arma::vec residOutlier);
  double getProbSpaG(const arma::vec MAF_outlier,
                       const arma::vec residOutlier,
                       double s,
                       bool lower_tail,
                       double mean_nonOutlier,
                       double var_nonOutlier);
  arma::vec simulate_uniform(int n, double lower, double upper);
  arma::vec fit_lm(const arma::vec& g, arma::vec& pvalues);
  arma::vec logistic_regression(const arma::mat& X, const arma::vec& y);
  arma::vec getMafEst(arma::vec g,
                      double altFreq,
                      double MAC_cutoff = 20,
                      double PCs_pvalue_cutoff = 0.05,
                      double MAF_est_negative_ratio_cutoff = 0.1);
  double getMarkerPval(arma::vec GVec, double altFreq);

};

}

#endif
