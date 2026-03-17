#ifndef SPACOX_H
#define SPACOX_H

// SPACox.h -- Saddlepoint approximation for Cox proportional hazards model

#include <RcppArmadillo.h>
#include <stdexcept>
#include "UTIL.h"


namespace SPACox{

class SPACoxClass {
private:

  approxfunClass m_K_0_emp;
  approxfunClass m_K_1_emp;
  approxfunClass m_K_2_emp;

  arma::vec m_mresid;
  double m_varResid;
  arma::mat m_XinvXX, m_tX;
  int m_N;
  double m_pVal_covaAdj_Cutoff;
  double m_SPA_Cutoff;

public:

  SPACoxClass(arma::mat cumul,
              arma::vec mresid,
              arma::mat XinvXX,
              arma::mat tX,
              int N,
              double pVal_covaAdj_Cutoff,
              double SPA_Cutoff);

  double K_0(double t, int N0, double adjG0, arma::vec adjG1);
  double K_1(double t, int N0, double adjG0, arma::vec adjG1, double q2);
  double K_2(double t, int N0, double adjG0, arma::vec adjG1);

  struct RootResult {
    double root;
    int iter;
    bool converge;
    double K2;
  };

  RootResult fastGetRootK1(double initX, int N0, double adjG0, arma::vec adjG1, double q2);
  double getProbSpa(double adjG0, arma::vec adjG1, int N0, double q2, bool lowerTail);
  double getMarkerPval(arma::vec GVec, double MAF, double& zScore);
  void getRegionPVec(arma::vec GVec, double& zScore, double& pval0, double& pval1, arma::vec& P1Vec, arma::vec& P2Vec);

};

}

#endif
