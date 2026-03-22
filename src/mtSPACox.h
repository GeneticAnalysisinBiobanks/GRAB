#ifndef SPACOX_H
#define SPACOX_H

// SPACox.h -- Saddlepoint approximation for Cox proportional hazards model

#include <RcppArmadillo.h>

// Piecewise-linear interpolation (port of R stats::approxfun)
class approxfunClass {
private:
  arma::vec m_xVec, m_yVec;
  double m_ylow, m_yhigh;
  int m_n;
  arma::vec m_slopeVec;
  void setApproxFun(arma::vec xVec, arma::vec yVec);
public:
  approxfunClass(arma::vec xVec, arma::vec yVec);
  double getValue(double v) const;
  arma::vec getVector(arma::vec vVec) const;
};

class mtSPACoxClass {
private:
  const approxfunClass m_K_0_emp;
  const approxfunClass m_K_1_emp;
  const approxfunClass m_K_2_emp;
  const arma::vec m_mresid;
  const double m_varResid;
  const arma::mat m_XinvXX, m_tX;
  const int m_N;
  const double m_pVal_covaAdj_Cutoff;
  const double m_SPA_Cutoff;

public:
  mtSPACoxClass(
    arma::mat cumul,
    arma::vec mresid,
    arma::mat XinvXX,
    arma::mat tX,
    int N,
    double pVal_covaAdj_Cutoff,
    double SPA_Cutoff
  );

  double getMarkerPval(arma::vec GVec, double MAF, double& zScore);
  void getRegionPVec(arma::vec GVec, double& zScore, double& pval0, double& pval1, arma::vec& P1Vec, arma::vec& P2Vec);

  // Returns [pval, zScore]
  std::vector<double> getResultVec(arma::vec GVec, double altFreq) {
    double zScore;
    double pval = getMarkerPval(std::move(GVec), altFreq, zScore);
    return {pval, zScore};
  }

  static int resultSize() { return 2; }

  std::string getHeaderColumns() const {
    return "\tPvalue\tzScore";
  }
};

#endif
