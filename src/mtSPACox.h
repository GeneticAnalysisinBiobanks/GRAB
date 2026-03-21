#ifndef SPACOX_H
#define SPACOX_H

// SPACox.h -- Saddlepoint approximation for Cox proportional hazards model

#include <RcppArmadillo.h>
#include <stdexcept>


namespace SPACox{

// Piecewise-linear interpolation (port of R stats::approxfun)
class approxfunClass {
  
private:

  arma::vec m_xVec, m_yVec;
  double m_ylow, m_yhigh;
  int m_n;
  arma::vec m_slopeVec;

public:

  approxfunClass() = default;

  approxfunClass(arma::vec xVec, arma::vec yVec) {
    setApproxFun(std::move(xVec), std::move(yVec));
  }

  void setApproxFun(arma::vec xVec, arma::vec yVec) {
    m_xVec = std::move(xVec);
    m_yVec = std::move(yVec);
    m_n = m_xVec.size();
    m_ylow = m_yVec(0);
    m_yhigh = m_yVec(m_n - 1);
    m_slopeVec.zeros(m_n - 1);

    for (int i = 0; i < m_n - 1; i ++)
      if (m_xVec(i+1) <= m_xVec(i)) throw std::runtime_error("xVec(i+1) should be greater than xVec(i).");

    for (int i = 0; i < m_n - 1; i ++)
      m_slopeVec(i) = (m_yVec(i+1) - m_yVec(i)) / (m_xVec(i+1) - m_xVec(i));
  }

  // Evaluate the interpolant at a single point via bisection lookup.
  double getValue(double v) const {
    int i, j, ij;
    i = 0;
    j = m_n - 1;

    if (v < m_xVec(i)) return m_ylow;
    if (v > m_xVec(j)) return m_yhigh;

    while (i < j - 1) {
      ij = (i + j)/2;
      if (v < m_xVec(ij)) j = ij; else i = ij;
    }

    if (v == m_xVec(j)) return m_yVec(j);
    if (v == m_xVec(i)) return m_yVec(i);

    return m_yVec(i) + (v - m_xVec(i)) * m_slopeVec(i);
  }

  // Evaluate the interpolant at each element of a vector.
  arma::vec getVector(arma::vec vVec) const {
    int p = vVec.size();
    arma::vec outVec(p);
    for (int i = 0; i < p; i++) {
      outVec(i) = getValue(vVec(i));
    }
    return outVec;
  }
};

class SPACoxClass {
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

  SPACoxClass(
    arma::mat cumul,
    arma::vec mresid,
    arma::mat XinvXX,
    arma::mat tX,
    int N,
    double pVal_covaAdj_Cutoff,
    double SPA_Cutoff
  );

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
