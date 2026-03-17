#ifndef UTIL_H
#define UTIL_H

// UTIL.h -- Shared utility functions

#include <stdexcept>
#include <limits>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/beta.hpp>
#include <RcppArmadillo.h>

double hwe_exact(int obs_hets, int obs_hom1, int obs_hom2);

void gethwepval(arma::vec GVec, double& hwepval, double hwepvalCutoff);

// Piecewise-linear interpolation (port of R stats::approxfun)
class approxfunClass {
private:

  arma::vec m_xVec, m_yVec;
  double m_ylow, m_yhigh;
  int m_n;
  arma::vec m_slopeVec;

public:

  void setApproxFun(arma::vec xVec, arma::vec yVec) {
    m_xVec = xVec;
    m_yVec = yVec;
    m_n = xVec.size();
    m_ylow = yVec(0);
    m_yhigh = yVec(m_n - 1);
    m_slopeVec.zeros(m_n - 1);

    for (int i = 0; i < m_n - 1; i ++)
      if (xVec(i+1) <= xVec(i)) throw std::runtime_error("xVec(i+1) should be greater than xVec(i).");

    for (int i = 0; i < m_n - 1; i ++)
      m_slopeVec(i) = (yVec(i+1) - yVec(i)) / (xVec(i+1) - xVec(i));
  }

  // Evaluate the interpolant at a single point via bisection lookup.
  double getValue(double v) {
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
  arma::vec getVector(arma::vec vVec) {
    int p = vVec.size();
    arma::vec outVec(p);
    for (int i = 0; i < p; i++) {
      outVec(i) = getValue(vVec(i));
    }
    return outVec;
  }
};

#endif
