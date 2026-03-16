// SPACox.cpp -- SPACoxClass method implementations

#include <RcppArmadillo.h>

#include "SPACox.h"
#include "approxfun.h"

namespace SPACox {

SPACoxClass::SPACoxClass(arma::mat cumul,
                         arma::vec mresid,
                         arma::mat XinvXX,
                         arma::mat tX,
                         int N,
                         double pVal_covaAdj_Cutoff,
                         double SPA_Cutoff) {
  m_mresid = mresid;
  m_varResid = var(m_mresid);
  m_XinvXX = XinvXX;
  m_tX = tX;
  m_N = N;
  m_pVal_covaAdj_Cutoff = pVal_covaAdj_Cutoff;
  m_SPA_Cutoff = SPA_Cutoff;

  m_K_0_emp.setApproxFun(cumul.col(0), cumul.col(1));
  m_K_1_emp.setApproxFun(cumul.col(0), cumul.col(2));
  m_K_2_emp.setApproxFun(cumul.col(0), cumul.col(3));
}

double SPACoxClass::K_0(double t, int N0, double adjG0, arma::vec adjG1) {
  double sG0 = t * adjG0;
  arma::vec sG1 = t * adjG1;
  double out = N0 * m_K_0_emp.getValue(sG0) + arma::sum(m_K_0_emp.getVector(sG1));
  return out;
}

double SPACoxClass::K_1(double t, int N0, double adjG0, arma::vec adjG1, double q2) {
  double sG0 = t * adjG0;
  arma::vec sG1 = t * adjG1;
  double out = N0 * sG0 * m_K_1_emp.getValue(sG0) + arma::sum(sG1 % m_K_1_emp.getVector(sG1)) - q2;
  return out;
}

double SPACoxClass::K_2(double t, int N0, double adjG0, arma::vec adjG1) {
  double sG0 = t * adjG0;
  arma::vec sG1 = t * adjG1;
  double out = N0 * pow(sG0, 2) * m_K_2_emp.getValue(sG0) + arma::sum(pow(sG1, 2) % m_K_2_emp.getVector(sG1));
  return out;
}

SPACoxClass::RootResult SPACoxClass::fastGetRootK1(double initX, int N0, double adjG0, arma::vec adjG1, double q2) {
  double x = initX, oldX;
  double K1 = 0, K2 = 0, oldK1;
  double diffX = arma::datum::inf, oldDiffX;
  bool converge = true;
  double tol = 0.001;
  int maxiter = 100;
  int iter = 0;

  for (iter = 0; iter < maxiter; iter++){
    oldX = x;
    oldDiffX = diffX;
    oldK1 = K1;

    K1 = K_1(x, N0, adjG0, adjG1, q2);
    K2 = K_2(x, N0, adjG0, adjG1);

    diffX = -1 * K1 / K2;

    if (!std::isfinite(K1)){
      x = arma::datum::inf;
      K2 = 0;
      break;
    }

    if (arma::sign(K1) != arma::sign(oldK1)){
      while (std::abs(diffX) > std::abs(oldDiffX) - tol){
        diffX = diffX / 2;
      }
    }

    if (std::abs(diffX) < tol) break;

    x = oldX + diffX;
  }

  if (iter == maxiter)
    converge = false;

  return {x, iter, converge, K2};
}

double SPACoxClass::getProbSpa(double adjG0, arma::vec adjG1, int N0, double q2, bool lowerTail) {
  double initX = 0;

  if (q2 > 0) initX = 3;
  if (q2 <= 0) initX = -3;

  RootResult rootRes = fastGetRootK1(initX, N0, adjG0, adjG1, q2);
  double zeta = rootRes.root;

  double k1 = K_0(zeta, N0, adjG0, adjG1);
  double k2 = K_2(zeta, N0, adjG0, adjG1);
  double temp1 = zeta * q2 - k1;

  double w = arma::sign(zeta) * sqrt(2 * temp1);
  double v = zeta * sqrt(k2);

  double pval = arma::normcdf(arma::sign(lowerTail - 0.5) * (w + 1/w * log(v/w)));
  return pval;
}

double SPACoxClass::getMarkerPval(arma::vec GVec, double MAF, double& zScore) {
  double S = sum(GVec % m_mresid);
  arma::vec adjGVec = GVec - 2 * MAF;
  arma::vec adjGVec2 = pow(adjGVec, 2);
  double VarS = m_varResid * sum(adjGVec2);
  zScore = S / sqrt(VarS);

  if (std::abs(zScore) < m_SPA_Cutoff){
    double pval = arma::normcdf(-1 * std::abs(zScore)) * 2;
    return pval;
  }

  arma::uvec N1set = arma::find(GVec != 0);
  int N0 = m_N - N1set.size();

  arma::vec adjGVecNorm = adjGVec / sqrt(VarS);

  arma::vec adjG1 = adjGVecNorm.elem(N1set);
  double adjG0 = -2 * MAF / sqrt(VarS);

  double pval1 = getProbSpa(adjG0, adjG1, N0, std::abs(zScore), false);
  double pval2 = getProbSpa(adjG0, adjG1, N0, -1 * std::abs(zScore), true);
  double pval = pval1 + pval2;

  if (pval > m_pVal_covaAdj_Cutoff){
    return pval;
  }

  adjGVec = GVec - m_XinvXX * m_tX.cols(N1set) * GVec.elem(N1set);
  adjGVec2 = pow(adjGVec, 2);
  VarS = m_varResid * sum(adjGVec2);
  zScore = S / sqrt(VarS);

  adjGVecNorm = adjGVec / sqrt(VarS);

  N0 = 0;
  adjG1 = adjGVecNorm;
  adjG0 = 0;

  pval1 = getProbSpa(adjG0, adjG1, N0, std::abs(zScore), false);
  pval2 = getProbSpa(adjG0, adjG1, N0, -1 * std::abs(zScore), true);
  pval = pval1 + pval2;

  return pval;
}

void SPACoxClass::getRegionPVec(arma::vec GVec,
                                double& zScore,
                                double& pval0,
                                double& pval1,
                                arma::vec& P1Vec,
                                arma::vec& P2Vec) {
  double S = sum(GVec % m_mresid);

  arma::uvec N1set = arma::find(GVec != 0);
  arma::vec adjGVec = GVec - m_XinvXX * m_tX.cols(N1set) * GVec.elem(N1set);
  arma::vec varR_adjGVec = m_varResid * adjGVec;
  double VarS = sum(adjGVec % varR_adjGVec);
  zScore = S / sqrt(VarS);

  pval0 = arma::normcdf(-1 * std::abs(zScore)) * 2;

  arma::vec adjGVecNorm = adjGVec / sqrt(VarS);

  if (std::abs(zScore) < m_SPA_Cutoff){
    pval1 = pval0;
  }else{
    int N0 = 0;
    double adjG0 = 0;
    double pval1 = getProbSpa(adjG0, adjGVecNorm, N0, std::abs(zScore), false);
    double pval2 = getProbSpa(adjG0, adjGVecNorm, N0, -1 * std::abs(zScore), true);
    pval1 = pval1 + pval2;
  }

  P1Vec = adjGVecNorm;
  P2Vec = varR_adjGVec / sqrt(VarS);
}

}
