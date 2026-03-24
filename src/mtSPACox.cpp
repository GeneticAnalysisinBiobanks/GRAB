// SPACox.cpp -- mtSPACoxClass method implementations

#include <RcppArmadillo.h>
#include <stdexcept>

#include "mtSPACox.h"

// ---- approxfunClass --------------------------------------------------------

approxfunClass::approxfunClass(arma::vec xVec, arma::vec yVec) {
  setApproxFun(std::move(xVec), std::move(yVec));
}

void approxfunClass::setApproxFun(arma::vec xVec, arma::vec yVec) {
  m_xVec = std::move(xVec);
  m_yVec = std::move(yVec);
  m_n = m_xVec.size();
  m_ylow = m_yVec(0);
  m_yhigh = m_yVec(m_n - 1);
  m_slopeVec.zeros(m_n - 1);

  for (int i = 0; i < m_n - 1; i++)
    if (m_xVec(i+1) <= m_xVec(i)) throw std::runtime_error("xVec(i+1) should be greater than xVec(i).");

  for (int i = 0; i < m_n - 1; i++)
    m_slopeVec(i) = (m_yVec(i+1) - m_yVec(i)) / (m_xVec(i+1) - m_xVec(i));
}

double approxfunClass::getValue(double v) const {
  int i = 0, j = m_n - 1, ij;

  if (v < m_xVec(i)) return m_ylow;
  if (v > m_xVec(j)) return m_yhigh;

  while (i < j - 1) {
    ij = (i + j) / 2;
    if (v < m_xVec(ij)) j = ij; else i = ij;
  }

  if (v == m_xVec(j)) return m_yVec(j);
  if (v == m_xVec(i)) return m_yVec(i);

  return m_yVec(i) + (v - m_xVec(i)) * m_slopeVec(i);
}

arma::vec approxfunClass::getVector(arma::vec vVec) const {
  int p = vVec.size();
  arma::vec outVec(p);
  for (int i = 0; i < p; i++)
    outVec(i) = getValue(vVec(i));
  return outVec;
}

// ---- Internal computation helpers (not part of public API) -----------------

namespace {

struct RootResult {
  double root;
  int iter;
  bool converge;
  double K2;
};

double K_0(double t, int N0, double adjG0, const arma::vec& adjG1, const approxfunClass& emp) {
  double sG0 = t * adjG0;
  arma::vec sG1 = t * adjG1;
  return N0 * emp.getValue(sG0) + arma::sum(emp.getVector(sG1));
}

double K_1(double t, int N0, double adjG0, const arma::vec& adjG1, double q2, const approxfunClass& emp) {
  double sG0 = t * adjG0;
  arma::vec sG1 = t * adjG1;
  return N0 * adjG0 * emp.getValue(sG0) + arma::sum(adjG1 % emp.getVector(sG1)) - q2;
}

double K_2(double t, int N0, double adjG0, const arma::vec& adjG1, const approxfunClass& emp) {
  double sG0 = t * adjG0;
  arma::vec sG1 = t * adjG1;
  return N0 * pow(adjG0, 2) * emp.getValue(sG0) + arma::sum(pow(adjG1, 2) % emp.getVector(sG1));
}

RootResult fastGetRootK1(
  double initX, int N0, double adjG0, const arma::vec& adjG1, double q2,
  const approxfunClass& k1_emp, const approxfunClass& k2_emp
) {
  double x = initX, oldX;
  double K1 = 0, K2 = 0, oldK1;
  double diffX = arma::datum::inf, oldDiffX;
  bool converge = true;
  const double tol = 0.001;
  const int maxiter = 100;
  int iter = 0;

  for (iter = 0; iter < maxiter; iter++) {
    oldX     = x;
    oldDiffX = diffX;
    oldK1    = K1;

    K1 = K_1(x, N0, adjG0, adjG1, q2, k1_emp);
    K2 = K_2(x, N0, adjG0, adjG1, k2_emp);

    diffX = -1 * K1 / K2;

    if (!std::isfinite(K1)) {
      x  = arma::datum::inf;
      K2 = 0;
      break;
    }

    if (arma::sign(K1) != arma::sign(oldK1)) {
      while (std::abs(diffX) > std::abs(oldDiffX) - tol)
        diffX = diffX / 2;
    }

    if (std::abs(diffX) < tol) break;

    x = oldX + diffX;
  }

  if (iter == maxiter) converge = false;

  return {x, iter, converge, K2};
}

double getProbSpa(
  double adjG0, const arma::vec& adjG1, int N0, double q2, bool lowerTail,
  const approxfunClass& k0_emp, const approxfunClass& k1_emp, const approxfunClass& k2_emp
) {
  double initX = (q2 > 0) ? 3.0 : -3.0;

  RootResult rootRes = fastGetRootK1(initX, N0, adjG0, adjG1, q2, k1_emp, k2_emp);
  double zeta = rootRes.root;

  double k1    = K_0(zeta, N0, adjG0, adjG1, k0_emp);
  double k2    = K_2(zeta, N0, adjG0, adjG1, k2_emp);
  double temp1 = zeta * q2 - k1;

  double w = arma::sign(zeta) * sqrt(2 * temp1);
  double v = zeta * sqrt(k2);

  return arma::normcdf(arma::sign(lowerTail - 0.5) * (w + 1.0 / w * log(v / w)));
}

} // namespace

// ---- mtSPACoxClass ---------------------------------------------------------

mtSPACoxClass::mtSPACoxClass(
  arma::mat cumul,
  arma::vec mresid,
  arma::mat XinvXX,
  arma::mat tX,
  int N,
  double pVal_covaAdj_Cutoff,
  double SPA_Cutoff
)
  : m_K_0_emp(cumul.col(0), cumul.col(1)),
    m_K_1_emp(cumul.col(0), cumul.col(2)),
    m_K_2_emp(cumul.col(0), cumul.col(3)),
    m_mresid(std::move(mresid)),
    m_varResid(arma::var(m_mresid)),
    m_XinvXX(std::move(XinvXX)),
    m_tX(std::move(tX)),
    m_N(N),
    m_pVal_covaAdj_Cutoff(pVal_covaAdj_Cutoff),
    m_SPA_Cutoff(SPA_Cutoff)
{}

double mtSPACoxClass::getMarkerPval(const arma::vec& GVec, double MAF, double& zScore) {
  double S = arma::sum(GVec % m_mresid);
  arma::vec adjGVec = GVec - 2 * MAF;
  double VarS = m_varResid * arma::sum(arma::pow(adjGVec, 2));
  zScore = S / sqrt(VarS);

  if (std::abs(zScore) < m_SPA_Cutoff)
    return arma::normcdf(-1 * std::abs(zScore)) * 2;

  arma::uvec N1set    = arma::find(GVec != 0);
  int        N0       = m_N - (int)N1set.size();
  arma::vec  adjGNorm = adjGVec / sqrt(VarS);
  arma::vec  adjG1    = adjGNorm.elem(N1set);
  double     adjG0    = -2 * MAF / sqrt(VarS);

  double pval = getProbSpa(adjG0, adjG1, N0,  std::abs(zScore), false, m_K_0_emp, m_K_1_emp, m_K_2_emp)
              + getProbSpa(adjG0, adjG1, N0, -std::abs(zScore), true,  m_K_0_emp, m_K_1_emp, m_K_2_emp);

  if (pval > m_pVal_covaAdj_Cutoff)
    return pval;

  adjGVec = GVec - m_XinvXX * m_tX.cols(N1set) * GVec.elem(N1set);
  VarS    = m_varResid * arma::sum(arma::pow(adjGVec, 2));
  zScore  = S / sqrt(VarS);
  adjGNorm = adjGVec / sqrt(VarS);

  pval = getProbSpa(0.0, adjGNorm, 0,  std::abs(zScore), false, m_K_0_emp, m_K_1_emp, m_K_2_emp)
       + getProbSpa(0.0, adjGNorm, 0, -std::abs(zScore), true,  m_K_0_emp, m_K_1_emp, m_K_2_emp);

  return pval;
}

void mtSPACoxClass::getRegionPVec(
  const arma::vec& GVec,
  double& zScore,
  double& pval0,
  double& pval1,
  arma::vec& P1Vec,
  arma::vec& P2Vec
) {
  double    S      = arma::sum(GVec % m_mresid);
  arma::uvec N1set = arma::find(GVec != 0);
  arma::vec adjGVec      = GVec - m_XinvXX * m_tX.cols(N1set) * GVec.elem(N1set);
  arma::vec varR_adjGVec = m_varResid * adjGVec;
  double    VarS         = arma::sum(adjGVec % varR_adjGVec);
  zScore = S / sqrt(VarS);

  pval0 = arma::normcdf(-1 * std::abs(zScore)) * 2;

  arma::vec adjGNorm = adjGVec / sqrt(VarS);

  if (std::abs(zScore) < m_SPA_Cutoff) {
    pval1 = pval0;
  } else {
    pval1 = getProbSpa(0.0, adjGNorm, 0,  std::abs(zScore), false, m_K_0_emp, m_K_1_emp, m_K_2_emp)
          + getProbSpa(0.0, adjGNorm, 0, -std::abs(zScore), true,  m_K_0_emp, m_K_1_emp, m_K_2_emp);
  }

  P1Vec = adjGNorm;
  P2Vec = varR_adjGVec / sqrt(VarS);
}

