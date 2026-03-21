// POLMM.cpp -- mtPOLMMClass marker-test implementation

#include <RcppArmadillo.h>
#include "mtPOLMM.h"

// ---- Anonymous namespace: internal types and helpers ----
namespace {

struct RootResult {
  double root;
  int    iter;
  bool   converge;
  double K2;
};

struct SaddleResult {
  double    pval;
  bool      converge;
  arma::vec K1roots;
};

// Build an (n*(J-1)) x p matrix repeating each row of Cova (J-1) times
arma::mat getCovaMat(const arma::mat& Cova, int J) {
  int n = static_cast<int>(Cova.n_rows);
  int p = static_cast<int>(Cova.n_cols);
  arma::mat out(n * (J-1), p);
  int idx = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < J-1; ++j)
      out.row(idx++) = Cova.row(i);
  return out;
}

// Column-wise sum: fold n*(J-1) columns into n columns
arma::mat sumCols(const arma::mat& xMat, int J) {
  int p = static_cast<int>(xMat.n_rows);
  int n = static_cast<int>(xMat.n_cols) / (J-1);
  arma::mat out(p, n, arma::fill::zeros);
  int idx = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < J-1; ++j)
      out.col(i) += xMat.col(idx++);
  return out;
}

// Unroll n x (J-1) matrix into n*(J-1) vector (row-major)
arma::vec Mat2Vec(const arma::mat& xMat, int n, int J) {
  arma::vec out(n * (J-1));
  int idx = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < J-1; ++j)
      out(idx++) = xMat(i, j);
  return out;
}

// Roll n*(J-1) vector back into n x (J-1) matrix (row-major)
arma::mat Vec2Mat(const arma::vec& xVec, int n, int J) {
  arma::mat out(n, J-1);
  int idx = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < J-1; ++j)
      out(i, j) = xVec(idx++);
  return out;
}

// ---- SPA cumulant functions ----

double K0(double x, const arma::mat& muMat, const arma::mat& cMat, double m1) {
  arma::mat temp = -muMat + muMat % arma::exp(cMat * x);
  return arma::accu(arma::log(1 + arma::sum(temp, 1))) - m1 * x;
}

arma::vec K12(double x, const arma::mat& muMat, const arma::mat& cMat, double m1) {
  arma::mat t0 = muMat % arma::exp(cMat * x);
  arma::mat t1 = -muMat + t0;
  arma::mat t2 = t0 % cMat;
  arma::mat t3 = t2 % cMat;
  arma::vec v1 = 1 + arma::sum(t1, 1);
  arma::vec v2 = arma::sum(t2, 1);
  arma::vec v3 = arma::sum(t3, 1);
  arma::vec y(2);
  y(0) = arma::accu(v2 / v1) - m1;
  y(1) = arma::accu((v3 % v1 - arma::square(v2)) / arma::square(v1));
  return y;
}

RootResult fastgetroot_K1(double Stat, double initX, double Ratio0,
                           const arma::mat& muMat, const arma::mat& cMat, double m1) {
  double x = initX, K1 = 0, K2 = 0, diffX = arma::datum::inf;
  const double tol = 0.0001;
  int iter = 0;
  for (; iter < 100; ++iter) {
    double oldX = x, oldDiffX = diffX, oldK1 = K1;
    arma::vec kv = K12(x, muMat, cMat, m1);
    K1 = kv(0) - Stat + Ratio0 * x;
    K2 = kv(1) + Ratio0;
    diffX = -K1 / K2;
    if (!std::isfinite(K1)) { x = arma::sign(Stat) * arma::datum::inf; K2 = 0; break; }
    if (arma::sign(K1) != arma::sign(oldK1))
      while (std::abs(diffX) > std::abs(oldDiffX) - tol) diffX /= 2;
    if (std::abs(diffX) < tol) break;
    x = oldX + diffX;
  }
  return {x, iter, iter < 100, K2};
}

double fastGet_Saddle_Prob(double Stat, double zeta, double K2, double Ratio0,
                           const arma::mat& muMat, const arma::mat& cMat,
                           double m1, bool lowerTail) {
  double k1 = K0(zeta, muMat, cMat, m1) + 0.5 * zeta * zeta * Ratio0;
  if (!std::isfinite(k1) || !std::isfinite(K2)) return 0.0;
  double w = arma::sign(zeta) * std::sqrt(2.0 * (zeta * Stat - k1));
  double v = zeta * std::sqrt(K2);
  double Z = w + std::log(v / w) / w;
  return arma::normcdf(arma::sign(lowerTail - 0.5) * Z);
}

SaddleResult fastSaddle_Prob(double Stat, double VarP, double VarW, double Ratio0,
                              arma::vec K1roots,
                              const arma::vec& adjGVec1,
                              const arma::mat& muMat1,
                              const arma::mat& iRMat1) {
  int J  = static_cast<int>(muMat1.n_cols);
  int N1 = static_cast<int>(muMat1.n_rows);
  arma::mat mu1 = muMat1.cols(0, J-2);
  double adjStat  = Stat / std::sqrt(VarP);
  double sqrtVarW = std::sqrt(VarW);

  arma::mat cMat(N1, J-1);
  for (int i = 0; i < N1; ++i)
    for (int j = 0; j < J-1; ++j)
      cMat(i, j) = adjGVec1(i) / iRMat1(i, j) / sqrtVarW;

  double m1 = arma::accu(mu1 % cMat);

  RootResult r1 = fastgetroot_K1( std::abs(adjStat), std::min(K1roots(0),  5.0), Ratio0, mu1, cMat, m1);
  RootResult r2 = fastgetroot_K1(-std::abs(adjStat), std::max(K1roots(1), -5.0), Ratio0, mu1, cMat, m1);

  if (r1.converge && r2.converge) {
    double p1 = fastGet_Saddle_Prob( std::abs(adjStat), r1.root, r1.K2, Ratio0, mu1, cMat, m1, false);
    double p2 = fastGet_Saddle_Prob(-std::abs(adjStat), r2.root, r2.K2, Ratio0, mu1, cMat, m1, true);
    double pval = p1 + p2;
    if (std::isfinite(pval) && pval > 0)
      return {pval, true, {r1.root, r2.root}};
  }
  return {2.0 * arma::normcdf(-std::abs(adjStat)), false, K1roots};
}

} // anonymous namespace

// ---- POLMM class implementations ----

mtPOLMMClass::mtPOLMMClass(
  arma::mat muMat,
  arma::mat iRMat,
  arma::mat Cova,
  arma::uvec yVec,
  double varRatio,
  double SPA_Cutoff
)
  : m_muMat(std::move(muMat)),
    m_iRMat(std::move(iRMat)),
    m_n(static_cast<int>(m_muMat.n_rows)),
    m_J(static_cast<int>(m_muMat.n_cols)),
    m_p(static_cast<int>(Cova.n_cols)),
    m_varRatio(varRatio),
    m_SPA_Cutoff(SPA_Cutoff)
{
  arma::mat CovaMat = getCovaMat(Cova, m_J);

  arma::mat XR_Psi_R(m_p, m_n * (m_J-1));
  for (int k = 0; k < m_p; ++k) {
    arma::mat xMat = Vec2Mat(CovaMat.col(k), m_n, m_J);
    XR_Psi_R.row(k) = Mat2Vec(getPsixMat(xMat / m_iRMat) / m_iRMat, m_n, m_J).t();
  }
  m_XXR_Psi_RX = Cova * arma::inv(XR_Psi_R * CovaMat);
  m_XR_Psi_R   = sumCols(XR_Psi_R, m_J);

  arma::mat yMat(m_n, m_J, arma::fill::zeros);
  for (int i = 0; i < m_n; ++i) yMat(i, static_cast<int>(yVec(i))) = 1;
  arma::mat temp = yMat - m_muMat;
  arma::mat RymuMat = temp.cols(0, m_J-2) / m_iRMat;
  m_RymuVec = sumCols(RymuMat, m_J);

  arma::vec RPsiRVec(m_n, arma::fill::zeros);
  arma::mat muRMat = m_muMat.cols(0, m_J-2) / m_iRMat;
  for (int i = 0; i < m_n; ++i) {
    for (int j1 = 0; j1 < m_J-1; ++j1) {
      RPsiRVec(i) += muRMat(i, j1) / m_iRMat(i, j1) - muRMat(i, j1) * muRMat(i, j1);
      for (int j2 = j1+1; j2 < m_J-1; ++j2)
        RPsiRVec(i) -= 2.0 * muRMat(i, j1) * muRMat(i, j2);
    }
  }
  m_RPsiR = std::move(RPsiRVec);
}

arma::mat mtPOLMMClass::getPsixMat(const arma::mat& xMat) const {
  arma::mat Psi_xMat(m_n, m_J-1);
  for (int i = 0; i < m_n; ++i) {
    arma::rowvec muVec(m_J-1);
    for (int j = 0; j < m_J-1; ++j) {
      Psi_xMat(i, j) = m_muMat(i, j) * xMat(i, j);
      muVec(j) = m_muMat(i, j);
    }
    Psi_xMat.row(i) -= muVec * arma::accu(Psi_xMat.row(i));
  }
  return Psi_xMat;
}

arma::vec mtPOLMMClass::getadjGFast(const arma::vec& GVec) const {
  arma::vec XR_Psi_RG(m_p, arma::fill::zeros);
  for (int i = 0; i < m_n; ++i)
    if (GVec(i) != 0) XR_Psi_RG += m_XR_Psi_R.col(i) * GVec(i);
  return GVec - m_XXR_Psi_RX * XR_Psi_RG;
}

double mtPOLMMClass::getStatFast(const arma::vec& adjGVec) const {
  double Stat = 0;
  for (int i = 0; i < m_n; ++i)
    if (adjGVec(i) != 0) Stat += adjGVec(i) * m_RymuVec(i);
  return Stat;
}

arma::vec mtPOLMMClass::getVarWVec(const arma::vec& adjGVec) const {
  return m_RPsiR % arma::square(adjGVec);
}

void mtPOLMMClass::getMarkerPval(
  arma::vec GVec,
  double& Beta,
  double& seBeta,
  double& pval,
  double altFreq,
  double& zScore
) const {
  arma::vec adjGVec = getadjGFast(GVec);
  double    statVal = getStatFast(adjGVec);
  arma::vec VarWVec = getVarWVec(adjGVec);
  double    VarW    = arma::accu(VarWVec);
  double    VarS    = VarW * m_varRatio;
  double    StdStat = std::abs(statVal) / std::sqrt(VarS);

  pval = 2.0 * arma::normcdf(-StdStat);

  if (StdStat > m_SPA_Cutoff) {
    arma::uvec posG1 = arma::find(GVec != 0);
    double VarW1  = arma::accu(VarWVec(posG1));
    double Ratio0 = (VarW - VarW1) / VarW;
    arma::vec k1roots = {3.0, -3.0};
    SaddleResult res = fastSaddle_Prob(statVal, VarS, VarW, Ratio0, k1roots,
                                       adjGVec.elem(posG1),
                                       m_muMat.rows(posG1),
                                       m_iRMat.rows(posG1));
    pval = res.pval;
  }

  Beta   = statVal / VarS;
  seBeta = std::abs(Beta) / StdStat;
  zScore = statVal / std::sqrt(VarS);
}
