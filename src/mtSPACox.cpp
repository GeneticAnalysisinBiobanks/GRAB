// SPACox.cpp -- mtSPACoxClass method implementations

#include <RcppArmadillo.h>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <limits>

#include "mtSPACox.h"

// ---- File-scope helper for slope computation --------------------------------

namespace {

arma::vec computeSlopes(const arma::vec& x, const arma::vec& y) {
  arma::uword n = x.n_elem;
  arma::vec s(n - 1);
  const double* xp = x.memptr();
  const double* yp = y.memptr();
  for (arma::uword i = 0; i < n - 1; ++i)
    s[i] = (yp[i + 1] - yp[i]) / (xp[i + 1] - xp[i]);
  return s;
}

} // namespace

// ---- Constructor ------------------------------------------------------------

mtSPACoxClass::mtSPACoxClass(
  arma::mat cumul,
  arma::vec mresid,
  arma::mat XinvXX,
  arma::mat tX,
  int N,
  double pVal_covaAdj_Cutoff,
  double SPA_Cutoff
)
  : m_xGrid(cumul.col(0)),
    m_nGrid(static_cast<int>(m_xGrid.n_elem)),
    m_yK0(cumul.col(1)),
    m_slopeK0([this]() {
      // Validate x-grid is strictly increasing (once)
      const double* xp = m_xGrid.memptr();
      for (int i = 0; i < m_nGrid - 1; ++i)
        if (xp[i + 1] <= xp[i])
          throw std::runtime_error("cumul x-grid must be strictly increasing");
      return computeSlopes(m_xGrid, m_yK0);
    }()),
    m_yK1(cumul.col(2)),
    m_slopeK1(computeSlopes(m_xGrid, m_yK1)),
    m_yK2(cumul.col(3)),
    m_slopeK2(computeSlopes(m_xGrid, m_yK2)),
    m_mresid(std::move(mresid)),
    m_varResid(arma::var(m_mresid)),
    m_XinvXX(std::move(XinvXX)),
    m_tX(std::move(tX)),
    m_N(N),
    m_pVal_covaAdj_Cutoff(pVal_covaAdj_Cutoff),
    m_SPA_Cutoff(SPA_Cutoff)
{}

// ---- Interpolation ----------------------------------------------------------

double mtSPACoxClass::interp(const double* yp, const double* sp, double v) const {
  const double* xp = m_xGrid.memptr();
  const int n = m_nGrid;

  if (v < xp[0])     return yp[0];
  if (v > xp[n - 1]) return yp[n - 1];

  int i = 0, j = n - 1;
  while (i < j - 1) {
    int ij = (i + j) / 2;
    if (v < xp[ij]) j = ij; else i = ij;
  }

  if (v == xp[j]) return yp[j];
  if (v == xp[i]) return yp[i];
  return yp[i] + (v - xp[i]) * sp[i];
}

// ---- Cumulant evaluation (scalar loops, no heap alloc) ----------------------

double mtSPACoxClass::evalK0(
  double t, int N0, double adjG0,
  const double* adjG, const arma::uword* idx, int n
) const {
  double sum = N0 * interpK0(t * adjG0);
  if (idx) {
    for (int k = 0; k < n; ++k) sum += interpK0(t * adjG[idx[k]]);
  } else {
    for (int k = 0; k < n; ++k) sum += interpK0(t * adjG[k]);
  }
  return sum;
}

double mtSPACoxClass::evalK1(
  double t, int N0, double adjG0,
  const double* adjG, const arma::uword* idx, int n, double q2
) const {
  double sum = N0 * adjG0 * interpK1(t * adjG0);
  if (idx) {
    for (int k = 0; k < n; ++k) {
      double a = adjG[idx[k]];
      sum += a * interpK1(t * a);
    }
  } else {
    for (int k = 0; k < n; ++k) {
      double a = adjG[k];
      sum += a * interpK1(t * a);
    }
  }
  return sum - q2;
}

double mtSPACoxClass::evalK2(
  double t, int N0, double adjG0,
  const double* adjG, const arma::uword* idx, int n
) const {
  double sum = N0 * adjG0 * adjG0 * interpK2(t * adjG0);
  if (idx) {
    for (int k = 0; k < n; ++k) {
      double a = adjG[idx[k]];
      sum += a * a * interpK2(t * a);
    }
  } else {
    for (int k = 0; k < n; ++k) {
      double a = adjG[k];
      sum += a * a * interpK2(t * a);
    }
  }
  return sum;
}

// ---- SPA root-finding -------------------------------------------------------

mtSPACoxClass::RootResult mtSPACoxClass::fastGetRootK1(
  double initX, int N0, double adjG0,
  const double* adjG, const arma::uword* idx, int n, double q2
) const {
  double x = initX, oldX;
  double K1val = 0.0, K2val = 0.0, oldK1;
  double diffX = std::numeric_limits<double>::infinity(), oldDiffX;
  bool converge = true;
  const double tol = 0.001;
  const int maxiter = 100;

  for (int iter = 0; iter < maxiter; ++iter) {
    oldX     = x;
    oldDiffX = diffX;
    oldK1    = K1val;

    K1val = evalK1(x, N0, adjG0, adjG, idx, n, q2);
    K2val = evalK2(x, N0, adjG0, adjG, idx, n);

    diffX = -K1val / K2val;

    if (!std::isfinite(K1val)) {
      x    = std::numeric_limits<double>::infinity();
      K2val = 0.0;
      break;
    }

    if ((K1val > 0) != (oldK1 > 0)) {
      while (std::abs(diffX) > std::abs(oldDiffX) - tol)
        diffX *= 0.5;
    }

    if (std::abs(diffX) < tol) break;

    x = oldX + diffX;

    if (iter == maxiter - 1) converge = false;
  }

  return {x, converge, K2val};
}

// ---- SPA probability --------------------------------------------------------

double mtSPACoxClass::getProbSpa(
  double adjG0, const double* adjG, const arma::uword* idx, int n,
  int N0, double q2, bool lowerTail
) const {
  double initX = (q2 > 0) ? 3.0 : -3.0;

  RootResult rootRes = fastGetRootK1(initX, N0, adjG0, adjG, idx, n, q2);
  double zeta = rootRes.root;

  double k0val = evalK0(zeta, N0, adjG0, adjG, idx, n);
  double k2val = evalK2(zeta, N0, adjG0, adjG, idx, n);
  double temp1 = zeta * q2 - k0val;

  if (!std::isfinite(zeta) || temp1 < 0.0 || k2val <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();

  double w = std::copysign(std::sqrt(2.0 * temp1), zeta);
  double v = zeta * std::sqrt(k2val);

  if (w == 0.0 || v == 0.0 || (v / w) <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();

  double sign = lowerTail ? 1.0 : -1.0;
  return arma::normcdf(sign * (w + 1.0 / w * std::log(v / w)));
}

// ---- getMarkerPval ----------------------------------------------------------

double mtSPACoxClass::getMarkerPval(const arma::vec& GVec, double MAF, double& zScore) {
  double S = arma::dot(GVec, m_mresid);
  double twoMAF = 2.0 * MAF;
  const double* gp = GVec.memptr();

  // Compute VarS = m_varResid * sum((g - 2*MAF)^2), no temp alloc
  double sumAdjG2 = 0.0;
  for (int i = 0; i < m_N; ++i) {
    double adj = gp[i] - twoMAF;
    sumAdjG2 += adj * adj;
  }
  double VarS = m_varResid * sumAdjG2;
  zScore = S / std::sqrt(VarS);

  if (std::abs(zScore) < m_SPA_Cutoff)
    return arma::normcdf(-std::abs(zScore)) * 2.0;

  // Build N1set (non-zero indices) and adjGNorm in one pass
  double sqrtVarS = std::sqrt(VarS);
  double adjG0 = -twoMAF / sqrtVarS;

  arma::vec adjGNorm(m_N);
  double* anp = adjGNorm.memptr();
  std::vector<arma::uword> N1set;
  N1set.reserve(m_N);
  for (int i = 0; i < m_N; ++i) {
    anp[i] = (gp[i] - twoMAF) / sqrtVarS;
    if (gp[i] != 0.0) N1set.push_back(static_cast<arma::uword>(i));
  }
  int nN1 = static_cast<int>(N1set.size());
  int N0  = m_N - nN1;

  // First SPA (indexed access — only non-zero subjects)
  double absZ = std::abs(zScore);
  double pval = getProbSpa(adjG0, anp, N1set.data(), nN1, N0,  absZ, false)
              + getProbSpa(adjG0, anp, N1set.data(), nN1, N0, -absZ, true);

  if (pval > m_pVal_covaAdj_Cutoff) return pval;

  // Covariate adjustment: adjGVec = GVec - XinvXX * (tX * g)_nonzero
  // Column accumulation avoids .cols() and .elem() heap allocs
  int nCov = static_cast<int>(m_tX.n_rows);
  arma::vec tX_g(nCov, arma::fill::zeros);
  double* tp = tX_g.memptr();
  for (int k = 0; k < nN1; ++k) {
    arma::uword j = N1set[k];
    const double* col = m_tX.colptr(j);
    double gj = gp[j];
    for (int r = 0; r < nCov; ++r)
      tp[r] += col[r] * gj;
  }
  arma::vec adjGVec = GVec - m_XinvXX * tX_g;

  VarS = m_varResid * arma::dot(adjGVec, adjGVec);
  zScore = S / std::sqrt(VarS);
  sqrtVarS = std::sqrt(VarS);

  // Recompute adjGNorm for full vector
  const double* avp = adjGVec.memptr();
  for (int i = 0; i < m_N; ++i)
    anp[i] = avp[i] / sqrtVarS;

  // Second SPA (full vector, N0=0)
  absZ = std::abs(zScore);
  pval = getProbSpa(0.0, anp, nullptr, m_N, 0,  absZ, false)
       + getProbSpa(0.0, anp, nullptr, m_N, 0, -absZ, true);

  return pval;
}

// ---- getRegionPVec ----------------------------------------------------------

void mtSPACoxClass::getRegionPVec(
  const arma::vec& GVec,
  double& zScore,
  double& pval0,
  double& pval1,
  arma::vec& P1Vec,
  arma::vec& P2Vec
) {
  double S = arma::dot(GVec, m_mresid);
  const double* gp = GVec.memptr();

  // Covariate adjustment via column accumulation (no .cols()/.elem())
  int nCov = static_cast<int>(m_tX.n_rows);
  arma::vec tX_g(nCov, arma::fill::zeros);
  double* tp = tX_g.memptr();
  for (int i = 0; i < m_N; ++i) {
    if (gp[i] == 0.0) continue;
    const double* col = m_tX.colptr(i);
    double gi = gp[i];
    for (int r = 0; r < nCov; ++r)
      tp[r] += col[r] * gi;
  }
  arma::vec adjGVec = GVec - m_XinvXX * tX_g;

  double VarS = m_varResid * arma::dot(adjGVec, adjGVec);
  zScore = S / std::sqrt(VarS);
  pval0 = arma::normcdf(-std::abs(zScore)) * 2.0;

  double sqrtVarS = std::sqrt(VarS);
  arma::vec adjGNorm = adjGVec / sqrtVarS;

  if (std::abs(zScore) < m_SPA_Cutoff) {
    pval1 = pval0;
  } else {
    const double* anp = adjGNorm.memptr();
    double absZ = std::abs(zScore);
    pval1 = getProbSpa(0.0, anp, nullptr, m_N, 0,  absZ, false)
          + getProbSpa(0.0, anp, nullptr, m_N, 0, -absZ, true);
  }

  P1Vec = adjGNorm;
  P2Vec = (m_varResid / sqrtVarS) * adjGVec;
}

