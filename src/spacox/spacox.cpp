// spacox.cpp — SPACox full implementation

#include "spacox/spacox.hpp"
#include "io/plink.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ======================================================================
// DesignMatrix
// ======================================================================

DesignMatrix::DesignMatrix(const Eigen::MatrixXd& X)
  : m_X(X)
{
  const int nCols = static_cast<int>(m_X.cols());
  m_tX = m_X.transpose();
  Eigen::MatrixXd XtX = m_tX * m_X;
  m_XinvXX = m_X * XtX.ldlt().solve(
      Eigen::MatrixXd::Identity(nCols, nCols));
}

void DesignMatrix::adjustGenotype(
    const double* G, const uint32_t* nzIdx, int nNz,
    Eigen::Ref<Eigen::VectorXd> adjG) const {

  const int p = nCols();
  const int N = nRows();

  Eigen::VectorXd tX_g = Eigen::VectorXd::Zero(p);
  for (int k = 0; k < nNz; ++k) {
    uint32_t j = nzIdx[k];
    double gj = G[j];
    tX_g.noalias() += gj * m_tX.col(j);
  }

  Eigen::Map<const Eigen::VectorXd> gVec(G, N);
  adjG.noalias() = gVec - m_XinvXX * tX_g;
}

// ======================================================================
// Build empirical CGF interpolation table
// ======================================================================

CumulantTable buildCumulantTable(const Eigen::VectorXd& residuals) {
  constexpr double rangeMax = 100.0;
  constexpr int    L        = 10000;

  // Cauchy-quantile spacing: idx0 = tan(pi*(k/(L+1) - 0.5))
  // then rescale so max(idx1) = rangeMax.
  Eigen::VectorXd xGrid(L);
  {
    double maxAbs = 0.0;
    for (int k = 0; k < L; ++k) {
      double p = static_cast<double>(k + 1) / static_cast<double>(L + 1);
      xGrid[k] = std::tan(M_PI * (p - 0.5));
      double a = std::fabs(xGrid[k]);
      if (a > maxAbs) maxAbs = a;
    }
    double scale = rangeMax / maxAbs;
    xGrid *= scale;
  }

  const int N = static_cast<int>(residuals.size());
  const double* rp = residuals.data();

  Eigen::VectorXd yK0(L), yK1(L), yK2(L);

  for (int i = 0; i < L; ++i) {
    double t = xGrid[i];
    double M0 = 0.0, M1 = 0.0, M2 = 0.0;
    for (int j = 0; j < N; ++j) {
      double r = rp[j];
      double e = std::exp(r * t);
      M0 += e;
      double re = r * e;
      M1 += re;
      M2 += r * re;
    }
    double invN = 1.0 / static_cast<double>(N);
    M0 *= invN;
    M1 *= invN;
    M2 *= invN;
    yK0[i] = std::log(M0);
    yK1[i] = M1 / M0;
    yK2[i] = (M0 * M2 - M1 * M1) / (M0 * M0);
  }

  // Pre-compute slopes for piecewise-linear interpolation.
  auto computeSlopes = [](const Eigen::VectorXd& x,
                          const Eigen::VectorXd& y) -> Eigen::VectorXd {
    const int n = static_cast<int>(x.size());
    Eigen::VectorXd s(n - 1);
    for (int i = 0; i < n - 1; ++i)
      s[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
    return s;
  };

  CumulantTable ct;
  ct.xGrid  = std::move(xGrid);
  ct.nGrid  = L;
  ct.yK0    = std::move(yK0);
  ct.slopeK0 = computeSlopes(ct.xGrid, ct.yK0);
  ct.yK1    = std::move(yK1);
  ct.slopeK1 = computeSlopes(ct.xGrid, ct.yK1);
  ct.yK2    = std::move(yK2);
  ct.slopeK2 = computeSlopes(ct.xGrid, ct.yK2);
  return ct;
}

// ======================================================================
// SPACoxMethod — construction & clone
// ======================================================================

SPACoxMethod::SPACoxMethod(
    const Eigen::VectorXd& residuals,
    double varResid,
    const CumulantTable& cumul,
    const DesignMatrix& design,
    double pvalCovAdjCut,
    double spaCutoff)
  : m_resid(residuals),
    m_varResid(varResid),
    m_cumul(cumul),
    m_design(design),
    m_N(static_cast<int>(residuals.size())),
    m_pvalCovAdjCut(pvalCovAdjCut),
    m_spaCutoff(spaCutoff),
    m_adjGNorm(residuals.size()),
    m_adjGVec(residuals.size())
{
  m_nzSet.reserve(m_N);
}

std::unique_ptr<MethodBase> SPACoxMethod::clone() const {
  // Shared const refs stay the same; scratch buffers are freshly allocated.
  return std::make_unique<SPACoxMethod>(
      m_resid, m_varResid, m_cumul, m_design,
      m_pvalCovAdjCut, m_spaCutoff);
}

// ======================================================================
// CGF interpolation (binary-search piecewise linear)
// ======================================================================

double SPACoxMethod::interp(const double* yp, const double* sp, double v) const {
  const double* xp = m_cumul.xGrid.data();
  const int n = m_cumul.nGrid;

  if (v <= xp[0])     return yp[0];
  if (v >= xp[n - 1]) return yp[n - 1];

  int lo = 0, hi = n - 1;
  while (lo < hi - 1) {
    int mid = (lo + hi) >> 1;
    if (v < xp[mid]) hi = mid; else lo = mid;
  }
  return yp[lo] + (v - xp[lo]) * sp[lo];
}

double SPACoxMethod::interpK0(double v) const {
  return interp(m_cumul.yK0.data(), m_cumul.slopeK0.data(), v);
}
double SPACoxMethod::interpK1(double v) const {
  return interp(m_cumul.yK1.data(), m_cumul.slopeK1.data(), v);
}
double SPACoxMethod::interpK2(double v) const {
  return interp(m_cumul.yK2.data(), m_cumul.slopeK2.data(), v);
}

// ======================================================================
// Cumulant evaluation — scalar loops, no heap alloc
// ======================================================================

double SPACoxMethod::evalK0(
    double t, int N0, double adjG0,
    const double* adjG, const uint32_t* idx, int n) const {
  double sum = N0 * interpK0(t * adjG0);
  if (idx) {
    for (int k = 0; k < n; ++k) sum += interpK0(t * adjG[idx[k]]);
  } else {
    for (int k = 0; k < n; ++k) sum += interpK0(t * adjG[k]);
  }
  return sum;
}

double SPACoxMethod::evalK1(
    double t, int N0, double adjG0,
    const double* adjG, const uint32_t* idx, int n, double q2) const {
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

double SPACoxMethod::evalK2(
    double t, int N0, double adjG0,
    const double* adjG, const uint32_t* idx, int n) const {
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

// ======================================================================
// Newton-Raphson root-finding on K1(t) - q2 = 0
// ======================================================================

SPACoxMethod::RootResult SPACoxMethod::fastGetRootK1(
    double initX, int N0, double adjG0,
    const double* adjG, const uint32_t* idx, int n, double q2) const {

  double x = initX, oldX;
  double K1val = 0.0, K2val = 0.0, oldK1;
  double diffX = std::numeric_limits<double>::infinity(), oldDiffX;
  bool converge = true;
  constexpr double tol = 0.001;
  constexpr int maxiter = 100;

  for (int iter = 0; iter < maxiter; ++iter) {
    oldX     = x;
    oldDiffX = diffX;
    oldK1    = K1val;

    K1val = evalK1(x, N0, adjG0, adjG, idx, n, q2);
    K2val = evalK2(x, N0, adjG0, adjG, idx, n);

    diffX = -K1val / K2val;

    if (!std::isfinite(K1val)) {
      x     = std::numeric_limits<double>::infinity();
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

// ======================================================================
// SPA tail probability (Lugannani-Rice)
// ======================================================================

double SPACoxMethod::getProbSpa(
    double adjG0, const double* adjG, const uint32_t* idx, int n,
    int N0, double q2, bool lowerTail) const {

  double initX = (q2 > 0) ? 3.0 : -3.0;

  RootResult rr = fastGetRootK1(initX, N0, adjG0, adjG, idx, n, q2);
  double zeta = rr.root;

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
  return math::pnorm(sign * (w + (1.0 / w) * std::log(v / w)));
}

// ======================================================================
// Per-marker p-value (two-stage SPA)
// ======================================================================

double SPACoxMethod::getMarkerPval(
    const Eigen::Ref<Eigen::VectorXd>& GVec,
    double MAF, double& zScore) {

  // Score statistic S = G' * resid
  double S = GVec.dot(m_resid);
  double twoMAF = 2.0 * MAF;
  const double* gp = GVec.data();

  // VarS = varResid * Σ(g_i - 2·MAF)²
  double sumAdjG2 = 0.0;
  for (int i = 0; i < m_N; ++i) {
    double adj = gp[i] - twoMAF;
    sumAdjG2 += adj * adj;
  }
  double VarS = m_varResid * sumAdjG2;
  if (VarS <= 0.0) {
    zScore = 0.0;
    return std::numeric_limits<double>::quiet_NaN();
  }
  zScore = S / std::sqrt(VarS);

  // Normal approximation if |z| below SPA cutoff
  if (std::abs(zScore) < m_spaCutoff)
    return 2.0 * math::pnorm(-std::abs(zScore));

  // ---- Stage 1: unadjusted SPA with non-zero index optimisation ----
  double sqrtVarS = std::sqrt(VarS);
  double adjG0 = -twoMAF / sqrtVarS;

  double* anp = m_adjGNorm.data();
  m_nzSet.clear();
  for (int i = 0; i < m_N; ++i) {
    anp[i] = (gp[i] - twoMAF) / sqrtVarS;
    if (gp[i] != 0.0)
      m_nzSet.push_back(static_cast<uint32_t>(i));
  }
  int nNz = static_cast<int>(m_nzSet.size());
  int N0  = m_N - nNz;

  double absZ = std::abs(zScore);
  double pval = getProbSpa(adjG0, anp, m_nzSet.data(), nNz, N0,  absZ, false)
              + getProbSpa(adjG0, anp, m_nzSet.data(), nNz, N0, -absZ, true);

  if (pval > m_pvalCovAdjCut) return pval;

  // ---- Stage 2: covariate-adjusted SPA ----
  m_design.adjustGenotype(gp, m_nzSet.data(), nNz, m_adjGVec);

  VarS = m_varResid * m_adjGVec.squaredNorm();
  if (VarS <= 0.0) {
    zScore = 0.0;
    return std::numeric_limits<double>::quiet_NaN();
  }
  zScore   = S / std::sqrt(VarS);
  sqrtVarS = std::sqrt(VarS);

  const double* avp = m_adjGVec.data();
  for (int i = 0; i < m_N; ++i)
    anp[i] = avp[i] / sqrtVarS;

  // Full-vector SPA (no sparsity shortcut after adjustment)
  absZ = std::abs(zScore);
  pval = getProbSpa(0.0, anp, nullptr, m_N, 0,  absZ, false)
       + getProbSpa(0.0, anp, nullptr, m_N, 0, -absZ, true);

  return pval;
}

// ======================================================================
// MethodBase::getResultVec — called by marker engine per marker
// ======================================================================

void SPACoxMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq, int /*markerInChunkIdx*/,
    bool flipped, std::vector<double>& result) {

  // altFreq is the alternate-allele frequency BEFORE any flip.
  // If the engine flipped GVec (altFreq > 0.5), MAF = 1 - altFreq.
  double MAF = flipped ? (1.0 - altFreq) : altFreq;

  double zScore;
  double pval = getMarkerPval(GVec, MAF, zScore);

  result.push_back(pval);
  result.push_back(zScore);
}

// ======================================================================
// runSPACox — orchestration
// ======================================================================

void runSPACox(
    const std::string& residFile,
    const std::vector<std::string>& covarNames,
    const std::string& phenoFile,
    const std::string& covarFile,
    const std::string& bfilePrefix,
    const std::string& outputFile,
    double pvalCovAdjCut,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff) {

  // ---- Load resid file and covariate data ----
  infoMsg("Loading resid file: %s", residFile.c_str());
  auto famIIDs = parseFamIIDs(bfilePrefix + ".fam");
  SubjectData sd(std::move(famIIDs));
  sd.loadResidOne(residFile);
  if (!phenoFile.empty())
    sd.loadPhenoFile(phenoFile);
  if (!covarFile.empty())
    sd.loadCovar(covarFile);
  sd.finalize();
  infoMsg("  %u subjects loaded", sd.nUsed());

  // ---- Build design-matrix projection (intercept + covariates) ----
  infoMsg("Building design matrix projection...");
  Eigen::MatrixXd X;
  if (!covarNames.empty()) {
    auto cov = sd.getColumns(covarNames);
    X.resize(sd.nUsed(), cov.cols() + 1);
    X.col(0).setOnes();
    X.rightCols(cov.cols()) = cov;
  } else if (sd.hasCovar()) {
    // Backward compat: no --covar-name but --covar was given → use all columns
    const auto& cov = sd.covar();
    X.resize(sd.nUsed(), cov.cols() + 1);
    X.col(0).setOnes();
    X.rightCols(cov.cols()) = cov;
  } else {
    X.resize(sd.nUsed(), 1);
    X.col(0).setOnes();
  }
  DesignMatrix design(X);
  infoMsg("  %d rows x %d cols", design.nRows(), design.nCols());

  // ---- Load PLINK data ----
  infoMsg("Loading PLINK data: %s", bfilePrefix.c_str());
  PlinkData plinkData(
      bfilePrefix + ".bed",
      bfilePrefix + ".bim",
      bfilePrefix + ".fam",
      sd.usedMask(),
      sd.nFam(),
      sd.nUsed(),
      "ref-first",
      {}, {}, {}, {},
      nSnpPerChunk);
  infoMsg("  %u subjects matched, %u markers",
          plinkData.nSubjUsed(), plinkData.nMarkers());

  // ---- Per-residual-column loop ----
  const int nRC = sd.residOneCols();
  if (nRC > 1)
    infoMsg("Multi-column residual file: %d columns", nRC);

  for (int rc = 0; rc < nRC; ++rc) {
    Eigen::VectorXd colBuf;
    if (nRC > 1) colBuf = sd.residMatrix().col(rc);
    const Eigen::VectorXd& resid = (nRC > 1) ? colBuf : sd.residuals();

    std::string outFile = (nRC == 1) ? outputFile
        : (outputFile + "." + std::to_string(rc + 1) + ".gz");
    if (nRC > 1)
      infoMsg("  Column %d/%d%s -> %s", rc + 1, nRC,
              (rc < static_cast<int>(sd.residColNames().size())
                   ? (" (" + sd.residColNames()[rc] + ")").c_str() : ""),
              outFile.c_str());

    // ---- Build empirical CGF table ----
    CumulantTable cumul = buildCumulantTable(resid);

    // ---- Variance of residuals ----
    double meanR = resid.mean();
    double varResid = (resid.array() - meanR).square().mean();
    double N = static_cast<double>(resid.size());
    varResid *= N / (N - 1.0);

    // ---- Construct method and run engine ----
    if (nRC == 1) infoMsg("Running SPACox marker tests (%d thread(s))...", nthread);
    SPACoxMethod method(resid, varResid, cumul, design,
                        pvalCovAdjCut, spaCutoff);

    markerEngine(plinkData, method, outFile,
                 nthread,
                 missingCutoff,
                 minMafCutoff,
                 minMacCutoff,
                 /*exactHwe=*/false);
  }
}
