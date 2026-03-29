// spamix.cpp — SPAmix method implementation
//
// Port of ref_code/src/mtSPAmix.cpp to pure C++17 / Eigen / Boost headers.
// Uses scalar kG0/kG1/kG2 loops from math_helper to avoid per-marker heap alloc.

#include "spamix/spamix.hpp"

#include <cmath>
#include <vector>

#include "io/data_matrix.hpp"
#include "io/plink.hpp"
#include "io/resid_file.hpp"
#include "spamix/indiv_af.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

// ======================================================================
// SPAmixMethod — construction / clone
// ======================================================================

SPAmixMethod::SPAmixMethod(
    const Eigen::VectorXd& residuals,
    const Eigen::VectorXd& resid2,
    const Eigen::MatrixXd& onePlusPCs,
    const Eigen::MatrixXd& XtX_inv_Xt,
    const Eigen::VectorXd& sqrt_XtX_inv_diag,
    const OutlierData& outlier,
    double spaCutoff)
  : m_resid(residuals),
    m_resid2(resid2),
    m_onePlusPCs(onePlusPCs),
    m_XtX_inv_Xt(XtX_inv_Xt),
    m_sqrt_XtX_inv_diag(sqrt_XtX_inv_diag),
    m_outlier(outlier),
    m_N(static_cast<int>(residuals.size())),
    m_nPC(static_cast<int>(onePlusPCs.cols()) - 1),
    m_spaCutoff(spaCutoff),
    m_AFVec(m_N),
    m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
    m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size()))
{}

std::unique_ptr<MethodBase> SPAmixMethod::clone() const {
  return std::make_unique<SPAmixMethod>(
      m_resid, m_resid2, m_onePlusPCs,
      m_XtX_inv_Xt, m_sqrt_XtX_inv_diag, m_outlier, m_spaCutoff);
}

// ======================================================================
// getMarkerPval — score test + SPA with outlier/non-outlier split
// ======================================================================

double SPAmixMethod::getMarkerPval(
    const Eigen::Ref<const Eigen::VectorXd>& GVec,
    double altFreq, double& zScore) {

  // Estimate per-individual AF → m_AFVec
  AFContext ctx{m_onePlusPCs, m_XtX_inv_Xt, m_sqrt_XtX_inv_diag,
               m_onePlusPCs.rightCols(m_nPC), m_N, m_nPC};
  computeAFVec(GVec, altFreq, ctx, m_AFVec);

  // Score statistic
  double S      = GVec.dot(m_resid);
  double S_mean = 2.0 * m_resid.dot(m_AFVec);

  // Variance: Σ resid² · 2 · AF · (1−AF)
  double VarS = 0.0;
  for (int i = 0; i < m_N; ++i)
    VarS += m_resid2[i] * 2.0 * m_AFVec[i] * (1.0 - m_AFVec[i]);

  if (VarS <= 0.0) {
    zScore = 0.0;
    return 1.0;
  }

  zScore = (S - S_mean) / std::sqrt(VarS);

  // Normal approximation when |z| < cutoff
  if (std::abs(zScore) < m_spaCutoff)
    return 2.0 * math::pnorm(-std::abs(zScore));

  // ---- SPA with outlier / non-outlier split ----

  const int nOut = static_cast<int>(m_outlier.posOutlier.size());
  const int nNon = static_cast<int>(m_outlier.posNonOutlier.size());

  // Gather per-individual MAF at outlier / non-outlier positions
  for (int i = 0; i < nOut; ++i)
    m_mafOutlier[i] = m_AFVec[m_outlier.posOutlier[i]];
  for (int i = 0; i < nNon; ++i)
    m_mafNonOutlier[i] = m_AFVec[m_outlier.posNonOutlier[i]];

  // Non-outlier normal approximation terms
  double mean_nonOutlier = 0.0, var_nonOutlier = 0.0;
  for (int i = 0; i < nNon; ++i) {
    double af = m_mafNonOutlier[i];
    mean_nonOutlier += m_outlier.residNonOutlier[i] * af;
    var_nonOutlier  += m_outlier.resid2NonOutlier[i] * af * (1.0 - af);
  }
  mean_nonOutlier *= 2.0;
  var_nonOutlier  *= 2.0;

  double absS = std::abs(S - S_mean);

  double pval1 = spa::getProbSpaG(
      m_mafOutlier.data(), m_outlier.residOutlier.data(), nOut,
      absS + S_mean, false,
      mean_nonOutlier, var_nonOutlier);

  double pval2 = spa::getProbSpaG(
      m_mafOutlier.data(), m_outlier.residOutlier.data(), nOut,
      -absS + S_mean, true,
      mean_nonOutlier, var_nonOutlier);

  return pval1 + pval2;
}

// ======================================================================
// getResultVec — MethodBase interface (handles flip internally)
// ======================================================================

void SPAmixMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq, int /*markerInChunkIdx*/,
    bool /*flipped*/, std::vector<double>& result) {

  // SPAmix handles flipping internally (skipFlip = true)
  double maf = altFreq;
  if (maf > 0.5) {
    GVec = 2.0 - GVec.array();
    maf = 1.0 - maf;
  }

  double zScore;
  double pval = getMarkerPval(GVec, maf, zScore);
  result.push_back(pval);
  result.push_back(zScore);
}

// ======================================================================
// runSPAmix — top-level orchestration
// ======================================================================

void runSPAmix(
    const std::string& residFile,
    const std::string& eigenVecsFile,
    const std::string& bfilePrefix,
    const std::string& outputFile,
    double spaCutoff,
    double outlierRatio,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff) {

  // ---- Load residual file ----
  infoMsg("Loading residual file: %s", residFile.c_str());
  ResidData resid = loadResidFile(residFile);
  infoMsg("  %zu subjects loaded", resid.subjects.size());

  // ---- Load eigenvectors (PCs) ----
  infoMsg("Loading eigenvector file: %s", eigenVecsFile.c_str());
  EigenVecData evd = loadEigenVecs(eigenVecsFile);
  Eigen::MatrixXd PCs = std::move(evd.PCs);
  infoMsg("  %zu subjects, %d PCs",
          evd.subjects.size(), static_cast<int>(PCs.cols()));

  if (static_cast<Eigen::Index>(evd.subjects.size()) !=
      static_cast<Eigen::Index>(resid.subjects.size()))
    throw std::runtime_error(
        "Eigenvector subject count (" + std::to_string(evd.subjects.size()) +
        ") does not match residual file (" +
        std::to_string(resid.subjects.size()) + ")");

  for (size_t i = 0; i < resid.subjects.size(); ++i) {
    if (evd.subjects[i] != resid.subjects[i])
      throw std::runtime_error(
          "Subject order mismatch at row " + std::to_string(i + 1) +
          ": eigen-vecs has '" + evd.subjects[i] +
          "', resid-file has '" + resid.subjects[i] + "'");
  }

  const int N   = static_cast<int>(resid.subjects.size());
  const int nPC = static_cast<int>(PCs.cols());

  // ---- Pre-compute OLS matrices for fit_lm ----
  Eigen::MatrixXd onePlusPCs(N, 1 + nPC);
  onePlusPCs.col(0).setOnes();
  onePlusPCs.rightCols(nPC) = PCs;
  PCs.resize(0, 0);  // free — data now lives only in onePlusPCs

  Eigen::MatrixXd XtX = onePlusPCs.transpose() * onePlusPCs;
  Eigen::MatrixXd XtX_inv =
      XtX.ldlt().solve(Eigen::MatrixXd::Identity(1 + nPC, 1 + nPC));
  Eigen::MatrixXd XtX_inv_Xt = XtX_inv * onePlusPCs.transpose();
  Eigen::VectorXd sqrt_XtX_inv_diag = XtX_inv.diagonal().cwiseSqrt();

  // ---- Squared residuals ----
  Eigen::VectorXd resid2 = resid.residuals.array().square();

  // ---- Outlier detection ----
  infoMsg("Detecting outlier residuals (ratio=%.2f)...", outlierRatio);
  OutlierData outlier = detectOutliers(resid.residuals, outlierRatio);
  infoMsg("  %zu outliers, %zu non-outliers",
          outlier.posOutlier.size(), outlier.posNonOutlier.size());

  // ---- Load PLINK data ----
  infoMsg("Loading PLINK data: %s", bfilePrefix.c_str());
  PlinkData plinkData(
      bfilePrefix + ".bed",
      bfilePrefix + ".bim",
      bfilePrefix + ".fam",
      resid.subjects,
      "ref-first",
      {}, {}, {}, {},
      nSnpPerChunk);
  infoMsg("  %u subjects matched, %u markers",
          plinkData.nSubjUsed(), plinkData.nMarkers());

  // ---- Construct method and run engine ----
  infoMsg("Running SPAmix marker tests (%d thread(s))...", nthread);
  SPAmixMethod method(resid.residuals, resid2, onePlusPCs,
                      XtX_inv_Xt, sqrt_XtX_inv_diag, outlier, spaCutoff);

  markerEngine(plinkData, method, outputFile,
               nthread,
               missingCutoff,
               minMafCutoff,
               minMacCutoff,
               /*exactHwe=*/false);
}
