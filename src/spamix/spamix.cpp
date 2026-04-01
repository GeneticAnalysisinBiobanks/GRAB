// spamix.cpp — SPAmix method implementation
//
// Port of ref_code/src/mtSPAmix.cpp to pure C++17 / Eigen / Boost headers.
// Uses scalar kG0/kG1/kG2 loops from math_helper to avoid per-marker heap alloc.

#include "spamix/spamix.hpp"

#include <cmath>
#include <vector>

#include "io/plink.hpp"
#include "io/subject_data.hpp"
#include "spamix/indiv_af.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

// ======================================================================
// SPAmixMethod — construction / clone
// ======================================================================

// On-the-fly AF constructor
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
    m_outlier(outlier),
    m_N(static_cast<int>(residuals.size())),
    m_nPC(static_cast<int>(onePlusPCs.cols()) - 1),
    m_spaCutoff(spaCutoff),
    m_XtX_inv_Xt(&XtX_inv_Xt),
    m_sqrt_XtX_inv_diag(&sqrt_XtX_inv_diag),
    m_afModels(nullptr),
    m_genoToFlat(nullptr),
    m_AFVec(m_N),
    m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
    m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size()))
{}

// Pre-computed AF constructor
SPAmixMethod::SPAmixMethod(
    const Eigen::VectorXd& residuals,
    const Eigen::VectorXd& resid2,
    const Eigen::MatrixXd& onePlusPCs,
    const OutlierData& outlier,
    double spaCutoff,
    const std::vector<AFModel>& afModels,
    const std::vector<uint32_t>& genoToFlat)
  : m_resid(residuals),
    m_resid2(resid2),
    m_onePlusPCs(onePlusPCs),
    m_outlier(outlier),
    m_N(static_cast<int>(residuals.size())),
    m_nPC(static_cast<int>(onePlusPCs.cols()) - 1),
    m_spaCutoff(spaCutoff),
    m_XtX_inv_Xt(nullptr),
    m_sqrt_XtX_inv_diag(nullptr),
    m_afModels(&afModels),
    m_genoToFlat(&genoToFlat),
    m_AFVec(m_N),
    m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
    m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size()))
{}

std::unique_ptr<MethodBase> SPAmixMethod::clone() const {
  if (m_XtX_inv_Xt) {
    return std::make_unique<SPAmixMethod>(
        m_resid, m_resid2, m_onePlusPCs,
        *m_XtX_inv_Xt, *m_sqrt_XtX_inv_diag, m_outlier, m_spaCutoff);
  }
  return std::make_unique<SPAmixMethod>(
      m_resid, m_resid2, m_onePlusPCs,
      m_outlier, m_spaCutoff,
      *m_afModels, *m_genoToFlat);
}

void SPAmixMethod::prepareChunk(const std::vector<uint64_t>& gIndices) {
  m_chunkGenoIndices = gIndices;
}

// ======================================================================
// getMarkerPval — score test + SPA with outlier/non-outlier split
// ======================================================================

double SPAmixMethod::getMarkerPval(
    const Eigen::Ref<const Eigen::VectorXd>& GVec,
    double altFreq, double& zScore) {

  // m_AFVec pre-filled by getResultVec

  // Score, mean, variance — single fused pass over N subjects
  double S = 0.0, S_mean = 0.0, VarS = 0.0;
  for (int i = 0; i < m_N; ++i) {
    const double af = m_AFVec[i];
    S      += GVec[i]    * m_resid[i];
    S_mean += m_resid[i] * af;
    VarS   += m_resid2[i] * af * (1.0 - af);
  }
  S_mean *= 2.0;
  VarS   *= 2.0;

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

  // Gather outlier MAFs
  for (int i = 0; i < nOut; ++i)
    m_mafOutlier[i] = m_AFVec[m_outlier.posOutlier[i]];

  // Gather non-outlier MAFs and accumulate normal-approx terms in one pass
  double mean_nonOutlier = 0.0, var_nonOutlier = 0.0;
  for (int i = 0; i < nNon; ++i) {
    const double af = m_AFVec[m_outlier.posNonOutlier[i]];
    m_mafNonOutlier[i]  = af;
    mean_nonOutlier    += m_outlier.residNonOutlier[i]  * af;
    var_nonOutlier     += m_outlier.resid2NonOutlier[i] * af * (1.0 - af);
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
    double altFreq, int markerInChunkIdx,
    bool /*flipped*/, std::vector<double>& result) {

  // SPAmix handles flipping internally (skipFlip = true)
  double maf = altFreq;
  bool didFlip = false;
  if (maf > 0.5) {
    GVec = 2.0 - GVec.array();
    maf = 1.0 - maf;
    didFlip = true;
  }

  // Fill m_AFVec: on-the-fly or from pre-computed model
  if (m_XtX_inv_Xt) {
    AFContext ctx{m_onePlusPCs, *m_XtX_inv_Xt, *m_sqrt_XtX_inv_diag,
                 m_onePlusPCs.rightCols(m_nPC), m_N, m_nPC};
    computeAFVec(GVec, maf, ctx, m_AFVec);
  } else {
    const uint64_t genoIdx = m_chunkGenoIndices[markerInChunkIdx];
    const uint32_t flatIdx = (*m_genoToFlat)[genoIdx];
    const AFModel& model   = (*m_afModels)[flatIdx];
    getAFVecFromModel(model, maf, m_onePlusPCs, m_N, m_AFVec);
    if (didFlip && model.status != 0)
      m_AFVec = 1.0 - m_AFVec.array();
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
    const std::string& afFile,
    const std::string& outputFile,
    double spaCutoff,
    double outlierRatio,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff) {

  // ---- Load residual file and eigenvectors ----
  infoMsg("Loading residual file: %s", residFile.c_str());
  auto famIIDs = parseFamIIDs(bfilePrefix + ".fam");
  SubjectData sd(std::move(famIIDs));
  sd.loadResidOne(residFile);
  sd.loadEigenVecs(eigenVecsFile);
  sd.finalize();
  infoMsg("  %u subjects loaded, %d PCs", sd.nUsed(), sd.nPC());

  const int N   = static_cast<int>(sd.nUsed());
  const int nPC = sd.nPC();

  // ---- Pre-compute OLS matrices for fit_lm ----
  Eigen::MatrixXd onePlusPCs(N, 1 + nPC);
  onePlusPCs.col(0).setOnes();
  onePlusPCs.rightCols(nPC) = sd.PCs();

  // OLS matrices — only needed for on-the-fly AF computation
  Eigen::MatrixXd XtX_inv_Xt;
  Eigen::VectorXd sqrt_XtX_inv_diag;
  if (afFile.empty()) {
    Eigen::MatrixXd XtX = onePlusPCs.transpose() * onePlusPCs;
    Eigen::MatrixXd XtX_inv =
        XtX.ldlt().solve(Eigen::MatrixXd::Identity(1 + nPC, 1 + nPC));
    XtX_inv_Xt        = XtX_inv * onePlusPCs.transpose();
    sqrt_XtX_inv_diag = XtX_inv.diagonal().cwiseSqrt();
  }

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

  const auto& markerInfo = plinkData.markerInfo();
  const uint32_t nMarkers = static_cast<uint32_t>(markerInfo.size());

  // ---- genoIndex → flat marker index mapping ----
  std::vector<uint32_t> genoToFlat(plinkData.nMarkers(), UINT32_MAX);
  std::vector<uint64_t> genoIndices(nMarkers);
  for (uint32_t fi = 0; fi < nMarkers; ++fi) {
    genoToFlat[markerInfo[fi].genoIndex] = fi;
    genoIndices[fi] = markerInfo[fi].genoIndex;
  }

  // ---- Load pre-computed AF models (if provided) ----
  std::vector<AFModel> afModels;
  if (!afFile.empty()) {
    infoMsg("Loading pre-computed AF models: %s", afFile.c_str());
    afModels = loadAFModels(afFile, nPC, nMarkers, genoIndices);
    infoMsg("  %zu AF models loaded", afModels.size());
  }

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

    // ---- Squared residuals ----
    Eigen::VectorXd resid2 = resid.array().square();

    // ---- Outlier detection ----
    OutlierData outlier = detectOutliers(resid, outlierRatio);
    if (nRC == 1)
      infoMsg("  %zu outliers, %zu non-outliers",
              outlier.posOutlier.size(), outlier.posNonOutlier.size());

    // ---- Construct method and run engine ----
    if (nRC == 1) infoMsg("Running SPAmix marker tests (%d thread(s))...", nthread);
    std::unique_ptr<SPAmixMethod> method;
    if (!afFile.empty()) {
      method = std::make_unique<SPAmixMethod>(
          resid, resid2, onePlusPCs,
          outlier, spaCutoff,
          afModels, genoToFlat);
    } else {
      method = std::make_unique<SPAmixMethod>(
          resid, resid2, onePlusPCs,
          XtX_inv_Xt, sqrt_XtX_inv_diag, outlier, spaCutoff);
    }

    markerEngine(plinkData, *method, outFile,
                 nthread,
                 missingCutoff,
                 minMafCutoff,
                 minMacCutoff,
                 /*exactHwe=*/false);
  }
}
