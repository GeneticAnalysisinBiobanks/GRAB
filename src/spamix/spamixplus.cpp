// spamixplus.cpp — SPAmixPlus method implementation
//
// Port of ref_code/src/mtSPAmixPlus.cpp to pure C++17 / Eigen / Boost headers.
// Key difference from SPAmix: sparse-GRM-based variance and ratio correction.

#include "spamix/spamixplus.hpp"

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <exception>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <vector>

#include <zlib.h>

#include "io/data_matrix.hpp"
#include "io/plink.hpp"
#include "io/resid_file.hpp"
#include "spamix/indiv_af.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

// ======================================================================
// SPAmixPlusMethod — construction / clone
// ======================================================================

// Pre-computed AF constructor
SPAmixPlusMethod::SPAmixPlusMethod(
    const Eigen::VectorXd& residuals,
    const Eigen::VectorXd& resid2,
    const Eigen::MatrixXd& onePlusPCs,
    const OutlierData& outlier,
    double spaCutoff,
    const SparseGRM& grm,
    const std::vector<AFModel>& afModels,
    const std::vector<uint32_t>& genoToFlat)
  : m_resid(residuals),
    m_resid2(resid2),
    m_onePlusPCs(onePlusPCs),
    m_outlier(outlier),
    m_spaCutoff(spaCutoff),
    m_grm(grm),
    m_N(static_cast<int>(residuals.size())),
    m_nPC(static_cast<int>(onePlusPCs.cols()) - 1),
    m_afModels(&afModels),
    m_genoToFlat(&genoToFlat),
    m_XtX_inv_Xt(nullptr),
    m_sqrt_XtX_inv_diag(nullptr),
    m_AFVec(m_N),
    m_R_new(m_N),
    m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
    m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size()))
{}

// On-the-fly AF constructor
SPAmixPlusMethod::SPAmixPlusMethod(
    const Eigen::VectorXd& residuals,
    const Eigen::VectorXd& resid2,
    const Eigen::MatrixXd& onePlusPCs,
    const OutlierData& outlier,
    double spaCutoff,
    const SparseGRM& grm,
    const Eigen::MatrixXd& XtX_inv_Xt,
    const Eigen::VectorXd& sqrt_XtX_inv_diag,
    int nPC)
  : m_resid(residuals),
    m_resid2(resid2),
    m_onePlusPCs(onePlusPCs),
    m_outlier(outlier),
    m_spaCutoff(spaCutoff),
    m_grm(grm),
    m_N(static_cast<int>(residuals.size())),
    m_nPC(nPC),
    m_afModels(nullptr),
    m_genoToFlat(nullptr),
    m_XtX_inv_Xt(&XtX_inv_Xt),
    m_sqrt_XtX_inv_diag(&sqrt_XtX_inv_diag),
    m_AFVec(m_N),
    m_R_new(m_N),
    m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
    m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size()))
{}

std::unique_ptr<MethodBase> SPAmixPlusMethod::clone() const {
  if (m_XtX_inv_Xt) {
    return std::make_unique<SPAmixPlusMethod>(
        m_resid, m_resid2, m_onePlusPCs, m_outlier, m_spaCutoff,
        m_grm, *m_XtX_inv_Xt, *m_sqrt_XtX_inv_diag, m_nPC);
  }
  return std::make_unique<SPAmixPlusMethod>(
      m_resid, m_resid2, m_onePlusPCs, m_outlier, m_spaCutoff,
      m_grm, *m_afModels, *m_genoToFlat);
}

void SPAmixPlusMethod::prepareChunk(const std::vector<uint64_t>& gIndices) {
  m_chunkGenoIndices = gIndices;
}

// ======================================================================
// getMarkerPval — score test with GRM variance + SPA
// ======================================================================

double SPAmixPlusMethod::getMarkerPval(
    const Eigen::Ref<const Eigen::VectorXd>& GVec,
    double altFreq, double& zScore) {

  // m_AFVec is pre-filled by caller (on-the-fly or from stored model)

  // Build R_new, S, S_mean, S_var_SPAmix in a single pass
  double S = 0.0, S_mean = 0.0, S_var_SPAmix = 0.0;
  for (int i = 0; i < m_N; ++i) {
    double af  = m_AFVec[i];
    double gvar = 2.0 * af * (1.0 - af);
    m_R_new[i]   = m_resid[i] * std::sqrt(gvar);
    S           += GVec[i] * m_resid[i];
    S_mean      += m_resid[i] * af;
    S_var_SPAmix += m_resid2[i] * gvar;
  }
  S_mean *= 2.0;

  // GRM-based variance
  double VarS = m_grm.spaVariance(m_R_new.data(),
                                  static_cast<uint32_t>(m_N));

  if (VarS <= 0.0) {
    zScore = 0.0;
    return 1.0;
  }

  zScore = (S - S_mean) / std::sqrt(VarS);

  // Normal approximation when |z| < cutoff
  if (std::abs(zScore) < m_spaCutoff)
    return 2.0 * math::pnorm(-std::abs(zScore));

  // ---- SPA with variance-ratio correction ----

  double Var_ratio  = S_var_SPAmix / VarS;
  double sqrt_ratio = std::sqrt(Var_ratio);
  double S_new      = S * sqrt_ratio;
  double S_mean_new = S_mean * sqrt_ratio;
  double S_upper    = std::max(S_new, 2.0 * S_mean_new - S_new);
  double S_lower    = std::min(S_new, 2.0 * S_mean_new - S_new);

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

  double pval1 = spa::getProbSpaG(
      m_mafOutlier.data(), m_outlier.residOutlier.data(), nOut,
      S_upper, false,
      mean_nonOutlier, var_nonOutlier);

  double pval2 = spa::getProbSpaG(
      m_mafOutlier.data(), m_outlier.residOutlier.data(), nOut,
      S_lower, true,
      mean_nonOutlier, var_nonOutlier);

  return pval1 + pval2;
}

// ======================================================================
// getResultVec — MethodBase interface (handles flip internally)
// ======================================================================

void SPAmixPlusMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq, int markerInChunkIdx,
    bool /*flipped*/, std::vector<double>& result) {

  // SPAmixPlus handles flipping internally (skipFlip = true)
  double maf = altFreq;
  bool didFlip = false;
  if (maf > 0.5) {
    GVec = 2.0 - GVec.array();
    maf = 1.0 - maf;
    didFlip = true;
  }

  // Fill m_AFVec: on-the-fly from genotype, or from pre-computed model
  if (m_XtX_inv_Xt) {
    AFContext ctx{m_onePlusPCs, *m_XtX_inv_Xt, *m_sqrt_XtX_inv_diag,
                 m_onePlusPCs.rightCols(m_nPC), m_N, m_nPC};
    computeAFVec(GVec, maf, ctx, m_AFVec);
  } else {
    const uint64_t genoIdx  = m_chunkGenoIndices[markerInChunkIdx];
    const uint32_t flatIdx  = (*m_genoToFlat)[genoIdx];
    const AFModel& model    = (*m_afModels)[flatIdx];
    getAFVecFromModel(model, maf, m_onePlusPCs, m_N, m_AFVec);
    // Model was fitted on original (unflipped) genotype.  For status 1/2
    // the betas encode AF for the original allele; flip to match test allele.
    if (didFlip && model.status != 0)
      m_AFVec = 1.0 - m_AFVec.array();
  }

  double zScore;
  double pval = getMarkerPval(GVec, maf, zScore);
  result.push_back(pval);
  result.push_back(zScore);
}

// ======================================================================
// AF model loading helpers
// ======================================================================

namespace {

// Load AF models from a binary file by genoIndex.
std::vector<AFModel> loadAFModelsBinary(
    const std::string& binFile,
    int nPC,
    const std::vector<PlinkData::MarkerInfo>& markerInfo) {

  const uint32_t nMarkers = static_cast<uint32_t>(markerInfo.size());
  std::vector<AFModel> models(nMarkers);
  IndivAFReader reader(binFile, nPC);
  for (uint32_t fi = 0; fi < nMarkers; ++fi)
    reader.read(markerInfo[fi].genoIndex, models[fi]);
  return models;
}

// Load AF models from a text or gzip-text file (rows in flat marker order).
std::vector<AFModel> loadAFModelsText(
    const std::string& path,
    int nPC,
    uint32_t nExpected) {

  // Detect gzip
  const bool isGz = path.size() > 3 &&
      path.compare(path.size() - 3, 3, ".gz") == 0;

  std::vector<AFModel> models;
  models.reserve(nExpected);

  auto parseLine = [&](const std::string& line) {
    if (line.empty() || line[0] == 'C') return;  // skip header (CHR...)
    std::istringstream iss(line);
    std::string chr;
    uint32_t bp;
    int status;
    iss >> chr >> bp >> status;
    AFModel m;
    m.status = static_cast<int8_t>(status);
    m.betas.resize(1 + nPC);
    for (int j = 0; j <= nPC; ++j)
      iss >> m.betas[j];
    models.push_back(std::move(m));
  };

  if (isGz) {
    gzFile gz = gzopen(path.c_str(), "rb");
    if (!gz)
      throw std::runtime_error("Cannot open gzip AF file: " + path);
    char buf[8192];
    std::string leftover;
    while (true) {
      int nRead = gzread(gz, buf, sizeof(buf));
      if (nRead <= 0) break;
      leftover.append(buf, nRead);
      size_t pos = 0;
      while (true) {
        size_t nl = leftover.find('\n', pos);
        if (nl == std::string::npos) {
          leftover = leftover.substr(pos);
          break;
        }
        parseLine(leftover.substr(pos, nl - pos));
        pos = nl + 1;
      }
    }
    if (!leftover.empty()) parseLine(leftover);
    gzclose(gz);
  } else {
    std::ifstream ifs(path);
    if (!ifs)
      throw std::runtime_error("Cannot open text AF file: " + path);
    std::string line;
    while (std::getline(ifs, line))
      parseLine(line);
  }

  return models;
}

} // anonymous namespace


// ======================================================================
// runSPAmixPlus — top-level orchestration
// ======================================================================

void runSPAmixPlus(
    const std::string& residFile,
    const std::string& eigenVecsFile,
    const std::string& bfilePrefix,
    const std::string& sparseGrmFile,
    const std::string& afFile,
    const std::string& outputFile,
    double spaCutoff,
    double outlierRatio,
    int    nthread,
    int    nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff)
{
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

  // ---- Design matrix: [1 | PCs] ----
  Eigen::MatrixXd onePlusPCs(N, 1 + nPC);
  onePlusPCs.col(0).setOnes();
  onePlusPCs.rightCols(nPC) = PCs;
  PCs.resize(0, 0);  // free — data now lives only in onePlusPCs

  // Squared residuals (shared across all threads)
  Eigen::VectorXd resid2 = resid.residuals.array().square();

  // OLS matrices — only needed for on-the-fly AF computation
  Eigen::MatrixXd XtX_inv_Xt;
  Eigen::VectorXd sqrt_XtX_inv_diag;
  if (afFile.empty()) {
    Eigen::MatrixXd XtX = onePlusPCs.transpose() * onePlusPCs;
    Eigen::MatrixXd XtX_inv =
        XtX.ldlt().solve(Eigen::MatrixXd::Identity(1 + nPC, 1 + nPC));
    XtX_inv_Xt    = XtX_inv * onePlusPCs.transpose();
    sqrt_XtX_inv_diag = XtX_inv.diagonal().cwiseSqrt();
  }

  // ---- Outlier detection ----
  infoMsg("Detecting outlier residuals (ratio=%.2f)...", outlierRatio);
  OutlierData outlier = detectOutliers(resid.residuals, outlierRatio);
  infoMsg("  %zu outliers, %zu non-outliers",
          outlier.posOutlier.size(), outlier.posNonOutlier.size());

  // ---- Load sparse GRM (raw mode — lower triangle + diagonal) ----
  infoMsg("Loading sparse GRM (raw): %s", sparseGrmFile.c_str());
  SparseGRM grm(sparseGrmFile, resid.subjects, /*symmetrize=*/false);
  infoMsg("  %u subjects, %zu entries", grm.nSubjects(), grm.nnz());

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

  const auto& markerInfo = plinkData.markerInfo();
  const uint32_t nMarkers = static_cast<uint32_t>(markerInfo.size());

  // ---- genoIndex → flat marker index mapping ----
  std::vector<uint32_t> genoToFlat(plinkData.nMarkers(), UINT32_MAX);
  for (uint32_t fi = 0; fi < nMarkers; ++fi)
    genoToFlat[markerInfo[fi].genoIndex] = fi;

  // ---- Load AF models (pre-computed) or prepare on-the-fly ----
  std::vector<AFModel> afModels;
  if (!afFile.empty()) {
    infoMsg("Loading pre-computed AF models: %s", afFile.c_str());
    const auto len = afFile.size();
    if (len > 4 && afFile.compare(len - 4, 4, ".bin") == 0) {
      afModels = loadAFModelsBinary(afFile, nPC, markerInfo);
    } else {
      afModels = loadAFModelsText(afFile, nPC, nMarkers);
    }
    infoMsg("  %zu AF models loaded", afModels.size());

    if (afModels.size() != nMarkers)
      throw std::runtime_error(
          "AF model count (" + std::to_string(afModels.size()) +
          ") does not match marker count (" +
          std::to_string(nMarkers) + ")");
  }
  // On-the-fly: AF computed per-marker inside the engine (single PLINK pass)

  // ---- Construct method and run engine ----
  infoMsg("Running SPAmixPlus marker tests (%d thread(s))...", nthread);
  std::unique_ptr<SPAmixPlusMethod> method;
  if (!afFile.empty()) {
    method = std::make_unique<SPAmixPlusMethod>(
        resid.residuals, resid2, onePlusPCs, outlier, spaCutoff,
        grm, afModels, genoToFlat);
  } else {
    infoMsg("AF models computed on-the-fly during marker testing.");
    method = std::make_unique<SPAmixPlusMethod>(
        resid.residuals, resid2, onePlusPCs, outlier, spaCutoff,
        grm, XtX_inv_Xt, sqrt_XtX_inv_diag, nPC);
  }

  markerEngine(plinkData, *method, outputFile,
               nthread,
               missingCutoff,
               minMafCutoff,
               minMacCutoff,
               /*exactHwe=*/false);
}
