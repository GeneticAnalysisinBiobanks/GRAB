// spasqr.cpp — SPAsqr: SPA-squared multi-tau marker association (pure C++17 / Eigen / Boost)

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "spasqr/spasqr.hpp"
#include "spasqr/conquer.hpp"
#include "spagrm/spagrm.hpp"
#include "engine/marker.hpp"
#include "io/subject_data.hpp"
#include "io/geno_data.hpp"
#include "io/sparse_grm.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <Eigen/Dense>


// ══════════════════════════════════════════════════════════════════════
// SPAsqrMethod — MethodBase adapter wrapping ntaus SPAGRMClass instances
// ══════════════════════════════════════════════════════════════════════

namespace {

class SPAsqrMethod : public MethodBase {
public:
  // Constructor with tau labels (used by pheno path):
  //   labels like "Tau0.1", "Tau0.3", ...
  //   Output order: P_CCT  P_tau0.1 ... P_tau0.9  Z_tau0.1 ... Z_tau0.9
  SPAsqrMethod(
      int ntaus,
      std::vector<SPAGRMClass> spagrm_vec,
      std::vector<std::string> tauLabels)
    : m_ntaus(ntaus),
      m_spagrm_vec(std::move(spagrm_vec)),
      m_tauLabels(std::move(tauLabels)),
      m_hasLabels(true)
  {}

  // Constructor without labels (used by legacy --null-resid path):
  //   Output order: Z_tau1 ... Z_tauK  P_tau1 ... P_tauK  P_CCT
  SPAsqrMethod(
      int ntaus,
      std::vector<SPAGRMClass> spagrm_vec)
    : m_ntaus(ntaus),
      m_spagrm_vec(std::move(spagrm_vec)),
      m_hasLabels(false)
  {}

  std::unique_ptr<MethodBase> clone() const override {
    return std::make_unique<SPAsqrMethod>(*this);
  }

  int resultSize() const override { return 2 * m_ntaus + 1; }

  std::string getHeaderColumns() const override {
    std::ostringstream oss;
    if (m_hasLabels) {
      // Pheno path: P_CCT first, then P_tau{val}..., then Z_tau{val}...
      oss << "\tP_CCT";
      for (int i = 0; i < m_ntaus; ++i)
        oss << "\tP_" << m_tauLabels[i];
      for (int i = 0; i < m_ntaus; ++i)
        oss << "\tZ_" << m_tauLabels[i];
    } else {
      // Legacy: Z_tau1..Z_tauK  P_tau1..P_tauK  P_CCT
      for (int i = 1; i <= m_ntaus; ++i)
        oss << "\tZ_tau" << i;
      for (int i = 1; i <= m_ntaus; ++i)
        oss << "\tP_tau" << i;
      oss << "\tP_CCT";
    }
    return oss.str();
  }

  void getResultVec(
      Eigen::Ref<Eigen::VectorXd> GVec,
      double altFreq,
      int /*markerInChunkIdx*/,
      std::vector<double>& result) override
  {
    result.clear();
    result.reserve(2 * m_ntaus + 1);

    std::vector<double> zScores(m_ntaus);
    std::vector<double> pvals(m_ntaus);

    for (int i = 0; i < m_ntaus; ++i) {
      double z;
      double p = m_spagrm_vec[i].getMarkerPval(GVec, altFreq, z);
      zScores[i] = z;
      pvals[i]   = p;
    }

    // CCT (Cauchy combination test) p-value
    std::vector<double> valid_p;
    valid_p.reserve(m_ntaus);
    for (int i = 0; i < m_ntaus; ++i)
      if (!std::isnan(pvals[i])) valid_p.push_back(pvals[i]);

    double pCCT = std::numeric_limits<double>::quiet_NaN();
    if (!valid_p.empty()) {
      bool hasZero = false;
      double tStat = 0.0;
      for (double p : valid_p) {
        if (p <= 0.0) { hasZero = true; break; }
        double pc = (p >= 1.0) ? 0.999 : p;
        tStat += std::tan((0.5 - pc) * M_PI);
      }
      if (hasZero) {
        pCCT = 0.0;
      } else {
        tStat /= static_cast<double>(valid_p.size());
        pCCT = (tStat > 1e15) ? (1.0 / tStat) / M_PI
                               : 0.5 - std::atan(tStat) / M_PI;
      }
    }

    if (m_hasLabels) {
      // P_CCT first, then P_tau{val}..., then Z_tau{val}...
      result.push_back(pCCT);
      for (int i = 0; i < m_ntaus; ++i) result.push_back(pvals[i]);
      for (int i = 0; i < m_ntaus; ++i) result.push_back(zScores[i]);
    } else {
      // Legacy: Z_tau1..Z_tauK  P_tau1..P_tauK  P_CCT
      for (int i = 0; i < m_ntaus; ++i) result.push_back(zScores[i]);
      for (int i = 0; i < m_ntaus; ++i) result.push_back(pvals[i]);
      result.push_back(pCCT);
    }
  }

private:
  int m_ntaus;
  std::vector<SPAGRMClass> m_spagrm_vec;
  std::vector<std::string> m_tauLabels;
  bool m_hasLabels;
};


// ══════════════════════════════════════════════════════════════════════
// Outlier detection (IQR-based, per column)
// ══════════════════════════════════════════════════════════════════════

// Returns an N × K boolean matrix (as std::vector<std::vector<bool>>).
// outlierIqrRatio  = multiplier for IQR (default 1.5)
// outlierAbsBound  = absolute clamp for cutoffs (default 0.55)
struct OutlierInfo {
  // per-column: indices of outlier subjects
  std::vector<std::vector<int>> outlierIdx;
  // per-column: boolean mask (size N)
  std::vector<std::vector<bool>> isOutlier;
};

OutlierInfo detectOutliers(
    const Eigen::MatrixXd& ResidMat,
    double outlierIqrRatio,
    double outlierAbsBound
) {
  const Eigen::Index N = ResidMat.rows();
  const Eigen::Index K = ResidMat.cols();

  OutlierInfo info;
  info.outlierIdx.resize(K);
  info.isOutlier.resize(K);

  // Scratch for sorting
  std::vector<double> scratch(N);

  for (Eigen::Index col = 0; col < K; ++col) {
    // Copy column for sorting
    for (Eigen::Index i = 0; i < N; ++i)
      scratch[i] = ResidMat(i, col);

    std::sort(scratch.begin(), scratch.end());

    // Q1, Q3 via linear interpolation (same as R type=7)
    auto quantile = [&](double prob) -> double {
      const double idx = prob * (N - 1);
      const Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
      const Eigen::Index hi = std::min(lo + 1, N - 1);
      const double frac = idx - lo;
      return scratch[lo] * (1.0 - frac) + scratch[hi] * frac;
    };

    const double Q1 = quantile(0.25);
    const double Q3 = quantile(0.75);
    const double IQR = Q3 - Q1;

    double cutLo = Q1 - outlierIqrRatio * IQR;
    double cutHi = Q3 + outlierIqrRatio * IQR;
    cutLo = std::max(cutLo, -outlierAbsBound);
    cutHi = std::min(cutHi,  outlierAbsBound);

    info.isOutlier[col].resize(N, false);
    int nOutlier = 0;
    for (Eigen::Index i = 0; i < N; ++i) {
      const double v = ResidMat(i, col);
      if (v < cutLo || v > cutHi) {
        info.isOutlier[col][i] = true;
        info.outlierIdx[col].push_back(static_cast<int>(i));
        ++nOutlier;
      }
    }

    infoMsg("  Column %d: outlier ratio = %.4f (%d / %lld)",
            static_cast<int>(col + 1),
            static_cast<double>(nOutlier) / N,
            nOutlier, static_cast<long long>(N));
  }
  return info;
}

} // anonymous namespace


// ══════════════════════════════════════════════════════════════════════
// Shared pipeline: outlier detection → GRM → SPAGRM → marker engine
// ══════════════════════════════════════════════════════════════════════

static void runSPAsqrPipeline(
    Eigen::MatrixXd& ResidMat,
    GenoMeta& plinkData,
    const std::vector<std::string>& subjOrder,
    const std::vector<std::string>& famIIDs,
    uint32_t nUsed,
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const std::string& outputFile,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    std::vector<std::string> tauLabels)  // empty → legacy column names
{
  const Eigen::Index K = ResidMat.cols();

  // ── 3. Outlier detection ───────────────────────────────────────────
  infoMsg("Detecting outliers (IQR ratio=%.2f, abs bound=%.2f)",
          outlierIqrRatio, outlierAbsBound);
  OutlierInfo outlierInfo = detectOutliers(ResidMat, outlierIqrRatio, outlierAbsBound);

  // ── 4. Load sparse GRM and compute variance terms ──────────────────

  struct GRMEntry {
    uint32_t row, col;
    double value;
    double factor;  // 1 for diagonal, 2 for off-diagonal
  };

  std::vector<GRMEntry> grmEntries;

  {
    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile,
                                     subjOrder, famIIDs);
    infoMsg("Sparse GRM: %zu entries after filtering", grm.nnz());
    grmEntries.reserve(grm.nnz());
    for (const auto& e : grm.entries())
      grmEntries.push_back({e.row, e.col, e.value,
                            (e.row == e.col) ? 1.0 : 2.0});
  }

  // ── 5. Compute per-column variance terms ───────────────────────────
  const int ntaus = static_cast<int>(K);
  std::vector<double> R_GRM_R_vec(ntaus, 0.0);
  std::vector<double> R_GRM_R_nonOutlier_vec(ntaus, 0.0);
  std::vector<double> sum_R_nonOutlier_vec(ntaus, 0.0);
  std::vector<Eigen::VectorXd> resid_outliers_lst(ntaus);

  for (int col = 0; col < ntaus; ++col) {
    const auto& isOut = outlierInfo.isOutlier[col];

    // R_GRM_R and R_GRM_R_nonOutlier
    double rgrm_r = 0.0;
    double rgrm_r_no = 0.0;
    for (const auto& e : grmEntries) {
      const double contrib = e.factor * e.value *
          ResidMat(e.row, col) * ResidMat(e.col, col);
      rgrm_r += contrib;
      if (!isOut[e.row] && !isOut[e.col])
        rgrm_r_no += contrib;
    }
    R_GRM_R_vec[col] = rgrm_r;
    R_GRM_R_nonOutlier_vec[col] = rgrm_r_no;

    // sum_R_nonOutlier and resid outlier values
    double sumNO = 0.0;
    std::vector<double> outVals;
    outVals.reserve(outlierInfo.outlierIdx[col].size());
    for (uint32_t i = 0; i < nUsed; ++i) {
      if (!isOut[i])
        sumNO += ResidMat(i, col);
      else
        outVals.push_back(ResidMat(i, col));
    }
    sum_R_nonOutlier_vec[col] = sumNO;
    resid_outliers_lst[col] = Eigen::Map<Eigen::VectorXd>(
        outVals.data(), static_cast<Eigen::Index>(outVals.size()));
  }

  // ── 6. Build SPAGRMClass instances (one per tau) ───────────────────
  static const std::vector<double> MAF_interval =
      {0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
  const double zeta = 0.01;
  const double tol  = 1e-6;

  std::vector<SPAGRMClass> spagrm_vec;
  spagrm_vec.reserve(ntaus);

  for (int col = 0; col < ntaus; ++col) {
    // SPAsqr: no families — empty FamilyData
    nsSPAGRM::FamilyData fam;
    fam.resid_unrelated_outliers = std::move(resid_outliers_lst[col]);

    spagrm_vec.emplace_back(
        ResidMat.col(col),             // resid
        sum_R_nonOutlier_vec[col],
        R_GRM_R_nonOutlier_vec[col],
        0.0,                           // R_GRM_R_TwoSubjOutlier = 0 (no families)
        R_GRM_R_vec[col],
        MAF_interval,
        std::move(fam),
        spaCutoff,
        zeta,
        tol);
  }

  // ── 7. Run marker engine ───────────────────────────────────────────
  std::unique_ptr<SPAsqrMethod> method;
  if (tauLabels.empty())
    method = std::make_unique<SPAsqrMethod>(ntaus, std::move(spagrm_vec));
  else
    method = std::make_unique<SPAsqrMethod>(ntaus, std::move(spagrm_vec), std::move(tauLabels));

  infoMsg("Starting marker-level association (SPAsqr, %d taus, %d threads)",
          ntaus, nthreads);
  markerEngine(plinkData, *method, outputFile,
               nthreads, missingCutoff, minMafCutoff, minMacCutoff,
               /*exactHwe=*/false);
}


// ══════════════════════════════════════════════════════════════════════
// runSPAsqr — legacy entry point (--null-resid)
// ══════════════════════════════════════════════════════════════════════

void runSPAsqr(
    const std::string& residFile,
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const GenoSpec& geno,
    const std::string& outputFile,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff
) {
  // ── 1. Load residual matrix ────────────────────────────────────────
  infoMsg("Loading residual matrix from %s", residFile.c_str());
  auto famIIDs = parseGenoIIDs(geno);
  SubjectData sd(std::move(famIIDs));
  sd.loadResidSPAsqr(residFile);
  sd.finalize();
  const Eigen::Index N = static_cast<Eigen::Index>(sd.nUsed());
  const Eigen::Index K = sd.residCols();
  infoMsg("Residual matrix: %lld subjects x %lld columns",
          static_cast<long long>(N), static_cast<long long>(K));

  // ── 2. Load genotype data ──────────────────────────────────────────
  auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(),
                               nSnpPerChunk);

  Eigen::MatrixXd ResidMat = sd.residMatrix();

  runSPAsqrPipeline(
      ResidMat, *genoData,
      sd.usedIIDs(), sd.famIIDs(),
      genoData->nSubjUsed(),
      spgrmGrabFile, spgrmGctaFile,
      outputFile, spaCutoff, outlierIqrRatio, outlierAbsBound,
      nthreads, missingCutoff, minMafCutoff, minMacCutoff,
      /*tauLabels=*/{});
}


// ══════════════════════════════════════════════════════════════════════
// runSPAsqrPheno — new entry point (--pheno + --pheno-quantitative)
// ══════════════════════════════════════════════════════════════════════

void runSPAsqrPheno(
    const std::string& phenoFile,
    const std::string& covarFile,
    const std::string& quantPhenoCol,
    const std::vector<std::string>& covarNames,
    const std::vector<double>& taus,
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const GenoSpec& geno,
    const std::string& outputFile,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff
) {
  const int ntaus = static_cast<int>(taus.size());

  // ── 1. Load phenotype/covariate data ────────────────────────────────────
  infoMsg("Loading phenotype and covariate data");
  auto famIIDs = parseGenoIIDs(geno);
  SubjectData sd(std::move(famIIDs));
  if (!phenoFile.empty())  sd.loadPhenoFile(phenoFile);
  if (!covarFile.empty())  sd.loadCovar(covarFile);
  sd.finalize();
  sd.dropNaInColumns({quantPhenoCol});  // remove subjects with missing phenotype

  const Eigen::Index N = static_cast<Eigen::Index>(sd.nUsed());
  infoMsg("Subjects after intersection: %lld", static_cast<long long>(N));

  // Extract Y and X
  Eigen::VectorXd Y = sd.getColumn(quantPhenoCol);
  Eigen::MatrixXd X = covarNames.empty()
      ? Eigen::MatrixXd(N, 0)
      : sd.getColumns(covarNames);
  const int p = static_cast<int>(X.cols());
  infoMsg("Quantitative phenotype: %s, %d covariates, %d tau levels",
          quantPhenoCol.c_str(), p, ntaus);

  // ── 2. Compute bandwidth h = IQR(Y) / 3 ──────────────────────────
  {
    std::vector<double> ysorted(N);
    Eigen::VectorXd::Map(ysorted.data(), N) = Y;
    std::sort(ysorted.begin(), ysorted.end());
    auto quantile = [&](double prob) -> double {
      double idx = prob * (N - 1);
      Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
      Eigen::Index hi = std::min(lo + 1, N - 1);
      double frac = idx - lo;
      return ysorted[lo] * (1.0 - frac) + ysorted[hi] * frac;
    };
    double iqr = quantile(0.75) - quantile(0.25);
    double h = iqr / 3.0;
    if (h <= 0.0)
      h = std::max(std::pow((std::log(N) + p) / static_cast<double>(N), 0.4), 0.05);
    infoMsg("Bandwidth h = %.6f", h);

  // ── 3. Run conquer for each tau → build residual matrix ───────────
    Eigen::MatrixXd ResidMat(N, ntaus);
    for (int t = 0; t < ntaus; ++t) {
      infoMsg("  conquer tau=%.4f ...", taus[t]);
      Eigen::VectorXd resid;
      Eigen::VectorXd beta = conquer::smqrGauss(X, Y, taus[t], h, &resid);
      infoMsg("    intercept=%.6f", beta(0));

      // ResidMat column = tau - Phi(-resid / h)   (smoothed quantile residual)
      const double h1 = 1.0 / h;
      for (Eigen::Index i = 0; i < N; ++i)
        ResidMat(i, t) = taus[t] - math::pnorm(-resid(i) * h1);
    }

  // ── 4. Build tau labels ────────────────────────────────────────────
    std::vector<std::string> tauLabels;
    tauLabels.reserve(ntaus);
    for (double tau : taus) {
      std::ostringstream oss;
      oss << "Tau" << tau;
      tauLabels.push_back(oss.str());
    }

  // ── 5. Load genotype data ────────────────────────────────────────────────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(),
                                 nSnpPerChunk);

  // ── 6. Run shared pipeline ───────────────────────────────────────────────
    runSPAsqrPipeline(
        ResidMat, *genoData,
        sd.usedIIDs(), sd.famIIDs(),
        genoData->nSubjUsed(),
        spgrmGrabFile, spgrmGctaFile,
        outputFile, spaCutoff, outlierIqrRatio, outlierAbsBound,
        nthreads, missingCutoff, minMafCutoff, minMacCutoff,
        std::move(tauLabels));
  }
}
