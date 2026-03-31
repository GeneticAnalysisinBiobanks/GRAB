// spasqr.cpp — SPAsqr: SPA-squared multi-tau marker association (pure C++17 / Eigen / Boost)

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "spasqr/spasqr.hpp"
#include "spagrm/spagrm.hpp"
#include "engine/marker.hpp"
#include "io/resid_file.hpp"
#include "io/plink.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>


// ══════════════════════════════════════════════════════════════════════
// SPAsqrMethod — MethodBase adapter wrapping ntaus SPAGRMClass instances
// ══════════════════════════════════════════════════════════════════════

namespace {

class SPAsqrMethod : public MethodBase {
public:
  SPAsqrMethod(
      int ntaus,
      std::vector<SPAGRMClass> spagrm_vec)
    : m_ntaus(ntaus),
      m_spagrm_vec(std::move(spagrm_vec))
  {}

  std::unique_ptr<MethodBase> clone() const override {
    return std::make_unique<SPAsqrMethod>(*this);
  }

  int resultSize() const override { return 2 * m_ntaus + 1; }

  std::string getHeaderColumns() const override {
    std::ostringstream oss;
    for (int i = 1; i <= m_ntaus; ++i)
      oss << "\tZ_tau" << i;
    for (int i = 1; i <= m_ntaus; ++i)
      oss << "\tP_tau" << i;
    oss << "\tP_CCT";
    return oss.str();
  }

  void getResultVec(
      Eigen::Ref<Eigen::VectorXd> GVec,
      double altFreq,
      int /*markerInChunkIdx*/,
      bool /*flipped*/,
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

    for (int i = 0; i < m_ntaus; ++i) result.push_back(zScores[i]);
    for (int i = 0; i < m_ntaus; ++i) result.push_back(pvals[i]);

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
    result.push_back(pCCT);
  }

private:
  int m_ntaus;
  std::vector<SPAGRMClass> m_spagrm_vec;
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
    double outlierAbsBound)
{
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
// runSPAsqr — entry point
// ══════════════════════════════════════════════════════════════════════

void runSPAsqr(
    const std::string& residFile,
    const std::string& sparseGrmFile,
    const std::string& bfilePrefix,
    const std::string& outputFile,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff)
{
  // ── 1. Load residual matrix ────────────────────────────────────────
  infoMsg("Loading residual matrix from %s", residFile.c_str());
  ResidMatrixData rmd = loadResidMatrix(residFile);
  const Eigen::Index N = rmd.residuals.rows();
  const Eigen::Index K = rmd.residuals.cols();
  infoMsg("Residual matrix: %lld subjects x %lld columns",
          static_cast<long long>(N), static_cast<long long>(K));

  // ── 2. Load PLINK data ────────────────────────────────────────────
  infoMsg("Loading PLINK files: %s", bfilePrefix.c_str());
  PlinkData plinkData(
      bfilePrefix + ".bed",
      bfilePrefix + ".bim",
      bfilePrefix + ".fam",
      rmd.subjects,
      "ref-first",
      {}, {}, {}, {},
      nSnpPerChunk);

  // The PlinkData genotype vectors follow the same order as rmd.subjects.
  // So the residual matrix rows are already aligned with genotype vectors.
  const uint32_t nUsed = plinkData.nSubjUsed();
  const auto& subjOrder = rmd.subjects;

  // The residual matrix is already in the correct order (same as rmd.subjects).
  Eigen::MatrixXd& ResidMat = rmd.residuals;

  // ── 3. Outlier detection ───────────────────────────────────────────
  infoMsg("Detecting outliers (IQR ratio=%.2f, abs bound=%.2f)",
          outlierIqrRatio, outlierAbsBound);
  OutlierInfo outlierInfo = detectOutliers(ResidMat, outlierIqrRatio, outlierAbsBound);

  // ── 4. Load sparse GRM and compute variance terms ──────────────────
  infoMsg("Loading sparse GRM from %s", sparseGrmFile.c_str());

  // Parse GRM file manually to get (id1, id2, value) with factor logic.
  // We need raw entries (not symmetrized) so we can apply factor=2 for off-diag.
  struct GRMEntry {
    uint32_t row, col;
    double value;
    double factor;  // 1 for diagonal, 2 for off-diagonal
  };

  std::unordered_map<std::string, uint32_t> subjIdMap;
  subjIdMap.reserve(nUsed);
  for (uint32_t i = 0; i < nUsed; ++i)
    subjIdMap.emplace(subjOrder[i], i);

  std::vector<GRMEntry> grmEntries;
  {
    std::ifstream ifs(sparseGrmFile);
    if (!ifs) throw std::runtime_error("Cannot open sparse GRM: " + sparseGrmFile);
    std::string line;
    while (std::getline(ifs, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty() || line[0] == '#') continue;
      const char*       p   = line.c_str();
      const char* const end = p + line.size();
      auto skipWS  = [&]() { while (p < end && (*p == ' ' || *p == '\t')) ++p; };
      auto nextTok = [&]() -> std::string {
        skipWS();
        const char* s = p;
        while (p < end && *p != ' ' && *p != '\t') ++p;
        return std::string(s, p);
      };
      skipWS();
      if (p >= end) continue;
      std::string id1 = nextTok();
      std::string id2 = nextTok();
      skipWS();
      if (p >= end) continue;
      char* endPtr;
      double val = std::strtod(p, &endPtr);
      if (endPtr == p) continue;

      auto it1 = subjIdMap.find(id1);
      auto it2 = subjIdMap.find(id2);
      if (it1 == subjIdMap.end() || it2 == subjIdMap.end()) continue;

      grmEntries.push_back({it1->second, it2->second, val,
                            (id1 == id2) ? 1.0 : 2.0});
    }
  }
  infoMsg("Sparse GRM: %zu entries after filtering", grmEntries.size());

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
  SPAsqrMethod method(ntaus, std::move(spagrm_vec));

  infoMsg("Starting marker-level association (SPAsqr, %d taus, %d threads)",
          ntaus, nthreads);
  markerEngine(plinkData, method, outputFile,
               nthreads, missingCutoff, minMafCutoff, minMacCutoff,
               /*exactHwe=*/false);
}
