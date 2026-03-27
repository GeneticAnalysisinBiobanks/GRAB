#ifndef WTCOXG_H
#define WTCOXG_H

// WtCoxG.h -- Weighted Cox-type G-test for batch-effect-aware association

#include <RcppArmadillo.h>
#include <vector>
#include <array>
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdexcept>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <unordered_map>
#include <memory>


class mtWtCoxGClass {
public:
  // Per-marker external-reference info (shared via refMap).
  struct RefInfo {
    double AF_ref, AN_ref, TPR, sigma2, pvalue_bat;
    double w_ext, var_ratio_w0, var_ratio_int, var_ratio_ext;
  };

private:
  // Per-marker metadata used during chunk computation.
  struct MarkerInfo {
    double AF_ref;
    double AN_ref;
    double TPR;
    double sigma2;
    double pvalue_bat;
    double w_ext;
    double var_ratio_w0;
    double var_ratio_int;
    double var_ratio_ext;

    MarkerInfo() = default;

    MarkerInfo(double af_ref, double an_ref, double tpr, double sig2,
        double pval_bat, double w_ext_val, double var_ratio_w0_val,
        double var_ratio_int_val, double var_ratio_ext_val)
      : AF_ref(af_ref), AN_ref(an_ref), TPR(tpr), sigma2(sig2),
      pvalue_bat(pval_bat), w_ext(w_ext_val), var_ratio_w0(var_ratio_w0_val),
      var_ratio_int(var_ratio_int_val), var_ratio_ext(var_ratio_ext_val) {}
  };

  const arma::vec m_R;
  const arma::vec m_w;
  const arma::vec m_w1;      // w / (2 * sum(w)), constant across markers
  const double m_meanR;      // mean(R), constant
  const double m_sumR;       // sum(R), constant
  const double m_cutoff;
  const double m_SPA_Cutoff;
  // Shared marker reference map (populated once, shared across thread copies)
  const std::shared_ptr<const std::unordered_map<uint64_t, RefInfo>> m_refMap;

  std::vector<MarkerInfo> m_markerInfoVec;
  // Per-marker scratch (avoid per-marker heap allocation)
  std::array<double,2> m_scoreArr;
  std::array<double,2> m_zScoreArr;
  // Per-chunk reference info (rebuilt each chunk)
  std::vector<RefInfo> m_chunkRefInfo;

public:

  mtWtCoxGClass(
    arma::vec R,
    arma::vec w,
    double cutoff,
    double SPA_Cutoff,
    std::unordered_map<uint64_t, RefInfo> refMap = {}
  );

  void updateMarkerInfo(
    const std::vector<double>& AF_ref,
    const std::vector<double>& AN_ref,
    const std::vector<double>& TPR,
    const std::vector<double>& sigma2,
    const std::vector<double>& pvalue_bat,
    const std::vector<double>& w_ext,
    const std::vector<double>& var_ratio_w0,
    const std::vector<double>& var_ratio_int,
    const std::vector<double>& var_ratio_ext
  );

  struct WtResult { double pval; double score; double zscore; };

  void getpvalVec(const arma::vec& GVec, const int i, double& pExt, double& pNoext);

  const std::array<double,2>& getScoreArr()  const { return m_scoreArr;  }
  const std::array<double,2>& getZScoreArr() const { return m_zScoreArr; }

  void prepareChunk(const std::vector<uint64_t>& genoIndices);

  const RefInfo& getChunkRefInfo(size_t i) const { return m_chunkRefInfo[i]; }

  // Fills rv with [pExt, pNoext, zExt, zNoext, AF_ref, AN_ref, pvalue_bat]
  void getResultVec(const arma::vec& GVec, int markerIdx, std::vector<double>& rv) {
    double pExt, pNoext;
    getpvalVec(GVec, markerIdx, pExt, pNoext);
    const auto& zArr = getZScoreArr();
    const RefInfo& ri = getChunkRefInfo(markerIdx);
    rv.clear();
    rv.push_back(pExt); rv.push_back(pNoext);
    rv.push_back(zArr[0]); rv.push_back(zArr[1]);
    rv.push_back(ri.AF_ref); rv.push_back(ri.AN_ref); rv.push_back(ri.pvalue_bat);
  }

  static int resultSize() { return 7; }

  std::string getHeaderColumns() const {
    return "\tWtCoxG.ext\tWtCoxG.noext\tzscore.ext\tzscore.noext\tAF_ref\tAN_ref\tpvalue_bat";
  }

private:

  WtResult wtCoxGTest(
    const arma::vec& g_input, const arma::vec& R, const arma::vec& w,
    double p_bat, double TPR, double sigma2, double b,
    double var_ratio_int, double var_ratio_w0, double var_ratio_w1,
    double var_ratio0, double var_ratio1, double mu_ext,
    double n_ext, double p_cut
  );
};

#endif
