#ifndef WTCOXG_H
#define WTCOXG_H

// WtCoxG.h -- Weighted Cox-type G-test for batch-effect-aware association

#include <RcppArmadillo.h>
#include <vector>
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
  const double m_cutoff;
  const double m_SPA_Cutoff;
  // Shared marker reference map (populated once, shared across thread copies)
  const std::shared_ptr<const std::unordered_map<uint64_t, RefInfo>> m_refMap;

  std::vector<MarkerInfo> m_markerInfoVec;
  arma::vec m_scoreVec;
  arma::vec m_zScoreVec;
  // Per-chunk reference info (rebuilt each chunk)
  std::vector<RefInfo> m_chunkRefInfo;

public:

  mtWtCoxGClass(
    arma::vec R,
    arma::vec w,
    double cutoff,
    double SPA_Cutoff,
    std::unordered_map<uint64_t, RefInfo> refMap = {}
  )
    : m_R(std::move(R)),
      m_w(std::move(w)),
      m_cutoff(cutoff),
      m_SPA_Cutoff(SPA_Cutoff),
      m_refMap(std::make_shared<const std::unordered_map<uint64_t, RefInfo>>(std::move(refMap)))
  {}

  void updateMarkerInfo(const std::vector<double>& AF_ref,
                          const std::vector<double>& AN_ref,
                          const std::vector<double>& TPR,
                          const std::vector<double>& sigma2,
                          const std::vector<double>& pvalue_bat,
                          const std::vector<double>& w_ext,
                          const std::vector<double>& var_ratio_w0,
                          const std::vector<double>& var_ratio_int,
                          const std::vector<double>& var_ratio_ext) {
    m_markerInfoVec.clear();

    size_t n = AF_ref.size();
    if (AN_ref.size() != n || TPR.size() != n || sigma2.size() != n ||
      pvalue_bat.size() != n || w_ext.size() != n || var_ratio_w0.size() != n ||
      var_ratio_int.size() != n || var_ratio_ext.size() != n) {
      throw std::runtime_error("WtCoxG::updateMarkerInfo received vectors with inconsistent lengths.");
    }

    m_markerInfoVec.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      m_markerInfoVec.emplace_back(
        AF_ref[i], AN_ref[i], TPR[i], sigma2[i], pvalue_bat[i],
        w_ext[i], var_ratio_w0[i], var_ratio_int[i], var_ratio_ext[i]
      );
    }
  }

  arma::vec getpvalVec(const arma::vec& GVec, const int i) {
    arma::vec result(2);
    const MarkerInfo& info = m_markerInfoVec[i];

    arma::vec res_ext = wtCoxGTest(
      GVec, m_R, m_w,
      info.pvalue_bat,
      info.TPR,
      info.sigma2,
      info.w_ext,
      info.var_ratio_int,
      info.var_ratio_w0,
      info.var_ratio_w0,
      info.var_ratio_ext,
      info.var_ratio_ext,
      info.AF_ref,
      info.AN_ref / 2.0,
      m_cutoff);
    result(0) = res_ext(0);

    arma::vec res_noext = wtCoxGTest(
      GVec, m_R, m_w,
      info.pvalue_bat,
      arma::datum::nan,
      arma::datum::nan,
      0.0,
      info.var_ratio_int,
      1.0,
      1.0,
      1.0,
      1.0,
      arma::datum::nan,
      arma::datum::nan,
      m_cutoff);
    result(1) = res_noext(0);

    if (m_scoreVec.n_elem < 2) {
      m_scoreVec.resize(2);
      m_zScoreVec.resize(2);
    }
    m_scoreVec(0) = res_ext(1);
    m_scoreVec(1) = res_noext(1);
    m_zScoreVec(0) = res_ext(2);
    m_zScoreVec(1) = res_noext(2);

    return result;
  }

  arma::vec getScoreVec() const {
    return m_scoreVec;
  }

  arma::vec getZScoreVec() const {
    return m_zScoreVec;
  }

  // Prepare per-chunk marker info by looking up the refMap
  void prepareChunk(const std::vector<uint64_t>& genoIndices) {
    size_t n = genoIndices.size();
    m_chunkRefInfo.resize(n);
    std::vector<double> AF_ref(n), AN_ref(n), TPR(n), sigma2(n);
    std::vector<double> pvalue_bat(n), w_ext(n), var_w0(n), var_int(n), var_ext(n);
    for (size_t i = 0; i < n; ++i) {
      RefInfo ri{NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
      if (m_refMap) {
        auto it = m_refMap->find(genoIndices[i]);
        if (it != m_refMap->end()) ri = it->second;
      }
      m_chunkRefInfo[i] = ri;
      AF_ref[i] = ri.AF_ref;   AN_ref[i] = ri.AN_ref;
      TPR[i] = ri.TPR;         sigma2[i] = ri.sigma2;
      pvalue_bat[i] = ri.pvalue_bat; w_ext[i] = ri.w_ext;
      var_w0[i] = ri.var_ratio_w0;   var_int[i] = ri.var_ratio_int;
      var_ext[i] = ri.var_ratio_ext;
    }
    updateMarkerInfo(AF_ref, AN_ref, TPR, sigma2, pvalue_bat,
                     w_ext, var_w0, var_int, var_ext);
  }

  const RefInfo& getChunkRefInfo(size_t i) const { return m_chunkRefInfo[i]; }

private:

  arma::vec wtCoxGTest(
    const arma::vec& g_input, const arma::vec& R, const arma::vec& w,
    double p_bat, double TPR, double sigma2, double b,
    double var_ratio_int, double var_ratio_w0, double var_ratio_w1,
    double var_ratio0, double var_ratio1, double mu_ext,
    double n_ext, double p_cut
  );
};

#endif
