#ifndef LEAF_H
#define LEAF_H

// LEAF.h -- Cluster-stratified LEAF meta-analysis delegating to mtWtCoxGClass per cluster

#include <RcppArmadillo.h>
#include <sstream>
#include <boost/math/distributions/chi_squared.hpp>
#include "mtWtCoxG.h"
#include <memory>
#include <unordered_map>

class mtLEAFClass {
  
private:

  const int m_Ncluster;
  std::vector<mtWtCoxGClass> m_WtCoxGobj_vec;
  std::vector<arma::uvec> m_clusterIdx_vec;
  const double m_cutoff;
  const double m_SPA_Cutoff;

  // Per-cluster shared reference maps
  const std::vector<std::shared_ptr<const std::unordered_map<uint64_t, mtWtCoxGClass::RefInfo>>> m_refMaps;
  // Per-chunk reference info per cluster
  std::vector<std::vector<mtWtCoxGClass::RefInfo>> m_chunkRefInfo;

  arma::uvec get_clusterIdx(int cluster_idx) const {
    if (cluster_idx >= 0 && cluster_idx < m_Ncluster) {
      return m_clusterIdx_vec[cluster_idx];
    }
    return arma::uvec();
  }

  void updateMarkerInfo(
    int cluster_idx,
    const std::vector<double>& AF_ref,
    const std::vector<double>& AN_ref,
    const std::vector<double>& TPR,
    const std::vector<double>& sigma2,
    const std::vector<double>& pvalue_bat,
    const std::vector<double>& w_ext,
    const std::vector<double>& var_ratio_w0,
    const std::vector<double>& var_ratio_int,
    const std::vector<double>& var_ratio_ext
  ) {
    if (cluster_idx >= 0 && cluster_idx < m_Ncluster) {
      m_WtCoxGobj_vec[cluster_idx].updateMarkerInfo(
        AF_ref, AN_ref, TPR, sigma2, pvalue_bat,
        w_ext, var_ratio_w0, var_ratio_int, var_ratio_ext
      );
    }
  }

public:

  mtLEAFClass(
    std::vector<arma::vec> residuals,
    std::vector<arma::vec> weights,
    std::vector<arma::uvec> clusterIdx,
    double cutoff,
    double SPA_Cutoff,
    std::vector<std::shared_ptr<const std::unordered_map<uint64_t, mtWtCoxGClass::RefInfo>>> refMaps
  )
    : m_Ncluster(static_cast<int>(residuals.size())),
      m_cutoff(cutoff),
      m_SPA_Cutoff(SPA_Cutoff),
      m_refMaps(std::move(refMaps))
  {
    m_WtCoxGobj_vec.reserve(m_Ncluster);
    m_clusterIdx_vec.resize(m_Ncluster);

    for (int i = 0; i < m_Ncluster; i++) {
      m_clusterIdx_vec[i] = std::move(clusterIdx[i]);
      m_WtCoxGobj_vec.emplace_back(
        std::move(residuals[i]),
        std::move(weights[i]),
        m_cutoff,
        m_SPA_Cutoff
      );
    }
  }

  int get_Ncluster() const { return m_Ncluster; }

  std::string getHeaderColumns() const {
    std::ostringstream oss;
    oss << "\tmeta.p_ext\tmeta.p_noext";
    for (int i = 1; i <= m_Ncluster; ++i)
      oss << "\tcl" << i << ".p_ext\tcl" << i << ".p_noext\tcl" << i << ".p_batch";
    return oss.str();
  }

  // Prepare all clusters for a chunk
  void prepareChunk(const std::vector<uint64_t>& genoIndices) {
    m_chunkRefInfo.resize(m_Ncluster);
    for (int c = 0; c < m_Ncluster; ++c) {
      size_t n = genoIndices.size();
      m_chunkRefInfo[c].resize(n);
      std::vector<double> AF_ref(n), AN_ref(n), TPR(n), sigma2(n);
      std::vector<double> pvalue_bat(n), w_ext(n), var_w0(n), var_int(n), var_ext(n);
      for (size_t i = 0; i < n; ++i) {
        mtWtCoxGClass::RefInfo ri{NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
        if (c < static_cast<int>(m_refMaps.size()) && m_refMaps[c]) {
          auto it = m_refMaps[c]->find(genoIndices[i]);
          if (it != m_refMaps[c]->end()) ri = it->second;
        }
        m_chunkRefInfo[c][i] = ri;
        AF_ref[i] = ri.AF_ref;   AN_ref[i] = ri.AN_ref;
        TPR[i] = ri.TPR;         sigma2[i] = ri.sigma2;
        pvalue_bat[i] = ri.pvalue_bat; w_ext[i] = ri.w_ext;
        var_w0[i] = ri.var_ratio_w0;   var_int[i] = ri.var_ratio_int;
        var_ext[i] = ri.var_ratio_ext;
      }
      updateMarkerInfo(c, AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext);
    }
  }

  const mtWtCoxGClass::RefInfo& getChunkRefInfo(int cluster, size_t markerIdx) const {
    return m_chunkRefInfo[cluster][markerIdx];
  }

  // Fills rv with [meta_p_ext, meta_p_noext, cl1_p_ext, cl1_p_noext, cl1_p_batch, ...]
  void getResultVec(const arma::vec& GVec, int markerIdx, std::vector<double>& rv) {
    // Per-cluster scores and pvals (stack-allocated for small cluster counts)
    std::vector<double> sExt(m_Ncluster), sNoext(m_Ncluster);
    std::vector<double> pExt(m_Ncluster), pNoext(m_Ncluster);

    for (int i = 0; i < m_Ncluster; i++) {
      arma::vec GVec_cluster = GVec.elem(m_clusterIdx_vec[i]);
      double pE, pN;
      m_WtCoxGobj_vec[i].getpvalVec(GVec_cluster, markerIdx, pE, pN);

      pExt[i] = pE;
      pNoext[i] = pN;

      const auto& scoreArr = m_WtCoxGobj_vec[i].getScoreArr();
      sExt[i] = scoreArr[0];
      sNoext[i] = scoreArr[1];
    }

    auto metaP = [](const std::vector<double>& scores, const std::vector<double>& pvals) {
      double sumScore = 0, sumVar = 0;
      for (size_t c = 0; c < scores.size(); ++c) {
        if (std::isnan(scores[c]) || std::isnan(pvals[c]) || pvals[c] <= 0 || pvals[c] >= 1) continue;
        double chisq = boost::math::quantile(boost::math::chi_squared(1.0), 1.0 - pvals[c]);
        if (chisq < 1e-30) chisq = 1e-30;
        double var = (scores[c] * scores[c]) / chisq;
        if (std::isnan(var)) var = 0;
        sumScore += scores[c]; sumVar += var;
      }
      if (sumVar <= 0) return std::numeric_limits<double>::quiet_NaN();
      double z = sumScore / std::sqrt(sumVar);
      return std::erfc(std::fabs(z) / std::sqrt(2.0));
    };

    rv.clear();
    rv.reserve(2 + 3 * m_Ncluster);
    rv.push_back(metaP(sExt, pExt));
    rv.push_back(metaP(sNoext, pNoext));
    for (int c = 0; c < m_Ncluster; ++c) {
      rv.push_back(pExt[c]);
      rv.push_back(pNoext[c]);
      rv.push_back(getChunkRefInfo(c, markerIdx).pvalue_bat);
    }
  }

  int resultSize() const { return 2 + 3 * m_Ncluster; }

};

#endif
