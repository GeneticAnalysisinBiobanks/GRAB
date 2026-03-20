#ifndef LEAF_H
#define LEAF_H

// LEAF.h -- Cluster-stratified LEAF meta-analysis delegating to WtCoxGClass per cluster

#include <RcppArmadillo.h>
#include "mtWtCoxG.h"
#include <memory>
#include <unordered_map>

namespace LEAF {

class LEAFClass {
  
private:

  int m_Ncluster;
  std::vector<WtCoxG::WtCoxGClass> m_WtCoxGobj_vec;
  std::vector<arma::uvec> m_clusterIdx_vec;
  double m_cutoff;
  double m_SPA_Cutoff;

  // Per-cluster shared reference maps
  std::vector<std::shared_ptr<const std::unordered_map<uint64_t, WtCoxG::WtCoxGClass::RefInfo>>> m_refMaps;
  // Per-chunk reference info per cluster
  std::vector<std::vector<WtCoxG::WtCoxGClass::RefInfo>> m_chunkRefInfo;

public:

  LEAFClass(
    const std::vector<arma::vec>& residuals,
    const std::vector<arma::vec>& weights,
    const std::vector<arma::uvec>& clusterIdx,
    double cutoff,
    double SPA_Cutoff
  ) {
    m_Ncluster = static_cast<int>(residuals.size());
    m_cutoff = cutoff;
    m_SPA_Cutoff = SPA_Cutoff;


    m_WtCoxGobj_vec.reserve(m_Ncluster);
    m_clusterIdx_vec.resize(m_Ncluster);

    for (int i = 0; i < m_Ncluster; i++) {
      m_clusterIdx_vec[i] = clusterIdx[i];
      m_WtCoxGobj_vec.emplace_back(
        residuals[i],
        weights[i],
        m_cutoff,
        m_SPA_Cutoff
      );
    }
  }

  int get_Ncluster() const { return m_Ncluster; }

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

  // Set the per-cluster reference maps (called once before threading)
  void setRefMaps(
      std::vector<std::shared_ptr<const std::unordered_map<uint64_t, WtCoxG::WtCoxGClass::RefInfo>>> maps) {
    m_refMaps = std::move(maps);
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
        WtCoxG::WtCoxGClass::RefInfo ri{NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN,NAN};
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

  const WtCoxG::WtCoxGClass::RefInfo& getChunkRefInfo(int cluster, size_t markerIdx) const {
    return m_chunkRefInfo[cluster][markerIdx];
  }

  arma::vec getMarkerZSP(const arma::vec& GVec, int marker_idx) {
    arma::vec zScoreVec(2 * m_Ncluster);
    arma::vec scoreVec(2 * m_Ncluster);
    arma::vec pvalVec(2 * m_Ncluster);

    for (int i = 0; i < m_Ncluster; i++) {

      arma::vec GVec_cluster = GVec.elem(m_clusterIdx_vec[i]);
      arma::vec pval_cluster = m_WtCoxGobj_vec[i].getpvalVec(GVec_cluster, marker_idx);

      pvalVec(2*i) = pval_cluster(0);
      pvalVec(2*i + 1) = pval_cluster(1);

      arma::vec score_cluster = m_WtCoxGobj_vec[i].getScoreVec();
      scoreVec(2*i) = score_cluster(0);
      scoreVec(2*i + 1) = score_cluster(1);

      arma::vec zscore_cluster = m_WtCoxGobj_vec[i].getZScoreVec();
      zScoreVec(2*i) = zscore_cluster(0);
      zScoreVec(2*i + 1) = zscore_cluster(1);
    }

    arma::vec result = arma::join_cols(arma::join_cols(zScoreVec, scoreVec), pvalVec);
    return result;
  }

};

}

#endif
