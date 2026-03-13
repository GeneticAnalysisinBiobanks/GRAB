#ifndef LEAF_H
#define LEAF_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "WtCoxG.h"

namespace LEAF {

class LEAFClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  int m_Ncluster;                                    // Number of clusters
  std::vector<WtCoxG::WtCoxGClass> m_WtCoxGobj_vec; // Vector of pre-built WtCoxG objects (one per cluster)
  std::vector<arma::uvec> m_clusterIdx_vec;          // Vector of cluster indices (0-based)
  double m_cutoff;                                   // Batch effect p-value cutoff
  double m_SPA_Cutoff;                               // SPA cutoff
  
public:
  
  LEAFClass(
    const std::vector<arma::vec>& t_residuals,  // Residual vectors (one per cluster)
    const std::vector<arma::vec>& t_weights,    // Weight vectors (one per cluster)
    const std::vector<arma::uvec>& t_clusterIdx, // Cluster index vectors (0-based)
    double t_cutoff,               // Batch effect p-value cutoff
    double t_SPA_Cutoff            // SPA cutoff
  ) {
    m_Ncluster = static_cast<int>(t_residuals.size());
    m_cutoff = t_cutoff;
    m_SPA_Cutoff = t_SPA_Cutoff;

    // Create Ncluster WtCoxG objects once
    m_WtCoxGobj_vec.reserve(m_Ncluster);
    m_clusterIdx_vec.resize(m_Ncluster);

    for (int i = 0; i < m_Ncluster; i++) {
      m_clusterIdx_vec[i] = t_clusterIdx[i];

      // Create WtCoxG object for this cluster and store it
      m_WtCoxGobj_vec.emplace_back(
        t_residuals[i],
        t_weights[i],
        m_cutoff,
        m_SPA_Cutoff
      );
    }
  }
  
  int get_Ncluster() const { return m_Ncluster; }
  
  // Get cluster indices for a specific cluster (0-based)
  arma::uvec get_clusterIdx(int cluster_idx) const {
    if (cluster_idx >= 0 && cluster_idx < m_Ncluster) {
      return m_clusterIdx_vec[cluster_idx];
    }
    return arma::uvec();
  }

  // Update marker info for a specific cluster
  void updateMarkerInfo(int cluster_idx, const Rcpp::DataFrame& t_mergeGenoInfo) {
    if (cluster_idx >= 0 && cluster_idx < m_Ncluster) {
      m_WtCoxGobj_vec[cluster_idx].updateMarkerInfo(t_mergeGenoInfo);
    }
  }

  // Thread-safe overload: update marker info for a specific cluster using native vectors.
  void updateMarkerInfo(int cluster_idx,
                        const std::vector<double>& AF_ref,
                        const std::vector<double>& AN_ref,
                        const std::vector<double>& TPR,
                        const std::vector<double>& sigma2,
                        const std::vector<double>& pvalue_bat,
                        const std::vector<double>& w_ext,
                        const std::vector<double>& var_ratio_w0,
                        const std::vector<double>& var_ratio_int,
                        const std::vector<double>& var_ratio_ext) {
    if (cluster_idx >= 0 && cluster_idx < m_Ncluster) {
      m_WtCoxGobj_vec[cluster_idx].updateMarkerInfo(
        AF_ref, AN_ref, TPR, sigma2, pvalue_bat,
        w_ext, var_ratio_w0, var_ratio_int, var_ratio_ext
      );
    }
  }

  // Get z-scores, scores, and p-values for all clusters in a single pass
  // Returns concatenated vectors: [zScores, Scores, Pvals]
  // Each sub-vector has length 2*Ncluster: [val.ext.clus1, val.noext.clus1, val.ext.clus2, val.noext.clus2, ...]
  arma::vec getMarkerZSP(const arma::vec& GVec, int marker_idx) {
    arma::vec zScoreVec(2 * m_Ncluster);
    arma::vec scoreVec(2 * m_Ncluster);
    arma::vec pvalVec(2 * m_Ncluster);
    
    for (int i = 0; i < m_Ncluster; i++) {
      // Extract cluster-specific genotypes using stored clusterIdx
      arma::vec GVec_cluster = GVec.elem(m_clusterIdx_vec[i]);
      
      // Get p-values (this computes and stores scores and z-scores internally)
      arma::vec pval_cluster = m_WtCoxGobj_vec[i].getpvalVec(GVec_cluster, marker_idx);
      pvalVec(2*i) = pval_cluster(0);     // p.ext for cluster i
      pvalVec(2*i + 1) = pval_cluster(1); // p.noext for cluster i
      
      // Retrieve pre-computed scores
      arma::vec score_cluster = m_WtCoxGobj_vec[i].getScoreVec();
      scoreVec(2*i) = score_cluster(0);     // score.ext for cluster i
      scoreVec(2*i + 1) = score_cluster(1); // score.noext for cluster i
      
      // Retrieve pre-computed z-scores
      arma::vec zscore_cluster = m_WtCoxGobj_vec[i].getZScoreVec();
      zScoreVec(2*i) = zscore_cluster(0);     // zscore.ext for cluster i
      zScoreVec(2*i + 1) = zscore_cluster(1); // zscore.noext for cluster i
    }
    
    // Concatenate: [zScores, Scores, Pvals]
    arma::vec result = arma::join_cols(arma::join_cols(zScoreVec, scoreVec), pvalVec);
    return result;
  }

};

} // namespace LEAF

#endif
