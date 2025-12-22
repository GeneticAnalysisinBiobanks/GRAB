#ifndef SPAsqr_H
#define SPAsqr_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "SPAGRM.h"

namespace SPAsqr{

class SPAsqrClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  arma::vec m_taus;                      // quantiles
  arma::mat m_Resid_mat;                  // residual matrix (N × ntaus)
  Rcpp::List m_Resid_unrelated_outliers_lst;  // list of unrelated outlier residual vectors (one per tau)
  arma::vec m_sum_R_nonOutlier_vec;      // sum of non-outlier residuals (length ntaus)
  arma::vec m_R_GRM_R_nonOutlier_vec;    // residuals x GRM x residuals for non-outlier families (length ntaus)
  arma::vec m_R_GRM_R_TwoSubjOutlier_vec;// residuals x GRM x residuals for outlier families n=2 (length ntaus)
  arma::vec m_R_GRM_R_vec;               // residuals x GRM x residuals (length ntaus)
  arma::vec m_MAF_interval;              // MAF interval divides the MAFs into several intervals
  Rcpp::List m_TwoSubj_list_lst;         // List of residuals (matrix 2 × ntaus) and IBD in outlier families (n = 2)
  Rcpp::List m_CLT_union_lst;            // Shared CLT cache (union across all taus)
  Rcpp::List m_ThreeSubj_family_idx_lst; // Family indices per tau (indices into CLT_union_lst)
  Rcpp::List m_ThreeSubj_stand_S_lst;    // stand.S values per tau
  
  double m_SPA_Cutoff;                   // cutoff of standardized score to use normal approximation or SPA
  double m_zeta;                         // initial saddle point for negative side, default is zero
  double m_tol;                          // accuracy of Newton's methods, default 1e-4 for beta; 1e-5 for tau 
  
public:
  
  SPAsqrClass(
    arma::vec t_taus,
    arma::mat t_Resid_mat,
    Rcpp::List t_Resid_unrelated_outliers_lst,
    arma::vec t_sum_R_nonOutlier_vec,
    arma::vec t_R_GRM_R_nonOutlier_vec,
    arma::vec t_R_GRM_R_TwoSubjOutlier_vec,
    arma::vec t_R_GRM_R_vec,
    arma::vec t_MAF_interval,
    Rcpp::List t_TwoSubj_list_lst,
    Rcpp::List t_CLT_union_lst,
    Rcpp::List t_ThreeSubj_family_idx_lst,
    Rcpp::List t_ThreeSubj_stand_S_lst,
    double t_SPA_Cutoff,
    double t_zeta,
    double t_tol
  ) {
    m_taus = t_taus;
    m_Resid_mat = t_Resid_mat;
    m_Resid_unrelated_outliers_lst = t_Resid_unrelated_outliers_lst;
    m_sum_R_nonOutlier_vec = t_sum_R_nonOutlier_vec;
    m_R_GRM_R_nonOutlier_vec = t_R_GRM_R_nonOutlier_vec;
    m_R_GRM_R_TwoSubjOutlier_vec = t_R_GRM_R_TwoSubjOutlier_vec;
    m_R_GRM_R_vec = t_R_GRM_R_vec;
    m_MAF_interval = t_MAF_interval;
    m_TwoSubj_list_lst = t_TwoSubj_list_lst;
    m_CLT_union_lst = t_CLT_union_lst;
    m_ThreeSubj_family_idx_lst = t_ThreeSubj_family_idx_lst;
    m_ThreeSubj_stand_S_lst = t_ThreeSubj_stand_S_lst;
    m_SPA_Cutoff = t_SPA_Cutoff;
    m_zeta = t_zeta;
    m_tol = t_tol;
  }
  
  int get_ntaus() const { return m_taus.n_elem; }

  arma::vec getMarkerPval(
    arma::vec GVec,
    double altFreq,
    arma::vec& zScoreVec,
    double& hwepval
  ) {
    int ntaus = m_taus.n_elem;
    arma::vec pvalVec(ntaus);
    zScoreVec.set_size(ntaus);
    
    for (int i = 0; i < ntaus; i++) {

      // Extract tau-specific data from vectors (i-th element)
      arma::vec resid_i = m_Resid_mat.col(i);
      
      arma::vec resid_unrelated_outliers_i;
      Rcpp::RObject resid_unrelated_outliers_obj = m_Resid_unrelated_outliers_lst[i];
      if (!resid_unrelated_outliers_obj.isNULL()) {
        resid_unrelated_outliers_i = Rcpp::as<arma::vec>(resid_unrelated_outliers_obj);
      }

      double sum_R_nonOutlier_i = m_sum_R_nonOutlier_vec(i);
      double R_GRM_R_nonOutlier_i = m_R_GRM_R_nonOutlier_vec(i);
      double R_GRM_R_TwoSubjOutlier_i = m_R_GRM_R_TwoSubjOutlier_vec(i);
      double R_GRM_R_i = m_R_GRM_R_vec(i);
      
      // Extract i-th element from TwoSubj_list_lst (which is a list of lists)
      Rcpp::List TwoSubj_list_i = Rcpp::as<Rcpp::List>(m_TwoSubj_list_lst[i]);
      
      // Reconstruct ThreeSubj_list_i from separated components
      Rcpp::List ThreeSubj_list_i;
      Rcpp::IntegerVector family_idx_i = Rcpp::as<Rcpp::IntegerVector>(m_ThreeSubj_family_idx_lst[i]);
      Rcpp::List stand_S_lst_i = Rcpp::as<Rcpp::List>(m_ThreeSubj_stand_S_lst[i]);
      
      int n_families_i = family_idx_i.size();
      ThreeSubj_list_i = Rcpp::List(n_families_i);
      
      for (int j = 0; j < n_families_i; j++) {
        int clt_idx = family_idx_i[j] - 1;  // R indices are 1-based, C++ is 0-based
        Rcpp::List family_j = Rcpp::List::create(
          Rcpp::Named("CLT") = m_CLT_union_lst[clt_idx],
          Rcpp::Named("stand.S") = stand_S_lst_i[j]
        );
        ThreeSubj_list_i[j] = family_j;
      }
      
      // Create temporary SPAGRM object for this tau
      SPAGRM::SPAGRMClass tempClass(
        resid_i,
        resid_unrelated_outliers_i,
        sum_R_nonOutlier_i,
        R_GRM_R_nonOutlier_i,
        R_GRM_R_TwoSubjOutlier_i,
        R_GRM_R_i,
        m_MAF_interval,
        TwoSubj_list_i,
        ThreeSubj_list_i,
        m_SPA_Cutoff,
        m_zeta,
        m_tol
      );
      
      // Get p-value for this tau
      double zScore_i, hwepval_i;
      double pval_i = tempClass.getMarkerPval(GVec, altFreq, zScore_i, hwepval_i);
      
      pvalVec(i) = pval_i;
      zScoreVec(i) = zScore_i;
      if (i == 0) hwepval = hwepval_i; // hwepval is same for all taus
    }
    
    return pvalVec;
  }

};
}

#endif