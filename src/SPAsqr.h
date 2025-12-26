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
  std::vector<SPAGRM::SPAGRMClass*> m_SPAGRMobj_vec;  // vector of pre-built SPAGRM objects (one per tau)
  arma::vec m_MAF_interval;              // MAF interval divides the MAFs into several intervals
  
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
    m_MAF_interval = t_MAF_interval;
    m_SPA_Cutoff = t_SPA_Cutoff;
    m_zeta = t_zeta;
    m_tol = t_tol;
    
    // Create ntaus SPAGRM objects once
    int ntaus = m_taus.n_elem;
    m_SPAGRMobj_vec.resize(ntaus);
    
    for (int i = 0; i < ntaus; i++) {
      // Extract tau-specific data from vectors (i-th element)
      arma::vec resid_i = t_Resid_mat.col(i);
      
      arma::vec resid_unrelated_outliers_i;
      Rcpp::RObject resid_unrelated_outliers_obj = t_Resid_unrelated_outliers_lst[i];
      if (!resid_unrelated_outliers_obj.isNULL()) {
        resid_unrelated_outliers_i = Rcpp::as<arma::vec>(resid_unrelated_outliers_obj);
      }

      double sum_R_nonOutlier_i = t_sum_R_nonOutlier_vec(i);
      double R_GRM_R_nonOutlier_i = t_R_GRM_R_nonOutlier_vec(i);
      double R_GRM_R_TwoSubjOutlier_i = t_R_GRM_R_TwoSubjOutlier_vec(i);
      double R_GRM_R_i = t_R_GRM_R_vec(i);
      
      // Extract i-th element from TwoSubj_list_lst (which is a list of lists)
      Rcpp::List TwoSubj_list_i = Rcpp::as<Rcpp::List>(t_TwoSubj_list_lst[i]);
      
      // Reconstruct ThreeSubj_list_i from separated components
      Rcpp::List ThreeSubj_list_i;
      Rcpp::IntegerVector family_idx_i = Rcpp::as<Rcpp::IntegerVector>(t_ThreeSubj_family_idx_lst[i]);
      Rcpp::List stand_S_lst_i = Rcpp::as<Rcpp::List>(t_ThreeSubj_stand_S_lst[i]);
      
      int n_families_i = family_idx_i.size();
      ThreeSubj_list_i = Rcpp::List(n_families_i);
      
      for (int j = 0; j < n_families_i; j++) {
        int clt_idx = family_idx_i[j] - 1;  // R indices are 1-based, C++ is 0-based
        Rcpp::List family_j = Rcpp::List::create(
          Rcpp::Named("CLT") = t_CLT_union_lst[clt_idx],
          Rcpp::Named("stand.S") = stand_S_lst_i[j]
        );
        ThreeSubj_list_i[j] = family_j;
      }
      
      // Create SPAGRM object for this tau and store it
      m_SPAGRMobj_vec[i] = new SPAGRM::SPAGRMClass(
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
    }
  }
  
  // Destructor to clean up SPAGRM objects
  ~SPAsqrClass() {
    for (size_t i = 0; i < m_SPAGRMobj_vec.size(); i++) {
      if (m_SPAGRMobj_vec[i] != nullptr) {
        delete m_SPAGRMobj_vec[i];
        m_SPAGRMobj_vec[i] = nullptr;
      }
    }
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
      // Use pre-built SPAGRM object for this tau
      double zScore_i, hwepval_i;
      double pval_i = m_SPAGRMobj_vec[i]->getMarkerPval(GVec, altFreq, zScore_i, hwepval_i);
      
      pvalVec(i) = pval_i;
      zScoreVec(i) = zScore_i;
      if (i == 0) hwepval = hwepval_i; // hwepval is same for all taus
    }
    
    return pvalVec;
  }

};
}

#endif