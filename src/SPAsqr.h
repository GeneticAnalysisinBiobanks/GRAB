#ifndef SPAsqr_H
#define SPAsqr_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "SPAGRM.h"

namespace SPAsqr{

struct TauFamilyNativeInput {
  arma::vec resid_unrelated_outliers;
  std::vector<arma::vec> twoSubj_resid;
  std::vector<arma::vec> twoSubj_rho;
  std::vector<arma::vec> threeSubj_standS;
  std::vector<arma::mat> threeSubj_CLT;
};

class SPAsqrClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  arma::vec m_taus;                      // quantiles
  std::vector<SPAGRM::SPAGRMClass> m_SPAGRMobj_vec;  // vector of pre-built SPAGRM objects (one per tau)
  arma::vec m_MAF_interval;              // MAF interval divides the MAFs into several intervals
  
  double m_SPA_Cutoff;                   // cutoff of standardized score to use normal approximation or SPA
  double m_zeta;                         // initial saddle point for negative side, default is zero
  double m_tol;                          // accuracy of Newton's methods, default 1e-4 for beta; 1e-5 for tau 
  
public:

  static std::vector<TauFamilyNativeInput> buildTauFamilyDataFromR(
    const Rcpp::List& t_Resid_unrelated_outliers_lst,
    const Rcpp::List& t_TwoSubj_list_lst,
    const Rcpp::List& t_CLT_union_lst,
    const Rcpp::List& t_ThreeSubj_family_idx_lst,
    const Rcpp::List& t_ThreeSubj_stand_S_lst,
    int ntaus
  ) {
    std::vector<TauFamilyNativeInput> out(ntaus);
    for (int i = 0; i < ntaus; ++i) {
      Rcpp::RObject resid_unrelated_outliers_obj = t_Resid_unrelated_outliers_lst[i];
      if (!resid_unrelated_outliers_obj.isNULL()) {
        out[i].resid_unrelated_outliers = Rcpp::as<arma::vec>(resid_unrelated_outliers_obj);
      }

      Rcpp::List TwoSubj_list_i = Rcpp::as<Rcpp::List>(t_TwoSubj_list_lst[i]);
      out[i].twoSubj_resid.reserve(TwoSubj_list_i.size());
      out[i].twoSubj_rho.reserve(TwoSubj_list_i.size());
      for (int j = 0; j < TwoSubj_list_i.size(); ++j) {
        Rcpp::List pair = Rcpp::as<Rcpp::List>(TwoSubj_list_i[j]);
        out[i].twoSubj_resid.push_back(Rcpp::as<arma::vec>(pair["Resid"]));
        out[i].twoSubj_rho.push_back(Rcpp::as<arma::vec>(pair["Rho"]));
      }

      Rcpp::IntegerVector family_idx_i = Rcpp::as<Rcpp::IntegerVector>(t_ThreeSubj_family_idx_lst[i]);
      Rcpp::List stand_S_lst_i = Rcpp::as<Rcpp::List>(t_ThreeSubj_stand_S_lst[i]);
      int n_families_i = family_idx_i.size();
      out[i].threeSubj_standS.reserve(n_families_i);
      out[i].threeSubj_CLT.reserve(n_families_i);
      for (int j = 0; j < n_families_i; ++j) {
        int clt_idx = family_idx_i[j] - 1;
        out[i].threeSubj_CLT.push_back(Rcpp::as<arma::mat>(t_CLT_union_lst[clt_idx]));
        out[i].threeSubj_standS.push_back(Rcpp::as<arma::vec>(stand_S_lst_i[j]));
      }
    }
    return out;
  }

  SPAsqrClass(
    arma::vec t_taus,
    arma::mat t_Resid_mat,
    const std::vector<TauFamilyNativeInput>& t_tauFamilyData,
    arma::vec t_sum_R_nonOutlier_vec,
    arma::vec t_R_GRM_R_nonOutlier_vec,
    arma::vec t_R_GRM_R_TwoSubjOutlier_vec,
    arma::vec t_R_GRM_R_vec,
    arma::vec t_MAF_interval,
    double t_SPA_Cutoff,
    double t_zeta,
    double t_tol
  ) {
    m_taus = t_taus;
    m_MAF_interval = t_MAF_interval;
    m_SPA_Cutoff = t_SPA_Cutoff;
    m_zeta = t_zeta;
    m_tol = t_tol;

    int ntaus = m_taus.n_elem;
    m_SPAGRMobj_vec.reserve(ntaus);

    for (int i = 0; i < ntaus; ++i) {
      arma::vec resid_i = t_Resid_mat.col(i);
      const TauFamilyNativeInput& fam = t_tauFamilyData[i];

      m_SPAGRMobj_vec.emplace_back(
        resid_i,
        fam.resid_unrelated_outliers,
        t_sum_R_nonOutlier_vec(i),
        t_R_GRM_R_nonOutlier_vec(i),
        t_R_GRM_R_TwoSubjOutlier_vec(i),
        t_R_GRM_R_vec(i),
        m_MAF_interval,
        fam.twoSubj_resid,
        fam.twoSubj_rho,
        fam.threeSubj_standS,
        fam.threeSubj_CLT,
        m_SPA_Cutoff,
        m_zeta,
        m_tol
      );
    }
  }
  
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
  ) : SPAsqrClass(
        t_taus,
        t_Resid_mat,
        buildTauFamilyDataFromR(
          t_Resid_unrelated_outliers_lst,
          t_TwoSubj_list_lst,
          t_CLT_union_lst,
          t_ThreeSubj_family_idx_lst,
          t_ThreeSubj_stand_S_lst,
          t_taus.n_elem
        ),
        t_sum_R_nonOutlier_vec,
        t_R_GRM_R_nonOutlier_vec,
        t_R_GRM_R_TwoSubjOutlier_vec,
        t_R_GRM_R_vec,
        t_MAF_interval,
        t_SPA_Cutoff,
        t_zeta,
        t_tol
      ) {}
  
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
      double pval_i = m_SPAGRMobj_vec[i].getMarkerPval(GVec, altFreq, zScore_i, hwepval_i);
      pvalVec(i) = pval_i;
      zScoreVec(i) = zScore_i;
      if (i == 0) hwepval = hwepval_i; // hwepval is same for all taus
    }
    
    return pvalVec;
  }

};
}

#endif