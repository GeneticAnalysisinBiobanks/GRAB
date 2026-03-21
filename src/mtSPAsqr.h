#ifndef SPAsqr_H
#define SPAsqr_H

// SPAsqr.h -- SPA-squared: multi-tau wrapper delegating to SPAGRMClass per tau

#include <RcppArmadillo.h>
#include "mtSPAGRM.h"

// Per-tau family data passed from R at null-model construction time.

class mtSPAsqrClass {
private:

  const arma::vec m_taus;
  std::vector<mtSPAGRMClass> m_SPAGRMobj_vec;
  const arma::vec m_MAF_interval;

  const double m_SPA_Cutoff;
  const double m_zeta;
  const double m_tol;

public:

  mtSPAsqrClass(
    arma::vec taus,
    arma::mat Resid_mat,
    std::vector<SPAGRMSpace::FamilyData> tauFamilyData,
    arma::vec sum_R_nonOutlier_vec,
    arma::vec R_GRM_R_nonOutlier_vec,
    arma::vec R_GRM_R_TwoSubjOutlier_vec,
    arma::vec R_GRM_R_vec,
    arma::vec MAF_interval,
    double SPA_Cutoff,
    double zeta,
    double tol
  )
    : m_taus(std::move(taus)),
      m_MAF_interval(std::move(MAF_interval)),
      m_SPA_Cutoff(SPA_Cutoff),
      m_zeta(zeta),
      m_tol(tol)
  {
    int ntaus = m_taus.n_elem;
    m_SPAGRMobj_vec.reserve(ntaus);

    for (int i = 0; i < ntaus; ++i) {
      m_SPAGRMobj_vec.emplace_back(
        Resid_mat.col(i),
        sum_R_nonOutlier_vec(i),
        R_GRM_R_nonOutlier_vec(i),
        R_GRM_R_TwoSubjOutlier_vec(i),
        R_GRM_R_vec(i),
        m_MAF_interval,
        std::move(tauFamilyData[i]),
        m_SPA_Cutoff,
        m_zeta,
        m_tol
      );
    }
  }

  int get_ntaus() const { return m_taus.n_elem; }
  std::vector<double> getTaus() const { return arma::conv_to<std::vector<double>>::from(m_taus); }

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

      double zScore_i, hwepval_i;
      double pval_i = m_SPAGRMobj_vec[i].getMarkerPval(GVec, altFreq, zScore_i, hwepval_i);
      pvalVec(i) = pval_i;
      zScoreVec(i) = zScore_i;
      if (i == 0) hwepval = hwepval_i;
    }

    return pvalVec;
  }
};

#endif
