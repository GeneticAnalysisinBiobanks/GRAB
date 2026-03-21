#ifndef SPAGRM_H
#define SPAGRM_H

// SPAGRM.h -- Saddlepoint approximation with genetic relationship matrix (GRM)

#include <RcppArmadillo.h>
#include <limits>
#include <boost/math/distributions/normal.hpp>


namespace SPAGRM{

class SPAGRMClass {
  
private:

  // ---- Members ----
  const arma::vec m_resid;
  const arma::vec m_resid_unrelated_outliers;
  const double m_sum_unrelated_outliers2;
  const double m_sum_R_nonOutlier;
  const double m_R_GRM_R_nonOutlier;
  const double m_R_GRM_R_TwoSubjOutlier;
  const double m_R_GRM_R;
  const arma::vec m_MAF_interval;
  const std::vector<arma::vec> m_TwoSubj_resid_list;
  const std::vector<arma::vec> m_TwoSubj_rho_list;
  const std::vector<arma::vec> m_ThreeSubj_standS_list;
  const std::vector<arma::mat> m_ThreeSubj_CLT_list;

  const double m_SPA_Cutoff;
  const double m_zeta;
  const double m_tol;

public:

  // ---- Public interface ----
  SPAGRMClass(
    arma::vec resid,
    arma::vec resid_unrelated_outliers,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double R_GRM_R_TwoSubjOutlier,
    double R_GRM_R,
    arma::vec MAF_interval,
    std::vector<arma::vec> TwoSubj_resid,
    std::vector<arma::vec> TwoSubj_rho,
    std::vector<arma::vec> ThreeSubj_standS,
    std::vector<arma::mat> ThreeSubj_CLT,
    double SPA_Cutoff,
    double zeta,
    double tol
  );

  arma::mat mgf(
    double t,
    const std::vector<arma::vec>& arr_prob_list,
    double MAF
  );

  double fastGetRoot(
    const std::vector<arma::vec>& arr_prob_list,
    double Score,
    double MAF,
    double init_t,
    double tol,
    int maxiter = 50
  );

  double getProbSpa(
    const std::vector<arma::vec>& arr_prob_list,
    double Score,
    double MAF,
    bool lower_tail,
    double zeta,
    double tol
  );

  double getMarkerPval(
    arma::vec GVec,
    double altFreq,
    double& zScore,
    double& hwepval,
    double hwepvalCutoff = 0.1
  );

};

}

#endif
