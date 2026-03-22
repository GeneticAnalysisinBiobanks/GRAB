#ifndef SPAGRM_H
#define SPAGRM_H

// SPAGRM.h -- Saddlepoint approximation with genetic relationship matrix (GRM)

#include <RcppArmadillo.h>
#include <limits>
#include <boost/math/distributions/normal.hpp>

namespace nsSPAGRM {

// Per-marker updated data for three-or-more-subject families (shared by SPAGRM and SAGELD)
struct UpdatedThreeSubj {
  arma::vec stand_S;
  arma::vec arr_prob;
};

// Family data passed from R at null-model construction time
struct FamilyData {
  arma::vec resid_unrelated_outliers;
  std::vector<arma::vec> twoSubj_resid;
  std::vector<arma::vec> twoSubj_rho;
  std::vector<arma::vec> threeSubj_standS;
  std::vector<arma::mat> threeSubj_CLT;
};

// MGF and its first two derivatives for the SPA approximation
arma::mat mgf(
  double t,
  const arma::vec& resid_unrelated_outliers,
  const std::vector<arma::vec>& TwoSubj_resid,
  const std::vector<arma::vec>& TwoSubj_rho,
  const std::vector<UpdatedThreeSubj>& threeSubj,
  double MAF
);

// Newton-Raphson root finder for the CGF equation
double fastGetRoot(
  const arma::vec& resid_unrelated_outliers,
  const std::vector<arma::vec>& TwoSubj_resid,
  const std::vector<arma::vec>& TwoSubj_rho,
  const std::vector<UpdatedThreeSubj>& threeSubj,
  double sum_R_nonOutlier,
  double R_GRM_R_nonOutlier,
  double Score, double MAF,
  double init_t, double tol,
  int maxiter = 50
);

// Saddlepoint approximation tail probability
double getProbSpa(
  const arma::vec& resid_unrelated_outliers,
  const std::vector<arma::vec>& TwoSubj_resid,
  const std::vector<arma::vec>& TwoSubj_rho,
  const std::vector<UpdatedThreeSubj>& threeSubj,
  double sum_R_nonOutlier,
  double R_GRM_R_nonOutlier,
  double Score, double MAF,
  bool lower_tail, double zeta, double tol
);

} // namespace nsSPAGRM

class mtSPAGRMClass {

private:

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

  mtSPAGRMClass(
    arma::vec resid,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double R_GRM_R_TwoSubjOutlier,
    double R_GRM_R,
    arma::vec MAF_interval,
    nsSPAGRM::FamilyData fam,
    double SPA_Cutoff,
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

  // Returns [zScore, pval, hwepval]
  std::vector<double> getResultVec(arma::vec GVec, double altFreq) {
    double zScore, hwepval;
    double pval = getMarkerPval(std::move(GVec), altFreq, zScore, hwepval);
    return {zScore, pval, hwepval};
  }

  static int resultSize() { return 3; }

  std::string getHeaderColumns() const {
    return "\tzScore\tPvalue\thwepval";
  }

};

#endif
