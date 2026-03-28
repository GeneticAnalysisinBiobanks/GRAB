#ifndef SPAGRM_H
#define SPAGRM_H

// SPAGRM.h -- Saddlepoint approximation with genetic relationship matrix (GRM)

#include <RcppArmadillo.h>
#include <limits>
#include <array>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>

namespace nsSPAGRM {

// Per-marker updated data for three-or-more-subject families (shared by SPAGRM and SAGELD)
struct UpdatedThreeSubj {
  std::vector<double> stand_S;
  std::vector<double> arr_prob;
};

// Family data passed from R at null-model construction time
struct FamilyData {
  arma::vec resid_unrelated_outliers;
  std::vector<std::array<double, 2>> twoSubj_resid;
  std::vector<std::vector<double>> twoSubj_rho;
  std::vector<std::vector<double>> threeSubj_standS;
  std::vector<arma::mat> threeSubj_CLT;
};

// Total number of elements in the mgf output vectors (constant for a given dataset).
inline size_t mgfOutputSize(
  size_t n_unrelated,
  const std::vector<std::vector<double>>& TwoSubj_rho,
  size_t nThreeSubj
) {
  size_t sz = n_unrelated;
  for (const auto& r : TwoSubj_rho) sz += r.size();
  sz += nThreeSubj;
  return sz;
}

// Pre-allocated scratch workspace for all mgf/NR computations.
// Constructed once per mtSPAGRMClass instance; reused across all marker calls.
// Safe because each worker thread owns its own mtSPAGRMClass copy.
struct MgfWorkspace {
  // mgf output vectors (size = mgfOutputSize)
  arma::vec MGF0, MGF1, MGF2;
  // NR working vector (same size)
  arma::vec temp;
  // Unrelated-outliers intermediate vectors (size = n_unrelated_outliers)
  arma::vec ul_lambda, ul_alpha, ul_alpha_1, ul_alpha_2;

  MgfWorkspace() = default;
  MgfWorkspace(arma::uword mgfSz, arma::uword n_unrelated)
    : MGF0(mgfSz), MGF1(mgfSz), MGF2(mgfSz), temp(mgfSz),
      ul_lambda(n_unrelated), ul_alpha(n_unrelated),
      ul_alpha_1(n_unrelated), ul_alpha_2(n_unrelated)
  {}
};

// MGF and its first two derivatives written into ws.MGF0/MGF1/MGF2.
// Uses ws.ul_* scratch for the unrelated-outliers block.
// Zero heap allocations inside this function.
void mgf(
  double t,
  const arma::vec& resid_unrelated_outliers,
  const std::vector<std::array<double, 2>>& TwoSubj_resid,
  const std::vector<std::vector<double>>& TwoSubj_rho,
  const std::vector<UpdatedThreeSubj>& threeSubj,
  double MAF,
  MgfWorkspace& ws
);

// Newton-Raphson root finder for the CGF equation.
// On return, ws.MGF0/MGF1/MGF2 hold the mgf values evaluated at the returned t.
double fastGetRoot(
  const arma::vec& resid_unrelated_outliers,
  const std::vector<std::array<double, 2>>& TwoSubj_resid,
  const std::vector<std::vector<double>>& TwoSubj_rho,
  const std::vector<UpdatedThreeSubj>& threeSubj,
  double sum_R_nonOutlier,
  double R_GRM_R_nonOutlier,
  double Score, double MAF,
  double init_t, double tol,
  MgfWorkspace& ws,
  int maxiter = 50
);

// Saddlepoint approximation tail probability (fast path: shared workspace, no extra mgf() call).
double getProbSpa(
  const arma::vec& resid_unrelated_outliers,
  const std::vector<std::array<double, 2>>& TwoSubj_resid,
  const std::vector<std::vector<double>>& TwoSubj_rho,
  const std::vector<UpdatedThreeSubj>& threeSubj,
  double sum_R_nonOutlier,
  double R_GRM_R_nonOutlier,
  double Score, double MAF,
  bool lower_tail, double zeta, double tol,
  MgfWorkspace& ws
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
  const std::vector<double> m_MAF_interval;
  const std::vector<std::array<double, 2>> m_TwoSubj_resid_list;
  const std::vector<std::vector<double>> m_TwoSubj_rho_list;
  const std::vector<std::vector<double>> m_ThreeSubj_standS_list;
  const std::vector<arma::mat> m_ThreeSubj_CLT_list;

  const double m_SPA_Cutoff;
  const double m_zeta;
  const double m_tol;

  // Per-instance scratch — pre-allocated at construction, reused across all marker calls.
  // Thread-safe: each worker thread owns its own mtSPAGRMClass copy (via makeThreadContext).
  nsSPAGRM::MgfWorkspace m_workspace;
  // Three-subject scratch: stand_S copied once at construction; arr_prob updated per marker.
  std::vector<nsSPAGRM::UpdatedThreeSubj> m_threeSubj_scratch;

public:

  mtSPAGRMClass(
    arma::vec resid,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double R_GRM_R_TwoSubjOutlier,
    double R_GRM_R,
    std::vector<double> MAF_interval,
    nsSPAGRM::FamilyData fam,
    double SPA_Cutoff,
    double zeta,
    double tol
  );

  // GVec is read-only; no copy is made.
  double getMarkerPval(
    const arma::vec& GVec,
    double altFreq,
    double& zScore
  );

  // Fills rv with [zScore, pval]
  void getResultVec(const arma::vec& GVec, double altFreq, std::vector<double>& rv) {
    double zScore;
    double pval = getMarkerPval(GVec, altFreq, zScore);
    rv.clear();
    rv.push_back(zScore);
    rv.push_back(pval);
  }

  static int resultSize() { return 2; }

  std::string getHeaderColumns() const {
    return "\tzScore\tPvalue";
  }

};

#endif
