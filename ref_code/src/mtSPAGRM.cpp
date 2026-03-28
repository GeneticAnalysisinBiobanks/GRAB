// SPAGRM.cpp -- nsSPAGRM free functions and mtSPAGRMClass method implementations

#include <RcppArmadillo.h>
#include <limits>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>

#include "mtSPAGRM.h"

namespace nsSPAGRM {

void mgf(
  double t,
  const arma::vec& resid_unrelated_outliers,
  const std::vector<std::array<double, 2>>& TwoSubj_resid,
  const std::vector<std::vector<double>>& TwoSubj_rho,
  const std::vector<UpdatedThreeSubj>& threeSubj,
  double MAF,
  MgfWorkspace& ws
) {
  arma::uword w = 0;
  const arma::uword nu = resid_unrelated_outliers.n_elem;

  // ---- Unrelated outliers: Armadillo SIMD lazy ops, pre-sized workspace ----
  // All assignments evaluate into pre-allocated buffers — zero heap allocation.
  if (nu > 0) {
    ws.ul_lambda  = arma::exp(t * resid_unrelated_outliers);
    ws.ul_alpha   = (1.0 - MAF) + MAF * ws.ul_lambda;
    ws.ul_alpha_1 = MAF * resid_unrelated_outliers % ws.ul_lambda;
    ws.ul_alpha_2 = resid_unrelated_outliers % ws.ul_alpha_1;
    ws.MGF0.subvec(0, nu - 1) = ws.ul_alpha % ws.ul_alpha;
    ws.MGF1.subvec(0, nu - 1) = 2.0 * ws.ul_alpha % ws.ul_alpha_1;
    ws.MGF2.subvec(0, nu - 1) = 2.0 * (ws.ul_alpha_1 % ws.ul_alpha_1
                                      + ws.ul_alpha   % ws.ul_alpha_2);
    w += nu;
  }

  // ---- Two-subject families: scalar inner loop, writes directly into ws.MGF* ----
  // All arithmetic is scalar — zero temporary arma::vec allocations.
  const int n1 = static_cast<int>(TwoSubj_resid.size());
  {
    double* g0 = ws.MGF0.memptr();
    double* g1 = ws.MGF1.memptr();
    double* g2 = ws.MGF2.memptr();
    const double maf01         = MAF * (1.0 - MAF);
    const double one_minus_MAF = 1.0 - MAF;

    for (int i = 0; i < n1; ++i) {
      const double* rho = TwoSubj_rho[i].data();
      const size_t ki = TwoSubj_rho[i].size();
      const double R1    = TwoSubj_resid[i][0];
      const double R2    = TwoSubj_resid[i][1];
      const double etR1  = std::exp(t * R1);
      const double etR2  = std::exp(t * R2);
      const double Rsum  = R1 + R2;
      const double R1sq  = R1   * R1;
      const double R2sq  = R2   * R2;
      const double Rssq  = Rsum * Rsum;
      const double etR12 = etR1 * etR2;

      for (size_t j = 0; j < ki; ++j) {
        const double tj  = (1.0 - rho[j]) * maf01;
        const double m1j = etR1  * tj;
        const double m2j = etR2  * tj;
        const double m3j = etR12 * (MAF - tj);
        const size_t idx = w + j;
        g0[idx] = m1j + m2j + m3j - tj + one_minus_MAF;
        g1[idx] = R1 * m1j + R2 * m2j + Rsum * m3j;
        g2[idx] = R1sq * m1j + R2sq * m2j + Rssq * m3j;
      }
      w += ki;
    }
  }

  // ---- Three-subject families: scalar reduction — zero temporary allocations ----
  const int n2 = static_cast<int>(threeSubj.size());
  for (int i = 0; i < n2; ++i) {
    const double* ss = threeSubj[i].stand_S.data();
    const double* ap = threeSubj[i].arr_prob.data();
    const size_t ns = threeSubj[i].stand_S.size();
    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (size_t j = 0; j < ns; ++j) {
      const double m0j = std::exp(t * ss[j]) * ap[j];
      const double m1j = ss[j] * m0j;
      s0 += m0j;
      s1 += m1j;
      s2 += ss[j] * m1j;
    }
    ws.MGF0[w] = s0;
    ws.MGF1[w] = s1;
    ws.MGF2[w] = s2;
    ++w;
  }
}

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
  int maxiter
) {
  double t = init_t;
  double CGF1 = 0.0, CGF2 = 0.0;
  double diff_t = std::numeric_limits<double>::infinity();

  const double mean = 2.0 * MAF * sum_R_nonOutlier;
  const double var  = 2.0 * MAF * (1.0 - MAF) * R_GRM_R_nonOutlier;

  // ws.MGF0/1/2 and ws.temp are pre-allocated — no heap inside the NR loop.
  for (int iter = 1; iter < maxiter; ++iter) {
    const double old_t      = t;
    const double old_diff_t = diff_t;
    const double old_CGF1   = CGF1;

    mgf(t, resid_unrelated_outliers, TwoSubj_resid, TwoSubj_rho, threeSubj, MAF, ws);

    ws.temp = ws.MGF1 / ws.MGF0;
    CGF1 = arma::accu(ws.temp) + mean + var * t - Score;
    CGF2 = arma::accu(ws.MGF2 / ws.MGF0) - arma::accu(ws.temp % ws.temp) + var;

    diff_t = -CGF1 / CGF2;

    if (std::isnan(diff_t) || std::isinf(CGF2)) {
      t = t / 2.0;
      diff_t = std::min(std::abs(t), 1.0) * arma::sign(Score);
      continue;
    }

    if (std::isnan(old_CGF1) || (arma::sign(old_CGF1) != 0 && arma::sign(CGF1) != arma::sign(old_CGF1))) {
      if (std::abs(diff_t) < tol) {
        t = old_t + diff_t;
        break;
      } else {
        while (std::abs(old_diff_t) > tol && std::abs(diff_t) > std::abs(old_diff_t) - tol)
          diff_t /= 2.0;
        t = old_t + diff_t;
        continue;
      }
    }

    if (arma::sign(Score) != arma::sign(old_t + diff_t) &&
        (arma::sign(old_CGF1) == 0 || arma::sign(CGF1) == arma::sign(old_CGF1))) {
      while (arma::sign(Score) != arma::sign(old_t + diff_t))
        diff_t /= 2.0;
      t = old_t + diff_t;
      continue;
    }

    t = old_t + diff_t;
    if (std::abs(diff_t) < tol) break;
  }

  // Final mgf evaluation at convergence so ws.MGF0/1/2 hold values at the returned t.
  // getProbSpa() uses these directly — eliminates a redundant mgf() call there.
  mgf(t, resid_unrelated_outliers, TwoSubj_resid, TwoSubj_rho, threeSubj, MAF, ws);
  return t;
}

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
) {
  // fastGetRoot fills ws.MGF0/1/2 at the returned zeta — no extra mgf() call needed.
  zeta = fastGetRoot(resid_unrelated_outliers, TwoSubj_resid, TwoSubj_rho, threeSubj,
                     sum_R_nonOutlier, R_GRM_R_nonOutlier, Score, MAF, zeta, tol, ws);

  const double mean = 2.0 * MAF * sum_R_nonOutlier;
  const double var  = 2.0 * MAF * (1.0 - MAF) * R_GRM_R_nonOutlier;

  // ws.temp may hold stale values from NR; recompute at zeta (no alloc, ws.temp is pre-sized).
  ws.temp = ws.MGF1 / ws.MGF0;
  const double CGF0     = arma::accu(arma::log(ws.MGF0)) + mean * zeta + 0.5 * var * zeta * zeta;
  const double CGF2_val = arma::accu(ws.MGF2 / ws.MGF0) - arma::accu(ws.temp % ws.temp) + var;

  const double w_val = arma::sign(zeta) * std::sqrt(2.0 * (zeta * Score - CGF0));
  const double v_val = zeta * std::sqrt(CGF2_val);
  const double u     = w_val + (1.0 / w_val) * std::log(v_val / w_val);

  boost::math::normal_distribution<double> ndist(0.0, 1.0);
  return lower_tail ? boost::math::cdf(ndist, u)
                    : boost::math::cdf(boost::math::complement(ndist, u));
}

} // namespace nsSPAGRM

mtSPAGRMClass::mtSPAGRMClass(
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
)
  : m_resid(std::move(resid)),
    m_resid_unrelated_outliers(std::move(fam.resid_unrelated_outliers)),
    m_sum_unrelated_outliers2(arma::dot(m_resid_unrelated_outliers, m_resid_unrelated_outliers)),
    m_sum_R_nonOutlier(sum_R_nonOutlier),
    m_R_GRM_R_nonOutlier(R_GRM_R_nonOutlier),
    m_R_GRM_R_TwoSubjOutlier(R_GRM_R_TwoSubjOutlier),
    m_R_GRM_R(R_GRM_R),
    m_MAF_interval(std::move(MAF_interval)),
    m_TwoSubj_resid_list(std::move(fam.twoSubj_resid)),
    m_TwoSubj_rho_list(std::move(fam.twoSubj_rho)),
    m_ThreeSubj_standS_list(std::move(fam.threeSubj_standS)),
    m_ThreeSubj_CLT_list(std::move(fam.threeSubj_CLT)),
    m_SPA_Cutoff(SPA_Cutoff),
    m_zeta(zeta),
    m_tol(tol)
{
  // Pre-allocate workspace (permanent size; reused across all marker calls)
  const arma::uword n_unrel = static_cast<arma::uword>(m_resid_unrelated_outliers.n_elem);
  const arma::uword mgfSz   = static_cast<arma::uword>(
      nsSPAGRM::mgfOutputSize(n_unrel, m_TwoSubj_rho_list, m_ThreeSubj_standS_list.size()));
  m_workspace = nsSPAGRM::MgfWorkspace(mgfSz, n_unrel);

  // Pre-allocate three-subject scratch:
  // stand_S copied once at construction (path 1 never changes it);
  // arr_prob sized for in-place update per marker.
  const int n3 = static_cast<int>(m_ThreeSubj_standS_list.size());
  m_threeSubj_scratch.resize(n3);
  for (int i = 0; i < n3; ++i) {
    m_threeSubj_scratch[i].stand_S = m_ThreeSubj_standS_list[i];
    m_threeSubj_scratch[i].arr_prob.resize(m_ThreeSubj_standS_list[i].size());
  }
}

double mtSPAGRMClass::getMarkerPval(
  const arma::vec& GVec,
  double altFreq,
  double& zScore
) {
  const double MAF       = std::min(altFreq, 1.0 - altFreq);
  const double Score     = arma::dot(GVec, m_resid) - arma::mean(GVec) * arma::accu(m_resid);
  const double G_var     = 2.0 * MAF * (1.0 - MAF);
  const double Score_var = G_var * m_R_GRM_R;
  zScore = Score / std::sqrt(Score_var);

  if (std::abs(zScore) <= m_SPA_Cutoff) {
    boost::math::normal_distribution<double> ndist(0.0, 1.0);
    return 2.0 * boost::math::cdf(boost::math::complement(ndist, std::abs(zScore)));
  }

  const int order2 = static_cast<int>(
      std::lower_bound(m_MAF_interval.begin(), m_MAF_interval.end(), MAF)
      - m_MAF_interval.begin());
  const int order1 = order2 - 1;
  const double MAF_ratio    = (m_MAF_interval[order2] - MAF)
                            / (m_MAF_interval[order2] - m_MAF_interval[order1]);
  const double one_minus_mr = 1.0 - MAF_ratio;

  // Update arr_prob in-place and compute Var_ThreeOutlier — zero allocation.
  // stand_S was pre-populated at construction and is constant across markers.
  double Var_ThreeOutlier = 0.0;
  const int n3 = static_cast<int>(m_threeSubj_scratch.size());
  for (int i = 0; i < n3; ++i) {
    const double* c1     = m_ThreeSubj_CLT_list[i].colptr(order1);
    const double* c2     = m_ThreeSubj_CLT_list[i].colptr(order2);
    const double* ss     = m_threeSubj_scratch[i].stand_S.data();
    double*       ap     = m_threeSubj_scratch[i].arr_prob.data();
    const size_t sz = m_threeSubj_scratch[i].stand_S.size();
    double s1 = 0.0, s2 = 0.0;
    for (size_t k = 0; k < sz; ++k) {
      ap[k] = MAF_ratio * c1[k] + one_minus_mr * c2[k];
      const double t = ss[k] * ap[k];
      s1 += t;
      s2 += ss[k] * t;
    }
    Var_ThreeOutlier += s2 - s1 * s1;
  }

  const double EmpVar   = G_var * (m_R_GRM_R_nonOutlier + m_sum_unrelated_outliers2
                                 + m_R_GRM_R_TwoSubjOutlier)
                        + Var_ThreeOutlier;
  const double Score_adj = Score * std::sqrt(EmpVar / Score_var);

  double zeta1 = std::abs(Score_adj) / Score_var; zeta1 = std::min(zeta1, 1.2);
  const double zeta2 = -std::abs(m_zeta);

  const double pval1 = nsSPAGRM::getProbSpa(
    m_resid_unrelated_outliers, m_TwoSubj_resid_list, m_TwoSubj_rho_list,
    m_threeSubj_scratch, m_sum_R_nonOutlier, m_R_GRM_R_nonOutlier,
    std::abs(Score_adj), MAF, false, zeta1, 1e-4, m_workspace);
  const double pval2 = nsSPAGRM::getProbSpa(
    m_resid_unrelated_outliers, m_TwoSubj_resid_list, m_TwoSubj_rho_list,
    m_threeSubj_scratch, m_sum_R_nonOutlier, m_R_GRM_R_nonOutlier,
    -std::abs(Score_adj), MAF, true, zeta2, m_tol, m_workspace);
  return pval1 + pval2;
}
