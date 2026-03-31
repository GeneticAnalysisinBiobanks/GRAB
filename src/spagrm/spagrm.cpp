// spagrm.cpp — nsSPAGRM free functions and SPAGRMClass (pure C++17 / Eigen / Boost)

#include "spagrm/spagrm.hpp"
#include "util/math_helper.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace nsSPAGRM {

// ──────────────────────────────────────────────────────────────────────
// MGF and its first two derivatives
// ──────────────────────────────────────────────────────────────────────

void mgf(
    double t,
    const Eigen::VectorXd& resid_unrelated_outliers,
    const std::vector<std::array<double, 2>>& TwoSubj_resid,
    const std::vector<std::vector<double>>& TwoSubj_rho,
    const std::vector<UpdatedThreeSubj>& threeSubj,
    double MAF,
    MgfWorkspace& ws)
{
  Eigen::Index w = 0;
  const Eigen::Index nu = resid_unrelated_outliers.size();

  // ── Unrelated outliers ─────────────────────────────────────────────
  if (nu > 0) {
    ws.ul_lambda  = (t * resid_unrelated_outliers).array().exp().matrix();
    ws.ul_alpha   = (Eigen::VectorXd::Constant(nu, 1.0 - MAF) +
                     MAF * ws.ul_lambda);
    ws.ul_alpha_1 = MAF * resid_unrelated_outliers.cwiseProduct(ws.ul_lambda);
    ws.ul_alpha_2 = resid_unrelated_outliers.cwiseProduct(ws.ul_alpha_1);

    ws.MGF0.segment(0, nu) = ws.ul_alpha.cwiseProduct(ws.ul_alpha);
    ws.MGF1.segment(0, nu) = 2.0 * ws.ul_alpha.cwiseProduct(ws.ul_alpha_1);
    ws.MGF2.segment(0, nu) = 2.0 * (ws.ul_alpha_1.cwiseProduct(ws.ul_alpha_1) +
                                     ws.ul_alpha.cwiseProduct(ws.ul_alpha_2));
    w += nu;
  }

  // ── Two-subject families ───────────────────────────────────────────
  {
    double* g0 = ws.MGF0.data();
    double* g1 = ws.MGF1.data();
    double* g2 = ws.MGF2.data();
    const double maf01         = MAF * (1.0 - MAF);
    const double one_minus_MAF = 1.0 - MAF;

    const int n1 = static_cast<int>(TwoSubj_resid.size());
    for (int i = 0; i < n1; ++i) {
      const double* rho = TwoSubj_rho[i].data();
      const size_t ki   = TwoSubj_rho[i].size();
      const double R1    = TwoSubj_resid[i][0];
      const double R2    = TwoSubj_resid[i][1];
      const double etR1  = std::exp(t * R1);
      const double etR2  = std::exp(t * R2);
      const double Rsum  = R1 + R2;
      const double R1sq  = R1 * R1;
      const double R2sq  = R2 * R2;
      const double Rssq  = Rsum * Rsum;
      const double etR12 = etR1 * etR2;

      for (size_t j = 0; j < ki; ++j) {
        const double tj  = (1.0 - rho[j]) * maf01;
        const double m1j = etR1  * tj;
        const double m2j = etR2  * tj;
        const double m3j = etR12 * (MAF - tj);
        const size_t idx = static_cast<size_t>(w) + j;
        g0[idx] = m1j + m2j + m3j - tj + one_minus_MAF;
        g1[idx] = R1 * m1j + R2 * m2j + Rsum * m3j;
        g2[idx] = R1sq * m1j + R2sq * m2j + Rssq * m3j;
      }
      w += static_cast<Eigen::Index>(ki);
    }
  }

  // ── Three-subject families ─────────────────────────────────────────
  const int n2 = static_cast<int>(threeSubj.size());
  for (int i = 0; i < n2; ++i) {
    const double* ss = threeSubj[i].stand_S.data();
    const double* ap = threeSubj[i].arr_prob.data();
    const size_t ns  = threeSubj[i].stand_S.size();
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


// ──────────────────────────────────────────────────────────────────────
// Newton-Raphson root finder
// ──────────────────────────────────────────────────────────────────────

static inline double signOf(double x) {
  return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);
}

double fastGetRoot(
    const Eigen::VectorXd& resid_unrelated_outliers,
    const std::vector<std::array<double, 2>>& TwoSubj_resid,
    const std::vector<std::vector<double>>& TwoSubj_rho,
    const std::vector<UpdatedThreeSubj>& threeSubj,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score, double MAF,
    double init_t, double tol,
    MgfWorkspace& ws,
    int maxiter)
{
  double t = init_t;
  double CGF1 = 0.0, CGF2 = 0.0;
  double diff_t = std::numeric_limits<double>::infinity();

  const double mean = 2.0 * MAF * sum_R_nonOutlier;
  const double var  = 2.0 * MAF * (1.0 - MAF) * R_GRM_R_nonOutlier;

  for (int iter = 1; iter < maxiter; ++iter) {
    const double old_t      = t;
    const double old_diff_t = diff_t;
    const double old_CGF1   = CGF1;

    mgf(t, resid_unrelated_outliers, TwoSubj_resid, TwoSubj_rho, threeSubj, MAF, ws);

    // temp = MGF1 ./ MGF0
    const Eigen::Index sz = ws.MGF0.size();
    for (Eigen::Index k = 0; k < sz; ++k)
      ws.temp[k] = ws.MGF1[k] / ws.MGF0[k];

    CGF1 = ws.temp.sum() + mean + var * t - Score;
    // CGF2 = sum(MGF2/MGF0) - sum(temp^2) + var
    double s1 = 0.0, s2 = 0.0;
    for (Eigen::Index k = 0; k < sz; ++k) {
      s1 += ws.MGF2[k] / ws.MGF0[k];
      s2 += ws.temp[k] * ws.temp[k];
    }
    CGF2 = s1 - s2 + var;

    diff_t = -CGF1 / CGF2;

    if (std::isnan(diff_t) || std::isinf(CGF2)) {
      t = t / 2.0;
      diff_t = std::min(std::abs(t), 1.0) * signOf(Score);
      continue;
    }

    if (std::isnan(old_CGF1) || (signOf(old_CGF1) != 0 && signOf(CGF1) != signOf(old_CGF1))) {
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

    if (signOf(Score) != signOf(old_t + diff_t) &&
        (signOf(old_CGF1) == 0 || signOf(CGF1) == signOf(old_CGF1))) {
      while (signOf(Score) != signOf(old_t + diff_t))
        diff_t /= 2.0;
      t = old_t + diff_t;
      continue;
    }

    t = old_t + diff_t;
    if (std::abs(diff_t) < tol) break;
  }

  // Final mgf evaluation at convergence.
  mgf(t, resid_unrelated_outliers, TwoSubj_resid, TwoSubj_rho, threeSubj, MAF, ws);
  return t;
}


// ──────────────────────────────────────────────────────────────────────
// Saddlepoint approximation tail probability
// ──────────────────────────────────────────────────────────────────────

double getProbSpa(
    const Eigen::VectorXd& resid_unrelated_outliers,
    const std::vector<std::array<double, 2>>& TwoSubj_resid,
    const std::vector<std::vector<double>>& TwoSubj_rho,
    const std::vector<UpdatedThreeSubj>& threeSubj,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score, double MAF,
    bool lower_tail, double zeta, double tol,
    MgfWorkspace& ws)
{
  zeta = fastGetRoot(resid_unrelated_outliers, TwoSubj_resid, TwoSubj_rho, threeSubj,
                     sum_R_nonOutlier, R_GRM_R_nonOutlier, Score, MAF, zeta, tol, ws);

  const double mean = 2.0 * MAF * sum_R_nonOutlier;
  const double var  = 2.0 * MAF * (1.0 - MAF) * R_GRM_R_nonOutlier;

  const Eigen::Index sz = ws.MGF0.size();
  for (Eigen::Index k = 0; k < sz; ++k)
    ws.temp[k] = ws.MGF1[k] / ws.MGF0[k];

  // CGF0 = sum(log(MGF0)) + mean*zeta + 0.5*var*zeta^2
  double logSum = 0.0;
  for (Eigen::Index k = 0; k < sz; ++k)
    logSum += std::log(ws.MGF0[k]);
  const double CGF0 = logSum + mean * zeta + 0.5 * var * zeta * zeta;

  // CGF2 = sum(MGF2/MGF0) - sum(temp^2) + var
  double s1 = 0.0, s2 = 0.0;
  for (Eigen::Index k = 0; k < sz; ++k) {
    s1 += ws.MGF2[k] / ws.MGF0[k];
    s2 += ws.temp[k] * ws.temp[k];
  }
  const double CGF2_val = s1 - s2 + var;

  const double w_val = signOf(zeta) * std::sqrt(2.0 * (zeta * Score - CGF0));
  const double v_val = zeta * std::sqrt(CGF2_val);
  const double u     = w_val + (1.0 / w_val) * std::log(v_val / w_val);

  return math::pnorm(u, 0.0, 1.0, lower_tail);
}

} // namespace nsSPAGRM


// ══════════════════════════════════════════════════════════════════════
// SPAGRMClass
// ══════════════════════════════════════════════════════════════════════

SPAGRMClass::SPAGRMClass(
    Eigen::VectorXd resid,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double R_GRM_R_TwoSubjOutlier,
    double R_GRM_R,
    std::vector<double> MAF_interval,
    nsSPAGRM::FamilyData fam,
    double SPA_Cutoff,
    double zeta,
    double tol)
  : m_resid(std::move(resid)),
    m_resid_unrelated_outliers(std::move(fam.resid_unrelated_outliers)),
    m_sum_unrelated_outliers2(m_resid_unrelated_outliers.squaredNorm()),
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
  const auto n_unrel = m_resid_unrelated_outliers.size();
  const auto mgfSz = static_cast<Eigen::Index>(
      nsSPAGRM::mgfOutputSize(
          static_cast<size_t>(n_unrel),
          m_TwoSubj_rho_list,
          m_ThreeSubj_standS_list.size()));
  m_workspace = nsSPAGRM::MgfWorkspace(mgfSz, n_unrel);

  const int n3 = static_cast<int>(m_ThreeSubj_standS_list.size());
  m_threeSubj_scratch.resize(n3);
  for (int i = 0; i < n3; ++i) {
    m_threeSubj_scratch[i].stand_S = m_ThreeSubj_standS_list[i];
    m_threeSubj_scratch[i].arr_prob.resize(m_ThreeSubj_standS_list[i].size());
  }
}

SPAGRMClass::SPAGRMClass(const SPAGRMClass& o)
  : m_resid(o.m_resid),
    m_resid_unrelated_outliers(o.m_resid_unrelated_outliers),
    m_sum_unrelated_outliers2(o.m_sum_unrelated_outliers2),
    m_sum_R_nonOutlier(o.m_sum_R_nonOutlier),
    m_R_GRM_R_nonOutlier(o.m_R_GRM_R_nonOutlier),
    m_R_GRM_R_TwoSubjOutlier(o.m_R_GRM_R_TwoSubjOutlier),
    m_R_GRM_R(o.m_R_GRM_R),
    m_MAF_interval(o.m_MAF_interval),
    m_TwoSubj_resid_list(o.m_TwoSubj_resid_list),
    m_TwoSubj_rho_list(o.m_TwoSubj_rho_list),
    m_ThreeSubj_standS_list(o.m_ThreeSubj_standS_list),
    m_ThreeSubj_CLT_list(o.m_ThreeSubj_CLT_list),
    m_SPA_Cutoff(o.m_SPA_Cutoff),
    m_zeta(o.m_zeta),
    m_tol(o.m_tol)
{
  // Rebuild scratch for thread safety
  const auto n_unrel = m_resid_unrelated_outliers.size();
  const auto mgfSz = static_cast<Eigen::Index>(
      nsSPAGRM::mgfOutputSize(
          static_cast<size_t>(n_unrel),
          m_TwoSubj_rho_list,
          m_ThreeSubj_standS_list.size()));
  m_workspace = nsSPAGRM::MgfWorkspace(mgfSz, n_unrel);

  const int n3 = static_cast<int>(m_ThreeSubj_standS_list.size());
  m_threeSubj_scratch.resize(n3);
  for (int i = 0; i < n3; ++i) {
    m_threeSubj_scratch[i].stand_S = m_ThreeSubj_standS_list[i];
    m_threeSubj_scratch[i].arr_prob.resize(m_ThreeSubj_standS_list[i].size());
  }
}


double SPAGRMClass::getMarkerPval(
    const Eigen::VectorXd& GVec,
    double altFreq,
    double& zScore)
{
  const double MAF       = std::min(altFreq, 1.0 - altFreq);
  const double Score     = GVec.dot(m_resid) - GVec.mean() * m_resid.sum();
  const double G_var     = 2.0 * MAF * (1.0 - MAF);
  const double Score_var = G_var * m_R_GRM_R;
  zScore = Score / std::sqrt(Score_var);

  if (std::abs(zScore) <= m_SPA_Cutoff) {
    return 2.0 * math::pnorm(std::abs(zScore), 0.0, 1.0, /*lower_tail=*/false);
  }

  const int order2 = static_cast<int>(
      std::lower_bound(m_MAF_interval.begin(), m_MAF_interval.end(), MAF)
      - m_MAF_interval.begin());
  const int order1 = order2 - 1;
  const double MAF_ratio    = (m_MAF_interval[order2] - MAF)
                            / (m_MAF_interval[order2] - m_MAF_interval[order1]);
  const double one_minus_mr = 1.0 - MAF_ratio;

  // Update arr_prob in-place for three-subject families.
  double Var_ThreeOutlier = 0.0;
  const int n3 = static_cast<int>(m_threeSubj_scratch.size());
  for (int i = 0; i < n3; ++i) {
    const double* c1 = m_ThreeSubj_CLT_list[i].col(order1).data();
    const double* c2 = m_ThreeSubj_CLT_list[i].col(order2).data();
    const double* ss = m_threeSubj_scratch[i].stand_S.data();
    double*       ap = m_threeSubj_scratch[i].arr_prob.data();
    const size_t  sz = m_threeSubj_scratch[i].stand_S.size();
    double s1 = 0.0, s2 = 0.0;
    for (size_t k = 0; k < sz; ++k) {
      ap[k] = MAF_ratio * c1[k] + one_minus_mr * c2[k];
      const double tmp = ss[k] * ap[k];
      s1 += tmp;
      s2 += ss[k] * tmp;
    }
    Var_ThreeOutlier += s2 - s1 * s1;
  }

  const double EmpVar = G_var * (m_R_GRM_R_nonOutlier + m_sum_unrelated_outliers2
                                 + m_R_GRM_R_TwoSubjOutlier)
                        + Var_ThreeOutlier;
  const double Score_adj = Score * std::sqrt(EmpVar / Score_var);

  double zeta1 = std::abs(Score_adj) / Score_var;
  zeta1 = std::min(zeta1, 1.2);
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
