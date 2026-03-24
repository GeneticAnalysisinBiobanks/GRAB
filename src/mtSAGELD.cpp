// SAGELD.cpp -- mtSAGELDClass method implementations

#include <RcppArmadillo.h>
#include <limits>
#include <algorithm>
#include <boost/math/distributions/normal.hpp>

#include "mtSAGELD.h"
#include "mtSPAGRM.h"

mtSAGELDClass::mtSAGELDClass(
  std::string Method,
  arma::mat XTs,
  arma::mat SS,
  arma::mat AtS,
  arma::mat Q,
  arma::mat A21,
  arma::mat TTs,
  arma::mat Tys,
  arma::vec sol,
  arma::vec blups,
  double sig,
  arma::vec resid,
  arma::vec resid_G,
  arma::vec resid_GxE,
  arma::vec resid_E,
  arma::vec resid_unrelated_outliers,
  arma::vec resid_unrelated_outliers_G,
  arma::vec resid_unrelated_outliers_GxE,
  double sum_R_nonOutlier,
  double sum_R_nonOutlier_G,
  double sum_R_nonOutlier_GxE,
  double R_GRM_R,
  double R_GRM_R_G,
  double R_GRM_R_GxE,
  double R_GRM_R_G_GxE,
  double R_GRM_R_E,
  double R_GRM_R_nonOutlier,
  double R_GRM_R_nonOutlier_G,
  double R_GRM_R_nonOutlier_GxE,
  double R_GRM_R_nonOutlier_G_GxE,
  double R_GRM_R_TwoSubjOutlier,
  double R_GRM_R_TwoSubjOutlier_G,
  double R_GRM_R_TwoSubjOutlier_GxE,
  double R_GRM_R_TwoSubjOutlier_G_GxE,
  std::vector<TwoSubjFamily> TwoSubj_list,
  std::vector<ThreeSubjFamily> ThreeSubj_list,
  std::vector<double> MAF_interval,
  double zScoreE_cutoff,
  double SPA_Cutoff,
  double zeta,
  double tol
)
  : m_Method(std::move(Method)),
    m_XTs(std::move(XTs)),
    m_SS(std::move(SS)),
    m_AtS(std::move(AtS)),
    m_Q(std::move(Q)),
    m_A21(std::move(A21)),
    m_TTs(std::move(TTs)),
    m_Tys(std::move(Tys)),
    m_sol(std::move(sol)),
    m_blups(std::move(blups)),
    m_sig(sig),
    m_ncov(static_cast<double>(m_sol.n_elem)),
    m_resid(std::move(resid)),
    m_resid_G(std::move(resid_G)),
    m_resid_GxE(std::move(resid_GxE)),
    m_resid_E(std::move(resid_E)),
    m_resid_unrelated_outliers(std::move(resid_unrelated_outliers)),
    m_resid_unrelated_outliers_G(std::move(resid_unrelated_outliers_G)),
    m_resid_unrelated_outliers_GxE(std::move(resid_unrelated_outliers_GxE)),
    m_sum_unrelated_outliers2(arma::dot(m_resid_unrelated_outliers, m_resid_unrelated_outliers)),
    m_sum_unrelated_outliers_G2(arma::dot(m_resid_unrelated_outliers_G, m_resid_unrelated_outliers_G)),
    m_sum_unrelated_outliers_GxE2(arma::dot(m_resid_unrelated_outliers_GxE, m_resid_unrelated_outliers_GxE)),
    m_sum_unrelated_outliers_G_GxE2(2.0 * arma::dot(m_resid_unrelated_outliers_G, m_resid_unrelated_outliers_GxE)),
    m_sum_R_nonOutlier(sum_R_nonOutlier),
    m_sum_R_nonOutlier_G(sum_R_nonOutlier_G),
    m_sum_R_nonOutlier_GxE(sum_R_nonOutlier_GxE),
    m_R_GRM_R(R_GRM_R),
    m_R_GRM_R_G(R_GRM_R_G),
    m_R_GRM_R_GxE(R_GRM_R_GxE),
    m_R_GRM_R_G_GxE(R_GRM_R_G_GxE),
    m_R_GRM_R_E(R_GRM_R_E),
    m_R_GRM_R_nonOutlier(R_GRM_R_nonOutlier),
    m_R_GRM_R_nonOutlier_G(R_GRM_R_nonOutlier_G),
    m_R_GRM_R_nonOutlier_GxE(R_GRM_R_nonOutlier_GxE),
    m_R_GRM_R_nonOutlier_G_GxE(R_GRM_R_nonOutlier_G_GxE),
    m_R_GRM_R_TwoSubjOutlier(R_GRM_R_TwoSubjOutlier),
    m_R_GRM_R_TwoSubjOutlier_G(R_GRM_R_TwoSubjOutlier_G),
    m_R_GRM_R_TwoSubjOutlier_GxE(R_GRM_R_TwoSubjOutlier_GxE),
    m_R_GRM_R_TwoSubjOutlier_G_GxE(R_GRM_R_TwoSubjOutlier_G_GxE),
    m_TwoSubj_list(std::move(TwoSubj_list)),
    m_ThreeSubj_list(std::move(ThreeSubj_list)),
    m_MAF_interval(std::move(MAF_interval)),
    m_zScoreE_cutoff(zScoreE_cutoff),
    m_SPA_Cutoff(SPA_Cutoff),
    m_zeta(zeta),
    m_tol(tol)
{
  int nTwo = static_cast<int>(m_TwoSubj_list.size());
  m_TwoSubj_resid_list.resize(nTwo);
  m_TwoSubj_rho_list.resize(nTwo);
  for (int i = 0; i < nTwo; ++i) {
    m_TwoSubj_resid_list[i] = m_TwoSubj_list[i].Resid;
    m_TwoSubj_rho_list[i] = m_TwoSubj_list[i].Rho;
  }

  m_pvalVec.resize(2);
  m_zScoreVec.resize(2);
  m_BetaVec.resize(2);
  m_seBetaVec.resize(2);

  // Pre-allocate SPA workspace (constant size; reused across all marker calls)
  const arma::uword n_unrel = static_cast<arma::uword>(m_resid_unrelated_outliers.n_elem);
  const arma::uword mgfSz   = static_cast<arma::uword>(
      nsSPAGRM::mgfOutputSize(n_unrel, m_TwoSubj_rho_list, m_ThreeSubj_list.size()));
  m_workspace = nsSPAGRM::MgfWorkspace(mgfSz, n_unrel);

  // Pre-allocate three-subject scratch:
  // stand_S pre-copied from ThreeSubjFamily.stand_S for path-1 (invariant across markers).
  // arr_prob sized for in-place update each marker.
  const int n3 = static_cast<int>(m_ThreeSubj_list.size());
  m_threeSubj_scratch.resize(n3);
  for (int i = 0; i < n3; ++i) {
    m_threeSubj_scratch[i].stand_S = m_ThreeSubj_list[i].stand_S;
    m_threeSubj_scratch[i].arr_prob.resize(m_ThreeSubj_list[i].stand_S.size());
  }

  // Pre-allocate path-2 (lambda_i) scratch: avoids heap allocation per marker
  m_resid_outliers_i.set_size(n_unrel);
  m_twoResid_i.resize(nTwo);
}

double mtSAGELDClass::getMarkerPval(
  const arma::vec& GVec,
  double altFreq
) {
  const double MAF = std::min(altFreq, 1.0 - altFreq);

  arma::mat GVec2, m_H1, m_H2, m_AtH, m_R, m_Cfix, m_Cran, m_GtG, m_V;

  if (m_Method == "GALLOP") {
    GVec2 = arma::repmat(GVec.t(), 2, 1);
    GVec2 = GVec2.reshape(GVec.n_elem * 2, 1);
    m_H1 = GVec.t() * m_XTs; m_H1 = m_H1.reshape(m_ncov, 2);
    m_H2 = m_SS.each_col() % GVec2;
    m_AtH = GVec.t() * m_AtS; m_AtH = m_AtH.reshape(m_ncov, 2);
    m_R = m_H1 - m_AtH;
    m_Cfix = arma::inv(m_Q) * m_R;
    m_Cran = m_H2 - m_A21 * m_Cfix;
    m_GtG = (GVec % GVec).t() * m_TTs; m_GtG = m_GtG.reshape(2, 2);
    arma::mat m_Gty = GVec.t() * m_Tys; m_Gty = m_Gty.reshape(2, 1);
    arma::mat m_V2  = m_GtG - m_H1.t() * m_Cfix - m_H2.t() * m_Cran;
    arma::mat m_v   = m_Gty - m_H1.t() * m_sol   - m_H2.t() * m_blups;
    arma::mat m_intV = arma::inv(m_V2);
    arma::vec Theta = m_intV * m_v;
    arma::vec SD    = m_intV.diag();
    arma::vec SE    = m_sig * arma::sqrt(SD);
    arma::vec pval  = 2.0 * arma::normcdf(-1.0 * arma::abs(Theta / SE));
    m_pvalVec[0]   = pval[0];   m_pvalVec[1]   = pval[1];
    m_BetaVec[0]   = Theta[0];  m_BetaVec[1]   = Theta[1];
    m_seBetaVec[0] = SE[0];     m_seBetaVec[1] = SE[1];
    return 0.0;
  }

  // ---- Non-GALLOP SPA branch ----
  const double G_var    = 2.0 * MAF * (1.0 - MAF);
  const double zScore_G = arma::dot(GVec, m_resid_G) / std::sqrt(G_var * m_R_GRM_R_G);
  m_zScoreVec[0] = zScore_G;
  m_pvalVec[0]   = 2.0 * arma::normcdf(-std::abs(zScore_G));

  const double zScore_E = arma::dot(GVec, m_resid_E) / std::sqrt(G_var * m_R_GRM_R_E);

  // Precompute MAF interpolation indices (same for both paths)
  const int order2 = static_cast<int>(
      std::lower_bound(m_MAF_interval.begin(), m_MAF_interval.end(), MAF)
      - m_MAF_interval.begin());
  const int order1 = order2 - 1;
  const double MAF_ratio    = (m_MAF_interval[order2] - MAF)
                            / (m_MAF_interval[order2] - m_MAF_interval[order1]);
  const double one_minus_mr = 1.0 - MAF_ratio;
  const int n3   = static_cast<int>(m_ThreeSubj_list.size());
  const int nTwo = static_cast<int>(m_TwoSubj_list.size());

  double zScore_GxE, pval_GxE;

  if (std::abs(zScore_E) < m_zScoreE_cutoff) {
    // ---- Path 1: standard SPA with fixed residuals ----
    const double Score     = arma::dot(GVec, m_resid);
    const double Score_var = G_var * m_R_GRM_R;
    zScore_GxE = Score / std::sqrt(Score_var);

    if (std::abs(zScore_GxE) <= m_SPA_Cutoff) {
      pval_GxE = 2.0 * arma::normcdf(-std::abs(zScore_GxE));
    } else {
      // Update arr_prob in-place; stand_S pre-set at construction — zero allocation.
      double Var_ThreeOutlier = 0.0;
      for (int i = 0; i < n3; ++i) {
        const double* c1     = m_ThreeSubj_list[i].CLT.colptr(order1);
        const double* c2     = m_ThreeSubj_list[i].CLT.colptr(order2);
        const double* ss     = m_threeSubj_scratch[i].stand_S.data();
        double*       ap     = m_threeSubj_scratch[i].arr_prob.data();
        const size_t sz = m_threeSubj_scratch[i].stand_S.size();
        double s1 = 0.0, s2 = 0.0;
        for (size_t k = 0; k < sz; ++k) {
          ap[k] = MAF_ratio * c1[k] + one_minus_mr * c2[k];
          const double t = ss[k] * ap[k];
          s1 += t; s2 += ss[k] * t;
        }
        Var_ThreeOutlier += s2 - s1 * s1;
      }
      const double EmpVar   = G_var * (m_R_GRM_R_nonOutlier + m_sum_unrelated_outliers2
                                     + m_R_GRM_R_TwoSubjOutlier)
                            + Var_ThreeOutlier;
      const double Score_adj = Score * std::sqrt(EmpVar / Score_var);
      double zeta1 = std::abs(Score_adj) / Score_var; zeta1 = std::min(zeta1, 1.2);
      const double zeta2 = -std::abs(m_zeta);
      pval_GxE =
        nsSPAGRM::getProbSpa(
          m_resid_unrelated_outliers, m_TwoSubj_resid_list, m_TwoSubj_rho_list,
          m_threeSubj_scratch, m_sum_R_nonOutlier, m_R_GRM_R_nonOutlier,
          std::abs(Score_adj), MAF, false, zeta1, 1e-4, m_workspace)
      + nsSPAGRM::getProbSpa(
          m_resid_unrelated_outliers, m_TwoSubj_resid_list, m_TwoSubj_rho_list,
          m_threeSubj_scratch, m_sum_R_nonOutlier, m_R_GRM_R_nonOutlier,
          -std::abs(Score_adj), MAF, true, zeta2, m_tol, m_workspace);
    }
  } else {
    // ---- Path 2: lambda_i-adjusted SPA ----
    GVec2 = arma::repmat(GVec.t(), 2, 1);
    GVec2 = GVec2.reshape(GVec.n_elem * 2, 1);
    m_H1  = GVec.t() * m_XTs;  m_H1  = m_H1.reshape(m_ncov, 2);
    m_H2  = m_SS.each_col() % GVec2;
    m_AtH = GVec.t() * m_AtS;  m_AtH = m_AtH.reshape(m_ncov, 2);
    m_R    = m_H1 - m_AtH;
    m_Cfix = arma::inv(m_Q) * m_R;
    m_Cran = m_H2 - m_A21 * m_Cfix;
    m_GtG  = (GVec % GVec).t() * m_TTs; m_GtG = m_GtG.reshape(2, 2);
    m_V    = m_GtG - m_H1.t() * m_Cfix - m_H2.t() * m_Cran;
    const double lambda_i        = m_V(0, 1) / m_V(0, 0);
    const double R_GRM_R_i       = m_R_GRM_R_GxE  + lambda_i * lambda_i * m_R_GRM_R_G
                                 - lambda_i * m_R_GRM_R_G_GxE;
    const double Score     = arma::dot(GVec, m_resid_GxE)
                           - lambda_i * arma::dot(GVec, m_resid_G);
    const double Score_var = G_var * R_GRM_R_i;
    zScore_GxE = Score / std::sqrt(Score_var);

    if (std::abs(zScore_GxE) <= m_SPA_Cutoff) {
      pval_GxE = 2.0 * arma::normcdf(-std::abs(zScore_GxE));
    } else {
      // Update stand_S and arr_prob in-place, accumulate Var_ThreeOutlier — one pass, zero alloc.
      double Var_ThreeOutlier = 0.0;
      for (int i = 0; i < n3; ++i) {
        const double* gxe_ss = m_ThreeSubj_list[i].stand_S_GxE.data();
        const double* g_ss   = m_ThreeSubj_list[i].stand_S_G.data();
        const double* c1     = m_ThreeSubj_list[i].CLT.colptr(order1);
        const double* c2     = m_ThreeSubj_list[i].CLT.colptr(order2);
        double*       ss     = m_threeSubj_scratch[i].stand_S.data();
        double*       ap     = m_threeSubj_scratch[i].arr_prob.data();
        const size_t sz = m_threeSubj_scratch[i].stand_S.size();
        double s1 = 0.0, s2 = 0.0;
        for (size_t k = 0; k < sz; ++k) {
          ss[k] = gxe_ss[k] - lambda_i * g_ss[k];
          ap[k] = MAF_ratio * c1[k] + one_minus_mr * c2[k];
          const double t = ss[k] * ap[k];
          s1 += t; s2 += ss[k] * t;
        }
        Var_ThreeOutlier += s2 - s1 * s1;
      }

      // Compute adjusted scalars
      const double R_GRM_R_nonOutlier_i =   m_R_GRM_R_nonOutlier_GxE
                                          + lambda_i * lambda_i * m_R_GRM_R_nonOutlier_G
                                          - lambda_i * m_R_GRM_R_nonOutlier_G_GxE;
      const double sum_unrel2_i =   m_sum_unrelated_outliers_GxE2
                                  + lambda_i * lambda_i * m_sum_unrelated_outliers_G2
                                  - lambda_i * m_sum_unrelated_outliers_G_GxE2;
      const double R_GRM_R_Two_i =   m_R_GRM_R_TwoSubjOutlier_GxE
                                   + lambda_i * lambda_i * m_R_GRM_R_TwoSubjOutlier_G
                                   - lambda_i * m_R_GRM_R_TwoSubjOutlier_G_GxE;
      const double sum_R_nonOutlier_i = m_sum_R_nonOutlier_GxE
                                      - lambda_i * m_sum_R_nonOutlier_G;
      const double EmpVar   = G_var * (R_GRM_R_nonOutlier_i + sum_unrel2_i + R_GRM_R_Two_i)
                            + Var_ThreeOutlier;
      const double Score_adj = Score * std::sqrt(EmpVar / Score_var);
      double zeta1 = std::abs(Score_adj) / Score_var; zeta1 = std::min(zeta1, 1.2);
      const double zeta2 = -std::abs(m_zeta);

      // Compute adjusted unrelated-outlier residuals in-place scratch — zero alloc.
      {
        const double* gxe = m_resid_unrelated_outliers_GxE.memptr();
        const double* g   = m_resid_unrelated_outliers_G.memptr();
        double*       ri  = m_resid_outliers_i.memptr();
        const arma::uword nu = m_resid_outliers_i.n_elem;
        for (arma::uword k = 0; k < nu; ++k) ri[k] = gxe[k] - lambda_i * g[k];
      }
      // Compute adjusted two-subject residuals in-place scratch — zero alloc.
      for (int j = 0; j < nTwo; ++j) {
        const double* gxe = m_TwoSubj_list[j].Resid_GxE.data();
        const double* g   = m_TwoSubj_list[j].Resid_G.data();
        double*       ri  = m_twoResid_i[j].data();
        for (int k = 0; k < 2; ++k) ri[k] = gxe[k] - lambda_i * g[k];
      }

      pval_GxE =
        nsSPAGRM::getProbSpa(
          m_resid_outliers_i, m_twoResid_i, m_TwoSubj_rho_list,
          m_threeSubj_scratch, sum_R_nonOutlier_i, R_GRM_R_nonOutlier_i,
          std::abs(Score_adj), MAF, false, zeta1, 1e-4, m_workspace)
      + nsSPAGRM::getProbSpa(
          m_resid_outliers_i, m_twoResid_i, m_TwoSubj_rho_list,
          m_threeSubj_scratch, sum_R_nonOutlier_i, R_GRM_R_nonOutlier_i,
          -std::abs(Score_adj), MAF, true, zeta2, m_tol, m_workspace);
    }
  }

  m_pvalVec[1]   = pval_GxE;
  m_zScoreVec[1] = zScore_GxE;
  return 0.0;
}

