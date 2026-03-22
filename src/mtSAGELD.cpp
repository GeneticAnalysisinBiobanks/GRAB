// SAGELD.cpp -- mtSAGELDClass method implementations

#include <RcppArmadillo.h>
#include <limits>
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
  arma::vec R_GRM_R,
  arma::vec R_GRM_R_nonOutlier,
  arma::vec R_GRM_R_TwoSubjOutlier,
  std::vector<TwoSubjFamily> TwoSubj_list,
  std::vector<ThreeSubjFamily> ThreeSubj_list,
  arma::vec MAF_interval,
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
    m_R_GRM_R(std::move(R_GRM_R)),
    m_R_GRM_R_nonOutlier(std::move(R_GRM_R_nonOutlier)),
    m_R_GRM_R_TwoSubjOutlier(std::move(R_GRM_R_TwoSubjOutlier)),
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
}

double mtSAGELDClass::getMarkerPval(
  arma::vec GVec,
  double altFreq
) {
  double MAF = std::min(altFreq, 1 - altFreq);

  arma::mat GVec2, GVecq2, m_H1, m_H2, m_AtH, m_R, m_Cfix, m_Cran, m_GtG, m_Gty, m_V, m_v, m_intV;

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
    m_Gty = GVec.t() * m_Tys; m_Gty = m_Gty.reshape(2, 1);
    m_V = m_GtG - m_H1.t() * m_Cfix - m_H2.t() * m_Cran;
    m_v = m_Gty - m_H1.t() * m_sol - m_H2.t() * m_blups;
    m_intV = arma::inv(m_V);
    arma::vec Theta = m_intV * m_v;
    arma::vec SD = m_intV.diag();
    arma::vec SE = m_sig * sqrt(SD);
    arma::vec pval = 2 * arma::normcdf(-1.0 * arma::abs(Theta / SE));
    m_pvalVec.at(0) = pval.at(0); m_pvalVec.at(1) = pval.at(1);
    m_BetaVec.at(0) = Theta.at(0); m_BetaVec.at(1) = Theta.at(1);
    m_seBetaVec.at(0) = SE.at(0); m_seBetaVec.at(1) = SE.at(1);
  } else {
    double zScore_G, zScore_GxE, pval_G, pval_GxE;
    double G_var = 2 * MAF * (1 - MAF);
    double Score_E = sum(GVec % m_resid_E);
    double zScore_E = Score_E / sqrt(G_var * m_R_GRM_R(4));
    zScore_G = sum(GVec % m_resid_G) / sqrt(G_var * m_R_GRM_R(1));
    pval_G = 2 * arma::normcdf(-1.0 * std::abs(zScore_G));

    if (std::abs(zScore_E) < m_zScoreE_cutoff) {
      double Score = sum(GVec % m_resid);
      double Score_var = G_var * m_R_GRM_R(0);
      zScore_GxE = Score / sqrt(Score_var);

      if (std::abs(zScore_GxE) <= m_SPA_Cutoff) {
        pval_GxE = 2 * arma::normcdf(-1.0 * std::abs(zScore_GxE));
      } else {
        int order2 = arma::index_max(m_MAF_interval >= MAF);
        int order1 = order2 - 1;
        double MAF_ratio = (m_MAF_interval[order2] - MAF) / (m_MAF_interval[order2] - m_MAF_interval[order1]);
        double Var_ThreeOutlier = 0;
        int n1 = static_cast<int>(m_ThreeSubj_list.size());
        std::vector<nsSPAGRM::UpdatedThreeSubj> update_ThreeSubj_list(n1);
        if (n1 != 0) {
          for (int i = 0; i < n1; i++) {
            const auto& tsf3 = m_ThreeSubj_list[i];
            const arma::mat& CLT_temp = tsf3.CLT;
            arma::vec stand_S = tsf3.stand_S;
            arma::vec CLT_temp1 = CLT_temp.col(order1);
            arma::vec CLT_temp2 = CLT_temp.col(order2);
            arma::vec arr_prob = MAF_ratio * CLT_temp1 + (1 - MAF_ratio) * CLT_temp2;
            update_ThreeSubj_list[i] = {stand_S, arr_prob};
            arma::vec temp1 = stand_S % arr_prob;
            double temp2 = arma::accu(temp1);
            double Var_ThreeOutlier_temp = arma::accu(stand_S % temp1) - temp2 * temp2;
            Var_ThreeOutlier = Var_ThreeOutlier + Var_ThreeOutlier_temp;
          }
        }
        double Var_nonOutlier = G_var * m_R_GRM_R_nonOutlier(0);
        double Var_unrelated_outliers = G_var * m_sum_unrelated_outliers2;
        double Var_TwoOutlier = G_var * m_R_GRM_R_TwoSubjOutlier(0);
        double EmpVar = Var_nonOutlier + Var_unrelated_outliers + Var_TwoOutlier + Var_ThreeOutlier;
        double Var_Ratio = Score_var / EmpVar;
        double Score_adj = Score / sqrt(Var_Ratio);
        double zeta1 = std::abs(Score_adj) / Score_var; zeta1 = std::min(zeta1, 1.2);
        double zeta2 = -std::abs(m_zeta);
        double pval1 = nsSPAGRM::getProbSpa(
          m_resid_unrelated_outliers, m_TwoSubj_resid_list, m_TwoSubj_rho_list,
          update_ThreeSubj_list, m_sum_R_nonOutlier, m_R_GRM_R_nonOutlier(0),
          std::abs(Score_adj), MAF, false, zeta1, 1e-4);
        double pval2 = nsSPAGRM::getProbSpa(
          m_resid_unrelated_outliers, m_TwoSubj_resid_list, m_TwoSubj_rho_list,
          update_ThreeSubj_list, m_sum_R_nonOutlier, m_R_GRM_R_nonOutlier(0),
          -std::abs(Score_adj), MAF, true, zeta2, m_tol);
        pval_GxE = pval1 + pval2;
      }
    } else {
      GVec2 = arma::repmat(GVec.t(), 2, 1);
      GVec2 = GVec2.reshape(GVec.n_elem * 2, 1);
      m_H1 = GVec.t() * m_XTs; m_H1 = m_H1.reshape(m_ncov, 2);
      m_H2 = m_SS.each_col() % GVec2;
      m_AtH = GVec.t() * m_AtS; m_AtH = m_AtH.reshape(m_ncov, 2);
      m_R = m_H1 - m_AtH;
      m_Cfix = arma::inv(m_Q) * m_R;
      m_Cran = m_H2 - m_A21 * m_Cfix;
      m_GtG = (GVec % GVec).t() * m_TTs; m_GtG = m_GtG.reshape(2, 2);
      m_V = m_GtG - m_H1.t() * m_Cfix - m_H2.t() * m_Cran;
      double lambda_i = m_V(0, 1) / m_V(0, 0);
      arma::vec m_resid_i = m_resid_GxE - lambda_i * m_resid_G;
      double m_R_GRM_R_i = m_R_GRM_R(2) + lambda_i * lambda_i * m_R_GRM_R(1) - lambda_i * m_R_GRM_R(3);
      double Score = sum(GVec % m_resid_i);
      double Score_var = G_var * m_R_GRM_R_i;
      zScore_GxE = Score / sqrt(Score_var);

      if (std::abs(zScore_GxE) <= m_SPA_Cutoff) {
        pval_GxE = 2 * arma::normcdf(-1.0 * std::abs(zScore_GxE));
      } else {
        int order2 = arma::index_max(m_MAF_interval >= MAF);
        int order1 = order2 - 1;
        double MAF_ratio = (m_MAF_interval[order2] - MAF) / (m_MAF_interval[order2] - m_MAF_interval[order1]);
        double Var_ThreeOutlier = 0;
        int n1 = static_cast<int>(m_ThreeSubj_list.size());
        std::vector<nsSPAGRM::UpdatedThreeSubj> update_ThreeSubj_list(n1);
        if (n1 != 0) {
          for (int i = 0; i < n1; i++) {
            const auto& tsf3 = m_ThreeSubj_list[i];
            const arma::mat& CLT_temp = tsf3.CLT;
            const arma::vec& stand_S_G = tsf3.stand_S_G;
            const arma::vec& stand_S_GxE = tsf3.stand_S_GxE;
            arma::vec stand_S = stand_S_GxE - lambda_i * stand_S_G;
            arma::vec CLT_temp1 = CLT_temp.col(order1);
            arma::vec CLT_temp2 = CLT_temp.col(order2);
            arma::vec arr_prob = MAF_ratio * CLT_temp1 + (1 - MAF_ratio) * CLT_temp2;
            update_ThreeSubj_list[i] = {stand_S, arr_prob};
            arma::vec temp1 = stand_S % arr_prob;
            double temp2 = arma::accu(temp1);
            double Var_ThreeOutlier_temp = arma::accu(stand_S % temp1) - temp2 * temp2;
            Var_ThreeOutlier = Var_ThreeOutlier + Var_ThreeOutlier_temp;
          }
        }
        double m_R_GRM_R_nonOutlier_i = m_R_GRM_R_nonOutlier(2) + lambda_i * lambda_i * m_R_GRM_R_nonOutlier(1) - lambda_i * m_R_GRM_R_nonOutlier(3);
        double Var_nonOutlier = G_var * m_R_GRM_R_nonOutlier_i;
        double m_sum_unrelated_outliers2_i = m_sum_unrelated_outliers_GxE2 + lambda_i * lambda_i * m_sum_unrelated_outliers_G2 - lambda_i * m_sum_unrelated_outliers_G_GxE2;
        double Var_unrelated_outliers = G_var * m_sum_unrelated_outliers2_i;
        double m_R_GRM_R_TwoSubjOutlier_i = m_R_GRM_R_TwoSubjOutlier(2) + lambda_i * lambda_i * m_R_GRM_R_TwoSubjOutlier(1) - lambda_i * m_R_GRM_R_TwoSubjOutlier(3);
        double Var_TwoOutlier = G_var * m_R_GRM_R_TwoSubjOutlier_i;
        double EmpVar = Var_nonOutlier + Var_unrelated_outliers + Var_TwoOutlier + Var_ThreeOutlier;
        double Var_Ratio = Score_var / EmpVar;
        double Score_adj = Score / sqrt(Var_Ratio);
        double zeta1 = std::abs(Score_adj) / Score_var; zeta1 = std::min(zeta1, 1.2);
        double zeta2 = -std::abs(m_zeta);

        // Construct modified data for lambda_i-adjusted SPA
        arma::vec resid_outliers_i = m_resid_unrelated_outliers_GxE - lambda_i * m_resid_unrelated_outliers_G;
        int nTwo = static_cast<int>(m_TwoSubj_list.size());
        std::vector<arma::vec> twoResid_i(nTwo);
        for (int j = 0; j < nTwo; ++j) {
          twoResid_i[j] = m_TwoSubj_list[j].Resid_GxE - lambda_i * m_TwoSubj_list[j].Resid_G;
        }
        double sum_R_nonOutlier_i = m_sum_R_nonOutlier_GxE - lambda_i * m_sum_R_nonOutlier_G;

        double pval1 = nsSPAGRM::getProbSpa(
          resid_outliers_i, twoResid_i, m_TwoSubj_rho_list,
          update_ThreeSubj_list, sum_R_nonOutlier_i, m_R_GRM_R_nonOutlier_i,
          std::abs(Score_adj), MAF, false, zeta1, 1e-4);
        double pval2 = nsSPAGRM::getProbSpa(
          resid_outliers_i, twoResid_i, m_TwoSubj_rho_list,
          update_ThreeSubj_list, sum_R_nonOutlier_i, m_R_GRM_R_nonOutlier_i,
          -std::abs(Score_adj), MAF, true, zeta2, m_tol);
        pval_GxE = pval1 + pval2;
      }
    }
    m_pvalVec.at(0) = pval_G; m_pvalVec.at(1) = pval_GxE;
    m_zScoreVec.at(0) = zScore_G; m_zScoreVec.at(1) = zScore_GxE;
  }
  return 0;
}

