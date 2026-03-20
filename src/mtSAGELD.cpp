// SAGELD.cpp -- SAGELDClass method implementations

#include <RcppArmadillo.h>
#include <limits>
#include <boost/math/distributions/normal.hpp>

#include "mtSAGELD.h"
#include "mtUTIL.h"


namespace SAGELD {

SAGELDClass::SAGELDClass(
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
) {
  m_Method = Method;

  m_XTs = XTs;
  m_SS = SS;
  m_AtS = AtS;
  m_Q = Q;
  m_A21 = A21;
  m_TTs = TTs;
  m_Tys = Tys;
  m_sol = sol;
  m_blups = blups;
  m_sig = sig;
  m_ncov = sol.n_elem;

  m_resid = resid;
  m_resid_G = resid_G;
  m_resid_GxE = resid_GxE;
  m_resid_E = resid_E;

  m_resid_unrelated_outliers = resid_unrelated_outliers;
  m_resid_unrelated_outliers_G = resid_unrelated_outliers_G;
  m_resid_unrelated_outliers_GxE = resid_unrelated_outliers_GxE;

  m_sum_unrelated_outliers2 = sum(resid_unrelated_outliers % resid_unrelated_outliers);
  m_sum_unrelated_outliers_G2 = sum(resid_unrelated_outliers_G % resid_unrelated_outliers_G);
  m_sum_unrelated_outliers_GxE2 = sum(resid_unrelated_outliers_GxE % resid_unrelated_outliers_GxE);
  m_sum_unrelated_outliers_G_GxE2 = 2 * sum(resid_unrelated_outliers_G % resid_unrelated_outliers_GxE);

  m_sum_R_nonOutlier = sum_R_nonOutlier;
  m_sum_R_nonOutlier_G = sum_R_nonOutlier_G;
  m_sum_R_nonOutlier_GxE = sum_R_nonOutlier_GxE;

  m_R_GRM_R = R_GRM_R;
  m_R_GRM_R_nonOutlier = R_GRM_R_nonOutlier;
  m_R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier;

  m_TwoSubj_list = std::move(TwoSubj_list);
  m_ThreeSubj_list = std::move(ThreeSubj_list);
  m_MAF_interval = MAF_interval;
  m_zScoreE_cutoff = zScoreE_cutoff;

  m_SPA_Cutoff = SPA_Cutoff;
  m_zeta = zeta;
  m_tol = tol;

  m_pvalVec.resize(2);
  m_zScoreVec.resize(2);
  m_BetaVec.resize(2);
  m_seBetaVec.resize(2);

}

arma::mat SAGELDClass::mgf(
  double t,
  const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
  double MAF
) {
  arma::vec lambda = arma::exp(t * m_resid_unrelated_outliers);
  arma::vec alpha = 1 - MAF + MAF * lambda;
  arma::vec alpha_1 = MAF * m_resid_unrelated_outliers % lambda;
  arma::vec alpha_2 = m_resid_unrelated_outliers % alpha_1;
  arma::vec M_G0_all = alpha % alpha;
  arma::vec M_G1_all = 2 * alpha % alpha_1;
  arma::vec M_G2_all = 2 * (alpha_1 % alpha_1 + alpha % alpha_2);

  int n1 = static_cast<int>(m_TwoSubj_list.size());
  if (n1 != 0) {
    for (int i = 0; i < n1; i++) {
      const auto& tsf = m_TwoSubj_list[i];
      const arma::vec& Resid = tsf.Resid;
      const arma::vec& Rho = tsf.Rho;
      arma::vec temp = (1 - Rho) * MAF * (1 - MAF);
      double R1 = Resid[0]; double etR1 = exp(t * R1);
      double R2 = Resid[1]; double etR2 = exp(t * R2);
      double Rsum = R1 + R2;
      arma::vec midterm1 = etR1 * temp;
      arma::vec midterm2 = etR2 * temp;
      arma::vec midterm3 = etR1 * etR2 * (MAF - temp);
      arma::vec M_G0 = midterm1 + midterm2 + midterm3 - temp + 1 - MAF;
      arma::vec M_G1 = R1 * midterm1 + R2 * midterm2 + Rsum * midterm3;
      arma::vec M_G2 = R1*R1 * midterm1 + R2*R2 * midterm2 + Rsum*Rsum * midterm3;
      M_G0_all = arma::join_cols(M_G0_all, M_G0);
      M_G1_all = arma::join_cols(M_G1_all, M_G1);
      M_G2_all = arma::join_cols(M_G2_all, M_G2);
    }
  }

  int n2 = static_cast<int>(update_ThreeSubj_list.size());
  if (n2 != 0) {
    for (int i = 0; i < n2; i++) {
      const auto& uts = update_ThreeSubj_list[i];
      const arma::vec& stand_S = uts.stand_S;
      const arma::vec& arr_prob = uts.arr_prob;
      arma::vec midterm0 = exp(t * stand_S) % arr_prob;
      arma::vec midterm1 = stand_S % midterm0;
      arma::vec midterm2 = stand_S % midterm1;
      M_G0_all = arma::join_cols(M_G0_all, arma::vec{arma::accu(midterm0)});
      M_G1_all = arma::join_cols(M_G1_all, arma::vec{arma::accu(midterm1)});
      M_G2_all = arma::join_cols(M_G2_all, arma::vec{arma::accu(midterm2)});
    }
  }

  return arma::join_rows(M_G0_all, M_G1_all, M_G2_all);
}

arma::mat SAGELDClass::mgf(
  double t,
  const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
  double MAF,
  double lambda_i
) {
  arma::vec m_resid_unrelated_outliers_i = m_resid_unrelated_outliers_GxE - lambda_i * m_resid_unrelated_outliers_G;
  arma::vec lambda = arma::exp(t * m_resid_unrelated_outliers_i);
  arma::vec alpha = 1 - MAF + MAF * lambda;
  arma::vec alpha_1 = MAF * m_resid_unrelated_outliers_i % lambda;
  arma::vec alpha_2 = m_resid_unrelated_outliers_i % alpha_1;
  arma::vec M_G0_all = alpha % alpha;
  arma::vec M_G1_all = 2 * alpha % alpha_1;
  arma::vec M_G2_all = 2 * (alpha_1 % alpha_1 + alpha % alpha_2);

  int n1 = static_cast<int>(m_TwoSubj_list.size());
  if (n1 != 0) {
    for (int i = 0; i < n1; i++) {
      const auto& tsf = m_TwoSubj_list[i];
      const arma::vec& Resid_G = tsf.Resid_G;
      const arma::vec& Resid_GxE = tsf.Resid_GxE;
      arma::vec Resid = Resid_GxE - lambda_i * Resid_G;
      const arma::vec& Rho = tsf.Rho;
      arma::vec temp = (1 - Rho) * MAF * (1 - MAF);
      double R1 = Resid[0]; double etR1 = exp(t * R1);
      double R2 = Resid[1]; double etR2 = exp(t * R2);
      double Rsum = R1 + R2;
      arma::vec midterm1 = etR1 * temp;
      arma::vec midterm2 = etR2 * temp;
      arma::vec midterm3 = etR1 * etR2 * (MAF - temp);
      arma::vec M_G0 = midterm1 + midterm2 + midterm3 - temp + 1 - MAF;
      arma::vec M_G1 = R1 * midterm1 + R2 * midterm2 + Rsum * midterm3;
      arma::vec M_G2 = R1*R1 * midterm1 + R2*R2 * midterm2 + Rsum*Rsum * midterm3;
      M_G0_all = arma::join_cols(M_G0_all, M_G0);
      M_G1_all = arma::join_cols(M_G1_all, M_G1);
      M_G2_all = arma::join_cols(M_G2_all, M_G2);
    }
  }

  int n2 = static_cast<int>(update_ThreeSubj_list.size());
  if (n2 != 0) {
    for (int i = 0; i < n2; i++) {
      const auto& uts = update_ThreeSubj_list[i];
      const arma::vec& stand_S = uts.stand_S;
      const arma::vec& arr_prob = uts.arr_prob;
      arma::vec midterm0 = exp(t * stand_S) % arr_prob;
      arma::vec midterm1 = stand_S % midterm0;
      arma::vec midterm2 = stand_S % midterm1;
      M_G0_all = arma::join_cols(M_G0_all, arma::vec{arma::accu(midterm0)});
      M_G1_all = arma::join_cols(M_G1_all, arma::vec{arma::accu(midterm1)});
      M_G2_all = arma::join_cols(M_G2_all, arma::vec{arma::accu(midterm2)});
    }
  }

  return arma::join_rows(M_G0_all, M_G1_all, M_G2_all);
}

double SAGELDClass::fastGetRoot(
  const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
  double Score,
  double MAF,
  double init_t,
  double tol,
  int maxiter
) {
  double t = init_t;
  arma::vec MGF0; arma::vec MGF1; arma::vec MGF2;
  double CGF1 = 0; double CGF2 = 0;
  double diff_t = std::numeric_limits<double>::infinity();

  double mean = 2 * MAF * m_sum_R_nonOutlier;
  double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier(0);

  for (int iter = 1; iter < maxiter; iter++) {
    double old_t = t;
    double old_diff_t = diff_t;
    double old_CGF1 = CGF1;

    arma::mat MGF_all = mgf(t, update_ThreeSubj_list, MAF);
    MGF0 = MGF_all.col(0);
    MGF1 = MGF_all.col(1);
    MGF2 = MGF_all.col(2);

    arma::vec temp = MGF1 / MGF0;
    CGF1 = arma::accu(temp) + mean + var * t - Score;
    CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
    diff_t = -CGF1 / CGF2;

    if (std::isnan(diff_t) || std::isinf(CGF2)) {
      t = t / 2;
      diff_t = std::min(std::abs(t), 1.0) * arma::sign(Score);
      continue;
    }
    if (std::isnan(old_CGF1) || (arma::sign(old_CGF1) != 0 && arma::sign(CGF1) != arma::sign(old_CGF1))) {
      if (std::abs(diff_t) < tol) {
        t = old_t + diff_t;
        break;
      } else {
        while (std::abs(old_diff_t) > tol && std::abs(diff_t) > std::abs(old_diff_t) - tol) {
          diff_t = diff_t / 2;
        }
        t = old_t + diff_t;
        continue;
      }
    }
    if (arma::sign(Score) != arma::sign(old_t + diff_t) &&
        (arma::sign(old_CGF1) == 0 || arma::sign(CGF1) == arma::sign(old_CGF1))) {
      while (arma::sign(Score) != arma::sign(old_t + diff_t)) {
        diff_t = diff_t / 2;
      }
      t = old_t + diff_t;
      continue;
    }
    t = old_t + diff_t;
    if (std::abs(diff_t) < tol) break;
  }
  return t;
}

double SAGELDClass::fastGetRoot(
  const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
  double m_R_GRM_R_nonOutlier_i,
  double lambda_i,
  double Score,
  double MAF,
  double init_t,
  double tol,
  int maxiter
) {
  double t = init_t;
  arma::vec MGF0; arma::vec MGF1; arma::vec MGF2;
  double CGF1 = 0; double CGF2 = 0;
  double diff_t = std::numeric_limits<double>::infinity();

  double mean = 2 * MAF * (m_sum_R_nonOutlier_GxE - lambda_i * m_sum_R_nonOutlier_G);
  double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier_i;

  for (int iter = 1; iter < maxiter; iter++) {
    double old_t = t;
    double old_diff_t = diff_t;
    double old_CGF1 = CGF1;

    arma::mat MGF_all = mgf(t, update_ThreeSubj_list, MAF, lambda_i);
    MGF0 = MGF_all.col(0);
    MGF1 = MGF_all.col(1);
    MGF2 = MGF_all.col(2);

    arma::vec temp = MGF1 / MGF0;
    CGF1 = arma::accu(temp) + mean + var * t - Score;
    CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
    diff_t = -CGF1 / CGF2;

    if (std::isnan(diff_t) || std::isinf(CGF2)) {
      t = t / 2;
      diff_t = std::min(std::abs(t), 1.0) * arma::sign(Score);
      continue;
    }
    if (std::isnan(old_CGF1) || (arma::sign(old_CGF1) != 0 && arma::sign(CGF1) != arma::sign(old_CGF1))) {
      if (std::abs(diff_t) < tol) {
        t = old_t + diff_t;
        break;
      } else {
        while (std::abs(old_diff_t) > tol && std::abs(diff_t) > std::abs(old_diff_t) - tol) {
          diff_t = diff_t / 2;
        }
        t = old_t + diff_t;
        continue;
      }
    }
    if (arma::sign(Score) != arma::sign(old_t + diff_t) &&
        (arma::sign(old_CGF1) == 0 || arma::sign(CGF1) == arma::sign(old_CGF1))) {
      while (arma::sign(Score) != arma::sign(old_t + diff_t)) {
        diff_t = diff_t / 2;
      }
      t = old_t + diff_t;
      continue;
    }
    t = old_t + diff_t;
    if (std::abs(diff_t) < tol) break;
  }
  return t;
}

double SAGELDClass::getProbSpa(
  const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
  double Score,
  double MAF,
  bool lower_tail,
  double zeta,
  double tol
) {
  zeta = fastGetRoot(update_ThreeSubj_list, Score, MAF, zeta, tol);
  arma::mat MGF_all = mgf(zeta, update_ThreeSubj_list, MAF);
  arma::vec MGF0 = MGF_all.col(0);
  arma::vec MGF1 = MGF_all.col(1);
  arma::vec MGF2 = MGF_all.col(2);

  double mean = 2 * MAF * m_sum_R_nonOutlier;
  double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier(0);

  arma::vec temp = MGF1 / MGF0;
  double CGF0 = arma::accu(log(MGF0)) + mean * zeta + 0.5 * var * zeta * zeta;
  double CGF2_val = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;

  double w = arma::sign(zeta) * sqrt(2 * (zeta * Score - CGF0));
  double v = zeta * sqrt(CGF2_val);
  double u = w + 1/w * log(v/w);

  boost::math::normal_distribution<double> ndist(0.0, 1.0);
  double pval = lower_tail ? boost::math::cdf(ndist, u)
                           : boost::math::cdf(boost::math::complement(ndist, u));
  return pval;
}

double SAGELDClass::getProbSpa(
  const std::vector<UpdatedThreeSubj>& update_ThreeSubj_list,
  double m_R_GRM_R_nonOutlier_i,
  double lambda_i,
  double Score,
  double MAF,
  bool lower_tail,
  double zeta,
  double tol
) {
  zeta = fastGetRoot(update_ThreeSubj_list, m_R_GRM_R_nonOutlier_i, lambda_i, Score, MAF, zeta, tol);
  arma::mat MGF_all = mgf(zeta, update_ThreeSubj_list, MAF, lambda_i);
  arma::vec MGF0 = MGF_all.col(0);
  arma::vec MGF1 = MGF_all.col(1);
  arma::vec MGF2 = MGF_all.col(2);

  double mean = 2 * MAF * (m_sum_R_nonOutlier_GxE - lambda_i * m_sum_R_nonOutlier_G);
  double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier_i;

  arma::vec temp = MGF1 / MGF0;
  double CGF0 = arma::accu(log(MGF0)) + mean * zeta + 0.5 * var * zeta * zeta;
  double CGF2_val = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;

  double w = arma::sign(zeta) * sqrt(2 * (zeta * Score - CGF0));
  double v = zeta * sqrt(CGF2_val);
  double u = w + 1/w * log(v/w);

  boost::math::normal_distribution<double> ndist(0.0, 1.0);
  double pval = lower_tail ? boost::math::cdf(ndist, u)
                           : boost::math::cdf(boost::math::complement(ndist, u));
  return pval;
}

double SAGELDClass::getMarkerPval(
  arma::vec GVec,
  double altFreq,
  double& hwepval,
  double hwepvalCutoff
) {
  gethwepval(GVec, hwepval, hwepvalCutoff);
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
        std::vector<UpdatedThreeSubj> update_ThreeSubj_list(n1);
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
        double pval1 = getProbSpa(update_ThreeSubj_list, std::abs(Score_adj), MAF, false, zeta1, 1e-4);
        double pval2 = getProbSpa(update_ThreeSubj_list, -std::abs(Score_adj), MAF, true, zeta2, m_tol);
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
        std::vector<UpdatedThreeSubj> update_ThreeSubj_list(n1);
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
        double pval1 = getProbSpa(update_ThreeSubj_list, m_R_GRM_R_nonOutlier_i, lambda_i, std::abs(Score_adj), MAF, false, zeta1, 1e-4);
        double pval2 = getProbSpa(update_ThreeSubj_list, m_R_GRM_R_nonOutlier_i, lambda_i, -std::abs(Score_adj), MAF, true, zeta2, m_tol);
        pval_GxE = pval1 + pval2;
      }
    }
    m_pvalVec.at(0) = pval_G; m_pvalVec.at(1) = pval_GxE;
    m_zScoreVec.at(0) = zScore_G; m_zScoreVec.at(1) = zScore_GxE;
  }
  return 0;
}

}
