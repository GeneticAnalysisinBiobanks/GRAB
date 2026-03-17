// SPAGRM.cpp -- SPAGRMClass method implementations

#include <RcppArmadillo.h>
#include <limits>
#include <boost/math/distributions/normal.hpp>

#include "SPAGRM.h"
#include "UTIL.h"

namespace SPAGRM {

SPAGRMClass::SPAGRMClass(arma::vec resid,
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
                         double tol) {
  m_resid = resid;
  m_resid_unrelated_outliers = resid_unrelated_outliers;
  m_sum_unrelated_outliers2 = sum(resid_unrelated_outliers % resid_unrelated_outliers);
  m_sum_R_nonOutlier = sum_R_nonOutlier;
  m_R_GRM_R_nonOutlier = R_GRM_R_nonOutlier;
  m_R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier;
  m_R_GRM_R = R_GRM_R;
  m_MAF_interval = MAF_interval;

  m_TwoSubj_resid_list = std::move(TwoSubj_resid);
  m_TwoSubj_rho_list   = std::move(TwoSubj_rho);
  m_ThreeSubj_standS_list = std::move(ThreeSubj_standS);
  m_ThreeSubj_CLT_list    = std::move(ThreeSubj_CLT);

  m_SPA_Cutoff = SPA_Cutoff;
  m_zeta = zeta;
  m_tol = tol;
}

arma::mat SPAGRMClass::mgf(double t,
                               const std::vector<arma::vec>& arr_prob_list,
                               double MAF) {
  arma::vec lambda = arma::exp(t * m_resid_unrelated_outliers);

  arma::vec alpha = 1 - MAF + MAF * lambda;
  arma::vec alpha_1 = MAF * m_resid_unrelated_outliers % lambda;
  arma::vec alpha_2 = m_resid_unrelated_outliers % alpha_1;

  arma::vec M_G0_all = alpha % alpha;
  arma::vec M_G1_all = 2 * alpha % alpha_1;
  arma::vec M_G2_all = 2 * (alpha_1 % alpha_1 + alpha % alpha_2);

  int n1 = static_cast<int>(m_TwoSubj_resid_list.size());
  if (n1 != 0) {
    for (int i = 0; i < n1; i++) {
      const arma::vec& Resid = m_TwoSubj_resid_list[i];
      const arma::vec& Rho = m_TwoSubj_rho_list[i];

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

  int n2 = static_cast<int>(arr_prob_list.size());
  if (n2 != 0) {
    for (int i = 0; i < n2; i++) {
      const arma::vec& stand_S = m_ThreeSubj_standS_list[i];
      const arma::vec& arr_prob = arr_prob_list[i];

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

double SPAGRMClass::fastGetRoot(const std::vector<arma::vec>& arr_prob_list,
                                    double Score,
                                    double MAF,
                                    double init_t,
                                    double tol,
                                    int maxiter) {
  double t = init_t;
  arma::vec MGF0; arma::vec MGF1; arma::vec MGF2;
  double CGF1 = 0; double CGF2 = 0;
  double diff_t = std::numeric_limits<double>::infinity();

  double mean = 2 * MAF * m_sum_R_nonOutlier;
  double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier;

  for (int iter = 1; iter < maxiter; iter++) {
    double old_t = t;
    double old_diff_t = diff_t;
    double old_CGF1 = CGF1;

    arma::mat MGF_all = mgf(t, arr_prob_list, MAF);

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

double SPAGRMClass::getProbSpa(const std::vector<arma::vec>& arr_prob_list,
                                double Score,
                                double MAF,
                                bool lower_tail,
                                double zeta,
                                double tol) {
  zeta = fastGetRoot(arr_prob_list, Score, MAF, zeta, tol);

  arma::mat MGF_all = mgf(zeta, arr_prob_list, MAF);

  arma::vec MGF0 = MGF_all.col(0);
  arma::vec MGF1 = MGF_all.col(1);
  arma::vec MGF2 = MGF_all.col(2);

  double mean = 2 * MAF * m_sum_R_nonOutlier;
  double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier;

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

double SPAGRMClass::getMarkerPval(arma::vec GVec,
                                  double altFreq,
                                  double& zScore,
                                  double& hwepval,
                                  double hwepvalCutoff) {
  gethwepval(GVec, hwepval, hwepvalCutoff);

  double MAF = std::min(altFreq, 1 - altFreq);

  double Score = sum(GVec % m_resid) - mean(GVec) * sum(m_resid);
  double G_var = 2 * MAF * (1 - MAF);
  double Score_var = G_var * m_R_GRM_R;
  zScore = Score / sqrt(Score_var);

  if (std::abs(zScore) <= m_SPA_Cutoff) {
    boost::math::normal_distribution<double> ndist2(0.0, 1.0);
    double pval = boost::math::cdf(boost::math::complement(ndist2, std::abs(zScore)));
    return 2 * pval;
  }

  int order2 = arma::index_max(m_MAF_interval >= MAF);
  int order1 = order2 - 1;

  double MAF_ratio = (m_MAF_interval[order2] - MAF) / (m_MAF_interval[order2] - m_MAF_interval[order1]);
  double Var_ThreeOutlier = 0;

  int n1 = static_cast<int>(m_ThreeSubj_standS_list.size());
  std::vector<arma::vec> arr_prob_list;
  arr_prob_list.reserve(n1);

  if (n1 != 0) {
    for (int i = 0; i < n1; i++) {
      const arma::mat& CLT_temp = m_ThreeSubj_CLT_list[i];
      const arma::vec& stand_S = m_ThreeSubj_standS_list[i];
      arma::vec CLT_temp1 = CLT_temp.col(order1);
      arma::vec CLT_temp2 = CLT_temp.col(order2);

      arma::vec arr_prob = MAF_ratio * CLT_temp1 + (1 - MAF_ratio) * CLT_temp2;
      arr_prob_list.push_back(arr_prob);

      arma::vec temp1 = stand_S % arr_prob;

      double temp2 = arma::accu(temp1);
      double Var_ThreeOutlier_temp = arma::accu(stand_S % temp1) - temp2 * temp2;
      Var_ThreeOutlier = Var_ThreeOutlier + Var_ThreeOutlier_temp;
    }
  }

  double Var_nonOutlier = G_var * m_R_GRM_R_nonOutlier;
  double Var_unrelated_outliers = G_var * m_sum_unrelated_outliers2;
  double Var_TwoOutlier = G_var * m_R_GRM_R_TwoSubjOutlier;

  double EmpVar = Var_nonOutlier + Var_unrelated_outliers + Var_TwoOutlier + Var_ThreeOutlier;
  double Var_Ratio = Score_var / EmpVar;
  double Score_adj = Score / sqrt(Var_Ratio);

  double zeta1 = std::abs(Score_adj) / Score_var; zeta1 = std::min(zeta1, 1.2);
  double zeta2 = -std::abs(m_zeta);

  double pval1 = getProbSpa(arr_prob_list, std::abs(Score_adj), MAF, false, zeta1, 1e-4);
  double pval2 = getProbSpa(arr_prob_list, -std::abs(Score_adj), MAF, true, zeta2, m_tol);
  double pval = pval1 + pval2;

  return pval;
}

}
