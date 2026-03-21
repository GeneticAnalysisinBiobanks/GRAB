// SPAmix.cpp -- SPAmixClass method implementations

#include <RcppArmadillo.h>

#include "mtSPAmix.h"

namespace SPAmix {

SPAmixClass::SPAmixClass(
  arma::mat resid,
  arma::mat PCs,
  int N,
  double SPA_Cutoff,
  std::vector<SPAmixClass::OutlierData> outlierVec
)
  : m_resid(std::move(resid)),
    m_onePlusPCs(arma::join_horiz(arma::ones(N), PCs)),
    m_N(N),
    m_Npheno(static_cast<int>(m_resid.n_cols)),
    m_SPA_Cutoff(SPA_Cutoff),
    m_PCs(std::move(PCs)),
    m_sqrt_XTX_inv_diag([&]() {
      arma::mat XTX_inv = arma::inv(arma::mat(m_onePlusPCs.t() * m_onePlusPCs));
      return arma::vec(arma::sqrt(XTX_inv.diag()));
    }()),
    m_diffTime1(arma::zeros(2)),
    m_diffTime2(arma::zeros(2)),
    m_outlierVec(std::move(outlierVec))
{
  m_pvalVec.resize(m_Npheno);
  m_zScoreVec.resize(m_Npheno);
}

arma::vec SPAmixClass::M_G0(arma::vec t, arma::vec MAF) {
  arma::vec re = pow((1 - MAF + MAF % arma::exp(t)), 2);
  return re;
}

arma::vec SPAmixClass::M_G1(arma::vec t, arma::vec MAF) {
  arma::vec re = 2 * (MAF % arma::exp(t)) % (1 - MAF + MAF % arma::exp(t));
  return re;
}

arma::vec SPAmixClass::M_G2(arma::vec t, arma::vec MAF) {
  arma::vec re = 2 * pow(MAF % arma::exp(t), 2) + 2 * (MAF % arma::exp(t)) % (1 - MAF + MAF % arma::exp(t));
  return re;
}

arma::vec SPAmixClass::K_G0(arma::vec t, arma::vec MAF) {
  arma::vec re = arma::log(M_G0(t, MAF));
  return re;
}

arma::vec SPAmixClass::K_G1(arma::vec t, arma::vec MAF) {
  arma::vec re = M_G1(t, MAF) / M_G0(t, MAF);
  return re;
}

arma::vec SPAmixClass::K_G2(arma::vec t, arma::vec MAF) {
  arma::vec re = (M_G0(t, MAF) % M_G2(t, MAF) - pow(M_G1(t, MAF), 2)) / pow(M_G0(t, MAF), 2);
  return re;
}

double SPAmixClass::H_org(double t, arma::vec R, const arma::vec& MAFVec) {
  double out = sum(K_G0(t * R, MAFVec));
  return out;
}

double SPAmixClass::H1_adj(double t, arma::vec R, const double& s, const arma::vec& MAFVec) {
  double out = sum(R % K_G1(t * R, MAFVec)) - s;
  return out;
}

double SPAmixClass::H2(double t, arma::vec R, const arma::vec& MAFVec) {
  double out = sum(pow(R, 2) % K_G2(t * R, MAFVec));
  return out;
}

arma::vec SPAmixClass::Horg_H2(double t, arma::vec R, const arma::vec MAFVec) {
  arma::vec Horg_H2_vec(2);
  arma::vec tR = t * R;
  arma::vec exp_tR = arma::exp(tR);
  arma::vec MAF_exp_tR = MAFVec % exp_tR;
  arma::vec M_G0_vec = pow((1 - MAFVec + MAF_exp_tR), 2);
  arma::vec M_G1_vec = 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
  arma::vec M_G2_vec = 2 * pow(MAF_exp_tR, 2) + 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
  arma::vec K_G0_vec = arma::log(M_G0_vec);
  arma::vec K_G2_vec = (M_G0_vec % M_G2_vec - pow(M_G1_vec, 2)) / pow(M_G0_vec, 2);
  double Horg = sum(K_G0_vec);
  double H2val = sum(pow(R, 2) % K_G2_vec);
  Horg_H2_vec.at(0) = Horg;
  Horg_H2_vec.at(1) = H2val;
  return Horg_H2_vec;
}

arma::vec SPAmixClass::H1_adj_H2(double t, arma::vec R, double s, const arma::vec MAFVec) {
  arma::vec H1_adj_H2_vec(2);
  arma::vec tR = t * R;
  arma::vec exp_tR = arma::exp(tR);
  arma::vec MAF_exp_tR = MAFVec % exp_tR;
  arma::vec M_G0_vec = pow((1 - MAFVec + MAF_exp_tR), 2);
  arma::vec M_G1_vec = 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
  arma::vec M_G2_vec = 2 * pow(MAF_exp_tR, 2) + 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
  arma::vec K_G1_vec = M_G1_vec / M_G0_vec;
  arma::vec K_G2_vec = (M_G0_vec % M_G2_vec - pow(M_G1_vec, 2)) / pow(M_G0_vec, 2);
  double H1_adj_val = sum(R % K_G1_vec) - s;
  double H2val = sum(pow(R, 2) % K_G2_vec);
  H1_adj_H2_vec.at(0) = H1_adj_val;
  H1_adj_H2_vec.at(1) = H2val;
  return H1_adj_H2_vec;
}

SPAmixClass::RootResult SPAmixClass::fastGetRootK1(
  double initX,
  const double& s,
  const arma::vec MAF_outlier,
  double mean_nonOutlier,
  double var_nonOutlier,
  const arma::vec residOutlier
) {
  double x = initX, oldX;
  double K1 = 0, K2 = 0, oldK1;
  double diffX = arma::datum::inf, oldDiffX;
  bool converge = true;
  double tol = 0.001;
  int maxiter = 100;
  int iter = 0;

  for (iter = 0; iter < maxiter; iter++){

    oldX = x;
    oldDiffX = diffX;
    oldK1 = K1;

    arma::vec H1_adj_H2_vec = H1_adj_H2(x, residOutlier, s, MAF_outlier);

    K1 = H1_adj_H2_vec.at(0) + mean_nonOutlier + var_nonOutlier * x;
    K2 = H1_adj_H2_vec.at(1) + var_nonOutlier;

    diffX = -1 * K1 / K2;

    if (!std::isfinite(K1)){
      x = arma::datum::inf;
      K2 = 0;
      break;
    }

    if (arma::sign(K1) != arma::sign(oldK1)){
      while (std::abs(diffX) > std::abs(oldDiffX) - tol){
        diffX = diffX / 2;
      }
    }

    if (std::abs(diffX) < tol) break;

    x = oldX + diffX;
  }

  if (iter == maxiter)
    converge = false;

  return {x, iter, converge, K2};
}

double SPAmixClass::getProbSpaG(const arma::vec MAF_outlier,
  const arma::vec residOutlier,
  double s,
  bool lower_tail,
  double mean_nonOutlier,
  double var_nonOutlier
) {
  double initX = 0;

  RootResult rootRes = fastGetRootK1(initX, s, MAF_outlier, mean_nonOutlier, var_nonOutlier, residOutlier);
  double zeta = rootRes.root;

  arma::vec k12 = Horg_H2(zeta, residOutlier, MAF_outlier);
  double k1 = k12.at(0) + mean_nonOutlier * zeta + 0.5 * var_nonOutlier * pow(zeta, 2);
  double k2 = k12.at(1) + var_nonOutlier;

  double temp1 = zeta * s - k1;

  double w = arma::sign(zeta) * sqrt(2 * temp1);
  double v = zeta * sqrt(k2);

  double pval = arma::normcdf(arma::sign(lower_tail - 0.5) * (w + log(v/w) / w));
  return pval;
}

arma::vec SPAmixClass::simulate_uniform(int n, double lower, double upper) {
  arma::vec vec(n);
  vec.randu();
  vec = lower + (upper - lower) * vec;
  return vec;
}

arma::vec SPAmixClass::fit_lm(const arma::vec& g, arma::vec& pvalues) {
  int n = m_N;
  int k = m_PCs.n_cols;

  arma::vec coef = arma::solve(m_onePlusPCs, g);
  arma::vec fittedValues = m_onePlusPCs * coef;

  double s2 = sum(square(g - fittedValues)) / (n - k - 1);

  arma::vec se = m_sqrt_XTX_inv_diag * sqrt(s2);
  arma::vec t = coef / se;

  for (int i = 0; i < k; i++){
    boost::math::students_t_distribution<double> tdist(n - k - 1);
    pvalues[i] = 2.0 * boost::math::cdf(boost::math::complement(tdist, std::abs(t[i+1])));
  }
  return fittedValues;
}

arma::vec SPAmixClass::logistic_regression(const arma::mat& X, const arma::vec& y) {
  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat WX_new(n, p + 1);
  arma::mat X_new = arma::join_horiz(arma::ones(n), X);

  arma::vec beta(p+1, arma::fill::zeros);
  double tol = 1e-6;
  int max_iter = 100;
  arma::vec mu(n);

  for (int i = 0; i < max_iter; i++) {
    mu = 1 / (1 + exp(-X_new * beta));

    arma::vec W = mu % (1 - mu);
    arma::vec z = X_new * beta + (y - mu) / W;

    for (int j = 0; j < p+1; j++){
      WX_new.col(j) = X_new.col(j) % W;
    }
    arma::vec beta_new = inv(X_new.t() * WX_new) * X_new.t() * (W % z);

    if (norm(beta_new - beta) < tol) {
      break;
    }
    beta = beta_new;
  }

  arma::vec MAFest = 1 - sqrt(1 - mu);

  return MAFest;
}

arma::vec SPAmixClass::getMafEst(
  arma::vec g,
  double altFreq,
  double MAC_cutoff,
  double PCs_pvalue_cutoff,
  double MAF_est_negative_ratio_cutoff
) {
  int N = g.n_elem;
  arma::vec g0(N, arma::fill::zeros);
  arma::vec MAF_all = arma::vec(N, arma::fill::value(altFreq));
  double MAC = altFreq * 2 * N;

  int PC_number = m_PCs.n_cols;
  double MAF0 = 0;

  arma::vec pvalues(PC_number);
  arma::vec MAF_est, topPCs_pvalueVec;

  if (MAC <= MAC_cutoff){
    MAF_est = MAF_all;
  }else{
    arma::vec fit = fit_lm(g, pvalues);
    fit = fit / 2;

    arma::uvec posZero = arma::find(fit < 0);
    arma::uvec posOne = arma::find(fit > 1);

    int nError = posZero.n_elem + posOne.n_elem;
    double propError = (double)nError / N;

    if (propError < MAF_est_negative_ratio_cutoff){
      fit.elem(posZero).fill(MAF0);
      fit.elem(posOne).fill(1-MAF0);
      MAF_est = fit;
    }else{
      arma::uvec posSigPCs = arma::find(pvalues < PCs_pvalue_cutoff);

      if (posSigPCs.n_elem == 0){
        MAF_est = MAF_all;
      }else{
        arma::mat sigPCs = m_PCs.cols(posSigPCs);
        arma::uvec posg12 = arma::find(g > 0.5);
        g0.elem(posg12).fill(1);

        double MAC_after = sum(g0);
        if (MAC_after <= MAC_cutoff){
          MAF_est = MAF_all;
        }else{
          MAF_est = logistic_regression(sigPCs, g0);
        }
      }
    }
  }

  return MAF_est;
}

double SPAmixClass::getMarkerPval(arma::vec GVec, double altFreq) {

  arma::vec AFVec = getMafEst(GVec, altFreq);
  arma::vec GVarVec = 2 * AFVec % (1 - AFVec);

  for (int i = 0; i < m_Npheno; i++){
    const OutlierData& tempOutlierList = m_outlierVec[i];
    const arma::uvec& posValue = tempOutlierList.posValue;
    const arma::uvec& posOutlier = tempOutlierList.posOutlier;
    const arma::uvec& posNonOutlier = tempOutlierList.posNonOutlier;
    const arma::vec& resid = tempOutlierList.resid;
    const arma::vec& resid2 = tempOutlierList.resid2;
    const arma::vec& residOutlier = tempOutlierList.residOutlier;
    const arma::vec& residNonOutlier = tempOutlierList.residNonOutlier;
    const arma::vec& resid2NonOutlier = tempOutlierList.resid2NonOutlier;

    double S = sum(GVec.elem(posValue) % resid);
    double VarS = sum(resid2 % GVarVec.elem(posValue));

    double S_mean = 2 * sum(resid % AFVec.elem(posValue));
    double zScore = (S-S_mean) / sqrt(VarS);

    m_zScoreVec.at(i) = zScore;

    if (std::abs(zScore) < m_SPA_Cutoff){
      m_pvalVec.at(i) = arma::normcdf(-1*std::abs(zScore))*2;
      continue;
    }

    arma::vec MAF_outlier = AFVec.elem(posOutlier);
    arma::vec MAF_nonOutlier = AFVec.elem(posNonOutlier);

    double mean_nonOutlier = sum(residNonOutlier % MAF_nonOutlier) * 2;
    double var_nonOutlier = sum(resid2NonOutlier % MAF_nonOutlier % (1-MAF_nonOutlier)) * 2;

    double pval1 = getProbSpaG(
      MAF_outlier,
      residOutlier,
      std::abs(S-S_mean)+S_mean,
      false,
      mean_nonOutlier,
      var_nonOutlier
    );

    double pval2 = getProbSpaG(
      MAF_outlier,
      residOutlier,
      -1*std::abs(S-S_mean)+S_mean,
      true,
      mean_nonOutlier,
      var_nonOutlier
    );

    m_pvalVec.at(i) = pval1 + pval2;
  }

  double pval = 0;
  return pval;
}

}
