// SPAmix.cpp -- mtSPAmixClass method implementations

#include <RcppArmadillo.h>

#include "mtSPAmix.h"

mtSPAmixClass::mtSPAmixClass(
  arma::vec resid,
  arma::mat PCs,
  int N,
  double SPA_Cutoff,
  OutlierData outlier
)
  : m_resid(std::move(resid)),
    m_onePlusPCs(arma::join_horiz(arma::ones(N), PCs)),
    m_N(N),
    m_SPA_Cutoff(SPA_Cutoff),
    m_PCs(std::move(PCs)),
    m_sqrt_XTX_inv_diag([&]() {
      arma::mat XTX_inv = arma::inv(arma::mat(m_onePlusPCs.t() * m_onePlusPCs));
      return arma::vec(arma::sqrt(XTX_inv.diag()));
    }()),
    m_outlier(std::move(outlier)),
    m_pval(0.0),
    m_zScore(0.0),
    m_scratch_posValue(m_outlier.posValue.n_elem),
    m_scratch_posOutlier(m_outlier.posOutlier.n_elem),
    m_scratch_posNonOutlier(m_outlier.posNonOutlier.n_elem)
{}

namespace {

struct RootResult {
  double root;
  int iter;
  bool converge;
  double K2;
};

std::pair<double,double> Horg_H2(double t, const arma::vec& R, const arma::vec& MAFVec) {
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
  return {Horg, H2val};
}

std::pair<double,double> H1_adj_H2(double t, const arma::vec& R, double s, const arma::vec& MAFVec) {
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
  return {H1_adj_val, H2val};
}

RootResult fastGetRootK1(
  double initX,
  double s,
  const arma::vec& MAF_outlier,
  double mean_nonOutlier,
  double var_nonOutlier,
  const arma::vec& residOutlier
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

    auto [h1_adj, h2_val] = H1_adj_H2(x, residOutlier, s, MAF_outlier);

    K1 = h1_adj + mean_nonOutlier + var_nonOutlier * x;
    K2 = h2_val + var_nonOutlier;

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

arma::vec logistic_regression(const arma::mat& X, const arma::vec& y) {
  arma::vec beta = nsSPAmix::logistic_regression_beta(X, y);
  arma::mat X_new = arma::join_horiz(arma::ones(X.n_rows), X);
  arma::vec mu = 1.0 / (1.0 + arma::exp(-X_new * beta));
  return 1.0 - arma::sqrt(1.0 - mu);
}

} // anonymous namespace

namespace nsSPAmix {

arma::vec logistic_regression_beta(const arma::mat& X, const arma::vec& y) {
  int n = X.n_rows;
  int p = X.n_cols;

  arma::mat WX_new(n, p + 1);
  arma::mat X_new = arma::join_horiz(arma::ones(n), X);

  arma::vec beta(p+1, arma::fill::zeros);
  double tol = 1e-6;
  int max_iter = 100;
  arma::vec mu(n);

  for (int i = 0; i < max_iter; i++) {
    mu = 1.0 / (1.0 + arma::exp(-X_new * beta));

    arma::vec W = mu % (1.0 - mu);
    arma::vec z = X_new * beta + (y - mu) / W;

    for (int j = 0; j < p+1; j++){
      WX_new.col(j) = X_new.col(j) % W;
    }
    arma::vec beta_new = arma::solve(X_new.t() * WX_new, X_new.t() * (W % z));

    if (arma::norm(beta_new - beta) < tol) {
      beta = beta_new;
      break;
    }
    beta = beta_new;
  }
  return beta;
}

double getProbSpaG(const arma::vec& MAF_outlier,
  const arma::vec& residOutlier,
  double s,
  bool lower_tail,
  double mean_nonOutlier,
  double var_nonOutlier
) {
  double initX = 0;

  RootResult rootRes = fastGetRootK1(initX, s, MAF_outlier, mean_nonOutlier, var_nonOutlier, residOutlier);
  double zeta = rootRes.root;

  auto [k0_val, k2_val] = Horg_H2(zeta, residOutlier, MAF_outlier);
  double k1 = k0_val + mean_nonOutlier * zeta + 0.5 * var_nonOutlier * pow(zeta, 2);
  double k2 = k2_val + var_nonOutlier;

  double temp1 = zeta * s - k1;

  if (!std::isfinite(zeta) || temp1 < 0.0 || k2 <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();

  double w = arma::sign(zeta) * sqrt(2 * temp1);
  double v = zeta * sqrt(k2);

  if (w == 0.0 || v == 0.0 || (v / w) <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();

  double pval = arma::normcdf(arma::sign(lower_tail - 0.5) * (w + log(v/w) / w));
  return pval;
}

} // namespace nsSPAmix

arma::vec mtSPAmixClass::fit_lm(const arma::vec& g, arma::vec& pvalues) {
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

arma::vec mtSPAmixClass::getMafEst(
  const arma::vec& g,
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

double mtSPAmixClass::getMarkerPval(const arma::vec& GVec, double altFreq) {

  arma::vec AFVec = getMafEst(GVec, altFreq);

  const arma::uvec& posValue = m_outlier.posValue;
  const arma::uvec& posOutlier = m_outlier.posOutlier;
  const arma::uvec& posNonOutlier = m_outlier.posNonOutlier;
  const arma::vec& resid = m_outlier.resid;
  const arma::vec& resid2 = m_outlier.resid2;
  const arma::vec& residOutlier = m_outlier.residOutlier;
  const arma::vec& residNonOutlier = m_outlier.residNonOutlier;
  const arma::vec& resid2NonOutlier = m_outlier.resid2NonOutlier;

  // Gather AFVec subsets into pre-allocated scratch (no heap alloc)
  {
    const double* af = AFVec.memptr();
    const arma::uword* idx;
    double* dst;

    idx = posValue.memptr(); dst = m_scratch_posValue.memptr();
    for (arma::uword k = 0; k < posValue.n_elem; ++k) dst[k] = af[idx[k]];
  }
  const arma::vec& AF_posValue = m_scratch_posValue;

  arma::vec GVarVec_posValue = 2.0 * AF_posValue % (1.0 - AF_posValue);

  double S;
  {
    double s = 0.0;
    const double* gp = GVec.memptr();
    const double* rp = resid.memptr();
    const arma::uword* idx = posValue.memptr();
    for (arma::uword k = 0; k < posValue.n_elem; ++k)
      s += gp[idx[k]] * rp[k];
    S = s;
  }
  double VarS = arma::dot(resid2, GVarVec_posValue);
  double S_mean = 2.0 * arma::dot(resid, AF_posValue);
  double zScore = (S - S_mean) / std::sqrt(VarS);

  m_zScore = zScore;

  if (std::abs(zScore) < m_SPA_Cutoff){
    m_pval = arma::normcdf(-1*std::abs(zScore))*2;
    return m_pval;
  }

  // Gather AFVec for outlier and non-outlier positions
  {
    const double* af = AFVec.memptr();
    const arma::uword* idx;
    double* dst;

    idx = posOutlier.memptr(); dst = m_scratch_posOutlier.memptr();
    for (arma::uword k = 0; k < posOutlier.n_elem; ++k) dst[k] = af[idx[k]];

    idx = posNonOutlier.memptr(); dst = m_scratch_posNonOutlier.memptr();
    for (arma::uword k = 0; k < posNonOutlier.n_elem; ++k) dst[k] = af[idx[k]];
  }
  const arma::vec& MAF_outlier = m_scratch_posOutlier;
  const arma::vec& MAF_nonOutlier = m_scratch_posNonOutlier;

  double mean_nonOutlier = arma::dot(residNonOutlier, MAF_nonOutlier) * 2.0;
  double var_nonOutlier = arma::sum(resid2NonOutlier % MAF_nonOutlier % (1.0 - MAF_nonOutlier)) * 2.0;

  double pval1 = nsSPAmix::getProbSpaG(
    MAF_outlier,
    residOutlier,
    std::abs(S-S_mean)+S_mean,
    false,
    mean_nonOutlier,
    var_nonOutlier
  );

  double pval2 = nsSPAmix::getProbSpaG(
    MAF_outlier,
    residOutlier,
    -1*std::abs(S-S_mean)+S_mean,
    true,
    mean_nonOutlier,
    var_nonOutlier
  );

  m_pval = pval1 + pval2;
  return m_pval;
}

