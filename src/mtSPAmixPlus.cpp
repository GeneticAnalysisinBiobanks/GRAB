// SPAmixPlus.cpp -- mtSPAmixPlusClass and namespace-scope SPA helpers

#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp>

#include "mtSPAmixPlus.h"

mtSPAmixPlusClass::mtSPAmixPlusClass(
  arma::vec resid,
  arma::mat PCs,
  int N,
  double SPA_Cutoff,
  OutlierData outlier,
  std::vector<std::tuple<int, int, double>> sparseTriplets,
  const std::string& afFilePath,
  const std::string& afFilePrecision
) {

  m_afFilePath = afFilePath;
  m_afFilePrecision = afFilePrecision;
  m_sparseTriplets = std::move(sparseTriplets);

  m_resid = std::move(resid);
  m_PCs = std::move(PCs);
  m_N = N;
  m_SPA_Cutoff = SPA_Cutoff;

  m_pval = 0.0;
  m_zScore = 0.0;
  m_Beta = 0.0;
  m_S = 0.0;
  m_Smean = 0.0;
  m_VarS = 0.0;

  m_outlier = std::move(outlier);
  m_onePlusPCs = arma::join_horiz(arma::ones(N), m_PCs);
  arma::mat X_t = m_onePlusPCs.t();
  arma::mat XTX = X_t * m_onePlusPCs;
  arma::mat XTX_inv = arma::inv(XTX);
  m_sqrt_XTX_inv_diag = arma::sqrt(XTX_inv.diag());

}

arma::vec mtSPAmixPlusClass::fit_lm_get_beta(const arma::vec& g, arma::vec& pvalues) {
  int n = m_N;
  int k = m_PCs.n_cols;

  arma::vec coef = arma::solve(m_onePlusPCs, g);
  arma::vec fittedValues = m_onePlusPCs * coef;
  double s2 = arma::sum(arma::square(g - fittedValues)) / (n - k - 1);
  arma::vec se = m_sqrt_XTX_inv_diag * std::sqrt(s2);
  arma::vec t = coef / se;

  for (int i = 0; i < k; i++){
    boost::math::students_t_distribution<double> tdist(n - k - 1);
    pvalues[i] = 2.0 * boost::math::cdf(boost::math::complement(tdist, std::abs(t[i+1])));
  }
  return coef;
}

arma::vec mtSPAmixPlusClass::fit_lm(const arma::vec& g, arma::vec& pvalues) {
    arma::vec coef = fit_lm_get_beta(g, pvalues);
    return m_onePlusPCs * coef;
}

mtSPAmixPlusClass::AFModelInfo mtSPAmixPlusClass::computeAFModel(arma::vec GVec, double altFreq) {
  AFModelInfo model;
  model.status = 0;
  int k = m_PCs.n_cols;
  model.betas = arma::vec(k + 1, arma::fill::zeros);

  double MAC_cutoff = 20;
  double PCs_pvalue_cutoff = 0.05;
  double MAF_est_negative_ratio_cutoff = 0.1;

  int N = m_N;
  arma::vec g = GVec;
  double MAC = altFreq * 2.0 * N;
  arma::vec pvalues(k);

  if (MAC <= MAC_cutoff){
      model.status = 0;

  }else{
    arma::vec coef_lm = fit_lm_get_beta(g, pvalues);
    arma::vec fit = m_onePlusPCs * coef_lm;
    fit = fit / 2.0;

    arma::uvec posZero = arma::find(fit < 0);
    arma::uvec posOne = arma::find(fit > 1);

    int nError = posZero.n_elem + posOne.n_elem;
    double propError = (double)nError / N;

    if (propError < MAF_est_negative_ratio_cutoff){
      model.status = 1;
      model.betas = coef_lm;
    }else{
      arma::uvec posSigPCs = arma::find(pvalues < PCs_pvalue_cutoff);
      if (posSigPCs.n_elem == 0){
          model.status = 0;
      } else {
        arma::mat sigPCs = m_PCs.cols(posSigPCs);
        arma::vec g0(N, arma::fill::zeros);
        arma::uvec posg12 = arma::find(g > 0.5);
        g0.elem(posg12).fill(1);

        double MAC_after = sum(g0);
        if (MAC_after <= MAC_cutoff){
          model.status = 0;
        }else{
          model.status = 2;
          arma::vec sub_beta = nsSPAmix::logistic_regression_beta(sigPCs, g0);
          model.betas(0) = sub_beta(0);
          for (unsigned int i=0; i < posSigPCs.n_elem; ++i){
            int pc_index = posSigPCs(i);
            model.betas(pc_index + 1) = sub_beta(i + 1);
          }
        }
      }
    }
  }
  return model;
}

arma::vec mtSPAmixPlusClass::getAFFromModel(AFModelInfo model, double altFreq) {
  if (model.status == 0){
      return arma::vec(m_N, arma::fill::value(altFreq));
  }
  if (model.status == 1){
      arma::vec fit = m_onePlusPCs * model.betas;
      fit = fit / 2.0;
      arma::uvec posZero = arma::find(fit < 0);
      arma::uvec posOne = arma::find(fit > 1);
      fit.elem(posZero).fill(0.0);
      fit.elem(posOne).fill(1.0);
      return fit;
  }
  if (model.status == 2){
      arma::vec linear_pred = m_onePlusPCs * model.betas;
      arma::vec mu = 1.0 / (1.0 + arma::exp(-linear_pred));
      return 1.0 - arma::sqrt(1.0 - mu);
  }
  return arma::vec(m_N, arma::fill::value(altFreq));
}

arma::vec mtSPAmixPlusClass::getMafEst(arma::vec GVec, double altFreq) {
  AFModelInfo model = computeAFModel(GVec, altFreq);
  return getAFFromModel(model, altFreq);
}

double mtSPAmixPlusClass::getMarkerPval(arma::vec GVec, double altFreq) {
  arma::vec AFVec = getMafEst(GVec, altFreq);
  m_MAFVec = AFVec;
  arma::vec GVarVec = 2.0 * AFVec % (1.0 - AFVec);

  const auto& od = m_outlier;
  const arma::uvec& posValue = od.posValue;
  const arma::uvec& posOutlier = od.posOutlier;
  const arma::uvec& posNonOutlier = od.posNonOutlier;

  const arma::vec& resid = od.resid;
  const arma::vec& residOutlier = od.residOutlier;
  const arma::vec& residNonOutlier = od.residNonOutlier;
  const arma::vec& resid2NonOutlier = od.resid2NonOutlier;

  const arma::vec& R_subset = resid;
  arma::vec GVar_subset = GVarVec.elem(posValue);
  arma::vec R_new = R_subset % arma::sqrt(GVar_subset);

  double VarS = calculateSparseVariance(R_new, posValue);
  double S = arma::sum(GVec.elem(posValue) % R_subset);
  double S_mean = 2.0 * arma::sum(R_subset % AFVec.elem(posValue));
  double zScore = (S - S_mean) / std::sqrt(VarS);

  m_zScore = zScore;
  m_Beta = (S - S_mean) / VarS;
  m_S = S;
  m_Smean = S_mean;
  m_VarS = VarS;

  if (std::abs(zScore) < m_SPA_Cutoff){
    m_pval = arma::normcdf(-1.0*std::abs(zScore))*2.0;
    return m_pval;
  }

  double S_var_SPAmix = arma::sum(arma::square(R_subset) % GVar_subset);
  double Var_ratio = S_var_SPAmix / VarS;
  double S_new = S * std::sqrt(Var_ratio);
  double S_mean_new = S_mean * std::sqrt(Var_ratio);
  double S_upper = std::max(S_new, 2.0*S_mean_new - S_new);
  double S_lower = std::min(S_new, 2.0*S_mean_new - S_new);

  arma::vec MAF_outlier = AFVec.elem(posOutlier);
  double mean_nonOutlier = arma::sum(residNonOutlier % AFVec.elem(posNonOutlier)) * 2.0;
  double var_nonOutlier = arma::sum(resid2NonOutlier % AFVec.elem(posNonOutlier) % (1.0 - AFVec.elem(posNonOutlier))) * 2.0;
  double pval1 = nsSPAmix::getProbSpaG(MAF_outlier, residOutlier, S_upper, false, mean_nonOutlier, var_nonOutlier);
  double pval2 = nsSPAmix::getProbSpaG(MAF_outlier, residOutlier, S_lower, true, mean_nonOutlier, var_nonOutlier);

  m_pval = pval1 + pval2;
  return m_pval;
}

double mtSPAmixPlusClass::getMarkerPvalFromModel(arma::vec GVec, AFModelInfo model, double altFreq) {
  arma::vec AFVec = getAFFromModel(model, altFreq);
  m_MAFVec = AFVec;
  arma::vec GVarVec = 2.0 * AFVec % (1.0 - AFVec);

  const auto& od = m_outlier;
  const arma::uvec& posValue = od.posValue;
  const arma::uvec& posOutlier = od.posOutlier;
  const arma::uvec& posNonOutlier = od.posNonOutlier;

  const arma::vec& resid = od.resid;
  const arma::vec& residOutlier = od.residOutlier;
  const arma::vec& residNonOutlier = od.residNonOutlier;
  const arma::vec& resid2NonOutlier = od.resid2NonOutlier;

  arma::vec R_subset = resid;
  arma::vec GVar_subset = GVarVec.elem(posValue);
  arma::vec R_new = R_subset % arma::sqrt(GVar_subset);
  double VarS = calculateSparseVariance(R_new, posValue);

  double S = arma::sum(GVec.elem(posValue) % R_subset);
  double S_mean = 2.0 * arma::sum(R_subset % AFVec.elem(posValue));
  double zScore = (S - S_mean) / std::sqrt(VarS);

  m_zScore = zScore;
  m_Beta = (S - S_mean) / VarS;
  m_S = S;
  m_Smean = S_mean;
  m_VarS = VarS;

  if (std::abs(zScore) < m_SPA_Cutoff){
    boost::math::normal_distribution<double> ndist(0.0, 1.0);
    m_pval = 2.0 * boost::math::cdf(boost::math::complement(ndist, std::abs(zScore)));
    return m_pval;
  }

  double S_var_SPAmix = arma::sum(arma::square(R_subset) % GVar_subset);
  double Var_ratio = S_var_SPAmix / VarS;
  double S_new = S * std::sqrt(Var_ratio);
  double S_mean_new = S_mean * std::sqrt(Var_ratio);

  double S_upper = std::max(S_new, 2.0*S_mean_new - S_new);
  double S_lower = std::min(S_new, 2.0*S_mean_new - S_new);
  arma::vec MAF_outlier = AFVec.elem(posOutlier);

  double mean_nonOutlier = arma::sum(residNonOutlier % AFVec.elem(posNonOutlier)) * 2.0;
  double var_nonOutlier = arma::sum(resid2NonOutlier % AFVec.elem(posNonOutlier) % (1.0 - AFVec.elem(posNonOutlier))) * 2.0;
  double pval1 = nsSPAmix::getProbSpaG(MAF_outlier, residOutlier, S_upper, false, mean_nonOutlier, var_nonOutlier);
  double pval2 = nsSPAmix::getProbSpaG(MAF_outlier, residOutlier, S_lower, true, mean_nonOutlier, var_nonOutlier);

  m_pval = pval1 + pval2;
  return m_pval;
}

double mtSPAmixPlusClass::calculateSparseVariance(const arma::vec& R_new, const arma::uvec& posValue) {
  double covSum = 0.0;
  std::unordered_set<int> validIndices;
  for (auto idx : posValue) {
      validIndices.insert(static_cast<int>(idx));
  }
  for (const auto& triplet : m_sparseTriplets) {
    int i = std::get<0>(triplet);
    int j = std::get<1>(triplet);

    if (validIndices.count(i) && validIndices.count(j)) {
      double grmValue = std::get<2>(triplet);
      covSum += grmValue * R_new(i) * R_new(j);
    }
  }
  return 2.0 * covSum - arma::sum(arma::square(R_new));
}
