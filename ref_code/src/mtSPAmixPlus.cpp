// SPAmixPlus.cpp -- mtSPAmixPlusClass and namespace-scope SPA helpers

#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
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
  std::string afFilePath,
  std::string afFilePrecision
)
  : m_N(N),
    m_SPA_Cutoff(SPA_Cutoff),
    m_resid(std::move(resid)),
    m_PCs(std::move(PCs)),
    m_onePlusPCs(arma::join_horiz(arma::ones(m_N), m_PCs)),
    m_sqrt_XTX_inv_diag([this]() {
      arma::mat XTX_inv = arma::inv(m_onePlusPCs.t() * m_onePlusPCs);
      return arma::vec(arma::sqrt(XTX_inv.diag()));
    }()),
    m_outlier(std::move(outlier)),
    m_sparseTriplets(std::move(sparseTriplets)),
    m_filteredTriplets([this]() {
      // Build global-to-local index map for posValue
      std::unordered_map<int, arma::uword> globalToLocal;
      const arma::uvec& pv = m_outlier.posValue;
      globalToLocal.reserve(pv.n_elem);
      for (arma::uword k = 0; k < pv.n_elem; ++k)
        globalToLocal[static_cast<int>(pv[k])] = k;
      // Filter triplets to only those with both indices in posValue
      std::vector<FilteredTriplet> ft;
      for (const auto& [i, j, grm] : m_sparseTriplets) {
        auto it_i = globalToLocal.find(i);
        if (it_i == globalToLocal.end()) continue;
        auto it_j = globalToLocal.find(j);
        if (it_j == globalToLocal.end()) continue;
        ft.push_back({it_i->second, it_j->second, grm});
      }
      return ft;
    }()),
    m_afFilePath(std::move(afFilePath)),
    m_afFilePrecision(std::move(afFilePrecision)),
    m_pval(0.0),
    m_zScore(0.0),
    m_Beta(0.0),
    m_S(0.0),
    m_Smean(0.0),
    m_VarS(0.0),
    m_scratch_posValue(m_outlier.posValue.n_elem),
    m_scratch_posOutlier(m_outlier.posOutlier.n_elem),
    m_scratch_posNonOutlier(m_outlier.posNonOutlier.n_elem)
{}

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

mtSPAmixPlusClass::AFModelInfo mtSPAmixPlusClass::computeAFModel(const arma::vec& GVec, double altFreq) {
  AFModelInfo model;
  model.status = 0;
  int k = m_PCs.n_cols;
  model.betas = arma::vec(k + 1, arma::fill::zeros);

  double MAC_cutoff = 20;
  double PCs_pvalue_cutoff = 0.05;
  double MAF_est_negative_ratio_cutoff = 0.1;

  int N = m_N;
  const arma::vec& g = GVec;
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

arma::vec mtSPAmixPlusClass::getAFFromModel(const AFModelInfo& model, double altFreq) {
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

arma::vec mtSPAmixPlusClass::getMafEst(const arma::vec& GVec, double altFreq) {
  AFModelInfo model = computeAFModel(GVec, altFreq);
  return getAFFromModel(model, altFreq);
}

double mtSPAmixPlusClass::getMarkerPval(const arma::vec& GVec, double altFreq) {
  arma::vec AFVec = getMafEst(GVec, altFreq);
  m_MAFVec = AFVec;

  const auto& od = m_outlier;
  const arma::uvec& posValue      = od.posValue;
  const arma::uvec& posOutlier    = od.posOutlier;
  const arma::uvec& posNonOutlier = od.posNonOutlier;
  const arma::uword nPV  = posValue.n_elem;
  const arma::uword nPO  = posOutlier.n_elem;
  const arma::uword nPNO = posNonOutlier.n_elem;

  const arma::vec& resid             = od.resid;
  const arma::vec& residOutlier      = od.residOutlier;
  const arma::vec& residNonOutlier   = od.residNonOutlier;
  const arma::vec& resid2NonOutlier  = od.resid2NonOutlier;

  // Gather AFVec for posValue into scratch
  const double* afp = AFVec.memptr();
  const double* gp  = GVec.memptr();
  const arma::uword* pvp  = posValue.memptr();
  const arma::uword* pop  = posOutlier.memptr();
  const arma::uword* pnop = posNonOutlier.memptr();
  double* sv  = m_scratch_posValue.memptr();
  double* so  = m_scratch_posOutlier.memptr();
  double* sno = m_scratch_posNonOutlier.memptr();

  for (arma::uword k = 0; k < nPV; ++k)  sv[k]  = afp[pvp[k]];
  for (arma::uword k = 0; k < nPO; ++k)  so[k]  = afp[pop[k]];
  for (arma::uword k = 0; k < nPNO; ++k) sno[k] = afp[pnop[k]];

  // Compute GVar_subset and R_new for posValue (local indices)
  const double* rp = resid.memptr();
  arma::vec R_new(nPV);
  double* rnp = R_new.memptr();
  double S = 0.0, S_mean = 0.0;
  double S_var_SPAmix = 0.0;
  for (arma::uword k = 0; k < nPV; ++k) {
    double af = sv[k];
    double gvar = 2.0 * af * (1.0 - af);
    rnp[k] = rp[k] * std::sqrt(gvar);
    S      += gp[pvp[k]] * rp[k];
    S_mean += rp[k] * af;
    S_var_SPAmix += rp[k] * rp[k] * gvar;
  }
  S_mean *= 2.0;

  double VarS = calculateSparseVariance(R_new);
  double zScore = (S - S_mean) / std::sqrt(VarS);

  m_zScore = zScore;
  m_Beta   = (S - S_mean) / VarS;
  m_S      = S;
  m_Smean  = S_mean;
  m_VarS   = VarS;

  if (std::abs(zScore) < m_SPA_Cutoff) {
    m_pval = arma::normcdf(-1.0 * std::abs(zScore)) * 2.0;
    return m_pval;
  }

  double Var_ratio    = S_var_SPAmix / VarS;
  double sqrt_ratio   = std::sqrt(Var_ratio);
  double S_new        = S * sqrt_ratio;
  double S_mean_new   = S_mean * sqrt_ratio;
  double S_upper      = std::max(S_new, 2.0 * S_mean_new - S_new);
  double S_lower      = std::min(S_new, 2.0 * S_mean_new - S_new);

  // MAF_outlier already gathered into m_scratch_posOutlier
  const arma::vec& MAF_outlier = m_scratch_posOutlier;

  // Compute mean_nonOutlier and var_nonOutlier using gathered scratch
  const double* rnop = residNonOutlier.memptr();
  const double* r2nop = resid2NonOutlier.memptr();
  double mean_nonOutlier = 0.0, var_nonOutlier = 0.0;
  for (arma::uword k = 0; k < nPNO; ++k) {
    double af = sno[k];
    mean_nonOutlier += rnop[k] * af;
    var_nonOutlier  += r2nop[k] * af * (1.0 - af);
  }
  mean_nonOutlier *= 2.0;
  var_nonOutlier  *= 2.0;

  double pval1 = nsSPAmix::getProbSpaG(MAF_outlier, residOutlier, S_upper, false, mean_nonOutlier, var_nonOutlier);
  double pval2 = nsSPAmix::getProbSpaG(MAF_outlier, residOutlier, S_lower, true,  mean_nonOutlier, var_nonOutlier);

  m_pval = pval1 + pval2;
  return m_pval;
}

double mtSPAmixPlusClass::getMarkerPvalFromModel(const arma::vec& GVec, const AFModelInfo& model, double altFreq) {
  arma::vec AFVec = getAFFromModel(model, altFreq);
  m_MAFVec = AFVec;

  const auto& od = m_outlier;
  const arma::uvec& posValue      = od.posValue;
  const arma::uvec& posOutlier    = od.posOutlier;
  const arma::uvec& posNonOutlier = od.posNonOutlier;
  const arma::uword nPV  = posValue.n_elem;
  const arma::uword nPO  = posOutlier.n_elem;
  const arma::uword nPNO = posNonOutlier.n_elem;

  const arma::vec& resid             = od.resid;
  const arma::vec& residOutlier      = od.residOutlier;
  const arma::vec& residNonOutlier   = od.residNonOutlier;
  const arma::vec& resid2NonOutlier  = od.resid2NonOutlier;

  // Gather AFVec at index positions into scratch buffers
  const double* afp  = AFVec.memptr();
  const double* gp   = GVec.memptr();
  const arma::uword* pvp  = posValue.memptr();
  const arma::uword* pop  = posOutlier.memptr();
  const arma::uword* pnop = posNonOutlier.memptr();
  double* sv  = m_scratch_posValue.memptr();
  double* so  = m_scratch_posOutlier.memptr();
  double* sno = m_scratch_posNonOutlier.memptr();

  for (arma::uword k = 0; k < nPV; ++k)  sv[k]  = afp[pvp[k]];
  for (arma::uword k = 0; k < nPO; ++k)  so[k]  = afp[pop[k]];
  for (arma::uword k = 0; k < nPNO; ++k) sno[k] = afp[pnop[k]];

  // Compute R_new, S, S_mean, S_var_SPAmix in a single pass
  const double* rp = resid.memptr();
  arma::vec R_new(nPV);
  double* rnp = R_new.memptr();
  double S = 0.0, S_mean = 0.0, S_var_SPAmix = 0.0;
  for (arma::uword k = 0; k < nPV; ++k) {
    double af = sv[k];
    double gvar = 2.0 * af * (1.0 - af);
    rnp[k] = rp[k] * std::sqrt(gvar);
    S      += gp[pvp[k]] * rp[k];
    S_mean += rp[k] * af;
    S_var_SPAmix += rp[k] * rp[k] * gvar;
  }
  S_mean *= 2.0;

  double VarS = calculateSparseVariance(R_new);
  double zScore = (S - S_mean) / std::sqrt(VarS);

  m_zScore = zScore;
  m_Beta   = (S - S_mean) / VarS;
  m_S      = S;
  m_Smean  = S_mean;
  m_VarS   = VarS;

  if (std::abs(zScore) < m_SPA_Cutoff) {
    boost::math::normal_distribution<double> ndist(0.0, 1.0);
    m_pval = 2.0 * boost::math::cdf(boost::math::complement(ndist, std::abs(zScore)));
    return m_pval;
  }

  double Var_ratio  = S_var_SPAmix / VarS;
  double sqrt_ratio = std::sqrt(Var_ratio);
  double S_new      = S * sqrt_ratio;
  double S_mean_new = S_mean * sqrt_ratio;
  double S_upper    = std::max(S_new, 2.0 * S_mean_new - S_new);
  double S_lower    = std::min(S_new, 2.0 * S_mean_new - S_new);

  const arma::vec& MAF_outlier = m_scratch_posOutlier;

  const double* rnop  = residNonOutlier.memptr();
  const double* r2nop = resid2NonOutlier.memptr();
  double mean_nonOutlier = 0.0, var_nonOutlier = 0.0;
  for (arma::uword k = 0; k < nPNO; ++k) {
    double af = sno[k];
    mean_nonOutlier += rnop[k] * af;
    var_nonOutlier  += r2nop[k] * af * (1.0 - af);
  }
  mean_nonOutlier *= 2.0;
  var_nonOutlier  *= 2.0;

  double pval1 = nsSPAmix::getProbSpaG(MAF_outlier, residOutlier, S_upper, false, mean_nonOutlier, var_nonOutlier);
  double pval2 = nsSPAmix::getProbSpaG(MAF_outlier, residOutlier, S_lower, true,  mean_nonOutlier, var_nonOutlier);

  m_pval = pval1 + pval2;
  return m_pval;
}

double mtSPAmixPlusClass::calculateSparseVariance(const arma::vec& R_new) {
  const double* rp = R_new.memptr();
  double covSum = 0.0;
  for (const auto& ft : m_filteredTriplets)
    covSum += ft.grm * rp[ft.li] * rp[ft.lj];
  return 2.0 * covSum - arma::dot(R_new, R_new);
}
