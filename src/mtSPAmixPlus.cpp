// SPAmixPlus.cpp -- SPAmixPlusClass and namespace-scope SPA helpers

#include <RcppArmadillo.h>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp>

#include "mtSPAmixPlus.h"


namespace SPAmixPlus {


  arma::vec hOrgH2(double t, arma::vec R, const arma::vec MAFVec) {
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
    double H2 = sum(pow(R, 2) % K_G2_vec);
    Horg_H2_vec.at(0) = Horg;
    Horg_H2_vec.at(1) = H2;
    return Horg_H2_vec;
  }

  arma::vec h1AdjH2(double t, arma::vec R, double s, const arma::vec MAFVec) {
    arma::vec H1_adj_H2_vec(2);
    arma::vec tR = t * R;
    arma::vec exp_tR = arma::exp(tR);
    arma::vec MAF_exp_tR = MAFVec % exp_tR;
    arma::vec M_G0_vec = pow((1 - MAFVec + MAF_exp_tR), 2);
    arma::vec M_G1_vec = 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
    arma::vec M_G2_vec = 2 * pow(MAF_exp_tR, 2) + 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
    arma::vec K_G1_vec = M_G1_vec / M_G0_vec;
    arma::vec K_G2_vec = (M_G0_vec % M_G2_vec - pow(M_G1_vec, 2)) / pow(M_G0_vec, 2);
    double H1_adj = sum(R % K_G1_vec) - s;
    double H2 = sum(pow(R, 2) % K_G2_vec);
    H1_adj_H2_vec.at(0) = H1_adj;
    H1_adj_H2_vec.at(1) = H2;
    return H1_adj_H2_vec;
  }

  struct RootResult {
    double root;
    int iter;
    bool converge;
    double K2;
  };

  RootResult fastGetRootK1(double initX,
                            const double& s,
                            const arma::vec MAF_outlier,
                            double mean_nonOutlier,
                            double var_nonOutlier,
                            const arma::vec residOutlier) {
    double x = initX, oldX;
    double K1 = 0, K2 = 0, oldK1;
    double diffX = arma::datum::inf, oldDiffX;
    bool converge = true;
    double tol = 0.001;
    int maxiter = 100;
    int iter = 0;

    for (iter = 0; iter < maxiter; iter ++){

      oldX = x;
      oldDiffX = diffX;
      oldK1 = K1;

      arma::vec H1_adj_H2_vec = h1AdjH2(x, residOutlier, s, MAF_outlier);

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

    RootResult yList = {x, iter, converge, K2};
    return yList;
  }

  double getProbSpaG(const arma::vec MAF_outlier,
                       const arma::vec residOutlier,
                       double s,
                       bool lower_tail,
                       double mean_nonOutlier,
                       double var_nonOutlier) {
    double initX = 0;

    RootResult rootList = fastGetRootK1(initX, s, MAF_outlier, mean_nonOutlier, var_nonOutlier, residOutlier);
    double zeta = rootList.root;

    arma::vec k12 = hOrgH2(zeta, residOutlier, MAF_outlier);
    double k1 = k12.at(0) + mean_nonOutlier * zeta + 0.5 * var_nonOutlier * pow(zeta, 2);
    double k2 = k12.at(1) + var_nonOutlier;

    double temp1 = zeta * s - k1;

    double w = arma::sign(zeta) * sqrt(2 * temp1);
    double v = zeta * sqrt(k2);

    double pval = arma::normcdf(arma::sign(lower_tail - 0.5) * (w + log(v/w) / w));
    return pval;
  }


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


SPAmixPlusClass::SPAmixPlusClass(
  const arma::mat& resid,
  const arma::mat& PCs,
  int N,
  double SPA_Cutoff,
  std::vector<PhenoOutlierData> outlierList,
  std::vector<std::tuple<int, int, double>> sparseTriplets,
  const std::string& afFilePath,
  const std::string& afFilePrecision
) {


  m_afFilePath = afFilePath;
  m_afFilePrecision = afFilePrecision;


  m_sparseTriplets = std::move(sparseTriplets);


  m_resid = resid;
  m_PCs = PCs;
  m_N = N;
  m_SPA_Cutoff = SPA_Cutoff;

  m_Npheno = m_resid.n_cols;

  m_pvalVec.zeros(m_Npheno);
  m_zScoreVec.zeros(m_Npheno);
  m_BetaVec.zeros(m_Npheno);
  m_SVec.zeros(m_Npheno);
  m_SmeanVec.zeros(m_Npheno);
  m_VarSVec.zeros(m_Npheno);

  m_outlierList = std::move(outlierList);


  m_onePlusPCs = arma::join_horiz(arma::ones(N), PCs);
  arma::mat X_t = m_onePlusPCs.t();
  arma::mat XTX = X_t * m_onePlusPCs;
  arma::mat XTX_inv = arma::inv(XTX);
  m_sqrt_XTX_inv_diag = arma::sqrt(XTX_inv.diag());

  m_diffTime1 = arma::vec(1, arma::fill::zeros);
  m_diffTime2 = arma::vec(1, arma::fill::zeros);
}


arma::vec SPAmixPlusClass::fit_lm_get_beta(const arma::vec& g, arma::vec& pvalues) {
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

arma::vec SPAmixPlusClass::fit_lm(const arma::vec& g, arma::vec& pvalues) {
    arma::vec coef = fit_lm_get_beta(g, pvalues);
    return m_onePlusPCs * coef;
}

SPAmixPlusClass::AFModelInfo SPAmixPlusClass::computeAFModel(arma::vec GVec, double altFreq) {
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
        }else{

           arma::mat sigPCs = m_PCs.cols(posSigPCs);


           arma::vec g0(N, arma::fill::zeros);
           arma::uvec posg12 = arma::find(g > 0.5);
           g0.elem(posg12).fill(1);

           double MAC_after = sum(g0);
           if (MAC_after <= MAC_cutoff){
             model.status = 0;
           }else{
             model.status = 2;

             arma::vec sub_beta = logistic_regression_beta(sigPCs, g0);


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

arma::vec SPAmixPlusClass::getAFFromModel(AFModelInfo model, double altFreq) {
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

arma::vec SPAmixPlusClass::getMafEst(arma::vec GVec, double altFreq) {

    AFModelInfo model = computeAFModel(GVec, altFreq);
    return getAFFromModel(model, altFreq);
}


double SPAmixPlusClass::getMarkerPval(arma::vec GVec, double altFreq) {

    arma::vec AFVec = getMafEst(GVec, altFreq);
    m_MAFVec = AFVec;


    arma::vec GVarVec = 2.0 * AFVec % (1.0 - AFVec);

    for (int i = 0; i < m_Npheno; i++){
        const auto& od = m_outlierList[i];

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

        m_zScoreVec.at(i) = zScore;
        m_BetaVec.at(i) = (S - S_mean) / VarS;
        m_SVec.at(i) = S;
        m_SmeanVec.at(i) = S_mean;
        m_VarSVec.at(i) = VarS;

        if (std::abs(zScore) < m_SPA_Cutoff){
            m_pvalVec.at(i) = arma::normcdf(-1.0*std::abs(zScore))*2.0;
            continue;
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

        double pval1 = getProbSpaG(MAF_outlier, residOutlier, S_upper, false, mean_nonOutlier, var_nonOutlier);
        double pval2 = getProbSpaG(MAF_outlier, residOutlier, S_lower, true, mean_nonOutlier, var_nonOutlier);

        m_pvalVec.at(i) = pval1 + pval2;
    }
    return 0.0;
}

double SPAmixPlusClass::getMarkerPvalFromModel(arma::vec GVec, AFModelInfo model, double altFreq) {
    arma::vec AFVec = getAFFromModel(model, altFreq);
    m_MAFVec = AFVec;


    arma::vec GVarVec = 2.0 * AFVec % (1.0 - AFVec);

    for (int i = 0; i < m_Npheno; i++){
        const auto& od = m_outlierList[i];

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


        m_zScoreVec.at(i) = zScore;
        m_BetaVec.at(i) = (S - S_mean) / VarS;
        m_SVec.at(i) = S;
        m_SmeanVec.at(i) = S_mean;
        m_VarSVec.at(i) = VarS;

        if (std::abs(zScore) < m_SPA_Cutoff){

            boost::math::normal_distribution<double> ndist(0.0, 1.0);
            m_pvalVec.at(i) = 2.0 * boost::math::cdf(boost::math::complement(ndist, std::abs(zScore)));
            continue;
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

        double pval1 = getProbSpaG(MAF_outlier, residOutlier, S_upper, false, mean_nonOutlier, var_nonOutlier);
        double pval2 = getProbSpaG(MAF_outlier, residOutlier, S_lower, true, mean_nonOutlier, var_nonOutlier);

        m_pvalVec.at(i) = pval1 + pval2;
    }
    return 0.0;
}

double SPAmixPlusClass::calculateSparseVariance(const arma::vec& R_new, const arma::uvec& posValue) {
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

}
