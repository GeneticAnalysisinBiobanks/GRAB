
#ifndef SPAMIX_HPP
#define SPAMIX_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "UTIL.hpp"

namespace SPAmix{

class SPAmixClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  arma::mat m_resid;  // residuals
  // arma::mat m_resid2; // residuals^2
  // arma::mat m_XinvXX, m_tX;
  arma::mat m_onePlusPCs;
  
  int m_N;             // sample size
  int m_Npheno;        // number of phenotypes (>= 1)
  
  double m_SPA_Cutoff;
  arma::mat m_PCs;     // SNP-derived principal components (PCs)
  arma::vec m_sqrt_XTX_inv_diag; // derived from PCs, check SPAmix.cpp for more details
  arma::vec m_diffTime1, m_diffTime2;
  
  Rcpp::List m_outlierList;
  
  arma::vec m_pvalVec;
  arma::vec m_zScoreVec;
  // arma::uvec m_posOutlier;
  // arma::uvec m_posNonOutlier;
  // 
  // arma::vec m_R_outlier;
  // arma::vec m_R_nonOutlier;
  
public:
  
  SPAmixClass(arma::mat t_resid,
              arma::mat t_PCs,
              int t_N,
              double t_SPA_Cutoff,
              Rcpp::List t_outlierList);
  
  arma::vec getTestTime1(){return m_diffTime1;}
  arma::vec getTestTime2(){return m_diffTime2;}
  int getNpheno(){return m_Npheno;}
  arma::vec getpvalVec(){return m_pvalVec;}
  arma::vec getzScoreVec(){return m_zScoreVec;}
  // a major update on 2023-04-23 to improve SPA process
  
  
  // The MGF of G (genotype)
  arma::vec M_G0(arma::vec t, arma::vec MAF){
    arma::vec re = pow((1 - MAF + MAF % arma::exp(t)), 2);
    return re;
  }
  
  // The first derivative of the MGF of G (genotype)
  arma::vec M_G1(arma::vec t, arma::vec MAF){
    arma::vec re = 2 * (MAF % arma::exp(t)) % (1 - MAF + MAF % arma::exp(t));
    return re;                           
  }
  
  // The second derivative of the MGF of G (genotype)
  arma::vec M_G2(arma::vec t, arma::vec MAF){
    arma::vec re = 2 * pow(MAF % arma::exp(t), 2) + 2 * (MAF % arma::exp(t)) % (1 - MAF + MAF % arma::exp(t));
    return re;
  }
  
  // The CGF of G (genotype)
  arma::vec K_G0(arma::vec t, arma::vec MAF){
    arma::vec re = arma::log(M_G0(t, MAF));
    return re;
  }
  
  // The first derivative of the CGF of G (genotype)
  arma::vec K_G1(arma::vec t, arma::vec MAF){
    arma::vec re = M_G1(t, MAF) / M_G0(t, MAF);
    return re;
  }
  
  // The second derivative of the CGF of G (genotype)
  arma::vec K_G2(arma::vec t, arma::vec MAF){
    arma::vec re = (M_G0(t, MAF) % M_G2(t, MAF) - pow(M_G1(t, MAF), 2)) / pow(M_G0(t, MAF), 2);
    return re;
  }
  
  // The CGF of score test statistic 
  double H_org(double t, arma::vec R, const arma::vec& MAFVec) {
    double out = sum(K_G0(t * R, MAFVec));
    return out;
  }
  
  // The first derivative of the CGF of score test statistic
  double H1_adj(double t, arma::vec R, const double& s, const arma::vec& MAFVec) {
    double out = sum(R % K_G1(t * R, MAFVec)) - s;
    return out;
  }
  
  // The second derivative of the CGF of score test statistic 
  double H2(double t, arma::vec R, const arma::vec& MAFVec) {
    double out = sum(pow(R, 2) % K_G2(t * R, MAFVec));
    return out;
  }
  
  // partial normal distribution approximation
  
  arma::vec Horg_H2(double t, arma::vec R, const arma::vec MAFVec)
  {
    arma::vec Horg_H2_vec(2);
    arma::vec t_R = t * R;
    arma::vec exp_tR = arma::exp(t_R);
    arma::vec MAF_exp_tR = MAFVec % exp_tR;
    arma::vec M_G0_vec = pow((1 - MAFVec + MAF_exp_tR), 2);;
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
  
  arma::vec H1_adj_H2(double t, arma::vec R, double s, const arma::vec MAFVec)
  {
    arma::vec H1_adj_H2_vec(2);
    arma::vec t_R = t * R;
    arma::vec exp_tR = arma::exp(t_R);
    arma::vec MAF_exp_tR = MAFVec % exp_tR;
    arma::vec M_G0_vec = pow((1 - MAFVec + MAF_exp_tR), 2);;
    arma::vec M_G1_vec = 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
    arma::vec M_G2_vec = 2 * pow(MAF_exp_tR, 2) + 2 * (MAF_exp_tR) % (1 - MAFVec + MAF_exp_tR);
    arma::vec K_G1_vec = M_G1_vec / M_G0_vec;
    arma::vec K_G2_vec = (M_G0_vec % M_G2_vec - pow(M_G1_vec, 2)) / pow(M_G0_vec, 2);
    double H1_adj = sum(R % K_G1_vec) - s;
    double H2 = sum(pow(R, 2) % K_G2_vec);
    // double H1_adj = sum(R % K_G1(t * R, MAFVec)) - s;
    // double H2 = sum(pow(R, 2) % K_G2(t * R, MAFVec));
    H1_adj_H2_vec.at(0) = H1_adj;
    H1_adj_H2_vec.at(1) = H2;
    return H1_adj_H2_vec;
  }
  
  // END of the major update on 2023-04-23 to improve SPA process
  
  // The below code is from SPACox.hpp
  Rcpp::List fastgetroot_K1(double t_initX,
                            const double& s,
                            const arma::vec MAF_outlier,
                            double mean_nonOutlier,
                            double var_nonOutlier,
                            const arma::vec residOutlier)
  {
    double x = t_initX, oldX;
    double K1 = 0, K2 = 0, oldK1;
    double diffX = arma::datum::inf, oldDiffX;
    bool converge = true;
    double tol = 0.001;
    int maxiter = 100;
    int iter = 0;
    
    for(iter = 0; iter < maxiter; iter ++){
      
      oldX = x;
      oldDiffX = diffX;
      oldK1 = K1;
      
      // K1 = H1_adj(x, R, s, MAFVec);
      // K2 = H2(x, R, MAFVec);
      
      arma::vec H1_adj_H2_vec = H1_adj_H2(x, residOutlier, s, MAF_outlier);
      
      K1 = H1_adj_H2_vec.at(0) + mean_nonOutlier + var_nonOutlier * x;
      K2 = H1_adj_H2_vec.at(1) + var_nonOutlier;
            
      diffX = -1 * K1 / K2;
      
      if(!std::isfinite(K1)){
        // checked it on 07/05:
        // if the solution 'x' tends to infinity, 'K2' tends to 0, and 'K1' tends to 0 very slowly.
        // then we can set the one sided p value as 0 (instead of setting converge = F)
        x = arma::datum::inf;
        K2 = 0;
        break;
      }
      
      if(arma::sign(K1) != arma::sign(oldK1)){
        while(std::abs(diffX) > std::abs(oldDiffX) - tol){
          diffX = diffX / 2;
        }
      }
      
      if(std::abs(diffX) < tol) break;
      
      x = oldX + diffX;
    }
    
    if(iter == maxiter) 
      converge = false;
    
    Rcpp::List yList = Rcpp::List::create(Rcpp::Named("root") = x,
                                          Rcpp::Named("iter") = iter,
                                          Rcpp::Named("converge") = converge,
                                          Rcpp::Named("K2") = K2);
    return yList;
  }
  
  double GetProb_SPA_G(const arma::vec MAF_outlier, 
                       const arma::vec residOutlier, 
                       double s, 
                       bool lower_tail,
                       double mean_nonOutlier,
                       double var_nonOutlier)
  {
    double initX = 0;
    
    // The following initial values are validated on 03/25/2021
    // if(q2 > 0) initX = 3;
    // if(q2 <= 0) initX = -3;
    
    Rcpp::List rootList = fastgetroot_K1(initX, s, MAF_outlier, mean_nonOutlier, var_nonOutlier, residOutlier);
    double zeta = rootList["root"];
    
    // std::cout << "zeta:\t" << zeta << std::endl;
    
    // double k1 = H_org(zeta, R, MAFVec);
    // double k2 = H2(zeta, R, MAFVec);
    arma::vec k12 = Horg_H2(zeta, residOutlier, MAF_outlier);
    double k1 = k12.at(0) + mean_nonOutlier * zeta + 0.5 * var_nonOutlier * pow(zeta, 2);
    double k2 = k12.at(1) + var_nonOutlier;
    
    double temp1 = zeta * s - k1;
    
    double w = arma::sign(zeta) * sqrt(2 * temp1);
    double v = zeta * sqrt(k2);
    
    double pval = arma::normcdf(arma::sign(lower_tail - 0.5) * (w + log(v/w) / w));
    return pval;
  }
  
  arma::vec simulate_uniform(int n, double lower, double upper) {
    arma::vec vec(n);
    vec.randu();
    vec = lower + (upper - lower) * vec;
    return vec;
  }
  
  // NOTE about the udpate (2023-04-20): PCs (and some objects calculated based on PCs) are the same for all genetic variants
  // arma::vec fit_lm(const arma::mat& PCs, const arma::vec& g, arma::vec& pvalues) 
  arma::vec fit_lm(const arma::vec& g, arma::vec& pvalues) 
  {
    // int n = PCs.n_rows;
    int n = m_N;
    int k = m_PCs.n_cols;
    
    // NOTE (2023-04-20): the below part has been moved to SPAmix.cpp since they are shared for all genetic variants
    // arma::mat X = arma::join_horiz(arma::ones(n), PCs);
    // arma::mat X_t = X.t();
    // arma::mat XTX = X_t * X;
    // arma::mat XTX_inv = arma::inv(XTX);
    // arma::vec XTX_inv_diag = XTX_inv.diag();
    
    // arma::vec coef = arma::solve(X, g);  // can be improved: 2023-04-20
    // arma::vec fittedValues = X * coef;
    
    arma::vec coef = arma::solve(m_onePlusPCs, g);  // can be improved: 2023-04-20: checked by BWJ, SVD idea and QR idea does not work
    arma::vec fittedValues = m_onePlusPCs * coef;
    
    double s2 = sum(square(g - fittedValues)) / (n - k - 1);
    
    // arma::vec se = arma::sqrt(m_XTX_inv_diag * s2);
    arma::vec se = m_sqrt_XTX_inv_diag * sqrt(s2);
    arma::vec t = coef / se;
    
    for(int i = 0; i < k; i++){
      pvalues[i] = 2 * R::pt(abs(t[i+1]), n-k-1, 0, 0);
    }
    return fittedValues;
  }
  
  arma::vec logistic_regression(const arma::mat& X, const arma::vec& y) {
    
    int n = X.n_rows;
    int p = X.n_cols;
    
    // arma::mat X_new(n, p + 1);
    arma::mat WX_new(n, p + 1);
    // X_new.col(0) = arma::ones<vec>(n);
    // X_new.cols(1, p) = X;
    arma::mat X_new = arma::join_horiz(arma::ones(n), X);
    
    arma::vec beta(p+1, arma::fill::zeros);
    double tol = 1e-6;
    int max_iter = 100;
    arma::vec mu(n);
    
    for (int i = 0; i < max_iter; i++) {
      mu = 1 / (1 + exp(-X_new * beta));
      
      // arma::mat W = diagmat(mu % (1 - mu));
      // arma::vec z = X_new * beta + inv(W) * (y - mu);
      // arma::vec beta_new = inv(X_new.t() * W * X_new) * X_new.t() * W * z;
      
      arma::vec W = mu % (1 - mu);
      arma::vec z = X_new * beta + (y - mu) / W;
      
      for(int j = 0; j < p+1; j++){
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
  
  arma::vec getMAFest(arma::vec g,
                      double t_altFreq,
                      double MAC_cutoff = 20,
                      double PCs_pvalue_cutoff = 0.05,
                      double MAF_est_negative_ratio_cutoff = 0.1) 
  {
    int N = g.n_elem;
    arma::vec g0(N, arma::fill::zeros);  // 0 if g < 0.5, 1 if g > 0.5.
    // arma::vec MAF_all = arma::vec(N, arma::fill::value(arma::mean(g)/2.0));  // check NA later
    // double MAC = arma::accu(g);
    arma::vec MAF_all = arma::vec(N, arma::fill::value(t_altFreq));  // check NA later
    double MAC = t_altFreq * 2 * N; // what if the alternative allele is not minor allele? 2023-04-20
    
    int PC_number = m_PCs.n_cols;
    double MAF0 = 0;
    
    arma::vec pvalues(PC_number);
    arma::vec MAF_est, topPCs_pvalueVec;
    
    // std::cout << "MAC:\t" << MAC << std::endl;
    // std::cout << "MAC_cutoff:\t" << MAC_cutoff << std::endl;
    
    if(MAC <= MAC_cutoff){
      MAF_est = MAF_all;
    }else{
      // arma::vec fit = fit_lm(m_PCs, g, pvalues);
      arma::vec fit = fit_lm(g, pvalues);  // PCs (and some objects calculated based on PCs) are the same for all genetic variants
      fit = fit / 2;
      
      arma::uvec posZero = arma::find(fit < 0);
      arma::uvec posOne = arma::find(fit > 1);
      
      // std::cout << "posZero:\t" << posZero.size() << std::endl;
      // std::cout << "posOne:\t" << posOne.size() << std::endl;
      
      int nError = posZero.n_elem + posOne.n_elem; 
      double propError = (double)nError / N;
      
      // std::cout << "propError:\t" << propError << std::endl;
      // std::cout << "MAF_est_negative_ratio_cutoff:\t" << MAF_est_negative_ratio_cutoff << std::endl;
      
      if(propError < MAF_est_negative_ratio_cutoff){
        fit.elem(posZero).fill(MAF0);
        fit.elem(posOne).fill(1-MAF0);
        MAF_est = fit;  // updated on 2023-04-23
      }else{
        arma::uvec posSigPCs = arma::find(pvalues < PCs_pvalue_cutoff);

        if(posSigPCs.n_elem == 0){
          MAF_est = MAF_all;
        }else{
          arma::mat sigPCs = m_PCs.cols(posSigPCs);
          arma::uvec posg12 = arma::find(g > 0.5);
          g0.elem(posg12).fill(1);
          
          // updated on 2023-04-23: caution! should be checked later
          double MAC_after = sum(g0);
          if(MAC_after <= MAC_cutoff){
            MAF_est = MAF_all;  // end of the update on 2023-04-23
          }else{
            MAF_est = logistic_regression(sigPCs, g0);
          }
        }
      }
    }
    
    return MAF_est;
  }
  
  double getMarkerPval(arma::vec t_GVec, 
                       double t_altFreq)  // later update score and variance here (2023-06-20)
  {
    // std::cout << "part1" << std::endl;
    
    // estimate allele frequency based on PC information and raw genotpe
    // arma::vec AFVec = getMAFest(m_PCs, t_GVec, t_altFreq);
    arma::vec time1 = getTime();
    
    arma::vec AFVec = getMAFest(t_GVec, t_altFreq);  // PCs are global variables and can be loaded when necessary
    
    arma::vec time2 = getTime();
    arma::vec diffTime = time2 - time1;
    // std::cout << "part2" << std::endl;
    
    // std::cout << "(MAF) diffTime:\t" << diffTime << std::endl;
    
    m_diffTime2 += diffTime;
    
    arma::vec GVarVec = 2 * AFVec % (1 - AFVec);
    
    // (BWJ) 2023-06-20: Support multiple phenotypes
    // outLierList[[i]] = list(posValue = posValue - 1,
    //                         posOutlier = posOutlier - 1,
    //                         posNonOutlier = posNonOutlier - 1,
    //                         resid = mresid.temp[posValue],
    //                         resid2 = mresid.temp[posValue]^2,
    //                         residOutleir = mresid.temp[posOutlier],
    //                         residNonOutlier = mresid.temp[posNonOutlier])
    for(int i = 0; i < m_Npheno; i++){
      Rcpp::List tempOutlierList = m_outlierList[i];
      arma::uvec posValue = tempOutlierList["posValue"];
      arma::uvec posOutlier = tempOutlierList["posOutlier"];
      arma::uvec posNonOutlier = tempOutlierList["posNonOutlier"];
      arma::vec resid = tempOutlierList["resid"];
      arma::vec resid2 = tempOutlierList["resid2"];
      arma::vec residOutlier = tempOutlierList["residOutlier"];
      arma::vec residNonOutlier = tempOutlierList["residNonOutlier"];
      arma::vec resid2NonOutlier = tempOutlierList["resid2NonOutlier"];
      
      double S = sum(t_GVec.elem(posValue) % resid);
      double VarS = sum(resid2 % GVarVec.elem(posValue));
      
      // updated on 2023-04-23
      double S_mean = 2 * sum(resid % AFVec.elem(posValue)); // NOTE: I think S_mean is somewhat weird and should be checked later (2023-04-22)
      double zScore = (S-S_mean) / sqrt(VarS);
      
      m_zScoreVec.at(i) = zScore;
      // std::cout << "part3" << std::endl;
      // std::cout << "S:\t" << S << std::endl;
      // std::cout << "VarS:\t" << VarS << std::endl;
      // std::cout << "t_zScore:\t" << t_zScore << std::endl;
      
      if(std::abs(zScore) < m_SPA_Cutoff){
        m_pvalVec.at(i) = arma::normcdf(-1*std::abs(zScore))*2;
        continue;
        // return pval;
      }
      
      // std::cout << "part4" << std::endl;
      
      // double pval1 = GetProb_SPA_G(AFVec, m_resid, std::abs(t_zScore), false);
      // double pval2 = GetProb_SPA_G(AFVec, m_resid, -1*std::abs(t_zScore), true);
      
      time1 = getTime();
      
      // put the below objects as a global object
      arma::vec MAF_outlier = AFVec.elem(posOutlier);
      arma::vec MAF_nonOutlier = AFVec.elem(posNonOutlier);
      
      double mean_nonOutlier = sum(residNonOutlier % MAF_nonOutlier) * 2;
      double var_nonOutlier = sum(resid2NonOutlier % MAF_nonOutlier % (1-MAF_nonOutlier)) * 2;
      
      double pval1 = GetProb_SPA_G(MAF_outlier, 
                                   residOutlier, 
                                   std::abs(S-S_mean)+S_mean, 
                                   false,
                                   mean_nonOutlier,
                                   var_nonOutlier);
      double pval2 = GetProb_SPA_G(MAF_outlier, 
                                   residOutlier, 
                                   -1*std::abs(S-S_mean)+S_mean, 
                                   true,
                                   mean_nonOutlier,
                                   var_nonOutlier);
      
      time2 = getTime();
      diffTime = time2 - time1;
      
      // std::cout << "(SPA_G) diffTime:\t" << diffTime << std::endl;
      
      m_diffTime1 += diffTime;

      m_pvalVec.at(i) = pval1 + pval2;
    }
    
    // std::cout << "part5" << std::endl;
    double pval = 0; // we modify the codes to save pval to m_pvalVec and thus the function does not output pvalue any more.
    return pval;
  }
  
};

}

#endif
