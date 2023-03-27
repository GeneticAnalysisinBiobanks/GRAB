
#ifndef SPAMIX_HPP
#define SPAMIX_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace SPAmix{

class SPAmixClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  arma::vec m_resid;  // residuals
  arma::vec m_resid2; // residuals^2
  arma::mat m_XinvXX, m_tX;
  int m_N;             // sample size
  double m_SPA_Cutoff;
  arma::mat m_PCs;     // SNP-derived principal components (PCs)
  
public:
  
  SPAmixClass(arma::vec t_resid,
              arma::mat t_XinvXX,
              arma::mat t_tX,
              arma::mat t_PCs,
              int t_N,
              double t_SPA_Cutoff);
  
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
  
  // The below code is from SPACox.hpp
  Rcpp::List fastgetroot_K1(double t_initX,
                            arma::vec R,
                            const double& s,
                            const arma::vec& MAFVec)
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
      
      K1 = H1_adj(x, R, s, MAFVec);
      K2 = H2(x, R, MAFVec);
      
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
  
  double GetProb_SPA_G(arma::vec MAFVec, arma::vec R, double s, bool lower_tail)
  {
    double initX = 0;
    
    // The following initial values are validated on 03/25/2021
    // if(q2 > 0) initX = 3;
    // if(q2 <= 0) initX = -3;
    
    Rcpp::List rootList = fastgetroot_K1(initX, R, s, MAFVec);
    double zeta = rootList["root"];
    
    std::cout << "zeta:\t" << zeta << std::endl;
    
    double k1 = H_org(zeta, R, MAFVec);
    double k2 = H2(zeta, R, MAFVec);
    
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
  
  arma::vec fit_lm(const arma::mat& PCs, const arma::vec& g, arma::vec& pvalues) 
  {
    int n = PCs.n_rows;
    int k = PCs.n_cols;
    
    arma::mat X = arma::join_horiz(arma::ones(n), PCs);
    arma::vec coef = arma::solve(X, g);
    
    arma::vec fittedValues = X * coef;
    
    double s2 = sum(square(g - fittedValues)) / (n - k - 1);
    
    arma::mat X_t = X.t();
    arma::mat XTX = X_t * X;
    arma::mat XTX_inv = arma::inv(XTX);
    arma::vec XTX_inv_diag = XTX_inv.diag();
    
    arma::vec se = arma::sqrt(XTX_inv_diag * s2);
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
  
  arma::vec getMAFest(arma::mat PCs,
                      arma::vec g,
                      double MAC_cutoff = 20,
                      double PCs_pvalue_cutoff = 0.05,
                      double MAF_est_negative_ratio_cutoff = 0.1) 
  {
    int N = g.n_elem;
    arma::vec g0(N, arma::fill::zeros);  // 0 if g < 0.5, 1 if g > 0.5.
    arma::vec MAF_all = arma::vec(N, arma::fill::value(arma::mean(g)/2.0));  // check NA later
    double MAC = arma::accu(g);
    int PC_number = PCs.n_cols;
    double MAF0 = 0;
    
    arma::vec pvalues(PC_number);
    arma::vec MAF_est, topPCs_pvalueVec;
    
    std::cout << "MAC:\t" << MAC << std::endl;
    std::cout << "MAC_cutoff:\t" << MAC_cutoff << std::endl;
    
    if(MAC <= MAC_cutoff){
      MAF_est = MAF_all;
    }else{
      arma::vec fit = fit_lm(PCs, g, pvalues);
      
      arma::uvec posZero = arma::find(fit < 0);
      arma::uvec posOne = arma::find(fit > 1);
      int nError = posZero.n_elem + posOne.n_elem; 
      double propError = nError / N;
      
      std::cout << "propError:\t" << propError << std::endl;
      std::cout << "MAF_est_negative_ratio_cutoff:\t" << MAF_est_negative_ratio_cutoff << std::endl;
      
      if(propError < MAF_est_negative_ratio_cutoff){
        fit.elem(posZero).fill(MAF0);
        fit.elem(posOne).fill(1-MAF0);
        MAF_est = fit;
      }else{
        arma::uvec posSigPCs = arma::find(pvalues < PCs_pvalue_cutoff);
        
        std::cout << "posSigPCs:\t" << posSigPCs << std::endl;
        
        if(posSigPCs.n_elem == 0){
          MAF_est = MAF_all;
        }else{
          arma::mat sigPCs = PCs.cols(posSigPCs);
          arma::uvec posg12 = arma::find(g > 0.5);
          g0.elem(posg12).fill(1);
          
          MAF_est = logistic_regression(sigPCs, g0);
          // Fit.logistic = glm(cbind(round(g), 2-round(g))~selected.PCs, family = binomial)
          // MAF.est = Fit.logistic$fitted.values
        }
      }
    }
    
    return MAF_est;
  }
  
  double getMarkerPval(arma::vec t_GVec)
  {
    arma::vec AFVec = getMAFest(m_PCs, t_GVec);
    
    arma::vec GVarVec = 2 * AFVec % (1 - AFVec);
    double S = sum(t_GVec % m_resid);
    double VarS = sum(m_resid2 % GVarVec);
    double z = S / sqrt(VarS);
    
    if(std::abs(z) < m_SPA_Cutoff){
      double pval = arma::normcdf(-1*std::abs(z))*2;
      return pval;
    }
    
    double pval1 = GetProb_SPA_G(AFVec, m_resid, std::abs(z), false);
    double pval2 = GetProb_SPA_G(AFVec, m_resid, -1*std::abs(z), true);
    double pval = pval1 + pval2;
    
    return pval;
  }
  
};

}

#endif
