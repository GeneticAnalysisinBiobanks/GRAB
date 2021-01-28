
#ifndef SPACOX_HPP
#define SPACOX_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "approxfun.hpp"

namespace SPACOX{

class SPACOXClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  approxfun::approxfunClass m_K_0_emp;
  approxfun::approxfunClass m_K_1_emp;
  approxfun::approxfunClass m_K_2_emp;

  arma::vec m_mresid;
  double m_varResid;
  arma::mat m_XinvXX, m_tX;
  int m_N;
  double m_pVal_covaAdj_Cutoff;
  double m_pVal_SPA_Cutoff;
  
public:
  
  SPACOXClass(arma::mat t_cumul,
              arma::vec t_mresid,
              arma::mat t_XinvXX,
              arma::mat t_tX,
              int t_N,
              double t_pVal_covaAdj_Cutoff,
              double t_pVal_SPA_Cutoff);
  
  double K_0(double t, 
             int N0, 
             double adjG0, 
             arma::vec adjG1)        // adjusted Genotype 
  {
    double t_adjG0 = t * adjG0;
    arma::vec t_adjG1 = t * adjG1;
    double out = N0 * m_K_0_emp.getValue(t_adjG0) + arma::sum(m_K_0_emp.getVector(t_adjG1));
    return out;
  }
  
  double K_1(double t,
             int N0, 
             double adjG0, 
             arma::vec adjG1,        // adjusted Genotype
             double q2)
  {
    double t_adjG0 = t * adjG0;
    arma::vec t_adjG1 = t * adjG1;
    double out = N0 * adjG0 * m_K_1_emp.getValue(t_adjG0) + arma::sum(m_K_1_emp.getVector(t_adjG1)) - q2;
    return out;
  }
  
  double K_2(double t, 
             int N0, 
             double adjG0, 
             arma::vec adjG1)       // adjusted Genotype
  {
    double t_adjG0 = t * adjG0;
    arma::vec t_adjG1 = t * adjG1;
    double out = N0 * pow(adjG0, 2) * m_K_2_emp.getValue(t_adjG0) + arma::sum(pow(t_adjG1, 2) * m_K_2_emp.getVector(t_adjG1));
    return out;
  }
  
  Rcpp::List fastgetroot_K1(double t_initX,
                            int N0, 
                            double adjG0, 
                            arma::vec adjG1,        // adjusted Genotype
                            double q2)
  {
    double x = t_initX;
    double K1 = 0;
    double K2 = 0;
    double diffX = arma::datum::inf;
    bool converge = true;
    double tol = 0.001;
    int maxiter = 100;
    int iter = 0;
    
    for(iter = 0; iter < maxiter; iter ++){
      double oldX = x;
      double oldDiffX = diffX;
      double oldK1 = K1;
      
      double K1 = K_1(x, N0, adjG0, adjG1, q2);
      double K2 = K_2(x, N0, adjG0, adjG1);
      
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
    
    if(iter == maxiter - 1) 
      converge = false;
    
    Rcpp::List yList = Rcpp::List::create(Rcpp::Named("root") = x,
                                          Rcpp::Named("iter") = iter,
                                          Rcpp::Named("converge") = converge,
                                          Rcpp::Named("K2") = K2);
    return yList;
  }
  
  double GetProb_SPA(double adjG0, 
                     arma::vec adjG1,
                     int N0, 
                     double q2, 
                     bool lowerTail)
  {
    double initX;
    if(q2 > 0) initX = 3;
    if(q2 <= 0) initX = -3;
    
    Rcpp::List rootList = fastgetroot_K1(initX, N0, adjG0, adjG1, q2);
    double zeta = rootList["root"];
    
    double k1 = K_0(zeta,  N0, adjG0, adjG1);
    double k2 = K_2(zeta,  N0, adjG0, adjG1);
    
    double temp1 = zeta * q2 - k1;
    
    double w = arma::sign(zeta) * pow(2 * temp1, 1/2);
    double v = zeta * pow(k2, 1/2);
    
    double pval = arma::normcdf(arma::sign(lowerTail-0.5) * w + 1/w * log(v/w));
    
    return pval;
  }
  
  double getMarkerPval(arma::vec t_GVec, 
                       double t_MAF)
  {
    double S = sum(t_GVec * m_mresid);
    arma::vec adjGVec = t_GVec - 2 * t_MAF;
    arma::vec adjGVec2 = pow(adjGVec, 2);
    double VarS = m_varResid * sum(adjGVec2);
    double z = S / sqrt(VarS);
    
    if(std::abs(z) < m_pVal_SPA_Cutoff){
      double pval = arma::normcdf(-1*std::abs(z))*2;
      return pval;
    }
      
    arma::uvec N1set = arma::find(t_GVec!=0);  // position of non-zero genotypes
    int N0 = m_N - N1set.size();
        
    arma::vec adjGVecNorm = adjGVec / sqrt(VarS); // normalized genotype (such that sd=1)
    
    arma::vec adjG1 = adjGVecNorm.elem(N1set);
    double adjG0 = -2 * t_MAF / sqrt(VarS);  // all subjects with g=0 share the same normlized genotype, this is to reduce computation time
        
    double pval1 = GetProb_SPA(adjG0, adjG1, N0, std::abs(z), false);
    double pval2 = GetProb_SPA(adjG0, adjG1, N0, -1*std::abs(z), true);
    double pval = pval1 + pval2;
    
    if(pval > m_pVal_covaAdj_Cutoff){
      return pval;
    }
          
    // estimated variance after adjusting for covariates
          
    adjGVec = t_GVec - m_XinvXX * m_tX.cols(N1set) * t_GVec.elem(N1set);
    adjGVec2 = pow(adjGVec, 2);
    VarS = m_varResid * sum(adjGVec2);
    z = S / sqrt(VarS);
            
    adjGVecNorm = adjGVec / sqrt(VarS);
            
    N0 = 0;
    adjG1 = adjGVecNorm;
    adjG0 = 0;   // since N0=0, this value actually does not matter

    pval1 = GetProb_SPA(adjG0, adjG1, N0, std::abs(z), false);
    pval2 = GetProb_SPA(adjG0, adjG1, N0, -1*std::abs(z), true);
    pval = pval1 + pval2;
    return pval;
  }
  
};

}

#endif
