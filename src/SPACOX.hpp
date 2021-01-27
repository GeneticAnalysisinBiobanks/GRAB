
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
  
  double GetProb_SPA(double adjG0, 
                     arma::vec adjG1,
                     int N0, 
                     double q2, 
                     bool lowerTail,
                     Rcpp::Function m_K_0_emp,
                     Rcpp::Function m_K_1_emp,
                     Rcpp::Function m_K_2_emp)
  {
    Rcpp::Environment pkg = Rcpp::Environment::namespace_env("stats");
    Rcpp::Function uniroot = pkg["uniroot"];
    arma::vec intervalVec = {-20, 20};
    
    // Rcpp::List rootList = Rcpp::as<Rcpp::List>(uniroot(K_1, 
    //                                                    Rcpp::Named("interval")=intervalVec, 
    //                                                    Rcpp::Named("extendInt")="upX", 
    //                                                    Rcpp::Named("N0")=N0, 
    //                                                    Rcpp::Named("adjG0")=adjG0,
    //                                                    Rcpp::Named("adjG1")=adjG1,
    //                                                    Rcpp::Named("q2")=q2,
    //                                                    Rcpp::Named("m_K_1_emp")=m_K_1_emp));
    // double zeta = rootList["root"];
    double zeta = 0;
    
    double k1 = K_0(zeta,  N0, adjG0, adjG1, m_K_0_emp);
    double k2 = K_2(zeta,  N0, adjG0, adjG1, m_K_2_emp);
    
    double temp1 = zeta * q2 - k1;
    
    double w = arma::sign(zeta) * pow(2 * temp1, 1/2);
    double v = zeta * pow(k2, 1/2);
    
    double pval = arma::normcdf(arma::sign(lowerTail-0.5) * w + 1/w * log(v/w));
    
    return pval;
  }
  
  double getMarkerPval(arma::vec t_GVec, 
                       double t_MAF,
                       Rcpp::Function m_K_0_emp,
                       Rcpp::Function m_K_1_emp,
                       Rcpp::Function m_K_2_emp)
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
        
    double pval1 = GetProb_SPA(adjG0, adjG1, N0, std::abs(z), false, m_K_0_emp, m_K_1_emp, m_K_2_emp);
    double pval2 = GetProb_SPA(adjG0, adjG1, N0, -1*std::abs(z), true, m_K_0_emp, m_K_1_emp, m_K_2_emp);
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

    pval1 = GetProb_SPA(adjG0, adjG1, N0, std::abs(z), false, m_K_0_emp, m_K_1_emp, m_K_2_emp);
    pval2 = GetProb_SPA(adjG0, adjG1, N0, -1*std::abs(z), true, m_K_0_emp, m_K_1_emp, m_K_2_emp);
    pval = pval1 + pval2;
    return pval;
  }
  
};

}

#endif
