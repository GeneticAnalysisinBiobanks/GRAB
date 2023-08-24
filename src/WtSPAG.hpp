
#ifndef WTSPAG_HPP
#define WTSPAG_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

namespace WtSPAG{

class WtSPAGClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  arma::vec m_mresid;
  arma::vec m_weight;
  int m_N;
  double m_SPA_Cutoff;
  arma::vec m_AF_ref;
  arma::vec m_AN_ref;
  arma::vec m_pvalue_bat;
  double m_pvalue_bat_cutoff;
  
public:
  
  WtSPAGClass(arma::mat t_mresid,
              arma::vec t_weight,
              int t_N,
              double t_SPA_Cutoff);
  
  void set_AF_ref(arma::vec t_AF_ref){m_AF_ref = t_AF_ref;}
  void set_AN_ref(arma::vec t_AN_ref){m_AN_ref = t_AN_ref;}
  void set_pvalue_bat(arma::vec t_pvalue_bat){m_pvalue_bat = t_pvalue_bat;}
  void set_pvalue_bat_cutoff(double t_pvalue_bat_cutoff){m_pvalue_bat_cutoff = t_pvalue_bat_cutoff;}
  
  double getMarkerPval(arma::vec t_GVec, 
                       double t_altFreq, 
                       double& t_zScore, 
                       bool t_flip,     // if true, genotype is flipped and AF_ref should be 1-AF_ref
                       int t_i)
  {
    double AF_ref = m_AF_ref.at(t_i);
    double AN_ref = m_AN_ref.at(t_i);
    double pvalue_bat = m_pvalue_bat.at(t_i);
    
    if(t_flip)
      AF_ref = 1 - AF_ref;
    
    // If pvalue for batch effect is less than the cutoff, then we do not use reference (external) allele frequency information
    if(pvalue_bat < m_pvalue_bat_cutoff)
      AN_ref = 0;
    
    // The below is from WtCoxG-2023-07-27-LY.R
    
    double AF = (m_N * t_altFreq + AN_ref * AF_ref) / (m_N + AN_ref);
    double S = sum(m_mresid % (t_GVec - 2 * AF));
    
    double G_var = 2 * AF * (1 - AF);
    
    double meanR = mean(m_mresid);
    double tildeR = meanR * m_N / (m_N + AN_ref);
    
    double S_var = (sum(pow(m_mresid - tildeR, 2)) + AN_ref * pow(tildeR, 2)) * G_var;
    t_zScore = S / sqrt(S_var);
    
    double pval = 0;
    
    // if(std::abs(t_zScore) < m_SPA_Cutoff){
    pval = arma::normcdf(-1*std::abs(t_zScore))*2;
    return pval;
    // }
    
    // Saddlepoint approximation (check SPAmix.hpp for more details)
    
    // return pval;
  }
  
};

}

#endif
