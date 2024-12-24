
#ifndef SAGELD_HPP
#define SAGELD_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "UTIL.hpp"

namespace SAGELD{

class SAGELDClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  std::string m_Method;                  // method: 'SAGELD' or 'GALLOP'
  arma::mat m_XTs;                       // used to calculate marker-level lambda
  arma::mat m_SS;                        // used to calculate marker-level lambda
  arma::mat m_AtS;                       // used to calculate marker-level lambda
  arma::mat m_Q;                         // used to calculate marker-level lambda
  arma::mat m_A21;                       // used to calculate marker-level lambda
  arma::mat m_TTs;                       // used to calculate marker-level lambda
  arma::mat m_Tys;                       // used to calculate marker-level lambda
  arma::vec m_sol;                       // used to calculate marker-level lambda
  arma::vec m_blups;                     // used to calculate marker-level lambda
  double m_sig;                          // used to calculate marker-level lambda
  double m_ncov;                         // used to calculate marker-level lambda
  
  arma::vec m_resid;                     // residuals
  arma::vec m_resid_G;                   // residuals for G
  arma::vec m_resid_GxE;                 // residuals for GxE
  arma::vec m_resid_E;                 // residuals for E
  
  arma::vec m_resid_unrelated_outliers;     // unrelated outlier residuals
  arma::vec m_resid_unrelated_outliers_G;   // unrelated outlier residuals for G
  arma::vec m_resid_unrelated_outliers_GxE; // unrelated outlier residuals for GxE
  
  double m_sum_unrelated_outliers2;      // sum of squares of unrelated outlier residuals
  double m_sum_unrelated_outliers_G2;    // sum of squares of unrelated outlier residuals for G
  double m_sum_unrelated_outliers_GxE2;  // sum of squares of unrelated outlier residuals for GxE
  double m_sum_unrelated_outliers_G_GxE2;// sum of squares of unrelated outlier residuals for G x GxE
  
  double m_sum_R_nonOutlier;             // sum of non-outlier residuals
  double m_sum_R_nonOutlier_G;           // sum of non-outlier residuals for G
  double m_sum_R_nonOutlier_GxE;         // sum of non-outlier residuals for GxE
  
  double m_R_GRM_R;                      // residuals x GRM x residuals
  double m_R_GRM_R_G;                    // residuals x GRM x residuals for G
  double m_R_GRM_R_GxE;                  // residuals x GRM x residuals for GxE
  double m_R_GRM_R_G_GxE;                // residuals x GRM x residuals for G x GxE
  double m_R_GRM_R_E;                    // residuals x GRM x residuals for E
  
  double m_R_GRM_R_nonOutlier;           // residuals x GRM x residuals for non-outlier families
  double m_R_GRM_R_nonOutlier_G;         // residuals x GRM x residuals for non-outlier families for G
  double m_R_GRM_R_nonOutlier_GxE;       // residuals x GRM x residuals for non-outlier families for GxE
  double m_R_GRM_R_nonOutlier_G_GxE;     // residuals x GRM x residuals for non-outlier families for G x GxE
  
  double m_R_GRM_R_TwoSubjOutlier;       // residuals x GRM x residuals for outlier families (n = 2)
  double m_R_GRM_R_TwoSubjOutlier_G;     // residuals x GRM x residuals for outlier families (n = 2) for G
  double m_R_GRM_R_TwoSubjOutlier_GxE;   // residuals x GRM x residuals for outlier families (n = 2) for GxE
  double m_R_GRM_R_TwoSubjOutlier_G_GxE; // residuals x GRM x residuals for outlier families (n = 2) for G x GxE
  
  Rcpp::List m_TwoSubj_list;             // List of residuals and IBD probabilities in outlier families (n = 2)
  Rcpp::List m_ThreeSubj_list;           // List of residuals and Chow-Liu tree in outlier families (n > 2)
  
  arma::vec m_MAF_interval;              // MAF interval divides the MAFs into several intervals
  double m_zScoreE_cutoff;               // cutoff of standardized score for G-E to use marker-level lambda or mean lambda
  double m_SPA_Cutoff;                   // cutoff of standardized score to use normal approximation or SPA
  double m_zeta;                         // initial saddle point for negative side, default is zero
  double m_tol;                          // accuracy of Newton's methods, default 1e-4
  
  arma::vec m_pvalVec;                   // used to save summary statistics of 'SAGELD' or 'GALLOP'
  arma::vec m_zScoreVec;                 // used to save summary statistics of 'SAGELD' or 'GALLOP'
  arma::vec m_BetaVec;                   // used to save summary statistics of 'SAGELD' or 'GALLOP'
  arma::vec m_seBetaVec;                 // used to save summary statistics of 'SAGELD' or 'GALLOP'
  
public:
  
  SAGELDClass(std::string t_Method,
              arma::mat t_XTs,
              arma::mat t_SS,
              arma::mat t_AtS,
              arma::mat t_Q,
              arma::mat t_A21,
              arma::mat t_TTs,
              arma::mat t_Tys,
              arma::vec t_sol,
              arma::vec t_blups,
              double t_sig,
              arma::vec t_resid,
              arma::vec t_resid_G,
              arma::vec t_resid_GxE,
              arma::vec t_resid_E,
              arma::vec t_resid_unrelated_outliers,
              arma::vec t_resid_unrelated_outliers_G,
              arma::vec t_resid_unrelated_outliers_GxE,
              double t_sum_R_nonOutlier,
              double t_sum_R_nonOutlier_G,
              double t_sum_R_nonOutlier_GxE,
              double t_R_GRM_R,
              double t_R_GRM_R_G,
              double t_R_GRM_R_GxE,
              double t_R_GRM_R_G_GxE,
              double t_R_GRM_R_E,
              double t_R_GRM_R_nonOutlier,
              double t_R_GRM_R_nonOutlier_G,
              double t_R_GRM_R_nonOutlier_GxE,
              double t_R_GRM_R_nonOutlier_G_GxE,
              double t_R_GRM_R_TwoSubjOutlier,
              double t_R_GRM_R_TwoSubjOutlier_G,
              double t_R_GRM_R_TwoSubjOutlier_GxE,
              double t_R_GRM_R_TwoSubjOutlier_G_GxE,
              Rcpp::List t_TwoSubj_list,
              Rcpp::List t_ThreeSubj_list,
              arma::vec t_MAF_interval,
              double t_zScoreE_cutoff,
              double t_SPA_Cutoff,
              double t_zeta,
              double t_tol);
  
  arma::vec getpvalVec(){return m_pvalVec;}
  arma::vec getzScoreVec(){return m_zScoreVec;}
  arma::vec getBetaVec(){return m_BetaVec;}
  arma::vec getseBetaVec(){return m_seBetaVec;}
  
  // The MGF and its first and second derivative MGF of G (genotype)
  arma::mat MGF_cpp(double t, 
                    const Rcpp::List update_ThreeSubj_list,
                    double MAF)
  {
    // Unrelated subjects.
    arma::vec lambda = arma::exp(t * m_resid_unrelated_outliers);
    
    arma::vec alpha = 1 - MAF + MAF * lambda; 
    arma::vec alpha_1 = MAF * m_resid_unrelated_outliers % lambda; 
    arma::vec alpha_2 = m_resid_unrelated_outliers % alpha_1;
    
    arma::vec M_G0_all = alpha % alpha;
    arma::vec M_G1_all = 2 * alpha % alpha_1;
    arma::vec M_G2_all = 2 * (alpha_1 % alpha_1 + alpha % alpha_2);
    
    // Two related subjects in a family.
    int n1 = m_TwoSubj_list.length();
    if (n1 != 0)
    {
      for (int i = 0; i < n1; i++)
      {
        Rcpp::List TwoSubj_list_temp = m_TwoSubj_list[i];
        // arma::vec Resid = Rcpp::as<arma::vec>(TwoSubj_list_temp["Resid"]);
        // arma::vec Rho = Rcpp::as<arma::vec>(TwoSubj_list_temp["Rho"]);
        arma::vec Resid = TwoSubj_list_temp["Resid"];
        arma::vec Rho = TwoSubj_list_temp["Rho"];
        
        arma::vec temp = (1 - Rho) * MAF * (1 - MAF);
        
        double R1 = Resid[0]; double etR1 = exp(t * R1);
        double R2 = Resid[1]; double etR2 = exp(t * R2);
        double Rsum = R1 + R2;
        
        arma::vec midterm1 = etR1 * temp;
        arma::vec midterm2 = etR2 * temp;
        arma::vec midterm3 = etR1 * etR2 * (MAF - temp);
        
        arma::vec M_G0 = midterm1 + midterm2 + midterm3 - temp + 1 - MAF;
        arma::vec M_G1 = R1 * midterm1 + R2 * midterm2 + Rsum * midterm3;
        arma::vec M_G2 = R1*R1 * midterm1 + R2*R2 * midterm2 + Rsum*Rsum * midterm3;
        
        M_G0_all = arma::join_cols(M_G0_all, M_G0);
        M_G1_all = arma::join_cols(M_G1_all, M_G1);
        M_G2_all = arma::join_cols(M_G2_all, M_G2);
      }
    }
    
    // Three above Related Subjects.
    int n2 = update_ThreeSubj_list.length();
    if (n2 != 0)
    {
      for (int i = 0; i < n2; i++)
      {
        Rcpp::List ThreeSubj_list_temp = update_ThreeSubj_list[i];
        // arma::vec stand_S = Rcpp::as<arma::vec>(ThreeSubj_list_temp["stand.S"]);
        // arma::vec arr_prob = Rcpp::as<arma::vec>(ThreeSubj_list_temp["arr.prob"]);
        arma::vec stand_S = ThreeSubj_list_temp["stand.S"];
        arma::vec arr_prob = ThreeSubj_list_temp["arr.prob"];
        
        arma::vec midterm0 = exp(t * stand_S) % arr_prob;
        arma::vec midterm1 = stand_S % midterm0;
        arma::vec midterm2 = stand_S % midterm1;
        
        M_G0_all = arma::join_cols(M_G0_all, arma::vec{arma::accu(midterm0)});
        M_G1_all = arma::join_cols(M_G1_all, arma::vec{arma::accu(midterm1)});
        M_G2_all = arma::join_cols(M_G2_all, arma::vec{arma::accu(midterm2)});
      }
    }
    
    return arma::join_rows(M_G0_all, M_G1_all, M_G2_all);
  }
  
  arma::mat MGF_cpp(double t, 
                    const Rcpp::List update_ThreeSubj_list,
                    double MAF,
                    double lambda_i)
  {
    // Unrelated subjects.
    arma::vec m_resid_unrelated_outliers_i = m_resid_unrelated_outliers_GxE - lambda_i * m_resid_unrelated_outliers_G;
    arma::vec lambda = arma::exp(t * m_resid_unrelated_outliers_i);
    
    arma::vec alpha = 1 - MAF + MAF * lambda; 
    arma::vec alpha_1 = MAF * m_resid_unrelated_outliers_i % lambda; 
    arma::vec alpha_2 = m_resid_unrelated_outliers_i % alpha_1;
    
    arma::vec M_G0_all = alpha % alpha;
    arma::vec M_G1_all = 2 * alpha % alpha_1;
    arma::vec M_G2_all = 2 * (alpha_1 % alpha_1 + alpha % alpha_2);
    
    // Two related subjects in a family.
    int n1 = m_TwoSubj_list.length();
    if (n1 != 0)
    {
      for (int i = 0; i < n1; i++)
      {
        Rcpp::List TwoSubj_list_temp = m_TwoSubj_list[i];
        // arma::vec Resid = Rcpp::as<arma::vec>(TwoSubj_list_temp["Resid"]);
        // arma::vec Rho = Rcpp::as<arma::vec>(TwoSubj_list_temp["Rho"]);
        arma::vec Resid_G = TwoSubj_list_temp["Resid_G"];
        arma::vec Resid_GxE = TwoSubj_list_temp["Resid_GxE"];
        
        arma::vec Resid = Resid_GxE - lambda_i * Resid_G;
        arma::vec Rho = TwoSubj_list_temp["Rho"];
        
        arma::vec temp = (1 - Rho) * MAF * (1 - MAF);
        
        double R1 = Resid[0]; double etR1 = exp(t * R1);
        double R2 = Resid[1]; double etR2 = exp(t * R2);
        double Rsum = R1 + R2;
        
        arma::vec midterm1 = etR1 * temp;
        arma::vec midterm2 = etR2 * temp;
        arma::vec midterm3 = etR1 * etR2 * (MAF - temp);
        
        arma::vec M_G0 = midterm1 + midterm2 + midterm3 - temp + 1 - MAF;
        arma::vec M_G1 = R1 * midterm1 + R2 * midterm2 + Rsum * midterm3;
        arma::vec M_G2 = R1*R1 * midterm1 + R2*R2 * midterm2 + Rsum*Rsum * midterm3;
        
        M_G0_all = arma::join_cols(M_G0_all, M_G0);
        M_G1_all = arma::join_cols(M_G1_all, M_G1);
        M_G2_all = arma::join_cols(M_G2_all, M_G2);
      }
    }
    
    // Three above Related Subjects.
    int n2 = update_ThreeSubj_list.length();
    if (n2 != 0)
    {
      for (int i = 0; i < n2; i++)
      {
        Rcpp::List ThreeSubj_list_temp = update_ThreeSubj_list[i];
        // arma::vec stand_S = Rcpp::as<arma::vec>(ThreeSubj_list_temp["stand.S"]);
        // arma::vec arr_prob = Rcpp::as<arma::vec>(ThreeSubj_list_temp["arr.prob"]);
        arma::vec stand_S = ThreeSubj_list_temp["stand.S"];
        arma::vec arr_prob = ThreeSubj_list_temp["arr.prob"];
        
        arma::vec midterm0 = exp(t * stand_S) % arr_prob;
        arma::vec midterm1 = stand_S % midterm0;
        arma::vec midterm2 = stand_S % midterm1;
        
        M_G0_all = arma::join_cols(M_G0_all, arma::vec{arma::accu(midterm0)});
        M_G1_all = arma::join_cols(M_G1_all, arma::vec{arma::accu(midterm1)});
        M_G2_all = arma::join_cols(M_G2_all, arma::vec{arma::accu(midterm2)});
      }
    }
    
    return arma::join_rows(M_G0_all, M_G1_all, M_G2_all);
  }
  
  // Newton's method to get the saddle point
  double fastgetroot_cpp(const Rcpp::List update_ThreeSubj_list,
                         double Score,
                         double MAF,
                         double init_t,
                         double tol,
                         int maxiter = 50)
  {
    double t = init_t;
    arma::vec MGF0; arma::vec MGF1; arma::vec MGF2;
    double CGF1 = 0; double CGF2 = 0;
    double diff_t = R_PosInf;
    int iter;
    
    double mean = 2 * MAF * m_sum_R_nonOutlier;
    double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier;
    
    for (iter = 1; iter < maxiter; iter++)
    {
      double old_t = t;
      double old_diff_t = diff_t;
      double old_CGF1 = CGF1;
      
      arma::mat MGF_all = MGF_cpp(t, update_ThreeSubj_list, MAF);
      
      MGF0 = MGF_all.col(0);
      MGF1 = MGF_all.col(1);
      MGF2 = MGF_all.col(2);
      
      arma::vec temp = MGF1 / MGF0;
      CGF1 = arma::accu(temp) + mean + var * t - Score;
      CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
      
      diff_t = - CGF1/CGF2;
      
      // std::cout << "iter:\t" << iter << std::endl;
      // std::cout << "t:\t" << t << std::endl;
      // std::cout << "CGF1:\t" << CGF1 << std::endl;
      // std::cout << "CGF2:\t" << CGF2 << std::endl;
      // std::cout << "diff_t:\t" << diff_t << std::endl;
      // std::cout << std::endl;
      
      if (std::isnan(diff_t) || std::isinf(CGF2))
      {
        t = t / 2;
        diff_t = std::min(std::abs(t), 1.0) * arma::sign(Score);
        continue;
      }
      
      if (std::isnan(old_CGF1) || (arma::sign(old_CGF1) != 0 && arma::sign(CGF1) != arma::sign(old_CGF1))) 
      {
        if (std::abs(diff_t) < tol) 
        {
          t = old_t + diff_t;
          break;
        } else {
          while (std::abs(old_diff_t) > tol && std::abs(diff_t) > std::abs(old_diff_t) - tol) 
          {
            diff_t = diff_t / 2;
          }
          t = old_t + diff_t;
          continue;
        }
      }
      
      if (arma::sign(Score) != arma::sign(old_t + diff_t) && 
          (arma::sign(old_CGF1) == 0 || arma::sign(CGF1) == arma::sign(old_CGF1))) 
      {
        while (arma::sign(Score) != arma::sign(old_t + diff_t)) 
        {
          diff_t = diff_t / 2;
        }
        t = old_t + diff_t;
        continue;
      }
      
      t = old_t + diff_t;
      if (std::abs(diff_t) < tol) break;
    }
    
    // return Rcpp::List::create(Named("root") = t,
    //                           Named("iter") = iter);
    
    return t;
  }
  
  double fastgetroot_cpp(const Rcpp::List update_ThreeSubj_list,
                         double m_R_GRM_R_nonOutlier_i,
                         double lambda_i,
                         double Score,
                         double MAF,
                         double init_t,
                         double tol,
                         int maxiter = 50)
  {
    double t = init_t;
    arma::vec MGF0; arma::vec MGF1; arma::vec MGF2;
    double CGF1 = 0; double CGF2 = 0;
    double diff_t = R_PosInf;
    int iter;
    
    double mean = 2 * MAF * (m_sum_R_nonOutlier_GxE - lambda_i * m_sum_R_nonOutlier_G);
    double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier_i;
    
    for (iter = 1; iter < maxiter; iter++)
    {
      double old_t = t;
      double old_diff_t = diff_t;
      double old_CGF1 = CGF1;
      
      arma::mat MGF_all = MGF_cpp(t, update_ThreeSubj_list, MAF, lambda_i);
      
      MGF0 = MGF_all.col(0);
      MGF1 = MGF_all.col(1);
      MGF2 = MGF_all.col(2);
      
      arma::vec temp = MGF1 / MGF0;
      CGF1 = arma::accu(temp) + mean + var * t - Score;
      CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
      
      diff_t = - CGF1/CGF2;
      
      if (std::isnan(diff_t) || std::isinf(CGF2))
      {
        t = t / 2;
        diff_t = std::min(std::abs(t), 1.0) * arma::sign(Score);
        continue;
      }
      
      if (std::isnan(old_CGF1) || (arma::sign(old_CGF1) != 0 && arma::sign(CGF1) != arma::sign(old_CGF1))) 
      {
        if (std::abs(diff_t) < tol) 
        {
          t = old_t + diff_t;
          break;
        } else {
          while (std::abs(old_diff_t) > tol && std::abs(diff_t) > std::abs(old_diff_t) - tol) 
          {
            diff_t = diff_t / 2;
          }
          t = old_t + diff_t;
          continue;
        }
      }
      
      if (arma::sign(Score) != arma::sign(old_t + diff_t) && 
          (arma::sign(old_CGF1) == 0 || arma::sign(CGF1) == arma::sign(old_CGF1))) 
      {
        while (arma::sign(Score) != arma::sign(old_t + diff_t)) 
        {
          diff_t = diff_t / 2;
        }
        t = old_t + diff_t;
        continue;
      }
      
      t = old_t + diff_t;
      if (std::abs(diff_t) < tol) break;
    }
    
    return t;
  }
  
  // function to get one side p value
  double GetProb_SPA(const Rcpp::List update_ThreeSubj_list,
                     double Score,
                     double MAF,
                     bool lower_tail,
                     double zeta,
                     double tol)
  {
    zeta = fastgetroot_cpp(update_ThreeSubj_list, Score, MAF, zeta, tol);
    
    arma::mat MGF_all = MGF_cpp(zeta, update_ThreeSubj_list, MAF);
    
    arma::vec MGF0 = MGF_all.col(0);
    arma::vec MGF1 = MGF_all.col(1);
    arma::vec MGF2 = MGF_all.col(2);
    
    double mean = 2 * MAF * m_sum_R_nonOutlier;
    double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier;
    
    arma::vec temp = MGF1 / MGF0;
    double CGF0 = arma::accu(log(MGF0)) + mean * zeta + 0.5 * var * zeta * zeta;
    double CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
    
    double w = arma::sign(zeta) * sqrt(2 * (zeta * Score - CGF0));
    double v = zeta * sqrt(CGF2);
    
    double u = w + 1/w * log(v/w);
    double pval = R::pnorm(u, 0, 1, lower_tail, false);
    
    // std::cout << "zeta:\t" << zeta << std::endl;
    // std::cout << "p value:\t" << pval << std::endl;
    // std::cout << std::endl;
    
    return pval;
  }
  
  double GetProb_SPA(const Rcpp::List update_ThreeSubj_list,
                     double m_R_GRM_R_nonOutlier_i,
                     double lambda_i,
                     double Score,
                     double MAF,
                     bool lower_tail,
                     double zeta,
                     double tol)
  {
    zeta = fastgetroot_cpp(update_ThreeSubj_list, m_R_GRM_R_nonOutlier_i, lambda_i, Score, MAF, zeta, tol);
    
    arma::mat MGF_all = MGF_cpp(zeta, update_ThreeSubj_list, MAF, lambda_i);
    
    arma::vec MGF0 = MGF_all.col(0);
    arma::vec MGF1 = MGF_all.col(1);
    arma::vec MGF2 = MGF_all.col(2);
    
    double mean = 2 * MAF * (m_sum_R_nonOutlier_GxE - lambda_i * m_sum_R_nonOutlier_G);
    double var = 2 * MAF * (1 - MAF) * m_R_GRM_R_nonOutlier_i;
    
    arma::vec temp = MGF1 / MGF0;
    double CGF0 = arma::accu(log(MGF0)) + mean * zeta + 0.5 * var * zeta * zeta;
    double CGF2 = arma::accu(MGF2 / MGF0) - arma::accu(temp % temp) + var;
    
    double w = arma::sign(zeta) * sqrt(2 * (zeta * Score - CGF0));
    double v = zeta * sqrt(CGF2);
    
    double u = w + 1/w * log(v/w);
    double pval = R::pnorm(u, 0, 1, lower_tail, false);
    
    // std::cout << "zeta:\t" << zeta << std::endl;
    // std::cout << "p value:\t" << pval << std::endl;
    // std::cout << std::endl;
    
    return pval;
  }
  
  // function to get two side p values, 2024-12-23
  double getMarkerPval(arma::vec t_GVec, 
                       double t_altFreq,
                       double& t_hwepval,
                       double t_hwepvalCutoff)
  {
    // updated on 2023-05-23 to get hwe pvalue
    gethwepval(t_GVec, t_hwepval, t_hwepvalCutoff);
    
    double MAF = std::min(t_altFreq, 1 - t_altFreq);
    
    arma::mat t_GVec2, t_GVecq2, m_H1, m_H2, m_AtH, m_R, m_Cfix, m_Cran, m_GtG, m_Gty, m_V, m_v, m_intV;
    
    if(m_Method == "GALLOP")
    {
      t_GVec2 = arma::repmat(t_GVec.t(), 2, 1);
      t_GVec2 = t_GVec2.reshape(t_GVec.n_elem * 2, 1);

      m_H1 = t_GVec.t() * m_XTs; m_H1 = m_H1.reshape(m_ncov, 2);
      m_H2 = m_SS.each_col() % t_GVec2;
      m_AtH = t_GVec.t() * m_AtS; m_AtH = m_AtH.reshape(m_ncov, 2);
      m_R = m_H1 - m_AtH;

      m_Cfix = arma::inv(m_Q) * m_R;
      m_Cran = m_H2 - m_A21 * m_Cfix;
      m_GtG = (t_GVec % t_GVec).t() * m_TTs; m_GtG = m_GtG.reshape(2, 2);
      m_Gty = t_GVec.t() * m_Tys; m_Gty = m_Gty.reshape(2, 1);

      m_V = m_GtG - m_H1.t() * m_Cfix - m_H2.t() * m_Cran;
      m_v = m_Gty - m_H1.t() * m_sol - m_H2.t() * m_blups;
      m_intV = arma::inv(m_V);

      arma::vec Theta = m_intV * m_v;
      arma::vec SD = m_intV.diag();
      arma::vec SE = m_sig * sqrt(SD);
      arma::vec pval = 2 * arma::normcdf(-abs(Theta / SE));

      m_pvalVec.at(0) = pval.at(0); m_pvalVec.at(1) = pval.at(1);
      m_BetaVec.at(0) = Theta.at(0); m_BetaVec.at(1) = Theta.at(1);
      m_seBetaVec.at(0) = SE.at(0); m_seBetaVec.at(1) = SE.at(1);
      
    }else
    {
      double zScore_G, zScore_GxE, pval_G, pval_GxE;

      double G_var = 2 * MAF * (1 - MAF);

      double Score_E = sum(t_GVec % m_resid_E);
      double zScore_E = Score_E/sqrt(G_var * m_R_GRM_R_E);

      zScore_G = sum(t_GVec % m_resid_G)/sqrt(G_var * m_R_GRM_R_G);
      pval_G = 2 * arma::normcdf(-abs(zScore_G));

      if(abs(zScore_E) < m_zScoreE_cutoff)
      {
        double Score = sum(t_GVec % m_resid);
        double Score_var = G_var * m_R_GRM_R;
        zScore_GxE = Score/sqrt(Score_var);

        if(abs(zScore_GxE) <= m_SPA_Cutoff)
        {
          pval_GxE = 2 * arma::normcdf(-abs(zScore_GxE));
        }else
        {
          int order2 = arma::index_max(m_MAF_interval >= MAF);
          int order1 = order2 - 1;

          double MAF_ratio = (m_MAF_interval[order2] - MAF)/(m_MAF_interval[order2] - m_MAF_interval[order1]);

          double Var_ThreeOutlier = 0;

          int n1 = m_ThreeSubj_list.length();
          Rcpp::List update_ThreeSubj_list(n1);

          if (n1 != 0)
          {
            for (int i = 0; i < n1; i++)
            {
              Rcpp::List ThreeSubj_list_temp = m_ThreeSubj_list[i];

              // arma::vec CLT_temp1(243, arma::fill::zeros);
              // arma::vec CLT_temp2(243, arma::fill::zeros);
              // arma::vec stand_S(243, arma::fill::zeros);

              arma::mat CLT_temp =  ThreeSubj_list_temp["CLT"];
              arma::vec stand_S = ThreeSubj_list_temp["stand.S"];
              arma::vec CLT_temp1 = CLT_temp.col(order1);
              arma::vec CLT_temp2 = CLT_temp.col(order2);

              arma::vec arr_prob = MAF_ratio * CLT_temp1 + (1 - MAF_ratio) * CLT_temp2;

              update_ThreeSubj_list[i] = Rcpp::List::create(Rcpp::Named("stand.S") = stand_S,
                                                            Rcpp::Named("arr.prob") = arr_prob);

              arma::vec temp1 = stand_S % arr_prob;

              double temp2 = arma::accu(temp1);
              double Var_ThreeOutlier_temp = arma::accu(stand_S % temp1) - temp2 * temp2;
              Var_ThreeOutlier = Var_ThreeOutlier + Var_ThreeOutlier_temp;
            }
          }

          double Var_nonOutlier = G_var * m_R_GRM_R_nonOutlier;
          double Var_unrelated_outliers = G_var * m_sum_unrelated_outliers2;
          double Var_TwoOutlier = G_var * m_R_GRM_R_TwoSubjOutlier;

          double EmpVar = Var_nonOutlier + Var_unrelated_outliers + Var_TwoOutlier + Var_ThreeOutlier;
          double Var_Ratio = Score_var / EmpVar;
          double Score_adj = Score / sqrt(Var_Ratio);

          double zeta1 = std::abs(Score_adj) / Score_var; zeta1 = std::min(zeta1, 1.2);
          double zeta2 = - std::abs(m_zeta);

          double pval1 = GetProb_SPA(update_ThreeSubj_list, std::abs(Score_adj), MAF, false, zeta1, 1e-4);
          double pval2 = GetProb_SPA(update_ThreeSubj_list, -std::abs(Score_adj), MAF, true, zeta2, m_tol);
          pval_GxE = pval1 + pval2;
        }
      }else
      {
        t_GVec2 = arma::repmat(t_GVec.t(), 2, 1);
        t_GVec2 = t_GVec2.reshape(t_GVec.n_elem * 2, 1);

        m_H1 = t_GVec.t() * m_XTs; m_H1 = m_H1.reshape(m_ncov, 2);
        m_H2 = m_SS.each_col() % t_GVec2;
        m_AtH = t_GVec.t() * m_AtS; m_AtH = m_AtH.reshape(m_ncov, 2);
        m_R = m_H1 - m_AtH;

        m_Cfix = arma::inv(m_Q) * m_R;
        m_Cran = m_H2 - m_A21 * m_Cfix;
        m_GtG = (t_GVec % t_GVec).t() * m_TTs; m_GtG = m_GtG.reshape(2, 2);

        m_V = m_GtG - m_H1.t() * m_Cfix - m_H2.t() * m_Cran;

        double lambda_i = m_V(0,1)/m_V(0,0);

        arma::vec m_resid_i = m_resid_GxE - lambda_i * m_resid_G;
        double m_R_GRM_R_i = m_R_GRM_R_GxE + lambda_i * lambda_i * m_R_GRM_R_G - lambda_i * m_R_GRM_R_G_GxE;

        double Score = sum(t_GVec % m_resid_i);
        double Score_var = G_var * m_R_GRM_R_i;
        zScore_GxE = Score/sqrt(Score_var);

        if(abs(zScore_GxE) <= m_SPA_Cutoff)
        {
          pval_GxE = 2 * arma::normcdf(-abs(zScore_GxE));
        }else
        {
          int order2 = arma::index_max(m_MAF_interval >= MAF);
          int order1 = order2 - 1;

          double MAF_ratio = (m_MAF_interval[order2] - MAF)/(m_MAF_interval[order2] - m_MAF_interval[order1]);

          double Var_ThreeOutlier = 0;

          int n1 = m_ThreeSubj_list.length();
          Rcpp::List update_ThreeSubj_list(n1);

          if (n1 != 0)
          {
            for (int i = 0; i < n1; i++)
            {
              Rcpp::List ThreeSubj_list_temp = m_ThreeSubj_list[i];

              arma::mat CLT_temp =  ThreeSubj_list_temp["CLT"];
              arma::vec stand_S_G = ThreeSubj_list_temp["stand.S_G"];
              arma::vec stand_S_GxE = ThreeSubj_list_temp["stand.S_GxE"];
              arma::vec stand_S = stand_S_GxE - lambda_i * stand_S_G;

              arma::vec CLT_temp1 = CLT_temp.col(order1);
              arma::vec CLT_temp2 = CLT_temp.col(order2);

              arma::vec arr_prob = MAF_ratio * CLT_temp1 + (1 - MAF_ratio) * CLT_temp2;

              update_ThreeSubj_list[i] = Rcpp::List::create(Rcpp::Named("stand.S") = stand_S,
                                                            Rcpp::Named("arr.prob") = arr_prob);

              arma::vec temp1 = stand_S % arr_prob;

              double temp2 = arma::accu(temp1);
              double Var_ThreeOutlier_temp = arma::accu(stand_S % temp1) - temp2 * temp2;
              Var_ThreeOutlier = Var_ThreeOutlier + Var_ThreeOutlier_temp;
            }
          }

          double m_R_GRM_R_nonOutlier_i = m_R_GRM_R_nonOutlier_GxE + lambda_i * lambda_i * m_R_GRM_R_nonOutlier_G - lambda_i * m_R_GRM_R_nonOutlier_G_GxE;
          double Var_nonOutlier = G_var * m_R_GRM_R_nonOutlier_i;

          double m_sum_unrelated_outliers2_i = m_sum_unrelated_outliers_GxE2 + lambda_i * lambda_i * m_sum_unrelated_outliers_G2 - lambda_i * m_sum_unrelated_outliers_G_GxE2;
          double Var_unrelated_outliers = G_var * m_sum_unrelated_outliers2_i;

          double m_R_GRM_R_TwoSubjOutlier_i = m_R_GRM_R_TwoSubjOutlier_GxE + lambda_i * lambda_i * m_R_GRM_R_TwoSubjOutlier_G - lambda_i * m_R_GRM_R_TwoSubjOutlier_G_GxE;
          double Var_TwoOutlier = G_var * m_R_GRM_R_TwoSubjOutlier_i;

          double EmpVar = Var_nonOutlier + Var_unrelated_outliers + Var_TwoOutlier + Var_ThreeOutlier;
          double Var_Ratio = Score_var / EmpVar;
          double Score_adj = Score / sqrt(Var_Ratio);

          double zeta1 = std::abs(Score_adj) / Score_var; zeta1 = std::min(zeta1, 1.2);
          double zeta2 = - std::abs(m_zeta);

          double pval1 = GetProb_SPA(update_ThreeSubj_list, m_R_GRM_R_nonOutlier_i, lambda_i, std::abs(Score_adj), MAF, false, zeta1, 1e-4);
          double pval2 = GetProb_SPA(update_ThreeSubj_list, m_R_GRM_R_nonOutlier_i, lambda_i, -std::abs(Score_adj), MAF, true, zeta2, m_tol);
          pval_GxE = pval1 + pval2;
        }
      }
      m_pvalVec.at(0) = pval_G; m_pvalVec.at(1) = pval_GxE;
      m_zScoreVec.at(0) = zScore_G; m_zScoreVec.at(1) = zScore_GxE;
    }
    
    return 0;
  }
  
};

}

#endif
