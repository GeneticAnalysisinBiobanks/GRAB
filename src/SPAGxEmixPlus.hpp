
#ifndef SPAGxEMIXPLUS_HPP
#define SPAGxEMIXPLUS_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>  // 必须包含以使用Rcpp环境
#include <RcppArmadillo.h>
#include "UTIL.hpp"
#include <unordered_map>   // 新增
#include <vector>          // 新增
#include <cmath>  // 需要包含此头文件以使用 std::isnan
#include <Rcpp/Formula.h>  // 必须包含Formula支持
#include <set>



using namespace Rcpp;           // 声明命名空间

// #define ARMA_DONT_USE_WRAPPER // 提升Armadillo性能
// #include <RcppArmadillo.h>
// using namespace arma;


namespace SPAGxEmixPlus{

class SPAGxEmixPlusClass
{
private:
  
  std::string m_ResidTraitType; // 表征类型成员变量
  arma::mat m_PhenoMat;
  arma::mat m_Covariates;
  
  
  // 修改后（正确）
  Environment m_statsEnv;
  Function m_glm;
  
  // Environment m_statsEnv = Environment::namespace_env("stats");
  // Function m_glm = m_statsEnv["glm"];
  
  // ==== 新增环境因子成员 ====
  arma::vec m_E;  // 必须与.cpp中的使用保持一致
  
  // ==== 修改1：新增成员 ====
  // std::unordered_map<std::string, int> m_idMap;  // ID映射表
  // arma::mat m_ResidMat;                          // 多表型残差矩阵（N×M）
  // std::vector<std::tuple<int, int, double>> m_sparseTriplets; // 稀疏矩阵三元组
  arma::vec m_MAFVec;                            // MAF估计值
  
  // ==== 新增成员 ====
  arma::ivec m_subjIndices;  // 存储SubjID_Index（整数索引）
  arma::mat m_ResidMat;      // 存储Resid_*列（双精度矩阵）
  std::vector<std::tuple<int, int, double>> m_sparseTriplets; // 稀疏矩阵三元组
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  arma::mat m_resid;  // residuals
  
  arma::mat m_resid_by_E;  // update by Yuzhuo Ma
  
  
  
  // arma::mat m_resid2; // residuals^2
  // arma::mat m_XinvXX, m_tX;
  arma::mat m_onePlusPCs;
  
  int m_N;             // sample size
  int m_Npheno;        // number of phenotypes (>= 1)
  
  double m_SPA_Cutoff;
  arma::mat m_PCs;     // SNP-derived principal components (PCs)
  arma::vec m_sqrt_XTX_inv_diag; // derived from PCs, check SPAGxEmixPlus.cpp for more details
  arma::vec m_diffTime1, m_diffTime2;
  
  Rcpp::List m_outlierList;
  
  arma::vec m_pvalVec;
  arma::vec m_zScoreVec;
  // arma::uvec m_posOutlier;
  // arma::uvec m_posNonOutlier;
  // 
  // arma::vec m_R_outlier;
  // arma::vec m_R_nonOutlier;
  
  
  // 其他成员变量...
  // Rcpp::DataFrame m_sparseGRM;   // 新增：稀疏GRM数据框
  // Rcpp::DataFrame m_ResidMat;    // 新增：残差矩阵数据框
  
  
public:
  
  SPAGxEmixPlusClass(arma::mat t_resid,
                     arma::mat t_resid_by_E,  // 新增参数：
                     arma::mat t_PCs,
                     int t_N,
                     double t_SPA_Cutoff,
                     Rcpp::List t_outlierList,
                     Rcpp::DataFrame t_sparseGRM, // update by Yuzhuo Ma
                     Rcpp::DataFrame t_ResidMat,  // update by Yuzhuo Ma
                     arma::vec t_E,               // update by Yuzhuo Ma
                     std::string t_ResidTraitType,
                     arma::mat t_PhenoMat,
                     arma::mat t_Covariates);              // update by Yuzhuo Ma
  
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
      
      // 新增：显式检查 NaN
      if (std::isnan(diffX)) {
        diffX = 5;
      }
      
      
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
      
      
      // 新增：显式检查 NaN
      if (std::isnan(diffX)) {
        diffX = 5;
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
    //std::cout << "SPAmixPlus_GetProb_SPA_G:\t" << std::endl;
    
    
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
    
    // double pval = arma::normcdf(arma::sign(lower_tail - 0.5) * (w + log(v/w) / w));
    // return pval;
    
    // 新增：处理数学异常
    double pval = 0.0;
    if (w == 0 || v == 0 || (v/w) <= 0) { // 防止除零和 log(负数)
      pval = 0.0;
    } else {
      double term = w + log(v/w) / w;
      pval = arma::normcdf(arma::sign(lower_tail - 0.5) * term);
    }
    
    // 新增：显式检查 NaN
    if (std::isnan(pval)) {
      pval = 0.0;
    }
    
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
    
    // NOTE (2023-04-20): the below part has been moved to SPAGxEmixPlus.cpp since they are shared for all genetic variants
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
  
  
  
  
  
  
  
  // // update  SPA  ////////////////////////////////////////////////////////////////
  // 
  // double getMarkerPval(arma::vec t_GVec,
  //                      double t_altFreq)  // later update score and variance here (2023-06-20)
  // {
  // 
  //   //std::cout << "SPAmixPlus_getMarkerPval" << std::endl;
  // 
  // 
  //   arma::vec time1 = getTime();
  // 
  //   arma::vec AFVec = getMAFest(t_GVec, t_altFreq);  // PCs are global variables and can be loaded when necessary
  // 
  //   m_MAFVec = AFVec;  // 存储到成员变量
  // 
  // 
  //   arma::vec time2 = getTime();
  //   arma::vec diffTime = time2 - time1;
  //   // std::cout << "part2" << std::endl;
  // 
  //   // std::cout << "(MAF) diffTime:\t" << diffTime << std::endl;
  // 
  //   m_diffTime2 += diffTime;
  // 
  //   arma::vec GVarVec = 2 * AFVec % (1 - AFVec);
  // 
  //   // (BWJ) 2023-06-20: Support multiple phenotypes
  //   // outLierList[[i]] = list(posValue = posValue - 1,
  //   //                         posOutlier = posOutlier - 1,
  //   //                         posNonOutlier = posNonOutlier - 1,
  //   //                         resid = mresid.temp[posValue],
  //   //                         resid2 = mresid.temp[posValue]^2,
  //   //                         residOutleir = mresid.temp[posOutlier],
  //   //                         residNonOutlier = mresid.temp[posNonOutlier])
  //   for(int i = 0; i < m_Npheno; i++){
  //     Rcpp::List tempOutlierList = m_outlierList[i];
  //     arma::uvec posValue = tempOutlierList["posValue"];
  //     arma::uvec posOutlier = tempOutlierList["posOutlier"];
  //     arma::uvec posNonOutlier = tempOutlierList["posNonOutlier"];
  //     arma::vec resid = tempOutlierList["resid"];
  //     arma::vec resid2 = tempOutlierList["resid2"];
  //     arma::vec residOutlier = tempOutlierList["residOutlier"];
  //     arma::vec residNonOutlier = tempOutlierList["residNonOutlier"];
  //     arma::vec resid2NonOutlier = tempOutlierList["resid2NonOutlier"];
  // 
  // 
  //     // 修改5：严格限制在posValue范围内计算
  //     arma::vec R_subset = resid.elem(posValue);       // 残差子集
  //     arma::vec GVar_subset = GVarVec.elem(posValue);  // 基因型方差子集
  //     // arma::vec AF_subset = AFVec.elem(posValue);      // MAF子集
  // 
  //     // 计算调整后的残差（仅考虑posValue个体）
  //     // arma::vec R_new = R_subset % sqrt(GVar_subset);
  //     // 计算调整残差 R_new = R * sqrt(g.var.est)
  //     arma::vec R_new = R_subset % arma::sqrt(GVar_subset);
  // 
  // 
  //     // 修改6：传入posValue和GVar子集
  //     double VarS = calculateSparseVariance(R_new, posValue);
  // 
  // 
  // 
  // 
  //     // ==== 5. 后续p值计算逻辑保持不变 ====
  // 
  //     // ===== 5. 统计量计算 =====
  //     double S = sum(t_GVec.elem(posValue) % R_subset);
  //     double S_mean = 2 * sum(R_subset % AFVec.elem(posValue));
  //     double zScore = (S - S_mean) / sqrt(VarS);
  //     // double zScore = -1 + (S - S_mean) / sqrt(VarS);  //  test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // 
  //     m_zScoreVec.at(i) = zScore;
  // 
  // 
  // 
  // 
  // 
  //     // double S = sum(t_GVec.elem(posValue) % resid);
  // 
  //     // double VarS = sum(resid2 % GVarVec.elem(posValue));
  // 
  //     // // updated on 2023-04-23
  //     // double S_mean = 2 * sum(resid % AFVec.elem(posValue)); // NOTE: I think S_mean is somewhat weird and should be checked later (2023-04-22)
  //     // double zScore = (S-S_mean) / sqrt(VarS);
  //     //
  //     // m_zScoreVec.at(i) = zScore;
  // 
  // 
  //     if(std::abs(zScore) < m_SPA_Cutoff){
  //       m_pvalVec.at(i) = arma::normcdf(-1*std::abs(zScore))*2;
  //       // m_pvalVec.at(i) = arma::normcdf(-1*std::abs(zScore));     //  test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // 
  //       continue;
  //       // return pval;
  //     }
  // 
  // 
  // 
  // 
  // 
  //     // ===== SPA-G调整（核心修改）=====
  //     time1 = getTime();
  // 
  //     // 计算方差比率
  //     double S_var_SPAmix = arma::sum(arma::square(R_subset) % GVar_subset);
  //     double Var_ratio = S_var_SPAmix / VarS;
  // 
  //     // 调整统计量
  //     double S_new = S * std::sqrt(Var_ratio);
  //     double S_mean_new = S_mean * std::sqrt(Var_ratio);
  // 
  //     // 对称统计量计算
  //     double S_upper = std::max(S_new, 2*S_mean_new - S_new);
  //     double S_lower = std::min(S_new, 2*S_mean_new - S_new);
  // 
  // 
  // 
  // 
  //     // std::cout << "part4" << std::endl;
  // 
  //     // double pval1 = GetProb_SPA_G(AFVec, m_resid, std::abs(t_zScore), false);
  //     // double pval2 = GetProb_SPA_G(AFVec, m_resid, -1*std::abs(t_zScore), true);
  // 
  //     // time1 = getTime();
  // 
  //     // put the below objects as a global object
  //     arma::vec MAF_outlier = AFVec.elem(posOutlier);
  //     arma::vec MAF_nonOutlier = AFVec.elem(posNonOutlier);
  // 
  //     double mean_nonOutlier = sum(residNonOutlier % MAF_nonOutlier) * 2;
  //     double var_nonOutlier = sum(resid2NonOutlier % MAF_nonOutlier % (1-MAF_nonOutlier)) * 2;
  // 
  // 
  // 
  // 
  // 
  //     // 调用SPA-G函数
  //     double pval1 = GetProb_SPA_G(MAF_outlier,
  //                                  residOutlier,
  //                                  S_upper,
  //                                  false,
  //                                  mean_nonOutlier,
  //                                  var_nonOutlier);
  // 
  //     double pval2 = GetProb_SPA_G(MAF_outlier,
  //                                  residOutlier,
  //                                  S_lower,
  //                                  true,
  //                                  mean_nonOutlier,
  //                                  var_nonOutlier);
  // 
  // 
  //     // double pval1 = GetProb_SPA_G(MAF_outlier,
  //     //                              residOutlier,
  //     //                              // std::abs(S-S_mean)+S_mean,
  //     //                              std::max(S_new, 2*S_mean_new - S_new),
  //     //                              false,
  //     //                              mean_nonOutlier,
  //     //                              var_nonOutlier);
  //     // double pval2 = GetProb_SPA_G(MAF_outlier,
  //     //                              residOutlier,
  //     //                              // -1*std::abs(S-S_mean)+S_mean,
  //     //                              std::min(S_new, 2*S_mean_new - S_new),
  //     //                              true,
  //     //                              mean_nonOutlier,
  //     //                              var_nonOutlier);
  // 
  //     time2 = getTime();
  //     diffTime = time2 - time1;
  // 
  //     // std::cout << "(SPA_G) diffTime:\t" << diffTime << std::endl;
  // 
  //     m_diffTime1 += diffTime;
  // 
  //     m_pvalVec.at(i) = pval1 + pval2;
  //   }
  // 
  //   // std::cout << "part5" << std::endl;
  //   double pval = 0; // we modify the codes to save pval to m_pvalVec and thus the function does not output pvalue any more.
  //   return pval;
  // }
  // 

  
  
  
  
  // update  SPAGxEmixPlus  ////////////////////////////////////////////////////////////////
  
  double getMarkerPval(arma::vec t_GVec,
                       double t_altFreq,  // later update score and variance here (2023-06-20)
                       std::string t_ResidTraitType,
                       arma::mat t_PhenoMat,
                       arma::mat t_Covariates,
                       double epsilon = 0.001)              
  {
    
    // Rcpp::Rcout << "Head t_GVec: " << t_GVec.head(5).t() << std::endl;
    Rcpp::Rcout << "Head t_altFreq: " << t_altFreq << std::endl;
    
    
    
    
    
    // std::cout << "ResidTraitType:\t" << t_ResidTraitType << std::endl; // update by Yuzhuo Ma
    // ==== 新增调试输出 ====
    // 检查矩阵是否为空
    // if(t_PhenoMat.n_elem == 0) {
    //   Rcpp::Rcout << "WARNING: t_PhenoMat is empty!" << std::endl;
    // } else {
    //   // Rcpp::Rcout << "=== Debug: Head of t_PhenoMat (first 5 rows) ===" << std::endl;
    //   // for(int i=0; i < 5 && i < t_PhenoMat.n_rows; ++i) {
    //     // Rcpp::Rcout << "Row " << i+1 << ": " << t_PhenoMat.row(i) << std::endl;
    //   }
    // }
    // 
    // if(t_Covariates.n_elem == 0) {
    //   Rcpp::Rcout << "WARNING: t_Covariates is empty!" << std::endl;
    // } else {
    //   // Rcpp::Rcout << "=== Debug: Head of t_Covariates (first 5 rows) ===" << std::endl;
    //   // for(int i=0; i < 5 && i < t_Covariates.n_rows; ++i) {
    //     // Rcpp::Rcout << "Row " << i+1 << ": " << t_Covariates.row(i) << std::endl;
    //   }
    // }
    
    
    // 新增输入校验
    if(t_ResidTraitType == "Binary") {
      if(t_PhenoMat.n_elem == 0 || t_Covariates.n_elem == 0)
        Rcpp::stop("Binary性状必须提供PhenoMat和Covariates矩阵");
    }
    
    // std::cout << "SPAGxEmixPlus_getMarkerPval part1" << std::endl; // update by Yuzhuo Ma
    
    
    arma::vec time1 = getTime();
    
    arma::vec AFVec = getMAFest(t_GVec, t_altFreq);  // PCs are global variables and can be loaded when necessary
    
    m_MAFVec = AFVec;  // 存储到成员变量
    
    
    arma::vec time2 = getTime();
    arma::vec diffTime = time2 - time1;
    // std::cout << "part2" << std::endl;
    
    // std::cout << "(MAF) diffTime:\t" << diffTime << std::endl;
    
    m_diffTime2 += diffTime;
    
    arma::vec GVarVec = 2 * AFVec % (1 - AFVec);
    
    // std::cout << "SPAGxEmixPlus_getMarkerPval part2" << std::endl; // update by Yuzhuo Ma
    
    
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
      
      
      // ==== 新增检查：posValue是否为空 ====
      if(posValue.n_elem == 0) {
        Rcpp::stop("No valid samples (posValue is empty) for phenotype %d.", i);
      }
      
      
      // // +++ 修改数据提取方式（使用成员变量中的完整数据）+++
      // arma::vec y_full = m_PhenoMat.col(i); // 使用类成员变量
      // arma::vec y = y_full.elem(posValue);  // 用全局索引取子集
      // 
      // arma::mat cov_full = m_Covariates;    // 使用类成员变量
      // arma::mat cov_subset = cov_full.rows(posValue);
      
      ///////////////////////////////////////////////////////////////////////////////////
      
      // 修改5：严格限制在posValue范围内计算
      arma::vec R_subset = resid.elem(posValue);       // 残差子集
      arma::vec GVar_subset = GVarVec.elem(posValue);  // 基因型方差子集
      // arma::vec AF_subset = AFVec.elem(posValue);      // MAF子集
      
      // 计算调整后的残差（仅考虑posValue个体）
      // arma::vec R_new = R_subset % sqrt(GVar_subset);
      // 计算调整残差 R_new = R * sqrt(g.var.est)
      arma::vec R_sqrtMAF = R_subset % arma::sqrt(GVar_subset);
      
      
      // 修改6：传入posValue和GVar子集
      double VarS1 = calculateSparseVariance(R_sqrtMAF, posValue);
      // std::cout << "VarS1:\t" << VarS1 << std::endl; // update by Yuzhuo Ma
      
      
      
      
      // ==== 5. 后续p值计算逻辑保持不变 ====
      
      // ===== 5. 统计量计算 =====
      double S1 = sum(t_GVec.elem(posValue) % R_subset);
      double S1_mean = 2 * sum(R_subset % AFVec.elem(posValue));
      double zScore = (S1 - S1_mean) / sqrt(VarS1);
      // double zScore = -1 + (S - S_mean) / sqrt(VarS);  //  test !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
      m_zScoreVec.at(i) = zScore;
      
      
      // std::cout << "S1:\t" << S1 << std::endl; // update by Yuzhuo Ma
      // std::cout << "S1_mean:\t" << S1_mean << std::endl; // update by Yuzhuo Ma
      // std::cout << "zScore:\t" << zScore << std::endl; // update by Yuzhuo Ma
      
      
      
      // update /////////////////////////////////////////////////////////////////
      double pval_norm1 = arma::normcdf(-1*std::abs(zScore))*2;

      // std::cout << "pval_norm1:\t" << pval_norm1 << std::endl; // update by Yuzhuo Ma
      
      
      
      // 在getMarkerPval或类似函数中添加以下逻辑
      if(pval_norm1 > epsilon){
        
        arma::vec E_subset = m_E.elem(posValue);
        arma::vec RE_subset = m_resid_by_E.elem(posValue);
        
        // ===== GxE效应计算 =====
        // 步骤1：计算基础统计量
        // arma::vec RE = E % R; // E为环境因子，R为残差
        double S2 = arma::sum(t_GVec.elem(posValue) % RE_subset);
        
        // 步骤2：计算调整残差
        arma::vec RE_sqrtMAF = RE_subset % arma::sqrt(GVar_subset);
        // arma::vec R_sqrtMAF = R % arma::sqrt(gVarEst);
        
        // 步骤3：交叉协方差计算（优化实现）
        double Cov_S1_S2 = calculateCrossCovariance(R_sqrtMAF, RE_sqrtMAF, posValue) 
          - arma::sum(GVar_subset % E_subset % arma::square(R_subset));
        
        // std::cout << "Cov_S1_S2:\t" << Cov_S1_S2 << std::endl; // update by Yuzhuo Ma
        
        
        // 步骤4：计算lambda和新残差
        double lambda = Cov_S1_S2 / VarS1;
        arma::vec R_new = (E_subset - lambda) % R_subset;
        
        // std::cout << "lambda:\t" << lambda << std::endl; // update by Yuzhuo Ma
        
        
        // 步骤5：计算新残差的调整方差
        arma::vec R_new_sqrtMAF = R_new % arma::sqrt(GVar_subset);
        double S_GxE_var = calculateSparseVariance(R_new_sqrtMAF, posValue);
        
        // 步骤6：标准化统计量
        double S_GxE = S2 - lambda * S1;
        double S_GxE_mean = 2 * arma::sum(R_new % AFVec.elem(posValue));
        double z_GxE = (S_GxE - S_GxE_mean) / std::sqrt(S_GxE_var);
        
        
        // std::cout << "S_GxE:\t" << S_GxE << std::endl; // update by Yuzhuo Ma
        // std::cout << "S_GxE_mean:\t" << S_GxE_mean << std::endl; // update by Yuzhuo Ma
        // std::cout << "z_GxE:\t" << z_GxE << std::endl; // update by Yuzhuo Ma
        
        
        // ===== p值计算 =====
        if(std::abs(z_GxE) < m_SPA_Cutoff){
          // ... 正态近似 ...
          double pval_norm_GxE = arma::normcdf(-1*std::abs(z_GxE))*2;
          m_pvalVec.at(i) = pval_norm_GxE;
        }else{
          // 方差比率调整
          double S_GxE_var_SPA = arma::sum(arma::square(R_new) % GVar_subset);
          double Var_ratio = S_GxE_var_SPA / S_GxE_var;
          
          // 对称统计量边界
          double S_upper = std::max(S_GxE, 2*S_GxE_mean - S_GxE);
          double S_lower = std::min(S_GxE, 2*S_GxE_mean - S_GxE);
          
          // // SPA-G计算（含异常处理）
          // double pval1 = GetProb_SPA_G(MAF_outlier, residOutlier, 
          //                              S_upper * std::sqrt(Var_ratio), false,
          //                              mean_nonOutlier, var_nonOutlier);
          // double pval2 = GetProb_SPA_G(MAF_outlier, residOutlier, 
          //                              S_lower * std::sqrt(Var_ratio), true,
          //                              mean_nonOutlier, var_nonOutlier);
          
          // SPA-G计算（含异常处理）
          double pval1 = GetProb_SPA_G(AFVec.elem(posValue),
                                       R_new, 
                                       S_upper * std::sqrt(Var_ratio), 
                                       false,
                                       0,
                                       0);
          
          double pval2 = GetProb_SPA_G(AFVec.elem(posValue),
                                       R_new, 
                                       S_lower * std::sqrt(Var_ratio),
                                       true,
                                       0,
                                       0);
          
          // 处理NaN
          pval1 = std::isnan(pval1) ? 0.0 : pval1;
          pval2 = std::isnan(pval2) ? 0.0 : pval2;
          m_pvalVec.at(i) = pval1 + pval2;
          
          // std::cout << "pval1:\t" << pval1 << std::endl; // update by Yuzhuo Ma
          // std::cout << "pval2:\t" << pval2 << std::endl; // update by Yuzhuo Ma
        }
      }else{
        
        // arma::vec R0; // 在父作用域声明
        arma::vec R0 = arma::vec(); // 显式初始化为空向量
        
        
        
        // quantitative phenotype
        
        if(t_ResidTraitType == "Quantitative"){
          // 1. 构造W矩阵
          arma::mat W = arma::join_rows(arma::ones(posValue.n_elem), t_GVec.elem(posValue));
          
          // 2. 投影矩阵计算 (R0 = R - W(W^TW)^{-1}W^TR)
          arma::mat WtW = W.t() * W;
          arma::mat inv_WtW = arma::inv(WtW);
          arma::vec R_subset = resid.elem(posValue);  // 从原有代码获取
          // ===== 修改点1：移除内部变量声明 =====
          R0 = R_subset - W * inv_WtW * W.t() * R_subset;
        }else if(t_ResidTraitType == "Binary") {
          // 检查输入完整性
          if(t_PhenoMat.n_elem == 0 || t_Covariates.n_elem == 0) {
            Rcpp::stop("Binary性状需要非空表型矩阵和协变量矩阵");
          }
          
          
          // +++ 修改数据提取方式（使用成员变量中的完整数据）+++
          arma::vec y_full = m_PhenoMat.col(i); // 使用类成员变量
          // 调试输出实际表型值
          Rcpp::Rcout << "Unique y_full values: " 
                      << arma::unique(y_full).t() << std::endl;
          
          
          arma::vec y = arma::round(y_full.elem(posValue));  // 四舍五入消除浮点误差
          
          // 调试输出实际表型值
          Rcpp::Rcout << "Unique y values: " 
                      << arma::unique(y).t() << std::endl;
          
          y = arma::clamp(y, 0, 1);  // 确保严格0/1
          
          // 调试输出实际表型值
          Rcpp::Rcout << "Unique y values: " 
                      << arma::unique(y).t() << std::endl;
          
          arma::mat cov_full = m_Covariates;    // 使用类成员变量
          // ==== 关键修复1：强制协变量为数值矩阵 ====
          arma::mat cov_subset = arma::conv_to<arma::mat>::from(
            cov_full.rows(posValue).eval() // 确保内存连续性
          );      
          

          
          
          // // 获取当前表型列（第i列）
          // // arma::vec y = t_PhenoMat.col(i).elem(posValue); 
          // // 修复后代码
          // arma::vec y_col = t_PhenoMat.col(i).eval();
          // arma::vec y = y_col.elem(posValue);
          // 
          // // 提取协变量子集（需与posValue对齐）
          // // arma::mat cov_subset = t_Covariates.rows(posValue);
          // // 修复后代码
          // arma::mat cov_subset = t_Covariates.rows(posValue).eval(); // 强制转换为连续内存矩阵
          

          
          // // 原错误代码：
          // // arma::vec y_col = t_PhenoMat.col(i).eval();
          // // 修改后：
          // arma::vec y_full = m_PhenoMat.col(i); // 使用成员变量获取完整表型
          // arma::vec y = y_full.elem(posValue);  // 用全局索引取子集
          // 
          // arma::mat cov_full = m_Covariates;    // 使用成员变量获取完整协变量
          // arma::mat cov_subset = cov_full.rows(posValue);
          
          // +++ 添加调试输出 +++
          Rcpp::Rcout << "Checking binary phenotype values:\n";
          Rcpp::Rcout << "Min: " << arma::min(y) 
                      << " Max: " << arma::max(y) 
                      << " N(1): " << arma::sum(y) << "\n";
          
          // +++ 新增响应变量有效性检查 +++
          if(arma::any(y < 0) || arma::any(y > 1)) {
            Rcpp::stop("Invalid phenotype values (must be 0 or 1)");
          }
          
          

          
          R0 = calculateGLMResidual_R(y, 
                                      t_GVec.elem(posValue),
                                      cov_subset);
          
     

          // 调试输出（与图片中的Rcout位置对应）
          Rcpp::Rcout << "Head R0: " << R0.head(5).t() << std::endl;
          
          // +++ 新增残差检查 +++
          if(arma::var(R0) < 1e-6) {
            Rcpp::Rcout << "WARNING: Constant residuals detected.\n";
          }
          
        }
        
        
        // ==== 后续统一使用父作用域的R0 ====
        if(R0.is_empty()) {
          Rcpp::stop("残差计算失败，R0为空");
        }
        
        // 在getMarkerPval函数中添加状态日志
        Rcpp::Rcout << "当前表型类型: " << t_ResidTraitType 
                    << " | R0维度: " << R0.n_elem << std::endl;
        
        Rcpp::Rcout << "Head R0: " << R0.head(5).t() << std::endl;
        
        
        //////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // 3. 基础统计量计算
        arma::vec E_subset = m_E.elem(posValue);  // 假设E是成员变量
        
        Rcpp::Rcout << "E_subset: " << E_subset.head(5).t() << std::endl;
        
        
        double S_GxE0 = arma::sum(t_GVec.elem(posValue) % E_subset % R0);
        
        std::cout << "S_GxE0:\t" << S_GxE0 << std::endl; // update by Yuzhuo Ma
        
        
        // 4. 调整残差计算
        arma::vec R_new0 = E_subset % R0;
        
        Rcpp::Rcout << "R_new0: " << R_new0.head(5).t() << std::endl;
        
        
        arma::vec R_new0_sqrtMAF = R_new0 % arma::sqrt(GVar_subset);
        
        Rcpp::Rcout << "R_new0_sqrtMAF: " << R_new0_sqrtMAF.head(5).t() << std::endl;
        
        
        
        // 5. 协方差计算（需要已实现的函数）
        // double Cov_sum = calculateSparseGRMCovariance(R_new0_sqrtMAF); 
        double S_GxE0_var = calculateSparseVariance(R_new0_sqrtMAF, posValue);
        
        std::cout << "S_GxE0_var:\t" << S_GxE0_var << std::endl; // update by Yuzhuo Ma
        
        
        // // 步骤3：交叉协方差计算（优化实现）
        // double Cov_S1_S2 = calculateCrossCovariance(R_sqrtMAF, RE_sqrtMAF, posValue)
        //   - arma::sum(GVar_subset % E_subset % arma::square(R_subset));
        
        
        // 6. 标准化统计量
        double S_GxE0_mean = 2 * arma::sum(R_new0 % AFVec.elem(posValue));
        double z_GxE0 = (S_GxE0 - S_GxE0_mean) / std::sqrt(S_GxE0_var);
        
        std::cout << "S_GxE0_mean:\t" << S_GxE0_mean << std::endl; // update by Yuzhuo Ma
        std::cout << "z_GxE0:\t" << z_GxE0 << std::endl; // update by Yuzhuo Ma
        
        
        // 7. p值计算分支
        if(std::abs(z_GxE0) < m_SPA_Cutoff) {
          m_pvalVec.at(i) = 2 * arma::normcdf(-std::abs(z_GxE0));
        } else {
          double S_GxE0_var_SPA = arma::sum(arma::square(R_new0) % GVar_subset);
          double Var_ratio = S_GxE0_var_SPA / S_GxE0_var;
          
          double S_upper = std::max(S_GxE0, 2*S_GxE0_mean - S_GxE0);
          double S_lower = std::min(S_GxE0, 2*S_GxE0_mean - S_GxE0);
          
          double pval1 = GetProb_SPA_G(AFVec.elem(posValue),
                                       R_new0,
                                       S_upper * std::sqrt(Var_ratio),
                                       false, 0, 0);
          
          double pval2 = GetProb_SPA_G(AFVec.elem(posValue),
                                       R_new0,
                                       S_lower * std::sqrt(Var_ratio),
                                       true, 0, 0);
          
          // 异常值处理
          pval1 = std::isnan(pval1) ? 0.0 : pval1;
          pval2 = std::isnan(pval2) ? 0.0 : pval2;
          m_pvalVec.at(i) = pval1 + pval2;
          
          std::cout << "pval1:\t" << pval1 << std::endl; // update by Yuzhuo Ma
          std::cout << "pval2:\t" << pval2 << std::endl; // update by Yuzhuo Ma
          
        }
        

        
        
      }
      
      
      
      

    }
    
    // std::cout << "part5" << std::endl;
    double pval = 0; // we modify the codes to save pval to m_pvalVec and thus the function does not output pvalue any more.
    return pval;
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  
  
  
  // // ==== 新增4：稀疏矩阵方差计算 ====
  // double calculateSparseVariance(const arma::vec& R_new) {
  //   double covSum = 0.0;
  //   for(const auto& triplet : m_sparseTriplets){
  //     int i = std::get<0>(triplet);
  //     int j = std::get<1>(triplet);
  //     double val = std::get<2>(triplet);
  //     covSum += val * R_new(i) * R_new(j);
  //   }
  //   return 2 * covSum - sum(square(R_new) % (2 * m_MAFVec % (1 - m_MAFVec)));
  // }
  
  
  
  //  优化方差计算函数 2025-04-09
  
  // ==== 最终正确的方差计算函数 ====
  double calculateSparseVariance(const arma::vec& R_new, 
                                 const arma::uvec& posValue) 
  {
    double covSum = 0.0;
    
    // 快速索引查询优化
    std::unordered_set<int> validIndices;
    for (auto idx : posValue) {
      validIndices.insert(static_cast<int>(idx));
    }
    
    // 遍历稀疏GRM三元组
    for (const auto& triplet : m_sparseTriplets) {
      int i = std::get<0>(triplet);
      int j = std::get<1>(triplet);
      
      // 仅处理双方都在posValue中的个体
      if (validIndices.count(i) && validIndices.count(j)) {
        double grmValue = std::get<2>(triplet);
        covSum += grmValue * R_new(i) * R_new(j);
      }
    }
    
    // 严格遵循R代码公式
    return 2.0 * covSum - arma::sum(arma::square(R_new));
  }
  
  
  
  /////////////////////////////////////////////////////////////////////////////////////////
  
  // 在cpp文件中实现核心逻辑
  double calculateCrossCovariance(const arma::vec& R_sqrtMAF,
                                  const arma::vec& RE_sqrtMAF,
                                  const arma::uvec& posValue)
  {
    double covSum_part1 = 0.0, covSum_part2 = 0.0;
    
    // 创建快速索引查询表
    std::unordered_map<int, int> indexMap;
    for(size_t i=0; i<posValue.n_elem; ++i){
      indexMap[posValue[i]] = i;
    }
    
    // 遍历稀疏GRM三元组
    for(const auto& triplet : m_sparseTriplets){
      int i = std::get<0>(triplet);
      int j = std::get<1>(triplet);
      
      // 检查双方是否都在有效样本中
      if(indexMap.count(i) && indexMap.count(j)){
        double value = std::get<2>(triplet);
        int idx_i = indexMap[i];
        int idx_j = indexMap[j];
        
        // 计算两个交叉项
        covSum_part1 += value * R_sqrtMAF[idx_i] * RE_sqrtMAF[idx_j];
        covSum_part2 += value * R_sqrtMAF[idx_j] * RE_sqrtMAF[idx_i];
      }
    }
    
    return covSum_part1 + covSum_part2;
  }
  
  
  
  
  
  // // 修改后的逻辑回归原始残差计算函数
  // arma::vec calculateLogisticRawResidual(const arma::vec& y, const arma::vec& g) {
  //   int n = y.n_elem;
  //   
  //   // 1. 构造设计矩阵（包含截距项）
  //   arma::mat X = arma::join_horiz(arma::ones<arma::vec>(n), g);
  //   
  //   // 2. 初始化参数（beta[0]截距，beta[1]基因型效应）
  //   arma::vec beta = arma::zeros<arma::vec>(2);
  //   
  //   // 3. 迭代加权最小二乘法（IRLS）
  //   double tol = 1e-6;
  //   int max_iter = 25;
  //   for (int iter = 0; iter < max_iter; ++iter) {
  //     arma::vec eta = X * beta;
  //     arma::vec mu = 1.0 / (1.0 + arma::exp(-eta));
  //     
  //     // 计算权重矩阵（对角阵）
  //     arma::vec weights = mu % (1 - mu);
  //     arma::mat W = arma::diagmat(weights);
  //     
  //     // 工作响应量计算（带数值稳定处理）
  //     arma::vec z = eta + (y - mu) / (weights + 1e-16);
  //     
  //     // 更新参数（优化矩阵运算）
  //     arma::mat XtW = X.t() * W;
  //     arma::vec beta_new = arma::solve(XtW * X, XtW * z);
  //     
  //     // 收敛判断（L2范数变化）
  //     if (arma::norm(beta_new - beta, 2) < tol) break;
  //     beta = beta_new;
  //   }
  //   
  //   // 4. 计算原始残差（观测值 - 预测概率）
  //   arma::vec mu = 1.0 / (1.0 + arma::exp(X * beta));
  //   arma::vec R0 = y - mu;  // 原始残差
  //   
  //   return R0;
  // }
  // 
  
  
  
  
  
  
  // 修改后的calculateGLMResidual函数
  // 新增函数：带协变量的逻辑回归残差计算
  // =============== 在calculateGLMResidual函数中修改 ===============
  // 修改位置：图片第906-915行附近
  // 原错误代码构造数据框部分：
  // Rcpp::DataFrame data = Rcpp::DataFrame::create(...)
  
  // // 修改后代码（完整函数）：
  // arma::vec calculateGLMResidual(const arma::vec& y, 
  //                                const arma::vec& g,
  //                                const arma::mat& covariates) 
  // {
  //   try {
  //     // ==== 输入校验增强 ====
  //     if (y.n_elem != g.n_elem || y.n_elem != covariates.n_rows) 
  //       Rcpp::stop("维度不匹配: y(%d), g(%d), cov(%dx%d)", 
  //                  y.n_elem, g.n_elem, covariates.n_rows, covariates.n_cols);
  //     
  //     // ==== 显式类型转换 ====
  //     Rcpp::NumericVector y_vec = Rcpp::wrap(y);  // Armadillo转Rcpp向量
  //     Rcpp::NumericVector g_vec = Rcpp::wrap(g);
  //     
  //     // ==== 构造数据框 ====
  //     Rcpp::DataFrame data = Rcpp::DataFrame::create(
  //       Rcpp::Named("response") = y_vec,
  //       Rcpp::Named("genotype") = g_vec,
  //       Rcpp::_["stringsAsFactors"] = false
  //     );
  //     
  //     // ==== 协变量列添加优化 ====
  //     for(int i=0; i<covariates.n_cols; i++) {
  //       std::string colname = "cov" + std::to_string(i+1);
  //       data[colname] = Rcpp::NumericVector(covariates.colptr(i), 
  //                                           covariates.colptr(i)+covariates.n_rows);
  //     }
  //     
  //     // ==== 动态公式构造 ====
  //     std::string formula_str = "response ~ genotype";
  //     for(int i=0; i<covariates.n_cols; i++) {
  //       formula_str += " + cov" + std::to_string(i+1);
  //     }
  //     Rcpp::Formula formula(formula_str);
  //     
  //     // ==== 调用R的glm ====
  //     Rcpp::Function glm = Rcpp::Environment::namespace_env("stats")["glm"];
  //     Rcpp::List model = glm(
  //       formula,
  //       Rcpp::_["data"] = data,
  //       Rcpp::_["family"] = Rcpp::Function("binomial")()
  //     );
  //     
  //     // ==== 残差提取 ====
  //     Rcpp::NumericVector mu = model["fitted.values"];
  //     return y_vec - mu;  // 返回原始残差
  //     
  //   } catch(...) {
  //     return arma::vec(); // 错误处理
  //   }
  // }
  
  // 
  // arma::vec calculateGLMResidual(const arma::vec& y, 
  //                                const arma::vec& g,
  //                                const arma::mat& covariates) 
  // {
  //   try {
  //     // ==== 输入校验 ====
  //     if (y.n_elem != g.n_elem || y.n_elem != covariates.n_rows) 
  //       Rcpp::stop("维度不匹配: y(%d), g(%d), cov(%dx%d)", 
  //                  y.n_elem, g.n_elem, covariates.n_rows, covariates.n_cols);
  //     
  //     // ==== 构造设计矩阵（包含截距项和协变量）====
  //     arma::mat X = arma::join_horiz(arma::ones<arma::vec>(y.n_elem), g, covariates);
  //     
  //     // ==== 迭代加权最小二乘法（IRLS）====
  //     int max_iter = 25;
  //     double tol = 1e-6;
  //     arma::vec beta(X.n_cols, arma::fill::zeros); // 初始化系数
  //     
  //     for(int iter = 0; iter < max_iter; ++iter) {
  //       arma::vec eta = X * beta;
  //       arma::vec mu = 1.0 / (1.0 + arma::exp(-eta));
  //       
  //       // 计算权重和对角矩阵
  //       arma::vec W_vec = mu % (1 - mu);
  //       arma::mat W = arma::diagmat(W_vec);
  //       
  //       // 工作响应量（带数值稳定处理）
  //       arma::vec z = eta + (y - mu) / (W_vec + 1e-16);
  //       
  //       // 更新系数
  //       arma::mat XtW = X.t() * W;
  //       arma::vec beta_new = arma::solve(XtW * X, XtW * z, arma::solve_opts::equilibrate);
  //       
  //       // 收敛检查
  //       if (arma::norm(beta_new - beta, 2) < tol) break;
  //       beta = beta_new;
  //     }
  //     
  //     // ==== 计算残差（观测值 - 预测概率）====
  //     arma::vec mu = 1.0 / (1.0 + arma::exp(X * beta));
  //     return y - mu;
  //     
  //   } catch(...) {
  //     Rcpp::Rcerr << "逻辑回归计算失败" << std::endl;
  //     return arma::vec();
  //   }
  // }
  
  // arma::vec calculateGLMResidual(const arma::vec& y, 
  //                                const arma::vec& g,
  //                                const arma::mat& covariates) 
  // {
  //   try {
  //     // 输入维度校验
  //     if (y.n_elem != g.n_elem || y.n_elem != covariates.n_rows) {
  //       Rcpp::stop("Dimension mismatch: y(%d), g(%d), covariates(%d x %d).", 
  //                  y.n_elem, g.n_elem, covariates.n_rows, covariates.n_cols);
  //     }
  //     
  //     // 构造设计矩阵（包含截距项、基因型和协变量）
  //     arma::mat X = arma::join_horiz(arma::ones<arma::vec>(y.n_elem), g, covariates);
  //     
  //     // IRLS拟合逻辑回归模型
  //     int max_iter = 25;
  //     double tol = 1e-6;
  //     arma::vec beta(X.n_cols, arma::fill::zeros); // 初始化系数
  //     
  //     for(int iter = 0; iter < max_iter; ++iter) {
  //       arma::vec eta = X * beta;
  //       arma::vec mu = 1.0 / (1.0 + arma::exp(-eta));
  //       arma::vec W_vec = mu % (1 - mu);
  //       arma::mat W = arma::diagmat(W_vec);
  //       arma::vec z = eta + (y - mu) / (W_vec + 1e-16); // 添加小数避免除零
  //       
  //       arma::mat XtW = X.t() * W;
  //       arma::vec beta_new = arma::solve(XtW * X, XtW * z, arma::solve_opts::equilibrate);
  //       
  //       if (arma::norm(beta_new - beta, "fro") < tol) 
  //         break;
  //       beta = beta_new;
  //     }
  //     
  //     // 计算残差：y - 预测概率
  //     arma::vec mu = 1.0 / (1.0 + arma::exp(X * beta));
  //     return y - mu;
  //     
  //   } catch(const std::exception& e) {
  //     Rcpp::Rcerr << "Error in logistic regression: " << e.what() << std::endl;
  //     return arma::vec(); // 返回空向量表示错误
  //   }
  // }
  
  
  
  
  

  
  // // [[Rcpp::export]]
  // arma::vec calculateGLMResidual_R(const arma::vec& y, 
  //                                  const arma::vec& g,
  //                                  const arma::mat& covariates) {
  //   try {
  //     // =============== 输入校验增强 ===============
  //     if(y.n_elem != g.n_elem || y.n_elem != covariates.n_rows) {
  //       Rcpp::stop("Dimension mismatch: y(%d), g(%d), covariates(%d x %d)",
  //                  y.n_elem, g.n_elem, covariates.n_rows, covariates.n_cols);
  //     }
  //     
  //     // =============== 数据转换优化 ===============
  //     Rcpp::NumericVector r_y = Rcpp::wrap(y);
  //     Rcpp::NumericVector r_g = Rcpp::wrap(g);
  //     
  //     // =============== 协变量处理关键修正 ===============
  //     Rcpp::List cov_list;
  //     for(int i=0; i<covariates.n_cols; ++i) {
  //       std::string colname = "V" + std::to_string(i+1);
  //       cov_list[colname] = Rcpp::NumericVector(
  //         covariates.colptr(i),
  //         covariates.colptr(i) + covariates.n_rows
  //       );
  //     }
  //     Rcpp::DataFrame cov_df(cov_list);
  //     
  //     // =============== 合并数据集 ===============
  //     Rcpp::DataFrame df = Rcpp::DataFrame::create(
  //       Rcpp::Named("response") = r_y,
  //       Rcpp::Named("genotype") = r_g,
  //       cov_df,
  //       Rcpp::_["stringsAsFactors"] = false
  //     );
  //     
  //     // =============== 动态公式构建 ===============
  //     std::string formula_str = "response ~ genotype";
  //     for(int i=0; i<covariates.n_cols; ++i) {
  //       formula_str += " + V" + std::to_string(i+1);
  //     }
  //     
  //     // =============== R环境调用优化 ===============
  //     Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  //     Rcpp::Function glm = stats["glm"];
  //     
  //     // =============== 模型拟合 ===============
  //     Rcpp::List model;
  //     try {
  //       model = glm(
  //         Rcpp::Formula(formula_str),
  //         Rcpp::_["data"]    = df,
  //         Rcpp::_["family"]  = Rcpp::Function("binomial")()
  //       );
  //     } catch(const std::exception& e) {
  //       Rcpp::Rcerr << "GLM Fitting Error: " << e.what() << std::endl;
  //       return arma::vec(y.n_elem, arma::fill::zeros);
  //     }
  //     
  //     // =============== 残差提取与校验 ===============
  //     // 修正1：使用正确的元素存在性检查方法
  //     if(!model.containsElementNamed("fitted.values")) {
  //       Rcpp::Rcerr << "Model object missing fitted.values" << std::endl;
  //       return arma::vec(y.n_elem, arma::fill::zeros);
  //     }
  //     
  //     Rcpp::NumericVector fitted = model["fitted.values"];
  //     
  //     // 维度二次校验
  //     if(fitted.size() != y.n_elem) {
  //       Rcpp::Rcerr << "Fitted values dimension mismatch: " 
  //                   << fitted.size() << " vs " << y.n_elem << std::endl;
  //       return arma::vec(y.n_elem, arma::fill::zeros);
  //     }
  //     
  //     // 修正2：正确的类型转换
  //     return Rcpp::as<arma::vec>(r_y) - Rcpp::as<arma::vec>(fitted);
  //     
  //   } catch(const std::exception& e) {
  //     Rcpp::Rcerr << "GLM Calculation Error: " << e.what() << std::endl;
  //     return arma::vec(y.n_elem, arma::fill::zeros);
  //   } catch(...) {
  //     Rcpp::Rcerr << "Unknown error occurred in GLM calculation" << std::endl;
  //     return arma::vec(y.n_elem, arma::fill::zeros);
  //   }
  // }
  

  
  
  
  
  // 
  // // [[Rcpp::export]]
  // arma::vec calculateGLMResidual_R(const arma::vec& y, 
  //                                  const arma::vec& g,
  //                                  const arma::mat& covariates) {
  //   try {
  //     // =============== 输入校验增强 ===============
  //     if(y.n_elem != g.n_elem || y.n_elem != covariates.n_rows) {
  //       Rcpp::stop("Dimension mismatch: y(%d), g(%d), covariates(%d x %d)",
  //                  y.n_elem, g.n_elem, covariates.n_rows, covariates.n_cols);
  //     }
  //     
  //     // =============== 协变量列名唯一性 ===============
  //     Rcpp::List cov_list;
  //     std::set<std::string> used_names{"response", "genotype"};
  //     for(int i=0; i<covariates.n_cols; ++i) {
  //       std::string base = "V";
  //       int counter = 1;
  //       std::string colname;
  //       do {
  //         colname = base + std::to_string(counter++);
  //       } while(used_names.count(colname));
  //       
  //       used_names.insert(colname);
  //       cov_list[colname] = Rcpp::NumericVector(
  //         covariates.colptr(i),
  //         covariates.colptr(i) + covariates.n_rows
  //       );
  //     }
  //     
  //     // =============== 合并数据集 ===============
  //     Rcpp::DataFrame df = Rcpp::DataFrame::create(
  //       Rcpp::Named("response") = Rcpp::wrap(y),
  //       Rcpp::Named("genotype") = Rcpp::wrap(g),
  //       cov_list,
  //       Rcpp::_["stringsAsFactors"] = false
  //     );
  //     
  //     // =============== 动态公式构建 ===============
  //     std::string formula_str = "response ~ genotype";
  //     for(int i=0; i<covariates.n_cols; ++i) {
  //       formula_str += " + V" + std::to_string(i+1);
  //     }
  //     
  //     // =============== 模型拟合 ===============
  //     Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
  //     Rcpp::Function glm = stats["glm"];
  //     Rcpp::List model = glm(
  //       Rcpp::Formula(formula_str),
  //       Rcpp::_["data"]    = df,
  //       Rcpp::_["family"]  = Rcpp::Function("binomial")()
  //     );
  //     
  //     // =============== 残差严格校验 ===============
  //     if(!model.containsElementNamed("fitted.values")) {
  //       Rcpp::Rcerr << "Missing fitted values in model object" << std::endl;
  //       return arma::vec(y.n_elem, arma::fill::zeros);
  //     }
  //     
  //     Rcpp::NumericVector fitted = model["fitted.values"];
  //     arma::vec residuals = Rcpp::as<arma::vec>(Rcpp::wrap(y)) - Rcpp::as<arma::vec>(fitted);
  //     
  //     if(residuals.n_elem != y.n_elem) {
  //       Rcpp::Rcerr << "Residual dimension corrupted! Actual: " 
  //                   << residuals.n_elem << " vs Expected: " << y.n_elem << std::endl;
  //       return arma::vec(y.n_elem, arma::fill::zeros);
  //     }
  //     
  //     return residuals;
  //     
  //   } catch(const std::exception& e) {
  //     Rcpp::Rcerr << "GLM Error: " << e.what() << std::endl;
  //     return arma::vec(y.n_elem, arma::fill::zeros);
  //   } catch(...) {
  //     Rcpp::Rcerr << "Unknown GLM error" << std::endl;
  //     return arma::vec(y.n_elem, arma::fill::zeros);
  //   }
  // }
  
  
  
  
  // [[Rcpp::export]]
  arma::vec calculateGLMResidual_R(const arma::vec& y, 
                                   const arma::vec& g,
                                   const arma::mat& covariates) 
  {
    try {
      // =============== 初始化与基本校验 ===============
      const int n = y.n_elem;
      arma::vec residuals(n, arma::fill::zeros);
      
      // 调试信息：输出关键统计量
      Rcpp::Rcout << "\n=== DEBUG: calculateGLMResidual_R ===" << std::endl;
      Rcpp::Rcout << "Sample size: " << n << std::endl;
      Rcpp::Rcout << "Genotype summary: "
                  << "min=" << arma::min(g) << ", "
                  << "max=" << arma::max(g) << ", "
                  << "mean=" << arma::mean(g) << std::endl;
      Rcpp::Rcout << "Phenotype distribution: "
                  << "N0=" << arma::sum(y == 0) << ", "
                  << "N1=" << arma::sum(y == 1) << std::endl;
      
      // =============== 异常情况处理 ===============
      // 情况1：全零基因型
      if(arma::all(g < 1e-6)) {
        Rcpp::Rcout << "WARNING: All genotypes are zero. Returning zero residuals." << std::endl;
        return residuals;
      }
      
      // 情况2：无效表型值
      if(arma::any(y < 0) || arma::any(y > 1)) {
        Rcpp::Rcerr << "ERROR: Invalid phenotype values (must be 0 or 1)" << std::endl;
        return y - arma::mean(y); // 返回中心化残差
        
        // 调试输出实际表型值
        Rcpp::Rcout << "Unique y values: " 
                    << arma::unique(y).t() << std::endl;
        
      }
      
      
      
      // 情况3：表型无变异
      if(arma::var(y) < 1e-6) {
        Rcpp::Rcout << "WARNING: Constant phenotype. Returning centered residuals." << std::endl;
        return y - arma::mean(y);
      }
      
      // =============== 数据准备 ===============
      // 构造数据框
      Rcpp::DataFrame df = Rcpp::DataFrame::create(
        Rcpp::Named("response") = Rcpp::wrap(y),
        Rcpp::Named("genotype") = Rcpp::wrap(g),
        Rcpp::_["stringsAsFactors"] = false
      );
      
      // 动态构造公式
      std::string formula_str = "response ~ genotype";
      Rcpp::List cov_list;
      
      // 处理协变量（带唯一性检查）
      std::set<std::string> used_names{"response", "genotype"};
      for(int i = 0; i < covariates.n_cols; ++i) {
        std::string base = "Cov";
        int counter = 1;
        std::string colname;
        do {
          colname = base + std::to_string(counter++);
        } while(used_names.count(colname));
        
        used_names.insert(colname);
        cov_list[colname] = Rcpp::NumericVector(
          covariates.colptr(i),
          covariates.colptr(i) + covariates.n_rows
        );
        
        formula_str += " + " + colname;
        Rcpp::Rcout << "Added covariate: " << colname << std::endl;
      }
      
      // 添加协变量到数据框
      if(covariates.n_cols > 0) {
        df = Rcpp::DataFrame::create(df, cov_list);
      }
      
      // =============== 模型拟合 ===============
      Rcpp::Environment stats = Rcpp::Environment::namespace_env("stats");
      Rcpp::Function glm = stats["glm"];
      Rcpp::List model;
      
      try {
        model = glm(
          Rcpp::Formula(formula_str),
          Rcpp::_["data"]    = df,
          Rcpp::_["family"]  = Rcpp::Function("binomial")()
        );
      } catch(const std::exception& e) {
        Rcpp::Rcerr << "\nGLM Fitting Error: " << e.what() 
                    << "\nFormula: " << formula_str
                    << "\nReturning mean-centered residuals." << std::endl;
        return y - arma::mean(y);
      }
      
      // =============== 残差提取与后处理 ===============
      // 检查模型有效性
      if(!model.inherits("glm")) {
        Rcpp::Rcerr << "ERROR: Model object is invalid" << std::endl;
        return y - arma::mean(y);
      }
      
      // 提取预测值
      Rcpp::NumericVector fitted = model["fitted.values"];
      residuals = Rcpp::as<arma::vec>(Rcpp::wrap(y)) - Rcpp::as<arma::vec>(fitted);
      
      // 后处理：限制精度
      residuals.transform( [](double val) { 
        return std::abs(val) < 1e-8 ? 0.0 : val; 
      });
      
      // =============== 最终校验 ===============
      Rcpp::Rcout << "Residual summary: "
                  << "min=" << arma::min(residuals) << ", "
                  << "max=" << arma::max(residuals) << ", "
                  << "sd=" << arma::stddev(residuals) << std::endl;
      
      return residuals;
      
    } catch(const std::exception& e) {
      Rcpp::Rcerr << "Fatal Error: " << e.what() << std::endl;
      return arma::vec(y.n_elem, arma::fill::zeros);
    } catch(...) {
      Rcpp::Rcerr << "Unknown Error" << std::endl;
      return arma::vec(y.n_elem, arma::fill::zeros);
    }
  }
  
  
  
  
  
};

}

#endif
