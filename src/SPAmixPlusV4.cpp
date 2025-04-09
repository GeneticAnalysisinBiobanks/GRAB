
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "SPAmixPlusV4.hpp"

#include <unordered_map> // 新增头文件
#include <vector>        // 新增头文件


namespace SPAmixPlusV4 {

SPAmixPlusV4Class::SPAmixPlusV4Class(arma::mat t_resid,
                                     arma::mat t_PCs,
                                     int t_N,
                                     double t_SPA_Cutoff,
                                     Rcpp::List t_outlierList,
                                     Rcpp::DataFrame t_sparseGRM,    // 新增参数：稀疏GRM数据
                                     Rcpp::DataFrame t_ResidMat)     // 新增参数：残差矩阵
{
  
  // // ==== 新增1：ID映射系统 ====
  // Rcpp::CharacterVector subjIDs = t_ResidMat["SubjID"];
  // Rcpp::NumericVector indices = t_ResidMat["SubjID_Index"];
  // std::unordered_map<std::string, int> idMap;
  // for(int i=0; i<subjIDs.size(); ++i){
  //   idMap[Rcpp::as<std::string>(subjIDs[i])] = indices[i];
  // }
  // 
  // // ==== 新增2：稀疏矩阵预处理 ====
  // Rcpp::CharacterVector id1 = t_sparseGRM["ID1"];
  // Rcpp::CharacterVector id2 = t_sparseGRM["ID2"];
  // // Rcpp::NumericVector values = t_sparseGRM["Value"];
  // // 修改后（显式转换）：
  // Rcpp::NumericVector values = Rcpp::as<Rcpp::NumericVector>(t_sparseGRM["Value"]);
  // 
  // std::vector<std::tuple<int, int, double>> sparseTriplets;
  // for(int i=0; i<values.size(); ++i){
  //   auto it1 = idMap.find(Rcpp::as<std::string>(id1[i]));
  //   auto it2 = idMap.find(Rcpp::as<std::string>(id2[i]));
  //   if(it1 != idMap.end() && it2 != idMap.end()){
  //     sparseTriplets.emplace_back(it1->second, it2->second, values[i]);
  //   }
  // }
  // m_sparseTriplets = sparseTriplets; // 存储预处理结果
  // 
  // // // ==== 修改3：多表型残差矩阵 ====
  // // Rcpp::NumericMatrix residR = t_ResidMat["Resid"];
  // // m_ResidMat = Rcpp::as<arma::mat>(residR); // 转换为N×M矩阵(N样本数,M表型数)
  // 
  // // 新代码：处理Resid列为列表（每个元素为表型的残差向量）
  // Rcpp::List residList = t_ResidMat["Resid"];
  // int numPheno = residList.size();
  // int numSamples = Rcpp::as<arma::vec>(residList[0]).n_elem;
  // arma::mat residMat(numSamples, numPheno);
  // 
  // // 逐列提取残差向量
  // for(int i = 0; i < numPheno; ++i) {
  //   residMat.col(i) = Rcpp::as<arma::vec>(residList[i]);
  // }
  // m_ResidMat = residMat;
  
  
  // ==== 处理ResidMat ====
  // 获取所有列名
  Rcpp::CharacterVector colNames = t_ResidMat.names();
  
  // 动态识别残差列（例如Resid_1, Resid_2）
  std::vector<std::string> residCols;
  for(int i=0; i<colNames.size(); ++i){
    std::string colName = Rcpp::as<std::string>(colNames[i]);
    if(colName.find("Resid_") != std::string::npos) {
      residCols.push_back(colName);
    }
  }
  
  // 提取残差列并构建矩阵
  int numSamples = t_ResidMat.nrows();
  int numPheno = residCols.size();
  arma::mat residMat(numSamples, numPheno);
  
  for(int i=0; i<numPheno; ++i){
    Rcpp::NumericVector residVec = t_ResidMat[residCols[i]];
    residMat.col(i) = Rcpp::as<arma::vec>(residVec);
  }
  m_ResidMat = residMat;
  
  // ==== 处理稀疏GRM ====
  Rcpp::CharacterVector id1 = t_sparseGRM["ID1"];
  Rcpp::CharacterVector id2 = t_sparseGRM["ID2"];
  Rcpp::NumericVector values = Rcpp::as<Rcpp::NumericVector>(t_sparseGRM["Value"]); // 显式转换
  
  // 构建ID到索引的映射表
  Rcpp::CharacterVector subjIDs = t_ResidMat["SubjID"];
  Rcpp::NumericVector indices = t_ResidMat["SubjID_Index"];
  std::unordered_map<std::string, int> idMap;
  for(int i=0; i<subjIDs.size(); ++i){
    idMap[Rcpp::as<std::string>(subjIDs[i])] = indices[i];
  }
  
  // 构建稀疏矩阵三元组
  std::vector<std::tuple<int, int, double>> sparseTriplets;
  for(int i=0; i<values.size(); ++i){
    auto it1 = idMap.find(Rcpp::as<std::string>(id1[i]));
    auto it2 = idMap.find(Rcpp::as<std::string>(id2[i]));
    if(it1 != idMap.end() && it2 != idMap.end()){
      sparseTriplets.emplace_back(it1->second, it2->second, values[i]);
    }
  }
  m_sparseTriplets = sparseTriplets;
  
  
  
  // ==== 原有初始化逻辑保持不变 ====
  //////////////////////////////////////////////////////////////////
  
  m_resid = t_resid;
  // m_resid2 = pow(m_resid, 2);
  m_Npheno = t_resid.n_cols;
  
  m_pvalVec.resize(m_Npheno);
  m_zScoreVec.resize(m_Npheno);
  
  m_N = t_N;
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_PCs = t_PCs;
  m_outlierList = t_outlierList;
  
  m_sparseGRM = t_sparseGRM; // update by Yuzhuo Ma
  m_ResidMat = t_ResidMat;   // update by Yuzhuo Ma
  
  // m_posOutlier = t_posOutlier;
  // m_posNonOutlier = t_posNonOutlier;
  
  // m_R_outlier = m_resid(m_posOutlier);  // what if no residuals are outlier?? check later (2023-04-23)
  // m_R_nonOutlier = m_resid(m_posNonOutlier);  
  
  arma::mat X = arma::join_horiz(arma::ones(t_N), t_PCs);
  arma::mat X_t = X.t();
  arma::mat XTX = X_t * X;
  arma::mat XTX_inv = arma::inv(XTX);
  arma::vec XTX_inv_diag = XTX_inv.diag();
  
  m_onePlusPCs = X;
  m_sqrt_XTX_inv_diag = arma::sqrt(XTX_inv_diag);
  
  m_diffTime1.zeros(2);
  m_diffTime2.zeros(2);
}
}



