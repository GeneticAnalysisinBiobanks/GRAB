
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "SPAGxEmixPlus.hpp"

#include <unordered_map> // 新增头文件
#include <vector>        // 新增头文件


namespace SPAGxEmixPlus {

SPAGxEmixPlusClass::SPAGxEmixPlusClass(arma::mat t_resid,
                                       arma::mat t_resid_by_E,  // 新增参数：
                                       arma::mat t_PCs,
                                       int t_N,
                                       double t_SPA_Cutoff,
                                       Rcpp::List t_outlierList,
                                       Rcpp::DataFrame t_sparseGRM,    // 新增参数：稀疏GRM数据
                                       Rcpp::DataFrame t_ResidMat,      // 新增参数：残差矩阵
                                       arma::vec t_E,                  // 新增环境因子参数
                                       std::string t_ResidTraitType,
                                       arma::mat t_PhenoMat,
                                       arma::mat t_Covariates
                                       ):     

// // 关键初始化必须在此处（成员初始化列表）
// m_statsEnv(Environment::namespace_env("stats")),
// m_glm(m_statsEnv["glm"])

// 成员初始化列表
m_statsEnv(Environment::namespace_env("stats")),
m_glm(m_statsEnv["glm"])
// m_E(t_E),
// m_resid(t_resid),
// m_resid_by_E(t_resid_by_E),
// m_PCs(t_PCs),
// m_N(t_N),
// m_SPA_Cutoff(t_SPA_Cutoff),
// m_outlierList(t_outlierList),
// m_ResidTraitType(t_ResidTraitType),  // 新增
// m_PhenoMat(t_PhenoMat),             // 新增
// m_Covariates(t_Covariates)         // 新增

{
  
  m_ResidTraitType = t_ResidTraitType;
  m_PhenoMat = t_PhenoMat;
  m_Covariates = t_Covariates;
  
  
  // ==== 新增环境因子处理 ====
  // 检查环境因子维度
  if(t_E.n_elem != static_cast<arma::uword>(t_N)) {  // 修正符号问题
    Rcpp::stop("Environment factor dimension mismatch...");
  }
  m_E = t_E;  // 存储环境因子
  
  // // ==== 新增：处理环境因子交互残差 ====
  // Rcpp::NumericVector residByEnvi = t_ResidByEnviMat["ResidByEnvi"];
  // m_ResidByEnviMat = arma::vec(residByEnvi.begin(), residByEnvi.size());
  
  // ==== 处理ResidMat ====
  Rcpp::IntegerVector subjID_Index = t_ResidMat["SubjID_Index"];
  arma::ivec subjIndices = Rcpp::as<arma::ivec>(subjID_Index);
  
  // 提取Resid_*列（双精度）
  Rcpp::CharacterVector colNames = t_ResidMat.names();
  std::vector<std::string> residCols;
  for(int i=0; i<colNames.size(); ++i){
    std::string colName = Rcpp::as<std::string>(colNames[i]);
    if(colName.find("Resid_") != std::string::npos) {
      residCols.push_back(colName);
    }
  }
  
  int numSamples = t_ResidMat.nrows();
  int numPheno = residCols.size();
  arma::mat residMat(numSamples, numPheno);
  for(int i=0; i<numPheno; ++i){
    Rcpp::NumericVector residVec = Rcpp::as<Rcpp::NumericVector>(t_ResidMat[residCols[i]]);
    residMat.col(i) = arma::vec(residVec.begin(), residVec.size(), false);
  }
  m_ResidMat = residMat;
  
  // ==== 处理sparseGRM ====
  Rcpp::IntegerVector id1_indices = t_sparseGRM["ID1_Index"];
  Rcpp::IntegerVector id2_indices = t_sparseGRM["ID2_Index"];
  Rcpp::NumericVector values = Rcpp::as<Rcpp::NumericVector>(t_sparseGRM["Value"]);
  
  std::vector<std::tuple<int, int, double>> sparseTriplets;
  for(int i=0; i<values.size(); ++i){
    int id1 = id1_indices[i];
    int id2 = id2_indices[i];
    double val = values[i];
    sparseTriplets.emplace_back(id1, id2, val);
  }
  m_sparseTriplets = sparseTriplets;
  
  
  
  // ==== 原有初始化逻辑保持不变 ====
  //////////////////////////////////////////////////////////////////
  
  m_resid = t_resid;
  
  m_resid_by_E = t_resid_by_E;    // update by Yuzhuo Ma
  
  
  // m_resid2 = pow(m_resid, 2);
  m_Npheno = t_resid.n_cols;
  
  m_pvalVec.resize(m_Npheno);
  m_zScoreVec.resize(m_Npheno);
  
  m_N = t_N;
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_PCs = t_PCs;
  m_outlierList = t_outlierList;
  
  // m_sparseGRM = t_sparseGRM; // update by Yuzhuo Ma
  // m_ResidMat = t_ResidMat;   // update by Yuzhuo Ma
  //
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



