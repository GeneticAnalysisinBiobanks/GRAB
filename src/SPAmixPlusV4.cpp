
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>             // 改用vector提升性能
#include "SPAmixPlusV4.hpp"

using namespace arma;         // update by Yuzhuo Ma
using namespace std;          // update by Yuzhuo Ma


namespace SPAmixPlusV4 {

// 修改1：构造函数改用整数ID处理
SPAmixPlusV4Class::SPAmixPlusV4Class(arma::mat t_resid,
                         arma::mat t_PCs,
                         int t_N,
                         double t_SPA_Cutoff,
                         Rcpp::List t_outlierList,
                         
                         const arma::mat& t_sparseGRM,  // update by Yuzhuo Ma      // [M x 3] (整数ID1, 整数ID2, Value)
                         const arma::mat& t_ResidMat)   // update by Yuzhuo Ma   // [N x (1 + nPheno)] (整数SubjID, Resid_1, ...)
// arma::uvec t_posOutlier,
// arma::uvec t_posNonOutlier)
{
  
  
  ////////////////////////////////////////////////////////////////////
  
  m_resid = t_resid;
  // m_resid2 = pow(m_resid, 2);
  m_Npheno = t_resid.n_cols;
  
  m_pvalVec.resize(m_Npheno);
  m_zScoreVec.resize(m_Npheno);
  
  m_N = t_N;
  m_SPA_Cutoff = t_SPA_Cutoff;
  m_PCs = t_PCs;
  m_outlierList = t_outlierList;
  
 
  
  // 修改2：直接提取整数ID并构建索引向量
  // ----------------------------------------------------------
  // 处理ResidMat的SubjID（第一列）
  arma::vec subjID_vec = t_ResidMat.col(0);
  arma::uvec subjID = arma::conv_to<arma::uvec>::from(subjID_vec);
  
  // 构建索引向量（O(1)访问）
  int max_id = arma::max(subjID);
  m_SubjID_to_Index.resize(max_id + 1, -1); // -1表示无效ID
  for (uword i = 0; i < subjID.n_elem; ++i) {
    int id = subjID(i);
    m_SubjID_to_Index[id] = i;
  }
  
  // 存储残差数据（去除SubjID列）
  m_ResidMat_Residuals = t_ResidMat.cols(1, t_ResidMat.n_cols-1);
  
  // 处理sparseGRM（直接存储为整数矩阵）
  // ----------------------------------------------------------
  m_sparseGRM_ID1 = arma::conv_to<arma::uvec>::from(t_sparseGRM.col(0));
  m_sparseGRM_ID2 = arma::conv_to<arma::uvec>::from(t_sparseGRM.col(1));
  m_sparseGRM_Value = t_sparseGRM.col(2);

  
  
  
  
  
  
  
  
  
  
  
  
  // 其他原有初始化逻辑保持不变
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








// // 计算基于稀疏GRM的方差（整数ID优化版）
// double SPAmixPlusV4Class::ComputeSparseGRMVariance(const arma::vec& MAF_est) {
//   // 1. 调整残差
//   arma::vec R = m_ResidMat_Residuals.col(0);  // 使用第一个表型的残差
//   arma::vec R_new = R % sqrt(2 * MAF_est % (1 - MAF_est));
//   
//   // 2. 计算协方差项
//   double cov_sum = 0.0;
//   for (uword i = 0; i < m_sparseGRM_Value.n_elem; ++i) {
//     int id1 = m_sparseGRM_ID1(i);
//     int id2 = m_sparseGRM_ID2(i);
//     double value = m_sparseGRM_Value(i);
//     
//     // 直接通过索引向量访问（无哈希开销）
//     if (id1 < m_SubjID_to_Index.size() && 
//         id2 < m_SubjID_to_Index.size()) {
//       int idx1 = m_SubjID_to_Index[id1];
//       int idx2 = m_SubjID_to_Index[id2];
//       if (idx1 != -1 && idx2 != -1) {
//         cov_sum += value * R_new(idx1) * R_new(idx2);
//       }
//     }
//   }
//   
//   // 3. 计算最终方差
//   arma::vec g_var_est = 2 * MAF_est % (1 - MAF_est);
//   return 2 * cov_sum - sum(R % R % g_var_est);
// }






} // namespace SPAmixPlusV4


