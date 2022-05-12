
#ifndef POLMM_HPP
#define POLMM_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "DenseGRM.hpp"
#include "BGEN.hpp"

namespace POLMM{

class POLMMClass
{
private:
  
  ////////////////////// -------------------- members ---------------------------------- //////////////////////
  
  arma::vec m_diffTimePOLMM1;
  arma::vec m_diffTimePOLMM2;
  arma::vec m_diffTimePOLMM3;
  arma::vec m_diffTimePOLMM4;
  arma::vec m_diffTimePOLMM5;
  arma::vec m_diffTimePOLMM6;
  arma::vec m_diffTimePOLMM7;
  arma::vec m_diffTimePOLMM8;
  
  // MAIN Dimensions: sample size, number of categories, number of covariates
  int m_n, m_J, m_p, m_M;
  
  // Input data:
  // m_yVec: a vector (m_n x 1): ordinal categorical data: 0, 1, 2, ..., m_J-1
  // m_yMat: a matrix (m_n x m_J): the equivalent representation of m_yVec
  // m_Cova: a matrix (m_n x m_p): the covariate matrix
  // m_CovaMat: a matrix (m_n(m_J-1) x m_p): the equivalent representation of m_Cova 
  
  arma::uvec m_yVec;
  arma::mat m_yMat;
  arma::mat m_Cova, m_CovaMat;
  
  // Parameters in the model
  arma::vec m_beta, m_bVec, m_eps;
  double m_tau;
  
  // Control parameters
  unsigned int m_iter, m_maxiter, m_maxiterPCG, m_maxiterEps, m_tracenrun, m_seed, m_nSNPsVarRatio, m_grainSize; 
  double m_tolBeta, m_tolTau, m_tolPCG, m_tolEps, m_CVcutoff, m_minMafVarRatio, m_maxMissingVarRatio, m_memoryChunk, m_minMafGRM, m_maxMissingGRM; 
  bool m_LOCO, m_showInfo, m_printPCGInfo, m_flagSparseGRM;
  Rcpp::List m_LOCOList;
  
  // working vectors/matrix
  arma::mat m_WMat, m_muMat, m_mMat, m_nuMat, m_iRMat, m_YMat, m_iSigma_CovaMat, m_iSigmaX_XSigmaX;
  arma::vec m_eta, m_iSigma_YVec, m_iSigma_VPYVec;
  
  arma::mat m_TraceRandMat, m_V_TRM, m_iSigma_V_TRM;
  
  arma::mat m_XXR_Psi_RX;  // XXR_Psi_RX ( n x p )
  arma::mat m_XR_Psi_R;    // XR_Psi_R ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject
  arma::vec m_RymuVec;     // n x 1: row sum of the n x (J-1) matrix R %*% (yMat - muMat)
  arma::vec m_RPsiR;
  
  PLINK::PlinkClass* m_ptrPlinkObj;
  BGEN::BgenClass* m_ptrBgenObj;
  
  DenseGRM::DenseGRMClass* m_ptrDenseGRMObj;
  
  arma::cube m_InvBlockDiagSigma;
  
  // SparseGRM 
  // Rcpp::List m_SparseGRM;
  // arma::sp_mat m_SparseGRM_all;
  arma::sp_mat m_SparseGRM;
  arma::sp_mat m_ZMat_sp;
  arma::sp_mat m_SigmaMat_sp;
  
  // for Efficient Resampling (ER)
  arma::Mat<uint8_t> m_SeqMat;
  
  void setControlList(Rcpp::List t_controlList)
  {
    m_memoryChunk = t_controlList["memoryChunk"];
    m_minMafGRM = t_controlList["minMafGRM"];
    m_maxMissingGRM = t_controlList["maxMissingGRM"];
    m_maxiter = t_controlList["maxiter"];
    m_maxiterPCG = t_controlList["maxiterPCG"]; 
    m_maxiterEps = t_controlList["maxiterEps"];
    m_tolBeta = t_controlList["tolBeta"]; 
    m_tolTau = t_controlList["tolTau"]; 
    m_tolPCG = t_controlList["tolPCG"]; 
    m_tolEps = t_controlList["tolEps"];
    m_tracenrun = t_controlList["tracenrun"]; 
    m_seed = t_controlList["seed"];
    m_minMafVarRatio = t_controlList["minMafVarRatio"];
    m_maxMissingVarRatio = t_controlList["maxMissingVarRatio"];
    m_nSNPsVarRatio = t_controlList["nSNPsVarRatio"];
    m_CVcutoff = t_controlList["CVcutoff"];
    m_LOCO = t_controlList["LOCO"];
    m_grainSize = t_controlList["grainSize"];
    m_showInfo = t_controlList["showInfo"];
  }
  
  void setArray()
  {
    m_WMat.zeros(m_n, m_J);    
    m_muMat.zeros(m_n, m_J);   
    m_mMat.zeros(m_n, m_J);
    m_nuMat.zeros(m_n, m_J);
    m_iRMat.zeros(m_n, m_J-1);
    m_YMat.zeros(m_n, m_J-1);
    m_iSigma_CovaMat.zeros(m_n * (m_J-1), m_p);
    //
    m_eta.zeros(m_n);
    m_iSigma_YVec.zeros(m_n * (m_J-1));
    m_iSigma_VPYVec.zeros(m_n * (m_J-1));
    //
    m_TraceRandMat.zeros(m_n * (m_J-1), m_tracenrun);
    m_V_TRM.zeros(m_n * (m_J-1), m_tracenrun);
    m_iSigma_V_TRM.zeros(m_n * (m_J-1), m_tracenrun);
  }
  
  // yMat: matrix with dim of n x J
  arma::mat getyMat(arma::uvec t_yVec) 
  {
    unsigned int n = t_yVec.size();
    unsigned int J = arma::max(t_yVec) + 1;
    
    arma::mat yMat(n, J, arma::fill::zeros);
    for(unsigned int i = 0; i < n; i++)
      yMat(i, t_yVec.at(i)) = 1;
    
    return yMat;
  };
  
  void setPOLMMInner(arma::mat t_Cova,
                     arma::uvec t_yVec,     // should be from 0 to J-1
                     arma::vec t_beta,
                     arma::vec t_bVec,
                     arma::vec t_eps,           // 
                     double t_tau);
  
  arma::sp_mat setZMat_sp()
  {
    arma::umat locations(2, m_n * (m_J-1));
    arma::urowvec l1 = arma::linspace<arma::urowvec>(0, m_n * (m_J-1) - 1, m_n * (m_J-1));
    arma::urowvec l2 = l1 / (m_J-1);
    locations.row(0) = l1;
    locations.row(1) = l2;
    arma::vec values(m_n * (m_J-1));
    values.fill(1);
    arma::sp_mat ZMat_sp(locations, values);
    return ZMat_sp;
  }
  
  // sum up each (J-1) elements: n(J-1) x 1 -> n x 1
  arma::vec ZMat(arma::vec t_xVec)
  {
    arma::vec y1Vec(m_n, arma::fill::zeros);
    int index = 0;
    for(int i = 0; i < m_n; i ++){
      for(int j = 0; j < m_J-1; j ++){
        y1Vec(i) += t_xVec(index);
        index++;
      }
    }
    return(y1Vec);
  }
  
  // duplicate each element for (J-1) times: n x 1 -> n(J-1) x 1 
  arma::vec tZMat(arma::vec t_xVec)
  {
    arma::vec y1Vec(m_n * (m_J-1));
    int index = 0;
    for(int i = 0; i < m_n; i ++){
      for(int j = 0; j < m_J-1; j ++){
        y1Vec(index) = t_xVec(i);
        index++;
      }
    }
    return(y1Vec);
  }
  
  // set up m_TraceRandMat (TRM) and m_V_TRM, only used once at ()
  void getTraceRandMat();
  
  arma::vec getKinbVecPOLMM(arma::vec t_bVec, std::string t_excludeChr);
  
  ////////////////////// -------------------- functions ---------------------------------- //////////////////////
  
public:
  
  arma::vec getTestTime1(){return m_diffTimePOLMM1;}
  arma::vec getTestTime2(){return m_diffTimePOLMM2;}
  arma::vec getTestTime3(){return m_diffTimePOLMM3;}
  arma::vec getTestTime4(){return m_diffTimePOLMM4;}
  arma::vec getTestTime5(){return m_diffTimePOLMM5;}
  arma::vec getTestTime6(){return m_diffTimePOLMM6;}
  arma::vec getTestTime7(){return m_diffTimePOLMM7;}
  arma::vec getTestTime8(){return m_diffTimePOLMM8;}
  
  void fitPOLMM();
  void estVarRatio(arma::mat GenoMat);
  void updateMats();
  void updateParaConv(std::string t_excludechr);
  void updateTau();
  void updatePara(std::string t_excludechr);
  void updateEps();
  void updateEpsOneStep();
  
  Rcpp::List getPOLMM();
  
  arma::mat getVarRatio(arma::mat t_GMatRatio, std::string t_excludechr);
  arma::rowvec getVarOneSNP(arma::vec GVec,
                                        std::string excludechr,
                                        Rcpp::List objP);
  void getPCGofSigmaAndVector(arma::vec t_y1Vec,    // vector with length of n(J-1)
                              arma::vec& t_xVec,    // vector with length of n(J-1)
                              std::string t_excludechr);
  arma::mat getSigmaxMat(arma::mat t_xMat,   // matrix: n x (J-1) 
                         std::string t_excludechr);
  arma::mat solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,   // (J-1) x (J-1) x n
                                 arma::mat& xMat);                 // n x (J-1)
  
  void getPCGofSigmaAndCovaMat(arma::mat t_xMat,              // matrix with dim of n(J-1) x p
                               arma::mat& t_iSigma_xMat,      // matrix with dim of n(J-1) x p
                               std::string t_excludechr);
  double getVarP(arma::vec t_adjGVec,
                             std::string t_excludechr);
  
  POLMMClass(bool t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
             DenseGRM::DenseGRMClass* t_ptrDenseGRMObj,
             PLINK::PlinkClass* t_ptrPlinkObj,
             BGEN::BgenClass* t_ptrBgenObj,
             arma::mat t_Cova,
             arma::uvec t_yVec,     // should be from 0 to J-1
             arma::vec t_beta,
             arma::vec t_bVec,
             arma::vec t_eps,           // 
             double t_tau,
             arma::sp_mat t_SparseGRM,    // results of function getKinMatList()
             Rcpp::List t_controlList)
  {
    setPOLMMObj(t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
                t_ptrDenseGRMObj,
                t_ptrPlinkObj,
                t_ptrBgenObj,
                t_Cova,
                t_yVec,     // should be from 0 to J-1
                t_beta,
                t_bVec,
                t_eps,           // 
                t_tau,
                t_SparseGRM,    // results of function getKinMatList()
                t_controlList); 
  }
  
  POLMMClass(arma::mat t_muMat,
             arma::mat t_iRMat,
             arma::mat t_Cova,
             arma::uvec t_yVec,
             arma::sp_mat t_SparseGRM,
             double t_tau,
             bool t_printPCGInfo,
             double t_tolPCG,
             int t_maxiterPCG,
             double t_varRatio, 
             double t_SPA_cutoff,
             bool t_flagSparseGRM);
  
  
  void setPOLMMObj(bool t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
                   DenseGRM::DenseGRMClass* t_ptrDenseGRMObj,
                   PLINK::PlinkClass* t_ptrPlinkObj,
                   BGEN::BgenClass* t_ptrBgenObj,
                   arma::mat t_Cova,
                   arma::uvec t_yVec,     // should be from 1 to J
                   arma::vec t_beta,
                   arma::vec t_bVec,
                   arma::vec t_eps,           // 
                   double t_tau,
                   arma::sp_mat t_SparseGRM,    // results of function getKinMatList()
                   Rcpp::List t_controlList);
  
  void getMarkerPval(arma::vec t_GVec, 
                     double& t_Beta, 
                     double& t_seBeta, 
                     double& t_pval, 
                     double t_altFreq,
                     double& t_zScore);
  
  void getRegionPVec(arma::vec t_GVec, 
                     double& t_Stat,
                     double& t_Beta, 
                     double& t_seBeta, 
                     double& t_pval0, 
                     double& t_pval1,
                     arma::vec& t_P1Vec, 
                     arma::vec& t_P2Vec);
  
  double m_varRatio, m_SPA_Cutoff;
  
  arma::vec getadjGFast(arma::vec t_GVec);
  double getStatFast(arma::vec t_adjGVec);
  arma::vec get_ZPZ_adjGVec(arma::vec t_adjGVec);
  void getPCGofSigmaAndCovaMat(arma::mat t_xMat,             // matrix with dim of n(J-1) x p
                               arma::mat& t_iSigma_xMat);    // matrix with dim of n(J-1) x p
  void getPCGofSigmaAndVector(arma::vec t_y1Vec,    // vector with length of n(J-1)
                              arma::vec& t_xVec);    // vector with length of n(J-1)
  arma::mat getSigmaxMat(arma::mat& t_xMat);
  
  arma::mat getiPsixMat(arma::mat t_xMat);
  arma::mat getPsixMat(arma::mat t_xMat);
  
  arma::cube getInvBlockDiagSigma();
  arma::mat solverBlockDiagSigma(arma::mat& t_xMat);
  
  void setRPsiR();
  arma::vec getVarWVec(arma::vec adjGVec);
  
  Rcpp::List MAIN_SPA(double t_Stat,
                      arma::vec t_adjGVec,
                      arma::vec t_K1roots,
                      double t_VarP,
                      double t_VarW,
                      double t_Ratio0,
                      arma::uvec t_posG1);
  
  double MAIN_ER(arma::vec t_GVec,
                 arma::uvec t_posG1);
  
  void setSeqMat(int t_NonZero_cutoff);
  
};

// duplicate each row for (J-1) times: n x p -> n(J-1) x p
arma::mat getCovaMat(arma::mat t_Cova, unsigned int t_J);

// outMat = PsiMat %*% xMat, PsiMat is determined by muMat
arma::mat getPsixMat(arma::mat t_xMat,    // matrix: n x (J-1)
                     arma::mat t_muMat);   // matrix: n x J

// sum each (J-1) cols to 1 col: p x n(J-1) -> p x n (OR) p x (J-1) -> p x 1
arma::mat sumCols(arma::mat t_xMat,
                  int J);

double calCV(arma::vec t_xVec);

// get RPsiP: (J-1) x (J-1) x n 
// Only used in getVarWFast(): 
arma::vec getRPsiR(arma::mat t_muMat,
                   arma::mat t_iRMat,
                   int t_n, int t_J, int t_p);

double K0(double t_x,
          arma::mat t_muMat,     // N x (J-1)
          arma::mat t_cMat,      // N x (J-1)
          double t_m1);           // sum(muMat * cMat)

arma::vec K12(double t_x,
              arma::mat t_muMat,
              arma::mat t_cMat,
              double t_m1);

Rcpp::List fastgetroot_K1(double t_Stat,
                          double t_initX,
                          double t_Ratio0,
                          arma::mat t_muMat,
                          arma::mat t_cMat,
                          double t_m1);

double fastGet_Saddle_Prob(double t_Stat,
                           double t_zeta,
                           double t_K2,
                           double t_Ratio0,
                           arma::mat t_muMat,
                           arma::mat t_cMat,
                           double t_m1,          // sum(muMat * cMat)
                           bool t_lowerTail);

// add partial normal approximation to speed up the SPA
Rcpp::List fastSaddle_Prob(double t_Stat,
                           double t_VarP,
                           double t_VarW,
                           double t_Ratio0,      // Ratio of variance (G==0)
                           arma::vec t_K1roots,
                           arma::vec t_adjGVec1, // N1 x 1, where N1 is length(G!=0)
                           arma::mat t_muMat1,   // N1 x (J-1)
                           arma::mat t_iRMat1);  // N1 x (J-1)

double getPvalER(arma::uvec t_yVec,     // N1 x 1 vector, from 0 to J-1
                 arma::vec t_GVec,      // N1 x 1 vector,
                 arma::mat t_muMat,     // N1 x J matrix,
                 arma::mat t_iRMat,     // N1 x (J-1) matrix
                 arma::Mat<uint8_t> t_SeqMat);   // N1 x nER

arma::vec getStatVec(arma::Mat<uint8_t> t_SeqMat,   // n x J^n matrix
                     arma::vec t_GVec,         // n x 1 vector, where n is number of subjects with Geno != 0
                     arma::mat t_muMat,        // n x J matrix, where n is number of subjects with Geno != 0
                     arma::mat t_iRMat);       // n x (J-1) matrix

double getProbOne(arma::Col<uint8_t> t_SeqVec,  // n x 1
                  arma::mat t_muMat);           // n x J

double getProb(arma::Mat<uint8_t> t_SeqMat,  // n x m matrix, where m \leq J^n is the number of resampling with abs(stat) > stat_obs
               arma::mat t_muMat);           // n x J matrix

// convert: n x (J-1) -> n(J-1) x 1
arma::vec convert1(arma::mat xMat, // matrix: n x (J-1)
                   int n, int J);

// convert: n(J-1) x 1 -> n x (J-1)
arma::mat convert2(arma::vec xVec, // n(J-1) x 1 
                   int n, int J);

// get a list for p value calculation in step 2
Rcpp::List getobjP(arma::mat t_Cova,     // matrix: n x p
                   arma::mat t_yMat,
                   arma::mat t_muMat,    // matrix: n x J
                   arma::mat t_iRMat);    // matrix: n x (J-1)

arma::vec getadjGFast(arma::vec GVec,
                      arma::mat XXR_Psi_RX_new,   // XXR_Psi_RX_new ( n x p )
                      arma::mat XR_Psi_R_new,     // XR_Psi_R_new ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject 
                      int n, int p);

double getStatFast(arma::vec GVec,         // n x 1
                   arma::vec RymuVec);      // n x 1: row sum of the n x (J-1) matrix R %*% (yMat - muMat)

double getVarWFast(arma::vec adjGVec,  // n x 1
                   arma::vec RPsiRVec); // n x 1


Rcpp::List outputadjGFast(arma::vec GVec,
                          Rcpp::List objP);

// sum up each row: n1 x n2 matrix -> n1 x 1 vector
arma::vec getRowSums(arma::mat t_xMat);

}

#endif
