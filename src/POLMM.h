#ifndef POLMM_H
#define POLMM_H

// POLMM.h -- Proportional-odds logistic mixed model for ordinal categorical traits

#include <RcppArmadillo.h>
#include <sys/time.h>

namespace POLMM{


// ---- Result structs ----

struct RootResult {
  double root;
  int iter;
  bool converge;
  double K2;
};

struct SaddleResult {
  double pval;
  bool converge;
  arma::vec K1roots;
};

struct ObjP {
  int n, J, p;
  arma::mat XXR_Psi_RX_new;
  arma::mat XR_Psi_R_new;
  arma::vec RymuVec;
  arma::vec RPsiR;
  arma::mat muMat;
  arma::mat iRMat;
};

struct AdjGResult {
  arma::vec adjGVec;
  double Stat;
  double VarW;
};


class POLMMClass {
private:


  arma::vec m_diffTimePOLMM1;
  arma::vec m_diffTimePOLMM2;
  arma::vec m_diffTimePOLMM3;
  arma::vec m_diffTimePOLMM4;
  arma::vec m_diffTimePOLMM5;
  arma::vec m_diffTimePOLMM6;
  arma::vec m_diffTimePOLMM7;
  arma::vec m_diffTimePOLMM8;


  int m_n, m_J, m_p, m_M;


  arma::uvec m_yVec;
  arma::mat m_yMat;
  arma::mat m_Cova, m_CovaMat;


  arma::vec m_beta, m_bVec, m_eps;
  double m_tau;


  unsigned int m_iter, m_maxiter, m_maxiterPCG, m_maxiterEps, m_tracenrun, m_nSNPsVarRatio, m_grainSize;
  double m_tolBeta, m_tolTau, m_tolPCG, m_tolEps, m_CVcutoff, m_minMafVarRatio, m_maxMissingVarRatio, m_memoryChunk, m_minMafGRM, m_maxMissingGRM;
  bool m_LOCO, m_showInfo, m_printPCGInfo, m_flagSparseGRM;
  int m_seed;


  arma::mat m_WMat, m_muMat, m_mMat, m_nuMat, m_iRMat, m_YMat, m_iSigma_CovaMat, m_iSigmaX_XSigmaX;
  arma::vec m_eta, m_iSigma_YVec, m_iSigma_VPYVec;

  arma::mat m_TraceRandMat, m_V_TRM, m_iSigma_V_TRM;

  arma::mat m_XXR_Psi_RX;
  arma::mat m_XR_Psi_R;
  arma::vec m_RymuVec;
  arma::vec m_RPsiR;

  arma::cube m_InvBlockDiagSigma;


  arma::sp_mat m_SparseGRM;
  arma::sp_mat m_ZMat_sp;
  arma::sp_mat m_SigmaMat_sp;


  arma::Mat<uint8_t> m_SeqMat;

  void setArray() {
    m_WMat.zeros(m_n, m_J);
    m_muMat.zeros(m_n, m_J);
    m_mMat.zeros(m_n, m_J);
    m_nuMat.zeros(m_n, m_J);
    m_iRMat.zeros(m_n, m_J-1);
    m_YMat.zeros(m_n, m_J-1);
    m_iSigma_CovaMat.zeros(m_n * (m_J-1), m_p);

    m_eta.zeros(m_n);
    m_iSigma_YVec.zeros(m_n * (m_J-1));
    m_iSigma_VPYVec.zeros(m_n * (m_J-1));

    m_TraceRandMat.zeros(m_n * (m_J-1), m_tracenrun);
    m_V_TRM.zeros(m_n * (m_J-1), m_tracenrun);
    m_iSigma_V_TRM.zeros(m_n * (m_J-1), m_tracenrun);
  }


  arma::mat getyMat(arma::uvec yVec) {
    unsigned int n = yVec.size();
    unsigned int J = arma::max(yVec) + 1;

    arma::mat yMat(n, J, arma::fill::zeros);
    for (unsigned int i = 0; i < n; i++)
      yMat(i, yVec.at(i)) = 1;

    return yMat;
  };

  void setPOLMMInner(arma::mat Cova,
                     arma::uvec yVec,
                     arma::vec beta,
                     arma::vec bVec,
                     arma::vec eps,
                     double tau);

  arma::sp_mat setZMat_sp() {
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


  arma::vec ZMat(arma::vec xVec) {
    arma::vec y1Vec(m_n, arma::fill::zeros);
    int index = 0;
    for (int i = 0; i < m_n; i ++){
      for (int j = 0; j < m_J-1; j ++){
        y1Vec(i) += xVec(index);
        index++;
      }
    }
    return y1Vec;
  }


  arma::vec tZMat(arma::vec xVec) {
    arma::vec y1Vec(m_n * (m_J-1));
    int index = 0;
    for (int i = 0; i < m_n; i ++){
      for (int j = 0; j < m_J-1; j ++){
        y1Vec(index) = xVec(i);
        index++;
      }
    }
    return y1Vec;
  }


  void getTraceRandMat();

  arma::vec getKinbVecPOLMM(arma::vec bVec, std::string excludeChr);


public:

  arma::vec getTestTime1(){return m_diffTimePOLMM1;}
  arma::vec getTestTime2(){return m_diffTimePOLMM2;}
  arma::vec getTestTime3(){return m_diffTimePOLMM3;}
  arma::vec getTestTime4(){return m_diffTimePOLMM4;}
  arma::vec getTestTime5(){return m_diffTimePOLMM5;}
  arma::vec getTestTime6(){return m_diffTimePOLMM6;}
  arma::vec getTestTime7(){return m_diffTimePOLMM7;}
  arma::vec getTestTime8(){return m_diffTimePOLMM8;}

  void updateMats();
  void updateTau();
  void updatePara(std::string excludechr);
  void updateEps();
  void updateEpsOneStep();

  void getPCGofSigmaAndVector(arma::vec y1Vec,
                              arma::vec& xVec,
                              std::string excludechr);
  arma::mat getSigmaxMat(arma::mat xMat,
                         std::string excludechr);
  arma::mat solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,
                                 arma::mat& xMat);

  void getPCGofSigmaAndCovaMat(arma::mat xMat,
                               arma::mat& iSigma_xMat,
                               std::string excludechr);
  double getVarP(arma::vec adjGVec,
                 std::string excludechr);

  POLMMClass(arma::mat muMat,
             arma::mat iRMat,
             arma::mat Cova,
             arma::uvec yVec,
             arma::sp_mat SparseGRM,
             double tau,
             bool printPCGInfo,
             double tolPCG,
             int maxiterPCG,
             double varRatio,
             double SPA_cutoff,
             bool flagSparseGRM);

  void getMarkerPval(arma::vec GVec,
                     double& Beta,
                     double& seBeta,
                     double& pval,
                     double altFreq,
                     double& zScore);

  void getRegionPVec(arma::vec GVec,
                     double& Stat,
                     double& Beta,
                     double& seBeta,
                     double& pval0,
                     double& pval1,
                     arma::vec& P1Vec,
                     arma::vec& P2Vec);

  double m_varRatio, m_SPA_Cutoff;

  arma::vec getadjGFast(arma::vec GVec);
  double getStatFast(arma::vec adjGVec);
  arma::vec get_ZPZ_adjGVec(arma::vec adjGVec);
  void getPCGofSigmaAndCovaMat(arma::mat xMat,
                               arma::mat& iSigma_xMat);
  void getPCGofSigmaAndVector(arma::vec y1Vec,
                              arma::vec& xVec);
  arma::mat getSigmaxMat(arma::mat& xMat);

  arma::mat getiPsixMat(arma::mat xMat);
  arma::mat getPsixMat(arma::mat xMat);

  arma::cube getInvBlockDiagSigma();
  arma::mat solverBlockDiagSigma(arma::mat& xMat);

  void setRPsiR();
  arma::vec getVarWVec(arma::vec adjGVec);

  SaddleResult MAIN_SPA(double Stat,
                        arma::vec adjGVec,
                        arma::vec K1roots,
                        double VarP,
                        double VarW,
                        double Ratio0,
                        arma::uvec posG1);

  double MAIN_ER(arma::vec GVec,
                 arma::uvec posG1);

  void setSeqMat(int NonZero_cutoff);

};


// ---- Free functions used by POLMMClass ----
arma::mat getCovaMat(arma::mat Cova, unsigned int J);

arma::mat getPsixMat(arma::mat xMat,
                     arma::mat muMat);

arma::mat sumCols(arma::mat xMat,
                  int J);

double calCV(arma::vec xVec);

arma::vec getRPsiR(arma::mat muMat,
                   arma::mat iRMat,
                   int n, int J, int p);

double K0(double x,
          arma::mat muMat,
          arma::mat cMat,
          double m1);

arma::vec K12(double x,
              arma::mat muMat,
              arma::mat cMat,
              double m1);

RootResult fastgetroot_K1(double Stat,
                          double initX,
                          double Ratio0,
                          arma::mat muMat,
                          arma::mat cMat,
                          double m1);

double fastGet_Saddle_Prob(double Stat,
                           double zeta,
                           double K2,
                           double Ratio0,
                           arma::mat muMat,
                           arma::mat cMat,
                           double m1,
                           bool lowerTail);

SaddleResult fastSaddle_Prob(double Stat,
                           double VarP,
                           double VarW,
                           double Ratio0,
                           arma::vec K1roots,
                           arma::vec adjGVec1,
                           arma::mat muMat1,
                           arma::mat iRMat1);

double getPvalER(arma::uvec yVec,
                 arma::vec GVec,
                 arma::mat muMat,
                 arma::mat iRMat,
                 arma::Mat<uint8_t> SeqMat);

arma::vec getStatVec(arma::Mat<uint8_t> SeqMat,
                     arma::vec GVec,
                     arma::mat muMat,
                     arma::mat iRMat);

double getProbOne(arma::Col<uint8_t> SeqVec,
                  arma::mat muMat);

double getProb(arma::Mat<uint8_t> SeqMat,
               arma::mat muMat);

arma::vec convert1(arma::mat xMat,
                   int n, int J);

arma::mat convert2(arma::vec xVec,
                   int n, int J);

ObjP getobjP(arma::mat Cova,
             arma::mat yMat,
             arma::mat muMat,
             arma::mat iRMat);

arma::vec getadjGFast(arma::vec GVec,
                      arma::mat XXR_Psi_RX_new,
                      arma::mat XR_Psi_R_new,
                      int n, int p);

double getStatFast(arma::vec GVec,
                   arma::vec RymuVec);

double getVarWFast(arma::vec adjGVec,
                   arma::vec RPsiRVec);

AdjGResult outputadjGFast(arma::vec GVec,
                          const ObjP& objP);

arma::vec getRowSums(arma::mat xMat);

arma::vec getTime();

void printTime(arma::vec t1, arma::vec t2, std::string message);

arma::vec nb (unsigned int n);

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat);

arma::vec Vec2LongVec(arma::vec xVec, int n, int J);

arma::vec LongVec2Vec(arma::vec xVec, int n, int J);

arma::mat Vec2Mat(arma::vec xVec, int n, int J);

arma::vec Mat2Vec(arma::mat xMat, int n, int J);

}

#endif
