// POLMM.cpp -- POLMMClass and free-function implementations

#include <RcppArmadillo.h>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <chrono>

#include "mtPOLMM.h"

namespace POLMM {

POLMMClass::POLMMClass(arma::mat muMat,
                       arma::mat iRMat,
                       arma::mat Cova,
                       arma::uvec yVec,
                       arma::sp_mat SparseGRM,
                       double tau,
                       bool printPCGInfo,
                       double tolPCG,
                       int maxiterPCG,
                       double varRatio,
                       double SPA_Cutoff,
                       bool flagSparseGRM) {
  m_diffTimePOLMM1.zeros(2);
  m_diffTimePOLMM2.zeros(2);
  m_diffTimePOLMM3.zeros(2);
  m_diffTimePOLMM4.zeros(2);
  m_diffTimePOLMM5.zeros(2);
  m_diffTimePOLMM6.zeros(2);
  m_diffTimePOLMM7.zeros(2);
  m_diffTimePOLMM8.zeros(2);

  m_muMat = muMat;
  m_iRMat = iRMat;
  m_Cova = Cova;
  m_varRatio = varRatio;
  m_SPA_Cutoff = SPA_Cutoff;
  m_flagSparseGRM = flagSparseGRM;

  m_SparseGRM = SparseGRM;
  m_printPCGInfo = printPCGInfo;
  m_tolPCG = tolPCG;
  m_maxiterPCG = maxiterPCG;

  m_n = m_muMat.n_rows;
  m_J = m_muMat.n_cols;
  m_p = m_Cova.n_cols;

  m_CovaMat = getCovaMat(m_Cova, m_J);

  m_tau = tau;

  if (m_flagSparseGRM == true){
    m_InvBlockDiagSigma = getInvBlockDiagSigma();
  }


  arma::mat XR_Psi_R(m_p, m_n * (m_J-1));
  for (int k = 0; k < m_p; k++){
    arma::mat xMat = Vec2Mat(m_CovaMat.col(k), m_n, m_J);
    arma::vec temp = Mat2Vec(getPsixMat(xMat / m_iRMat) / m_iRMat, m_n, m_J);
    XR_Psi_R.row(k) = temp.t();
  }

  m_XXR_Psi_RX = m_Cova * inv(XR_Psi_R * m_CovaMat);


  m_XR_Psi_R = sumCols(XR_Psi_R, m_J);

  m_yVec = yVec;
  arma::mat yMat(m_n, m_J, arma::fill::zeros);
  for (int i = 0; i < m_n; i++)
    yMat(i, m_yVec(i)) = 1;

  arma::mat ymuMat = yMat - m_muMat;
  arma::mat RymuMat = ymuMat.cols(0, m_J-2) / iRMat;
  m_RymuVec = sumCols(RymuMat, m_J);

  if (m_flagSparseGRM == true){
    arma::mat iSigma_CovaMat(m_n * (m_J-1), m_p);
    getPCGofSigmaAndCovaMat(m_CovaMat, iSigma_CovaMat);
    arma::mat XSigmaX = inv(m_CovaMat.t() * iSigma_CovaMat);
    m_iSigmaX_XSigmaX = iSigma_CovaMat * XSigmaX;
  }

  setRPsiR();
}


void set_seed (unsigned int seed) {
  arma::arma_rng::set_seed(seed);
}

void POLMMClass::getTraceRandMat() {
  arma::vec t1  = getTime();
  for (unsigned int itrace = 0; itrace < m_tracenrun; itrace++) {
    arma::vec uVec = nb(m_n * (m_J-1));
    uVec = uVec * 2 - 1;
    m_TraceRandMat.col(itrace) = uVec;
    arma::vec ZuVec = ZMat(uVec);


    arma::vec tempVec = getKinbVecPOLMM(ZuVec, "none");
    m_V_TRM.col(itrace) = tZMat(tempVec);
  }

  arma::vec t2  = getTime();
  std::string info = "calculate " + std::to_string(m_tracenrun) + " genKinbVec()";
  printTime(t1, t2, info);
}

void POLMMClass::setPOLMMInner(arma::mat Cova,
                               arma::uvec yVec,
                               arma::vec beta,
                               arma::vec bVec,
                               arma::vec eps,
                               double tau) {
  m_n = Cova.n_rows;
  m_p = Cova.n_cols;
  m_J = arma::max(yVec) + 1;

  m_CovaMat = getCovaMat(Cova, m_J);
  m_yVec = yVec;
  m_yMat = getyMat(yVec);

  m_Cova = Cova;
  m_beta = beta;
  m_bVec = bVec;
  m_eps = eps;
  m_tau = tau;

  setArray();


  if (m_seed != -1){
    set_seed(static_cast<unsigned int>(m_seed));
  }
}


void POLMMClass::getMarkerPval(arma::vec GVec,
                               double& Beta,
                               double& seBeta,
                               double& pval,
                               double altFreq,
                               double& zScore) {
  arma::vec adjGVec = getadjGFast(GVec);

  double statVal = getStatFast(adjGVec);
  arma::vec VarWVec = getVarWVec(adjGVec);
  double VarW = sum(VarWVec);
  double VarS = VarW * m_varRatio;

  double StdStat = std::abs(statVal) / sqrt(VarS);
  double pvalNorm = 2 * arma::normcdf(-1*StdStat);
  pval = pvalNorm;

  arma::vec k1rootsLocal = {3, -3};
  if (StdStat > m_SPA_Cutoff){

    arma::uvec posG1 = arma::find(GVec != 0);
    double VarW1 = sum(VarWVec(posG1));
    double VarW0 = VarW - VarW1;
    double Ratio0 = VarW0 / VarW;

    SaddleResult resSPA = MAIN_SPA(statVal, adjGVec, k1rootsLocal, VarS, VarW, Ratio0, posG1);
    pval = resSPA.pval;
  }

  Beta = statVal / VarS;
  seBeta = std::abs(Beta) / StdStat;
  zScore = statVal / sqrt(VarS);
}


void POLMMClass::getRegionPVec(arma::vec GVec,
                               double& Stat,
                               double& Beta,
                               double& seBeta,
                               double& pval0,
                               double& pval1,
                               arma::vec& P1Vec,
                               arma::vec& P2Vec) {
  arma::vec test11 = getTime();

  arma::vec test21 = getTime();
  arma::vec adjGVec = getadjGFast(GVec);

  arma::vec test22 = getTime();
  arma::uvec testposG1 = arma::find(GVec != 0);

  arma::vec test23 = getTime();
  Stat = getStatFast(adjGVec);

  arma::vec ZPZ_adjGVec = get_ZPZ_adjGVec(adjGVec);

  arma::vec test24 = getTime();
  double VarS = as_scalar(adjGVec.t() * ZPZ_adjGVec);
  double StdStat = std::abs(Stat) / sqrt(VarS);
  double pvalNorm = 2 * arma::normcdf(-1*StdStat);
  double pvalLocal = pvalNorm;

  arma::vec k1rootsLocal = {3, -3};
  if (StdStat > m_SPA_Cutoff){

    arma::vec VarWVec = getVarWVec(adjGVec);
    double VarW = sum(VarWVec);


    arma::uvec posG1 = arma::find(GVec != 0);
    double VarW1 = sum(VarWVec(posG1));
    double VarW0 = VarW - VarW1;
    double Ratio0 = VarW0 / VarW;

    SaddleResult resSPA = MAIN_SPA(Stat, adjGVec, k1rootsLocal, VarS, VarW, Ratio0, posG1);
    pvalLocal = resSPA.pval;
  }

  pval0 = pvalNorm;
  pval1 = pvalLocal;
  Beta = Stat / VarS;
  seBeta = Beta / StdStat;

  P1Vec = adjGVec;
  P2Vec = ZPZ_adjGVec;


  arma::vec test12 = getTime();

  m_diffTimePOLMM1 += (test12 - test11);


}

arma::vec POLMMClass::getadjGFast(arma::vec GVec) {

  arma::vec test111 = getTime();

  arma::vec XR_Psi_RG(m_p, arma::fill::zeros);

  for (int i = 0; i < m_n; i++){
    if (GVec(i) != 0){
      XR_Psi_RG += m_XR_Psi_R.col(i) * GVec(i);
    }
  }

  arma::vec test112 = getTime();
  arma::vec adjGVec = GVec - m_XXR_Psi_RX * XR_Psi_RG;

  arma::vec test113 = getTime();

  m_diffTimePOLMM5 += (test112 - test111);
  m_diffTimePOLMM6 += (test113 - test112);

  return adjGVec;
}

double POLMMClass::getStatFast(arma::vec adjGVec) {
  double Stat = 0;
  for (int i = 0; i < m_n; i++){
    if (adjGVec(i) != 0){
      Stat += adjGVec(i) * m_RymuVec(i);
    }
  }
  return Stat;
}

arma::vec POLMMClass::get_ZPZ_adjGVec(arma::vec adjGVec) {
  arma::vec test111 = getTime();

  arma::vec adjGVecLong = Vec2LongVec(adjGVec, m_n, m_J);

  arma::vec iSigmaGVec(m_n * (m_J-1), arma::fill::zeros);

  getPCGofSigmaAndVector(adjGVecLong, iSigmaGVec);

  arma::vec test112 = getTime();

  arma::vec PZ_adjGVec = iSigmaGVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigmaGVec);
  arma::vec ZPZ_adjGVec = LongVec2Vec(PZ_adjGVec, m_n, m_J);

  arma::vec test113 = getTime();

  m_diffTimePOLMM7 += (test112 - test111);
  m_diffTimePOLMM8 += (test113 - test112);

  return ZPZ_adjGVec;
}


void POLMMClass::getPCGofSigmaAndCovaMat(arma::mat xMat,
                                         arma::mat& iSigma_xMat) {
  int p1 = xMat.n_cols;
  for (int i = 0; i < p1; i++){
    arma::vec y1Vec = xMat.col(i);
    arma::vec iSigma_y1Vec = iSigma_xMat.col(i);
    getPCGofSigmaAndVector(y1Vec, iSigma_y1Vec);
    iSigma_xMat.col(i) = iSigma_y1Vec;
  }
}


void POLMMClass::getPCGofSigmaAndVector(arma::vec y1Vec,
                                        arma::vec& xVec,
                                        std::string excludechr) {
  arma::mat xMat = convert2(xVec, m_n, m_J);
  arma::mat y1Mat = convert2(y1Vec, m_n, m_J);

  unsigned int iter = 0;
  arma::mat r2Mat = y1Mat - getSigmaxMat(xMat, excludechr);
  double meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));

  if (meanL2 <= m_tolPCG){

  }else{
    iter++;
    arma::cube InvBlockDiagSigma = getInvBlockDiagSigma();
    arma::mat z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);

    arma::mat z1Mat, r1Mat;
    double beta1 = 0;
    arma::mat pMat = z2Mat;
    arma::mat ApMat = getSigmaxMat(pMat, excludechr);
    double alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
    xMat = xMat + alpha * pMat;
    r1Mat = r2Mat;
    z1Mat = z2Mat;
    r2Mat = r1Mat - alpha * ApMat;

    meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));

    while (meanL2 > m_tolPCG && iter < m_maxiterPCG){
      iter++;


      z2Mat = solverBlockDiagSigma(InvBlockDiagSigma, r2Mat);

      beta1 = getInnerProd(z2Mat, r2Mat) / getInnerProd(z1Mat, r1Mat);
      pMat = z2Mat + beta1 * pMat;
      ApMat = getSigmaxMat(pMat, excludechr);
      alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);

      xMat = xMat + alpha * pMat;
      r1Mat = r2Mat;
      z1Mat = z2Mat;
      r2Mat = r1Mat - alpha * ApMat;
      meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));
    }
  }

  xVec = convert1(xMat, m_n, m_J);
  if (iter >= m_maxiterPCG){
    std::cout << "    PCG did not converge (increase maxiter)" << std::endl;
  }
  if (m_showInfo)
    std::cout << "    PCG iterations: " << iter << std::endl;

}


void POLMMClass::getPCGofSigmaAndVector(arma::vec y1Vec,
                                        arma::vec& xVec) {
  arma::vec test11 = getTime();
  arma::mat xMat = Vec2Mat(xVec, m_n, m_J);
  arma::mat y1Mat = Vec2Mat(y1Vec, m_n, m_J);
  arma::vec test12 = getTime();
  m_diffTimePOLMM4 += (test12 - test11);

  unsigned int iter = 0;
  arma::mat r2Mat = y1Mat - getSigmaxMat(xMat);
  double meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));

  if (meanL2 <= m_tolPCG){

  }else{

    iter++;
    arma::mat z2Mat = solverBlockDiagSigma(r2Mat);

    arma::mat z1Mat, r1Mat;
    double beta1 = 0;
    arma::mat pMat = z2Mat;
    arma::mat ApMat = getSigmaxMat(pMat);
    double alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
    xMat = xMat + alpha * pMat;
    r1Mat = r2Mat;
    z1Mat = z2Mat;
    r2Mat = r1Mat - alpha * ApMat;

    meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));

    while (meanL2 > m_tolPCG && iter < m_maxiterPCG){

      iter++;

      z2Mat = solverBlockDiagSigma(r2Mat);

      beta1 = getInnerProd(z2Mat, r2Mat) / getInnerProd(z1Mat, r1Mat);
      pMat = z2Mat + beta1 * pMat;
      ApMat = getSigmaxMat(pMat);
      alpha = getInnerProd(z2Mat, r2Mat) / getInnerProd(pMat, ApMat);
      xMat = xMat + alpha * pMat;
      r1Mat = r2Mat;
      z1Mat = z2Mat;
      r2Mat = r1Mat - alpha * ApMat;
      meanL2 = sqrt(getInnerProd(r2Mat, r2Mat)) / sqrt(m_n * (m_J-1));


    }
  }

  xVec = Mat2Vec(xMat, m_n, m_J);
  if (iter >= m_maxiterPCG){
    std::cout << "    PCG did not converge (increase maxiter)" << std::endl;
  }
  if (m_printPCGInfo)
    std::cout << "    PCG iterations: " << iter << std::endl;

}


arma::mat POLMMClass::getSigmaxMat(arma::mat& xMat) {
  arma::mat iR_xMat = m_iRMat % xMat;

  arma::vec test11 = getTime();

  arma::mat iPsi_iR_xMat = getiPsixMat(iR_xMat);

  arma::vec test12 = getTime();
  m_diffTimePOLMM2 += (test12 - test11);

  arma::mat yMat = m_iRMat % iPsi_iR_xMat;

  if (m_tau == 0){}
  else{
    arma::vec tZ_xMat = arma::sum(xMat, 1);
    arma::vec V_tZ_xMat = m_SparseGRM * tZ_xMat;
    yMat.each_col() += m_tau * V_tZ_xMat;
  }
  return yMat;
}


arma::mat POLMMClass::getiPsixMat(arma::mat xMat) {
  arma::mat iPsi_xMat(m_n, m_J-1);
  for (int i = 0; i < m_n; i ++){
    double sumx = arma::sum(xMat.row(i));
    double sumx_divided_by_mu = sumx / m_muMat(i, m_J-1);
    for (int j = 0; j < m_J-1; j++){

      iPsi_xMat(i,j) = sumx_divided_by_mu + xMat(i,j) / m_muMat(i,j);
    }
  }
  return iPsi_xMat;
}


arma::mat POLMMClass::getPsixMat(arma::mat xMat) {
  arma::mat Psi_xMat(m_n, m_J-1);

  for (int i = 0; i < m_n; i++){
    arma::rowvec muVec(m_J-1);
    for (int j = 0; j < m_J-1; j++){
      Psi_xMat(i,j) = m_muMat(i,j) * xMat(i,j);
      muVec(j) = m_muMat(i,j);
    }
    double sum_mu_x = sum(Psi_xMat.row(i));
    Psi_xMat.row(i) -= muVec * sum_mu_x;
  }
  return Psi_xMat;
}


arma::cube POLMMClass::getInvBlockDiagSigma() {
  arma::vec DiagGRM;
  if (m_flagSparseGRM){

    DiagGRM = m_tau * m_SparseGRM.diag();

  }else{
    throw std::runtime_error("DenseGRM not available in mtMarker build");
  }

  double temp;

  arma::cube InvBlockDiagSigma(m_J-1, m_J-1, m_n, arma::fill::zeros);
  for (int i = 0; i < m_n; i++){
    for (int j2 = 0; j2 < m_J-1; j2++){
      for (int j1 = 0; j1 < m_J-1; j1++){
        temp = m_iRMat(i,j2) * (1 / m_muMat(i, m_J-1)) * m_iRMat(i,j1) + DiagGRM(i);
        if (j2 == j1){
          temp += m_iRMat(i,j2) * (1 / m_muMat(i,j2)) * m_iRMat(i,j1);
        }
        InvBlockDiagSigma(j2, j1, i) = temp;
      }
    }
    InvBlockDiagSigma.slice(i) = inv(InvBlockDiagSigma.slice(i));
  }
  return InvBlockDiagSigma;
}

arma::mat POLMMClass::solverBlockDiagSigma(arma::mat& xMat) {
  arma::vec test11 = getTime();
  arma::mat outMat(m_n, m_J-1);
  for (int i = 0; i < m_n; i++){
    outMat.row(i) = xMat.row(i) * m_InvBlockDiagSigma.slice(i);
  }
  arma::vec test12 = getTime();
  m_diffTimePOLMM3 += (test12 - test11);
  return outMat;
}

SaddleResult POLMMClass::MAIN_SPA(double Stat,
                                arma::vec adjGVec,
                                arma::vec K1roots,
                                double VarP,
                                double VarW,
                                double Ratio0,
                                arma::uvec posG1) {
  SaddleResult resSPA = fastSaddle_Prob(Stat, VarP, VarW, Ratio0, K1roots,
                                      adjGVec.elem(posG1), m_muMat.rows(posG1), m_iRMat.rows(posG1));
  return resSPA;
}

void POLMMClass::setRPsiR() {

  arma::vec RPsiRVec(m_n, arma::fill::zeros);
  arma::mat muRMat = m_muMat.cols(0, m_J-2) / m_iRMat;
  for (int i = 0; i < m_n; i++){
    for (int j1 = 0; j1 < m_J-1; j1++){

      RPsiRVec(i) += muRMat(i,j1) / m_iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for (int j2 = j1+1; j2 < m_J-1; j2++){

        RPsiRVec(i) -= 2* muRMat(i,j1) * muRMat(i,j2);
      }
    }
  }
  m_RPsiR = RPsiRVec;
}

arma::vec POLMMClass::getVarWVec(arma::vec adjGVec) {
  arma::vec VarWVec = m_RPsiR % pow(adjGVec, 2);
  return VarWVec;
}

double K0(double x,
          arma::mat muMat,
          arma::mat cMat,
          double m1) {
  arma::mat temp1Mat = - muMat + muMat % exp(cMat * x);
  arma::vec temp1Vec = log(1 + arma::sum(temp1Mat, 1));
  double y = sum(temp1Vec) - m1 * x;


  return y;
}

arma::vec K12(double x,
              arma::mat muMat,
              arma::mat cMat,
              double m1) {
  arma::mat temp0Mat = muMat % exp(cMat * x);
  arma::mat temp1Mat = - muMat + temp0Mat;
  arma::mat temp2Mat = temp0Mat % cMat;
  arma::mat temp3Mat = temp2Mat % cMat;

  arma::vec temp1Vec = 1 + arma::sum(temp1Mat, 1);
  arma::vec temp2Vec = arma::sum(temp2Mat, 1);
  arma::vec temp3Vec = arma::sum(temp3Mat, 1);

  arma::vec yVec(2);
  yVec(0) = sum(temp2Vec / temp1Vec) - m1;

  yVec(1) = sum((temp3Vec % temp1Vec - pow(temp2Vec, 2)) / pow(temp1Vec, 2));

  return yVec;
}

RootResult fastgetroot_K1(double Stat,
                          double initX,
                          double Ratio0,
                          arma::mat muMat,
                          arma::mat cMat,
                          double m1) {
  double x = initX;
  double K1 = 0;
  double K2 = 0;
  double diffX = arma::datum::inf;
  bool converge = true;
  double tol = 0.0001;
  int maxiter = 100;
  int iter = 0;

  for (iter = 0; iter < maxiter; iter ++){
    double oldX = x;
    double oldDiffX = diffX;
    double oldK1 = K1;

    arma::vec K12Vec = K12(x, muMat, cMat, m1);

    K1 = K12Vec(0) - Stat + Ratio0 * x;
    K2 = K12Vec(1) + Ratio0;

    diffX = -1 * K1 / K2;


    if (!std::isfinite(K1)){


      x = arma::sign(Stat) * arma::datum::inf;
      K2 = 0;
      break;
    }

    if (arma::sign(K1) != arma::sign(oldK1)){
      while (std::abs(diffX) > std::abs(oldDiffX) - tol){
        diffX = diffX / 2;
      }
    }
    if (std::abs(diffX) < tol) break;

    x = oldX + diffX;

  }

  if (iter == maxiter - 1)
    converge = false;

  return {x, iter, converge, K2};
}

double fastGet_Saddle_Prob(double Stat,
                           double zeta,
                           double K2,
                           double Ratio0,
                           arma::mat muMat,
                           arma::mat cMat,
                           double m1,
                           bool lowerTail) {
  double k1 = K0(zeta, muMat, cMat, m1) + 0.5 * pow(zeta, 2) * Ratio0;
  double k2 = K2;
  double pval = 0;
  if (std::isfinite(k1) && std::isfinite(k2)) {
    double w = arma::sign(zeta) * sqrt(2 * (zeta * Stat - k1));
    double v = zeta * sqrt(K2);

    double Z = w + 1/w * log(v/w);
    pval = arma::normcdf(arma::sign(lowerTail-0.5) * Z);
  }

  return pval;
}


SaddleResult fastSaddle_Prob(double Stat,
                           double VarP,
                           double VarW,
                           double Ratio0,
                           arma::vec K1roots,
                           arma::vec adjGVec1,
                           arma::mat muMat1,
                           arma::mat iRMat1) {
  int J = muMat1.n_cols;
  int N1 = muMat1.n_rows;

  muMat1 = muMat1.cols(0, J-2);
  double adjStat = Stat / sqrt(VarP);

  double sqrtVarW = sqrt(VarW);

  arma::mat cMat(N1, J-1);
  for (int i = 0; i < N1; i ++){
    for (int j = 0; j < J-1; j ++){
      cMat(i,j) = adjGVec1(i) / iRMat1(i,j) / sqrtVarW;
    }
  }

  double m1 = arma::accu(muMat1 % cMat);

  RootResult outUni1 = fastgetroot_K1(std::abs(adjStat), std::min(K1roots(0), 5.0),
                                      Ratio0, muMat1, cMat, m1);
  RootResult outUni2 = fastgetroot_K1(-1 * std::abs(adjStat), std::max(K1roots(1), -5.0),
                                      Ratio0, muMat1, cMat, m1);

  bool converge = false;
  double pval = 0;
  arma::vec k1rootsOut;

  if (outUni1.converge && outUni2.converge){

    double p1 = fastGet_Saddle_Prob(std::abs(adjStat), outUni1.root,
                                    outUni1.K2, Ratio0, muMat1, cMat, m1, false);

    double p2 = fastGet_Saddle_Prob(-1 * std::abs(adjStat), outUni2.root,
                                    outUni2.K2, Ratio0, muMat1, cMat, m1, true);

    pval = p1 + p2;

    converge = true;
    k1rootsOut = {outUni1.root, outUni2.root};
  }else{
    std::cout << "    SPA did not converge, using normal approximation" << std::endl;
    pval = 2 * arma::normcdf(-1 * std::abs(adjStat));
    k1rootsOut = K1roots;
  }

  if ((!std::isfinite(pval)) || (pval == 0)){
    std::cout << "    SPA gave invalid p-value, using normal approximation" << std::endl;
    pval = 2 * arma::normcdf(-1 * std::abs(adjStat));
    k1rootsOut = K1roots;
  }

  return {pval, converge, k1rootsOut};
}

void POLMMClass::setSeqMat(int NonZero_cutoff) {
  int n = NonZero_cutoff;
  std::cout << "    Setting up Efficient Resampling (ER) ..." << std::endl;
  uint32_t nER = pow(m_J, n);
  arma::Col<uint32_t> y = arma::linspace<arma::Col<uint32_t>>(0, nER-1, nER);
  m_SeqMat.resize(n, nER);
  uint32_t powJ = nER / m_J;
  arma::Col<uint8_t> SeqVec(nER);
  for (int i = 0; i < n; i++){
    int pos_row = n - 1 - i;
    for (uint32_t j = 0; j < nER; j++){
      SeqVec(j) = y(j) / powJ;
      y(j) = y(j) - SeqVec(j) * powJ;
    }
    powJ = powJ / m_J;
    m_SeqMat.row(pos_row) = SeqVec.t();
  }
}


double POLMMClass::MAIN_ER(arma::vec GVec,
                           arma::uvec posG1) {
  int N1 = posG1.size();
  uint32_t nER = pow(m_J, N1);
  arma::Mat<uint8_t> SeqMat = m_SeqMat.submat(0, 0, N1-1, nER-1);
  double pvalER = getPvalER(m_yVec.elem(posG1), GVec.elem(posG1), m_muMat.rows(posG1), m_iRMat.rows(posG1), SeqMat);

  return pvalER;
}


double getPvalER(arma::uvec yVec,
                 arma::vec GVec,
                 arma::mat muMat,
                 arma::mat iRMat,
                 arma::Mat<uint8_t> SeqMat) {
  uint32_t nER = SeqMat.n_cols;
  arma::vec StatVec = getStatVec(SeqMat, GVec, muMat, iRMat);

  arma::Col<uint8_t> yVecU8 = arma::conv_to<arma::Col<uint8_t>>::from(yVec);
  double StatObs = arma::as_scalar(getStatVec(yVecU8, GVec, muMat, iRMat));

  double eps = 1e-10;

  double pvalER_pos = 0;
  double pvalER_neg = 0;
  double absStatObs = std::abs(StatObs);
  for (uint32_t i = 0; i < nER; i++){


    double StatTmp = StatVec(i);

    if (StatTmp > absStatObs + eps){
      pvalER_pos += getProb(SeqMat.col(i), muMat);
    }else if (StatTmp > absStatObs - eps){
      pvalER_pos += 0.5 * getProb(SeqMat.col(i), muMat);
    }

    if (StatTmp < -1 * absStatObs - eps){
      pvalER_neg += getProb(SeqMat.col(i), muMat);
    }else if (StatTmp < -1 * absStatObs + eps){
      pvalER_neg += 0.5 * getProb(SeqMat.col(i), muMat);
    }
  }


  double pvalER = pvalER_pos + pvalER_neg;
  return pvalER;
}

arma::vec getStatVec(arma::Mat<uint8_t> SeqMat,
                     arma::vec GVec,
                     arma::mat muMat,
                     arma::mat iRMat) {
  int n = muMat.n_rows;
  int J = muMat.n_cols;
  int nER = SeqMat.n_cols;

  arma::vec StatVec(nER);

  arma::mat A(n, J-1);
  for (int i = 0; i < J-1; i++){
    A.col(i) = GVec / iRMat.col(i);
  }

  double a1 = arma::accu(A % muMat.cols(0, J-2));

  for (int i = 0; i < nER; i++){
    double a2 = 0;
    for (int j = 0; j < n; j++){
      int idxL = SeqMat(j, i);
      if (idxL != J-1){
        a2 += A(j,idxL);
      }
    }

    StatVec(i) = a2 - a1;
  }

  return StatVec;
}

double getProbOne(arma::Col<uint8_t> SeqVec,
                  arma::mat muMat) {
  int n = muMat.n_rows;
  double tempProb = 1;
  for (int j = 0; j < n; j++){
    tempProb *= muMat(j, SeqVec(j));
  }
  return tempProb;
}


double getProb(arma::Mat<uint8_t> SeqMat,
               arma::mat muMat) {
  int nER = SeqMat.n_cols;

  double prob = 0;

  for (int i = 0; i < nER; i++){
    arma::Col<uint8_t> SeqVec = SeqMat.col(i);
    double tempProb = getProbOne(SeqVec, muMat);
    prob += tempProb;
  }

  return prob;
}


void POLMMClass::updateMats() {

  m_eta = m_Cova * m_beta + m_bVec;


  double tmpExp, tmpnu0, tmpnu1;
  for (int i = 0; i < m_n; i ++){
    tmpnu0 = 0;
    for (int j = 0; j < m_J-1; j ++){
      tmpExp = exp(m_eps(j) - m_eta(i));
      tmpnu1 = tmpExp / (1 + tmpExp);
      m_muMat(i,j) = tmpnu1 - tmpnu0;
      m_WMat(i,j) = tmpnu1 * (1 - tmpnu1);
      m_mMat(i,j) = tmpnu1 + tmpnu0 - 1;
      m_nuMat(i,j) = tmpnu1;
      tmpnu0 = tmpnu1;
    }
    int j = m_J-1;
    tmpnu1 = 1;
    m_muMat(i,j) = tmpnu1 - tmpnu0;
    m_WMat(i,j) = tmpnu1 * (1 - tmpnu1);
    m_mMat(i,j) = tmpnu1 + tmpnu0 - 1;
    m_nuMat(i,j) = tmpnu1;
  }


  for (int i = 0; i < m_n; i ++){
    for (int j = 0; j < m_J-1; j ++){
      m_iRMat(i,j) = 1 / (m_mMat(i,j) - m_mMat(i, m_J-1));
    }
  }


  arma::mat xMat = m_yMat.cols(0, m_J-2) - m_muMat.cols(0, m_J-2);
  arma::mat iPsi_xMat = getiPsixMat(xMat);
  for (int i = 0; i < m_n; i++){
    for (int j = 0; j < m_J-1; j++)
      m_YMat(i,j) = m_eta(i) + (m_iRMat(i,j) * iPsi_xMat(i,j));
  }
}

void POLMMClass::updateEpsOneStep() {

  arma::vec d1eps(m_J-2, arma::fill::zeros);
  arma::mat d2eps(m_J-2, m_J-2, arma::fill::zeros);
  double temp1, temp2, temp3;

  for (int k = 1; k < m_J-1; k++){
    for (int i = 0; i < m_n; i++){
      temp1 = m_yMat(i, k) / m_muMat(i, k) - m_yMat(i, k+1) / m_muMat(i, k+1);
      temp2 = - m_yMat(i, k) / m_muMat(i, k) / m_muMat(i, k) - m_yMat(i, k+1) / m_muMat(i, k+1) / m_muMat(i, k+1);

      d1eps(k - 1) += m_WMat(i, k) * temp1;
      d2eps(k - 1, k - 1) += m_WMat(i, k) * (1 - 2 * m_nuMat(i, k)) * temp1 + m_WMat(i, k) * m_WMat(i, k) * temp2;
      if (k < m_J-2){
        temp3 = m_WMat(i,k) * m_WMat(i, k+1) * m_yMat(i, k+1) / m_muMat(i, k+1) / m_muMat(i, k+1);
        d2eps(k-1, k) += temp3;
        d2eps(k, k-1) += temp3;
      }
    }
  }

  arma::vec deps = -1 * inv(d2eps) * d1eps;

  for (int k = 1; k < m_J-1; k ++){
    m_eps(k) += deps(k-1);
  }
}

void POLMMClass::updateEps() {
  for (unsigned int iter = 0; iter < m_maxiterEps; iter ++){

    arma::vec eps0 = m_eps;

    updateEpsOneStep();
    updateMats();

    double diffeps = arma::max(arma::abs(m_eps - eps0)/(arma::abs(m_eps) + arma::abs(eps0) + m_tolEps));

    if (diffeps < m_tolEps){

      std::cout << "    Eps converged at iteration " << iter << std::endl;
      break;
    }
  }
}

void POLMMClass::updatePara(std::string excludechr) {
  getPCGofSigmaAndCovaMat(m_CovaMat, m_iSigma_CovaMat, excludechr);
  arma::vec YVec = convert1(m_YMat, m_n, m_J);
  getPCGofSigmaAndVector(YVec, m_iSigma_YVec, excludechr);


  arma::mat XSigmaX = inv(m_CovaMat.t() * m_iSigma_CovaMat);
  arma::vec Cova_iSigma_YVec = m_CovaMat.t() * m_iSigma_YVec;
  m_beta = XSigmaX * Cova_iSigma_YVec;
  m_iSigmaX_XSigmaX = m_iSigma_CovaMat * XSigmaX;


  arma::vec Z_iSigma_YVec = ZMat(m_iSigma_YVec);
  arma::vec Z_iSigma_Xbeta = ZMat(m_iSigma_CovaMat * m_beta);
  arma::vec tempVec = Z_iSigma_YVec - Z_iSigma_Xbeta;
  m_bVec = m_tau * getKinbVecPOLMM(tempVec, excludechr);
}


arma::mat POLMMClass::solverBlockDiagSigma(arma::cube& InvBlockDiagSigma,
                                           arma::mat& xMat) {
  arma::mat outMat(m_n, m_J-1);
  for (int i = 0; i < m_n; i++){
    outMat.row(i) = xMat.row(i) * InvBlockDiagSigma.slice(i);
  }
  return outMat;
}


arma::mat POLMMClass::getSigmaxMat(arma::mat xMat,
                                   std::string excludechr) {
  arma::mat iR_xMat = m_iRMat % xMat;
  arma::mat iPsi_iR_xMat = getiPsixMat(iR_xMat);
  arma::mat yMat = m_iRMat % iPsi_iR_xMat;
  if (m_tau == 0){}
  else{
    arma::vec tZ_xMat = getRowSums(xMat);
    arma::vec V_tZ_xMat = getKinbVecPOLMM(tZ_xMat, excludechr);
    yMat.each_col() += m_tau * V_tZ_xMat;
  }
  return yMat;
}

arma::vec POLMMClass::getKinbVecPOLMM(arma::vec bVec,
                                      std::string excludeChr) {
  arma::vec KinbVec = m_SparseGRM * bVec;
  return KinbVec;
}


void POLMMClass::getPCGofSigmaAndCovaMat(arma::mat xMat,
                                         arma::mat& iSigma_xMat,
                                         std::string excludechr) {
  int p1 = xMat.n_cols;
  for (int i = 0; i < p1; i++){

    arma::vec y1Vec = xMat.col(i);
    arma::vec iSigma_y1Vec = iSigma_xMat.col(i);
    getPCGofSigmaAndVector(y1Vec, iSigma_y1Vec, excludechr);

    iSigma_xMat.col(i) = iSigma_y1Vec;
  }
}

double POLMMClass::getVarP(arma::vec adjGVec,
                           std::string excludechr) {
  arma::vec adjGVecLong = tZMat(adjGVec);
  arma::vec iSigmaGVec(m_n * (m_J-1), arma::fill::zeros);
  getPCGofSigmaAndVector(adjGVecLong, iSigmaGVec, excludechr);
  double VarP = as_scalar(adjGVecLong.t() * (iSigmaGVec - m_iSigmaX_XSigmaX * (m_CovaMat.t() * iSigmaGVec)));
  return VarP;
}

arma::vec convert1(arma::mat xMat,
                   int n, int J) {
  arma::vec xVec(n*(J-1));
  int index = 0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < J-1; j++){
      xVec(index) = xMat(i,j);
      index++;
    }
  }
  return xVec;
}

arma::mat convert2(arma::vec xVec,
                   int n, int J) {
  arma::mat xMat(n,(J-1));
  int index = 0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < J-1; j++){
      xMat(i,j) = xVec(index);
      index++;
    }
  }
  return xMat;
}

double getVarWFast(arma::vec adjGVec,
                   arma::vec RPsiRVec) {
  int n = adjGVec.size();
  double VarW = 0;
  for (int i = 0; i < n; i++){
    VarW += RPsiRVec(i) * adjGVec(i) * adjGVec(i);
  }
  return VarW;
}


ObjP getobjP(arma::mat Cova,
             arma::mat yMat,
             arma::mat muMat,
             arma::mat iRMat) {
  int n = muMat.n_rows;
  int J = muMat.n_cols;
  int p = Cova.n_cols;


  arma::mat XR_Psi_R(p, n*(J-1));
  arma::mat CovaMat = getCovaMat(Cova, J);
  for (int k = 0; k < p; k++){
    arma::mat xMat = convert2(CovaMat.col(k), n, J);
    arma::vec temp = convert1(getPsixMat(xMat / iRMat, muMat) / iRMat, n, J);
    XR_Psi_R.row(k) = temp.t();
  }
  arma::mat XXR_Psi_RX_new = Cova * inv(XR_Psi_R * CovaMat);


  arma::mat XR_Psi_R_new = sumCols(XR_Psi_R, J);
  arma::mat ymuMat = yMat - muMat;
  arma::mat RymuMat = ymuMat.cols(0, J-2) / iRMat;
  arma::vec RymuVec = arma::vectorise(sumCols(RymuMat, J));
  arma::vec RPsiR = getRPsiR(muMat, iRMat, n, J, p);

  ObjP objP;
  objP.n = n;
  objP.J = J;
  objP.p = p;
  objP.XXR_Psi_RX_new = XXR_Psi_RX_new;
  objP.XR_Psi_R_new = XR_Psi_R_new;
  objP.RymuVec = RymuVec;
  objP.RPsiR = RPsiR;
  objP.muMat = muMat;
  objP.iRMat = iRMat;
  return objP;
}

arma::vec getadjGFast(arma::vec GVec,
                      arma::mat XXR_Psi_RX_new,
                      arma::mat XR_Psi_R_new,
                      int n, int p) {

  arma::vec XR_Psi_RG1(p, arma::fill::zeros);
  for (int i = 0; i < n; i++){

    if (GVec(i) != 0){

      XR_Psi_RG1 += XR_Psi_R_new.col(i) * GVec(i);
    }

  }

  arma::vec adjGVec = GVec - XXR_Psi_RX_new * XR_Psi_RG1;
  return adjGVec;
}

double getStatFast(arma::vec GVec,
                   arma::vec RymuVec) {
  int n = GVec.size();
  double Stat = 0;
  for (int i = 0; i < n; i++){
    if (GVec(i) != 0){
      Stat += GVec(i) * RymuVec(i);
    }
  }
  return Stat;
}

AdjGResult outputadjGFast(arma::vec GVec,
                          const ObjP& objP) {
  arma::vec adjGVec = getadjGFast(GVec, objP.XXR_Psi_RX_new, objP.XR_Psi_R_new, objP.n, objP.p);
  double Stat = getStatFast(adjGVec, objP.RymuVec);
  double VarW = getVarWFast(adjGVec, objP.RPsiR);
  return {adjGVec, Stat, VarW};
}


arma::vec getRowSums(arma::mat xMat) {
  int n1 = xMat.n_rows;
  int n2 = xMat.n_cols;
  arma::vec y1Vec(n1, arma::fill::zeros);
  for (int i = 0; i < n1; i++){
    for (int j = 0; j < n2; j++){
      y1Vec(i) += xMat(i,j);
    }
  }
  return y1Vec;
}


arma::vec getRPsiR(arma::mat muMat,
                   arma::mat iRMat,
                   int n, int J, int p) {

  arma::vec RPsiRVec(n, arma::fill::zeros);
  arma::mat muRMat = muMat.cols(0, J-2) / iRMat;
  for (int i = 0; i < n; i++){
    for (int j1 = 0; j1 < J-1; j1++){

      RPsiRVec(i) += muRMat(i,j1) / iRMat(i,j1) - muRMat(i,j1) * muRMat(i,j1);
      for (int j2 = j1+1; j2 < J-1; j2++){

        RPsiRVec(i) -= 2* muRMat(i,j1) * muRMat(i,j2);
      }
    }
  }

  return RPsiRVec;
}

double calCV(arma::vec xVec){
  unsigned int n = xVec.size();
  double Mean = arma::mean(xVec);
  double Sd = arma::stddev(xVec);
  double CV = (Sd/Mean)/n;
  return CV;
}


arma::mat sumCols(arma::mat xMat,
                  int J) {
  int n = xMat.n_cols / (J-1);
  int p = xMat.n_rows;
  arma::mat outMat(p, n, arma::fill::zeros);
  int index = 0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < J-1; j++){
      outMat.col(i) += xMat.col(index);
      index++;
    }
  }
  return outMat;
}


arma::mat getPsixMat(arma::mat xMat,
                     arma::mat muMat) {
  int n = muMat.n_rows;
  int J = muMat.n_cols;

  arma::mat Psi_xMat(n, J-1);

  for (int i = 0; i < n; i++){
    arma::rowvec muVec(J-1);
    for (int j = 0; j < J-1; j++){
      Psi_xMat(i,j) = muMat(i,j) * xMat(i,j);
      muVec(j) = muMat(i,j);
    }
    double sum_mu_x = sum(Psi_xMat.row(i));
    Psi_xMat.row(i) -= muVec * sum_mu_x;
  }
  return Psi_xMat;
}


arma::mat getCovaMat(arma::mat Cova, unsigned int J) {
  unsigned int n = Cova.n_rows;
  unsigned int p = Cova.n_cols;

  arma::mat CovaMat(n * (J-1), p);
  int index = 0;
  for (unsigned int i = 0; i < n; i++){
    for (unsigned int j = 0; j < J-1; j++){
      CovaMat.row(index) = Cova.row(i);
      index++;
    }
  }
  return CovaMat;
}

arma::vec getTime(){
  arma::vec Time(2, arma::fill::zeros);
  struct timeval time;
  Time(0) = 0;
  if (!gettimeofday(&time,NULL))
    Time(0) = (double)time.tv_sec + (double)time.tv_usec * .000001;
  Time(1) = (double)clock() / CLOCKS_PER_SEC;
  return Time;
}

void printTime(arma::vec t1, arma::vec t2, std::string message){
  double wallTime = t2(0) - t1(0);
  double cpuTime = t2(1) - t1(1);
  if (wallTime < 60){
    Rprintf ("    It took %.2f seconds (%.2f CPU seconds) to %s.\n",
             wallTime, cpuTime, message.c_str());
  }else if (wallTime < 3600){
    Rprintf ("    It took %.2f minutes (%.2f CPU minutes) to %s.\n",
             wallTime/60, cpuTime/60, message.c_str());
  }else{
    Rprintf ("    It took %.2f hours (%.2f CPU hours) to %s.\n",
             wallTime/3600, cpuTime/3600, message.c_str());
  }
}

arma::vec nb (unsigned int n) {
  arma::vec v(n);
  for (unsigned int i = 0; i < n; ++i)
    v(i) = (arma::randu() < 0.5) ? 0.0 : 1.0;
  return v;
}

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat) {
  double innerProd = arma::accu(x1Mat % x2Mat);
  return innerProd;
}

arma::vec Vec2LongVec(arma::vec xVec, int n, int J) {
  arma::vec yVec(n * (J-1));
  int index = 0;
  for (int i = 0; i < n; i ++){
    for (int j = 0; j < J-1; j ++){
      yVec(index) = xVec(i);
      index++;
    }
  }
  return yVec;
}

arma::vec LongVec2Vec(arma::vec xVec, int n, int J) {
  arma::vec yVec(n, arma::fill::zeros);
  int index = 0;
  for (int i = 0; i < n; i ++){
    for (int j = 0; j < J-1; j ++){
      yVec(i) += xVec(index);
      index++;
    }
  }
  return yVec;
}

arma::mat Vec2Mat(arma::vec xVec, int n, int J) {
  arma::mat xMat(n, (J-1));
  int index = 0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < J-1; j++){
      xMat(i, j) = xVec(index);
      index++;
    }
  }
  return xMat;
}

arma::vec Mat2Vec(arma::mat xMat, int n, int J) {
  arma::vec xVec(n * (J-1));
  int index = 0;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < J-1; j++){
      xVec(index) = xMat(i,j);
      index++;
    }
  }
  return xVec;
}

}
