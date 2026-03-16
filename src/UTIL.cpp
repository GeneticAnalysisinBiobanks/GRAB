// UTIL.cpp -- Implementations of shared utility functions

#include <RcppArmadillo.h>
#include "UTIL.h"
#include <sys/time.h>
#include <stdexcept>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/beta.hpp>

arma::vec nb (unsigned int n) {
  arma::vec v(n);
  for (unsigned int i = 0; i < n; ++i)
    v(i) = (arma::randu() < 0.5) ? 0.0 : 1.0;
  return v;
}

double getWeights (std::string kernel,
                   double freq,
                   arma::vec wBeta) {
  if (wBeta.size() != 2)
    throw std::runtime_error("The size of argument wBeta should be 2.");

  double weights = 1.0;
  if (kernel == "linear.weighted") {
    boost::math::beta_distribution<double> bdist(wBeta(0), wBeta(1));
    weights = boost::math::pdf(bdist, freq);
  }
  return weights;
}

void imputeGeno(arma::vec& GVec,
                double altFreq,
                std::vector<uint32_t> indexForMissing,
                std::string imputeMethod) {
  int nMissing = indexForMissing.size();

  double imputeG = 0;

  if (imputeMethod == "mean")
    imputeG = 2 * altFreq;

  if (imputeMethod == "none")
    imputeG = arma::datum::nan;

  if (imputeMethod == "bestguess")
    imputeG = std::round(2 * altFreq);

  for (int i = 0; i < nMissing; i++){
    uint32_t index = indexForMissing.at(i);
    GVec.at(index) = imputeG;
  }
}


bool imputeGenoAndFlip(arma::vec& GVec,
                       double altFreq,
                       std::vector<uint32_t> indexForMissing,
                       std::string impute_method) {
  int nMissing = indexForMissing.size();

  double imputeG = 0;
  if (impute_method == "mean"){
    imputeG = 2 * altFreq;
  }

  if (impute_method == "minor"){
    if (altFreq > 0.5){
      imputeG = 2;
    }else{
      imputeG = 0;
    }
  }

  for (int i = 0; i < nMissing; i++){
    uint32_t index = indexForMissing.at(i);
    GVec.at(index) = imputeG;
  }

  bool flip = false;
  if (altFreq > 0.5){
    GVec = 2 - GVec;
    flip = true;
  }

  return flip;
}


bool imputeGenoAndFlip(arma::vec& GVec,
                       double altFreq,
                       std::vector<uint32_t> indexForMissing,
                       double missingRate,
                       std::string impute_method,
                       const std::string& method) {
  int nMissing = indexForMissing.size();

  double imputeG = 0;
  if (impute_method == "mean"){
    imputeG = 2 * altFreq;
  }

  if (impute_method == "minor"){
    if (altFreq > 0.5){
      imputeG = 2;
    }else{
      imputeG = 0;
    }
  }

  for (int i = 0; i < nMissing; i++){
    uint32_t index = indexForMissing.at(i);
    GVec.at(index) = imputeG;
  }

  bool flip = false;

  if (method != "WtCoxG") {
    if (altFreq > 0.5){
      GVec = 2 - GVec;
      flip = true;
    }
  }
  return flip;
}

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat) {
  double innerProd = arma::accu(x1Mat % x2Mat);
  return innerProd;
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

void printTimeDiff(arma::vec timeDiff,
                   std::string message) {
  double wallTime = timeDiff(0);
  double cpuTime = timeDiff(1);
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

double getinvStd(double freq) {
  double Std = sqrt(2 * freq * (1-freq));
  if (Std == 0)
    return 0;
  else
    return 1/Std;
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

void gethwepval(arma::vec GVec,
                double& hwepval,
                double hwepvalCutoff) {
  int n = GVec.size();
  double n0 = 0;
  double n1 = 0;
  double n2 = 0;
  for (int i = 0; i < n; i++) {
    double g = GVec.at(i);
    if (g < hwepvalCutoff){
      n0 += 1;
      continue;
    }

    if (g > 2 - hwepvalCutoff) {
      n2 += 1;
      continue;
    }


    if ((g > 1 - hwepvalCutoff) & (g < 1 + hwepvalCutoff)) {
      n1 += 1;
      continue;
    }
  }


  double nsum = n0 + n1 + n2;

  if (nsum == 0){
    hwepval = 0;
  }else{
    double Gfreq = (0.5 * n1 + n2) / nsum;

    arma::vec exp_GVec = {(1-Gfreq)*(1-Gfreq) * nsum, 2*Gfreq*(1-Gfreq) * nsum, Gfreq*Gfreq * nsum};

    double hwechisq = arma::accu(arma::square(arma::vec{n0, n1, n2} - exp_GVec) / exp_GVec);

    boost::math::chi_squared_distribution<double> chisq_dist(1.0);
    hwepval = boost::math::cdf(boost::math::complement(chisq_dist, hwechisq));
  }
}
