#ifndef UTIL_H
#define UTIL_H

// UTIL.h -- Shared utility functions (imputation, genotype helpers, timing)

#include <RcppArmadillo.h>
#include <sys/time.h>

double getWeights(std::string kernel,
                  double freq,
                  arma::vec wBeta);

void imputeGeno(arma::vec& GVec,
                double freq,
                std::vector<uint32_t> posMissingGeno);

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat);


arma::vec Vec2LongVec(arma::vec xVec, int n, int J);


arma::vec LongVec2Vec(arma::vec xVec, int n, int J);


arma::mat Vec2Mat(arma::vec xVec, int n, int J);


arma::vec Mat2Vec(arma::mat xMat, int n, int J);

arma::mat sumCols(arma::mat xMat, int J);

arma::vec getRPsiR(arma::mat muMat, arma::mat iRMat, int n, int J, int p);

bool imputeGenoAndFlip(arma::vec& GVec,
                       double altFreq,
                       std::vector<uint32_t> indexForMissing,
                       std::string impute_method);

bool imputeGenoAndFlip(arma::vec& GVec,
                       double altFreq,
                       std::vector<uint32_t> indexForMissing,
                       double missingRate,
                       std::string impute_method,
                       const std::string& method = "foo");

arma::vec getTime();

void printTime(arma::vec t1, arma::vec t2, std::string message);
void printTimeDiff(arma::vec timeDiff, std::string message);

double getinvStd(double freq);

arma::vec nb(unsigned int n);

void imputeGeno(arma::vec& GVec,
                double altFreq,
                std::vector<uint32_t> indexForMissing,
                std::string imputeMethod);

void gethwepval(arma::vec GVec,
                double& hwepval,
                double hwepvalCutoff);

#endif
