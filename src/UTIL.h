
#ifndef UTIL_H
#define UTIL_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <sys/time.h>

double getWeights(std::string t_kernel, 
                  double t_freq, 
                  arma::vec t_wBeta);

void imputeGeno(arma::vec& GVec, 
                double freq, 
                std::vector<uint32_t> posMissingGeno);

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat);

// duplicate each element for (J-1) times: n x 1 -> n(J-1) x 1 
arma::vec Vec2LongVec(arma::vec t_xVec, int n, int J);

// sum up each (J-1) elements: n(J-1) x 1 -> n x 1
arma::vec LongVec2Vec(arma::vec t_xVec, int n, int J);

// convert: n(J-1) x 1 -> n x (J-1) 
arma::mat Vec2Mat(arma::vec xVec, int n, int J);

// convert: n x (J-1) -> n(J-1) x 1
arma::vec Mat2Vec(arma::mat xMat, int n, int J);

arma::mat sumCols(arma::mat t_xMat, int J);

arma::vec getRPsiR(arma::mat t_muMat, arma::mat t_iRMat, int t_n, int t_J, int t_p); 

bool imputeGenoAndFlip(arma::vec& t_GVec, 
                       double t_altFreq, 
                       std::vector<uint32_t> t_indexForMissing,
                       std::string t_impute_method);

bool imputeGenoAndFlip(arma::vec& t_GVec, 
                       double t_altFreq, 
                       std::vector<uint32_t> t_indexForMissing,
                       double missingRate,
                       std::string t_impute_method,
                       const std::string& t_method = "foo");

arma::vec getTime();

void printTime(arma::vec t1, arma::vec t2, std::string message);
void printTimeDiff(arma::vec t_timeDiff, std::string t_message);

double getinvStd(double t_freq);

// http://thecoatlessprofessor.com/programming/set_rs_seed_in_rcpp_sequential_case/
// void set_seed(unsigned int seed) {
//   Rcpp::Environment base_env("package:base");
//   Rcpp::Function set_seed_r = base_env["set.seed"];
//   set_seed_r(seed);  
// };

// 
// arma::vec nb(int n){
//   return(Rcpp::rbinom(n,1,0.5));
// }

arma::vec nb(unsigned int n);

void imputeGeno(arma::vec& t_GVec, 
                double t_altFreq, 
                std::vector<uint32_t> t_indexForMissing,
                std::string t_imputeMethod);

void gethwepval(arma::vec t_GVec,
                double& t_hwepval,
                double t_hwepvalCutoff);

#endif

