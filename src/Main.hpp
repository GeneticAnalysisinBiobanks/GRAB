
#ifndef MAIN_HPP
#define MAIN_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

arma::fmat getSymmMat(arma::fmat& xMat1,    // n x p
                      arma::fmat& xMat2,    // n x p
                      int p);

arma::fmat getfmatMulti(arma::fmat& xMat1,  // n x p
                        arma::fmat& xMat2); // n x p

Rcpp::List MAIN_REGION(std::vector<std::string> t_MarkerReqstd,
                       double t_NonZero_cutoff,
                       double t_StdStat_cutoff,
                       int t_maxMarkers,
                       std::string t_outputFile,
                       double t_missingRate_cutoff,
                       double t_maxMAF_cutoff,
                       std::string t_kernel,
                       arma::vec t_wBeta);

Rcpp::List MAIN_MARKER(std::vector<std::string> t_MarkerReqstd,
                       double t_StdStat_cutoff,
                       double t_missingRate_cutoff,
                       double t_minMAF_cutoff,
                       int t_minMAC_cutoff,
                       double t_varRatio);


// a unified function to get single marker from genotype file
arma::vec Unified_getOneMarker(std::string t_genoType,   // "PLINK", "BGEN"
                               uint64_t t_gIndex,        // different meanings for different genoType
                               std::string& t_ref,       // REF allele
                               std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                               std::string& t_marker,    // marker ID extracted from genotype file
                               uint32_t& t_pd,           // base position
                               std::string& t_chr,       // chromosome
                               double& t_altFreq,        // frequency of ALT allele
                               double& t_altCounts,      // counts of ALT allele
                               double& t_missingRate,    // missing rate
                               double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                               bool t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                               std::vector<uint32_t>& t_indexForMissing,     // index of missing genotype data
                               bool t_isOnlyOutputNonZero,                   // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                               std::vector<uint32_t>& t_indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.

// a unified function to get marker-level p-value
void Unified_getMarkerPval(std::string t_method,   // "POLMM", "SPACox", "SAIGE"
                           arma::vec t_GVec,
                           bool t_isOnlyOutputNonZero,
                           std::vector<uint32_t> t_indexForNonZero,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
                           double& t_zScore,
                           double t_altFreq);

// updated on 2023-05-22: add option to output hwepvalue, mainly for SPAGRM
void Unified_getMarkerPval(std::string t_method,   // "POLMM", "SPACox", "SAIGE"
                           arma::vec t_GVec,
                           bool t_isOnlyOutputNonZero,
                           std::vector<uint32_t> t_indexForNonZero,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
                           double& t_zScore,
                           double t_altFreq,
                           double& t_hwepval,
                           double t_hwepvalCutoff); // default value = 0.1

void Unified_getRegionPVec(std::string t_method, 
                           arma::vec t_GVec, 
                           double& t_Stat,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval0, 
                           double& t_pval1,
                           arma::vec& t_P1Vec, 
                           arma::vec& t_P2Vec);
  
#endif
