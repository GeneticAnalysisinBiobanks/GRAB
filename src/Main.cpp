
// This includes the main codes to connect C++ and R

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
// std::this_thread::sleep_for (std::chrono::seconds(1));
#include <cstdio>         // std::remove

// Currently, omp does not work well, will check it later
// error: SET_VECTOR_ELT() can only be applied to a 'list', not a 'character'
// remove all Rcpp::List to check if it works
// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]]

#include "Main.hpp"
#include "PLINK.hpp"
#include "BGEN.hpp"
#include "VCF.hpp"
#include "POLMM.hpp"
#include "UTIL.hpp"
#include "SPACox.hpp"
#include "DenseGRM.hpp"

// global objects for different genotype formats

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static BGEN::BgenClass* ptr_gBGENobj = NULL;
// static VCF::VcfClass* ptr_gVCFobj = NULL;

// global objects for different analysis methods
static POLMM::POLMMClass* ptr_gPOLMMobj = NULL;
static SPACox::SPACoxClass* ptr_gSPACoxobj = NULL;
static DenseGRM::DenseGRMClass* ptr_gDenseGRMobj = NULL;

// global variables for analysis
std::string g_impute_method;      // "mean", "minor", or "drop"
double g_missingRate_cutoff;
unsigned int g_omp_num_threads;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
double g_region_minMAC_cutoff;    // for Rare Variants (RVs) whose MAC < this value, we aggregate these variants like SAIGE-GENE+ 
double g_region_maxMAF_cutoff;
unsigned int g_region_maxMarkers_cutoff;   // maximal number of markers in one chunk, only used for region-based analysis to reduce memory usage


// global variables for sparse GRM
arma::sp_mat g_SparseGRM;

// set up global variables for analysis

// arma::sp_mat getKinMatList(Rcpp::List KinMatListR)
// {
  // int nKin = KinMatListR.size();
  // Rcpp::CharacterVector NameKin = KinMatListR.names();
  // Rcpp::List KinMatList_sp;
  // for(int i = 0; i < nKin; i ++){
  //   string excludeChr = string(NameKin[i]);
  //   Rcpp::List KinMatTemp = KinMatListR[excludeChr];
  //   arma::umat locations = KinMatTemp["locations"];
  //   arma::vec values = KinMatTemp["values"];
  //   int n = KinMatTemp["nSubj"];
  //   // make a sparse matrix
  //   arma::sp_mat KinMat(locations, values, n, n);
  //   KinMatList_sp[excludeChr] = KinMat;
  // }
  // arma::umat locations = KinMatListR["locations"];
  // arma::vec values = KinMatListR["values"];
  // int n = KinMatListR["nSubj"];
  // // make a sparse matrix
  // arma::sp_mat KinMat(locations, values, n, n);
  // return KinMat;
// }

// [[Rcpp::export]]
void setSparseGRMInCPP(Rcpp::List t_KinMatListR)
{
  arma::umat locations = t_KinMatListR["locations"];
  arma::vec values = t_KinMatListR["values"];
  int n = t_KinMatListR["nSubj"];
  // make a sparse matrix
  arma::sp_mat KinMat(locations, values, n, n);
  g_SparseGRM = KinMat;
}

// [[Rcpp::export]]
void setDenseGRMInCPP(double t_memoryChunk,
                      double t_minMafGRM,
                      double t_maxMissingGRM)
{
  ptr_gDenseGRMobj = new DenseGRM::DenseGRMClass(ptr_gPLINKobj, t_memoryChunk, t_minMafGRM, t_maxMissingGRM);
}

// [[Rcpp::export]]
arma::vec getDenseGRMInCPP(arma::vec t_bVec,
                           std::string t_excludeChr, 
                           int t_grainSize)
{
  arma::vec yVec = DenseGRM::getKinbVec(t_bVec, ptr_gDenseGRMobj, t_excludeChr, t_grainSize);
  return yVec;
}

// [[Rcpp::export]]
void setMarker_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               unsigned int t_omp_num_threads)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
  g_omp_num_threads = t_omp_num_threads;
}

// [[Rcpp::export]]
void setRegion_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_max_maf_region,
                               double t_min_mac_region,
                               unsigned int t_max_markers_region,
                               unsigned int t_omp_num_threads)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_region_minMAC_cutoff = t_min_mac_region;
  g_region_maxMAF_cutoff = t_max_maf_region;
  g_region_maxMarkers_cutoff = t_max_markers_region;
  g_omp_num_threads = t_omp_num_threads;
}

//////// ---------- Main function for marker-level analysis --------- ////////////

// [[Rcpp::export]]
Rcpp::List mainMarkerInCPP(std::string t_method,       // "POLMM", "SPACox", "SAIGE"
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint32_t> t_genoIndex)  
{
  int q = t_genoIndex.size();  // number of markers
  
  // set up output
  std::vector<std::string> markerVec(q);  // marker IDs
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q);    // allele counts of ALT allele.
  std::vector<double> missingRateVec(q);  // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> BetaVec(q);         // beta value for ALT allele
  std::vector<double> seBetaVec(q);       
  std::vector<double> pvalVec(q);
  std::vector<double> zScoreVec(q);

  // std::cout << "Totally " << g_omp_num_threads << " thread(s) were used for parallel computation." << std::endl;
  
  // loop for all markers
//   omp_set_dynamic(0);     // Explicitly disable dynamic teams
//   omp_set_num_threads(g_omp_num_threads); // Use 4 threads for all consecutive parallel regions
//   
// #pragma omp parallel
// {
  for(int i = 0; i < q; i++){
    
    if(i % 1000 == 0){
      std::cout << "Completed " << i << "/" << q << " markers in the chunk." << std::endl;
    }
    
    // information of marker
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    
    uint32_t gIndex = t_genoIndex.at(i);
    
    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
    
    int n = GVec.size();
    
    // e.g. 21:1000234:A:T
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;
    
    // record basic information for the marker
    markerVec.at(i) = marker;               // marker IDs
    infoVec.at(i) = info;    // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;         // allele frequencies of ALT allele, this is not always < 0.5.
    altCountsVec.at(i) = altCounts;         // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    
    // MAF and MAC are for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * n * (1 - missingRate);
    
    // Quality Control (QC) based on missing rate, MAF, and MAC
    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff)){
      BetaVec.at(i) = arma::datum::nan;         
      seBetaVec.at(i) = arma::datum::nan;       
      pvalVec.at(i) = arma::datum::nan;
      continue;
    }
    
    // if(missingRate != 0){
      // Function imputeGenoAndFlip is in UTIL.cpp
      // flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, g_impute_method);  // in UTIL.cpp
    flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, g_impute_method);
    // }
    
    // analysis results for single-marker
    double Beta, seBeta, pval, zScore;
    
    Unified_getMarkerPval(t_method, GVec, 
                          false, // bool t_isOnlyOutputNonZero, 
                          indexForNonZero, Beta, seBeta, pval, zScore, altFreq);
    
    BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true       
    seBetaVec.at(i) = seBeta;       
    pvalVec.at(i) = pval;
    zScoreVec.at(i) = zScore;
  }
// }

  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("altCountsVec") = altCountsVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                          Rcpp::Named("BetaVec") = BetaVec,
                                          Rcpp::Named("seBetaVec") = seBetaVec,
                                          Rcpp::Named("pvalVec") = pvalVec,
                                          Rcpp::Named("zScoreVec") = zScoreVec);
  
  return OutList;  
}

//////// ---------- Main function for region-level analysis --------- ////////////

// [[Rcpp::export]]
Rcpp::List mainRegionInCPP(std::string t_method,       // "POLMM", "SAIGE"
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint32_t> t_genoIndex,
                           std::string t_outputFile,
                           unsigned int t_n)           // sample size  
{
  unsigned int q = t_genoIndex.size();                 // number of markers (before QC) in one region
  
  // set up output (Ultra-Rare Variants, URV)
  std::vector<std::string> markerVec(q), markerURVVec(q);      // marker IDs
  std::vector<std::string> infoVec(q), infoURVVec(q);          // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q), altFreqURVVec(q);         // allele frequencies of the ALT allele, this is not always < 0.5.
  std::vector<double> MACVec(q), MACURVVec(q);
  std::vector<double> MAFVec(q), MAFURVVec(q);
  std::vector<double> missingRateVec(q), missingRateURVVec(q); // missing rate
  std::vector<double> BetaVec(q);            // beta value for ALT allele
  std::vector<double> seBetaVec(q);          // seBeta value
  std::vector<double> pval0Vec(q);           // p values from normal distribution approximation  // might be confused, is this needed?
  std::vector<double> pval1Vec(q);           // p values from more accurate methods including SPA and ER
  
  // std::vector<bool> passQCVec(q, false);     // false: does not pass QC; true: pass QC
  // std::vector<unsigned int> passRVVec(q, 0); // 0: does not pass QC; 1: pass QC, common variants; 2: pass QC, rare variants
  
  std::vector<double> StatVec(q);            // score statistics
  std::vector<double> adjPVec(q);            // adjusted p-values
  
  // example #1: (q = 999, m1 = 10) -> (nchunks = 100, m2 = 9)
  // example #2: (q = 1000, m1 = 10) -> (nchunks = 100, m2 = 10)
  // example #3: (q = 1001, m1 = 10) -> (nchunks = 101, m2 = 1)
  
  unsigned int m1 = g_region_maxMarkers_cutoff;  // number of markers in all chunks expect for the last chunk
  unsigned int nchunks = 0;
  // unsigned int nchunks = (q-1) / m1 + 1;         // number of chunks.  e.g. 42 / 3 = 14, 2 / 3 = 0.
  // unsigned int m2 = (q-1) % m1 + 1;              // number of markers in the last chunk. e.g. 42 % 3 = 0, 2 % 3 = 2.
  // unsigned int m3 = m1;                          // number of markers in the current chunk
  
  // std::cout << "In region-based analysis, all " << q << " markers are splitted into " << nchunks << " chunks." << std::endl;
  // std::cout << "Each chunk includes no more than " << m1 << " markers." << std::endl;
  // std::cout << "The last chunk includes " << m2 << " markers." << std::endl;
  
  // Suppose that 
  // n is the sample size in analysis 
  // m (<q) is the number of markers that pass the marker-level QC (e.g., g_missingRate_cutoff and g_region_maxMAF_cutoff)
  // added on 08-10-2021: we collapse all ultra-rare variants to get one "fake" marker (+1)
  // VarMat (m+1 x m+1) is the variance matrix of these m markers
  // VarMat = P1Mat %*% P2Mat, where P1Mat is of (m+1 x n) and P2Mat is of (n x m+1)
  
  // std::vector<unsigned int> mPassQCVec, mPassRVVec, mPassCVVec;
  // unsigned int mPassQCTot, mPassRVTot, mPassCVTot;
  std::vector<unsigned int> mPassCVVec;
    
  // arma::vec GVecBurden(t_n, arma::fill::zeros);
  // arma::vec GVecURV(t_n, arma::fill::zeros);
  
  // arma::sp_mat P1Mat, P2Mat;
  // arma::mat P1Mat_DNS, P2Mat_DNS;
  arma::mat P1Mat, P2Mat;
  arma::vec GVecURV(t_n, arma::fill::zeros);    // aggregate ultra-rare variants (URV) whose MAC less than cutoff (g_region_minMAC_cutoff)
  
  // conduct marker-level analysis
  double Stat, Beta, seBeta, pval0, pval1;
  arma::vec P1Vec(t_n), P2Vec(t_n);
  
  // cycle for multiple chunks
  // for(unsigned int ichunk = 0; ichunk < nchunks; ichunk++)
  unsigned int ichunk = 0;
  unsigned int indexInChunk = 0;
  unsigned int i1 = 0;    // index of Markers (non-URV)
  unsigned int i2 = 0;    // index of Markers (Ultra-Rare Variants, URV)
  
  for(unsigned int i = 0; i < q; i++)
  {
    // std::cout << "Start analyzing chunk " << ichunk << "/" << nchunks - 1 << "." << std::endl;
    // if(ichunk == nchunks - 1) m3 = m2;  // number of markers in the last chunk
    // P1Mat.resize(m3, t_n);
    // P1Mat.resize(t_n, m3);
    // P2Mat.resize(t_n, m3);
    // arma::mat GMat(t_n, m3, arma::fill::zeros);
    // arma::mat GMatRV(t_n, m3, arma::fill::zeros);
    
    // std::vector<bool> passQCVecInChunk(m3, false);     // false: does not pass QC; true: pass QC
    // std::vector<unsigned int> passRVVecInChunk(m3, 0); // 0: does not pass QC; 1: pass QC, not ultra-rare variants; 2: pass QC, ultra-rare variants
    
    // cycle for markers in chunk
    // for(unsigned int i = 0; i < m3; i++)  // index of marker in the current chunk
    // {
    //   unsigned int i1 = i + m1 * ichunk;  // index of marker in all markers
      
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    
    // uint32_t gIndex = t_genoIndex.at(i1);
    uint32_t gIndex = t_genoIndex.at(i);
    
    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
    
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;
    
      // if(missingRate != 0){
        // Function imputeGenoAndFlip is in UTIL.cpp
        // flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, g_impute_method);  // in UTIL.cpp
    flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, g_impute_method);
      // }
      
      // MAF for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
    
    if((missingRate > g_missingRate_cutoff) || (MAF > g_region_maxMAF_cutoff) || MAF == 0){
      continue;  // does not pass QC
    }

    if(MAC > g_region_minMAC_cutoff){  // not Ultra-Rare Variants
      
      if(indexInChunk == 0){
        std::cout << "Start analyzing chunk " << ichunk << "....." << std::endl;
        // P1Mat.resize(t_n, m1);
        P1Mat.resize(m1, t_n);
        P2Mat.resize(t_n, m1);
        // arma::mat GMat(t_n, m1, arma::fill::zeros);
        // GMatURV.resize(t_n, m1);
      }
      
      markerVec.at(i1) = marker;             // marker IDs
      infoVec.at(i1) = info;                 // marker information: CHR:POS:REF:ALT
      altFreqVec.at(i1) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
      missingRateVec.at(i1) = missingRate;
      MACVec.at(i1) = MAC;
      MAFVec.at(i1) = MAF;
      
      Unified_getRegionPVec(t_method, GVec, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
      
      std::cout << "MAC:\t" << MAC << std::endl;
      std::cout << "pval0:\t" << pval0 << std::endl;
      std::cout << "pval1:\t" << pval1 << std::endl;
      
      // insert results to pre-setup vectors and matrix
      StatVec.at(i1) = Stat;        
      BetaVec.at(i1) = Beta * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true       
      seBetaVec.at(i1) = seBeta;       
      pval0Vec.at(i1) = pval0;
      pval1Vec.at(i1) = pval1;
      adjPVec.at(i1) = pval1;
      // P1Mat.row(i) = P1Vec.t();
      // P1Mat.col(indexInChunk) = P1Vec;
      P1Mat.row(indexInChunk) = P1Vec.t();
      P2Mat.col(indexInChunk) = P2Vec;
      
      i1 += 1;
      indexInChunk += 1;
    }else{   // Ultra-Rare Variants (URV)
      markerURVVec.at(i2) = marker;             // marker IDs
      infoURVVec.at(i2) = info;                 // marker information: CHR:POS:REF:ALT
      altFreqURVVec.at(i2) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
      missingRateURVVec.at(i2) = missingRate;
      MACURVVec.at(i2) = MAC;
      MAFURVVec.at(i2) = MAF;
      
      if(altFreq < 0.5){
        GVecURV += GVec;
      }else{
        GVecURV += 2 - GVec;
      }
        
      i2 += 1;
    }
    
    if(indexInChunk == m1){
      std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      std::cout << "P1Mat.n_rows:\t" << P1Mat.n_rows << std::endl;
      std::cout << "P1Mat.n_cols:\t" << P1Mat.n_cols << std::endl;
      std::cout << "P2Mat.n_rows:\t" << P2Mat.n_rows << std::endl;
      std::cout << "P2Mat.n_cols:\t" << P2Mat.n_cols << std::endl;
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      mPassCVVec.push_back(m1);
      ichunk += 1;
      indexInChunk = 0;
    }
  }
  
  if(i1 == 0){
    std::cout << "Only ultra-rare variants are found. This region will be skipped." << std::endl;
    Rcpp::List OutList = Rcpp::List::create();
    return OutList;
  }
    
  
  nchunks = ichunk + 1;
  arma::mat VarMat(i1, i1);
  
  // non Ultra Rare Variants
  markerVec.resize(i1);
  infoVec.resize(i1);                 // marker information: CHR:POS:REF:ALT
  altFreqVec.resize(i1);           // allele frequencies of ALT allele, this is not always < 0.5.
  missingRateVec.resize(i1);
  MACVec.resize(i1);
  MAFVec.resize(i1);
  StatVec.resize(i1);        
  BetaVec.resize(i1);  // Beta if flip = false, -1 * Beta is flip = true       
  seBetaVec.resize(i1);       
  pval0Vec.resize(i1);
  pval1Vec.resize(i1);
  adjPVec.resize(i1);
  
  // Ultra Rare Variants
  markerURVVec.resize(i2);             // marker IDs
  infoURVVec.resize(i2);                 // marker information: CHR:POS:REF:ALT
  altFreqURVVec.resize(i2);           // allele frequencies of ALT allele, this is not always < 0.5.
  missingRateURVVec.resize(i2);
  MACURVVec.resize(i2);
  MAFURVVec.resize(i2);
  
  if(i2 != 0){
    Unified_getRegionPVec(t_method, GVecURV, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
    StatVec.push_back(Stat);
    adjPVec.push_back(pval1);
    double minMAF_cutoff = g_region_minMAC_cutoff / (2 * t_n);
    MAFVec.push_back(minMAF_cutoff);
    
    // P1Mat.insert_rows(mPassCVVec.back(), P1Vec.t());
    // P2Mat.insert_cols(mPassCVVec.back(), P2Vec);
    P1Mat.row(indexInChunk) = P1Vec.t();
    P2Mat.col(indexInChunk) = P2Vec;
    
    indexInChunk += 1;
    VarMat.resize(i1 + 1, i1 + 1);
    // mPassCVTot += 1;
    // mPassCVVec.at(nchunks-1) += 1;
  }
  
  mPassCVVec.push_back(indexInChunk);

  if(indexInChunk != 0){
    P1Mat = P1Mat.rows(0, indexInChunk - 1);
    P2Mat = P2Mat.cols(0, indexInChunk - 1);
    if(nchunks != 1){
      std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      std::cout << "P1Mat.n_rows:\t" << P1Mat.n_rows << std::endl;
      std::cout << "P1Mat.n_cols:\t" << P1Mat.n_cols << std::endl;
      std::cout << "P2Mat.n_rows:\t" << P2Mat.n_rows << std::endl;
      std::cout << "P2Mat.n_cols:\t" << P2Mat.n_cols << std::endl;
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
    }
  }
  
  
  
    // std::cout << "MAC:\t" << MAC << std::endl;
    // std::cout << "altCounts:\t" << altCounts << std::endl;
    // std::cout << "g_region_minMAC_cutoff:\t" << g_region_minMAC_cutoff << std::endl;
    
    // Quality Control (QC) based on missing rate, MAF, and MAC
    // if((missingRate > g_missingRate_cutoff) || (MAF > g_region_maxMAF_cutoff) || MAF == 0){
    //   continue;
    // }
    
    // passQCVec.at(i1) = true;
    // passQCVecInChunk.at(i) = true;
    
    // GMat.col(i) = GVec;
    // if(MAC > g_region_minMAC_cutoff){
    //   passRVVec.at(i1) = 1;
    //   passRVVecInChunk.at(i) = 1;
    // }else{
    //   // std::cout << "i1:\t" << i1 << std::endl; 
    //   passRVVec.at(i1) = 2;
    //   passRVVecInChunk.at(i) = 2;
    //   continue;
    // }
      
    // The below function should be the most important one in region-based analysis
      
      
      
    // }
    
    // index vector for markers passing QC
    // arma::uvec indexQCVecInChunk(m3);
    // arma::uvec indexRVVecInChunk(m3);
    // arma::uvec indexCVVecInChunk(m3);
    // 
    // unsigned int mPassQCInChunk = 0;
    // unsigned int mPassRVInChunk = 0;
    // unsigned int mPassCVInChunk = 0;
    // 
    // for(unsigned int i = 0; i < m3; i++)
    // {
    //   if(passQCVecInChunk.at(i)){   // index of markers that pass QC: output marker-level information
    //     indexQCVecInChunk.at(mPassQCInChunk) = i;
    //     mPassQCInChunk ++;
    //   }
    //   
    //   if(passRVVecInChunk.at(i) == 2){   // index of markers that pass QC & MAC < cutoff (to be aggregated to one genotype)
    //     indexRVVecInChunk.at(mPassRVInChunk) = i;
    //     mPassRVInChunk ++;
    //   }
    //   
    //   if(passRVVecInChunk.at(i) == 1){   // index of markers that pass QC & MAC > cutoff (to be used for gene-based analysis)
    //     indexCVVecInChunk.at(mPassCVInChunk) = i;
    //     mPassCVInChunk ++;
    //   }
    // }
    // 
    // indexQCVecInChunk.resize(mPassQCInChunk);
    // indexRVVecInChunk.resize(mPassRVInChunk);
    // indexCVVecInChunk.resize(mPassCVInChunk);
    // 
    // mPassQCVec.push_back(mPassQCInChunk);
    // mPassRVVec.push_back(mPassRVInChunk);
    // mPassCVVec.push_back(mPassCVInChunk);
    
    // std::cout << "In chunk " << ichunk << ", totally " << mPassQCInChunk << " markers pass QC." << std::endl;
    // std::cout << "Of which, " << mPassRVInChunk << " variants are ultra-rare variants with MAC < " << g_region_minMAC_cutoff << std::endl;
    // std::cout << "and " << mPassCVInChunk << " variants are not ultra-rare variants." << std::endl;
    
    // P1Mat = P1Mat.rows(indexQCVecInChunk);
    // P2Mat = P2Mat.cols(indexQCVecInChunk);
    // StatVec = StatVec.rows(indexCVVecInChunk);
    
    // std::cout << "before submitting, P1Mat.n_rows:\t" << P1Mat.n_rows << std::endl;
    // std::cout << "before submitting, P1Mat.n_cols:\t" << P1Mat.n_cols << std::endl;
    // 
    // // P1Mat = P1Mat.rows(indexCVVecInChunk);
    // P1Mat = P1Mat.cols(indexCVVecInChunk);
    // P2Mat = P2Mat.cols(indexCVVecInChunk);
    // 
    // P1Mat_DNS = (arma::mat)P1Mat;
    // P1Mat_DNS = P1Mat_DNS.t();
    // P2Mat_DNS = (arma::mat)P2Mat;
    // 
    // std::cout << "after submitting, P1Mat.n_rows:\t" << P1Mat.n_rows << std::endl;
    // std::cout << "after submitting, P1Mat.n_cols:\t" << P1Mat.n_cols << std::endl;
    // 
    // arma::mat GMatRV = GMat.cols(indexRVVecInChunk);
    // arma::vec GVecURVInChunk = arma::sum(GMatRV, 1);
    // GVecURV += GVecURVInChunk;
    // GVecBurden += arma::sum(GMat, 1);
    
  //   if(ichunk == nchunks - 1){ // in the last chunk
  //     
  //     // remove markers that did not pass QC
  //     unsigned int tempIndex = 0;
  //     // unsigned int tempIndex1 = 0;
  //     unsigned int tempIndex2 = 0;
  //     for(unsigned int i = 0; i < q; i++)
  //     {
  //       // if(!passQCVec.at(i)){
  //       if(passRVVec.at(i) != 1){   // edited on 08-10-2021
  //         markerVec.erase(markerVec.begin()+tempIndex);
  //         infoVec.erase(infoVec.begin()+tempIndex);
  //         altFreqVec.erase(altFreqVec.begin()+tempIndex);
  //         missingRateVec.erase(missingRateVec.begin()+tempIndex);
  //         // StatVec.erase(StatVec.begin()+tempIndex);
  //         BetaVec.erase(BetaVec.begin()+tempIndex);
  //         seBetaVec.erase(seBetaVec.begin()+tempIndex);
  //         pval0Vec.erase(pval0Vec.begin()+tempIndex);
  //         pval1Vec.erase(pval1Vec.begin()+tempIndex);
  //         // added on 08-10-2021
  //         MACVec.erase(MACVec.begin()+tempIndex);
  //         MAFVec.erase(MAFVec.begin()+tempIndex);
  //         
  //         StatVec.erase(StatVec.begin()+tempIndex);
  //         adjPVec.erase(adjPVec.begin()+tempIndex);
  //       }else{
  //         tempIndex++;
  //       }
  //       
  //       if(passRVVec.at(i) != 2){   // edited on 08-10-2021
  //         markerURVVec.erase(markerURVVec.begin()+tempIndex2);
  //         infoURVVec.erase(infoURVVec.begin()+tempIndex2);
  //         altFreqURVVec.erase(altFreqURVVec.begin()+tempIndex2);
  //         missingRateURVVec.erase(missingRateURVVec.begin()+tempIndex2);
  //         MACURVVec.erase(MACURVVec.begin()+tempIndex2);
  //         MAFURVVec.erase(MAFURVVec.begin()+tempIndex2);
  //       }else{
  //         tempIndex2++;
  //       }
  //       // std::cout << "StatVec:\t" << StatVec.size() << std::endl;
  //       // std::cout << "passRVVec.at(i):\t" << passRVVec.at(i) << std::endl; 
  //       // std::cout << "passRVVec:\t" << passRVVec.size() << std::endl;
  //       
  //       // if(passRVVec.at(i) != 1){
  //       //   StatVec.erase(StatVec.begin()+tempIndex1);
  //       //   adjPVec.erase(adjPVec.begin()+tempIndex1);
  //       // }else{
  //       //   tempIndex1++;
  //       // }
  //       
  //       // std::cout << "StatVec:\t" << StatVec.size() << std::endl;
  //     }
  //     
  //     mPassQCTot = std::accumulate(mPassQCVec.begin(), mPassQCVec.end(), 0);
  //     mPassCVTot = std::accumulate(mPassCVVec.begin(), mPassCVVec.end(), 0);
  //     mPassRVTot = std::accumulate(mPassRVVec.begin(), mPassRVVec.end(), 0);
  //     
  //     std::cout << "mPassQCTot:\t" << mPassQCTot << std::endl;
  //     std::cout << "mPassCVTot:\t" << mPassCVTot << std::endl;
  //     std::cout << "mPassRVTot:\t" << mPassRVTot << std::endl;
  //     // std::cout << "mPassCVVec:\t" << mPassCVVec << std::endl;
  //     std::cout << "mPassCVVec.back():\t" << mPassCVVec.back() << std::endl;
  //     if(mPassRVTot != 0){
  //       Unified_getRegionPVec(t_method, GVecURV, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
  //       StatVec.push_back(Stat);
  //       adjPVec.push_back(pval1);
  //       double minMAF_cutoff = g_region_minMAC_cutoff / (2 * t_n);
  //       MAFVec.push_back(minMAF_cutoff);
  //       
  //       P1Mat_DNS.insert_rows(mPassCVVec.back(), P1Vec.t());
  //       P2Mat_DNS.insert_cols(mPassCVVec.back(), P2Vec);
  //       
  //       mPassCVTot += 1;
  //       mPassCVVec.at(nchunks-1) += 1;
  //     }
  //   }
  //   
  //   // save information to hard drive to avoid high memory usage
  //   // if((nchunks > 1) & (mPassQCInChunk != 0)){ 
  //   if((nchunks > 1) & (P1Mat_DNS.n_rows != 0)){ 
  //     // std::cout << "P1Mat.n_rows:\t" << P1Mat.n_rows << std::endl;
  //     // std::cout << "P1Mat.n_cols:\t" << P1Mat.n_cols << std::endl;
  //     // std::cout << "P2Mat.n_rows:\t" << P2Mat.n_rows << std::endl;
  //     // std::cout << "P2Mat.n_cols:\t" << P2Mat.n_cols << std::endl;
  //     
  //     P1Mat_DNS.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
  //     P2Mat_DNS.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
  //     // GVecURV.save(t_outputFile + "_GVecURV_Chunk_" + std::to_string(ichunk) + ".bin");
  //   }
  // }
  
  // std::cout << "P1Mat.n_rows:\t" << P1Mat.n_rows << std::endl;
  // std::cout << "P1Mat.n_cols:\t" << P1Mat.n_cols << std::endl;
  // std::cout << "P2Mat.n_rows:\t" << P2Mat.n_rows << std::endl;
  // std::cout << "P2Mat.n_cols:\t" << P2Mat.n_cols << std::endl;
  
  // calculate variance-covariance matrix VarMat = P1Mat %*% P2Mat
  // arma::mat VarMat(mPassCVTot+1, mPassCVTot+1);    // variance-covariance matrix (after QC)
  // arma::mat VarMat(mPassCVTot, mPassCVTot);
  // if(mPassRVTot == 0){
  //   VarMat.resize(mPassCVTot, mPassCVTot);
  // }else{
  //   Unified_getRegionPVec(t_method, GVecURV, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
  //   StatVec.push_back(Stat);
  //   adjPVec.push_back(pval1);
    

    // 
    // 
    // 
    // std::cout << "P1Mat.n_rows:\t" << P1Mat.n_rows << std::endl;
    // std::cout << "P1Mat.n_cols:\t" << P1Mat.n_cols << std::endl;
    // std::cout << "P2Mat.n_rows:\t" << P2Mat.n_rows << std::endl;
    // std::cout << "P2Mat.n_cols:\t" << P2Mat.n_cols << std::endl;
    
    // if(nchunks > 1){ 
    //   P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(nchunks-1) + ".bin");
    //   P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(nchunks-1) + ".bin");
    //   
    //   std::cout << "mPassCVVec.size():\t" << mPassCVVec.size() << std::endl;
    //   std::cout << "nchunks:\t" << nchunks << std::endl;
    //   
    //   // mPassCVVec.at(nchunks-1) += 1;
    //   // GVecURV.save(t_outputFile + "_GVecURV_Chunk_" + std::to_string(ichunk) + ".bin");
    // }
  // }
  
  // not so many markers in the region, so all matrix is in memory
  if(nchunks == 1){
    VarMat = P1Mat * P2Mat;
  }
    

  // the region includes more markers than limitation, so P1Mat and P2Mat have been put in hard drive
  if(nchunks > 1)
  {
    int first_row = 0, first_col = 0, last_row = 0, last_col = 0;
    
    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      // last_row = first_row + mPassQCVec.at(index1) - 1;
      last_row = first_row + mPassCVVec.at(index1) - 1;
      
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      
      // P1Mat_DNS.load(P1MatFile);
      P1Mat.load(P1MatFile);
      
      // std::cout << "P1Mat.n_rows: " << P1Mat.n_rows << std::endl;
      // std::cout << "P1Mat.n_cols: " << P1Mat.n_cols << std::endl;
      
      if(P1Mat.n_cols == 0) continue;
      
      // off-diagonal sub-matrix
      for(unsigned int index2 = 0; index2 < index1; index2++)
      {
        std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index2 << "/" << nchunks - 1 << ")........" << std::endl;
        
        P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index2) + ".bin");
        
        // std::cout << "P2Mat.n_rows: " << P2Mat.n_rows << std::endl;
        // std::cout << "P2Mat.n_cols: " << P2Mat.n_cols << std::endl;
        
        if(P2Mat.n_cols == 0) continue;
        
        arma::mat offVarMat = P1Mat * P2Mat;
        
        // last_col = first_col + mPassQCVec.at(index2) - 1;
        last_col = first_col + mPassCVVec.at(index2) - 1;
        
        VarMat.submat(first_row, first_col, last_row, last_col) = offVarMat;
        VarMat.submat(first_col, first_row, last_col, last_row) = offVarMat.t();
        
        first_col = last_col + 1;
      }
      
      // diagonal sub-matrix
      // last_col = first_col + mPassQCVec.at(index1) - 1;
      last_col = first_col + mPassCVVec.at(index1) - 1;
      std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index1 << "/" << nchunks - 1 << ")........" << std::endl;
      P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin");
      
      // std::cout << "P2Mat.n_rows: " << P2Mat.n_rows << std::endl;
      // std::cout << "P2Mat.n_cols: " << P2Mat.n_cols << std::endl;
      
      arma::mat diagVarMat = P1Mat * P2Mat;
      
      VarMat.submat(first_row, first_col, last_row, last_col) = diagVarMat;
      
      first_row = last_row + 1;
      first_col = 0;
      Rcpp::checkUserInterrupt();
      
    }
    
    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      std::string P2MatFile = t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin";
      const char* File1 = P1MatFile.c_str();
      const char* File2 = P2MatFile.c_str();
      // int temp1 = std::remove(File1);
      // int temp2 = std::remove(File2);
      std::remove(File1);
      std::remove(File2);
    }
  }
  
  // calculate p-values for the burden test
  
  // To be added later
  // double pval0Burden, pval1Burden;
  // Unified_getRegionPVec(t_method, GVec, Beta, seBeta, pval, P1Vec, P2Vec);
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("VarMat") = VarMat,
                                          Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("markerURVVec") = markerURVVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("infoURVVec") = infoURVVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("altFreqURVVec") = altFreqURVVec,
                                          Rcpp::Named("MAFVec") = MAFVec,
                                          Rcpp::Named("MAFURVVec") = MAFURVVec,
                                          Rcpp::Named("MACVec") = MACVec,
                                          Rcpp::Named("MACURVVec") = MACURVVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                          Rcpp::Named("missingRateURVVec") = missingRateURVVec,
                                          Rcpp::Named("StatVec") = StatVec,
                                          Rcpp::Named("BetaVec") = BetaVec,
                                          Rcpp::Named("seBetaVec") = seBetaVec,
                                          Rcpp::Named("pval0Vec") = pval0Vec,
                                          Rcpp::Named("pval1Vec") = pval1Vec,
                                          Rcpp::Named("adjPVec") = adjPVec);
                                          // Rcpp::Named("pval0Burden") = pval0Burden,
                                          // Rcpp::Named("pval1Burden") = pval1Burden);
  
  return OutList;
}


//////// ---------- Main function for genotype extraction --------- ////////////

// [[Rcpp::export]]
arma::mat getGenoInCPP(std::string t_genoType,
                       Rcpp::DataFrame t_markerInfo,
                       int n)
{
  int q = t_markerInfo.nrow();         // number of markers requested
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];
  arma::mat GMat(n, q);
  
  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;
  
  for(int i = 0; i < q; i++){
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType,    // "PLINK", "BGEN"
                                          gIndex,        // different meanings for different genoType
                                          ref,           // REF allele
                                          alt,           // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                          marker,        // marker ID extracted from genotype file
                                          pd,            // base position
                                          chr,           // chromosome
                                          altFreq,       // frequency of ALT allele
                                          altCounts,     // counts of ALT allele
                                          missingRate,   // missing rate
                                          imputeInfo,    // imputation information score, i.e., R2 (all 1 for PLINK)
                                          false,         // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,               // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    GMat.col(i) = GVec;
  }
  
  return GMat;
}

// [[Rcpp::export]]
arma::sp_mat getSpGenoInCPP(std::string t_genoType,
                            Rcpp::DataFrame t_markerInfo,
                            int n)
{
  int q = t_markerInfo.nrow();         // number of markers requested
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];
  arma::sp_mat GMat(n, q);             // change #1 compared to getGenoInCPP()
  
  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;
  
  for(int i = 0; i < q; i++){
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType,    // "PLINK", "BGEN"
                                          gIndex,        // different meanings for different genoType
                                          ref,           // REF allele
                                          alt,           // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                          marker,        // marker ID extracted from genotype file
                                          pd,            // base position
                                          chr,           // chromosome
                                          altFreq,       // frequency of ALT allele
                                          altCounts,     // counts of ALT allele
                                          missingRate,   // missing rate
                                          imputeInfo,    // imputation information score, i.e., R2 (all 1 for PLINK)
                                          false,         // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,               // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    GMat.col(i) = arma::sp_mat(GVec); // change #2 compared to getGenoInCPP()
  }
  
  return GMat;
}

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
                               std::vector<uint32_t>& t_indexForNonZero)     // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
{
  arma::vec GVec;
  if(t_genoType == "PLINK"){
    GVec = ptr_gPLINKobj->getOneMarker(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                       t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                       true);   // t_isTrueGenotype, only used for PLINK format.
  }
  
  if(t_genoType == "BGEN"){
    bool isBoolRead;
    GVec = ptr_gBGENobj->getOneMarker(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                      isBoolRead);
  }
  
  return GVec;
}

// a unified function to get marker-level p-value
void Unified_getMarkerPval(std::string t_method,   // "POLMM", "SPACox", "SAIGE"
                           arma::vec t_GVec,
                           bool t_isOnlyOutputNonZero,
                           std::vector<uint32_t> t_indexForNonZero,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
                           double& t_zScore,
                           double t_altFreq)
{
  if(t_method == "POLMM"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using POLMM method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
    
    ptr_gPOLMMobj->getMarkerPval(t_GVec, t_Beta, t_seBeta, t_pval, t_altFreq);
  }
  
  if(t_method == "SPACox"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SPACox method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
    t_pval = ptr_gSPACoxobj->getMarkerPval(t_GVec, t_altFreq, t_zScore);
  }
  
}

// a unified function to get marker-level information for region-level analysis

// Unified_getRegionPVec(t_method, GVec, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
void Unified_getRegionPVec(std::string t_method, 
                           arma::vec t_GVec, 
                           double& t_Stat,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval0, 
                           double& t_pval1,
                           arma::vec& t_P1Vec, 
                           arma::vec& t_P2Vec)
{
  // something to add
  if(t_method == "POLMM")
  {
    ptr_gPOLMMobj->getRegionPVec(t_GVec, t_Stat, t_Beta, t_seBeta, t_pval0, t_pval1, t_P1Vec, t_P2Vec);
  }
}

//////// ---------- Main functions to set objects for different genotype format --------- ////////////

// [[Rcpp::export]]
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> t_SampleInModel,
                      std::string t_AlleleOrder)
{
  ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
                                        t_famFile,
                                        t_bedFile,
                                        t_SampleInModel,
                                        t_AlleleOrder);
  
  int n = ptr_gPLINKobj->getN();
  std::cout << "Number of samples:\t" << n << std::endl;
  
}

// [[Rcpp::export]]
void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> t_SampleInBgen,
                     std::vector<std::string> t_SampleInModel,
                     bool t_isSparseDosageInBgen,
                     bool t_isDropmissingdosagesInBgen,
                     std::string t_AlleleOrder)
{
  ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName,
                                     t_bgenFileIndex,
                                     t_SampleInBgen,
                                     t_SampleInModel,
                                     t_isSparseDosageInBgen,
                                     t_isDropmissingdosagesInBgen,
                                     t_AlleleOrder);
  int n = ptr_gBGENobj->getN();
  std::cout << "Number of samples:\t" << n << std::endl;
}


//////// ---------- Main functions to set objects for different analysis methods --------- ////////////

// [[Rcpp::export]]
void setPOLMMobjInCPP(arma::mat t_muMat,
                      arma::mat t_iRMat,
                      arma::mat t_Cova,
                      arma::uvec t_yVec,
                      double t_tau,
                      bool t_printPCGInfo,
                      double t_tolPCG,
                      int t_maxiterPCG,
                      double t_varRatio, 
                      double t_SPA_cutoff,
                      bool t_flagSparseGRM)
{
  // check POLMM.cpp
  ptr_gPOLMMobj = new POLMM::POLMMClass(t_muMat,
                                        t_iRMat,
                                        t_Cova,
                                        t_yVec,
                                        g_SparseGRM,
                                        t_tau,
                                        t_printPCGInfo,
                                        t_tolPCG,
                                        t_maxiterPCG,
                                        t_varRatio, 
                                        t_SPA_cutoff,
                                        t_flagSparseGRM);
}

// [[Rcpp::export]]
Rcpp::List setPOLMMobjInCPP_NULL(bool t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
                                 arma::mat t_Cova,
                                 arma::uvec t_yVec,     // should be from 0 to J-1
                                 arma::vec t_beta,
                                 arma::vec t_bVec,
                                 arma::vec t_eps,           // 
                                 double t_tau,
                                 Rcpp::List t_SPmatR,    // output of makeSPmatR()
                                 Rcpp::List t_controlList)
{
  // arma::umat locations = t_SPmatR["locations"];
  // arma::vec values = t_SPmatR["values"];
  // std::cout << "Setting Sparse GRM...." << std::endl;
  // arma::sp_mat SparseGRM = arma::sp_mat(locations, values);
  
  // std::cout << "test, t_flagSparseGRM" << t_flagSparseGRM << std::endl;
  
  // The following function is in POLMM.cpp
  ptr_gPOLMMobj = new POLMM::POLMMClass(t_flagSparseGRM,       // if 1, then use Sparse GRM, otherwise, use Dense GRM
                                        ptr_gDenseGRMobj,
                                        ptr_gPLINKobj,
                                        t_Cova,
                                        t_yVec,     // should be from 0 to J-1
                                        t_beta,
                                        t_bVec,
                                        t_eps,           // 
                                        t_tau,
                                        g_SparseGRM,    // results of function setSparseGRMInCPP()
                                        t_controlList);
  
  ptr_gPOLMMobj->fitPOLMM();
  
  Rcpp::List outList = ptr_gPOLMMobj->getPOLMM();
  
  // ptr_gDenseGRMobj->closeDenseGRMObj();
  return outList;
}


// [[Rcpp::export]]
void setSPACoxobjInCPP(arma::mat t_cumul,
                       arma::vec t_mresid,
                       arma::mat t_XinvXX,
                       arma::mat t_tX,
                       int t_N,
                       double t_pVal_covaAdj_Cutoff,
                       double t_SPA_Cutoff)
{
  ptr_gSPACoxobj = new SPACox::SPACoxClass(t_cumul,
                                           t_mresid,
                                           t_XinvXX,
                                           t_tX,
                                           t_N,
                                           t_pVal_covaAdj_Cutoff,
                                           t_SPA_Cutoff);
}

