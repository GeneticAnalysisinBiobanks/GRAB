
// This includes the main codes to connect C++ and R

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "Main.hpp"
#include "PLINK.hpp"
#include "BGEN.hpp"
#include "VCF.hpp"
#include "POLMM.hpp"
#include "UTIL.hpp"
#include "SPACox.hpp"

// global objects for different genotype formats

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static BGEN::BgenClass* ptr_gBGENobj = NULL;
// static VCF::VcfClass* ptr_gVCFobj = NULL;

// global objects for different analysis methods
static POLMM::POLMMClass* ptr_gPOLMMobj = NULL;
static SPACox::SPACoxClass* ptr_gSPACoxobj = NULL;

// global varables for marker-level analysis
uint8_t g_indexForImpute;   // 0: "mean"; 1: "minor"; 2: "drop"
double g_missingRate_cutoff;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
double g_region_maxMAF_cutoff;
unsigned int g_region_maxMarkers_cutoff;   // maximal number of markers in one chunk, only used for region-based analysis

void setMarker_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker)
{
  if(t_impute_method == "mean")
    g_indexForImpute = 0;
  
  if(t_impute_method == "minor")
    g_indexForImpute = 1;
  
  if(t_impute_method == "drop")
    g_indexForImpute = 2;
  
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
}

void setRegion_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_max_maf_region,
                               unsigned int t_max_markers_region)
{
  if(t_impute_method == "mean")
    g_indexForImpute = 0;
  
  if(t_impute_method == "minor")
    g_indexForImpute = 1;
  
  if(t_impute_method == "drop")
    g_indexForImpute = 2;
  
  g_missingRate_cutoff = t_missing_cutoff;
  g_region_maxMAF_cutoff = t_max_maf_region;
  g_region_maxMarkers_cutoff = t_max_markers_region;
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
                           double t_altFreq)
{
  if(t_method == "POLMM"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using POLMM method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
      
    ptr_gPOLMMobj->getMarkerPval(t_GVec, t_Beta, t_seBeta, t_pval, t_altFreq);
  }
    
  if(t_method == "SPACox"){
    // ptr_gSPACoxobj->getMarkerPval(t_GVec, t_Beta, t_seBeta, t_pval, t_altFreq);
  }
  
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
  std::vector<double> altFreqVec(q);         // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> missingRateVec(q);         // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> BetaVec(q);         // beta value for ALT allele
  std::vector<double> seBetaVec(q);       
  std::vector<double> pvalVec(q);
  
  
  // loop for all markers
  int indexPassingQC = 0;
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
    
    if(missingRate != 0){
      flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, g_indexForImpute);  // in UTIL.cpp
    }
    
    // analysis results for single-marker
    double Beta, seBeta, pval;
    
    Unified_getMarkerPval(t_method, GVec, 
                          false, // bool t_isOnlyOutputNonZero, 
                          indexForNonZero, Beta, seBeta, pval, altFreq);
    
    BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true       
    seBetaVec.at(i) = seBeta;       
    pvalVec.at(i) = pval;
  }
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                          Rcpp::Named("BetaVec") = BetaVec,
                                          Rcpp::Named("seBetaVec") = seBetaVec,
                                          Rcpp::Named("pvalVec") = pvalVec);
  
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
    arma::vec GVec = Unified_getOneMarker(t_genoType,   // "PLINK", "BGEN"
                                          gIndex,        // different meanings for different genoType
                                          ref,       // REF allele
                                          alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                          marker,    // marker ID extracted from genotype file
                                          pd,           // base position
                                          chr,       // chromosome
                                          altFreq,        // frequency of ALT allele
                                          altCounts,      // counts of ALT allele
                                          missingRate,    // missing rate
                                          imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                                          false,               // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,                   // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);     // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    GMat.col(i) = GVec;
  }
  
  return GMat;
}

//////// ---------- Main functions to set objects for different genotype format --------- ////////////

// [[Rcpp::export]]
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> t_SampleInModel)
{
  ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile,
                                        t_famFile,
                                        t_bedFile,
                                        t_SampleInModel);
  
  int n = ptr_gPLINKobj->getN();
  std::cout << "Number of samples:\t" << n << std::endl;
  
}

// [[Rcpp::export]]
void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> t_SampleInBgen,
                     std::vector<std::string> t_SampleInModel,
                     bool t_isSparseDosageInBgen,
                     bool t_isDropmissingdosagesInBgen)
{
  ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName,
                                     t_bgenFileIndex,
                                     t_SampleInBgen,
                                     t_SampleInModel,
                                     t_isSparseDosageInBgen,
                                     t_isDropmissingdosagesInBgen);
  int n = ptr_gBGENobj->getN();
  std::cout << "Number of samples:\t" << n << std::endl;
}


//////// ---------- Main functions to set objects for different analysis methods --------- ////////////

// [[Rcpp::export]]
void setPOLMMobjInCPP(arma::mat t_muMat,
                      arma::mat t_iRMat,
                      arma::mat t_Cova,
                      arma::uvec t_yVec,
                      Rcpp::List t_SPmatR,    // output of makeSPmatR()
                      double t_tau,
                      bool t_printPCGInfo,
                      double t_tolPCG,
                      int t_maxiterPCG,
                      double t_varRatio, 
                      double t_StdStat_cutoff)
{
  arma::umat locations = t_SPmatR["locations"];
  arma::vec values = t_SPmatR["values"];
  std::cout << "Setting Sparse GRM...." << std::endl;
  arma::sp_mat SparseGRM = arma::sp_mat(locations, values);
  ptr_gPOLMMobj = new POLMM::POLMMClass(t_muMat,
                                        t_iRMat,
                                        t_Cova,
                                        t_yVec,
                                        SparseGRM,
                                        t_tau,
                                        t_printPCGInfo,
                                        t_tolPCG,
                                        t_maxiterPCG,
                                        t_varRatio, 
                                        t_StdStat_cutoff);
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

//////// ---------- Main function for region-level analysis --------- ////////////

// [[Rcpp::export]]
Rcpp::List mainRegionInCPP(std::string t_method,       // "POLMM", "SAIGE"
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint32_t> t_genoIndex,
                           std::string t_outputFile)  
{
  int q = t_genoIndex.size();  // number of markers in one region
  
  // set up output
  std::vector<std::string> markerVec(q);    // marker IDs
  std::vector<std::string> infoVec(q);      // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q);        // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> missingRateVec(q);    // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> BetaVec(q);           // beta value for ALT allele
  std::vector<double> seBetaVec(q);         // seBeta value
  std::vector<double> pvalVec(q);           // p values
  
  unsigned int m1 = g_region_maxMarkers_cutoff;  // number of markers in all chunks expect for the last chunk (1, 2, 3, ...)
  unsigned int nchunks = (q-1) / m1;             // number of chunks (0, 1, 2, ...).  e.g. 42 / 3 = 14, 2 / 3 = 0.
  unsigned int m2 = (q-1) % m1;                  // number of markers in the last chunk, (0, 1, 2, ...). e.g. 42 % 3 = 0, 2 % 3 = 2.
  unsigned int m3 = m1 - 1;                      // use "-1" to start from 0
  
  std::cout << "In region-based analysis, all " << q << "markers are splitted into ", nchunks+1 << " chunks." << std::endl;
  
  // Suppose that 
  // n is sample size in analysis 
  // m (<q) is the number of markers that pass the marker-level filter (such as g_missingRate_cutoff and g_region_maxMAF_cutoff)
  // VarMat (m x m) is the variance matrix of these m markers
  // VarMat = P1Mat %*% P2Mat, where P1Mat is of (m x n) and P2Mat is of (n x m)
  
  // cycle for multiple chunks
  for(unsigned int ichunk = 0; ichunk < nchunks; ichunk++) // pay special attention to nchunks = 0
  {
    std::cout << "Start analyzing chunk " << ichunk << "/" << nchunks << "." << std::endl;
    
    if(ichunk == nchunks) m3 = m2;  // number of markers in one chunk (0, 1, 2, ...)
    arma::mat P1Mat(m3, n);
    arma::mat P2Mat(n, m3);
    std::vector<bool> passQCVecInChunk(m3, true);
    
    std::vector<unsigned int> mPassQCVec;
    // cycle for markers in each chunk
    for(unsigned int i = 0; i < m3; i++)
    {
      unsigned int i1 = i + m1 * ichunk;  // index in all markers
      
      // information of marker
      double altFreq, altCounts, missingRate, imputeInfo;
      std::vector<uint32_t> indexForMissing, indexForNonZero;
      std::string chr, ref, alt, marker;
      uint32_t pd;
      bool flip = false;
      
      uint32_t gIndex = t_genoIndex.at(i1);
      
      arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                            true, // bool t_isOutputIndexForMissing,
                                            indexForMissing,
                                            false, // bool t_isOnlyOutputNonZero,
                                            indexForNonZero);
      
      std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;
      
      if(missingRate != 0){
        flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, g_indexForImpute);  // in UTIL.cpp
      }
      
      // record basic information for the marker
      markerVec.at(i1) = marker;             // marker IDs
      infoVec.at(i1) = info;                 // marker information: CHR:POS:REF:ALT
      altFreqVec.at(i1) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
      missingRateVec.at(i1) = missingRate;
      
      // MAF and MAC are for Quality Control (QC)
      double MAF = std::min(altFreq, 1 - altFreq);
      double MAC = MAF * n * (1 - missingRate);
      
      // Quality Control (QC) based on missing rate, MAF, and MAC
      if((missingRate > g_missingRate_cutoff) || (MAF > g_region_maxMAF_cutoff)){
        BetaVec.at(i1) = arma::datum::nan;         
        seBetaVec.at(i1) = arma::datum::nan;       
        pvalVec.at(i1) = arma::datum::nan;
        passQCVecInChunk.at(i) = false;
        continue;
      }
      
      // analysis results for single-marker
      double Beta, seBeta, pval;
      
      Unified_getMarkerPval(t_method, GVec, 
                            false, // bool t_isOnlyOutputNonZero, 
                            indexForNonZero, Beta, seBeta, pval, altFreq);
      
      BetaVec.at(i1) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true       
      seBetaVec.at(i1) = seBeta;       
      pvalVec.at(i1) = pval;
      
      arma::vec P1Vec(n), P2Vec(n);
      
      // add one unified function to calculate P1Vec and P2Vec
      // Unified_getPVec(t_method, xxx, P1Vec, P2Vec);
      
      P1Mat.row(i) = P1Vec.t();
      P2Mat.col(i) = P2Vec;
      
      // index vector for markers passing QC
      std::vector<unsigned int> indexQCVecInChunk;
      for(unsigned int i = 0; i < m3; i++)
      {
        if(passQCVecInChunk.at(i)) 
          indexQCVecInChunk.push_back(i);
      }
      mPassQC = indexQCVecInChunk.size();
      mPassQCVec.push_back(mPassQC);
      
      P1Mat = P1Mat.rows(indexQCVecInChunk);
      P2Mat = P2Mat.cols(indexQCVecInChunk);
      
      // save information to hard-drive to avoid very high memory usage
      if((nchunks > 0) & (mPassQC != 0)){ 
        P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
        P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
        std::cout << "Completed chunk "<< ichunk << "!" << std::endl;
      }
      
    }
    
    // calculate variance-covariance matrix
    unsigned int mPassQCTot = std::accumulate(mPassQCVec.begin(), mPassQCVec.end(), 0);
    
    arma::mat VarMat(mPassQCTot, mPassQCTot);    // variance matrix (after adjusting for relatedness)
    
    // not so many markers in the region, so everything is in memory
    if(nchunks == 0)
      VarMat = P1Mat * P2Mat;
    
    // the region includes more markers than limitation, so everything is in hard-drive storage
    if(nchunks > 0)
    {
      int first_row = 0;
      int first_col = 0;
      int last_row = 0;
      int last_col = 0;
      // int last_row = t_maxMarkers - 1;
      // int last_col = t_maxMarkers - 1;
      
      for(int index1 = 0; index1 < nchunks; index1++)
      {
        last_row = first_row + mPassQCVec.at(index1) - 1;
        P1Mat.load(t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin");
        
        // off-diagonal sub-matrix
        for(int index2 = 0; index2 < index1; index2++)
        {
          std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index2 << "/" << nchunks - 1 << ")........" << std::endl;
          P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index2) + ".bin");
          arma::mat offVarMat = P1Mat * P2Mat;
          
          last_col = first_col + mPassQCVec.at(index2) - 1;
          
          VarMat.submat(first_row, first_col, last_row, last_col) = offVarMat;
          VarMat.submat(first_col, first_row, last_col, last_row) = offVarMat.t();
          
          first_col = last_col + 1;
        }
        
        // // diagonal sub-matrix
        last_col = first_col + mPassQCVec.at(index1) - 1;
        std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index1 << "/" << nchunks - 1 << ")........" << std::endl;
        P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin");
        arma::mat diagVarMat = P1Mat * P2Mat;
        
        VarMat.submat(first_row, first_col, last_row, last_col) = diagVarMat;
        
        first_row += mPassQCVec.at(index1);
        first_col = 0;
        Rcpp::checkUserInterrupt();
      }
    }
    
    Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("StatVec") = StatVec,
                                            Rcpp::Named("VarSMat") = VarSMat,
                                            Rcpp::Named("markerVec") = markerVec,
                                            Rcpp::Named("freqVec") = freqVec,
                                            Rcpp::Named("weightVec") = weightVec,
                                            Rcpp::Named("a1Vec") = a1Vec,
                                            Rcpp::Named("a2Vec") = a2Vec,
                                            Rcpp::Named("pdVec") = pdVec,
                                            Rcpp::Named("flipVec") = flipVec,
                                            Rcpp::Named("pvalNormVec") = pvalNormVec,
                                            Rcpp::Named("pvalVec") = pvalVec,
                                            Rcpp::Named("posVec") = posVec,    // starting from 0, not 1
                                            Rcpp::Named("rBT") = rBT);
    
  }
}







