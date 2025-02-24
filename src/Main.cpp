
// This file includes the main codes to connect C++ and R

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <boost/math/distributions/beta.hpp>

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
#include "SPAmix.hpp"
#include "SPAGRM.hpp"
#include "SAGELD.hpp"
#include "WtSPAG.hpp"

// global objects for different genotype formats

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static BGEN::BgenClass* ptr_gBGENobj = NULL;
// static VCF::VcfClass* ptr_gVCFobj = NULL;

static DenseGRM::DenseGRMClass* ptr_gDenseGRMobj = NULL;

// global objects for different analysis methods
static POLMM::POLMMClass* ptr_gPOLMMobj = NULL;
static SPACox::SPACoxClass* ptr_gSPACoxobj = NULL;
static SPAmix::SPAmixClass* ptr_gSPAmixobj = NULL;
static SPAGRM::SPAGRMClass* ptr_gSPAGRMobj = NULL;
static SAGELD::SAGELDClass* ptr_gSAGELDobj = NULL;
static WtSPAG::WtSPAGClass* ptr_gWtSPAGobj = NULL;

// global variables for analysis
std::string g_impute_method;      // "mean", "minor", or "drop"
double g_missingRate_cutoff;
unsigned int g_omp_num_threads;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
double g_region_minMAC_cutoff;    // for Rare Variants (RVs) whose MAC < this value, we aggregate these variants like SAIGE-GENE+ 
double g_region_maxMAF_cutoff;
unsigned int g_region_maxMarkers_cutoff;   // maximal number of markers in one chunk, only used for region-based analysis to reduce memory usage
arma::vec g_region_weight_beta;
arma::vec g_region_max_maf_vec;

// the below is only valid for POLMM
arma::uvec g_group;
bool g_ifOutGroup;
unsigned int g_nGroup;

// global variables for sparse GRM
arma::sp_mat g_SparseGRM;

// set up global variables for analysis
arma::vec g_compTime1(2, arma::fill::zeros);   // Unified_getOneMarker
arma::vec g_compTime2(2, arma::fill::zeros);   // Unified_getRegionPVec
arma::vec g_compTime3(2, arma::fill::zeros);

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
  if(ptr_gDenseGRMobj)
    delete ptr_gDenseGRMobj;
  
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
                               unsigned int t_omp_num_threads,
                               arma::uvec t_group,
                               bool t_ifOutGroup,
                               unsigned int t_nGroup)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
  g_omp_num_threads = t_omp_num_threads;
  g_group = t_group;
  g_ifOutGroup = t_ifOutGroup;
  g_nGroup = t_nGroup;
}

// [[Rcpp::export]]
void setRegion_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_max_maf_region,
                               double t_min_mac_region,
                               unsigned int t_max_markers_region,
                               unsigned int t_omp_num_threads,
                               arma::vec t_region_weight_beta,
                               arma::vec t_region_max_maf_vec)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_region_minMAC_cutoff = t_min_mac_region;
  g_region_maxMAF_cutoff = t_max_maf_region;
  g_region_maxMarkers_cutoff = t_max_markers_region;
  g_omp_num_threads = t_omp_num_threads;
  g_region_weight_beta = t_region_weight_beta;
  g_region_max_maf_vec = t_region_max_maf_vec;
}

void updateGroupInfo(arma::vec t_GVec,
                     std::vector<uint32_t> t_indexForMissing,
                     arma::vec& nSamplesInGroupVec,
                     arma::vec& AltCountsInGroupVec,
                     arma::vec& AltFreqInGroupVec)
{
  unsigned int n1 = t_GVec.size();
  nSamplesInGroupVec.zeros();
  AltCountsInGroupVec.zeros();
  
  // std::cout << "test1.21" << std::endl;
  // std::cout << "n1:\t" << n1 << std::endl;
  // std::cout << "t_indexForMissing.size()\t" << t_indexForMissing.size() << std::endl;
  
  // debug on 2023-04-21
  if(t_indexForMissing.size() == 0)
    t_indexForMissing.push_back(n1);
  
  // std::cout << "t_indexForMissing.size()\t" << t_indexForMissing.size() << std::endl;
    
  unsigned int i1 = 0;
  for(unsigned int i = 0; i < n1; i++){
    if(i == t_indexForMissing.at(i1)){
      if(i1 < t_indexForMissing.size() - 1)
        i1 ++;
    }else{
      unsigned int grp = g_group.at(i);
      
      // std::cout << "grp:\t" << grp << std::endl;
      // std::cout << "i:\t" << i << std::endl;
      // std::cout << "i1:\t" << i1 << std::endl;
      
      nSamplesInGroupVec.at(grp) += 1;
      AltCountsInGroupVec.at(grp) += t_GVec.at(i);
    }
  }
  
  // std::cout << "test1.22" << std::endl;
  
  AltFreqInGroupVec = AltCountsInGroupVec / nSamplesInGroupVec / 2;
}

//////// ---------- Main function for marker-level analysis --------- ////////////

// [[Rcpp::export]]
Rcpp::List mainMarkerInCPP(std::string t_method,       // "POLMM", "SPACox", "SAIGE", "SPAmix", "SPAGRM"
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint64_t> t_genoIndex)  
{
  int q = t_genoIndex.size();  // number of markers

  // set up output
  std::vector<std::string> markerVec(q);  // marker IDs
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q);    // allele counts of ALT allele.
  std::vector<double> missingRateVec(q);  // missing rate of markers
  // std::vector<double> BetaVec(q, arma::datum::nan);         // beta value for ALT allele
  // std::vector<double> seBetaVec(q, arma::datum::nan);
  std::vector<double> hwepvalVec(q, arma::datum::nan);
  
  int Npheno = 1;
  if(t_method == "SPAmix")
    Npheno = ptr_gSPAmixobj->getNpheno();
  
  if(t_method == "SAGELD")
    Npheno = 2;

  std::vector<double> pvalVec(q*Npheno, arma::datum::nan);
  std::vector<double> zScoreVec(q*Npheno, arma::datum::nan);
  std::vector<double> BetaVec(q*Npheno, arma::datum::nan);         // beta value for ALT allele
  std::vector<double> seBetaVec(q*Npheno, arma::datum::nan);  
  
  arma::mat nSamplesInGroup;
  arma::mat AltCountsInGroup;
  arma::mat AltFreqInGroup;
  
  if(g_ifOutGroup){
    nSamplesInGroup.resize(q, g_nGroup);
    AltCountsInGroup.resize(q, g_nGroup);
    AltFreqInGroup.resize(q, g_nGroup);
  }
  
  // std::cout << "Totally " << g_omp_num_threads << " thread(s) were used for parallel computation." << std::endl;
  
  // loop for all markers
  //   omp_set_dynamic(0);     // Explicitly disable dynamic teams
  //   omp_set_num_threads(g_omp_num_threads); // Use 4 threads for all consecutive parallel regions
  //   
  // #pragma omp parallel
  // {
  for(int i = 0; i < q; i++)
  {
    if(i % 1000 == 0)
    std::cout << "Completed " << i << "/" << q << " markers in the chunk." << std::endl;
    
    // information of marker
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    
    uint64_t gIndex = t_genoIndex.at(i);

    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
    int n = GVec.size();
    
    // std::cout << "marker:\t" << marker << std::endl;
    // std::cout << "test1.1" << std::endl;
    // std::cout << "GVec.size():\t" << n << std::endl;
    // std::cout << "g_ifOutGroup:\t" << g_ifOutGroup << std::endl;
    // std::cout << "g_nGroup:\t" << g_nGroup << std::endl;
    
    if(g_ifOutGroup){
      arma::vec nSamplesInGroupVec(g_nGroup);
      arma::vec AltCountsInGroupVec(g_nGroup);
      arma::vec AltFreqInGroupVec(g_nGroup);
      
      // std::cout << "test1.2" << std::endl;
      // std::cout << "g_nGroup:\t" << g_nGroup << std::endl;
      
      updateGroupInfo(GVec, indexForMissing, nSamplesInGroupVec, AltCountsInGroupVec, AltFreqInGroupVec);
      
      // std::cout << "test1.3" << std::endl;
      
      nSamplesInGroup.row(i) = nSamplesInGroupVec.t();
      AltCountsInGroup.row(i) = AltCountsInGroupVec.t();
      AltFreqInGroup.row(i) = AltFreqInGroupVec.t();
    }
    
    // std::cout << "test1.2" << std::endl;
    
    // e.g. 21:1000234:A:T
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;
    
    // record basic information for the marker
    markerVec.at(i) = marker;               // marker IDs
    infoVec.at(i) = info;                   // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;             // allele frequencies of ALT allele, this is not always < 0.5
    altCountsVec.at(i) = altCounts;         // allele counts of ALT allele
    missingRateVec.at(i) = missingRate;     // missing rate
    
    // c(MAF, MAC) for Quality Control (QC)
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = 2 * MAF * n * (1 - missingRate);
    
    // Quality Control (QC) based on missing rate, MAF, and MAC
    if((missingRate > g_missingRate_cutoff) || (MAF < g_marker_minMAF_cutoff) || (MAC < g_marker_minMAC_cutoff))
      continue;
    
    // Check UTIL.cpp
    flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, g_impute_method);
    
    // analysis results for single-marker
    double Beta, seBeta, pval, zScore;
    
    // std::cout << "test1.4" << std::endl;
    
    double hwepvalCutoff = 0.1;  // to be changed to a option, instead of a default value, later
    double hwepval = 0;
    
    if(t_method != "WtSPAG"){
      Unified_getMarkerPval(t_method, GVec,
                            false, // bool t_isOnlyOutputNonZero,
                            indexForNonZero, Beta, seBeta, pval, zScore, altFreq, hwepval, hwepvalCutoff);
    }
    
    if(t_method == "WtSPAG"){
      pval = ptr_gWtSPAGobj->getMarkerPval(GVec, altFreq, zScore, flip, i);
    }
    
    
    // std::cout << "test1.5" << std::endl;
    
    // BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true       
    // std::cout << "test1.51" << std::endl;
    // seBetaVec.at(i) = seBeta;       
    // std::cout << "test1.52" << std::endl;
    
    if(t_method == "SPAmix"){
      arma::vec pvalVecTemp = ptr_gSPAmixobj->getpvalVec();
      arma::vec zScoreVecTemp = ptr_gSPAmixobj->getzScoreVec();
      
      for(int j = 0; j < Npheno; j++){
        pvalVec.at(i*Npheno+j) = pvalVecTemp.at(j);
        zScoreVec.at(i*Npheno+j) = zScoreVecTemp.at(j);
      }
    }else if(t_method == "SAGELD"){
      arma::vec pvalVecTemp = ptr_gSAGELDobj->getpvalVec();
      arma::vec zScoreVecTemp = ptr_gSAGELDobj->getzScoreVec();
      arma::vec BetaVecTemp = ptr_gSAGELDobj->getBetaVec();
      arma::vec seBetaVecTemp = ptr_gSAGELDobj->getseBetaVec();
      
      for(int j = 0; j < 2; j++){
        pvalVec.at(2*i+j) = pvalVecTemp.at(j);
        zScoreVec.at(2*i+j) = zScoreVecTemp.at(j);
        BetaVec.at(2*i+j) = BetaVecTemp.at(j) * (1 - 2*flip);
        seBetaVec.at(2*i+j) = seBetaVecTemp.at(j);
      }
    }else{
      pvalVec.at(i) = pval;
      zScoreVec.at(i) = zScore; 
      BetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1*Beta is flip = true   
      seBetaVec.at(i) = seBeta;       
    }
    
    hwepvalVec.at(i) = hwepval;
  }
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("altCountsVec") = altCountsVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                          Rcpp::Named("pvalVec") = pvalVec,
                                          Rcpp::Named("beta") = BetaVec,
                                          Rcpp::Named("seBeta") = seBetaVec,
                                          Rcpp::Named("zScore") = zScoreVec,
                                          Rcpp::Named("nSamplesInGroup") = nSamplesInGroup,
                                          Rcpp::Named("AltCountsInGroup") = AltCountsInGroup,
                                          Rcpp::Named("AltFreqInGroup") = AltFreqInGroup,
                                          Rcpp::Named("hwepvalVec") = hwepvalVec);
  
  return OutList;  
}

//////// ---------- Main function for region-level analysis --------- ////////////

// [[Rcpp::export]]
Rcpp::List mainRegionURVInCPP(std::string t_method,       // "POLMM", "SPACox", "SAIGE" (to be continued)
                              std::string t_genoType,     // "PLINK", "BGEN"
                              std::vector<uint64_t> t_genoIndex,
                              unsigned int t_n)           // sample size
{
  unsigned int q = t_genoIndex.size();                 // number of Ultra-Rare Variants (URV) markers (after QC) in one region
  double Stat, Beta, seBeta, pval0, pval1;
  arma::vec P1Vec(t_n), P2Vec(t_n);
  
  arma::vec GVecURV(t_n, arma::fill::zeros);    // aggregate ultra-rare variants (URV) whose MAC less than cutoff (g_region_minMAC_cutoff)
  
  for(unsigned int i = 0; i < q; i++)
  {
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    
    uint64_t gIndex = t_genoIndex.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
    
    imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, g_impute_method);
    
    if(altFreq < 0.5){
      GVecURV = arma::max(GVecURV, GVec);  // edited on 2021-08-20
    }else{
      GVecURV = arma::max(GVecURV, 2 - GVec); // edited on 2021-08-20
    }
  }
  
  Unified_getRegionPVec(t_method, GVecURV, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("Stat") = Stat,
                                          Rcpp::Named("pval1") = pval1);
  
  return OutList;  
}

void getLabelInfo(arma::vec t_GVec,    // genotype vector after imputation
                  std::vector<unsigned int> t_labelVec,
                  unsigned int t_rowIndex,
                  // std::vector<uint32_t> t_indexForMissing,
                  arma::mat& t_MACLabelMat,
                  arma::mat& t_MAFLabelMat)
{
  unsigned int n = t_labelVec.size();
  unsigned int nLabel = t_MACLabelMat.n_cols;
  arma::vec MACLabelVec(nLabel, arma::fill::zeros);
  arma::vec counttLabelVec(nLabel, arma::fill::zeros);
  for(unsigned int i = 0; i < n; i++)
  {
    counttLabelVec.at(t_labelVec.at(i) - 1) += 1;
    MACLabelVec.at(t_labelVec.at(i) - 1) += t_GVec.at(i);
  }
  arma::vec MAFLabelVec = MACLabelVec / (counttLabelVec * 2);
  
  t_MACLabelMat.row(t_rowIndex) = MACLabelVec.t();
  t_MAFLabelMat.row(t_rowIndex) = MAFLabelVec.t();
}

// [[Rcpp::export]]
Rcpp::List mainRegionInCPP(std::string t_method,       // "POLMM", "SPACox", "SAIGE" (to be continued)
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint64_t> t_genoIndex,
                           std::vector<double> t_weightVec,
                           std::string t_outputFile,
                           std::vector<unsigned int> t_labelVec,
                           unsigned int t_nLabel,           // # 2022-04-27: give labels to each subject (e.g. 0 for control and 1 for case), to be extended later. Start from 0.
                           arma::mat t_annoMat,             // # 2022-05-01: matrix to indicate if the marker is in annotation list (0 or 1)
                           std::vector<std::string> t_annoVec)
{
  // arma::mat P1Mat,            // edited on 2021-08-19: to avoid repeated memory allocation of P1Mat and P2Mat
  // arma::mat P2Mat)
  unsigned int n = t_labelVec.size();
  unsigned int q = t_genoIndex.size();                 // number of markers (before QC) in one region
  unsigned int nAnno = t_annoMat.n_cols;               // 2022-05-01: number of annotations
  
  // +1 corresponds to "collapsing ultra-rare variants"
  arma::uvec indicatorVec(q+nAnno, arma::fill::zeros);       // 0: does not pass QC, 1: non-URV, 2: URV
  Rcpp::StringVector markerVec(q+nAnno);
  Rcpp::StringVector infoVec(q+nAnno);
  arma::vec altFreqVec(q+nAnno);         // allele frequencies of the ALT allele, this is not always < 0.5.
  arma::vec MACVec(q+nAnno);
  arma::vec MAFVec(q+nAnno);
  arma::vec missingRateVec(q+nAnno);     // missing rate
  
  std::vector<double> altBetaVec(q+nAnno);         // beta value for ALT allele
  std::vector<double> seBetaVec(q+nAnno);          // seBeta value
  std::vector<double> pval0Vec(q+nAnno);           // p values from normal distribution approximation  // might be confused, is this needed?
  std::vector<double> pval1Vec(q+nAnno);           // p values from more accurate methods including SPA and ER (ER is not used any more)
  
  std::vector<double> StatVec(q+nAnno);            // score statistics
  
  // updated on 2022-05-09: related to label
  arma::mat MACLabelMat(q+nAnno, t_nLabel);
  arma::mat MAFLabelMat(q+nAnno, t_nLabel);
  
  unsigned int m1 = g_region_maxMarkers_cutoff;     // number of markers in each chunk, the only exception is the last chunk.
  
  arma::mat P1Mat(m1, n);
  arma::mat P2Mat(n, m1); 
  
  std::vector<unsigned int> mPassCVVec;
  
  // conduct marker-level analysis
  double Stat, Beta, seBeta, pval0, pval1;
  arma::vec P1Vec(n), P2Vec(n);
  
  // initiate chunk information
  unsigned int nchunks = 0;    // total number of chunks
  unsigned int ichunk = 0;     // index of chunk
  unsigned int i1InChunk = 0;  // index of URV markers in chunk
  unsigned int i1 = 0;         // index of non-URV markers ()
  unsigned int i2 = 0;         // index of URV markers (Ultra-Rare Variants, URV)
  
  // arma::vec GVecURV(n, arma::fill::zeros);
  arma::mat GMatURV(n, nAnno, arma::fill::zeros);
  
  // updated on 2022-06-24 (save sum of genotype to conduct burden test and adjust p-values using SPA)
  unsigned int n_max_maf = g_region_max_maf_vec.size();
  arma::mat GMatBurden(n, nAnno * n_max_maf, arma::fill::zeros);
  arma::mat pvalBurden(nAnno * n_max_maf, 2, arma::fill::zeros);   // (only for GMatBurden) column 1: normal distribution; column 2: SPA.
  
  // updated on 2023-02-06
  arma::mat GMatBurdenNoWeight(n, nAnno * n_max_maf, arma::fill::zeros);       // no weights
    // column 1: Anno; 
    // column 2: MAF;
    // column 3: Sum; 
    // column 4: Beta; 
    // column 5: se.beta; 
    // column 6: normal distribution p-value; 
    // column 7: SPA p-value
  arma::mat infoBurdenNoWeight(nAnno * n_max_maf, 7, arma::fill::zeros);   
  
  
  
  boost::math::beta_distribution<> beta_dist(g_region_weight_beta[0], g_region_weight_beta[1]);
  
  // cycle for q markers
  for(unsigned int i = 0; i < q; i++)
  {
    double weight = t_weightVec.at(i);
    
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing, indexForNonZero;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    
    uint64_t gIndex = t_genoIndex.at(i);
    
    arma::vec test11 = getTime();
    
    arma::vec GVec = Unified_getOneMarker(t_genoType, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
    
    arma::vec test12 = getTime();
    g_compTime1 += test12 - test11;
    
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;
    
    flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, g_impute_method);
    
    double MAF = std::min(altFreq, 1 - altFreq);
    double MAC = MAF * 2 * n * (1 - missingRate);   // checked on 08-10-2021
    
    markerVec.at(i) = marker;             // marker IDs
    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT
    altFreqVec.at(i) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    MACVec.at(i) = MAC;
    MAFVec.at(i) = MAF;
    
    // Quality Control (QC)
    if((missingRate > g_missingRate_cutoff) || (MAF > g_region_maxMAF_cutoff) || MAF == 0){
      continue;  // does not pass QC
    }
    
    if(MAC > g_region_minMAC_cutoff){  // not Ultra-Rare Variants
      
      indicatorVec.at(i) = 1;
      
      if(i1InChunk == 0){
        std::cout << "Start analyzing chunk " << ichunk << "....." << std::endl;
      }
      
      arma::vec test21 = getTime();
      Unified_getRegionPVec(t_method, GVec, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
      
      arma::vec test22 = getTime();
      g_compTime2 += test22 - test21;
      
      // insert results to pre-setup vectors and matrix
      StatVec.at(i) = Stat;        
      altBetaVec.at(i) = Beta * (1 - 2*flip);  // Beta if flip = false, -1 * Beta is flip = true (hence, beta is only for alt allele)       
      seBetaVec.at(i) = seBeta;       
      pval0Vec.at(i) = pval0;
      pval1Vec.at(i) = pval1;
      
      P1Mat.row(i1InChunk) = P1Vec.t();
      P2Mat.col(i1InChunk) = P2Vec;
      
      if(t_nLabel != 1)
        getLabelInfo(GVec, t_labelVec, i, MACLabelMat, MAFLabelMat);
      
      i1 += 1;
      i1InChunk += 1;
      
      // updated on 2022-06-24 (save sum of genotype to conduct burden test and adjust p-values using SPA)
      double w0 = boost::math::pdf(beta_dist, MAF);
      for(unsigned int j = 0; j < n; j++)
      {
        if(GVec.at(j) != 0)
        {
          for(unsigned int iAnno = 0; iAnno < nAnno; iAnno++)
          {
            if(t_annoMat(i, iAnno) == 1)  // 1 or 0
            {
              for(unsigned int i_max_maf = 0; i_max_maf < n_max_maf; i_max_maf++)
              {
                double max_maf = g_region_max_maf_vec.at(i_max_maf);
                if(MAF < max_maf)
                {
                  GMatBurden(j, iAnno*n_max_maf+i_max_maf) += w0 * GVec.at(j);    // anno0,maf0; anno0,maf1; anno0,maf2; anno1,maf0; anno1,maf1, anno1,maf2; ...
                  // GMatBurden(j, iAnno) += w0 * weight * GVec.at(j); check it later (2022-06-24)
                  GMatBurdenNoWeight(j, iAnno*n_max_maf+i_max_maf) += GVec.at(j);
                }
              }
            }
          }
        }
      }
    }else{  // Ultra-Rare Variants
      
      indicatorVec.at(i) = 2;
      
      // double weight = t_weightVec.at(i);
      for(unsigned int j = 0; j < n; j++)
      {
        if(GVec.at(j) != 0)
        {
          // GVecURV.at(j) = std::max(GVecURV.at(j), weight*GVec.at(j)); 
          for(unsigned iAnno = 0; iAnno < nAnno; iAnno++)
          {
            if(t_annoMat(i, iAnno) == 1)  // 1 or 0
            {
              GMatURV(j, iAnno) = std::max(GMatURV(j, iAnno), weight * GVec.at(j)); 
              for(unsigned int i_max_maf = 0; i_max_maf < n_max_maf; i_max_maf++)
              {
                GMatBurdenNoWeight(j, iAnno*n_max_maf+i_max_maf) += GVec.at(j);
              }
            }
          }
        }
      }
      i2 += 1;
    }
    
    if(i1InChunk == m1)
    {
      std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      
      mPassCVVec.push_back(m1);
      ichunk += 1;
      i1InChunk = 0;
    }

    Rcpp::checkUserInterrupt();
  }
  
  printTimeDiff(g_compTime1, "Unified_getOneMarker");
  printTimeDiff(g_compTime2, "Unified_getRegionPVec");
  
  for(unsigned int iAnno = 0; iAnno < nAnno; iAnno++)
  {
    arma::vec GVecURV = GMatURV.col(iAnno);
    
    // updated on 2022-06-24 (save sum of genotype to conduct burden test and adjust p-values using SPA)
    double MAFURV = mean(GVecURV) / 2;
    double w0URV = boost::math::pdf(beta_dist, MAFURV);
    for(unsigned int i_max_maf = 0; i_max_maf < n_max_maf; i_max_maf++)
    {
      unsigned int i_pos = iAnno*n_max_maf+i_max_maf;
      GMatBurden.col(i_pos) += w0URV * GVecURV; 
      Unified_getRegionPVec(t_method, GMatBurden.col(i_pos), Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec); // use Unified_getMarkerPval to replace it later
      pvalBurden.at(i_pos, 0) = pval0;
      pvalBurden.at(i_pos, 1) = pval1;
      // std::cout << "iAnno:\t" << iAnno << std::endl;
      // std::cout << "i_max_maf:\t" << i_max_maf << std::endl;
      // std::cout << "pval0:\t" << pval0 << std::endl;
      // std::cout << "pval1:\t" << pval1 << std::endl;
      
      Unified_getRegionPVec(t_method, GMatBurdenNoWeight.col(i_pos), Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
      infoBurdenNoWeight.at(i_pos, 0) = iAnno;
      infoBurdenNoWeight.at(i_pos, 1) = i_max_maf;
      infoBurdenNoWeight.at(i_pos, 2) = sum(GMatBurdenNoWeight.col(i_pos));
      infoBurdenNoWeight.at(i_pos, 3) = Stat;
      infoBurdenNoWeight.at(i_pos, 4) = Beta;
      infoBurdenNoWeight.at(i_pos, 5) = seBeta;
      infoBurdenNoWeight.at(i_pos, 6) = pval1;
    }
    
    indicatorVec.at(q+iAnno) = 3;
    markerVec.at(q+iAnno) = t_annoVec.at(iAnno);
    infoVec.at(q+iAnno) = "Ultra-Rare Variants";
    altFreqVec.at(q+iAnno) = MAFVec.at(q+iAnno) = MAFURV;
    MACVec.at(q+iAnno) = sum(GVecURV);
    missingRateVec.at(q+iAnno) = 0;
    
    if(t_nLabel != 1)
      getLabelInfo(GVecURV, t_labelVec, q+iAnno, MACLabelMat, MAFLabelMat);
    
    // 2022-05-01: check if GVecURV are zero-vector later
    
    Unified_getRegionPVec(t_method, GVecURV, Stat, Beta, seBeta, pval0, pval1, P1Vec, P2Vec);
    
    StatVec.at(q+iAnno) = Stat;        
    altBetaVec.at(q+iAnno) = Beta;  // Beta if flip = false, -1 * Beta is flip = true (hence, beta is only for alt allele)       
    seBetaVec.at(q+iAnno) = seBeta;       
    pval0Vec.at(q+iAnno) = pval0;
    pval1Vec.at(q+iAnno) = pval1;
    
    if(i1InChunk >= m1)
    {
      P1Mat.resize(i1InChunk+1, n);
      P2Mat.resize(n, i1InChunk+1);
    }
    
    P1Mat.row(i1InChunk) = P1Vec.t();
    P2Mat.col(i1InChunk) = P2Vec;
    i1 += 1;
    i1InChunk += 1;
  }
  
  if(i2 == 0)
    std::cout << "i2 == 0." << std::endl;
   
  if((i1 == 0) & (i2 == 0))
    std::cout << "Cannot find any valid rare variants. This region will be skipped." << std::endl;

  mPassCVVec.push_back(i1InChunk);
  nchunks = ichunk + 1;
  
  if(i1InChunk != 0){
    P1Mat = P1Mat.rows(0, i1InChunk - 1);
    P2Mat = P2Mat.cols(0, i1InChunk - 1);
    if(nchunks != 1){
      std::cout << "In chunks 0-" << ichunk << ", " << i2 << " markers are ultra-rare and " << i1 << " markers are not ultra-rare." << std::endl;
      P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
      P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");
    }
  }
  
  arma::mat VarMat(i1, i1);
  
  // the number of non-ultra-rare markers in the region is less than limitation, so all data is in memory
  if(nchunks == 1){
    VarMat = P1Mat * P2Mat;
  }
  
  // the number of non-ultra-rare markers in the region is greater than limitation, so P1Mat and P2Mat are put in hard drive
  if(nchunks > 1)
  {
    int first_row = 0, first_col = 0, last_row = 0, last_col = 0;
    
    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      last_row = first_row + mPassCVVec.at(index1) - 1;
      
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      
      P1Mat.load(P1MatFile);
      
      if(P1Mat.n_cols == 0) continue;
      
      // off-diagonal sub-matrix
      for(unsigned int index2 = 0; index2 < index1; index2++)
      {
        std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index2 << "/" << nchunks - 1 << ")........" << std::endl;
        
        P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index2) + ".bin");
        
        if(P2Mat.n_cols == 0) continue;
        
        arma::mat offVarMat = P1Mat * P2Mat;
        
        // last_col = first_col + mPassQCVec.at(index2) - 1;
        last_col = first_col + mPassCVVec.at(index2) - 1;
        
        VarMat.submat(first_row, first_col, last_row, last_col) = offVarMat;
        
        // NOTE on 2022-04-28: here we assume that VarMat is symmetric, which is slightly different from (P1Mat * P2Mat)
        VarMat.submat(first_col, first_row, last_col, last_row) = offVarMat.t();
        
        first_col = last_col + 1;
      }
      
      // diagonal sub-matrix
      last_col = first_col + mPassCVVec.at(index1) - 1;
      std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index1 << "/" << nchunks - 1 << ")........" << std::endl;
      P2Mat.load(t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin");
      
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
      std::remove(File1);
      std::remove(File2);
    }
  }
  
  std::cout << "m1:\t" << m1 << std::endl;
  
  // Rcpp::Named("GVecURV") = GVecURV
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("MACVec") = MACVec,
                                          Rcpp::Named("MAFVec") = MAFVec,
                                          Rcpp::Named("MACLabelMat") = MACLabelMat,
                                          Rcpp::Named("MAFLabelMat") = MAFLabelMat,
                                          Rcpp::Named("StatVec") = StatVec,
                                          Rcpp::Named("altBetaVec") = altBetaVec,
                                          Rcpp::Named("seBetaVec") = seBetaVec,
                                          Rcpp::Named("pval0Vec") = pval0Vec,
                                          Rcpp::Named("pval1Vec") = pval1Vec,
                                          Rcpp::Named("indicatorVec") = indicatorVec,
                                          Rcpp::Named("VarMat") = VarMat,
                                          Rcpp::Named("pvalBurden") = pvalBurden,
                                          Rcpp::Named("infoBurdenNoWeight") = infoBurdenNoWeight);
  
  return OutList;
}

// [[Rcpp::export]]
void printTimeDiffInCPP()
{
  printTimeDiff(ptr_gPOLMMobj->getTestTime1(), "getRegionPVec");
  printTimeDiff(ptr_gPOLMMobj->getTestTime2(), "getSigmaxMat");
  printTimeDiff(ptr_gPOLMMobj->getTestTime3(), "solverBlockDiagSigma");
  printTimeDiff(ptr_gPOLMMobj->getTestTime4(), "get_ZPZ_adjGVec");
  printTimeDiff(ptr_gPOLMMobj->getTestTime5(), "getadjGFast Step 1");
  printTimeDiff(ptr_gPOLMMobj->getTestTime6(), "getadjGFast Step 2");
  printTimeDiff(ptr_gPOLMMobj->getTestTime7(), "get_ZPZ_adjGVec Step 1");
  printTimeDiff(ptr_gPOLMMobj->getTestTime8(), "get_ZPZ_adjGVec Step 2");
}

// [[Rcpp::export]]
void printTimeDiffSPAmixInCPP()
{
  printTimeDiff(ptr_gSPAmixobj->getTestTime1(), "SPAmix_SPA");
  printTimeDiff(ptr_gSPAmixobj->getTestTime2(), "SPAmix_MAF");
}


//////// ---------- Main function for genotype extraction --------- ////////////

// (2023-05-03): output two columns: column #1 is allele frequency and column #2 is missing rate

// [[Rcpp::export]]
arma::mat getGenoInfoInCPP(std::string t_genoType,
                           Rcpp::DataFrame t_markerInfo,
                           std::string t_imputeMethod) // 0: "mean"; 1: "minor"; 2: "drop" (to be continued)
{
  int q = t_markerInfo.nrow();         // number of markers requested
  arma::mat genoInfoMat(q, 2);
  
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];
  
  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;
  
  for(int i = 0; i < q; i++){
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType,          // "PLINK", "BGEN"
                                          gIndex,              // different meanings for different genoType
                                          ref,                 // REF allele
                                          alt,                 // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                          marker,              // marker ID extracted from genotype file
                                          pd,                  // base position
                                          chr,                 // chromosome
                                          altFreq,             // frequency of ALT allele
                                          altCounts,           // counts of ALT allele
                                          missingRate,         // missing rate
                                          imputeInfo,          // imputation information score, i.e., R2 (all 1 for PLINK)
                                          true,                // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,               // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    
    if((i+1) % 1000 == 0)
      std::cout << "Completed " << (i+1) << "/" << q << " genetic variants." << std::endl;
    
    // The below is not needed for information extraction
    // imputeGeno(GVec, altFreq, indexForMissing, t_imputeMethod);  // check UTIL.cpp
    genoInfoMat.at(i, 0) = altFreq;
    genoInfoMat.at(i, 1) = missingRate;
  }
  
  return genoInfoMat;
}
                         
                                          
// [[Rcpp::export]]
arma::mat getGenoInCPP(std::string t_genoType,
                       Rcpp::DataFrame t_markerInfo,
                       int n,
                       std::string t_imputeMethod) // 0: "mean"; 1: "minor"; 2: "drop" (to be continued)
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
    arma::vec GVec = Unified_getOneMarker(t_genoType,          // "PLINK", "BGEN"
                                          gIndex,              // different meanings for different genoType
                                          ref,                 // REF allele
                                          alt,                 // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                          marker,              // marker ID extracted from genotype file
                                          pd,                  // base position
                                          chr,                 // chromosome
                                          altFreq,             // frequency of ALT allele
                                          altCounts,           // counts of ALT allele
                                          missingRate,         // missing rate
                                          imputeInfo,          // imputation information score, i.e., R2 (all 1 for PLINK)
                                          true,                // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,               // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    
    imputeGeno(GVec, altFreq, indexForMissing, t_imputeMethod);  // check UTIL.cpp
    GMat.col(i) = GVec;
  }
  
  return GMat;
}

// This function will replace the above function (2022-01-28)

// [[Rcpp::export]]
arma::mat getGenoInCPP_fixedNumber(std::string t_genoType,
                                   Rcpp::DataFrame t_markerInfo,
                                   int n,
                                   std::string t_imputeMethod, // 0: "mean"; 1: "minor"; 2: "drop" (to be continued)
                                   int m,  // number of selected markers
                                   double missingRateCutoff,
                                   double minMAFCutoff)
{
  int q = t_markerInfo.nrow();         // number of markers requested
  std::vector<uint64_t> gIndexVec = t_markerInfo["genoIndex"];
  
  arma::mat GMat(n, m);
  int index = 0;
  
  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;
  
  for(int i = 0; i < q; i++){
    uint64_t gIndex = gIndexVec.at(i);
    std::cout << "gIndex:\t" << gIndex << std::endl;
    arma::vec GVec = Unified_getOneMarker(t_genoType,          // "PLINK", "BGEN"
                                          gIndex,              // different meanings for different genoType
                                          ref,                 // REF allele
                                          alt,                 // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                          marker,              // marker ID extracted from genotype file
                                          pd,                  // base position
                                          chr,                 // chromosome
                                          altFreq,             // frequency of ALT allele
                                          altCounts,           // counts of ALT allele
                                          missingRate,         // missing rate
                                          imputeInfo,          // imputation information score, i.e., R2 (all 1 for PLINK)
                                          true,                // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,               // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    
    if((altFreq < minMAFCutoff) | (altFreq > 1-minMAFCutoff) | (missingRate > missingRateCutoff))
      continue;
    
    imputeGeno(GVec, altFreq, indexForMissing, t_imputeMethod);  // check UTIL.cpp
    GMat.col(index) = GVec;
    index ++;
    if(index == m) break;
  }
  
  if(index < m)
    Rcpp::stop("No enough variants are for variance ratio estimation.");
  
  return GMat;
}

// [[Rcpp::export]]
arma::sp_mat getSpGenoInCPP(std::string t_genoType,
                            Rcpp::DataFrame t_markerInfo,
                            int n,
                            std::string t_imputeMethod)
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
                                          true,         // if true, output index of missing genotype data
                                          indexForMissing,     // index of missing genotype data
                                          false,               // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                          indexForNonZero);    // the index of non-zero genotype in the all subjects. Only valid if t_isOnlyOutputNonZero == true.
    
    imputeGeno(GVec, altFreq, indexForMissing, t_imputeMethod);  // check UTIL.hpp
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
                               bool t_isOnlyOutputNonZero,                   // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
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
void Unified_getMarkerPval(std::string t_method,   // "POLMM", "SPACox", "SAIGE", "SPAmix", and "SPAGRM"
                           arma::vec t_GVec,
                           bool t_isOnlyOutputNonZero,
                           std::vector<uint32_t> t_indexForNonZero,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
                           double& t_zScore,
                           double t_altFreq)
{
  // (BWJ) updated on 2023-04-20: I forgot what "t_isOnlyOutputNonZero" means, please let me know if anyone knows it,
  // (BWJ) it seems the "t_isOnlyOutputNonZero" can further save storage by removing the genotype == 0,
  // (BWJ) which is especially useful for region-based testing.
  if(t_method == "POLMM"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using POLMM method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
    
    ptr_gPOLMMobj->getMarkerPval(t_GVec, t_Beta, t_seBeta, t_pval, t_altFreq, t_zScore);
  }
  
  if(t_method == "SPACox"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SPACox method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
    t_pval = ptr_gSPACoxobj->getMarkerPval(t_GVec, t_altFreq, t_zScore);
  }
  
  if(t_method == "SPAmix"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SPAmix method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
    t_pval = ptr_gSPAmixobj->getMarkerPval(t_GVec, t_altFreq);
  }
}

// a unified function to get marker-level p-value
void Unified_getMarkerPval(std::string t_method,   // "POLMM", "SPACox", "SAIGE", "SPAmix", and "SPAGRM"
                           arma::vec t_GVec,
                           bool t_isOnlyOutputNonZero,
                           std::vector<uint32_t> t_indexForNonZero,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
                           double& t_zScore,
                           double t_altFreq,
                           double& t_hwepval,
                           double t_hwepvalCutoff)
{
  if(t_method == "SPAGRM"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SPAGRM method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
    t_pval = ptr_gSPAGRMobj->getMarkerPval(t_GVec, t_altFreq, t_zScore, t_hwepval, t_hwepvalCutoff);
  }else if(t_method == "SAGELD"){
    if(t_isOnlyOutputNonZero == true)
      Rcpp::stop("When using SAGELD method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
    t_pval = ptr_gSAGELDobj->getMarkerPval(t_GVec, t_altFreq, t_hwepval, t_hwepvalCutoff);
  }else{
    Unified_getMarkerPval(t_method, t_GVec,
                          false, // bool t_isOnlyOutputNonZero,
                          t_indexForNonZero, t_Beta, t_seBeta, t_pval, t_zScore, t_altFreq);
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
  // Check src/POLMM.cpp and src/POLMM.hpp
  if(t_method == "POLMM")
  {
    ptr_gPOLMMobj->getRegionPVec(t_GVec, t_Stat, t_Beta, t_seBeta, t_pval0, t_pval1, t_P1Vec, t_P2Vec);
  }
  
  // Check src/SPACox.cpp and src/SPACox.hpp
  if(t_method == "SPACox")
  {
    ptr_gSPACoxobj->getRegionPVec(t_GVec, t_Stat, t_pval0, t_pval1, t_P1Vec, t_P2Vec);
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
  if(ptr_gPLINKobj)
    delete ptr_gPLINKobj;
  
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
  if(ptr_gBGENobj){
    std::cout << "Deleting `ptr_gBGENobj`...." << std::endl;
    delete ptr_gBGENobj;
    ptr_gBGENobj = NULL;
  }
  
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

// [[Rcpp::export]]
void closeGenoInputInCPP(std::string t_genoType)  // "PLINK" or "BGEN"
{
  if(t_genoType == "PLINK")
  {
    delete ptr_gPLINKobj;
    ptr_gPLINKobj = NULL;
  }
  if(t_genoType == "BGEN")
  {
    delete ptr_gBGENobj;
    ptr_gBGENobj = NULL;
  }
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
  if(ptr_gPOLMMobj)
    delete ptr_gPOLMMobj;
  
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
                                 Rcpp::List t_controlList,
                                 arma::mat GenoMat)
{
  // arma::umat locations = t_SPmatR["locations"];
  // arma::vec values = t_SPmatR["values"];
  // std::cout << "Setting Sparse GRM...." << std::endl;
  // arma::sp_mat SparseGRM = arma::sp_mat(locations, values);
  
  // std::cout << "test, t_flagSparseGRM" << t_flagSparseGRM << std::endl;
  
  // The following function is in POLMM.cpp
  if(ptr_gPOLMMobj)
    delete ptr_gPOLMMobj;
  
  ptr_gPOLMMobj = new POLMM::POLMMClass(t_flagSparseGRM,       // if 1, then use Sparse GRM, otherwise, use Dense GRM
                                        ptr_gDenseGRMobj,
                                        ptr_gPLINKobj,
                                        ptr_gBGENobj,
                                        t_Cova,
                                        t_yVec,     // should be from 0 to J-1
                                        t_beta,
                                        t_bVec,
                                        t_eps,           // 
                                        t_tau,
                                        g_SparseGRM,    // results of function setSparseGRMInCPP()
                                        t_controlList);
  
  ptr_gPOLMMobj->fitPOLMM();
  
  // bool LOCO = t_controlList["LOCO"];
  
  // if(LOCO){
  //   
  //   // turn on LOCO option
  //   Rcpp::StringVector chrVec = m_ptrPlinkObj->getChrVec();
  //   Rcpp::StringVector uniqchr = unique(chrVec);
  //   
  //   std::cout << "uniqchr is " << uniqchr << std::endl;
  //   
  //   for(int i = 0; i < uniqchr.size(); i ++){
  //     
  //     std::string excludechr = std::string(uniqchr(i));
  //     std::cout << std::endl << "Leave One Chromosome Out: Chr " << excludechr << std::endl;
  //     
  //     updateParaConv(excludechr);
  //     
  //     arma::mat GMatRatio = m_ptrPlinkObj->getGMat(100, excludechr, m_minMafVarRatio, m_maxMissingVarRatio);
  //     arma::mat VarRatioMat = getVarRatio(GMatRatio, excludechr);
  //     double VarRatio = arma::mean(VarRatioMat.col(4));
  //     
  //     Rcpp::List temp = Rcpp::List::create(Rcpp::Named("muMat") = m_muMat,
  //                                          Rcpp::Named("iRMat") = m_iRMat,
  //                                          Rcpp::Named("VarRatioMat") = VarRatioMat,
  //                                          Rcpp::Named("VarRatio") = VarRatio);
  //     
  //     m_LOCOList[excludechr] = temp;
  //   }
  //   
  // }else{
  //   
  //   // turn off LOCO option
  //   arma::mat GMatRatio = m_ptrPlinkObj->getGMat(100, "none", m_minMafVarRatio, m_maxMissingVarRatio);
  //   arma::mat VarRatioMat = getVarRatio(GMatRatio, "none");
  //   double VarRatio = arma::mean(VarRatioMat.col(4));
  //   
  //   Rcpp::List temp = Rcpp::List::create(Rcpp::Named("muMat") = m_muMat,
  //                                        Rcpp::Named("iRMat") = m_iRMat,
  //                                        Rcpp::Named("VarRatioMat") = VarRatioMat,
  //                                        Rcpp::Named("VarRatio") = VarRatio);
  //   m_LOCOList["LOCO=F"] = temp;
  // }
  
  ptr_gPOLMMobj->estVarRatio(GenoMat);
  
  Rcpp::List outList = ptr_gPOLMMobj->getPOLMM();
  
  // ptr_gDenseGRMobj->closeDenseGRMObj();
  return outList;
}

// [[Rcpp::export]]
void setSPAGRMobjInCPP(arma::vec t_resid,
                       arma::vec t_resid_unrelated_outliers,
                       double t_sum_R_nonOutlier,
                       double t_R_GRM_R_nonOutlier,
                       double t_R_GRM_R_TwoSubjOutlier,
                       double t_R_GRM_R,
                       arma::vec t_MAF_interval,
                       Rcpp::List t_TwoSubj_list,
                       Rcpp::List t_ThreeSubj_list,
                       double t_SPA_Cutoff,
                       double t_zeta,
                       double t_tol)
{
  if(ptr_gSPAGRMobj)
    delete ptr_gSPAGRMobj;
  
  ptr_gSPAGRMobj = new SPAGRM::SPAGRMClass(t_resid,
                                           t_resid_unrelated_outliers,
                                           t_sum_R_nonOutlier,
                                           t_R_GRM_R_nonOutlier,
                                           t_R_GRM_R_TwoSubjOutlier,
                                           t_R_GRM_R,
                                           t_MAF_interval,
                                           t_TwoSubj_list,
                                           t_ThreeSubj_list,
                                           t_SPA_Cutoff,
                                           t_zeta,
                                           t_tol);
}

// [[Rcpp::export]]
void setSAGELDobjInCPP(std::string t_Method,
                       arma::mat t_XTs,
                       arma::mat t_SS,
                       arma::mat t_AtS,
                       arma::mat t_Q,
                       arma::mat t_A21,
                       arma::mat t_TTs,
                       arma::mat t_Tys,
                       arma::vec t_sol,
                       arma::vec t_blups,
                       double t_sig,
                       arma::vec t_resid,
                       arma::vec t_resid_G,
                       arma::vec t_resid_GxE,
                       arma::vec t_resid_E,
                       arma::vec t_resid_unrelated_outliers,
                       arma::vec t_resid_unrelated_outliers_G,
                       arma::vec t_resid_unrelated_outliers_GxE,
                       double t_sum_R_nonOutlier,
                       double t_sum_R_nonOutlier_G,
                       double t_sum_R_nonOutlier_GxE,
                       double t_R_GRM_R,
                       double t_R_GRM_R_G,
                       double t_R_GRM_R_GxE,
                       double t_R_GRM_R_G_GxE,
                       double t_R_GRM_R_E,
                       double t_R_GRM_R_nonOutlier,
                       double t_R_GRM_R_nonOutlier_G,
                       double t_R_GRM_R_nonOutlier_GxE,
                       double t_R_GRM_R_nonOutlier_G_GxE,
                       double t_R_GRM_R_TwoSubjOutlier,
                       double t_R_GRM_R_TwoSubjOutlier_G,
                       double t_R_GRM_R_TwoSubjOutlier_GxE,
                       double t_R_GRM_R_TwoSubjOutlier_G_GxE,
                       Rcpp::List t_TwoSubj_list,
                       Rcpp::List t_ThreeSubj_list,
                       arma::vec t_MAF_interval,
                       double t_zScoreE_cutoff,
                       double t_SPA_Cutoff,
                       double t_zeta,
                       double t_tol)
{
  if(ptr_gSAGELDobj)
    delete ptr_gSAGELDobj;
  
  ptr_gSAGELDobj = new SAGELD::SAGELDClass(t_Method, t_XTs, t_SS, t_AtS, t_Q, t_A21, t_TTs, t_Tys, t_sol, t_blups, t_sig, 
                                           t_resid, t_resid_G, t_resid_GxE, t_resid_E, t_resid_unrelated_outliers, t_resid_unrelated_outliers_G, t_resid_unrelated_outliers_GxE, 
                                           t_sum_R_nonOutlier, t_sum_R_nonOutlier_G, t_sum_R_nonOutlier_GxE, t_R_GRM_R, t_R_GRM_R_G, t_R_GRM_R_GxE, t_R_GRM_R_G_GxE, t_R_GRM_R_E, 
                                           t_R_GRM_R_nonOutlier, t_R_GRM_R_nonOutlier_G, t_R_GRM_R_nonOutlier_GxE, t_R_GRM_R_nonOutlier_G_GxE, t_R_GRM_R_TwoSubjOutlier, 
                                           t_R_GRM_R_TwoSubjOutlier_G, t_R_GRM_R_TwoSubjOutlier_GxE, t_R_GRM_R_TwoSubjOutlier_G_GxE, t_TwoSubj_list, t_ThreeSubj_list, 
                                           t_MAF_interval, t_zScoreE_cutoff, t_SPA_Cutoff, t_zeta, t_tol);
}

// [[Rcpp::export]]
void setSPAmixobjInCPP(arma::mat t_resid,
                       arma::mat t_PCs,
                       int t_N,
                       double t_SPA_Cutoff,
                       Rcpp::List t_outlierList)
{
  if(ptr_gSPAmixobj)
    delete ptr_gSPAmixobj;
  
  ptr_gSPAmixobj = new SPAmix::SPAmixClass(t_resid,
                                           // t_XinvXX,
                                           // t_tX,
                                           t_PCs,
                                           t_N,
                                           t_SPA_Cutoff,
                                           t_outlierList);
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
  if(ptr_gSPACoxobj)
    delete ptr_gSPACoxobj;
  
  ptr_gSPACoxobj = new SPACox::SPACoxClass(t_cumul,
                                           t_mresid,
                                           t_XinvXX,
                                           t_tX,
                                           t_N,
                                           t_pVal_covaAdj_Cutoff,
                                           t_SPA_Cutoff);
}

// [[Rcpp::export]]
void setWtSPAGobjInCPP(arma::vec t_mresid,
                       int t_N,
                       double t_SPA_Cutoff,
                       Rcpp::List t_outlierList)
{
  if(ptr_gWtSPAGobj)
    delete ptr_gWtSPAGobj;
  
  ptr_gWtSPAGobj = new WtSPAG::WtSPAGClass(t_mresid,
                                           t_N,
                                           t_SPA_Cutoff,
                                           t_outlierList);
}

// [[Rcpp::export]]
void updateQCInCPP(arma::vec t_AF_ref,
                   arma::vec t_AN_ref,
                   arma::vec t_pvalue_bat,
                   double t_pvalue_bat_cutoff)
{
  ptr_gWtSPAGobj->set_AF_ref(t_AF_ref);
  ptr_gWtSPAGobj->set_AN_ref(t_AN_ref);
  ptr_gWtSPAGobj->set_pvalue_bat(t_pvalue_bat);
  ptr_gWtSPAGobj->set_pvalue_bat_cutoff(t_pvalue_bat_cutoff);
}


