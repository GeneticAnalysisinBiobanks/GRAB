
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
double g_missingRate_cutoff;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;
std::string g_impute_method;

void setMarker_GlobalVarsInCPP(std::string t_impute_method,
                               double t_missing_cutoff,
                               double t_min_maf_marker,
                               double t_min_mac_marker)
{
  g_impute_method = t_impute_method;
  g_missingRate_cutoff = t_missing_cutoff;
  g_marker_minMAF_cutoff = t_min_maf_marker;
  g_marker_minMAC_cutoff = t_min_mac_marker;
}

// global varables for region-level analysis
// double g_marker_minMAF_cutoff;
// double g_marker_minMAC_cutoff;
// std::string g_impute_method;
// 
// void setMarker_GlobalVarsInCPP(std::string t_impute_method,
//                                double t_missing_cutoff,
//                                double t_min_maf_marker,
//                                double t_min_mac_marker)
// {
//   g_impute_method = t_impute_method;
//   g_missingRate_cutoff = t_missing_cutoff;
//   g_marker_minMAF_cutoff = t_min_maf_marker;
//   g_marker_minMAC_cutoff = t_min_mac_marker;
// }

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
      flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing);  // in UTIL.cpp
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
Rcpp::List MAIN_REGION(std::vector<std::string> t_MarkerReqstd,
                       double t_NonZero_cutoff,
                       double t_StdStat_cutoff,
                       int t_maxMarkers,
                       std::string t_outputFile,
                       double t_missingRate_cutoff,
                       double t_maxMAF_cutoff,
                       std::string t_kernel,
                       arma::vec t_wBeta)
{
  // extract information from global variable ptr_gPLINKobj
  int n = ptr_gPLINKobj->getN();
  std::vector<uint32_t> posMarkerInPlink = ptr_gPLINKobj->getPosMarkerInPlink(t_MarkerReqstd);
  int q = posMarkerInPlink.size();         // number of markers in the region
  
  // set up output
  std::vector<std::string> a1Vec; 
  std::vector<std::string> a2Vec; 
  std::vector<std::string> markerVec;
  std::vector<uint32_t> pdVec;
  std::vector<double> freqVec; 
  std::vector<double> weightVec; 
  std::vector<double> StatVec;
  std::vector<bool> flipVec;
  std::vector<double> pvalNormVec;
  std::vector<double> pvalVec;
  std::vector<int> posVec;
  arma::mat VarSMat;
  
  //
  arma::mat adjGMat(t_maxMarkers, n);       // adjusted genotype vector 
  arma::mat ZPZ_adjGMat(n, t_maxMarkers);   // t(Z) %*% P %*% Z %*% adjGMat
  ptr_gPOLMMobj->setSeqMat(t_NonZero_cutoff);
  
  int indexPassingQC = 0;
  int indexChunkSave = 0;
  
  arma::vec K1roots(2);
  K1roots(0) = 3;
  K1roots(1) = -3;
  
  arma::vec wGVecBT(n, arma::fill::zeros);
  
  // loop for all markers
  for(int i = 0; i < q; i++){
    
    if(i % 1000 == 0){
      std::cout << "Completed " << i << "/" << q << " markers in the region." << std::endl;
      std::cout << "indexPassingQC:\t" << indexPassingQC << std::endl;
    }
    
    uint32_t posMarker = posMarkerInPlink.at(i);
    double freq, missingRate;
    std::vector<uint32_t> posMissingGeno;
    std::string a1, a2, marker;
    uint32_t pd;
    uint8_t chr;
    bool flip = false;
    bool flagTrueGeno = true;
    
    arma::vec GVec;  // to be added later
    // arma::vec GVec = ptr_gPLINKobj->getOneMarker(posMarker, freq, missingRate, posMissingGeno,
    //                                              a1, a2, marker, pd, chr, flagTrueGeno);
    double MAF = std::min(freq, 1 - freq);
    
    // Quality Control (QC) based on missing rate and allele frequency
    if((missingRate > t_missingRate_cutoff) || (MAF > t_maxMAF_cutoff) || (MAF == 0))
      continue;
    
    // push back to output
    double weight = getWeights(t_kernel, MAF, t_wBeta);
    weightVec.push_back(weight);
    
    a1Vec.push_back(a1);
    a2Vec.push_back(a2);
    markerVec.push_back(marker);
    pdVec.push_back(pd);
    posVec.push_back(i);
    
    if(missingRate != 0){
      // imputeGeno(GVec, freq, posMissingGeno);
      // flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing);  // in UTIL.cpp
    }
      
    
    if(freq > 0.5){
      GVec = 2 - GVec;
      flip = true;
    }
    
    freqVec.push_back(MAF);
    flipVec.push_back(flip);
    wGVecBT += GVec * weight;
    
    arma::vec adjGVec = ptr_gPOLMMobj->getadjGFast(GVec);
    double Stat = ptr_gPOLMMobj->getStatFast(adjGVec);
    StatVec.push_back(Stat);
    
    // get t(Z) %*% P %*% Z %*% adjGVec for each marker
    arma::vec ZPZ_adjGVec = ptr_gPOLMMobj->get_ZPZ_adjGVec(adjGVec);
    double VarS = as_scalar(adjGVec.t() * ZPZ_adjGVec);
    double StdStat = std::abs(Stat) / sqrt(VarS);
    double pvalNorm = 2 * arma::normcdf(-1*StdStat);
    double pval = pvalNorm;
    
    if(StdStat > t_StdStat_cutoff){ // then use SPA/ER to correct p value
      // functions of SPA or ER
      arma::uvec posG1 = arma::find(GVec != 0);
      std::cout << "posG1.size():\t" << posG1.size() << std::endl;
      int nG1 = posG1.size();
      
      if(nG1 <= t_NonZero_cutoff){
        double pvalER = ptr_gPOLMMobj->MAIN_ER(GVec, posG1);
        pval = pvalER;
        // std::cout << "pvalER:\t" << pvalER << std::endl;
      }else{
        arma::vec VarWVec = ptr_gPOLMMobj->getVarWVec(adjGVec);
        double VarW = sum(VarWVec);
        double VarW1 = sum(VarWVec(posG1));
        double VarW0 = VarW - VarW1;
        double Ratio0 = VarW0 / VarW;
        
        Rcpp::List resSPA = ptr_gPOLMMobj->MAIN_SPA(Stat, adjGVec, K1roots, VarS, VarW, Ratio0, posG1);
        pval = resSPA["pval"];
      }
      // if(nG1 > t_NonZero_cutoff){
      
      
      // double pvalSPA = resSPA["pval"];
      // std::cout << "pvalSPA:\t" << pvalSPA << std::endl;
      
      
      
      // std::cout << resSPA << std::endl;
      // std::cout << "pvalNorm:\t" << 2 * arma::normcdf(-1*StdStat) << std::endl;
      
      // }else{
      // something to add for Efficient Resampling (ER)
      //  Rcpp::List resSPA = ptr_gPOLMMobj->MAIN_SPA(Stat, GVec, adjGVec, K1roots, VarP, VarW, Ratio0, posG1);
      // }
    }
    
    // insert adjGVec and ZPZ_adjGVec into a pre-defined matrix
    pvalNormVec.push_back(pvalNorm);
    pvalVec.push_back(pval);
    adjGMat.row(indexPassingQC) = adjGVec.t();
    ZPZ_adjGMat.col(indexPassingQC) = ZPZ_adjGVec;
    
    indexPassingQC++;
    
    if(indexPassingQC % t_maxMarkers == 0){
      adjGMat.save(t_outputFile + "_adjGMat" + std::to_string(indexChunkSave) + ".bin");
      ZPZ_adjGMat.save(t_outputFile + "_ZPZ_adjGMat" + std::to_string(indexChunkSave) + ".bin");
      indexPassingQC = 0;
      std::cout << "Completed chunk "<< indexChunkSave << "!" << std::endl;
      indexChunkSave++;
    }
    Rcpp::checkUserInterrupt();
  }
  
  // additional burden test to further adjust for the variance
  arma::vec wadjGVecBT = ptr_gPOLMMobj->getadjGFast(wGVecBT);
  double wStatBT = ptr_gPOLMMobj->getStatFast(wadjGVecBT);
  arma::vec ZPZ_wadjGVecBT = ptr_gPOLMMobj->get_ZPZ_adjGVec(wadjGVecBT);
  double wVarSBT = as_scalar(wadjGVecBT.t() * ZPZ_wadjGVecBT);
  double wStdStatBT = std::abs(wStatBT) / sqrt(wVarSBT);
  double rBT = 1;
  
  if(wStdStatBT > t_StdStat_cutoff){
    arma::uvec poswG1BT = arma::find(wGVecBT != 0);
    double wVarSBT = as_scalar(wadjGVecBT.t() * ZPZ_wadjGVecBT);
    arma::vec wVarWVecBT = ptr_gPOLMMobj->getVarWVec(wadjGVecBT);
    double wVarWBT = sum(wVarWVecBT);
    double wVarW1BT = sum(wVarWVecBT(poswG1BT));
    double wVarW0BT = wVarWBT - wVarW1BT;
    double wRatio0BT = wVarW0BT / wVarWBT;
    Rcpp::List resSPA = ptr_gPOLMMobj->MAIN_SPA(wStatBT, wadjGVecBT, K1roots, wVarSBT, wVarWBT, wRatio0BT, poswG1BT);
    Rcpp::NumericVector wadjPvalBT = {resSPA["pval"]};
    Rcpp::NumericVector temp = Rcpp::qchisq(wadjPvalBT, 1, false, false);
    double wadjVarSBT = pow(wStatBT,2) / temp(0);
    rBT = wadjVarSBT / wVarSBT;
    std::cout << "rBT:\t" << rBT << std::endl;
  }
  
  // total number of markers that pass QC from MAF and missing rate
  int nPassingQC = indexPassingQC + indexChunkSave * t_maxMarkers;
  
  if(indexPassingQC != 0){
    adjGMat = adjGMat.rows(0, indexPassingQC - 1);
    ZPZ_adjGMat = ZPZ_adjGMat.cols(0, indexPassingQC - 1);
  }
  
  // the region includes more markers than memory requested
  if((indexChunkSave > 0) & (indexPassingQC != 0)){ 
    adjGMat.save(t_outputFile + "_adjGMat" + std::to_string(indexChunkSave) + ".bin");
    ZPZ_adjGMat.save(t_outputFile + "_ZPZ_adjGMat" + std::to_string(indexChunkSave) + ".bin");
    std::cout << "Completed chunk "<< indexChunkSave << "!" << std::endl;
    indexChunkSave++;
  }
  
  // calculate variance-covariance matrix
  VarSMat.resize(nPassingQC, nPassingQC);    // variance matrix (after adjusting for relatedness)
  
  // not so many markers in the region, so everything is in memory
  if(indexChunkSave == 0)
  {
    VarSMat = adjGMat * ZPZ_adjGMat;
  }
  
  // the region includes more markers than limitation, so everything is in hard-drive storage
  if(indexChunkSave > 0)
  {
    int first_row = 0;
    int first_col = 0;
    int last_row = t_maxMarkers - 1;
    int last_col = t_maxMarkers - 1;
    
    for(int index1 = 0; index1 < indexChunkSave; index1++)
    {
      adjGMat.load(t_outputFile + "_adjGMat" + std::to_string(index1) + ".bin");
      
      // off-diagonal sub-matrix
      for(int index2 = 0; index2 < index1; index2++)
      {
        std::cout << "Analyzing chunks (" << index1 << "/" << indexChunkSave - 1 << ", " << index2 << "/" << indexChunkSave - 1 << ")........" << std::endl;
        ZPZ_adjGMat.load(t_outputFile + "_ZPZ_adjGMat" + std::to_string(index2) + ".bin");
        
        arma::mat offVarSMat = adjGMat * ZPZ_adjGMat;
        VarSMat.submat(first_row, first_col, 
                       std::min(last_row, nPassingQC - 1), 
                       std::min(last_col, nPassingQC - 1)) = offVarSMat;
        VarSMat.submat(first_col, first_row, 
                       std::min(last_col, nPassingQC - 1), 
                       std::min(last_row, nPassingQC - 1)) = offVarSMat.t();
        first_col += t_maxMarkers;
        last_col += t_maxMarkers;
      }
      
      // // diagonal sub-matrix
      std::cout << "Analyzing chunks (" << index1 << "/" << indexChunkSave - 1 << ", " << index1 << "/" << indexChunkSave - 1 << ")........" << std::endl;
      ZPZ_adjGMat.load(t_outputFile + "_ZPZ_adjGMat" + std::to_string(index1) + ".bin");
      // arma::fmat diagVarSMat = getSymmMat(adjGMat, ZPZ_adjGMat, adjGMat.n_cols);
      arma::mat diagVarSMat = adjGMat * ZPZ_adjGMat;
      VarSMat.submat(first_row, first_col, 
                     std::min(last_row, nPassingQC - 1), 
                     std::min(last_col, nPassingQC - 1)) = diagVarSMat;
      first_row += t_maxMarkers;
      last_row += t_maxMarkers;
      first_col = 0;
      last_col = t_maxMarkers - 1;
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
  return OutList;
}






