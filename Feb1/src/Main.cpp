
// This file includes the main codes to connect C++ and R

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <boost/math/distributions/beta.hpp>

#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
#include <cstdio>         // std::remove

#include "Main.hpp"
#include "PLINK.hpp"
#include "BGEN.hpp"
#include "UTIL.hpp"
#include "SPAmixPlus.hpp"

// global objects for different genotype formats

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static BGEN::BgenClass* ptr_gBGENobj = NULL;

// global objects for analysis method
static SPAmixPlus::SPAmixPlusClass* ptr_gSPAmixPlusobj = NULL;

// global variables for analysis
std::string g_impute_method;      // "mean", "minor", or "drop"
double g_missingRate_cutoff;
unsigned int g_omp_num_threads;
double g_marker_minMAF_cutoff;
double g_marker_minMAC_cutoff;

// the below is only valid for Group output
arma::uvec g_group;
bool g_ifOutGroup;
unsigned int g_nGroup;


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

void updateGroupInfo(arma::vec t_GVec,
                     std::vector<uint32_t> t_indexForMissing,
                     arma::vec& nSamplesInGroupVec,
                     arma::vec& AltCountsInGroupVec,
                     arma::vec& AltFreqInGroupVec)
{
  unsigned int n1 = t_GVec.size();
  nSamplesInGroupVec.zeros();
  AltCountsInGroupVec.zeros();
  
  if(t_indexForMissing.size() == 0)
    t_indexForMissing.push_back(n1);
  
  unsigned int i1 = 0;
  for(unsigned int i = 0; i < n1; i++){
    if(i == t_indexForMissing.at(i1)){
      if(i1 < t_indexForMissing.size() - 1)
        i1 ++;
    }else{
      unsigned int grp = g_group.at(i);
      nSamplesInGroupVec.at(grp) += 1;
      AltCountsInGroupVec.at(grp) += t_GVec.at(i);
    }
  }
  
  AltFreqInGroupVec = AltCountsInGroupVec / nSamplesInGroupVec / 2;
}

//////// ---------- Main function for marker-level analysis --------- ////////////

// [[Rcpp::export]]
Rcpp::List mainMarkerInCPP(std::string t_method,       // "SPAmixPlus"
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<uint64_t> t_genoIndex,
                           std::vector<std::string> t_outputColumns = std::vector<std::string>(),
                           std::string t_afFile = "") // New arg
{
  int q = t_genoIndex.size();  // number of markers

  std::ifstream afStream;
  bool useAFFile = false;
  int nPCs = 0;
  size_t recordSize = 0;
  
  if(!t_afFile.empty() && t_method == "SPAmixPlus" && ptr_gSPAmixPlusobj){
      afStream.open(t_afFile, std::ios::binary);
      if(!afStream){
         Rcpp::stop("Cannot open AF file: " + t_afFile);
      }
      useAFFile = true;
      nPCs = ptr_gSPAmixPlusobj->getNPCs();
      recordSize = sizeof(int) + (nPCs + 1) * sizeof(double);
  }

  // set up output
  std::vector<std::string> markerVec(q);  // marker IDs
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> altFreqVec(q);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> altCountsVec(q);    // allele counts of ALT allele.
  std::vector<double> missingRateVec(q);  // missing rate of markers
  
  int Npheno = 1;
  if(t_method == "SPAmixPlus" && ptr_gSPAmixPlusobj)
    Npheno = ptr_gSPAmixPlusobj->getNpheno();
  
  std::vector<double> pvalVec(q*Npheno, arma::datum::nan);
  std::vector<double> zScoreVec(q*Npheno, arma::datum::nan);
  std::vector<double> SVec(q*Npheno, arma::datum::nan);
  std::vector<double> SmeanVec(q*Npheno, arma::datum::nan);
  std::vector<double> VarSVec(q*Npheno, arma::datum::nan);
  std::vector<double> BetaVec(q*Npheno, arma::datum::nan);         
  std::vector<double> seBetaVec(q*Npheno, arma::datum::nan);  
  
  arma::mat nSamplesInGroup;
  arma::mat AltCountsInGroup;
  arma::mat AltFreqInGroup;
  
  if(g_ifOutGroup){
    nSamplesInGroup.resize(q, g_nGroup);
    AltCountsInGroup.resize(q, g_nGroup);
    AltFreqInGroup.resize(q, g_nGroup);
  }
  
  // loop for all markers
  for(int i = 0; i < q; i++)
  {
    if(i % 1000 == 0)
      Rcpp::Rcout << "Completed " << i << "/" << q << " markers in the chunk." << std::endl;
    
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
    
    if(g_ifOutGroup){
      arma::vec nSamplesInGroupVec(g_nGroup);
      arma::vec AltCountsInGroupVec(g_nGroup);
      arma::vec AltFreqInGroupVec(g_nGroup);
      
      updateGroupInfo(GVec, indexForMissing, nSamplesInGroupVec, AltCountsInGroupVec, AltFreqInGroupVec);
      
      nSamplesInGroup.row(i) = nSamplesInGroupVec.t();
      AltCountsInGroup.row(i) = AltCountsInGroupVec.t();
      AltFreqInGroup.row(i) = AltFreqInGroupVec.t();
    }
    
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
    
    SPAmixPlus::SPAmixPlusClass::AFModelInfo model;
    SPAmixPlus::SPAmixPlusClass::AFModelInfo* pModel = NULL;

    if(useAFFile){
        // Read from file using random access
        afStream.seekg(gIndex * recordSize, std::ios::beg);
        if(!afStream.fail()){
            afStream.read(reinterpret_cast<char*>(&model.status), sizeof(int));
            model.betas.set_size(nPCs + 1);
            afStream.read(reinterpret_cast<char*>(model.betas.memptr()), (nPCs + 1) * sizeof(double));
            if(afStream){
                 pModel = &model;
            }
        }
    }
    
    // analysis results for single-marker
    double Beta, seBeta, pval, zScore;
    // double hwepval = 0;
    // double hwepvalCutoff = 0.1;

    Unified_getMarkerPval(t_method, GVec, false, indexForNonZero, Beta, seBeta, pval, zScore, altFreq, pModel);
    
    // Unified_getMarkerPval call handled above with model support
    
    if(t_method == "SPAmixPlus" && ptr_gSPAmixPlusobj){
      arma::vec pvalVecTemp = ptr_gSPAmixPlusobj->getpvalVec();
      arma::vec zScoreVecTemp = ptr_gSPAmixPlusobj->getzScoreVec();
      arma::vec BetaVecTemp = ptr_gSPAmixPlusobj->getBetaVec();
      arma::vec SVecTemp = ptr_gSPAmixPlusobj->getSVec();
      arma::vec SmeanVecTemp = ptr_gSPAmixPlusobj->getSmeanVec();
      arma::vec VarSVecTemp = ptr_gSPAmixPlusobj->getVarSVec();
      
      for(int j = 0; j < Npheno; j++){
        pvalVec.at(i*Npheno+j) = pvalVecTemp.at(j);
        zScoreVec.at(i*Npheno+j) = zScoreVecTemp.at(j);
        BetaVec.at(i*Npheno+j) = BetaVecTemp.at(j);
        SVec.at(i*Npheno+j) = SVecTemp.at(j);
        SmeanVec.at(i*Npheno+j) = SmeanVecTemp.at(j);
        VarSVec.at(i*Npheno+j) = VarSVecTemp.at(j);
      }
    }else{
      // Fallback if not specifically handled above or Npheno=1
       pvalVec.at(i) = pval;
       zScoreVec.at(i) = zScore; 
       BetaVec.at(i) = Beta * (1 - 2*flip);  
       seBetaVec.at(i) = seBeta;       
    }
  }
  
  Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("markerVec") = markerVec,
                                          Rcpp::Named("infoVec") = infoVec,
                                          Rcpp::Named("altFreqVec") = altFreqVec,
                                          Rcpp::Named("altCountsVec") = altCountsVec,
                                          Rcpp::Named("missingRateVec") = missingRateVec,
                                          Rcpp::Named("Stat") = SVec,
                                          Rcpp::Named("StatMean") = SmeanVec,
                                          Rcpp::Named("StatVar") = VarSVec,
                                          Rcpp::Named("zScore") = zScoreVec,
                                          Rcpp::Named("pvalVec") = pvalVec,
                                          Rcpp::Named("BetaG") = BetaVec,
                                          Rcpp::Named("seBeta") = seBetaVec,
                                          Rcpp::Named("nSamplesInGroup") = nSamplesInGroup,
                                          Rcpp::Named("AltCountsInGroup") = AltCountsInGroup,
                                          Rcpp::Named("AltFreqInGroup") = AltFreqInGroup);
  
  return OutList;  
}



// a unified function to get marker-level p-value
void Unified_getMarkerPval(std::string t_method, 
                           arma::vec t_GVec,
                           bool t_isOnlyOutputNonZero,
                           std::vector<uint32_t> t_indexForNonZero,
                           double& t_Beta, 
                           double& t_seBeta, 
                           double& t_pval,
                           double& t_zScore,
                           double t_altFreq,
                           SPAmixPlus::SPAmixPlusClass::AFModelInfo* t_afModel) // Optional model
{
  if(t_method == "SPAmixPlus"){
    if(t_isOnlyOutputNonZero == true)
       Rcpp::stop("When using SPAmixPlus method to calculate marker-level p-values, 't_isOnlyOutputNonZero' shold be false.");
     
     if(ptr_gSPAmixPlusobj){
        if(t_afModel){
            // Use pre-calculated model
            t_pval = ptr_gSPAmixPlusobj->getMarkerPvalFromModel(t_GVec, *t_afModel, t_altFreq);
        } else {
            // Use on-the-fly estimation
            t_pval = ptr_gSPAmixPlusobj->getMarkerPval(t_GVec, t_altFreq);
        }
     }
     else
        t_pval = arma::datum::nan;
  }
}

// [[Rcpp::export]]
void exportAFModelInCPP(std::string t_method,
                        std::string t_genoType,
                        std::vector<uint64_t> t_genoIndex,
                        std::string t_outputFile)
{
    if (t_method != "SPAmixPlus" || !ptr_gSPAmixPlusobj) {
        Rcpp::stop("exportAFModelInCPP only supports SPAmixPlus method and initialized object.");
    }
    
    // Create file if not exists, to check/ensure we can open in random access mode
    {
        std::ifstream testOpen(t_outputFile);
        if(!testOpen.good()){
            std::ofstream create(t_outputFile, std::ios::binary);
            create.close();
        }
    }

    // Open binary file in update mode
    std::fstream outFile(t_outputFile, std::ios::binary | std::ios::in | std::ios::out);
    if (!outFile) Rcpp::stop("Failed to open output file " + t_outputFile);
    
    int nPCs = ptr_gSPAmixPlusobj->getNPCs();
    long long recordSize = sizeof(int) + (long long)(nPCs + 1) * sizeof(double);
    // Write header? No, just raw records to allow seek.
    // Assuming t_genoIndex matches the file structure or we append.
    // If we append, we can't seek easily unless we know order.
    // We assume the user manages the file creation properly (e.g. per chunk output).
    
    int q = t_genoIndex.size();
    
    // Progress
    int nPercent = q / 100;
    if(nPercent == 0) nPercent = 1;
    
    for(int i = 0; i < q; i++){
        if(i % nPercent == 0){
             Rcpp::checkUserInterrupt();
             Rprintf("Processed %d / %d markers (%.1f%%) \r", i, q, (100.0 * i) / q); 
        }

        uint64_t genoIndex = t_genoIndex[i];
        
        std::string marker;
        std::string chr;
        uint32_t pd;
        std::string ref;
        std::string alt;
        double altFreq;
        double altCounts; // Added
        double missingRate;
        std::vector<uint32_t> indexForMissing;
        std::vector<uint32_t> indexForNonZero;
        arma::vec GVec;
        
        double imputeInfo; // QC metrics

        // Use Unified getter
        GVec = Unified_getOneMarker(t_genoType, genoIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          true, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          false, // bool t_isOnlyOutputNonZero,
                                          indexForNonZero);
        
        // Impute
        // Note: fit_lm/logistic etc handle imputation internally in `computeAFModel` for logistic?
        // Wait, `fit_lm` expects imputed?
        // Let's check `getMAFest`. It uses `t_GVec`. `getMAFest` calls `fit_lm(g)`.
        // `g` is `t_GVec`.
        // Prior to passing `t_GVec` to `getMAFest`, it is imputed in `Main.cpp` loop.
        imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, g_impute_method);
        // Note: If flipped, GVec is 2-G. `computeAFModel` works on G.
        // The model coefficients will be for the flipped G.
        // In Step 2, if we encounter the same marker, `imputeGenoAndFlip` will flip it again.
        // Then we apply model.
        // The model expects the specific G vector distribution.
        
        // Compute Model
        SPAmixPlus::SPAmixPlusClass::AFModelInfo model = ptr_gSPAmixPlusobj->computeAFModel(GVec, altFreq);
        
        // Write to file with random access
        long long pos = (long long)genoIndex * recordSize;
        outFile.seekp(pos, std::ios::beg);
        
        // Format: [Status: int] [Betas: double array of size K+1]
        outFile.write(reinterpret_cast<const char*>(&model.status), sizeof(int));
        outFile.write(reinterpret_cast<const char*>(model.betas.memptr()), (nPCs + 1) * sizeof(double));
    }
    outFile.close();
}


void Unified_getMarkerPval(std::string t_method,   
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
    // Delegate to simpler function
    Unified_getMarkerPval(t_method, t_GVec, t_isOnlyOutputNonZero, t_indexForNonZero, t_Beta, t_seBeta, t_pval, t_zScore, t_altFreq);
}


// [[Rcpp::export]]
arma::mat getGenoInCPP(std::string t_genoType,
                       Rcpp::DataFrame t_markerInfo,
                       int n,
                       std::string t_imputeMethod)
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
    arma::vec GVec = Unified_getOneMarker(t_genoType,          
                                          gIndex,              
                                          ref,                 
                                          alt,                
                                          marker,              
                                          pd,                  
                                          chr,                 
                                          altFreq,             
                                          altCounts,           
                                          missingRate,         
                                          imputeInfo,          
                                          true,                
                                          indexForMissing,     
                                          false,               
                                          indexForNonZero);    
    
    imputeGeno(GVec, altFreq, indexForMissing, t_imputeMethod);  
    GMat.col(i) = GVec;
  }
  
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
  arma::sp_mat GMat(n, q);             
  
  std::string ref, alt, marker, chr;
  uint32_t pd;
  double altFreq, altCounts, missingRate, imputeInfo;
  std::vector<uint32_t> indexForMissing, indexForNonZero;
  
  for(int i = 0; i < q; i++){
    uint64_t gIndex = gIndexVec.at(i);
    arma::vec GVec = Unified_getOneMarker(t_genoType,    
                                          gIndex,        
                                          ref,           
                                          alt,           
                                          marker,        
                                          pd,            
                                          chr,           
                                          altFreq,       
                                          altCounts,     
                                          missingRate,   
                                          imputeInfo,    
                                          true,         
                                          indexForMissing,     
                                          false,               
                                          indexForNonZero);    
    
    imputeGeno(GVec, altFreq, indexForMissing, t_imputeMethod);  
    GMat.col(i) = arma::sp_mat(GVec); 
  }
  
  return GMat;
}

// a unified function to get single marker from genotype file
arma::vec Unified_getOneMarker(std::string t_genoType,   // "PLINK", "BGEN"
                               uint64_t t_gIndex,        // different meanings for different genoType
                               std::string& t_ref,       // REF allele
                               std::string& t_alt,       // ALT allele 
                               std::string& t_marker,    // marker ID 
                               uint32_t& t_pd,           // base position
                               std::string& t_chr,       // chromosome
                               double& t_altFreq,        // frequency of ALT allele
                               double& t_altCounts,      // counts of ALT allele
                               double& t_missingRate,    // missing rate
                               double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                               bool t_isOutputIndexForMissing,               
                               std::vector<uint32_t>& t_indexForMissing,     
                               bool t_isOnlyOutputNonZero,                   
                               std::vector<uint32_t>& t_indexForNonZero)     
{
  arma::vec GVec;
  if(t_genoType == "PLINK" && ptr_gPLINKobj){
    GVec = ptr_gPLINKobj->getOneMarker(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                       t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                       true);   
  }
  
  if(t_genoType == "BGEN" && ptr_gBGENobj){
    bool isBoolRead;
    GVec = ptr_gBGENobj->getOneMarker(t_gIndex, t_ref, t_alt, t_marker, t_pd, t_chr, t_altFreq, t_altCounts, t_missingRate, t_imputeInfo, 
                                      t_isOutputIndexForMissing, t_indexForMissing, t_isOnlyOutputNonZero, t_indexForNonZero,
                                      isBoolRead);
  }
  
  return GVec;
}

// [[Rcpp::export]]
void setPLINKobjInCPP(std::string t_bimFile,
                      std::string t_famFile,
                      std::string t_bedFile,
                      std::vector<std::string> t_SampleInModel,
                      std::string t_AlleleOrder)
{
  if(ptr_gPLINKobj) delete ptr_gPLINKobj;
  ptr_gPLINKobj = new PLINK::PlinkClass(t_bimFile, t_famFile, t_bedFile, t_SampleInModel, t_AlleleOrder);
}

// [[Rcpp::export]]
void setBGENobjInCPP(std::string t_bgenFileName,
                     std::string t_bgiFileName,
                     std::vector<std::string> t_SampleInBgen,
                     std::vector<std::string> t_SampleInModel,
                     bool t_isSparseDosageInBgen,
                     bool t_isDropVarNotWork, 
                     std::string t_AlleleOrder) 
{
  if(ptr_gBGENobj) delete ptr_gBGENobj;
  ptr_gBGENobj = new BGEN::BgenClass(t_bgenFileName, t_bgiFileName, t_SampleInBgen, t_SampleInModel, t_isSparseDosageInBgen, t_isDropVarNotWork, t_AlleleOrder);
}

// [[Rcpp::export]]
void closeGenoInputInCPP(std::string t_genoType)
{
  if(t_genoType == "PLINK" && ptr_gPLINKobj){
    delete ptr_gPLINKobj;
    ptr_gPLINKobj = NULL;
  }

  if(t_genoType == "BGEN" && ptr_gBGENobj){
    delete ptr_gBGENobj; 
    ptr_gBGENobj = NULL;
  }
}

// [[Rcpp::export]]
void setSPAmixPlusobjInCPP(arma::mat t_resid,
                           arma::mat t_PCs,
                           int t_N,
                           double t_SPA_Cutoff,
                           Rcpp::List t_outlierList,
                           Rcpp::DataFrame t_sparseGRM,
                           Rcpp::DataFrame t_ResidMat
                           )
{
  if(ptr_gSPAmixPlusobj)
    delete ptr_gSPAmixPlusobj;

  ptr_gSPAmixPlusobj = new SPAmixPlus::SPAmixPlusClass(t_resid,
                                                      t_PCs,
                                                      t_N,
                                                      t_SPA_Cutoff,
                                                      t_outlierList,
                                                      t_sparseGRM,
                                                      t_ResidMat
                                                      );
}
