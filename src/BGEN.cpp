
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(BH)]]

#include <cstring>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <numeric>
#include <utility>
#include <stdexcept>
#include <time.h>
#include <stdint.h>
#include <zlib.h>
#include "BGEN.h"

namespace BGEN {

BgenClass::BgenClass(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> t_SampleInBgen,
                     std::vector<std::string> t_SampleInModel,
                     bool t_isSparseDosageInBgen,
                     bool t_isDropmissingdosagesInBgen,
                     std::string t_AlleleOrder)      // added by Wenjian Bi on 03/14/2021: "ref-first" or "alt-first"
{
  setBgenObj(t_bgenFileName,
             t_bgenFileIndex,
             t_SampleInBgen);
  
  setPosSampleInBgen(t_SampleInModel);
  
  setIsDropMissingDosagesInBgen(t_isDropmissingdosagesInBgen);
  
  setIsSparseDosageInBgen(t_isSparseDosageInBgen);
  m_AlleleOrder = t_AlleleOrder;
}


void BgenClass::setBgenObj(const std::string t_bgenFileName,
                           const std::string t_bgenFileIndex,
                           std::vector<std::string> & t_SampleInBgen)
{
  m_isQuery = false;
  // if(t_bgenFileIndex == ""){
  //   m_isQuery = false;
  //   Rcpp::Rcout << "    no index file for bgen is provided" << std::endl;
  // }  
  
  /****code from BOLT-LMM v2.3.4***/
  
  /********** READ HEADER v1.2**********/
  m_fin = fopen(t_bgenFileName.c_str(), "rb");
  uint32_t offset; size_t bytesRead = fread(&offset, 4, 1, m_fin); (void)bytesRead; //cout << "offset: " << offset << endl;
  uint32_t L_H; bytesRead = fread(&L_H, 4, 1, m_fin); (void)bytesRead; //cout << "L_H: " << L_H << endl;
  bytesRead = fread(&m_M0, 4, 1, m_fin); (void)bytesRead; //Rcpp::Rcout << "    snpBlocks (Mbgen): " << m_M0 << std::endl;
  assert(m_M0 != 0);
  //unsigned int Nbgen; fread(&Nbgen, 4, 1, m_fin); Rcpp::Rcout << "    samples (Nbgen): " << Nbgen << std::endl;
  bytesRead = fread(&m_N0, 4, 1, m_fin); (void)bytesRead; //Rcpp::Rcout << "    samples (Nbgen): " << m_N0 << std::endl;
  unsigned int m_Nsample = t_SampleInBgen.size();
  m_SampleInBgen = t_SampleInBgen;
  if (m_N0 != m_Nsample) {
    Rcpp::stop("ERROR: Number of samples in BGEN header does not match sample file");
  }
  char magic[5]; bytesRead = fread(magic, 1, 4, m_fin); (void)bytesRead; magic[4] = '\0'; //cout << "magic bytes: " << string(magic) << endl;
  fseek(m_fin, L_H-20, SEEK_CUR); //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
  uint32_t flags; bytesRead = fread(&flags, 4, 1, m_fin); (void)bytesRead; //cout << "flags: " << flags << endl;
  uint32_t CompressedSNPBlocks = flags&3; //Rcpp::Rcout << "    CompressedSNPBlocks: " << CompressedSNPBlocks << std::endl;
  assert(CompressedSNPBlocks==1); // REQUIRE CompressedSNPBlocks==1
  uint32_t Layout = (flags>>2)&0xf; //Rcpp::Rcout << "    Layout: " << Layout << std::endl;
  assert(Layout==1 || Layout==2); // REQUIRE Layout==1 or Layout==2
  fseek(m_fin, offset+4, SEEK_SET);
}


void BgenClass::setPosSampleInBgen(std::vector<std::string> & t_SampleInModel)
{
  // Rcpp::Rcout << "    Setting positions in BGEN files ..." << std::endl;	  
  m_N = t_SampleInModel.size();
  
  // updated by BWJ on 03/14/2021
  
  Rcpp::CharacterVector SampleInBgen(m_N0);
  for(uint32_t i = 0; i < m_N0; i++)
    SampleInBgen(i) = m_SampleInBgen.at(i);
  
  Rcpp::CharacterVector SampleInModel(m_N);
  for(uint32_t i = 0; i < m_N; i++)
    SampleInModel(i) = t_SampleInModel.at(i);
  
  Rcpp::IntegerVector posSampleInBgen = Rcpp::match(SampleInModel, SampleInBgen);
  for(uint32_t i = 0; i < m_N; i++){
    if(Rcpp::IntegerVector::is_na(posSampleInBgen.at(i)))
      Rcpp::stop("At least one subject requested is not in BGEN file.");
  }
  
  Rcpp::IntegerVector posSampleInModel = Rcpp::match(SampleInBgen, SampleInModel);
  m_posSampleInModel.resize(m_N0);
  for(uint32_t i = 0; i < m_N0; i++){
    if(Rcpp::IntegerVector::is_na(posSampleInModel.at(i))){
      m_posSampleInModel.at(i) = -1;
    }else{
      m_posSampleInModel.at(i) = posSampleInModel.at(i) - 1;   // convert "starting from 1" to "starting from 0"
    }
  }
  
  // end of the update on 03/14/2021
}


void BgenClass::Parse2(unsigned char *buf, 
                       unsigned int bufLen, 
                       const unsigned char *zBuf, 
                       unsigned int zBufLen,
                       std::string & snpName,
                       std::vector< double > & dosages, 
                       double & AC, 
                       double & AF, 
                       std::vector<uint32_t> & indexforMissing, 
                       double & info, 
                       std::vector<unsigned int> & indexNonZero)
{
  uLong destLen = bufLen;
  if (uncompress(buf, &destLen, zBuf, zBufLen) != Z_OK || destLen != bufLen) {
    if(uncompress(buf, &destLen, zBuf, zBufLen) == Z_BUF_ERROR)
      Rcpp::stop("ERROR: uncompress() failed: The buffer dest was not large enough to hold the uncompressed data.");
    
    if(uncompress(buf, &destLen, zBuf, zBufLen) == Z_MEM_ERROR)
      Rcpp::stop("ERROR: uncompress() failed: Insufficient memory.");
    
    if(uncompress(buf, &destLen, zBuf, zBufLen) == Z_DATA_ERROR)
      Rcpp::stop("ERROR: uncompress() failed: The compressed data (referenced by source) was corrupted.");
    
    Rcpp::stop("ERROR: uncompress() failed");
  }
  
  unsigned char *bufAt = buf;
  unsigned int N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;
  
  if (N != m_N0) {
    Rcpp::stop("ERROR: " + snpName + " has N = " + std::to_string(N) + " (mismatch with header block)");
  }
  unsigned int K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
  if (K != 2U) {
    Rcpp::stop("ERROR: " + snpName + " has K = " + std::to_string(K) + " (non-bi-allelic)");
  }
  unsigned int Pmin = *bufAt; bufAt++;
  if (Pmin != 2U) {
    Rcpp::stop("ERROR: " + snpName + " has minimum ploidy = " + std::to_string(Pmin) + " (not 2)");
  }
  unsigned int Pmax = *bufAt; bufAt++;
  if (Pmax != 2U) {
    Rcpp::stop("ERROR: " + snpName + " has maximum ploidy = " + std::to_string(Pmax) + " (not 2)");
  }
  
  const unsigned char *ploidyMissBytes = bufAt;
  for (unsigned int i = 0; i < N; i++) {
    unsigned int ploidyMiss = *bufAt; bufAt++;
    if (ploidyMiss != 2U && ploidyMiss != 130U) {
      Rcpp::stop("ERROR: " + snpName + " has ploidy/missingness byte = " + std::to_string(ploidyMiss) + " (not 2 or 130)");
    }
  }
  unsigned int Phased = *bufAt; bufAt++;
  if (Phased != 0U) {
    Rcpp::stop("ERROR: " + snpName + " has Phased = " + std::to_string(Pmax) + " (not 0)");
  }
  unsigned int B = *bufAt; bufAt++;
  if (B != 8U) {
    Rcpp::stop("ERROR: " + snpName + " has B = " + std::to_string(B) + " (not 8)");
  }
  double lut[256];
  for (int i = 0; i <= 255; i++)
    lut[i] = i/255.0;
  
  double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0; // sum_fij_minus_eij2_sub = 0; // for INFO
  double p11,p10,dosage,eij,fij;  // p00, eijsub, fijsub
  dosages.clear();
  dosages.reserve(m_N);
  if(!m_isSparseDosagesInBgen){
    dosages.resize(m_N);
  }
  std::size_t missing_cnt = 0;
  
  for (unsigned int i = 0; i < N; i++) {
    // if(i == 1){Rcpp::Rcout << "    ploidyMissBytes[i] " << ploidyMissBytes[i] << std::endl;}
    if (ploidyMissBytes[i] != 130U){
      // bufAt += 2;
      p11 = lut[*bufAt]; bufAt++;
      p10 = lut[*bufAt]; bufAt++;
      // p00 = 1 - p11 - p10; //can remove
      dosage = 2*p11 + p10;
      
      eij = dosage;
      fij = 4*p11 + p10;
      sum_eij += eij;
      sum_fij_minus_eij2 += fij - eij*eij;
      if(m_posSampleInModel[i] >= 0){
        if(!m_isSparseDosagesInBgen){
          // dosages[m_posSampleInModel[i]] = 2 - dosage;
          // updated by BWJ on 2021-02-28
          // dosages[m_posSampleInModel[i]] = dosage;
          // changed back by BWJ on 2021-03-14 (change default setting to "ref-first")
          dosages[m_posSampleInModel[i]] = 2 - dosage;
        }else{
          if(2 - dosage > 0){
            // dosages.push_back(2 - dosage);
            // updated by BWJ on 2021-02-28
            dosages.push_back(dosage);
            indexNonZero.push_back(m_posSampleInModel[i]+1);
          }
        }
        sum_eij_sub += eij;
      }
    }else if(ploidyMissBytes[i] == 130U){
      bufAt += 2;
      if(m_posSampleInModel[i] >= 0){
        indexforMissing.push_back(m_posSampleInModel[i]);
        ++missing_cnt;
        if(!m_isSparseDosagesInBgen){
          dosages[m_posSampleInModel[i]] = -1;
        }
      }
    }
  }
  // Rcpp::Rcout << "    sum_eij_sub: " << sum_eij_sub << std::endl;
  AC = 2* ((double) (m_N - missing_cnt)) - sum_eij_sub;
  if(m_N == missing_cnt){
    AF = 0;
  }else{
    AF = AC/ 2/ ((double) (m_N - missing_cnt)) ;
  }
  
  double thetaHat = sum_eij / (2* (m_N - missing_cnt));
  // Rcpp::Rcout << "    sum_eij " << sum_eij << std::endl;
  // Rcpp::Rcout << "    missing_cnt " << sum_eij << std::endl;
  info = thetaHat==0 || thetaHat==1 ? 1 :
    1 - sum_fij_minus_eij2 / (2*(m_N - missing_cnt)*thetaHat*(1-thetaHat));
  
  // updated on 2021-08-27: the imputation part was put in UTIL.cpp
  // if(missing_cnt > 0){
  //   // Rcpp::Rcout << "    AC: " << AC << std::endl;
  //   // Rcpp::Rcout << "    sample index with missing dosages for snpName " << snpName << " :";
  //   
  //   if(!m_isDropMissingDosagesInBgen){
  //     double imputeDosage = 2*AF;
  //     for (unsigned int i = 0; i < indexforMissing.size(); i++)
  //     {
  //       // Rcpp::Rcout << indexforMissing[i]+1 << ",";
  //       if(!m_isSparseDosagesInBgen){
  //         dosages[indexforMissing[i]] = imputeDosage;
  //       }else{
  //         dosages.push_back(imputeDosage);
  //         indexNonZero.push_back(indexforMissing[i]+1);
  //       }
  //       AC = AC + imputeDosage;
  //     }
  //     // Rcpp::Rcout << "    AC new: " << AC << std::endl;
  //   }
  // }
}

arma::vec BgenClass::getOneMarker(uint64_t t_gIndex,        // different meanings for different genoType
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
                                  std::vector<uint32_t>& t_indexForNonZero,
                                  bool& t_isBoolRead)        // only used in BGEN, Wei, if you want, you can add details here.
{
  if(t_gIndex > 0){fseek(m_fin, t_gIndex, SEEK_SET);}
  std::string SNPID, RSID, chromosome, first_allele, second_allele;
  uint32_t position;
  std::vector< std::string > alleles ;
  std::vector< double > dosages;
  double AC, AF, info;
  std::vector<uint32_t> indexforMissing;
  std::vector< unsigned int > indexNonZero;
  char snpID[65536], rsID[65536], chrStr[65536];
  unsigned int maxLA = 65536, maxLB = 65536;
  char *allele1, *allele0;
  allele1 = (char *) malloc(maxLA+1);
  allele0 = (char *) malloc(maxLB+1);
  uint16_t LS; size_t numBoolRead = fread(&LS, 2, 1, m_fin); // cout << "LS: " << LS << " " << std::flush;
  // bool isBoolRead;  // BWJ (2021-02-28): I think it will not be used since we currently use t_gIndex to specify bytes position. This bool value can still be outputted through a reference. If we are sure it is not useful any more, we can remove it.
  // Rcpp::List result;
  if ( numBoolRead > 0 ) {
    // isBoolRead = true;
    t_isBoolRead = true;
    size_t bytesRead = fread(snpID, 1, LS, m_fin); (void)bytesRead; snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;
    uint16_t LR; bytesRead = fread(&LR, 2, 1, m_fin); (void)bytesRead; // cout << "LR: " << LR << " " << std::flush;
    bytesRead = fread(rsID, 1, LR, m_fin); (void)bytesRead; rsID[LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;
    RSID = std::string(rsID)=="." ? snpID : rsID;
    //std::string SNPID = string(snpID);
    
    uint16_t LC; bytesRead = fread(&LC, 2, 1, m_fin); (void)bytesRead; // cout << "LC: " << LC << " " << std::flush;
    bytesRead = fread(chrStr, 1, LC, m_fin); (void)bytesRead; chrStr[LC] = '\0';
    chromosome  = std::string(chrStr);
    
    uint32_t physpos; bytesRead = fread(&physpos, 4, 1, m_fin); (void)bytesRead; // cout << "physpos: " << physpos << " " << std::flush;
    position = physpos;
    uint16_t K; bytesRead = fread(&K, 2, 1, m_fin); (void)bytesRead; //cout << "K: " << K << endl;
    if (K != 2) {
      Rcpp::stop("ERROR: Non-bi-allelic variant found: " + std::to_string(K) + " alleles");
    }
    uint32_t LA; bytesRead = fread(&LA, 4, 1, m_fin); (void)bytesRead; // cout << "LA: " << LA << " " << std::flush;
    if (LA > maxLA) {
      maxLA = 2*LA;
      free(allele1);
      allele1 = (char *) malloc(maxLA+1);
    }
    bytesRead = fread(allele1, 1, LA, m_fin); (void)bytesRead; allele1[LA] = '\0';
    first_allele = std::string(allele1);
    free(allele1);
    
    uint32_t LB; bytesRead = fread(&LB, 4, 1, m_fin); (void)bytesRead; // cout << "LB: " << LB << " " << std::flush;
    if (LB > maxLB) {
      maxLB = 2*LB;
      free(allele0);
      allele0 = (char *) malloc(maxLB+1);
    }
    bytesRead = fread(allele0, 1, LB, m_fin); (void)bytesRead; allele0[LB] = '\0';
    second_allele = std::string(allele0);
    free(allele0);
    
    uint32_t C; bytesRead = fread(&C, 4, 1, m_fin); (void)bytesRead; //cout << "C: " << C << endl;
    if (C > m_zBuf.size()) m_zBuf.resize(C-4);
    //Rcpp::Rcout << "    m_zBuf.size() " << m_zBuf.size() << std::endl;
    uint32_t D; bytesRead = fread(&D, 4, 1, m_fin); (void)bytesRead; //cout << "D: " << D << endl;
    m_zBufLens = C-4; m_bufLens = D;
    bytesRead = fread(&m_zBuf[0], 1, C-4, m_fin); (void)bytesRead;
    AC = 0;
    AF = 0;
    info = 0;
    if (m_bufLens > m_buf.size()) m_buf.resize(m_bufLens); //fix the length
    
    Parse2(&m_buf[0], m_bufLens, &m_zBuf[0], m_zBufLens, RSID, dosages, AC, AF, indexforMissing, info, indexNonZero);
    
    // output
    // t_alt = first_allele;       // ALT allele (usually minor allele)
    // t_ref = second_allele;       // REF allele (usually major allele)
    t_alt = second_allele;  // default setting is "ref-first" (03-14-2021)
    t_ref = first_allele;
    t_marker = RSID;    // marker ID extracted from genotype file
    t_pd = position;           // base position
    t_chr = chromosome;       // chromosome
    t_altFreq = AF;        // frequency of ALT allele
    t_altCounts = AC;      // counts of ALT allele
    t_imputeInfo = info;     // imputation information score, i.e., R2 (all 1 for PLINK)
    t_indexForMissing = indexforMissing;     // index of missing genotype data
    t_missingRate = (double) t_indexForMissing.size() / (double) m_N;    // missing rate
    t_indexForNonZero = indexNonZero;
    
    if(m_AlleleOrder == "alt-first"){  // added by Wenjian Bi on 03/14/2021
      t_alt = first_allele;
      t_ref = second_allele;
      t_altFreq = 1 - t_altFreq;
      t_altCounts = t_altFreq * 2 * ((double)m_N - (double)t_indexForMissing.size());
      for(unsigned int i = 0; i < dosages.size(); i++)
        dosages.at(i) = 2 - dosages.at(i);
    }
    
  }else{
    // isBoolRead = false;
    t_isBoolRead = false;
    // Rcpp::DataFrame variants = NULL;
    // result["isBoolRead"] = isBoolRead;
  }
  // dosages.clear();
  // indexforMissing.clear();
  return dosages;
}


void BgenClass::setIsSparseDosageInBgen (bool t_isSparseDosageInBgen){
  m_isSparseDosagesInBgen = t_isSparseDosageInBgen;
}

void BgenClass::setMarkerIndicesToIncludeInBgen (std::vector< int > & t_markerIndicesToInclude){
  m_markerIndicesToInclude = t_markerIndicesToInclude;
}

void BgenClass::setIsDropMissingDosagesInBgen (bool t_isDropmissingdosagesInBgen){
  m_isDropMissingDosagesInBgen = t_isDropmissingdosagesInBgen;
}


}  
