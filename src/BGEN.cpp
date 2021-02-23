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
#include <memory>
#include <time.h>
#include <stdint.h>
#include <zlib.h>

// #include <boost/iostreams/filter/zstd.hpp>
// #include "zstd.h"
#include <boost/date_time.hpp>

#include <Rcpp.h>

#include "BGEN.hpp"

namespace BGEN {

BgenClass::BgenClass(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string> t_SampleInBgen,
                     std::vector<std::string> t_SampleInModel,
                     bool t_isSparseDosageInBgen,
                     bool t_isDropmissingdosagesInBgen)
{
  setBgenObj(t_bgenFileName,
             t_bgenFileIndex,
             t_SampleInBgen);
  
  setPosSampleInBgen(t_SampleInModel);
  
  setIsDropMissingDosagesInBgen(t_isDropmissingdosagesInBgen);
  
  setIsSparseDosageInBgen(t_isSparseDosageInBgen);
  
}


void BgenClass::setBgenObj(const std::string t_bgenFileName,
                           const std::string t_bgenFileIndex,
                           std::vector<std::string> & t_SampleInBgen)
{
  m_isQuery = false;
  if(t_bgenFileIndex == ""){
    m_isQuery = false;
    std::cout << "no index file for bgen is provided" << std::endl;
  }  
  
  /****code from BOLT-LMM v2.3.4***/
  
  /********** READ HEADER v1.2**********/
  m_fin = fopen(t_bgenFileName.c_str(), "rb");
  uint32_t offset; fread(&offset, 4, 1, m_fin); //cout << "offset: " << offset << endl;
  uint32_t L_H; fread(&L_H, 4, 1, m_fin); //cout << "L_H: " << L_H << endl;
  fread(&m_M0, 4, 1, m_fin); std::cout << "snpBlocks (Mbgen): " << m_M0 << std::endl;
  assert(m_M0 != 0);
  //unsigned int Nbgen; fread(&Nbgen, 4, 1, m_fin); std::cout << "samples (Nbgen): " << Nbgen << std::endl;
  fread(&m_N0, 4, 1, m_fin); std::cout << "samples (Nbgen): " << m_N0 << std::endl;
  unsigned int m_Nsample = t_SampleInBgen.size();
  m_SampleInBgen = t_SampleInBgen;
  if (m_N0 != m_Nsample) {
    std::cerr << "ERROR: Number of samples in BGEN header does not match sample file" << std::endl;
    exit(1);
  }
  char magic[5]; fread(magic, 1, 4, m_fin); magic[4] = '\0'; //cout << "magic bytes: " << string(magic) << endl;
  fseek(m_fin, L_H-20, SEEK_CUR); //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
  uint32_t flags; fread(&flags, 4, 1, m_fin); //cout << "flags: " << flags << endl;
  uint32_t CompressedSNPBlocks = flags&3; std::cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << std::endl;
  assert(CompressedSNPBlocks==1); // REQUIRE CompressedSNPBlocks==1
  uint32_t Layout = (flags>>2)&0xf; std::cout << "Layout: " << Layout << std::endl;
  assert(Layout==1 || Layout==2); // REQUIRE Layout==1 or Layout==2
  fseek(m_fin, offset+4, SEEK_SET);
}


void BgenClass::setPosSampleInBgen(std::vector<std::string> & t_SampleInModel)
{
  std::cout << "Setting position of samples in Bgen files...." << std::endl;	  
  m_N = t_SampleInModel.size();
  
  for(uint32_t i = 0; i < m_N; i++){
    std::string sample = t_SampleInModel.at(i);
    auto pos = std::find(m_SampleInBgen.begin(), m_SampleInBgen.end(), sample);
    if(pos == m_SampleInBgen.end()){
      Rcpp::stop("At least one subject requested is not in Bgen file.");
    }
  }
  
  m_posSampleInModel.clear();
  for(uint32_t i = 0; i < m_N0; i++){
    std::string sample = m_SampleInBgen.at(i);
    auto pos = std::find(t_SampleInModel.begin(), t_SampleInModel.end(), sample);
    if(pos != t_SampleInModel.end()){
      m_posSampleInModel.push_back(pos - t_SampleInModel.begin());
    }else{
      m_posSampleInModel.push_back(-1);      
    }
  }
  
}


void BgenClass::Parse2(unsigned char *buf, unsigned int bufLen, const unsigned char *zBuf, unsigned int zBufLen,std::string & snpName,std::vector< double > & dosages, double & AC, double & AF, std::vector<int> & indexforMissing, double & info, std::vector<unsigned int> & iIndex) {
  
  uLong destLen = bufLen;
  if (uncompress(buf, &destLen, zBuf, zBufLen) != Z_OK || destLen != bufLen) {
    std::cerr << "ERROR: uncompress() failed" << std::endl;
    exit(1);
  }
  
  unsigned char *bufAt = buf;
  unsigned int N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;
  
  if (N != m_N0) {
    std::cerr << "ERROR: " << snpName << " has N = " << N << " (mismatch with header block)" << std::endl;
    exit(1);
  }
  unsigned int K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
  if (K != 2U) {
    std::cerr << "ERROR: " << snpName << " has K = " << K << " (non-bi-allelic)" << std::endl;
    exit(1);
  }
  unsigned int Pmin = *bufAt; bufAt++;
  if (Pmin != 2U) {
    std::cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin << " (not 2)" << std::endl;
    exit(1);
  }
  unsigned int Pmax = *bufAt; bufAt++;
  if (Pmax != 2U) {
    std::cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax << " (not 2)" << std::endl;
    exit(1);
  }
  
  const unsigned char *ploidyMissBytes = bufAt;
  for (unsigned int i = 0; i < N; i++) {
    unsigned int ploidyMiss = *bufAt; bufAt++;
    if (ploidyMiss != 2U && ploidyMiss != 130U) {
      std::cerr << "ERROR: " << snpName << " has ploidy/missingness byte = " << ploidyMiss
                << " (not 2 or 130)" << std::endl;
      exit(1);
    }
  }
  unsigned int Phased = *bufAt; bufAt++;
  if (Phased != 0U) {
    std::cerr << "ERROR: " << snpName << " has Phased = " << Pmax << " (not 0)" << std::endl;
    exit(1);
  }
  unsigned int B = *bufAt; bufAt++;
  if (B != 8U) {
    std::cerr << "ERROR: " << snpName << " has B = " << B << " (not 8)" << std::endl;
    exit(1);
  }
  double lut[256];
  for (int i = 0; i <= 255; i++)
    lut[i] = i/255.0;
  
  double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0, sum_fij_minus_eij2_sub = 0; // for INFO
  double p11,p10,p00,dosage,eij,fij, eijsub, fijsub;
  dosages.clear();
  dosages.reserve(m_N);
  if(!m_isSparseDosagesInBgen){
    dosages.resize(m_N);
  }
  std::size_t missing_cnt = 0;
  
  for (unsigned int i = 0; i < N; i++) {
    //if(i == 1){std::cout << "ploidyMissBytes[i] " << ploidyMissBytes[i] << std::endl;}
    if (ploidyMissBytes[i] != 130U){
      //bufAt += 2;
      p11 = lut[*bufAt]; bufAt++;
      p10 = lut[*bufAt]; bufAt++;
      p00 = 1 - p11 - p10; //can remove
      dosage = 2*p11 + p10;
      
      eij = dosage;
      fij = 4*p11 + p10;
      sum_eij += eij;
      sum_fij_minus_eij2 += fij - eij*eij;
      if(m_posSampleInModel[i] >= 0){
        if(!m_isSparseDosagesInBgen){
          dosages[m_posSampleInModel[i]] = 2 - dosage;
        }else{
          if(2 - dosage > 0){
            dosages.push_back(2 - dosage);
            iIndex.push_back(m_posSampleInModel[i]+1);
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
  //std::cout << "sum_eij_sub: " << sum_eij_sub << std::endl;
  AC = 2* ((double) (m_N - missing_cnt)) - sum_eij_sub;
  if(m_N == missing_cnt){
    AF = 0;
  }else{
    AF = AC/ 2/ ((double) (m_N - missing_cnt)) ;
  }
  
  double thetaHat = sum_eij / (2* (m_N - missing_cnt));
  //std::cout << "sum_eij " << sum_eij << std::endl;
  //std::cout << "missing_cnt " << sum_eij << std::endl;
  info = thetaHat==0 || thetaHat==1 ? 1 :
    1 - sum_fij_minus_eij2 / (2*(m_N - missing_cnt)*thetaHat*(1-thetaHat));
  
  if(missing_cnt > 0){
    std::cout << "AC: " << AC << std::endl;
    std::cout << "sample index with missing dosages for snpName " << snpName << " :";
    
    if(!m_isDropMissingDosagesInBgen){
      double imputeDosage = 2*AF;
      for (unsigned int i = 0; i < indexforMissing.size(); i++)
      {
        std::cout << indexforMissing[i]+1 << ",";
        if(!m_isSparseDosagesInBgen){
          dosages[indexforMissing[i]] = imputeDosage;
        }else{
          dosages.push_back(imputeDosage);
          iIndex.push_back(indexforMissing[i]+1);
        }
        AC = AC + imputeDosage;
      }
      std::cout << "AC new: " << AC << std::endl;
    }
  }
  
}


// get dosages/genotypes of one marker
Rcpp::List BgenClass::getOneMarker(int t_fileStartPos)
{
  if(t_fileStartPos > 0){fseek(m_fin, t_fileStartPos, SEEK_SET);}
  std::string SNPID, RSID, chromosome, first_allele,second_allele ;
  unsigned int position;
  std::vector< std::string > alleles ;
  std::vector< double > dosages;
  double AC, AF, info;
  std::vector< int > indexforMissing;
  std::vector< unsigned int > iIndex;
  char snpID[65536], rsID[65536], chrStr[65536];
  unsigned int maxLA = 65536, maxLB = 65536;
  char *allele1, *allele0;
  allele1 = (char *) malloc(maxLA+1);
  allele0 = (char *) malloc(maxLB+1);
  uint16_t LS; size_t numBoolRead = fread(&LS, 2, 1, m_fin); // cout << "LS: " << LS << " " << std::flush;
  bool isBoolRead;
  Rcpp::List result ;
  if ( numBoolRead > 0 ) {
    isBoolRead = true;
    fread(snpID, 1, LS, m_fin); snpID[LS] = '\0'; // cout << "snpID: " << string(snpID) << " " << std::flush;
    uint16_t LR; fread(&LR, 2, 1, m_fin); // cout << "LR: " << LR << " " << std::flush;
    fread(rsID, 1, LR, m_fin); rsID[LR] = '\0'; // cout << "rsID: " << string(rsID) << " " << std::flush;
    RSID = std::string(rsID)=="." ? snpID : rsID;
    //std::string SNPID = string(snpID);
    
    uint32_t LC; fread(&LC, 2, 1, m_fin); // cout << "LC: " << LC << " " << std::flush;
    std::cout << "11" << std::endl;
    std::cout << LC << std::endl;
    fread(chrStr, 1, LC, m_fin); chrStr[LC] = '\0';
    std::cout << "22" << std::endl;
    chromosome  = std::string(chrStr);
    
    uint32_t physpos; fread(&physpos, 4, 1, m_fin); // cout << "physpos: " << physpos << " " << std::flush;
    position = physpos;
    uint16_t K; fread(&K, 2, 1, m_fin); //cout << "K: " << K << endl;
    if (K != 2) {
      std::cerr << "ERROR: Non-bi-allelic variant found: " << K << " alleles" << std::endl;
      exit(1);
    }
    uint32_t LA; fread(&LA, 4, 1, m_fin); // cout << "LA: " << LA << " " << std::flush;
    if (LA > maxLA) {
      maxLA = 2*LA;
      free(allele1);
      allele1 = (char *) malloc(maxLA+1);
    }
    fread(allele1, 1, LA, m_fin); allele1[LA] = '\0';
    second_allele = std::string(allele1);
    uint32_t LB; fread(&LB, 4, 1, m_fin); // cout << "LB: " << LB << " " << std::flush;
    if (LB > maxLB) {
      maxLB = 2*LB;
      free(allele0);
      allele0 = (char *) malloc(maxLB+1);
    }
    fread(allele0, 1, LB, m_fin); allele0[LB] = '\0';
    first_allele = std::string(allele0);
    
    uint32_t C; fread(&C, 4, 1, m_fin); //cout << "C: " << C << endl;
    if (C > m_zBuf.size()) m_zBuf.resize(C-4);
    //std::cout << "m_zBuf.size() " << m_zBuf.size() << std::endl;
    uint32_t D; fread(&D, 4, 1, m_fin); //cout << "D: " << D << endl;
    m_zBufLens = C-4; m_bufLens = D;
    fread(&m_zBuf[0], 1, C-4, m_fin);
    AC = 0;
    AF = 0;
    info = 0;
    if (m_bufLens > m_buf.size()) m_buf.resize(m_bufLens); //fix the length
    Parse2(&m_buf[0], m_bufLens, &m_zBuf[0], m_zBufLens, RSID, dosages, AC, AF, indexforMissing, info, iIndex);
    Rcpp::DataFrame variants = Rcpp::DataFrame::create(
      Rcpp::Named("chromosome") = chromosome,
      Rcpp::Named("position") = position,
      Rcpp::Named("rsid") = RSID,
      //        Named("number_of_alleles") = number_of_allele,
      Rcpp::Named("allele0") = first_allele,
      Rcpp::Named("allele1") = second_allele,
      //       _["stringsAsFactors"] = false,
      Rcpp::Named("AC") = AC,
      Rcpp::Named("AF") = AF
      // Named("info") = info
      //Named("homN_cases") = homN_cases,
      //Named("hetN_cases") = hetN_cases,
      //Named("homN_ctrls") = homN_ctrls,
      //Named("hetN_ctrls") = hetN_ctrls
    );
    result[ "variants" ] = variants ;
    result["info"] = info;
    result[ "dosages" ] = dosages ;
    result["iIndex"] = iIndex;
    result["indexforMissing"] = indexforMissing;
    
  }else{
    isBoolRead = false;
    Rcpp::DataFrame variants = NULL;
    result["isBoolRead"] = isBoolRead;
  }
  dosages.clear();
  indexforMissing.clear();
  return(result);
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
