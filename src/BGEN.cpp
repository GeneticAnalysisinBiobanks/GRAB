#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "../thirdParty/bgen/genfile/include/genfile/bgen/bgen.hpp"
#include "../thirdParty/bgen/genfile/include/genfile/bgen/View.hpp"
#include "../thirdParty/bgen/genfile/include/genfile/bgen/IndexQuery.hpp"
#include <sstream>
#include <time.h>
#include <stdint.h>

#include <Rcpp.h>


#include "BGEN.hpp"


namespace BGEN {

  BgenClass::BgenClass(std::string t_bgenFileName,
            std::string t_bgenFileIndex,
            Rcpp::DataFrame t_rangesInclude,
            Rcpp::DataFrame t_rangesExclude,
            std::vector< std::string > t_markerIDsToExclude,
            std::vector< std::string > t_markerIDsToInclude,
            std::vector< std::string > t_SampleInBgen,
            bool t_isDropMissingDosages,
            double t_minMAF,
            double t_minInfo,
            double t_maxMAF,
            std::vector<std::string> t_SampleInModel)
  {
    std::vector< double > t_dosages;	  
    setBgenObj(t_bgenFileName,
	       t_bgenFileIndex,
	       t_rangesInclude,
	       t_rangesExclude,
	       t_markerIDsToExclude,
	       t_markerIDsToInclude,
	       t_SampleInBgen,
	       t_isDropMissingDosages,
	       t_minMAF,
	       t_minInfo,
	       t_maxMAF,
	       t_dosages);

    setPosSampleInBgen(t_SampleInModel);

  }


  void BgenClass::setBgenObj(const std::string t_bgenFileName,
                 const std::string t_bgenFileIndex,
		 Rcpp::DataFrame t_rangesInclude,
  		 Rcpp::DataFrame t_rangesExclude,
                 std::vector< std::string > t_markerIDsToExclude,
                 std::vector< std::string > t_markerIDsToInclude,
		 std::vector< std::string > t_SampleInBgen,
		 bool t_isDropMissingDosages,
		 double t_minMAF,
                 double t_minInfo,
		 double t_maxMAF,
		 std::vector< double > & t_dosages)
  {

    m_minMAF = t_minMAF;
    m_minInfo = t_minInfo;
    m_maxMAF = t_maxMAF;
    m_SampleInBgen = t_SampleInBgen;
    m_isDropMissingDosages = t_isDropMissingDosages;
    m_bgenFileName = t_bgenFileName;
    m_bgenFileIndex = t_bgenFileIndex;

    m_isQuery = false;
    if(t_bgenFileIndex == ""){
     m_isQuery = false;
     std::cout << "no index file for bgen is provided" << std::endl;
    }else if(t_rangesInclude.nrow() == 0 &&
     t_rangesExclude.nrow() == 0 &&
     t_markerIDsToExclude.size() == 0 &&
     t_markerIDsToInclude.size() == 0){
     std::cout << "no query list is provided" << std::endl;
     m_isQuery = false;
    }else{
     m_isQuery = true;
    }


    if(m_isQuery){
      m_genoToTestBgen = genfile::bgen::View::create( t_bgenFileName ) ;
      genfile::bgen::IndexQuery::UniquePtr m_query = genfile::bgen::IndexQuery::create( t_bgenFileIndex ) ;
      std::cout << "t_rangesInclude.nrow() " << t_rangesInclude.nrow() << std::endl;
      std::cout << "t_rangesExclude.nrow() " << t_rangesExclude.nrow() << std::endl;
      std::cout << "t_markerIDsToExclude.size() " << t_markerIDsToExclude.size() << std::endl;
      std::cout << "t_markerIDsToInclude.size() " << t_markerIDsToInclude.size() << std::endl;

      if (t_rangesInclude.nrow() > 0){
	Rcpp::StringVector const& m_chromosome = t_rangesInclude["chromosome"] ;
        Rcpp::IntegerVector const& m_start = t_rangesInclude["start"] ;
        Rcpp::IntegerVector const& m_end = t_rangesInclude["end"] ;
        for( int i = 0; i < t_rangesInclude.nrow(); ++i){
        
	  //if( m_end[i] < m_start[i] ) {
          //  throw std::invalid_argument( "Range (" + m_chromosome[i] + ":" + atoi( m_start[i] ) + "-" + std::atoi( m_end[i] ) + ") is malformed." ) ;
          //}
          m_query->include_range( genfile::bgen::IndexQuery::GenomicRange( std::string( m_chromosome[i] ), m_start[i], m_end[i] )) ;
        }
      }

      if (t_rangesExclude.nrow() > 0){

        Rcpp::StringVector const& m_chromosome_exclude = t_rangesExclude["chromosome"] ;
        Rcpp::IntegerVector const& m_start_exclude = t_rangesExclude["start"] ;
        Rcpp::IntegerVector const& m_end_exclude = t_rangesExclude["end"] ;
        for( int i = 0; i < t_rangesExclude.nrows(); ++i ) {
          //if( m_end_exclude[i] < m_start_exclude[i] ) {
          //  throw std::invalid_argument( "Range (" + m_chromosome_exclude[i] + ":" + std::atoi( m_start_exclude[i] ) + "-" + std::atoi( m_end_exclude[i] ) + ") is malformed." ) ;
          //}
          m_query->exclude_range( genfile::bgen::IndexQuery::GenomicRange( std::string( m_chromosome_exclude[i] ), m_start_exclude[i], m_end_exclude[i] )) ;
        }
      }

      if (t_markerIDsToInclude.size() != 0){
        m_query->include_rsids(t_markerIDsToInclude);
      }

      if (t_markerIDsToInclude.size() != 0){
        m_query->exclude_rsids(t_markerIDsToInclude);
      }

      m_query->initialise() ;

      if(m_query->number_of_variants() > 0){
        m_genoToTestBgen->set_query( m_query ) ;
        m_M0 = m_genoToTestBgen->number_of_variants() ;
        m_N0 = m_genoToTestBgen->number_of_samples();
        std::cout << m_M0 << " markers will be analyzed " << std::endl;
        //return m_M0 ;
      }else{
        std::cout << "No queried variant is found in the bgen file! All variants bgen file will be analyzed" << std::endl;
        m_isQuery = false;
      }
    }


    if(!m_isQuery){

      m_gm_stream.reset(
        new std::ifstream( t_bgenFileName.c_str(), std::ifstream::binary )
      ) ;

      if( !*m_gm_stream ) {
        throw std::invalid_argument( t_bgenFileName ) ;
      }

      m_gm_stream->seekg( 0, std::ios::beg ) ;
      genfile::bgen::read_offset( *m_gm_stream, &m_gm_offset ) ;
      genfile::bgen::read_header_block( *m_gm_stream, &m_gm_context ) ;
      uint m_N0_temp = m_gm_context.number_of_samples;
      m_N0 = int(m_N0_temp);
      std::cout << m_N0 << " samples are found in the bgen file" << std::endl;

      // Jump to the first variant data block.
      m_gm_stream->seekg( m_gm_offset + 4 ) ;
      //printf("4\n");fflush(NULL);
      uint m_M0_temp = m_gm_context.number_of_variants;
      m_M0 = m_M0_temp;
      //std::cout << "All " << numMarkers << " markers will be analyzed " << std::endl;
      std::cout << m_M0 << " markers are found in the bgen file " << std::endl;
      //return m_M0 ;
    } 
  }



  void BgenClass::setPosSampleInBgen(std::vector<std::string> & t_SampleInModel)
  {	  
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

  // get dosages/genotypes of one marker
  bool BgenClass::getOneMarker(uint32_t t_posMarker,
		         double& t_freq,
			 double& t_ac,
			 double& t_info,
                         double& t_missingRate,
                         std::vector<uint32_t>& t_posMissingGeno,
                         std::string& t_a1,
                         std::string& t_a2,
                         std::string& t_marker,
			 std::string& t_rsid,
                         genfile::bgen::uint32_t& t_pd,
                         std::string& t_chrMarker,
			 std::vector<double> & t_dosage)
  {
 
    std::vector< genfile::byte_t > m_buffer2;	  

    if(!m_isQuery){	  
      std::vector< genfile::byte_t > m_buffer1;
      m_isReadVariant = genfile::bgen::read_snp_identifying_data(
                        *m_gm_stream,
                        m_gm_context,
                        &t_marker,
                        &t_rsid,
                        &t_chrMarker,
                        &t_pd,
                        &t_a1,
                        &t_a2);

      genfile::bgen::read_genotype_data_block(
                        *m_gm_stream,
                        m_gm_context,
                        &m_buffer1);

      genfile::bgen::uncompress_probability_data(
                        m_gm_context,
                        m_buffer1,
                        &m_buffer2);

    }else{
      std::vector< std::string > m_alleles;	    
      m_isReadVariant = m_genoToTestBgen->read_variant(&t_marker, &t_rsid, &t_chrMarker, &t_pd, &m_alleles) ;
      m_buffer2 = m_genoToTestBgen->read_and_uncompress_genotype_data_block();      
      t_a1 = m_alleles[0];
      t_a2 = m_alleles[1];
    }	   

    unsigned char * buf  = (unsigned char *) m_buffer2.data();
    uint Nbgen = m_gm_context.number_of_samples;
    t_ac = 0;
    t_freq = 0;
    t_info = Parse(buf, m_buffer2.size(), t_marker, Nbgen, t_dosage, t_ac, t_freq, t_posMissingGeno);
    return(m_isReadVariant);    	  
  }



  double  BgenClass::Parse(unsigned char * buf, size_t bufLen,  std::string & snpName, uint Nbgen,std::vector< double > & dosages, double & AC, double & AF, std::vector<unsigned int> & indexforMissing){

    size_t destLen = bufLen;
    unsigned char * bufAt = buf;
    uint N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;

    if (N != Nbgen) {
      std::cerr << "ERROR: " << snpName << " has N = " << N << " (mismatch with header block)" << std::endl;
      exit(1);
    }
    uint K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
    if (K != 2U) {
      std::cerr << "ERROR: " << snpName << " has K = " << K << " (non-bi-allelic)" << std::endl;
      exit(1);
    }
    uint Pmin = *bufAt; bufAt++;
    if (Pmin != 2U) {
      std::cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin << " (not 2)" << std::endl;
      exit(1);
    }
    uint Pmax = *bufAt; bufAt++;
    if (Pmax != 2U) {
      std::cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax << " (not 2)" << std::endl;
      exit(1);
    }

    //deal with missing dosages
    std::vector <bool> missingIdxVec;
    missingIdxVec.clear();
    missingIdxVec.reserve(N);
    missingIdxVec.resize(N);
    int missingSamplesize = 0;

    for (uint i = 0; i < N; i++) {
      uint ploidyMiss = *bufAt; bufAt++;
      bool const missing = (ploidyMiss & 0x80) ;
      missingIdxVec[i] = missing;
        if(missing){
          missingSamplesize = missingSamplesize + 1;
        }
    }
    uint Phased = *bufAt; bufAt++;
    if (Phased != 0U) {
      std::cerr << "ERROR: " << snpName << " has Phased = " << Pmax << " (not 0)" << std::endl;
      exit(1);
    }
    uint B = *bufAt; bufAt++;
    if (B != 8U) {
      std::cerr << "ERROR: " << snpName << " has B = " << B << " (not 8)" << std::endl;
      exit(1);
    }

        // Parse
    double lut[256];
    for (int i = 0; i <= 255; i++)
      lut[i] = i/255.0;

    double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0, sum_fij_minus_eij2_sub = 0; // for INFO
    double p11,p10,p00,dosage,eij,fij, eijsub, fijsub;
    dosages.clear();
    dosages.reserve(m_N);
    dosages.resize(m_N);
    std::size_t missing_cnt = 0;

    for (uint i = 0; i < N; i++) {
      p11 = lut[*bufAt]; bufAt++;
      p10 = lut[*bufAt]; bufAt++;
      p00 = 1 - p11 - p10;
      dosage = 2*p11 + p10;

      if(!missingIdxVec[i]){
        eij = dosage;
        fij = 4*p11 + p10;
        sum_eij += eij;
        sum_fij_minus_eij2 += fij - eij*eij;
        if(m_posSampleInModel[i] >= 0){
          dosages[m_posSampleInModel[i]] = 2 - dosage;
          sum_eij_sub += eij;
        }
     }else{
        if(m_posSampleInModel[i] >= 0){
          indexforMissing.push_back(m_posSampleInModel[i]);
          ++missing_cnt;
          dosages[m_posSampleInModel[i]] = -1;
        }
     }
    }

    AC = 2* ((double) (m_N - missing_cnt)) - sum_eij_sub;

    double thetaHat = sum_eij / (2* (N - missingSamplesize));
    double info = thetaHat==0 || thetaHat==1 ? 1 :
    1 - sum_fij_minus_eij2 / (2*(N-missingSamplesize)*thetaHat*(1-thetaHat));

    if(missing_cnt > 0){
      std::cout << "sample index with missing dosages for snpName " << snpName << " :";
      if(m_N == missing_cnt){
        AF = 0;
      }else{
        AF = AC/ 2/ ((double) (m_N - missing_cnt)) ;
      }

      if(!m_isDropMissingDosages){
        double imputeDosage = 2*AF;
        for (unsigned int i = 0; i < indexforMissing.size(); i++)
        {
          dosages[indexforMissing[i]] = imputeDosage;
            //std::cout << indexforMissing[i]+1 << ",";
	  AC = AC + imputeDosage;
        }
      }
    }   

    return(info);
  }


  double  BgenClass::Parse_Sparse(unsigned char * buf, size_t bufLen,  std::string & snpName, uint Nbgen,std::vector< double > & dosages, double & AC, double & AF, std::vector<unsigned int> & indexforMissing, std::vector< int > & iIndex){

    size_t destLen = bufLen;
    unsigned char * bufAt = buf;
    uint N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;

    if (N != Nbgen) {
      std::cerr << "ERROR: " << snpName << " has N = " << N << " (mismatch with header block)" << std::endl;
      exit(1);
    }
    uint K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
    if (K != 2U) {
      std::cerr << "ERROR: " << snpName << " has K = " << K << " (non-bi-allelic)" << std::endl;
      exit(1);
    }
    uint Pmin = *bufAt; bufAt++;
    if (Pmin != 2U) {
      std::cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin << " (not 2)" << std::endl;
      exit(1);
    }
    uint Pmax = *bufAt; bufAt++;
    if (Pmax != 2U) {
      std::cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax << " (not 2)" << std::endl;
      exit(1);
    }

    //deal with missing dosages
    std::vector <bool> missingIdxVec;
    missingIdxVec.clear();
    missingIdxVec.reserve(N);
    missingIdxVec.resize(N);
    int missingSamplesize = 0;

    for (uint i = 0; i < N; i++) {
      uint ploidyMiss = *bufAt; bufAt++;
      bool const missing = (ploidyMiss & 0x80) ;
      missingIdxVec[i] = missing;
        if(missing){
          missingSamplesize = missingSamplesize + 1;
        }
    }
    uint Phased = *bufAt; bufAt++;
    if (Phased != 0U) {
      std::cerr << "ERROR: " << snpName << " has Phased = " << Pmax << " (not 0)" << std::endl;
      exit(1);
    }
    uint B = *bufAt; bufAt++;
    if (B != 8U) {
      std::cerr << "ERROR: " << snpName << " has B = " << B << " (not 8)" << std::endl;
      exit(1);
    }

        // Parse
    double lut[256];
    for (int i = 0; i <= 255; i++)
      lut[i] = i/255.0;

    double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0, sum_fij_minus_eij2_sub = 0; // for INFO
    double p11,p10,p00,dosage,eij,fij, eijsub, fijsub;
    dosages.clear();
    dosages.reserve(m_N);
    dosages.resize(m_N);
    iIndex.clear();
    std::size_t missing_cnt = 0;

    for (uint i = 0; i < N; i++) {
      p11 = lut[*bufAt]; bufAt++;
      p10 = lut[*bufAt]; bufAt++;
      p00 = 1 - p11 - p10;
      dosage = 2*p11 + p10;

      if(!missingIdxVec[i]){
        eij = dosage;
        fij = 4*p11 + p10;
        sum_eij += eij;
        sum_fij_minus_eij2 += fij - eij*eij;
        if(m_posSampleInModel[i] >= 0){
          //dosages[m_posSampleInModel[i]] = 2 - dosage;
	  if(dosage < 2){
            dosages.push_back(2 - dosage);
            iIndex.push_back(m_posSampleInModel[i]+1);
          }
     	  sum_eij_sub += eij;
        }
     }else{
        if(m_posSampleInModel[i] >= 0){
          indexforMissing.push_back(m_posSampleInModel[i]);
          ++missing_cnt;
          //dosages[m_posSampleInModel[i]] = -1;
        }
     }
    }

    AC = 2* ((double) (m_N - missing_cnt)) - sum_eij_sub;

    double thetaHat = sum_eij / (2* (N - missingSamplesize));
    double info = thetaHat==0 || thetaHat==1 ? 1 :
    1 - sum_fij_minus_eij2 / (2*(N-missingSamplesize)*thetaHat*(1-thetaHat));

    if(missing_cnt > 0){
      std::cout << "sample index with missing dosages for snpName " << snpName << " :";
      if(m_N == missing_cnt){
        AF = 0;
      }else{
        AF = AC/ 2/ ((double) (m_N - missing_cnt)) ;
      }

      if(!m_isDropMissingDosages){
        double imputeDosage = 2*AF;
        for (unsigned int i = 0; i < indexforMissing.size(); i++)
        {
          //dosages[indexforMissing[i]] = imputeDosage;
            //std::cout << indexforMissing[i]+1 << ",";
	  AC = AC + imputeDosage;
	  iIndex.push_back(indexforMissing[i]+1);
	  dosages.push_back(imputeDosage);
        }
      }
    }   

    return(info);
  }




  bool BgenClass::getMultiMarker(std::string& t_markerline,
                  std::vector<uint32_t>& t_posMissingGeno,
                  std::vector<double> & t_freqMulti,
                  std::vector<double> & t_acMulti,
                  std::vector<std::string> & t_a1Multi,
                  std::vector<std::string> & t_a2Multi,
                  std::vector<std::string> & t_markerMulti,
                  std::vector<uint32_t> & t_pdMulti,
                  std::vector<uint8_t> & t_chrMarkerMulti,
                  std::vector<double> & t_dosageMulti,
                  std::vector<uint32_t> & t_iIndexMulti,
                  std::vector<uint32_t> & t_jIndexMulti
                  ){

      
   std::vector<std::string> m_markerIDsToInclude_multivar;
   std::istringstream iss(t_markerline); 
   int i = 0;
   std::string s;
   while(iss >> s){
     i = i + 1;	   
     if(i > 0){	   
       m_markerIDsToInclude_multivar.push_back(s); 
     }else{
       
     }	     
   }
   bool m_isMultiMarker;
   m_genoToTestBgen = genfile::bgen::View::create( m_bgenFileName ) ;
   genfile::bgen::IndexQuery::UniquePtr m_query = genfile::bgen::IndexQuery::create( m_bgenFileIndex ) ;
   m_query->include_rsids(m_markerIDsToInclude_multivar);   
   m_query->initialise() ;
   

   if(m_query->number_of_variants() > 0){
     m_genoToTestBgen->set_query( m_query ) ;
     m_M0 = m_genoToTestBgen->number_of_variants() ;
     m_N0 = m_genoToTestBgen->number_of_samples();
     std::cout << m_M0 << " markers are found " << std::endl;
     std::vector< genfile::byte_t > m_buffer2;
     std::vector< std::string > m_alleles;
     uint32_t m_posMarker;
     double m_freq, m_ac, m_info, m_missingRate;
     std::vector<uint32_t>  m_posMissingGeno;
     std::string m_a1,m_a2,m_marker,m_rsid,m_chrMarker;
     genfile::bgen::uint32_t  m_pd;
     std::vector<double> m_dosage;
     std::vector<int> m_iIndex;
     double maf;
     int missing_cnt;
     int cnt = 0;

     for (int i = 0; i < m_M0; i++) {	     
       m_isReadVariant = m_genoToTestBgen->read_variant(& m_marker, &m_rsid, &m_chrMarker, &m_pd, &m_alleles) ;
       if(m_isReadVariant){
         m_buffer2 = m_genoToTestBgen->read_and_uncompress_genotype_data_block();
         m_a1 = m_alleles[0];
         m_a2 = m_alleles[1];
         unsigned char * buf  = (unsigned char *) m_buffer2.data();
         uint Nbgen = m_gm_context.number_of_samples;
         m_ac = 0;
         m_freq = 0;
         m_info = Parse_Sparse(buf, m_buffer2.size(), m_marker, Nbgen, m_dosage, m_ac, m_freq, m_posMissingGeno, m_iIndex);
         missing_cnt = m_posMissingGeno.size();

         maf = m_freq;
         if(m_freq > 0.5){maf = 1-maf;}
	 std::vector<int> jIndexforOneMarker(m_iIndex.size(), i+1); 

         if(maf >= m_minMAF && maf <= m_maxMAF && m_info >= m_info){
           t_dosageMulti.insert(std::end(t_dosageMulti), std::begin(m_dosage), std::end(m_dosage));
           t_iIndexMulti.insert(std::end(t_iIndexMulti), std::begin(m_iIndex), std::end(m_iIndex));
           t_jIndexMulti.insert(std::end(t_jIndexMulti), std::begin(jIndexforOneMarker), std::end(jIndexforOneMarker));

           cnt = cnt + 1;
           t_markerMulti.push_back(m_marker);
           t_freqMulti.push_back(m_freq);
           t_acMulti.push_back(m_ac);
           t_pdMulti.push_back(m_pd);
         }
       }else{
         m_isMultiMarker=false;
       }
     }  

   }else{
     std::cout << "No queried variant is found in the bgen file! All variants bgen file will be analyzed" << std::endl;
     m_isMultiMarker=false;
   }

    return(m_isMultiMarker);
  }

}
