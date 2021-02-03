#ifndef BGEN_HPP
#define BGEN_HPP

#include "../thirdParty/bgen/genfile/include/genfile/bgen/bgen.hpp"
#include "../thirdParty/bgen/genfile/include/genfile/bgen/View.hpp"
#include "../thirdParty/bgen/genfile/include/genfile/bgen/IndexQuery.hpp"

#include <Rcpp.h>

namespace BGEN {

class BgenClass{
private:

  // vcf file and the index file (.tbi)
  std::string m_bgenFileName, m_bgenFileIndex;
  std::string m_markerIDsToExclude;
  std::string m_markerIDsToInclude; 

  std::vector<uint32_t> m_posSampleInBgen;
  std::vector<uint32_t> m_posSampleInModel;
  bool m_isDropMissingDosages;

  std::auto_ptr< std::istream > m_gm_stream;
  genfile::bgen::View::UniquePtr m_genoToTestBgen;
  uint32_t  m_gm_offset ;
  genfile::bgen::Context m_gm_context ;
  bool gm_have_sample_ids ;
  int gmtest_samplesize;  


  bool m_isQuery;
  //field in bgen to test. GT or DS 
  uint32_t m_M0, m_M;
  uint32_t m_N0, m_N;
  std::vector<std::string> m_MarkerInBgen;     // Variant identifier  
  std::vector<std::string> m_SampleInBgen;

  bool m_isBgenOpen; 
  bool m_isReadVariant;

  double m_minMAF, m_maxMAF;
  double m_minInfo;
  
public:

  BgenClass(std::string t_bgenFileName,
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
	    std::vector<std::string> t_SampleInModel);



  void setBgenObj(const std::string t_bgenFileName,
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
		 std::vector< double > & t_dosages);



  void setPosSampleInBgen(std::vector<std::string> & t_SampleInModel);  

  // get dosages/genotypes of one marker
  bool getOneMarker(uint32_t t_posMarker,
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
			 std::vector<double> & t_dosage);

  double Parse(unsigned char * buf, size_t bufLen,  std::string & snpName, uint Nbgen,std::vector< double > & dosages, double & AC, double & AF, std::vector<unsigned int> & indexforMissing);

  double Parse_Sparse(unsigned char * buf, size_t bufLen,  std::string & snpName, uint Nbgen,std::vector< double > & dosages, double & AC, double & AF, std::vector<unsigned int> & indexforMissing, std::vector< int > & iIndex);


  bool getMultiMarker(std::string& t_markerline,
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
                  );


  //uint32_t getN0(){return m_N0;}
  //uint32_t getN(){return m_N;}
  //uint32_t getM0(){return m_M0;}
  //uint32_t getM(){return m_M;}
};

}

#endif
