#ifndef BGEN_HPP
#define BGEN_HPP

#include <Rcpp.h>

namespace BGEN {

class BgenClass{
private:

  // vcf file and the index file (.tbi)
  std::string m_bgenFileName, m_bgenFileIndex;

  //variables to reading in bgen
  FILE *m_fin;
  std::vector <unsigned char> m_buf;
  std::vector <unsigned char> m_zBuf;
  unsigned int m_zBufLens;
  unsigned int m_bufLens; 


  std::vector<uint32_t> m_posSampleInBgen;
  std::vector<uint32_t> m_posSampleInModel;
  bool m_isDropMissingDosagesInBgen;
  bool m_isQuery;
  bool m_isSparseDosagesInBgen;
  std::vector< int > m_markerIndicesToInclude;

  uint32_t m_M0, m_M;
  uint32_t m_N0, m_N;
  std::vector<std::string> m_MarkerInBgen;     // Variant identifier  
  std::vector<std::string> m_SampleInBgen;

  double m_minMAF, m_maxMAF;
  double m_minInfo;
  
public:

  BgenClass(std::string t_bgenFileName,
            std::string t_bgenFileIndex,
            std::vector<std::string> t_SampleInBgen,
            std::vector<std::string> t_SampleInModel,
            bool t_isSparseDosageInBgen,
            bool t_isDropmissingdosagesInBgen);
	  
  void setBgenObj(const std::string t_bgenFileName,
                 const std::string t_bgenFileIndex,
                 std::vector<std::string> & t_SampleInBgen);

  void setPosSampleInBgen(std::vector<std::string> & t_SampleInModel);

  void Parse2(unsigned char *buf, unsigned int bufLen, const unsigned char *zBuf, unsigned int zBufLen,std::string & snpName,std::vector< double > & dosages, double & AC, double & AF, std::vector<int> & indexforMissing, double & info, std::vector<unsigned int> & iIndex);

  Rcpp::List getOneMarker(int t_fileStartPos);
  // get dosages/genotypes of one marker
  /*
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
  */

  void setIsSparseDosageInBgen (bool t_isSparseDosageInBgen);
  void setMarkerIndicesToIncludeInBgen (std::vector< int > & t_markerIndicesToInclude);
  void setIsDropMissingDosagesInBgen (bool t_isDropmissingdosagesInBgen);

  uint32_t getN0(){return m_N0;};
  uint32_t getN(){return m_N;};
  uint32_t getM0(){return m_M0;};
  uint32_t getM(){return m_M;};
};

}

#endif
