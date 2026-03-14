#ifndef PLINK4_H
#define PLINK4_H

#include <RcppArmadillo.h>

#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <vector>

namespace PLINK4 {

class PlinkReader {
private:
  std::string m_alleleOrder;

  uint32_t m_m0 = 0;
  uint32_t m_m = 0;
  std::vector<std::string> m_chr;
  std::vector<std::string> m_markerInPlink;
  std::vector<float> m_gd;
  std::vector<uint32_t> m_pd;
  std::vector<std::string> m_alt;
  std::vector<std::string> m_ref;

  std::vector<std::string> m_sampleInPlink;
  uint32_t m_n0 = 0;
  uint32_t m_n = 0;
  unsigned long long m_numBytesOfEachMarker0 = 0;
  unsigned long long m_numBytesOfEachMarker = 0;

  std::ifstream m_bedStream;
  std::string m_bimFile;
  std::string m_famFile;
  std::string m_bedFile;
  std::vector<uint32_t> m_posSampleInPlink;

  static constexpr unsigned char HOM_REF = 0x3;
  static constexpr unsigned char HET = 0x2;
  static constexpr unsigned char HOM_ALT = 0x0;
  static constexpr unsigned char MISSING = 0x1;

  std::map<int8_t, int8_t> m_genoMapsAltFirst = {{3, 0}, {2, 1}, {0, 2}, {1, -1}};
  std::map<int8_t, int8_t> m_genoMapsRefFirst = {{3, 2}, {2, 1}, {0, 0}, {1, -1}};
  std::vector<unsigned char> m_oneMarkerG4;

  // Worker-local input buffer for multiple consecutive markers.
  uint64_t m_blockMarkerCapacity = 64;
  std::vector<unsigned char> m_markerBlockBytes;
  uint64_t m_blockStartMarker = 0;
  uint64_t m_blockEndMarker = 0;
  bool m_hasBufferedBlock = false;

  bool m_hasSequentialCursor = false;
  uint64_t m_nextMarkerIndex = 0;

  void readBimFile();
  void readFamFile();
  void setPlinkObject(const std::string &bimFile,
                      const std::string &famFile,
                      const std::string &bedFile);
  void setPosSampleInPlink(const std::vector<std::string> &sampleInModel);
  void loadSequentialMarkerBlock(uint64_t startMarkerIndex);
  void copyBufferedMarkerToCurrent(uint64_t markerIndex);
  void readCurrentMarkerBytes(uint64_t gIndex);

  static void getGenotype(unsigned char *c, int pos, int &geno) {
    geno = ((*c) >> (pos << 1)) & 0x3;
  }

public:
  PlinkReader(const std::string &bimFile,
              const std::string &famFile,
              const std::string &bedFile,
              const std::vector<std::string> &sampleInModel,
              const std::string &alleleOrder);

  void beginSequentialBlock(uint64_t firstMarkerIndex);
  void resetSequentialCursor();

  arma::vec getOneMarker(uint64_t gIndex,
                         std::string &ref,
                         std::string &alt,
                         std::string &marker,
                         uint32_t &pd,
                         std::string &chr,
                         double &altFreq,
                         double &altCounts,
                         double &missingRate,
                         double &imputeInfo,
                         bool isOutputIndexForMissing,
                         std::vector<uint32_t> &indexForMissing,
                         bool isOnlyOutputNonZero,
                         std::vector<uint32_t> &indexForNonZero,
                         bool isTrueGenotype);
};

}

#endif