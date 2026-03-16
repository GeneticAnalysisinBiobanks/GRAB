// PLINK4.cpp -- PlinkReader method implementations

#include "PLINK4.h"

#include <boost/algorithm/string.hpp>

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

namespace PLINK4 {

PlinkReader::PlinkReader(const std::string &bimFile,
                         const std::string &famFile,
                         const std::string &bedFile,
                         const std::vector<std::string> &sampleInModel,
                         const std::string &alleleOrder)
  : m_alleleOrder(alleleOrder) {
  setPlinkObject(bimFile, famFile, bedFile);
  setPosSampleInPlink(sampleInModel);
}

void PlinkReader::setPlinkObject(const std::string &bimFile,
                                 const std::string &famFile,
                                 const std::string &bedFile) {
  m_bimFile = bimFile;
  m_famFile = famFile;
  m_bedFile = bedFile;

  readBimFile();
  readFamFile();

  m_bedStream.open(m_bedFile.c_str(), std::ios::binary);
  if (!m_bedStream.is_open()) {
    throw std::runtime_error("Cannot open PLINK bed file: " + m_bedFile);
  }

  m_bedStream.seekg(2);
  char magicNumber3 = 0;
  m_bedStream.read(&magicNumber3, 1);
  if (!m_bedStream.good()) {
    throw std::runtime_error("Failed to read PLINK bed header: " + m_bedFile);
  }
  if (magicNumber3 != 1) {
    throw std::runtime_error("The third magic number of the plink bed file is not 00000001. Please use SNP-major plink (plink version >= 1.9) files.");
  }

  m_hasSequentialCursor = false;
}

void PlinkReader::readBimFile() {
  std::ifstream bim(m_bimFile);
  if (!bim.is_open()) {
    throw std::runtime_error("Cannot open PLINK bim file: " + m_bimFile);
  }

  m_m0 = 0;
  m_chr.clear();
  m_markerInPlink.clear();
  m_gd.clear();
  m_pd.clear();
  m_alt.clear();
  m_ref.clear();

  std::string line;
  while (getline(bim, line)) {
    ++m_m0;
    std::vector<std::string> lineElements;
    boost::split(lineElements, line, boost::is_any_of("\t "));
    boost::replace_all(lineElements.back(), "\r", "");

    m_chr.push_back(lineElements[0]);
    m_markerInPlink.push_back(lineElements[1]);
    m_gd.push_back(std::stof(lineElements[2]));
    m_pd.push_back(std::stoi(lineElements[3]));
    std::transform(lineElements[4].begin(), lineElements[4].end(), lineElements[4].begin(), toupper);
    std::transform(lineElements[5].begin(), lineElements[5].end(), lineElements[5].begin(), toupper);
    m_alt.push_back(lineElements[4]);
    m_ref.push_back(lineElements[5]);
  }

  m_m = m_m0;
}

void PlinkReader::readFamFile() {
  std::ifstream fam(m_famFile);
  if (!fam.is_open()) {
    throw std::runtime_error("Cannot open PLINK fam file: " + m_famFile);
  }

  m_n0 = 0;
  m_sampleInPlink.clear();

  std::string line;
  while (getline(fam, line)) {
    ++m_n0;
    std::vector<std::string> lineElements;
    boost::split(lineElements, line, boost::is_any_of("\t "));
    m_sampleInPlink.push_back(lineElements[1]);
  }

  m_n = m_n0;
  m_numBytesOfEachMarker0 = (m_n0 + 3) / 4;
  m_oneMarkerG4.resize(m_numBytesOfEachMarker0);
}

void PlinkReader::setPosSampleInPlink(const std::vector<std::string> &sampleInModel) {
  m_n = static_cast<uint32_t>(sampleInModel.size());
  m_numBytesOfEachMarker = (m_n + 3) / 4;

  std::unordered_map<std::string, uint32_t> plinkPos;
  plinkPos.reserve(m_sampleInPlink.size());
  for (uint32_t i = 0; i < m_sampleInPlink.size(); ++i) {
    const std::string &id = m_sampleInPlink[i];
    if (plinkPos.find(id) == plinkPos.end()) {
      plinkPos.emplace(id, i);
    }
  }

  m_posSampleInPlink.resize(m_n);
  for (uint32_t i = 0; i < m_n; ++i) {
    auto it = plinkPos.find(sampleInModel[i]);
    if (it == plinkPos.end()) {
      throw std::runtime_error("At least one subject requested is not in Plink file.");
    }
    m_posSampleInPlink[i] = it->second;
  }
}

void PlinkReader::beginSequentialBlock(uint64_t firstMarkerIndex) {
  if (firstMarkerIndex >= m_m0) {
    throw std::runtime_error("PLINK marker index out of range in beginSequentialBlock.");
  }

  const uint64_t posSeek = 3 + m_numBytesOfEachMarker0 * firstMarkerIndex;
  m_bedStream.clear();
  m_bedStream.seekg(posSeek);
  if (!m_bedStream.good()) {
    throw std::runtime_error("Failed to seek PLINK bed file for sequential block.");
  }
  m_hasSequentialCursor = true;
  m_nextMarkerIndex = firstMarkerIndex;
  m_hasBufferedBlock = false;
}

void PlinkReader::resetSequentialCursor() {
  m_hasSequentialCursor = false;
  m_hasBufferedBlock = false;
}

void PlinkReader::loadSequentialMarkerBlock(uint64_t startMarkerIndex) {
  if (startMarkerIndex >= m_m0) {
    throw std::runtime_error("PLINK marker index out of range when loading marker block.");
  }

  const uint64_t nRemain = static_cast<uint64_t>(m_m0) - startMarkerIndex;
  const uint64_t nMarkersToRead = std::min(m_blockMarkerCapacity, nRemain);
  const uint64_t nBytesToRead = nMarkersToRead * m_numBytesOfEachMarker0;

  if (m_markerBlockBytes.size() != nBytesToRead) {
    m_markerBlockBytes.resize(static_cast<size_t>(nBytesToRead));
  }

  if (!m_hasSequentialCursor || m_nextMarkerIndex != startMarkerIndex) {
    const uint64_t posSeek = 3 + m_numBytesOfEachMarker0 * startMarkerIndex;
    m_bedStream.clear();
    m_bedStream.seekg(posSeek);
    if (!m_bedStream.good()) {
      throw std::runtime_error("Failed to seek PLINK bed file for block read.");
    }
  }

  m_bedStream.read(reinterpret_cast<char *>(&m_markerBlockBytes[0]), static_cast<std::streamsize>(nBytesToRead));
  if (!m_bedStream.good()) {
    throw std::runtime_error("Failed to read PLINK marker block.");
  }

  m_blockStartMarker = startMarkerIndex;
  m_blockEndMarker = startMarkerIndex + nMarkersToRead;
  m_hasBufferedBlock = true;

  m_hasSequentialCursor = true;
  m_nextMarkerIndex = m_blockEndMarker;
}

void PlinkReader::copyBufferedMarkerToCurrent(uint64_t markerIndex) {
  if (!m_hasBufferedBlock || markerIndex < m_blockStartMarker || markerIndex >= m_blockEndMarker) {
    throw std::runtime_error("Marker index is outside of buffered PLINK marker block.");
  }

  const uint64_t markerOffset = markerIndex - m_blockStartMarker;
  const uint64_t byteOffset = markerOffset * m_numBytesOfEachMarker0;

  std::copy(
    m_markerBlockBytes.begin() + static_cast<std::ptrdiff_t>(byteOffset),
    m_markerBlockBytes.begin() + static_cast<std::ptrdiff_t>(byteOffset + m_numBytesOfEachMarker0),
    m_oneMarkerG4.begin()
  );
}

void PlinkReader::readCurrentMarkerBytes(uint64_t gIndex) {
  if (gIndex >= m_m0) {
    throw std::runtime_error("PLINK marker index out of range.");
  }

  if (m_hasBufferedBlock && gIndex >= m_blockStartMarker && gIndex < m_blockEndMarker) {
    copyBufferedMarkerToCurrent(gIndex);
    return;
  }

  if (m_hasSequentialCursor && gIndex == m_nextMarkerIndex) {
    loadSequentialMarkerBlock(gIndex);
    copyBufferedMarkerToCurrent(gIndex);
    return;
  }

  if (!m_hasSequentialCursor || gIndex != m_nextMarkerIndex) {
    const uint64_t posSeek = 3 + m_numBytesOfEachMarker0 * gIndex;
    m_bedStream.clear();
    m_bedStream.seekg(posSeek);
    if (!m_bedStream.good()) {
      throw std::runtime_error("Failed to seek PLINK bed file.");
    }
  }

  m_bedStream.read(reinterpret_cast<char *>(&m_oneMarkerG4[0]), m_numBytesOfEachMarker0);
  if (!m_bedStream.good()) {
    throw std::runtime_error("Failed to read PLINK marker bytes.");
  }

  m_hasBufferedBlock = false;
  m_hasSequentialCursor = true;
  m_nextMarkerIndex = gIndex + 1;
}

arma::vec PlinkReader::getOneMarker(uint64_t gIndex,
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
                                    bool isTrueGenotype) {
  int sum = 0;
  int numMissing = 0;
  std::vector<double> oneMarkerG1;

  if (!isTrueGenotype) {
    if (isOutputIndexForMissing) {
      throw std::runtime_error("Check PlinkReader::getOneMarker, if isTrueGenotype = FALSE, then isOutputIndexForMissing should be FALSE.");
    }
    if (isOnlyOutputNonZero) {
      throw std::runtime_error("Check PlinkReader::getOneMarker, if isTrueGenotype = FALSE, then isOnlyOutputNonZero should be FALSE.");
    }
  }

  if (!isOnlyOutputNonZero) {
    oneMarkerG1.resize(m_n);
  }

  readCurrentMarkerBytes(gIndex);

  indexForMissing.clear();
  indexForNonZero.clear();

  marker = m_markerInPlink[gIndex];
  pd = m_pd[gIndex];
  chr = m_chr[gIndex];

  const std::map<int8_t, int8_t> *genoMaps = nullptr;
  if (m_alleleOrder == "alt-first") {
    ref = m_ref[gIndex];
    alt = m_alt[gIndex];
    genoMaps = &m_genoMapsAltFirst;
  } else if (m_alleleOrder == "ref-first") {
    ref = m_alt[gIndex];
    alt = m_ref[gIndex];
    genoMaps = &m_genoMapsRefFirst;
  } else {
    throw std::runtime_error("Unsupported PLINK allele order: " + m_alleleOrder);
  }

  for (uint32_t i = 0; i < m_n; ++i) {
    const uint32_t ind = m_posSampleInPlink[i];
    unsigned char bufferG4 = m_oneMarkerG4[ind / 4];
    int bufferG1 = 0;
    getGenotype(&bufferG4, ind % 4, bufferG1);

    switch (bufferG1) {
      case HOM_REF: break;
      case HET: sum += 1; break;
      case HOM_ALT: sum += 2; break;
      case MISSING:
        ++numMissing;
        if (isOutputIndexForMissing) {
          indexForMissing.push_back(i);
        }
        break;
    }

    if (isTrueGenotype) {
      bufferG1 = genoMaps->at(bufferG1);
    }

    if (isOnlyOutputNonZero) {
      if (bufferG1 > 0) {
        indexForNonZero.push_back(i);
        oneMarkerG1.push_back(bufferG1);
      }
    } else {
      oneMarkerG1[i] = bufferG1;
    }
  }

  const int count = static_cast<int>(m_n) - numMissing;
  missingRate = static_cast<double>(numMissing) / static_cast<double>(m_n);
  imputeInfo = 1.0;
  altCounts = static_cast<double>(sum);
  altFreq = altCounts / static_cast<double>(count) / 2.0;

  if (m_alleleOrder == "ref-first") {
    altFreq = 1.0 - altFreq;
    altCounts = 2.0 * static_cast<double>(count) * altFreq;
  }

  return oneMarkerG1;
}

}
