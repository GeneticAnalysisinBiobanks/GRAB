// Main4.cpp - chunk-level multi-threaded marker analysis and file writing
//
// Runtime workflow (aligned with GRAB.Marker single-thread flow, but chunk-parallel):
// 1) prepareAndRunMarkerInCPP4
//    - validate model/genotype mode and choose method
//    - build method object from null model (set*objInCPP)
//    - build PLINK reader config and marker list
//    - apply include/exclude marker filters
//    - split markers into chromosome-aware chunks
//    - delegates to mainMarkerChunksCore
// 2) unified cleanup and error propagation

#include <RcppArmadillo.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include <limits>
#include <cmath>
#include <exception>
#include <cstdio>

#include <boost/math/distributions/chi_squared.hpp>

#include "PLINK4.h"
#include "POLMM.h"
#include "SPACox.h"
#include "SPAmix.h"
#include "SPAGRM.h"
#include "SAGELD.h"
#include "WtCoxG.h"
#include "SPAsqr.h"
#include "LEAF.h"
#include "SPAmixPlus.h"
#include "UTIL.h"
#include "Main.h"

namespace {

struct Main4ExtraParams;

void mainMarkerChunksCore(
  const std::string &method,
  const std::vector<std::vector<uint64_t>> &chunkMarkers,
  const std::string &outputFile,
  unsigned int nThreads,
  const std::string &imputeMethod,
  double missingCutoff,
  double minMafMarker,
  double minMacMarker,
  const Main4ExtraParams &extraParams
);

// ===== Types used by native PLINK preparation =====

struct RangeFilter {
  std::string chrom;
  uint32_t start;
  uint32_t end;
};

struct PlinkMarkerInfo {
  std::string chrom;
  uint32_t pos;
  std::string id;
  std::string ref;
  std::string alt;
  uint64_t genoIndex;
};

struct ReaderConfig {
  std::string genoType;
  std::string bimFile;
  std::string famFile;
  std::string bedFile;
  std::vector<std::string> sampleInModel;
  std::string alleleOrder;
};

struct MarkerFilterConfig {
  std::string idsToIncludeFile;
  std::string rangesToIncludeFile;
  std::string idsToExcludeFile;
  std::string rangesToExcludeFile;
};

struct WtRow {
  double AF_ref;
  double AN_ref;
  double TPR;
  double sigma2;
  double pvalue_bat;
  double w_ext;
  double var_ratio_w0;
  double var_ratio_int;
  double var_ratio_ext;
};

struct LeafRow {
  double AF_ref;
  double AN_ref;
  double TPR;
  double sigma2;
  double pvalue_bat;
  double w_ext;
  double var_ratio_w0;
  double var_ratio_int;
  double var_ratio_ext;
};

// Native configuration passed from the R boundary into the pure C++ marker core.
// After prepareAndRunMarkerInCPP4 builds this struct, downstream code no longer
// needs to inspect Rcpp::List payloads to decide reader or method-specific state.
struct Main4ExtraParams {
  ReaderConfig reader;
  std::string sageldMethod = "SAGELD";
  std::vector<double> spasqrTaus;
  std::unordered_map<uint64_t, WtRow> wtMap;
  std::vector<std::unordered_map<uint64_t, LeafRow>> leafMaps;
  int leafNcluster = 0;
};

std::string toLowerCopy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

std::vector<std::string> splitWhitespace(const std::string &line) {
  std::istringstream iss(line);
  std::vector<std::string> tokens;
  std::string token;
  while (iss >> token) {
    tokens.push_back(token);
  }
  return tokens;
}

std::string getMethodFromNullModelClass(const std::string &nullModelClass) {
  if (nullModelClass == "POLMM_NULL_Model") return "POLMM";
  if (nullModelClass == "SPACox_NULL_Model") return "SPACox";
  if (nullModelClass == "SPAmix_NULL_Model") return "SPAmix";
  if (nullModelClass == "SPAmixPlus_NULL_Model") return "SPAmixPlus";
  if (nullModelClass == "SPAGRM_NULL_Model") return "SPAGRM";
  if (nullModelClass == "SAGELD_NULL_Model") return "SAGELD";
  if (nullModelClass == "WtCoxG_NULL_Model") return "WtCoxG";
  if (nullModelClass == "SPAsqr_NULL_Model") return "SPAsqr";
  if (nullModelClass == "LEAF_NULL_Model") return "LEAF";
  throw std::runtime_error("Unsupported null model class in Marker4: " + nullModelClass);
}

bool fileExists(const std::string &path) {
  if (path.empty()) return false;
  std::ifstream in(path.c_str(), std::ios::binary);
  return in.good();
}

std::string replaceExtension(const std::string &path, const std::string &newExtension) {
  const size_t pos = path.find_last_of('.');
  if (pos == std::string::npos) {
    return path + newExtension;
  }
  return path.substr(0, pos) + newExtension;
}

std::vector<std::string> readSingleColumnFile(const std::string &path) {
  std::ifstream in(path.c_str());
  if (!in.is_open()) {
    throw std::runtime_error("Cannot open marker ID filter file: " + path);
  }

  std::vector<std::string> values;
  std::string line;
  while (std::getline(in, line)) {
    std::vector<std::string> tokens = splitWhitespace(line);
    if (!tokens.empty()) {
      values.push_back(tokens[0]);
    }
  }
  return values;
}

std::vector<RangeFilter> readRangeFile(const std::string &path) {
  std::ifstream in(path.c_str());
  if (!in.is_open()) {
    throw std::runtime_error("Cannot open range filter file: " + path);
  }

  std::vector<RangeFilter> ranges;
  std::string line;
  while (std::getline(in, line)) {
    std::vector<std::string> tokens = splitWhitespace(line);
    if (tokens.empty()) {
      continue;
    }
    if (tokens.size() != 3) {
      throw std::runtime_error("Range filter file should include exactly three whitespace-separated columns: " + path);
    }
    ranges.push_back(RangeFilter{
      tokens[0],
      static_cast<uint32_t>(std::stoul(tokens[1])),
      static_cast<uint32_t>(std::stoul(tokens[2]))
    });
  }
  return ranges;
}

bool markerInRanges(const PlinkMarkerInfo &marker, const std::vector<RangeFilter> &ranges) {
  for (const RangeFilter &range : ranges) {
    if (marker.chrom == range.chrom && marker.pos >= range.start && marker.pos <= range.end) {
      return true;
    }
  }
  return false;
}

std::vector<PlinkMarkerInfo> readPlinkMarkerInfo(const std::string &bimFile, const std::string &alleleOrder) {
  std::ifstream in(bimFile.c_str());
  if (!in.is_open()) {
    throw std::runtime_error("Cannot open .bim file: " + bimFile);
  }

  std::vector<PlinkMarkerInfo> markers;
  std::string line;
  uint64_t markerIndex = 0;
  while (std::getline(in, line)) {
    std::vector<std::string> tokens = splitWhitespace(line);
    if (tokens.empty()) {
      continue;
    }
    if (tokens.size() != 6) {
      throw std::runtime_error("PLINK .bim file should include 6 whitespace-separated columns: " + bimFile);
    }

    std::string ref = tokens[5];
    std::string alt = tokens[4];
    if (alleleOrder == "ref-first") {
      ref = tokens[4];
      alt = tokens[5];
    }
    std::transform(ref.begin(), ref.end(), ref.begin(), ::toupper);
    std::transform(alt.begin(), alt.end(), alt.begin(), ::toupper);

    markers.push_back(PlinkMarkerInfo{
      tokens[0],
      static_cast<uint32_t>(std::stoul(tokens[3])),
      tokens[1],
      ref,
      alt,
      markerIndex
    });
    ++markerIndex;
  }
  return markers;
}

std::vector<std::string> readPlinkSampleIds(const std::string &famFile) {
  std::ifstream in(famFile.c_str());
  if (!in.is_open()) {
    throw std::runtime_error("Cannot open .fam file: " + famFile);
  }

  std::vector<std::string> sampleIds;
  std::string line;
  while (std::getline(in, line)) {
    std::vector<std::string> tokens = splitWhitespace(line);
    if (tokens.empty()) {
      continue;
    }
    if (tokens.size() != 6) {
      throw std::runtime_error("PLINK .fam file should include 6 whitespace-separated columns: " + famFile);
    }
    sampleIds.push_back(tokens[1]);
  }
  return sampleIds;
}

void validateRequestedSamples(const std::vector<std::string> &requested,
                              const std::vector<std::string> &available) {
  std::unordered_set<std::string> availableSet(available.begin(), available.end());
  for (const std::string &sampleId : requested) {
    if (availableSet.find(sampleId) == availableSet.end()) {
      throw std::runtime_error("At least one subject requested is not in PLINK file.");
    }
  }
}

std::vector<PlinkMarkerInfo> applyPlinkMarkerFilters(const std::vector<PlinkMarkerInfo> &markers,
                                                     const MarkerFilterConfig &filterConfig) {
  // Include/exclude semantics intentionally mirror setGenoInput() in R/Geno.R.
  std::unordered_set<std::string> includeIds;
  std::unordered_set<std::string> excludeIds;
  std::vector<RangeFilter> includeRanges;
  std::vector<RangeFilter> excludeRanges;

  if (!filterConfig.idsToIncludeFile.empty()) {
    std::vector<std::string> ids = readSingleColumnFile(filterConfig.idsToIncludeFile);
    includeIds.insert(ids.begin(), ids.end());
  }
  if (!filterConfig.rangesToIncludeFile.empty()) {
    includeRanges = readRangeFile(filterConfig.rangesToIncludeFile);
  }
  if (!filterConfig.idsToExcludeFile.empty()) {
    std::vector<std::string> ids = readSingleColumnFile(filterConfig.idsToExcludeFile);
    excludeIds.insert(ids.begin(), ids.end());
  }
  if (!filterConfig.rangesToExcludeFile.empty()) {
    excludeRanges = readRangeFile(filterConfig.rangesToExcludeFile);
  }

  const bool anyInclude = !includeIds.empty() || !includeRanges.empty();
  const bool anyExclude = !excludeIds.empty() || !excludeRanges.empty();

  std::vector<PlinkMarkerInfo> filtered;
  filtered.reserve(markers.size());
  for (const PlinkMarkerInfo &marker : markers) {
    bool includeMarker = !anyInclude;
    if (anyInclude) {
      includeMarker = (includeIds.find(marker.id) != includeIds.end()) || markerInRanges(marker, includeRanges);
    }
    if (!includeMarker) {
      continue;
    }

    bool excludeMarker = false;
    if (anyExclude) {
      excludeMarker = (excludeIds.find(marker.id) != excludeIds.end()) || markerInRanges(marker, excludeRanges);
    }
    if (!excludeMarker) {
      filtered.push_back(marker);
    }
  }
  return filtered;
}

std::vector<std::vector<uint64_t>> buildPlinkChunkIndexList(const std::vector<PlinkMarkerInfo> &markers,
                                                            const int nMarkersEachChunk) {
  std::vector<std::vector<uint64_t>> chunks;
  size_t start = 0;
  while (start < markers.size()) {
    const std::string chrom = markers[start].chrom;
    size_t chromEnd = start;
    while (chromEnd < markers.size() && markers[chromEnd].chrom == chrom) {
      ++chromEnd;
    }
    for (size_t chunkStart = start; chunkStart < chromEnd; chunkStart += static_cast<size_t>(nMarkersEachChunk)) {
      const size_t chunkEnd = std::min(chunkStart + static_cast<size_t>(nMarkersEachChunk), chromEnd);
      std::vector<uint64_t> chunk;
      chunk.reserve(chunkEnd - chunkStart);
      for (size_t i = chunkStart; i < chunkEnd; ++i) {
        chunk.push_back(markers[i].genoIndex);
      }
      chunks.push_back(std::move(chunk));
    }
    start = chromEnd;
  }
  return chunks;
}

ReaderConfig buildPlinkReaderConfig(const std::string &genoFile,
                                    const std::vector<std::string> &genoFileIndex,
                                    const std::vector<std::string> &sampleInModel,
                                    const std::string &alleleOrder) {
  // This native setup replaces the old setGenoInput() callback path for PLINK.
  const size_t dotPos = genoFile.find_last_of('.');
  const std::string ext = dotPos == std::string::npos ? std::string() : toLowerCopy(genoFile.substr(dotPos + 1));
  if (ext != "bed") {
    throw std::runtime_error("Pure C++ GRAB.Marker4 currently supports PLINK .bed input only.");
  }

  ReaderConfig reader;
  reader.genoType = "PLINK";
  reader.bedFile = genoFile;
  reader.sampleInModel = sampleInModel;
  reader.alleleOrder = alleleOrder;

  if (!genoFileIndex.empty()) {
    if (genoFileIndex.size() < 2) {
      throw std::runtime_error("For PLINK input, GenoFileIndex should contain .bim and .fam paths.");
    }
    reader.bimFile = genoFileIndex[0];
    reader.famFile = genoFileIndex[1];
  } else {
    reader.bimFile = replaceExtension(genoFile, ".bim");
    reader.famFile = replaceExtension(genoFile, ".fam");
  }

  if (!fileExists(reader.bedFile) || !fileExists(reader.bimFile) || !fileExists(reader.famFile)) {
    throw std::runtime_error("One or more PLINK files are missing.");
  }

  return reader;
}

void setMarkerObjectForMethod(const std::string &nullModelClass,
                              const Rcpp::List &objNull,
                              const Rcpp::List &control) {
  // Match the single-thread dispatcher pattern: initialize exactly one method object.
  if (nullModelClass == "POLMM_NULL_Model") {
    Rcpp::List objCHR = Rcpp::as<Rcpp::List>(Rcpp::as<Rcpp::List>(objNull["LOCOList"])["LOCO=F"]);
    arma::uvec group = Rcpp::as<arma::uvec>(objNull["yVec"]);
    std::unordered_set<int> uniqueGroups;
    for (arma::uword value : group) {
      uniqueGroups.insert(static_cast<int>(value));
    }

    setPOLMMobjInCPP(
      Rcpp::as<arma::mat>(objCHR["muMat"]),
      Rcpp::as<arma::mat>(objCHR["iRMat"]),
      Rcpp::as<arma::mat>(objNull["Cova"]),
      Rcpp::as<arma::uvec>(objNull["yVec"]),
      Rcpp::as<double>(objNull["tau"]),
      false,
      0.001,
      100,
      Rcpp::as<double>(objCHR["VarRatio"]),
      Rcpp::as<double>(control["SPA_Cutoff"]),
      false,
      group,
      Rcpp::as<bool>(control["ifOutGroup"]),
      static_cast<int>(uniqueGroups.size())
    );
    return;
  }

  if (nullModelClass == "SPACox_NULL_Model") {
    setSPACoxobjInCPP(
      Rcpp::as<arma::mat>(objNull["cumul"]),
      Rcpp::as<arma::vec>(objNull["mresid"]),
      Rcpp::as<arma::mat>(objNull["X.invXX"]),
      Rcpp::as<arma::mat>(objNull["tX"]),
      Rf_length(objNull["mresid"]),
      Rcpp::as<double>(control["pVal_covaAdj_Cutoff"]),
      Rcpp::as<double>(control["SPA_Cutoff"])
    );
    return;
  }

  if (nullModelClass == "SPAmix_NULL_Model") {
    setSPAmixobjInCPP(
      Rcpp::as<arma::mat>(objNull["resid"]),
      Rcpp::as<arma::mat>(objNull["PCs"]),
      Rcpp::as<int>(objNull["N"]),
      Rcpp::as<double>(control["SPA_Cutoff"]),
      Rcpp::as<Rcpp::List>(objNull["outLierList"])
    );
    return;
  }

  if (nullModelClass == "SPAmixPlus_NULL_Model") {
    setSPAmixPlusobjInCPP(
      Rcpp::as<arma::mat>(objNull["resid"]),
      Rcpp::as<arma::mat>(objNull["PCs"]),
      Rcpp::as<int>(objNull["N"]),
      Rcpp::as<double>(control["SPA_Cutoff"]),
      Rcpp::as<Rcpp::List>(objNull["outLierList"]),
      Rcpp::as<Rcpp::DataFrame>(objNull["sparseGRM"]),
      Rcpp::as<std::string>(control["afFilePath"]),
      Rcpp::as<std::string>(control["afFilePrecision"])
    );
    return;
  }

  if (nullModelClass == "SPAGRM_NULL_Model") {
    setSPAGRMobjInCPP(
      Rcpp::as<arma::vec>(objNull["Resid"]),
      Rcpp::as<arma::vec>(objNull["Resid.unrelated.outliers"]),
      Rcpp::as<double>(objNull["sum_R_nonOutlier"]),
      Rcpp::as<double>(objNull["R_GRM_R_nonOutlier"]),
      Rcpp::as<double>(objNull["R_GRM_R_TwoSubjOutlier"]),
      Rcpp::as<double>(objNull["R_GRM_R"]),
      Rcpp::as<arma::vec>(objNull["MAF_interval"]),
      Rcpp::as<Rcpp::List>(objNull["TwoSubj_list"]),
      Rcpp::as<Rcpp::List>(objNull["ThreeSubj_list"]),
      Rcpp::as<double>(control["SPA_Cutoff"]),
      Rcpp::as<double>(control["zeta"]),
      Rcpp::as<double>(control["tol"])
    );
    return;
  }

  if (nullModelClass == "SAGELD_NULL_Model") {
    setSAGELDobjInCPP(
      Rcpp::as<std::string>(objNull["Method"]),
      Rcpp::as<arma::mat>(objNull["XTs"]),
      Rcpp::as<arma::mat>(objNull["SS"]),
      Rcpp::as<arma::mat>(objNull["AtS"]),
      Rcpp::as<arma::mat>(objNull["Q"]),
      Rcpp::as<arma::mat>(objNull["A21"]),
      Rcpp::as<arma::mat>(objNull["TTs"]),
      Rcpp::as<arma::mat>(objNull["Tys"]),
      Rcpp::as<arma::vec>(objNull["sol"]),
      Rcpp::as<arma::vec>(objNull["blups"]),
      Rcpp::as<double>(objNull["sig"]),
      Rcpp::as<arma::vec>(objNull["Resid"]),
      Rcpp::as<arma::vec>(objNull["Resid_G"]),
      Rcpp::as<arma::vec>(objNull["Resid_GxE"]),
      Rcpp::as<arma::vec>(objNull["Resid_E"]),
      Rcpp::as<arma::vec>(objNull["Resid.unrelated.outliers"]),
      Rcpp::as<arma::vec>(objNull["Resid.unrelated.outliers_G"]),
      Rcpp::as<arma::vec>(objNull["Resid.unrelated.outliers_GxE"]),
      Rcpp::as<double>(objNull["sum_R_nonOutlier"]),
      Rcpp::as<double>(objNull["sum_R_nonOutlier_G"]),
      Rcpp::as<double>(objNull["sum_R_nonOutlier_GxE"]),
      Rcpp::as<double>(objNull["R_GRM_R"]),
      Rcpp::as<double>(objNull["R_GRM_R_G"]),
      Rcpp::as<double>(objNull["R_GRM_R_GxE"]),
      Rcpp::as<double>(objNull["R_GRM_R_G_GxE"]),
      Rcpp::as<double>(objNull["R_GRM_R_E"]),
      Rcpp::as<double>(objNull["R_GRM_R_nonOutlier"]),
      Rcpp::as<double>(objNull["R_GRM_R_nonOutlier_G"]),
      Rcpp::as<double>(objNull["R_GRM_R_nonOutlier_GxE"]),
      Rcpp::as<double>(objNull["R_GRM_R_nonOutlier_G_GxE"]),
      Rcpp::as<double>(objNull["R_GRM_R_TwoSubjOutlier"]),
      Rcpp::as<double>(objNull["R_GRM_R_TwoSubjOutlier_G"]),
      Rcpp::as<double>(objNull["R_GRM_R_TwoSubjOutlier_GxE"]),
      Rcpp::as<double>(objNull["R_GRM_R_TwoSubjOutlier_G_GxE"]),
      Rcpp::as<Rcpp::List>(objNull["TwoSubj_list"]),
      Rcpp::as<Rcpp::List>(objNull["ThreeSubj_list"]),
      Rcpp::as<arma::vec>(objNull["MAF_interval"]),
      Rcpp::as<double>(objNull["zScoreE_cutoff"]),
      Rcpp::as<double>(control["SPA_Cutoff"]),
      Rcpp::as<double>(control["zeta"]),
      Rcpp::as<double>(control["tol"])
    );
    return;
  }

  if (nullModelClass == "WtCoxG_NULL_Model") {
    setWtCoxGobjInCPP(
      Rcpp::as<arma::vec>(objNull["mresid"]),
      Rcpp::as<arma::vec>(objNull["weight"]),
      Rcpp::as<double>(control["cutoff"]),
      Rcpp::as<double>(control["SPA_Cutoff"])
    );
    return;
  }

  if (nullModelClass == "SPAsqr_NULL_Model") {
    if (ptr_gSPAsqrobj) {
      delete ptr_gSPAsqrobj;
      ptr_gSPAsqrobj = nullptr;
    }
    const arma::vec taus = Rcpp::as<arma::vec>(objNull["taus"]);
    const arma::mat residMat = Rcpp::as<arma::mat>(objNull["Resid_mat"]);
    const Rcpp::List residOutlierLst = Rcpp::as<Rcpp::List>(objNull["Resid.unrelated.outliers_lst"]);
    const arma::vec sumRnonOutlier = Rcpp::as<arma::vec>(objNull["sum_R_nonOutlier_vec"]);
    const arma::vec rGrmRnonOutlier = Rcpp::as<arma::vec>(objNull["R_GRM_R_nonOutlier_vec"]);
    const arma::vec rGrmRtwoSubj = Rcpp::as<arma::vec>(objNull["R_GRM_R_TwoSubjOutlier_vec"]);
    const arma::vec rGrmR = Rcpp::as<arma::vec>(objNull["R_GRM_R_vec"]);
    const arma::vec mafInterval = Rcpp::as<arma::vec>(objNull["MAF_interval"]);
    const Rcpp::List twoSubjLst = Rcpp::as<Rcpp::List>(objNull["TwoSubj_list_lst"]);
    const Rcpp::List cltUnionLst = Rcpp::as<Rcpp::List>(objNull["CLT_union_lst"]);
    const Rcpp::List threeSubjFamilyIdx = Rcpp::as<Rcpp::List>(objNull["ThreeSubj_family_idx_lst"]);
    const Rcpp::List threeSubjStandS = Rcpp::as<Rcpp::List>(objNull["ThreeSubj_stand_S_lst"]);

    std::vector<SPAsqr::TauFamilyNativeInput> tauData = SPAsqr::SPAsqrClass::buildTauFamilyDataFromR(
      residOutlierLst,
      twoSubjLst,
      cltUnionLst,
      threeSubjFamilyIdx,
      threeSubjStandS,
      taus.n_elem
    );

    ptr_gSPAsqrobj = new SPAsqr::SPAsqrClass(
      taus,
      residMat,
      tauData,
      sumRnonOutlier,
      rGrmRnonOutlier,
      rGrmRtwoSubj,
      rGrmR,
      mafInterval,
      Rcpp::as<double>(control["SPA_Cutoff"]),
      Rcpp::as<double>(control["zeta"]),
      Rcpp::as<double>(control["tol"])
    );
    return;
  }

  if (nullModelClass == "LEAF_NULL_Model") {
    const int nCluster = Rcpp::as<int>(objNull["Ncluster"]);
    Rcpp::IntegerVector clusterIdx = objNull["clusterIdx"];
    std::vector<std::vector<int>> clusterPositions(nCluster);
    for (int i = 0; i < clusterIdx.size(); ++i) {
      const int cluster = clusterIdx[i];
      if (cluster >= 1 && cluster <= nCluster) {
        clusterPositions[cluster - 1].push_back(i + 1);
      }
    }

    Rcpp::List clusterIdxList(nCluster);
    for (int i = 0; i < nCluster; ++i) {
      clusterIdxList[i] = Rcpp::wrap(clusterPositions[i]);
    }

    setLEAFobjInCPP(
      Rcpp::as<Rcpp::List>(objNull["residuals_list"]),
      Rcpp::as<Rcpp::List>(objNull["weights_list"]),
      clusterIdxList,
      Rcpp::as<double>(control["cutoff"]),
      Rcpp::as<double>(control["SPA_Cutoff"])
    );
    return;
  }

  throw std::runtime_error("Unsupported null model class in Marker4 setter dispatch: " + nullModelClass);
}

Main4ExtraParams buildNativeExtraParams(const std::string &nullModelClass,
                                        const Rcpp::List &objNull,
                                        const ReaderConfig &readerConfig) {
  Main4ExtraParams extra;
  extra.reader = readerConfig;

  if (nullModelClass == "SAGELD_NULL_Model") {
    extra.sageldMethod = Rcpp::as<std::string>(objNull["Method"]);
  }

  if (nullModelClass == "SPAsqr_NULL_Model") {
    extra.spasqrTaus = Rcpp::as<std::vector<double>>(objNull["taus"]);
  }

  if (nullModelClass == "WtCoxG_NULL_Model") {
    Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(objNull["mergeGenoInfo"]);
    std::vector<int64_t> genoIndex = Rcpp::as<std::vector<int64_t>>(df["genoIndex"]);
    std::vector<double> AF_ref = Rcpp::as<std::vector<double>>(df["AF_ref"]);
    std::vector<double> AN_ref = Rcpp::as<std::vector<double>>(df["AN_ref"]);
    std::vector<double> TPR = Rcpp::as<std::vector<double>>(df["TPR"]);
    std::vector<double> sigma2 = Rcpp::as<std::vector<double>>(df["sigma2"]);
    std::vector<double> pvalue_bat = Rcpp::as<std::vector<double>>(df["pvalue_bat"]);
    std::vector<double> w_ext = Rcpp::as<std::vector<double>>(df["w.ext"]);
    std::vector<double> var_ratio_w0 = Rcpp::as<std::vector<double>>(df["var.ratio.w0"]);
    std::vector<double> var_ratio_int = Rcpp::as<std::vector<double>>(df["var.ratio.int"]);
    std::vector<double> var_ratio_ext = Rcpp::as<std::vector<double>>(df["var.ratio.ext"]);

    for (size_t i = 0; i < genoIndex.size(); ++i) {
      extra.wtMap[static_cast<uint64_t>(genoIndex[i])] = WtRow{
        AF_ref[i], AN_ref[i], TPR[i], sigma2[i], pvalue_bat[i],
        w_ext[i], var_ratio_w0[i], var_ratio_int[i], var_ratio_ext[i]
      };
    }
  }

  if (nullModelClass == "LEAF_NULL_Model") {
    Rcpp::List subList = Rcpp::as<Rcpp::List>(objNull["subGenoInfo"]);
    extra.leafNcluster = Rcpp::as<int>(objNull["Ncluster"]);
    extra.leafMaps.resize(subList.size());

    for (int c = 0; c < subList.size(); ++c) {
      Rcpp::DataFrame df = Rcpp::as<Rcpp::DataFrame>(subList[c]);
      std::vector<int64_t> genoIndex = Rcpp::as<std::vector<int64_t>>(df["genoIndex"]);
      std::vector<double> AF_ref = Rcpp::as<std::vector<double>>(df["AF_ref"]);
      std::vector<double> AN_ref = Rcpp::as<std::vector<double>>(df["AN_ref"]);
      std::vector<double> TPR = Rcpp::as<std::vector<double>>(df["TPR"]);
      std::vector<double> sigma2 = Rcpp::as<std::vector<double>>(df["sigma2"]);
      std::vector<double> pvalue_bat = Rcpp::as<std::vector<double>>(df["pvalue_bat"]);
      std::vector<double> w_ext = Rcpp::as<std::vector<double>>(df["w.ext"]);
      std::vector<double> var_ratio_w0 = Rcpp::as<std::vector<double>>(df["var.ratio.w0"]);
      std::vector<double> var_ratio_int = Rcpp::as<std::vector<double>>(df["var.ratio.int"]);
      std::vector<double> var_ratio_ext = Rcpp::as<std::vector<double>>(df["var.ratio.ext"]);

      for (size_t i = 0; i < genoIndex.size(); ++i) {
        extra.leafMaps[c][static_cast<uint64_t>(genoIndex[i])] = LeafRow{
          AF_ref[i], AN_ref[i], TPR[i], sigma2[i], pvalue_bat[i],
          w_ext[i], var_ratio_w0[i], var_ratio_int[i], var_ratio_ext[i]
        };
      }
    }
  }

  return extra;
}

std::string numToStr(double x) {
  if (std::isnan(x)) return "NA";
  if (std::isinf(x)) return (x > 0 ? "Inf" : "-Inf");
  std::ostringstream oss;
  oss << std::setprecision(15) << x;
  return oss.str();
}

void appendCell(std::ostringstream &oss, const std::string &v, bool first = false) {
  if (!first) oss << '\t';
  oss << v;
}

void appendCell(std::ostringstream &oss, double v, bool first = false) {
  appendCell(oss, numToStr(v), first);
}

double cctPvalue(const std::vector<double> &pvals) {
  std::vector<double> p;
  p.reserve(pvals.size());
  for (double x : pvals) {
    if (!std::isnan(x)) p.push_back(x);
  }
  if (p.empty()) return std::numeric_limits<double>::quiet_NaN();

  double tStat = 0.0;
  for (double x : p) {
    if (x <= 0.0) return 0.0;
    if (x >= 1.0) x = 0.999;
    tStat += std::tan((0.5 - x) * M_PI);
  }
  tStat /= static_cast<double>(p.size());

  if (tStat > 1e15) {
    return (1.0 / tStat) / M_PI;
  }
  return 0.5 - std::atan(tStat) / M_PI;
}

std::string getHeader(const std::string &method,
                      const std::string &sageldMethod,
                      const std::vector<double> &taus,
                      int nPheno,
                      int nCluster) {
  std::ostringstream oss;

  if (method == "POLMM") {
    oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tPvalue\tbeta\tseBeta\tzScore";
  } else if (method == "SPACox") {
    oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tPvalue\tzScore";
  } else if (method == "SPAmix" || method == "SPAmixPlus") {
    oss << "Pheno\tMarker\tInfo\tAltFreq\tAltCounts\tMissingRate\tPvalue\tzScore";
  } else if (method == "SPAGRM") {
    oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tzScore\tPvalue\thwepval";
  } else if (method == "SAGELD") {
    if (sageldMethod == "GALLOP") {
      oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tMethod\tBeta_G\tBeta_GxE\tSE_G\tSE_GxE\tPvalue_G\tPvalue_GxE\thwepval";
    } else {
      oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\tMethod\tzScore_G\tzScore_GxE\tPvalue_G\tPvalue_GxE\thwepval";
    }
  } else if (method == "WtCoxG") {
    oss << "WtCoxG.ext\tWtCoxG.noext\tscore.ext\tscore.noext\tzscore.ext\tzscore.noext"
        << "\tAF_ref\tAN_ref\tTPR\tsigma2\tpvalue_bat\tw.ext\tvar.ratio.w0\tvar.ratio.int\tvar.ratio.ext";
  } else if (method == "SPAsqr") {
    oss << "Marker\tInfo\tAltFreq\tAltCounts\tMissingRate\thwepval";
    for (double tau : taus) {
      oss << "\tZ_tau" << numToStr(tau);
    }
    for (double tau : taus) {
      oss << "\tP_tau" << numToStr(tau);
    }
    oss << "\tP_CCT";
  } else if (method == "LEAF") {
    oss << "Marker\tInfo\tAltFreq\tMissingRate\tmeta.p_ext\tmeta.p_noext";
    for (int i = 1; i <= nCluster; ++i) {
      oss << "\tcl" << i << ".p_ext\tcl" << i << ".p_noext\tcl" << i << ".p_batch";
    }
  } else {
    throw std::runtime_error("Unsupported method in mainMarkerChunksInCPP4: " + method);
  }

  return oss.str();
}

struct ThreadContext {
  std::unique_ptr<POLMM::POLMMClass> polmm;
  std::unique_ptr<SPACox::SPACoxClass> spacox;
  std::unique_ptr<SPAmix::SPAmixClass> spamix;
  std::unique_ptr<SPAmixPlus::SPAmixPlusClass> spamixPlus;
  std::unique_ptr<SPAGRM::SPAGRMClass> spagrm;
  std::unique_ptr<SAGELD::SAGELDClass> sageld;
  std::unique_ptr<WtCoxG::WtCoxGClass> wtcoxg;
  std::unique_ptr<SPAsqr::SPAsqrClass> spasqr;
  std::unique_ptr<LEAF::LEAFClass> leaf;
};

ThreadContext makeThreadContext(const std::string &method) {
  ThreadContext ctx;
  // Each worker gets an independent method copy to avoid cross-thread mutation.
  if (method == "POLMM") {
    if (!ptr_gPOLMMobj) throw std::runtime_error("POLMM object is not initialized.");
    ctx.polmm.reset(new POLMM::POLMMClass(*ptr_gPOLMMobj));
  } else if (method == "SPACox") {
    if (!ptr_gSPACoxobj) throw std::runtime_error("SPACox object is not initialized.");
    ctx.spacox.reset(new SPACox::SPACoxClass(*ptr_gSPACoxobj));
  } else if (method == "SPAmix") {
    if (!ptr_gSPAmixobj) throw std::runtime_error("SPAmix object is not initialized.");
    ctx.spamix.reset(new SPAmix::SPAmixClass(*ptr_gSPAmixobj));
  } else if (method == "SPAmixPlus") {
    if (!ptr_gSPAmixPlusobj) throw std::runtime_error("SPAmixPlus object is not initialized.");
    ctx.spamixPlus.reset(new SPAmixPlus::SPAmixPlusClass(*ptr_gSPAmixPlusobj));
  } else if (method == "SPAGRM") {
    if (!ptr_gSPAGRMobj) throw std::runtime_error("SPAGRM object is not initialized.");
    ctx.spagrm.reset(new SPAGRM::SPAGRMClass(*ptr_gSPAGRMobj));
  } else if (method == "SAGELD") {
    if (!ptr_gSAGELDobj) throw std::runtime_error("SAGELD object is not initialized.");
    ctx.sageld.reset(new SAGELD::SAGELDClass(*ptr_gSAGELDobj));
  } else if (method == "WtCoxG") {
    if (!ptr_gWtCoxGobj) throw std::runtime_error("WtCoxG object is not initialized.");
    ctx.wtcoxg.reset(new WtCoxG::WtCoxGClass(*ptr_gWtCoxGobj));
  } else if (method == "SPAsqr") {
    if (!ptr_gSPAsqrobj) throw std::runtime_error("SPAsqr object is not initialized.");
    ctx.spasqr.reset(new SPAsqr::SPAsqrClass(*ptr_gSPAsqrobj));
  } else if (method == "LEAF") {
    if (!ptr_gLEAFobj) throw std::runtime_error("LEAF object is not initialized.");
    ctx.leaf.reset(new LEAF::LEAFClass(*ptr_gLEAFobj));
  }
  return ctx;
}

std::string getHeaderNative(const std::string &method,
                            const std::string &sageldMethod,
                            const std::vector<double> &taus,
                            int nPheno,
                            int nCluster) {
  return getHeader(method, sageldMethod, taus, nPheno, nCluster);
}

void mainMarkerChunksCore(
  const std::string &method,
  const std::vector<std::vector<uint64_t>> &chunkMarkers,
  const std::string &outputFile,
  unsigned int nThreads,
  const std::string &imputeMethod,
  double missingCutoff,
  double minMafMarker,
  double minMacMarker,
  const Main4ExtraParams &extra
) {
  if (chunkMarkers.empty()) {
    throw std::runtime_error("chunkMarkers is empty.");
  }
  if (nThreads < 1) {
    throw std::runtime_error("nThreads should be >= 1.");
  }

  const ReaderConfig &reader = extra.reader;
  if (reader.genoType.empty()) {
    throw std::runtime_error("Reader configuration is required in Main4 core.");
  }

  int nPheno = 1;
  if (method == "SPAmix") nPheno = ptr_gSPAmixobj->getNpheno();
  if (method == "SPAmixPlus") nPheno = ptr_gSPAmixPlusobj->getNpheno();
  if (method == "SPAsqr") nPheno = ptr_gSPAsqrobj->get_ntaus();

  const std::string header = getHeaderNative(method, extra.sageldMethod, extra.spasqrTaus, nPheno, extra.leafNcluster);

  std::vector<std::string> chunkOutput(chunkMarkers.size());
  std::vector<char> chunkReady(chunkMarkers.size(), 0);
  std::atomic<size_t> nextChunk(0);

  std::exception_ptr workerError = nullptr;
  std::mutex errorMutex;
  std::mutex writeMutex;
  std::mutex logMutex;
  std::condition_variable writeCv;
  bool stopWriter = false;

  std::thread writerThread([&]() {
    std::ofstream out(outputFile.c_str());
    if (!out.is_open()) {
      std::lock_guard<std::mutex> lock(errorMutex);
      workerError = std::make_exception_ptr(std::runtime_error("Cannot open output file: " + outputFile));
      stopWriter = true;
      writeCv.notify_all();
      return;
    }

    out << header << '\n';

    for (size_t i = 0; i < chunkOutput.size(); ++i) {
      std::unique_lock<std::mutex> lk(writeMutex);
      writeCv.wait(lk, [&]() { return chunkReady[i] || stopWriter; });
      if (!chunkReady[i]) {
        break;
      }
      out << chunkOutput[i];
      {
        std::lock_guard<std::mutex> lk(logMutex);
        std::fprintf(stderr, "[INFO] Writing finished: chunk %zu/%zu\n", i + 1, chunkOutput.size());
        std::fflush(stderr);
      }
    }
    out.close();
  });

  auto workerFn = [&]() {
    try {
      ThreadContext ctx = makeThreadContext(method);

      std::unique_ptr<PLINK4::PlinkReader> plinkReader;

      if (reader.genoType != "PLINK") {
        throw std::runtime_error("Unsupported reader genoType in Main4 core: " + reader.genoType);
      }
      plinkReader.reset(new PLINK4::PlinkReader(
        reader.bimFile,
        reader.famFile,
        reader.bedFile,
        reader.sampleInModel,
        reader.alleleOrder
      ));

      while (true) {
        size_t cidx = nextChunk.fetch_add(1);
        if (cidx >= chunkMarkers.size()) break;

        const auto &gIdx = chunkMarkers[cidx];
        std::ostringstream out;

        if (!gIdx.empty()) {
          plinkReader->beginSequentialBlock(gIdx.front());
        }

        std::vector<WtRow> wtChunkRows;
        std::vector<std::vector<LeafRow>> leafChunkRows;

        if (method == "WtCoxG") {
          wtChunkRows.reserve(gIdx.size());
          std::vector<double> AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext;
          AF_ref.reserve(gIdx.size());
          AN_ref.reserve(gIdx.size());
          TPR.reserve(gIdx.size());
          sigma2.reserve(gIdx.size());
          pvalue_bat.reserve(gIdx.size());
          w_ext.reserve(gIdx.size());
          var_w0.reserve(gIdx.size());
          var_int.reserve(gIdx.size());
          var_ext.reserve(gIdx.size());

          for (auto gi : gIdx) {
            auto it = extra.wtMap.find(gi);
            if (it == extra.wtMap.end()) {
              wtChunkRows.push_back(WtRow{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
            } else {
              wtChunkRows.push_back(it->second);
            }
          }

          for (const auto &r : wtChunkRows) {
            AF_ref.push_back(r.AF_ref);
            AN_ref.push_back(r.AN_ref);
            TPR.push_back(r.TPR);
            sigma2.push_back(r.sigma2);
            pvalue_bat.push_back(r.pvalue_bat);
            w_ext.push_back(r.w_ext);
            var_w0.push_back(r.var_ratio_w0);
            var_int.push_back(r.var_ratio_int);
            var_ext.push_back(r.var_ratio_ext);
          }

          ctx.wtcoxg->updateMarkerInfo(AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext);
        }

        if (method == "LEAF") {
          leafChunkRows.resize(extra.leafNcluster);
          for (int c = 0; c < extra.leafNcluster; ++c) {
            auto &rows = leafChunkRows[c];
            rows.reserve(gIdx.size());

            std::vector<double> AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext;
            AF_ref.reserve(gIdx.size());
            AN_ref.reserve(gIdx.size());
            TPR.reserve(gIdx.size());
            sigma2.reserve(gIdx.size());
            pvalue_bat.reserve(gIdx.size());
            w_ext.reserve(gIdx.size());
            var_w0.reserve(gIdx.size());
            var_int.reserve(gIdx.size());
            var_ext.reserve(gIdx.size());

            for (auto gi : gIdx) {
              auto it = extra.leafMaps[c].find(gi);
              if (it == extra.leafMaps[c].end()) {
                rows.push_back(LeafRow{NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN});
              } else {
                rows.push_back(it->second);
              }
            }

            for (const auto &r : rows) {
              AF_ref.push_back(r.AF_ref);
              AN_ref.push_back(r.AN_ref);
              TPR.push_back(r.TPR);
              sigma2.push_back(r.sigma2);
              pvalue_bat.push_back(r.pvalue_bat);
              w_ext.push_back(r.w_ext);
              var_w0.push_back(r.var_ratio_w0);
              var_int.push_back(r.var_ratio_int);
              var_ext.push_back(r.var_ratio_ext);
            }

            ctx.leaf->updateMarkerInfo(c, AF_ref, AN_ref, TPR, sigma2, pvalue_bat, w_ext, var_w0, var_int, var_ext);
          }
        }

        for (size_t i = 0; i < gIdx.size(); ++i) {
          double altFreq = NAN, altCounts = NAN, missingRate = NAN, imputeInfo = NAN;
          std::vector<uint32_t> indexForMissing, indexForNonZero;
          std::string chr, ref, alt, marker;
          uint32_t pd = 0;

          arma::vec GVec = plinkReader->getOneMarker(gIdx[i], ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo, true, indexForMissing, false, indexForNonZero, true);

          const std::string info = chr + ":" + std::to_string(pd) + ":" + ref + ":" + alt;
          const double n = static_cast<double>(GVec.size());
          const double maf = std::min(altFreq, 1.0 - altFreq);
          const double mac = 2.0 * maf * n * (1.0 - missingRate);
          bool passQC = !((missingRate > missingCutoff) || (maf < minMafMarker) || (mac < minMacMarker));

          bool flip = false;
          if (passQC) {
            flip = imputeGenoAndFlip(GVec, altFreq, indexForMissing, missingRate, imputeMethod, method);
          }

          if (method == "SPAmix" || method == "SPAmixPlus") {
            std::vector<double> p(nPheno, NAN), z(nPheno, NAN);
            if (passQC) {
              if (method == "SPAmix") {
                ctx.spamix->getMarkerPval(GVec, altFreq);
                arma::vec pTmp = ctx.spamix->getpvalVec();
                arma::vec zTmp = ctx.spamix->getzScoreVec();
                for (int j = 0; j < nPheno; ++j) {
                  p[j] = pTmp[j];
                  z[j] = zTmp[j];
                }
              } else {
                ctx.spamixPlus->getMarkerPval(GVec, altFreq);
                arma::vec pTmp = ctx.spamixPlus->getpvalVec();
                arma::vec zTmp = ctx.spamixPlus->getzScoreVec();
                for (int j = 0; j < nPheno; ++j) {
                  p[j] = pTmp[j];
                  z[j] = zTmp[j];
                }
              }
            }

            for (int j = 0; j < nPheno; ++j) {
              std::ostringstream row;
              appendCell(row, "pheno_" + std::to_string(j + 1), true);
              appendCell(row, marker);
              appendCell(row, info);
              appendCell(row, altFreq);
              appendCell(row, altCounts);
              appendCell(row, missingRate);
              appendCell(row, p[j]);
              appendCell(row, z[j]);
              out << row.str() << '\n';
            }
            continue;
          }

          if (method == "POLMM") {
            double Beta = NAN, seBeta = NAN, pval = NAN, zScore = NAN;
            if (passQC) {
              ctx.polmm->getMarkerPval(GVec, Beta, seBeta, pval, altFreq, zScore);
              Beta = Beta * (1.0 - 2.0 * static_cast<double>(flip));
            }
            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, pval);
            appendCell(row, Beta);
            appendCell(row, seBeta);
            appendCell(row, zScore);
            out << row.str() << '\n';
            continue;
          }

          if (method == "SPACox") {
            double pval = NAN, zScore = NAN;
            if (passQC) {
              pval = ctx.spacox->getMarkerPval(GVec, altFreq, zScore);
            }
            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, pval);
            appendCell(row, zScore);
            out << row.str() << '\n';
            continue;
          }

          if (method == "SPAGRM") {
            double pval = NAN, zScore = NAN, hwepval = NAN;
            if (passQC) {
              pval = ctx.spagrm->getMarkerPval(GVec, altFreq, zScore, hwepval);
            }
            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, zScore);
            appendCell(row, pval);
            appendCell(row, hwepval);
            out << row.str() << '\n';
            continue;
          }

          if (method == "SAGELD") {
            double hwepval = NAN;
            double zG = NAN, zGxE = NAN, pG = NAN, pGxE = NAN;
            double bG = NAN, bGxE = NAN, seG = NAN, seGxE = NAN;
            if (passQC) {
              ctx.sageld->getMarkerPval(GVec, altFreq, hwepval);
              arma::vec pTmp = ctx.sageld->getpvalVec();
              arma::vec zTmp = ctx.sageld->getzScoreVec();
              arma::vec bTmp = ctx.sageld->getBetaVec();
              arma::vec seTmp = ctx.sageld->getseBetaVec();
              zG = zTmp[0];
              zGxE = zTmp[1];
              pG = pTmp[0];
              pGxE = pTmp[1];
              bG = bTmp[0] * (1.0 - 2.0 * static_cast<double>(flip));
              bGxE = bTmp[1] * (1.0 - 2.0 * static_cast<double>(flip));
              seG = seTmp[0];
              seGxE = seTmp[1];
            }

            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, extra.sageldMethod);
            if (extra.sageldMethod == "GALLOP") {
              appendCell(row, bG);
              appendCell(row, bGxE);
              appendCell(row, seG);
              appendCell(row, seGxE);
              appendCell(row, pG);
              appendCell(row, pGxE);
            } else {
              appendCell(row, zG);
              appendCell(row, zGxE);
              appendCell(row, pG);
              appendCell(row, pGxE);
            }
            appendCell(row, hwepval);
            out << row.str() << '\n';
            continue;
          }

          if (method == "WtCoxG") {
            double pExt = NAN, pNoext = NAN, sExt = NAN, sNoext = NAN, zExt = NAN, zNoext = NAN;
            if (passQC) {
              arma::vec pTmp = ctx.wtcoxg->getpvalVec(GVec, static_cast<int>(i));
              arma::vec sTmp = ctx.wtcoxg->getScoreVec();
              arma::vec zTmp = ctx.wtcoxg->getZScoreVec();
              pExt = pTmp[0];
              pNoext = pTmp[1];
              sExt = sTmp[0];
              sNoext = sTmp[1];
              zExt = zTmp[0];
              zNoext = zTmp[1];
            }

            const WtRow wr = wtChunkRows[i];
            std::ostringstream row;
            appendCell(row, pExt, true);
            appendCell(row, pNoext);
            appendCell(row, sExt);
            appendCell(row, sNoext);
            appendCell(row, zExt);
            appendCell(row, zNoext);
            appendCell(row, wr.AF_ref);
            appendCell(row, wr.AN_ref);
            appendCell(row, wr.TPR);
            appendCell(row, wr.sigma2);
            appendCell(row, wr.pvalue_bat);
            appendCell(row, wr.w_ext);
            appendCell(row, wr.var_ratio_w0);
            appendCell(row, wr.var_ratio_int);
            appendCell(row, wr.var_ratio_ext);
            out << row.str() << '\n';
            continue;
          }

          if (method == "SPAsqr") {
            double hwepval = NAN;
            std::vector<double> z(extra.spasqrTaus.size(), NAN), p(extra.spasqrTaus.size(), NAN);
            if (passQC) {
              arma::vec zTmp;
              arma::vec pTmp = ctx.spasqr->getMarkerPval(GVec, altFreq, zTmp, hwepval);
              for (size_t j = 0; j < p.size(); ++j) {
                p[j] = pTmp[j];
                z[j] = zTmp[j];
              }
            }
            double pCCT = cctPvalue(p);

            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, altCounts);
            appendCell(row, missingRate);
            appendCell(row, hwepval);
            for (double v : z) appendCell(row, v);
            for (double v : p) appendCell(row, v);
            appendCell(row, pCCT);
            out << row.str() << '\n';
            continue;
          }

          if (method == "LEAF") {
            std::vector<double> pExt(extra.leafNcluster, NAN), pNoext(extra.leafNcluster, NAN);
            std::vector<double> sExt(extra.leafNcluster, NAN), sNoext(extra.leafNcluster, NAN);

            if (passQC) {
              arma::vec all = ctx.leaf->getMarkerZSP(GVec, static_cast<int>(i));
              const int nOut = 2 * extra.leafNcluster;
              for (int c = 0; c < extra.leafNcluster; ++c) {
                int extIdx = 2 * c;
                int noextIdx = 2 * c + 1;
                sExt[c] = all[nOut + extIdx];
                sNoext[c] = all[nOut + noextIdx];
                pExt[c] = all[2 * nOut + extIdx];
                pNoext[c] = all[2 * nOut + noextIdx];
              }
            }

            auto metaP = [&](const std::vector<double> &scores, const std::vector<double> &pvals) {
              double sumScore = 0.0;
              double sumVar = 0.0;
              boost::math::chi_squared dist(1.0);

              for (size_t c = 0; c < scores.size(); ++c) {
                double sc = scores[c];
                double pv = pvals[c];
                if (std::isnan(sc) || std::isnan(pv) || pv <= 0.0 || pv >= 1.0) continue;
                double chisq = boost::math::quantile(dist, 1.0 - pv);
                if (chisq < 1e-30) chisq = 1e-30;
                double var = (sc * sc) / chisq;
                if (std::isnan(var)) var = 0.0;
                sumScore += sc;
                sumVar += var;
              }

              if (sumVar <= 0.0) return std::numeric_limits<double>::quiet_NaN();
              double z = sumScore / std::sqrt(sumVar);
              return std::erfc(std::fabs(z) / std::sqrt(2.0));
            };

            double metaExt = metaP(sExt, pExt);
            double metaNoext = metaP(sNoext, pNoext);

            std::ostringstream row;
            appendCell(row, marker, true);
            appendCell(row, info);
            appendCell(row, altFreq);
            appendCell(row, missingRate);
            appendCell(row, metaExt);
            appendCell(row, metaNoext);

            for (int c = 0; c < extra.leafNcluster; ++c) {
              appendCell(row, pExt[c]);
              appendCell(row, pNoext[c]);
              appendCell(row, leafChunkRows[c][i].pvalue_bat);
            }
            out << row.str() << '\n';
            continue;
          }
        }

        {
          std::lock_guard<std::mutex> lk(writeMutex);
          chunkOutput[cidx] = out.str();
          chunkReady[cidx] = 1;
        }
        {
          std::lock_guard<std::mutex> lk(logMutex);
          std::fprintf(stderr, "[INFO] Calculation finished: chunk %zu/%zu\n", cidx + 1, chunkMarkers.size());
          std::fflush(stderr);
        }
        writeCv.notify_all();
      }
    } catch (...) {
      {
        std::lock_guard<std::mutex> lock(errorMutex);
        if (!workerError) workerError = std::current_exception();
      }
      {
        std::lock_guard<std::mutex> lk(writeMutex);
        stopWriter = true;
      }
      writeCv.notify_all();
    }
  };

  nThreads = std::min<unsigned int>(nThreads, static_cast<unsigned int>(chunkMarkers.size()));
  const bool runWorkerInline = (nThreads == 1 && (method == "WtCoxG" || method == "LEAF"));
  if (runWorkerInline) {
    // WtCoxG/LEAF use WtCoxG internals that call R-level mvtnorm; keep execution
    // on the main thread when single-threaded to avoid R API calls in a worker thread.
    workerFn();
  } else {
    std::vector<std::thread> workers;
    workers.reserve(nThreads);
    for (unsigned int t = 0; t < nThreads; ++t) {
      workers.emplace_back(workerFn);
    }
    for (auto &th : workers) {
      th.join();
    }
  }

  {
    std::lock_guard<std::mutex> lk(writeMutex);
    stopWriter = true;
  }
  writeCv.notify_all();
  writerThread.join();

  if (workerError) {
    std::rethrow_exception(workerError);
  }
}

} // namespace

// [[Rcpp::export]]
void prepareAndRunMarkerInCPP4(
  const Rcpp::List t_objNull,
  const std::string t_GenoFile,
  const std::string t_OutputFile,
  const Rcpp::Nullable<Rcpp::CharacterVector> t_GenoFileIndex,
  const Rcpp::List t_control,
  const unsigned int t_nThreads
) {
  try {
    // Stage A: decode method and ensure this pure-native path receives PLINK input.
    if (!t_objNull.hasAttribute("class"))
      throw std::runtime_error("objNull must have an R class attribute.");
    const Rcpp::CharacterVector _nullClasses = t_objNull.attr("class");
    if (_nullClasses.size() < 1)
      throw std::runtime_error("objNull must have a non-empty R class attribute.");
    const std::string nullModelClass = Rcpp::as<std::string>(_nullClasses[0]);
    const std::string method = getMethodFromNullModelClass(nullModelClass);
    const int nMarkersEachChunk = Rcpp::as<int>(t_control["nMarkersEachChunk"]);
    const std::string genoFileLower = toLowerCopy(t_GenoFile);
    if (genoFileLower.size() < 4 || genoFileLower.substr(genoFileLower.size() - 4) != ".bed") {
      throw std::runtime_error("Pure C++ Marker4 currently supports only PLINK .bed input. BGEN is not supported in this path yet.");
    }

    // Parse R control/config at boundary, then pass native C++ config downstream.
    auto hasNonNull = [](const Rcpp::List& lst, const char* name) -> bool {
      return lst.containsElementNamed(name) && !Rcpp::as<Rcpp::RObject>(lst[name]).isNULL();
    };
    std::vector<std::string> genoFileIndex;
    if (t_GenoFileIndex.isNotNull()) {
      Rcpp::CharacterVector indexVec(t_GenoFileIndex);
      genoFileIndex.reserve(indexVec.size());
      for (int i = 0; i < indexVec.size(); ++i) {
        genoFileIndex.push_back(Rcpp::as<std::string>(indexVec[i]));
      }
    }
    const std::string alleleOrder =
      hasNonNull(t_control, "AlleleOrder") ?
      Rcpp::as<std::string>(t_control["AlleleOrder"]) :
      "alt-first";
    MarkerFilterConfig filterConfig;
    if (hasNonNull(t_control, "IDsToIncludeFile")) filterConfig.idsToIncludeFile = Rcpp::as<std::string>(t_control["IDsToIncludeFile"]);
    if (hasNonNull(t_control, "RangesToIncludeFile")) filterConfig.rangesToIncludeFile = Rcpp::as<std::string>(t_control["RangesToIncludeFile"]);
    if (hasNonNull(t_control, "IDsToExcludeFile")) filterConfig.idsToExcludeFile = Rcpp::as<std::string>(t_control["IDsToExcludeFile"]);
    if (hasNonNull(t_control, "RangesToExcludeFile")) filterConfig.rangesToExcludeFile = Rcpp::as<std::string>(t_control["RangesToExcludeFile"]);

    // Stage B: initialize method object exactly once from objNull/control.
    setMarkerObjectForMethod(nullModelClass, t_objNull, t_control);

    // Stage C: native genotype preparation (no R callback) and marker chunking.
    const std::vector<std::string> requestedSamples =
      Rcpp::as<std::vector<std::string>>(Rcpp::as<Rcpp::CharacterVector>(t_objNull["subjData"]));
    ReaderConfig readerConfig = buildPlinkReaderConfig(
      t_GenoFile,
      genoFileIndex,
      requestedSamples,
      alleleOrder
    );

    std::vector<PlinkMarkerInfo> markerInfo = readPlinkMarkerInfo(readerConfig.bimFile, readerConfig.alleleOrder);
    const std::vector<std::string> famSamples = readPlinkSampleIds(readerConfig.famFile);
    validateRequestedSamples(requestedSamples, famSamples);
    markerInfo = applyPlinkMarkerFilters(markerInfo, filterConfig);

    if (markerInfo.empty()) {
      throw std::runtime_error("No markers remain after PLINK marker filtering.");
    }
    std::vector<std::vector<uint64_t>> chunkIndices = buildPlinkChunkIndexList(markerInfo, nMarkersEachChunk);

    Rcpp::Rcout << "Number of markers to test: " << markerInfo.size() << std::endl;
    Rcpp::Rcout << "Genotype type: PLINK" << std::endl;
    Rcpp::Rcout << "Number of markers in each chunk: " << nMarkersEachChunk << std::endl;
    Rcpp::Rcout << "Number of chunks for all markers: " << chunkIndices.size() << std::endl;
    Rcpp::Rcout << "Number of threads: " << t_nThreads << std::endl;

    // Stage D: convert any remaining method metadata once, then run the pure C++ core.
    Main4ExtraParams extraParams = buildNativeExtraParams(nullModelClass, t_objNull, readerConfig);
    mainMarkerChunksCore(
      method,
      chunkIndices,
      t_OutputFile,
      t_nThreads,
      Rcpp::as<std::string>(t_control["impute_method"]),
      Rcpp::as<double>(t_control["missing_cutoff"]),
      Rcpp::as<double>(t_control["min_maf_marker"]),
      Rcpp::as<double>(t_control["min_mac_marker"]),
      extraParams
    );
  } catch (const std::exception &e) {
    Rcpp::stop(e.what());
  }
}
