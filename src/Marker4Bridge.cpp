
// Marker4Bridge.cpp -- R/Rcpp boundary: per-method Rcpp entry points
//
// This header contains the full definitions of all runMarkerInCPP_FOO functions.
//
// Functions in this file:
//
//   [Array splitters]  (static, file-local)
//     splitVec        — split a concatenated arma::vec into a vector of arma::vec
//     splitUvec       — split a concatenated arma::uvec into a vector of arma::uvec
//     splitMat        — split a concatenated arma::mat (row-stacked) into a vector of arma::mat
//     matRowsToVecs   — convert each row of a matrix into a separate arma::vec
//
//   [Internal engine runner]  (static, file-local)
//     runEngineCore   — build MarkerFilterConfig/ReaderConfig, read+filter .bim/.fam,
//                       build chunks, print summary, and invoke mainMarkerChunksCore
//
//   [Rcpp entry points]
//     runMarkerInCPP_POLMM      — construct POLMMClass
//     runMarkerInCPP_WtCoxG     — construct WtCoxGClass
//     runMarkerInCPP_LEAF       — construct LEAFClass
//     runMarkerInCPP_SPAGRM     — construct SPAGRMClass
//     runMarkerInCPP_SAGELD     — construct SAGELDClass
//     runMarkerInCPP_SPAsqr     — construct SPAsqrClass
//     runMarkerInCPP_SPACox     — construct SPACoxClass
//     runMarkerInCPP_SPAmix     — construct SPAmixClass
//     runMarkerInCPP_SPAmixPlus — construct SPAmixPlusClass

#define MARKER4BRIDGE_H

#include <RcppArmadillo.h>

#include "Marker4.h"
#include "PLINK4.h"
#include "POLMM.h"
#include "WtCoxG.h"
#include "LEAF.h"
#include "SPAGRM.h"
#include "SAGELD.h"
#include "SPAsqr.h"
#include "SPACox.h"
#include "SPAmix.h"
#include "SPAmixPlus.h"

// Global method pointers — defined in Marker4Engine.cpp.
extern POLMM::POLMMClass*            ptr_gPOLMMobj;
extern WtCoxG::WtCoxGClass*          ptr_gWtCoxGobj;
extern LEAF::LEAFClass*              ptr_gLEAFobj;
extern SPACox::SPACoxClass*          ptr_gSPACoxobj;
extern SPAmix::SPAmixClass*          ptr_gSPAmixobj;
extern SPAGRM::SPAGRMClass*          ptr_gSPAGRMobj;
extern SAGELD::SAGELDClass*          ptr_gSAGELDobj;
extern SPAsqr::SPAsqrClass*          ptr_gSPAsqrobj;
extern SPAmixPlus::SPAmixPlusClass*  ptr_gSPAmixPlusobj;

namespace {

// =========================================================================
// Array splitters: reconstruct R-concatenated arma containers
// =========================================================================

static std::vector<arma::vec> splitVec(
  const arma::vec& data,
  const arma::uvec& lens
) {
  std::vector<arma::vec> out(lens.n_elem);
  arma::uword off = 0;
  for (arma::uword i = 0; i < lens.n_elem; ++i) {
    if (lens[i] > 0)
      out[i] = data.subvec(off, off + lens[i] - 1);
    off += lens[i];
  }
  return out;
}

static std::vector<arma::uvec> splitUvec(
  const arma::uvec& data,
  const arma::uvec& lens
) {
  std::vector<arma::uvec> out(lens.n_elem);
  arma::uword off = 0;
  for (arma::uword i = 0; i < lens.n_elem; ++i) {
    if (lens[i] > 0)
      out[i] = data.subvec(off, off + lens[i] - 1);
    off += lens[i];
  }
  return out;
}

static std::vector<arma::mat> splitMat(
  const arma::mat& data,
  const arma::uvec& nrows
) {
  std::vector<arma::mat> out(nrows.n_elem);
  arma::uword off = 0;
  for (arma::uword i = 0; i < nrows.n_elem; ++i) {
    if (nrows[i] > 0)
      out[i] = data.rows(off, off + nrows[i] - 1);
    off += nrows[i];
  }
  return out;
}

// Each row of the matrix becomes a column vector in the returned vector.
static std::vector<arma::vec> matRowsToVecs(const arma::mat& m) {
  std::vector<arma::vec> out(m.n_rows);
  for (arma::uword i = 0; i < m.n_rows; ++i)
    out[i] = m.row(i).t();
  return out;
}

// =========================================================================
// Internal engine runner (not exported to R)
// =========================================================================

static void runEngineCore(
  const std::string& method,
  const std::string& bedFile,
  const std::string& bimFile,
  const std::string& famFile,
  const std::string& outputFile,
  const std::vector<std::string>& subjData,
  const std::string& AlleleOrder,
  int nMarkersEachChunk,
  unsigned int nthreads,
  const std::string& impute_method,
  double missing_cutoff,
  double min_maf_marker,
  double min_mac_marker,
  const std::string& IDsToIncludeFile,
  const std::string& RangesToIncludeFile,
  const std::string& IDsToExcludeFile,
  const std::string& RangesToExcludeFile
) {

  PLINK4::PlinkData plinkData(bedFile, bimFile, famFile, subjData, AlleleOrder);

  PLINK4::MarkerFilterConfig filterConfig;
  filterConfig.IDsToIncludeFile    = IDsToIncludeFile;
  filterConfig.RangesToIncludeFile = RangesToIncludeFile;
  filterConfig.IDsToExcludeFile    = IDsToExcludeFile;
  filterConfig.RangesToExcludeFile = RangesToExcludeFile;

  auto markerInfo = plinkData.getFilteredMarkers(filterConfig); // bim lines for included markers

  if (markerInfo.empty())
    throw std::runtime_error("No markers remain after PLINK marker filtering.");

  std::vector<std::vector<uint64_t>> chunkIndices = PLINK4::PlinkData::buildChunks(markerInfo, nMarkersEachChunk);

  Rprintf("Number of subjects in the input file: %u\n", plinkData.nSamplesInFile());
  Rprintf("Number of subjects to test: %u\n", plinkData.nSamplesUsed());
  Rprintf("Number of markers in the input file: %u\n", plinkData.nMarkers());
  Rprintf("Number of markers to test: %zu\n", markerInfo.size());
  Rprintf("Number of markers in each chunk: %d\n", nMarkersEachChunk);
  Rprintf("Number of chunks for all markers: %zu\n", chunkIndices.size());
  Rprintf("Number of threads: %u\n", nthreads);

  Marker4::mainMarkerChunksCore(
    method, chunkIndices, outputFile, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker, plinkData
  );
}

} // anonymous namespace

// =========================================================================
// Per-method Rcpp entry points
// =========================================================================

// [[Rcpp::export("runMarkerInCPP.POLMM")]]
void runMarkerInCPP_POLMM(
    arma::mat muMat, arma::mat iRMat, arma::mat Cova, arma::uvec yVec,
    double tau, double varRatio, double SPA_Cutoff,
    std::string bedFile, std::string bimFile, std::string famFile,
    std::string outputFile,
    std::vector<std::string> subjData, std::string AlleleOrder,
    int nMarkersEachChunk, unsigned int nthreads,
    std::string impute_method, double missing_cutoff,
    double min_maf_marker, double min_mac_marker,
    std::string IDsToIncludeFile, std::string RangesToIncludeFile,
    std::string IDsToExcludeFile, std::string RangesToExcludeFile
  ) {
  if (ptr_gPOLMMobj) {
    delete ptr_gPOLMMobj;
  }
  ptr_gPOLMMobj = new POLMM::POLMMClass(
    muMat, iRMat, Cova, yVec,
    arma::sp_mat(), tau, false, 0.001, 100,
    varRatio, SPA_Cutoff, false
  );

  runEngineCore(
    "POLMM", bedFile, bimFile, famFile, outputFile,
    subjData, AlleleOrder, nMarkersEachChunk, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker,
    IDsToIncludeFile, RangesToIncludeFile,
    IDsToExcludeFile, RangesToExcludeFile
  );
}

// [[Rcpp::export("runMarkerInCPP.WtCoxG")]]
void runMarkerInCPP_WtCoxG(
  arma::vec R, arma::vec w, double cutoff, double SPA_Cutoff,
  arma::vec wt_genoIndex, arma::vec wt_AF_ref, arma::vec wt_AN_ref,
  arma::vec wt_TPR, arma::vec wt_sigma2, arma::vec wt_pvalue_bat,
  arma::vec wt_w_ext, arma::vec wt_var_w0, arma::vec wt_var_int,
  arma::vec wt_var_ext,
  std::string bedFile, std::string bimFile, std::string famFile,
  std::string outputFile,
  std::vector<std::string> subjData, std::string AlleleOrder,
  int nMarkersEachChunk, unsigned int nthreads,
  std::string impute_method, double missing_cutoff,
  double min_maf_marker, double min_mac_marker,
  std::string IDsToIncludeFile, std::string RangesToIncludeFile,
  std::string IDsToExcludeFile, std::string RangesToExcludeFile
) {
    std::unordered_map<uint64_t, WtCoxG::WtCoxGClass::RefInfo> refMap;
    refMap.reserve(wt_genoIndex.n_elem);
    for (arma::uword i = 0; i < wt_genoIndex.n_elem; ++i) {
      refMap[static_cast<uint64_t>(wt_genoIndex[i])] = WtCoxG::WtCoxGClass::RefInfo{
        wt_AF_ref[i], wt_AN_ref[i], wt_TPR[i], wt_sigma2[i], wt_pvalue_bat[i],
        wt_w_ext[i], wt_var_w0[i], wt_var_int[i], wt_var_ext[i]
      };
    }

  if (ptr_gWtCoxGobj) {
    delete ptr_gWtCoxGobj;
  }
  ptr_gWtCoxGobj = new WtCoxG::WtCoxGClass(R, w, cutoff, SPA_Cutoff);
  ptr_gWtCoxGobj->setRefMap(std::move(refMap));

  runEngineCore(
    "WtCoxG", bedFile, bimFile, famFile, outputFile,
    subjData, AlleleOrder, nMarkersEachChunk, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker,
    IDsToIncludeFile, RangesToIncludeFile,
    IDsToExcludeFile, RangesToExcludeFile
  );
}

// [[Rcpp::export("runMarkerInCPP.LEAF")]]
void runMarkerInCPP_LEAF(
  arma::vec residuals_all, arma::uvec residuals_lens,
  arma::vec weights_all, arma::uvec weights_lens,
  arma::uvec clusterIdx_all, arma::uvec clusterIdx_lens,
  double cutoff, double SPA_Cutoff,
  int nCluster,
  arma::vec leaf_genoIndex, arma::vec leaf_AF_ref, arma::vec leaf_AN_ref,
  arma::vec leaf_TPR, arma::vec leaf_sigma2, arma::vec leaf_pvalue_bat,
  arma::vec leaf_w_ext, arma::vec leaf_var_w0, arma::vec leaf_var_int,
  arma::vec leaf_var_ext, arma::uvec leaf_nrows,
  std::string bedFile, std::string bimFile, std::string famFile,
  std::string outputFile,
  std::vector<std::string> subjData, std::string AlleleOrder,
  int nMarkersEachChunk, unsigned int nthreads,
  std::string impute_method, double missing_cutoff,
  double min_maf_marker, double min_mac_marker,
  std::string IDsToIncludeFile, std::string RangesToIncludeFile,
  std::string IDsToExcludeFile, std::string RangesToExcludeFile
) {
  auto residuals      = splitVec(residuals_all, residuals_lens);
  auto weights        = splitVec(weights_all, weights_lens);
  auto clusterIdxVecs = splitUvec(clusterIdx_all, clusterIdx_lens);
  for (auto& v : clusterIdxVecs) v -= 1;  // R 1-based -> C++ 0-based

  std::vector<std::shared_ptr<const std::unordered_map<uint64_t, WtCoxG::WtCoxGClass::RefInfo>>> refMaps(nCluster);
  arma::uword off = 0;
  for (int c = 0; c < nCluster; ++c) {
    auto map = std::make_shared<std::unordered_map<uint64_t, WtCoxG::WtCoxGClass::RefInfo>>();
    arma::uword nr = leaf_nrows[c];
    map->reserve(nr);
    for (arma::uword i = 0; i < nr; ++i) {
      arma::uword idx = off + i;
      (*map)[static_cast<uint64_t>(leaf_genoIndex[idx])] = WtCoxG::WtCoxGClass::RefInfo{
        leaf_AF_ref[idx], leaf_AN_ref[idx], leaf_TPR[idx],
        leaf_sigma2[idx], leaf_pvalue_bat[idx], leaf_w_ext[idx],
        leaf_var_w0[idx], leaf_var_int[idx], leaf_var_ext[idx]
      };
    }
    off += nr;
    refMaps[c] = map;
  }

  if (ptr_gLEAFobj) {
    delete ptr_gLEAFobj;
  }
  ptr_gLEAFobj = new LEAF::LEAFClass(residuals, weights, clusterIdxVecs, cutoff, SPA_Cutoff);
  ptr_gLEAFobj->setRefMaps(std::move(refMaps));

  runEngineCore(
    "LEAF", bedFile, bimFile, famFile, outputFile,
    subjData, AlleleOrder, nMarkersEachChunk, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker,
    IDsToIncludeFile, RangesToIncludeFile,
    IDsToExcludeFile, RangesToExcludeFile
  );
}

// [[Rcpp::export("runMarkerInCPP.SPAGRM")]]
void runMarkerInCPP_SPAGRM(
  arma::vec resid, arma::vec resid_unrelated_outliers,
  double sum_R_nonOutlier, double R_GRM_R_nonOutlier,
  double R_GRM_R_TwoSubjOutlier, double R_GRM_R,
  arma::vec MAF_interval,
  arma::mat twoSubj_resid, arma::mat twoSubj_rho,
  arma::vec threeSubj_standS_all, arma::uvec threeSubj_standS_lens,
  arma::mat threeSubj_CLT_all, arma::uvec threeSubj_CLT_nrows,
  double SPA_Cutoff, double zeta, double tol,
  std::string bedFile, std::string bimFile, std::string famFile,
  std::string outputFile,
  std::vector<std::string> subjData, std::string AlleleOrder,
  int nMarkersEachChunk, unsigned int nthreads,
  std::string impute_method, double missing_cutoff,
  double min_maf_marker, double min_mac_marker,
  std::string IDsToIncludeFile, std::string RangesToIncludeFile,
  std::string IDsToExcludeFile, std::string RangesToExcludeFile
) {
  auto twoResid    = matRowsToVecs(twoSubj_resid);
  auto twoRho      = matRowsToVecs(twoSubj_rho);
  auto threeStandS = splitVec(threeSubj_standS_all, threeSubj_standS_lens);
  auto threeCLT    = splitMat(threeSubj_CLT_all, threeSubj_CLT_nrows);

  if (ptr_gSPAGRMobj) {
    delete ptr_gSPAGRMobj;
  }
  ptr_gSPAGRMobj = new SPAGRM::SPAGRMClass(
    resid, resid_unrelated_outliers,
    sum_R_nonOutlier, R_GRM_R_nonOutlier, R_GRM_R_TwoSubjOutlier, R_GRM_R,
    MAF_interval,
    std::move(twoResid), std::move(twoRho),
    std::move(threeStandS), std::move(threeCLT),
    SPA_Cutoff, zeta, tol
  );

  runEngineCore(
    "SPAGRM", bedFile, bimFile, famFile, outputFile,
    subjData, AlleleOrder, nMarkersEachChunk, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker,
    IDsToIncludeFile, RangesToIncludeFile,
    IDsToExcludeFile, RangesToExcludeFile
  );
}

// [[Rcpp::export("runMarkerInCPP.SAGELD")]]
void runMarkerInCPP_SAGELD(
  std::string Method, arma::mat XTs, arma::mat SS, arma::mat AtS,
  arma::mat Q, arma::mat A21, arma::mat TTs, arma::mat Tys,
  arma::vec sol, arma::vec blups, double sig,
  arma::vec resid, arma::vec resid_G,
  arma::vec resid_GxE, arma::vec resid_E,
  arma::vec resid_unrelated_outliers,
  arma::vec resid_unrelated_outliers_G,
  arma::vec resid_unrelated_outliers_GxE,
  double sum_R_nonOutlier, double sum_R_nonOutlier_G, double sum_R_nonOutlier_GxE,
  arma::vec R_GRM_R, arma::vec R_GRM_R_nonOutlier, arma::vec R_GRM_R_TwoSubjOutlier,
  arma::mat twoSubj_Resid, arma::mat twoSubj_Rho,
  arma::mat twoSubj_Resid_G, arma::mat twoSubj_Resid_GxE,
  arma::vec threeSubj_standS_all, arma::uvec threeSubj_standS_lens,
  arma::vec threeSubj_standS_G_all, arma::uvec threeSubj_standS_G_lens,
  arma::vec threeSubj_standS_GxE_all, arma::uvec threeSubj_standS_GxE_lens,
  arma::mat threeSubj_CLT_all, arma::uvec threeSubj_CLT_nrows, arma::vec MAF_interval,
  double zScoreE_cutoff, double SPA_Cutoff, double zeta, double tol,
  std::string bedFile, std::string bimFile, std::string famFile, std::string outputFile,
  std::vector<std::string> subjData, std::string AlleleOrder,
  int nMarkersEachChunk, unsigned int nthreads, std::string impute_method,
  double missing_cutoff, double min_maf_marker, double min_mac_marker,
  std::string IDsToIncludeFile, std::string RangesToIncludeFile,
  std::string IDsToExcludeFile, std::string RangesToExcludeFile
) {

  int nTwo = twoSubj_Resid.n_rows;
  std::vector<SAGELD::SAGELDClass::TwoSubjFamily> twoSubj(nTwo);
  for (int i = 0; i < nTwo; ++i) {
    twoSubj[i].Resid     = twoSubj_Resid.row(i).t();
    twoSubj[i].Rho       = twoSubj_Rho.row(i).t();
    twoSubj[i].Resid_G   = twoSubj_Resid_G.row(i).t();
    twoSubj[i].Resid_GxE = twoSubj_Resid_GxE.row(i).t();
  }
  auto standS     = splitVec(threeSubj_standS_all, threeSubj_standS_lens);
  auto standS_G   = splitVec(threeSubj_standS_G_all, threeSubj_standS_G_lens);
  auto standS_GxE = splitVec(threeSubj_standS_GxE_all, threeSubj_standS_GxE_lens);
  auto CLTs       = splitMat(threeSubj_CLT_all, threeSubj_CLT_nrows);
  int nThree = standS.size();
  std::vector<SAGELD::SAGELDClass::ThreeSubjFamily> threeSubj(nThree);
  for (int i = 0; i < nThree; ++i) {
    threeSubj[i].CLT         = std::move(CLTs[i]);
    threeSubj[i].stand_S     = std::move(standS[i]);
    threeSubj[i].stand_S_G   = std::move(standS_G[i]);
    threeSubj[i].stand_S_GxE = std::move(standS_GxE[i]);
  }

  if (ptr_gSAGELDobj) {
    delete ptr_gSAGELDobj;
  }
  ptr_gSAGELDobj = new SAGELD::SAGELDClass(
    Method, XTs, SS, AtS, Q, A21, TTs, Tys,
    sol, blups, sig,
    resid, resid_G, resid_GxE, resid_E,
    resid_unrelated_outliers, resid_unrelated_outliers_G, resid_unrelated_outliers_GxE,
    sum_R_nonOutlier, sum_R_nonOutlier_G, sum_R_nonOutlier_GxE,
    R_GRM_R, R_GRM_R_nonOutlier, R_GRM_R_TwoSubjOutlier,
    std::move(twoSubj), std::move(threeSubj),
    MAF_interval, zScoreE_cutoff, SPA_Cutoff, zeta, tol
  );

  runEngineCore(
    "SAGELD", bedFile, bimFile, famFile, outputFile,
    subjData, AlleleOrder, nMarkersEachChunk, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker,
    IDsToIncludeFile, RangesToIncludeFile,
    IDsToExcludeFile, RangesToExcludeFile
  );
}

// [[Rcpp::export("runMarkerInCPP.SPAsqr")]]
void runMarkerInCPP_SPAsqr(
  arma::vec taus, arma::mat Resid_mat,
  arma::vec resid_outlier_all, arma::uvec resid_outlier_lens,
  arma::vec sum_R_nonOutlier_vec, arma::vec R_GRM_R_nonOutlier_vec,
  arma::vec R_GRM_R_TwoSubjOutlier_vec, arma::vec R_GRM_R_vec,
  arma::vec MAF_interval,
  arma::mat twoSubj_resid_all, arma::mat twoSubj_rho_all,
  arma::uvec twoSubj_perTau,
  arma::vec threeSubj_standS_all, arma::uvec threeSubj_standS_lens,
  arma::mat threeSubj_CLT_all, arma::uvec threeSubj_CLT_nrows,
  arma::uvec threeSubj_perTau,
  double SPA_Cutoff, double zeta, double tol,
  std::string bedFile, std::string bimFile, std::string famFile,
  std::string outputFile,
  std::vector<std::string> subjData, std::string AlleleOrder,
  int nMarkersEachChunk, unsigned int nthreads,
  std::string impute_method, double missing_cutoff,
  double min_maf_marker, double min_mac_marker,
  std::string IDsToIncludeFile, std::string RangesToIncludeFile,
  std::string IDsToExcludeFile, std::string RangesToExcludeFile
) {
  int ntaus            = static_cast<int>(taus.n_elem);
  auto residOutliers  = splitVec(resid_outlier_all, resid_outlier_lens);
  auto allTwoResid    = matRowsToVecs(twoSubj_resid_all);
  auto allTwoRho      = matRowsToVecs(twoSubj_rho_all);
  auto allThreeStandS = splitVec(threeSubj_standS_all, threeSubj_standS_lens);
  auto allThreeCLT    = splitMat(threeSubj_CLT_all, threeSubj_CLT_nrows);

  std::vector<SPAsqr::TauFamilyNativeInput> tauData(ntaus);
  arma::uword twoOff = 0, threeOff = 0;
  for (int i = 0; i < ntaus; ++i) {
    if (residOutliers[i].n_elem > 0)
      tauData[i].resid_unrelated_outliers = std::move(residOutliers[i]);

    arma::uword nTwo = twoSubj_perTau[i];
    tauData[i].twoSubj_resid.assign(
      allTwoResid.begin() + twoOff, allTwoResid.begin() + twoOff + nTwo);
    tauData[i].twoSubj_rho.assign(
      allTwoRho.begin() + twoOff, allTwoRho.begin() + twoOff + nTwo);
    twoOff += nTwo;

    arma::uword nThree = threeSubj_perTau[i];
    tauData[i].threeSubj_standS.assign(
      allThreeStandS.begin() + threeOff, allThreeStandS.begin() + threeOff + nThree);
    tauData[i].threeSubj_CLT.assign(
      allThreeCLT.begin() + threeOff, allThreeCLT.begin() + threeOff + nThree);
    threeOff += nThree;
  }

  if (ptr_gSPAsqrobj) {
    delete ptr_gSPAsqrobj;
  }
  ptr_gSPAsqrobj = new SPAsqr::SPAsqrClass(
    taus, Resid_mat, tauData,
    sum_R_nonOutlier_vec, R_GRM_R_nonOutlier_vec,
    R_GRM_R_TwoSubjOutlier_vec, R_GRM_R_vec,
    MAF_interval, SPA_Cutoff, zeta, tol
  );

  runEngineCore(
    "SPAsqr", bedFile, bimFile, famFile, outputFile,
    subjData, AlleleOrder, nMarkersEachChunk, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker,
    IDsToIncludeFile, RangesToIncludeFile,
    IDsToExcludeFile, RangesToExcludeFile
  );
}

// [[Rcpp::export("runMarkerInCPP.SPACox")]]
void runMarkerInCPP_SPACox(
    arma::mat cumul, arma::vec mresid, arma::mat XinvXX, arma::mat tX,
    int N, double pVal_covaAdj_Cutoff, double SPA_Cutoff,
    std::string bedFile, std::string bimFile, std::string famFile,
    std::string outputFile,
    std::vector<std::string> subjData, std::string AlleleOrder,
    int nMarkersEachChunk, unsigned int nthreads,
    std::string impute_method, double missing_cutoff,
    double min_maf_marker, double min_mac_marker,
    std::string IDsToIncludeFile, std::string RangesToIncludeFile,
    std::string IDsToExcludeFile, std::string RangesToExcludeFile
) {
  if (ptr_gSPACoxobj) {
    delete ptr_gSPACoxobj;
  }
  ptr_gSPACoxobj = new SPACox::SPACoxClass(
    cumul, mresid, XinvXX, tX, N,
    pVal_covaAdj_Cutoff, SPA_Cutoff
  );

  runEngineCore(
    "SPACox", bedFile, bimFile, famFile, outputFile,
    subjData, AlleleOrder, nMarkersEachChunk, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker,
    IDsToIncludeFile, RangesToIncludeFile,
    IDsToExcludeFile, RangesToExcludeFile
  );
}

// [[Rcpp::export("runMarkerInCPP.SPAmix")]]
void runMarkerInCPP_SPAmix(
  arma::mat resid, arma::mat PCs, int N, double SPA_Cutoff, int nPheno,
  arma::uvec posValue_all, arma::uvec posValue_lens,
  arma::uvec posOutlier_all, arma::uvec posOutlier_lens,
  arma::uvec posNonOutlier_all, arma::uvec posNonOutlier_lens,
  std::string bedFile, std::string bimFile, std::string famFile,
  std::string outputFile,
  std::vector<std::string> subjData, std::string AlleleOrder,
  int nMarkersEachChunk, unsigned int nthreads,
  std::string impute_method, double missing_cutoff,
  double min_maf_marker, double min_mac_marker,
  std::string IDsToIncludeFile, std::string RangesToIncludeFile,
  std::string IDsToExcludeFile, std::string RangesToExcludeFile
) {
  auto posVals = splitUvec(posValue_all, posValue_lens);
  auto posOuts = splitUvec(posOutlier_all, posOutlier_lens);
  auto posNons = splitUvec(posNonOutlier_all, posNonOutlier_lens);
  std::vector<SPAmix::SPAmixClass::OutlierData> outlierVec(nPheno);

  for (int i = 0; i < nPheno; ++i) {
    outlierVec[i].posValue         = posVals[i];
    outlierVec[i].posOutlier       = posOuts[i];
    outlierVec[i].posNonOutlier    = posNons[i];
    arma::vec col_i                = resid.col(i);
    outlierVec[i].resid            = col_i.elem(posVals[i]);
    outlierVec[i].resid2           = arma::square(col_i.elem(posVals[i]));
    outlierVec[i].residOutlier     = col_i.elem(posOuts[i]);
    outlierVec[i].residNonOutlier  = col_i.elem(posNons[i]);
    outlierVec[i].resid2NonOutlier = arma::square(col_i.elem(posNons[i]));
  }

  if (ptr_gSPAmixobj) {
    delete ptr_gSPAmixobj;
  }
  ptr_gSPAmixobj = new SPAmix::SPAmixClass(
    resid, PCs, N, SPA_Cutoff, std::move(outlierVec)
  );

  runEngineCore(
    "SPAmix", bedFile, bimFile, famFile, outputFile,
    subjData, AlleleOrder, nMarkersEachChunk, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker,
    IDsToIncludeFile, RangesToIncludeFile,
    IDsToExcludeFile, RangesToExcludeFile
  );
}

// [[Rcpp::export("runMarkerInCPP.SPAmixPlus")]]
void runMarkerInCPP_SPAmixPlus(
  arma::mat resid, arma::mat PCs, int N, double SPA_Cutoff, int nPheno,
  arma::uvec posValue_all, arma::uvec posValue_lens,
  arma::uvec posOutlier_all, arma::uvec posOutlier_lens,
  arma::uvec posNonOutlier_all, arma::uvec posNonOutlier_lens,
  arma::ivec sparseId1, arma::ivec sparseId2, arma::vec sparseVal,
  std::string afFilePath, std::string afFilePrecision,
  std::string bedFile, std::string bimFile, std::string famFile,
  std::string outputFile,
  std::vector<std::string> subjData, std::string AlleleOrder,
  int nMarkersEachChunk, unsigned int nthreads,
  std::string impute_method, double missing_cutoff,
  double min_maf_marker, double min_mac_marker,
  std::string IDsToIncludeFile, std::string RangesToIncludeFile,
  std::string IDsToExcludeFile, std::string RangesToExcludeFile
) {
  auto posVals = splitUvec(posValue_all, posValue_lens);
  auto posOuts = splitUvec(posOutlier_all, posOutlier_lens);
  auto posNons = splitUvec(posNonOutlier_all, posNonOutlier_lens);
  std::vector<SPAmixPlus::SPAmixPlusClass::PhenoOutlierData> outlierVec(nPheno);

  for (int i = 0; i < nPheno; ++i) {
    outlierVec[i].posValue         = posVals[i];
    outlierVec[i].posOutlier       = posOuts[i];
    outlierVec[i].posNonOutlier    = posNons[i];
    arma::vec col_i                = resid.col(i);
    outlierVec[i].resid            = col_i.elem(posVals[i]);
    outlierVec[i].residOutlier     = col_i.elem(posOuts[i]);
    outlierVec[i].residNonOutlier  = col_i.elem(posNons[i]);
    outlierVec[i].resid2NonOutlier = arma::square(col_i.elem(posNons[i]));
  }
  std::vector<std::tuple<int, int, double>> triplets;
  triplets.reserve(sparseId1.n_elem);
  for (arma::uword i = 0; i < sparseId1.n_elem; ++i)
    triplets.emplace_back(sparseId1[i], sparseId2[i], sparseVal[i]);

  if (ptr_gSPAmixPlusobj) {
    delete ptr_gSPAmixPlusobj;
  }
  ptr_gSPAmixPlusobj = new SPAmixPlus::SPAmixPlusClass(
    resid, PCs, N, SPA_Cutoff,
    std::move(outlierVec), std::move(triplets),
    afFilePath, afFilePrecision
  );

  runEngineCore(
    "SPAmixPlus", bedFile, bimFile, famFile, outputFile,
    subjData, AlleleOrder, nMarkersEachChunk, nthreads,
    impute_method, missing_cutoff, min_maf_marker, min_mac_marker,
    IDsToIncludeFile, RangesToIncludeFile,
    IDsToExcludeFile, RangesToExcludeFile
  );
}
