// mtMarkerBridge.cpp -- Unified R/Rcpp boundary for marker-level association

#include <RcppArmadillo.h>
#include "mtPLINK.h"
#include "mtPOLMM.h"
#include "mtWtCoxG.h"
#include "mtLEAF.h"
#include "mtSPAGRM.h"
#include "mtSAGELD.h"
#include "mtSPAsqr.h"
#include "mtSPACox.h"
#include "mtSPAmix.h"
#include "mtSPAmixPlus.h"

// Forward declarations — defined in mtMarkerEngine.cpp.
extern std::unique_ptr<mtPOLMMClass>      ptr_gPOLMMobj;
extern std::unique_ptr<mtWtCoxGClass>     ptr_gWtCoxGobj;
extern std::unique_ptr<mtLEAFClass>       ptr_gLEAFobj;
extern std::unique_ptr<mtSPACoxClass>     ptr_gSPACoxobj;
extern std::unique_ptr<mtSPAmixClass>     ptr_gSPAmixobj;
extern std::unique_ptr<mtSPAGRMClass>     ptr_gSPAGRMobj;
extern std::unique_ptr<mtSAGELDClass>     ptr_gSAGELDobj;
extern std::unique_ptr<mtSPAsqrClass>     ptr_gSPAsqrobj;
extern std::unique_ptr<mtSPAmixPlusClass> ptr_gSPAmixPlusobj;
extern std::unique_ptr<const PlinkData>   ptr_gPlinkDataObj;

void mtMarkerEngine(
  const std::string method,
  const std::string outputFile,
  const int nthreads,
  const std::string impute_method,
  const double missing_cutoff,
  const double min_maf_marker,
  const double min_mac_marker,
  const bool exactHwe
);


namespace {

// ---- Helper: extract string from R list, return "" if NULL or missing ----

std::string strOrEmpty(const Rcpp::List& lst, const char* key) {
  if (lst.containsElementNamed(key) && !Rf_isNull(lst[key]))
    return Rcpp::as<std::string>(lst[key]);
  return "";
}

// ---- Helper: build RefInfo map from an R data.frame ----

std::unordered_map<uint64_t, mtWtCoxGClass::RefInfo>
buildRefMap(const Rcpp::DataFrame& df) {
  Rcpp::NumericVector genoIndex  = df["genoIndex"];
  Rcpp::NumericVector AF_ref     = df["AF_ref"];
  Rcpp::NumericVector AN_ref     = df["AN_ref"];
  Rcpp::NumericVector TPR        = df["TPR"];
  Rcpp::NumericVector sigma2     = df["sigma2"];
  Rcpp::NumericVector pvalue_bat = df["pvalue_bat"];
  Rcpp::NumericVector w_ext      = df["w.ext"];
  Rcpp::NumericVector var_w0     = df["var.ratio.w0"];
  Rcpp::NumericVector var_int    = df["var.ratio.int"];
  Rcpp::NumericVector var_ext    = df["var.ratio.ext"];
  std::unordered_map<uint64_t, mtWtCoxGClass::RefInfo> refMap;
  refMap.reserve(genoIndex.size());
  for (R_xlen_t i = 0; i < genoIndex.size(); ++i) {
    refMap[static_cast<uint64_t>(genoIndex[i])] = mtWtCoxGClass::RefInfo{
      AF_ref[i], AN_ref[i], TPR[i], sigma2[i], pvalue_bat[i],
      w_ext[i], var_w0[i], var_int[i], var_ext[i]
    };
  }
  return refMap;
}

// ---- Per-method setup functions ----

void setupPOLMM(const Rcpp::List& obj, const Rcpp::List& ctl) {
  Rcpp::List locoList = obj["LOCOList"];
  Rcpp::List objCHR = locoList["LOCO=F"];
  ptr_gPOLMMobj.reset(new mtPOLMMClass(
    Rcpp::as<arma::mat>(objCHR["muMat"]),
    Rcpp::as<arma::mat>(objCHR["iRMat"]),
    Rcpp::as<arma::mat>(obj["Cova"]),
    Rcpp::as<arma::uvec>(obj["yVec"]),
    Rcpp::as<double>(objCHR["VarRatio"]),
    Rcpp::as<double>(ctl["SPA_Cutoff"])
  ));
}

void setupWtCoxG(const Rcpp::List& obj, const Rcpp::List& ctl) {
  auto refMap = buildRefMap(Rcpp::as<Rcpp::DataFrame>(obj["mergeGenoInfo"]));
  ptr_gWtCoxGobj.reset(new mtWtCoxGClass(
    Rcpp::as<arma::vec>(obj["mresid"]),
    Rcpp::as<arma::vec>(obj["weight"]),
    Rcpp::as<double>(ctl["cutoff"]),
    Rcpp::as<double>(ctl["SPA_Cutoff"]),
    std::move(refMap)
  ));
}

void setupLEAF(const Rcpp::List& obj, const Rcpp::List& ctl) {
  int nCluster = Rcpp::as<int>(obj["Ncluster"]);
  arma::ivec clusterIdx = Rcpp::as<arma::ivec>(obj["clusterIdx"]);

  // Build per-cluster 0-based index vectors directly
  std::vector<arma::uvec> clusterIdxVecs(nCluster);
  for (int k = 0; k < nCluster; ++k)
    clusterIdxVecs[k] = arma::find(clusterIdx == (k + 1));

  Rcpp::List resList = obj["residuals_list"];
  Rcpp::List wgtList = obj["weights_list"];
  std::vector<arma::vec> residuals(nCluster), weights(nCluster);
  for (int k = 0; k < nCluster; ++k) {
    residuals[k] = Rcpp::as<arma::vec>(resList[k]);
    weights[k]   = Rcpp::as<arma::vec>(wgtList[k]);
  }

  // Build per-cluster RefMaps from subGenoInfo data frames
  Rcpp::List sgi = obj["subGenoInfo"];
  std::vector<std::shared_ptr<
    const std::unordered_map<uint64_t, mtWtCoxGClass::RefInfo>>> refMaps(nCluster);
  for (int c = 0; c < nCluster; ++c)
    refMaps[c] = std::make_shared<
      const std::unordered_map<uint64_t, mtWtCoxGClass::RefInfo>>(
        buildRefMap(Rcpp::as<Rcpp::DataFrame>(sgi[c])));

  ptr_gLEAFobj.reset(new mtLEAFClass(
    std::move(residuals), std::move(weights), std::move(clusterIdxVecs),
    Rcpp::as<double>(ctl["cutoff"]),
    Rcpp::as<double>(ctl["SPA_Cutoff"]),
    std::move(refMaps)
  ));
}

void setupSPAGRM(const Rcpp::List& obj, const Rcpp::List& ctl) {
  Rcpp::List twoList = obj["TwoSubj_list"];
  int nTwo = twoList.size();
  std::vector<std::array<double, 2>> twoResid(nTwo);
  std::vector<std::vector<double>> twoRho(nTwo);
  for (int i = 0; i < nTwo; ++i) {
    Rcpp::List item = twoList[i];
    Rcpp::NumericVector tmpR = item["Resid"];
    twoResid[i] = {tmpR[0], tmpR[1]};
    Rcpp::NumericVector tmpRho = item["Rho"];
    twoRho[i].assign(tmpRho.begin(), tmpRho.end());
  }

  Rcpp::List threeList = obj["ThreeSubj_list"];
  int nThree = threeList.size();
  std::vector<std::vector<double>> threeStandS(nThree);
  std::vector<arma::mat> threeCLT(nThree);
  for (int i = 0; i < nThree; ++i) {
    Rcpp::List item = threeList[i];
    Rcpp::NumericVector tmpSS = item["stand.S"];
    threeStandS[i].assign(tmpSS.begin(), tmpSS.end());
    threeCLT[i]    = Rcpp::as<arma::mat>(item["CLT"]);
  }

  nsSPAGRM::FamilyData fam{
    Rcpp::as<arma::vec>(obj["Resid.unrelated.outliers"]),
    std::move(twoResid), std::move(twoRho),
    std::move(threeStandS), std::move(threeCLT)
  };

  ptr_gSPAGRMobj.reset(new mtSPAGRMClass(
    Rcpp::as<arma::vec>(obj["Resid"]),
    Rcpp::as<double>(obj["sum_R_nonOutlier"]),
    Rcpp::as<double>(obj["R_GRM_R_nonOutlier"]),
    Rcpp::as<double>(obj["R_GRM_R_TwoSubjOutlier"]),
    Rcpp::as<double>(obj["R_GRM_R"]),
    Rcpp::as<std::vector<double>>(obj["MAF_interval"]),
    std::move(fam),
    Rcpp::as<double>(ctl["SPA_Cutoff"]),
    Rcpp::as<double>(ctl["zeta"]),
    Rcpp::as<double>(ctl["tol"])
  ));
}

void setupSAGELD(const Rcpp::List& obj, const Rcpp::List& ctl) {
  Rcpp::List twoList = obj["TwoSubj_list"];
  int nTwo = twoList.size();
  std::vector<mtSAGELDClass::TwoSubjFamily> twoSubj(nTwo);
  for (int i = 0; i < nTwo; ++i) {
    Rcpp::List item = twoList[i];
    {
      Rcpp::NumericVector tmp = item["Resid"];
      twoSubj[i].Resid = {tmp[0], tmp[1]};
    }
    {
      Rcpp::NumericVector tmp = item["Rho"];
      twoSubj[i].Rho.assign(tmp.begin(), tmp.end());
    }
    {
      Rcpp::NumericVector tmp = item["Resid_G"];
      twoSubj[i].Resid_G = {tmp[0], tmp[1]};
    }
    {
      Rcpp::NumericVector tmp = item["Resid_GxE"];
      twoSubj[i].Resid_GxE = {tmp[0], tmp[1]};
    }
  }

  Rcpp::List threeList = obj["ThreeSubj_list"];
  int nThree = threeList.size();
  std::vector<mtSAGELDClass::ThreeSubjFamily> threeSubj(nThree);
  for (int i = 0; i < nThree; ++i) {
    Rcpp::List item = threeList[i];
    threeSubj[i].CLT         = Rcpp::as<arma::mat>(item["CLT"]);
    {
      Rcpp::NumericVector tmp = item["stand.S"];
      threeSubj[i].stand_S.assign(tmp.begin(), tmp.end());
    }
    {
      Rcpp::NumericVector tmp = item["stand.S_G"];
      threeSubj[i].stand_S_G.assign(tmp.begin(), tmp.end());
    }
    {
      Rcpp::NumericVector tmp = item["stand.S_GxE"];
      threeSubj[i].stand_S_GxE.assign(tmp.begin(), tmp.end());
    }
  }

  // Pass all R_GRM_R values directly as scalars — no temporary arma::vec needed.
  ptr_gSAGELDobj.reset(new mtSAGELDClass(
    Rcpp::as<std::string>(obj["Method"]),
    Rcpp::as<arma::mat>(obj["XTs"]),
    Rcpp::as<arma::mat>(obj["SS"]),
    Rcpp::as<arma::mat>(obj["AtS"]),
    Rcpp::as<arma::mat>(obj["Q"]),
    Rcpp::as<arma::mat>(obj["A21"]),
    Rcpp::as<arma::mat>(obj["TTs"]),
    Rcpp::as<arma::mat>(obj["Tys"]),
    Rcpp::as<arma::vec>(obj["sol"]),
    Rcpp::as<arma::vec>(obj["blups"]),
    Rcpp::as<double>(obj["sig"]),
    Rcpp::as<arma::vec>(obj["Resid"]),
    Rcpp::as<arma::vec>(obj["Resid_G"]),
    Rcpp::as<arma::vec>(obj["Resid_GxE"]),
    Rcpp::as<arma::vec>(obj["Resid_E"]),
    Rcpp::as<arma::vec>(obj["Resid.unrelated.outliers"]),
    Rcpp::as<arma::vec>(obj["Resid.unrelated.outliers_G"]),
    Rcpp::as<arma::vec>(obj["Resid.unrelated.outliers_GxE"]),
    Rcpp::as<double>(obj["sum_R_nonOutlier"]),
    Rcpp::as<double>(obj["sum_R_nonOutlier_G"]),
    Rcpp::as<double>(obj["sum_R_nonOutlier_GxE"]),
    Rcpp::as<double>(obj["R_GRM_R"]),
    Rcpp::as<double>(obj["R_GRM_R_G"]),
    Rcpp::as<double>(obj["R_GRM_R_GxE"]),
    Rcpp::as<double>(obj["R_GRM_R_G_GxE"]),
    Rcpp::as<double>(obj["R_GRM_R_E"]),
    Rcpp::as<double>(obj["R_GRM_R_nonOutlier"]),
    Rcpp::as<double>(obj["R_GRM_R_nonOutlier_G"]),
    Rcpp::as<double>(obj["R_GRM_R_nonOutlier_GxE"]),
    Rcpp::as<double>(obj["R_GRM_R_nonOutlier_G_GxE"]),
    Rcpp::as<double>(obj["R_GRM_R_TwoSubjOutlier"]),
    Rcpp::as<double>(obj["R_GRM_R_TwoSubjOutlier_G"]),
    Rcpp::as<double>(obj["R_GRM_R_TwoSubjOutlier_GxE"]),
    Rcpp::as<double>(obj["R_GRM_R_TwoSubjOutlier_G_GxE"]),
    std::move(twoSubj), std::move(threeSubj),
    Rcpp::as<std::vector<double>>(obj["MAF_interval"]),
    Rcpp::as<double>(obj["zScoreE_cutoff"]),
    Rcpp::as<double>(ctl["SPA_Cutoff"]),
    Rcpp::as<double>(ctl["zeta"]),
    Rcpp::as<double>(ctl["tol"])
  ));
}

void setupSPAsqr(const Rcpp::List& obj, const Rcpp::List& ctl) {
  arma::vec taus = Rcpp::as<arma::vec>(obj["taus"]);
  int ntaus = static_cast<int>(taus.n_elem);

  // Wrap without copying; extract each column once directly from R memory.
  Rcpp::NumericMatrix rmat(obj["Resid_mat"]);
  int N = rmat.nrow();
  std::vector<arma::vec> residCols(ntaus);
  for (int i = 0; i < ntaus; ++i)
    residCols[i] = arma::vec(rmat.column(i).begin(), N);

  Rcpp::List residOutlierLst     = obj["Resid.unrelated.outliers_lst"];
  Rcpp::List twoSubjLstLst       = obj["TwoSubj_list_lst"];
  Rcpp::List cltUnionLst         = obj["CLT_union_lst"];
  Rcpp::List threeSubjFamIdxLst  = obj["ThreeSubj_family_idx_lst"];
  Rcpp::List threeSubjStandSLst  = obj["ThreeSubj_stand_S_lst"];

  std::vector<nsSPAGRM::FamilyData> tauData(ntaus);
  for (int i = 0; i < ntaus; ++i) {
    // Outlier residuals (may be NULL for some taus)
    SEXP outlierSexp = residOutlierLst[i];
    if (!Rf_isNull(outlierSexp))
      tauData[i].resid_unrelated_outliers = Rcpp::as<arma::vec>(outlierSexp);

    // Two-subject families for this tau
    Rcpp::List pairs = twoSubjLstLst[i];
    int nTwo = pairs.size();
    tauData[i].twoSubj_resid.resize(nTwo);
    tauData[i].twoSubj_rho.resize(nTwo);
    for (int j = 0; j < nTwo; ++j) {
      Rcpp::List p = pairs[j];
      tauData[i].twoSubj_resid[j] = [&]{
        Rcpp::NumericVector tmp = p["Resid"];
        return std::array<double, 2>{tmp[0], tmp[1]};
      }();
      {
        Rcpp::NumericVector tmp = p["Rho"];
        tauData[i].twoSubj_rho[j].assign(tmp.begin(), tmp.end());
      }
    }

    // Three-subject families: standS from per-tau list, CLT from shared union
    Rcpp::IntegerVector famIdx = Rcpp::as<Rcpp::IntegerVector>(threeSubjFamIdxLst[i]);
    Rcpp::List standSList = threeSubjStandSLst[i];
    int nThree = famIdx.size();
    tauData[i].threeSubj_standS.resize(nThree);
    tauData[i].threeSubj_CLT.resize(nThree);
    for (int j = 0; j < nThree; ++j) {
      {
        Rcpp::NumericVector tmp = standSList[j];
        tauData[i].threeSubj_standS[j].assign(tmp.begin(), tmp.end());
      }
      tauData[i].threeSubj_CLT[j]    = Rcpp::as<arma::mat>(cltUnionLst[famIdx[j] - 1]);
    }
  }

  // Zero-copy wrappers for the per-tau scalar vectors.
  Rcpp::NumericVector rvSum  = obj["sum_R_nonOutlier_vec"];
  Rcpp::NumericVector rvNonO = obj["R_GRM_R_nonOutlier_vec"];
  Rcpp::NumericVector rvTwoS = obj["R_GRM_R_TwoSubjOutlier_vec"];
  Rcpp::NumericVector rvAll  = obj["R_GRM_R_vec"];

  ptr_gSPAsqrobj.reset(new mtSPAsqrClass(
    std::move(taus),
    std::move(residCols),
    std::move(tauData),
    std::vector<double>(rvSum.begin(),  rvSum.end()),
    std::vector<double>(rvNonO.begin(), rvNonO.end()),
    std::vector<double>(rvTwoS.begin(), rvTwoS.end()),
    std::vector<double>(rvAll.begin(),  rvAll.end()),
    Rcpp::as<std::vector<double>>(obj["MAF_interval"]),
    Rcpp::as<double>(ctl["SPA_Cutoff"]),
    Rcpp::as<double>(ctl["zeta"]),
    Rcpp::as<double>(ctl["tol"])
  ));
}

void setupSPACox(const Rcpp::List& obj, const Rcpp::List& ctl) {
  arma::vec mresid = Rcpp::as<arma::vec>(obj["mresid"]);
  int N = static_cast<int>(mresid.n_elem);
  ptr_gSPACoxobj.reset(new mtSPACoxClass(
    Rcpp::as<arma::mat>(obj["cumul"]),
    std::move(mresid),
    Rcpp::as<arma::mat>(obj["X.invXX"]),
    Rcpp::as<arma::mat>(obj["tX"]),
    N,
    Rcpp::as<double>(ctl["pVal_covaAdj_Cutoff"]),
    Rcpp::as<double>(ctl["SPA_Cutoff"])
  ));
}

void setupSPAmix(const Rcpp::List& obj, const Rcpp::List& ctl) {
  arma::vec resid = Rcpp::as<arma::vec>(obj["resid"]);

  Rcpp::List ol = Rcpp::as<Rcpp::List>(obj["outLierList"]);
  Rcpp::List ph = ol[0];
  arma::uvec posVal = Rcpp::as<arma::uvec>(ph["posValue"]);
  arma::uvec posOut = Rcpp::as<arma::uvec>(ph["posOutlier"]);
  arma::uvec posNon = Rcpp::as<arma::uvec>(ph["posNonOutlier"]);

  mtSPAmixClass::OutlierData outlier;
  outlier.posValue         = posVal;
  outlier.posOutlier       = posOut;
  outlier.posNonOutlier    = posNon;
  outlier.resid            = resid.elem(posVal);
  outlier.resid2           = arma::square(outlier.resid);
  outlier.residOutlier     = resid.elem(posOut);
  outlier.residNonOutlier  = resid.elem(posNon);
  outlier.resid2NonOutlier = arma::square(outlier.residNonOutlier);

  ptr_gSPAmixobj.reset(new mtSPAmixClass(
    std::move(resid),
    Rcpp::as<arma::mat>(obj["PCs"]),
    Rcpp::as<int>(obj["N"]),
    Rcpp::as<double>(ctl["SPA_Cutoff"]),
    std::move(outlier)
  ));
}

void setupSPAmixPlus(const Rcpp::List& obj, const Rcpp::List& ctl) {
  arma::vec resid = Rcpp::as<arma::vec>(obj["resid"]);

  Rcpp::List ol = Rcpp::as<Rcpp::List>(obj["outLierList"]);
  Rcpp::List ph = ol[0];
  arma::uvec posVal = Rcpp::as<arma::uvec>(ph["posValue"]);
  arma::uvec posOut = Rcpp::as<arma::uvec>(ph["posOutlier"]);
  arma::uvec posNon = Rcpp::as<arma::uvec>(ph["posNonOutlier"]);

  mtSPAmixPlusClass::OutlierData outlier;
  outlier.posValue         = posVal;
  outlier.posOutlier       = posOut;
  outlier.posNonOutlier    = posNon;
  outlier.resid            = resid.elem(posVal);
  outlier.residOutlier     = resid.elem(posOut);
  outlier.residNonOutlier  = resid.elem(posNon);
  outlier.resid2NonOutlier = arma::square(outlier.residNonOutlier);

  Rcpp::List sparseGRM = obj["sparseGRM"];
  Rcpp::IntegerVector id1 = sparseGRM["id1_index"];
  Rcpp::IntegerVector id2 = sparseGRM["id2_index"];
  Rcpp::NumericVector val = sparseGRM["value"];
  std::vector<std::tuple<int, int, double>> triplets;
  triplets.reserve(id1.size());
  for (R_xlen_t i = 0; i < id1.size(); ++i)
    triplets.emplace_back(id1[i], id2[i], val[i]);

  ptr_gSPAmixPlusobj.reset(new mtSPAmixPlusClass(
    std::move(resid),
    Rcpp::as<arma::mat>(obj["PCs"]),
    Rcpp::as<int>(obj["N"]),
    Rcpp::as<double>(ctl["SPA_Cutoff"]),
    std::move(outlier), std::move(triplets),
    strOrEmpty(ctl, "afFilePath"),
    Rcpp::as<std::string>(ctl["afFilePrecision"])
  ));
}

} // namespace


// [[Rcpp::export("mtMarkerBridgeInCPP")]]
void mtMarkerBridgeInCPP(
    const Rcpp::List& objNull,
    const std::string& OutputFile,
    const Rcpp::List& control,
    const std::string& bedFile,
    const std::string& bimFile,
    const std::string& famFile
) {
  std::string nullClass = Rcpp::as<std::string>(objNull.attr("class"));
  std::string method;

  if (nullClass == "POLMM_NULL_Model") {
    method = "POLMM";
    setupPOLMM(objNull, control);
  } else if (nullClass == "WtCoxG_NULL_Model") {
    method = "WtCoxG";
    setupWtCoxG(objNull, control);
  } else if (nullClass == "LEAF_NULL_Model") {
    method = "LEAF";
    setupLEAF(objNull, control);
  } else if (nullClass == "SPAGRM_NULL_Model") {
    method = "SPAGRM";
    setupSPAGRM(objNull, control);
  } else if (nullClass == "SAGELD_NULL_Model") {
    method = "SAGELD";
    setupSAGELD(objNull, control);
  } else if (nullClass == "SPAsqr_NULL_Model") {
    method = "SPAsqr";
    setupSPAsqr(objNull, control);
  } else if (nullClass == "SPACox_NULL_Model") {
    method = "SPACox";
    setupSPACox(objNull, control);
  } else if (nullClass == "SPAmix_NULL_Model") {
    method = "SPAmix";
    setupSPAmix(objNull, control);
  } else if (nullClass == "SPAmixPlus_NULL_Model") {
    method = "SPAmixPlus";
    setupSPAmixPlus(objNull, control);
  } else {
    Rcpp::stop("Unsupported null model class: " + nullClass);
  }

  std::vector<std::string> subjData =
    Rcpp::as<std::vector<std::string>>(objNull["subjData"]);

  // Create PlinkData — shared read-only across all worker threads
  ptr_gPlinkDataObj = std::make_unique<const PlinkData>(
    bedFile, bimFile, famFile, subjData,
    Rcpp::as<std::string>(control["AlleleOrder"]),
    strOrEmpty(control, "IDsToIncludeFile"),
    strOrEmpty(control, "RangesToIncludeFile"),
    strOrEmpty(control, "IDsToExcludeFile"),
    strOrEmpty(control, "RangesToExcludeFile"),
    Rcpp::as<int>(control["nMarkersEachChunk"])
  );
  // PlinkCursor is created per-worker by the ThreadContext constructor

  mtMarkerEngine(
    method, OutputFile,
    Rcpp::as<int>(control["nthreads"]),
    Rcpp::as<std::string>(control["impute_method"]),
    Rcpp::as<double>(control["missing_cutoff"]),
    Rcpp::as<double>(control["min_maf_marker"]),
    Rcpp::as<double>(control["min_mac_marker"]),
    Rcpp::as<std::string>(control["hwe"]) == "exact"
  );
}
