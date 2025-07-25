// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// setSparseGRMInCPP
void setSparseGRMInCPP(Rcpp::List t_KinMatListR);
RcppExport SEXP _GRAB_setSparseGRMInCPP(SEXP t_KinMatListRSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type t_KinMatListR(t_KinMatListRSEXP);
    setSparseGRMInCPP(t_KinMatListR);
    return R_NilValue;
END_RCPP
}
// setDenseGRMInCPP
void setDenseGRMInCPP(double t_memoryChunk, double t_minMafGRM, double t_maxMissingGRM);
RcppExport SEXP _GRAB_setDenseGRMInCPP(SEXP t_memoryChunkSEXP, SEXP t_minMafGRMSEXP, SEXP t_maxMissingGRMSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type t_memoryChunk(t_memoryChunkSEXP);
    Rcpp::traits::input_parameter< double >::type t_minMafGRM(t_minMafGRMSEXP);
    Rcpp::traits::input_parameter< double >::type t_maxMissingGRM(t_maxMissingGRMSEXP);
    setDenseGRMInCPP(t_memoryChunk, t_minMafGRM, t_maxMissingGRM);
    return R_NilValue;
END_RCPP
}
// getDenseGRMInCPP
arma::vec getDenseGRMInCPP(arma::vec t_bVec, std::string t_excludeChr, int t_grainSize);
RcppExport SEXP _GRAB_getDenseGRMInCPP(SEXP t_bVecSEXP, SEXP t_excludeChrSEXP, SEXP t_grainSizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t_bVec(t_bVecSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_excludeChr(t_excludeChrSEXP);
    Rcpp::traits::input_parameter< int >::type t_grainSize(t_grainSizeSEXP);
    rcpp_result_gen = Rcpp::wrap(getDenseGRMInCPP(t_bVec, t_excludeChr, t_grainSize));
    return rcpp_result_gen;
END_RCPP
}
// setMarker_GlobalVarsInCPP
void setMarker_GlobalVarsInCPP(std::string t_impute_method, double t_missing_cutoff, double t_min_maf_marker, double t_min_mac_marker, unsigned int t_omp_num_threads, arma::uvec t_group, bool t_ifOutGroup, unsigned int t_nGroup);
RcppExport SEXP _GRAB_setMarker_GlobalVarsInCPP(SEXP t_impute_methodSEXP, SEXP t_missing_cutoffSEXP, SEXP t_min_maf_markerSEXP, SEXP t_min_mac_markerSEXP, SEXP t_omp_num_threadsSEXP, SEXP t_groupSEXP, SEXP t_ifOutGroupSEXP, SEXP t_nGroupSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_impute_method(t_impute_methodSEXP);
    Rcpp::traits::input_parameter< double >::type t_missing_cutoff(t_missing_cutoffSEXP);
    Rcpp::traits::input_parameter< double >::type t_min_maf_marker(t_min_maf_markerSEXP);
    Rcpp::traits::input_parameter< double >::type t_min_mac_marker(t_min_mac_markerSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_omp_num_threads(t_omp_num_threadsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type t_group(t_groupSEXP);
    Rcpp::traits::input_parameter< bool >::type t_ifOutGroup(t_ifOutGroupSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_nGroup(t_nGroupSEXP);
    setMarker_GlobalVarsInCPP(t_impute_method, t_missing_cutoff, t_min_maf_marker, t_min_mac_marker, t_omp_num_threads, t_group, t_ifOutGroup, t_nGroup);
    return R_NilValue;
END_RCPP
}
// setRegion_GlobalVarsInCPP
void setRegion_GlobalVarsInCPP(std::string t_impute_method, double t_missing_cutoff, double t_max_maf_region, double t_min_mac_region, unsigned int t_max_markers_region, unsigned int t_omp_num_threads, arma::vec t_region_weight_beta, arma::vec t_region_max_maf_vec);
RcppExport SEXP _GRAB_setRegion_GlobalVarsInCPP(SEXP t_impute_methodSEXP, SEXP t_missing_cutoffSEXP, SEXP t_max_maf_regionSEXP, SEXP t_min_mac_regionSEXP, SEXP t_max_markers_regionSEXP, SEXP t_omp_num_threadsSEXP, SEXP t_region_weight_betaSEXP, SEXP t_region_max_maf_vecSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_impute_method(t_impute_methodSEXP);
    Rcpp::traits::input_parameter< double >::type t_missing_cutoff(t_missing_cutoffSEXP);
    Rcpp::traits::input_parameter< double >::type t_max_maf_region(t_max_maf_regionSEXP);
    Rcpp::traits::input_parameter< double >::type t_min_mac_region(t_min_mac_regionSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_max_markers_region(t_max_markers_regionSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_omp_num_threads(t_omp_num_threadsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_region_weight_beta(t_region_weight_betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_region_max_maf_vec(t_region_max_maf_vecSEXP);
    setRegion_GlobalVarsInCPP(t_impute_method, t_missing_cutoff, t_max_maf_region, t_min_mac_region, t_max_markers_region, t_omp_num_threads, t_region_weight_beta, t_region_max_maf_vec);
    return R_NilValue;
END_RCPP
}
// mainMarkerInCPP
Rcpp::List mainMarkerInCPP(std::string t_method, std::string t_genoType, std::vector<uint64_t> t_genoIndex);
RcppExport SEXP _GRAB_mainMarkerInCPP(SEXP t_methodSEXP, SEXP t_genoTypeSEXP, SEXP t_genoIndexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_method(t_methodSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_genoType(t_genoTypeSEXP);
    Rcpp::traits::input_parameter< std::vector<uint64_t> >::type t_genoIndex(t_genoIndexSEXP);
    rcpp_result_gen = Rcpp::wrap(mainMarkerInCPP(t_method, t_genoType, t_genoIndex));
    return rcpp_result_gen;
END_RCPP
}
// mainRegionURVInCPP
Rcpp::List mainRegionURVInCPP(std::string t_method, std::string t_genoType, std::vector<uint64_t> t_genoIndex, unsigned int t_n);
RcppExport SEXP _GRAB_mainRegionURVInCPP(SEXP t_methodSEXP, SEXP t_genoTypeSEXP, SEXP t_genoIndexSEXP, SEXP t_nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_method(t_methodSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_genoType(t_genoTypeSEXP);
    Rcpp::traits::input_parameter< std::vector<uint64_t> >::type t_genoIndex(t_genoIndexSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_n(t_nSEXP);
    rcpp_result_gen = Rcpp::wrap(mainRegionURVInCPP(t_method, t_genoType, t_genoIndex, t_n));
    return rcpp_result_gen;
END_RCPP
}
// mainRegionInCPP
Rcpp::List mainRegionInCPP(std::string t_method, std::string t_genoType, std::vector<uint64_t> t_genoIndex, std::vector<double> t_weightVec, std::string t_outputFile, std::vector<unsigned int> t_labelVec, unsigned int t_nLabel, arma::mat t_annoMat, std::vector<std::string> t_annoVec);
RcppExport SEXP _GRAB_mainRegionInCPP(SEXP t_methodSEXP, SEXP t_genoTypeSEXP, SEXP t_genoIndexSEXP, SEXP t_weightVecSEXP, SEXP t_outputFileSEXP, SEXP t_labelVecSEXP, SEXP t_nLabelSEXP, SEXP t_annoMatSEXP, SEXP t_annoVecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_method(t_methodSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_genoType(t_genoTypeSEXP);
    Rcpp::traits::input_parameter< std::vector<uint64_t> >::type t_genoIndex(t_genoIndexSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type t_weightVec(t_weightVecSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_outputFile(t_outputFileSEXP);
    Rcpp::traits::input_parameter< std::vector<unsigned int> >::type t_labelVec(t_labelVecSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type t_nLabel(t_nLabelSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_annoMat(t_annoMatSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type t_annoVec(t_annoVecSEXP);
    rcpp_result_gen = Rcpp::wrap(mainRegionInCPP(t_method, t_genoType, t_genoIndex, t_weightVec, t_outputFile, t_labelVec, t_nLabel, t_annoMat, t_annoVec));
    return rcpp_result_gen;
END_RCPP
}
// printTimeDiffInCPP
void printTimeDiffInCPP();
RcppExport SEXP _GRAB_printTimeDiffInCPP() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    printTimeDiffInCPP();
    return R_NilValue;
END_RCPP
}
// printTimeDiffSPAmixInCPP
void printTimeDiffSPAmixInCPP();
RcppExport SEXP _GRAB_printTimeDiffSPAmixInCPP() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    printTimeDiffSPAmixInCPP();
    return R_NilValue;
END_RCPP
}
// getGenoInfoInCPP
arma::mat getGenoInfoInCPP(std::string t_genoType, Rcpp::DataFrame t_markerInfo, std::string t_imputeMethod);
RcppExport SEXP _GRAB_getGenoInfoInCPP(SEXP t_genoTypeSEXP, SEXP t_markerInfoSEXP, SEXP t_imputeMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_genoType(t_genoTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type t_markerInfo(t_markerInfoSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_imputeMethod(t_imputeMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(getGenoInfoInCPP(t_genoType, t_markerInfo, t_imputeMethod));
    return rcpp_result_gen;
END_RCPP
}
// getGenoInCPP
arma::mat getGenoInCPP(std::string t_genoType, Rcpp::DataFrame t_markerInfo, int n, std::string t_imputeMethod);
RcppExport SEXP _GRAB_getGenoInCPP(SEXP t_genoTypeSEXP, SEXP t_markerInfoSEXP, SEXP nSEXP, SEXP t_imputeMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_genoType(t_genoTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type t_markerInfo(t_markerInfoSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_imputeMethod(t_imputeMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(getGenoInCPP(t_genoType, t_markerInfo, n, t_imputeMethod));
    return rcpp_result_gen;
END_RCPP
}
// getGenoInCPP_fixedNumber
arma::mat getGenoInCPP_fixedNumber(std::string t_genoType, Rcpp::DataFrame t_markerInfo, int n, std::string t_imputeMethod, int m, double missingRateCutoff, double minMAFCutoff);
RcppExport SEXP _GRAB_getGenoInCPP_fixedNumber(SEXP t_genoTypeSEXP, SEXP t_markerInfoSEXP, SEXP nSEXP, SEXP t_imputeMethodSEXP, SEXP mSEXP, SEXP missingRateCutoffSEXP, SEXP minMAFCutoffSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_genoType(t_genoTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type t_markerInfo(t_markerInfoSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_imputeMethod(t_imputeMethodSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< double >::type missingRateCutoff(missingRateCutoffSEXP);
    Rcpp::traits::input_parameter< double >::type minMAFCutoff(minMAFCutoffSEXP);
    rcpp_result_gen = Rcpp::wrap(getGenoInCPP_fixedNumber(t_genoType, t_markerInfo, n, t_imputeMethod, m, missingRateCutoff, minMAFCutoff));
    return rcpp_result_gen;
END_RCPP
}
// getSpGenoInCPP
arma::sp_mat getSpGenoInCPP(std::string t_genoType, Rcpp::DataFrame t_markerInfo, int n, std::string t_imputeMethod);
RcppExport SEXP _GRAB_getSpGenoInCPP(SEXP t_genoTypeSEXP, SEXP t_markerInfoSEXP, SEXP nSEXP, SEXP t_imputeMethodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_genoType(t_genoTypeSEXP);
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type t_markerInfo(t_markerInfoSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_imputeMethod(t_imputeMethodSEXP);
    rcpp_result_gen = Rcpp::wrap(getSpGenoInCPP(t_genoType, t_markerInfo, n, t_imputeMethod));
    return rcpp_result_gen;
END_RCPP
}
// setPLINKobjInCPP
void setPLINKobjInCPP(std::string t_bimFile, std::string t_famFile, std::string t_bedFile, std::vector<std::string> t_SampleInModel, std::string t_AlleleOrder);
RcppExport SEXP _GRAB_setPLINKobjInCPP(SEXP t_bimFileSEXP, SEXP t_famFileSEXP, SEXP t_bedFileSEXP, SEXP t_SampleInModelSEXP, SEXP t_AlleleOrderSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_bimFile(t_bimFileSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_famFile(t_famFileSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_bedFile(t_bedFileSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type t_SampleInModel(t_SampleInModelSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_AlleleOrder(t_AlleleOrderSEXP);
    setPLINKobjInCPP(t_bimFile, t_famFile, t_bedFile, t_SampleInModel, t_AlleleOrder);
    return R_NilValue;
END_RCPP
}
// setBGENobjInCPP
void setBGENobjInCPP(std::string t_bgenFileName, std::string t_bgenFileIndex, std::vector<std::string> t_SampleInBgen, std::vector<std::string> t_SampleInModel, bool t_isSparseDosageInBgen, bool t_isDropmissingdosagesInBgen, std::string t_AlleleOrder);
RcppExport SEXP _GRAB_setBGENobjInCPP(SEXP t_bgenFileNameSEXP, SEXP t_bgenFileIndexSEXP, SEXP t_SampleInBgenSEXP, SEXP t_SampleInModelSEXP, SEXP t_isSparseDosageInBgenSEXP, SEXP t_isDropmissingdosagesInBgenSEXP, SEXP t_AlleleOrderSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_bgenFileName(t_bgenFileNameSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_bgenFileIndex(t_bgenFileIndexSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type t_SampleInBgen(t_SampleInBgenSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type t_SampleInModel(t_SampleInModelSEXP);
    Rcpp::traits::input_parameter< bool >::type t_isSparseDosageInBgen(t_isSparseDosageInBgenSEXP);
    Rcpp::traits::input_parameter< bool >::type t_isDropmissingdosagesInBgen(t_isDropmissingdosagesInBgenSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_AlleleOrder(t_AlleleOrderSEXP);
    setBGENobjInCPP(t_bgenFileName, t_bgenFileIndex, t_SampleInBgen, t_SampleInModel, t_isSparseDosageInBgen, t_isDropmissingdosagesInBgen, t_AlleleOrder);
    return R_NilValue;
END_RCPP
}
// closeGenoInputInCPP
void closeGenoInputInCPP(std::string t_genoType);
RcppExport SEXP _GRAB_closeGenoInputInCPP(SEXP t_genoTypeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_genoType(t_genoTypeSEXP);
    closeGenoInputInCPP(t_genoType);
    return R_NilValue;
END_RCPP
}
// setPOLMMobjInCPP
void setPOLMMobjInCPP(arma::mat t_muMat, arma::mat t_iRMat, arma::mat t_Cova, arma::uvec t_yVec, double t_tau, bool t_printPCGInfo, double t_tolPCG, int t_maxiterPCG, double t_varRatio, double t_SPA_cutoff, bool t_flagSparseGRM);
RcppExport SEXP _GRAB_setPOLMMobjInCPP(SEXP t_muMatSEXP, SEXP t_iRMatSEXP, SEXP t_CovaSEXP, SEXP t_yVecSEXP, SEXP t_tauSEXP, SEXP t_printPCGInfoSEXP, SEXP t_tolPCGSEXP, SEXP t_maxiterPCGSEXP, SEXP t_varRatioSEXP, SEXP t_SPA_cutoffSEXP, SEXP t_flagSparseGRMSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type t_muMat(t_muMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_iRMat(t_iRMatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_Cova(t_CovaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type t_yVec(t_yVecSEXP);
    Rcpp::traits::input_parameter< double >::type t_tau(t_tauSEXP);
    Rcpp::traits::input_parameter< bool >::type t_printPCGInfo(t_printPCGInfoSEXP);
    Rcpp::traits::input_parameter< double >::type t_tolPCG(t_tolPCGSEXP);
    Rcpp::traits::input_parameter< int >::type t_maxiterPCG(t_maxiterPCGSEXP);
    Rcpp::traits::input_parameter< double >::type t_varRatio(t_varRatioSEXP);
    Rcpp::traits::input_parameter< double >::type t_SPA_cutoff(t_SPA_cutoffSEXP);
    Rcpp::traits::input_parameter< bool >::type t_flagSparseGRM(t_flagSparseGRMSEXP);
    setPOLMMobjInCPP(t_muMat, t_iRMat, t_Cova, t_yVec, t_tau, t_printPCGInfo, t_tolPCG, t_maxiterPCG, t_varRatio, t_SPA_cutoff, t_flagSparseGRM);
    return R_NilValue;
END_RCPP
}
// setPOLMMobjInCPP_NULL
Rcpp::List setPOLMMobjInCPP_NULL(bool t_flagSparseGRM, arma::mat t_Cova, arma::uvec t_yVec, arma::vec t_beta, arma::vec t_bVec, arma::vec t_eps, double t_tau, Rcpp::List t_SPmatR, Rcpp::List t_controlList, arma::mat GenoMat);
RcppExport SEXP _GRAB_setPOLMMobjInCPP_NULL(SEXP t_flagSparseGRMSEXP, SEXP t_CovaSEXP, SEXP t_yVecSEXP, SEXP t_betaSEXP, SEXP t_bVecSEXP, SEXP t_epsSEXP, SEXP t_tauSEXP, SEXP t_SPmatRSEXP, SEXP t_controlListSEXP, SEXP GenoMatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< bool >::type t_flagSparseGRM(t_flagSparseGRMSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_Cova(t_CovaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type t_yVec(t_yVecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_beta(t_betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_bVec(t_bVecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_eps(t_epsSEXP);
    Rcpp::traits::input_parameter< double >::type t_tau(t_tauSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type t_SPmatR(t_SPmatRSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type t_controlList(t_controlListSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type GenoMat(GenoMatSEXP);
    rcpp_result_gen = Rcpp::wrap(setPOLMMobjInCPP_NULL(t_flagSparseGRM, t_Cova, t_yVec, t_beta, t_bVec, t_eps, t_tau, t_SPmatR, t_controlList, GenoMat));
    return rcpp_result_gen;
END_RCPP
}
// setSPAGRMobjInCPP
void setSPAGRMobjInCPP(arma::vec t_resid, arma::vec t_resid_unrelated_outliers, double t_sum_R_nonOutlier, double t_R_GRM_R_nonOutlier, double t_R_GRM_R_TwoSubjOutlier, double t_R_GRM_R, arma::vec t_MAF_interval, Rcpp::List t_TwoSubj_list, Rcpp::List t_ThreeSubj_list, double t_SPA_Cutoff, double t_zeta, double t_tol);
RcppExport SEXP _GRAB_setSPAGRMobjInCPP(SEXP t_residSEXP, SEXP t_resid_unrelated_outliersSEXP, SEXP t_sum_R_nonOutlierSEXP, SEXP t_R_GRM_R_nonOutlierSEXP, SEXP t_R_GRM_R_TwoSubjOutlierSEXP, SEXP t_R_GRM_RSEXP, SEXP t_MAF_intervalSEXP, SEXP t_TwoSubj_listSEXP, SEXP t_ThreeSubj_listSEXP, SEXP t_SPA_CutoffSEXP, SEXP t_zetaSEXP, SEXP t_tolSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t_resid(t_residSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_resid_unrelated_outliers(t_resid_unrelated_outliersSEXP);
    Rcpp::traits::input_parameter< double >::type t_sum_R_nonOutlier(t_sum_R_nonOutlierSEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_nonOutlier(t_R_GRM_R_nonOutlierSEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_TwoSubjOutlier(t_R_GRM_R_TwoSubjOutlierSEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R(t_R_GRM_RSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_MAF_interval(t_MAF_intervalSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type t_TwoSubj_list(t_TwoSubj_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type t_ThreeSubj_list(t_ThreeSubj_listSEXP);
    Rcpp::traits::input_parameter< double >::type t_SPA_Cutoff(t_SPA_CutoffSEXP);
    Rcpp::traits::input_parameter< double >::type t_zeta(t_zetaSEXP);
    Rcpp::traits::input_parameter< double >::type t_tol(t_tolSEXP);
    setSPAGRMobjInCPP(t_resid, t_resid_unrelated_outliers, t_sum_R_nonOutlier, t_R_GRM_R_nonOutlier, t_R_GRM_R_TwoSubjOutlier, t_R_GRM_R, t_MAF_interval, t_TwoSubj_list, t_ThreeSubj_list, t_SPA_Cutoff, t_zeta, t_tol);
    return R_NilValue;
END_RCPP
}
// setSAGELDobjInCPP
void setSAGELDobjInCPP(std::string t_Method, arma::mat t_XTs, arma::mat t_SS, arma::mat t_AtS, arma::mat t_Q, arma::mat t_A21, arma::mat t_TTs, arma::mat t_Tys, arma::vec t_sol, arma::vec t_blups, double t_sig, arma::vec t_resid, arma::vec t_resid_G, arma::vec t_resid_GxE, arma::vec t_resid_E, arma::vec t_resid_unrelated_outliers, arma::vec t_resid_unrelated_outliers_G, arma::vec t_resid_unrelated_outliers_GxE, double t_sum_R_nonOutlier, double t_sum_R_nonOutlier_G, double t_sum_R_nonOutlier_GxE, double t_R_GRM_R, double t_R_GRM_R_G, double t_R_GRM_R_GxE, double t_R_GRM_R_G_GxE, double t_R_GRM_R_E, double t_R_GRM_R_nonOutlier, double t_R_GRM_R_nonOutlier_G, double t_R_GRM_R_nonOutlier_GxE, double t_R_GRM_R_nonOutlier_G_GxE, double t_R_GRM_R_TwoSubjOutlier, double t_R_GRM_R_TwoSubjOutlier_G, double t_R_GRM_R_TwoSubjOutlier_GxE, double t_R_GRM_R_TwoSubjOutlier_G_GxE, Rcpp::List t_TwoSubj_list, Rcpp::List t_ThreeSubj_list, arma::vec t_MAF_interval, double t_zScoreE_cutoff, double t_SPA_Cutoff, double t_zeta, double t_tol);
RcppExport SEXP _GRAB_setSAGELDobjInCPP(SEXP t_MethodSEXP, SEXP t_XTsSEXP, SEXP t_SSSEXP, SEXP t_AtSSEXP, SEXP t_QSEXP, SEXP t_A21SEXP, SEXP t_TTsSEXP, SEXP t_TysSEXP, SEXP t_solSEXP, SEXP t_blupsSEXP, SEXP t_sigSEXP, SEXP t_residSEXP, SEXP t_resid_GSEXP, SEXP t_resid_GxESEXP, SEXP t_resid_ESEXP, SEXP t_resid_unrelated_outliersSEXP, SEXP t_resid_unrelated_outliers_GSEXP, SEXP t_resid_unrelated_outliers_GxESEXP, SEXP t_sum_R_nonOutlierSEXP, SEXP t_sum_R_nonOutlier_GSEXP, SEXP t_sum_R_nonOutlier_GxESEXP, SEXP t_R_GRM_RSEXP, SEXP t_R_GRM_R_GSEXP, SEXP t_R_GRM_R_GxESEXP, SEXP t_R_GRM_R_G_GxESEXP, SEXP t_R_GRM_R_ESEXP, SEXP t_R_GRM_R_nonOutlierSEXP, SEXP t_R_GRM_R_nonOutlier_GSEXP, SEXP t_R_GRM_R_nonOutlier_GxESEXP, SEXP t_R_GRM_R_nonOutlier_G_GxESEXP, SEXP t_R_GRM_R_TwoSubjOutlierSEXP, SEXP t_R_GRM_R_TwoSubjOutlier_GSEXP, SEXP t_R_GRM_R_TwoSubjOutlier_GxESEXP, SEXP t_R_GRM_R_TwoSubjOutlier_G_GxESEXP, SEXP t_TwoSubj_listSEXP, SEXP t_ThreeSubj_listSEXP, SEXP t_MAF_intervalSEXP, SEXP t_zScoreE_cutoffSEXP, SEXP t_SPA_CutoffSEXP, SEXP t_zetaSEXP, SEXP t_tolSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type t_Method(t_MethodSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_XTs(t_XTsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_SS(t_SSSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_AtS(t_AtSSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_Q(t_QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_A21(t_A21SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_TTs(t_TTsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_Tys(t_TysSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_sol(t_solSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_blups(t_blupsSEXP);
    Rcpp::traits::input_parameter< double >::type t_sig(t_sigSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_resid(t_residSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_resid_G(t_resid_GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_resid_GxE(t_resid_GxESEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_resid_E(t_resid_ESEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_resid_unrelated_outliers(t_resid_unrelated_outliersSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_resid_unrelated_outliers_G(t_resid_unrelated_outliers_GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_resid_unrelated_outliers_GxE(t_resid_unrelated_outliers_GxESEXP);
    Rcpp::traits::input_parameter< double >::type t_sum_R_nonOutlier(t_sum_R_nonOutlierSEXP);
    Rcpp::traits::input_parameter< double >::type t_sum_R_nonOutlier_G(t_sum_R_nonOutlier_GSEXP);
    Rcpp::traits::input_parameter< double >::type t_sum_R_nonOutlier_GxE(t_sum_R_nonOutlier_GxESEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R(t_R_GRM_RSEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_G(t_R_GRM_R_GSEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_GxE(t_R_GRM_R_GxESEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_G_GxE(t_R_GRM_R_G_GxESEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_E(t_R_GRM_R_ESEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_nonOutlier(t_R_GRM_R_nonOutlierSEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_nonOutlier_G(t_R_GRM_R_nonOutlier_GSEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_nonOutlier_GxE(t_R_GRM_R_nonOutlier_GxESEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_nonOutlier_G_GxE(t_R_GRM_R_nonOutlier_G_GxESEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_TwoSubjOutlier(t_R_GRM_R_TwoSubjOutlierSEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_TwoSubjOutlier_G(t_R_GRM_R_TwoSubjOutlier_GSEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_TwoSubjOutlier_GxE(t_R_GRM_R_TwoSubjOutlier_GxESEXP);
    Rcpp::traits::input_parameter< double >::type t_R_GRM_R_TwoSubjOutlier_G_GxE(t_R_GRM_R_TwoSubjOutlier_G_GxESEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type t_TwoSubj_list(t_TwoSubj_listSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type t_ThreeSubj_list(t_ThreeSubj_listSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_MAF_interval(t_MAF_intervalSEXP);
    Rcpp::traits::input_parameter< double >::type t_zScoreE_cutoff(t_zScoreE_cutoffSEXP);
    Rcpp::traits::input_parameter< double >::type t_SPA_Cutoff(t_SPA_CutoffSEXP);
    Rcpp::traits::input_parameter< double >::type t_zeta(t_zetaSEXP);
    Rcpp::traits::input_parameter< double >::type t_tol(t_tolSEXP);
    setSAGELDobjInCPP(t_Method, t_XTs, t_SS, t_AtS, t_Q, t_A21, t_TTs, t_Tys, t_sol, t_blups, t_sig, t_resid, t_resid_G, t_resid_GxE, t_resid_E, t_resid_unrelated_outliers, t_resid_unrelated_outliers_G, t_resid_unrelated_outliers_GxE, t_sum_R_nonOutlier, t_sum_R_nonOutlier_G, t_sum_R_nonOutlier_GxE, t_R_GRM_R, t_R_GRM_R_G, t_R_GRM_R_GxE, t_R_GRM_R_G_GxE, t_R_GRM_R_E, t_R_GRM_R_nonOutlier, t_R_GRM_R_nonOutlier_G, t_R_GRM_R_nonOutlier_GxE, t_R_GRM_R_nonOutlier_G_GxE, t_R_GRM_R_TwoSubjOutlier, t_R_GRM_R_TwoSubjOutlier_G, t_R_GRM_R_TwoSubjOutlier_GxE, t_R_GRM_R_TwoSubjOutlier_G_GxE, t_TwoSubj_list, t_ThreeSubj_list, t_MAF_interval, t_zScoreE_cutoff, t_SPA_Cutoff, t_zeta, t_tol);
    return R_NilValue;
END_RCPP
}
// setSPAmixobjInCPP
void setSPAmixobjInCPP(arma::mat t_resid, arma::mat t_PCs, int t_N, double t_SPA_Cutoff, Rcpp::List t_outlierList);
RcppExport SEXP _GRAB_setSPAmixobjInCPP(SEXP t_residSEXP, SEXP t_PCsSEXP, SEXP t_NSEXP, SEXP t_SPA_CutoffSEXP, SEXP t_outlierListSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type t_resid(t_residSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_PCs(t_PCsSEXP);
    Rcpp::traits::input_parameter< int >::type t_N(t_NSEXP);
    Rcpp::traits::input_parameter< double >::type t_SPA_Cutoff(t_SPA_CutoffSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type t_outlierList(t_outlierListSEXP);
    setSPAmixobjInCPP(t_resid, t_PCs, t_N, t_SPA_Cutoff, t_outlierList);
    return R_NilValue;
END_RCPP
}
// setSPACoxobjInCPP
void setSPACoxobjInCPP(arma::mat t_cumul, arma::vec t_mresid, arma::mat t_XinvXX, arma::mat t_tX, int t_N, double t_pVal_covaAdj_Cutoff, double t_SPA_Cutoff);
RcppExport SEXP _GRAB_setSPACoxobjInCPP(SEXP t_cumulSEXP, SEXP t_mresidSEXP, SEXP t_XinvXXSEXP, SEXP t_tXSEXP, SEXP t_NSEXP, SEXP t_pVal_covaAdj_CutoffSEXP, SEXP t_SPA_CutoffSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type t_cumul(t_cumulSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_mresid(t_mresidSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_XinvXX(t_XinvXXSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type t_tX(t_tXSEXP);
    Rcpp::traits::input_parameter< int >::type t_N(t_NSEXP);
    Rcpp::traits::input_parameter< double >::type t_pVal_covaAdj_Cutoff(t_pVal_covaAdj_CutoffSEXP);
    Rcpp::traits::input_parameter< double >::type t_SPA_Cutoff(t_SPA_CutoffSEXP);
    setSPACoxobjInCPP(t_cumul, t_mresid, t_XinvXX, t_tX, t_N, t_pVal_covaAdj_Cutoff, t_SPA_Cutoff);
    return R_NilValue;
END_RCPP
}
// setWtCoxGobjInCPP
void setWtCoxGobjInCPP(arma::vec t_mresid, arma::vec t_weight, std::string t_imputeMethod, double t_cutoff);
RcppExport SEXP _GRAB_setWtCoxGobjInCPP(SEXP t_mresidSEXP, SEXP t_weightSEXP, SEXP t_imputeMethodSEXP, SEXP t_cutoffSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type t_mresid(t_mresidSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type t_weight(t_weightSEXP);
    Rcpp::traits::input_parameter< std::string >::type t_imputeMethod(t_imputeMethodSEXP);
    Rcpp::traits::input_parameter< double >::type t_cutoff(t_cutoffSEXP);
    setWtCoxGobjInCPP(t_mresid, t_weight, t_imputeMethod, t_cutoff);
    return R_NilValue;
END_RCPP
}
// updateWtCoxGChunkInCPP
void updateWtCoxGChunkInCPP(DataFrame t_mergeGenoInfo_subset);
RcppExport SEXP _GRAB_updateWtCoxGChunkInCPP(SEXP t_mergeGenoInfo_subsetSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< DataFrame >::type t_mergeGenoInfo_subset(t_mergeGenoInfo_subsetSEXP);
    updateWtCoxGChunkInCPP(t_mergeGenoInfo_subset);
    return R_NilValue;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_GRAB_setSparseGRMInCPP", (DL_FUNC) &_GRAB_setSparseGRMInCPP, 1},
    {"_GRAB_setDenseGRMInCPP", (DL_FUNC) &_GRAB_setDenseGRMInCPP, 3},
    {"_GRAB_getDenseGRMInCPP", (DL_FUNC) &_GRAB_getDenseGRMInCPP, 3},
    {"_GRAB_setMarker_GlobalVarsInCPP", (DL_FUNC) &_GRAB_setMarker_GlobalVarsInCPP, 8},
    {"_GRAB_setRegion_GlobalVarsInCPP", (DL_FUNC) &_GRAB_setRegion_GlobalVarsInCPP, 8},
    {"_GRAB_mainMarkerInCPP", (DL_FUNC) &_GRAB_mainMarkerInCPP, 3},
    {"_GRAB_mainRegionURVInCPP", (DL_FUNC) &_GRAB_mainRegionURVInCPP, 4},
    {"_GRAB_mainRegionInCPP", (DL_FUNC) &_GRAB_mainRegionInCPP, 9},
    {"_GRAB_printTimeDiffInCPP", (DL_FUNC) &_GRAB_printTimeDiffInCPP, 0},
    {"_GRAB_printTimeDiffSPAmixInCPP", (DL_FUNC) &_GRAB_printTimeDiffSPAmixInCPP, 0},
    {"_GRAB_getGenoInfoInCPP", (DL_FUNC) &_GRAB_getGenoInfoInCPP, 3},
    {"_GRAB_getGenoInCPP", (DL_FUNC) &_GRAB_getGenoInCPP, 4},
    {"_GRAB_getGenoInCPP_fixedNumber", (DL_FUNC) &_GRAB_getGenoInCPP_fixedNumber, 7},
    {"_GRAB_getSpGenoInCPP", (DL_FUNC) &_GRAB_getSpGenoInCPP, 4},
    {"_GRAB_setPLINKobjInCPP", (DL_FUNC) &_GRAB_setPLINKobjInCPP, 5},
    {"_GRAB_setBGENobjInCPP", (DL_FUNC) &_GRAB_setBGENobjInCPP, 7},
    {"_GRAB_closeGenoInputInCPP", (DL_FUNC) &_GRAB_closeGenoInputInCPP, 1},
    {"_GRAB_setPOLMMobjInCPP", (DL_FUNC) &_GRAB_setPOLMMobjInCPP, 11},
    {"_GRAB_setPOLMMobjInCPP_NULL", (DL_FUNC) &_GRAB_setPOLMMobjInCPP_NULL, 10},
    {"_GRAB_setSPAGRMobjInCPP", (DL_FUNC) &_GRAB_setSPAGRMobjInCPP, 12},
    {"_GRAB_setSAGELDobjInCPP", (DL_FUNC) &_GRAB_setSAGELDobjInCPP, 41},
    {"_GRAB_setSPAmixobjInCPP", (DL_FUNC) &_GRAB_setSPAmixobjInCPP, 5},
    {"_GRAB_setSPACoxobjInCPP", (DL_FUNC) &_GRAB_setSPACoxobjInCPP, 7},
    {"_GRAB_setWtCoxGobjInCPP", (DL_FUNC) &_GRAB_setWtCoxGobjInCPP, 4},
    {"_GRAB_updateWtCoxGChunkInCPP", (DL_FUNC) &_GRAB_updateWtCoxGChunkInCPP, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_GRAB(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
