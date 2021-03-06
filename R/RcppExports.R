# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

setSparseGRMInCPP <- function(t_KinMatListR) {
    invisible(.Call(`_GRAB_setSparseGRMInCPP`, t_KinMatListR))
}

setDenseGRMInCPP <- function(t_memoryChunk, t_minMafGRM, t_maxMissingGRM) {
    invisible(.Call(`_GRAB_setDenseGRMInCPP`, t_memoryChunk, t_minMafGRM, t_maxMissingGRM))
}

getDenseGRMInCPP <- function(t_bVec, t_excludeChr, t_grainSize) {
    .Call(`_GRAB_getDenseGRMInCPP`, t_bVec, t_excludeChr, t_grainSize)
}

setMarker_GlobalVarsInCPP <- function(t_impute_method, t_missing_cutoff, t_min_maf_marker, t_min_mac_marker, t_omp_num_threads) {
    invisible(.Call(`_GRAB_setMarker_GlobalVarsInCPP`, t_impute_method, t_missing_cutoff, t_min_maf_marker, t_min_mac_marker, t_omp_num_threads))
}

setRegion_GlobalVarsInCPP <- function(t_impute_method, t_missing_cutoff, t_max_maf_region, t_max_markers_region, t_omp_num_threads) {
    invisible(.Call(`_GRAB_setRegion_GlobalVarsInCPP`, t_impute_method, t_missing_cutoff, t_max_maf_region, t_max_markers_region, t_omp_num_threads))
}

mainMarkerInCPP <- function(t_method, t_genoType, t_genoIndex) {
    .Call(`_GRAB_mainMarkerInCPP`, t_method, t_genoType, t_genoIndex)
}

mainRegionInCPP <- function(t_method, t_genoType, t_genoIndex, t_outputFile, t_n) {
    .Call(`_GRAB_mainRegionInCPP`, t_method, t_genoType, t_genoIndex, t_outputFile, t_n)
}

getGenoInCPP <- function(t_genoType, t_markerInfo, n) {
    .Call(`_GRAB_getGenoInCPP`, t_genoType, t_markerInfo, n)
}

getSpGenoInCPP <- function(t_genoType, t_markerInfo, n) {
    .Call(`_GRAB_getSpGenoInCPP`, t_genoType, t_markerInfo, n)
}

setPLINKobjInCPP <- function(t_bimFile, t_famFile, t_bedFile, t_SampleInModel, t_AlleleOrder) {
    invisible(.Call(`_GRAB_setPLINKobjInCPP`, t_bimFile, t_famFile, t_bedFile, t_SampleInModel, t_AlleleOrder))
}

setBGENobjInCPP <- function(t_bgenFileName, t_bgenFileIndex, t_SampleInBgen, t_SampleInModel, t_isSparseDosageInBgen, t_isDropmissingdosagesInBgen, t_AlleleOrder) {
    invisible(.Call(`_GRAB_setBGENobjInCPP`, t_bgenFileName, t_bgenFileIndex, t_SampleInBgen, t_SampleInModel, t_isSparseDosageInBgen, t_isDropmissingdosagesInBgen, t_AlleleOrder))
}

setPOLMMobjInCPP <- function(t_muMat, t_iRMat, t_Cova, t_yVec, t_SPmatR, t_tau, t_printPCGInfo, t_tolPCG, t_maxiterPCG, t_varRatio, t_SPA_cutoff, t_flagSparseGRM) {
    invisible(.Call(`_GRAB_setPOLMMobjInCPP`, t_muMat, t_iRMat, t_Cova, t_yVec, t_SPmatR, t_tau, t_printPCGInfo, t_tolPCG, t_maxiterPCG, t_varRatio, t_SPA_cutoff, t_flagSparseGRM))
}

setPOLMMobjInCPP_NULL <- function(t_flagSparseGRM, t_Cova, t_yVec, t_beta, t_bVec, t_eps, t_tau, t_SPmatR, t_controlList) {
    .Call(`_GRAB_setPOLMMobjInCPP_NULL`, t_flagSparseGRM, t_Cova, t_yVec, t_beta, t_bVec, t_eps, t_tau, t_SPmatR, t_controlList)
}

setSPACoxobjInCPP <- function(t_cumul, t_mresid, t_XinvXX, t_tX, t_N, t_pVal_covaAdj_Cutoff, t_SPA_Cutoff) {
    invisible(.Call(`_GRAB_setSPACoxobjInCPP`, t_cumul, t_mresid, t_XinvXX, t_tX, t_N, t_pVal_covaAdj_Cutoff, t_SPA_Cutoff))
}

squares <- function(data) {
    .Call(`_GRAB_squares`, data)
}

