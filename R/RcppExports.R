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

setMarker_GlobalVarsInCPP <- function(t_impute_method, t_missing_cutoff, t_min_maf_marker, t_min_mac_marker, t_omp_num_threads, t_group, t_ifOutGroup, t_nGroup) {
    invisible(.Call(`_GRAB_setMarker_GlobalVarsInCPP`, t_impute_method, t_missing_cutoff, t_min_maf_marker, t_min_mac_marker, t_omp_num_threads, t_group, t_ifOutGroup, t_nGroup))
}

setRegion_GlobalVarsInCPP <- function(t_impute_method, t_missing_cutoff, t_max_maf_region, t_min_mac_region, t_max_markers_region, t_omp_num_threads, t_region_weight_beta, t_region_max_maf_vec) {
    invisible(.Call(`_GRAB_setRegion_GlobalVarsInCPP`, t_impute_method, t_missing_cutoff, t_max_maf_region, t_min_mac_region, t_max_markers_region, t_omp_num_threads, t_region_weight_beta, t_region_max_maf_vec))
}

mainMarkerInCPP <- function(t_method, t_genoType, t_genoIndex) {
    .Call(`_GRAB_mainMarkerInCPP`, t_method, t_genoType, t_genoIndex)
}

mainRegionURVInCPP <- function(t_method, t_genoType, t_genoIndex, t_n) {
    .Call(`_GRAB_mainRegionURVInCPP`, t_method, t_genoType, t_genoIndex, t_n)
}

mainRegionInCPP <- function(t_method, t_genoType, t_genoIndex, t_weightVec, t_outputFile, t_labelVec, t_nLabel, t_annoMat, t_annoVec) {
    .Call(`_GRAB_mainRegionInCPP`, t_method, t_genoType, t_genoIndex, t_weightVec, t_outputFile, t_labelVec, t_nLabel, t_annoMat, t_annoVec)
}

printTimeDiffInCPP <- function() {
    invisible(.Call(`_GRAB_printTimeDiffInCPP`))
}

printTimeDiffSPAmixInCPP <- function() {
    invisible(.Call(`_GRAB_printTimeDiffSPAmixInCPP`))
}

getGenoInfoInCPP <- function(t_genoType, t_markerInfo, t_imputeMethod) {
    .Call(`_GRAB_getGenoInfoInCPP`, t_genoType, t_markerInfo, t_imputeMethod)
}

getGenoInCPP <- function(t_genoType, t_markerInfo, n, t_imputeMethod) {
    .Call(`_GRAB_getGenoInCPP`, t_genoType, t_markerInfo, n, t_imputeMethod)
}

getGenoInCPP_fixedNumber <- function(t_genoType, t_markerInfo, n, t_imputeMethod, m, missingRateCutoff, minMAFCutoff) {
    .Call(`_GRAB_getGenoInCPP_fixedNumber`, t_genoType, t_markerInfo, n, t_imputeMethod, m, missingRateCutoff, minMAFCutoff)
}

getSpGenoInCPP <- function(t_genoType, t_markerInfo, n, t_imputeMethod) {
    .Call(`_GRAB_getSpGenoInCPP`, t_genoType, t_markerInfo, n, t_imputeMethod)
}

setPLINKobjInCPP <- function(t_bimFile, t_famFile, t_bedFile, t_SampleInModel, t_AlleleOrder) {
    invisible(.Call(`_GRAB_setPLINKobjInCPP`, t_bimFile, t_famFile, t_bedFile, t_SampleInModel, t_AlleleOrder))
}

setBGENobjInCPP <- function(t_bgenFileName, t_bgenFileIndex, t_SampleInBgen, t_SampleInModel, t_isSparseDosageInBgen, t_isDropmissingdosagesInBgen, t_AlleleOrder) {
    invisible(.Call(`_GRAB_setBGENobjInCPP`, t_bgenFileName, t_bgenFileIndex, t_SampleInBgen, t_SampleInModel, t_isSparseDosageInBgen, t_isDropmissingdosagesInBgen, t_AlleleOrder))
}

closeGenoInputInCPP <- function(t_genoType) {
    invisible(.Call(`_GRAB_closeGenoInputInCPP`, t_genoType))
}

setPOLMMobjInCPP <- function(t_muMat, t_iRMat, t_Cova, t_yVec, t_tau, t_printPCGInfo, t_tolPCG, t_maxiterPCG, t_varRatio, t_SPA_cutoff, t_flagSparseGRM) {
    invisible(.Call(`_GRAB_setPOLMMobjInCPP`, t_muMat, t_iRMat, t_Cova, t_yVec, t_tau, t_printPCGInfo, t_tolPCG, t_maxiterPCG, t_varRatio, t_SPA_cutoff, t_flagSparseGRM))
}

setPOLMMobjInCPP_NULL <- function(t_flagSparseGRM, t_Cova, t_yVec, t_beta, t_bVec, t_eps, t_tau, t_SPmatR, t_controlList, GenoMat) {
    .Call(`_GRAB_setPOLMMobjInCPP_NULL`, t_flagSparseGRM, t_Cova, t_yVec, t_beta, t_bVec, t_eps, t_tau, t_SPmatR, t_controlList, GenoMat)
}

setSPAGRMobjInCPP <- function(t_resid, t_resid_unrelated_outliers, t_sum_R_nonOutlier, t_R_GRM_R_nonOutlier, t_R_GRM_R_TwoSubjOutlier, t_R_GRM_R, t_MAF_interval, t_TwoSubj_list, t_ThreeSubj_list, t_SPA_Cutoff, t_zeta, t_tol) {
    invisible(.Call(`_GRAB_setSPAGRMobjInCPP`, t_resid, t_resid_unrelated_outliers, t_sum_R_nonOutlier, t_R_GRM_R_nonOutlier, t_R_GRM_R_TwoSubjOutlier, t_R_GRM_R, t_MAF_interval, t_TwoSubj_list, t_ThreeSubj_list, t_SPA_Cutoff, t_zeta, t_tol))
}

setSAGELDobjInCPP <- function(t_Method, t_XTs, t_SS, t_AtS, t_Q, t_A21, t_TTs, t_Tys, t_sol, t_blups, t_sig, t_resid, t_resid_G, t_resid_GxE, t_resid_E, t_resid_unrelated_outliers, t_resid_unrelated_outliers_G, t_resid_unrelated_outliers_GxE, t_sum_R_nonOutlier, t_sum_R_nonOutlier_G, t_sum_R_nonOutlier_GxE, t_R_GRM_R, t_R_GRM_R_G, t_R_GRM_R_GxE, t_R_GRM_R_G_GxE, t_R_GRM_R_E, t_R_GRM_R_nonOutlier, t_R_GRM_R_nonOutlier_G, t_R_GRM_R_nonOutlier_GxE, t_R_GRM_R_nonOutlier_G_GxE, t_R_GRM_R_TwoSubjOutlier, t_R_GRM_R_TwoSubjOutlier_G, t_R_GRM_R_TwoSubjOutlier_GxE, t_R_GRM_R_TwoSubjOutlier_G_GxE, t_TwoSubj_list, t_ThreeSubj_list, t_MAF_interval, t_zScoreE_cutoff, t_SPA_Cutoff, t_zeta, t_tol) {
    invisible(.Call(`_GRAB_setSAGELDobjInCPP`, t_Method, t_XTs, t_SS, t_AtS, t_Q, t_A21, t_TTs, t_Tys, t_sol, t_blups, t_sig, t_resid, t_resid_G, t_resid_GxE, t_resid_E, t_resid_unrelated_outliers, t_resid_unrelated_outliers_G, t_resid_unrelated_outliers_GxE, t_sum_R_nonOutlier, t_sum_R_nonOutlier_G, t_sum_R_nonOutlier_GxE, t_R_GRM_R, t_R_GRM_R_G, t_R_GRM_R_GxE, t_R_GRM_R_G_GxE, t_R_GRM_R_E, t_R_GRM_R_nonOutlier, t_R_GRM_R_nonOutlier_G, t_R_GRM_R_nonOutlier_GxE, t_R_GRM_R_nonOutlier_G_GxE, t_R_GRM_R_TwoSubjOutlier, t_R_GRM_R_TwoSubjOutlier_G, t_R_GRM_R_TwoSubjOutlier_GxE, t_R_GRM_R_TwoSubjOutlier_G_GxE, t_TwoSubj_list, t_ThreeSubj_list, t_MAF_interval, t_zScoreE_cutoff, t_SPA_Cutoff, t_zeta, t_tol))
}

setSPAmixobjInCPP <- function(t_resid, t_PCs, t_N, t_SPA_Cutoff, t_outlierList) {
    invisible(.Call(`_GRAB_setSPAmixobjInCPP`, t_resid, t_PCs, t_N, t_SPA_Cutoff, t_outlierList))
}

setSPACoxobjInCPP <- function(t_cumul, t_mresid, t_XinvXX, t_tX, t_N, t_pVal_covaAdj_Cutoff, t_SPA_Cutoff) {
    invisible(.Call(`_GRAB_setSPACoxobjInCPP`, t_cumul, t_mresid, t_XinvXX, t_tX, t_N, t_pVal_covaAdj_Cutoff, t_SPA_Cutoff))
}

setWtCoxGobjInCPP <- function(t_mresid, t_weight, t_imputeMethod, t_cutoff) {
    invisible(.Call(`_GRAB_setWtCoxGobjInCPP`, t_mresid, t_weight, t_imputeMethod, t_cutoff))
}

updateWtCoxGChunkInCPP <- function(t_mergeGenoInfo_subset) {
    invisible(.Call(`_GRAB_updateWtCoxGChunkInCPP`, t_mergeGenoInfo_subset))
}

