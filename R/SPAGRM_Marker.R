# SPAGRM_Marker.R -- Unwrap SPAGRM null model and run marker engine

runMarker.SPAGRM <- function(objNull, OutputFile, control, bedFile, bimFile, famFile) {

  # Pack 2-element vectors into matrices (nTwo x 2)
  twoSubjList   <- objNull$TwoSubj_list
  twoSubj_resid <- if (length(twoSubjList) > 0) do.call(rbind, lapply(twoSubjList, function(x) x$Resid)) else matrix(numeric(0), nrow = 0, ncol = 2)
  twoSubj_rho   <- if (length(twoSubjList) > 0) do.call(rbind, lapply(twoSubjList, function(x) x$Rho))   else matrix(numeric(0), nrow = 0, ncol = 2)

  # Concatenate variable-length vectors + stacked matrices
  threeSubjList         <- objNull$ThreeSubj_list
  standS_list           <- lapply(threeSubjList, function(x) x[["stand.S"]])
  CLT_list              <- lapply(threeSubjList, function(x) x$CLT)
  threeSubj_standS_all  <- if (length(standS_list) > 0) unlist(standS_list) else numeric(0)
  threeSubj_standS_lens <- as.integer(sapply(standS_list, length))
  threeSubj_CLT_all     <- if (length(CLT_list) > 0) do.call(rbind, CLT_list) else matrix(numeric(0), nrow = 0, ncol = 1)
  threeSubj_CLT_nrows   <- as.integer(sapply(CLT_list, nrow))

  runMarkerInCPP.SPAGRM(
    resid                    = objNull$Resid,
    resid_unrelated_outliers = objNull[["Resid.unrelated.outliers"]],
    sum_R_nonOutlier         = objNull$sum_R_nonOutlier,
    R_GRM_R_nonOutlier       = objNull$R_GRM_R_nonOutlier,
    R_GRM_R_TwoSubjOutlier   = objNull$R_GRM_R_TwoSubjOutlier,
    R_GRM_R                  = objNull$R_GRM_R,
    MAF_interval             = objNull$MAF_interval,
    twoSubj_resid            = twoSubj_resid,
    twoSubj_rho              = twoSubj_rho,
    threeSubj_standS_all     = threeSubj_standS_all,
    threeSubj_standS_lens    = threeSubj_standS_lens,
    threeSubj_CLT_all        = threeSubj_CLT_all,
    threeSubj_CLT_nrows      = threeSubj_CLT_nrows,
    SPA_Cutoff               = control$SPA_Cutoff,
    zeta                     = control$zeta,
    tol                      = control$tol,
    bedFile             = bedFile,
    bimFile             = bimFile,
    famFile             = famFile,
    outputFile          = OutputFile,
    subjData            = objNull$subjData,
    AlleleOrder         = control$AlleleOrder,
    nMarkersEachChunk   = as.integer(control$nMarkersEachChunk),
    nthreads            = as.integer(control$nthreads),
    impute_method       = control$impute_method,
    missing_cutoff      = as.numeric(control$missing_cutoff),
    min_maf_marker      = as.numeric(control$min_maf_marker),
    min_mac_marker      = as.numeric(control$min_mac_marker),
    IDsToIncludeFile    = if (is.null(control$IDsToIncludeFile)) "" else control$IDsToIncludeFile,
    RangesToIncludeFile = if (is.null(control$RangesToIncludeFile)) "" else control$RangesToIncludeFile,
    IDsToExcludeFile    = if (is.null(control$IDsToExcludeFile)) "" else control$IDsToExcludeFile,
    RangesToExcludeFile = if (is.null(control$RangesToExcludeFile)) "" else control$RangesToExcludeFile
  )
}
