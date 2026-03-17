# SAGELD_Marker.R -- Unwrap SAGELD null model and run marker engine

runMarker.SAGELD <- function(objNull, control, bedFile, bimFile, famFile, OutputFile, nThreads) {

  twoSubjList <- objNull$TwoSubj_list
  # Pack 2-element vectors into matrices (nTwo x 2)
  twoSubj_Resid     <- if (length(twoSubjList) > 0) do.call(rbind, lapply(twoSubjList, function(x) x$Resid))     else matrix(numeric(0), nrow = 0, ncol = 2)
  twoSubj_Rho       <- if (length(twoSubjList) > 0) do.call(rbind, lapply(twoSubjList, function(x) x$Rho))       else matrix(numeric(0), nrow = 0, ncol = 2)
  twoSubj_Resid_G   <- if (length(twoSubjList) > 0) do.call(rbind, lapply(twoSubjList, function(x) x$Resid_G))   else matrix(numeric(0), nrow = 0, ncol = 2)
  twoSubj_Resid_GxE <- if (length(twoSubjList) > 0) do.call(rbind, lapply(twoSubjList, function(x) x$Resid_GxE)) else matrix(numeric(0), nrow = 0, ncol = 2)

  threeSubjList <- objNull$ThreeSubj_list
  CLT_list        <- lapply(threeSubjList, function(x) x$CLT)
  standS_list     <- lapply(threeSubjList, function(x) x[["stand.S"]])
  standS_G_list   <- lapply(threeSubjList, function(x) x[["stand.S_G"]])
  standS_GxE_list <- lapply(threeSubjList, function(x) x[["stand.S_GxE"]])
  # Concatenate variable-length vectors + stacked matrices
  threeSubj_standS_all      <- if (length(standS_list) > 0) unlist(standS_list) else numeric(0)
  threeSubj_standS_lens     <- as.integer(sapply(standS_list, length))
  threeSubj_standS_G_all    <- if (length(standS_G_list) > 0) unlist(standS_G_list) else numeric(0)
  threeSubj_standS_G_lens   <- as.integer(sapply(standS_G_list, length))
  threeSubj_standS_GxE_all  <- if (length(standS_GxE_list) > 0) unlist(standS_GxE_list) else numeric(0)
  threeSubj_standS_GxE_lens <- as.integer(sapply(standS_GxE_list, length))
  threeSubj_CLT_all         <- if (length(CLT_list) > 0) do.call(rbind, CLT_list) else matrix(numeric(0), nrow = 0, ncol = 1)
  threeSubj_CLT_nrows       <- as.integer(sapply(CLT_list, nrow))

  runMarkerInCPP.SAGELD(
    Method                       = objNull$Method,
    XTs                          = objNull$XTs,
    SS                           = objNull$SS,
    AtS                          = objNull$AtS,
    Q                            = objNull$Q,
    A21                          = objNull$A21,
    TTs                          = objNull$TTs,
    Tys                          = objNull$Tys,
    sol                          = objNull$sol,
    blups                        = objNull$blups,
    sig                          = objNull$sig,
    resid                        = objNull$Resid,
    resid_G                      = objNull$Resid_G,
    resid_GxE                    = objNull$Resid_GxE,
    resid_E                      = objNull$Resid_E,
    resid_unrelated_outliers     = objNull[["Resid.unrelated.outliers"]],
    resid_unrelated_outliers_G   = objNull[["Resid.unrelated.outliers_G"]],
    resid_unrelated_outliers_GxE = objNull[["Resid.unrelated.outliers_GxE"]],
    sum_R_nonOutlier             = objNull$sum_R_nonOutlier,
    sum_R_nonOutlier_G           = objNull$sum_R_nonOutlier_G,
    sum_R_nonOutlier_GxE         = objNull$sum_R_nonOutlier_GxE,
    R_GRM_R                      = c(objNull$R_GRM_R, objNull$R_GRM_R_G, objNull$R_GRM_R_GxE, objNull$R_GRM_R_G_GxE, objNull$R_GRM_R_E),
    R_GRM_R_nonOutlier           = c(objNull$R_GRM_R_nonOutlier, objNull$R_GRM_R_nonOutlier_G, objNull$R_GRM_R_nonOutlier_GxE, objNull$R_GRM_R_nonOutlier_G_GxE),
    R_GRM_R_TwoSubjOutlier       = c(objNull$R_GRM_R_TwoSubjOutlier, objNull$R_GRM_R_TwoSubjOutlier_G, objNull$R_GRM_R_TwoSubjOutlier_GxE, objNull$R_GRM_R_TwoSubjOutlier_G_GxE),
    twoSubj_Resid                = twoSubj_Resid,
    twoSubj_Rho                  = twoSubj_Rho,
    twoSubj_Resid_G              = twoSubj_Resid_G,
    twoSubj_Resid_GxE            = twoSubj_Resid_GxE,
    threeSubj_standS_all         = threeSubj_standS_all,
    threeSubj_standS_lens        = threeSubj_standS_lens,
    threeSubj_standS_G_all       = threeSubj_standS_G_all,
    threeSubj_standS_G_lens      = threeSubj_standS_G_lens,
    threeSubj_standS_GxE_all     = threeSubj_standS_GxE_all,
    threeSubj_standS_GxE_lens    = threeSubj_standS_GxE_lens,
    threeSubj_CLT_all            = threeSubj_CLT_all,
    threeSubj_CLT_nrows          = threeSubj_CLT_nrows,
    MAF_interval                 = objNull$MAF_interval,
    zScoreE_cutoff               = objNull$zScoreE_cutoff,
    SPA_Cutoff                   = control$SPA_Cutoff,
    zeta                         = control$zeta,
    tol                          = control$tol,
    bedFile             = bedFile,
    bimFile             = bimFile,
    famFile             = famFile,
    outputFile          = OutputFile,
    sampleInModel       = objNull$subjData,
    alleleOrder         = if (is.null(control$AlleleOrder)) "alt-first" else control$AlleleOrder,
    nMarkersEachChunk   = as.integer(control$nMarkersEachChunk),
    nThreads            = nThreads,
    imputeMethod        = control$impute_method,
    missingCutoff       = control$missing_cutoff,
    minMafMarker        = control$min_maf_marker,
    minMacMarker        = control$min_mac_marker,
    idsToIncludeFile    = if (is.null(control$IDsToIncludeFile)) "" else control$IDsToIncludeFile,
    rangesToIncludeFile = if (is.null(control$RangesToIncludeFile)) "" else control$RangesToIncludeFile,
    idsToExcludeFile    = if (is.null(control$IDsToExcludeFile)) "" else control$IDsToExcludeFile,
    rangesToExcludeFile = if (is.null(control$RangesToExcludeFile)) "" else control$RangesToExcludeFile
  )
}
