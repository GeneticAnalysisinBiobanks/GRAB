# SPAsqr_Marker.R -- Unwrap SPAsqr null model and run marker engine

runMarker.SPAsqr <- function(objNull, control, bedFile, bimFile, famFile, OutputFile, nThreads) {
  
  taus  <- objNull$taus
  ntaus <- length(taus)

  residOutlierLst    <- objNull[["Resid.unrelated.outliers_lst"]]
  twoSubjLstLst      <- objNull$TwoSubj_list_lst
  cltUnionLst        <- objNull$CLT_union_lst
  threeSubjFamIdxLst <- objNull$ThreeSubj_family_idx_lst
  threeSubjStandSLst <- objNull$ThreeSubj_stand_S_lst

  twoSubj_resid_lists <- vector("list", ntaus)
  twoSubj_rho_lists   <- vector("list", ntaus)
  for (i in seq_len(ntaus)) {
    pairs <- twoSubjLstLst[[i]]
    twoSubj_resid_lists[[i]] <- lapply(pairs, function(p) p$Resid)
    twoSubj_rho_lists[[i]]   <- lapply(pairs, function(p) p$Rho)
  }

  threeSubj_standS_lists <- vector("list", ntaus)
  threeSubj_CLT_lists    <- vector("list", ntaus)
  for (i in seq_len(ntaus)) {
    famIdx  <- threeSubjFamIdxLst[[i]]
    standS  <- threeSubjStandSLst[[i]]
    nFam    <- length(famIdx)
    threeSubj_standS_lists[[i]] <- standS
    threeSubj_CLT_lists[[i]]    <- lapply(seq_len(nFam), function(j) {
      cltUnionLst[[ famIdx[j] ]]
    })
  }

  resid_outlier_list <- lapply(seq_len(ntaus), function(i) {
    v <- residOutlierLst[[i]]
    if (is.null(v)) numeric(0) else v
  })

  # Flatten resid_outlier: concat across taus
  resid_outlier_all  <- unlist(resid_outlier_list)
  if (is.null(resid_outlier_all)) resid_outlier_all <- numeric(0)
  resid_outlier_lens <- as.integer(sapply(resid_outlier_list, length))

  # Flatten twoSubj: stack all 2-elem vecs across taus into (totalTwo x 2) matrix
  all_two_resid <- do.call(c, twoSubj_resid_lists)
  all_two_rho   <- do.call(c, twoSubj_rho_lists)
  twoSubj_resid_all <- if (length(all_two_resid) > 0) do.call(rbind, all_two_resid) else matrix(numeric(0), nrow = 0, ncol = 2)
  twoSubj_rho_all   <- if (length(all_two_rho) > 0)   do.call(rbind, all_two_rho)   else matrix(numeric(0), nrow = 0, ncol = 2)
  twoSubj_perTau <- as.integer(sapply(twoSubj_resid_lists, length))

  # Flatten threeSubj: concat all families across taus
  all_three_standS <- do.call(c, threeSubj_standS_lists)
  threeSubj_standS_all  <- if (length(all_three_standS) > 0) unlist(all_three_standS) else numeric(0)
  threeSubj_standS_lens <- if (length(all_three_standS) > 0) as.integer(sapply(all_three_standS, length)) else integer(0)
  all_three_CLT <- do.call(c, threeSubj_CLT_lists)
  threeSubj_CLT_all   <- if (length(all_three_CLT) > 0) do.call(rbind, all_three_CLT) else matrix(numeric(0), nrow = 0, ncol = 1)
  threeSubj_CLT_nrows <- if (length(all_three_CLT) > 0) as.integer(sapply(all_three_CLT, nrow)) else integer(0)
  threeSubj_perTau <- as.integer(sapply(threeSubj_standS_lists, length))

  runMarkerInCPP.SPAsqr(
    taus                       = taus,
    Resid_mat                  = objNull$Resid_mat,
    resid_outlier_all          = resid_outlier_all,
    resid_outlier_lens         = resid_outlier_lens,
    sum_R_nonOutlier_vec       = objNull$sum_R_nonOutlier_vec,
    R_GRM_R_nonOutlier_vec     = objNull$R_GRM_R_nonOutlier_vec,
    R_GRM_R_TwoSubjOutlier_vec = objNull$R_GRM_R_TwoSubjOutlier_vec,
    R_GRM_R_vec                = objNull$R_GRM_R_vec,
    MAF_interval               = objNull$MAF_interval,
    twoSubj_resid_all          = twoSubj_resid_all,
    twoSubj_rho_all            = twoSubj_rho_all,
    twoSubj_perTau             = twoSubj_perTau,
    threeSubj_standS_all       = threeSubj_standS_all,
    threeSubj_standS_lens      = threeSubj_standS_lens,
    threeSubj_CLT_all          = threeSubj_CLT_all,
    threeSubj_CLT_nrows        = threeSubj_CLT_nrows,
    threeSubj_perTau           = threeSubj_perTau,
    SPA_Cutoff                 = control$SPA_Cutoff,
    zeta                       = control$zeta,
    tol                        = control$tol,
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
