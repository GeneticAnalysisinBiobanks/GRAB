# LEAF_Marker.R -- Unwrap LEAF null model and run marker engine

runMarker.LEAF <- function(objNull, control, bedFile, bimFile, famFile,
                           OutputFile, nThreads) {
  nCluster   <- objNull$Ncluster
  clusterIdx <- objNull$clusterIdx

  clusterIdx_list <- lapply(seq_len(nCluster), function(k) {
    which(clusterIdx == k)
  })

  # Concat per-cluster vectors
  residuals_all    <- unlist(objNull$residuals_list)
  residuals_lens   <- as.integer(sapply(objNull$residuals_list, length))
  weights_all      <- unlist(objNull$weights_list)
  weights_lens     <- as.integer(sapply(objNull$weights_list, length))
  clusterIdx_all   <- as.integer(unlist(clusterIdx_list))
  clusterIdx_lens  <- as.integer(sapply(clusterIdx_list, length))

  # Flatten leafGenoInfo columns across clusters
  sgi <- objNull$subGenoInfo
  leaf_genoIndex  <- as.numeric(do.call(c, lapply(sgi, function(df) df$genoIndex)))
  leaf_AF_ref     <- as.numeric(do.call(c, lapply(sgi, function(df) df$AF_ref)))
  leaf_AN_ref     <- as.numeric(do.call(c, lapply(sgi, function(df) df$AN_ref)))
  leaf_TPR        <- as.numeric(do.call(c, lapply(sgi, function(df) df$TPR)))
  leaf_sigma2     <- as.numeric(do.call(c, lapply(sgi, function(df) df$sigma2)))
  leaf_pvalue_bat <- as.numeric(do.call(c, lapply(sgi, function(df) df$pvalue_bat)))
  leaf_w_ext      <- as.numeric(do.call(c, lapply(sgi, function(df) df[["w.ext"]])))
  leaf_var_w0     <- as.numeric(do.call(c, lapply(sgi, function(df) df[["var.ratio.w0"]])))
  leaf_var_int    <- as.numeric(do.call(c, lapply(sgi, function(df) df[["var.ratio.int"]])))
  leaf_var_ext    <- as.numeric(do.call(c, lapply(sgi, function(df) df[["var.ratio.ext"]])))
  leaf_nrows      <- as.integer(sapply(sgi, nrow))

  runMarkerInCPP.LEAF(
    residuals_all   = residuals_all,
    residuals_lens  = residuals_lens,
    weights_all     = weights_all,
    weights_lens    = weights_lens,
    clusterIdx_all  = clusterIdx_all,
    clusterIdx_lens = clusterIdx_lens,
    cutoff          = control$cutoff,
    SPA_Cutoff      = control$SPA_Cutoff,
    nCluster        = nCluster,
    leaf_genoIndex  = leaf_genoIndex,
    leaf_AF_ref     = leaf_AF_ref,
    leaf_AN_ref     = leaf_AN_ref,
    leaf_TPR        = leaf_TPR,
    leaf_sigma2     = leaf_sigma2,
    leaf_pvalue_bat = leaf_pvalue_bat,
    leaf_w_ext      = leaf_w_ext,
    leaf_var_w0     = leaf_var_w0,
    leaf_var_int    = leaf_var_int,
    leaf_var_ext    = leaf_var_ext,
    leaf_nrows      = leaf_nrows,
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
