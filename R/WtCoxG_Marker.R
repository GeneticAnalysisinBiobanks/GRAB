# WtCoxG_Marker.R -- Unwrap WtCoxG null model and run marker engine

runMarker.WtCoxG <- function(objNull, control, bedFile, bimFile, famFile, OutputFile, nThreads) {

  gi <- objNull$mergeGenoInfo
  
  runMarkerInCPP.WtCoxG(
    R             = objNull$mresid,
    w             = objNull$weight,
    cutoff        = control$cutoff,
    SPA_Cutoff    = control$SPA_Cutoff,
    wt_genoIndex  = as.numeric(gi$genoIndex),
    wt_AF_ref     = as.numeric(gi$AF_ref),
    wt_AN_ref     = as.numeric(gi$AN_ref),
    wt_TPR        = as.numeric(gi$TPR),
    wt_sigma2     = as.numeric(gi$sigma2),
    wt_pvalue_bat = as.numeric(gi$pvalue_bat),
    wt_w_ext      = as.numeric(gi[["w.ext"]]),
    wt_var_w0     = as.numeric(gi[["var.ratio.w0"]]),
    wt_var_int    = as.numeric(gi[["var.ratio.int"]]),
    wt_var_ext    = as.numeric(gi[["var.ratio.ext"]]),
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
