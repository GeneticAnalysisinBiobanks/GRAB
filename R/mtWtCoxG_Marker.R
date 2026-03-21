# WtCoxG_Marker.R -- Unwrap WtCoxG null model and run marker engine

mtMarker.WtCoxG <- function(objNull, OutputFile, control, bedFile, bimFile, famFile) {

  gi <- objNull$mergeGenoInfo
  
  mtMarkerInCPP.WtCoxG(
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
