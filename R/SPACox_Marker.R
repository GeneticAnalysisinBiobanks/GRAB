# SPACox_Marker.R -- Unwrap SPACox null model and run marker engine

runMarker.SPACox <- function(objNull, OutputFile, control, bedFile, bimFile, famFile) {
  
  runMarkerInCPP.SPACox(
    cumul               = objNull$cumul,
    mresid              = objNull$mresid,
    XinvXX              = objNull[["X.invXX"]],
    tX                  = objNull$tX,
    N                   = length(objNull$mresid),
    pVal_covaAdj_Cutoff = control$pVal_covaAdj_Cutoff,
    SPA_Cutoff          = control$SPA_Cutoff,
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
