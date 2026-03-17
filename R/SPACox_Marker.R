# SPACox_Marker.R -- Unwrap SPACox null model and run marker engine

runMarker.SPACox <- function(objNull, control, bedFile, bimFile, famFile, OutputFile, nThreads) {
  
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
