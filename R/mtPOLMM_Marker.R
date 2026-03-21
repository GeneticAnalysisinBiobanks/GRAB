# POLMM_Marker.R -- Unwrap POLMM null model and run marker engine

mtMarker.POLMM <- function(objNull, OutputFile, control, bedFile, bimFile, famFile) {
  
  objCHR <- objNull$LOCOList[["LOCO=F"]]

  mtMarkerInCPP.POLMM(
    muMat      = objCHR$muMat,
    iRMat      = objCHR$iRMat,
    Cova       = objNull$Cova,
    yVec       = objNull$yVec,
    varRatio   = objCHR$VarRatio,
    SPA_Cutoff = control$SPA_Cutoff,
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
