# POLMM_Marker.R -- Unwrap POLMM null model and run marker engine

runMarker.POLMM <- function(objNull, control, bedFile, bimFile, famFile,
                            OutputFile, nThreads) {
  objCHR <- objNull$LOCOList[["LOCO=F"]]

  runMarkerInCPP.POLMM(
    muMat       = objCHR$muMat,
    iRMat       = objCHR$iRMat,
    Cova        = objNull$Cova,
    yVec        = objNull$yVec,
    tau         = objNull$tau,
    varRatio    = objCHR$VarRatio,
    SPA_Cutoff  = control$SPA_Cutoff,
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
