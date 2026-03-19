# SPAmix_Marker.R -- Unwrap SPAmix null model and run marker engine

runMarker.SPAmix <- function(objNull, OutputFile, control, bedFile, bimFile, famFile) {
  
  resid  <- objNull$resid
  nPheno <- ncol(resid)
  ol     <- objNull$outLierList

  # Concatenate per-pheno position vectors; C++ derives resid values from full matrix
  posValue_all       <- as.integer(unlist(lapply(seq_len(nPheno), function(i) ol[[i]]$posValue)))
  posValue_lens      <- as.integer(sapply(seq_len(nPheno), function(i) length(ol[[i]]$posValue)))
  posOutlier_all     <- as.integer(unlist(lapply(seq_len(nPheno), function(i) ol[[i]]$posOutlier)))
  posOutlier_lens    <- as.integer(sapply(seq_len(nPheno), function(i) length(ol[[i]]$posOutlier)))
  posNonOutlier_all  <- as.integer(unlist(lapply(seq_len(nPheno), function(i) ol[[i]]$posNonOutlier)))
  posNonOutlier_lens <- as.integer(sapply(seq_len(nPheno), function(i) length(ol[[i]]$posNonOutlier)))

  runMarkerInCPP.SPAmix(
    resid               = resid,
    PCs                 = objNull$PCs,
    N                   = objNull$N,
    SPA_Cutoff          = control$SPA_Cutoff,
    nPheno              = nPheno,
    posValue_all        = posValue_all,
    posValue_lens       = posValue_lens,
    posOutlier_all      = posOutlier_all,
    posOutlier_lens     = posOutlier_lens,
    posNonOutlier_all   = posNonOutlier_all,
    posNonOutlier_lens  = posNonOutlier_lens,
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
