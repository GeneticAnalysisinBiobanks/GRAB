## ------------------------------------------------------------------------------
## GRAB_Marker4.R
##
## Chunk-parallel marker analysis in C++ threads.
## ------------------------------------------------------------------------------

#' Perform chunk-parallel single-marker association tests
#'
#' Similar to \code{GRAB.Marker}, but chunk dispatch, buffering, and file writing are
#' executed fully in C++ with multi-threading.
#'
#' @inheritParams GRAB.Marker
#' @param nthreads Integer or NULL. Number of worker threads for chunk-level C++
#'   multithreading. If \code{NULL}, uses \code{data.table::getDTthreads()}.
#' @param overwrite Logical. If TRUE, overwrite existing \code{OutputFile}.
#'
#' @return Invisible \code{NULL}. Results are written to \code{OutputFile}.
#' @export
GRAB.Marker4 <- function(
  objNull,
  GenoFile,
  OutputFile,
  GenoFileIndex = NULL,
  OutputFileIndex = NULL,
  control = NULL,
  nthreads = NULL,
  overwrite = FALSE
) {

  supported_classes <- c(
    "POLMM_NULL_Model",
    "SPACox_NULL_Model",
    "SPAmix_NULL_Model",
    "SPAmixPlus_NULL_Model",
    "SPAGRM_NULL_Model",
    "SAGELD_NULL_Model",
    "WtCoxG_NULL_Model",
    "SPAsqr_NULL_Model",
    "LEAF_NULL_Model"
  )

  NullModelClass <- class(objNull)
  if (!NullModelClass %in% supported_classes) {
    stop(
      "class(objNull) should be one of: ",
      paste(paste0('"', supported_classes, '"'), collapse = ", ")
    )
  }

  if (any(!c("subjData") %in% names(objNull))) {
    stop("c('subjData') should be in names(objNull).")
  }

  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be an R list.")
  }

  control <- checkControl.ReadGeno(control)

  default.marker.control <- list(
    impute_method = "mean",
    missing_cutoff = 0.15,
    min_maf_marker = 0.001,
    min_mac_marker = 20,
    nMarkersEachChunk = 10000,
    SPA_Cutoff = 2
  )

  control <- updateControl(control, default.marker.control)

  if (!is.null(control$omp_num_threads)) {
    stop("control$omp_num_threads is not used in GRAB.Marker4. Use argument 'nthreads' instead.")
  }

  nThreads <- if (is.null(nthreads)) {
    as.integer(data.table::getDTthreads())
  } else {
    if (!is.numeric(nthreads) || nthreads < 1 || (nthreads %% 1) != 0) {
      stop("Argument 'nthreads' should be a positive integer.")
    }
    as.integer(nthreads)
  }

  if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
    stop("Argument 'overwrite' should be TRUE or FALSE.")
  }

  if (!control$impute_method %in% c("mean", "minor", "drop")) {
    stop("control$impute_method should be 'mean', 'minor', or 'drop'.")
  }

  if (!is.numeric(control$missing_cutoff) ||
      control$missing_cutoff < 0 || control$missing_cutoff > 0.5) {
    stop("control$missing_cutoff should be numeric in [0, 0.5].")
  }

  if (!is.numeric(control$min_maf_marker) ||
      control$min_maf_marker < 0 || control$min_maf_marker > 0.1) {
    stop("control$min_maf_marker should be numeric in [0, 0.1].")
  }

  if (!is.numeric(control$min_mac_marker) ||
      control$min_mac_marker < 0 || control$min_mac_marker > 100) {
    stop("control$min_mac_marker should be numeric in [0, 100].")
  }

  if (!is.numeric(control$nMarkersEachChunk) ||
      control$nMarkersEachChunk < 1e2 || control$nMarkersEachChunk > 1e5) {
    stop("control$nMarkersEachChunk should be numeric in [100, 100000].")
  }

  if (!is.numeric(control$SPA_Cutoff) || control$SPA_Cutoff <= 0) {
    stop("control$SPA_Cutoff should be a numeric value > 0.")
  }

  control <- switch(
    NullModelClass,
    POLMM_NULL_Model = checkControl.Marker.POLMM(control),
    SPACox_NULL_Model = checkControl.Marker.SPACox(control),
    SPAmix_NULL_Model = checkControl.Marker.SPAmix(control),
    SPAmixPlus_NULL_Model = checkControl.Marker.SPAmixPlus(control),
    SPAGRM_NULL_Model = checkControl.Marker.SPAGRM(control, objNull$MAF_interval),
    SAGELD_NULL_Model = checkControl.Marker.SAGELD(control, objNull$MAF_interval),
    WtCoxG_NULL_Model = checkControl.Marker.WtCoxG(control),
    SPAsqr_NULL_Model = checkControl.Marker.SPAsqr(control),
    LEAF_NULL_Model = checkControl.Marker.LEAF(control)
  )

  # Some methods currently rely on components that are not thread-safe in the
  # Marker4 chunk-parallel backend. Fall back to one thread to avoid crashes.
  unsafe_parallel_classes <- c(
    "WtCoxG_NULL_Model",
    "LEAF_NULL_Model"
  )
  if (NullModelClass %in% unsafe_parallel_classes && nThreads > 1) {
    .message(
      "Method %s currently uses thread-unsafe internals in GRAB.Marker4; forcing nthreads=1.",
      NullModelClass
    )
    nThreads <- 1L
  }

  if (file.exists(OutputFile)) {
    if (overwrite) {
      ok <- file.remove(OutputFile)
      if (!ok) {
        stop("Failed to remove existing OutputFile: ", OutputFile)
      }
      .message("Existing output file removed because overwrite = TRUE.")
    } else {
      stop("'OutputFile' exists. Set overwrite = TRUE to replace it, or use another output file.")
    }
  }

  if (!is.null(OutputFileIndex)) {
    .message("OutputFileIndex is ignored in GRAB.Marker4.")
  }

  params <- list(
    Method = NullModelClass,
    `Genotype file` = GenoFile,
    `Genotype index file` = ifelse(is.null(GenoFileIndex), "Default", GenoFileIndex),
    `Output file` = OutputFile
  )
  .printParameters("Parameters for Marker-Level Tests (GRAB.Marker4)", params, control)

  subjData <- as.character(objNull$subjData)
  objGeno <- setGenoInput(GenoFile, GenoFileIndex, subjData, control)
  genoType <- objGeno$genoType
  markerInfo <- objGeno$markerInfo

  CHROM <- markerInfo$CHROM
  uCHROM <- unique(CHROM)
  genoIndex <- markerInfo$genoIndex
  nMarkersEachChunk <- control$nMarkersEachChunk

  iTot <- 1
  genoIndexList <- list()
  for (chrom in uCHROM) {
    pos <- which(CHROM == chrom)
    gIdx <- genoIndex[pos]
    M <- length(gIdx)

    idxStart <- seq(1, M, nMarkersEachChunk)
    idxEnd <- idxStart + nMarkersEachChunk - 1
    nChunks <- length(idxStart)
    idxEnd[nChunks] <- M

    for (i in 1:nChunks) {
      idxMarker <- idxStart[i]:idxEnd[i]
      genoIndexList[[iTot]] <- as.integer(gIdx[idxMarker])
      iTot <- iTot + 1
    }
  }

  nChunks <- length(genoIndexList)
  .message("Number of markers to test: %d", nrow(markerInfo))
  .message("Genotype type: %s", genoType)
  .message("Number of markers in each chunk: %d", nMarkersEachChunk)
  .message("Number of chunks for all markers: %d", nChunks)
  .message("Number of threads: %d", nThreads)

  switch(
    NullModelClass,
    POLMM_NULL_Model = setMarker.POLMM(objNull, control),
    SPACox_NULL_Model = setMarker.SPACox(objNull, control),
    SPAmix_NULL_Model = setMarker.SPAmix(objNull, control),
    SPAmixPlus_NULL_Model = setMarker.SPAmixPlus(objNull, control),
    SPAGRM_NULL_Model = setMarker.SPAGRM(objNull, control),
    SAGELD_NULL_Model = setMarker.SAGELD(objNull, control),
    WtCoxG_NULL_Model = setMarker.WtCoxG(objNull, control),
    SPAsqr_NULL_Model = setMarker.SPAsqr(objNull, control),
    LEAF_NULL_Model = setMarker.LEAF(objNull, control)
  )

  readerConfig <- if (genoType == "PLINK") {
    list(
      genoType = "PLINK",
      bimFile = objGeno$GenoFileIndex[1],
      famFile = objGeno$GenoFileIndex[2],
      bedFile = objGeno$GenoFile,
      sampleInModel = as.character(objGeno$SampleIDs),
      alleleOrder = objGeno$AlleleOrder
    )
  } else {
    sampleFile <- objGeno$GenoFileIndex[2]
    if (!is.na(sampleFile) && file.exists(sampleFile)) {
      sampleData <- data.table::fread(sampleFile, header = TRUE)
      samplesInBgen <- as.character(sampleData$ID_2[-1])
    } else {
      samplesInBgen <- getSampleIDsFromBGEN(objGeno$GenoFile)
    }

    list(
      genoType = "BGEN",
      bgenFile = objGeno$GenoFile,
      bgiFile = objGeno$GenoFileIndex[1],
      sampleInBgen = as.character(samplesInBgen),
      sampleInModel = as.character(objGeno$SampleIDs),
      isSparseDosageInBgen = FALSE,
      isDropmissingdosagesInBgen = FALSE,
      alleleOrder = objGeno$AlleleOrder
    )
  }

  method <- switch(
    NullModelClass,
    POLMM_NULL_Model = "POLMM",
    SPACox_NULL_Model = "SPACox",
    SPAmix_NULL_Model = "SPAmix",
    SPAmixPlus_NULL_Model = "SPAmixPlus",
    SPAGRM_NULL_Model = "SPAGRM",
    SAGELD_NULL_Model = "SAGELD",
    WtCoxG_NULL_Model = "WtCoxG",
    SPAsqr_NULL_Model = "SPAsqr",
    LEAF_NULL_Model = "LEAF"
  )

  extraParams <- list(
    reader_config = readerConfig,
    sageld_method = if (NullModelClass == "SAGELD_NULL_Model") objNull$Method else NA_character_,
    spasqr_taus = if (NullModelClass == "SPAsqr_NULL_Model") objNull$taus else numeric(0),
    wtcoxg_merge = if (NullModelClass == "WtCoxG_NULL_Model") objNull$mergeGenoInfo else NULL,
    leaf_subgeno = if (NullModelClass == "LEAF_NULL_Model") objNull$subGenoInfo else NULL,
    leaf_ncluster = if (NullModelClass == "LEAF_NULL_Model") objNull$Ncluster else 0
  )

  mainMarkerChunksInCPP4(
    t_method = method,
    t_chunkIndexList = genoIndexList,
    t_outputFile = OutputFile,
    t_nThreads = nThreads,
    t_impute_method = control$impute_method,
    t_missing_cutoff = control$missing_cutoff,
    t_min_maf_marker = control$min_maf_marker,
    t_min_mac_marker = control$min_mac_marker,
    t_extraParams = extraParams
  )

  .message("Analysis complete! Results saved to '%s'", OutputFile)
  invisible(NULL)
}
