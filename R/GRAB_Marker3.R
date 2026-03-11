## ------------------------------------------------------------------------------
## GRAB_Marker.R
##
## Functions:
##   GRAB.Marker: High-level wrapper to perform marker-level tests given a
##                fitted null model and genotype source.
## ------------------------------------------------------------------------------

#' Perform single-marker association tests using a fitted null model
#'
#' Conducts single-marker association tests between genetic variants and phenotypes using
#' various statistical methods supported by GRAB.
#'
#' @param objNull (S3 object) Null model object from \code{\link{GRAB.NullModel}},
#'   \code{\link{SPAGRM.NullModel}} or \code{\link{SAGELD.NullModel}}. Supported classes:
#'   \itemize{
#'     \item \code{POLMM_NULL_Model}: See \code{?\link{GRAB.POLMM}}.
#'     \item \code{SPACox_NULL_Model}: See \code{?\link{GRAB.SPACox}}.
#'     \item \code{SPAmix_NULL_Model}: See \code{?\link{GRAB.SPAmix}}.
#'     \item \code{WtCoxG_NULL_Model}: See \code{?\link{GRAB.WtCoxG}}.
#'     \item \code{SPAGRM_NULL_Model}: See \code{?\link{GRAB.SPAGRM}}.
#'     \item \code{SAGELD_NULL_Model}: See \code{?\link{GRAB.SAGELD}}.
#'   }
#' @param GenoFile Path to genotype file. Supported formats determined by extension:
#'   \itemize{
#'     \item PLINK: "prefix.bed"
#'     \item BGEN: "prefix.bgen" (version 1.2 with 8-bit compression)
#'   }
#' @param OutputFile (character) Path for saving association test results.
#' @param GenoFileIndex (character vector or NULL) Associated files for the genotype file (auto-detected if NULL):
#'   \itemize{
#'     \item PLINK: c("prefix.bim", "prefix.fam")
#'     \item BGEN: c("prefix.bgen.bgi", "prefix.sample") or c("prefix.bgen.bgi")
#'   }
#' @param OutputFileIndex (character or NULL) #' Path to the progress tracking file from a previous unfinished run.
#'   Enables analysis to restart if interrupted. If \code{NULL} (default), uses \code{paste0(OutputFile, ".index")}.
#' @param control (list or NULL) List of control parameters with the following elements:
#' \itemize{
#'   \item \code{AlleleOrder} (character or NULL): Allele order in genotype file. Options: "ref-first",
#'     "alt-first", or NULL (default: "alt-first" for BGEN, "ref-first" for PLINK).
#'   \item \strong{Marker Selection:}
#'   \itemize{
#'     \item \code{AllMarkers} (logical): Set to TRUE (default) to analyze all markers. Automatically
#'       set to FALSE if any include/exclude files are provided.
#'     \item \code{IDsToIncludeFile} (character or NULL): Path to file with marker IDs to include.
#'     \item \code{RangesToIncludeFile} (character or NULL): Path to file with genomic ranges to include.
#'       Can be used with IDsToIncludeFile (union will be used).
#'     \item \code{IDsToExcludeFile} (character or NULL): Path to file with marker IDs to exclude.
#'     \item \code{RangesToExcludeFile} (character or NULL): Path to file with genomic ranges to exclude.
#'       Can be used with IDsToExcludeFile (union will be excluded).
#'     \item Note: Cannot use both include and exclude files simultaneously.
#'   }
#'   \item \code{impute_method} (character): Imputation method for handling missing genotypes
#'     during analysis in C++ backend. Applies to all genotype formats. Options: "mean" (default), "minor", "drop".
#'   \item \code{missing_cutoff} (numeric): Exclude markers with missing rate above this threshold.
#'     Range: 0 to 0.5. Default: 0.15.
#'   \item \code{min_maf_marker} (numeric): Exclude markers with MAF below this threshold.
#'     Range: 0 to 0.1. Default: 0.001.
#'   \item \code{min_mac_marker} (numeric): Exclude markers with MAC below this threshold.
#'     Range: 0 to 100. Default: 20.
#'   \item \code{nMarkersEachChunk} (integer): Number of markers processed per chunk.
#'     Range: 1000 to 100000. Default: 10000.
#'   \item \code{parallel_ncores} (integer): Number of worker processes for chunk-level
#'     parallelization. Default: 1 (serial).
#'   \item \code{write_buffer_size} (integer or NULL): Pre-allocation hint for the ordered
#'     write buffer. The buffer always grows as needed to hold any number of out-of-order
#'     completed chunks while waiting for the next in-order chunk to arrive. With dynamic
#'     dispatch, a minimum of \code{parallel_ncores - 1} slots can be out-of-order at any
#'     instant, but a single stalled chunk lets later ones accumulate without bound. Default:
#'     \code{NULL}, which auto-derives to \code{parallel_ncores}.
#'   \item \code{omp_num_threads} (deprecated): Ignored for marker analysis.
#'   \item \code{SPA_Cutoff} (numeric): Z-score cutoff for saddlepoint approximation. When the absolute
#'     value of the test statistic exceeds this cutoff, SPA is used to calculate more accurate p-values. Default: 2.
#' }
#'
#' @return
#' The function returns \code{NULL} invisibly. Results are written to \code{OutputFile}.
#' For method-specific examples and output columns and format, see:
#' \itemize{
#'   \item POLMM method: \code{\link{GRAB.POLMM}}
#'   \item SPACox method: \code{\link{GRAB.SPACox}}
#'   \item SPAmix method: \code{\link{GRAB.SPAmix}}
#'   \item WtCoxG method: \code{\link{GRAB.WtCoxG}}
#'   \item SPAGRM method: \code{\link{GRAB.SPAGRM}}
#'   \item SAGELD method: \code{\link{GRAB.SAGELD}}
#' }
#'
GRAB.Marker3 <- function(
  objNull,
  GenoFile,
  OutputFile,
  GenoFileIndex = NULL,
  OutputFileIndex = NULL,
  control = NULL
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

  # ========== Validate and configure parameters ==========

  # Validate objNull
  NullModelClass <- class(objNull)                            # character

  if (!NullModelClass %in% supported_classes) {
    stop(
      "class(objNull) should be one of: ",
      paste(paste0('"', supported_classes, '"'), collapse = ", ")
    )
  }

  if (any(!c("subjData") %in% names(objNull))) {
    stop("c('subjData') should be in names(objNull).")
  }

  # OutputFile and OutputFileIndex will be further validated in checkOutputFile()
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  # ========== Validate and configure the control list ==========

  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be an R list.")
  }

  # Validate genotype-reading parameters: impute_method, AlleleOrder, AllMarkers,
  # IDsToIncludeFile, IDsToExcludeFile, RangesToIncludeFile, RangesToExcludeFile
  control <- checkControl.ReadGeno(control)                   # list

  default.marker.control <- list(                             # list
    impute_method = "mean",
    missing_cutoff = 0.15,
    min_maf_marker = 0.001,
    min_mac_marker = 20,
    nMarkersEachChunk = 10000,
    parallel_ncores = 1,
    write_buffer_size = NULL,               # auto-derived from parallel_ncores below
    SPA_Cutoff = 2
  )

  control <- updateControl(control, default.marker.control)  # list

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
        control$nMarkersEachChunk < 100 || control$nMarkersEachChunk > 1e5) {
    stop("control$nMarkersEachChunk should be numeric in [100, 100000].")
  }

  if (!is.numeric(control$parallel_ncores) ||
        control$parallel_ncores < 1 ||
        (control$parallel_ncores %% 1) != 0) {
    stop("control$parallel_ncores should be a positive integer.")
  }

  if (is.null(control$write_buffer_size)) {
    # auto-derive: with dynamic dispatch, at most (parallel_ncores - 1) chunks
    # can be out-of-order at any instant, but a stalled chunk allows later
    # completed chunks to accumulate unboundedly. The buffer always grows as
    # needed; this is only a pre-allocation hint.
    control$write_buffer_size <- control$parallel_ncores
  } else if (!is.numeric(control$write_buffer_size) ||
               control$write_buffer_size < 1 ||
               (control$write_buffer_size %% 1) != 0) {
    stop("control$write_buffer_size should be a positive integer.")
  }

  if (!is.null(control$omp_num_threads)) {
    .message("control$omp_num_threads is deprecated and ignored in GRAB.Marker.")
  }

  if (!is.numeric(control$SPA_Cutoff) || control$SPA_Cutoff <= 0) {
    stop("control$SPA_Cutoff should be a numeric value > 0.")
  }

  # Validate method-specific control parameters
  control <- switch(                                          # list
    NullModelClass,
    POLMM_NULL_Model  = checkControl.Marker.POLMM(control),
    SPACox_NULL_Model = checkControl.Marker.SPACox(control),
    SPAmix_NULL_Model = checkControl.Marker.SPAmix(control),
    SPAmixPlus_NULL_Model = checkControl.Marker.SPAmixPlus(control),
    SPAGRM_NULL_Model = checkControl.Marker.SPAGRM(control, objNull$MAF_interval),
    SAGELD_NULL_Model = checkControl.Marker.SAGELD(control, objNull$MAF_interval),
    WtCoxG_NULL_Model = checkControl.Marker.WtCoxG(control),
    SPAsqr_NULL_Model = checkControl.Marker.SPAsqr(control),
    LEAF_NULL_Model = checkControl.Marker.LEAF(control)
  )

  # ========== Check output file status and determine restart point ==========

  nMarkersEachChunk <- control$nMarkersEachChunk              # numeric

  indexChunk <- checkOutputFile(                              # integer
    OutputFile, OutputFileIndex, "Marker",
    format(nMarkersEachChunk, scientific = FALSE)
  )

  # ========== Print all parameters ==========

  params <- list(
    Method = NullModelClass,
    `Genotype file` = GenoFile,
    `Genotype index file` = ifelse(is.null(GenoFileIndex), "Default", GenoFileIndex),
    `Output file` = OutputFile,
    `Output index file` = ifelse(is.null(OutputFileIndex), "Default", OutputFileIndex)
  )
  .printParameters("Parameters for Marker-Level Tests", params, control)

  # ========== Initialize genotype reader and create chunks ==========

  # Initialize genotype reader with file paths and subject filtering options
  subjData <- as.character(objNull$subjData)                  # character vector
  objGeno <- setGenoInput(GenoFile, GenoFileIndex, subjData, control) # list:
    # genoType, markerInfo, SampleIDs, AlleleOrder, GenoFile, GenoFileIndex, anyQueue
  genoType <- objGeno$genoType                                # character
  markerInfo <- objGeno$markerInfo                            # data.frame
  CHROM <- markerInfo$CHROM                                   # character vector
  uCHROM <- unique(CHROM)                                     # character vector
  genoIndex <- markerInfo$genoIndex                           # integer vector

  # Iterate over chromosomes to create chunks
  iTot <- 1                                                   # integer
  genoIndexList <- list()                                     # list: to hold chunks of marker indices
  for (chrom in uCHROM) {
    # Extract all markers belonging to current chromosome
    pos <- which(CHROM == chrom)                              # integer vector
    gIdx <- genoIndex[pos]                                    # integer vector
    M <- length(gIdx)                                         # integer

    # Calculate chunk boundaries within this chromosome
    idxStart <- seq(1, M, nMarkersEachChunk)                  # integer vector
    idxEnd <- idxStart + nMarkersEachChunk - 1                # integer vector

    nChunks <- length(idxStart)                               # integer
    idxEnd[nChunks] <- M                                      # Ensure last chunk includes all remaining markers

    # Create individual chunks for this chromosome
    for (i in 1:nChunks) {
      idxMarker <- idxStart[i]:idxEnd[i]                      # integer vector
      genoIndexList[[iTot]] <- list(                          # list
        chrom = chrom,                                        # Chromosome identifier
        genoIndex = gIdx[idxMarker]                           # Marker indices for this chunk
      )
      iTot <- iTot + 1
    }
  }

  nChunks <- length(genoIndexList)                            # integer

  .message("Number of markers to test: %d", nrow(markerInfo))
  .message("Number of markers in each chunk: %d", nMarkersEachChunk)
  .message("Number of chunks for all markers: %d", nChunks)

  if (control$parallel_ncores > 1) {
    .runMarkerChunkAnalysisParallel(
      NullModelClass = NullModelClass,
      objNull = objNull,
      control = control,
      GenoFile = GenoFile,
      GenoFileIndex = GenoFileIndex,
      subjData = subjData,
      genoIndexList = genoIndexList,
      indexChunk = indexChunk,
      nChunks = nChunks,
      OutputFile = OutputFile,
      OutputFileIndex = OutputFileIndex,
      nMarkersEachChunk = nMarkersEachChunk
    )
  } else {
    .runMarkerChunkAnalysis(
      NullModelClass = NullModelClass,
      objNull = objNull,
      control = control,
      genoType = genoType,
      genoIndexList = genoIndexList,
      indexChunk = indexChunk,
      nChunks = nChunks,
      OutputFile = OutputFile,
      OutputFileIndex = OutputFileIndex,
      nMarkersEachChunk = nMarkersEachChunk
    )
  }

  .message("Analysis complete! Results saved to '%s'", OutputFile)
  return(invisible(NULL))
}

# Internal helper to run marker-level chunk analysis in C++ backend.
.runMarkerChunkAnalysis <- function(
  NullModelClass,
  objNull,
  control,
  genoType,
  genoIndexList,
  indexChunk,
  nChunks,
  OutputFile,
  OutputFileIndex,
  nMarkersEachChunk
) {

  # Set global objects in C++
  setMarker_GlobalVarsInCPP(
    t_impute_method = control$impute_method,      # character: "mean", "minor", or "drop"
    t_missing_cutoff = control$missing_cutoff,    # numeric: Max missing rate for markers
    t_min_maf_marker = control$min_maf_marker,    # numeric: Min MAF threshold
    t_min_mac_marker = control$min_mac_marker,    # numeric: Min MAC threshold
    t_omp_num_threads = 1                          # deprecated in marker path; use process-level parallelism
  )

  # Set method-specific objects in C++
  switch(
    NullModelClass,
    POLMM_NULL_Model  = setMarker.POLMM(objNull, control),
    SPACox_NULL_Model = setMarker.SPACox(objNull, control),
    SPAmix_NULL_Model = setMarker.SPAmix(objNull, control),
    SPAmixPlus_NULL_Model = setMarker.SPAmixPlus(objNull, control),
    SPAGRM_NULL_Model = setMarker.SPAGRM(objNull, control),
    SAGELD_NULL_Model = setMarker.SAGELD(objNull, control),
    WtCoxG_NULL_Model = setMarker.WtCoxG(objNull, control),
    SPAsqr_NULL_Model = setMarker.SPAsqr(objNull, control),
    LEAF_NULL_Model = setMarker.LEAF(objNull, control)
  )

  for (i in seq.int(indexChunk + 1, nChunks)) {
    tempList <- genoIndexList[[i]]                            # list: chrom, genoIndex
    genoIndex <- tempList$genoIndex                           # integer vector
    tempChrom <- tempList$chrom                               # character

    .message("---- Analyzing Chunk %d/%d: chrom %s ----", i, nChunks, tempChrom)

    # Test one chunk in C++ backend
    resMarker <- switch(                                      # data.frame
      NullModelClass,
      POLMM_NULL_Model  = mainMarker.POLMM(genoType, genoIndex, control),
      SPACox_NULL_Model = mainMarker.SPACox(genoType, genoIndex),
      SPAmix_NULL_Model = mainMarker.SPAmix(genoType, genoIndex, objNull),
      SPAmixPlus_NULL_Model = mainMarker.SPAmixPlus(genoType, genoIndex, objNull, control),
      SPAGRM_NULL_Model = mainMarker.SPAGRM(genoType, genoIndex),
      SAGELD_NULL_Model = mainMarker.SAGELD(genoType, genoIndex, objNull),
      WtCoxG_NULL_Model = mainMarker.WtCoxG(genoType, genoIndex, objNull),
      SPAsqr_NULL_Model = mainMarker.SPAsqr(genoType, genoIndex, objNull),
      LEAF_NULL_Model = mainMarker.LEAF(genoType, genoIndex, objNull)
    )

    # Write chunk results to output file and update progress tracking
    writeOutputFile(
      Output = list(resMarker),
      OutputFile = list(OutputFile),
      OutputFileIndex = OutputFileIndex,
      AnalysisType = "Marker",
      nEachChunk = format(nMarkersEachChunk, scientific = FALSE),
      indexChunk = i,
      Start = (i == 1),
      End = (i == nChunks)
    )
  }

  return(invisible(NULL))
}


.createMarkerWriteBuffer <- function(indexChunk, nChunks, OutputFile, OutputFileIndex, nMarkersEachChunk) {
  buffer <- new.env(parent = emptyenv())
  buffer$pending <- list()                  # named list: chunkId -> resMarker; grows unboundedly
  buffer$nextChunk <- indexChunk + 1        # next chunk id expected for in-order writing
  buffer$nChunks <- nChunks
  buffer$OutputFile <- OutputFile
  buffer$OutputFileIndex <- OutputFileIndex
  buffer$nMarkersEachChunk <- nMarkersEachChunk
  buffer
}


.flushMarkerWriteBuffer <- function(buffer) {
  # Drain all consecutive completed chunks starting at buffer$nextChunk.
  # Stops as soon as the next expected chunk is not yet in the buffer.
  repeat {
    key <- as.character(buffer$nextChunk)
    if (is.null(buffer$pending[[key]])) break

    resMarker <- buffer$pending[[key]]
    buffer$pending[[key]] <- NULL

    writeOutputFile(
      Output = list(resMarker),
      OutputFile = list(buffer$OutputFile),
      OutputFileIndex = buffer$OutputFileIndex,
      AnalysisType = "Marker",
      nEachChunk = format(buffer$nMarkersEachChunk, scientific = FALSE),
      indexChunk = buffer$nextChunk,
      Start = (buffer$nextChunk == 1),
      End = (buffer$nextChunk == buffer$nChunks)
    )

    buffer$nextChunk <- buffer$nextChunk + 1
  }

  invisible(NULL)
}


.pushMarkerWriteBuffer <- function(buffer, chunkId, resMarker) {
  buffer$pending[[as.character(chunkId)]] <- resMarker
  .flushMarkerWriteBuffer(buffer)
  invisible(NULL)
}


.runMarkerChunkWorkerInit <- function(
  NullModelClass,
  objNull,
  control,
  GenoFile,
  GenoFileIndex,
  subjData
) {
  objGeno <- setGenoInput(
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    SampleIDs = subjData,
    control = control
  )

  genoType <- objGeno$genoType

  setMarker_GlobalVarsInCPP(
    t_impute_method = control$impute_method,
    t_missing_cutoff = control$missing_cutoff,
    t_min_maf_marker = control$min_maf_marker,
    t_min_mac_marker = control$min_mac_marker,
    t_omp_num_threads = 1
  )

  switch(
    NullModelClass,
    POLMM_NULL_Model  = setMarker.POLMM(objNull, control),
    SPACox_NULL_Model = setMarker.SPACox(objNull, control),
    SPAmix_NULL_Model = setMarker.SPAmix(objNull, control),
    SPAmixPlus_NULL_Model = setMarker.SPAmixPlus(objNull, control),
    SPAGRM_NULL_Model = setMarker.SPAGRM(objNull, control),
    SAGELD_NULL_Model = setMarker.SAGELD(objNull, control),
    WtCoxG_NULL_Model = setMarker.WtCoxG(objNull, control),
    SPAsqr_NULL_Model = setMarker.SPAsqr(objNull, control),
    LEAF_NULL_Model = setMarker.LEAF(objNull, control)
  )

  assign(
    x = ".markerWorkerState",
    value = list(
      NullModelClass = NullModelClass,
      objNull = objNull,
      control = control,
      genoType = genoType
    ),
    envir = .GlobalEnv
  )

  invisible(TRUE)
}


.runMarkerChunkWorker <- function(chunkId, chunkInfo) {
  state <- get(".markerWorkerState", envir = .GlobalEnv, inherits = FALSE)
  genoIndex <- chunkInfo$genoIndex

  resMarker <- switch(
    state$NullModelClass,
    POLMM_NULL_Model  = mainMarker.POLMM(state$genoType, genoIndex, state$control),
    SPACox_NULL_Model = mainMarker.SPACox(state$genoType, genoIndex),
    SPAmix_NULL_Model = mainMarker.SPAmix(state$genoType, genoIndex, state$objNull),
    SPAmixPlus_NULL_Model = mainMarker.SPAmixPlus(state$genoType, genoIndex, state$objNull, state$control),
    SPAGRM_NULL_Model = mainMarker.SPAGRM(state$genoType, genoIndex),
    SAGELD_NULL_Model = mainMarker.SAGELD(state$genoType, genoIndex, state$objNull),
    WtCoxG_NULL_Model = mainMarker.WtCoxG(state$genoType, genoIndex, state$objNull),
    SPAsqr_NULL_Model = mainMarker.SPAsqr(state$genoType, genoIndex, state$objNull),
    LEAF_NULL_Model = mainMarker.LEAF(state$genoType, genoIndex, state$objNull)
  )

  list(
    chunkId = chunkId,
    chrom = chunkInfo$chrom,
    resMarker = resMarker
  )
}


.runMarkerChunkAnalysisParallel <- function(
  NullModelClass,
  objNull,
  control,
  GenoFile,
  GenoFileIndex,
  subjData,
  genoIndexList,
  indexChunk,
  nChunks,
  OutputFile,
  OutputFileIndex,
  nMarkersEachChunk
) {
  if (indexChunk >= nChunks) {
    return(invisible(NULL))
  }

  remainingChunkIds <- seq.int(indexChunk + 1, nChunks)
  nWorkers <- min(control$parallel_ncores, length(remainingChunkIds))

  .message("Running chunk-level parallel analysis with %d workers", nWorkers)

  cl <- parallel::makeCluster(nWorkers)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(
    cl,
    varlist = c(
      "setGenoInput",
      "setMarker_GlobalVarsInCPP",
      "setMarker.POLMM", "setMarker.SPACox", "setMarker.SPAmix", "setMarker.SPAmixPlus",
      "setMarker.SPAGRM", "setMarker.SAGELD", "setMarker.WtCoxG", "setMarker.SPAsqr", "setMarker.LEAF",
      "mainMarker.POLMM", "mainMarker.SPACox", "mainMarker.SPAmix", "mainMarker.SPAmixPlus",
      "mainMarker.SPAGRM", "mainMarker.SAGELD", "mainMarker.WtCoxG", "mainMarker.SPAsqr", "mainMarker.LEAF",
      ".runMarkerChunkWorkerInit", ".runMarkerChunkWorker"
    ),
    envir = environment()
  )

  parallel::clusterCall(
    cl,
    fun = .runMarkerChunkWorkerInit,
    NullModelClass = NullModelClass,
    objNull = objNull,
    control = control,
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    subjData = subjData
  )

  writeBuffer <- .createMarkerWriteBuffer(
    indexChunk = indexChunk,
    nChunks = nChunks,
    OutputFile = OutputFile,
    OutputFileIndex = OutputFileIndex,
    nMarkersEachChunk = nMarkersEachChunk
  )

  submitted <- 0L
  received <- 0L
  nextPos <- 1L
  nRemain <- length(remainingChunkIds)

  while (nextPos <= nRemain && submitted < nWorkers) {
    chunkId <- remainingChunkIds[nextPos]
    workerId <- submitted + 1L
    parallel:::sendCall(
      cl[[workerId]],
      fun = .runMarkerChunkWorker,
      args = list(chunkId = chunkId, chunkInfo = genoIndexList[[chunkId]]),
      tag = as.integer(chunkId)
    )
    submitted <- submitted + 1L
    nextPos <- nextPos + 1L
  }

  while (received < nRemain) {
    recv <- parallel:::recvOneResult(cl)
    value <- recv$value
    chunkId <- value$chunkId

    .message("---- Completed Chunk %d/%d: chrom %s ----", chunkId, nChunks, value$chrom)

    .pushMarkerWriteBuffer(writeBuffer, chunkId, value$resMarker)
    received <- received + 1L

    if (nextPos <= nRemain) {
      nextChunkId <- remainingChunkIds[nextPos]
      parallel:::sendCall(
        cl[[recv$node]],
        fun = .runMarkerChunkWorker,
        args = list(chunkId = nextChunkId, chunkInfo = genoIndexList[[nextChunkId]]),
        tag = as.integer(nextChunkId)
      )
      nextPos <- nextPos + 1L
    }
  }

  .flushMarkerWriteBuffer(writeBuffer)

  invisible(NULL)
}