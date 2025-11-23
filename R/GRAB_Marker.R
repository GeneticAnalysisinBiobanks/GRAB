## ------------------------------------------------------------------------------
## GRAB_Marker.R
## Single-variant (marker-level) association testing orchestrator. Handles:
##   * Control validation & defaults
##   * Output/restart file management
##   * Chunked genotype processing for scalability
##   * Optional group-wise allele frequency/count summaries
##   * Dispatch to method-specific C++ backends (POLMM, SPACox, SPAmix, etc.)
##
## Functions:
##   GRAB.Marker: High-level wrapper to perform marker-level tests given a
##                fitted null model and genotype source.
## ------------------------------------------------------------------------------

#' Perform single-marker association tests
#'
#' Conducts single-marker association tests between genetic variants and phenotypes using
#' various statistical methods supported by GRAB.
#'
#' @param objNull (S3 object) Null model object from \code{\link{GRAB.NullModel}}.
#'   Must be one of: POLMM_NULL_Model, SPACox_NULL_Model, SPAmix_NULL_Model,
#'   SPAGRM_NULL_Model, SAGELD_NULL_Model, or WtCoxG_NULL_Model.
#' @param GenoFile (character) Path to genotype file. Supports PLINK (.bed/.bim/.fam)
#'   and BGEN formats. See \code{\link{GRAB.ReadGeno}} for format details.
#' @param OutputFile (character) Path for saving association test results.
#' @param GenoFileIndex (character or NULL) Index files for the genotype file. If
#'   \code{NULL} (default), uses the same prefix as \code{GenoFile}. See
#'   \code{\link{GRAB.ReadGeno}} for details.
#' @param OutputFileIndex (character or NULL) Path for progress tracking file. Enables
#'   analysis restart if interrupted. If \code{NULL} (default), uses
#'   \code{paste0(OutputFile, ".index")}.
#' @param control (list or NULL) List of control parameters with the following elements:
#' \itemize{
#'   \item \code{BgenImputeMethod} (character): Imputation method for BGEN-specific probabilistic
#'     genotype data during file reading. Options: "none" (default), "bestguess", "mean".
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
#'   \item \code{omp_num_threads} (integer): Number of OpenMP threads for parallel processing.
#'     Default: \code{data.table::getDTthreads()}.
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
GRAB.Marker <- function(
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
    "SPAGRM_NULL_Model",
    "SAGELD_NULL_Model",
    "WtCoxG_NULL_Model"
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

  if (any(!c("subjData", "N") %in% names(objNull))) {
    stop("c('subjData', 'N') should be in names(objNull).")
  }

  # OutputFile and OutputFileIndex will be further validated in checkOutputFile()
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  # ========== Validate and configure the control list ==========

  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be an R list.")
  }

  # Validate genotype-reading parameters: BgenImputeMethod, AlleleOrder, AllMarkers,
  # IDsToIncludeFile, IDsToExcludeFile, RangesToIncludeFile, RangesToExcludeFile
  control <- checkControl.ReadGeno(control)                   # list

  default.marker.control <- list(                             # list
    impute_method = "mean",
    missing_cutoff = 0.15,
    min_maf_marker = 0.001,
    min_mac_marker = 20,
    nMarkersEachChunk = 10000,
    omp_num_threads = data.table::getDTthreads(),
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
        control$nMarkersEachChunk < 1e3 || control$nMarkersEachChunk > 1e5) {
    stop("control$nMarkersEachChunk should be numeric in [1000, 100000].")
  }

  if (control$omp_num_threads < 0) {
    stop("control$omp_num_threads should be a positive integer.")
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
    SPAGRM_NULL_Model = checkControl.Marker.SPAGRM(control, objNull$MAF_interval),
    SAGELD_NULL_Model = checkControl.Marker.SAGELD(control, objNull$MAF_interval),
    WtCoxG_NULL_Model = checkControl.Marker.WtCoxG(control)
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

  # ========== Iterate over chunks to perform tests and write to OutputFile ==========

  # Set global objects in C++
  setMarker_GlobalVarsInCPP(
    t_impute_method = control$impute_method,      # character: "mean", "minor", or "drop"
    t_missing_cutoff = control$missing_cutoff,    # numeric: Max missing rate for markers
    t_min_maf_marker = control$min_maf_marker,    # numeric: Min MAF threshold
    t_min_mac_marker = control$min_mac_marker,    # numeric: Min MAC threshold
    t_omp_num_threads = control$omp_num_threads   # integer: Number of OpenMP threads
  )

  # Set method-specific objects in C++
  switch(
    NullModelClass,
    POLMM_NULL_Model  = setMarker.POLMM(objNull, control),
    SPACox_NULL_Model = setMarker.SPACox(objNull, control),
    SPAmix_NULL_Model = setMarker.SPAmix(objNull, control),
    SPAGRM_NULL_Model = setMarker.SPAGRM(objNull, control),
    SAGELD_NULL_Model = setMarker.SAGELD(objNull, control),
    WtCoxG_NULL_Model = setMarker.WtCoxG(objNull, control)
  )

  for (i in (indexChunk + 1):nChunks) {
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
      SPAGRM_NULL_Model = mainMarker.SPAGRM(genoType, genoIndex),
      SAGELD_NULL_Model = mainMarker.SAGELD(genoType, genoIndex, objNull),
      WtCoxG_NULL_Model = mainMarker.WtCoxG(genoType, genoIndex, objNull)
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

  .message("Analysis complete! Results saved to '%s'", OutputFile)
  return(invisible(NULL))
}