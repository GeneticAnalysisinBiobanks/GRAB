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
#' @param GenoFileIndex (character or NULL) Index files for the genotype file. If 
#'   \code{NULL} (default), uses the same prefix as \code{GenoFile}. See 
#'   \code{\link{GRAB.ReadGeno}} for details.
#' @param OutputFile (character) Path for saving association test results.
#' @param OutputFileIndex (character or NULL) Path for progress tracking file. Enables 
#'   analysis restart if interrupted. If \code{NULL} (default), uses 
#'   \code{paste0(OutputFile, ".index")}.
#' @param control (list or NULL) List of control parameters with the following elements:
#' \itemize{
#'   \item \code{ImputeMethod} (character): Imputation method for genotype data.
#'     Options: "none" (default), "bestguess", "mean".
#'   \item \code{AlleleOrder} (character or NULL): Allele order in genotype file. Options: "ref-first",
#'     "alt-first", or NULL (default, uses file-type default).
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
#'   \item \code{impute_method} (character): Imputation method in C++ backend.
#'     Options: "mean" (default), "minor", "drop".
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
#'   \item \code{outputColumns} (character vector or NULL): Additional columns to include in output.
#'     Available options depend on the method. Default: NULL.
#' }
#'
#' @return
#' Results are written to \code{OutputFile}. The function returns \code{NULL} invisibly.
#' For details on output columns and format, see the documentation for the specific method:
#' \code{\link{GRAB.POLMM}}, \code{\link{GRAB.SPACox}}, \code{\link{GRAB.SPAmix}},
#' \code{\link{GRAB.SPAGRM}}, \code{\link{GRAB.SAGELD}}, or \code{\link{GRAB.WtCoxG}}.
#'
#' @examples
#' objNullFile <- system.file("extdata", "objPOLMMnull.RData", package = "GRAB")
#' load(objNullFile) # load a an example object, obj.POLMM, from step 1
#'
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "simuOUTPUT.txt")
#' outputColumns <- c(
#'   "beta", "seBeta", "zScore",
#'   "nSamplesInGroup", "AltCountsInGroup", "AltFreqInGroup"
#' )
#'
#' GRAB.Marker(
#'   objNull = obj.POLMM,
#'   GenoFile = GenoFile,
#'   OutputFile = OutputFile,
#'   control = list(outputColumns = outputColumns)
#' )
#'
#' data.table::fread(OutputFile)
#'
GRAB.Marker <- function(
  objNull,
  GenoFile,
  GenoFileIndex = NULL,
  OutputFile,
  OutputFileIndex = NULL,
  control = NULL
) {

# ========== Validate and configure parameters ==========

  supported_classes <- c(                                      # character vector
    "POLMM_NULL_Model",
    "SPACox_NULL_Model",
    "SPAmix_NULL_Model",
    "SPAGRM_NULL_Model",
    "SAGELD_NULL_Model",
    "WtCoxG_NULL_Model"
  )

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

  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  # ========== Validate and configure the control list ==========

  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be an R list.")
  }

  # Validate genotype-reading parameters: ImputeMethod, AlleleOrder, AllMarkers,
  # IDsToIncludeFile, IDsToExcludeFile, RangesToIncludeFile, RangesToExcludeFile
  control <- checkControl.ReadGeno(control)                   # list

  default.marker.control <- list(                             # list
    impute_method = "mean",
    missing_cutoff = 0.15,
    min_maf_marker = 0.001,
    min_mac_marker = 20,
    nMarkersEachChunk = 10000,
    omp_num_threads = data.table::getDTthreads(),
    outputColumns = NULL
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

  # Validate method-specific parameters
  control <- switch(                                          # list
    NullModelClass,
    POLMM_NULL_Model  = checkControl.Marker.POLMM(control),
    SPACox_NULL_Model = checkControl.Marker.SPACox(control),
    SPAmix_NULL_Model = checkControl.Marker.SPAmix(control),
    SPAGRM_NULL_Model = checkControl.Marker.SPAGRM(control, objNull$MAF_interval),
    SAGELD_NULL_Model = checkControl.Marker.SAGELD(control, objNull$MAF_interval),
    WtCoxG_NULL_Model = checkControl.Marker.WtCoxG(control)
  )

  # Pretty-print final control list
  .message("Control parameters for marker-level association tests:")
  tmp <- capture.output(str(control))
  for (line in tmp[startsWith(tmp, " $")]) {
    message(sub("^ \\$", strrep(" ", 8), line))
  }

  # ========== Check output file status and determine restart point ==========

  nMarkersEachChunk <- control$nMarkersEachChunk              # numeric

  indexChunk <- checkOutputFile(                              # integer
    OutputFile, OutputFileIndex, "Marker",
    format(nMarkersEachChunk, scientific = FALSE)
  )

  # ========== Grouping phenotypic values ==========

  yVec <- objNull$yVec                                        # numeric vector
  m1 <- length(unique(yVec))                                  # integer
  if (m1 <= 10) {
    Group <- as.numeric(as.factor(yVec)) - 1                  # numeric vector: Groups 0 to (m1-1)
  } else {
    Group <- floor((rank(yVec, ties.method = "max") - 1) / length(yVec) * 10) # numeric vector: Groups 0 to 9
  }

  nGroup <- length(unique(Group))                             # integer
  ifOutGroup <- any(c("AltFreqInGroup", "AltCountsInGroup") %in% control$outputColumns) # logical

  # ========== Initialize genotype reader and create chunks ==========

  # Initialize genotype reader with file paths and subject filtering options
  subjData <- as.character(objNull$subjData)                  # character vector
  objGeno <- setGenoInput(GenoFile, GenoFileIndex, subjData, control) # list
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

  # Display analysis information to user
  .message("Number of markers to test: %d", nrow(markerInfo))
  .message("Number of markers in each chunk: %d", nMarkersEachChunk)
  .message("Number of chunks for all markers: %d", nChunks)
  
  # ========== Iterate over chunks to perform tests ==========

  chrom <- "InitialChunk"                                     # character
  for (i in (indexChunk + 1):nChunks) {
    tempList <- genoIndexList[[i]]                            # list
    genoIndex <- tempList$genoIndex                           # integer vector
    tempChrom <- tempList$chrom                               # character
    
    if (tempChrom != chrom) {

      # Set global objects in C++ backend when tempChrom changes
      setMarker_GlobalVarsInCPP(
        control$impute_method,
        control$missing_cutoff,
        control$min_maf_marker,
        control$min_mac_marker,
        control$omp_num_threads,
        Group, ifOutGroup, nGroup
      )

      # Set method-specific objects in C++ backend when tempChrom changes
      switch(
        NullModelClass,
        POLMM_NULL_Model  = setMarker.POLMM(objNull, control, chrom),
        SPACox_NULL_Model = setMarker.SPACox(objNull, control),
        SPAmix_NULL_Model = setMarker.SPAmix(objNull, control),
        SPAGRM_NULL_Model = setMarker.SPAGRM(objNull, control),
        SAGELD_NULL_Model = setMarker.SAGELD(objNull, control),
        WtCoxG_NULL_Model = setMarker.WtCoxG(objNull, control)
      )

      chrom <- tempChrom
    }

    .message("---- Analyzing Chunk %d/%d: chrom %s ----", i, nChunks, chrom)

    # Test one SNP in C++ backend
    resMarker <- switch(                                      # data.frame
      NullModelClass,
      POLMM_NULL_Model  = mainMarker.POLMM(genoType, genoIndex, control$outputColumns),
      SPACox_NULL_Model = mainMarker.SPACox(genoType, genoIndex, control$outputColumns),
      SPAmix_NULL_Model = mainMarker.SPAmix(genoType, genoIndex, control$outputColumns, objNull),
      SPAGRM_NULL_Model = mainMarker.SPAGRM(genoType, genoIndex),
      SAGELD_NULL_Model = mainMarker.SAGELD(genoType, genoIndex, objNull),
      WtCoxG_NULL_Model = mainMarker.WtCoxG(genoType, genoIndex, control, objNull)
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