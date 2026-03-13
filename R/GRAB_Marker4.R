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
    nMarkersEachChunk = 1000,
    SPA_Cutoff = 2
  )

  control <- updateControl(control, default.marker.control)

  # Backward-compatible thread control: accept control$nthreads when the
  # dedicated function argument is not provided.
  if (!is.null(control$nthreads)) {
    if (is.null(nthreads)) {
      nthreads <- control$nthreads
    } else if (!identical(as.integer(nthreads), as.integer(control$nthreads))) {
      .message("Both argument 'nthreads' and control$nthreads are set; using argument 'nthreads'.")
    }
    # Remove to avoid confusion in printed control parameters.
    control$nthreads <- NULL
  }

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

  # Methods listed here are forced to run with one thread.
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

  prepareAndRunMarkerInCPP4(
    t_objNull = objNull,
    t_GenoFile = GenoFile,
    t_OutputFile = OutputFile,
    t_GenoFileIndex = GenoFileIndex,
    t_control = control,
    t_nThreads = nThreads
  )

  .message("Analysis complete! Results saved to '%s'", OutputFile)
  invisible(NULL)
}
