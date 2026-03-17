# GRAB_Marker4.R
# Single-entry-point for chunk-parallel marker association testing.

.log <- function (msg, ...) {
  message(
    sprintf("[INFO] %s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(msg, ...))
  )
}

.print_params <- function (title, params, control) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[INFO] %s %s ...", ts, title))
  for (nm in names(params)) {
    if (!is.null(params[[nm]])) {
      message(sprintf("    %s: %s", nm, paste(params[[nm]], collapse = ", ")))
    }
  }
  if (length(control) > 0) {
    message("    Control parameters:")
    for (nm in names(control)) {
      v <- control[[nm]]
      if (is.numeric(v)) {
        message(sprintf("      %s: %g", nm, v))
      } else if (is.logical(v)) {
        message(sprintf("      %s: %s", nm, as.character(v)))
      } else if (is.character(v)) {
        message(sprintf("      %s: %s", nm, v))
      } else {
        message(sprintf("      %s: %s", nm, paste(v, collapse = ", ")))
      }
    }
  }
}

.merge_defaults <- function (control, defaults) {
  if (is.null(defaults)) return(control)
  if (!is.null(control)) {
    for (nm in names(control)) defaults[[nm]] <- control[[nm]]
  }
  defaults
}

.check_geno_control <- function (control) {
  defaults <- list(
    AlleleOrder     = NULL,
    AllMarkers      = TRUE,
    IDsToIncludeFile    = NULL,
    IDsToExcludeFile    = NULL,
    RangesToIncludeFile = NULL,
    RangesToExcludeFile = NULL
  )
  control <- if (is.null(control)) defaults else .merge_defaults(control, defaults)

  if (!is.null(control$AlleleOrder) &&
      !control$AlleleOrder %in% c("ref-first", "alt-first")) {
    stop("control$AlleleOrder must be 'ref-first' or 'alt-first'.")
  }

  inc <- c("IDsToIncludeFile", "RangesToIncludeFile")
  exc <- c("IDsToExcludeFile", "RangesToExcludeFile")
  has_inc <- any(vapply(inc, function(x) !is.null(control[[x]]), FALSE))
  has_exc <- any(vapply(exc, function(x) !is.null(control[[x]]), FALSE))
  if (has_inc && has_exc) {
    stop("Cannot provide both include and exclude marker-selection files.")
  }

  for (ft in c(inc, exc)) {
    f <- control[[ft]]
    if (!is.null(f) && !file.exists(f)) stop("File not found: ", f)
  }

  if (has_inc || has_exc) control$AllMarkers <- FALSE
  control
}

.check_marker_control <- function (null_class, control, obj_null) {
  method_defaults <- switch(
    null_class,
    POLMM_NULL_Model    = list(ifOutGroup = FALSE),
    SPACox_NULL_Model   = list(pVal_covaAdj_Cutoff = 5e-05),
    SPAmix_NULL_Model   = list(dosage_option = "rounding_first"),
    SPAmixPlus_NULL_Model = list(
      dosage_option   = "rounding_first",
      afFilePrecision = "double"
    ),
    SPAGRM_NULL_Model   = list(zeta = 0, tol = 1e-5),
    SAGELD_NULL_Model   = list(SPA_Cutoff = 2, zeta = 0, tol = 1e-4),
    WtCoxG_NULL_Model   = list(cutoff = 0.05),
    SPAsqr_NULL_Model   = list(zeta = 0, tol = 1e-5),
    LEAF_NULL_Model     = list(cutoff = 0.05),
    list()
  )
  control <- .merge_defaults(control, method_defaults)

  if (null_class == "SPACox_NULL_Model") {
    if (!is.numeric(control$pVal_covaAdj_Cutoff) || control$pVal_covaAdj_Cutoff <= 0) {
      stop("control$pVal_covaAdj_Cutoff must be numeric > 0.")
    }
  }

  if (null_class %in% c("SPAmix_NULL_Model", "SPAmixPlus_NULL_Model")) {
    if (!control$dosage_option %in% c("rounding_first", "rounding_last")) {
      stop("control$dosage_option must be 'rounding_first' or 'rounding_last'.")
    }
  }

  if (null_class == "SPAmixPlus_NULL_Model") {
    if (!control$afFilePrecision %in% c("double", "single", "text")) {
      stop("control$afFilePrecision must be 'double', 'single', or 'text'.")
    }
    if (is.null(control$afFilePath)) {
      stop("control$afFilePath must be provided for SPAmixPlus.")
    }
  }

  if (null_class %in% c("SPAGRM_NULL_Model", "SAGELD_NULL_Model")) {
    maf_iv <- obj_null$MAF_interval
    if (length(maf_iv) > 1 && control$min_maf_marker <= min(maf_iv)) {
      stop("min_maf_marker is out of MAF_interval. Please reset.")
    }
  }

  if (null_class %in% c("WtCoxG_NULL_Model", "LEAF_NULL_Model")) {
    lo <- if (null_class == "LEAF_NULL_Model") 0 else 0
    hi <- if (null_class == "LEAF_NULL_Model") 1 else 1
    strict_hi <- null_class == "LEAF_NULL_Model"
    if (!is.numeric(control$cutoff) || control$cutoff <= lo ||
        (strict_hi && control$cutoff >= hi) ||
        (!strict_hi && control$cutoff > hi)) {
      rng <- if (strict_hi) "(0, 1)" else "(0, 1]"
      stop(sprintf("control$cutoff must be numeric in %s.", rng))
    }
  }

  control
}

#' @export
GRAB.Marker4 <- function (
    objNull,
    GenoFile,
    OutputFile,
    GenoFileIndex = NULL,
    OutputFileIndex = NULL,
    control = NULL,
    nthreads = NULL,
    overwrite = FALSE
) {
  supported <- c(
    "POLMM_NULL_Model", "SPACox_NULL_Model", "SPAmix_NULL_Model",
    "SPAmixPlus_NULL_Model", "SPAGRM_NULL_Model", "SAGELD_NULL_Model",
    "WtCoxG_NULL_Model", "SPAsqr_NULL_Model", "LEAF_NULL_Model"
  )

  null_class <- class(objNull)
  if (!null_class %in% supported) {
    stop("class(objNull) must be one of: ", paste0('"', supported, '"', collapse = ", "))
  }
  if (!"subjData" %in% names(objNull)) stop("objNull must contain 'subjData'.")
  if (!is.null(control) && !is.list(control)) stop("'control' must be a list.")

  control <- .check_geno_control(control)

  marker_defaults <- list(
    impute_method      = "mean",
    missing_cutoff     = 0.15,
    min_maf_marker     = 0.001,
    min_mac_marker     = 20,
    nMarkersEachChunk  = 100,
    SPA_Cutoff         = 2,
    gzip_output        = FALSE
  )
  control <- .merge_defaults(control, marker_defaults)

  if (!is.null(control$nthreads)) {
    if (is.null(nthreads)) {
      nthreads <- control$nthreads
    } else if (!identical(as.integer(nthreads), as.integer(control$nthreads))) {
      .log("Both argument 'nthreads' and control$nthreads set; using argument.")
    }
    control$nthreads <- NULL
  }
  if (!is.null(control$omp_num_threads)) {
    stop("control$omp_num_threads is not used. Use argument 'nthreads'.")
  }

  n_threads <- if (is.null(nthreads)) {
    as.integer(data.table::getDTthreads())
  } else {
    if (!is.numeric(nthreads) || nthreads < 1 || (nthreads %% 1) != 0) {
      stop("'nthreads' must be a positive integer.")
    }
    as.integer(nthreads)
  }

  if (!is.logical(overwrite) || length(overwrite) != 1 || is.na(overwrite)) {
    stop("'overwrite' must be TRUE or FALSE.")
  }
  if (!control$impute_method %in% c("mean", "minor", "drop")) {
    stop("control$impute_method must be 'mean', 'minor', or 'drop'.")
  }
  if (!is.numeric(control$missing_cutoff) ||
      control$missing_cutoff < 0 || control$missing_cutoff > 0.5) {
    stop("control$missing_cutoff must be numeric in [0, 0.5].")
  }
  if (!is.numeric(control$min_maf_marker) ||
      control$min_maf_marker < 0 || control$min_maf_marker > 0.1) {
    stop("control$min_maf_marker must be numeric in [0, 0.1].")
  }
  if (!is.numeric(control$min_mac_marker) ||
      control$min_mac_marker < 0 || control$min_mac_marker > 100) {
    stop("control$min_mac_marker must be numeric in [0, 100].")
  }
  if (!is.numeric(control$nMarkersEachChunk) ||
      control$nMarkersEachChunk < 1e2 || control$nMarkersEachChunk > 1e5) {
    stop("control$nMarkersEachChunk must be numeric in [100, 100000].")
  }
  if (!is.numeric(control$SPA_Cutoff) || control$SPA_Cutoff <= 0) {
    stop("control$SPA_Cutoff must be numeric > 0.")
  }

  control <- .check_marker_control(null_class, control, objNull)

  single_only <- c(
    "SPACox_NULL_Model", "SPAmixPlus_NULL_Model",
    "SAGELD_NULL_Model", "WtCoxG_NULL_Model", "LEAF_NULL_Model"
  )
  if (null_class %in% single_only && n_threads > 1) {
    .log("%s currently supports only nthreads=1; forcing nthreads=1.", null_class)
    n_threads <- 1L
  }

  if (isTRUE(control$gzip_output) && !grepl("\\.gz$", OutputFile)) {
    OutputFile <- paste0(OutputFile, ".gz")
  }

  if (file.exists(OutputFile)) {
    if (overwrite) {
      if (!file.remove(OutputFile)) stop("Cannot remove existing OutputFile: ", OutputFile)
      .log("Existing output file removed (overwrite = TRUE).")
    } else {
      stop("OutputFile exists. Set overwrite=TRUE or use another path.")
    }
  }

  .print_params(
    "Parameters for Marker-Level Tests (GRAB.Marker4)",
    list(
      Method               = null_class,
      `Genotype file`      = GenoFile,
      `Genotype index file` = if (is.null(GenoFileIndex)) "Default" else GenoFileIndex,
      `Output file`        = OutputFile
    ),
    control
  )

  # Resolve PLINK file paths from GenoFile + optional GenoFileIndex
  bedFile <- GenoFile
  if (!is.null(GenoFileIndex) && length(GenoFileIndex) >= 2) {
    bimFile <- GenoFileIndex[1]
    famFile <- GenoFileIndex[2]
  } else {
    bimFile <- sub("\\.[^.]+$", ".bim", GenoFile)
    famFile <- sub("\\.[^.]+$", ".fam", GenoFile)
  }
  for (f in c(bedFile, bimFile, famFile)) {
    if (!file.exists(f)) stop("PLINK file not found: ", f)
  }

  # Unwrap null model and run marker analysis in C++
  runFn <- switch(null_class,
    POLMM_NULL_Model      = runMarker.POLMM,
    SPACox_NULL_Model     = runMarker.SPACox,
    SPAmix_NULL_Model     = runMarker.SPAmix,
    SPAmixPlus_NULL_Model = runMarker.SPAmixPlus,
    SPAGRM_NULL_Model     = runMarker.SPAGRM,
    SAGELD_NULL_Model     = runMarker.SAGELD,
    WtCoxG_NULL_Model     = runMarker.WtCoxG,
    SPAsqr_NULL_Model     = runMarker.SPAsqr,
    LEAF_NULL_Model       = runMarker.LEAF,
    stop("Unsupported null model class: ", null_class)
  )
  runFn(objNull, control, bedFile, bimFile, famFile, OutputFile, n_threads)

  .log("Analysis complete. Results saved to '%s'.", OutputFile)
  invisible(NULL)
}
