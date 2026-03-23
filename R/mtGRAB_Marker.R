# mtGRAB_Marker.R
# Single-entry-point for chunk-parallel marker association testing.

GRAB.Marker6 <- function (
    objNull,
    GenoFile,
    OutputFile,
    GenoFileIndex = NULL,
    overwrite = FALSE,
    control = list()
) {

  null_class <- class(objNull)
  supported_classes <- c(
    "POLMM_NULL_Model", "WtCoxG_NULL_Model", "LEAF_NULL_Model",
    "SPAGRM_NULL_Model", "SAGELD_NULL_Model", "SPAsqr_NULL_Model",
    "SPACox_NULL_Model", "SPAmix_NULL_Model", "SPAmixPlus_NULL_Model"
  )
  if (!null_class %in% supported_classes) {
    stop("Unsupported null model class: ", null_class)
  }

  common_defaults <- list(
    SPA_Cutoff          = 2,
    impute_method       = "mean",
    missing_cutoff      = 0.15,
    min_maf_marker      = 0.001,
    min_mac_marker      = 20,
    nMarkersEachChunk   = 1024,
    AlleleOrder         = "alt-first",
    hwe                 = "exact",
    AllMarkers          = TRUE,
    IDsToIncludeFile    = NULL,
    IDsToExcludeFile    = NULL,
    RangesToIncludeFile = NULL,
    RangesToExcludeFile = NULL,
    nthreads = min(as.numeric(system("nproc", intern = TRUE)) - 1, 8)
  )

  method_defaults <- switch(
    null_class,
    POLMM_NULL_Model      = list(),
    WtCoxG_NULL_Model     = list(cutoff = 0.05),
    LEAF_NULL_Model       = list(cutoff = 0.05),
    SPAGRM_NULL_Model     = list(zeta = 0, tol = 1e-5),
    SAGELD_NULL_Model     = list(zeta = 0, tol = 1e-4),
    SPAsqr_NULL_Model     = list(zeta = 0, tol = 1e-5),
    SPACox_NULL_Model     = list(pVal_covaAdj_Cutoff = 5e-05),
    SPAmix_NULL_Model     = list(dosage_option = "rounding_first"),
    SPAmixPlus_NULL_Model = list(dosage_option = "rounding_first", afFilePrecision = "double"),
    list()
  )

  # validate genotype file
  if (!is.character(GenoFile) || length(GenoFile) != 1) {
    stop("'GenoFile' must be a single string (PLINK bed file).")
  }
  ext <- tools::file_ext(GenoFile)
  if (ext != "bed") {
    stop("GenoFile must have a .bed extension (PLINK bed file).")
  }

  bedFile <- GenoFile
  if (!is.null(GenoFileIndex)) {
    if (length(GenoFileIndex) == 2) {
      bimFile <- GenoFileIndex[1]
      famFile <- GenoFileIndex[2]
    } else {
      stop("GenoFileIndex must be a character vector of length 2 (PLINK bim and fam files).")
    }
  } else {
    bimFile <- sub("\\.bed$", ".bim", GenoFile)
    famFile <- sub("\\.bed$", ".fam", GenoFile)
  }
  for (f in c(bedFile, bimFile, famFile)) {
    if (!file.exists(f)) stop("PLINK file not found: ", f)
  }

  # validate output file
  if (!is.character(OutputFile) || length(OutputFile) != 1) {
    stop("'OutputFile' must be a single string.")
  }
  if (!is.logical(overwrite) || length(overwrite) != 1) {
    stop("'overwrite' must be TRUE or FALSE.")
  }

  if (file.exists(OutputFile)) {
    if (overwrite) {
      if (!file.remove(OutputFile)) stop("Cannot remove existing OutputFile: ", OutputFile)
      .message("Existing output file removed (overwrite = TRUE).")
    } else {
      stop("OutputFile exists. Set overwrite=TRUE or use another path.")
    }
  }

  # validate control parameters
  if (!is.null(control) && !is.list(control)) stop("'control' must be a list.")
  control <- modifyList(common_defaults, control)
  control <- .check_common_control(control)
  control <- modifyList(method_defaults, control)
  control <- .check_marker_control(control, objNull)

  # print parameters
  .message("Parameters for Marker-Level Tests:")

  params <- list(
    Method = null_class, OutputFile = OutputFile, overwrite = overwrite,
    bedFile = bedFile, bimFile = bimFile, famFile = famFile
  )
  for (nm in names(params)) {
    val <- params[[nm]]
    msg_val <- if (is.null(val)) "NULL" else as.character(val)
    message(sprintf("    %s: %s", nm, msg_val))
  }

  message("    Control parameters:")
  for (nm in names(control)) {
    val <- control[[nm]]
    msg_val <- if (is.null(val)) "NULL" else as.character(val)
    message(sprintf("        %s: %s", nm, msg_val))
  }

  if (!control$AllMarkers) {
    .message("Marker filters specified. The final set is: (union of includes) minus (union of excludes).")
  }

  # Prevent BLAS oversubscription: each of the nthreads worker threads may call
  # Armadillo, which in turn calls BLAS. Without pinning BLAS to 1 thread, a
  # machine with 8 cores running 4 workers could spin up 32 BLAS threads.
  if (control$nthreads > 1L) .enforce_single_blas_thread(control$nthreads)

  # dispatch to unified C++ bridge
  mtMarkerBridgeInCPP(objNull, OutputFile, control, bedFile, bimFile, famFile)

  .message("Analysis complete. Results saved to '%s'.", OutputFile)
  invisible(NULL)
}


.check_common_control <- function (control) {

  if (!is.numeric(control$SPA_Cutoff) || control$SPA_Cutoff <= 0) {
    stop("control$SPA_Cutoff must be numeric > 0.")
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
  if (!is.numeric(control$nMarkersEachChunk) || control$nMarkersEachChunk < 256) {
    stop("control$nMarkersEachChunk must be numeric >= 256.")
  }
  if (!is.null(control$AlleleOrder) &&
      !control$AlleleOrder %in% c("ref-first", "alt-first")) {
    stop("control$AlleleOrder must be 'ref-first' or 'alt-first'.")
  }
  if (!control$hwe %in% c("chi2", "exact")) {
    stop("control$hwe must be 'chi2' or 'exact'.")
  }

  marker_files <- c(
    control$IDsToIncludeFile, control$RangesToIncludeFile,
    control$IDsToExcludeFile, control$RangesToExcludeFile
  )
  marker_files <- marker_files[!is.null(marker_files)]
  if (length(marker_files) > 0) {
    for (f in marker_files) {
      if (!file.exists(f)) stop("File not found: ", f)
    }
    control$AllMarkers <- FALSE
  }

  if (!is.numeric(control$nthreads) || control$nthreads < 1 || (control$nthreads %% 1) != 0) {
      stop("'nthreads' must be a positive integer.")
  }

  control
}

.check_marker_control <- function (control, obj_null) {
  null_class <- class(obj_null)

  if (null_class == "POLMM_NULL_Model") {
    return(control)
  } else if (null_class %in% c("WtCoxG_NULL_Model", "LEAF_NULL_Model")) {
    if (!is.numeric(control$cutoff) || control$cutoff <= 0 || control$cutoff >= 1) {
      stop("control$cutoff must be numeric in (0, 1).")
    }
  } else if (null_class == "SPACox_NULL_Model") {
    if (!is.numeric(control$pVal_covaAdj_Cutoff) || control$pVal_covaAdj_Cutoff <= 0) {
      stop("control$pVal_covaAdj_Cutoff must be numeric > 0.")
    }
  } else if (null_class == "SPAmix_NULL_Model") {
    if (!control$dosage_option %in% c("rounding_first", "rounding_last")) {
      stop("control$dosage_option must be 'rounding_first' or 'rounding_last'.")
    }
  } else if (null_class == "SPAmixPlus_NULL_Model") {
    if (!control$dosage_option %in% c("rounding_first", "rounding_last")) {
      stop("control$dosage_option must be 'rounding_first' or 'rounding_last'.")
    }
    if (!control$afFilePrecision %in% c("double", "single", "text")) {
      stop("control$afFilePrecision must be 'double', 'single', or 'text'.")
    }
    if (is.null(control$afFilePath)) {
      stop("control$afFilePath must be provided for SPAmixPlus.")
    }
  } else if (null_class %in% c("SPAGRM_NULL_Model", "SAGELD_NULL_Model", "SPAsqr_NULL_Model")) {
    if (length(obj_null$MAF_interval) > 1 && control$min_maf_marker < min(obj_null$MAF_interval)) {
      stop("min_maf_marker is out of MAF_interval. Please reset.")
    }
  } else {
    stop("Unsupported null model class: ", null_class)
  }

  control
}

.enforce_single_blas_thread <- function(nthreads) {
  # Preferred: RhpcBLASctl calls the library API at runtime (OpenBLAS, MKL, etc.)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    cur <- RhpcBLASctl::blas_get_num_procs()
    if (cur != 1L) {
      RhpcBLASctl::blas_set_num_threads(1L)
      .message("BLAS threads set to 1 (was %d) to avoid oversubscription with %d workers.", cur, nthreads)
    }
    return(invisible(NULL))
  } else {
     .message("RhpcBLASctl not available. Install RhpcBLASctl or set BLAS threads via environment variables.")
  }
}

.message <- function (msg, ...) {
  message(
    sprintf("[INFO] %s %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(msg, ...))
  )
}
