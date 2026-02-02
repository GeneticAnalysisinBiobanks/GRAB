## ------------------------------------------------------------------------------
## GRAB_Marker2.R
##
## Simplified marker-level analysis with async I/O and multithreading.
## All processing happens in C++ - no chunking, no checkpointing.
## ------------------------------------------------------------------------------

#' Marker-level association analysis with async I/O and multithreading
#'
#' @param objNull Null model object from GRAB.NullModel()
#' @param GenoFile Path to genotype file (PLINK .bed or BGEN)
#' @param OutputFile Path to output file
#' @param GenoFileIndex Path to .bim (PLINK) or .bgi (BGEN) file
#' @param control List of control parameters (nWorkers, inputBufferSize, outputBufferSize, method-specific params)
#'
#' @return NULL (results written to OutputFile)
#' @export
GRAB.Marker2 <- function(
  objNull,
  GenoFile,
  OutputFile,
  GenoFileIndex = NULL,
  control = list()
) {

  method <- class(objNull)[1]
  method <- gsub("_NULL_Model", "", method)
  bfile <- tools::file_path_sans_ext(GenoFile)
  bedfile <- paste0(bfile, ".bed")
  bimfile <- paste0(bfile, ".bim")
  famfile <- paste0(bfile, ".fam")

  # Set default control parameters
  nCores <- parallel::detectCores()
  alleleOrderDefault <- if (genoType == "BGEN") "alt-first" else "ref-first"
  default.control <- list(
    nWorkers = min(10, max(1, nCores - 2)),
    inputBufferSize = 50,
    outputBufferSize = 200
  )
  control <- updateControl(control, default.control)


  mainMarkerInCPP2(
    outputFile = OutputFile,
    inputBufferSize = control$inputBufferSize,
    outputBufferSize = control$outputBufferSize,
    nWorkers = control$nWorkers,

    t_bimFile = bimfile,
    t_famFile = famfile,
    t_bedFile = GenoFile,
    t_SampleInModel = objNull$subjData,
    t_AlleleOrder = alleleOrderDefault,

    t_impute_method = "drop",
    t_missing_cutoff = 0.05,
    t_min_maf_marker = 1e-4,
    t_min_mac_marker = 20,

    t_resid = objNull$resid,
    t_PCs = objNull$PCs,
    t_N = objNull$N,
    t_SPA_Cutoff = 2,
    t_outlierList = objNull$outLierList
  )

  invisible(NULL)
}
