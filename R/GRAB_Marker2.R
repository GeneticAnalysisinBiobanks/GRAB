## ------------------------------------------------------------------------------
## GRAB_Marker2.R
##
## Simplified marker-level analysis with async I/O and multithreading.
## All processing happens in C++ - no chunking, no checkpointing.
## ------------------------------------------------------------------------------

# Marker-level association analysis with async I/O and multithreading
GRAB_Marker2 <- function(
  objNull,
  GenoFile,
  OutputFile,
  control = list()
) {

  method <- gsub("_NULL_Model", "", class(objNull)[1])
  control$method <- method
  control$outputFile <- OutputFile

  # Set multithreading defaults
  nCores <- parallel::detectCores()
  multithread.default <- list(
    nWorkers = min(10, max(1, nCores - 2)),
    inputBufferSize = 50,
    outputBufferSize = 200
  )
  control <- updateControl(control, multithread.default)

  # Set analysis method defaults
  control <- switch(
    method,
    POLMM = checkControl.Marker.POLMM(control),
    SPACox = checkControl.Marker.SPACox(control),
    SPAmix = checkControl.Marker.SPAmix(control),
    SPAGRM = checkControl.Marker.SPAGRM(control, objNull$MAF_interval),
    SAGELD = checkControl.Marker.SAGELD(control, objNull$MAF_interval),
    WtCoxG = checkControl.Marker.WtCoxG(control),
    SPAsqr = checkControl.Marker.SPAsqr(control),
    LEAF = checkControl.Marker.LEAF(control)
  )

  # Set genotype file parameters
  geno_prefix <- tools::file_path_sans_ext(GenoFile)
  geno_suffix <- tools::file_ext(GenoFile)

  genotype.default <- list(
      AlleleOrder = "",
      bedfile = "",
      bimfile = "",
      famfile = "",
      bgenFileName = "",
      bgenFileIndex = "",
      t_SampleInBgen = c(),
      t_isSparseDosageInBgen = FALSE,
      t_isDropmissingdosagesInBgen = FALSE,
      impute_method = "mean",
      missing_cutoff = 0.15,
      min_maf_marker = 0.001,
      min_mac_marker = 20,
      SPA_Cutoff = 2
  )
  control <- updateControl(control, genotype.default)  

  if (geno_suffix == "bed") {  
    control$genoType <- "PLINK"
    control$AlleleOrder = "alt-first"
    control$bedfile = paste0(geno_prefix, ".bed")
    control$bimfile = paste0(geno_prefix, ".bim")
    control$famfile = paste0(geno_prefix, ".fam")
  } else if (geno_suffix == "bgen") {
    control$genoType <- "BGEN"
    control$AlleleOrder = "ref-first"
    control$bgenFileName = paste0(geno_prefix, ".bgen")
    control$bgenFileIndex = paste0(geno_prefix, ".bgen.bgi")
  }

  mainMarkerInCPP2(control = control, objNull = objNull)

  return(invisible(NULL))
}
