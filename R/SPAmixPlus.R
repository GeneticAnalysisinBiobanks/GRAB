#' Instruction of SPAmixPlus method
#'
#' SPAmixPlus extends SPAmix with improved individual-specific allele frequency estimation using 
#' principal components, accounting for sample relatedness through a sparse genetic relationship 
#' matrix (GRM). It performs retrospective single-variant association tests using 
#' genotypes and residuals from null models of any complex trait in large-scale biobanks. 
#' SPAmixPlus supports complex population structures and sample relatedness.
#' 
#' @return NULL
#'
#' @examples
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultSPAmixPlus.txt")
#' afFileOutput <- file.path(tempdir(), "afModels.bin")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' 
#' # Step 0: Pre-calculate individual-specific allele frequencies
#' SPAmixPlus.AF(
#'   GenoFile = GenoFile,
#'   PCs = as.matrix(PhenoData[, c("PC1", "PC2", "PC3")]),
#'   subjData = PhenoData$IID,
#'   afFileOutput = afFileOutput,
#'   control = list(afFilePrecision = "double")
#' )
#'
#' # Step 1: Fit a null model and obtain residuals
#' residuals <- survival::coxph(
#'   survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
#'   data = PhenoData
#' )$residuals
#'
#' obj.SPAmixPlus <- GRAB.NullModel(
#'   residuals ~ AGE + GENDER + PC1 + PC2,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   SparseGRMFile = SparseGRMFile,
#'   method = "SPAmixPlus",
#'   traitType = "Residual",
#'   control = list(PC_columns = "PC1,PC2")
#' )
#'
#' # Step 2: Run marker-level association tests using SPAmixPlus
#' OutputFile <- file.path(tempdir(), "resultSPAmixPlus1.txt")
#' GRAB.Marker(
#'   obj.SPAmixPlus,
#'   GenoFile = GenoFile,
#'   OutputFile = OutputFile,
#'   control = list(
#'     afFilePath = afFileOutput,
#'     afFilePrecision = "double"
#'   )
#' )
#'
#' head(data.table::fread(OutputFile))
#'
#' @details
#'
#' \strong{Additional Control Parameters for GRAB.NullModel()}:
#' \itemize{
#'   \item \code{PC_columns} (character, required): Comma-separated column names
#'      of principal components (e.g., \code{"PC1,PC2"}).
#'   \item \code{OutlierRatio} (numeric, default: 1.5): IQR multiplier for outlier detection.
#'      Outliers are defined as values outside \eqn{[Q1 - r \times IQR, Q3 + r \times IQR]}, where \eqn{r} is the multiplier.
#' }
#'
#' \strong{Method-specific elements in the \code{SPAmixPlus_NULL_Model} object returned by \code{GRAB.NullModel()}:}:
#' \itemize{
#'   \item \code{resid}: Residuals from mixed model (matrix or "Residual" class).
#'   \item \code{yVec}: Phenotype vector (numeric or "Residual" class).
#'   \item \code{PCs}: Principal components for dimension reduction (matrix).
#'   \item \code{nPheno}: Number of phenotypes analyzed (integer).
#'   \item \code{outLierList}: List identifying outlier subjects for SPA adjustment.
#'   \item \code{sparseGRM}: Sparse genetic relationship matrix for relatedness adjustment.
#' }
#' 
#' \strong{Additional Control Parameters for GRAB.Marker()}:
#' \itemize{
#'   \item \code{afFilePath} (character, required): Path to pre-computed AF model file from 
#'      \code{SPAmixPlus.AF()} (Step 0).
#'   \item \code{afFilePrecision} (character, default: "double"): Precision level for AF model file.
#'      Must be either "double" (float64), "single" (float32), or "text" (TSV format).
#'   \item \code{dosage_option} (character, default: "rounding_first"):
#'      Dosage handling option. Must be either "rounding_first" or "rounding_last".
#' }
#'
#' \strong{Output file columns}:
#' \describe{
#'   \item{Pheno}{Phenotype identifier (for multi-trait analysis).}
#'   \item{Marker}{Marker identifier (rsID or CHR:POS:REF:ALT).}
#'   \item{Info}{Marker information in format CHR:POS:REF:ALT.}
#'   \item{AltFreq}{Alternative allele frequency in the sample.}
#'   \item{AltCounts}{Total count of alternative alleles.}
#'   \item{MissingRate}{Proportion of missing genotypes.}
#'   \item{Pvalue}{P-value from the score test.}
#'   \item{zScore}{Z-score from the score test.}
#' }
#'
#' @references
#' Ma et al. (2025). SPAmix: a scalable, accurate, and universal analysis framework for large‑scale 
#' genetic association studies in admixed populations. \doi{10.1186/s13059-025-03827-9}
#'
#' @export
GRAB.SPAmixPlus <- function() {
  .message("?GRAB.SPAmixPlus for instructions")
}



#' Check control parameters for SPAmixPlus null model
#'
#' @keywords internal
checkControl.NullModel.SPAmixPlus <- function(traitType, GenoFile, SparseGRMFile, control) {
  
  # Use SPAmix control checking as base
  control <- checkControl.NullModel.SPAmix(traitType, GenoFile, SparseGRMFile, control)

  if (is.null(SparseGRMFile) || !file.exists(SparseGRMFile)) {
    stop("SparseGRMFile must be provided and exist for SPAmixPlus method.")
  }
  
  return(control)
}


#' Fit null model for SPAmixPlus
#'
#' @keywords internal
fitNullModel.SPAmixPlus <- function(
  response,
  designMat,
  subjData,
  control,
  SparseGRMFile
) {
  
  if (!inherits(response, "Residual")) {
    stop("For SPAmixPlus, the response variable should be 'Residual'.")
  }

  # Call SPAmix null model fitting first
  objNull <- fitNullModel.SPAmix(response, designMat, subjData, control)
  
  # Load sparse GRM and create DataFrame with proper column names
  sparseGRM_df <- data.table::fread(SparseGRMFile, header = TRUE)
  
  # Create index mapping from subjData to 0-based indices
  id_map <- setNames(seq_along(subjData) - 1, subjData)
  
  # Map ID columns to indices
  id1_index <- id_map[as.character(sparseGRM_df[[1]])]
  id2_index <- id_map[as.character(sparseGRM_df[[2]])]
  value <- as.numeric(sparseGRM_df[[3]])
  
  # Remove NA indices
  valid_idx <- !is.na(id1_index) & !is.na(id2_index)
  
  # Create DataFrame with correct column names for C++
  objNull$sparseGRM <- data.frame(
    id1_index = as.integer(id1_index[valid_idx]),
    id2_index = as.integer(id2_index[valid_idx]),
    value = value[valid_idx]
  )
  
  # Update class to SPAmixPlus_NULL_Model
  class(objNull) <- "SPAmixPlus_NULL_Model"
  
  return(objNull)
}


#' Check control parameters for SPAmixPlus marker testing
#'
#' @keywords internal
checkControl.Marker.SPAmixPlus <- function(control) {
  # Set default control parameters
  default.control <- list(
    dosage_option = "rounding_first",  # Options: "rounding_first" or "rounding_last"
    afFilePrecision = "double"  # Options: "double", "single", or "text"
  )
  
  control <- updateControl(control, default.control)
  
  # Validate dosage option
  if (!control$dosage_option %in% c("rounding_first", "rounding_last")) {
    stop("control$dosage_option should be 'rounding_first' or 'rounding_last'.")
  }

  if (!control$afFilePrecision %in% c("double", "single", "text")) {
    stop("control$afFilePrecision should be 'double', 'single', or 'text'.")
  }

  if (is.null(control$afFilePath)) {
    stop("control$afFilePath must be provided for SPAmixPlus marker testing.")
  }
  
  return(control)
}


setMarker.SPAmixPlus <- function(objNull, control) {

  setSPAmixPlusobjInCPP(
    t_resid = objNull$resid,
    t_PCs = objNull$PCs,
    t_N = objNull$N,
    t_SPA_Cutoff = control$SPA_Cutoff,
    t_outlierList = objNull$outLierList,
    t_sparseGRM = objNull$sparseGRM,
    t_afFilePath = control$afFilePath,
    t_afFilePrecision = control$afFilePrecision
  )

}



#' Run SPAmixPlus marker-level tests (used by GRAB.Marker framework)
#'
#' @param genoType Character, genotype file type ("PLINK" or "BGEN").
#' @param genoIndex Integer vector of marker indices to test.
#' @param objNull SPAmixPlus_NULL_Model object from fitNullModel.SPAmixPlus.
#' @param control List of control parameters (optional).
#'
#' @return Data frame with marker-level association statistics.
#'
#' @keywords internal
mainMarker.SPAmixPlus <- function(genoType, genoIndex, objNull, control = NULL) {

  # Call C++ backend for marker testing
  OutList <- mainMarkerInCPP(
    t_method = "SPAmixPlus",
    t_genoType = genoType,
    t_genoIndex = genoIndex
  )
  
  nPheno <- objNull$nPheno
  
  # Format output data frame with marker-level statistics
  obj.mainMarker <- data.frame(
    Pheno = paste0("pheno_", 1:nPheno),
    Marker = rep(OutList$markerVec, each = nPheno),
    Info = rep(OutList$infoVec, each = nPheno),
    AltFreq = rep(OutList$altFreqVec, each = nPheno),
    AltCounts = rep(OutList$altCountsVec, each = nPheno),
    MissingRate = rep(OutList$missingRateVec, each = nPheno),
    Pvalue = OutList$pvalVec,
    zScore = OutList$zScore
  )
  
  return(obj.mainMarker)
}


#' Pre-compute allele frequency models for SPAmixPlus
#'
#' @description
#' Optional step to pre-compute allele frequency estimation models for all markers.
#' Uses only principal components and sample information (no null model required).
#' These models can be reused in Step 2 to speed up analysis.
#'
#' @param GenoFile Path to genotype file (PLINK or BGEN format)
#' @param PCs Matrix of principal components (N x nPCs)
#' @param subjData Vector of subject IDs aligned with PC rows
#' @param afFileOutput Path to output file for storing AF models
#' @param GenoFileIndex Path to genotype file index (optional)
#' @param control List of control parameters
#'
#' @details
#' The function fits allele frequency estimation models using principal components:
#' \itemize{
#'   \item Status 0: Mean-based (low MAC markers)
#'   \item Status 1: Linear regression AF ~ intercept + PC1 + PC2 + ...
#'   \item Status 2: Logistic regression (when linear model has boundary violations)
#' }
#' 
#' Precision levels (set via control$afFilePrecision):
#' \itemize{
#'   \item double (float64): ~15-16 decimal digits, best precision, largest files
#'   \item single (float32): ~7-8 decimal digits, good balance of precision and size (50\\% smaller than double)
#'   \item text: 6 significant digits, human-readable TSV format, good for inspection and debugging but slow to read in step 2
#' }
#'
#' @return Invisible NULL. AF models saved to outputFile.
#' @export
SPAmixPlus.AF <- function(
  GenoFile,
  PCs,
  subjData,
  afFileOutput,
  GenoFileIndex = NULL,
  control = list()
) {
  
  # Validate inputs
  if (!is.matrix(PCs) && !is.data.frame(PCs)) {
    stop("PCs must be a matrix or data.frame.")
  }

  if (nrow(PCs) != length(subjData)) {
    stop("Number of rows in PCs must match length of subjData.")
  }
  
  # Set default control parameters
  if (is.null(control)) control <- list()
  default.control <- list(
    impute_method = "mean",
    missing_cutoff = 0.15,
    min_maf_marker = 0.0001,
    min_mac_marker = 1,
    omp_num_threads = 1,
    afFilePrecision = "double"  # Options: "double", "single", "text"
  )
  
  control <- updateControl(control, default.control)
  
  # Validate afFilePrecision
  if (!control$afFilePrecision %in% c("double", "single", "text")) {
    stop("control$afFilePrecision must be 'double', 'single', or 'text'.")
  }
  
  # Ensure all markers are included (no filtering at this stage)
  control$AllMarkers <- TRUE
  
  # Initialize genotype input
  objGeno <- setGenoInput(
    GenoFile, 
    GenoFileIndex = GenoFileIndex, 
    SampleIDs = as.character(subjData), 
    control = control
  )
  
  cat("Pre-computing AF models for ", length(objGeno$markerInfo$genoIndex), 
      " markers...\n", sep = "")
  
  # Set global C++ variables
  setMarker_GlobalVarsInCPP(
    t_impute_method = control$impute_method,
    t_missing_cutoff = control$missing_cutoff,
    t_min_maf_marker = control$min_maf_marker,
    t_min_mac_marker = control$min_mac_marker,
    t_omp_num_threads = control$omp_num_threads
  )

  # Export AF models to file (PCs passed directly)
  exportAFModelInCPP(
    genoType = objGeno$genoType,
    genoIndex = objGeno$markerInfo$genoIndex,
    afFileOutput = afFileOutput,
    t_pcs = as.matrix(PCs),
    afFilePrecision = control$afFilePrecision
  )
  
  cat("AF models saved to: ", afFileOutput, " (precision: ", control$afFilePrecision, ")\n", sep = "")
  
  invisible(NULL)
}

