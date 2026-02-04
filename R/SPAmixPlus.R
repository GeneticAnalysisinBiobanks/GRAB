#' SPAmixPlus: Scalable and accurate genetic association testing with mixed models
#'
#' @description
#' SPAmixPlus extends SPAmix with improved allele frequency estimation using 
#' principal components. It follows a 3-step workflow:
#' 
#' Step 0: Pre-compute AF models (optional, for efficiency)
#' Step 1: Fit null model using SPAmix
#' Step 2: Perform marker tests using SPAmixPlus with optional pre-computed AF models
#'
#' @examples
#' # Step 0 (Optional): Pre-compute AF models for efficiency
#' SPAmixPlus.AF(GenoFile = GenoFile,
#'               PCs = PCmatrix,
#'               subjData = IID,
#'               outputFile = "af_model.db")
#'
#' # Step 1: Fit null model using SPAmix
#' objNull <- GRAB.NullModel(resid ~ AGE + GENDER + PC1 + PC2,
#'                           data = PhenoData,
#'                           subjData = IID,
#'                           method = "SPAmix",
#'                           traitType = "Residual",
#'                           control = list(PC_columns = "PC1,PC2"))
#'
#' # Step 2: Perform marker tests with SPAmixPlus
#' GRAB.Marker(objNull,
#'             GenoFile = GenoFile,
#'             OutputFile = OutputFile,
#'             control = list(afFile = "af_model.db"))  # Use pre-computed AF models
#'
#' @export
GRAB.SPAmixPlus <- function() {
  .message("?GRAB.SPAmixPlus for instructions")
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
    t_genoIndex = genoIndex,
    t_extraParams = list(afFile = afFile)
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
#' @keywords internal
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
  if (!is.matrix(PCs)) {
    stop("PCs must be a matrix.")
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
    t_genoType = objGeno$genoType,
    t_genoIndex = objGeno$markerInfo$genoIndex,
    t_afFileOutput = afFileOutput,
    t_pcs = PCs,
    t_afFilePrecision = control$afFilePrecision
  )
  
  cat("AF models saved to: ", afFileOutput, " (precision: ", control$afFilePrecision, ")\n", sep = "")
  
  invisible(NULL)
}

