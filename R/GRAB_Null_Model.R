## ------------------------------------------------------------------------------
## GRAB_Null_Model.R
## Build method-specific null models used by downstream association tests. Handles:
##   * Parse formula/data and construct design matrices
##   * Handle intercepts, missingness, and multi-response residual inputs
##   * Configure GRM options (dense vs sparse) and genotype metadata
##   * Validate/merge method-specific control parameters
##   * Dispatch to method-specific null-model fitters
##
## Functions:
##   GRAB.NullModel: Top-level API to fit a null model object used by
##                   GRAB.Marker and GRAB.Region.
## ------------------------------------------------------------------------------

#' Fit a null model for genetic association testing
#'
#' Estimates model parameters and residuals from a null model containing phenotypes,
#' covariates, and optionally a genetic relationship matrix (GRM) for subsequent
#' association testing.
#'
#' @param formula (formula) Formula with response variable(s) on the left and covariates
#'   on the right. Do not include an intercept (added automatically). Missing values
#'   should be coded as \code{NA}. Other values (e.g., -9, -999) are treated as numeric.
#'   For SPAmix with traitType "Residual", multiple response variables are supported.
#' @param data (data.frame) Data frame containing response variables and covariates in the formula.
#' @param subjIDcol (character) Column name in \code{data} containing subject IDs.
#' @param method (character) Statistical method. Supported methods:
#'   \itemize{
#'     \item "POLMM": Ordinal traits. See \code{\link{GRAB.POLMM}}.
#'     \item "SPACox": Time-to-event or residual traits. See \code{\link{GRAB.SPACox}}.
#'     \item "SPAmix": Time-to-event or residual traits. See \code{\link{GRAB.SPAmix}}.
#'     \item "WtCoxG": Time-to-event traits. See \code{\link{GRAB.WtCoxG}}.
#'   }
#' @param traitType (character) Trait type: "ordinal", "time-to-event", or "Residual".
#' @param GenoFile (character or NULL) Path to genotype file (PLINK or BGEN format).
#'   Required for dense GRM construction. GRAB supports both dense and sparse GRM for
#'   relatedness adjustment. Dense GRM is constructed from GenoFile (PLINK/BGEN format).
#'   See \code{\link{GRAB.ReadGeno}} and \code{\link{getTempFilesFullGRM}} for details.
#' @param GenoFileIndex (character or NULL) Index files for the genotype file. If \code{NULL}
#'   (default), uses same prefix as \code{GenoFile}. See \code{\link{GRAB.ReadGeno}}
#'   for details.
#' @param SparseGRMFile (character or NULL) Path to sparse GRM file. Alternative to
#'   dense GRM construction. Pre-computed sparse GRM matrix.
#'   See \code{\link{getSparseGRM}} for details.
#' @param control (list or NULL) List of method-specific control parameters.
#'   See the corresponding method documentation for available options and defaults.
#' @param ... Additional arguments for method-specific functions.
#'
#' @return
#' S3 object with class "\{method\}_NULL_Model" containing:
#' \describe{
#'   \item{N}{Sample size in analysis.}
#'   \item{yVec}{Phenotype vector used in analysis.}
#'   \item{beta}{Estimated covariate effect coefficients.}
#'   \item{subjData}{Subject IDs included in analysis.}
#'   \item{sessionInfo}{R session and package information.}
#'   \item{Call}{Function call with all arguments.}
#'   \item{time}{Analysis completion timestamp.}
#'   \item{control}{Control parameters used in fitting.}
#'   \item{tau}{Variance component estimates (for mixed models).}
#'   \item{SparseGRM}{Sparse GRM (if used).}
#' }
#'
#' This object serves as input for \code{\link{GRAB.Marker}} and \code{\link{GRAB.Region}}.
#'
#' For method-specific examples, see:
#' \itemize{
#'   \item POLMM method: \code{\link{GRAB.POLMM}}
#'   \item SPACox method: \code{\link{GRAB.SPACox}}
#'   \item SPAmix method: \code{\link{GRAB.SPAmix}}
#'   \item WtCoxG method: \code{\link{GRAB.WtCoxG}}
#' }
#'
GRAB.NullModel <- function(
  formula,
  data,
  subjIDcol,
  method,
  traitType,
  GenoFile = NULL,
  GenoFileIndex = NULL,
  SparseGRMFile = NULL,
  control = NULL,
  ...
) {

  supported_traitTypes <- c("ordinal", "time-to-event", "Residual")
  supported_methods <- c("POLMM", "SPACox", "SPAmix", "WtCoxG")

  # ========== Validate and configure parameters ==========

  # Validate formula (required, no default) and extract variable names
  if (!inherits(formula, "formula")) {
    stop("Argument 'formula' should be a formula object.")
  }

  responseVars <- all.vars(formula[[2]])                     # character vector: left side
  covariateVars <- all.vars(formula[[3]])                    # character vector: right side

  if (length(responseVars) == 0) {
    stop("formula must include at least one response variable on the left side.")
  }
  if (length(covariateVars) == 0) {
    stop("Formula must include at least one covariate on the right side.")
  }
  
  # Validate data (required, no default)
  if (!is.data.frame(data)) {
    stop("Argument 'data' should be a data.frame.")
  }
  if (nrow(data) < 1) {
    stop("Argument 'data' should have at least one row.")
  }
  
  missingVars <- setdiff(covariateVars, colnames(data))      # character vector
  if (length(missingVars) > 0) {
    stop("Variables not found in data: ", paste(missingVars, collapse = ", "))
  }

  # Validate subjIDcol (required, no default)
  if (!is.character(subjIDcol) || length(subjIDcol) != 1) {
    stop("Argument 'subjIDcol' should be a character string (column name).")
  }
  if (!subjIDcol %in% colnames(data)) {
    stop("Column '", subjIDcol, "' not found in data.")
  }

  # Extract subjData and validate subject IDs
  subjData <- data[[subjIDcol]]                              # character vector
  if (!is.character(subjData) && !is.numeric(subjData)) {
    stop("Column '", subjIDcol, "' should contain character or numeric subject IDs.")
  }
  subjData <- as.character(subjData)                         # character vector

  if (any(duplicated(subjData))) {
    stop("Column '", subjIDcol, "' contains duplicated subject IDs, which are not supported.")
  }

  # Validate method (required, no default)
  if (!is.character(method) || length(method) != 1) {
    stop("Argument 'method' is required and should be a character string.")
  }

  if (!method %in% supported_methods) {
    stop(
      "Argument 'method' should be one of: ",
      paste(paste0('"', supported_methods, '"'), collapse = ", ")
    )
  }

  # Validate traitType (required, no default)
  if (!is.character(traitType) || length(traitType) != 1) {
    stop("Argument 'traitType' is required and should be a character string.")
  }

  if (!traitType %in% supported_traitTypes) {
    stop(
      "Argument 'traitType' should be one of: ",
      paste(paste0('"', supported_traitTypes, '"'), collapse = ", ")
    )
  }

  # Validate control parameter (optional, default NULL)
  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be a list of control parameters.")
  }

  # Method-specific validation
  checkResult <- switch(method,
    POLMM = checkControl.NullModel.POLMM(traitType, GenoFile, SparseGRMFile, control),
    SPACox = checkControl.NullModel.SPACox(traitType, GenoFile, SparseGRMFile, control),
    SPAmix = checkControl.NullModel.SPAmix(traitType, GenoFile, SparseGRMFile, control),
    WtCoxG = checkControl.NullModel.WtCoxG(traitType, GenoFile, SparseGRMFile, control, ...)
  )
  
  control <- checkResult$control
  optionGRM <- checkResult$optionGRM

  # ========== Print all parameters ==========

  params <- list(
    Method = method,
    `Trait type` = traitType,
    `Formula` = deparse(formula),
    `Subject ID column` = subjIDcol,
    `Sample size of input` = nrow(data),
    `Genotype file` = ifelse(is.null(GenoFile), "Not provided", GenoFile),
    `Genotype index file` = ifelse(is.null(GenoFileIndex), "Default", GenoFileIndex),
    `Sparse GRM file` = ifelse(is.null(SparseGRMFile), "Not provided", SparseGRMFile),
    `GRM option` = ifelse(!is.null(optionGRM), optionGRM, "Not applicable")
  )
  .printParameters("Parameters for Null Model Fitting", params, control)

  # ========== Extract and validate designMat, response, subjData ==========

  # Extract designMat
  designMat <- as.matrix(data[, covariateVars, drop = FALSE])

  # Extract response based on traitType
  if (traitType == "ordinal") {

    if (length(responseVars) != 1) {
      stop("For traitType 'ordinal', the left side of the formula must be a single column name.")
    }

    if (!responseVars[1] %in% colnames(data)) {
      stop("Response variable '", responseVars[1], "' not found in data.")
    }

    responseData <- data[[responseVars[1]]]                  # vector
    if (!is.factor(responseData)) {
      stop("For traitType 'ordinal', response variable must be factor. ",
           "The class of the current response variable is'", class(responseData), "'.")
    }
    response <- ordered(responseData)                        # ordered factor
    anyNA <- is.na(response)                                 # logical vector

  } else if (traitType == "time-to-event") {

    if (length(responseVars) != 2) {
      stop("For traitType 'time-to-event', the left side of the formula must have exactly two items. ",
           "The first is time, the second is event. They are extracted by all.vars(formula[[2]])".)
    }

    missingResponseVars <- setdiff(responseVars, colnames(data))
    if (length(missingResponseVars) > 0) {
      stop("Response variables not found in data: ", paste(missingResponseVars, collapse = ", "))
    }
    
    timeVar <- data[[responseVars[1]]]                               # numeric vector
    eventVar <- data[[responseVars[2]]]                              # numeric/logical vector
    response <- survival::Surv(time = timeVar, event = eventVar)     # Surv object
    anyNA <- is.na(response[, "time"]) | is.na(response[, "status"]) # logical vector

  } else if (traitType == "Residual") {

    if (length(responseVars) == 1) {

      responseData <- eval(as.name(responseVars), envir = parent.frame()) # numeric vector

      if (!is.numeric(responseData)) {
        stop("For traitType 'Residual', response variable should be numeric.")
      }

      response <- matrix(responseData, ncol = 1)             # single column matrix 
      class(response) <- "Residual"
      anyNA <- is.na(response[, 1])                          # logical vector
      
    } else {
      # Multiple residual variables (for SPAmix multi-trait analysis)
      if (!method %in% c("SPAmix")) {
        stop("Multiple response variables for traitType 'Residual' are only supported for method 'SPAmix'.")
      }
      .message("SPAmix analysis will use residuals from %d models.", length(responseVars))

      # Extract all residual variables as matrix from environment
      response <- sapply(responseVars, function(varName) {   # matrix
        eval(as.name(varName), envir = parent.frame(2))      # numeric vector
      })
      
      if (!is.numeric(response)) {
        stop("For traitType 'Residual', all response variables should be numeric.")
      }
                        
      class(response) <- "Residual"
      allNA <- apply(response, 1, function(x) all(is.na(x))) # logical vector
    }
  }

  # ========== Remove subjects with missing phenotype data ==========
  
  # Determine which subjects have missing phenotype data
  if (exists("allNA")) {
    # For SPAmix multi-trait: use allNA (subjects where all residuals are NA)
    naSubjects <- allNA
  } else if (exists("anyNA")) {
    # For other methods: use anyNA (subjects with any NA in response)
    naSubjects <- anyNA
  } else {
    stop("Internal error: could not determine missingness in phenotype data.")
  }
  
  # Remove subjects with missing data
  if (any(naSubjects)) {
    .message("Removing %d subjects with missing phenotype data", sum(naSubjects))
    
    if (is.matrix(response)) {
      response <- response[!naSubjects, , drop = FALSE]      # matrix
    } else if (inherits(response, "Surv")) {
      response <- response[!naSubjects, ]                    # Surv object
    } else {
      response <- response[!naSubjects]                      # ordered factor or vector
    }
    
    designMat <- designMat[!naSubjects, , drop = FALSE]      # matrix
    subjData <- subjData[!naSubjects]                        # character vector
  }

  nData <- length(subjData)                                   # integer
  .message("Number of subjects with phenotype: %d", nData)

  # ========== Fit null model ==========

  objNull <- switch(method,
    POLMM = fitNullModel.POLMM(response, designMat, subjData, control, optionGRM,
                               GenoFile, GenoFileIndex, SparseGRMFile),
    SPACox = fitNullModel.SPACox(response, designMat, subjData, control),
    SPAmix = fitNullModel.SPAmix(response, designMat, subjData, control),
    WtCoxG = fitNullModel.WtCoxG(response, designMat, subjData, control, data,
                                 GenoFile, GenoFileIndex, SparseGRMFile,
                                 responseVars[1], responseVars[2], ...)
  )

  # Add metadata to the null model object
  objNull$subjData <- subjData
  objNull$Call <- match.call()
  objNull$sessionInfo <- sessionInfo()
  objNull$time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  objNull$control <- control

  .message("Successfully finished fitting the null model")

  return(objNull)
}
