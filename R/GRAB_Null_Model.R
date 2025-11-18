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
#'   Parameter "subset" is deprecated. All subjects with phenotype data will be used.
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

  # Validate data (required, no default)
  if (!is.data.frame(data)) {
    stop("Argument 'data' should be a data.frame.")
  }
  if (nrow(data) < 1) {
    stop("Argument 'data' should have at least one row.")
  }

  # Validate subjIDcol (required, no default)
  if (!is.character(subjIDcol) || length(subjIDcol) != 1) {
    stop("Argument 'subjIDcol' should be a character string (column name).")
  }

  # Validate that columns in formula exist in data
  responseVars <- all.vars(formula[[2]])                     # character vector: left side
  covariateVars <- all.vars(formula[[3]])                    # character vector: right side

  # For Residual trait type, response variables come from environment, not data
  if (traitType == "Residual") {
    neededVars <- c(covariateVars, subjIDcol)
  } else if (is.symbol(formula[[2]])) {
    neededVars <- c(responseVars, covariateVars, subjIDcol)
  } else {
    neededVars <- c(covariateVars, subjIDcol)
  }

  dataCols <- colnames(data)                                 # character vector
  missingVars <- setdiff(neededVars, dataCols)               # character vector

  if (length(missingVars) > 0) {
    stop("Variables not found in data: ", paste(missingVars, collapse = ", "))
  }

  if (!subjIDcol %in% dataCols) {
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
  designMat <- as.matrix(as.data.frame(data)[, covariateVars, drop = FALSE])

  # Extract response
  LeftInFormula <- deparse(formula[[2]])                       # character string
  LeftIncludesAdd <- grepl("\\+", LeftInFormula)               # logical
  if (LeftIncludesAdd) {
    
    if (method %in% c("SPAmix") && traitType == "Residual") {
      .message("SPAmix analysis will use residuals from %d models.", length(responseVars))

      # Evaluate all variable names on the left side of the formula
      response <- sapply(responseVars, function(varName) {      # matrix
        eval(as.name(varName), envir = parent.frame(2)) 
      })
    } else {
      stop("Only a single response variable is supported for method '", method,
           "' and trait type '", traitType, "'.")
    }
    
    class(response) <- "Residual"
    naSubjects <- apply(response, 1, function(x) all(is.na(x)))   # logical vector
  } else {

    mf <- stats::model.frame(formula, data, na.action = na.pass)
    response <- model.response(mf)                                # vector or matrix

    if (traitType == "time-to-event") {
      if (inherits(response, "Surv")) {
        naSubjects <- is.na(response[, 1]) | is.na(response[, 2])
      } else {
        stop("For time-to-event traits, the response variable must be a Surv object.")
      } 
    } else if (traitType == "ordinal") {
      if (is.factor(response)) {
        response <- droplevels(response)
        naSubjects <- is.na(response)
      } else {
        stop("For POLMM method, the response variable must be a factor (ordinal trait).")
      }
    } else if (traitType == "Residual") {
      class(response) <- "Residual"
      naSubjects <- is.na(response)
    } else {
      stop("Internal error: '", traitType, "' for method '", method, "'.")
    }
  }

  # ========== Remove subjects with missing phenotype data ==========

  nRemoved <- sum(naSubjects)
  if (nRemoved > 0) {
    .message("Removing %d subjects with missing phenotype data", nRemoved)

    if (LeftIncludesAdd) {
      response <- response[!naSubjects, , drop = FALSE]
    } else {
      response <- response[!naSubjects]
    }

    designMat <- designMat[!naSubjects, , drop = FALSE]
    subjData <- subjData[!naSubjects]
  } else {
    .message("All subjects have phenotype data.")
  }
      
  .message("Number of subjects included for subsequent analysis: %d", length(subjData))

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
