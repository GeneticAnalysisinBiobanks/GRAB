## ------------------------------------------------------------------------------
## GRAB_Null_Model.R
##
## Functions:
##   GRAB.NullModel: Top-level API to fit a null model object used by
##                   GRAB.Marker and GRAB.Region.
## ------------------------------------------------------------------------------

#' @title Top-level API for generating a null model object used by GRAB.Marker and GRAB.Region
#'
#' @description GRAB performs two-step genetic association testing. This function 
#' implements the first step: fitting a null model and preparing the dataset required for 
#' downstream marker-level (\code{GRAB.Marker}) and region-level (\code{GRAB.Region}) analyses.
#'
#' @param formula (formula) Formula with response variable(s) on the left and covariates
#'   on the right. Do not include an intercept (added automatically). 
#'   For SPAmix with traitType "Residual", multiple response variables (separated by "+") are supported.
#' @param data (data.frame) Data frame containing response variables and covariates in the formula.
#'   Parameter "subset" is deprecated. All subjects with phenotype data will be used. Missing values
#'   should be coded as \code{NA}. Other values (e.g., -9, -999) are treated as numeric.
#' @param subjIDcol (character or NULL) Column name in \code{data} containing subject IDs.
#' @param subjData (character vector or NULL) Subject IDs aligned with rows of \code{data}.
#'   Exactly one of \code{subjIDcol} or \code{subjData} must be provided.
#' @param method (character) Supported methods:
#'   \itemize{
#'     \item "POLMM": Ordinal traits. See \code{?\link{GRAB.POLMM}}.
#'     \item "SPACox": Time-to-event or Residual. See \code{?\link{GRAB.SPACox}}.
#'     \item "SPAmix": Time-to-event or Residual. See \code{?\link{GRAB.SPAmix}}.
#'     \item "WtCoxG": Time-to-event traits. See \code{?\link{GRAB.WtCoxG}}.
#'   }
#' @param traitType (character) Supported: "ordinal", "time-to-event", and "Residual".
#' @param GenoFile (character or NULL) Path to genotype file. Supported formats determined by extension:
#'   \itemize{
#'     \item PLINK: "prefix.bed"
#'     \item BGEN: "prefix.bgen" (version 1.2 with 8-bit compression)
#'   }
#' @param GenoFileIndex (character vector or NULL) Associated files for the genotype file (auto-detected if NULL):
#'   \itemize{
#'     \item PLINK: c("prefix.bim", "prefix.fam")
#'     \item BGEN: c("prefix.bgen.bgi", "prefix.sample") or c("prefix.bgen.bgi")
#'   }
#' @param SparseGRMFile (character or NULL) Path to a sparse GRM file. 
#'   The file must be whitespace-delimited with three columns in the order:
#'   \itemize{
#'     \item Column 1: Subject ID 1
#'     \item Column 2: Subject ID 2
#'     \item Column 3: Genetic correlation between the two subjects
#'   }
#'   See \code{system.file("extdata", "SparseGRM.txt", package = "GRAB")} for an example.
#'   See \code{?\link{getSparseGRM}} for details on generating a sparse GRM.
#' @param control (list or NULL) List of additional, less commonly used parameters.
#'   See the corresponding method documentation for available options and defaults.
#' @param ... Additional method-specific parameters.
#'
#' @return
#' An S3 object with class "\{method\}_NULL_Model". All returned objects contain the following elements:
#' \describe{
#'   \item{N}{Sample size (integer).}
#'   \item{subjData}{Character vector of subject IDs included in analysis.}
#'   \item{Call}{Original function call.}
#'   \item{sessionInfo}{R session and package information.}
#'   \item{time}{Analysis completion timestamp (character).}
#'   \item{control}{List of control parameters used in fitting.}
#' }
#'
#' This object serves as input for \code{\link{GRAB.Marker}} or \code{\link{GRAB.Region}}.
#' 
#' See method-specific documentation for additional elements included in the returned object.
#'
GRAB.NullModel <- function(
  formula,
  data,
  subjIDcol = NULL,
  subjData = NULL,
  method,
  traitType,
  GenoFile = NULL,
  GenoFileIndex = NULL,
  SparseGRMFile = NULL,
  control = NULL,
  PairwiseIBDFile = NULL,
  ...
) {

  supported_traitTypes <- c("ordinal", "time-to-event", "Residual", "quantitative")
  supported_methods <- c("POLMM", "SPACox", "SPAmix", "WtCoxG", "SPAsqr")

  # ========== Validate and configure parameters ==========

  # Validate formula (required, no default)
  if (!inherits(formula, "formula")) {
    stop("Argument 'formula' should be a formula object.")
  }
  responseVars <- all.vars(formula[[2]])                     # character vector: left side
  covariateVars <- all.vars(formula[[3]])                    # character vector: right side

  # Validate data (required, no default)
  if (!is.data.frame(data)) {
    stop("Argument 'data' should be a data.frame.")
  }
  if (nrow(data) < 1) {
    stop("Argument 'data' should have at least one row.")
  }
  header <- colnames(data)                                   # character vector

  # Validate subjIDcol and subjData: exactly one must be provided
  if (is.null(subjIDcol) && is.null(subjData)) {
    stop("Exactly one of 'subjIDcol' or 'subjData' must be provided.")
  }
  if (!is.null(subjIDcol) && !is.null(subjData)) {
    stop("Exactly one of 'subjIDcol' or 'subjData' must be provided, not both.")
  }
  
  if (!is.null(subjIDcol)) {
    if (!is.character(subjIDcol) || length(subjIDcol) != 1) {
      stop("Argument 'subjIDcol' should be a character string (column name).")
    }
  }
  
  if (!is.null(subjData)) {
    if (length(subjData) != nrow(data)) {
      stop("Length of 'subjData' (", length(subjData), ") must match number of rows in 'data' (", nrow(data), ").")
    }
    if (!is.character(subjData) && !is.numeric(subjData)) {
      stop("Argument 'subjData' should be a character or numeric vector of subject IDs.")
    }
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

  # Validate that columns exist in data
  missingVars <- setdiff(covariateVars, header)              # character vector
  
  if ((traitType %in% c("ordinal")) && is.symbol(formula[[2]]) && !(responseVars %in% header)) {
    # The left side of the formula is a single column name, but is absent in header.
    missingVars <- c(responseVars, missingVars)
  }

  if (!is.null(subjIDcol) && !(subjIDcol %in% header)) {
    missingVars <- c(subjIDcol, missingVars)
  }

  if (length(missingVars) > 0) {
    stop("Columns not found in data: ", paste(missingVars, collapse = ", "))
  }

  # Extract and further validate subjData
  if (!is.null(subjIDcol)) {
    subjData <- data[[subjIDcol]]
    if (!is.character(subjData) && !is.numeric(subjData)) {
      stop("Column '", subjIDcol, "' should contain character or numeric subject IDs.")
    }
  }

  subjData <- as.character(subjData)
  if (any(duplicated(subjData))) {
    stop("Subject IDs contain duplicates, which are not supported.")
  }

  # Validate SparseGRMFile format (if provided)
  if (!is.null(SparseGRMFile)) {
    if (!file.exists(SparseGRMFile)) {
      stop("SparseGRMFile does not exist: ", SparseGRMFile)
    }
    
    firstLines <- readLines(SparseGRMFile, n = 2)
    if (length(firstLines) < 2) {
      stop("SparseGRMFile has fewer than 2 lines.")
    }

    firstRow <- strsplit(firstLines[1], "\\s+")[[1]]
    secondRow <- strsplit(firstLines[2], "\\s+")[[1]]

    if (length(firstRow) != 3 || length(secondRow) != 3) {
      stop("SparseGRMFile should have exactly 3 columns.")
    }
    
    if (is.na(suppressWarnings(as.numeric(secondRow[3])))) {
      stop("SparseGRMFile third column (genetic correlation) must be numeric.")
    }
  }

  # Validate control parameter (optional, default NULL)
  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be a list of control parameters.")
  }

  checkResult <- switch(method,
    POLMM = checkControl.NullModel.POLMM(traitType, GenoFile, SparseGRMFile, control),
    SPACox = checkControl.NullModel.SPACox(traitType, GenoFile, SparseGRMFile, control),
    SPAmix = checkControl.NullModel.SPAmix(traitType, GenoFile, SparseGRMFile, control),
    WtCoxG = checkControl.NullModel.WtCoxG(traitType, GenoFile, SparseGRMFile, control, ...),
    SPAsqr = checkControl.NullModel.SPAsqr(traitType, GenoFile, SparseGRMFile, control, PairwiseIBDFile)
  )
  
  control <- checkResult$control
  optionGRM <- checkResult$optionGRM

  # ========== Print all parameters ==========

  params <- list(
    Method = method,
    `Trait type` = traitType,
    `Formula` = deparse(formula),
    `Subject ID source` = ifelse(!is.null(subjIDcol), paste0("Column '", subjIDcol, "'"), "Provided vector"),
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
    } else if (traitType == "quantitative") {
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
                                 responseVars[1], responseVars[2], ...),
    SPAsqr = fitNullModel.SPAsqr(response, designMat, subjData, control, 
                                 GenoFile, SparseGRMFile, PairwiseIBDFile)
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
