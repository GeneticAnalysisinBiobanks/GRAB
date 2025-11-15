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
#' @param subset (vector or NULL) Row specification for subsetting subjects. Default: NULL for all.
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
  subset = NULL,
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

  # formula and data will be validated during model frame extraction

  # Validate subset (optional, default NULL)
  if (!is.null(subset) && !is.vector(subset) && !is.logical(subset)) {
    stop("Argument 'subset' should be a vector or logical expression for subsetting subjects")
  }

  # Validate subjIDcol (required, no default)
  if (!is.character(subjIDcol) || length(subjIDcol) != 1) {
    stop("Argument 'subjIDcol' should be a character string (column name).")
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

  # Validate GenoFile and SparseGRMFile and assign optionGRM for POLMM method
  if (method == "POLMM") {
    # GenoFile validation (required)
    if (is.null(GenoFile)) {
      stop("Argument 'GenoFile' is required for method 'POLMM'.")
    }
    if (!is.character(GenoFile) || length(GenoFile) != 1) {
      stop("Argument 'GenoFile' should be a character string (file path).")
    }
    if (!file.exists(GenoFile)) {
        stop("Cannot find GenoFile: ", GenoFile)
    }

    # SparseGRMFile validation (optional)
    if (!is.null(SparseGRMFile)) {
      if (!is.character(SparseGRMFile) || length(SparseGRMFile) != 1) {
        stop("Argument 'SparseGRMFile' should be a character string (file path).")
      }
      if (!file.exists(SparseGRMFile)) {
        stop("Cannot find SparseGRMFile: ", SparseGRMFile)
      }
      optionGRM <- "SparseGRM"
    } else {
      optionGRM <- "DenseGRM"
    }
  }

  # Validate GenoFile and SparseGRMFile for WtCoxG method
  if (method == "WtCoxG") {
    # GenoFile validation (required)
    if (is.null(GenoFile)) {
      stop("Argument 'GenoFile' is required for method 'WtCoxG'.")
    }
    if (!is.character(GenoFile) || length(GenoFile) != 1) {
      stop("Argument 'GenoFile' should be a character string (file path).")
    }
    if (!file.exists(GenoFile)) {
        stop("Cannot find GenoFile: ", GenoFile)
    }

    # SparseGRMFile validation (optional)
    if (!is.null(SparseGRMFile)) {
      if (!is.character(SparseGRMFile) || length(SparseGRMFile) != 1) {
        stop("Argument 'SparseGRMFile' should be a character string (file path).")
      }
      if (!file.exists(SparseGRMFile)) {
        stop("Cannot find SparseGRMFile: ", SparseGRMFile)
      }
      optionGRM <- "SparseGRM"
    }
  }

  # ========== Validate and configure the control list ==========

  # Validate control parameter (optional, default NULL)
  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be a list of control parameters.")
  }

  # Method-specific control parameter validation and default setting
  control <- switch(method,
    POLMM = checkControl.NullModel.POLMM(control, traitType, optionGRM),
    SPACox = checkControl.NullModel.SPACox(control, traitType),
    SPAmix = checkControl.NullModel.SPAmix(control, traitType),
    WtCoxG = checkControl.NullModel.WtCoxG(control, traitType),
    stop("Unsupported method: ", method)
  )

  # ========== Print all parameters ==========
  
  params <- list(
    Method = method,
    `Trait type` = traitType,
    `Formula` = deparse(formula),
    `Subject ID column` = subjIDcol,
    `Sample size of input` = nrow(data),
    `Sample size in subset` = ifelse(is.null(subset), "All subjects", length(subset)),
    `Genotype file` = ifelse(is.null(GenoFile), "Not provided", GenoFile),
    `Genotype index file` = ifelse(is.null(GenoFileIndex), "Default", GenoFileIndex),
    `Sparse GRM file` = ifelse(is.null(SparseGRMFile), "Not provided", SparseGRMFile),
    `GRM option` = ifelse(exists("optionGRM"), optionGRM, "Not applicable")
  )
  .printParameters("Parameters for Null Model Fitting", params, control)

  # ========== Process formula and data to extract response and design matrix ==========

  # Handle multiple response variables for SPAmix/SPACox with residual input
  mf <- match.call(expand.dots = FALSE)                       # call object

  # Check if formula left side contains multiple variables (only for SPAmix/SPACox residual analysis)
  if (method %in% c("SPAmix", "SPACox") && traitType == "Residual") {
    LeftInFormula <- deparse(formula[[2]])                    # character
    LeftIncludesAdd <- grepl("\\+", LeftInFormula)            # logical

    if (LeftIncludesAdd) {
      nInLeft <- length(strsplit(LeftInFormula, "\\+")[[1]]) # integer
      .message("%s method supports multiple response variables of model residuals", method)

      # Parse and reconstruct formula for multiple responses
      RightInFormula <- deparse(formula[[3]])                 # character
      NewLeftInFormula <- paste0("paste(", gsub("\\+", ",", LeftInFormula), ")") # character
      NewRightInFormula <- paste0(RightInFormula, collapse = " ") # character
      mf$formula <- as.formula(paste(NewLeftInFormula, "~", NewRightInFormula)) # formula
    }
  }

  # Extract model frame components (formula, data, subset)
  m <- match(                                                 # integer vector
    x = c("formula", "data", "subset"),
    table = names(mf), nomatch = 0L
  )

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())             # data.frame

  mt <- attr(x = mf, which = "terms")                        # terms object

  # Extract response variable(s) and design matrix from model frame
  response <- model.response(mf)                              # vector or matrix
  designMat <- model.matrix(object = mt, data = mf)          # matrix
  
  # Extract subject IDs from the original data using subjIDcol
  if (is.null(subset)) {
    subjData <- data[[subjIDcol]]                            # character vector
  } else {
    subjData <- data[[subjIDcol]][subset]                    # character vector
  }

  # Validate subject IDs
  if (!is.character(subjData) || length(subjData) == 0) {
    stop("Column '", subjIDcol, "' should contain character subject IDs.")
  }

  if (any(duplicated(subjData))) {
    stop("Column '", subjIDcol, "' contains duplicated subject IDs, which are not supported.")
  }

  # Handle multiple response variables for SPAmix/SPACox residual analysis
  if (method %in% c("SPAmix", "SPACox") && traitType == "Residual") {
    LeftInFormula <- deparse(formula[[2]])                    # character
    LeftIncludesAdd <- grepl("\\+", LeftInFormula)            # logical
    
    if (LeftIncludesAdd) {
      nInLeft <- length(strsplit(LeftInFormula, "\\+")[[1]]) # integer
      
      # Create pattern for missing values across all phenotypes
      noValueInAnyPheno <- paste(rep(NA, nInLeft), collapse = " ") # character
      posNoValue <- which(response == noValueInAnyPheno)      # integer vector
      response.temp <- response                               # character vector

      # Remove subjects without any phenotype data
      if (length(posNoValue) > 0) {
        .message("Removing %d subjects with no phenotype", length(posNoValue))
        response.temp <- response[-posNoValue]                # character vector
        designMat <- designMat[-posNoValue, , drop = FALSE]  # matrix
        subjData <- subjData[-posNoValue]                    # character vector
      }

      # Convert response to matrix format for multiple residuals
      nRes <- length(response.temp)                           # integer
      response <- matrix(NA, nRes, nInLeft)                  # matrix
      for (i in 1:nRes) {
        response[i, ] <- as.numeric(unlist(strsplit(response.temp[i], split = " "))) # numeric vector
      }
      
      # Set response class for residual analysis
      class(response) <- "Residual"
    } else {
      # Single response variable converted to matrix format
      response <- matrix(response, ncol = 1)                  # matrix
      class(response) <- "Residual"
    }
  } else if (traitType == "Residual") {
    # Single response variable for other methods with Residual trait type
    response <- matrix(response, ncol = 1)                    # matrix
    class(response) <- "Residual"
  }

  # Remove intercept column if present (automatically added during model fitting)
  if (colnames(designMat)[1] == "(Intercept)") {
    designMat <- designMat[, -1, drop = FALSE]               # matrix
  }

  nData <- length(subjData)                                   # integer
  .message("Number of subjects with phenotype: %d", nData)

  # ========== Configure GRM ==========

  # Only certain methods require genotype information for GRM construction
  if (method %in% c("POLMM")) {
    # Setup GRM options and extract genotype file information
    if (is.null(GenoFile)) {
      stop("Argument of 'GenoFile' is required to estimate variance ratio.")
    }

    genoList <- setGenoInput(GenoFile, GenoFileIndex, subjData) # list

    if (!is.null(SparseGRMFile)) {
      .message("Using sparse GRM for null model")
      SparseGRM <- data.table::fread(SparseGRMFile)          # data.frame
      SparseGRM <- as.data.frame(SparseGRM)                  # data.frame

      KinMatListR <- updateSparseGRM(SparseGRM, subjData)    # list

      # Configure sparse GRM in C++ backend (see Main.cpp)
      setSparseGRMInCPP(
        t_KinMatListR = KinMatListR  # list: Sparse kinship matrix (locations, values, nSubj)
      )
    } else {
      .message("Using dense GRM for null model")
      if (genoList$genoType != "PLINK") {
        stop(
          "If DenseGRM is used when fitting a null model, ",
          "then only PLINK format is supported."
        )
      }

      memoryChunk <- 2                                        # numeric: memory in GB
      minMafGRM <- 0.01                                       # numeric
      maxMissingGRM <- 0.1                                    # numeric

      # Configure dense GRM in C++ backend (see Main.cpp)
      setDenseGRMInCPP(
        t_memoryChunk = memoryChunk,      # numeric: Memory allocation in GB for GRM
        t_minMafGRM = minMafGRM,          # numeric: Min MAF for variants in GRM
        t_maxMissingGRM = maxMissingGRM   # numeric: Max missing rate for GRM variants
      )
    }

    genoType <- genoList$genoType                             # character: "PLINK" or "BGEN"
    markerInfo <- genoList$markerInfo                         # data.frame
  }

  # ========== Fit null model ==========

  objNull <- switch(method,
    POLMM = fitNullModel.POLMM(response, designMat, subjData, control,
                               optionGRM, genoType, markerInfo),
    SPACox = fitNullModel.SPACox(response, designMat, subjData, control),
    SPAmix = fitNullModel.SPAmix(response, designMat, subjData, control),
    WtCoxG = fitNullModel.WtCoxG(response, designMat, subjData, control, data,
                                 GenoFile, GenoFileIndex, SparseGRMFile, ...)
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
