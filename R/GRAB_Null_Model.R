#' Fit a null model to estimate parameters and residuals
#'
#' Fits a null model that includes response variables, covariates, and optionally a Genetic 
#' Relationship Matrix (GRM) to estimate model parameters and residuals for subsequent 
#' association testing.
#'
#' @param formula A formula object with the response on the left of a ~ operator and 
#'   covariates on the right. Do not include an intercept term (i.e., a vector of ones) 
#'   on the right side. Missing values should be coded as NA and corresponding samples 
#'   will be excluded from analysis. Other values (e.g., -9, -999) will be treated as 
#'   ordinary numeric values.
#' @param data A data.frame, list, or environment (or object coercible by 
#'   \code{\link{as.data.frame}} to a data.frame) containing the variables in the formula. 
#'   Neither a matrix nor an array will be accepted.
#' @param subset A specification of the rows to be used; defaults to all rows. This can be 
#'   any valid indexing vector for the rows of data, or if data is not supplied, a data frame 
#'   made up of the variables used in the formula.
#' @param subjData A character vector of subject IDs. The order should match the subject 
#'   order in the formula and data (before any subset processing).
#' @param method A character string specifying the statistical method: "POLMM" 
#'   (see \code{\link{GRAB.POLMM}}), "SPACox" (see \code{\link{GRAB.SPACox}}), 
#'   "SPAmix" (see \code{\link{GRAB.SPAmix}}), or "WtCoxG" (see \code{\link{GRAB.WtCoxG}}).
#' @param traitType A character string specifying the trait type: "binary", "ordinal", 
#'   "quantitative", or "time-to-event".
#' @param GenoFile A character string specifying the genotype file path. Currently, two 
#'   genotype formats are supported: PLINK and BGEN. See \code{\link{GRAB.ReadGeno}} for details.
#' @param GenoFileIndex Additional index files corresponding to \code{GenoFile}. 
#'   If \code{NULL} (default), the same prefix as GenoFile is used. 
#'   See \code{\link{GRAB.ReadGeno}} for details.
#' @param SparseGRMFile A character string specifying the sparse GRM file path. An example is 
#'   \code{system.file("SparseGRM","SparseGRM.txt",package="GRAB")}.
#' @param control A list of parameters for controlling the model fitting process. 
#'   See the \code{Details} section for comprehensive information.
#' @param ... Additional arguments passed to or from other methods.
#' @return A list object with class "XXXXX_NULL_Model" where XXXXX is the specified \code{method}.
#'   The returned object contains the following components:
#'   \describe{
#'     \item{\code{N}}{Sample size in analysis}
#'     \item{\code{yVec}}{Phenotype data vector}
#'     \item{\code{beta}}{Coefficient parameters corresponding to covariates}
#'     \item{\code{subjData}}{Subject IDs included in analysis}
#'     \item{\code{sessionInfo}}{Version information about R, OS, and attached packages}
#'     \item{\code{Call}}{Function call with all specified arguments by their full names}
#'     \item{\code{time}}{Timestamp when analysis was completed}
#'     \item{\code{control}}{Control parameters used for null model fitting}
#'     \item{\code{tau}}{Estimated variance components (if using mixed models)}
#'     \item{\code{SparseGRM}}{Sparse genetic relationship matrix (if specified)}
#'   }
#'   This object serves as input for downstream association testing with 
#'   \code{\link{GRAB.Marker}} and \code{\link{GRAB.Region}}.
#'
#' @details The \code{GRAB} package uses score testing which consists of two steps:
#' \enumerate{
#'   \item \code{GRAB.NullModel} fits a null model including response variable, covariates, and
#'     Genetic Relationship Matrix (GRM) if needed
#'   \item \code{\link{GRAB.Marker}} and \code{\link{GRAB.Region}} perform genome-wide marker-level
#'     analysis and region-level analysis, respectively
#' }
#'
#' Step 1 fits a null model to get an R object, which is passed to Step 2 for association testing.
#' Functions \code{\link{save}} and \code{\link{load}} can save and load this object.
#'
#' \code{GRAB} package includes multiple methods which support a wide variety of phenotypes as follows.
#' \itemize{
#'   \item \code{POLMM}: Support \code{traitType} = \code{"ordinal"}.
#'     Check \code{\link{GRAB.POLMM}} for more details.
#'   \item \code{SPACox}: Support \code{traitType} = \code{"time-to-event"} or
#'     \code{"Residual"}. Check \code{\link{GRAB.SPACox}} for more details.
#'   \item \code{SPAmix}: Support \code{traitType} = \code{"time-to-event"} or
#'     \code{"Residual"}. Check \code{\link{GRAB.SPAmix}} for more details.
#'   \item \code{WtCoxG}: Support \code{traitType} = \code{"time-to-event"}. 
#'     Check \code{\link{GRAB.WtCoxG}} for more details.
#' }
#'
#' The \code{GRAB} package supports both Dense and Sparse GRM to adjust for sample relatedness.
#' If Dense GRM is used, then \code{GenoFile} is required to construct GRM.
#' If Sparse GRM is used, then \code{SparseGRMFile} is required. See
#' \code{\link{getTempFilesFullGRM}} and \code{\link{getSparseGRM}} for details.
#'
#' ## Control Parameters
#' The \code{control} argument includes a list of parameters for controlling the null model
#' fitting process:
#'
#' **Basic Parameters:**
#' \describe{
#'   \item{\code{maxiter}}{Maximum number of iterations used to fit the null model (default: 100)}
#'   \item{\code{seed}}{Random seed for reproducible results (default: 12345678)}
#'   \item{\code{tolBeta}}{Tolerance for fixed effects convergence:
#'     |beta - beta_old| / (|beta| + |beta_old| + tolBeta) < tolBeta (default: 0.001)}
#'   \item{\code{showInfo}}{Whether to show detailed information for troubleshooting (default: FALSE)}
#' }
#'
#' **Variance Component Parameters:**
#' \describe{
#'   \item{\code{tau}}{Initial value of the variance component (default: 0.2)}
#'   \item{\code{tolTau}}{Tolerance for variance component convergence:
#'     |tau - tau_old| / (|tau| + |tau_old| + tolTau) < tolTau (default: 0.002)}
#' }
#'
#' **Dense GRM Parameters (when using PLINK files):**
#' \describe{
#'   \item{\code{maxiterPCG}}{Maximum iterations for Preconditioned Conjugate Gradient 
#'     (default: 100)}
#'   \item{\code{tolEps}}{Tolerance for PCG convergence (default: 1e-6)}
#'   \item{\code{minMafVarRatio}}{Minimum MAF for markers used in variance ratio estimation 
#'     (default: 0.1)}
#'   \item{\code{maxMissingVarRatio}}{Maximum missing rate for markers used in variance ratio 
#'     estimation (default: 0.1)}
#'   \item{\code{nSNPsVarRatio}}{Initial number of markers for variance ratio estimation (default: 20)}
#'   \item{\code{CVcutoff}}{Maximum coefficient of variation for variance ratio estimation (default: 0.0025)}
#'   \item{\code{LOCO}}{Whether to apply leave-one-chromosome-out approach (default: TRUE)}
#'   \item{\code{stackSize}}{Stack size (bytes) for worker threads (default: "auto")}
#'   \item{\code{grainSize}}{Minimum chunk size for parallelization (default: 1)}
#'   \item{\code{minMafGRM}}{Minimum MAF for markers used in dense GRM construction (default: 0.01)}
#'   \item{\code{memoryChunk}}{Memory chunk size (GB) when reading PLINK files (default: 2)}
#'   \item{\code{tracenrun}}{Number of runs for trace estimator (default: 30)}
#'   \item{\code{maxMissingGRM}}{Maximum missing rate for markers used in dense GRM construction (default: 0.1)}
#'   \item{\code{onlyCheckTime}}{Only check computation time without fitting model (default: FALSE)}
#' }
#'
#' @examples
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- read.table(PhenoFile, header = TRUE)
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' 
#' # Fit a null model using POLMM with a dense GRM constructed from PLINK files.
#' Sys.setenv(RCPP_PARALLEL_NUM_THREADS = 2) # Limit threads for CRAN checks (optional for users).
#' 
#' obj.POLMM <- GRAB.NullModel(
#'   formula = factor(OrdinalPheno) ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjData = IID,
#'   method = "POLMM",
#'   traitType = "ordinal",
#'   GenoFile = GenoFile,
#'   control = list(showInfo = FALSE, LOCO = FALSE, tolTau = 0.2, tolBeta = 0.1)
#' )
#'
#' names(obj.POLMM)
#'
#' # Fit a null model using POLMM with a sparse GRM pre-calculated by getSparseGRM()
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#'
#' obj.POLMM <- GRAB.NullModel(
#'   formula = factor(OrdinalPheno) ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjData = IID,
#'   method = "POLMM",
#'   traitType = "ordinal",
#'   GenoFile = GenoFile,
#'   SparseGRMFile = SparseGRMFile,
#'   control = list(showInfo = FALSE, LOCO = FALSE, tolTau = 0.2, tolBeta = 0.1)
#' )
#' 
#' # Save the null model object for downstream analysis
#' OutputFile <- file.path(tempdir(), "objPOLMMnull.RData")
#' save(obj.POLMM, file = OutputFile)
#' 
#' # For SPACox method, check ?GRAB.SPACox
#' # For SPAmix method, check ?GRAB.SPAmix
#' # For SPAGRM method, check ?GRAB.SPAGRM
#' # For WtCoxG method, check ?GRAB.WtCoxG
#'
GRAB.NullModel <- function(formula,
                           data = NULL,
  subset = NULL,
  subjData,
  method,
  traitType, # "binary", "ordinal", "quantitative", "time-to-event"
  GenoFile = NULL,
  GenoFileIndex = NULL,
  SparseGRMFile = NULL,
  control = NULL,
  ...
) {
  # Validate required arguments
  if (missing(subjData)) {
    stop("Argument 'subjData' is required to specify the subject IDs in 'formula' and/or 'data'.")
  }
  
  # Validate method-specific requirements for wtCoxG
  if (method == "wtCoxG") {
    required_args <- c("RefAfFile", "OutputFile", "SampleIDColumn", "SurvTimeColumn", "IndicatorColumn")
    dots <- list(...)
    missing_args <- required_args[!required_args %in% names(dots)]
    if (length(missing_args) > 0) {
      stop("The following arguments are missing for wtCoxG method: ", 
           paste(missing_args, collapse = ", "))
    }

    required_columns <- c(SampleIDColumn, SurvTimeColumn, IndicatorColumn)
    missing_columns <- required_columns[!required_columns %in% colnames(data)]
    if (length(missing_columns) > 0) {
      stop("The following columns are missing in 'data': ", 
           paste(missing_columns, collapse = ", "))
    }
  }

  # Store function call for output
  Call <- match.call()

  #### START: Handle formula and data processing ####
  # Input: formula, data, subset, subjData
  # Output: response, designMat, subjData

  mf <- match.call(expand.dots = FALSE)

  # Support multiple response variables for SPAmix with residuals as input
  LeftInFormula <- deparse(formula[[2]])
  LeftIncludesAdd <- grepl("\\+", LeftInFormula)

  if (LeftIncludesAdd) {
    if (method != "SPAmix" || traitType != "Residual") {
      stop("Only 'SPAmix' method with traitType of 'Residual' supports multiple response variables in 'formula'.")
    }

    nInLeft <- length(strsplit(LeftInFormula, "\\+")[[1]])
    message("SPAmix method supports multiple response variables of model residuals.")
    
    # Parse and reconstruct formula for multiple responses
    RightInFormula <- deparse(formula[[3]])
    NewLeftInFormla <- paste0("paste(", gsub("\\+", ",", LeftInFormula), ")")
    NewRightInFormula <- paste0(RightInFormula, collapse = " ")
    # Convert "cov1 + cov2 +" and "cov3" to "cov1 + cov2 + cov3"
    mf$formula <- as.formula(paste(NewLeftInFormla, "~", NewRightInFormula))
  }

  # Match and extract components from model frame
  m <- match(
    x = c("formula", "data", "subset", "subjData"),
    table = names(mf), nomatch = 0L
  )

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())

  mt <- attr(x = mf, which = "terms")

  # Extract response variable and design matrix
  response <- model.response(mf)
  designMat <- model.matrix(object = mt, data = mf)
  subjData <- model.extract(mf, "subjData")

  # Handle multiple response variables for SPAmix with residuals
  if (traitType == "Residual") {
    if (LeftIncludesAdd) {
      # Create pattern for missing values across all phenotypes
      noValueInAnyPheno <- paste(rep(NA, nInLeft), collapse = " ")
      posNoValue <- which(response == noValueInAnyPheno)
      response.temp <- response

      # Remove individuals without any phenotype data
      if (length(posNoValue) > 0) {
        message("Removing ", length(posNoValue), 
                " individuals without any phenotype data in analysis.")
        response.temp <- response[-posNoValue]
        designMat <- designMat[-posNoValue, , drop = FALSE]
        subjData <- subjData[-posNoValue]
      }

      # Convert response to matrix format
      nRes <- length(response.temp)
      response <- matrix(NA, nRes, nInLeft)
      for (i in 1:nRes) {
        response[i, ] <- as.numeric(unlist(strsplit(response.temp[i], split = " ")))
      }
    } else {
      response <- matrix(response, ncol = 1)
    }
    # Set response class for residual analysis
    class(response) <- "Residual"
  }

  # Remove intercept column if present (will be added automatically in model fitting)
  if (colnames(designMat)[1] == "(Intercept)") {
    designMat <- designMat[, -1, drop = FALSE]
  }

  nData <- length(subjData)
  message("Number of subjects in 'formula':\t", nData)

  # Check for duplicate subject IDs
  if (any(duplicated(subjData))) {
    stop("Duplicated subject IDs in 'subjData' are not supported!")
  }

  #### END: Handle formula and data processing ####

  #### START: Setup genetic relationship matrix (GRM) ####
  # Only certain methods require 'GenoFile' information to adjust for sample relatedness
  optionGRM <- NULL
  if (method %in% c("POLMM")) {
    # Setup GRM options and extract genotype information
    objGRM <- setGRM(GenoFile, GenoFileIndex, SparseGRMFile, subjData) # Check 'SparseGRM.R'
    optionGRM <- objGRM$optionGRM
    genoType <- objGRM$genoType # "PLINK" or "BGEN"
    markerInfo <- objGRM$markerInfo # Columns: "CHROM", "POS", "ID", "REF", "ALT", "genoIndex"

    # Validate and set control parameters (see 'control.R')
    control <- checkControl.NullModel(control, method, traitType, optionGRM)
    
    # Prepare function call for null model fitting with GRM
    textToParse <- paste0("objNull = fitNullModel.", method,
                          "(response, designMat, subjData, control, optionGRM, genoType, markerInfo)")
  } else {
    # Validate and set control parameters for methods without GRM (see 'control.R')
    control <- checkControl.NullModel(control, method, traitType)
    
    # Prepare function call for null model fitting without GRM
    textToParse <- paste0("objNull = fitNullModel.", method, 
                          "(response, designMat, subjData, control)")
  }

  #### START: Fit the null model ####
  # Dynamically call the appropriate method-specific fitting function
  # e.g., if method == "POLMM", then call fitNullModel.POLMM()
  eval(parse(text = textToParse))

  # Add additional information to the fitted null model object
  objNull$subjData <- subjData
  objNull$Call <- Call
  objNull$sessionInfo <- sessionInfo()
  objNull$time <- paste0("Analysis completed at: ", Sys.time())
  objNull$control <- control

  # Special handling for WtCoxG method
  if (method == "WtCoxG") {
    objNull$mergeGenoInfo <- TestforBatchEffect(
      objNull = objNull,
      data = data,
      GenoFile = GenoFile,
      GenoFileIndex = GenoFileIndex,
      SparseGRMFile = SparseGRMFile,
      ...
    )
  }

  message("Successfully completed null model fitting in GRAB package:\t", objNull$time)
  return(objNull)
}
