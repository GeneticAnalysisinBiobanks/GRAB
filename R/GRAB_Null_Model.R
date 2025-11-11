#' Fit a null model for genetic association testing
#'
#' Estimates model parameters and residuals from a null model containing phenotypes,
#' covariates, and optionally a genetic relationship matrix (GRM) for subsequent
#' association testing.
#'
#' @param formula Formula with response variable(s) on the left and covariates on the right.
#'   Do not include an intercept (added automatically). Missing values should be coded as
#'   \code{NA}. Other values (e.g., -9, -999) are treated as numeric.
#' @param data Data frame containing variables in the formula.
#' @param subset Row specification for subsetting data (default: all rows).
#' @param subjData Character vector of subject IDs matching the order in formula/data
#'   before subsetting.
#' @param method Statistical method: "POLMM", "SPACox", "SPAmix", "SPAGRM", "SAGELD",
#'   or "WtCoxG".
#' @param traitType Trait type: "binary", "ordinal", "quantitative", "time-to-event",
#'   or "Residual".
#' @param GenoFile Path to genotype file (PLINK or BGEN format). Required for dense GRM
#'   construction. See \code{\link{GRAB.ReadGeno}} for details.
#' @param GenoFileIndex Index files for genotype file. If \code{NULL} (default), uses
#'   same prefix as \code{GenoFile}. See \code{\link{GRAB.ReadGeno}} for details.
#' @param SparseGRMFile Path to sparse GRM file. Alternative to dense GRM construction.
#' @param control List of control parameters. See \code{Details} for options.
#' @param ... Additional arguments for method-specific functions.

#' @return
#' List object with class "\{method\}_NULL_Model" containing:
#' \describe{
#'   \item{\code{N}}{Sample size in analysis.}
#'   \item{\code{yVec}}{Phenotype vector used in analysis.}
#'   \item{\code{beta}}{Estimated covariate effect coefficients.}
#'   \item{\code{subjData}}{Subject IDs included in analysis.}
#'   \item{\code{sessionInfo}}{R session and package information.}
#'   \item{\code{Call}}{Function call with all arguments.}
#'   \item{\code{time}}{Analysis completion timestamp.}
#'   \item{\code{control}}{Control parameters used in fitting.}
#'   \item{\code{tau}}{Variance component estimates (for mixed models).}
#'   \item{\code{SparseGRM}}{Sparse GRM (if used).}
#' }
#'
#' This object serves as input for \code{\link{GRAB.Marker}} and \code{\link{GRAB.Region}}.
#'   \code{\link{GRAB.Marker}} and \code{\link{GRAB.Region}}.
#'
#' @details
#' GRAB uses score testing in two steps:
#' \enumerate{
#'   \item \code{GRAB.NullModel} fits a null model with phenotypes, covariates, and GRM.
#'   \item \code{\link{GRAB.Marker}} and \code{\link{GRAB.Region}} perform association testing.
#' }
#'
#' The null model object can be saved using \code{\link{save}} and loaded with \code{\link{load}}.
#'
#' **Supported Methods and Trait Types:**
#' \itemize{
#'   \item \code{POLMM}: Ordinal traits. See \code{\link{GRAB.POLMM}}.
#'   \item \code{SPACox}: Time-to-event or residual traits. See \code{\link{GRAB.SPACox}}.
#'   \item \code{SPAmix}: Time-to-event or residual traits. See \code{\link{GRAB.SPAmix}}.
#'   \item \code{SPAGRM}: Various trait types with sparse GRM. See \code{\link{GRAB.SPAGRM}}.
#'   \item \code{SAGELD}: Various trait types with SAGELD algorithm. See \code{\link{GRAB.SAGELD}}.
#'   \item \code{WtCoxG}: Time-to-event traits. See \code{\link{GRAB.WtCoxG}}.
#' }
#'
#' **Genetic Relationship Matrix (GRM):**
#' GRAB supports both dense and sparse GRM for relatedness adjustment:
#' \itemize{
#'   \item Dense GRM: Constructed from \code{GenoFile} (PLINK/BGEN format).
#'   \item Sparse GRM: Pre-computed matrix from \code{SparseGRMFile}.
#' }
#' See \code{\link{getTempFilesFullGRM}} and \code{\link{getSparseGRM}} for details.
#'
#' ## Control Parameters
#'
#' **Convergence and Iteration:**
#' \describe{
#'   \item{\code{maxiter}}{Maximum iterations for null model fitting (default: 100).}
#'   \item{\code{seed}}{Random seed for reproducibility (default: 12345678).}
#'   \item{\code{tolBeta}}{Convergence tolerance for fixed effects (default: 0.001).}
#'   \item{\code{showInfo}}{Show detailed fitting information (default: \code{FALSE}).}
#' }
#'
#' **Variance Components:**
#' \describe{
#'   \item{\code{tau}}{Initial variance component value (default: 0.2).}
#'   \item{\code{tolTau}}{Convergence tolerance for variance components (default: 0.002).}
#' }
#'
#' **Dense GRM Parameters:**
#' \describe{
#'   \item{\code{maxiterPCG}}{Maximum PCG iterations (default: 100).}
#'   \item{\code{tolEps}}{PCG convergence tolerance (default: 1e-6).}
#'   \item{\code{minMafVarRatio}}{Minimum MAF for variance ratio estimation (default: 0.1).}
#'   \item{\code{maxMissingVarRatio}}{Maximum missing rate for variance ratio (default: 0.1).}
#'   \item{\code{nSNPsVarRatio}}{Initial markers for variance ratio (default: 20).}
#'   \item{\code{CVcutoff}}{CV cutoff for variance ratio estimation (default: 0.0025).}
#'   \item{\code{LOCO}}{Use leave-one-chromosome-out (default: \code{TRUE}).}
#'   \item{\code{minMafGRM}}{Minimum MAF for GRM construction (default: 0.01).}
#'   \item{\code{memoryChunk}}{Memory chunk size in GB for PLINK reading (default: 2).}
#'   \item{\code{tracenrun}}{Trace estimator runs (default: 30).}
#'   \item{\code{maxMissingGRM}}{Maximum missing rate for GRM (default: 0.1).}
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
GRAB.NullModel <- function(
  formula,
  data,
  subset = NULL,
  subjData,
  method,
  traitType, # "binary", "ordinal", "quantitative", "time-to-event", "Residual"
  GenoFile = NULL,
  GenoFileIndex = NULL,
  SparseGRMFile = NULL,
  control = NULL,
  ...
) {
  # Validate required arguments
  if (missing(subjData)) {
    stop(
      "Argument 'subjData' is required to specify subject IDs in 'formula' and/or 'data'."
    )
  }

  # Store function call for reproducibility and debugging
  Call <- match.call()

  #### Process formula and data to extract response and design matrix ####
  # Handle multiple response variables for SPAmix with residual input
  mf <- match.call(expand.dots = FALSE)

  # Check if formula left side contains multiple variables (for SPAmix residual analysis)
  LeftInFormula <- deparse(formula[[2]])
  LeftIncludesAdd <- grepl("\\+", LeftInFormula)

  if (LeftIncludesAdd) {
    if (method != "SPAmix" || traitType != "Residual") {
      stop(
        "Only 'SPAmix' method with traitType 'Residual' supports ",
        "multiple response variables in formula."
      )
    }

    nInLeft <- length(strsplit(LeftInFormula, "\\+")[[1]])
    .message("SPAmix method supports multiple response variables of model residuals")

    # Parse and reconstruct formula for multiple responses
    RightInFormula <- deparse(formula[[3]])
    NewLeftInFormula <- paste0("paste(", gsub("\\+", ",", LeftInFormula), ")")
    NewRightInFormula <- paste0(RightInFormula, collapse = " ")
    mf$formula <- as.formula(paste(NewLeftInFormula, "~", NewRightInFormula))
  }

  # Extract model frame components (formula, data, subset, subjData)
  m <- match(
    x = c("formula", "data", "subset", "subjData"),
    table = names(mf), nomatch = 0L
  )

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())

  mt <- attr(x = mf, which = "terms")

  # Extract response variable(s) and design matrix from model frame
  response <- model.response(mf)
  designMat <- model.matrix(object = mt, data = mf)
  subjData <- model.extract(mf, "subjData")

  # Handle multiple response variables for SPAmix residual analysis
  if (traitType == "Residual") {
    if (LeftIncludesAdd) {
      # Create pattern for missing values across all phenotypes
      noValueInAnyPheno <- paste(rep(NA, nInLeft), collapse = " ")
      posNoValue <- which(response == noValueInAnyPheno)
      response.temp <- response

      # Remove subjects without any phenotype data
      if (length(posNoValue) > 0) {
        .message("Removing %d subjects with no phenotype", length(posNoValue))
        response.temp <- response[-posNoValue]
        designMat <- designMat[-posNoValue, , drop = FALSE]
        subjData <- subjData[-posNoValue]
      }

      # Convert response to matrix format for multiple residuals
      nRes <- length(response.temp)
      response <- matrix(NA, nRes, nInLeft)
      for (i in 1:nRes) {
        response[i, ] <- as.numeric(unlist(strsplit(response.temp[i], split = " ")))
      }
    } else {
      # Single response variable converted to matrix format
      response <- matrix(response, ncol = 1)
    }
    # Set response class for residual analysis
    class(response) <- "Residual"
  }

  # Remove intercept column if present (automatically added during model fitting)
  if (colnames(designMat)[1] == "(Intercept)") {
    designMat <- designMat[, -1, drop = FALSE]
  }

  nData <- length(subjData)
  .message("Number of subjects with phenotype: %d", nData)

  # Validate subject ID uniqueness
  if (any(duplicated(subjData))) {
    stop("Duplicated subject IDs in 'subjData' are not supported!")
  }

  #### Configure genetic relationship matrix (GRM) for relatedness adjustment ####
  # Only certain methods require genotype information for GRM construction
  optionGRM <- NULL
  if (method %in% c("POLMM")) {
    # Setup GRM options and extract genotype file information
    # ---- BEGIN inlined: setGRM ----
    if (is.null(GenoFile)) {
      stop("Argument of 'GenoFile' is required to estimate variance ratio.")
    }

    genoList <- setGenoInput(GenoFile, GenoFileIndex, subjData) # check Geno.R for more details

    if (!is.null(SparseGRMFile)) {
      .message("Using sparse GRM for null model")
      SparseGRM <- data.table::fread(SparseGRMFile)
      SparseGRM <- as.data.frame(SparseGRM)

      KinMatListR <- updateSparseGRM(SparseGRM, subjData)

      # the following function is in Main.cpp
      setSparseGRMInCPP(KinMatListR)
      optionGRM_local <- "SparseGRM"
    } else {
      .message("Using dense GRM for null model")
      if (genoList$genoType != "PLINK") {
        stop(
          "If DenseGRM is used when fitting a null model, ",
          "then only PLINK format is supported."
        )
      }

      memoryChunk <- 2 # (GB)
      minMafGRM <- 0.01
      maxMissingGRM <- 0.1

      # the following function is in Main.cpp
      setDenseGRMInCPP(memoryChunk, minMafGRM, maxMissingGRM)
      optionGRM_local <- "DenseGRM"
    }

    objGRM <- list(
      optionGRM = optionGRM_local,
      genoType = genoList$genoType,
      markerInfo = genoList$markerInfo
    )
    # ---- END inlined: setGRM ----
    optionGRM <- objGRM$optionGRM
    genoType <- objGRM$genoType # "PLINK" or "BGEN"
    markerInfo <- objGRM$markerInfo # Columns: "CHROM", "POS", "ID", "REF", "ALT", "genoIndex"
  }

  #### Validate and set control parameters for null model fitting ####
  if (!is.null(control) && !is.list(control)) {
    stop("If specified, the argument 'control' should be an R list.")
  }

  # Method-specific control parameter validation and default setting
  if (method == "POLMM") {
    control <- checkControl.NullModel.POLMM(control, traitType, optionGRM)
  } else if (method == "SPACox") {
    control <- checkControl.NullModel.SPACox(control, traitType)
  } else if (method == "SPAmix") {
    control <- checkControl.NullModel.SPAmix(control, traitType)
  } else if (method == "WtCoxG") {
    control <- checkControl.NullModel.WtCoxG(control, traitType)
  }

  # Display control parameters for user verification
  .message("Control parameters for null model fitting:")
  tmp <- capture.output(str(control))
  for (line in tmp) {
    if (startsWith(line, " $")) {
      message(sub("^ \\$", strrep(" ", 8), line))
    }
  }

  #### Fit the null model using the appropriate method ####
  # Dispatch to method-specific fitting functions based on selected method
  if (method == "POLMM") {
    objNull <- fitNullModel.POLMM(response, designMat, subjData, control,
                                  optionGRM, genoType, markerInfo)
  } else if (method == "SPACox") {
    objNull <- fitNullModel.SPACox(response, designMat, subjData, control)
  } else if (method == "SPAmix") {
    objNull <- fitNullModel.SPAmix(response, designMat, subjData, control)
  } else if (method == "WtCoxG") {
    objNull <- fitNullModel.WtCoxG(
      response, designMat, subjData, control, data,
      GenoFile, GenoFileIndex, SparseGRMFile, ...
    )
  }

  # Add metadata to the fitted null model object for downstream analysis
  objNull$subjData <- subjData
  objNull$Call <- Call
  objNull$sessionInfo <- sessionInfo()
  objNull$time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  objNull$control <- control

  .message("Successfully finished fitting the null model")
  return(objNull)
}
