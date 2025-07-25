#' Fit a null model to estimate parameters and residuals
#'
#' We fit a null model including response variable, covariates, and Genetic Relationship Matrix (GRM, if needed) to estimate parameters and residuals.
#'
#' @param formula a formula object, with the response on the left of a ~ operator and the covariates on the right. Do not add a column of intercept (i.e. a vector of ones) on the right. Missing values should be denoted by NA and the corresponding samples will be removed from analysis. Other values (e.g. -9, -999) will be treated as ordinary numeric values in analysis.
#' @param data a data.frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data.frame), containing the variables in formula. Neither a matrix nor an array will be accepted.
#' @param subset a specification of the rows to be used: defaults to all rows. This can be any valid indexing vector for the rows of data or if that is not supplied, a data frame made up of the variables used in formula.
#' @param subjData a character vector of subject IDs. Its order should be the same as the subject order in the formula and data (before subset process).
#' @param method a character: "SPACox" (check \code{\link{GRAB.SPACox}}), "POLMM" (check \code{\link{GRAB.POLMM}}), "SPAGE" (will be supported later), or "GATE" (will be supported later).
#' @param traitType a character: "binary", "ordinal" (check \code{\link{GRAB.POLMM}}), "quantitative", or "time-to-event" (check \code{\link{GRAB.SPACox}}).
#' @param GenoFile a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param GenoFileIndex additional index files corresponding to the \code{GenoFile}. If \code{NULL} (default), the same prefix as GenoFile is used. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param SparseGRMFile a character of sparseGRM file. An example is \code{system.file("SparseGRM","SparseGRM.txt",package="GRAB")}
#' @param control a list of parameters for controlling the model fitting process. For more details, please check \code{Details} section.
#' @param ... other arguments passed to or from other methods.
#' @return an R object with a class of "XXXXX_NULL_Model" in which XXXXX is the 'method' used in analysis. The following elements are required for all methods.
#' \itemize{
#'   \item{N}: Sample size in analysis
#'   \item{yVec}: Phenotype data
#'   \item{beta}: Coefficient parameters corresponding to covariates
#'   \item{subjData}: Subject IDs in analysis
#'   \item{sessionInfo}: Version information about R, the OS and attached or loaded packages.
#'   \item{Call}: A call in which all of the specified arguments are specified by their full names.
#'   \item{time}: The time when analysis is finished
#'   \item{control}: The R list of control in null model fitting
#' }
#' @details
#' \code{GRAB} package uses score testing which consists of two steps. In Step 1, function \code{GRAB.NullModel} fits a null model including response variable, covariates, and Genetic Relationship Matrix (GRM) if needed.
#' In Step 2, functions \code{\link{GRAB.Marker}} and \code{\link{GRAB.Region}} perform genome-wide marker-level analysis and region-level analysis, respectively.
#' Step 1 fits a null model to get an R object, which is passed to Step 2 for association testing. Functions of \code{\link{save}} and \code{\link{load}} can save and load the object.
#'
#' \code{GRAB} package includes multiple methods which support a wide variety of phenotypes as follows.
#' \itemize{
#'   \item \code{POLMM}: Support \code{traitType} = \code{"ordinal"}. Check \code{\link{GRAB.POLMM}} for more details.
#'   \item \code{SPACox}: Support \code{traitType} = \code{"time-to-event"} or \code{"Residual"}. Check \code{\link{GRAB.SPACox}} for more details.
#'   \item \code{SPAmix}: Support \code{traitType} = \code{"time-to-event"} or \code{"Residual"}. Check \code{\link{GRAB.SPAmix}} for more details.
#'   \item \code{SPAGRM}: Support \code{traitType} = \code{"time-to-event"} or \code{"Residual"}. Check \code{\link{GRAB.SPAGRM}} for more details.
#' }
#'
#' \code{GRAB} package supports both Dense and Sparse GRM to adjust for sample relatedness.
#' If Dense GRM is used, then \code{GenoFile} is required to construct GRM.
#' If Sparse GRM is used, then \code{SparseGRMFile} is required, whose details can be seen in \code{\link{getTempFilesFullGRM}} and \code{\link{getSparseGRM}}.
#'
#' ## The following details are about argument \code{control}
#' Argument \code{control} includes a list of parameters for controlling the null model fitting process.
#' \itemize{
#'     \item \code{maxiter}: Maximum number of iterations used to fit the null model. *(default=100)*
#'     \item \code{seed}: An integer as a random seed. Used when random process is involved. *(default=12345678)*
#'     \item \code{tolBeta}: Positive tolerance: the iterations converge when |beta - beta_old| / (|beta| + |beta_old| + tolBeta) < tolBeta. *(default=0.001)*
#'     \item \code{showInfo}: Whether to show more detailed information for trouble shooting. *(default=FALSE)*
#' }
#'
#' To adjust for sample relatedness, mixed effect model incorporates a random effect with a variance component.
#' Argument \code{control} includes additional parameters to estimate the variance component.
#' \itemize{
#'   \item \code{tau}: Initial value of the variance component (tau). *(default=0.2)*.
#'   \item \code{tolTau}: Positive tolerance: the iterations converge when |tau - tau_old| / (|tau| + |tau_old| + tolTau) < tolTau. *(default=0.002)*
#' }
#'
#' If dense GRM is used to adjust for sample relatedness, \code{GenoFile} should be PLINK files and argument \code{control} includes additional parameters as follows.
#' \itemize{
#'   \item \code{maxiterPCG}: Maximum number of iterations for PCG to converge. *(default=100)*
#'   \item \code{tolEps}: Positive tolerance for PCG to converge. *(default=1e-6)*
#'   \item \code{minMafVarRatio}: Minimal value of MAF cutoff to select markers (from PLINK files) to estimate variance ratio. *(default=0.1)*
#'   \item \code{maxMissingVarRatio}: Maximal value of missing rate cutoff to select markers (from PLINK files) to estimate variance ratio. *(default=0.1)*
#'   \item \code{nSNPsVarRatio}: Initial number of the selected markers to estimate variance ratio *(default=20)* the number will be automatically added by 10 until the coefficient of variantion (CV) of the variance ratio estimate is below CVcutoff.
#'   \item \code{CVcutoff}: Minimal cutoff of coefficient of variantion (CV) to estimate variance ratio *(default=0.0025)*
#'   \item \code{LOCO}: Whether to apply the leave-one-chromosome-out (LOCO) approach. *(default=TRUE)*
#'   \item \code{numThreads}: Number of threads (CPUs) to use. Only valid if dense GRM is used, check \code{\link[RcppParallel]{defaultNumThreads}}. *(default="auto")*
#'   \item \code{stackSize}: Stack size (in bytes) to use for worker threads. For more details, check \code{\link[RcppParallel]{setThreadOptions}}. *(default="auto")*
#'   \item \code{grainSize}: Grain size of a parallel algorithm sets a minimum chunk size for parallelization. In other words, at what point to stop processing input on separate threads. *(default=1)*
#'   \item \code{minMafGRM}: Minimal value of MAF cutoff to select markers (from PLINK files) to construct dense GRM. *(default=0.01)*
#'   \item \code{memoryChunk}: Size (Gb) for each memory chunk when reading in PLINK files. *(default=2)*
#'   \item \code{tracenrun}: Number of runs for trace estimator. *(default=30)*
#'   \item \code{maxMissingGRM}: Maximal value of missing rate to select markers (from PLINK files) to construct dense GRM. *(default=0.1)*
#'   \item \code{onlyCheckTime}: Not fit the null model, only check the computation time of reading PLINK files and running 30 KinbVec() functions. *(default=FALSE)*
#' }
#'
#' @examples
#' # For POLMM method (ordinal categorical data analysis while adjusting for sample relatedness)
#' 
#' # Step 1: fit a null model using a sparse GRM
#' # A sparse GRM file can be obtained by getSparseGRM().
#' # If SparseGRMFile isn't provided, GRAB.NullModel() will calculate dense GRM from GenoFile.
#' 
#' PhenoData <- read.table(system.file("extdata", "simuPHENO.txt", package = "GRAB"), header = TRUE)
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' SparseGRMFile <- system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#' obj.POLMM <- GRAB.NullModel(factor(OrdinalPheno) ~ AGE + GENDER,
#'   data = PhenoData, subjData = IID,
#'   method = "POLMM", traitType = "ordinal",
#'   GenoFile = GenoFile,
#'   SparseGRMFile = SparseGRMFile,
#'   control = list(showInfo = FALSE, LOCO = FALSE, tolTau = 0.2, tolBeta = 0.1)
#' )
#'
#' names(obj.POLMM)
#' obj.POLMM$tau
#'
#' # save(obj.POLMM, "obj.POLMM.RData")  # save the object for analysis in step 2
#'
#' # For SPACox method, check ?GRAB.SPACox.
#' # For SPAmix method, check ?GRAB.SPAmix.
#' # For SPAGRM method, check ?GRAB.SPAGRM
#' # For WtCoxG method, check ?GRAB.WtCoxG
#'
#' @export
GRAB.NullModel <- function(formula,
                           data = NULL,
                           subset = NULL,
                           subjData,
                           method = "SPACox",
                           traitType = "time-to-event", # "binary", "ordinal", "quantitative", "time-to-event"
                           GenoFile = NULL,
                           GenoFileIndex = NULL,
                           SparseGRMFile = NULL,
                           control = NULL,
                           ...) {
  if (missing(subjData)) {
    stop("Argument 'subjData' is required to specify the subjects IDs in 'formula' and/or 'data'.")
  }

  if (method == "wtCoxG") {
    required_args <- c("RefAfFile", "OutputFile", "SampleIDColumn", "SurvTimeColumn", "IndicatorColumn")
    dots <- list(...)
    missing_args <- required_args[!required_args %in% names(dots)]
    if (length(missing_args) > 0) {
      stop(paste("The following arguments are missing:", paste(missing_args, collapse = ", ")))
    }

    required_columns <- c(SampleIDColumn, SurvTimeColumn, IndicatorColumn)
    missing_columns <- required_columns[!required_columns %in% colnames(data)]
    if (length(missing_columns) > 0) {
      stop(paste("The following columns are missing in `data`:", paste(missing_columns, collapse = ", ")))
    }
  }

  Call <- match.call()

  #### START: formula.R

  #### input: formula, data, subset, subjData
  #### output: response, designMat, subjData

  mf <- match.call(expand.dots = FALSE)

  ### The below is to support multiple response variables for SPAmix with residuals as input

  LeftInFormula <- deparse(formula[[2]])
  LeftIncludesAdd <- grepl("\\+", LeftInFormula)

  if (LeftIncludesAdd) {
    if (method != "SPAmix" | traitType != "Residual") {
      stop("Only 'SPAmix' method with traitType of 'Residual' supports multiple responses variables in 'formula'.")
    }

    nInLeft <- length(strsplit(LeftInFormula, "\\+")[[1]])
    cat("SPAmix method supports multiple response variables of model residuals.\n")
    RightInFormula <- deparse(formula[[3]])
    NewLeftInFormla <- paste0("paste(", gsub("\\+", ",", LeftInFormula), ")")
    NewRightInFormula <- paste0(RightInFormula, collapse = " ") # c("cov1 + cov2 +", "cov3") -> "cov1 + cov2 + cov3"
    # mf$formula = as.formula(paste(NewLeftInFormla, "~", RightInFormula))
    mf$formula <- as.formula(paste(NewLeftInFormla, "~", NewRightInFormula))
  }
  ##

  m <- match(
    x = c("formula", "data", "subset", "subjData"),
    table = names(mf), nomatch = 0L
  )

  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())

  mt <- attr(x = mf, which = "terms")

  response <- model.response(mf)
  designMat <- model.matrix(object = mt, data = mf)
  subjData <- model.extract(mf, "subjData")

  ### The below is to support multiple response variables for SPAmix with residuals as input
  if (traitType == "Residual") {
    if (LeftIncludesAdd) {
      noValueInAnyPheno <- paste(rep(NA, nInLeft), collapse = " ")
      posNoValue <- which(response == noValueInAnyPheno)
      response.temp <- response

      if (length(posNoValue) > 0) {
        cat("We remove", length(posNoValue), "individuals without any phenotyeps in analysis.\n")
        response.temp <- response[-1 * posNoValue]
        designMat <- designMat[-1 * posNoValue, , drop = F]
        subjData <- subjData[-1 * posNoValue]
      }

      nRes <- length(response.temp)
      response <- matrix(NA, nRes, nInLeft)
      for (i in 1:nRes) {
        response[i, ] <- as.numeric(unlist(strsplit(response.temp[i], split = " ")))
      }
    } else {
      response <- matrix(response, ncol = 1)
    }
    class(response) <- "Residual"
  }

  if (colnames(designMat)[1] == "(Intercept)") {
    designMat <- designMat[, -1, drop = F]
  }

  nData <- length(subjData)
  cat("Number of subjects in 'formula':\t", nData, "\n")

  if (any(duplicated(subjData))) {
    stop("Duplicated subject IDs in 'subjData' is not supported!")
  }

  #### END: formula.R

  ## Only the below methods requires 'GenoFile' related information to adjust for sample relatedness
  optionGRM <- NULL
  if (method %in% c("POLMM", "GATE")) {
    objGRM <- setGRM(GenoFile, GenoFileIndex, SparseGRMFile, subjData) # Check 'SparseGRM.R'
    optionGRM <- objGRM$optionGRM
    genoType <- objGRM$genoType # "PLINK" or "BGEN"
    markerInfo <- objGRM$markerInfo # Columns: "CHROM", "POS", "ID", "REF", "ALT", "genoIndex"

    # Check 'control.R'
    control <- checkControl.NullModel(control, method, traitType, optionGRM)
    textToParse <- paste0("objNull = fitNullModel.", method, "(response, designMat, subjData, control, optionGRM, genoType, markerInfo)")
  } else {
    # Check 'control.R'
    control <- checkControl.NullModel(control, method, traitType)
    textToParse <- paste0("objNull = fitNullModel.", method, "(response, designMat, subjData, control)")
  }

  # e.g. if(method == "POLMM"){objNull = fitNullModel.POLMM(...)}  # fitNullModel.POLMM() function is in POLMM.R
  eval(parse(text = textToParse))

  objNull$subjData <- subjData

  # (BWJ: 2023-08-09): not sure if the below works?
  # objNull$N = length(subjData)

  objNull$Call <- Call
  objNull$sessionInfo <- sessionInfo()
  objNull$time <- paste0("Complete Time: ", Sys.time())
  objNull$control <- control

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

  cat("Complete the null model fitting in package GRAB:\t", objNull$time, "\n")
  return(objNull)
}
