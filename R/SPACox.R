#' SPACox method in GRAB package
#'
#' SPACox method is an empirical approach to analyzing complex traits
#' (including but not limited to time-to-event trait) for unrelated samples
#' in a large-scale biobank.
#'
#' @details
#' Additional list of \code{control} in \code{GRAB.NullModel()} function.
#'
#' Additional list of \code{control} in \code{GRAB.Marker()} function.
#'
#' @return No return value, called for side effects.
#'
#' @examples
#' # Step 1: fit a null model
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' obj.SPACox <- GRAB.NullModel(
#'   survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjData = IID,
#'   method = "SPACox",
#'   traitType = "time-to-event"
#' )
#'
#' # Using model residuals performs exactly the same as the above. Note that
#' # confounding factors are still required in the right of the formula.
#' obj.coxph <- survival::coxph(
#'   survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
#'   data = PhenoData,
#'   x = TRUE
#' )
#'
#' obj.SPACox <- GRAB.NullModel(
#'   obj.coxph$residuals ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjData = IID,
#'   method = "SPACox",
#'   traitType = "Residual"
#' )
#'
#' # Step 2: conduct score test
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputDir <- tempdir()
#' OutputFile <- file.path(OutputDir, "Results_SPACox.txt")
#' GRAB.Marker(
#'   objNull = obj.SPACox,
#'   GenoFile = GenoFile,
#'   OutputFile = OutputFile,
#'   control = list(outputColumns = "zScore")
#' )
#' data.table::fread(OutputFile)
#'
GRAB.SPACox <- function() {
  .message("Using SPACox method - see ?GRAB.SPACox for details")
}


#' @title Validate control parameters for SPACox null model
#' @param control List of control parameters for null model fitting.
#' @param traitType Character string specifying the trait type.
#' @return Updated control list with validated parameters and defaults.
#' @keywords internal
checkControl.NullModel.SPACox <- function(
  control,
  traitType
) {
  if (!traitType %in% c("time-to-event", "Residual")) {
    stop("For 'SPACox' method, only traitType of 'time-to-event' or 'Residual' is supported.")
  }

  default.control <- list(
    range = c(-100, 100),
    length.out = 10000
  )

  control <- updateControl(control, default.control)

  # check the parameters
  range <- control$range
  length.out <- control$length.out

  if (range[1] >= -50 || range[2] <= 50 || length.out <= 1000) {
    stop("We suggest setting argument 'control$range=c(-100,100)' and 'control$length.out=10000'.")
  }

  if (range[2] != -1 * range[1]) {
    stop("range[2] should be -1*range[1]")
  }

  return(control)
}


#' @title Fit null model for SPACox method
#' @param response Response variable (Surv or Residual object).
#' @param designMat Design matrix containing covariates.
#' @param subjData Vector of subject identifiers.
#' @param control List of control parameters.
#' @param ... Additional arguments passed to internal functions.
#' @return List containing fitted null model components.
#' @keywords internal
fitNullModel.SPACox <- function(
  response,
  designMat,
  subjData,
  control,
  ...
) {
  if (!(inherits(response, "Surv") || inherits(response, "Residual"))) {
    stop("For SPAcox, the response variable should be of class 'Surv' or 'Residual'.")
  }

  if (inherits(response, "Surv")) {
    formula <- response ~ designMat
    obj.coxph <- survival::coxph(formula, x = TRUE, ...)

    y <- obj.coxph$y
    yVec <- y[, ncol(y)]

    mresid <- obj.coxph$residuals
    Cova <- designMat
  }

  if (inherits(response, "Residual")) {
    yVec <- mresid <- response
    Cova <- designMat
  }

  range <- control$range
  length.out <- control$length.out

  if (length(mresid) != length(subjData)) {
    stop("Please check the consistency between 'formula' and 'subjData'.")
  }

  ### Get the covariate matrix to adjust for genotype
  X <- cbind(1, Cova)
  X.invXX <- X %*% solve(t(X) %*% X)
  tX <- t(X)

  ### calculate empirical CGF for martingale residuals
  idx0 <- qcauchy(1:length.out / (length.out + 1))
  idx1 <- idx0 * max(range) / max(idx0)

  cumul <- NULL
  .message("Calculating empirical CGF for martingale residuals")
  c <- 0
  for (i in idx1) {
    c <- c + 1
    t <- i
    e_resid <- exp(mresid * t)
    M0 <- mean(e_resid)
    M1 <- mean(mresid * e_resid)
    M2 <- mean(mresid^2 * e_resid)
    K0 <- log(M0)
    K1 <- M1 / M0
    K2 <- (M0 * M2 - M1^2) / M0^2
    cumul <- rbind(cumul, c(t, K0, K1, K2))
    if (c %% 1000 == 0) .message("CGF progress: %d/%d", c, length.out)
  }

  re <- list(
    N = length(mresid),
    mresid = mresid,
    cumul = cumul,
    tX = tX,
    yVec = yVec,
    X.invXX = X.invXX,
    subjData = subjData
  )

  class(re) <- "SPACox_NULL_Model"
  return(re)
}


#' @title Validate control parameters for SPACox marker testing
#' @param control List of control parameters for marker-level analysis.
#' @return Updated control list with validated parameters and defaults.
#' @keywords internal
checkControl.Marker.SPACox <- function(control) {
  default.control <- list(
    pVal_covaAdj_Cutoff = 5e-05,
    SPA_Cutoff = 2
  )

  control <- updateControl(control, default.control)

  return(control)
}


#' @title Set up marker-level testing for SPACox method
#' @param objNull Null model object from GRAB.NullModel().
#' @param control List of control parameters.
#' @return List containing setup parameters for marker testing.
#' @keywords internal
setMarker.SPACox <- function(
  objNull,
  control
) {
  cumul <- objNull$cumul
  mresid <- objNull$mresid
  XinvXX <- objNull$X.invXX
  tX <- objNull$tX
  N <- length(mresid)
  pVal_covaAdj_Cutoff <- control$pVal_covaAdj_Cutoff
  SPA_Cutoff <- control$SPA_Cutoff

  setSPACoxobjInCPP(
    cumul,
    mresid,
    XinvXX,
    tX,
    N,
    pVal_covaAdj_Cutoff,
    SPA_Cutoff
  )
}


#' @title Perform marker-level analysis for SPACox method
#' @param genoType Character string specifying genotype file format.
#' @param genoIndex Integer vector of genotype indices to analyze.
#' @param outputColumns Character vector specifying output columns to include.
#' @return Data frame containing analysis results.
#' @keywords internal
mainMarker.SPACox <- function(
  genoType,
  genoIndex,
  outputColumns
) {
  OutList <- mainMarkerInCPP("SPACox", genoType, genoIndex)
  obj.mainMarker <- data.frame(
    Marker = OutList$markerVec, # marker IDs
    Info = OutList$infoVec, # marker information: CHR:POS:REF:ALT
    AltFreq = OutList$altFreqVec, # alternative allele frequencies
    AltCounts = OutList$altCountsVec, # alternative allele counts
    MissingRate = OutList$missingRateVec, # alternative allele counts
    Pvalue = OutList$pvalVec # marker-level p-values
  ) 

  optionalColumns <- c("zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  additionalColumns <- intersect(optionalColumns, outputColumns)

  if (length(additionalColumns) > 0) {
    obj.mainMarker <- cbind.data.frame(
      obj.mainMarker,
      as.data.frame(OutList[additionalColumns])
    )
  }

  return(obj.mainMarker)
}
