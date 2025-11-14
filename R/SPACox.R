## ------------------------------------------------------------------------------
## SPACox.R
## Saddlepoint approximation–based Cox model for time-to-event data (unrelated).
## Provides null model fitting and single-variant testing utilities.
##
## Functions:
##   GRAB.SPACox                  : Print brief method information.
##   checkControl.NullModel.SPACox: Validate and populate null-model controls.
##   fitNullModel.SPACox          : Fit SPACox null model.
##   checkControl.Marker.SPACox   : Validate marker-level controls.
##   setMarker.SPACox             : Initialize marker-level analysis objects.
##   mainMarker.SPACox            : Run marker-level SPACox tests.
## ------------------------------------------------------------------------------

#' SPACox method in GRAB package
#'
#' SPACox method is an empirical approach to analyzing complex traits
#' (including but not limited to time-to-event trait) for unrelated samples
#' in a large-scale biobank.
#'
#' @details
#' \strong{Additional Control Parameters for GRAB.NullModel()}:
#' \itemize{
#'   \item \code{range} (numeric vector, default: c(-100, 100)): Range for saddlepoint approximation grid. Must be symmetric (range\[2\] = -range\[1\]).
#'   \item \code{length.out} (integer, default: 10000): Number of grid points for saddlepoint approximation.
#' }
#'
#' \strong{Additional Control Parameters for GRAB.Marker()}:
#' \itemize{
#'   \item \code{pVal_covaAdj_Cutoff} (numeric, default: 5e-05): P-value cutoff for covariate adjustment.
#'   \item \code{SPA_Cutoff} (numeric, default: 2): Cutoff for saddlepoint approximation.
#' }
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
#'   subjIDcol = "IID",
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
#'   subjIDcol = "IID",
#'   method = "SPACox",
#'   traitType = "Residual"
#' )
#'
#' # Step 2: conduct score test
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputDir <- tempdir()
#' OutputFile <- file.path(OutputDir, "resultSPACox.txt")
#' GRAB.Marker(obj.SPACox, GenoFile, OutputFile)
#' data.table::fread(OutputFile)
#'
GRAB.SPACox <- function() {
  .message("Using SPACox method - see ?GRAB.SPACox for details")
}


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


#' Fit SPACox null model from survival outcomes or residuals
#'
#' Computes martingale residuals (or uses provided residuals) and an empirical
#' cumulant generating function (CGF) for SPA-based single-variant tests.
#'
#' @param response Either a \code{survival::Surv} object (time-to-event) or a
#'   numeric residual vector with class \code{"Residual"}.
#' @param designMat Numeric design matrix (n x p) of covariates.
#' @param subjIDcol Vector of subject IDs aligned with rows of \code{designMat}.
#' @param control List with fields such as \code{range} and \code{length.out}
#'   for the CGF grid.
#' @param ... Extra arguments passed to \code{survival::coxph} when
#'   \code{response} is \code{Surv}.
#'
#' @return A list of class \code{"SPACox_NULL_Model"} with elements:
#'   \describe{
#'     \item{N}{Number of subjects.}
#'     \item{mresid}{Martingale residuals (numeric vector).}
#'     \item{cumul}{CGF grid as a matrix with columns t, K0, K1, K2.}
#'     \item{tX}{Transpose of design matrix with intercept (p+1 x n).}
#'     \item{yVec}{Status/event indicator or residual-based response.}
#'     \item{X.invXX}{Projection helper: X %*% solve(t(X) %*% X).}
#'     \item{subjData}{Character vector of subject IDs.}
#'   }
#'
#' @keywords internal
fitNullModel.SPACox <- function(
  response,
  designMat,
  subjIDcol,
  control,
  ...
) {
  subjData <- subjIDcol  # Use subjData internally for compatibility
  
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


checkControl.Marker.SPACox <- function(control) {
  default.control <- list(
    pVal_covaAdj_Cutoff = 5e-05,
    SPA_Cutoff = 2
  )

  control <- updateControl(control, default.control)

  # Validate parameters
  if (!is.numeric(control$pVal_covaAdj_Cutoff) || control$pVal_covaAdj_Cutoff <= 0) {
    stop("control$pVal_covaAdj_Cutoff should be a numeric value > 0.")
  }

  if (!is.numeric(control$SPA_Cutoff) || control$SPA_Cutoff <= 0) {
    stop("control$SPA_Cutoff should be a numeric value > 0.")
  }

  return(control)
}


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


mainMarker.SPACox <- function(
  genoType,
  genoIndex
) {
  OutList <- mainMarkerInCPP("SPACox", genoType, genoIndex)
  
  obj.mainMarker <- data.frame(
    Marker = OutList$markerVec, # marker IDs
    Info = OutList$infoVec, # marker information: CHR:POS:REF:ALT
    AltFreq = OutList$altFreqVec, # alternative allele frequencies
    AltCounts = OutList$altCountsVec, # alternative allele counts
    MissingRate = OutList$missingRateVec, # alternative allele counts
    Pvalue = OutList$pvalVec, # marker-level p-values
    zScore = OutList$zScore
  ) 

  return(obj.mainMarker)
}
