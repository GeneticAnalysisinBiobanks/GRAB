## ------------------------------------------------------------------------------
## SPAmix.R
## SPAmix: saddlepoint approximation for mixture populations. Provides null
## model fitting and single-variant testing with ancestry-mixing adjustments.
##
## Functions:
##   GRAB.SPAmix                  : Print brief method information.
##   checkControl.NullModel.SPAmix: Validate and populate null-model controls.
##   fitNullModel.SPAmix          : Fit SPAmix null model.
##   checkControl.Marker.SPAmix   : Validate marker-level controls.
##   setMarker.SPAmix             : Initialize marker-level analysis objects.
##   mainMarker.SPAmix            : Run marker-level SPAmix tests.
## ------------------------------------------------------------------------------

#' SPAmix method in GRAB package
#'
#' SPAmix method is an empirical approach to analyzing complex traits
#' (including but not limited to time-to-event trait) for unrelated samples
#' in a large-scale biobank. SPAmix extends SPACox to support an admixture
#' population or multiple populations.
#'
#' @return NULL
#'
#' @examples
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultSPAmix.txt")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#'
#' # Step 1 option 1
#' obj.SPAmix <- GRAB.NullModel(
#'   survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   method = "SPAmix",
#'   traitType = "time-to-event",
#'   control = list(PC_columns = "PC1,PC2")
#' )
#'
#' # Step 1 option 2
#' residuals <- survival::coxph(
#'   survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
#'   data = PhenoData
#' )$residuals
#'
#' obj.SPAmix <- GRAB.NullModel(
#'   residuals ~ AGE + GENDER + PC1 + PC2,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   method = "SPAmix",
#'   traitType = "Residual",
#'   control = list(PC_columns = "PC1,PC2")
#' )
#'
#' # Step 1 option 2: analyze multiple traits at once
#' res_cox <- survival::coxph(
#'   survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
#'   data = PhenoData
#' )$residuals
#'
#' res_lm <- lm(QuantPheno ~ AGE + GENDER + PC1 + PC2, data = PhenoData)$residuals
#'
#' obj.SPAmix <- GRAB.NullModel(
#'   res_cox + res_lm ~ AGE + GENDER + PC1 + PC2,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   method = "SPAmix",
#'   traitType = "Residual",
#'   control = list(PC_columns = "PC1,PC2")
#' )
#'
#' # Step 2
#' GRAB.Marker(obj.SPAmix, GenoFile, OutputFile)
#'
#' head(data.table::fread(OutputFile))
#'
#' @details
#'
#' \strong{Additional Control Parameters for GRAB.NullModel()}:
#' \itemize{
#'   \item \code{PC_columns} (character, required): Comma-separated column names
#'      of principal components (e.g., "PC1,PC2").
#'   \item \code{OutlierRatio} (numeric, default: 1.5): IQR multiplier for outlier detection.
#'      Outliers are defined as values outside \[Q1 - r*IQR, Q3 + r*IQR\].
#' }
#'
#' \strong{Additional Control Parameters for GRAB.Marker()}:
#' \itemize{
#'   \item \code{dosage_option} (character, default: "rounding_first"):
#'    Dosage handling option. Must be either "rounding_first" or "rounding_last".
#' }
#'
#' \strong{Output file columns}:
#' \describe{
#'   \item{Pheno}{Phenotype identifier (for multi-trait analysis).}
#'   \item{Marker}{Marker identifier (rsID or CHR:POS:REF:ALT).}
#'   \item{Info}{Marker information in format CHR:POS:REF:ALT.}
#'   \item{AltFreq}{Alternative allele frequency in the sample.}
#'   \item{AltCounts}{Total count of alternative alleles.}
#'   \item{MissingRate}{Proportion of missing genotypes.}
#'   \item{Pvalue}{P-value from the score test.}
#'   \item{zScore}{Z-score from the score test.}
#' }
#'
GRAB.SPAmix <- function() {
  .message("?GRAB.SPAmix for instructions")
}


checkControl.NullModel.SPAmix <- function(traitType, GenoFile, SparseGRMFile, control) {

  if (!traitType %in% c("time-to-event", "Residual")) {
    stop("For 'SPAmix' method, only traitType of 'time-to-event' or 'Residual' is supported.")
  }

  if (!is.null(GenoFile)) {
    warning("Argument 'GenoFile' is ignored for method 'SPAmix'.")
  }
  
  if (!is.null(SparseGRMFile)) {
    warning("Argument 'SparseGRMFile' is ignored for method 'SPAmix'.")
  }

  default.control <- list(
    OutlierRatio = 1.5
  )
  control <- updateControl(control, default.control)

  if (is.null(control$PC_columns)) {
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') is required for 'SPAmix' method.")
  }

  if (length(control$PC_columns) != 1) {
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') should be a character, not a character vector.")
  }

  control$PC_columns <- unlist(strsplit(control$PC_columns, split = ","))
  if (length(control$PC_columns) == 1) {
    warning("We detected that only one PC column exists, is that what you want? ",
            "Note that control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') should be ",
            "a character string split using ','.")
  }

  return(list(control = control, optionGRM = NULL))
}


#' Fit a SPAmix null model from a survival response (\code{Surv}) with
#' covariates or from precomputed residuals. Principal components (PCs) named
#' in \code{control$PC_columns} are extracted; residual outliers are detected
#' using an IQR rule with adjustable multiplier and stored for SPA testing.
#'
#' @param response Either a \code{survival::Surv} object (time-to-event data)
#'   or a numeric residual vector/matrix with class \code{"Residual"}.
#' @param designMat Numeric matrix (n x p) of covariates; must include the PC
#'   columns specified in \code{control$PC_columns}.
#' @param subjData Vector of subject IDs aligned with rows of \code{designMat}
#'   and \code{response}.
#' @param control List of options. Required element: \code{PC_columns}, a
#'   single comma-separated string of PC column names (e.g.
#'   \code{"PC1,PC2,PC3,PC4"}). \code{OutlierRatio}.
#' @param ... Extra arguments passed to \code{survival::coxph} when
#'   \code{response} is \code{Surv}.
#'
#' @details If \code{response} is \code{Surv}, a Cox model is fit and martingale
#' residuals are used. If \code{response} is \code{Residual}, its values are
#' used directly. Outliers per phenotype are defined by
#' \code{[Q1 - r*IQR, Q3 + r*IQR]} with \code{r = OutlierRatio}; if none are
#' found, \code{r} is iteratively reduced by 20% until at least one appears.
#'
#' @return A list of class \code{"SPAmix_NULL_Model"} containing:
#'   \describe{
#'     \item{resid}{Residual matrix (n x k)}
#'     \item{N}{Number of subjects}
#'     \item{yVec}{Response vector (event indicator for survival models)}
#'     \item{PCs}{Selected principal component columns}
#'     \item{nPheno}{Number of phenotypes (columns of residuals)}
#'     \item{outLierList}{List of per-phenotype indices (0-based) and residual
#'       subsets for outlier/non-outlier strata}
#'   }
#'
#' @keywords internal
fitNullModel.SPAmix <- function(
  response,
  designMat,
  subjData,
  control,
  ...
) {

  if (!(inherits(response, "Surv") || inherits(response, "Residual"))) {
    stop("For SPAmix, the response variable should be of class 'Surv' or 'Residual'.")
  }

  if (inherits(response, "Surv")) {
    formula <- response ~ designMat

    obj.coxph <- survival::coxph(formula, x = TRUE, ...)

    ### Check input arguments

    y <- obj.coxph$y
    yVec <- y[, ncol(y)]

    mresid <- obj.coxph$residuals
    Cova <- obj.coxph$x

    if (length(mresid) != length(subjData)) {
      stop("Please check the consistency between 'formula' and 'subjData'.")
    }

    mresid <- matrix(mresid, ncol = 1)
  }

  if (inherits(response, "Residual")) {
    yVec <- mresid <- response
    Cova <- designMat

    .message("Design matrix: %d samples x %d covariates", nrow(designMat), ncol(designMat))
    if (nrow(mresid) != length(subjData)) {
      stop("Please check the consistency between 'formula' and 'subjData'.")
    }
  }

  PC_columns <- control$PC_columns

  # Debug info for covariate structure
  if (length(PC_columns) > 0) {
    .message("PC columns: %s", paste(PC_columns, collapse = ", "))
  }

  if (any(!PC_columns %in% colnames(designMat))) {
    stop("PC columns specified in 'control$PC_columns' should be in 'formula'.")
  }

  pos_col <- match(PC_columns, colnames(designMat))

  PCs <- Cova[, pos_col, drop = FALSE]
  outLierList <- list()
  nPheno <- ncol(mresid)
  for (i in 1:nPheno) {
    mresid.temp <- mresid[, i]

    ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
    q25 <- quantile(mresid.temp, 0.25, na.rm = TRUE)
    q75 <- quantile(mresid.temp, 0.75, na.rm = TRUE)
    IQR <- q75 - q25
    # outlier ratio with fallback if not specified in control
    r.outlier <- control$OutlierRatio
    # put this to the control argument later
    cutoff <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier <- which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])

    while (length(posOutlier) == 0) {
      r.outlier <- r.outlier * 0.8
      cutoff <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
      posOutlier <- which(mresid < cutoff[1] | mresid > cutoff[2])
      .message("Outlier ratio adjusted to: %.2f (%d outliers)", r.outlier, length(posOutlier))
    }


    posValue <- which(!is.na(mresid.temp))
    posNonOutlier <- setdiff(posValue, posOutlier)

    .message("Outlier cutoffs: [%.2f, %.2f]", cutoff[1], cutoff[2])
    .message("Outliers for SPA analysis: %d/%d (%.1f%%)",
             length(posOutlier), length(posValue),
             100 * length(posOutlier) / length(posValue))

    if (length(posOutlier) == 0) {
      stop("No outlier is observed. SPA is not required in this case.")
    }

    # "-1" is to convert R style (index starting from 1) to C++ style (index starting from 0)
    outLierList[[i]] <- list(
      posValue = posValue - 1,
      posOutlier = posOutlier - 1,
      posNonOutlier = posNonOutlier - 1,
      resid = mresid.temp[posValue],
      resid2 = mresid.temp[posValue]^2,
      residOutlier = mresid.temp[posOutlier],
      residNonOutlier = mresid.temp[posNonOutlier],
      resid2NonOutlier = mresid.temp[posNonOutlier]^2
    )
  }

  objNull <- list(
    resid = mresid,
    N = nrow(Cova),
    yVec = yVec, # event variable: 0 or 1
    PCs = PCs,
    nPheno = nPheno,
    outLierList = outLierList
  )

  class(objNull) <- "SPAmix_NULL_Model"
  return(objNull)
}


checkControl.Marker.SPAmix <- function(control) {

  default.control <- list(
    dosage_option = "rounding_first"
  )
  control <- updateControl(control, default.control)

  if (!control$dosage_option %in% c("rounding_first", "rounding_last")) {
    stop("control$dosage_option should be 'rounding_first' or 'rounding_last'.")
  }

  return(control)
}


setMarker.SPAmix <- function(objNull, control) {

  setSPAmixobjInCPP(
    t_resid = objNull$resid,              # matrix: Residuals from null model
    t_PCs = objNull$PCs,                  # matrix: Principal components for population structure
    t_N = objNull$N,                      # integer: Sample size
    t_SPA_Cutoff = control$SPA_Cutoff,    # numeric: P-value cutoff for SPA correction
    t_outlierList = objNull$outLierList   # list: Outlier subject information
  )
}


mainMarker.SPAmix <- function(
  genoType,
  genoIndex,
  objNull
) {
  
  OutList <- mainMarkerInCPP(
    t_method = "SPAmix",      # character: Statistical method name
    t_genoType = genoType,    # character: "PLINK" or "BGEN"
    t_genoIndex = genoIndex   # integer vector: Genotype indices to analyze
  )

  nPheno <- objNull$nPheno
  obj.mainMarker <- data.frame(
    Pheno = paste0("pheno_", 1:nPheno),
    Marker = rep(OutList$markerVec, each = nPheno), # marker IDs
    Info = rep(OutList$infoVec, each = nPheno), # marker information: CHR:POS:REF:ALT
    AltFreq = rep(OutList$altFreqVec, each = nPheno), # alternative allele frequencies
    AltCounts = rep(OutList$altCountsVec, each = nPheno), # alternative allele counts
    MissingRate = rep(OutList$missingRateVec, each = nPheno), # alternative allele counts
    Pvalue = OutList$pvalVec, # marker-level p-values
    zScore = OutList$zScore
  )

  return(obj.mainMarker)
}
