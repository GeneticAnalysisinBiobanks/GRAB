## ------------------------------------------------------------------------------
## WtCoxG.R
## Weighted Cox regression with case-ascertainment bias correction (two-step):
##   1) Null model fitting + batch effect assessment using reference AFs
##   2) Genome-wide association testing with external reference integration
##
## Functions:
##   GRAB.WtCoxG                  : Print brief method information.
##   checkControl.NullModel.WtCoxG: Validate/populate null-model controls.
##   fitNullModel.WtCoxG          : Fit weighted Cox null model & test batch effects.
##   setMarker.WtCoxG             : Initialize marker-level analysis objects.
##   mainMarker.WtCoxG            : Marker-level association testing.
##   TestforBatchEffect           : QC + parameter estimation (TPR, sigma2, weights).
## ------------------------------------------------------------------------------

#' Weighted Cox regression for genetic association analysis with case ascertainment
#'
#' WtCoxG provides an accurate, powerful, and computationally efficient Cox-based
#' approach for genome-wide time-to-event analyses in study cohorts with case
#' ascertainment.
#'
#' @details
#' \strong{Two-Step Analysis Process:}
#' 1. **Step 1**: Fit null model and test for batch effects using sample SNPs
#' 2. **Step 2**: Conduct GWAS with allele frequencies of a reference sample
#'
#' \strong{Additional Parameters for \code{GRAB.NullModel()}:}
#' \itemize{
#'   \item \code{RefAfFile} (character, required): Reference allele frequency file path. 
#'     File must contain columns: CHROM, POS, ID, REF, ALT, AF_ref, AN_ref
#'   \item \code{SurvTimeColumn} (character, default: "SurvTime"): Column name in 
#'     \code{data} containing survival times
#'   \item \code{IndicatorColumn} (character, default: "Indicator"): Column name in 
#'     \code{data} for case-control status (0 = control, 1 = case)
#' }
#'
#' \strong{Additional Control Parameters for GRAB.NullModel()}:
#' \itemize{
#'   \item \code{RefPrevalence} (numeric, required): Population-level disease prevalence 
#'     for weighting. Must be in range (0, 0.5)
#'   \item \code{OutlierRatio} (numeric, default: 1.5): IQR multiplier for outlier detection
#'   \item \code{SNPnum} (numeric, default: 1e4): Minimum number of SNPs for batch effect testing
#' }
#'
#' \strong{Additional Control Parameters for GRAB.Marker()}:
#' \itemize{
#'   \item \code{cutoff} (numeric, default: 0.1): Cutoff of batch effect test p-value for 
#'     association testing. Variants with batch effect p-value below this cutoff
#'     will be excluded from association testing.
#' }
#'
#' @return No return value. Called for informational side effects.
#'
#'
#' @examples
#' # Step0&1: fit a null model and estimate parameters according to batch effect p-values
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#'
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' RefAfFile <- system.file("extdata", "simuRefAf.txt", package = "GRAB")
#' RefPrevalence <- 0.1 # population-level disease prevalence
#'
#' OutputDir <- tempdir()
#' OutputFile <- file.path(OutputDir, "WtCoxG_step2_out.txt")
#'
#' obj.WtCoxG <- GRAB.NullModel(
#'   formula = survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   method = "WtCoxG",
#'   traitType = "time-to-event",
#'   GenoFile = GenoFile,
#'   SparseGRMFile = SparseGRMFile,
#'   control = list(
#'     AlleleOrder = "ref-first",
#'     AllMarkers = TRUE,
#'     RefPrevalence = RefPrevalence,
#'     SNPnum = 1000
#'   ),
#'   RefAfFile = RefAfFile,
#'   SurvTimeColumn = "SurvTime",
#'   IndicatorColumn = "SurvEvent"
#' )
#'
#' obj.WtCoxG$mergeGenoInfo[, c("CHROM", "POS", "pvalue_bat")]
#'
#' # Step2: conduct association testing
#' OutputFile <- file.path(tempdir(), "resultWtCoxG.txt")
#' GRAB.Marker(obj.WtCoxG, GenoFile, OutputFile,
#'   control = list(
#'     AlleleOrder = "ref-first",
#'     AllMarkers = TRUE,
#'     cutoff = 0.1,
#'     nMarkersEachChunk = 5000
#'   )
#' )
#'
#' data.table::fread(OutputFile)[, c("CHROM", "POS", "WtCoxG.noext", "WtCoxG.ext")]
#'
GRAB.WtCoxG <- function() {
  .message("Using WtCoxG method - see ?GRAB.WtCoxG for details")
}


checkControl.NullModel.WtCoxG <- function(control, traitType) {
  # Ensure only time-to-event traits are supported
  if (!traitType %in% c("time-to-event")) {
    stop("For 'WtCoxG' method, only traitType of 'time-to-event' is supported.")
  }

  # Set default control parameters
  default.control <- list(RefPrevalence = 0, OutlierRatio = 1.5, SNPnum = 1e4)
  control <- updateControl(control, default.control)

  # Validate reference prevalence parameter
  if (control$RefPrevalence <= 0 || control$RefPrevalence >= 0.5) {
    stop("control$RefPrevalence is required and should be between (0, 0.5).")
  }

  return(control)
}


#' Fit weighted Cox null model with outlier handling and batch-effect QC
#'
#' Fits a weighted Cox model using case/control weighting based on reference
#' prevalence and identifies residual outliers for SPA testing. Optionally
#' performs batch-effect QC by cross-referencing external allele frequencies.
#'
#' @param response \code{survival::Surv} response (time-to-event). Residuals are
#'   not supported here.
#' @param designMat Numeric matrix (n x p) of covariates.
#' @param subjData Character vector of subject IDs aligned with rows of \code{designMat}.
#' @param control List with fields such as \code{RefPrevalence} (0, 0.5),
#'   \code{OutlierRatio}, \code{SNPnum}.
#' @param data Data frame used for optional batch-effect QC.
#' @param GenoFile Character. PLINK prefix (without extension) used when
#'   sampling markers for QC.
#' @param GenoFileIndex Character. Path to an index file used in QC workflow.
#' @param SparseGRMFile Character. Path to sparse GRM used in QC workflow.
#' @param ... Optional named parameters forwarded to QC (e.g.,
#'   \code{RefAfFile}, \code{IndicatorColumn},
#'   \code{SurvTimeColumn}).
#'
#' @return A list of class \code{"WtCoxG_NULL_Model"} with elements:
#'   \describe{
#'     \item{mresid}{Martingale residuals from weighted Cox model.}
#'     \item{Cova}{Design matrix used in the null model (n x p).}
#'     \item{yVec}{Event indicator extracted from \code{Surv}.}
#'     \item{weight}{Observation weights derived from reference prevalence.}
#'     \item{RefPrevalence}{Reference prevalence used to define weights.}
#'     \item{N}{Number of subjects.}
#'     \item{outLierList}{Lists indices (0-based) and residual subsets for SPA.}
#'     \item{control}{Copy of control options used.}
#'     \item{mergeGenoInfo}{QC-derived marker metadata for batch-effect testing (if run).}
#'     \item{subjData}{Character vector of subject IDs.}
#'   }
#'
#' @keywords internal
fitNullModel.WtCoxG <- function(
  response, designMat, subjData, control, data,
  GenoFile, GenoFileIndex, SparseGRMFile, ...
) {

  # ========== Validate additional parameters from dots ==========
  
  dots <- list(...)
  
  # Validate RefAfFile (required)
  if (is.null(dots$RefAfFile)) {
    stop("Argument 'RefAfFile' is required for WtCoxG method.")
  }
  if (!is.character(dots$RefAfFile) || length(dots$RefAfFile) != 1) {
    stop("Argument 'RefAfFile' should be a character string (file path).")
  }
  if (!file.exists(dots$RefAfFile)) {
    stop("Cannot find RefAfFile: ", dots$RefAfFile)
  }
  
  # Set defaults and validate optional column name parameters
  SurvTimeColumn <- if (is.null(dots$SurvTimeColumn)) "SurvTime" else dots$SurvTimeColumn
  if (!is.character(SurvTimeColumn) || length(SurvTimeColumn) != 1) {
    stop("Argument 'SurvTimeColumn' should be a character string (column name).")
  }
  
  IndicatorColumn <- if (is.null(dots$IndicatorColumn)) "Indicator" else dots$IndicatorColumn
  if (!is.character(IndicatorColumn) || length(IndicatorColumn) != 1) {
    stop("Argument 'IndicatorColumn' should be a character string (column name).")
  }
  
  # Validate that columns exist in data
  if (!SurvTimeColumn %in% colnames(data)) {
    stop("Column '", SurvTimeColumn, "' not found in data.")
  }
  if (!IndicatorColumn %in% colnames(data)) {
    stop("Column '", IndicatorColumn, "' not found in data.")
  }

  # ========== Validate response type ==========

  if (!(inherits(response, "Surv") || inherits(response, "Residual"))) {
    stop("For WtCoxG, the response variable should be of class 'Surv' or 'Residual'.")
  }

  if (inherits(response, "Surv")) {
    formula <- response ~ designMat

    Indicator <- as.matrix(response)[, "status"]
    RefPrevalence <- control$RefPrevalence

    # ---- BEGIN inlined: getWeight.WtCoxG ----
    if (any(!unique(Indicator) %in% c(0, 1))) {
      stop("The value of Indicator should be 0 or 1.")
    }

    sumOnes <- sum(Indicator)
    sumZeros <- sum(1 - Indicator)
    ratio <- sumOnes / sumZeros

    weight <- ifelse(
      Indicator == 1,
      1,
      (1 - RefPrevalence) / RefPrevalence * ratio
    )
    # ---- END inlined: getWeight.WtCoxG ----

    obj.coxph <- survival::coxph(formula, x = TRUE, weight = weight, robust = TRUE)

    y <- obj.coxph$y
    yVec <- y[, ncol(y)] # status

    mresid <- obj.coxph$residuals
    Cova <- designMat
  } else {
    stop("We only support 'time-to-event' trait for WtCoxG by 2023-08-08.")
  }

  ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
  q25 <- quantile(mresid, 0.25, na.rm = TRUE)
  q75 <- quantile(mresid, 0.75, na.rm = TRUE)
  IQR <- q75 - q25

  r.outlier <- ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)
  cutoff <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
  posOutlier <- which(mresid < cutoff[1] | mresid > cutoff[2])

  while (length(posOutlier) == 0) {
    r.outlier <- r.outlier * 0.8
    cutoff <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier <- which(mresid < cutoff[1] | mresid > cutoff[2])
    .message("Adjusted outlier ratio: %.2f (%d outliers)", r.outlier, length(posOutlier))
  }

  # The original code is from SPAmix in which multiple residuals were analysis simultaneously
  posValue <- seq_along(mresid)
  posNonOutlier <- setdiff(posValue, posOutlier)

  .message("Outliers detected: %d/%d (%.1f%%)", length(posOutlier), length(posValue),
           100 * length(posOutlier) / length(posValue))

  if (length(posOutlier) == 0) {
    stop("No outlier is observed. SPA is not required in this case.")
  }

  # "-1" is to convert R style (index starting from 1) to C++ style (index starting from 0)
  outLierList <- list(
    posOutlier = posOutlier - 1,
    posNonOutlier = posNonOutlier - 1,
    resid = mresid,
    resid2 = mresid^2,
    residOutlier = mresid[posOutlier],
    residNonOutlier = mresid[posNonOutlier],
    resid2NonOutlier = mresid[posNonOutlier]^2
  )

  re <- list(
    mresid = mresid, Cova = Cova, yVec = yVec, weight = weight,
    RefPrevalence = RefPrevalence,
    N = length(mresid),
    outLierList = outLierList,
    control = control,
    subjData = subjData
  )

  class(re) <- "WtCoxG_NULL_Model"

  # Test for batch effect using validated parameters
  re$mergeGenoInfo <- TestforBatchEffect(
    objNull = re,
    data = data,
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    SparseGRMFile = SparseGRMFile,
    RefAfFile = dots$RefAfFile,
    IndicatorColumn = IndicatorColumn,
    SurvTimeColumn = SurvTimeColumn
  )

  return(re)
}


checkControl.Marker.WtCoxG <- function(control) {
  default.control <- list(
    cutoff = 0.1
  )

  control <- updateControl(control, default.control)

  # Validate parameters
  if (!is.numeric(control$cutoff) || control$cutoff < 0 || control$cutoff > 1) {
    stop("control$cutoff should be a numeric value in [0, 1].")
  }

  return(control)
}


setMarker.WtCoxG <- function(objNull, control) {
  ImputeMethod <- if (is.null(control$ImputeMethod)) "none" else control$ImputeMethod
  cutoff <- if (is.null(control$cutoff)) 0.1 else control$cutoff
  setWtCoxGobjInCPP(
    objNull$mresid,
    objNull$weight,
    ImputeMethod,
    cutoff
  )
}


mainMarker.WtCoxG <- function(genoType, genoIndex, control, objNull) {
  mergeGenoInfo <- objNull$mergeGenoInfo
  mergeGenoInfo_subset <- mergeGenoInfo[mergeGenoInfo$genoIndex %in% genoIndex, ]

  # Use match to reorder, but check for missing values
  match_indices <- match(genoIndex, mergeGenoInfo_subset$genoIndex)
  if (any(is.na(match_indices))) {
    missing_indices <- genoIndex[is.na(match_indices)]
    stop(paste0(
      "Missing marker info for genoIndex values: ",
      paste(missing_indices[seq_len(min(10, length(missing_indices)))], collapse = ", "),
      if (length(missing_indices) > 10) " ..." else ""
    ))
  }

  mergeGenoInfo_subset <- mergeGenoInfo_subset[match_indices, ]

  # Safety check: ensure sizes match
  if (nrow(mergeGenoInfo_subset) != length(genoIndex)) {
    stop(paste0(
      "Size mismatch: genoIndex has ", length(genoIndex),
      " elements, but mergeGenoInfo_subset has ", nrow(mergeGenoInfo_subset),
      " rows. Some markers in genoIndex may not be found in mergeGenoInfo."
    ))
  }

  # Update WtCoxG object with marker information for current chunk
  updateWtCoxGChunkInCPP(mergeGenoInfo_subset)

  OutList <- mainMarkerInCPP("WtCoxG", genoType, genoIndex)
  pvals <- data.frame(matrix(OutList$pvalVec, ncol = 2, byrow = TRUE))
  colnames(pvals) <- c("WtCoxG.ext", "WtCoxG.noext")

  obj.mainMarker <- cbind(pvals, mergeGenoInfo_subset)
  return(obj.mainMarker)
}


#' Quality control to check batch effect between study cohort and reference population.
#'
#' This function performs quality control to test for the batch effect between a study
#' cohort and a reference population. And fit a weighted null model.
#'
#' @param objNull a \code{WtCoxG_NULL_Model} object, which is the output of \code{\link{GRAB.NullModel}}.
#' @param data a data.frame, list or environment (or object coercible by \code{\link{as.data.frame}}
#'   to a data.frame), containing the variables in formula. Neither a matrix nor an array will be accepted.
#' @param GenoFile A character string of the genotype file. See Details section for more details.
#' @param GenoFileIndex Additional index file(s) corresponding to GenoFile. See Details section for more details.
#' @param Geno.mtx A matrix of genotype data. If provided, it will be used instead of GenoFile.
#'   The matrix should have samples in rows and markers in columns.
#' @param SparseGRMFile a path to file of output to be passed to \code{\link{GRAB.NullModel}}.
#' @param RefAfFile A character string of the reference file. The reference file must be a \code{txt}
#'   file (header required) including at least 7 columns: \code{CHROM}, \code{POS}, \code{ID},
#'   \code{REF}, \code{ALT}, \code{AF_ref}, \code{AN_ref}.
#' @param IndicatorColumn A character string of the column name in \code{data} that indicates
#'   the case-control status. The value should be 0 for controls and 1 for cases.
#' @param SurvTimeColumn A character string of the column name in \code{data} that indicates the survival time.
#' @return A dataframe of marker info and reference MAF.
#'
#' @keywords internal
TestforBatchEffect <- function(
  objNull,
  data,
  GenoFile = NULL, # a character of file names of genotype files
  GenoFileIndex = NULL, # additional index file(s) corresponding to GenoFile
  Geno.mtx = NULL, # genotype matrix, if provided, will be used instead of GenoFile
  SparseGRMFile = NULL, # sparse genotype relatedness matrix
  RefAfFile, # header should include c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref","AN_ref")
  IndicatorColumn,
  SurvTimeColumn
) {
  .message("Testing for batch effect ...")

  if (!is.null(SparseGRMFile)) {
    sparseGRM <- data.table::fread(SparseGRMFile)
  } else {
    sparseGRM <- NULL
  }

  control <- objNull$control
  RefPrevalence <- control$RefPrevalence # refernce population prevalence, the proportion of indicator == 1.
  SNPnum <- control$SNPnum

  posCol <- which(colnames(data) == IndicatorColumn)
  colnames(data)[posCol] <- "Indicator"

  posCol <- which(colnames(data) == SurvTimeColumn)
  colnames(data)[posCol] <- "SurvTime"

  # Add SampleID column from objNull$subjData
  data$SampleID <- objNull$subjData


  # step1: quality control--------------------------------------------------------
  ## reference genoInfo----------------------------------------------------------
  refGenoInfo <- data.table::fread(RefAfFile) %>% as_tibble()

  # check if there are 7 columns in RefAfFile
  for (colname in c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref", "AN_ref")) {
    if (!colname %in% colnames(refGenoInfo)) {
      stop(paste0(colname, " is missing in RefAfFile!"))
    }
  }

  ## merge sample genoInfo and ref genoInfo--------------------------------------
  if (is.null(Geno.mtx)) {
    .message("Getting genotype info for controls (Indicator=0) ...")
    GenoInfo.ctrl <- GRAB.getGenoInfo(
      GenoFile = GenoFile,
      GenoFileIndex = GenoFileIndex,
      SampleIDs = with(data, SampleID[Indicator == 0]),
      control = control
    ) %>%
      rename(mu0 = altFreq, mr0 = missingRate) %>%
      select(mu0, mr0)

    if (nrow(GenoInfo.ctrl) < SNPnum) {
      stop("The number of genetic variants < ", SNPnum)
    }

    .message("Getting genotype info for cases (Indicator=1) ...")
    GenoInfo <- GRAB.getGenoInfo(
      GenoFile = GenoFile,
      GenoFileIndex = GenoFileIndex,
      SampleIDs = with(data, SampleID[Indicator == 1]),
      control = control
    ) %>%
      rename(mu1 = altFreq, mr1 = missingRate) %>%
      cbind(., GenoInfo.ctrl) %>%
      as_tibble() %>%
      mutate(RA = paste0(pmin(REF, ALT), pmax(REF, ALT))) %>%
      mutate(index = seq_len(n()))

    mergeGenoInfo <- refGenoInfo %>%
      mutate(RA = paste0(pmin(REF, ALT), pmax(REF, ALT))) %>%
      merge(., GenoInfo, by = c("CHROM", "POS", "RA"), all.y = TRUE, sort = FALSE) %>%
      rename(REF = REF.y, ALT = ALT.y, ID = ID.y) %>%
      mutate(AF_ref = ifelse(REF == REF.x, AF_ref, 1 - AF_ref)) %>%
      select(-REF.x, -ALT.x, -ID.x, -RA) %>%
      filter(!duplicated(index)) %>%
      mutate(
        n1 = sum(data$Indicator) * (1 - mr1),
        n0 = sum(1 - data$Indicator) * (1 - mr0),
        mu.int = 0.5 * mu1 + 0.5 * mu0,
        mu.int = ifelse(mu.int > 0.5, 1 - mu.int, mu.int)
      )
  } else {
    GenoInfo <- data.frame(
      ID = colnames(Geno.mtx),
      mu0 = apply(Geno.mtx[data$Indicator == 0, ], 2, function(x) {
        mean(na.omit(x) / 2)
      }),
      mu1 = apply(Geno.mtx[data$Indicator == 1, ], 2, function(x) {
        mean(na.omit(x) / 2)
      }),
      n0 = apply(Geno.mtx[data$Indicator == 0, ], 2, function(x) {
        sum(!is.na(x))
      }),
      n1 = apply(Geno.mtx[data$Indicator == 1, ], 2, function(x) {
        sum(!is.na(x))
      })
    ) %>%
      mutate(
        mu.int = 0.5 * mu1 + 0.5 * mu0,
        mu.int = ifelse(mu.int > 0.5, 1 - mu.int, mu.int),
        index = seq_len(n())
      )
    mergeGenoInfo <- merge(GenoInfo, refGenoInfo, by = "ID", all.x = TRUE, sort = FALSE) %>%
      filter(!duplicated(index))
  }

  #### calculate batch effect p-value for each genetic variant------------------------------------
  w1 <- objNull$weight / (2 * sum(objNull$weight))
  names(w1) <- data$SampleID
  meanR <- mean(objNull$mresid)
  R_tilde <- objNull$mresid - meanR
  names(R_tilde) <- data$SampleID

  if (!is.null(sparseGRM)) {
    sparseGRM <- sparseGRM %>% mutate(
      cov = Value * w1[as.character(ID1)] * w1[as.character(ID2)],
      cov_R = Value * R_tilde[as.character(ID1)] * R_tilde[as.character(ID2)]
    )
    var.ratio.w0 <- (sum(sparseGRM$cov) + 1 / (2 * mergeGenoInfo$AN_ref)) / 
      (sum(w1^2) + 1 / (2 * mergeGenoInfo$AN_ref))
    var.ratio.int <- sum(sparseGRM$cov_R) / sum(R_tilde^2)
  } else {
    var.ratio.w0 <- var.ratio.int <- 1
  }

  mergeGenoInfo <- mergeGenoInfo %>%
    mutate(
      var.ratio.w0 = var.ratio.w0,
      var.ratio.int = var.ratio.int
    )

  pvalue_bat <- lapply(seq_len(nrow(mergeGenoInfo)), function(ind) {
    # ---- BEGIN inlined: Batcheffect.TestOneMarker ----
    n0 <- mergeGenoInfo$n0[ind]
    n1 <- mergeGenoInfo$n1[ind]
    n.ext <- mergeGenoInfo$AN_ref[ind] / 2
    maf0 <- mergeGenoInfo$mu0[ind]
    maf1 <- mergeGenoInfo$mu1[ind]
    maf.ext <- mergeGenoInfo$AF_ref[ind]
    pop.prev <- RefPrevalence
    var.ratio <- mergeGenoInfo$var.ratio.w0[ind]
    
    er <- n1 / (n1 + n0)
    w0 <- (1 - pop.prev) / pop.prev / ((1 - er) / er)
    w1 <- 1

    ## weighted mean of genotypes
    weight.maf <- sum(maf0 * w0 * n0 + maf1 * w1 * n1) / sum(w0 * n0 + w1 * n1)
    ## MAF estimates
    est.maf <- sum(maf0 * w0 * n0 + maf1 * w1 * n1 + maf.ext * n.ext * w0) /
      sum(n1 * w1 + n0 * w0 + n.ext * w0)

    ## variance of test statistics
    v <- ((n1 * w1^2 + n0 * w0^2) / (2 * (n1 * w1 + n0 * w0)^2) + 1 / (2 * n.ext)) *
      est.maf * (1 - est.maf)
    z <- (weight.maf - maf.ext) / sqrt(v) ## standardized statistics
    z.adj <- z / sqrt(var.ratio) ## adjusted statistics by variance ratio
    p <- 2 * pnorm(-abs(z.adj), lower.tail = TRUE)
    # ---- END inlined: Batcheffect.TestOneMarker ----
    
    return(p)
  }) %>% unlist()

  mergeGenoInfo <- mergeGenoInfo %>% mutate(pvalue_bat)
  rm(pvalue_bat)

  #### estimate unknown parameters according to batch effect p-values---------------------------------
  .message("Estimating TPR and sigma2 ...")
  maf.group <- c(seq(-1e-4, 0.4, 0.05), max(mergeGenoInfo$mu.int))
  mergeGenoInfo <- lapply(1:(length(maf.group) - 1), function(i) {
    # Progress indicator for MAF group processing

    ## assume that genotypes with MAF in [ maf.group[i] , maf.group[i+1]] have the same mixture distribution
    mergeGenoInfo_1 <- mergeGenoInfo %>% filter(mu.int > maf.group[i] & mu.int <= maf.group[i + 1])

    ## using batcheffect p-values with MAF in [maf.group[i]-0.1 , maf.group[i+1]+0.1] to estimate parameters
    mergeGenoInfo_2 <- mergeGenoInfo %>%
      filter(mu.int >= max(maf.group[i] - 0.1, 0) & mu.int < min(1, maf.group[i + 1] + 0.1))

    mu <- (maf.group[i] + maf.group[i + 1]) / 2

    n.ext <- mean(na.omit(mergeGenoInfo_1$AN_ref)[1]) / 2
    var_mu_ext <- mu * (1 - mu) / (2 * n.ext)

    var_Sbat <- ifelse(is.null(sparseGRM), sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext,
      na.omit(mergeGenoInfo$var.ratio.w0)[1] * (sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext)
    )

    # ---- BEGIN inlined: fun.est.param ----
    vec_p_bat <- mergeGenoInfo_2$pvalue_bat
    vec_var_Sbat <- var_Sbat
    vec_cutoff <- seq(0.01, 0.4, 0.1)

    vec_p_deno <- sapply(vec_cutoff, function(p_cut) {
      mean(na.omit(vec_p_bat > p_cut))
    })

    opti_fun <- function(var_Sbat, vec_p_deno, par) {
      diff <- sum(sapply(seq_along(vec_cutoff), function(j) {
        p_cut <- vec_cutoff[j]
        lb <- -qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
        ub <- qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
        p_deno <- vec_p_deno[j]

        c <- pnorm(ub, 0, sqrt(var_Sbat + par[2]), log.p = TRUE)
        d <- pnorm(lb, 0, sqrt(var_Sbat + par[2]), log.p = TRUE)

        pro.cut <- par[1] * (exp(d) * (exp(c - d) - 1)) + (1 - par[1]) * (1 - p_cut)
        ((p_deno - pro.cut) / p_deno)^2
      }))
      diff
    }

    obj_par <- optim(
      par = c(0.01, 0.01),
      fn = opti_fun,
      vec_p_deno = vec_p_deno,
      var_Sbat = vec_var_Sbat
    )$par

    obj <- list(
      TPR = min(1, max(0, obj_par[1])),
      sigma2 = min(1, max(0, obj_par[2]))
    )

    TPR <- obj$TPR
    sigma2 <- obj$sigma2
    # ---- END inlined: fun.est.param ----

    # Function of optimal external weight contribution. Needed by the call of optim().
    fun.optimalWeight <- function(par, pop.prev, R, y, mu1, w, mu, N, n.ext, sigma2, TPR) {
      b <- par[1]

      p.fun <- function(b, pop.prev, R, y, mu1, mu, w, N, n.ext, sigma2, TPR) {
        meanR <- mean(R)
        sumR <- sum(R)

        mu0 <- mu
        mu.pop <- mu1 * pop.prev + mu0 * (1 - pop.prev)

        mu.i <- ifelse(y == 1, 2 * mu1, 2 * mu0)

        S <- sum((R - (1 - b) * meanR) * mu.i) - sumR * 2 * b * mu.pop

        w1 <- w / (2 * sum(w))
        mu <- mean(mu.i) / 2

        var_mu_ext <- mu * (1 - mu) / (2 * n.ext)
        var_Sbat <- sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext

        p_cut <- 0.1
        lb <- -qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
        ub <- qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
        c <- pnorm(ub, 0, sqrt(var_Sbat + sigma2), log.p = TRUE)
        d <- pnorm(lb, 0, sqrt(var_Sbat + sigma2), log.p = TRUE)
        p_deno <- TPR * (exp(d) * (exp(c - d) - 1)) + (1 - TPR) * (1 - p_cut)

        var.int <- sum((R - (1 - b) * meanR)^2) * 2 * mu * (1 - mu)
        var_S <- var.int + 4 * b^2 * sumR^2 * var_mu_ext
        cov_Sbat_S <- sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) + 2 * b * sumR * var_mu_ext
        VAR <- matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat), nrow = 2)
        p0 <- max(0, mvtnorm::pmvnorm(lower = c(-Inf, lb), upper = c(-abs(S), ub), mean = c(0, 0), sigma = VAR))

        var_S1 <- var.int + 4 * b^2 * sumR^2 * (var_mu_ext + sigma2)
        cov_Sbat_S1 <- sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) + 2 * b * sumR * (var_mu_ext + sigma2)
        var_Sbat1 <- var_Sbat + sigma2
        VAR1 <- matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1), nrow = 2)
        p1 <- max(0, mvtnorm::pmvnorm(lower = c(-Inf, lb), upper = c(-abs(S), ub), mean = c(0, 0), sigma = VAR1))

        p.con <- 2 * (TPR * p1 + (1 - TPR) * p0) / p_deno
        diff <- -log10(p.con / 5e-8)

        return(diff)
      }

      mu1 <- uniroot(p.fun,
        lower = mu, upper = 1,
        b = b, pop.prev = pop.prev, mu = mu,
        R = R, y = y, w = w, N = N, n.ext = n.ext, sigma2 = sigma2, TPR = TPR
      )$root

      return(mu1)
    }

    w.ext <- optim(
      par = 0.5, method = "L-BFGS-B", lower = 0, upper = 1,
      fn = fun.optimalWeight,
      pop.prev = RefPrevalence,
      y = data$Indicator,
      R = objNull$mresid,
      w = objNull$weight,
      mu = mu,
      N = nrow(data),
      n.ext = n.ext,
      sigma2 = obj$sigma2,
      TPR = obj$TPR
    )$par[1]

    if (is.null(sparseGRM)) {
      var.ratio.ext <- 1
    } else {
      R_tilde_w <- objNull$mresid - mean(objNull$mresid) * w.ext
      names(R_tilde_w) <- data$SampleID
      sparseGRM <- sparseGRM %>%
        mutate(cov_Rext = Value * R_tilde_w[as.character(ID1)] * R_tilde_w[as.character(ID2)])
      numerator <- sum(sparseGRM$cov_Rext) + w.ext^2 * sum(objNull$mresid)^2 / n.ext
      denominator <- sum(R_tilde_w^2) + w.ext^2 * sum(objNull$mresid)^2 / n.ext
      var.ratio.ext <- numerator / denominator
    }

    mergeGenoInfo_1 <- mergeGenoInfo_1 %>% cbind(., TPR, sigma2, w.ext, var.ratio.ext)
  }) %>%
    do.call("rbind", .) %>%
    as_tibble() %>%
    arrange(index) %>%
    select(-index)

  return(mergeGenoInfo)
}
