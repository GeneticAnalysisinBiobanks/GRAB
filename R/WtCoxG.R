## ------------------------------------------------------------------------------
## WtCoxG.R
##
## Functions:
##   GRAB.WtCoxG                  : Print brief method information.
##   checkControl.NullModel.WtCoxG: Validate/populate null-model controls.
##   fitNullModel.WtCoxG          : Fit weighted Cox null model & test batch effects.
##   setMarker.WtCoxG             : Initialize marker-level analysis objects.
##   mainMarker.WtCoxG            : Marker-level association testing.
##   TestforBatchEffect           : QC + parameter estimation (TPR, sigma2, weights).
## ------------------------------------------------------------------------------

#' Instruction of WtCoxG method
#'
#' WtCoxG is a Cox-based association test method for time-to-event traits. It effectively 
#' addresses case ascertainment and rare variant analysis. By leveraging external minor 
#' allele frequencies from public resources, WtCoxG can further boost statistical power.
#'
#' @return NULL
#'
#' @examples
#' # Step 1: fit null model and test batch effect
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' RefAfFile <- system.file("extdata", "simuRefAf.txt", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultWtCoxG.txt")
#'
#' obj.WtCoxG <- GRAB.NullModel(
#'   survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   method = "WtCoxG",
#'   traitType = "time-to-event",
#'   GenoFile = GenoFile,
#'   SparseGRMFile = SparseGRMFile,
#'   RefAfFile = RefAfFile,
#'   RefPrevalence = 0.1
#' )
#'
#' # Step2
#' GRAB.Marker(obj.WtCoxG, GenoFile, OutputFile)
#'
#' head(data.table::fread(OutputFile))
#'
#' @details
#'
#' \strong{Additional Parameters for \code{GRAB.NullModel()}:}
#' \itemize{
#'   \item \code{RefAfFile} (character, required): Reference allele frequency file path.
#'     File must contain columns: CHROM, POS, ID, REF, ALT, AF_ref, AN_ref
#'   \item \code{RefPrevalence} (numeric, required): Population-level disease prevalence
#'     for weighting. Must be in range (0, 0.5)
#' }
#'
#' \strong{Additional Control Parameters for GRAB.NullModel()}:
#' \itemize{
#'   \item \code{OutlierRatio} (numeric, default: 1.5): IQR multiplier for outlier detection
#' }
#'
#' \strong{Method-specific elements in the \code{WtCoxG_NULL_Model} object returned by \code{GRAB.NullModel()}:}:
#' \itemize{
#'   \item \code{mresid}: Martingale residuals from weighted Cox model (numeric).
#'   \item \code{Cova}: Design matrix of covariates (matrix).
#'   \item \code{yVec}: Event indicator (numeric).
#'   \item \code{weight}: Observation weights based on reference prevalence (numeric).
#'   \item \code{RefPrevalence}: Reference population prevalence used for weighting (numeric).
#'   \item \code{outLierList}: List identifying outlier subjects for SPA adjustment.
#'   \item \code{mergeGenoInfo}: Data frame with batch effect QC results and external reference data.
#' }
#' 
#' \strong{Additional Control Parameters for GRAB.Marker()}:
#' \itemize{
#'   \item \code{cutoff} (numeric, default: 0.1): Cutoff of batch effect test p-value for
#'     association testing. Variants with batch effect p-value below this cutoff
#'     will be excluded from association testing.
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
#' @references
#' Li et al. (2025). Applying weighted Cox regression to genome-wide association studies of 
#' time-to-event phenotypes. \doi{10.1038/s43588-025-00864-z}
#'
GRAB.WtCoxG <- function() {
  .message("?GRAB.WtCoxG for instructions")
}


checkControl.NullModel.WtCoxG <- function(traitType, GenoFile, SparseGRMFile, control, ...) {

  if (!traitType %in% c("time-to-event")) {
    stop("For 'WtCoxG' method, only traitType of 'time-to-event' is supported.")
  }

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
  } else {
    optionGRM <- NULL
  }

  # Validate RefPrevalence (required)
  RefPrevalence <- list(...)$RefPrevalence
  if (is.null(RefPrevalence)) {
    stop("Argument 'RefPrevalence' is required for WtCoxG method.")
  }
  if (!is.numeric(RefPrevalence) || length(RefPrevalence) != 1) {
    stop("Argument 'RefPrevalence' should be a single numeric value.")
  }
  if (RefPrevalence <= 0 || RefPrevalence >= 0.5) {
    stop("Argument 'RefPrevalence' should be between (0, 0.5).")
  }

  # Validate RefAfFile (required)
  RefAfFile = list(...)$RefAfFile
  if (is.null(RefAfFile)) {
    stop("Argument 'RefAfFile' is required for WtCoxG method.")
  }
  if (!is.character(RefAfFile) || length(RefAfFile) != 1) {
    stop("Argument 'RefAfFile' should be a character string (file path).")
  }
  if (!file.exists(RefAfFile)) {
    stop("Cannot find RefAfFile: ", RefAfFile)
  }

  # Validate RefAfFile columns (read header only)
  refGenoHeader <- colnames(data.table::fread(RefAfFile, nrows = 0))
  requiredCols <- c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref", "AN_ref")
  missingCols <- setdiff(requiredCols, refGenoHeader)
  if (length(missingCols) > 0) {
    stop("RefAfFile must contain columns: ", paste(requiredCols, collapse = ", "),
         ". Missing columns: ", paste(missingCols, collapse = ", "))
  }

  default.control <- list(
    OutlierRatio = 1.5
  )
  control <- updateControl(control, default.control)

  return(list(control = control, optionGRM = optionGRM))
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
#' @param control List with fields such as \code{OutlierRatio}.
#' @param data Data frame used for optional batch-effect QC.
#' @param GenoFile Character. PLINK prefix (without extension) used when
#'   sampling markers for QC.
#' @param GenoFileIndex Character. Path to an index file used in QC workflow.
#' @param SparseGRMFile Character. Path to sparse GRM used in QC workflow.
#' @param SurvTimeColumn Character. Column name in \code{data} containing survival time.
#' @param IndicatorColumn Character. Column name in \code{data} containing event indicator (0/1).
#' @param ... Optional named parameters forwarded to QC (e.g., \code{RefAfFile}).
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
  GenoFile, GenoFileIndex, SparseGRMFile, 
  SurvTimeColumn, IndicatorColumn, ...
) {

  dots <- list(...)
  RefAfFile <- dots$RefAfFile
  RefPrevalence <- dots$RefPrevalence
  
  Indicator <- as.matrix(response)[, "status"]
  formula <- response ~ designMat

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

  ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
  q25 <- quantile(mresid, 0.25, na.rm = TRUE)
  q75 <- quantile(mresid, 0.75, na.rm = TRUE)
  IQR <- q75 - q25

  r.outlier <- ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)
  outlier_bounds <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR) # 
  posOutlier <- which(mresid < outlier_bounds[1] | mresid > outlier_bounds[2])

  while (length(posOutlier) == 0) {
    r.outlier <- r.outlier * 0.8
    outlier_bounds <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier <- which(mresid < outlier_bounds[1] | mresid > outlier_bounds[2])
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

  # Test for batch effect using validated parameters
  re$mergeGenoInfo <- TestforBatchEffect(
    objNull = re,
    data = data,
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    SparseGRMFile = SparseGRMFile,
    RefAfFile = RefAfFile,
    IndicatorColumn = IndicatorColumn,
    SurvTimeColumn = SurvTimeColumn,
    RefPrevalence = RefPrevalence
  )

  class(re) <- "WtCoxG_NULL_Model"
  return(re)
}


checkControl.Marker.WtCoxG <- function(control) {

  default.control <- list(
    cutoff = 0.1
  )
  control <- updateControl(control, default.control)

  if (!is.numeric(control$cutoff) || control$cutoff < 0 || control$cutoff > 1) {
    stop("control$cutoff should be a numeric value in [0, 1].")
  }

  return(control)
}


setMarker.WtCoxG <- function(objNull, control) {

  setWtCoxGobjInCPP(
    t_mresid = objNull$mresid,               # numeric vector: Martingale residuals from Cox model
    t_weight = objNull$weight,               # numeric vector: Weight vector for analysis
    t_cutoff = control$cutoff,               # numeric: batch effect p-value cutoff for association testing
    t_SPA_Cutoff = control$SPA_Cutoff        # numeric: P-value cutoff for SPA
  )
}


mainMarker.WtCoxG <- function(genoType, genoIndex, objNull) {

  mergeGenoInfo <- objNull$mergeGenoInfo
  mergeGenoInfo_chunk <- mergeGenoInfo[mergeGenoInfo$genoIndex %in% genoIndex, ]
  cols <- c("AF_ref", "AN_ref", "TPR", "sigma2", "pvalue_bat", "w.ext", "var.ratio.w0", "var.ratio.int", "var.ratio.ext")

  OutList <- mainMarkerInCPP(
    t_method = "WtCoxG",      # character: Statistical method name
    t_genoType = genoType,    # character: "PLINK" or "BGEN"
    t_genoIndex = genoIndex,   # numeric vector: Genotype indices to analyze
    t_extraParams = list("mergeGenoInfo_chunk" = mergeGenoInfo_chunk[, cols]) # list: additional parameters for the method
  )

  pvals <- data.frame(matrix(OutList$pvalVec, ncol = 2, byrow = TRUE))
  colnames(pvals) <- c("WtCoxG.ext", "WtCoxG.noext")

  obj.mainMarker <- cbind(pvals, mergeGenoInfo_chunk)
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
  SurvTimeColumn,
  RefPrevalence # refernce population prevalence

) {
  .message("Testing for batch effect ...")

  if (!is.null(SparseGRMFile)) {
    sparseGRM <- data.table::fread(SparseGRMFile)
  } else {
    sparseGRM <- NULL
  }

  control <- objNull$control

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

    # All rows from GenoInfo will be kept. The result will have as many rows as GenoInfo.
    # Only matching rows from refGenoInfo (based on CHROM, POS, RA) will be included.
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
  
  maf.group <- sort(unique(
    c(seq(0, 0.4, 0.05), max(mergeGenoInfo$mu.int, na.rm = TRUE))
  ))
  
  mergeGenoInfo.est <- lapply(seq_len(length(maf.group) - 1), function(i) {
    
    mergeGenoInfo_1 <- mergeGenoInfo %>%
      filter(mu.int > maf.group[i], mu.int <= maf.group[i + 1])
    
    ## No SNP in this group
    if (nrow(mergeGenoInfo_1) = 0) {
      .message(" No in this MAF group, skipping")
      return(NULL)
    }
    
    mergeGenoInfo_2 <- mergeGenoInfo %>%
      filter(
        mu.int >= max(maf.group[i] - 0.1, 0.01),
        mu.int <  min(1, maf.group[i + 1] + 0.1)
      )
    
    mu    <- (maf.group[i] + maf.group[i + 1]) / 2
    n.ext <- mean(na.omit(mergeGenoInfo_1$AN_ref)) / 2
    
    if (!is.finite(n.ext) || n.ext < 5) {
      return(
        mergeGenoInfo_1 %>%
          dplyr::mutate(TPR = NA, sigma2 = NA, w.ext = 0, var.ratio.ext = NA)
      )
    }
    
    ## ---- variance parts ----
    var_mu_ext <- mu * (1 - mu) / (2 * n.ext)
    
    if (is.null(sparseGRM)) {
      var_Sbat <- sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext
    } else {
      vr <- na.omit(mergeGenoInfo_2$var.ratio.w0)
      if (length(vr) == 0) return(NULL)
      var_Sbat <- mean(vr) * (sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext)
    }
    
    ## ---- estimate TPR & sigma2 ----
    vec_p_bat <- na.omit(mergeGenoInfo_2$pvalue_bat)
    if (length(vec_p_bat) < 50) return(NULL)
    
    vec_cutoff <- c(0.01,  0.1, 0.2)
    vec_p_deno <- sapply(vec_cutoff, function(p) mean(vec_p_bat > p))
    
    ## —— guard 2: remove degenerate cutoffs
    keep <- which(is.finite(vec_p_deno) & vec_p_deno > 0)
    if (length(keep) < 2) return(NULL)
    
    vec_cutoff  <- vec_cutoff[keep]
    vec_p_deno  <- vec_p_deno[keep]
    
    opti_fun <- function(par) {
      TPR    <- par[1]
      sigma2 <- par[2]
      if (TPR < 0 || TPR > 1 || sigma2 < 0) return(1e6)
      
      sum(sapply(seq_along(vec_cutoff), function(j) {
        p_cut <- vec_cutoff[j]
        lb <- -qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
        ub <-  qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
        
        p_deno <- vec_p_deno[j]
        
        v2 <- var_Sbat + sigma2
        c <- pnorm(ub, 0, sqrt(v2), log.p = TRUE)
        d <- pnorm(lb, 0, sqrt(v2), log.p = TRUE)
        
        pro.cut <- TPR * (exp(d) * (exp(c - d) - 1)) +
          (1 - TPR) * (1 - p_cut)
        
        ((p_deno - pro.cut) / p_deno)^2
      }))
    }
    
    opt <- tryCatch(
      optim(
        par = c(0.05, 0.01),
        fn = opti_fun,
        method = "L-BFGS-B",
        lower = c(0, 0),
        upper = c(1, Inf)
      ),
      error = function(e) NULL
    )
    
    if (is.null(opt) || opt$convergence != 0) {
      TPR <- NA
      sigma2 <- NA
    } else {
      TPR    <- opt$par[1]
      sigma2 <- opt$par[2]
    }
    
    ## ---- estimate w.ext (never crash) ----
    w.ext <- tryCatch(
      optim(
        par = 0.5,
        method = "L-BFGS-B",
        lower = 0, upper = 1,
        fn = fun.optimalWeight,
        pop.prev = RefPrevalence,
        y = data$Indicator,
        R = objNull$mresid,
        w = objNull$weight,
        mu = mu,
        N = nrow(data),
        n.ext = n.ext,
        sigma2 = sigma2,
        TPR = TPR
      )$par[1],
      error = function(e) 0
    )
    
    if (is.null(sparseGRM)) {
      var.ratio.ext <- 1
    } else {
      R_tilde_w <- objNull$mresid - mean(objNull$mresid) * w.ext
      names(R_tilde_w) <- data$SampleID
      
      sparseGRM$tmp <- sparseGRM$Value *
        R_tilde_w[as.character(sparseGRM$ID1)] *
        R_tilde_w[as.character(sparseGRM$ID2)]
      
      var.ratio.ext <- sum(sparseGRM$tmp, na.rm = TRUE) /
        sum(R_tilde_w^2, na.rm = TRUE)
    }
    
    cbind(
      mergeGenoInfo_1,
      TPR = TPR,
      sigma2 = sigma2,
      w.ext = w.ext,
      var.ratio.ext = var.ratio.ext
    )
  }) %>%
    Filter(Negate(is.null), .) %>%
    dplyr::bind_rows() %>%
    dplyr::arrange(index) %>%
    dplyr::select(-index)
  
  mergeGenoInfo <- mergeGenoInfo.est
  return(as.data.frame(mergeGenoInfo))
}
