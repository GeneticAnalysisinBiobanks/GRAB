################################################################################
# LEAF.R - LEAF Method Implementation
################################################################################
#
# This file contains the complete implementation of the LEAF (Logistic regression
# with External Ancestry estimation and adjustment Framework) method for binary 
# and time-to-event trait association testing in admixed or multi-ancestry populations.
#
# MAIN WORKFLOW:
# 1. checkControl.NullModel.LEAF: Validates input parameters and reference files
# 2. fitNullModel.LEAF: Performs k-means clustering and computes null model per cluster
# 3. checkControl.Marker.LEAF: Validates marker-level control parameters
# 4. setMarker.LEAF: Sets up C++ objects for marker testing
# 5. mainMarker.LEAF: Performs association testing and meta-analysis
#
# HELPER FUNCTIONS:
# - LEAF.Cluster: K-means clustering and ancestry proportion estimation via Summix
# - LEAF.fitWtGLM: Weighted GLM fitting for binary traits
# - LEAF.fitWtCox: Weighted Cox regression for time-to-event traits
# - LEAF.BatchEffect: Tests for batch effects between internal and external AFs
# - summix/summix_calc/calc_scaledObj: Ancestry proportion estimation (from Summix package)
#
# KEY CONCEPTS:
# - Ancestry-informed clustering partitions samples by genetic ancestry
# - External reference AFs are ancestry-matched for each cluster
# - Batch effect testing identifies markers with systematic AF differences
# - Meta-analysis combines evidence across clusters (ext/noext)
#
################################################################################

#' Instruction of LEAF method
#'
#' LEAF is a method for binary trait association testing in admixed or multi-ancestry populations. 
#' It partitions samples into ancestry-informed clusters, estimates ancestry-specific external 
#' allele frequencies, and performs cluster-specific association tests to improve power and 
#' control for population stratification.
#'
#' @return NULL
#'
#' @examples
#' # Step 1: fit null model and test batch effect
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' RefAfFile <- system.file("extdata", "simuRefAf_2pop.txt", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultLEAF.txt")
#'
#' obj.LEAF <- GRAB.NullModel(
#'   BinaryPheno ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   method = "LEAF",
#'   traitType = "binary",
#'   GenoFile = GenoFile,
#'   SparseGRMFile = SparseGRMFile,
#'   RefAfFile = RefAfFile,
#'   RefPrevalence = 0.1,
#'   Ncluster = 2,
#'   PCmatrix = PhenoData[, c("PC1", "PC2", "PC3", "PC4")]
#' )
#'
#' # Step 2: Perform marker-level association testing
#' GRAB.Marker(obj.LEAF, GenoFile, OutputFile)
#'
#' head(data.table::fread(OutputFile))
#'
#' @details
#'
#' \strong{Additional Parameters for \code{GRAB.NullModel()}:}
#' \itemize{
#'   \item \code{RefAfFile} (character, required): Reference allele frequency file path.
#'     File must contain columns: CHROM, POS, ID, REF, ALT, and population-specific AF/AN column
#'     pairs (e.g., EUR_AF & EUR_AN, AFR_AF & AFR_AN). Column naming is flexible: both underscore 
#'     ("_") and dot (".") separators are supported, and matching is case-insensitive 
#'     (e.g., EUR_AF/EUR_AN or eur.af/eur.an). Each AF column must have a corresponding AN column.
#'   \item \code{RefPrevalence} (numeric, required): Population-level disease prevalence
#'     for weighting. Must be in range (0, 0.5)
#'   \item \code{Ncluster} (integer, required): Number of ancestry clusters to create using
#'     k-means clustering on the first 4 PCs. Typically 2-5 depending on ancestry diversity
#'   \item \code{PCmatrix} (matrix or data.frame, required): Containing the principal components
#'     for ancestry clustering. Recommended to include the first 4 PCs.
#' }
#'
#' \strong{Additional Control Parameters for GRAB.Marker()}:
#' \itemize{
#'   \item \code{cutoff} (numeric, default: 0.05): P-value cutoff for batch effect testing.
#'     Markers with batch effect p-value < cutoff will be excluded from analysis with external AFs.
#' }
#'
#' \strong{Output Columns from GRAB.Marker()}:
#'
#' The output contains meta analysis p-values, and for each cluster, the output contains:
#' \itemize{
#'   \item \code{p.ext.clus[i]} (numeric): P-value with external reference allele frequencies
#'   \item \code{p.noext.clus[i]} (numeric): P-value without external reference (internal-only)
#'   \item \code{score.ext.clus[i]} (numeric): Score statistic with external reference
#'   \item \code{score.noext.clus[i]} (numeric): Score statistic without external reference
#'   \item \code{zscore.ext.clus[i]} (numeric): Z-score statistic with external reference
#'   \item \code{zscore.noext.clus[i]} (numeric): Z-score statistic without external reference
#' }
#'
#' @references
#' Ying Li et al. (in prep). Leveraging external allele frequency to boost powers of genome-wide association studies.
#'
GRAB.LEAF <- function() {
  .message("?GRAB.LEAF for instructions")
}


checkControl.NullModel.LEAF <- function(traitType, GenoFile, SparseGRMFile, control, ...) {

  if (!traitType %in% c("binary", "time-to-event")) {
    stop("For 'LEAF' method, only traitType of 'binary' or 'time-to-event' is supported.")
  }

  # GenoFile validation (required)
  if (is.null(GenoFile)) {
    stop("Argument 'GenoFile' is required for method 'LEAF'.")
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

  # Validate Ncluster (required)
  dots <- list(...)
  Ncluster <- dots$Ncluster
  if (is.null(Ncluster)) {
    stop("Argument 'Ncluster' is required for LEAF method.")
  }
  if (!is.numeric(Ncluster) || length(Ncluster) != 1 || Ncluster < 1) {
    stop("Argument 'Ncluster' should be a positive integer.")
  }

  # Validate RefPrevalence (required)
  RefPrevalence <- dots$RefPrevalence
  if (is.null(RefPrevalence)) {
    stop("Argument 'RefPrevalence' is required for LEAF method.")
  }
  if (!is.numeric(RefPrevalence) || length(RefPrevalence) != 1) {
    stop("Argument 'RefPrevalence' should be a single numeric value.")
  }
  if (RefPrevalence <= 0 || RefPrevalence >= 0.5) {
    stop("Argument 'RefPrevalence' should be between (0, 0.5).")
  }

  # Validate RefAfFile (required)
  RefAfFile <- dots$RefAfFile
  if (is.null(RefAfFile)) {
    stop("Argument 'RefAfFile' is required for LEAF method.")
  }
  if (!is.character(RefAfFile) || length(RefAfFile) != 1) {
    stop("Argument 'RefAfFile' should be a character string (file path).")
  }
  if (!file.exists(RefAfFile)) {
    stop("Cannot find RefAfFile: ", RefAfFile)
  }

  ## Validate RefAfFile columns (read header only)
  # Uppercase all column names and replace dots with underscores
  refGenoHeader <- colnames(data.table::fread(RefAfFile, nrows = 0))
  refGenoHeader <- toupper(gsub("\\.", "_", refGenoHeader))
  requiredCols <- c("CHROM", "POS", "ID", "REF", "ALT")
  
  # Check for required columns
  missingCols <- setdiff(requiredCols, refGenoHeader)
  if (length(missingCols) > 0) {
    stop("RefAfFile must contain columns: ", paste(requiredCols, collapse = ", "),
         ". Missing columns: ", paste(missingCols, collapse = ", "))
  }
  
  # Check AF columns
  afCols <- grep("_AF$", refGenoHeader, value = TRUE)
  if (length(afCols) == 0) {
    stop("RefAfFile must contain at least one AF column (e.g., 'EUR_AF', 'AFR_AF').")
  }
  
  if (any(duplicated(afCols))) {
    stop("RefAfFile contains duplicated AF columns: ", paste(afCols[duplicated(afCols)], collapse = ", "))
  }
  
  # For each AF column, check corresponding AN column
  for (afCol in afCols) {
    prefix <- sub("_AF$", "", afCol)
    expectedAN <- paste0(prefix, "_AN")
    
    if (!expectedAN %in% refGenoHeader) {
      stop("RefAfFile has AF column '", afCol, "' but missing corresponding AN column '", expectedAN, "'.")
    }
    
    anCols <- grep(paste0("^", prefix, "_AN$"), refGenoHeader, value = TRUE)
    if (length(anCols) > 1) {
      stop("RefAfFile has duplicated AN columns for population '", prefix, "': ", paste(anCols, collapse = ", "))
    }
  }

  ## Validate PCmatrix (required)
  PCmatrix <- dots$PCmatrix
  if (is.null(PCmatrix)) {
    stop("PCmatrix is required for kmeans clustering in LEAF method.")
  }
  if (!is.matrix(PCmatrix) && !is.data.frame(PCmatrix)) {
    stop("PCmatrix must be a matrix or data.frame.")
  }

  # default.control <- list(
  #   OutlierRatio = 1.5
  # )
  # control <- updateControl(control, default.control)

  return(list(control = control, optionGRM = optionGRM))
}


# kmeans clustering based on PCs
# Summix: for each subcohort k, estimate proportion of ancestry-specific external AFs
# calculate ancestry-matched external AFs according to output of Summix
fitNullModel.LEAF <- function(
  response,
  designMat, 
  subjData, 
  control, 
  GenoFile = NULL,
  GenoFileIndex = NULL,
  SparseGRMFile = NULL,
  ...
) {
  dots <- list(...)
  Ncluster = dots$Ncluster
  RefAfFile = dots$RefAfFile
  RefPrevalence = dots$RefPrevalence
  PCmatrix = dots$PCmatrix

  colsForBatchEffect <- c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref", "AN_ref")
  colsForStep2 <- c("genoIndex", "AF_ref", "AN_ref", "TPR", "sigma2", "pvalue_bat", 
                    "w.ext", "var.ratio.w0", "var.ratio.int", "var.ratio.ext")

  SparseGRM <- if (!is.null(SparseGRMFile)) data.table::fread(SparseGRMFile) else NULL

  # Get cluster assignments and ancestry-matched reference AFs for each cluster
  .message("Performing k-means clustering with K=%d ...", Ncluster)
  clusterResult <- LEAF.Cluster(
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    control = control,
    Ncluster = Ncluster,
    subjData = subjData,
    PCmatrix = PCmatrix,
    RefAfFile = RefAfFile
  )
  
  clusterIdx <- clusterResult$clusterIdx  # vector of cluster assignments for each subject
  subGenoInfo_lst <- clusterResult$subGenoInfo_lst

  .message("Fitting null models and testing for batch effects in each cluster ...")
  residuals_list <- vector("list", Ncluster)
  weights_list <- vector("list", Ncluster)
  subGenoInfo <- vector("list", Ncluster)
  for (i in 1:Ncluster) {
    # Subset data for this cluster
    clusterIdx_i <- which(clusterIdx == i)
    subjData_i <- subjData[clusterIdx_i]
    subGenoInfo_i <- subGenoInfo_lst[[i]]

    # Subset sparse GRM for this cluster if available
    if (!is.null(SparseGRM)) {
      SparseGRM_i <- SparseGRM[SparseGRM$ID1 %in% subjData_i & SparseGRM$ID2 %in% subjData_i, ]
    } else {
      SparseGRM_i <- NULL
    }

    # Calculate weights and residuals
    if (inherits(response, "Surv")) {
      weightsResids <- LEAF.fitWtCox(
        response = response[clusterIdx_i, ],
        designMat = designMat[clusterIdx_i, ],
        RefPrevalence = RefPrevalence
      )

      Indicator <- response[clusterIdx_i, "status"]
    } else if (is.vector(response) || is.factor(response)) {
      weightsResids <- LEAF.fitWtGLM(
        response = response[clusterIdx_i],
        designMat = designMat[clusterIdx_i, ],
        RefPrevalence = RefPrevalence
      )

      Indicator <- response[clusterIdx_i]
    } else {
      stop("Internal error: unsupported response class, ", class(response))
    }

    weights <- weightsResids$weights
    residuals <- weightsResids$residuals
    residuals_list[[i]] <- residuals
    weights_list[[i]] <- weights

    .message("Cluster %d: Finished fitting null model. Performing batch effect testing...", i)
    mergeGenoInfo <- LEAF.BatchEffect(
      Indicator = Indicator,
      SampleID = subjData[clusterIdx_i],
      control = control,
      GenoFile = GenoFile,
      GenoFileIndex = GenoFileIndex,
      SparseGRM = SparseGRM_i,
      refGenoInfo = subGenoInfo_i[, colsForBatchEffect],
      RefPrevalence = RefPrevalence,
      weights = weights,
      residuals = residuals
    )

    # Store results in lists
    subGenoInfo[[i]] <- mergeGenoInfo[, colsForStep2]

    .message("Cluster %d: Completed batch effect testing.", i)
    # outLierList <- getOutlierList(residuals, control$OutlierRatio)

  }
  .message("Finished fitting null models and batch effect testing for all clusters.")

  objNull <- list(
    residuals_list = residuals_list,
    weights_list = weights_list,
    subjData = subjData,
    clusterIdx = clusterIdx,
    subGenoInfo = subGenoInfo,
    Ncluster = Ncluster
  )

  class(objNull) <- "LEAF_NULL_Model"
  return(objNull)
}

checkControl.Marker.LEAF <- function(control) {

  default.control <- list(
    cutoff = 0.05
  )
  control <- updateControl(control, default.control)

  if (!is.numeric(control$cutoff) || control$cutoff <= 0 || control$cutoff >= 1) {
    stop("control$cutoff should be a numeric value in (0, 1).")
  }

  return(control)
}


setMarker.LEAF <- function(objNull, control) {

  # Convert clusterIdx to list of index vectors (one per cluster)
  clusterIdx_list <- lapply(1:objNull$Ncluster, function(i) {
    which(objNull$clusterIdx == i)
  })

  setLEAFobjInCPP(
    t_residuals = objNull$residuals_list,
    t_weight = objNull$weights_list,
    t_clusterIdx = clusterIdx_list,
    t_cutoff = control$cutoff,             # numeric: batch effect p-value cutoff for association testing
    t_SPA_Cutoff = control$SPA_Cutoff      # numeric: P-value cutoff for SPA
  )

}


mainMarker.LEAF <- function(genoType, genoIndex, objNull) {

  Ncluster <- objNull$Ncluster
  cols <- c("AF_ref", "AN_ref", "TPR", "sigma2", "pvalue_bat", "w.ext", "var.ratio.w0", "var.ratio.int", "var.ratio.ext")

  # Prepare extraParams with subGenoInfo for each cluster (index-based list)
  extraParams <- vector("list", Ncluster)
  for (i in 1:Ncluster) {
    subGenoInfo <- objNull$subGenoInfo[[i]]
    extraParams[[i]] <- subGenoInfo[subGenoInfo$genoIndex %in% genoIndex, cols]
  }

  OutList <- mainMarkerInCPP(
    t_method = "LEAF",         # character: Statistical method name
    t_genoType = genoType,     # character: "PLINK" or "BGEN"
    t_genoIndex = genoIndex,   # numeric vector: Genotype indices to analyze
    t_extraParams = extraParams # list: subGenoInfo and clusterIdx for each cluster
  )

  obj.mainMarker <- data.frame(
    Marker = OutList$markerVec,
    Info = OutList$infoVec,
    AltFreq = OutList$altFreqVec,
    # AltCounts = OutList$altCountsVec,
    MissingRate = OutList$missingRateVec
  )
 
  ## Compute meta-analysis statistics (ext and noext separately)
  compute_meta <- function(scores, pvals) {
    # input matrices: rows=markers, cols=clusters
    chisq <- qchisq(pvals, df = 1, lower.tail = FALSE)
    chisq[chisq < 1e-30] <- 1e-30
    var <- (scores^2) / chisq
    var[is.na(var)] <- 0
    sum_scores <- rowSums(scores, na.rm = TRUE)
    sum_var <- rowSums(var, na.rm = TRUE)
    Z_meta <- sum_scores / sqrt(sum_var)
    Z_meta[sum_var == 0] <- 0
    P_meta <- 2 * pnorm(-abs(Z_meta))
    return(list(score = sum_scores, var = sum_var, z = Z_meta, p = P_meta))
  }
  
  # Extract matrices and convert to data frame
  pvals_mat <- matrix(OutList$pvalVec, ncol = 2 * Ncluster, byrow = TRUE)
  scores_mat <- matrix(OutList$score, ncol = 2 * Ncluster, byrow = TRUE)
  zscores_mat <- matrix(OutList$zScore, ncol = 2 * Ncluster, byrow = TRUE)

  # Meta-analysis for ext and noext
  ext_idx <- seq(1, 2 * Ncluster, by = 2)   # odd indices (1,3,5,...)
  noext_idx <- seq(2, 2 * Ncluster, by = 2) # even indices (2,4,6,...)

  meta_ext <- compute_meta(scores_mat[, ext_idx, drop = FALSE], pvals_mat[, ext_idx, drop = FALSE])
  meta_noext <- compute_meta(scores_mat[, noext_idx, drop = FALSE], pvals_mat[, noext_idx, drop = FALSE])
  
  # Add meta-analysis p-values
  obj.mainMarker$meta.p_ext <- meta_ext$p
  obj.mainMarker$meta.p_noext <- meta_noext$p
  
  # Add cluster-specific p-values: p.ext, p.noext, p.batch for each cluster
  for (i in 1:Ncluster) {
    cluster_prefix <- paste0("cl", i)
    
    # Extract p.ext and p.noext from pvals_mat
    obj.mainMarker[[paste0(cluster_prefix, ".p_ext")]] <- pvals_mat[, 2*i - 1]
    obj.mainMarker[[paste0(cluster_prefix, ".p_noext")]] <- pvals_mat[, 2*i]
    
    # Extract p.batch from subGenoInfo
    subGenoInfo <- objNull$subGenoInfo[[i]]
    # Match markers by genoIndex
    batch_pvals <- subGenoInfo$pvalue_bat[match(genoIndex, subGenoInfo$genoIndex)]
    obj.mainMarker[[paste0(cluster_prefix, ".p_batch")]] <- batch_pvals
  }

  return(obj.mainMarker)
}

# ======================================================================================

#' Fit Weighted GLM for Binary Traits
#'
#' Performs weighted logistic regression for binary traits using case-control
#' weighting based on reference population prevalence.
#'
#' @param response Numeric vector of binary outcomes (0/1)
#' @param designMat Design matrix including covariates
#' @param RefPrevalence Population-level disease prevalence (0 < x < 0.5)
#'
#' @return List with two components:
#'   \item{weights}{Regression weights for each sample}
#'   \item{residuals}{Model residuals}
#'
#' @keywords internal
LEAF.fitWtGLM <- function(response, designMat, RefPrevalence) {

  weights <- calRegrWeight(response, RefPrevalence)

  residuals <- glm.fit(
    designMat, response, weights = weights,
    family = binomial(), intercept = TRUE
  )$residuals

  return(list(
    weights= weights,
    residuals = residuals
  ))
}

#' Fit Weighted Cox Regression for Time-to-Event Traits
#'
#' Performs weighted Cox proportional hazards regression for survival outcomes
#' using case-control weighting based on reference population prevalence.
#'
#' @param response Surv object containing time and status
#' @param designMat Design matrix including covariates
#' @param RefPrevalence Population-level disease prevalence (0 < x < 0.5)
#'
#' @return List with two components:
#'   \item{weights}{Regression weights for each sample}
#'   \item{residuals}{Cox model residuals}
#'
#' @keywords internal
LEAF.fitWtCox <- function(response, designMat, RefPrevalence) {
  
  Indicator <- response[, "status"]
  weights <- calRegrWeight(Indicator, RefPrevalence)

  residuals <- obj.coxph <- survival::coxph(
    response ~ designMat, x = TRUE, 
    weight = weights, robust = TRUE
  )$residuals

  return(list(
    weights = weights,
    residuals = residuals
  ))
}

calRegrWeight <- function(Indicator, RefPrevalence) {
  if (any(!unique(Indicator) %in% c(0, 1))) {
    stop("The value of Indicator should be 0 or 1.")
  }

  sumOnes <- sum(Indicator)
  sumZeros <- sum(1 - Indicator)
  ratio <- sumOnes / sumZeros

  weights<- ifelse(
    Indicator == 1,
    1,
    (1 - RefPrevalence) / RefPrevalence * ratio
  )

  return(weights)
}


#' Test for Batch Effects Between Internal and External Allele Frequencies
#'
#' Compares internal sample AFs with external reference AFs to detect systematic
#' differences that could indicate batch effects. Estimates true positive rate (TPR),
#' variance inflation (sigma2), and optimal external weight contribution for each
#' marker grouped by MAF.
#'
#' @param Indicator Binary outcome vector (0/1 or case/control status)
#' @param SampleID Character vector of sample identifiers
#' @param control Control parameters list
#' @param GenoFile Path to genotype file (PLINK/BGEN)
#' @param GenoFileIndex Optional genotype file index
#' @param SparseGRM Sparse genetic relationship matrix (data.frame with ID1, ID2, Value)
#' @param refGenoInfo Reference genotype info with columns: CHROM, POS, ID, REF, ALT, AF_ref, AN_ref
#' @param RefPrevalence Population-level disease prevalence
#' @param weights Regression weights from null model
#' @param residuals Model residuals from null model
#'
#' @return Data.frame with merged genotype info including batch effect test results:
#'   TPR, sigma2, pvalue_bat, w.ext, var.ratio.w0, var.ratio.int, var.ratio.ext
#'
#' @keywords internal
LEAF.BatchEffect <- function(
  Indicator,
  SampleID,
  control = NULL,
  GenoFile,
  GenoFileIndex = NULL,
  SparseGRM = NULL,
  refGenoInfo,
  RefPrevalence,
  weights,
  residuals
) {

  refGenoInfo <- as.data.frame(refGenoInfo)

  # check if there are 7 columns in RefAfFile
  for (colname in c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref", "AN_ref")) {
    if (!colname %in% colnames(refGenoInfo)) {
      stop(paste0(colname, " is missing in RefAfFile!"))
    }
  }

  ## merge sample genoInfo and ref genoInfo
  # .message("Getting genotype info for controls (Indicator=0) ...")
  GenoInfo.ctrl <- GRAB.getGenoInfo(
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    SampleIDs = SampleID[Indicator == 0],
    control = control
  )
  GenoInfo.ctrl <- GenoInfo.ctrl[, c("altFreq", "missingRate")]
  colnames(GenoInfo.ctrl) <- c("mu0", "mr0")

  # .message("Getting genotype info for cases (Indicator=1) ...")
  GenoInfo <- GRAB.getGenoInfo(
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    SampleIDs = SampleID[Indicator == 1],
    control = control
  )
  colnames(GenoInfo)[colnames(GenoInfo) == "altFreq"] <- "mu1"
  colnames(GenoInfo)[colnames(GenoInfo) == "missingRate"] <- "mr1"
  GenoInfo <- cbind(GenoInfo, GenoInfo.ctrl)
  GenoInfo$RA <- paste0(pmin(GenoInfo$REF, GenoInfo$ALT), pmax(GenoInfo$REF, GenoInfo$ALT))
  GenoInfo$index <- seq_len(nrow(GenoInfo))

  # All rows from GenoInfo will be kept. The result will have as many rows as GenoInfo.
  # Only matching rows from refGenoInfo (based on CHROM, POS, RA) will be included.
  refGenoInfo$RA <- paste0(pmin(refGenoInfo$REF, refGenoInfo$ALT), pmax(refGenoInfo$REF, refGenoInfo$ALT))
  mergeGenoInfo <- merge(refGenoInfo, GenoInfo, by = c("CHROM", "POS", "RA"), all.y = TRUE, sort = FALSE)
  mergeGenoInfo$REF <- mergeGenoInfo$REF.y
  mergeGenoInfo$ALT <- mergeGenoInfo$ALT.y
  mergeGenoInfo$ID <- mergeGenoInfo$ID.y
  mergeGenoInfo$AF_ref <- ifelse(mergeGenoInfo$REF == mergeGenoInfo$REF.x, 
                                   mergeGenoInfo$AF_ref, 
                                   1 - mergeGenoInfo$AF_ref)
  mergeGenoInfo <- mergeGenoInfo[, !(colnames(mergeGenoInfo) %in% c("REF.x", "ALT.x", "ID.x", "REF.y", "ALT.y", "ID.y", "RA"))]
  mergeGenoInfo <- mergeGenoInfo[!duplicated(mergeGenoInfo$index), ]
  mergeGenoInfo$n1 <- sum(Indicator) * (1 - mergeGenoInfo$mr1)
  mergeGenoInfo$n0 <- sum(1 - Indicator) * (1 - mergeGenoInfo$mr0)
  mergeGenoInfo$mu.int <- 0.5 * mergeGenoInfo$mu1 + 0.5 * mergeGenoInfo$mu0
  mergeGenoInfo$mu.int <- ifelse(mergeGenoInfo$mu.int > 0.5, 1 - mergeGenoInfo$mu.int, mergeGenoInfo$mu.int)

  ## calculate variance ratios
  w1 <- weights / (2 * sum(weights))
  names(w1) <- SampleID
  meanR <- mean(residuals)
  R_tilde <- residuals - meanR
  names(R_tilde) <- SampleID

  if (!is.null(SparseGRM)) {
    SparseGRM$cov <- SparseGRM$Value * w1[as.character(SparseGRM$ID1)] * w1[as.character(SparseGRM$ID2)]
    SparseGRM$cov_R <- SparseGRM$Value * R_tilde[as.character(SparseGRM$ID1)] * R_tilde[as.character(SparseGRM$ID2)]
    var.ratio.w0 <- (sum(SparseGRM$cov) + 1 / (2 * mergeGenoInfo$AN_ref)) /
      (sum(w1^2) + 1 / (2 * mergeGenoInfo$AN_ref))
    var.ratio.int <- sum(SparseGRM$cov_R) / sum(R_tilde^2)
  } else {
    var.ratio.w0 <- var.ratio.int <- 1
  }

  mergeGenoInfo$var.ratio.w0 <- var.ratio.w0
  mergeGenoInfo$var.ratio.int <- var.ratio.int

  ## batch effect p-value calculation
  pvalue_bat <- unlist(lapply(seq_len(nrow(mergeGenoInfo)), function(ind) {

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

    return(p)
  }))

  mergeGenoInfo$pvalue_bat <- pvalue_bat
  rm(pvalue_bat)

  ## estimate statistics according to batch effect p-values
  # .message("Estimating TPR and sigma2 ...")
  maf.group <- c(seq(-1e-4, 0.4, 0.05), max(mergeGenoInfo$mu.int))
  mergeGenoInfo <- do.call("rbind", lapply(1:(length(maf.group) - 1), function(i) {
    # Progress indicator for MAF group processing

    ## assume that genotypes with MAF in [ maf.group[i] , maf.group[i+1]] have the same mixture distribution
    mergeGenoInfo_1 <- mergeGenoInfo[mergeGenoInfo$mu.int > maf.group[i] & mergeGenoInfo$mu.int <= maf.group[i + 1], ]

    ## using batcheffect p-values with MAF in [maf.group[i]-0.1 , maf.group[i+1]+0.1] to estimate parameters
    mergeGenoInfo_2 <- mergeGenoInfo[mergeGenoInfo$mu.int >= max(maf.group[i] - 0.1, 0) & mergeGenoInfo$mu.int < min(1, maf.group[i + 1] + 0.1), ]

    mu <- (maf.group[i] + maf.group[i + 1]) / 2

    n.ext <- mean(na.omit(mergeGenoInfo_1$AN_ref)[1]) / 2
    var_mu_ext <- mu * (1 - mu) / (2 * n.ext)

    var_Sbat <- ifelse(is.null(SparseGRM), sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext,
      na.omit(mergeGenoInfo$var.ratio.w0)[1] * (sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext)
    )

    ## BEGIN inlined: fun.est.param
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

        var_Sbat_par2 <- var_Sbat + par[2]
        if (var_Sbat_par2 >=0) {
          c <- pnorm(ub, 0, sqrt(var_Sbat_par2), log.p = TRUE)
          d <- pnorm(lb, 0, sqrt(var_Sbat_par2), log.p = TRUE)
        } else {
          c <- NaN
          d <- NaN
        }

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
    ## END inlined: fun.est.param

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
      y = Indicator,
      R = residuals,
      w = weights,
      mu = mu,
      N = length(Indicator),
      n.ext = n.ext,
      sigma2 = obj$sigma2,
      TPR = obj$TPR
    )$par[1]

    if (is.null(SparseGRM)) {
      var.ratio.ext <- 1
    } else {
      R_tilde_w <- residuals - mean(residuals) * w.ext
      names(R_tilde_w) <- SampleID
      SparseGRM$cov_Rext <- SparseGRM$Value * R_tilde_w[as.character(SparseGRM$ID1)] * R_tilde_w[as.character(SparseGRM$ID2)]
      var.ratio.ext <- (sum(SparseGRM$cov_Rext) + w.ext^2 * sum(residuals)^2 / n.ext) /
        (sum(R_tilde_w^2) + w.ext^2 * sum(residuals)^2 / n.ext)
    }
 
    return(cbind(mergeGenoInfo_1, TPR, sigma2, w.ext, var.ratio.ext))
  }))
  
  mergeGenoInfo <- mergeGenoInfo[order(mergeGenoInfo$index), ]
  mergeGenoInfo$index <- NULL

  return(as.data.frame(mergeGenoInfo))
}


#' Perform K-means Clustering and Ancestry-Matched Reference AF Synthesis
#'
#' Uses k-means clustering on principal components to partition samples into
#' ancestry-informed clusters. For each cluster, estimates ancestry proportions
#' using Summix and synthesizes ancestry-matched reference allele frequencies
#' and sample sizes from multi-population reference data.
#'
#' @param GenoFile Path to genotype file (PLINK/BGEN)
#' @param GenoFileIndex Optional genotype file index
#' @param Geno.mtx Optional genotype matrix (alternative to GenoFile)
#' @param control Control parameters list (e.g., AlleleOrder)
#' @param Ncluster Number of ancestry clusters (typically 2-5)
#' @param subjData Character vector of subject IDs
#' @param PCmatrix Matrix or data.frame of principal components for clustering (recommend first 4 PCs)
#' @param RefAfFile Path to reference allele frequency file with population-specific AF/AN columns
#'
#' @return List with two components:
#'   \item{clusterIdx}{Integer vector of cluster assignments for each subject}
#'   \item{subGenoInfo_lst}{List of data.frames (one per cluster) with ancestry-matched AF_ref and AN_ref}
#'
#' @details
#' Column naming in RefAfFile is flexible: supports both underscore ("_") and dot (".")
#' separators, case-insensitive matching (e.g., EUR_AF/eur.af both work).
#' Ancestry-matched AF_ref is computed as weighted sum: Sum(Proportion_i * AF_i).
#' Effective sample size AN_ref is computed as: 1 / Sum(Proportion_i^2 / AN_i).
#'
#' @keywords internal
LEAF.Cluster <- function(
  GenoFile = NULL,
  GenoFileIndex = NULL,
  Geno.mtx = NULL,
  control = list(AlleleOrder = "ref-first"),
  Ncluster = 1,
  subjData,
  PCmatrix,
  RefAfFile
) {

  # Detect reference populations
  refGenoInfo <- data.table::fread(RefAfFile)
  colnames(refGenoInfo) <- toupper(gsub("\\.", "_", colnames(refGenoInfo)))
  afCols <- grep("_AF$", colnames(refGenoInfo), value = TRUE)
  pops <- sub("_AF$", "", afCols)
  
  # Population Stratification via K-means
  km_cluster <- kmeans(PCmatrix, centers = Ncluster, nstart = 25)$cluster
  
  # Get marker info for allele alignment (no need to load full genotype matrix)
  # Use GRAB.getGenoInfo to get marker info without loading all genotypes
  # Get marker information without loading full matrix
  .message("Aligning alleles between study data and reference data ...")
  .message("Getting REF and ALT alleles using 1 subject ...")
  GenoInfo <- GRAB.getGenoInfo(
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    SampleIDs = subjData[1],  # Only need 1 subject to get marker info
    control = control
  )[, c("CHROM", "POS", "ID", "REF", "ALT")]
  
  GenoInfo$RA <- paste0(pmin(GenoInfo$REF, GenoInfo$ALT), pmax(GenoInfo$REF, GenoInfo$ALT))
  GenoInfo$index <- seq_len(nrow(GenoInfo))
  
  # Prepare reference info for merge
  refGenoInfo$RA <- paste0(pmin(refGenoInfo$REF, refGenoInfo$ALT), pmax(refGenoInfo$REF, refGenoInfo$ALT))
  
  # Merge by position and allele combination
  mergeGenoInfo <- merge(refGenoInfo, GenoInfo, by = c("CHROM", "POS", "RA"), all.y = TRUE, sort = FALSE)
  
  # Rename columns and identify flips
  mergeGenoInfo$REF <- mergeGenoInfo$REF.y
  mergeGenoInfo$ALT <- mergeGenoInfo$ALT.y
  mergeGenoInfo$ID <- mergeGenoInfo$ID.y
  mergeGenoInfo$flip <- ifelse(mergeGenoInfo$REF == mergeGenoInfo$REF.x, FALSE, TRUE)
  
  # Remove duplicates and clean up
  mergeGenoInfo <- mergeGenoInfo[!duplicated(mergeGenoInfo$index), ]
  mergeGenoInfo <- mergeGenoInfo[, !(colnames(mergeGenoInfo) %in% c("REF.x", "ALT.x", "ID.x", "REF.y", "ALT.y", "ID.y", "RA"))]

  # Apply Allele Flipping for External AF
  for (pop in pops) {
    af_col <- paste0(pop, "_AF")
    if (af_col %in% colnames(mergeGenoInfo)) {
      mergeGenoInfo[[af_col]] <- ifelse(
        mergeGenoInfo$flip, 
        1 - mergeGenoInfo[[af_col]], 
        mergeGenoInfo[[af_col]])
    }
  }
  
  # Cluster-specific Analysis (Summix and AF Synthesis)
  subGenoInfo_lst <- vector("list", Ncluster)

  for (clu in 1:Ncluster) {
    .message("Processing cluster %d of %d: estimating ancestry proportions and ancestry-matched reference AFs...", clu, Ncluster)
    samples_in_clu <- subjData[km_cluster == clu]
    if (length(samples_in_clu) == 0) {
      message(paste("Warning: Cluster", clu, "is empty. Skipping..."))
      next
    }
    
    # Calculate Internal AF for current cluster
    # Memory-efficient: get AF directly without loading full matrix
    cluster_genoInfo <- GRAB.getGenoInfo(
      GenoFile = GenoFile,
      GenoFileIndex = GenoFileIndex,
      SampleIDs = samples_in_clu,
      control = control
    )
    int_af_vec <- setNames(cluster_genoInfo$altFreq, cluster_genoInfo$ID)
  
    # Merge cluster AF with reference info
    subGenoInfo <- mergeGenoInfo
    subGenoInfo$int_AF <- int_af_vec[match(subGenoInfo$ID, names(int_af_vec))]
    
    # Estimate Ancestry Proportions using Summix
    pop_af_cols <- paste0(pops, "_AF")
    summix_cols <- c("int_AF", pop_af_cols)
    summix_input <- subGenoInfo[, summix_cols]
    summix_input <- summix_input[complete.cases(summix_input), ]
    
    if (nrow(summix_input) > 0) {
      anc_prop_res <- summix(
        data = summix_input,
        reference = pop_af_cols,
        observed = "int_AF",
        pi.start = rep(1/length(pops), length(pops))
      )
      anc_prop <- unlist(anc_prop_res)[pop_af_cols]
      anc_prop[is.na(anc_prop)] <- 0
    } else {
      anc_prop <- setNames(rep(0, length(pops)), pop_af_cols)
    }
    
    # Synthesize Ancestry-Matched Reference AF and AN
    # AF_ref = Sum( Proportion_i * AF_i )
    af_mat <- as.matrix(subGenoInfo[, pop_af_cols])
    subGenoInfo$AF_ref <- as.numeric(af_mat %*% anc_prop)
    
    # AN_ref = 1 / Sum( Proportion_i^2 / AN_i )
    pop_an_cols <- paste0(pops, "_AN")
    n_mat <- as.matrix(subGenoInfo[, pop_an_cols])
    n_mat[n_mat <= 0] <- NA
    w2_mat <- matrix(rep(anc_prop^2, each = nrow(n_mat)), nrow = nrow(n_mat))
    denom <- rowSums(w2_mat / n_mat, na.rm = TRUE)
    subGenoInfo$AN_ref <- ifelse(denom > 0, 1 / denom, 0)

    subGenoInfo_lst[[clu]] <- subGenoInfo
  }

  return(list(
    clusterIdx = km_cluster,
    subGenoInfo_lst = subGenoInfo_lst
  ))
}


# ==== Functions from Summix ====

summix <- function(data, 
                   reference, 
                   observed, 
                   pi.start = NA, 
                   goodness.of.fit = TRUE, 
                   override_removeSmallRef = FALSE,
                   network = FALSE,
                   N_reference = NA,
                   reference_colors=NA) {
  start_time = Sys.time()
  if(network) {
    if (length(N_reference) != length(reference)){
      stop("ERROR: Please make sure N_reference is the same length as reference")
    }
  }
  if (length(reference_colors) != 1){
    if (length(reference_colors) != length(reference)){
      stop("ERROR: Please make sure reference_colors is the same length as reference")
    }
  }
  
  if(!override_removeSmallRef) {
    globTest <- invisible(summix_calc(data = data, reference = reference,
                                      observed = observed, pi.start = pi.start))
    if(sum(globTest[1,5:ncol(globTest)] < 0.01) > 0) {
      toremove <- which(globTest[1,5:ncol(globTest)] < 0.01)
      reference <- reference[-toremove]
    }
    pi.start = NA
    sum_res <- summix_calc(data = data, reference = reference,
                           observed = observed)
  }
  sum_res <- summix_calc(data = data, reference = reference,
                         observed = observed, pi.start = pi.start)
  if(goodness.of.fit) {
    new_obj <- calc_scaledObj(data = data, observed = observed, reference = reference, pi.start= pi.start)
    sum_res[1] <- new_obj
  }
  end_time = Sys.time()
  ttime = end_time - start_time
  sum_res[3] <- ttime
  if(goodness.of.fit) {
    if(sum_res[1] >= 0.5 & sum_res[1] < 1.5) {
      print(paste0("CAUTION: Objective/1000SNPs = ", round(sum_res[1], 4),
                   " which is within the moderate fit range"))
    } else if (sum_res[1] >= 1.5) {
      print(paste0("WARNING: Objective/1000SNPs = ", round(sum_res[1], 4),
                   " which is above the poor fit threshold"))
    }
  } else {
    if(sum_res[1]/nrow(data)/1000 >= 0.5 & sum_res[1]/nrow(data)/1000 < 1.5) {
      print(paste0("CAUTION: Objective/1000SNPs = ", round(sum_res[1]/nrow(data)/1000, 4),
                   " which is within the moderate fit range"))
    } else if (sum_res[1]/nrow(data)/1000 >= 1.5) {
      print(paste0("WARNING: Objective/1000SNPs = ", round(sum_res[1]/nrow(data)/1000, 4),
                   " which is above the poor fit threshold"))
    }
  }
  if(network) {
    return(list(sum_res, summix_network(data = data, sum_res = sum_res, reference = reference, N_reference = N_reference, reference_colors=reference_colors)))
  }else{
    return(sum_res)
  }
}

summix_calc = function(data, reference, observed, pi.start=NA){
  if(!is(object = data, class2 = "data.frame")){
    stop("ERROR: data must be a data.frame as described in the vignette")
  }
  if(typeof(observed)!="character"){
    stop("ERROR: 'observed' must be a character string for the column name of the observed group in data")
  }
  if(!(observed %in% names(data))){
    stop("ERROR: 'observed' must be the column name of the observed group in data")
  }
  if(typeof(reference)!="character"){
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }
  if(all(reference %in% names(data))==FALSE){
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  }

  data <- as.data.frame(data)
  filteredNA <- length(which(is.na(data[[observed]]==TRUE)))
  observed.b  <- as.data.frame(data[which(is.na(data[[observed]])==FALSE), observed, drop = FALSE])
  refmatrix <- as.data.frame(data[which(is.na(data[[observed]])==FALSE), reference, drop = FALSE])

  if(!is.na(pi.start)[1]){
    if (is.numeric(pi.start)==FALSE){
      stop("ERROR: Please make sure pi.start is a numeric vector")
    }
    if (length(pi.start) != length(reference)){
      stop("ERROR: Please make sure pi.start is the same length as reference")
    }
    if (all(pi.start>0) == FALSE){
      stop("ERROR: Please make sure pi.start is a positive numeric vector")
    }
    if (sum(pi.start)!=1){
      stop("ERROR: Please make sure pi.start sums to one")
    }
    starting = pi.start
  } else{
    starting = rep( 1/ncol(refmatrix), ncol(refmatrix) )
  }
  fn.refmix = function(x){
    expected=x%*%t(as.matrix(refmatrix))
    minfunc = sum((expected - observed.b)**2)
    return(minfunc)
  }
  gr.refmix <- function(x){
    gradfunc = x%*%t(as.matrix(refmatrix)) - observed.b
    gradvec <- apply(2*refmatrix*t(gradfunc), 2, sum)
    return(gradvec)
  }
  heq.refmix = function(x){
    equality = sum(x)
    return(equality - 1)
  }

  start_time = Sys.time()
  
  # hin.refmix <- function(x){
  #   h = numeric(ncol(refmatrix))
  #   h=x
  #   return(h)
  # }
  #
  # S = suppressMessages( nloptr::slsqp(starting,
  #                                     fn = fn.refmix,
  #                                     gr = gr.refmix,
  #                                     hin = hin.refmix,
  #                                     heq = heq.refmix)
  # )

  # Use built-in lower bounds instead of manual hin constraints to avoid deprecation warnings
  S = nloptr::slsqp(starting,
                    fn = fn.refmix,
                    gr = gr.refmix,
                    lower = rep(0, ncol(refmatrix)),
                    heq = heq.refmix)

  end_time = Sys.time()
  ttime = end_time - start_time

  d <- data.frame(matrix(ncol = length(reference)+4, nrow = 1))
  colnames(d) <- c("goodness.of.fit", "iterations", "time",
                   "filtered", colnames(refmatrix))
  
  d[1] <- S$value
  d[2] <- S$iter
  d[3] <- ttime
  d[4] <- filteredNA
  d[5:(length(reference)+4)] <- round(S$par[1:length(reference)], 6)
  return(d)
}

calc_scaledObj <- function(data, reference, observed, pi.start) {
  start_time <- Sys.time()
  
  multiplier <- c(5, 1.5, 1)
  data$obs <- data[[observed]]
  
  data$maf_obs <- ifelse(data$obs > 0.5, 1-data$obs, data$obs)
  
  breaks <- c(0, 0.1, 0.3, 0.5)
  data$bin <- cut(data$maf_obs, breaks = breaks, include.lowest = T, right = FALSE,
                  labels = c(1, 2, 3))
  
  for(b in 1:3) {
    subdata <- data[data$bin == b,]
    subdata <- subdata[complete.cases(subdata), ]
    if(nrow(subdata) > 0) {
      res <- summix_calc(subdata, reference = reference, observed = observed, pi.start)
      res$obj_adj <-  res$goodness.of.fit/(nrow(subdata)/1000)
      res$nSNPs <- nrow(subdata)
      res$bin <- b
      if(b == 1) {
        sum_res <- res
      } else {
        sum_res <- rbind(sum_res, res)
      }
    } else {
      res <- summix_calc(data, reference = reference,
                         observed = observed, pi.start)
      res$obj_adj <- NA
      res$nSNPs <- 0
      res$bin <- b
      if(b == 1) {
        sum_res <- res
      } else {
        sum_res <- rbind(sum_res, res)
      }
    }
  }
  objective_scaled <- 0
  nonNA <- 0
  
  for(b in 1:3) {
    if(!is.na(sum_res[b, ]$obj_adj)) {
      
      nonNA <- nonNA + 1
      objective_scaled <- objective_scaled +
        (sum_res[b, ]$obj_adj*multiplier[b])
    }
  }
  objective_scaled <- objective_scaled/nonNA
  end_time <- Sys.time()
  return(objective_scaled)
}
