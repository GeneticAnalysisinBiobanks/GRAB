

#' Instruction of SPAsqr method
#'
#' SPAsqr is a smoothed quantile regression-based association test method for 
#' quantitative traits. It accounts for sample relatedness using the SPAGRM framework, 
#' performs association testing across multiple quantiles, and combines p-values 
#' using the Cauchy combination test.
#'
#' @return NULL
#'
#' @examples
#' # Step 1: fit null model and calculate joint distribution of genotypes
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultSPAsqr.txt")
#'
#' obj.SPAsqr <- GRAB.NullModel(
#'   QuantPheno ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   method = "SPAsqr",
#'   traitType = "quantitative",
#'   GenoFile = GenoFile,
#'   SparseGRMFile = SparseGRMFile,
#'   control = list(
#'     taus = c(0.05, 0.2, 0.5, 0.8, 0.95),
#'     h = 0
#'   )
#' )
#'
#' # Step 2: perform association tests
#' GRAB.Marker(obj.SPAsqr, GenoFile, OutputFile)
#'
#' head(data.table::fread(OutputFile))
#'
#' @details
#'
#' \strong{Additional Control Parameters for \code{GRAB.NullModel()}:}
#' \itemize{
#'   \item \code{taus} (numeric vector, default: c(0.05, 0.2, 0.5, 0.8, 0.95)): Quantiles 
#'     to examine for association testing. All values must be between 0 and 1. P-values 
#'     across quantiles are combined using Cauchy combination test.
#'   \item \code{h} (numeric, default: 0): Bandwidth parameter for smooth quantile regression. 
#'     If h = 0, bandwidth is automatically selected as IQR(y)/3.
#'   \item \code{frqFile} (character, default: NULL): Path to allele frequency file. If NULL,
#'     uses PlinkPrefix.frq for pairwise IBD calculation.
#'   \item \code{tempDir} (character, default: NULL): Directory for temporary files during
#'     pairwise IBD calculation. If NULL, uses tempdir().
#'   \item \code{maxSampleNums} (integer, default: 2500): Maximum number of subjects' genotypes
#'     to read for pairwise IBD analysis.
#'   \item \code{minMafIBD} (numeric, default: 0.01): Minimum MAF cutoff to select markers
#'     for pairwise IBD calculation.
#'   \item \code{rm.tempFile} (logical, default: FALSE): Whether to delete temporary files
#'     after pairwise IBD calculation.
#'   \item \code{MaxQuantile} (numeric, default: 0.75): Upper quantile for outlier detection.
#'     Must be greater than MinQuantile.
#'   \item \code{MinQuantile} (numeric, default: 0.25): Lower quantile for outlier detection.
#'     Must be less than MaxQuantile.
#'   \item \code{OutlierRatio} (numeric, default: 1.5): IQR multiplier for outlier cutoff 
#'     calculation. Must be greater than 0.
#'   \item \code{ControlOutlier} (logical, default: TRUE): Whether to automatically adjust 
#'     outlier ratio to keep outliers below 5\% of sample size.
#'   \item \code{MaxNuminFam} (integer, default: 5): Maximum family size for Chow-Liu tree 
#'     construction in related samples.
#'   \item \code{MAF_interval} (numeric vector, default: c(0.0001, 0.0005, 0.001, 0.005, 
#'     0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)): MAF breakpoints for genotype distribution 
#'     approximation in families.
#'   \item \code{sqr_tol} (numeric, default: 1e-7): Tolerance level for convergence in SQR fitting.
#' }
#'
#' \strong{Method-specific elements in the \code{SPAsqr_NULL_Model} object returned by \code{GRAB.NullModel()}:}
#' \itemize{
#'   \item All elements from \code{SPAGRM_NULL_Model} (see \code{\link{SPAGRM.NullModel}}).
#'   \item \code{Resid}: Numeric vector (or matrix for multiple quantiles) of residuals from 
#'     quantile regression model.
#'   \item \code{subjData}: Character vector of subject IDs.
#'   \item \code{N}: Number of subjects in analysis.
#' }
#'
#' \strong{Output file columns}:
#' \describe{
#'   \item{Marker}{Marker identifier (rsID or CHR:POS:REF:ALT).}
#'   \item{Info}{Marker information in format CHR:POS:REF:ALT.}
#'   \item{AltFreq}{Alternative allele frequency in the sample.}
#'   \item{AltCounts}{Total count of alternative alleles.}
#'   \item{MissingRate}{Proportion of missing genotypes.}
#'   \item{zScore}{Z-score from the score test.}
#'   \item{Pvalue}{P-value from the score test (Cauchy combined across quantiles).}
#'   \item{hwepval}{Hardy-Weinberg equilibrium p-value.}
#' }
#'
#' @references
#' Xu et al. (2025). SPAGRM: effectively controlling for sample relatedness in large-scale 
#' genome-wide association studies of longitudinal traits. Nature Communications, 16:1018.
#' 
#' Liu, Yaowu and Jun Xie (2020). Cauchy combination test: a powerful test with analytic p-value 
#' calculation under arbitrary dependency structures. Journal of the American Statistical 
#' Association, 115.529:393-402.
#'
GRAB.SPAsqr <- function() {
  .message("?GRAB.SPAsqr for instructions")
}


checkControl.NullModel.SPAsqr <- function(traitType, GenoFile, SparseGRMFile, control) {

  if (!traitType %in% c("quantitative")) {
    stop("For 'SPAsqr' method, only traitType of 'quantitative' is supported.")
  }

  if (is.null(GenoFile)) {
    stop("Argument 'GenoFile' is required for method 'SPAsqr'.")
  }
  if (!is.character(GenoFile) || length(GenoFile) != 1) {
    stop("Argument 'GenoFile' should be a character string (file path).")
  }
  if (!file.exists(GenoFile)) {
    stop("Cannot find GenoFile: ", GenoFile)
  }

  if (is.null(SparseGRMFile)) {
    stop("Argument 'SparseGRMFile' is required for method 'SPAsqr'.")
  }
  if (!is.character(SparseGRMFile) || length(SparseGRMFile) != 1) {
    stop("Argument 'SparseGRMFile' should be a character string (file path).")
  }
  if (!file.exists(SparseGRMFile)) {
    stop("Cannot find SparseGRMFile: ", SparseGRMFile)
  }

  default.control <- list(
    # Parameters for null model fitting
    taus = c(0.05,0.2,0.5,0.8,0.95),
    h = 0,

    # Parameters for getPairwiseIBD
    frqFile = NULL,
    tempDir = NULL,
    maxSampleNums = 2500,
    minMafIBD = 0.01,
    rm.tempFile = FALSE,

    # Parameters for SPAGRM.NullModel
    MaxQuantile = 0.75,
    MinQuantile = 0.25,
    OutlierRatio = 1.5,
    ControlOutlier = TRUE,
    MaxNuminFam = 5,
    MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
    sqr_tol=1e-7
  )
  control <- updateControl(control, default.control)

  # Validate taus
  if (!is.numeric(control$taus) || length(control$taus) == 0 ||
      any(control$taus <= 0) || any(control$taus >= 1)) {
    stop("'control$taus' should be a numeric vector with each element between 0 and 1.")
  }

  # Validate h
  if (!is.numeric(control$h) || length(control$h) != 1 || control$h < 0) {
    stop("'control$h' should be a single non-negative numeric value (default is 0 for automatic selection).")
  }

  # Validate frqFile if provided
  if (!is.null(control$frqFile)) {
    if (!is.character(control$frqFile) || length(control$frqFile) != 1) {
      stop("'control$frqFile' should be a character string (file path).")
    }
    if (!file.exists(control$frqFile)) {
      stop("Cannot find frqFile: ", control$frqFile)
    }
  }

  # Validate tempDir if provided
  if (!is.null(control$tempDir)) {
    if (!is.character(control$tempDir) || length(control$tempDir) != 1) {
      stop("'control$tempDir' should be a character string (directory path).")
    }
  }

  # Validate maxSampleNums
  if (!is.numeric(control$maxSampleNums) || length(control$maxSampleNums) != 1 || control$maxSampleNums <= 0) {
    stop("'control$maxSampleNums' should be a single numeric value greater than 0 (default is 2500).")
  }

  # Validate minMafIBD
  if (!is.numeric(control$minMafIBD) || length(control$minMafIBD) != 1 ||
      control$minMafIBD < 0 || control$minMafIBD > 0.5) {
    stop("'control$minMafIBD' should be a single numeric value between 0 and 0.5 (default is 0.01).")
  }

  # Validate rm.tempFile
  if (!is.logical(control$rm.tempFile) || length(control$rm.tempFile) != 1) {
    stop("'control$rm.tempFile' should be a single logical value (TRUE or FALSE).")
  }

  # Validate MaxQuantile and MinQuantile
  if (!is.numeric(control$MaxQuantile) || length(control$MaxQuantile) != 1 || 
      !is.numeric(control$MinQuantile) || length(control$MinQuantile) != 1 ||
      control$MaxQuantile <= 0 || control$MaxQuantile >= 1 || 
      control$MinQuantile <= 0 || control$MinQuantile >= 1) {
    stop("'control$MaxQuantile' and 'control$MinQuantile' should each be a single numeric value between 0 and 1.")
  }
  if (control$MaxQuantile <= control$MinQuantile) {
    stop("'control$MaxQuantile' (default is 0.75) should be larger than 'control$MinQuantile' (default is 0.25).")
  }
 
  # Validate OutlierRatio
  if (!is.numeric(control$OutlierRatio) || length(control$OutlierRatio) != 1 || control$OutlierRatio <= 0) {
    stop("'control$OutlierRatio' should be a single numeric value greater than 0 (default is 1.5).")
  }

  # Validate ControlOutlier
  if (!is.logical(control$ControlOutlier) || length(control$ControlOutlier) != 1) {
    stop("'control$ControlOutlier' should be a single logical value (TRUE or FALSE).")
  }

  # Validate MaxNuminFam
  if (!is.numeric(control$MaxNuminFam) || length(control$MaxNuminFam) != 1 || 
       control$MaxNuminFam %% 1 != 0 || control$MaxNuminFam < 1) {
    stop("'control$MaxNuminFam' should be an integer greater than or equal to 1.")
  }

  # Validate MAF_interval
  if (!is.numeric(control$MAF_interval) || any(control$MAF_interval < 0) || 
      any(control$MAF_interval > 0.5) || is.unsorted(control$MAF_interval)) {
    stop("'control$MAF_interval' should be a numeric vector in ascending order, with each element between 0 and 0.5.")
  }

  return(list(control = control, optionGRM = "SparseGRM"))
}


fitNullModel.SPAsqr <- function(
  response,
  designMat,
  subjData,
  control,
  GenoFile, 
  SparseGRMFile
) {
  # ========== Fit SQR null model ==========
  X <- designMat
  y <- response
  taus <- control$taus
  ntaus <- length(taus)
  h <- control$h
  if (h == 0) h <- IQR(y)/3

  ResidMat <- matrix(0, nrow = length(y), ncol = ntaus)
  for (i in seq_along(taus)) {
    current_tau <- taus[i]
    current_fit <- conquer::conquer(X, y, tau = current_tau, kernel = "Gaussian", h = h, tol = control$sqr_tol)
    current_resid <- as.numeric(y - current_fit$coeff[1] - X %*% current_fit$coeff[2:(ncol(X)+1)])
    ResidMat[, i] <- current_tau - pnorm((-current_resid)/h)
  }

  # ========= # Read SparseGRM and make ID set to be same as subjData ==========
  SparseGRM <- data.table::fread(SparseGRMFile)
  data.table::setnames(SparseGRM, c("ID1", "ID2", "Value"))
  data.table::set(SparseGRM, j = "ID1", value = as.character(SparseGRM$ID1))
  data.table::set(SparseGRM, j = "ID2", value = as.character(SparseGRM$ID2))
  SparseGRM <- SparseGRM[SparseGRM$ID1 %in% subjData & SparseGRM$ID2 %in% subjData, ]

  missSubjInGRM <- subjData[!subjData %in% unique(c(SparseGRM$ID1, SparseGRM$ID2))]
  num_missing <- length(missSubjInGRM)

  if (num_missing > 0) {
    # Append one line for each missing subject: ID ID 1
    new_rows <- data.table::data.table(
      ID1 = missSubjInGRM,
      ID2 = missSubjInGRM,
      Value = 1
    )
    SparseGRM <- rbind(SparseGRM, new_rows)

    # Prepare message (show up to 5 missing IDs)
    show_n <- min(5, num_missing)
    msg <- paste0(
      num_missing, " subjects do not have GRM info, impute as:\n",
      paste0(missSubjInGRM[1:show_n], " ", missSubjInGRM[1:show_n], " 1", collapse = "\n"),
      if (num_missing > show_n) "\n..." else ""
    )
    .message(msg)
  }

  # ========== SPAGRM Workflow ==========
  PairwiseIBD <- getPairwiseIBD(
    PlinkPrefix = sub("\\.bed$", "", GenoFile, ignore.case = TRUE),
    SparseGRMFile = SparseGRM,
    PairwiseIBDOutput = NULL,
    frqFile = control$frqFile,
    tempDir = control$tempDir,
    maxSampleNums = control$maxSampleNums,
    minMafIBD = control$minMafIBD,
    rm.tempFile = control$rm.tempFile
  )

  obj <- SPAGRM.NullModel.Multi(subjData, ResidMat, SparseGRM, PairwiseIBD, control)
  class(obj) <- "SPAsqr_NULL_Model"

  return(obj)
}


checkControl.Marker.SPAsqr <- function(control) {

  default.control <- list(
    zeta = 0,
    tol = 1e-5
  )
  control <- updateControl(control, default.control)

  return(control)
}


setMarker.SPAsqr <- function(objNull, control) {
  # Pass all tau data directly to C++ backend
  # C++ will handle indexing and processing for each tau
  
  setSPAsqrobjInCPP(
    t_taus = objNull$taus,                                # numeric vector: Quantiles (length ntaus)
    t_Resid_mat = objNull$Resid_mat,                      # numeric matrix: Residuals (N × ntaus)
    t_Resid_unrelated_outliers_lst = objNull$Resid.unrelated.outliers_lst,  # list (ntaus elements)
    t_sum_R_nonOutlier_vec = objNull$sum_R_nonOutlier_vec,     # numeric vector (length ntaus)
    t_R_GRM_R_nonOutlier_vec = objNull$R_GRM_R_nonOutlier_vec, # numeric vector (length ntaus)
    t_R_GRM_R_TwoSubjOutlier_vec = objNull$R_GRM_R_TwoSubjOutlier_vec,  # numeric vector (length ntaus)
    t_R_GRM_R_vec = objNull$R_GRM_R_vec,                  # numeric vector: Full R'*GRM*R (length ntaus)
    t_MAF_interval = objNull$MAF_interval,                # numeric vector: MAF intervals for binning
    t_TwoSubj_list_lst = objNull$TwoSubj_list_lst,                # list: Two-subject outlier pairs
    t_ThreeSubj_list_lst = objNull$ThreeSubj_list_lst,            # list: Three+ subject outlier families
    t_SPA_Cutoff = control$SPA_Cutoff,                    # numeric: P-value cutoff for SPA
    t_zeta = control$zeta,                                # numeric: SPA moment approximation parameter
    t_tol = control$tol                                   # numeric: Numerical tolerance for SPA
  )
}


mainMarker.SPAsqr <- function(genoType, genoIndex, objNull) {
  # Perform main marker analysis for SPAGRM method
  OutList <- mainMarkerInCPP(
    t_method = "SPAsqr",      # character: Statistical method name
    t_genoType = genoType,    # character: "PLINK" or "BGEN"
    t_genoIndex = genoIndex   # integer vector: Genotype indices to analyze
  )

  # Format results into output data frame
  obj.mainMarker <- data.frame(
    Marker = OutList$markerVec, # marker IDs
    Info = OutList$infoVec, # marker information: CHR:POS:REF:ALT
    AltFreq = OutList$altFreqVec, # alternative allele frequencies
    AltCounts = OutList$altCountsVec, # alternative allele counts
    MissingRate = OutList$missingRateVec, # missing data rate
    hwepval = OutList$hwepvalVec # Hardy-Weinberg equilibrium p-value
  )

  taus <- objNull$taus
  ntaus <- length(taus)

  z_df <- data.frame(matrix(OutList$zScore, ncol = ntaus, byrow = TRUE))
  colnames(z_df) <- paste0("Z_tau", taus)

  p_df <- data.frame(matrix(OutList$pvalVec, ncol = ntaus, byrow = TRUE))
  colnames(p_df) <- paste0("P_tau", taus)

  p_cct <- apply(p_df, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA)
    CCT(x)
  })
  obj.mainMarker <- cbind(obj.mainMarker, z_df, p_df, P_CCT = p_cct)

  return(obj.mainMarker)
}


SPAGRM.NullModel.Multi <- function(
  subjData, 
  ResidMat,
  SparseGRM,
  PairwiseIBD,
  control
) {
  taus <- control$taus
  ntaus <- length(taus)
  MaxQuantile <- control$MaxQuantile
  MinQuantile <- control$MinQuantile
  OutlierRatio <- control$OutlierRatio
  ControlOutlier <- control$ControlOutlier
  MaxNuminFam <- control$MaxNuminFam
  MAF_interval <- control$MAF_interval

  #### Identify outliers based on quantiles
  Quants <- apply(ResidMat, 2, function(x) quantile(x, probs = c(MinQuantile, MaxQuantile)))
  Ranges <- Quants[2, ] - Quants[1, ]
  cutoffLower <- Quants[1, ] - OutlierRatio * Ranges
  cutoffUpper <- Quants[2, ] + OutlierRatio * Ranges

  tooSmall <- sweep(ResidMat, 2, cutoffLower, "<")
  tooLarge <- sweep(ResidMat, 2, cutoffUpper, ">")
  Outlier <- tooSmall | tooLarge

  #### Pre-compute tau-independent graph structure
  edges <- t(SparseGRM[, c("ID1", "ID2")])
  graph_GRM <- igraph::make_graph(edges, directed = FALSE)
  graph_list_all <- igraph::decompose(graph_GRM)
  graph_vertex_names <- lapply(graph_list_all, function(g) igraph::V(g)$name)
  graph_length <- sapply(graph_list_all, length)

  graph_list_1 <- graph_list_all[graph_length == 1]
  SubjID.unrelated <- lapply(graph_list_1, igraph::get.vertex.attribute) %>% unlist(use.names = FALSE)
  graph_list <- graph_list_all[graph_length > 1]
  nGraph <- length(graph_list)

  if (nGraph == 0) {
    stop("No family found in SparseGRM. Please check SparseGRMFile.")
  }

  #### Step 1: Calculate contributions from UNRELATED subjects (vectorized across all taus)
  
  # Extract unrelated subject data
  unrelated_idx <- which(subjData %in% SubjID.unrelated)
  ResidMat_unrelated <- ResidMat[unrelated_idx, , drop = FALSE]
  Outlier_unrelated <- Outlier[unrelated_idx, , drop = FALSE]
  
  # Filter GRM for unrelated subjects only
  SparseGRM_unrelated <- SparseGRM %>% filter(ID1 %in% SubjID.unrelated)
  pos1_idx <- match(SparseGRM_unrelated$ID1, subjData)
  pos2_idx <- match(SparseGRM_unrelated$ID2, subjData)
  ResidMat_pos1 <- ResidMat[pos1_idx, , drop = FALSE]
  ResidMat_pos2 <- ResidMat[pos2_idx, , drop = FALSE]
  
  # Calculate R'*GRM*R for unrelated subjects (all taus simultaneously)
  cov_matrix <- abs(SparseGRM_unrelated$Value) * ResidMat_pos1 * ResidMat_pos2
  R_GRM_R_vec <- colSums(cov_matrix)
  
  # Calculate R'*GRM*R for non-outlier unrelated subjects (all taus simultaneously)
  id1_in_unrelated_idx <- match(SparseGRM_unrelated$ID1, subjData[unrelated_idx])
  mask_matrix <- !Outlier_unrelated[id1_in_unrelated_idx, , drop = FALSE]
  cov_masked <- cov_matrix * mask_matrix
  R_GRM_R_nonOutlier_vec <- colSums(cov_masked)
  
  # Calculate sum of non-outlier residuals for unrelated subjects (all taus simultaneously)
  sum_R_nonOutlier_vec <- colSums(ResidMat_unrelated * (!Outlier_unrelated))
  
  # Extract outlier residuals for unrelated subjects (per tau)
  Resid.unrelated.outliers_lst <- lapply(seq_along(taus), function(i) {
    ResidMat_unrelated[Outlier_unrelated[, i], i]
  })
  
  # Initialize output structures for RELATED subjects
  R_GRM_R_TwoSubjOutlier_vec <- numeric(ntaus)
  TwoSubj_list_lst <- lapply(1:ntaus, function(x) list())
  ThreeSubj_list_lst <- lapply(1:ntaus, function(x) list())

  arr.index <- list()
  for (n in seq_len(MaxNuminFam)) {
    temp <- c()
    for (i_idx in seq_len(n)) {
      indexString <- rep("c(1, 1, 1)", n)
      indexString[i_idx] <- "0:2"
      indexString <- paste0(indexString, collapse = "%o%")
      cmd <- paste0("temp = c(temp, list(arr.index", i_idx, "=", indexString, "))")
      eval(parse(text = cmd))
    }
    arr.index[[n]] <- temp
  }

  #### Step 2: First tau loop - Build graph_list_updated_lst and add RELATED family contributions

  # graph_list_updated_lst is a list of length ntaus, where each element i contains:
  # A list of igraph objects representing families with outliers for that specific tau value
  # These families have been decomposed/adjusted so that no family exceeds MaxNuminFam members
  graph_list_updated_lst <- vector("list", ntaus)
  for (i in seq_along(taus)) {

    ResidMat_df <- data.frame(
      SubjID = subjData,
      Resid = ResidMat[, i],
      Outlier = Outlier[, i]
    )

    # Calculate tau-specific covariances for edge weighting
    SparseGRM1 <- SparseGRM
    SparseGRM1$pos1 <- ResidMat_df$Resid[match(SparseGRM$ID1, ResidMat_df$SubjID)]
    SparseGRM1$pos2 <- ResidMat_df$Resid[match(SparseGRM$ID2, ResidMat_df$SubjID)]
    SparseGRM1 <- SparseGRM1 %>% mutate(Cov = abs(Value * pos1 * pos2))

    # Initialize for this tau
    graph_list_updated <- list()
    index.outlier <- 1

    .message("Processing %d family groups for tau %g", nGraph, taus[i])

    for (i_fam in seq_len(nGraph)) {
      if (i_fam %% 1000 == 0) {
        .message("Processing family group %d/%d", i_fam, nGraph)
      }

      comp1 <- graph_list[[i_fam]]
      comp3 <- graph_vertex_names[[which(graph_length > 1)[i_fam]]]
      pos1 <- match(comp3, subjData)
      outlierInFam <- any(ResidMat_df$Outlier[pos1])

      # Add related family variance to unrelated contribution
      block_GRM <- make.block.GRM(comp1, SparseGRM)
      R_GRM_R.temp <- as.numeric(t(ResidMat_df$Resid[pos1]) %*% block_GRM %*% ResidMat_df$Resid[pos1])
      R_GRM_R_vec[i] <- R_GRM_R_vec[i] + R_GRM_R.temp

      if (!outlierInFam) {
        # Add non-outlier family contributions to unrelated contributions
        sum_R_nonOutlier_vec[i] <- sum_R_nonOutlier_vec[i] + sum(ResidMat_df$Resid[pos1])
        R_GRM_R_nonOutlier_vec[i] <- R_GRM_R_nonOutlier_vec[i] + R_GRM_R.temp
        next
      }

      # Use pre-computed graph length instead of igraph::vcount
      vcount <- graph_length[which(graph_length > 1)[i_fam]]

      if (vcount <= MaxNuminFam) {
        graph_list_updated[[index.outlier]] <- comp1
        index.outlier <- index.outlier + 1
        next
      }

      # Step 1: remove the edges until the largest family size is <= MaxNuminFam, default is 5.
      comp1.temp <- comp1
      tempGRM1 <- SparseGRM1 %>%
        filter(ID1 %in% comp3 | ID2 %in% comp3) %>%
        arrange(Cov)
      for (j in seq_len(nrow(tempGRM1))) {
        # Remove edge and calculate vertex counts for new graph components
        edgesToRemove <- paste0(tempGRM1$ID1[j], "|", tempGRM1$ID2[j])
        comp1.temp <- igraph::delete.edges(comp1.temp, edgesToRemove)
        # Get vertex count for each component after edge removal
        vcount <- igraph::decompose(comp1.temp) %>% sapply(igraph::vcount)
        if (max(vcount) <= MaxNuminFam) {
          break
        }
      }

      # Step 2: add the (removed) edges while keeping the largest family size <= MaxNuminFam, default is 5.
      tempGRM1 <- tempGRM1[seq_len(j), ] %>% arrange(desc(Cov))
      comp1 <- comp1.temp
      for (k in seq_len(nrow(tempGRM1))) {
        # Add edge and calculate vertex counts for new graph components
        edgesToAdd <- c(tempGRM1$ID1[k], tempGRM1$ID2[k])
        comp1.temp <- igraph::add.edges(comp1, edgesToAdd)

        # Get vertex count for each component after edge addition
        vcount <- igraph::decompose(comp1.temp) %>% sapply(igraph::vcount)

        if (max(vcount) <= MaxNuminFam) {
          comp1 <- comp1.temp
        }
      }

      comp1 <- igraph::decompose(comp1)

      for (k in seq_len(length(comp1))) {
        comp11 <- comp1[[k]]
        comp13 <- igraph::V(comp11)$name

        pos2 <- match(comp13, subjData)
        outlierInFam <- any(ResidMat_df$Outlier[pos2])

        block_GRM <- make.block.GRM(comp11, SparseGRM)
        R_GRM_R.temp <- as.numeric(t(ResidMat_df$Resid[pos2]) %*% block_GRM %*% ResidMat_df$Resid[pos2])

        if (!outlierInFam) {
          sum_R_nonOutlier_vec[i] <- sum_R_nonOutlier_vec[i] + sum(ResidMat_df$Resid[pos2])
          R_GRM_R_nonOutlier_vec[i] <- R_GRM_R_nonOutlier_vec[i] + R_GRM_R.temp
        } else {
          graph_list_updated[[index.outlier]] <- comp11
          index.outlier <- index.outlier + 1
        }
      }
    }

    graph_list_updated_lst[[i]] <- graph_list_updated
  }

  #### Step 3: Build CLT cache for union of all outlier families (across all taus)

  # Collect all unique outlier families across all taus
  all_unique_families <- list()
  family_id_map <- new.env(hash = TRUE)
  
  for (i in seq_along(taus)) {
    graph_list_updated <- graph_list_updated_lst[[i]]
    if (length(graph_list_updated) == 0) next
    
    for (comp1 in graph_list_updated) {
      comp3 <- igraph::V(comp1)$name
      family_key <- paste(sort(comp3), collapse = "_")
      
      if (!exists(family_key, envir = family_id_map)) {
        assign(family_key, length(all_unique_families) + 1, envir = family_id_map)
        all_unique_families[[length(all_unique_families) + 1]] <- comp1
      }
    }
  }

  # Pre-compute CLT for all unique families
  CLT_cache <- vector("list", length(all_unique_families))
  
  if (length(all_unique_families) > 0) {
    .message("Computing Chow-Liu trees for cache %d outlier unique families", length(all_unique_families))
    
    for (fam_idx in seq_along(all_unique_families)) {

      if (fam_idx %% 100 == 0) {
        .message("Processing family group %d/%d", fam_idx, length(all_unique_families))
      }
      comp1 <- all_unique_families[[fam_idx]]
      comp3 <- igraph::V(comp1)$name
      n1 <- length(comp3)
      
      if (n1 <= 2) next  # CLT only needed for n >= 3
      
      tempIBD <- PairwiseIBD %>% filter(ID1 %in% comp3 & ID2 %in% comp3)
      
      CLT_cache[[fam_idx]] <- chow.liu.tree(
        N = n1,
        IBD = tempIBD,
        IDs = comp3,
        MAF_interval = MAF_interval
      )
    }
  }

  #### Step 4: Fill ThreeSubj_list_lst using cached CLTs
  for (i in seq_along(taus)) {

    ResidMat_df <- data.frame(
      SubjID = subjData,
      Resid = ResidMat[, i],
      Outlier = Outlier[, i]
    )

    graph_list_updated <- graph_list_updated_lst[[i]]
    R_GRM_R_TwoSubjOutlier <- 0
    TwoSubj_list <- ThreeSubj_list <- list()

    if (length(graph_list_updated) != 0) {
      .message("Building standardized scores for %d outlier families for tau %g", length(graph_list_updated), taus[i])

      # build chou-liu-tree.
      n.outliers <- length(graph_list_updated)
      ## The below values are only used in chou.liu.tree
      TwofamID.index <- ThreefamID.index <- 0
      for (index.outlier in seq_len(n.outliers)) {

        comp1 <- graph_list_updated[[index.outlier]]
        comp3 <- igraph::V(comp1)$name
        n1 <- length(comp3)
        pos3 <- match(comp3, subjData)

        Resid.temp <- ResidMat_df$Resid[pos3]

        if (n1 == 1) {
          Resid.unrelated.outliers_lst[[i]] <- c(Resid.unrelated.outliers_lst[[i]], Resid.temp)
          next
        }

        block_GRM <- make.block.GRM(comp1, SparseGRM)

        tempIBD <- PairwiseIBD %>% filter(ID1 %in% comp3 & ID2 %in% comp3)

        if (n1 == 2) {
          TwofamID.index <- TwofamID.index + 1

          R_GRM_R_TwoSubjOutlier.temp <- as.numeric(t(Resid.temp) %*% block_GRM %*% Resid.temp)
          R_GRM_R_TwoSubjOutlier <- R_GRM_R_TwoSubjOutlier + R_GRM_R_TwoSubjOutlier.temp

          Rho.temp <- tempIBD$pa + 0.5 * tempIBD$pb
          midterm <- sqrt(Rho.temp^2 - tempIBD$pa)

          TwoSubj_list[[TwofamID.index]] <- list(
            Resid = Resid.temp,
            Rho = c(Rho.temp + midterm, Rho.temp - midterm)
          )
          next
        }

        ThreefamID.index <- ThreefamID.index + 1

        # Use pre-computed CLT from cache
        family_key <- paste(sort(comp3), collapse = "_")
        fam_cache_idx <- get(family_key, envir = family_id_map)
        CLT <- CLT_cache[[fam_cache_idx]]

        # Calculate standardized score using array operations and mapply
        stand.S.temp <- array(
          rowSums(mapply(function(x, y) x * y, arr.index[[n1]], Resid.temp)),
          rep(3, n1)
        )

        ThreeSubj_list[[ThreefamID.index]] <- list(
          CLT = CLT,
          stand.S = c(stand.S.temp)
        )
      }
    }

    R_GRM_R_TwoSubjOutlier_vec[i] <- R_GRM_R_TwoSubjOutlier
    TwoSubj_list_lst[[i]] <- TwoSubj_list
    ThreeSubj_list_lst[[i]] <- ThreeSubj_list
  }

  # Return single list with vectors/matrices for all taus
  obj <- list(
    taus = taus,
    Resid_mat = ResidMat,
    subjData = subjData,
    N = length(subjData),
    R_GRM_R_vec = R_GRM_R_vec,
    R_GRM_R_TwoSubjOutlier_vec = R_GRM_R_TwoSubjOutlier_vec,
    sum_R_nonOutlier_vec = sum_R_nonOutlier_vec,
    R_GRM_R_nonOutlier_vec = R_GRM_R_nonOutlier_vec,
    Resid.unrelated.outliers_lst = Resid.unrelated.outliers_lst,
    TwoSubj_list_lst = TwoSubj_list_lst,
    ThreeSubj_list_lst = ThreeSubj_list_lst,
    MAF_interval = MAF_interval
  )

  return(obj)
}
