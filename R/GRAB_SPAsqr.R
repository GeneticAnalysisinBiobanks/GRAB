#' Instructions for SPAsqr method (Smoothed Quantile Regression with Saddlepoint Approximation)
#'
#' SPAsqr is a smoothed quantile regression-based association test method for 
#' quantitative traits. It accounts for sample relatedness using the SPAGRM framework, 
#' performs association testing across multiple quantiles, and combines p-values 
#' using the Cauchy combination test (CCT).
#' 
#' @return NULL
#'
#' @examples
#' # Step 1: Fit null model and prepare genotype distribution cache
#' # See ?getSparseGRM for details on generating a sparse GRM.
#' # See ?getPairwiseIBD for details on computing pairwise IBD estimates.
#'
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' PairwiseIBDFile <- system.file("extdata", "PairwiseIBD.txt", package = "GRAB")
#' 
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' 
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultSPAsqr.txt")
#'
#' obj.SPAsqr <- GRAB.NullModel(
#'   QuantPheno ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   method = "SPAsqr",
#'   traitType = "quantitative",
#'   SparseGRMFile = SparseGRMFile,
#'   control = list(taus = c(0.2, 0.5, 0.8)),
#'   PairwiseIBDFile = PairwiseIBDFile
#' )
#'
#' # Step 2: Perform single-marker association tests
#' GRAB.Marker(obj.SPAsqr, GenoFile, OutputFile)
#'
#' # View results
#' head(data.table::fread(OutputFile))
#'
#' @details
#'
#' \strong{Usage with \code{GRAB.NullModel()}:}
#' 
#' Set \code{method = "SPAsqr"} and \code{traitType = "quantitative"}. Two additional required arguments:
#' \itemize{
#'   \item \code{SparseGRMFile}: Path to sparse GRM file (whitespace-delimited: ID1 ID2 Value).
#'   \item \code{PairwiseIBDFile}: Path to pairwise IBD file (see \code{?\link{getPairwiseIBD}}).
#' }
#'
#' \strong{Control parameters for \code{GRAB.NullModel(control = list(...))}:}
#' \itemize{
#'   \item \code{taus} (numeric vector, default: c(0.05, 0.2, 0.5, 0.8, 0.95)): Quantiles 
#'     to examine for association testing. All values must be between 0 and 1 (exclusive). 
#'     P-values across quantiles are combined using Cauchy combination test.
#'   \item \code{smooth} (logical, default: TRUE): Whether to use smoothed quantile regression.
#'     If FALSE, nonsmooth quantile regression is used; \code{h} and \code{sqr_tol} are ignored.
#'   \item \code{h} (numeric, default: 0): Bandwidth parameter for smooth quantile regression. 
#'     If h = 0, bandwidth is automatically selected as IQR(y)/3.
#'   \item \code{MaxNuminFam} (integer, default: 5): Maximum family size for Chow-Liu tree 
#'     construction in related samples. Larger families are decomposed into smaller components.
#'   \item \code{MAF_interval} (numeric vector, default: c(0.0001, 0.0005, 0.001, 0.005, 
#'     0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)): MAF breakpoints for genotype distribution 
#'     approximation in families.
#'   \item \code{sqr_tol} (numeric, default: 1e-7): Tolerance level for convergence in smoothed quantile regression fitting.
#' }
#'
#' \strong{Control parameters for \code{GRAB.Marker(control = list(...))}:}
#' \itemize{
#'   \item \code{SPA_Cutoff} (numeric, default: 2): Z-score cutoff for applying saddlepoint approximation.
#'   \item \code{zeta} (numeric, default: 0): SPA moment approximation parameter.
#'   \item \code{tol} (numeric, default: 1e-5): Numerical tolerance for SPA calculations.
#'   \item See \code{?\link{GRAB.Marker}} for additional quality control parameters (MAF, MAC, missing rate).
#' }
#'
#' \strong{Elements in the \code{SPAsqr_NULL_Model} object returned by \code{GRAB.NullModel()}:}
#' \itemize{
#'   \item \code{taus}: Numeric vector of quantiles analyzed.
#'   \item \code{Resid_mat}: Numeric matrix of residuals from smoothed quantile regression (N × ntaus).
#'   \item \code{subjData}: Character vector of subject IDs included in analysis.
#'   \item \code{N}: Number of subjects in analysis (integer).
#'   \item \code{R_GRM_R_vec}: Numeric vector of R'*GRM*R values for each quantile (length ntaus).
#'   \item \code{R_GRM_R_TwoSubjOutlier_vec}: Numeric vector of R'*GRM*R contributions from two-subject outlier pairs (length ntaus).
#'   \item \code{sum_R_nonOutlier_vec}: Numeric vector of summed residuals for non-outlier subjects (length ntaus).
#'   \item \code{R_GRM_R_nonOutlier_vec}: Numeric vector of R'*GRM*R contributions from non-outlier subjects (length ntaus).
#'   \item \code{Resid.unrelated.outliers_lst}: List of residual vectors for unrelated outliers (ntaus elements).
#'   \item \code{TwoSubj_list_lst}: List of two-subject outlier family information (ntaus elements).
#'   \item \code{CLT_union_lst}: List of Chow-Liu tree structures for related families (shared across all quantiles).
#'   \item \code{ThreeSubj_family_idx_lst}: List of family indices into CLT cache for 3+ member families (ntaus elements).
#'   \item \code{ThreeSubj_stand_S_lst}: List of standardized score arrays for 3+ member families (ntaus elements).
#'   \item \code{MAF_interval}: Numeric vector of MAF breakpoints used for genotype approximation.
#'   \item \code{Call}: Original function call.
#'   \item \code{sessionInfo}: R session and package information.
#'   \item \code{time}: Analysis completion timestamp (character).
#'   \item \code{control}: List of control parameters used in fitting.
#' }
#'
#' \strong{Output file columns from \code{GRAB.Marker()}:}
#' \describe{
#'   \item{Marker}{Marker identifier (rsID or CHR:POS:REF:ALT).}
#'   \item{Info}{Marker information in format CHR:POS:REF:ALT.}
#'   \item{AltFreq}{Alternative allele frequency in the sample.}
#'   \item{AltCounts}{Total count of alternative alleles (integer).}
#'   \item{MissingRate}{Proportion of missing genotypes.}
#'   \item{hwepval}{Hardy-Weinberg equilibrium p-value.}
#'   \item{Z_tau[tau]}{Z-score from score test at quantile [tau] (one column per tau).}
#'   \item{P_tau[tau]}{P-value from score test at quantile [tau] (one column per tau).}
#'   \item{P_CCT}{Cauchy combination test p-value aggregating across all quantiles.}
#' }
#'
#' @references
#' Heng et al. (in prep). Discovering Heterogeneous Associations with Smoothed Quantile Regression and Saddle Point Approximation.
#'
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


checkControl.NullModel.SPAsqr <- function(traitType, GenoFile, SparseGRMFile, control, ...) {

  PairwiseIBDFile <- list(...)$PairwiseIBDFile

  if (!traitType %in% c("quantitative")) {
    stop("For 'SPAsqr' method, only traitType of 'quantitative' is supported.")
  }

  if (!is.null(GenoFile)) {
    warning("Argument 'GenoFile' is ignored for method 'SPACox'.")
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
  
  if (is.null(PairwiseIBDFile)) {
    stop("Argument 'PairwiseIBDFile' is required for method 'SPAsqr'.")
  }
  if (!is.character(PairwiseIBDFile) || length(PairwiseIBDFile) != 1) {
    stop("Argument 'PairwiseIBDFile' should be a character string (file path).")
  }
  if (!file.exists(PairwiseIBDFile)) {
    stop("Cannot find PairwiseIBDFile: ", PairwiseIBDFile)
  }

  default.control <- list(
    # Parameters for null model fitting
    taus = c(0.05,0.2,0.5,0.8,0.95),
    smooth=TRUE,
    h = 0,
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

  # Validate smooth
  if (!is.logical(control$smooth) || length(control$smooth) != 1) {
    stop("'control$smooth' should be a single logical value (TRUE or FALSE).")
  }

  # Validate h
  if (control$smooth) {
    if (!is.numeric(control$h) || length(control$h) != 1 || control$h < 0) {
      stop("'control$h' should be a single non-negative numeric value (default is 0 for automatic selection).")
    }
    # Validate sqr_tol
    if (!is.numeric(control$sqr_tol) || length(control$sqr_tol) != 1 || control$sqr_tol <= 0) {
      stop("'control$sqr_tol' should be a single positive numeric value.")
    }
  } else {
    control$h <- NULL  # Not used if smooth = FALSE
    control$sqr_tol <- NULL  # Not used if smooth = FALSE
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
  SparseGRMFile,
  ...
) {
  PairwiseIBDFile <- list(...)$PairwiseIBDFile

  # ========== Fit null model ==========
  y <- response
  taus <- control$taus
  ntaus <- length(taus)
  ResidMat <- matrix(0, nrow = length(y), ncol = ntaus)

  if (control$smooth) {
    # Smoothed quantile regression
    X <- designMat
    h <- control$h
    if (h == 0) h <- IQR(y)/3

    for (i in seq_along(taus)) {
      current_tau <- taus[i]
      current_fit <- conquer::conquer(X, y, tau = current_tau, kernel = "Gaussian", h = h, tol = control$sqr_tol)
      current_resid <- as.numeric(y - current_fit$coeff[1] - X %*% current_fit$coeff[2:(ncol(X)+1)])
      ResidMat[, i] <- current_tau - pnorm((-current_resid)/h)
    }
  } else {
    # Nonsmooth quantile regression
    X <- cbind(1, designMat)

    for (i in seq_along(taus)) {
      current_tau <- taus[i]
      current_fit <- quantreg::rq.fit(X, y, tau = current_tau)
      ResidMat[, i] <- current_tau - ifelse(current_fit$residuals < 0, 1, 0)
    }
  }

  # ========= # Read SparseGRM and make ID set to be same as subjData ==========
  SparseGRM <- data.table::fread(SparseGRMFile)
  data.table::setnames(SparseGRM, c("ID1", "ID2", "Value"))
  data.table::set(SparseGRM, j = "ID1", value = as.character(SparseGRM$ID1))
  data.table::set(SparseGRM, j = "ID2", value = as.character(SparseGRM$ID2))
  SparseGRM <- SparseGRM[SparseGRM$ID1 %in% subjData & SparseGRM$ID2 %in% subjData, ]

  # ========== SPAGRM Workflow ==========
  PairwiseIBD <- data.table::fread(PairwiseIBDFile)
  data.table::set(PairwiseIBD, j = "ID1", value = as.character(PairwiseIBD$ID1))
  data.table::set(PairwiseIBD, j = "ID2", value = as.character(PairwiseIBD$ID2))
  PairwiseIBD <- PairwiseIBD[PairwiseIBD$ID1 %in% subjData & PairwiseIBD$ID2 %in% subjData, ]

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
    t_CLT_union_lst = objNull$CLT_union_lst,                      # list: Shared CLT cache (union across all taus)
    t_ThreeSubj_family_idx_lst = objNull$ThreeSubj_family_idx_lst, # list: Family indices per tau
    t_ThreeSubj_stand_S_lst = objNull$ThreeSubj_stand_S_lst,       # list: stand.S values per tau
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
  MaxNuminFam <- control$MaxNuminFam
  MAF_interval <- control$MAF_interval
  
  #### Identify outliers based on quantiles
  tooSmall <- sweep(ResidMat, 2, -0.8, "<")
  tooLarge <- sweep(ResidMat, 2, 0.8, ">")
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
    .message("No family found in SparseGRM. Treating all individuals as unrelated.")
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
  
  .message("Pre-computing block GRMs and decomposing outlier families")
  
  # Pre-compute block GRMs for ALL ORIGINAL families (before any decomposition)
  # These are needed for correct variance calculation
  block_GRM_cache_original <- vector("list", nGraph)
  .message("Pre-computing block GRMs for %d original families", nGraph)
  
  for (i_fam in seq_len(nGraph)) {
    comp1 <- graph_list[[i_fam]]
    block_GRM_cache_original[[i_fam]] <- make.block.GRM(comp1, SparseGRM)
  }
  
  # Identify which families have outliers (union across all taus)
  families_with_outliers <- logical(nGraph)
  
  .message("Identifying families with outliers across all taus")
  
  for (i_fam in seq_len(nGraph)) {
    comp3 <- graph_vertex_names[[which(graph_length > 1)[i_fam]]]
    pos1 <- match(comp3, subjData)
    
    # Check if family has outliers in ANY tau
    has_outlier <- any(Outlier[pos1, ])
    families_with_outliers[i_fam] <- has_outlier
  }
  
  n_outlier_families <- sum(families_with_outliers)
  .message("Found %d families with outliers (out of %d total)", n_outlier_families, nGraph)
  
  # Decompose large outlier families and compute block GRMs for components
  decomposed_families_lst <- vector("list", nGraph)
  decomposed_block_GRM_lst <- vector("list", nGraph)
  
  families_to_decompose <- which(families_with_outliers & graph_length[graph_length > 1] > MaxNuminFam)
  
  if (length(families_to_decompose) > 0) {
    .message("Decomposing %d large outlier families (size > MaxNuminFam=%d)", 
             length(families_to_decompose), MaxNuminFam)
    
    for (idx in seq_along(families_to_decompose)) {
      i_fam <- families_to_decompose[idx]
      
      if (idx %% 100 == 0) {
        .message("Decomposing family %d/%d", idx, length(families_to_decompose))
      }
      
      comp1 <- graph_list[[i_fam]]
      comp3 <- graph_vertex_names[[which(graph_length > 1)[i_fam]]]
      
      # Decompose using kinship values (tau-independent)
      comp1.temp <- comp1
      tempGRM1 <- SparseGRM %>%
        filter(ID1 %in% comp3 | ID2 %in% comp3) %>%
        mutate(Weight = abs(Value)) %>%
        arrange(Weight)
      
      # Remove edges until largest component <= MaxNuminFam
      for (j in seq_len(nrow(tempGRM1))) {
        edgesToRemove <- paste0(tempGRM1$ID1[j], "|", tempGRM1$ID2[j])
        comp1.temp <- igraph::delete_edges(comp1.temp, edgesToRemove)
        vcount_decomp <- igraph::decompose(comp1.temp) %>% sapply(igraph::vcount)
        if (max(vcount_decomp) <= MaxNuminFam) {
          break
        }
      }
      
      # Add back edges (strongest first) while maintaining size constraint
      tempGRM1 <- tempGRM1[seq_len(j), ] %>% arrange(desc(Weight))
      comp1 <- comp1.temp
      for (k in seq_len(nrow(tempGRM1))) {
        edgesToAdd <- c(tempGRM1$ID1[k], tempGRM1$ID2[k])
        comp1.temp <- igraph::add_edges(comp1, edgesToAdd)
        vcount_decomp <- igraph::decompose(comp1.temp) %>% sapply(igraph::vcount)
        if (max(vcount_decomp) <= MaxNuminFam) {
          comp1 <- comp1.temp
        }
      }
      
      # Decompose into separate components
      comp1_decomposed <- igraph::decompose(comp1)
      
      # Pre-compute block GRM for each decomposed component
      block_GRMs_decomposed <- lapply(comp1_decomposed, function(comp11) {
        make.block.GRM(comp11, SparseGRM)
      })
      
      decomposed_families_lst[[i_fam]] <- comp1_decomposed
      decomposed_block_GRM_lst[[i_fam]] <- block_GRMs_decomposed
    }
    
    .message("Family decomposition complete")
  }
  
  .message("Computing tau-specific variance components")
  
  graph_list_updated_lst <- vector("list", ntaus)
  
  for (i in seq_along(taus)) {
    
    ResidMat_df <- data.frame(
      SubjID = subjData,
      Resid = ResidMat[, i],
      Outlier = Outlier[, i]
    )
    
    # Initialize for this tau
    graph_list_updated <- list()
    index.outlier <- 1
    
    .message("Processing %d families for tau %g", nGraph, taus[i])
    
    for (i_fam in seq_len(nGraph)) {
      if (i_fam %% 1000 == 0) {
        .message("Processing family %d/%d", i_fam, nGraph)
      }
      
      # Get ORIGINAL family info for variance calculation
      comp3_original <- graph_vertex_names[[which(graph_length > 1)[i_fam]]]
      pos_original <- match(comp3_original, subjData)
      
      # ALWAYS use ORIGINAL block GRM for variance calculation (correct total variance)
      block_GRM_original <- block_GRM_cache_original[[i_fam]]
      R_GRM_R_original <- as.numeric(t(ResidMat_df$Resid[pos_original]) %*% block_GRM_original %*% ResidMat_df$Resid[pos_original])
      R_GRM_R_vec[i] <- R_GRM_R_vec[i] + R_GRM_R_original
      
      # Check for outliers in this tau
      outlierInFam_thisTau <- any(ResidMat_df$Outlier[pos_original])
      
      if (!outlierInFam_thisTau) {
        # No outliers in this tau - add to non-outlier contributions
        sum_R_nonOutlier_vec[i] <- sum_R_nonOutlier_vec[i] + sum(ResidMat_df$Resid[pos_original])
        R_GRM_R_nonOutlier_vec[i] <- R_GRM_R_nonOutlier_vec[i] + R_GRM_R_original
        next
      }
      
      # Family has outliers in this tau
      vcount <- graph_length[which(graph_length > 1)[i_fam]]
      
      if (vcount <= MaxNuminFam) {
        # Family is small enough - use original family for outlier processing
        comp1 <- graph_list[[i_fam]]
        graph_list_updated[[index.outlier]] <- comp1
        index.outlier <- index.outlier + 1
        next
      }
      
      # Family is large - use decomposed components for outlier processing
      decomposed_families <- decomposed_families_lst[[i_fam]]
      decomposed_block_GRMs <- decomposed_block_GRM_lst[[i_fam]]
      
      for (k in seq_along(decomposed_families)) {
        comp11 <- decomposed_families[[k]]
        comp13 <- igraph::V(comp11)$name
        pos2 <- match(comp13, subjData)
        
        # Check for outliers in this component for this tau
        outlierInComponent <- any(ResidMat_df$Outlier[pos2])
        
        if (!outlierInComponent) {
          # Component has no outliers in this tau
          # Add to non-outlier contributions using DECOMPOSED block GRM
          block_GRM_component <- decomposed_block_GRMs[[k]]
          R_GRM_R_component <- as.numeric(t(ResidMat_df$Resid[pos2]) %*% block_GRM_component %*% ResidMat_df$Resid[pos2])
          
          sum_R_nonOutlier_vec[i] <- sum_R_nonOutlier_vec[i] + sum(ResidMat_df$Resid[pos2])
          R_GRM_R_nonOutlier_vec[i] <- R_GRM_R_nonOutlier_vec[i] + R_GRM_R_component
        } else {
          # Component has outliers - save for later processing
          graph_list_updated[[index.outlier]] <- comp11
          index.outlier <- index.outlier + 1
        }
      }
    }
    
    graph_list_updated_lst[[i]] <- graph_list_updated
  }
  
  .message("Family structures updated and variance components computed")
  
  # Build lookup for block GRMs of outlier families (for Step 4)
  block_GRM_lookup <- new.env(hash = TRUE)
  
  # Add small outlier families (original, not decomposed)
  for (i_fam in seq_len(nGraph)) {
    vcount <- graph_length[which(graph_length > 1)[i_fam]]
    has_outlier <- families_with_outliers[i_fam]
    
    if (has_outlier && vcount <= MaxNuminFam) {
      comp3 <- graph_vertex_names[[which(graph_length > 1)[i_fam]]]
      family_key <- paste(sort(comp3), collapse = "_")
      assign(family_key, block_GRM_cache_original[[i_fam]], envir = block_GRM_lookup)
    }
  }
  
  # Add decomposed outlier family components
  for (i_fam in seq_len(nGraph)) {
    if (!is.null(decomposed_families_lst[[i_fam]])) {
      decomposed_families <- decomposed_families_lst[[i_fam]]
      decomposed_block_GRMs <- decomposed_block_GRM_lst[[i_fam]]
      
      for (k in seq_along(decomposed_families)) {
        comp11 <- decomposed_families[[k]]
        comp13 <- igraph::V(comp11)$name
        family_key <- paste(sort(comp13), collapse = "_")
        assign(family_key, decomposed_block_GRMs[[k]], envir = block_GRM_lookup)
      }
    }
  }
  
  .message("Built block GRM lookup with %d entries", length(ls(envir = block_GRM_lookup)))
  
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
  
  #### Step 4: Calculate stand.S for each tau, stored separately from CLT
  # Initialize storage for 3+ member families
  ThreeSubj_family_idx_lst <- vector("list", ntaus)  # Family indices into CLT_cache
  ThreeSubj_stand_S_lst <- vector("list", ntaus)     # stand.S values per tau
  
  for (i in seq_along(taus)) {
    
    ResidMat_df <- data.frame(
      SubjID = subjData,
      Resid = ResidMat[, i],
      Outlier = Outlier[, i]
    )
    
    graph_list_updated <- graph_list_updated_lst[[i]]
    R_GRM_R_TwoSubjOutlier <- 0
    TwoSubj_list <- list()
    
    # Storage for 3+ member families (this tau)
    ThreeSubj_family_idx <- integer(0)
    ThreeSubj_stand_S <- list()
    
    if (length(graph_list_updated) != 0) {
      .message("Building standardized scores for %d outlier families for tau %g", length(graph_list_updated), taus[i])
      
      n.outliers <- length(graph_list_updated)
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
        
        family_key <- paste(sort(comp3), collapse = "_")
        block_GRM <- get(family_key, envir = block_GRM_lookup)
        
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
        
        # For families with 3+ members: store family index and stand.S separately
        ThreefamID.index <- ThreefamID.index + 1
        
        # Get family index in CLT_cache
        fam_cache_idx <- get(family_key, envir = family_id_map)
        
        # Calculate standardized score using array operations and mapply
        stand.S.temp <- array(
          rowSums(mapply(function(x, y) x * y, arr.index[[n1]], Resid.temp)),
          rep(3, n1)
        )
        
        # Store family index and stand.S separately
        ThreeSubj_family_idx[ThreefamID.index] <- fam_cache_idx
        ThreeSubj_stand_S[[ThreefamID.index]] <- c(stand.S.temp)
      }
    }
    
    R_GRM_R_TwoSubjOutlier_vec[i] <- R_GRM_R_TwoSubjOutlier
    TwoSubj_list_lst[[i]] <- TwoSubj_list
    ThreeSubj_family_idx_lst[[i]] <- ThreeSubj_family_idx
    ThreeSubj_stand_S_lst[[i]] <- ThreeSubj_stand_S
  }
  
  # Return single list with CLT and stand.S stored separately
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
    # New separated structure for 3+ member families:
    CLT_union_lst = CLT_cache,                           # Shared CLT cache (union across all taus)
    ThreeSubj_family_idx_lst = ThreeSubj_family_idx_lst,  # Family indices per tau
    ThreeSubj_stand_S_lst = ThreeSubj_stand_S_lst,        # stand.S values per tau
    MAF_interval = MAF_interval
  )
  
  return(obj)
}

# =======================================================================

# When sample size is large and many traits are analyzed or loco is used,
# it is more efficient to separate the GRM cache building (Step 1a) from
# the null model fitting (Step 1). In this file, SPAsqr.Step1a builds 
# the GRM cache, and SPAsqr.Step1b computes the null model object using 
# the pre-computed GRM cache. The returned object is for GRAB.Marker input.


# Step 1a: Build GRM Cache (ResidMat-independent)
SPAsqr.Step1a <- function(
  subjData,
  SparseGRMFile,
  PairwiseIBDFile,
  MaxNuminFam = 5,
  MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5),
  nsubprocess = 1
) {

  .message("Loading Sparse GRM and Pairwise IBD data")

  SparseGRM <- data.table::fread(SparseGRMFile)
  data.table::setnames(SparseGRM, c("ID1", "ID2", "Value"))
  data.table::set(SparseGRM, j = "ID1", value = as.character(SparseGRM$ID1))
  data.table::set(SparseGRM, j = "ID2", value = as.character(SparseGRM$ID2))
  SparseGRM <- SparseGRM[SparseGRM$ID1 %in% subjData & SparseGRM$ID2 %in% subjData, ]

  PairwiseIBD <- data.table::fread(PairwiseIBDFile)
  data.table::set(PairwiseIBD, j = "ID1", value = as.character(PairwiseIBD$ID1))
  data.table::set(PairwiseIBD, j = "ID2", value = as.character(PairwiseIBD$ID2))
  PairwiseIBD <- PairwiseIBD[PairwiseIBD$ID1 %in% subjData & PairwiseIBD$ID2 %in% subjData, ]
  

  .message("Building GRM cache (ResidMat-independent)")
  
  edges <- t(SparseGRM[, c("ID1", "ID2")])
  graph_GRM <- igraph::make_graph(edges, directed = FALSE)
  graph_list_all <- igraph::decompose(graph_GRM)
  graph_vertex_names <- parallel::mclapply(graph_list_all, function(g) igraph::V(g)$name, mc.cores = nsubprocess)
  graph_length <- sapply(graph_list_all, length)
  
  graph_list_1 <- graph_list_all[graph_length == 1]
  SubjID.unrelated <- unlist(lapply(graph_list_1, igraph::get.vertex.attribute), use.names = FALSE)
  graph_list <- graph_list_all[graph_length > 1]
  nGraph <- length(graph_list)
  
  if (nGraph == 0) {
    stop("No related subjects found in the SparseGRM. Please check the input GRM file.")
  }
  
  #### Build arr.index structure
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
  
  #### Pre-compute block GRMs for ALL ORIGINAL families (before any decomposition)
  .message("Pre-computing block GRMs for %d original families", nGraph)
  
  # Parallelize block GRM computation
  block_GRM_cache_original <- parallel::mclapply(seq_len(nGraph), function(i_fam) {
    if (i_fam %% 10000 == 0) {
      .message("Computing block GRM %d/%d", i_fam, nGraph)
    }
    comp1 <- graph_list[[i_fam]]
    make.block.GRM(comp1, SparseGRM)
  }, mc.cores = nsubprocess)
  
  #### Identify families that need decomposition (all families > MaxNuminFam)
  families_needing_decomposition <- which(graph_length[graph_length > 1] > MaxNuminFam)
  
  #### Decompose large families and compute block GRMs for components
  decomposed_families_lst <- vector("list", nGraph)
  decomposed_block_GRM_lst <- vector("list", nGraph)
  
  if (length(families_needing_decomposition) > 0) {
    .message("Pre-decomposing %d large families (size > MaxNuminFam=%d)", 
             length(families_needing_decomposition), MaxNuminFam)
    
    # Parallelize decomposition
    decomposition_results <- parallel::mclapply(seq_along(families_needing_decomposition), function(idx) {
      i_fam <- families_needing_decomposition[idx]
      
      if (idx %% 1000 == 0) {
        .message("Decomposing family %d/%d", idx, length(families_needing_decomposition))
      }
      
      comp1 <- graph_list[[i_fam]]
      comp3 <- graph_vertex_names[[which(graph_length > 1)[i_fam]]]
      
      # Decompose using kinship values (tau-independent)
      comp1.temp <- comp1
      # Use data.table operations with proper allocation
      temp_filter <- SparseGRM$ID1 %in% comp3 | SparseGRM$ID2 %in% comp3
      tempGRM1 <- data.table::copy(SparseGRM[temp_filter, ])
      data.table::setalloccol(tempGRM1)
      data.table::set(tempGRM1, j = "Weight", value = abs(tempGRM1$Value))
      data.table::setorder(tempGRM1, Weight)
      
      # Remove edges until largest component <= MaxNuminFam
      for (j in seq_len(nrow(tempGRM1))) {
        edgesToRemove <- paste0(tempGRM1$ID1[j], "|", tempGRM1$ID2[j])
        comp1.temp <- igraph::delete_edges(comp1.temp, edgesToRemove)
        vcount_decomp <- sapply(igraph::decompose(comp1.temp), igraph::vcount)
        if (max(vcount_decomp) <= MaxNuminFam) {
          break
        }
      }
      
      # Add back edges (strongest first) while maintaining size constraint
      tempGRM1 <- tempGRM1[seq_len(j), ]
      data.table::setorder(tempGRM1, -Weight)
      comp1 <- comp1.temp
      for (k in seq_len(nrow(tempGRM1))) {
        edgesToAdd <- c(tempGRM1$ID1[k], tempGRM1$ID2[k])
        comp1.temp <- igraph::add_edges(comp1, edgesToAdd)
        vcount_decomp <- sapply(igraph::decompose(comp1.temp), igraph::vcount)
        if (max(vcount_decomp) <= MaxNuminFam) {
          comp1 <- comp1.temp
        }
      }
      
      # Decompose into separate components
      comp1_decomposed <- igraph::decompose(comp1)
      
      # Pre-compute block GRM for each decomposed component (parallel)
      block_GRMs_decomposed <- parallel::mclapply(comp1_decomposed, function(comp11) {
        make.block.GRM(comp11, SparseGRM)
      }, mc.cores = nsubprocess)
      
      list(families = comp1_decomposed, block_GRMs = block_GRMs_decomposed, i_fam = i_fam)
    }, mc.cores = nsubprocess)
    
    # Unpack results into original structure
    for (result in decomposition_results) {
      i_fam <- result$i_fam
      decomposed_families_lst[[i_fam]] <- result$families
      decomposed_block_GRM_lst[[i_fam]] <- result$block_GRMs
    }
    
    # Memory cleanup: remove temporary decomposition results
    rm(decomposition_results)
    gc()
    
    .message("Family decomposition complete")
  }
  
  #### Combine operations: Build block GRM lookup and collect unique families
  # This includes both small families and decomposed large families
  .message("Building block GRM lookup and collecting unique families")
  
  # Pre-compute which families are related (avoid repeated indexing)
  related_idx <- which(graph_length > 1)
  vcounts <- graph_length[related_idx]
  
  # Parallelize family processing to collect data
  family_data <- parallel::mclapply(seq_len(nGraph), function(i_fam) {
    vcount <- vcounts[i_fam]
    result <- list()
    
    if (vcount <= MaxNuminFam) {
      # Small family - process original
      comp3 <- graph_vertex_names[[related_idx[i_fam]]]
      family_key <- paste(sort(comp3), collapse = "_")
      
      result[[1]] <- list(
        family_key = family_key,
        block_GRM = block_GRM_cache_original[[i_fam]],
        comp_graph = graph_list[[i_fam]]
      )
    } else if (!is.null(decomposed_families_lst[[i_fam]])) {
      # Large family - process decomposed components
      decomposed_families <- decomposed_families_lst[[i_fam]]
      decomposed_block_GRMs <- decomposed_block_GRM_lst[[i_fam]]
      
      for (k in seq_along(decomposed_families)) {
        comp11 <- decomposed_families[[k]]
        comp13 <- igraph::V(comp11)$name
        family_key <- paste(sort(comp13), collapse = "_")
        
        result[[length(result) + 1]] <- list(
          family_key = family_key,
          block_GRM = decomposed_block_GRMs[[k]],
          comp_graph = comp11
        )
      }
    }
    return(result)
  }, mc.cores = nsubprocess)
  
  # Flatten results and build lookup tables sequentially (must be sequential for uniqueness)
  # Store references to original data instead of duplicating block GRMs
  # Use named lists instead of environments for proper serialization with save/saveRDS
  family_key_to_source <- list()  # family_key -> list(source, index)
  all_unique_families <- list()
  family_id_map <- list()
  
  current_ifam <- 1
  for (i in seq_len(nGraph)) {
    vcount <- vcounts[i]
    
    if (vcount <= MaxNuminFam) {
      # Small family - reference original cache
      comp3 <- graph_vertex_names[[related_idx[i]]]
      family_key <- paste(sort(comp3), collapse = "_")
      
      # Store reference to original cache
      family_key_to_source[[family_key]] <- list(source = "original", index = i)
      
      # Add to unique families if not seen before
      if (is.null(family_id_map[[family_key]])) {
        all_unique_families[[length(all_unique_families) + 1]] <- graph_list[[i]]
        family_id_map[[family_key]] <- length(all_unique_families)
      }
    } else if (!is.null(decomposed_families_lst[[i]])) {
      # Large family - reference decomposed cache
      decomposed_families <- decomposed_families_lst[[i]]
      
      for (k in seq_along(decomposed_families)) {
        comp11 <- decomposed_families[[k]]
        comp13 <- igraph::V(comp11)$name
        family_key <- paste(sort(comp13), collapse = "_")
        
        # Store reference to decomposed cache
        family_key_to_source[[family_key]] <- list(source = "decomposed", i_fam = i, k = k)
        
        # Add to unique families if not seen before
        if (is.null(family_id_map[[family_key]])) {
          all_unique_families[[length(all_unique_families) + 1]] <- comp11
          family_id_map[[family_key]] <- length(all_unique_families)
        }
      }
    }
  }

  .message("Built block GRM reference map with %d entries and %d unique families", 
           length(family_key_to_source), length(all_unique_families))
  
  #### Pre-compute CLT for all unique families (only for families with 3+ members)
  n_unique <- length(all_unique_families)
  CLT_cache <- vector("list", n_unique)
  
  if (n_unique > 0) {
    .message("Computing Chow-Liu trees for %d unique families", n_unique)
    
    # Progress reporting thresholds
    progress_step <- max(1, floor(n_unique / 20))  # Report every 5%
    
    # Parallelize CLT computation
    CLT_cache <- parallel::mclapply(seq_len(n_unique), function(fam_idx) {
      if (fam_idx %% 1000 == 0) {
        .message("Processing family %d/%d", fam_idx, n_unique)
      }

      comp1 <- all_unique_families[[fam_idx]]
      comp3 <- igraph::V(comp1)$name
      n1 <- length(comp3)
      
      if (n1 <= 2) return(NULL)  # CLT only needed for n >= 3
      
      # Use base R subsetting for reliability
      tempIBD <- PairwiseIBD[PairwiseIBD$ID1 %in% comp3 & PairwiseIBD$ID2 %in% comp3, ]
      
      chow.liu.tree(
        N = n1,
        IBD = tempIBD,
        IDs = comp3,
        MAF_interval = MAF_interval
      )
    }, mc.cores = nsubprocess)
  }
  
  .message("GRM cache built successfully")
  
  # Estimated memory usage for 400k subjects, 62k families:
  # - SparseGRM: 2GB (depending on sparsity)
  # - block_GRM_cache_original: 5GB (62k matrices, avg 3-5 subjects each)
  # - CLT_cache: 3GB (depends on family sizes)
  # - Total peak: 15GB
  
  # Return pre-computed structures
  GRMCache <- list(
    graph_list = graph_list,
    graph_vertex_names = graph_vertex_names,
    graph_length = graph_length,
    SubjID.unrelated = SubjID.unrelated,
    nGraph = nGraph,
    block_GRM_cache_original = block_GRM_cache_original,
    decomposed_families_lst = decomposed_families_lst,
    decomposed_block_GRM_lst = decomposed_block_GRM_lst,
    family_key_to_source = family_key_to_source,
    family_id_map = family_id_map,
    CLT_cache = CLT_cache,
    arr.index = arr.index,
    SparseGRM = SparseGRM,
    PairwiseIBD = PairwiseIBD,
    MaxNuminFam = MaxNuminFam,
    MAF_interval = MAF_interval
  )
  
  return(GRMCache)
}


# Step 1b: Fit Null Model using pre-computed GRM Cache
SPAsqr.Step1b <- function(
  y,
  X,
  subjData,
  GRMCache,
  taus = c(0.05, 0.2, 0.5, 0.8, 0.95),
  smooth = TRUE,
  h = 0,
  sqr_tol = 1e-7
) {
  
  # ========== Fit null model ==========
  .message("Fitting null model and computing residuals for taus: %s", paste(taus, collapse = ", "))

  y[is.na(y)] <- median(y, na.rm = TRUE)
  ntaus <- length(taus)
  ResidMat <- matrix(0, nrow = length(y), ncol = ntaus)

  if (smooth) {
    # Smoothed quantile regression
    if (h == 0) h <- IQR(y)/3

    for (i in seq_along(taus)) {
      current_tau <- taus[i]
      current_fit <- conquer::conquer(X, y, tau = current_tau, kernel = "Gaussian", h = h, tol = sqr_tol)
      current_resid <- as.numeric(y - current_fit$coeff[1] - X %*% current_fit$coeff[2:(ncol(X)+1)])
      ResidMat[, i] <- current_tau - pnorm((-current_resid)/h)
    }
  } else {
    # Nonsmooth quantile regression
    X <- cbind(1, X)

    for (i in seq_along(taus)) {
      current_tau <- taus[i]
      current_fit <- quantreg::rq.fit(X, y, tau = current_tau)
      ResidMat[, i] <- current_tau - ifelse(current_fit$residuals < 0, 1, 0)
    }
  }

  # ========== Compute variance components using GRM cache ==========
  .message("Computing variance components using pre-built cache")
  
  # Extract from cache
  graph_list <- GRMCache$graph_list
  graph_vertex_names <- GRMCache$graph_vertex_names
  graph_length <- GRMCache$graph_length
  SubjID.unrelated <- GRMCache$SubjID.unrelated
  nGraph <- GRMCache$nGraph
  block_GRM_cache_original <- GRMCache$block_GRM_cache_original
  decomposed_families_lst <- GRMCache$decomposed_families_lst
  decomposed_block_GRM_lst <- GRMCache$decomposed_block_GRM_lst
  family_key_to_source <- GRMCache$family_key_to_source
  family_id_map <- GRMCache$family_id_map
  CLT_cache <- GRMCache$CLT_cache
  arr.index <- GRMCache$arr.index
  SparseGRM <- GRMCache$SparseGRM
  PairwiseIBD <- GRMCache$PairwiseIBD
  MaxNuminFam <- GRMCache$MaxNuminFam
  MAF_interval <- GRMCache$MAF_interval
  
  # Helper function to retrieve block GRM from original sources (no duplication)
  get_block_GRM <- function(family_key) {
    ref <- family_key_to_source[[family_key]]
    if (ref$source == "original") {
      block_GRM_cache_original[[ref$index]]
    } else {
      decomposed_block_GRM_lst[[ref$i_fam]][[ref$k]]
    }
  }
  
  #### Identify outliers based on quantiles
  tooSmall <- sweep(ResidMat, 2, -0.8, "<")
  tooLarge <- sweep(ResidMat, 2, 0.8, ">")
  Outlier <- tooSmall | tooLarge
  
  #### Step 1: Calculate contributions from UNRELATED subjects (vectorized across all taus)
  
  # Extract unrelated subject data
  unrelated_idx <- which(subjData %in% SubjID.unrelated)
  ResidMat_unrelated <- ResidMat[unrelated_idx, , drop = FALSE]
  Outlier_unrelated <- Outlier[unrelated_idx, , drop = FALSE]
  
  # Filter GRM for unrelated subjects only
  SparseGRM_unrelated <- SparseGRM[SparseGRM$ID1 %in% SubjID.unrelated, ]
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
  
  # Extract outlier residuals for unrelated subjects (per tau) - parallel
  Resid.unrelated.outliers_lst <- lapply(seq_along(taus), function(i) {
    ResidMat_unrelated[Outlier_unrelated[, i], i]
  })
  
  # Initialize output structures for RELATED subjects
  R_GRM_R_TwoSubjOutlier_vec <- numeric(ntaus)
  TwoSubj_list_lst <- lapply(1:ntaus, function(x) list())
  
  #### Step 2: First tau loop - Build graph_list_updated_lst and add RELATED family contributions
  
  .message("Identifying outlier families and computing variance components")
  
  graph_list_updated_lst <- vector("list", ntaus)
  
  # Pre-compute related family indices to avoid repeated calculations
  related_idx <- which(graph_length > 1)
  vcounts <- graph_length[related_idx]
  
  for (i in seq_along(taus)) {
    .message("Processing %d families for tau %g (%d/%d)", nGraph, taus[i], i, ntaus)
    
    # Pre-allocate data structures
    ResidMat_col <- ResidMat[, i]
    Outlier_col <- Outlier[, i]
    
    # Initialize for this tau
    graph_list_updated <- vector("list", nGraph)  # Pre-allocate max size
    index.outlier <- 1
    
    for (i_fam in seq_len(nGraph)) {
      # Get ORIGINAL family info for variance calculation
      comp3_original <- graph_vertex_names[[related_idx[i_fam]]]
      pos_original <- match(comp3_original, subjData)
      
      # ALWAYS use ORIGINAL block GRM for variance calculation (correct total variance)
      block_GRM_original <- block_GRM_cache_original[[i_fam]]
      R_GRM_R_original <- as.numeric(crossprod(ResidMat_col[pos_original], block_GRM_original %*% ResidMat_col[pos_original]))
      R_GRM_R_vec[i] <- R_GRM_R_vec[i] + R_GRM_R_original
      
      # Check for outliers in this tau
      outlierInFam_thisTau <- any(Outlier_col[pos_original])
      
      if (!outlierInFam_thisTau) {
        # No outliers in this tau - add to non-outlier contributions
        sum_R_nonOutlier_vec[i] <- sum_R_nonOutlier_vec[i] + sum(ResidMat_col[pos_original])
        R_GRM_R_nonOutlier_vec[i] <- R_GRM_R_nonOutlier_vec[i] + R_GRM_R_original
        next
      }
      
      # Family has outliers in this tau
      vcount <- vcounts[i_fam]
      
      if (vcount <= MaxNuminFam) {
        # Family is small enough - use original family for outlier processing
        comp1 <- graph_list[[i_fam]]
        graph_list_updated[[index.outlier]] <- comp1
        index.outlier <- index.outlier + 1
        next
      }
      
      # Family is large - use decomposed components for outlier processing
      decomposed_families <- decomposed_families_lst[[i_fam]]
      decomposed_block_GRMs <- decomposed_block_GRM_lst[[i_fam]]
      
      for (k in seq_along(decomposed_families)) {
        comp11 <- decomposed_families[[k]]
        comp13 <- igraph::V(comp11)$name
        pos2 <- match(comp13, subjData)
        
        # Check for outliers in this component for this tau
        outlierInComponent <- any(Outlier_col[pos2])
        
        if (!outlierInComponent) {
          # Component has no outliers in this tau
          # Add to non-outlier contributions using DECOMPOSED block GRM
          block_GRM_component <- decomposed_block_GRMs[[k]]
          R_GRM_R_component <- as.numeric(crossprod(ResidMat_col[pos2], block_GRM_component %*% ResidMat_col[pos2]))
          
          sum_R_nonOutlier_vec[i] <- sum_R_nonOutlier_vec[i] + sum(ResidMat_col[pos2])
          R_GRM_R_nonOutlier_vec[i] <- R_GRM_R_nonOutlier_vec[i] + R_GRM_R_component
        } else {
          # Component has outliers - save for later processing
          graph_list_updated[[index.outlier]] <- comp11
          index.outlier <- index.outlier + 1
        }
      }
    }
    
    # Trim to actual size
    graph_list_updated <- graph_list_updated[seq_len(index.outlier - 1)]
    graph_list_updated_lst[[i]] <- graph_list_updated
  }
  
  .message("Family structures updated and variance components computed")
  
  #### Build filtered CLT cache for outlier families only (to match fitNullModel.SPAsqr)
  
  # Collect all unique outlier families across all taus
  all_outlier_families <- list()
  outlier_family_id_map <- list()
  
  for (i in seq_along(taus)) {
    graph_list_updated <- graph_list_updated_lst[[i]]
    if (length(graph_list_updated) == 0) next
    
    for (comp1 in graph_list_updated) {
      comp3 <- igraph::V(comp1)$name
      family_key <- paste(sort(comp3), collapse = "_")
      
      if (is.null(outlier_family_id_map[[family_key]])) {
        # Map to the new filtered index
        new_idx <- length(all_outlier_families) + 1
        outlier_family_id_map[[family_key]] <- new_idx
        
        # Get the original index from full cache
        original_idx <- family_id_map[[family_key]]
        
        # Copy CLT from full cache to filtered cache
        all_outlier_families[[new_idx]] <- CLT_cache[[original_idx]]
      }
    }
  }
  
  .message("Extracted %d outlier families from %d total families", length(all_outlier_families), length(CLT_cache))
  
  #### Step 3: Calculate stand.S for each tau, stored separately from CLT
  # Initialize storage for 3+ member families
  ThreeSubj_family_idx_lst <- vector("list", ntaus)  # Family indices into filtered CLT_cache
  ThreeSubj_stand_S_lst <- vector("list", ntaus)     # stand.S values per tau
  
  for (i in seq_along(taus)) {
    
    # Pre-allocate column vectors
    ResidMat_col <- ResidMat[, i]
    Outlier_col <- Outlier[, i]
    
    graph_list_updated <- graph_list_updated_lst[[i]]
    R_GRM_R_TwoSubjOutlier <- 0
    TwoSubj_list <- list()
    
    # Storage for 3+ member families (this tau)
    ThreeSubj_family_idx <- integer(0)
    ThreeSubj_stand_S <- list()
    
    if (length(graph_list_updated) != 0) {
      if (i %% 5 == 1 || i == length(taus)) {
        .message("Building standardized scores for %d outlier families for tau %g (%d/%d)", 
                 length(graph_list_updated), taus[i], i, ntaus)
      }
      
      n.outliers <- length(graph_list_updated)
      TwofamID.index <- ThreefamID.index <- 0
      
      for (index.outlier in seq_len(n.outliers)) {
        
        comp1 <- graph_list_updated[[index.outlier]]
        comp3 <- igraph::V(comp1)$name
        n1 <- length(comp3)
        pos3 <- match(comp3, subjData)
        
        Resid.temp <- ResidMat_col[pos3]
        
        if (n1 == 1) {
          Resid.unrelated.outliers_lst[[i]] <- c(Resid.unrelated.outliers_lst[[i]], Resid.temp)
          next
        }
        
        family_key <- paste(sort(comp3), collapse = "_")
        block_GRM <- get_block_GRM(family_key)
        
        tempIBD <- PairwiseIBD[PairwiseIBD$ID1 %in% comp3 & PairwiseIBD$ID2 %in% comp3, ]
        
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
        
        # For families with 3+ members: store family index and stand.S separately
        ThreefamID.index <- ThreefamID.index + 1
        
        # Get family index in filtered outlier CLT_cache
        fam_cache_idx <- outlier_family_id_map[[family_key]]
        
        # Calculate standardized score using array operations and mapply
        stand.S.temp <- array(
          rowSums(mapply(function(x, y) x * y, arr.index[[n1]], Resid.temp)),
          rep(3, n1)
        )
        
        # Store family index and stand.S separately
        ThreeSubj_family_idx[ThreefamID.index] <- fam_cache_idx
        ThreeSubj_stand_S[[ThreefamID.index]] <- c(stand.S.temp)
      }
    }
    
    R_GRM_R_TwoSubjOutlier_vec[i] <- R_GRM_R_TwoSubjOutlier
    TwoSubj_list_lst[[i]] <- TwoSubj_list
    ThreeSubj_family_idx_lst[[i]] <- ThreeSubj_family_idx
    ThreeSubj_stand_S_lst[[i]] <- ThreeSubj_stand_S
  }
  
  .message("Variance computation complete")
  
  # Return single list with filtered CLT (outliers only) and stand.S stored separately
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
    CLT_union_lst = all_outlier_families,
    ThreeSubj_family_idx_lst = ThreeSubj_family_idx_lst,
    ThreeSubj_stand_S_lst = ThreeSubj_stand_S_lst,
    MAF_interval = MAF_interval
  )

  class(obj) <- "SPAsqr_NULL_Model"

  return(obj)
}
