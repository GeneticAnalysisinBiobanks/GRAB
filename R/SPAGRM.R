## ------------------------------------------------------------------------------
## SPAGRM.R
##
## Functions:
##   GRAB.SPAGRM                  : Print brief method information.
##   SPAGRM.NullModel             : Fit SPAGRM null model and residual handling.
##   checkControl.SPAGRM.NullModel: Validate/populate null-model controls.
##   checkControl.Marker.SPAGRM   : Validate marker-level controls.
##   setMarker.SPAGRM             : Initialize marker-level analysis objects.
##   mainMarker.SPAGRM            : Run marker-level SPAGRM tests.
##   make.block.GRM               : Construct GRM blocks for subgraphs/families.
##   chow.liu.tree                : Build Chow–Liu trees for family structures.
## ------------------------------------------------------------------------------

#' Instruction of SPAGRM method
#'
#' SPAGRM is a scalable and accurate framework for retrospective association tests. 
#' It treats genetic loci as random vectors and uses a precise approximation of their 
#' joint distribution. This approach enables SPAGRM to handle any type of complex trait, 
#' including longitudinal and unbalanced phenotypes. SPAGRM extends SPACox to support 
#' sample relatedness.
#'
#' @return NULL
#'
#' @examples
#' ResidMatFile <- system.file("extdata", "ResidMat.txt", package = "GRAB")
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' PairwiseIBDFile <- system.file("extdata", "PairwiseIBD.txt", package = "GRAB")
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultSPAGRM.txt")
#'
#' # Step 2a: pre-calculate genotype distributions
#' obj.SPAGRM <- SPAGRM.NullModel(
#'   ResidMatFile = ResidMatFile,
#'   SparseGRMFile = SparseGRMFile,
#'   PairwiseIBDFile = PairwiseIBDFile,
#'   control = list(ControlOutlier = FALSE)
#' )
#'
#' # Step 2b: perform association tests
#' GRAB.Marker(obj.SPAGRM, GenoFile, OutputFile)
#'
#' head(data.table::fread(OutputFile))
#'
#' @details
#' See \code{\link{SPAGRM.NullModel}} for detailed instructions
#' on preparing a SPAGRM_NULL_Model object required for GRAB.Marker().
#'
#' \strong{Additional Control Parameters for GRAB.Marker()}:
#' \itemize{
#'   \item \code{zeta} (numeric, default: 0): SPA moment approximation parameter.
#'   \item \code{tol} (numeric, default: 1e-5): Numerical tolerance for SPA convergence.
#' }
#' 
#' \strong{Marker-level results} (\code{OutputFile}) columns:
#' \describe{
#'   \item{Marker}{Marker identifier (rsID or CHR:POS:REF:ALT).}
#'   \item{Info}{Marker information in format CHR:POS:REF:ALT.}
#'   \item{AltFreq}{Alternative allele frequency in the sample.}
#'   \item{AltCounts}{Total count of alternative alleles.}
#'   \item{MissingRate}{Proportion of missing genotypes.}
#'   \item{zScore}{Z-score from the score test.}
#'   \item{Pvalue}{P-value from the score test.}
#'   \item{hwepval}{Hardy-Weinberg equilibrium p-value.}
#' }
#'
#' @references
#' Xu et al. (2025). SPAGRM: effectively controlling for sample relatedness in large-scale 
#' genome-wide association studies of longitudinal traits. \doi{10.1038/s41467-025-56669-1}
#'
GRAB.SPAGRM <- function() {
  .message("?SPAGRM for instructions")
}


#' Fit SPAGRM null model from residuals and relatedness inputs
#'
#' Builds the SPAGRM null model object using subject residuals, sparse GRM,
#' and pairwise IBD estimates, detecting residual outliers and constructing
#' family-level graph structures for downstream saddlepoint marker tests.
#'
#' @param ResidMatFile Data frame or file path with columns \code{SubjID, Resid}.
#' @param SparseGRMFile File path to sparse GRM (tab-delimited: ID1, ID2, Value).
#' @param PairwiseIBDFile File path to pairwise IBD table (ID1, ID2, pa, pb, pc).
#' @param control List of options controlling outlier handling and family
#'   decomposition (see \code{checkControl.SPAGRM.NullModel}).
#'
#' @return A list of class \code{"SPAGRM_NULL_Model"} with elements:
#'   \describe{
#'     \item{Resid}{Numeric vector of residuals used in analysis.}
#'     \item{subjData}{Character vector of subject IDs (length = N).}
#'     \item{N}{Number of subjects.}
#'     \item{Resid.unrelated.outliers}{Residuals of unrelated outlier subjects.}
#'     \item{R_GRM_R}{Sum of quadratic form Resid' * GRM * Resid for all subjects.}
#'     \item{R_GRM_R_TwoSubjOutlier}{Aggregate contribution from two-subject outlier families.}
#'     \item{sum_R_nonOutlier}{Sum of residuals for non-outlier unrelated subjects.}
#'     \item{R_GRM_R_nonOutlier}{Quadratic form contribution for non-outlier unrelated subjects.}
#'     \item{TwoSubj_list}{List with per two-member family residual/Rho info.}
#'     \item{ThreeSubj_list}{List with Chow–Liu tree structures and standardized scores for larger families.}
#'     \item{MAF_interval}{Vector of MAF breakpoints used in tree construction.}
#'   }
#'
SPAGRM.NullModel <- function(
  ResidMatFile, # two columns: column 1 is subjID, column 2 is Resid
  SparseGRMFile, # a path of SparseGRMFile get from getSparseGRM() function.
  PairwiseIBDFile, # a path of PairwiseIBDFile get from getPairwiseIBD() function.
  control = list(
    MaxQuantile = 0.75,
    MinQuantile = 0.25,
    OutlierRatio = 1.5,
    ControlOutlier = TRUE,
    MaxNuminFam = 5,
    MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  )
) {
  # Read input files if they are file paths rather than data frames
  if (is.data.frame(ResidMatFile)) {
    ResidMat <- ResidMatFile
  } else {
    ResidMat <- data.table::fread(ResidMatFile)
  }
  SparseGRM <- data.table::fread(SparseGRMFile)
  PairwiseIBD <- data.table::fread(PairwiseIBDFile)

  # Ensure all ID columns are character type for consistent matching
  ResidMat$SubjID <- as.character(ResidMat$SubjID)
  SparseGRM$ID1 <- as.character(SparseGRM$ID1)
  SparseGRM$ID2 <- as.character(SparseGRM$ID2)
  PairwiseIBD$ID1 <- as.character(PairwiseIBD$ID1)
  PairwiseIBD$ID2 <- as.character(PairwiseIBD$ID2)

  control <- checkControl.SPAGRM.NullModel(control, ResidMat, SparseGRM, PairwiseIBD)

  MaxQuantile <- control$MaxQuantile
  MinQuantile <- control$MinQuantile
  OutlierRatio <- control$OutlierRatio
  ControlOutlier <- control$ControlOutlier
  MaxNuminFam <- control$MaxNuminFam
  MAF_interval <- control$MAF_interval

  SubjID <- ResidMat$SubjID
  SparseGRM <- SparseGRM %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
  PairwiseIBD <- PairwiseIBD %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)

  # Use residual information to define outliers / non-outliers
  Resid <- ResidMat$Resid
  Quant <- quantile(Resid, probs = c(MinQuantile, MaxQuantile))
  Range <- max(Quant) - min(Quant)
  cutoffVec <- c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)

  .message("Outlier cutoffs: [%.3f, %.3f]", cutoffVec[1], cutoffVec[2])
  ResidMat$Outlier <- ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
    TRUE, FALSE
  )

  if (ControlOutlier) {
    .message("ControlOutlier=TRUE: keeping outliers <5%% (set FALSE for higher accuracy)")

    while (sum(ResidMat$Outlier) == 0) {
      OutlierRatio <- OutlierRatio * 0.8
      cutoffVec <- c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
      .message("Adjusted cutoffs: [%.3f, %.3f]", cutoffVec[1], cutoffVec[2])
      ResidMat$Outlier <- ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
        TRUE, FALSE
      )
      .message("Outliers found: %d", sum(ResidMat$Outlier))
    }

    while (sum(ResidMat$Outlier) / nrow(ResidMat) > 0.05) {
      OutlierRatio <- OutlierRatio + 0.5
      cutoffVec <- c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
      .message("Reducing outliers: cutoffs [%.3f, %.3f]", cutoffVec[1], cutoffVec[2])
      ResidMat$Outlier <- ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
        TRUE, FALSE
      )
      .message("Outliers: %d (%.1f%%)", sum(ResidMat$Outlier),
               100 * sum(ResidMat$Outlier) / nrow(ResidMat))
    }
  }

  .message("Outlier summary:")
  # Only show outlier info in debug mode to avoid cluttering console
  if (nrow(ResidMat %>% filter(Outlier == TRUE)) > 0) {
    outlier_summary <- ResidMat %>%
      filter(Outlier == TRUE) %>%
      dplyr::select(SubjID, Resid, Outlier) %>%
      arrange(Resid)
    .message("Found %d outliers (range: %.3f to %.3f)",
             nrow(outlier_summary), min(outlier_summary$Resid), max(outlier_summary$Resid))
  }

  # Decompose the subjects based on family structure and use a greedy algorithm to reduce family size if needed
  SparseGRM1 <- SparseGRM
  SparseGRM1$pos1 <- ResidMat$Resid[match(SparseGRM$ID1, ResidMat$SubjID)]
  SparseGRM1$pos2 <- ResidMat$Resid[match(SparseGRM$ID2, ResidMat$SubjID)]
  SparseGRM1 <- SparseGRM1 %>% mutate(Cov = abs(Value * pos1 * pos2))

  edges <- t(SparseGRM1[, c("ID1", "ID2")])
  graph_GRM <- igraph::make_graph(edges, directed = FALSE)
  graph_list_all <- graph_GRM %>% igraph::decompose()
  graph_length <- lapply(graph_list_all, length)

  graph_list_1 <- graph_list_all[graph_length == 1]
  SubjID.unrelated <- lapply(graph_list_1, igraph::get.vertex.attribute) %>% unlist(use.names = FALSE)
  ResidMat.unrelated <- ResidMat %>% filter(SubjID %in% SubjID.unrelated)
  SubjID.unrelated.nonOutlier <- ResidMat.unrelated %>%
    filter(Outlier == FALSE) %>%
    select(SubjID) %>%
    unlist(use.names = FALSE)

  # Values used in association analysys
  R_GRM_R <- SparseGRM1 %>%
    filter(ID1 %in% SubjID.unrelated) %>%
    select(Cov) %>%
    sum()
  sum_R_nonOutlier <- ResidMat.unrelated %>%
    filter(Outlier == FALSE) %>%
    select(Resid) %>%
    sum()
  R_GRM_R_nonOutlier <- SparseGRM1 %>%
    filter(ID1 %in% SubjID.unrelated.nonOutlier) %>%
    select(Cov) %>%
    sum()
  Resid.unrelated.outliers <- ResidMat.unrelated %>%
    filter(Outlier == TRUE) %>%
    select(Resid) %>%
    unlist(use.names = FALSE)
  R_GRM_R_TwoSubjOutlier <- 0
  TwoSubj_list <- ThreeSubj_list <- list()

  # initialize parameters
  graph_list_updated <- list()
  graph_list <- graph_list_all[graph_length > 1]
  nGraph <- length(graph_list)
  index.outlier <- 1

  if (nGraph != 0) {
    .message("Processing %d family groups with related residuals", nGraph)

    for (i in seq_len(nGraph)) {
      if (i %% 1000 == 0) {
        .message("Processing family group %d/%d", i, nGraph)
      }

      comp1 <- graph_list[[i]]
      comp3 <- igraph::V(comp1)$name

      # Step 0: calculate variance for the family
      pos1 <- match(comp3, SubjID)
      outlierInFam <- any(ResidMat$Outlier[pos1])

      block_GRM <- make.block.GRM(comp1, SparseGRM)

      R_GRM_R.temp <- as.numeric(t(Resid[pos1]) %*% block_GRM %*% Resid[pos1])
      R_GRM_R <- R_GRM_R + R_GRM_R.temp

      if (!outlierInFam) {
        sum_R_nonOutlier <- sum_R_nonOutlier + sum(ResidMat$Resid[pos1])
        R_GRM_R_nonOutlier <- R_GRM_R_nonOutlier + R_GRM_R.temp
        next
      }

      vcount <- igraph::vcount(comp1) # number of vertices

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

        pos2 <- match(comp13, SubjID)
        outlierInFam <- any(ResidMat$Outlier[pos2])

        block_GRM <- make.block.GRM(comp11, SparseGRM)

        R_GRM_R.temp <- as.numeric(t(Resid[pos2]) %*% block_GRM %*% Resid[pos2])

        if (!outlierInFam) {
          sum_R_nonOutlier <- sum_R_nonOutlier + sum(ResidMat$Resid[pos2])
          R_GRM_R_nonOutlier <- R_GRM_R_nonOutlier + R_GRM_R.temp
        } else {
          graph_list_updated[[index.outlier]] <- comp11
          index.outlier <- index.outlier + 1
        }
      }
    }

    .message("Building Chow-Liu tree for family outliers")

    # Make a list of array index.
    arr.index <- list()
    for (n in seq_len(MaxNuminFam)) {
      temp <- c()
      for (i in seq_len(n)) {
        indexString <- rep("c(1, 1, 1)", n)
        indexString[i] <- "0:2"
        indexString <- paste0(indexString, collapse = "%o%")
        cmd <- paste0("temp = c(temp, list(arr.index", i, "=", indexString, "))")
        eval(parse(text = cmd))
      }
      arr.index[[n]] <- temp
    }

    # build chou-liu-tree.
    n.outliers <- length(graph_list_updated)
    if (n.outliers != 0) {
      ## The below values are only used in chou.liu.tree
      TwofamID.index <- ThreefamID.index <- 0
      for (index.outlier in seq_len(n.outliers)) {
        if (index.outlier %% 1000 == 0) {
          .message("Processing CLT for outlier families: %d, %d/%d", TwofamID.index, ThreefamID.index, nGraph)
        }

        comp1 <- graph_list_updated[[index.outlier]]
        comp3 <- igraph::V(comp1)$name
        n1 <- length(comp3)
        pos3 <- match(comp3, SubjID)

        Resid.temp <- ResidMat$Resid[pos3]

        if (n1 == 1) {
          Resid.unrelated.outliers <- c(Resid.unrelated.outliers, Resid.temp)
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

        CLT <- chow.liu.tree(
          N = n1,
          IBD = tempIBD,
          IDs = comp3,
          MAF_interval = MAF_interval
        )

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
      .message(
        "Completed CLT processing for %d families (%d, %d/%d)",
        n.outliers, TwofamID.index, ThreefamID.index, nGraph
      )
    }
  }

  obj <- list(
    Resid = Resid, subjData = SubjID, N = length(SubjID), Resid.unrelated.outliers = Resid.unrelated.outliers,
    R_GRM_R = R_GRM_R, R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier,
    sum_R_nonOutlier = sum_R_nonOutlier, R_GRM_R_nonOutlier = R_GRM_R_nonOutlier,
    TwoSubj_list = TwoSubj_list, ThreeSubj_list = ThreeSubj_list,
    MAF_interval = MAF_interval
  )

  class(obj) <- "SPAGRM_NULL_Model"

  return(obj)
}


checkControl.SPAGRM.NullModel <- function(
  control,
  ResidMat,
  SparseGRM,
  PairwiseIBD
) {
  default.control <- list(
    MaxQuantile = 0.75,
    MinQuantile = 0.25,
    OutlierRatio = 1.5,
    ControlOutlier = TRUE,
    MaxNuminFam = 5,
    MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  )

  control <- updateControl(control, default.control) # This file is in 'control.R'

  if (control$MaxQuantile < control$MinQuantile) {
    stop("MaxQuantile(default is 0.75) should be larger than MinQuantile(default is 0.25).")
  }

  if (control$OutlierRatio < 0) {
    stop("OutlierRatio should be larger than or equal 0 (default is 1.5).")
  }

  if (any(colnames(ResidMat) != c("SubjID", "Resid"))) {
    stop("The column names of ResidMat should be ['SubjID', 'Resid'].")
  }

  if (any(colnames(SparseGRM) != c("ID1", "ID2", "Value"))) {
    stop("The column names of SparseGRM should be ['ID1', 'ID2', 'Value'].")
  }

  if (any(colnames(PairwiseIBD) != c("ID1", "ID2", "pa", "pb", "pc"))) {
    stop("The column names of PairwiseIBD should be ['ID1', 'ID2', 'pa', 'pb', 'pc'].")
  }

  SubjID.In.Resid <- ResidMat$SubjID
  SubjID.In.GRM <- unique(c(SparseGRM$ID1, SparseGRM$ID2))
  SubjID.In.IBD <- unique(c(PairwiseIBD$ID1, PairwiseIBD$ID2))

  if (any(!SubjID.In.Resid %in% SubjID.In.GRM)) {
    stop("At least one subject in residual matrix does not have GRM information.")
  }

  if (any(!SubjID.In.IBD %in% SubjID.In.GRM)) {
    stop("At least one subject has IBD information but does not have GRM information.")
  }

  return(control)
}


checkControl.Marker.SPAGRM <- function(control, MAF_interval) {

  # Validate MAF interval constraints specific to SPAGRM
  if (length(MAF_interval) > 1) {
    if (control$min_maf_marker <= min(MAF_interval)) {
      stop(
        "min_maf_marker is out of MAF_interval. ",
        "Please reset min_maf_marker or check MAF_interval."
      )
    }
  }

  default.control <- list(
    zeta = 0,
    tol = 1e-5
  )

  control <- updateControl(control, default.control)

  return(control)
}


setMarker.SPAGRM <- function(objNull, control) {
  # Initialize marker-level analysis in C++ for SPAGRM method
  setSPAGRMobjInCPP(
    t_resid = objNull$Resid,                              # numeric vector: Residuals from null model
    t_resid_unrelated_outliers = objNull$Resid.unrelated.outliers,  # numeric vector: Outlier residuals
    t_sum_R_nonOutlier = objNull$sum_R_nonOutlier,        # numeric: Sum of non-outlier residuals
    t_R_GRM_R_nonOutlier = objNull$R_GRM_R_nonOutlier,    # numeric: R'*GRM*R for non-outliers
    t_R_GRM_R_TwoSubjOutlier = objNull$R_GRM_R_TwoSubjOutlier,  # numeric: Two-subject outlier quad form
    t_R_GRM_R = objNull$R_GRM_R,                          # numeric: Full R'*GRM*R quadratic form
    t_MAF_interval = objNull$MAF_interval,                # numeric vector: MAF intervals for binning
    t_TwoSubj_list = objNull$TwoSubj_list,                # list: Two-subject outlier pair info
    t_ThreeSubj_list = objNull$ThreeSubj_list,            # list: Three-subject outlier combinations
    t_SPA_Cutoff = control$SPA_Cutoff,                    # numeric: P-value cutoff for SPA
    t_zeta = control$zeta,                                # numeric: SPA moment approximation parameter
    t_tol = control$tol                                   # numeric: Numerical tolerance for SPA
  )
}


mainMarker.SPAGRM <- function(genoType, genoIndex) {
  # Perform main marker analysis for SPAGRM method
  OutList <- mainMarkerInCPP(
    t_method = "SPAGRM",      # character: Statistical method name
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
    zScore = OutList$zScore, # standardized score statistics
    Pvalue = OutList$pvalVec, # marker-level p-value
    hwepval = OutList$hwepvalVec # Hardy-Weinberg equilibrium p-value
  )

  return(obj.mainMarker)
}


make.block.GRM <- function(
  graph,
  GRM # three columns: "ID1", "ID2", and "Value"
) {
  comp2 <- igraph::get.data.frame(graph)

  # igraph gives an unexpected additional loop, which may change the block GRM
  # the below is to remove the additional loop
  comp2 <- comp2[!duplicated(comp2), ]

  comp3 <- igraph::V(graph)$name

  colnames(GRM) <- c("to", "from", "Value")

  n1 <- nrow(comp2)
  comp2 <- merge(comp2, GRM)
  n2 <- nrow(comp2)

  if (n1 != n2) {
    stop("Internal error, 'n1 != n2'.")
  }

  block_GRM <- Matrix::sparseMatrix(
    i = match(comp2$from, comp3),
    j = match(comp2$to, comp3),
    x = comp2$Value,
    symmetric = TRUE
  )
  return(block_GRM)
}


chow.liu.tree <- function(
  N,
  IBD,
  IDs,
  MAF_interval
) {
  # Build Chow-Liu tree for modeling genetic dependencies in family structures
  CLT <- c()

  # Iterate through different minor allele frequency intervals
  for (index in seq_len(length(MAF_interval))) {
    mu <- MAF_interval[index]

    # Calculate baseline genotype probabilities p = c(G0, G1, G2)
    p0 <- c((1 - mu)^2, 2 * mu * (1 - mu), mu^2)

    # Calculate allele frequency probabilities for different IBD states
    pa.allele2 <- c((1 - mu)^2, 0, 0, 0, 2 * mu * (1 - mu), 0, 0, 0, mu^2)

    # IBD probability for allele 1 sharing
    pb.allele1 <- c(
      (1 - mu)^3, mu * (1 - mu)^2, 0,
      mu * (1 - mu)^2, mu * (1 - mu), mu^2 * (1 - mu),
      0, mu^2 * (1 - mu), mu^3
    )

    pc.allele0 <- c(
      (1 - mu)^4, 2 * mu * (1 - mu)^3, mu^2 * (1 - mu)^2, 2 * mu * (1 - mu)^3, 4 * mu^2 * (1 - mu)^2,
      2 * mu^3 * (1 - mu), mu^2 * (1 - mu)^2, 2 * mu^3 * (1 - mu), mu^4
    )

    # calculate entropy I(Gi, Gj). Noting that entropy of unrelated pairs is zero.
    for (j in seq_len(nrow(IBD))) {
      pro <- IBD$pa[j] * pa.allele2 + IBD$pb[j] * pb.allele1 + IBD$pc[j] * pc.allele0

      entropy <- sum(pro * log(pro / pc.allele0), na.rm = TRUE)
      IBD$entropy[j] <- entropy
    }

    # use the "prim" lgorithm to bulid a maximum spanning tree.
    Max_span_tree <- IBD %>%
      igraph::graph_from_data_frame(directed = TRUE) %>%
      igraph::mst(weights = -IBD$entropy, algorithm = "prim") %>%
      igraph::get.edgelist() %>%
      data.table::as.data.table() %>%
      rename(ID1 = V1, ID2 = V2)

    mst.IBD <- merge(Max_span_tree, IBD, all.x = TRUE) %>%
      mutate(idxID1 = match(ID1, IDs), idxID2 = match(ID2, IDs))

    arr.prob <- array(1, dim = rep(3, N))
    for (i in seq_len(N)) {
      dimnames(arr.prob)[[i]] <- paste0("ID", i, ":", 0:2)
    }

    vec <- c(mst.IBD$idxID1, mst.IBD$idxID2)
    vec <- vec[duplicated(vec)]

    for (k in seq_len(N - 1)) {
      pro <- mst.IBD$pa[k] * pa.allele2 + mst.IBD$pb[k] * pb.allele1 + mst.IBD$pc[k] * pc.allele0

      matrix.prob <- matrix(pro, 3, 3)  # Used in eval() below
      matrix.index1 <- mst.IBD$idxID1[k]
      matrix.index2 <- mst.IBD$idxID2[k]
      for (i in seq_len(3)) {
        for (j in seq_len(3)) {
          indexString <- rep("", N)
          indexString[matrix.index1] <- i
          indexString[matrix.index2] <- j
          indexString <- paste0(indexString, collapse = ",")
          cmd <- paste0("arr.prob[", indexString, "] = arr.prob[", indexString, "] * matrix.prob[", i, ",", j, "]")
          eval(parse(text = cmd))
        }
      }
    }

    for (k in seq_len(N - 2)) {
      vector.prob <- p0  # Used in eval() below
      vector.index <- vec[k]
      for (i in seq_len(3)) {
        indexString <- rep("", N)
        indexString[vector.index] <- i
        indexString <- paste0(indexString, collapse = ",")
        cmd <- paste0("arr.prob[", indexString, "] = arr.prob[", indexString, "] / vector.prob[", i, "]")
        eval(parse(text = cmd))
      }
    }

    CLT <- cbind(CLT, c(arr.prob))
  }

  return(CLT)
}
