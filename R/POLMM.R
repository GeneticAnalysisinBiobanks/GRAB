## ------------------------------------------------------------------------------
## POLMM.R
## Core implementation and helpers for the POLMM method (ordinal mixed model)
## covering null model fitting, marker-level and region-level (POLMM-GENE) tests.
##
## Functions:
##   GRAB.POLMM                  : Print brief method information.
##   checkControl.NullModel.POLMM: Validate and populate null model control list.
##   fitNullModel.POLMM          : Fit the POLMM null model (C++ backend setup).
##   checkControl.Marker.POLMM   : Validate marker-level control parameters.
##   setMarker.POLMM             : Initialize marker-level analysis objects.
##   checkControl.Region.POLMM   : Validate region-level control parameters.
##   setRegion.POLMM             : Prepare region (gene/set) analysis context.
##   mainRegion.POLMM            : Run region-based association tests.
## ------------------------------------------------------------------------------

#' POLMM method in GRAB package
#'
#' POLMM method is to analyze ordinal categorical data for related samples in a large-scale biobank.
#'
#' @details
#' \strong{Additional Control Parameters for GRAB.NullModel()}:
#' \itemize{
#'   \item \code{memoryChunk} (numeric, default: 2): Memory chunk size for computation.
#'   \item \code{seed} (integer, default: -1): Random seed (-1 means no seed is set).
#'   \item \code{tracenrun} (integer, default: 30): Number of runs for trace calculation.
#'   \item \code{maxiter} (integer, default: 100): Maximum number of iterations for model fitting.
#'   \item \code{tolBeta} (numeric, default: 0.001): Convergence tolerance for beta estimates.
#'   \item \code{tolTau} (numeric, default: 0.002): Convergence tolerance for tau estimates.
#'   \item \code{tau} (numeric, default: 0.2): Initial variance component value.
#'   \item \code{maxiterPCG} (integer, default: 100): Maximum iterations for preconditioned conjugate gradient.
#'   \item \code{tolPCG} (numeric, default: 1e-6): Tolerance for preconditioned conjugate gradient.
#'   \item \code{maxiterEps} (integer, default: 100): Maximum iterations for epsilon estimation.
#'   \item \code{tolEps} (numeric, default: 1e-10): Tolerance for epsilon estimation.
#'   \item \code{minMafVarRatio} (numeric, default: 0.1): Minimum MAF for variance ratio estimation.
#'   \item \code{maxMissingVarRatio} (numeric, default: 0.1): Maximum missing rate for variance ratio estimation.
#'   \item \code{nSNPsVarRatio} (integer, default: 20): Number of SNPs used for variance ratio estimation.
#'   \item \code{CVcutoff} (numeric, default: 0.0025): Coefficient of variation cutoff.
#'   \item \code{stackSize} (character, default: "auto"): Stack size for parallel processing.
#'   \item \code{grainSize} (integer, default: 1): Grain size for parallel processing.
#'   \item \code{minMafGRM} (numeric, default: 0.01): Minimum MAF for GRM construction.
#'   \item \code{maxMissingGRM} (numeric, default: 0.1): Maximum missing rate for GRM construction.
#'   \item \code{showInfo} (logical, default: TRUE): Whether to show progress information.
#'   \item \code{onlyCheckTime} (logical, default: FALSE): Whether to only check computation time.
#'   \item \code{LOCO} (logical, default: TRUE for DenseGRM, FALSE for SparseGRM): Leave-one-chromosome-out analysis.
#' }
#'
#' \strong{Additional Control Parameters for GRAB.Marker()}:
#' \itemize{
#'   \item \code{SPA_Cutoff} (numeric, default: 2): Cutoff for saddlepoint approximation.
#'   \item \code{outputColumns} (character vector, default: c("beta", "seBeta")): Columns to include in output. 
#'     Optional columns: "beta", "seBeta", "zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup".
#' }
#'
#' \strong{Additional Control Parameters for GRAB.Region()}:
#' \itemize{
#'   \item \code{SPA_Cutoff} (numeric, default: 2): Cutoff for saddlepoint approximation.
#' }
#'
#' @return No return value, called for side effects (prints information about the POLMM method to the console).
#'
#' @examples
#' ### First, read phenotype data and convert to a factor
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' PhenoData$OrdinalPheno <- factor(PhenoData$OrdinalPheno, levels = c(0, 1, 2))
#'
#' ### Step 1: Fit a null model
#' # If SparseGRMFile is provided, the sparse GRM will be used in model fitting.
#' # If SparseGRMFile isn't provided, GRAB.NullModel() will calculate dense GRM from GenoFile.
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#'
#' obj.POLMM <- GRAB.NullModel(
#'   formula = OrdinalPheno ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjIDcol = "IID",
#'   method = "POLMM",
#'   traitType = "ordinal",
#'   GenoFile = GenoFile,
#'   SparseGRMFile = SparseGRMFile,
#'   control = list(
#'     showInfo = FALSE,
#'     LOCO = FALSE,
#'     tolTau = 0.2,
#'     tolBeta = 0.1
#'   )
#' )
#'
#' ### Step 2(a): Single-variant tests using POLMM
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultPOLMMmarker.txt")
#' outputColumns <- c(
#'   "beta", "seBeta", "zScore",
#'   "nSamplesInGroup", "AltCountsInGroup", "AltFreqInGroup"
#' )
#'
#' GRAB.Marker(obj.POLMM, GenoFile, OutputFile,
#'   control = list(outputColumns = outputColumns)
#' )
#'
#' data.table::fread(OutputFile)
#'
#' ### Step 2(b): Set-based tests using POLMM-GENE
#' GenoFile <- system.file("extdata", "simuPLINK_RV.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultPOLMMregion.txt")
#' GroupFile <- system.file("extdata", "simuPLINK_RV.group", package = "GRAB")
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#'
#' GRAB.Region(obj.POLMM, GenoFile, OutputFile,
#'   GroupFile = GroupFile,
#'   SparseGRMFile = SparseGRMFile,
#'   MaxMAFVec = "0.01,0.005"
#' )
#'
#' data.table::fread(OutputFile)
#'
GRAB.POLMM <- function() {
  .message("Using POLMM method - see ?GRAB.POLMM for details")
}


checkControl.NullModel.POLMM <- function(
  control,
  traitType,
  optionGRM
) {
  if (traitType != "ordinal") {
    stop("For method of 'POLMM', only traitType of 'ordinal' is supported.")
  }

  # default setting of control for POLMM method
  default.control <- list(
    memoryChunk = 2,
    seed = -1, # use -1 to indicate no seed should be set
    tracenrun = 30,
    maxiter = 100,
    tolBeta = 0.001,
    tolTau = 0.002,
    tau = 0.2,
    maxiterPCG = 100,
    tolPCG = 1e-6,
    maxiterEps = 100,
    tolEps = 1e-10,
    minMafVarRatio = 0.1,
    maxMissingVarRatio = 0.1,
    nSNPsVarRatio = 20,
    CVcutoff = 0.0025,
    stackSize = "auto",
    grainSize = 1,
    minMafGRM = 0.01,
    maxMissingGRM = 0.1,
    showInfo = TRUE,
    onlyCheckTime = FALSE
  )

  if (optionGRM == "DenseGRM") {
    default.control$LOCO <- TRUE
  }

  if (optionGRM == "SparseGRM") {
    default.control$LOCO <- FALSE
  }

  control <- updateControl(control, default.control)
  return(control)
}


#' Fit POLMM null model for ordinal outcomes
#'
#' Initializes the POLMM null model from an ordered categorical response and
#' covariate matrix, preparing C++ state for subsequent marker/region tests.
#'
#' @param response Ordered factor response (lowest level coded as 0 internally).
#' @param designMat Numeric covariate matrix or data.frame (n x p).
#' @param subjIDcol Character vector of subject IDs aligned with rows of
#'   \code{designMat} and \code{response}.
#' @param control List of POLMM options (e.g., \code{tau}, \code{LOCO},
#'   \code{maxMissingVarRatio}, \code{minMafVarRatio}).
#' @param optionGRM Character, either \code{"DenseGRM"} or \code{"SparseGRM"}.
#' @param genoType Character, genotype backend: \code{"PLINK"} or \code{"BGEN"}.
#' @param markerInfo Data frame of marker metadata with columns
#'   \code{CHROM, POS, ID, REF, ALT, genoIndex} used to fetch genotypes.
#'
#' @return An object of class \code{"POLMM_NULL_Model"} representing the
#'   initialized null model; state is stored in C++ and not intended for direct
#'   element-wise access from R.
#'
#' @keywords internal
fitNullModel.POLMM <- function(
  response, designMat, subjIDcol, control, optionGRM,
  genoType, # "PLINK" or "BGEN"
  markerInfo # colnames: CHROM, POS, ID, REF, ALT, genoIndex
) {

  subjData <- subjIDcol  # Use subjData internally for compatibility

  ######## -------------- first set up the object in C++ -------- ########

  if (!is.factor(response)) {
    stop(
      "The response variable in POLMM method should be a factor. The class of the current response variable is '",
      class(response), "'."
    )
  }

  obj.clm <- summary(ordinal::clm(response ~ designMat))
  beta <- c(-1 * obj.clm$alpha[1], obj.clm$beta)
  eps <- c(0, obj.clm$alpha[-1] - obj.clm$alpha[1])
  bVec <- rep(0, length(response)) # initiate random effect of 0

  yVec <- as.numeric(response) - 1 # "-1" means change from R to C++
  Cova <- cbind(1, designMat)
  tau <- control$tau

  # This value is not used any more, remove it later. 03/31/2021
  SPmatR <- list(
    locations = matrix(c(0, 0), 2, 1),
    values = rep(0, 1)
  )

  LOCO <- control$LOCO

  if (LOCO) {
    stop("Option of LOCO will be supported later. (2022-01-28)")
  }

  m <- nrow(markerInfo)
  markerInfo <- markerInfo[sample(m), ]

  # Main.cpp
  GenoMat <- getGenoInCPP_fixedNumber(
    genoType, markerInfo, length(yVec), "mean", 100,
    control$maxMissingVarRatio,
    control$minMafVarRatio
  )

  controlList <- control
  flagSparseGRM <- ifelse(optionGRM == "SparseGRM", TRUE, FALSE)

  # The following function is in 'Main.cpp'
  objNull <- setPOLMMobjInCPP_NULL(
    flagSparseGRM,
    Cova,
    yVec,
    beta,
    bVec,
    eps,
    tau,
    SPmatR,
    controlList,
    GenoMat
  )

  class(objNull) <- "POLMM_NULL_Model"
  return(objNull)
}


checkControl.Marker.POLMM <- function(control) {
  default.control <- list(
    SPA_Cutoff = 2,
    outputColumns = c("beta", "seBeta")
  )

  control <- updateControl(control, default.control)

  # Validate parameters
  if (!is.numeric(control$SPA_Cutoff) || control$SPA_Cutoff <= 0) {
    stop("control$SPA_Cutoff should be a numeric value > 0.")
  }

  return(control)
}


setMarker.POLMM <- function(
  objNull,
  control,
  chrom
) {
  if (objNull$control$LOCO) {
    if (!chrom %in% names(objNull$LOCOList)) {
      stop("If control$LOCO == TRUE, then 'chrom' should be in names(objNull$LOCOList).")
    }
    objCHR <- objNull$LOCOList[[chrom]]
  } else {
    objCHR <- objNull$LOCOList[["LOCO=F"]]
  }

  # marker-level analysis does not require the following parameters
  # Note: it might be not so accurate if min_mac_marker is very low
  flagSparseGRM <- FALSE
  printPCGInfo <- FALSE
  tolPCG <- 0.001
  maxiterPCG <- 100

  # Check 'Main.cpp'
  setPOLMMobjInCPP(
    objCHR$muMat,
    objCHR$iRMat,
    objNull$Cova,
    objNull$yVec, # 0 to J-1
    objNull$tau,
    printPCGInfo,
    tolPCG,
    maxiterPCG,
    objCHR$VarRatio,
    control$SPA_Cutoff,
    flagSparseGRM
  )
}

mainMarker.POLMM <- function(
  genoType,
  genoIndex,
  outputColumns
) {
  # Call optimized C++ implementation for computational efficiency
  OutList <- mainMarkerInCPP("POLMM", genoType, genoIndex)

  # Construct base output data frame with required columns
  obj.mainMarker <- data.frame(
    Marker = OutList$markerVec,        # Marker IDs
    Info = OutList$infoVec,            # Marker info: CHR:POS:REF:ALT
    AltFreq = OutList$altFreqVec,      # Alternative allele frequencies
    AltCounts = OutList$altCountsVec,  # Alternative allele counts
    MissingRate = OutList$missingRateVec, # Missing rates per marker
    Pvalue = OutList$pvalVec           # Association test p-values
  )

  # Add optional columns if requested by user
  optionalColumns <- c("beta", "seBeta", "zScore", "PvalueNorm",
                       "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  additionalColumns <- intersect(optionalColumns, outputColumns)

  if (length(additionalColumns) > 0) {
    obj.mainMarker <- cbind.data.frame(
      obj.mainMarker,
      as.data.frame(OutList[additionalColumns])
    )
  }

  return(obj.mainMarker)
}


checkControl.Region.POLMM <- function(control) {
  default.control <- list(
    SPA_Cutoff = 2
  )

  control <- updateControl(control, default.control)

  # Validate parameters
  if (!is.numeric(control$SPA_Cutoff) || control$SPA_Cutoff <= 0) {
    stop("control$SPA_Cutoff should be a numeric value > 0.")
  }

  return(control)
}


setRegion.POLMM <- function(
  objNull,
  control,
  SparseGRMFile
) {

  objCHR <- objNull$LOCOList[["LOCO=F"]]

  # Since region-level analysis mainly focuses on rare variants, we use sparse GRM for all markers

  .message("Using sparse GRM for POLMM-GENE analysis")

  # ---- BEGIN inlined: setSparseGRMInStep2 ----
  SparseGRM <- data.table::fread(SparseGRMFile)
  SparseGRM <- as.data.frame(SparseGRM)
  KinMatListR <- updateSparseGRM(SparseGRM, objNull$subjData)
  setSparseGRMInCPP(KinMatListR) # check Main.cpp
  # ---- END inlined: setSparseGRMInStep2 ----

  # The following parameters are not used any more
  flagSparseGRM <- TRUE
  printPCGInfo <- FALSE
  tolPCG <- 0.001
  maxiterPCG <- 100
  VarRatio <- 1

  # Check 'Main.cpp'
  setPOLMMobjInCPP(
    objCHR$muMat,
    objCHR$iRMat,
    objNull$Cova,
    objNull$yVec, # 0 to J-1
    objNull$tau,
    printPCGInfo,
    tolPCG,
    maxiterPCG,
    VarRatio,
    control$SPA_Cutoff,
    flagSparseGRM
  )
}


mainRegion.POLMM <- function(
  genoType,
  genoIndex,
  OutputFile,
  n,
  obj.setRegion,
  obj.mainRegionInCPP,
  nLabel
) {
  ## required columns for all methods
  info.Region <- with(obj.mainRegionInCPP, data.frame(
    ID = markerVec,
    Info = infoVec,
    Anno = AnnoVec,
    AltFreq = altFreqVec,
    MAC = MACVec,
    MAF = MAFVec,
    MissingRate = missingRateVec,
    IndicatorVec = indicatorVec,
    StatVec = StatVec,
    altBetaVec = altBetaVec,
    seBetaVec = seBetaVec,
    pval0Vec = pval0Vec,
    pval1Vec = pval1Vec,
    stringsAsFactors = FALSE
  ))

  if (nLabel != 1) {
    info.Region <- with(
      obj.mainRegionInCPP,
      cbind(info.Region, MACLabelMat, MAFLabelMat)
    )
  }

  RV.Markers <- info.Region %>%
    filter(IndicatorVec == 1 | IndicatorVec == 3)

  RV.Markers <- RV.Markers %>%
    mutate(posRow = seq_len(nrow(RV.Markers)))

  Other.Markers <- info.Region %>%
    filter(IndicatorVec == 2 | IndicatorVec == 0) %>%
    select(-(StatVec:pval1Vec))

  return(list(
    RV.Markers = RV.Markers,
    Other.Markers = Other.Markers,
    VarMat = obj.mainRegionInCPP$VarMat
  ))
}
