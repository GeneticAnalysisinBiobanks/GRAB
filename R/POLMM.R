#' POLMM method in GRAB package
#'
#' POLMM method is to analyze ordinal categorical data for related samples in a large-scale biobank.
#'
#' @details
#' Please check \code{?GRAB.control} for the generic list of \code{control} in
#' \code{GRAB.NullModel()} and \code{GRAB.Marker()}.
#'
#' Additional list of \code{control} in \code{GRAB.NullModel()} function
#' Additional list of \code{control} in \code{GRAB.Marker()} function
#' Additional list of \code{control} in \code{GRAB.Region()} function
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
#'   subjData = PhenoData$IID,
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
#' OutputFile <- file.path(tempdir(), "simuMarkerOutput.txt")
#'
#' GRAB.Marker(
#'   objNull = obj.POLMM,
#'   GenoFile = GenoFile,
#'   OutputFile = OutputFile
#' )
#'
#' data.table::fread(OutputFile)
#'
#' ### Step 2(b): Set-based tests using POLMM-GENE
#' GenoFile <- system.file("extdata", "simuPLINK_RV.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "simuRegionOutput.txt")
#' GroupFile <- system.file("extdata", "simuPLINK_RV.group", package = "GRAB")
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#'
#' GRAB.Region(
#'   objNull = obj.POLMM,
#'   GenoFile = GenoFile,
#'   GenoFileIndex = NULL,
#'   OutputFile = OutputFile,
#'   OutputFileIndex = NULL,
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


# check the control list in null model fitting for POLMM method
#' @title Validate control parameters for POLMM null model
#' @param control List of control parameters for null model fitting.
#' @param traitType Character string specifying the trait type.
#' @param optionGRM Character string specifying GRM option.
#' @return Updated control list with validated parameters and defaults.
#' @keywords internal
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


# fit null model using POLMM method
fitNullModel.POLMM <- function(
  response, designMat, subjData, control, optionGRM,
  genoType, # "PLINK" or "BGEN"
  markerInfo # colnames: CHROM, POS, ID, REF, ALT, genoIndex
) {

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


#' @title Validate control parameters for POLMM marker testing
#' @param control List of control parameters for marker-level analysis.
#' @return Updated control list with validated parameters and defaults.
#' @keywords internal
checkControl.Marker.POLMM <- function(control) {
  default.control <- list(
    SPA_Cutoff = 2,
    outputColumns = c("beta", "seBeta")
  )

  control <- updateControl(control, default.control) # This file is in 'Util.R'
  return(control)
}


#' @title Set up marker-level testing for POLMM method
#' @param objNull Null model object from GRAB.NullModel().
#' @param control List of control parameters.
#' @param chrom Character string specifying chromosome for LOCO analysis.
#' @return List containing setup parameters for marker testing.
#' @keywords internal
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


# check the control list in region-level testing
#' @title Validate control parameters for POLMM region testing
#' @param control List of control parameters for region-based analysis.
#' @return Updated control list with validated parameters and defaults.
#' @keywords internal
checkControl.Region.POLMM <- function(control) {
  default.control <- list(
    SPA_Cutoff = 2,
    outputColumns = c("beta", "seBeta")
  )

  control <- updateControl(control, default.control) # This file is in 'Util.R'

  # check the parameter
  if (!is.numeric(control$SPA_Cutoff) || control$SPA_Cutoff <= 0) {
    stop("control$SPA_Cutoff should be a numeric value > 0")
  }

  return(control)
}


# Used in setRegion() function in GRAB_Region.R
#' @title Set up region-based testing for POLMM method
#' @param objNull Null model object from GRAB.NullModel().
#' @param control List of control parameters.
#' @param chrom Character string specifying chromosome for LOCO analysis.
#' @param SparseGRMFile Character string specifying sparse GRM file path.
#' @return List containing setup parameters for region testing.
#' @keywords internal
setRegion.POLMM <- function(
  objNull,
  control,
  chrom,
  SparseGRMFile
) {
  if (chrom != "LOCO=F") {
    .message("Chromosome: %s", chrom)
    if (!"LOCOList" %in% names(objNull)) {
      stop("If argument 'chrom' is not 'LOCO=FALSE', then objNull should includes element of 'LOCOList'.")
    }
  }

  objCHR <- objNull$LOCOList[[chrom]]

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


#' @title Perform region-based analysis for POLMM method
#' @param genoType Character string specifying genotype file format.
#' @param genoIndex Integer vector of genotype indices to analyze.
#' @param OutputFile Character string specifying output file path.
#' @param control List of control parameters.
#' @param n Integer specifying number of subjects.
#' @param obj.setRegion List containing region setup parameters.
#' @param obj.mainRegionInCPP List containing C++ analysis results.
#' @param nLabel Integer specifying number of labels/categories.
#' @return Data frame containing analysis results.
#' @keywords internal
mainRegion.POLMM <- function(
  genoType,
  genoIndex,
  OutputFile,
  control,
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
