## ------------------------------------------------------------------------------
## POLMM.R
## Core implementation and helpers for the POLMM method (ordinal mixed model)
## covering null model fitting, marker-level and region-level (POLMM-GENE) tests.
##
## Functions:
##   GRAB.POLMM                  : Print brief method information for marker-level analysis.
##   GRAB.POLMM.Region           : Print brief method information for region-level analysis.
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
#' For region-based analysis using POLMM-GENE, see \code{\link{GRAB.POLMM.Region}}.
#'
#' @return NULL
#'
#' @examples
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultPOLMMmarker.txt")
#'
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' PhenoData$OrdinalPheno <- factor(PhenoData$OrdinalPheno, levels = c(0, 1, 2))
#
#' # Step 1
#' obj.POLMM <- GRAB.NullModel(
#'  OrdinalPheno ~ AGE + GENDER,
#'  data = PhenoData,
#'  subjIDcol = "IID",
#'  method = "POLMM",
#'  traitType = "ordinal",
#'  GenoFile = GenoFile,
#'  SparseGRMFile = SparseGRMFile,
#'  control = list(tolTau = 0.2, tolBeta = 0.1)
#' )
#'
#' # Step 2
#' GRAB.Marker(obj.POLMM, GenoFile, OutputFile,
#'   control = list(ifOutGroup = TRUE))
#'
#' head(data.table::fread(OutputFile))
#'
#' @details
#' \strong{Genotype file}: \code{GenoFile} is mandatory for \code{GRAB.NullModel()} when using POLMM.
#' It is required for estimating the variance ratio parameter, which is essential for calibrating
#' the test statistics in subsequent association tests.
#'
#' \strong{Genetic Relationship Matrix (GRM) Options}:
#' POLMM supports both sparse and dense GRM for modeling genetic relatedness:
#' \itemize{
#'   \item If \code{SparseGRMFile} is provided to \code{GRAB.NullModel()},
#'      the sparse GRM will be used in model fitting.
#'   \item If \code{SparseGRMFile} is not provided, \code{GRAB.NullModel()}
#'      will calculate a dense GRM from \code{GenoFile}.
#' }
#'
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
#'   \item \code{showInfo} (logical, default: FALSE): Whether to print PCG iteration information for debugging.
#'   \item \code{maxiterEps} (integer, default: 100): Maximum iterations for epsilon estimation.
#'   \item \code{tolEps} (numeric, default: 1e-10): Tolerance for epsilon estimation.
#'   \item \code{minMafVarRatio} (numeric, default: 0.1): Minimum MAF for variance ratio estimation.
#'   \item \code{maxMissingVarRatio} (numeric, default: 0.1): Maximum missing rate for variance ratio estimation.
#'   \item \code{nSNPsVarRatio} (integer, default: 20): Number of SNPs used for variance ratio estimation.
#'   \item \code{CVcutoff} (numeric, default: 0.0025): Coefficient of variation cutoff.
#'   \item \code{grainSize} (integer, default: 1): Grain size for parallel processing.
#'   \item \code{minMafGRM} (numeric, default: 0.01): Minimum MAF for GRM construction.
#'   \item \code{maxMissingGRM} (numeric, default: 0.1): Maximum missing rate for GRM construction.
#' }
#'
#' \strong{Additional Control Parameters for GRAB.Marker()}:
#' \itemize{
#'   \item \code{ifOutGroup} (logical, default: FALSE): Whether to output group-specific statistics
#'      (alternative allele frequency, counts, and sample size for each ordinal category).
#'      When TRUE, adds columns AltFreqInGroup.*, AltCountsInGroup.*, and nSamplesInGroup.* to the output file.
#' }
#'
#' \strong{Output file columns}:
#' \describe{
#'   \item{Marker}{Marker identifier (rsID or CHR:POS:REF:ALT).}
#'   \item{Info}{Marker information in format CHR:POS:REF:ALT.}
#'   \item{AltFreq}{Alternative allele frequency in the overall sample.}
#'   \item{AltCounts}{Total count of alternative alleles.}
#'   \item{MissingRate}{Proportion of missing genotypes.}
#'   \item{Pvalue}{P-value from the score test.}
#'   \item{beta}{Effect size estimate (log-odds scale).}
#'   \item{seBeta}{Standard error of beta.}
#'   \item{zScore}{Z-score from the score test.}
#'   \item{AltFreqInGroup.1, AltFreqInGroup.2, ...}{(Only if \code{ifOutGroup = TRUE})
#'      Alternative allele frequency in each ordinal category.}
#'   \item{AltCountsInGroup.1, AltCountsInGroup.2, ...}{(Only if \code{ifOutGroup = TRUE})
#'      Alternative allele counts in each ordinal category.}
#'   \item{nSamplesInGroup.1, nSamplesInGroup.2, ...}{(Only if \code{ifOutGroup = TRUE}) Sample size in each ordinal category.}
#' }
#'
GRAB.POLMM <- function() {
  .message("?GRAB.POLMM for instructions on step 1, and step 2 of marker-level analysis")
  .message("?GRAB.POLMM.Region for instructions on step 2 of region-based analysis")
}


#' Region-based analysis using POLMM-GENE method
#'
#' POLMM-GENE is an extension of POLMM for region-based (gene-based or set-based) analysis
#' of ordinal categorical traits in related samples. It uses a variant-set mixed model framework
#' with SKAT-O, SKAT, and Burden tests.
#' For single-variant tests, see \code{\link{GRAB.POLMM}}.
#'
#' @return NULL
#'
#' @examples
#' GenoFileStep1 <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' GenoFileStep2 <- system.file("extdata", "simuPLINK_RV.bed", package = "GRAB")
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' GroupFile <- system.file("extdata", "simuPLINK_RV.group", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "resultPOLMMregion.txt")
#'
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' PhenoData$OrdinalPheno <- factor(PhenoData$OrdinalPheno, levels = c(0, 1, 2))
#
#' # Step 1
#' obj.POLMM <- GRAB.NullModel(
#'  OrdinalPheno ~ AGE + GENDER,
#'  data = PhenoData,
#'  subjIDcol = "IID",
#'  method = "POLMM",
#'  traitType = "ordinal",
#'  GenoFile = GenoFileStep1,
#'  SparseGRMFile = SparseGRMFile,
#'  control = list(tolTau = 0.2, tolBeta = 0.1)
#' )
#'
#' # Step 2
#' GRAB.Region(obj.POLMM, GenoFileStep2, OutputFile,
#'   GroupFile = GroupFile,
#'   SparseGRMFile = SparseGRMFile,
#'   MaxMAFVec = "0.01,0.005"
#' )
#'
#' head(data.table::fread(OutputFile))
#' head(data.table::fread(paste0(OutputFile, ".markerInfo")))
#' head(data.table::fread(paste0(OutputFile, ".otherMarkerInfo")))
#' head(data.table::fread(paste0(OutputFile, ".infoBurdenNoWeight")))
#'
#' @details
#'
#' See \code{\link{GRAB.POLMM}} for details on step 1.
#'
#' \strong{Additional Control Parameters for GRAB.Region() with POLMM}:
#' \itemize{
#'   \item \code{showInfo} (logical, default: FALSE): Whether to print PCG iteration information for debugging.
#'   \item \code{tolPCG} (numeric, default: 0.001): Tolerance for PCG in region testing.
#'   \item \code{maxiterPCG} (integer, default: 100): Maximum PCG iterations in region testing.
#' }
#'
#' \strong{Results are saved to four files}:
#' \enumerate{
#'   \item \code{OutputFile}: Region-based test results (SKAT-O, SKAT, Burden p-values).
#'   \item \code{paste0(OutputFile, ".markerInfo")}: Marker-level results for rare variants
#'     (MAC >= \code{min_mac_region}) included in region tests.
#'   \item \code{paste0(OutputFile, ".otherMarkerInfo")}: Information for excluded markers
#'     (ultra-rare variants or failed QC).
#'   \item \code{paste0(OutputFile, ".infoBurdenNoWeight")}: Summary statistics for burden
#'     tests without weights.
#' }
#'
#' **Region-level results** (\code{OutputFile}) columns:
#' \describe{
#'   \item{Region}{Region identifier from \code{GroupFile}.}
#'   \item{nMarkers}{Number of rare variants with MAF < cutoff and MAC >= \code{min_mac_region}.}
#'   \item{nMarkersURV}{Number of ultra-rare variants with MAC < \code{min_mac_region}.}
#'   \item{Anno.Type}{Annotation type from \code{GroupFile}.}
#'   \item{MaxMAF.Cutoff}{Maximum MAF cutoff used for variant selection.}
#'   \item{pval.SKATO}{SKAT-O test p-value.}
#'   \item{pval.SKAT}{SKAT test p-value.}
#'   \item{pval.Burden}{Burden test p-value.}
#' }
#'
#' **Marker-level results** (\code{paste0(OutputFile, ".markerInfo")}) columns:
#' \describe{
#'   \item{Region}{Region identifier.}
#'   \item{ID}{Marker identifier.}
#'   \item{Info}{Marker information in format CHR:POS:REF:ALT.}
#'   \item{Anno}{Annotation from \code{GroupFile}.}
#'   \item{AltFreq}{Alternative allele frequency.}
#'   \item{MAC}{Minor allele count.}
#'   \item{MAF}{Minor allele frequency.}
#'   \item{MissingRate}{Proportion of missing genotypes.}
#'   \item{IndicatorVec}{Marker status indicator (1 = rare variant included, 3 = ultra-rare variant included).}
#'   \item{StatVec}{Score test statistic.}
#'   \item{altBetaVec}{Effect size estimate.}
#'   \item{seBetaVec}{Standard error of effect size estimate.}
#'   \item{pval0Vec}{Unadjusted p-value.}
#'   \item{pval1Vec}{SPA-adjusted p-value.}
#'   \item{posRow}{Position row index.}
#' }
#'
#' **Other marker info** (\code{paste0(OutputFile, ".otherMarkerInfo")}) columns:
#' \describe{
#'   \item{ID}{Marker identifier.}
#'   \item{Annos}{Annotation from \code{GroupFile}.}
#'   \item{Region}{Region identifier.}
#'   \item{Info}{Marker information in format CHR:POS:REF:ALT.}
#'   \item{Anno}{Annotation category.}
#'   \item{AltFreq}{Alternative allele frequency.}
#'   \item{MAC}{Minor allele count.}
#'   \item{MAF}{Minor allele frequency.}
#'   \item{MissingRate}{Proportion of missing genotypes.}
#'   \item{IndicatorVec}{Status indicator (0 or 2 for excluded markers).}
#' }
#'
#' **Burden test summary** (\code{paste0(OutputFile, ".infoBurdenNoWeight")}) columns:
#' \describe{
#'   \item{region}{Region identifier.}
#'   \item{anno}{Annotation type.}
#'   \item{max_maf}{Maximum MAF cutoff.}
#'   \item{sum}{Sum of genotypes.}
#'   \item{Stat}{Score test statistic.}
#'   \item{beta}{Effect size estimate.}
#'   \item{se.beta}{Standard error of effect size estimate.}
#'   \item{pvalue}{P-value for burden test.}
#' }
#'
GRAB.POLMM.Region <- function() {
  .message("?GRAB.POLMM for instructions on step 1")
  .message("?GRAB.POLMM.Region for instructions on step 2")
}


checkControl.NullModel.POLMM <- function(traitType, GenoFile, SparseGRMFile, control) {

  if (traitType != "ordinal") {
    stop("For method of 'POLMM', only traitType of 'ordinal' is supported.")
  }

  # GenoFile validation (required)
  if (is.null(GenoFile)) {
    stop("Argument 'GenoFile' is required for method 'POLMM'.")
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
    optionGRM <- "DenseGRM"
  }

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
    LOCO = FALSE,
    grainSize = 1,
    minMafGRM = 0.01,
    maxMissingGRM = 0.1,
    showInfo = FALSE
  )
  control <- updateControl(control, default.control)

  return(list(control = control, optionGRM = optionGRM))
}


#' Fit POLMM null model for ordinal outcomes
#'
#' Initializes the POLMM null model from an ordered categorical response and
#' covariate matrix, preparing C++ state for subsequent marker/region tests.
#'
#' @param response Ordered factor response (lowest level coded as 0 internally).
#' @param designMat Numeric covariate matrix or data.frame (n x p).
#' @param subjData Character vector of subject IDs aligned with rows of
#'   \code{designMat} and \code{response}.
#' @param control List of POLMM options (e.g., \code{tau},
#'   \code{maxMissingVarRatio}, \code{minMafVarRatio}).
#' @param optionGRM Character, either \code{"DenseGRM"} or \code{"SparseGRM"}.
#' @param GenoFile Character, path to genotype file (PLINK or BGEN format).
#' @param GenoFileIndex Character or NULL, path to genotype index files.
#'   If NULL, uses same prefix as \code{GenoFile}.
#' @param SparseGRMFile Character or NULL, path to sparse GRM file.
#'   If provided, sparse GRM is used; otherwise dense GRM is constructed.
#'
#' @return An object of class \code{"POLMM_NULL_Model"} representing the
#'   initialized null model; state is stored in C++ and not intended for direct
#'   element-wise access from R.
#'
#' @keywords internal
fitNullModel.POLMM <- function(
  response, designMat, subjData, control, optionGRM,
  GenoFile, GenoFileIndex, SparseGRMFile
) {
  
  ######## -------------- Setup GRM in C++ -------- ########

  genoList <- setGenoInput(GenoFile, GenoFileIndex, subjData, control) # list
  flagSparseGRM <- optionGRM == "SparseGRM"

  if (flagSparseGRM) {
    SparseGRM <- data.table::fread(SparseGRMFile)
    KinMatListR <- updateSparseGRM(as.data.frame(SparseGRM), subjData)
    setSparseGRMInCPP(t_KinMatListR = KinMatListR)         # C++ backend setup
  } else {
    if (genoList$genoType != "PLINK") {
      stop("If DenseGRM is used when fitting a null model, ",
        "then only PLINK format is supported.")
    }

    setDenseGRMInCPP(
      t_memoryChunk = control$memoryChunk,      # numeric: Memory allocation in GB for GRM
      t_minMafGRM = control$minMafGRM,          # numeric: Min MAF for variants in GRM
      t_maxMissingGRM = control$maxMissingGRM   # numeric: Max missing rate for GRM variants
    )
  }

  ######## -------------- Fit null model -------- ########

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

  # Set seed for reproducible marker sampling (if seed != -1)
  if (control$seed != -1) {
    set.seed(control$seed)
  }
  
  markerInfo <- genoList$markerInfo[sample(nrow(genoList$markerInfo)), ]

  # Main.cpp
  GenoMat <- getGenoInCPP_fixedNumber(
    t_genoType = genoList$genoType,               # character: "PLINK" or "BGEN"
    t_markerInfo = markerInfo,                    # data.frame: Marker info with genoIndex
    n = length(yVec),                             # integer: Sample size
    t_imputeMethod = "mean",                      # character: Imputation method
    m = 100,                                      # integer: Number of markers to select
    missingRateCutoff = control$maxMissingVarRatio, # numeric: Max missing rate cutoff
    minMAFCutoff = control$minMafVarRatio         # numeric: Min MAF cutoff
  )

  # The following function is in 'Main.cpp'
  objNull <- setPOLMMobjInCPP_NULL(
    t_flagSparseGRM = flagSparseGRM,   # logical: Use sparse (TRUE) or dense (FALSE) GRM
    t_Cova = Cova,                     # matrix: Covariate matrix (n x p) with intercept
    t_yVec = yVec,                     # integer vector: Response (0 to J-1 for J categories)
    t_beta = beta,                     # numeric vector: Fixed effect coefficients
    t_bVec = bVec,                     # numeric vector: Random effect coefficients
    t_eps = eps,                       # numeric vector: Threshold parameters
    t_tau = tau,                       # numeric: Variance component parameter
    t_SPmatR = SPmatR,                 # list: Sparse matrix representation (deprecated)
    t_controlList = control,           # list: Control parameters for optimization
    GenoMat = GenoMat                  # matrix: Genotype matrix for variance ratio estimation
  )

  class(objNull) <- "POLMM_NULL_Model"
  return(objNull)
}


checkControl.Marker.POLMM <- function(control) {

  default.control <- list(
    ifOutGroup = FALSE
  )
  control <- updateControl(control, default.control)

  return(control)
}


setMarker.POLMM <- function(objNull, control) {
  
  objCHR <- objNull$LOCOList[["LOCO=F"]]

  # Calculate grouping for phenotypic values
  Group <- objNull$yVec                                       # numeric vector
  nGroup <- length(unique(Group))                             # integer

  # Check 'Main.cpp'
  setPOLMMobjInCPP(
    t_muMat = objCHR$muMat,               # matrix: Mean probability matrix (n x J)
    t_iRMat = objCHR$iRMat,               # matrix: Inverse correlation matrix (n x (J-1))
    t_Cova = objNull$Cova,                # matrix: Covariate matrix (n x p) with intercept
    t_yVec = objNull$yVec,                # integer vector: Response (0 to J-1)
    t_tau = objNull$tau,                  # numeric: Variance component parameter
    t_printPCGInfo = FALSE,               # logical: (not used in marker)
    t_tolPCG = 0.001,                     # numeric: (not used in marker)
    t_maxiterPCG = 100,                   # integer: (not used in marker)
    t_varRatio = objCHR$VarRatio,         # numeric: Variance ratio from null model
    t_SPA_cutoff = control$SPA_Cutoff,    # numeric: Cutoff for saddlepoint approximation
    t_flagSparseGRM = FALSE,              # logical: Use sparse GRM (FALSE for marker)
    t_group = Group,                      # integer vector: Group assignments for each individual
    t_ifOutGroup = control$ifOutGroup,    # logical: Output group-specific statistics
    t_nGroup = nGroup                     # integer: Total number of groups
  )
}

mainMarker.POLMM <- function(
  genoType,
  genoIndex,
  control
) {

  OutList <- mainMarkerInCPP(
    t_method = "POLMM",       # character: Statistical method name
    t_genoType = genoType,    # character: "PLINK" or "BGEN"
    t_genoIndex = genoIndex   # integer vector: Genotype indices to analyze
  )

  obj.mainMarker <- data.frame(
    Marker = OutList$markerVec,               # Marker IDs
    Info = OutList$infoVec,                   # Marker info: CHR:POS:REF:ALT
    AltFreq = OutList$altFreqVec,             # Alternative allele frequencies
    AltCounts = OutList$altCountsVec,         # Alternative allele counts
    MissingRate = OutList$missingRateVec,     # Missing rates per marker
    Pvalue = OutList$pvalVec,                 # Association test p-values
    beta = OutList$beta,                      # Effect size estimates
    seBeta = OutList$seBeta,                  # Standard errors of beta
    zScore = OutList$zScore                   # Z-scores
  )

  if (control$ifOutGroup) {
    obj.mainMarker <- cbind(
      obj.mainMarker,
      AltFreqInGroup = OutList$AltFreqInGroup,      # Alt freq in analysis group
      AltCountsInGroup = OutList$AltCountsInGroup,  # Alt counts in analysis group
      nSamplesInGroup = OutList$nSamplesInGroup     # Sample size in analysis group
    )
  }

  return(obj.mainMarker)
}


checkControl.Region.POLMM <- function(control) {

  default.control <- list(
    showInfo = FALSE,
    tolPCG = 0.001,
    maxiterPCG = 100
  )
  control <- updateControl(control, default.control)

  return(control)
}


setRegion.POLMM <- function(
  objNull,
  control,
  SparseGRMFile
) {

  # Since region-level analysis mainly focuses on rare variants, we use sparse GRM for all markers

  .message("Using sparse GRM for POLMM-GENE analysis")

  # ---- BEGIN inlined: setSparseGRMInStep2 ----
  SparseGRM <- data.table::fread(SparseGRMFile)
  SparseGRM <- as.data.frame(SparseGRM)
  KinMatListR <- updateSparseGRM(SparseGRM, objNull$subjData)
  setSparseGRMInCPP(
    t_KinMatListR = KinMatListR  # list: Sparse kinship matrix (locations, values, nSubj)
  )
  # ---- END inlined: setSparseGRMInStep2 ----

  # Calculate grouping for phenotypic values (POLMM-specific)
  Group <- objNull$yVec                                       # numeric vector
  nGroup <- length(unique(Group))                             # integer
  ifOutGroup <- TRUE                                          # logical

  # region-level analysis uses sparse GRM
  objCHR <- objNull$LOCOList[["LOCO=F"]]
  flagSparseGRM <- TRUE

  # Check 'Main.cpp'
  setPOLMMobjInCPP(
    t_muMat = objCHR$muMat,                 # matrix: Mean probability matrix (n x J)
    t_iRMat = objCHR$iRMat,                 # matrix: Inverse correlation matrix (n x (J-1))
    t_Cova = objNull$Cova,                  # matrix: Covariate matrix (n x p) with intercept
    t_yVec = objNull$yVec,                  # integer vector: Response (0 to J-1)
    t_tau = objNull$tau,                    # numeric: Variance component parameter
    t_printPCGInfo = control$showInfo,      # logical: Print PCG iteration info for debugging
    t_tolPCG = control$tolPCG,              # numeric: PCG convergence tolerance
    t_maxiterPCG = control$maxiterPCG,      # integer: Max PCG iterations
    t_varRatio = objCHR$VarRatio,           # numeric: Variance ratio from null model
    t_SPA_cutoff = control$SPA_Cutoff,      # numeric: Cutoff for saddlepoint approximation
    t_flagSparseGRM = flagSparseGRM,        # logical: Use sparse GRM (TRUE for region)
    t_group = Group,                        # integer vector: Group assignments for each individual
    t_ifOutGroup = ifOutGroup,              # logical: Output group-specific statistics
    t_nGroup = nGroup                       # integer: Total number of groups
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
