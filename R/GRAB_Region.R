## ------------------------------------------------------------------------------
## GRAB_Region.R
## Region / gene / set-based rare variant analysis workflow. Handles:
##   * Control processing for region-level parameters & kernels
##   * Parsing of user-supplied group (region) definition files
##   * Chunked region iteration with restart (index) support
##   * Integration with SKAT / SKAT-O and Cauchy combination p-values
##   * Method-specific region setup (currently POLMM-GENE)
##
## Functions:
##   GRAB.Region         : High-level API to execute region-based tests and persist
##                         both region- and marker-level outputs.
##   processOneRegion    : Process a single region: genotype reading, variant filtering,
##                         SKAT/SKAT-O tests, Cauchy combination, and file output.
## ------------------------------------------------------------------------------

#' Perform region-based association tests
#'
#' Tests for association between phenotypes and genomic regions containing multiple
#' genetic variants, primarily low-frequency and rare variants.
#'
#' @param objNull (S3 object) Null model object from \code{\link{GRAB.NullModel}}. 
#'   Currently supports POLMM_NULL_Model.
#' @param GenoFile (character) Path to genotype file (PLINK or BGEN format). See
#'   \code{\link{GRAB.ReadGeno}} for details.
#' @param OutputFile (character) Path for saving region-based association results.
#' @param GenoFileIndex (character or NULL) Index files for the genotype file. If 
#'   \code{NULL} (default), uses same prefix as \code{GenoFile}. See 
#'   \code{\link{GRAB.ReadGeno}} for details.
#' @param OutputFileIndex (character or NULL) Path for progress tracking file. If 
#'   \code{NULL} (default), uses \code{paste0(OutputFile, ".index")}.
#' @param GroupFile (character) Path to region definition file specifying region-marker 
#'   mappings and annotation information. Tab-separated format with 2-3 columns per region.
#' @param SparseGRMFile (character or NULL) Path to sparse GRM file (optional).
#' @param MaxMAFVec (character) Comma-separated MAF cutoffs for including variants in 
#'   analysis (default: "0.01,0.001,0.0005").
#' @param annoVec (character) Comma-separated annotation groups for analysis
#'   (default: "lof,lof:missense,lof:missense:synonymous").
#' @param control (list or NULL) List of the following parameters:
#'   \itemize{
#'     \item \code{impute_method} (character): Method for imputing missing genotypes: "mean", "minor", or "drop". Default: "minor".
#'     \item \code{missing_cutoff} (numeric): Exclude markers with missing rate > this value. Range: 0 to 0.5. Default: 0.15.
#'     \item \code{min_mac_region} (numeric): Minimum MAC threshold; markers with MAC < this value are treated as ultra-rare variants. Default: 5.
#'     \item \code{max_markers_region} (integer): Maximum number of markers allowed per region. Default: 100.
#'     \item \code{r.corr} (numeric vector): Rho parameters for SKAT-O test. Range: 0 to 1. Default: c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1).
#'     \item \code{weights.beta} (numeric vector): Beta distribution parameters for variant weights (length 2). Default: c(1, 25).
#'     \item \code{omp_num_threads} (integer): Number of OpenMP threads for parallel computation. Default: data.table::getDTthreads().
#'     \item \code{min_nMarker} (integer): Minimum number of markers required for region analysis. Default: 3.
#'     \item \code{SPA_Cutoff} (numeric): Z-score cutoff for calculating p-value by SPA. Default: 2.
#'   }
#' 
#' @return
#' The function returns \code{NULL} invisibly. Results are saved to four files:
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
#'   \item{Info}{Marker status ("Rare Variants" or "Ultra-Rare Variants").}
#'   \item{Anno}{Annotation from \code{GroupFile}.}
#'   \item{AltFreq, MAC, MAF}{Allele frequency, minor allele count, and minor allele frequency.}
#'   \item{MissingRate}{Proportion of missing genotypes.}
#'   \item{StatVec}{Score test statistic.}
#'   \item{altBetaVec, seBetaVec}{Effect size estimate and standard error.}
#'   \item{pval0Vec, pval1Vec}{Unadjusted and SPA-adjusted p-values.}
#' }
#'
#' **Other marker info** (\code{paste0(OutputFile, ".otherMarkerInfo")}) columns:
#' \describe{
#'   \item{Region}{Region identifier.}
#'   \item{ID}{Marker identifier.}
#'   \item{Annos}{Annotation from \code{GroupFile}.}
#'   \item{Info}{Reason for exclusion (e.g., "Ultra-Rare Variants", "Missing rate > cutoff").}
#'   \item{AltFreq, MAC, MAF, MissingRate}{Allele frequency and QC metrics.}
#' }
#'
#' **Burden test summary** (\code{paste0(OutputFile, ".infoBurdenNoWeight")}) columns:
#' \describe{
#'   \item{region}{Region identifier.}
#'   \item{anno}{Annotation type.}
#'   \item{max_maf}{Maximum MAF cutoff.}
#'   \item{sum}{Sum of genotypes.}
#'   \item{Stat}{Score test statistic.}
#'   \item{beta, se.beta}{Effect size estimate and standard error.}
#'   \item{pvalue}{P-value for burden test.}
#' }
#' 
#' @examples
#' objNullFile <- system.file("extdata", "objPOLMMnull.RData", package = "GRAB")
#' load(objNullFile) # load a an example object, obj.POLMM, from step 1
#'
#' OutputDir <- tempdir()
#' OutputFile <- file.path(OutputDir, "resultPOLMMregion1.txt")
#' GenoFile <- system.file("extdata", "simuPLINK_RV.bed", package = "GRAB")
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
#' data.table::fread(paste0(OutputFile, ".markerInfo"))
#' data.table::fread(paste0(OutputFile, ".otherMarkerInfo"))
#' data.table::fread(paste0(OutputFile, ".index"), sep = "\t", header = FALSE)
#'
GRAB.Region <- function(
  objNull,
  GenoFile,
  OutputFile,
  GenoFileIndex = NULL,
  OutputFileIndex = NULL,
  GroupFile,
  SparseGRMFile = NULL,
  MaxMAFVec = "0.01,0.001,0.0005",
  annoVec = "lof,lof:missense,lof:missense:synonymous",
  control = NULL
) {

  supported_classes <- c("POLMM_NULL_Model")

  # ========== Validate and configure parameters ==========

  # Validate objNull
  NullModelClass <- class(objNull)                            # character

  if (!NullModelClass %in% supported_classes) {
    stop(
      "class(objNull) should be one of: ",
      paste(paste0('"', supported_classes, '"'), collapse = ", ")
    )
  }

  if (any(!c("subjData", "N") %in% names(objNull))) {
    stop("c('subjData', 'N') should be in names(objNull).")
  }

  method <- gsub("_NULL_Model", "", NullModelClass)          # character

  # GenoFile and GenoFileIndex will be validated in GRAB.ReadGeno()

  # OutputFile and OutputFileIndex will be further validated in checkOutputFile()
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  # Validate GroupFile
  if (!file.exists(GroupFile)) {
    stop("Cannot find GroupFile: ", GroupFile)
  }

  # Validate SparseGRMFile
  if (!is.null(SparseGRMFile) && !file.exists(SparseGRMFile)) {
      stop("Cannot find SparseGRMFile: ", SparseGRMFile)
  }
  
  # Parse and validate MAF cutoffs for variant selection
  MaxMAFVec <- MaxMAFVec %>%                                  # numeric vector
    strsplit(split = ",") %>%
    unlist() %>%
    as.numeric()

  if (any(is.na(MaxMAFVec))) {
    stop("MaxMAFVec contains invalid (NA) values. Please check your input.")
  }

  MaxMAF <- max(MaxMAFVec)                                    # numeric
  if (MaxMAF > 0.05) {
    stop("Maximum value of 'MaxMAFVec' should be <= 0.05.")
  }

  # Parse annotation groups for variant filtering
  annoVec <- annoVec %>%                                      # character vector
    strsplit(split = ",") %>%
    unlist()
  annoList <- annoVec %>% strsplit(split = ":")               # list
  allAnno <- annoList %>%                                     # character vector
    unlist() %>%
    unique()

  # ========== Validate and configure the control list ==========

  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be an R list.")
  }

  default.region.control <- list(                             # list
    impute_method = "minor",
    missing_cutoff = 0.15,
    min_mac_region = 5,
    max_markers_region = 100,
    r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),
    weights.beta = c(1, 25),
    omp_num_threads = data.table::getDTthreads(),
    min_nMarker = 3
  )

  control <- updateControl(control, default.region.control)  # list

  if (!control$impute_method %in% c("mean", "minor", "drop")) {
    stop("control$impute_method should be 'mean', 'minor', or 'drop'.")
  }

  if (!is.numeric(control$missing_cutoff) ||
      control$missing_cutoff < 0 || control$missing_cutoff > 0.5) {
    stop("control$missing_cutoff should be numeric in [0, 0.5].")
  }

  if (!is.numeric(control$min_mac_region) || control$min_mac_region < 0) {
    stop("control$min_mac_region should be numeric >= 0.")
  }

  if (!is.numeric(control$max_markers_region) || control$max_markers_region < 50) {
    stop("control$max_markers_region should be integer >= 50.")
  }

  if (!is.numeric(control$r.corr) ||
      min(control$r.corr) < 0 || max(control$r.corr) > 1) {
    stop("control$r.corr should be numeric vector with elements in [0, 1].")
  }

  if (!is.numeric(control$weights.beta) ||
      length(control$weights.beta) != 2 || min(control$weights.beta) < 0) {
    stop("control$weights.beta should be numeric vector with two non-negative elements.")
  }

  if (!is.numeric(control$min_nMarker) || control$min_nMarker <= 0) {
    stop("control$min_nMarker should be a positive integer.")
  }

  # Validate method-specific control parameters
  control <- switch(                                          # list
    NullModelClass,
    POLMM_NULL_Model = checkControl.Region.POLMM(control)
  )

  # ========== Print all parameters ==========

  params <- list(
    Method = NullModelClass,
    `Genotype file` = GenoFile,
    `Genotype index file` = ifelse(is.null(GenoFileIndex), "Default", GenoFileIndex),
    `Output file` = OutputFile,
    `Output index file` = ifelse(is.null(OutputFileIndex), "Default", OutputFileIndex),
    `Group file` = GroupFile,
    `Sparse GRM file` = ifelse(is.null(SparseGRMFile), "NULL", SparseGRMFile),
    `Max MAF cutoffs` = paste(MaxMAFVec, collapse = ", "),
    `Annotation groups` = paste(annoVec, collapse = ", ")
  )
  .printParameters("Parameters for Region-Level Tests", params, control)

  # ========== Check output file status and determine restart point ==========

  indexChunk <- checkOutputFile(                              # integer
    OutputFile, OutputFileIndex, "Region", nEachChunk = 1
  )
  
  # ========== Extract subject information and process sample grouping ==========

  # Extract subject information from null model
  subjData <- as.character(objNull$subjData)                  # character vector
  n <- length(subjData)                                       # integer

  # Process sample grouping information for stratified analysis
  # Default: single group for all samples
  SampleLabelNumber <- rep(1, n)                              # numeric vector
  SampleLabelLevels <- NULL                                   # NULL or character vector
  nLabel <- max(SampleLabelNumber)                            # integer

  # ========== Extract region information and process variant grouping ==========

  .message("Extracting marker information from GroupFile: %s", GroupFile)

  gf <- file(GroupFile, "r")                                  # connection
  regionList <- list()                                        # list
  nLine <- 1                                                  # integer

  previousType <- "first"                                     # character
  previousGene <- "first"                                     # character
  Weights <- NA                                               # numeric or NA
  nRegion <- 1                                                # integer

  while (TRUE) {
    markerGroupLine <- readLines(gf, n = 1)                   # character

    if (length(markerGroupLine) == 0) {
      if (nRegion == 1) {
        stop("Cannot find any region information in 'GroupFile'.")
      }
      regionList[[nRegion]] <- list(                          # list
        regionID = previousGene,
        regionInfo = data.frame(
          ID = Markers,
          Annos = Annos,
          Weights = Weights
        )
      )
      close(gf)
      break
    }

    # Parse group file line
    if (length(markerGroupLine) == 0) {
      stop("The line ", nLine, " in `groupFile` is empty.")
    }

    info <- strsplit(markerGroupLine, "\t")[[1]]              # character vector
    if (length(info) < 3) {
      stop("The line ", nLine, " in 'groupFile' includes < 3 elements, ",
        "please note that each line should be seperated by 'tab'.")
    }

    geneID <- info[1]                                         # character
    type <- info[2]                                           # character
    values <- info[-c(1, 2)]                                  # character vector

    grepTemp <- grep(" ", values, value = TRUE)               # character vector
    if (length(grepTemp) > 0) {
      stop("'GroupFile' cannot contain 'space':\n",
        paste0(unique(grepTemp), collapse = "\t"))
    }

    grepTemp <- grep(";", values, value = TRUE)               # character vector
    if (length(grepTemp) > 0) {
      stop("'GroupFile' cannot contain ';':\n",
        paste0(unique(grepTemp), collapse = "\t"))
    }

    if (type == "weight") {
      values <- as.numeric(values)                            # numeric vector
    }

    n <- length(values)                                       # integer
    nLine <- nLine + 1

    if (!type %in% c("var", "anno", "weight")) {
      stop("The second column of the groupFile (tab-seperated) should be one of 'var', 'anno', and 'weight'.\n",
        "         Please double check line ", nLine, ".")
    }

    if (type == "var") {
      if (previousType == "var") {
        stop("Cannot find 'anno' line for region ", previousGene, ".")
      }
      if (previousType != "first") {
        regionList[[nRegion]] <- list(                        # list
          regionID = previousGene,
            regionInfo = data.frame(
              ID = Markers,
              Annos = Annos,
              Weights = Weights
            )
        )
        nRegion <- nRegion + 1
      }

      Markers <- values                                       # character vector
      n1 <- n                                                 # integer
      Weights <- NA                                           # NA
    }

    if (type == "anno") {
      if (n != n1) {
        stop("The length of annotations for markers is not equal to the length of marker IDs")
      }
      if (previousType != "var") {
        stop("In the 'GroupFile', the 'anno' line should follow the 'var' line.")
      }
      Annos <- values                                         # character vector
    }

    if (type == "weight") {
      if (n != n1) {
        stop("The length of weights for markers is not equal to the length of marker IDs")
      }
      if (previousType != "anno") {
        stop("In the 'GroupFile', the 'weight' line should follow the 'anno' line.")
      }
      Weights <- values                                       # numeric vector
    }

    previousType <- type
    previousGene <- geneID
  }

  .message("Found %d groups in GroupFile", nRegion)
  RegionList <- regionList                                    # list
  nRegions <- length(RegionList)                              # integer

  # ========== Configure C++ backend ==========
 
  # Initialize genotype reader in C++ backend
  control$max_maf_region <- MaxMAF                            # numeric
  objGeno <- setGenoInput(GenoFile, GenoFileIndex, subjData, control) # list
  genoType <- objGeno$genoType                                # character
  markerInfo <- objGeno$markerInfo                            # data.frame

  # Configure global variables in C++ backend
  with(
    control,
    setRegion_GlobalVarsInCPP(
      impute_method,
      missing_cutoff,
      max_maf_region,
      min_mac_region,
      max_markers_region,
      omp_num_threads,
      weights.beta,
      MaxMAFVec
    )
  )

  # Set method-specific objects in C++ backend
  obj.setRegion <- switch(                                    # list
    NullModelClass,
    POLMM_NULL_Model = setRegion.POLMM(objNull, control, SparseGRMFile)
  )

  # Use SKAT.Met_SKAT_Get_Pvalue instead of SKAT:::Met_SKAT_Get_Pvalue to be CRAN-compliant
  SKAT.Met_SKAT_Get_Pvalue <- getFromNamespace("Met_SKAT_Get_Pvalue", "SKAT") # function

  # ========== Iterate over regions to perform tests ==========

  for (i in (indexChunk + 1):nRegions) {
    processOneRegion(
      i = i,
      RegionList = RegionList,
      markerInfo = markerInfo,
      allAnno = allAnno,
      annoList = annoList,
      annoVec = annoVec,
      MaxMAFVec = MaxMAFVec,
      method = method,
      genoType = genoType,
      SampleLabelNumber = SampleLabelNumber,
      nLabel = nLabel,
      SampleLabelLevels = SampleLabelLevels,
      NullModelClass = NullModelClass,
      control = control,
      n = n,
      obj.setRegion = obj.setRegion,
      SKAT.Met_SKAT_Get_Pvalue = SKAT.Met_SKAT_Get_Pvalue,
      OutputFile = OutputFile,
      OutputFileIndex = OutputFileIndex,
      nRegions = nRegions
    )
  }

  .message("Analysis complete! Results saved to '%s'", OutputFile)
  return(invisible(NULL))
}


# Internal function to analyze one region: read genotypes, filter variants,
# compute SKAT/SKAT-O/Burden tests, apply Cauchy combination, and write results.
processOneRegion <- function(
  i,                          # Integer index of the current region
  RegionList,                 # List of regions with region IDs and marker information
  markerInfo,                 # Data frame with marker information from genotype file
  allAnno,                    # Character vector of all annotation types
  annoList,                   # List of annotation groups for filtering
  annoVec,                    # Character vector of annotation names
  MaxMAFVec,                  # Numeric vector of MAF cutoffs
  method,                     # Character string specifying the method (e.g., "POLMM")
  genoType,                   # Character string specifying genotype file type
  SampleLabelNumber,          # Numeric vector of sample group labels
  nLabel,                     # Integer number of sample groups
  SampleLabelLevels,          # Character vector of sample group names (or NULL)
  NullModelClass,             # Character string specifying null model class
  control,                    # List of control parameters
  n,                          # Integer number of samples
  obj.setRegion,              # List with method-specific region setup
  SKAT.Met_SKAT_Get_Pvalue,   # Function to compute SKAT p-values
  OutputFile,                 # Character path to main output file
  OutputFileIndex,            # Character path to index file for restart
  nRegions                    # Integer total number of regions
) {
  region <- RegionList[[i]]                                 # list

  regionID <- region$regionID                               # character
  regionInfo <- region$regionInfo                           # data.frame

  regionInfo <- markerInfo %>%                              # data.frame
    select(ID, genoIndex) %>%
    merge(regionInfo, by = "ID") %>%
    arrange(genoIndex) %>%
    filter(Annos %in% allAnno)

  nMarkers <- nrow(regionInfo)                              # integer

  if (nMarkers == 0) {
    return(invisible(NULL))
  }

  nAnno <- length(annoList)                                 # integer
  annoMat <- matrix(0, nrow = nMarkers, ncol = nAnno)       # matrix
  colnames(annoMat) <- annoVec

  for (iAnno in 1:nAnno) {
    annoMat[, iAnno] <- ifelse(regionInfo$Annos %in% annoList[[iAnno]], 1, 0) # numeric vector
  }

  genoIndex <- regionInfo$genoIndex                         # integer vector
  weightVec <- regionInfo$Weights                           # numeric vector or NA

  if (all(is.na(weightVec))) {
    weightVec <- rep(1, nMarkers)                           # numeric vector
  } else {
    if (any(is.na(weightVec) || weightVec <= 0)) {
      stop("Marker weights cannot be non-positive (<= 0) or NA.")
    }
  }

  .message("Analyzing region %s (%d/%d)", regionID, i, nRegions)
  .message(
    "Region contains %d markers: %s",
    length(regionInfo$ID),
    paste0(head(regionInfo$ID, 6), collapse = ", ")
  )

  obj.mainRegionInCPP <- mainRegionInCPP(                   # list
    method, genoType, genoIndex, weightVec, OutputFile,
    SampleLabelNumber, nLabel,
    annoMat, annoVec
  )

  # Updated on 2022-06-24: save sum of genotype to conduct burden test and adjust p-values using SPA
  pvalBurden <- obj.mainRegionInCPP$pvalBurden              # matrix

  # Updated on 2023-02-06: record summary statistics for sum of genotype for a region
  infoBurdenNoWeight <- obj.mainRegionInCPP$infoBurdenNoWeight # matrix
  infoBurdenNoWeight <- as.data.frame(infoBurdenNoWeight)   # data.frame
  infoBurdenNoWeight <- cbind(regionID, infoBurdenNoWeight) # data.frame
  colnames(infoBurdenNoWeight) <- c("region", "anno", "max_maf", "sum", "Stat", "beta", "se.beta", "pvalue")

  infoBurdenNoWeight$anno <- annoVec[infoBurdenNoWeight$anno + 1]
  infoBurdenNoWeight$max_maf <- MaxMAFVec[infoBurdenNoWeight$max_maf + 1]

  # Add annotation information
  obj.mainRegionInCPP$AnnoVec <- c(regionInfo$Annos, annoVec) # character vector
  if (!is.null(SampleLabelLevels)) {
    colnames(obj.mainRegionInCPP$MACLabelMat) <- paste0("MAC_", SampleLabelLevels)
    colnames(obj.mainRegionInCPP$MAFLabelMat) <- paste0("MAF_", SampleLabelLevels)
  }

  # Perform method-specific region analysis
  obj.mainRegion <- switch(                                 # list
    NullModelClass,
    POLMM_NULL_Model = mainRegion.POLMM(genoType, genoIndex, OutputFile, n,
                                         obj.setRegion, obj.mainRegionInCPP, nLabel),
    stop("Unknown NullModelClass: ", NullModelClass)
  )

  Other.Markers <- obj.mainRegion$Other.Markers %>%         # data.frame
    mutate(Region = regionID, .before = ID)
  VarMat <- obj.mainRegion$VarMat                           # matrix
  RV.Markers0 <- obj.mainRegion$RV.Markers %>%              # data.frame
    mutate(Region = regionID, .before = ID)

  Other.Markers <- regionInfo %>%                           # data.frame
    select(ID, Annos) %>%
    merge(Other.Markers, by = "ID")

  if (nrow(VarMat) != nrow(RV.Markers0)) {
    stop("nrow(VarMat) != nrow(RV.Markers0)!")
  }

  RV.Markers <- RV.Markers0 %>%                             # data.frame
    mutate(
      betaWeights = dbeta(MAF, control$weights.beta[1], control$weights.beta[2]),
      adjVarSVec = StatVec^2 / qchisq(pval1Vec, df = 1, lower.tail = FALSE),
      r0 = pmax(adjVarSVec / diag(VarMat), 1),
      wr0 = sqrt(r0) * betaWeights,
      wStatVec = StatVec * betaWeights
    )

  wr0 <- RV.Markers$wr0                                     # numeric vector

  wadjVarSMat <- t(VarMat * wr0) * wr0                      # matrix

  RV.MarkersWithAnno <- regionInfo %>%                      # data.frame
    select(-genoIndex) %>%
    merge(RV.Markers %>% select(ID, MAF, posRow), by = "ID")

  Other.MarkersWithAnno <- regionInfo %>%                   # data.frame
    select(ID, Annos) %>%
    merge(Other.Markers %>% filter(IndicatorVec == 2) %>% select(ID), by = "ID")

  RV.MarkersURV <- RV.Markers %>%                           # data.frame
    filter(Info == "Ultra-Rare Variants") %>%
    select(ID, posRow)

  pval.Region <- data.frame()                               # data.frame
  iSPA <- 1                                                 # integer
  for (anno in annoVec) {
    annoTemp <- unlist(strsplit(anno, split = ":"))        # character vector

    posURV <- RV.MarkersURV %>%                             # integer vector
      filter(ID == anno) %>%
      select(posRow) %>%
      unlist()
    nMarkersURV <- Other.MarkersWithAnno %>%                # integer
      filter(Annos %in% annoTemp) %>%
      nrow()
    if (length(posURV) != 1) {
      stop("length(posURV) != 1")
    }

    for (MaxMAF in MaxMAFVec) {
      posRV <- RV.MarkersWithAnno %>%                       # integer vector
        filter(MAF < MaxMAF & Annos %in% annoTemp) %>%
        select(posRow) %>%
        unlist()
      pos <- c(posRV, posURV)                               # integer vector
      n1 <- length(pos)                                     # integer

      ScoreBurden <- sum(RV.Markers$wStatVec[pos])          # numeric
      VarBurden <- sum(wadjVarSMat[pos, pos])               # numeric
      pvalBurdenSPA <- pvalBurden[iSPA, 2]                  # numeric
      VarBurdenSPA <- ScoreBurden^2 / qchisq(pvalBurdenSPA, df = 1, lower.tail = FALSE) # numeric
      ratioBurdenSPA <- max(VarBurdenSPA / VarBurden, 1)    # numeric
      iSPA <- iSPA + 1

      out_SKAT_List <- with(RV.Markers, try(              # list or try-error
        SKAT.Met_SKAT_Get_Pvalue(
          Score = wStatVec[pos],
          Phi = ratioBurdenSPA * wadjVarSMat[pos, pos],
          r.corr = control$r.corr,
          method = "optimal.adj",
          Score.Resampling = NULL
        ),
        silent = TRUE
      ))

      if (inherits(out_SKAT_List, "try-error")) {
        Pvalue <- c(NA, NA, NA)                             # numeric vector
      } else if (!any(c(0, 1) %in% out_SKAT_List$param$rho)) {
        Pvalue <- c(NA, NA, NA)                             # numeric vector
      } else {
        pos00 <- which(out_SKAT_List$param$rho == 0)        # integer
        pos01 <- which(out_SKAT_List$param$rho == 1)        # integer
        Pvalue <- c(                                        # numeric vector
          out_SKAT_List$p.value, # SKAT-O
          out_SKAT_List$param$p.val.each[pos00], # SKAT
          out_SKAT_List$param$p.val.each[pos01] # Burden Test
        )
      }

      pval.Region <- rbind.data.frame(                        # data.frame
        pval.Region,
        data.frame(
          Region = regionID,
          nMarkers = length(posRV),
          nMarkersURV = nMarkersURV,
          Anno.Type = anno,
          MaxMAF.Cutoff = MaxMAF,
          pval.SKATO = Pvalue[1],
          pval.SKAT = Pvalue[2],
          pval.Burden = Pvalue[3]
        )
      )
    }
  }

  # Cauchy combination test
  pval.Cauchy.SKATO <- CCT(pval.Region$pval.SKATO)          # numeric
  pval.Cauchy.SKAT <- CCT(pval.Region$pval.SKAT)            # numeric
  pval.Cauchy.Burden <- CCT(pval.Region$pval.Burden)        # numeric

  pval.Region <- rbind.data.frame(                          # data.frame
    pval.Region,
    data.frame(
      Region = regionID,
      nMarkers = NA,
      nMarkersURV = NA,
      Anno.Type = "Cauchy",
      MaxMAF.Cutoff = NA,
      pval.SKATO = pval.Cauchy.SKATO,
      pval.SKAT = pval.Cauchy.SKAT,
      pval.Burden = pval.Cauchy.Burden
    )
  )

  # Write chunk results to output files and update progress tracking
  writeOutputFile(
    Output = list(
      pval.Region,
      RV.Markers0,
      Other.Markers,
      infoBurdenNoWeight
    ),
    OutputFile = list(
      OutputFile,
      paste0(OutputFile, ".markerInfo"),
      paste0(OutputFile, ".otherMarkerInfo"),
      paste0(OutputFile, ".infoBurdenNoWeight")
    ),
    OutputFileIndex = OutputFileIndex,
    AnalysisType = "Region",
    nEachChunk = 1,
    indexChunk = i,
    Start = (i == 1),
    End = (i == nRegions)
  )

  return(invisible(NULL))
}
