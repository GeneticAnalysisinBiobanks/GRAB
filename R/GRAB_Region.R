#' Perform region-based genetic association testing
#'
#' Tests for association between phenotypes and genomic regions containing multiple
#' genetic variants, primarily low-frequency and rare variants.
#'
#' @param objNull Null model object from \code{\link{GRAB.NullModel}}.
#' @param GenoFile Path to genotype file (PLINK or BGEN format). See
#'   \code{\link{GRAB.ReadGeno}} for details.
#' @param GenoFileIndex Index files for genotype file. If \code{NULL} (default), uses
#'   same prefix as \code{GenoFile}. See \code{\link{GRAB.ReadGeno}} for details.
#' @param OutputFile Path for saving region-based association results.
#' @param OutputFileIndex Path for progress tracking file. If \code{NULL} (default),
#'   uses \code{paste0(OutputFile, ".index")}.
#' @param GroupFile Path to region definition file specifying region-marker mappings
#'   and annotation information. Tab-separated format with 2-3 columns per region.
#' @param SparseGRMFile Path to sparse GRM file (optional).
#' @param SampleFile Path to sample information file with header (optional).
#' @param MaxMAFVec Comma-separated MAF cutoffs for including variants in analysis
#'   (default: "0.01,0.001,0.0005").
#' @param annoVec Comma-separated annotation groups for analysis
#'   (default: "lof,lof:missense,lof:missense:synonymous").
#' @param chrom Chromosome-specific options (default: "LOCO=F").
#' @param control List of control parameters. See \code{Details} for options.
#' @details
#' GRAB supports region-based testing using multiple methods: \code{POLMM}, \code{SPACox},
#' \code{SPAGRM}, \code{SPAmix}, and \code{WtCoxG}. The method is automatically detected
#' from \code{class(objNull)}. See \code{\link{GRAB.NullModel}} for method details.
#'
#' ## Control Parameters
#'
#' **Genotype Processing:**
#' \itemize{
#'   \item \code{AlleleOrder}: "alt-first" (PLINK default) or "ref-first" (BGEN default).
#'   \item \code{ImputeMethod}: "mean", "bestguess" (default), or "drop".
#'   \item \code{omp_num_threads}: Number of OpenMP threads for parallel computation.
#' }
#'
#' **Quality Control:**
#' \itemize{
#'   \item \code{MissingRateCutoff}: Exclude markers with missing rate > 0.15 (default).
#'   \item \code{MinMACCutoff}: Treat markers with MAC < 5 (default) as ultra-rare variants.
#'   \item \code{nRegionsEachChunk}: Number of regions per output chunk (default: 1).
#' }
#'
#' **Kernel-Based Testing (SKAT/SKAT-O):**
#' \itemize{
#'   \item \code{kernel}: Kernel type (default: "linear.weighted").
#'   \item \code{weights_beta}: Beta weight parameters (default: c(1, 25)).
#'   \item \code{weights}: Custom weight vector (overrides \code{weights_beta} if specified).
#'   \item \code{r.corr}: Rho parameters for SKAT-O (default: c(0, 0.1^2, 0.2^2, 0.3^2,
#'     0.4^2, 0.5^2, 0.5, 1)).
#' }
#'
#' **Output Customization:**
#' \itemize{
#'   \item \code{outputColumns}: Additional columns for marker-level output.
#'     \itemize{
#'       \item \code{POLMM}: Default: \code{beta}, \code{seBeta}; Optional: \code{zScore},
#'         \code{AltFreqInGroup}, \code{nSamplesInGroup}, \code{AltCountsInGroup}
#'       \item \code{SPACox}: Optional: \code{zScore}
#'     }
#' }
#' See \code{\link{GRAB.ReadGeno}} for genotype processing details.
#' @return
#' Results are saved to two files:
#' \enumerate{
#'   \item \code{OutputFile}: Region-based test results
#'   \item \code{paste0(OutputFile, ".markerInfo")}: Marker-level results (same format as
#'     \code{\link{GRAB.Marker}})
#' }
#'
#' **Region-level results** (\code{OutputFile}) contain:
#' \describe{
#'   \item{Region}{Region identifiers from \code{GroupFile}.}
#'   \item{Anno.Type}{Annotation type from \code{GroupFile}.}
#'   \item{maxMAF}{Maximum MAF cutoff used for variant selection.}
#'   \item{nSamples}{Number of samples in analysis.}
#'   \item{nMarkers}{Number of markers with MAF < cutoff and MAC > \code{MinMACCutoff}.}
#'   \item{nMarkersURV}{Number of ultra-rare variants with MAC < \code{MinMACCutoff}.}
#'   \item{pval.SKATO}{SKAT-O test p-values.}
#'   \item{pval.SKAT}{SKAT test p-values.}
#'   \item{pval.Burden}{Burden test p-values.}
#' }
#' @examples
#' # Load a precomputed example object to perform step 2 without repeating step 1
#' objNullFile <- system.file("extdata", "objPOLMMnull.RData", package = "GRAB")
#' load(objNullFile)
#' class(obj.POLMM) # "POLMM_NULL_Model" is an object from POLMM method.
#'
#' OutputDir <- tempdir()
#' OutputFile <- file.path(OutputDir, "simuRegionOutput.txt")
#' GenoFile <- system.file("extdata", "simuPLINK_RV.bed", package = "GRAB")
#' GroupFile <- system.file("extdata", "simuPLINK_RV.group", package = "GRAB")
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#'
#' GRAB.Region(
#'   objNull = obj.POLMM,
#'   GenoFile = GenoFile,
#'   OutputFile = OutputFile,
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
  GenoFileIndex = NULL,
  OutputFile,
  OutputFileIndex = NULL,
  GroupFile,
  SparseGRMFile = NULL,
  SampleFile = NULL,
  MaxMAFVec = "0.01,0.001,0.0005",
  annoVec = "lof,lof:missense,lof:missense:synonymous",
  chrom = "LOCO=F",
  control = NULL
) {
  # Validate null model object and extract method information
  # ---- BEGIN inlined: checkObjNull ----
  NullModelClass <- class(objNull)
  nm <- names(objNull)

  supported_classes <- c(
    "POLMM_NULL_Model"
  )

  if (!NullModelClass %in% supported_classes) {
    stop(
      "class(objNull) should be one of: ",
      paste(paste0('"', supported_classes, '"'), collapse = ", ")
    )
  }

  if (any(!c("subjData", "N") %in% nm)) {
    stop("c('subjData', 'N') should be in names(objNull).")
  }
  # ---- END inlined: checkObjNull ----
  method <- gsub("_NULL_Model", "", NullModelClass)

  # Set default output index file if not provided
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  # Check output file status and determine restart point if needed
  outList <- checkOutputFile(OutputFile, OutputFileIndex, "Region", nEachChunk = 1)
  indexChunk <- outList$indexChunk
  Start <- outList$Start
  End <- outList$End

  # Check if analysis has already been completed
  if (End) {
    msg_text <- paste0(
      "Analysis completed in earlier run. Results saved in '", OutputFile, "'. ",
      "Use a different 'OutputFile' to restart analysis."
    )
    .message("%s", msg_text)
    return(msg_text)
  }

  # Check if analysis was partially completed and needs restart
  if (!Start) {
    msg_text <- paste0(
      "Partial analysis completed based on index file: ", OutputFileIndex, "\n",
      "Restarting from chunk ", indexChunk + 1
    )
    .message("%s", msg_text)
  }

  # Validate and set control parameters with method-specific defaults
  # ---- Begin inline of checkControl.Region(control) ----
  if (!is.null(control)) {
    if (!is.list(control)) {
      stop("If specified, the argument of 'control' should be an R 'list'.")
    }
  }

  # uniform default control setting for region-level analysis
  default.region.control <- list(
    impute_method = "minor",
    missing_cutoff = 0.15,
    min_mac_region = 5,
    max_markers_region = 100,
    r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),
    weights.beta = c(1, 25),
    omp_num_threads = data.table::getDTthreads(),
    min_nMarker = 3
  )

  control <- updateControl(control, default.region.control)

  # check if 'control' is reasonable
  if (!control$impute_method %in% c("mean", "minor", "drop")) {
    stop("control$impute_method should be 'mean', 'minor', or 'drop'.")
  }

  if (!is.numeric(control$missing_cutoff) || control$missing_cutoff < 0 || control$missing_cutoff > 0.5) {
    stop("control$missing_cutoff should be a numeric value ranging from 0 to 0.5.")
  }

  if (!is.numeric(control$min_mac_region) || control$min_mac_region < 0) {
    stop("control$min_mac_region should be a numeric value >= 0.")
  }

  if (!is.numeric(control$max_markers_region) || control$max_markers_region < 50) {
    stop("control$max_markers_region should be a integer >= 50.")
  }

  if (!is.numeric(control$r.corr) || min(control$r.corr) < 0 || max(control$r.corr) > 1) {
    stop("control$r.corr should be a numeric vector whose elements are between 0 and 1.")
  }

  if (!is.numeric(control$weights.beta) || length(control$weights.beta) != 2 || min(control$weights.beta) < 0) {
    stop("control$weights.beta should be a numeric vector with two non-negative elements.")
  }

  if (!is.numeric(control$min_nMarker) || control$min_nMarker <= 0) {
    stop("control$min_nMarker should be a positive integer.")
  }
  # ---- End inline of checkControl.Region(control) ----

  if (NullModelClass == "POLMM_NULL_Model") {
    control <- checkControl.Region.POLMM(control)
  } else {
    stop("Unknown NullModelClass: ", NullModelClass)
  }

  # Display control parameters for user verification
  .message("Control parameters for region-level genetic association analysis:")
  tmp <- capture.output(str(control))
  for (line in tmp) {
    if (startsWith(line, " $")) {
      message(sub("^ \\$", strrep(" ", 8), line))
    }
  }

  # Parse and validate MAF cutoffs for variant selection
  MaxMAFVec <- MaxMAFVec %>%
    strsplit(split = ",") %>%
    unlist() %>%
    as.numeric()

  if (any(is.na(MaxMAFVec))) {
    stop("MaxMAFVec contains invalid (NA) values. Please check your input.")
  }

  MaxMAF <- max(MaxMAFVec)
  if (MaxMAF > 0.05) {
    stop("Maximum value of 'MaxMAFVec' should be <= 0.05.")
  }
  control$max_maf_region <- MaxMAF

  # Parse annotation groups for variant filtering
  annoVec <- annoVec %>%
    strsplit(split = ",") %>%
    unlist()
  annoList <- annoVec %>% strsplit(split = ":")
  allAnno <- annoList %>%
    unlist() %>%
    unique()

  # Extract subject information from null model
  subjData <- as.character(objNull$subjData)
  n <- length(subjData)

  # Process sample grouping information for stratified analysis
  # Default: single group for all samples
  SampleLabelNumber <- rep(1, n)
  SampleLabelLevels <- NULL

  if (!is.null(SampleFile)) {
    # Read sample information file with required 'IID' column
    SampleInfo <- data.table::fread(SampleFile)
    if (colnames(SampleInfo)[1] != "IID") {
      stop("The header of the first column in 'SampleFile' should be 'IID'.")
    }

    # Validate that all subjects in null model are present in sample file
    pos <- which(!subjData %in% SampleInfo$IID)
    if (length(pos) > 0) {
      stop(
        "At least one subject in null model fitting not found in 'SampleFile':\n",
        paste0(subjData[pos], collapse = "\t")
      )
    }

    # Extract sample group labels if specified
    if (!is.null(control$SampleLabelCol)) {
      SampleLabelColName <- control$SampleLabelCol
      if (!SampleLabelColName %in% colnames(SampleInfo)) {
        stop("'SampleFile' should include column: ", SampleLabelColName)
      }

      # Map sample labels to numeric codes for stratified analysis
      posInSampleInfo <- match(subjData, SampleInfo$IID)
      SampleLabel <- SampleInfo[[SampleLabelColName]][posInSampleInfo]
      SampleLabelFactor <- as.factor(SampleLabel)
      SampleLabelNumber <- as.numeric(SampleLabelFactor)
      SampleLabelLevels <- levels(SampleLabelFactor)
    }
  }
  nLabel <- max(SampleLabelNumber)

  # Initialize genotype reader with file paths and subject filtering options
  objGeno <- setGenoInput(GenoFile, GenoFileIndex, subjData, control) # Function in 'Geno.R'
  genoType <- objGeno$genoType
  markerInfo <- objGeno$markerInfo

  # Parse region definitions from group file
  # ---- BEGIN inlined: getInfoGroupFile ----
  .message("Extracting marker information from GroupFile: %s", GroupFile)

  if (!file.exists(GroupFile)) {
    stop("cannot find the below file:\n", GroupFile)
  }

  gf <- file(GroupFile, "r")
  regionList <- list()
  nLine <- 1

  previousType <- "first"
  previousGene <- "first"
  Weights <- NA
  nRegion <- 1

  while (TRUE) {
    markerGroupLine <- readLines(gf, n = 1)

    if (length(markerGroupLine) == 0) {
      if (nRegion == 1) {
        stop("Cannot find any region information in 'GroupFile'.")
      }
      regionList[[nRegion]] <- list(
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

    # ---- inline getInfoGroupLine() ----
    if (length(markerGroupLine) == 0) {
      stop("The line ", nLine, " in `groupFile` is empty.")
    }

    info <- strsplit(markerGroupLine, "\t")[[1]]
    if (length(info) < 3) {
      stop("The line ", nLine, " in 'groupFile' includes < 3 elements, ",
        "please note that each line should be seperated by 'tab'.")
    }

    geneID <- info[1]
    type <- info[2]
    values <- info[-c(1, 2)]

    grepTemp <- grep(" ", values, value = TRUE)
    if (length(grepTemp) > 0) {
      stop("'GroupFile' cannot contain 'space':\n",
        paste0(unique(grepTemp), collapse = "\t"))
    }

    grepTemp <- grep(";", values, value = TRUE)
    if (length(grepTemp) > 0) {
      stop("'GroupFile' cannot contain ';':\n",
        paste0(unique(grepTemp), collapse = "\t"))
    }

    if (type == "weight") {
      values <- as.numeric(values)
    }

    n <- length(values)
    # ---- end inline getInfoGroupLine() ----
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
        regionList[[nRegion]] <- list(
          regionID = previousGene,
            regionInfo = data.frame(
              ID = Markers,
              Annos = Annos,
              Weights = Weights
            )
        )
        nRegion <- nRegion + 1
      }

      Markers <- values
      n1 <- n
      Weights <- NA
    }

    if (type == "anno") {
      if (n != n1) {
        stop("The length of annotations for markers is not equal to the length of marker IDs")
      }
      if (previousType != "var") {
        stop("In the 'GroupFile', the 'anno' line should follow the 'var' line.")
      }
      Annos <- values
    }

    if (type == "weight") {
      if (n != n1) {
        stop("The length of weights for markers is not equal to the length of marker IDs")
      }
      if (previousType != "anno") {
        stop("In the 'GroupFile', the 'weight' line should follow the 'anno' line.")
      }
      Weights <- values
    }

    previousType <- type
    previousGene <- geneID
  }

  .message("Found %d groups in GroupFile", nRegion)
  RegionList <- regionList
  # ---- END inlined: getInfoGroupFile ----
  nRegions <- length(RegionList)

  # Configure global variables in C++ for efficient region-based analysis
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

  if (NullModelClass == "POLMM_NULL_Model") {
    obj.setRegion <- setRegion.POLMM(objNull, control, chrom, SparseGRMFile)
  } else {
    stop("Unknown NullModelClass: ", NullModelClass)
  }

  diffTime1 <- 0
  diffTime2 <- 0
  diffTime3 <- 0

  # Use SKAT.Met_SKAT_Get_Pvalue instead of SKAT:::Met_SKAT_Get_Pvalue to be CRAN-compliant
  SKAT.Met_SKAT_Get_Pvalue <- getFromNamespace("Met_SKAT_Get_Pvalue", "SKAT")

  for (i in (indexChunk + 1):nRegions) {
    region <- RegionList[[i]]

    regionID <- region$regionID
    regionInfo <- region$regionInfo

    regionInfo <- markerInfo %>%
      select(ID, genoIndex) %>%
      merge(regionInfo, by = "ID") %>%
      arrange(genoIndex) %>%
      filter(Annos %in% allAnno)

    nMarkers <- nrow(regionInfo)

    if (nMarkers == 0) {
      next
    }

    nAnno <- length(annoList)
    annoMat <- matrix(0, nrow = nMarkers, ncol = nAnno)
    colnames(annoMat) <- annoVec

    for (iAnno in 1:nAnno) {
      annoMat[, iAnno] <- ifelse(regionInfo$Annos %in% annoList[[iAnno]], 1, 0)
    }

    genoIndex <- regionInfo$genoIndex
    weightVec <- regionInfo$Weights

    if (all(is.na(weightVec))) {
      weightVec <- rep(1, nMarkers)
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

    t11 <- Sys.time()
    obj.mainRegionInCPP <- mainRegionInCPP(
      method, genoType, genoIndex, weightVec, OutputFile,
      SampleLabelNumber, nLabel,
      annoMat, annoVec
    )
    t12 <- Sys.time()
    diffTime1 <- diffTime1 + (t12 - t11)

    # updated on 2022-06-24 (save sum of genotype to conduct burden test and adjust p-values using SPA)
    pvalBurden <- obj.mainRegionInCPP$pvalBurden

    # updated on 2023-02-06 (record summary statistics for sum of genotype for a region)
    infoBurdenNoWeight <- obj.mainRegionInCPP$infoBurdenNoWeight
    infoBurdenNoWeight <- as.data.frame(infoBurdenNoWeight)
    infoBurdenNoWeight <- cbind(regionID, infoBurdenNoWeight)
    colnames(infoBurdenNoWeight) <- c("region", "anno", "max_maf", "sum", "Stat", "beta", "se.beta", "pvalue")

    infoBurdenNoWeight$anno <- annoVec[infoBurdenNoWeight$anno + 1]
    infoBurdenNoWeight$max_maf <- MaxMAFVec[infoBurdenNoWeight$max_maf + 1]

    ## add annotation information
    obj.mainRegionInCPP$AnnoVec <- c(regionInfo$Annos, annoVec)
    if (!is.null(SampleLabelLevels)) {
      colnames(obj.mainRegionInCPP$MACLabelMat) <- paste0("MAC_", SampleLabelLevels)
      colnames(obj.mainRegionInCPP$MAFLabelMat) <- paste0("MAF_", SampleLabelLevels)
    }

    if (NullModelClass == "POLMM_NULL_Model") {
      obj.mainRegion <- mainRegion.POLMM(genoType, genoIndex, OutputFile, control, n, 
                                         obj.setRegion, obj.mainRegionInCPP, nLabel)
    } else {
      stop("Unknown NullModelClass: ", NullModelClass)
    }

    Other.Markers <- obj.mainRegion$Other.Markers %>% mutate(Region = regionID, .before = ID)
    VarMat <- obj.mainRegion$VarMat
    RV.Markers0 <- obj.mainRegion$RV.Markers %>% mutate(Region = regionID, .before = ID)

    Other.Markers <- regionInfo %>%
      select(ID, Annos) %>%
      merge(Other.Markers, by = "ID")

    if (nrow(VarMat) != nrow(RV.Markers0)) {
      stop("nrow(VarMat) != nrow(RV.Markers0)!")
    }

    RV.Markers <- RV.Markers0 %>%
      mutate(
        betaWeights = dbeta(MAF, control$weights.beta[1], control$weights.beta[2]),
        adjVarSVec = StatVec^2 / qchisq(pval1Vec, df = 1, lower.tail = FALSE),
        # r0 = adjVarSVec / diag(VarMat),  # edited on 06/22/2022
        r0 = pmax(adjVarSVec / diag(VarMat), 1),
        wr0 = sqrt(r0) * betaWeights,
        wStatVec = StatVec * betaWeights
      )

    # check given weights version later: 2022-05-01

    wr0 <- RV.Markers$wr0

    wadjVarSMat <- t(VarMat * wr0) * wr0

    RV.MarkersWithAnno <- regionInfo %>%
      select(-genoIndex) %>%
      merge(RV.Markers %>% select(ID, MAF, posRow), by = "ID")

    Other.MarkersWithAnno <- regionInfo %>%
      select(ID, Annos) %>%
      merge(Other.Markers %>% filter(IndicatorVec == 2) %>% select(ID), by = "ID")

    RV.MarkersURV <- RV.Markers %>%
      filter(Info == "Ultra-Rare Variants") %>%
      select(ID, posRow)

    t21 <- Sys.time()
    pval.Region <- data.frame()
    iSPA <- 1
    for (anno in annoVec) {
      annoTemp <- unlist(strsplit(anno, split = ":"))

      posURV <- RV.MarkersURV %>%
        filter(ID == anno) %>%
        select(posRow) %>%
        unlist()
      nMarkersURV <- Other.MarkersWithAnno %>%
        filter(Annos %in% annoTemp) %>%
        nrow()
      if (length(posURV) != 1) {
        stop("length(posURV) != 1")
      }

      for (MaxMAF in MaxMAFVec) {
        posRV <- RV.MarkersWithAnno %>%
          filter(MAF < MaxMAF & Annos %in% annoTemp) %>%
          select(posRow) %>%
          unlist()
        pos <- c(posRV, posURV)
        n1 <- length(pos)

        ScoreBurden <- sum(RV.Markers$wStatVec[pos])
        VarBurden <- sum(wadjVarSMat[pos, pos])
        pvalBurdenSPA <- pvalBurden[iSPA, 2]
        VarBurdenSPA <- ScoreBurden^2 / qchisq(pvalBurdenSPA, df = 1, lower.tail = FALSE)
        ratioBurdenSPA <- max(VarBurdenSPA / VarBurden, 1)
        iSPA <- iSPA + 1

        t31 <- Sys.time()
        out_SKAT_List <- with(RV.Markers, try(
          SKAT.Met_SKAT_Get_Pvalue(
            Score = wStatVec[pos],
            Phi = ratioBurdenSPA * wadjVarSMat[pos, pos],
            r.corr = control$r.corr,
            method = "optimal.adj",
            Score.Resampling = NULL
          ),
          silent = TRUE
        ))

        t32 <- Sys.time()
        diffTime3 <- diffTime3 + (t32 - t31)

        if (inherits(out_SKAT_List, "try-error")) {
          Pvalue <- c(NA, NA, NA)
        } else if (!any(c(0, 1) %in% out_SKAT_List$param$rho)) {
          Pvalue <- c(NA, NA, NA)
        } else {
          pos00 <- which(out_SKAT_List$param$rho == 0)
          pos01 <- which(out_SKAT_List$param$rho == 1)
          Pvalue <- c(
            out_SKAT_List$p.value, # SKAT-O
            out_SKAT_List$param$p.val.each[pos00], # SKAT
            out_SKAT_List$param$p.val.each[pos01]
          ) # Burden Test
        }

        pval.Region <- rbind.data.frame(
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

    ## Cauchy Combination
    pval.Cauchy.SKATO <- CCT(pval.Region$pval.SKATO)
    pval.Cauchy.SKAT <- CCT(pval.Region$pval.SKAT)
    pval.Cauchy.Burden <- CCT(pval.Region$pval.Burden)

    pval.Region <- rbind.data.frame(
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

    t22 <- Sys.time()
    diffTime2 <- diffTime2 + (t22 - t21)

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
  }

  .message("Region analysis timing: mainRegionInCPP %.2f seconds, SKATO %.2f seconds", diffTime1, diffTime3)

  .message(
    "Analysis complete! Results saved to:\n    %s\n    %s\n    %s",
    OutputFile,
    paste0(OutputFile, ".markerInfo"),
    paste0(OutputFile, ".otherMarkerInfo")
  )

  return(invisible(NULL))
}
