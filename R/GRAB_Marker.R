#' Perform genome-wide association testing at the marker level
#'
#' Conducts single-marker association tests between genetic variants and phenotypes using
#' various statistical methods supported by GRAB.
#'
#' @param objNull Null model object from \code{\link{GRAB.NullModel}}.
#' @param GenoFile Path to genotype file. Supports PLINK (.bed/.bim/.fam) and BGEN formats.
#'   See \code{\link{GRAB.ReadGeno}} for format details.
#' @param GenoFileIndex Index files for the genotype file. If \code{NULL} (default), uses
#'   the same prefix as \code{GenoFile}. See \code{\link{GRAB.ReadGeno}} for details.
#' @param OutputFile Path for saving association test results.
#' @param OutputFileIndex Path for progress tracking file. Enables analysis restart if
#'   interrupted. If \code{NULL} (default), uses \code{paste0(OutputFile, ".index")}.
#' @param control List of control parameters. See \code{Details} for available options.
#' @details
#' GRAB supports multiple association testing methods: \code{POLMM}, \code{SPACox},
#' \code{SPAGRM}, \code{SPAmix}, and \code{WtCoxG}. The appropriate method is automatically
#' detected from \code{class(objNull)}. See \code{\link{GRAB.NullModel}} for method details.
#'
#' ## Control Parameters
#'
#' **Marker Selection**: Control which markers to analyze (default: all markers in file).
#' \itemize{
#'   \item \code{IDsToIncludeFile}, \code{IDsToExcludeFile}: Include/exclude specific markers.
#'   \item \code{RangesToIncludeFile}, \code{RangesToExcludeFile}: Include/exclude genomic ranges.
#'   \item \code{AlleleOrder}: "alt-first" (PLINK default) or "ref-first" (BGEN default).
#' }
#' See \code{\link{GRAB.ReadGeno}} for details on these parameters.
#'
#' **Quality Control**: Filter markers based on data quality.
#' \itemize{
#'   \item \code{ImputeMethod}: "mean" (default), "bestguess", or "drop".
#'   \item \code{MissingRateCutoff}: Exclude markers with missing rate > 0.15 (default).
#'   \item \code{MinMAFCutoff}: Exclude markers with MAF < 0.001 (default).
#'   \item \code{MinMACCutoff}: Exclude markers with MAC < 20 (default).
#'   \item \code{nMarkersEachChunk}: Process 10,000 markers per chunk (default).
#' }
#'
#' **Output Customization**: Control which columns appear in results.
#' \itemize{
#'   \item \code{outputColumns}: Additional columns beyond the standard set. Available options
#'     depend on the analysis method:
#'     \itemize{
#'       \item \code{POLMM}: Default: \code{beta}, \code{seBeta}; Optional: \code{zScore},
#'         \code{AltFreqInGroup}, \code{nSamplesInGroup}, \code{AltCountsInGroup}
#'       \item \code{SPACox}: Optional: \code{zScore}
#'     }
#' }
#' @return
#' Results are saved to \code{OutputFile} with the following columns:
#'
#' **Standard Columns** (included for all methods):
#' \describe{
#'   \item{Marker}{Marker identifiers from genotype files.}
#'   \item{Info}{Marker information as "CHR:POS:REF:ALT". REF/ALT order depends on
#'     \code{control$AlleleOrder}.}
#'   \item{AltFreq}{Alternative allele frequency before imputation. Values > 0.5 suggest
#'     reconsidering \code{control$AlleleOrder}.}
#'   \item{AltCounts}{Alternative allele counts before imputation.}
#'   \item{MissingRate}{Proportion of missing genotypes per marker.}
#'   \item{Pvalue}{Association test p-value.}
#' }
#'
#' **Optional Columns** (controlled by \code{control$outputColumns}):
#' \describe{
#'   \item{beta}{Effect size estimate for the alternative allele.}
#'   \item{seBeta}{Standard error of the effect size estimate.}
#'   \item{zScore}{Standardized test statistic (approximately standard normal).}
#'   \item{nSamplesInGroup}{Sample counts by phenotype group (may vary due to missing genotypes).}
#'   \item{AltCountsInGroup}{Alternative allele counts by phenotype group before imputation.}
#'   \item{AltFreqInGroup}{Alternative allele frequencies by phenotype group before imputation.}
#' }
#'
#' See \code{\link{makeGroup}} for phenotype grouping details used in group-specific statistics.
#'
#' @examples
#' # Load a precomputed POLMM_NULL_Model object to perform step 2 without repeating step 1
#' objNullFile <- system.file("extdata", "objPOLMMnull.RData", package = "GRAB")
#' load(objNullFile)
#' class(obj.POLMM)
#'
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputFile <- file.path(tempdir(), "simuOUTPUT.txt")
#' outputColumns <- c(
#'   "beta", "seBeta", "zScore",
#'   "nSamplesInGroup", "AltCountsInGroup", "AltFreqInGroup"
#' )
#'
#' GRAB.Marker(
#'   objNull = obj.POLMM,
#'   GenoFile = GenoFile,
#'   OutputFile = OutputFile,
#'   control = list(outputColumns = outputColumns)
#' )
#'
#' data.table::fread(OutputFile)
#'
GRAB.Marker <- function(
  objNull,
  GenoFile,
  GenoFileIndex = NULL,
  OutputFile,
  OutputFileIndex = NULL,
  control = NULL
) {
  # Check the validity of the null model object and extract its class (Util.R)
  # ---- BEGIN inlined: checkObjNull ----
  NullModelClass <- class(objNull)
  nm <- names(objNull)

  supported_classes <- c(
    "SPACox_NULL_Model",
    "POLMM_NULL_Model",
    "SPAmix_NULL_Model",
    "SPAGRM_NULL_Model",
    "SAGELD_NULL_Model",
    "WtCoxG_NULL_Model"
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

  # Set default output index file if not provided
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  # Validate and set control parameters with default values if not specified (control.R)
  checkControl.ReadGeno(control)
  # ---- start checkControl.Marker(control, NullModelClass) ----
  # check if control is an R list
  if (!is.null(control)) {
    if (!is.list(control)) {
      stop("If specified, the argument of 'control' should be an R 'list'.")
    }
  }

  # uniform default control setting for marker-level analysis
  default.marker.control <- list(
    impute_method = "mean",
    missing_cutoff = 0.15,
    min_maf_marker = 0.001,
    min_mac_marker = 20,
    nMarkersEachChunk = 10000,
    omp_num_threads = data.table::getDTthreads()
  ) # if 0, value is from omp_get_num_threads(). Not supported on 2022-02-07

  control <- updateControl(control, default.marker.control)

  # check if argument of 'control' is reasonable
  if (!control$impute_method %in% c("mean", "minor", "drop")) {
    stop("control$impute_method should be 'mean', 'minor', or 'drop'.")
  }

  if (!is.numeric(control$missing_cutoff) || control$missing_cutoff < 0 || control$missing_cutoff > 0.5) {
    stop("control$missing_cutoff should be a numeric value ranging from 0 to 0.5.")
  }

  if (!is.numeric(control$min_maf_marker) || control$min_maf_marker < 0 || control$min_maf_marker > 0.1) {
    stop("control$min_maf_marker should be a numeric value ranging from 0 to 0.1.")
  }

  if (!is.numeric(control$min_mac_marker) || control$min_mac_marker < 0 || control$min_mac_marker > 100) {
    stop("control$min_mac_marker should be a numeric value ranging from 0 to 100.")
  }

  if (!is.numeric(control$nMarkersEachChunk) || control$nMarkersEachChunk < 1e3 || 
      control$nMarkersEachChunk > 1e5) {
    stop("control$nMarkersEachChunk should be a numeric value ranging from 1e3 to 1e5.")
  }

  if (control$omp_num_threads < 0) {
    stop("control$omp_num_threads should be a positive integral value.")
  }

  if (NullModelClass == "POLMM_NULL_Model") {
    control <- checkControl.Marker.POLMM(control)
  } else if (NullModelClass == "SPACox_NULL_Model") {
    control <- checkControl.Marker.SPACox(control)
  } else if (NullModelClass == "SPAmix_NULL_Model") {
    control <- checkControl.Marker.SPAmix(control)
  } else if (NullModelClass == "SPAGRM_NULL_Model") {
    control <- checkControl.Marker.SPAGRM(control)
  } else if (NullModelClass == "SAGELD_NULL_Model") {
    control <- checkControl.Marker.SAGELD(control)
  } else if (NullModelClass == "WtCoxG_NULL_Model") {
    control <- checkControl.Marker.WtCoxG(control)
  } else {
    stop("Unknown NullModelClass: ", NullModelClass)
  }

  .message("Control parameters for marker-level association tests:")
  tmp <- capture.output(str(control))
  for (line in tmp) {
    if (startsWith(line, " $")) {
      message(sub("^ \\$", strrep(" ", 8), line))
    }
  }
  # ---- end checkControl.Marker(control, NullModelClass) ----

  nMarkersEachChunk <- control$nMarkersEachChunk

  # Check output file status and determine restart point if needed (Util.R)
  outList <- checkOutputFile(
    OutputFile, OutputFileIndex, "Marker",
    format(nMarkersEachChunk, scientific = FALSE)
  )

  # Special validation for SPAGRM and SAGELD methods: check MAF interval constraints
  # These methods require marker MAF to be within the interval used in null model fitting
  if (NullModelClass %in% c("SPAGRM_NULL_Model", "SAGELD_NULL_Model")) {
    if (length(objNull$MAF_interval) > 1) {
      if (control$min_maf_marker <= min(objNull$MAF_interval)) {
        stop(
          "min_maf_marker is out of MAF_interval. ",
          "Please reset min_maf_marker or check MAF_interval."
        )
      }
    }
  }

  # Extract restart information from output file checking
  indexChunk <- outList$indexChunk
  Start <- outList$Start
  End <- outList$End

  # Check if analysis has already been completed
  if (End) {
    stop(
      "Analysis completed in an earlier run. Results saved in '",
      OutputFile,
      "'. Use a different 'OutputFile' to restart analysis."
    )
  }

  # Check if analysis was partially completed and needs to restart
  if (!Start) {
    .message(
      "Part of analysis completed and saved in %s. Restarting from chunk %d",
      OutputFileIndex, indexChunk + 1
    )
  }

  # Extract subject IDs from null model object for genotype filtering
  subjData <- as.character(objNull$subjData)

  # Create phenotype groups for calculating group-specific allele statistics
  # ---- BEGIN inlined: makeGroup ----
  yVec <- objNull$yVec
  m1 <- length(unique(yVec))
  if (m1 <= 10) {
    Group <- as.numeric(as.factor(yVec)) - 1  # Groups 0 to (m1-1)
  } else {
    Group <- floor((rank(yVec, ties.method = "max") - 1) / length(yVec) * 10)  # Groups 0 to 9
  }
  # ---- END inlined: makeGroup ----
  ifOutGroup <- any(c("AltFreqInGroup", "AltCountsInGroup") %in% control$outputColumns)

  # Initialize genotype reader with file paths and subject filtering options
  objGeno <- setGenoInput(GenoFile, GenoFileIndex, subjData, control) # Function in 'Geno.R'
  genoType <- objGeno$genoType
  markerInfo <- objGeno$markerInfo
  CHROM <- markerInfo$CHROM
  genoIndex <- markerInfo$genoIndex

  # Split all markers into processing chunks by chromosome
  # Strategy: Group markers by chromosome, then divide each chromosome into fixed-size chunks
  # This maintains data locality and enables efficient chromosome-specific optimizations
  genoIndexList <- list()
  iTot <- 1

  # Process each chromosome separately to maintain genomic locality
  uCHROM <- unique(CHROM)
  for (chrom in uCHROM) {
    # Extract all markers belonging to current chromosome
    pos <- which(CHROM == chrom)
    gIdx <- genoIndex[pos]
    M <- length(gIdx)

    # Calculate chunk boundaries within this chromosome
    idxStart <- seq(1, M, nMarkersEachChunk)
    idxEnd <- idxStart + nMarkersEachChunk - 1

    nChunks <- length(idxStart)
    idxEnd[nChunks] <- M  # Ensure last chunk includes all remaining markers

    # Create individual chunks for this chromosome
    for (i in 1:nChunks) {
      idxMarker <- idxStart[i]:idxEnd[i]
      genoIndexList[[iTot]] <- list(
        chrom = chrom,                    # Chromosome identifier
        genoIndex = gIdx[idxMarker]      # Marker indices for this chunk
      )
      iTot <- iTot + 1
    }
  }

  nChunks <- length(genoIndexList)

  # Display analysis information to user
  .message("Number of markers to test: %d", nrow(markerInfo))
  .message("Number of markers in each chunk: %d", nMarkersEachChunk)
  .message("Number of chunks for all markers: %d", nChunks)

  # Initialize chromosome tracking for efficient method-specific setup
  chrom <- "InitialChunk"

  # Process each chunk of markers
  for (i in (indexChunk + 1):nChunks) {
    tempList <- genoIndexList[[i]]
    genoIndex <- tempList$genoIndex
    tempChrom <- tempList$chrom

    # Set up method-specific objects when chromosome changes
    # This optimization avoids redundant setup for markers within the same chromosome
    if (tempChrom != chrom) {
      # Configure global variables in C++ for efficient marker processing
      # See Main.cpp for implementation details
      nGroup <- length(unique(Group))
      setMarker_GlobalVarsInCPP(
        control$impute_method,
        control$missing_cutoff,
        control$min_maf_marker,
        control$min_mac_marker,
        control$omp_num_threads,
        Group, ifOutGroup, nGroup
      )

      # Initialize method-specific objects based on the null model class
      # Each method has its own setup requirements and computational optimizations
      if (NullModelClass == "POLMM_NULL_Model") {
        setMarker.POLMM(objNull, control, chrom)
      } else if (NullModelClass == "SPACox_NULL_Model") {
        setMarker.SPACox(objNull, control)
      } else if (NullModelClass == "SPAmix_NULL_Model") {
        setMarker.SPAmix(objNull, control)
      } else if (NullModelClass == "SPAGRM_NULL_Model") {
        setMarker.SPAGRM(objNull, control)
      } else if (NullModelClass == "SAGELD_NULL_Model") {
        setMarker.SAGELD(objNull, control)
      } else if (NullModelClass == "WtCoxG_NULL_Model") {
        setMarker.WtCoxG(objNull, control)
      }

      chrom <- tempChrom
    }

    .message("---- Analyzing Chunk %d/%d: chrom %s ----", i, nChunks, chrom)

    # Perform association testing for all markers in current chunk
    # Method selection based on null model class determines the appropriate algorithm
    if (NullModelClass == "POLMM_NULL_Model") {
      # POLMM: Proportional odds logistic mixed model for ordinal traits
      # Call optimized C++ implementation for computational efficiency
      OutList <- mainMarkerInCPP("POLMM", genoType, genoIndex)

      # Construct base output data frame with required columns
      resMarker <- data.frame(
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
      additionalColumns <- intersect(optionalColumns, control$outputColumns)

      if (length(additionalColumns) > 0) {
        resMarker <- cbind.data.frame(
          resMarker,
          as.data.frame(OutList[additionalColumns])
        )
      }
    } else if (NullModelClass == "SPACox_NULL_Model") {
      resMarker <- mainMarker.SPACox(genoType, genoIndex, control$outputColumns)
    } else if (NullModelClass == "SPAmix_NULL_Model") {
      resMarker <- mainMarker.SPAmix(genoType, genoIndex, control$outputColumns, objNull)
    } else if (NullModelClass == "SPAGRM_NULL_Model") {
      resMarker <- mainMarker.SPAGRM(genoType, genoIndex, control$outputColumns)
    } else if (NullModelClass == "SAGELD_NULL_Model") {
      resMarker <- mainMarker.SAGELD(genoType, genoIndex, control$outputColumns, objNull)
    } else if (NullModelClass == "WtCoxG_NULL_Model") {
      resMarker <- mainMarker.WtCoxG(genoType, genoIndex, control, objNull)
    }

    # Write chunk results to output file and update progress tracking
    writeOutputFile(
      Output = list(resMarker),
      OutputFile = list(OutputFile),
      OutputFileIndex = OutputFileIndex,
      AnalysisType = "Marker",
      nEachChunk = format(nMarkersEachChunk, scientific = FALSE),
      indexChunk = i,
      Start = (i == 1),
      End = (i == nChunks)
    )
  }

  .message("Analysis complete! Results saved to '%s'", OutputFile)
  return(invisible(NULL))
}