#' Conduct marker-level genetic association testing
#'
#' Performs GWAS between a trait and individual genetic markers.
#'
#' @param objNull The output object from function \code{\link{GRAB.NullModel}}.
#' @param GenoFile A character string specifying the genotype file path. Currently, two
#'   genotype formats are supported: PLINK and BGEN. See \code{\link{GRAB.ReadGeno}} for details.
#' @param GenoFileIndex Additional index files corresponding to \code{GenoFile}. If \code{NULL}
#'   (default), the same prefix as GenoFile is used. See \code{\link{GRAB.ReadGeno}} for details.
#' @param OutputFile A character string specifying the output file path to save analysis results.
#' @param OutputFileIndex A character string specifying the output index file to record the
#'   progress. If the program terminates unexpectedly, this helps \code{GRAB} understand where
#'   to restart the analysis. If \code{NULL} (default), 
#'   \code{OutputFileIndex = paste0(OutputFile, ".index")}.
#' @param control A list of parameters for controlling \code{GRAB.Marker} function behavior.
#'   See the \code{Details} section for more information.
#' @details
#' The \code{GRAB} package supports multiple statistical methods: \code{POLMM}, \code{SPACox},
#' \code{SPAGRM}, \code{SPAmix}, and \code{WtCoxG}.
#' Detailed information about these analysis methods is provided in the \code{Details} section
#' of \code{\link{GRAB.NullModel}}.
#' Users do not need to specify the method explicitly since \code{GRAB.Marker} and
#' \code{\link{GRAB.Region}} automatically detect it from \code{class(objNull)}.
#'
#' ## Control Parameters
#' The following parameters allow users to customize which markers to include in the analysis.
#' If these parameters are not specified, \code{GRAB} will analyze all markers in the file.
#' For PLINK files, the default is \code{control$AlleleOrder = "alt-first"};
#' for BGEN files, the default is \code{control$AlleleOrder = "ref-first"}.
#'   \itemize{
#'   \item \code{IDsToIncludeFile}: See the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{IDsToExcludeFile}: See the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{RangesToIncludeFile}: See the \code{Details} section of 
#'     \code{\link{GRAB.ReadGeno}}.
#'   \item \code{RangesToExcludeFile}: See the \code{Details} section of 
#'     \code{\link{GRAB.ReadGeno}}.
#'   \item \code{AlleleOrder}: See the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   }
#' 
#' The following parameters customize the quality control (QC) process:
#'   \itemize{
#'   \item \code{ImputeMethod}: A character string specifying imputation method: "mean"
#'     (default), "bestguess", or "drop". See the \code{Details} section of 
#'     \code{\link{GRAB.ReadGeno}}.
#'   \item \code{MissingRateCutoff}: A numeric value *(default=0.15)*. Markers with missing
#'     rate exceeding this value will be excluded from analysis.
#'   \item \code{MinMAFCutoff}: A numeric value *(default=0.001)*. Markers with minor allele
#'     frequency (MAF) below this value will be excluded from analysis.
#'   \item \code{MinMACCutoff}: A numeric value *(default=20)*. Markers with minor allele
#'     count (MAC) below this value will be excluded from analysis.
#'   \item \code{nMarkersEachChunk}: Number of markers *(default=10000)* processed in each
#'     output chunk.
#'   }
#'  
#' The following parameters customize the columns in the \code{OutputFile}.
#' The columns \code{Marker}, \code{Info}, \code{AltFreq}, \code{AltCounts}, 
#' \code{MissingRate}, and \code{Pvalue} are included for all methods.
#'  \itemize{
#'  \item \code{outputColumns}: Specifies additional columns to include in the output.
#'     For example, for the POLMM method, users can set 
#'     \code{control$outputColumns = c("beta", "seBeta", "AltFreqInGroup")}:
#'     \itemize{
#'     \item \code{POLMM}: Default: \code{beta}, \code{seBeta}; 
#'       Optional: \code{zScore}, \code{AltFreqInGroup}, \code{nSamplesInGroup}, 
#'       \code{AltCountsInGroup}
#'     \item \code{SPACox}: Optional: \code{zScore}
#'     }
#'  }
#' @return The analysis results are written to \code{OutputFile}, which includes the following
#'   columns:
#' \describe{
#' \item{Marker}{Marker IDs extracted from \code{GenoFile} and \code{GenoFileIndex}.}
#' \item{Info}{Marker information in format "CHR:POS:REF:ALT". The order of REF/ALT depends
#'   on \code{control$AlleleOrder}: "ref-first" or "alt-first".}
#' \item{AltFreq}{Alternative allele frequency (before genotype imputation, might be > 0.5).
#'   If most markers have \code{AltFreq} > 0.5, consider resetting \code{control$AlleleOrder}.}
#' \item{AltCounts}{Alternative allele counts (before genotype imputation).}
#' \item{MissingRate}{Missing rate for each marker.}
#' \item{Pvalue}{Association test p-value.}
#' }
#' 
#' The following columns can be customized using \code{control$outputColumns}. 
#' See \code{\link{makeGroup}} for details about phenotype grouping, which is used for
#' \code{nSamplesInGroup}, \code{AltCountsInGroup}, and \code{AltFreqInGroup}.
#' \describe{
#' \item{beta}{Estimated effect size of the ALT allele.}
#' \item{seBeta}{Estimated standard error of the effect size.}
#' \item{zScore}{Standardized score statistic, usually follows a standard normal distribution.}
#' \item{nSamplesInGroup}{Number of subjects in different phenotype groups. This may differ
#'   slightly from the original distribution due to missing genotypes.}
#' \item{AltCountsInGroup}{Alternative allele counts (before genotype imputation) in different
#'   phenotype groups.}
#' \item{AltFreqInGroup}{Alternative allele frequency (before genotype imputation) in different
#'   phenotype groups.}
#' }
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
#'   obj.POLMM,
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
  NullModelClass <- checkObjNull(objNull)

  # Set default output index file if not provided
  if (is.null(OutputFileIndex)) {
    OutputFileIndex <- paste0(OutputFile, ".index")
  }

  # Validate and set control parameters with default values if not specified (control.R)
  checkControl.ReadGeno(control)
  control <- checkControl.Marker(control, NullModelClass)
  nMarkersEachChunk <- control$nMarkersEachChunk

  # Check output file status and determine restart point if needed (Util.R)
  outList <- checkOutputFile(
    OutputFile, OutputFileIndex, "Marker",
    format(nMarkersEachChunk, scientific = FALSE)
  )

  # Special validation for SPAGRM and SAGELD methods: check MAF interval constraints
  # Added by XH-2023-05-09
  if (NullModelClass %in% c("SPAGRM_NULL_Model", "SAGELD_NULL_Model")) {
    if (length(objNull$MAF_interval) > 1) {
      if (control$min_maf_marker <= min(objNull$MAF_interval)) {
        stop("min_maf_marker is out of MAF_interval. Please reset min_maf_marker or check MAF_interval.")
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
      "The analysis has been completed in an earlier run. Results have been saved in '", 
      OutputFile, 
      "'. If you want to restart the analysis, please use a different 'OutputFile'."
    )
  }

  # Check if analysis was partially completed and needs to restart
  if (!Start) {
    .message(
      "Part of the analysis has been completed and saved in %s. Restarting the analysis from chunk %d",
      OutputFileIndex, indexChunk + 1
    )
  }

  # Extract subject IDs from null model object
  subjData <- as.character(objNull$subjData)

  # Categorizes subjects into groups and enables checking of AltFreq/AltCounts within each group
  Group <- makeGroup(objNull$yVec) # Function defined in Util.R
  ifOutGroup <- any(c("AltFreqInGroup", "AltCountsInGroup") %in% control$outputColumns)

  # Set up genotype reading object with file information and subject filtering
  objGeno <- setGenoInput(GenoFile, GenoFileIndex, subjData, control) # Function in 'Geno.R'
  genoType <- objGeno$genoType
  markerInfo <- objGeno$markerInfo
  CHROM <- markerInfo$CHROM
  genoIndex <- markerInfo$genoIndex

  # Split all markers into chunks for processing. Strategy: 
  # 1. SNPs in the same chromosome are grouped into chunks
  # 2. Chunks are ordered by chromosome
  # [INLINED FUNCTION: splitMarker()]
  genoIndexList <- list()
  iTot <- 1

  # Process each chromosome separately to maintain data locality
  uCHROM <- unique(CHROM)
  for (chrom in uCHROM) {
    # Extract markers belonging to current chromosome
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
  # [END INLINED FUNCTION: splitMarker()]

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
    # This optimization avoids redundant setup for markers in the same chromosome
    if (tempChrom != chrom) {
      # [INLINED FUNCTION: setMarker()]
      # Set global variables in C++ for efficient processing
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
      # Each method has its own setup requirements and optimizations
      if (NullModelClass == "POLMM_NULL_Model") {
        obj.setMarker <- setMarker.POLMM(objNull, control, chrom)
      } else if (NullModelClass == "SPACox_NULL_Model") {
        obj.setMarker <- setMarker.SPACox(objNull, control)
      } else if (NullModelClass == "SPAmix_NULL_Model") {
        obj.setMarker <- setMarker.SPAmix(objNull, control)
      } else if (NullModelClass == "SPAGRM_NULL_Model") {
        obj.setMarker <- setMarker.SPAGRM(objNull, control)
      } else if (NullModelClass == "SAGELD_NULL_Model") {
        obj.setMarker <- setMarker.SAGELD(objNull, control)
      } else if (NullModelClass == "WtCoxG_NULL_Model") {
        obj.setMarker <- setMarker.WtCoxG(objNull, control)
      }
      # [END INLINED FUNCTION: setMarker()]

      chrom <- tempChrom
    }

    .message("---- Analyzing Chunk %d/%d: chrom %s ----", 
             i, nChunks, chrom)

    # Calculate summary statistics for all markers in current chunk
    # [INLINED FUNCTION: mainMarker()]
    # POLMM: Process ordinal traits using proportional odds logistic mixed model
    if (NullModelClass == "POLMM_NULL_Model") {
      # Call C++ implementation for efficient computation
      # See 'Main.cpp' for the underlying algorithm
      OutList <- mainMarkerInCPP("POLMM", genoType, genoIndex)

      # Construct base output data frame with required columns
      resMarker <- data.frame(
        Marker = OutList$markerVec,        # Marker IDs
        Info = OutList$infoVec,            # Marker information: CHR:POS:REF:ALT
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
    # [END INLINED FUNCTION: mainMarker()]

    # Write results to output file and update progress
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
