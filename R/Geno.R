## ------------------------------------------------------------------------------
## Geno.R
## Genotype I/O utilities: reading PLINK/BGEN data, extracting marker summaries,
## validating control settings, and preparing subject-specific genotype subsets
## for downstream null model / marker / region analyses.
##
## Functions:
##   GRAB.ReadGeno        : Read genotype file(s) with optional filtering & imputation.
##   GRAB.getGenoInfo     : Retrieve marker-level allele frequency & missingness info.
##   checkControl.ReadGeno: Validate and populate control list for reading genotypes.
##   setGenoInput         : Orchestrate genotype reading; returns marker metadata + type.
##   updateSampleIDs      : Harmonize requested SampleIDs with those present in data.
##   getSampleIDsFromBGEN : Extract sample IDs from a BGEN file header.
##   getVersionFromBGEN   : Parse BGEN file version.
##   convert4BitsToNumber : Helper for bit-level genotype encoding.
##   checkIfSampleIDsExist: Sanity check for sample presence in BGEN header.
## ------------------------------------------------------------------------------

#' Read genotype data from multiple file formats
#'
#' Reads genotype data from PLINK or BGEN format files with flexible filtering
#' and processing options. Supports efficient memory usage and various
#' imputation methods for missing genotypes.
#'
#' @param GenoFile Path to genotype file. Supported formats determined by extension:
#'   \itemize{
#'     \item PLINK: "prefix.bed" (binary format)
#'     \item BGEN: "prefix.bgen" (version 1.2 with 8-bit compression)
#'   }
#' @param GenoFileIndex Associated index files for the genotype file:
#'   \itemize{
#'     \item PLINK: c("prefix.bim", "prefix.fam") (auto-detected if NULL)
#'     \item BGEN: "prefix.bgen.bgi" or c("prefix.bgen.bgi", "prefix.sample")
#'   }
#' @param SampleIDs Character vector of sample IDs to extract. If NULL,
#'   extracts all samples.
#' @param control List of control parameters with the following options:
#' \itemize{
#'   \item \code{imputeMethod}: Imputation method for genotype data.
#'     Options: "none" (default), "mean" (2 times allele frequency).
#'     "bestguess" (round mean to the nearest integer, 0, 1, or 2).
#'   \item \code{AlleleOrder}: Allele order in genotype file. Options: "ref-first",
#'     "alt-first", or NULL (default: "alt-first" for BGEN, "ref-first" for PLINK).
#'   \item \strong{Marker Selection:}
#'   \itemize{
#'     \item \code{AllMarkers}: Set to TRUE (default) to analyze all markers.
#'       Automatically set to FALSE if any include/exclude files are provided.
#'     \item \code{IDsToIncludeFile}: Path to file with marker IDs to include.
#'     \item \code{RangesToIncludeFile}: Path to file with genomic ranges to include.
#'       Can be used with IDsToIncludeFile (union will be used).
#'     \item \code{IDsToExcludeFile}: Path to file with marker IDs to exclude.
#'     \item \code{RangesToExcludeFile}: Path to file with genomic ranges to exclude.
#'       Can be used with IDsToExcludeFile (union will be excluded).
#'     \item Note: Cannot use both include and exclude files simultaneously.
#'   }
#' }
#' @param sparse Logical indicating whether to return sparse genotype matrix
#'   (default: FALSE).
#' @return List containing:
#' \describe{
#'   \item{GenoMat}{Genotype matrix (samples × markers) with values 0, 1, 2, or NA.}
#'   \item{markerInfo}{Data frame with columns CHROM, POS, ID, REF, ALT.}
#' }
#' @details
#' **File Format Support:**
#'
#' *PLINK Format:* Binary BED/BIM/FAM files. See
#' \url{https://www.cog-genomics.org/plink/2.0/} for specifications.
#'
#' *BGEN Format:* Version 1.2 with 8-bit compression. See
#' \url{https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html} for details.
#' Requires BGI index file created with bgenix tool.
#'
#' @examples
#' ## Raw genotype data
#' RawFile <- system.file("extdata", "simuRAW.raw.gz", package = "GRAB")
#' GenoMat <- data.table::fread(RawFile)
#' GenoMat[1:10, 1:10]
#'
#' ## PLINK files
#' PLINKFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' # If include/exclude files are not specified, then control$AllMarker should be TRUE
#' GenoList <- GRAB.ReadGeno(PLINKFile, control = list(AllMarkers = TRUE))
#' GenoMat <- GenoList$GenoMat
#' markerInfo <- GenoList$markerInfo
#' head(GenoMat[, 1:6])
#' head(markerInfo)
#'
#' ## BGEN files (Note the different REF/ALT order for BGEN and PLINK formats)
#' BGENFile <- system.file("extdata", "simuBGEN.bgen", package = "GRAB")
#' GenoList <- GRAB.ReadGeno(BGENFile, control = list(AllMarkers = TRUE))
#' GenoMat <- GenoList$GenoMat
#' markerInfo <- GenoList$markerInfo
#' head(GenoMat[, 1:6])
#' head(markerInfo)
#'
#' ## The below is to demonstrate parameters in control
#' PLINKFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' IDsToIncludeFile <- system.file("extdata", "simuGENO.IDsToInclude", package = "GRAB")
#' RangesToIncludeFile <- system.file("extdata", "RangesToInclude.txt", package = "GRAB")
#' GenoList <- GRAB.ReadGeno(PLINKFile,
#'   control = list(
#'     IDsToIncludeFile = IDsToIncludeFile,
#'     RangesToIncludeFile = RangesToIncludeFile,
#'     AlleleOrder = "ref-first"
#'   )
#' )
#' GenoMat <- GenoList$GenoMat
#' head(GenoMat)
#' markerInfo <- GenoList$markerInfo
#' head(markerInfo)
#'
#' ## The below is for PLINK/BGEN files with missing data
#' PLINKFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' GenoList <- GRAB.ReadGeno(PLINKFile, control = list(AllMarkers = TRUE))
#' head(GenoList$GenoMat)
#'
#' GenoList <- GRAB.ReadGeno(PLINKFile, control = list(AllMarkers = TRUE, imputeMethod = "mean"))
#' head(GenoList$GenoMat)
#'
#' BGENFile <- system.file("extdata", "simuBGEN.bgen", package = "GRAB")
#' GenoList <- GRAB.ReadGeno(BGENFile, control = list(AllMarkers = TRUE))
#' head(GenoList$GenoMat)
#'
GRAB.ReadGeno <- function(
  GenoFile,
  GenoFileIndex = NULL,
  SampleIDs = NULL,
  control = NULL,
  sparse = FALSE
) {
  # Validate control is a list if not NULL
  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be an R list.")
  }

  # Validate and set default control parameters
  control <- checkControl.ReadGeno(control)

  # Initialize genotype input object with file paths and parameters
  objGeno <- setGenoInput(GenoFile, GenoFileIndex, SampleIDs, control)
  genoType <- objGeno$genoType # "PLINK" or "BGEN"
  markerInfo <- objGeno$markerInfo
  SampleIDs <- objGeno$SampleIDs
  anyQueue <- objGeno$anyQueue # if FALSE, no include/exclude is specified

  # Extract marker and sample information
  MarkerIDs <- markerInfo$ID
  n <- length(SampleIDs)
  m <- length(MarkerIDs)

  .message("Reading genotypes: %d subjects, %d markers", n, m)

  # Read genotype data using appropriate C++ backend
  if (sparse == TRUE) {
    # Use sparse matrix representation for memory efficiency
    GenoMat <- getSpGenoInCPP(
      t_genoType = genoType,              # character: "PLINK" or "BGEN"
      t_markerInfo = markerInfo,          # data.frame: Marker info with genoIndex column
      n = n,                              # integer: Number of subjects
      t_imputeMethod = control$imputeMethod # character: Imputation method ("none", "mean", "bestguess")
    )
  } else {
    # Use standard dense matrix representation
    GenoMat <- getGenoInCPP(
      t_genoType = genoType,              # character: "PLINK" or "BGEN"
      t_markerInfo = markerInfo,          # data.frame: Marker info with genoIndex column
      n = n,                              # integer: Number of subjects
      t_imputeMethod = control$imputeMethod # character: Imputation method ("none", "mean", "bestguess")
    )
  }

  # Set matrix row and column names for identification
  colnames(GenoMat) <- MarkerIDs
  rownames(GenoMat) <- SampleIDs

  # Keep only essential marker information columns
  markerInfo <- markerInfo[, 1:5]

  .message("Genotype reading completed")

  closeGenoInputInCPP(
    t_genoType = genoType  # character: "PLINK" or "BGEN" - file type to close
  )

  return(list(
    GenoMat = GenoMat,
    markerInfo = markerInfo
  ))
}


#' Get allele frequency and missing rate information from genotype data
#'
#' This function shares input as in function \code{GRAB.ReadGeno}, please
#' check \code{?GRAB.ReadGeno} for more details.
#'
#' @param GenoFile a character of genotype file. See \code{Details} section
#'   for more details.
#' @param GenoFileIndex additional index file(s) corresponding to
#'   \code{GenoFile}. See \code{Details} section for more details.
#' @param SampleIDs a character vector of sample IDs to extract. The default
#'   is \code{NULL}, that is, all samples in \code{GenoFile} will be extracted.
#' @param control List of control parameters with the following options:
#' \itemize{
#'   \item \code{AlleleOrder}: Allele order in genotype file. Options: "ref-first",
#'     "alt-first", or NULL (default: "alt-first" for BGEN, "ref-first" for PLINK).
#'   \item \strong{Marker Selection:}
#'   \itemize{
#'     \item \code{AllMarkers}: Set to TRUE (default) to analyze all markers.
#'       Automatically set to FALSE if any include/exclude files are provided.
#'     \item \code{IDsToIncludeFile}: Path to file with marker IDs to include.
#'     \item \code{RangesToIncludeFile}: Path to file with genomic ranges to include.
#'       Can be used with IDsToIncludeFile (union will be used).
#'     \item \code{IDsToExcludeFile}: Path to file with marker IDs to exclude.
#'     \item \code{RangesToExcludeFile}: Path to file with genomic ranges to exclude.
#'       Can be used with IDsToExcludeFile (union will be excluded).
#'     \item Note: Cannot use both include and exclude files simultaneously.
#'   }
#' }
#' @return A data frame containing marker information with allele frequencies
#'   and missing rates. The data frame includes columns from marker information
#'   (CHROM, POS, ID, REF, ALT, etc.) plus additional columns:
#' \describe{
#'   \item{altFreq}{Alternative allele frequency (before genotype imputation)}
#'   \item{missingRate}{Missing rate for each marker}
#' }
#'
GRAB.getGenoInfo <- function(
  GenoFile,
  GenoFileIndex = NULL,
  SampleIDs = NULL,
  control = NULL
) {
  # Validate control is a list if not NULL
  if (!is.null(control) && !is.list(control)) {
    stop("Argument 'control' should be an R list.")
  }

  # Validate and set default control parameters
  control <- checkControl.ReadGeno(control)

  objGeno <- setGenoInput(GenoFile, GenoFileIndex, SampleIDs, control)

  genoType <- objGeno$genoType # "PLINK" or "BGEN"
  markerInfo <- objGeno$markerInfo
  SampleIDs <- objGeno$SampleIDs
  anyQueue <- objGeno$anyQueue # if FALSE, no include/exclude is specified

  if (!anyQueue) {
    if (!control$AllMarkers) {
      stop("If include/exclude files are not specified, control$AllMarkers should be TRUE.")
    }
  }

  MarkerIDs <- markerInfo$ID

  n <- length(SampleIDs)
  m <- length(MarkerIDs)

  .message("Getting genotype info: %d subjects, %d markers", n, m)

  GenoInfoMat <- getGenoInfoInCPP(
    t_genoType = genoType,              # character: "PLINK" or "BGEN"
    t_markerInfo = markerInfo           # data.frame: Marker info with genoIndex column
  )
  GenoInfoMat <- as.data.frame(GenoInfoMat)
  colnames(GenoInfoMat) <- c("altFreq", "missingRate")

  GenoInfoMat <- cbind(markerInfo, GenoInfoMat)
  return(GenoInfoMat)
}


# Validate and set default control parameters for SPACox marker analysis
checkControl.ReadGeno <- function(control) {

  # Merge user-provided control with defaults
  default.control <- list(
    imputeMethod = "none",
    AlleleOrder = NULL,
    AllMarkers = TRUE,
    IDsToIncludeFile = NULL,
    IDsToExcludeFile = NULL,
    RangesToIncludeFile = NULL,
    RangesToExcludeFile = NULL
  )

  if (is.null(control)) {
    control <- default.control
  } else {
    control <- updateControl(control, default.control)
  }

  # Validate parameters
  if (!is.null(control$AlleleOrder)) {
    if (control$AlleleOrder != "ref-first" && control$AlleleOrder != "alt-first") {
      stop("control$AlleleOrder should be 'ref-first' or 'alt-first'.")
    }
  }

  if (!control$imputeMethod %in% c("none", "bestguess", "mean")) {
    stop("control$imputeMethod should be 'none', 'bestguess', or 'mean'.")
  }

  # Check marker selection parameters
  include_files <- c("IDsToIncludeFile", "RangesToIncludeFile")
  exclude_files <- c("IDsToExcludeFile", "RangesToExcludeFile")

  # Check which files are provided
  include_provided <- sapply(include_files, function(x) !is.null(control[[x]]))
  exclude_provided <- sapply(exclude_files, function(x) !is.null(control[[x]]))

  has_include <- any(include_provided)
  has_exclude <- any(exclude_provided)

  # Validate: cannot provide both include and exclude files
  if (has_include && has_exclude) {
    stop(
      "Cannot provide both include and exclude files. ",
      "Either use include files (IDsToIncludeFile, RangesToIncludeFile) ",
      "or exclude files (IDsToExcludeFile, RangesToExcludeFile), but not both."
    )
  }

  # Check if the files specified exist
  FileType <- c("IDsToIncludeFile", "IDsToExcludeFile", "RangesToIncludeFile", "RangesToExcludeFile")
  for (ft in FileType) {
    if (ft %in% names(control)) {
      file <- control[[ft]]
      if (!is.null(file) && !file.exists(file)) {
        stop(paste0("Cannot find the file of ", file, "..."))
      }
    }
  }

  # Log which files are being used
  if (has_include || has_exclude) {
    control$AllMarkers <- FALSE

    if (has_include) {
      used_files <- include_files[include_provided]
      .message("Marker selection: Including markers from union of:")
      for (file_param in used_files) {
        .message("  - %s: %s", file_param, control[[file_param]])
      }
    } else {
      used_files <- exclude_files[exclude_provided]
      .message("Marker selection: Excluding markers from union of:")
      for (file_param in used_files) {
        .message("  - %s: %s", file_param, control[[file_param]])
      }
    }
  } else {
    # Default: analyze all markers
    if (is.null(control$AllMarkers)) {
      control$AllMarkers <- TRUE
    }
    if (isTRUE(control$AllMarkers)) {
      .message("Marker selection: Analyzing all markers in the genotype file.")
    }
  }

  return(control)
}


# Setup an object in C++ (Main.cpp)
# PLINK format: ptr_gPLINKobj;
# BGEN format: ptr_gBGENobj;
setGenoInput <- function(
  GenoFile,
  GenoFileIndex = NULL,
  SampleIDs,
  control
) {

  if (!file.exists(GenoFile)) {
    stop("Cannot find GenoFile: ", GenoFile, ".")
  }

  GenoFileExt <- tools::file_ext(GenoFile)
  AlleleOrder <- control$AlleleOrder

  if (GenoFileExt == "bed") {

    ########## ----------  PLINK format ---------- ##########
    genoType <- "PLINK"

    if (is.null(AlleleOrder)) AlleleOrder <- "alt-first"

    if (is.null(GenoFileIndex)) {
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bim' and 'fam' files
      GenoFileIndex <- c(
        gsub("bed$", "bim", GenoFile),
        gsub("bed$", "fam", GenoFile)
      )
    }

    bimFile <- GenoFileIndex[1]
    famFile <- GenoFileIndex[2]
    bedFile <- GenoFile

    if (!file.exists(bimFile) || !file.exists(famFile) || !file.exists(bedFile)) {
      stop("One or more genotype files are missing: ",
           paste(c(bimFile, famFile, bedFile)[!file.exists(c(bimFile, famFile, bedFile))], collapse = ", "))
    }

    # Read BIM file
    .message("Reading bim file: %s", basename(bimFile))
    markerInfo <- data.table::fread(bimFile, header = FALSE)
    markerInfo <- as.data.frame(markerInfo)

    if (ncol(markerInfo) != 6) {
      stop("bim file should include 6 columns seperated by whitespace.")
    }

    # https://www.cog-genomics.org/plink/2.0/formats#bim
    if (AlleleOrder == "alt-first") {
      markerInfo <- markerInfo[, c(1, 4, 2, 6, 5)]
    } 
    if (AlleleOrder == "ref-first") {
      markerInfo <- markerInfo[, c(1, 4, 2, 5, 6)]
    }

    colnames(markerInfo) <- c("CHROM", "POS", "ID", "REF", "ALT")
    markerInfo$genoIndex <- seq_len(nrow(markerInfo)) - 1 # -1 is to convert 'R' to 'C++'

    # Read FAM file
    .message("Reading fam file: %s", basename(famFile))
    sampleInfo <- data.table::fread(famFile, header = FALSE)

    if (ncol(sampleInfo) != 6) {
      stop("fam file should include 6 columns seperated by space or '\t'.")
    }

    samplesInGeno <- sampleInfo$V2
    SampleIDs <- updateSampleIDs(SampleIDs, samplesInGeno)

    .message("Setting up PLINK object in C++ ...")
    setPLINKobjInCPP(
      t_bimFile = bimFile,        # character: Path to .bim file (marker info)
      t_famFile = famFile,        # character: Path to .fam file (sample info)
      t_bedFile = bedFile,        # character: Path to .bed file (genotype data)
      t_SampleInModel = SampleIDs, # character vector: Sample IDs to include
      t_AlleleOrder = AlleleOrder  # character: "alt-first" or "ref-first"
    )
  } else if (GenoFileExt == "bgen") {

    ########## ----------  BGEN format ---------- ##########
    genoType <- "BGEN"
    bgenFile <- GenoFile

    if (getVersionFromBGEN(bgenFile) != "v1.2") {
      stop("Package GRAB currently only supports version 1.2 of BGEN format.")
    }

    if (is.null(AlleleOrder)) AlleleOrder <- "ref-first"

    if (is.null(GenoFileIndex)) {
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bgen.bgi' file
      GenoFileIndex <- c(
        gsub("bgen$", "bgen.bgi", GenoFile),
        gsub("bgen$", "sample", GenoFile)
      )
    }

    bgiFile <- GenoFileIndex[1]
    if (!file.exists(bgiFile)) {
      stop("Cannot find bgen.bgi file of ", GenoFileIndex[1])
    }

    if (!file.exists(GenoFileIndex[2])) {
      # No sample file provided, read sample IDs from BGEN header
      samplesInGeno <- getSampleIDsFromBGEN(bgenFile)
    } else {
      # Sample file provided, read sample IDs from sample file
      sampleFile <- GenoFileIndex[2]
      .message("Reading sample file: %s", basename(sampleFile))
      sampleData <- data.table::fread(sampleFile, header = TRUE)
      if (ncol(sampleData) < 4) stop("Column number of sample file should be >= 4.")

      expected_colnames <- c("ID_1", "ID_2", "missing", "sex")
      expected_first_row <- c(0, 0, 0, "D")
      if (any(colnames(sampleData)[1:4] != expected_colnames) ||
            any(sampleData[1, 1:4] != expected_first_row)) {
        stop("Column names of sample file should be c('ID_1', 'ID_2', 'missing', 'sex') and ",
              "the first row of sample file should be c(0, 0, 0, 'D')")
      }
      
      samplesInGeno <- as.character(sampleData$ID_2[-1])
    }

    .message("Reading bgi file: %s", basename(bgiFile))
    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgiFile)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
    bgiData <- dplyr::tbl(db_con, "Variant")
    bgiData <- as.data.frame(bgiData)

    # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    if (AlleleOrder == "alt-first") {
      markerInfo <- bgiData[, c(1, 2, 3, 6, 5, 7)]
    } 
    if (AlleleOrder == "ref-first") {
      markerInfo <- bgiData[, c(1, 2, 3, 5, 6, 7)]
    }
    colnames(markerInfo) <- c("CHROM", "POS", "ID", "REF", "ALT", "genoIndex")

    SampleIDs <- updateSampleIDs(SampleIDs, samplesInGeno)

    .message("Setting up BGEN object in C++ ...")
    setBGENobjInCPP(
      t_bgenFileName = bgenFile,              # character: Path to .bgen file
      t_bgenFileIndex = bgiFile,              # character: Path to .bgi index file
      t_SampleInBgen = samplesInGeno,         # character vector: Sample IDs in BGEN
      t_SampleInModel = SampleIDs,            # character vector: Sample IDs to include
      t_isSparseDosageInBgen = FALSE,         # logical: Use sparse dosage encoding
      t_isDropmissingdosagesInBgen = FALSE,   # logical: Drop missing dosages
      t_AlleleOrder = AlleleOrder             # character: "alt-first" or "ref-first"
    )
  } else {
    stop("The current version only supports genotype input of PLINK (filename extension is ",
         "'.bed') and BGEN (filename extension is '.bgen').")
  }

  anyInclude <- FALSE
  anyExclude <- FALSE

  markersInclude <- c()
  markersExclude <- c()

  if (!is.null(control$IDsToIncludeFile)) {
    IDsToInclude <- data.table::fread(control$IDsToIncludeFile, header = FALSE, colClasses = c("character"))
    if (ncol(IDsToInclude) != 1) {
      stop("'IDsToIncludeFile' of ", control$IDsToIncludeFile, " must include exactly one column.")
    }
    IDsToInclude <- IDsToInclude[, 1]

    posRows <- which(markerInfo$ID %in% IDsToInclude)
    if (length(posRows) != 0) {
      markersInclude <- c(markersInclude, markerInfo$ID[posRows])
    }
    anyInclude <- TRUE
  }

  if (!is.null(control$RangesToIncludeFile)) {
    col_classes <- c("character", "numeric", "numeric")
    RangesToInclude <- data.table::fread(control$RangesToIncludeFile, header = FALSE,
                                         colClasses = col_classes)
    if (ncol(RangesToInclude) != 3) {
      stop("RangesToIncludeFile must include exactly three columns.")
    }

    colnames(RangesToInclude) <- c("CHROM", "START", "END")

    for (i in seq_len(nrow(RangesToInclude))) {
      posRows <- which(
        markerInfo$CHROM == RangesToInclude$CHROM[i] &
          markerInfo$POS >= RangesToInclude$START[i] &
          markerInfo$POS <= RangesToInclude$END[i]
      )
      if (length(posRows) != 0) {
        markersInclude <- c(markersInclude, markerInfo$ID[posRows])
      }
    }
    anyInclude <- TRUE
  }

  if (!is.null(control$IDsToExcludeFile)) {
    if (anyInclude) {
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
    }
    IDsToExclude <- data.table::fread(control$IDsToExcludeFile, header = FALSE, colClasses = c("character"))
    if (ncol(IDsToExclude) != 1) {
      stop("IDsToExcludeFile should only include one column.")
    }
    IDsToExclude <- IDsToExclude[, 1]

    posRows <- which(markerInfo$ID %in% IDsToExclude)
    if (length(posRows) != 0) {
      markersExclude <- c(markersExclude, markerInfo$ID[posRows])
    }
    anyExclude <- TRUE
  }

  if (!is.null(control$RangesToExcludeFile)) {
    if (anyInclude) {
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
    }
    col_classes <- c("character", "numeric", "numeric")
    RangesToExclude <- data.table::fread(control$RangesToExcludeFile, header = FALSE,
                                         colClasses = col_classes)
    if (ncol(RangesToExclude) != 3) {
      stop("RangesToExcludeFile should only include three columns.")
    }

    colnames(RangesToExclude) <- c("CHROM", "START", "END")

    for (i in seq_len(nrow(RangesToExclude))) {
      posRows <- which(
        markerInfo$CHROM == RangesToExclude$CHROM[i] &
          markerInfo$POS >= RangesToExclude$START[i] &
          markerInfo$POS <= RangesToExclude$END[i]
      )
      if (length(posRows) != 0) {
        markersExclude <- c(markersExclude, markerInfo$ID[posRows])
      }
    }
    anyExclude <- TRUE
  }

  markersInclude <- unique(markersInclude)
  markersExclude <- unique(markersExclude)

  # return genotype
  if (anyInclude) {
    markerInfo <- subset(markerInfo, ID %in% markersInclude)
  }

  if (anyExclude) {
    markerInfo <- subset(markerInfo, !ID %in% markersExclude)
  }

  anyQueue <- anyInclude | anyExclude

  # convert integer64 to numeric, which is supported in c++
  markerInfo$genoIndex <- as.numeric(markerInfo$genoIndex)

  genoList <- list(
    genoType = genoType,            # character: "PLINK" or "BGEN"
    markerInfo = markerInfo,        # data.frame: CHROM, POS, ID, REF, ALT, genoIndex
    SampleIDs = SampleIDs,          # character vector: IDs also in genotype file
    AlleleOrder = AlleleOrder,      # character: "ref-first" or "alt-first"
    GenoFile = GenoFile,            # character: Genotype file path
    GenoFileIndex = GenoFileIndex,  # character vector: Index file path(s)
    anyQueue = anyQueue             # logical: if any include/exclude is specified
  )

  return(genoList)
}


# Update SampleIDs based on samples in genotype file
updateSampleIDs <- function(SampleIDs, samplesInGeno) {
  if (is.null(SampleIDs)) {
    .message("Using all samples from genotype file (%d samples)", length(samplesInGeno))
    SampleIDs <- samplesInGeno
  }

  if (any(!SampleIDs %in% samplesInGeno)) {
    missing_samples <- SampleIDs[!SampleIDs %in% samplesInGeno]
    .message("Warning: %d samples not found in genotype file", length(missing_samples))
    stop("The above samples from 'SampleIDs' are not in 'GenoFile' and 'GenoFileIndex'.")
  }

  SampleIDs <- as.character(SampleIDs)

  return(SampleIDs)
}


# Extract sample identifiers from BGEN file (only support BGEN v1.2)
# Check https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
getSampleIDsFromBGEN <- function(bgenFile) {
  if (!checkIfSampleIDsExist(bgenFile)) {
    stop("The BGEN file does not include subject IDs. Check ",
         "https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html for details.")
  }

  .message("Extracting sample information from BGEN file")
  con <- file(bgenFile, "rb")
  seek(con, 4)
  LH <- readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH + 4)
  N <- readBin(con, n = 1, what = "integer", size = 4) # number of samples
  samplesInGeno <- rep(0, N)

  # cycle for all samples to extract IDs
  for (i in 1:N) {
    LS <- readBin(con, n = 1, what = "integer", size = 2)
    sample <- readChar(con, nchars = LS, useBytes = TRUE)
    samplesInGeno[i] <- sample
  }

  # close connection
  close(con)

  return(samplesInGeno)
}


# Get version information from BGEN file
getVersionFromBGEN <- function(bgenFile) {
  con <- file(bgenFile, "rb")
  seek(con, 4)
  LH <- readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH - 4)
  header <- rawToBits(readBin(con, n = 4, what = "raw", size = 1, signed = FALSE))

  VersionNum <- convert4BitsToNumber(header[3:6])
  if (VersionNum == 0) {
    version <- paste0("Version Layout = 0, which is not supported. Please check ",
                      "https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html for more details.")
  }

  if (VersionNum == 1) {
    version <- "v1.1"
  }

  if (VersionNum == 2) {
    version <- "v1.2"
  }

  if (VersionNum > 2) {
    version <- paste0("Version Layout > 2, which is reserved for future use. Please check ",
                      "https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html for more details.")
  }

  close(con)
  return(version)
}


# Helper function to convert 4 bits (least-significant first) to a number
convert4BitsToNumber <- function(leastSignificantBit) {
  leastSignificantBit <- as.numeric(leastSignificantBit)
  if (length(leastSignificantBit) != 4) {
    stop("Input should be 4 bits in which least-significant first.")
  }

  Number <- 0
  for (i in 1:4) {
    Number <- Number + 2^(i - 1) * leastSignificantBit[i]
  }
  return(Number)
}


# Check if sample identifiers are stored in a BGEN file
checkIfSampleIDsExist <- function(bgenFile) {
  con <- file(bgenFile, "rb")
  seek(con, 4)
  LH <- readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH - 4)
  header <- rawToBits(readBin(con, n = 4, what = "raw", size = 1, signed = FALSE))
  close(con)
  return(header[32] == 01)
}
