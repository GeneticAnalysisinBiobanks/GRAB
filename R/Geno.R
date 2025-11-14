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
#'   \item \code{ImputeMethod}: Imputation method for genotype data.
#'     Options: "none" (default), "bestguess", "mean".
#'   \item \code{AlleleOrder}: Allele order in genotype file. Options: "ref-first",
#'     "alt-first", or NULL (default, uses file-type default).
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
#' GenoList <- GRAB.ReadGeno(PLINKFile, control = list(AllMarkers = TRUE, ImputeMethod = "mean"))
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

  # Validate that AllMarkers is TRUE when no filtering is specified
  if (!anyQueue) {
    if (!control$AllMarkers) {
      stop("If include/exclude files are not specified, control$AllMarkers should be TRUE.")
    }
  }

  # Extract marker and sample information
  MarkerIDs <- markerInfo$ID
  n <- length(SampleIDs)
  m <- length(MarkerIDs)

  .message("Reading genotypes: %d subjects, %d markers", n, m)

  # Read genotype data using appropriate C++ backend
  if (sparse == TRUE) {
    # Use sparse matrix representation for memory efficiency
    GenoMat <- getSpGenoInCPP(genoType, markerInfo, n, control$ImputeMethod) # check Main.cpp
  } else {
    # Use standard dense matrix representation
    GenoMat <- getGenoInCPP(genoType, markerInfo, n, control$ImputeMethod) # check Main.cpp
  }

  # Set matrix row and column names for identification
  colnames(GenoMat) <- MarkerIDs
  rownames(GenoMat) <- SampleIDs

  # Keep only essential marker information columns
  markerInfo <- markerInfo[, 1:5]

  .message("Genotype reading completed")

  closeGenoInputInCPP(genoType) # "PLINK" or "BGEN"

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
#'   \item \code{ImputeMethod}: Imputation method for genotype data.
#'     Options: "none" (default), "bestguess", "mean".
#'   \item \code{AlleleOrder}: Allele order in genotype file. Options: "ref-first",
#'     "alt-first", or NULL (default, uses file-type default).
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

  GenoInfoMat <- getGenoInfoInCPP(genoType, markerInfo, control$ImputeMethod) # check Main.cpp
  GenoInfoMat <- as.data.frame(GenoInfoMat)
  colnames(GenoInfoMat) <- c("altFreq", "missingRate")

  GenoInfoMat <- cbind(markerInfo, GenoInfoMat)
  return(GenoInfoMat)
}


# Validate and set default control parameters for SPACox marker analysis
checkControl.ReadGeno <- function(control) {

  # Merge user-provided control with defaults
  default.control <- list(
    ImputeMethod = "none",
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
    control <- utils::modifyList(default.control, control)
  }

  # Validate parameters
  if (!is.null(control$AlleleOrder)) {
    if (control$AlleleOrder != "ref-first" && control$AlleleOrder != "alt-first") {
      stop("control$AlleleOrder should be 'ref-first' or 'alt-first'.")
    }
  }

  if (!control$ImputeMethod %in% c("none", "bestguess", "mean")) {
    stop("control$ImputeMethod should be 'none', 'bestguess', or 'mean'.")
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
  SampleIDs = NULL,
  control = NULL
) {

  if (!file.exists(GenoFile)) {
    stop("Cannot find GenoFile: ", GenoFile, ".")
  }

  GenoFileExt <- tools::file_ext(GenoFile)

  # Currently, only support PLINK and BGEN

  if (GenoFileExt != "bed" && GenoFileExt != "bgen") {
    stop("The current version only supports genotype input of PLINK (filename extension is ",
         "'.bed') and BGEN (filename extension is '.bgen').")
  }

  AlleleOrder <- control$AlleleOrder

  ########## ----------  PLINK format ---------- ##########

  if (GenoFileExt == "bed") {
    genoType <- "PLINK"

    if (is.null(AlleleOrder)) AlleleOrder <- "alt-first"

    if (is.null(GenoFileIndex)) {
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bim' and 'fam' files
      GenoFileIndex <- c(
        gsub("bed$", "bim", GenoFile),
        gsub("bed$", "fam", GenoFile)
      )
    }

    if (length(GenoFileIndex) != 2) {
      stop("If PLINK format is used, argument 'GenoFileIndex' should be 'NULL' or a character ",
           "vector of c(bimFile, famFile).")
    }

    bimFile <- GenoFileIndex[1]
    famFile <- GenoFileIndex[2]
    bedFile <- GenoFile

    # Read in BIM file

    if (!file.exists(bimFile)) stop(paste("Cannot find bim file of", bimFile))

    .message("Reading bim file: %s", basename(bimFile))
    markerInfo <- data.table::fread(bimFile, header = FALSE, sep = "\t")
    markerInfo <- as.data.frame(markerInfo)

    if (ncol(markerInfo) != 6) {
      stop("bim file should include 6 columns seperated by '\t'.")
    }

    if (AlleleOrder == "alt-first") {
      markerInfo <- markerInfo[, c(1, 4, 2, 6, 5)]
    } # https://www.cog-genomics.org/plink/2.0/formats#bim
    if (AlleleOrder == "ref-first") {
      markerInfo <- markerInfo[, c(1, 4, 2, 5, 6)]
    } # https://www.cog-genomics.org/plink/2.0/formats#bim

    colnames(markerInfo) <- c("CHROM", "POS", "ID", "REF", "ALT")
    markerInfo$genoIndex <- seq_len(nrow(markerInfo)) - 1 # -1 is to convert 'R' to 'C++'

    # Read in FAM file

    if (!file.exists(famFile)) stop(paste("Cannot find fam file of", famFile))

    .message("Reading fam file: %s", basename(famFile))
    sampleInfo <- data.table::fread(famFile, header = FALSE, sep = " ")

    if (ncol(sampleInfo) == 1) {
      sampleInfo <- data.table::fread(famFile, header = FALSE, sep = "\t")
    }

    if (ncol(sampleInfo) != 6) {
      stop("fam file should include 6 columns seperated by space or '\t'.")
    }

    samplesInGeno <- sampleInfo$V2
    SampleIDs <- updateSampleIDs(SampleIDs, samplesInGeno)

    .message("Setting up PLINK object in C++ ...")
    setPLINKobjInCPP(bimFile, famFile, bedFile, SampleIDs, AlleleOrder)
  }

  ########## ----------  BGEN format ---------- ##########

  if (GenoFileExt == "bgen") {
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

    if (length(GenoFileIndex) != 1 && length(GenoFileIndex) != 2) {
      stop("For genotype input of BGEN format, 'GenoFileIndex' should be of length 1 or 2. ",
           "Check 'Details' section in '?GRAB.ReadGeno' for more details.")
    }

    if (length(GenoFileIndex) == 1) {
      samplesInGeno <- getSampleIDsFromBGEN(bgenFile)
    }

    if (length(GenoFileIndex) == 2) {
      sampleFile <- GenoFileIndex[2]
      if (!file.exists(sampleFile)) {
        if (!checkIfSampleIDsExist(bgenFile)) {
          stop("Cannot find bgen.samples file of", sampleFile)
        } else {
          samplesInGeno <- getSampleIDsFromBGEN(bgenFile)
        }
      } else {
        .message("Reading sample file: %s", basename(sampleFile))
        sampleData <- data.table::fread(sampleFile, header = TRUE, sep = " ")
        if (ncol(sampleData) < 4) {
          stop("Column number of sample file should be >= 4.")
        }

        expected_colnames <- c("ID_1", "ID_2", "missing", "sex")
        expected_first_row <- c(0, 0, 0, "D")
        if (any(colnames(sampleData)[1:4] != expected_colnames) ||
            any(sampleData[1, 1:4] != expected_first_row)) {
          stop("Column names of sample file should be c('ID_1', 'ID_2', 'missing', 'sex') and ",
               "the first row of sample file should be c(0,0,0,'D')")
        }

        samplesInGeno <- as.character(sampleData$ID_2[-1])
      }
    }

    bgiFile <- GenoFileIndex[1]

    if (!file.exists(bgiFile)) stop(paste("Cannot find bgi file of", bgiFile))

    .message("Reading bgi file: %s", basename(bgiFile))
    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgiFile)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
    bgiData <- dplyr::tbl(db_con, "Variant")
    bgiData <- as.data.frame(bgiData)

    if (AlleleOrder == "alt-first") {
      markerInfo <- bgiData[, c(1, 2, 3, 6, 5, 7)]
    } # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    if (AlleleOrder == "ref-first") {
      markerInfo <- bgiData[, c(1, 2, 3, 5, 6, 7)]
    } # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html

    colnames(markerInfo) <- c("CHROM", "POS", "ID", "REF", "ALT", "genoIndex")

    SampleIDs <- updateSampleIDs(SampleIDs, samplesInGeno)

    .message("Setting up BGEN object in C++ ...")
    setBGENobjInCPP(bgenFile, bgiFile, samplesInGeno, SampleIDs, FALSE, FALSE, AlleleOrder)
  }

  ########## ----------  More format such as VCF will be supported later ---------- ##########

  anyInclude <- FALSE
  anyExclude <- FALSE

  markersInclude <- c()
  markersExclude <- c()

  if (!is.null(control$IDsToIncludeFile)) {
    IDsToInclude <- data.table::fread(control$IDsToIncludeFile, header = FALSE, colClasses = c("character"))
    if (ncol(IDsToInclude) != 1) {
      stop("'IDsToIncludeFile' of ", control$IDsToIncludeFile, " should only include one column.")
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
      stop("RangesToIncludeFile should only include three columns.")
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

  # added on 2022-04-07: avoid potential error due to "integer64", which is not well
  # supported between C++ and R
  markerInfo$genoIndex <- as.numeric(markerInfo$genoIndex)

  genoList <- list(
    genoType = genoType,
    markerInfo = markerInfo,
    SampleIDs = SampleIDs,
    AlleleOrder = AlleleOrder,
    GenoFile = GenoFile,
    GenoFileIndex = GenoFileIndex,
    anyQueue = anyQueue
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
