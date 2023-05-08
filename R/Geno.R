
#' Read in genotype data
#' 
#' \code{GRAB} package provides functions to read in genotype data. Currently, we support genotype formats of PLINK and BGEN. Other formats such as VCF will be added later.
#' 
#' @param GenoFile a character of genotype file. See \code{Details} section for more details.
#' @param GenoFileIndex additional index file(s) corresponding to \code{GenoFile}. See \code{Details} section for more details.
#' @param SampleIDs a character vector of sample IDs to extract. The default is \code{NULL}, that is, all samples in \code{GenoFile} will be extracted.
#' @param control a list of parameters to decide which markers to extract. See \code{Details} section for more details.
#' @param sparse a logical value *(default: FALSE)* to indicate if the output of genotype matrix is sparse.
#' @return An R list including a genotype matrix and an information matrix. 
#' \itemize{
#'  \item \code{GenoMat}: Genotype matrix, each row is for one sample and each column is for one marker. 
#'  \item \code{markerInfo}: Information matrix including 5 columns of CHROM, POS, ID, REF, and ALT.
#' }
#' @details
#' ## Details about \code{GenoFile} and \code{GenoFileIndex}
#' Currently, we support two formats of genotype input including PLINK and BGEN. Other formats such as VCF will be added later. 
#' Users do not need to specify the genotype format, \code{GRAB} package will check the extension of the file name for that purpose.  
#' If \code{GenoFileIndex} is not specified, \code{GRAB} package assumes the prefix is the same as \code{GenoFile}.
#' \describe{
#'   PLINK format (Check [link](https://www.cog-genomics.org/plink/2.0/) for more details about this format)
#'   \itemize{
#'   \item \code{GenoFile}: "prefix.bed". The full file name (including the extension ".bed") of the PLINK binary \code{bed} file. 
#'   \item \code{GenoFileIndex}: c("prefix.bim", "prefix.fam"). If not specified, \code{GRAB} package assumes that \code{bim} and \code{fam} files have the same prefix as the \code{bed} file.
#'   }
#'   BGEN format (Check [link](https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html) for more details about this format. Currently, only version 1.2 with 8 bits suppression is supported)
#'     \itemize{
#'     \item \code{GenoFile}: "prefix.bgen". The full file name (including the extension ".bgen") of the BGEN binary \code{bgen} file. 
#'     \item \code{GenoFileIndex}: "prefix.bgen.bgi" or c("prefix.bgen.bgi", "prefix.sample"). If not specified, \code{GRAB} package assumes that \code{bgi} and \code{sample} files have the same prefix as the \code{bgen} file.
#'     If only one element is given for \code{GenoFileIndex}, then it should be a \code{bgi} file.  Check [link](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md) for more details about \code{bgi} file.
#'     \item If the \code{bgen} file does not include sample identifiers, then \code{sample} file is required, whose detailed description can ben seen in [link](https://www.cog-genomics.org/plink/2.0/formats#sample). 
#'     If you are not sure if sample identifiers are in BGEN file, please refer to \code{\link{checkIfSampleIDsExist}}.
#'     }
#'   VCF format will be supported later. \code{GenoFile}: "prefix.vcf"; \code{GenoFileIndex}: "prefix.vcf.tbi"
#' }
#' 
#' ## Details about argument \code{control}
#' Argument \code{control} is used to include and exclude markers for function \code{GRAB.ReadGeno}. 
#' The function supports two include files of (\code{IDsToIncludeFile}, \code{RangesToIncludeFile}) and two exclude files of (\code{IDsToExcludeFile}, \code{RangesToExcludeFile}), 
#' but does not support both include and exclude files at the same time.
#' \describe{
#'   \itemize{
#'   \item \code{IDsToIncludeFile}: a file of marker IDs to include, one column (no header). Check \code{system.file("extdata", "IDsToInclude.txt", package = "GRAB")} for an example. 
#'   \item \code{IDsToExcludeFile}: a file of marker IDs to exclude, one column (no header). 
#'   \item \code{RangesToIncludeFile}: a file of ranges to include, three columns (no headers): chromosome, start position, end position. Check \code{system.file("extdata", "RangesToInclude.txt", package = "GRAB")} for an example.
#'   \item \code{RangesToExcludeFile}: a file of ranges to exclude, three columns (no headers): chromosome, start position, end position.
#'   \item \code{AlleleOrder}: a character, "ref-first" or "alt-first", to determine whether the REF/major allele should appear first or second. Default is "alt-first" for PLINK and "ref-first" for BGEN. If the ALT allele frequencies of most markers are > 0.5, you should consider resetting this option. NOTE, if you use plink2 to convert PLINK file to BGEN file, then 'ref-first' modifier is to reset the order.
#'   \item \code{AllMarkers}: a logical value (default: FALSE) to indicate if all markers are extracted. It might take too much memory to put genotype of all markers in R. This parameter is to remind users.  
#'   \item \code{ImputeMethod}: a character, "none" (default), "bestguess", or "mean". By default, missing genotype is \code{NA}. Suppose alternative allele frequency is \code{p}, then missing genotype is imputed as \code{2p} (ImputeMethod = "mean") or \code{round(2p)} (ImputeMethod = "bestguess").
#'   }
#' }
#' 
#' @examples
#' 
#' ## Raw genotype data 
#' RawFile = system.file("extdata", "simuRAW.raw", package = "GRAB")
#' GenoMat = data.table::fread(RawFile)
#' GenoMat[1:10,1:10]
#' 
#' ## PLINK files
#' PLINKFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' GenoList = GRAB.ReadGeno(PLINKFile, control = list(AllMarkers = TRUE)) # If include/exclude files are not specified, then control$AllMarker should be TRUE
#' GenoMat = GenoList$GenoMat
#' markerInfo = GenoList$markerInfo
#' head(GenoMat[,1:6])
#' head(markerInfo)
#' 
#' ## BGEN files (Note the different REF/ALT order for BGEN and PLINK formats)
#' BGENFile = system.file("extdata", "simuBGEN.bgen", package = "GRAB")
#' GenoList = GRAB.ReadGeno(BGENFile, control = list(AllMarkers = TRUE))
#' GenoMat = GenoList$GenoMat
#' markerInfo = GenoList$markerInfo
#' head(GenoMat[,1:6])
#' head(markerInfo)
#' 
#' ## The below is to demonstrate parameters in control
#' PLINKFile = system.file("extdata", "example.bed", package = "GRAB")
#' IDsToIncludeFile = system.file("extdata", "example.IDsToIncludeFile.txt", package = "GRAB")
#' RangesToIncludeFile = system.file("extdata", "example.RangesToIncludeFile.txt", package = "GRAB")
#' GenoList = GRAB.ReadGeno(PLINKFile, 
#'                          control = list(IDsToIncludeFile = IDsToIncludeFile, 
#'                                         RangesToIncludeFile = RangesToIncludeFile,
#'                                         AlleleOrder = "ref-first"))
#' GenoMat = GenoList$GenoMat
#' head(GenoMat)
#' markerInfo = GenoList$markerInfo
#' head(markerInfo)
#' 
#' ## The below is for PLINK/BGEN files with missing data
#' PLINKFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' GenoList = GRAB.ReadGeno(PLINKFile, control = list(AllMarkers = TRUE))
#' head(GenoList$GenoMat)
#' 
#' GenoList = GRAB.ReadGeno(PLINKFile, control = list(AllMarkers = TRUE, ImputeMethod = "mean"))
#' head(GenoList$GenoMat)
#' 
#' BGENFile = system.file("extdata", "simuBGEN.bgen", package = "GRAB")
#' GenoList = GRAB.ReadGeno(BGENFile, control = list(AllMarkers = TRUE))
#' head(GenoList$GenoMat)
#' 
#' @export
#' @import data.table, tidyr, dbplyr, RSQLite, Matrix
GRAB.ReadGeno = function(GenoFile,
                         GenoFileIndex = NULL,
                         SampleIDs = NULL,
                         control = NULL,
                         sparse = FALSE)
{
  control = checkControl.ReadGeno(control)  # check 'control.R'
  
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs, control)
  
  genoType = objGeno$genoType   # "PLINK" or "BGEN"
  markerInfo = objGeno$markerInfo
  SampleIDs = objGeno$SampleIDs
  anyQueue = objGeno$anyQueue   # if FALSE, no include/exclude is specified
  
  if(!anyQueue)
    if(!control$AllMarkers)
      stop("If include/exclude files are not specified, control$AllMarkers should be TRUE.")

  MarkerIDs = markerInfo$ID
  
  n = length(SampleIDs)
  m = length(MarkerIDs)
  
  cat("Number of Samples:\t", n, "\n")
  cat("Number of Markers:\t", m, "\n")
  
  if(sparse == TRUE){
    GenoMat = getSpGenoInCPP(genoType, markerInfo, n, control$ImputeMethod)  # check Main.cpp
  }else{
    GenoMat = getGenoInCPP(genoType, markerInfo, n, control$ImputeMethod)  # check Main.cpp
  }
    
  colnames(GenoMat) = MarkerIDs;
  rownames(GenoMat) = SampleIDs;
  
  markerInfo = markerInfo[,1:5]
  
  cat("Complete the genotype reading.\n")
  
  closeGenoInputInCPP(genoType)  # "PLINK" or "BGEN"
  
  return(list(GenoMat = GenoMat,
              markerInfo = markerInfo))
}

# (2023-05-03) Get allele frequency and missing rate information from genotype data

#' Get allele frequency and missing rate information from genotype data
#' 
#' This function shares input as in function \code{GRAB.ReadGeno}, please check \code{?GRAB.ReadGeno} for more details.
#' 
#' @param GenoFile a character of genotype file. See \code{Details} section for more details.
#' @param GenoFileIndex additional index file(s) corresponding to \code{GenoFile}. See \code{Details} section for more details.
#' @param SampleIDs a character vector of sample IDs to extract. The default is \code{NULL}, that is, all samples in \code{GenoFile} will be extracted.
#' @param control a list of parameters to decide which markers to extract. See \code{Details} section for more details.
GRAB.getGenoInfo = function(GenoFile,
                            GenoFileIndex = NULL,
                            SampleIDs = NULL,
                            control = NULL)
{
  control = checkControl.ReadGeno(control)  # check 'control.R'
  
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs, control)
  
  genoType = objGeno$genoType   # "PLINK" or "BGEN"
  markerInfo = objGeno$markerInfo
  SampleIDs = objGeno$SampleIDs
  anyQueue = objGeno$anyQueue   # if FALSE, no include/exclude is specified
  
  if(!anyQueue)
    if(!control$AllMarkers)
      stop("If include/exclude files are not specified, control$AllMarkers should be TRUE.")
  
  MarkerIDs = markerInfo$ID
  
  n = length(SampleIDs)
  m = length(MarkerIDs)
  
  cat("Number of Samples:\t", n, "\n")
  cat("Number of Markers:\t", m, "\n")
  
  GenoInfoMat = getGenoInfoInCPP(genoType, markerInfo, control$ImputeMethod)  # check Main.cpp
  GenoInfoMat = as.data.frame(GenoInfoMat)
  colnames(GenoInfoMat) = c("altFreq", "missingRate")
  
  GenoInfoMat = cbind(markerInfo, GenoInfoMat)
  return(GenoInfoMat)
}


# setGenoInput() is to setup the following object in C++ (Main.cpp)
# PLINK format: ptr_gPLINKobj;
# BGEN format: ptr_gBGENobj;
setGenoInput = function(GenoFile, 
                        GenoFileIndex = NULL, 
                        SampleIDs = NULL,
                        control = NULL)
{
  if(missing(GenoFile))
    stop("Argument 'GenoFile' is required.")
  
  if(!file.exists(GenoFile))
    stop("Cannot find genotype file of ", GenoFile, ".")
  
  GenoFileExt = tools::file_ext(GenoFile);
  
  # Currently, only support PLINK and BGEN
  
  if(GenoFileExt != "bed" & GenoFileExt != "bgen")
    stop("The current version only supports genotype input of PLINK (filename extension is '.bed') and BGEN (filename extension is '.bgen').")  
  
  AlleleOrder = control$AlleleOrder
  
  ########## ----------  PLINK format ---------- ##########
  
  if(GenoFileExt == "bed"){
    
    genoType = "PLINK"
    
    if(is.null(AlleleOrder)) AlleleOrder = "alt-first"
    
    if(is.null(GenoFileIndex)){  
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bim' and 'fam' files
      GenoFileIndex = c(gsub("bed$", "bim", GenoFile),
                        gsub("bed$", "fam", GenoFile))
    }
    
    if(length(GenoFileIndex) != 2)
      stop("If PLINK format is used, argument 'GenoFileIndex' should be 'NULL' or a character vector of c(bimFile, famFile).")
    
    bimFile = GenoFileIndex[1]
    famFile = GenoFileIndex[2]
    bedFile = GenoFile
    
    # Read in BIM file
    
    if(!file.exists(bimFile)) stop(paste("Cannot find bim file of", bimFile))
    
    cat("Reading bim file:\t", bimFile, "\n")
    markerInfo = data.table::fread(bimFile, header = F, sep = "\t")
    markerInfo = as.data.frame(markerInfo)
    
    if(ncol(markerInfo) != 6)
      stop("bim file should include 6 columns seperated by '\t'.")
    
    if(AlleleOrder == "alt-first")
      markerInfo = markerInfo[,c(1,4,2,6,5)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    if(AlleleOrder == "ref-first")
      markerInfo = markerInfo[,c(1,4,2,5,6)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    
    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT")
    markerInfo$genoIndex = 1:nrow(markerInfo) - 1  # -1 is to convert 'R' to 'C++' 
    
    # Read in FAM file
    
    if(!file.exists(famFile)) stop(paste("Cannot find fam file of", famFile))
    
    cat("Reading fam file:\t", famFile, "\n")
    sampleInfo = data.table::fread(famFile, header = F, sep = " ")
    
    if(ncol(sampleInfo) == 1)
      sampleInfo = data.table::fread(famFile, header = F, sep = "\t")
      
    if(ncol(sampleInfo) != 6)
      stop("fam file should include 6 columns seperated by space or '\t'.")
    
    samplesInGeno = sampleInfo$V2
    SampleIDs = updateSampleIDs(SampleIDs, samplesInGeno)
      
    cat("setting PLINK object in CPP....\n")
    setPLINKobjInCPP(bimFile, famFile, bedFile, SampleIDs, AlleleOrder)
  }
  
  ########## ----------  BGEN format ---------- ##########
  
  if(GenoFileExt == "bgen"){
    
    genoType = "BGEN"
    bgenFile = GenoFile
    
    if(getVersionFromBGEN(bgenFile) != "v1.2")
      stop("Package GRAB currently only supports version 1.2 of BGEN format.")
    
    if(is.null(AlleleOrder)) AlleleOrder = "ref-first"
    
    if(is.null(GenoFileIndex)){  
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bgen.bgi' file
      GenoFileIndex = c(gsub("bgen$", "bgen.bgi", GenoFile),
                        gsub("bgen$", "sample", GenoFile))
    }
    
    if(length(GenoFileIndex) != 1 & length(GenoFileIndex) != 2)
      stop("For genotype input of BGEN format, 'GenoFileIndex' should be of length 1 or 2. Check 'Details' section in '?GRAB.ReadGeno' for more details.")
    
    if(length(GenoFileIndex) == 1){
      samplesInGeno = getSampleIDsFromBGEN(bgenFile)
    }
      
    if(length(GenoFileIndex) == 2){
      sampleFile = GenoFileIndex[2]
      if(!file.exists(sampleFile)){
        if(!checkIfSampleIDsExist(bgenFile)){
          stop("Cannot find bgen.samples file of", sampleFile)
        }else{
          samplesInGeno = getSampleIDsFromBGEN(bgenFile)
        }
      }else{
        cat("Reading sample file:\t", sampleFile, "\n")
        sampleData = data.table::fread(sampleFile, header=T, sep = " ")
        if(ncol(sampleData) < 4)
          stop("Column number of sample file should be >= 4.")
        
        if(any(colnames(sampleData)[1:4] != c("ID_1", "ID_2", "missing", "sex")) | any(sampleData[1,1:4] != c(0,0,0,"D")))
          stop("Column names of sample file should be c('ID_1', 'ID_2', 'missing', 'sex') and the first row of sample file should be c(0,0,0,'D')")
        
        samplesInGeno = as.character(sampleData$ID_2[-1])
      }
    }
    
    bgiFile = GenoFileIndex[1]
    
    if(!file.exists(bgiFile)) stop(paste("Cannot find bgi file of", bgiFile))
    
    cat("Reading bgi file:\t", bgiFile, "\n")
    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgiFile)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
    bgiData = dplyr::tbl(db_con, "Variant")
    bgiData = as.data.frame(bgiData)
    
    if(AlleleOrder == "alt-first")
      markerInfo = bgiData[,c(1,2,3,6,5,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    if(AlleleOrder == "ref-first")
      markerInfo = bgiData[,c(1,2,3,5,6,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    
    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT","genoIndex")
    
    SampleIDs = updateSampleIDs(SampleIDs, samplesInGeno)
    
    setBGENobjInCPP(bgenFile, bgiFile, samplesInGeno, SampleIDs, F, F, AlleleOrder)
  }
  
  ########## ----------  More format such as VCF will be supported later ---------- ##########

  Files = c("IDsToIncludeFile", "IDsToExcludeFile", "RangesToIncludeFile", "RangesToExcludeFile")
  
  anyInclude = FALSE
  anyExclude = FALSE
  
  markersInclude = c()
  markersExclude = c()
  
  if(!is.null(control$IDsToIncludeFile)){
    IDsToInclude = data.table::fread(control$IDsToIncludeFile, header = F, colClasses = c("character"))
    if(ncol(IDsToInclude) != 1)
      stop("'IDsToIncludeFile' of ", control$IDsToIncludeFile, " should only include one column.")
    IDsToInclude = IDsToInclude[,1]
    
    posRows = which(markerInfo$ID %in% IDsToInclude)
    if(length(posRows) != 0) 
      markersInclude = c(markersInclude, markerInfo$ID[posRows])
    anyInclude = TRUE
  }
  
  if(!is.null(control$RangesToIncludeFile)){
    RangesToInclude = data.table::fread(control$RangesToIncludeFile, header = F, colClasses = c("character", "numeric", "numeric"))
    if(ncol(RangesToInclude) != 3)
      stop("RangesToIncludeFile should only include three columns.")
    
    colnames(RangesToInclude) = c("CHROM", "START", "END")
    
    for(i in 1:nrow(RangesToInclude)){
      CHROM1 = RangesToInclude$CHROM[i]
      START = RangesToInclude$START[i]
      END = RangesToInclude$END[i]
      posRows = with(markerInfo, which(CHROM == CHROM1 & POS >= START & POS <= END))
      if(length(posRows) != 0) 
        markersInclude = c(markersInclude, markerInfo$ID[posRows])
    }
    anyInclude = TRUE
  }
  
  if(!is.null(control$IDsToExcludeFile)){
    if(anyInclude) 
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
    IDsToExclude = data.table::fread(control$IDsToExcludeFile, header = F, colClasses = c("character"))
    if(ncol(IDsToExclude) != 1)
      stop("IDsToExcludeFile should only include one column.")
    IDsToExclude = IDsToExclude[,1]
    
    posRows = which(markerInfo$ID %in% IDsToExclude)
    if(length(posRows) != 0) 
      markersExclude = c(markersExclude, markerInfo$ID[posRows])
    anyExclude = TRUE
  }
  
  if(!is.null(control$RangesToExcludeFile)){
    if(anyInclude) 
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile' at the same time.")
    RangesToExclude = data.table::fread(control$RangesToExcludeFile, header = F, colClasses = c("character", "numeric", "numeric"))
    if(ncol(RangesToExclude) != 3)
      stop("RangesToExcludeFile should only include three columns.")
    
    colnames(RangesToExclude) = c("CHROM", "START", "END")
    
    for(i in 1:nrow(RangesToExclude)){
      CHROM1 = RangesToExclude$CHROM[i]
      START = RangesToExclude$START[i]
      END = RangesToExclude$END[i]
      posRows = with(markerInfo, which(CHROM == CHROM1 & POS >= START & POS <= END))
      if(length(posRows) != 0) 
        markersExclude = c(markersExclude, markerInfo$ID[posRows])
    }
    anyExclude = TRUE
  }
  
  markersInclude = unique(markersInclude)
  markersExclude = unique(markersExclude)
  
  # return genotype
  cat("Based on the 'GenoFile' and 'GenoFileIndex',", genoType, "format is used for genotype data.\n")
  
  if(anyInclude)
    markerInfo = subset(markerInfo, ID %in% markersInclude)
  
  if(anyExclude)
    markerInfo = subset(markerInfo, !ID %in% markersExclude)
  
  anyQueue = anyInclude | anyExclude
  
  markerInfo$genoIndex = as.numeric(markerInfo$genoIndex) # added on 2022-04-07: avoid potential error due to "integer64", which is not well supported between C++ and R
  
  genoList = list(genoType = genoType, 
                  markerInfo = markerInfo, 
                  SampleIDs = SampleIDs, 
                  AlleleOrder = AlleleOrder, 
                  GenoFile = GenoFile, 
                  GenoFileIndex = GenoFileIndex, 
                  anyQueue = anyQueue)
  
  return(genoList)
}

updateSampleIDs = function(SampleIDs, samplesInGeno)
{
  if(is.null(SampleIDs)){
    cat("Since 'SampleIDs' is not specified, all samples in 'GenoFile' are used.\n")
    SampleIDs = samplesInGeno;
  }
  
  if(any(!SampleIDs %in% samplesInGeno)){
    print(SampleIDs[!SampleIDs %in% samplesInGeno])
    stop("The above samples from 'SampleIDs' are not in 'GenoFile' and 'GenoFileIndex'.")
  }
  
  SampleIDs = as.character(SampleIDs)
  
  return(SampleIDs)
}

# https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
#' Get sample identifiers from BGEN file
#' 
#' Extract sample identifiers from BGEN file (only support BGEN v1.2, check [link](https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html))
#' 
#' @param bgenFile a character of BGEN file. 
#' @examples
#' 
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "GRAB")
#' getSampleIDsFromBGEN(BGENFile)
#' @export
getSampleIDsFromBGEN = function(bgenFile)
{
  if(!checkIfSampleIDsExist(bgenFile))
    stop("The BGEN file does not include sample identifiers. Please refer to help(checkIfSampleIDsExist) for more details")
  
  cat("extracting sample information from bgen file\n")
  con = file(bgenFile, "rb")
  seek(con, 4)
  LH = readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH + 4)
  N = readBin(con, n = 1, what = "integer", size = 4)  # number of samples
  samplesInGeno = rep(0, N)
  
  # cycle for all samples to extract IDs
  for(i in 1:N){
    LS = readBin(con, n = 1, what = "integer", size = 2)
    sample = readChar(con, nchars = LS, useBytes = T)
    samplesInGeno[i] = sample
  }
  
  # close connection
  close(con)
  
  return(samplesInGeno)
}

# https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
#' Get version information from BGEN file
#' 
#' Get version information from BGEN file (check [link](https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html))
#' 
#' @param bgenFile a character of BGEN file. 
#' @examples
#' 
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "GRAB")
#' getVersionFromBGEN(BGENFile)
#' @export
getVersionFromBGEN = function(bgenFile)
{
  con = file(bgenFile, "rb")
  seek(con, 4)
  LH = readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH - 4)
  header = rawToBits(readBin(con, n = 4, what = "raw", size = 1, signed = FALSE))
  
  VersionNum = convert4BitsToNumber(header[3:6])
  if(VersionNum == 0)
    version = "Version Layout = 0, which is not supported. Please check https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html for more details."
  
  if(VersionNum == 1)
    version = "v1.1"
  
  if(VersionNum == 2)
    version = "v1.2"
  
  if(VersionNum > 2)
    version = "Version Layout > 2, which is reserved for future use. Please check https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html for more details."
  
  close(con)
  return(version)
}

# https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
convert4BitsToNumber = function(leastSignificantBit)
{
  leastSignificantBit = as.numeric(leastSignificantBit)
  if(length(leastSignificantBit) != 4)
    stop("Input should be 4 bits in which least-significant first.")
  
  Number  = 0;
  for(i in 1:4){
    Number = Number + 2^(i-1) * leastSignificantBit[i];
  }
  return(Number)
}

#' Check if sample identifiers are stored in a BGEN file
#' 
#' Check if sample identifiers are stored in a BGEN file, only support BGEN v1.2. Check [link](https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html) for more details.
#' 
#' @param bgenFile a character of BGEN file. Sometimes, BGEN file does not include sample IDs. This information can be extracted from BGEN file. Please refer to [link](https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html) for more details. 
#' @examples
#' 
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "GRAB")
#' checkIfSampleIDsExist(BGENFile)
#' @export
checkIfSampleIDsExist = function(bgenFile)
{
  con = file(bgenFile, "rb")
  seek(con, 4)
  LH = readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH - 4)
  header = rawToBits(readBin(con, n = 4, what = "raw", size = 1, signed = FALSE))
  close(con)
  return(header[32] == 01)
}


