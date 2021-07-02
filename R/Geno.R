
## R codes to read in genotype data
## The package supports the following genotype format
## 1. plink files: bedFile as GenoFile, c(bimFile, famFile) as GenoFileIndex

#' Read in genotype data
#' 
#' GRAB package provides functions to read in genotype data from multiple format (including Plink, BGEN, and VCF) into R
#' 
#' @param GenoFile a character of genotype file. See \code{Details} section for more information.
#' @param GenoFileIndex additional index file(s) corresponding to the \code{GenoFile}. See \code{Details} section for more information.
#' @param SampleIDs a character vector of sample IDs to extract. The default is NULL, that is, to use all samples in GenoFile.
#' @param control a list of parameters to decide which markers to extract. If not specified, the first 10 markers will be extracted. For more details, please check \code{?GRAB.control}.
#' @param sparse a logical value (default: FALSE) to indicate if sparse genotype matrix is outputted.
#' @return An R list include an R genotype matrix (each row is for one sample and each column is for one marker) and an R SNP information matrix.
#' @details
#' We support three genotype formats including Plink, BGEN, and VCF. 
#' Users do not need to specify the genotype format, GRAB package will check the filename extention for that purpose.  
#' If \code{GenoFileIndex} is not specified, then GRAB uses the same prefix as \code{GenoFile}.
#' \describe{
#'   \item{Plink}{
#'   \itemize{
#'   \item \code{GenoFile}: "prefix.bed". 
#'   \item \code{GenoFileIndex}: c("prefix.bim", "prefix.fam").
#'   }
#'   }
#'   \item{BGEN}{
#'     \itemize{
#'     \item \code{GenoFile}: "prefix.bgen"; 
#'     \item \code{GenoFileIndex}: "prefix.bgen.bgi" or c("prefix.bgen.bgi", "prefix.bgen.samples").
#'     \item IMPORTANT NOTE: If only one element is given for \code{GenoFileIndex}, then we assume it should be "prefix.bgen.bgi". 
#'     If BGEN file does not include sample identifiers, then "prefix.bgen.samples" is required, which should be a file with only one column whose column name is "GRAB_BGEN_SAMPLE" (case insensitive). 
#'     One example is \code{system.file("extdata", "example_bgen_1.2_8bits.bgen.samples", package = "GRAB")}.
#'     If you are not sure if sample identifiers are in BGEN file, you can try function \code{?checkIfSampleIDsExist}.
#'     }
#'   }
#'   \item{VCF}{Not available now. \code{GenoFile}: "prefix.vcf"; \code{GenoFileIndex}: "prefix.vcf.tbi"}
#' }
#' 
#' @examples
#' 
#' ## The below is raw data 
#' RawFile = system.file("extdata", "example.raw", package = "GRAB")
#' GenoMat = data.table::fread(RawFile)
#' head(GenoMat[,1:15])
#' 
#' ## The below is for Plink format
#' PlinkFile = system.file("extdata", "example.bed", package = "GRAB")
#' GenoList = GRAB.ReadGeno(PlinkFile)
#' GenoMat = GenoList$GenoMat
#' markerInfo = GenoList$markerInfo
#' head(GenoMat)
#' markerInfo
#' 
#' ## The below is for BGEN format (Note the different REF/ALT order for BGEN and Plink formats)
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "GRAB")
#' GenoList = GRAB.ReadGeno(BGENFile)
#' GenoMat = GenoList$GenoMat
#' markerInfo = GenoList$markerInfo
#' head(GenoMat)
#' markerInfo
#' 
#' ## The below is to demonstrate parameters in control
#' PlinkFile = system.file("extdata", "example.bed", package = "GRAB")
#' IDsToIncludeFile = system.file("extdata", "example.IDsToIncludeFile.txt", package = "GRAB")
#' RangesToIncludeFile = system.file("extdata", "example.RangesToIncludeFile.txt", package = "GRAB")
#' GenoList = GRAB.ReadGeno(PlinkFile, 
#'                          control = list(IDsToIncludeFile = IDsToIncludeFile, 
#'                                         RangesToIncludeFile = RangesToIncludeFile,
#'                                         AlleleOrder = "ref-first"))
#' GenoMat = GenoList$GenoMat
#' head(GenoMat)
#' markerInfo = GenoList$markerInfo
#' markerInfo
#' 
#' @export
#' @import data.table, tidyr, dbplyr, RSQLite, Matrix
GRAB.ReadGeno = function(GenoFile,
                         GenoFileIndex = NULL,
                         SampleIDs = NULL,
                         control = NULL,
                         sparse = FALSE)
{
  checkControl.ReadGeno(control)  # this function is in 'control.R': indexGeno can be 0, 1, 2, 3, 4 depending on which argument is given.
  
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs, control)
  
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  SampleIDs = objGeno$SampleIDs
  anyQueue = objGeno$anyQueue   # if FALSE, no include/exclude is specified
  
  if(!anyQueue){
    print("Since no markers or regions were selected, we use the first 10 markers in 'GenoFile'.")
    markerInfo = markerInfo[1:min(10,nrow(markerInfo)),]
  }
  
  MarkerIDs = markerInfo$ID
  
  n = length(SampleIDs)
  m = length(MarkerIDs)
  
  cat("Number of Samples:\t",n,"\nNumber of Markers:\t",m)
  
  if(sparse == TRUE){
    GenoMat = getSpGenoInCPP(genoType, markerInfo, n)
  }else{
    GenoMat = getGenoInCPP(genoType, markerInfo, n)  # This function is in Main.cpp
  }
    
  colnames(GenoMat) = MarkerIDs;
  rownames(GenoMat) = SampleIDs;
  
  return(list(GenoMat = GenoMat,
              markerInfo = markerInfo))
}

## add something in setGenoInput(.) to select specific markers requested by users
setGenoInput = function(GenoFile, 
                        GenoFileIndex = NULL, 
                        SampleIDs = NULL,
                        control = NULL)
{
  if(missing(GenoFile))
    stop("Argument 'GenoFile' is required.")
  
  if(!file.exists(GenoFile))
    stop(paste("Cannot find genotype file of", GenoFile))
  
  GenoFileExt = tools::file_ext(GenoFile);
  
  # support more genotype format later
  if(GenoFileExt != "bed" & GenoFileExt != "bgen")
    stop("Current version only supports genotype input of Plink (filename extension of '.bed') or BGEN (filename extension of '.bgen').")  
  
  AlleleOrder = control$AlleleOrder
  ########## ----------  Plink format ---------- ##########
  
  if(GenoFileExt == "bed"){
    
    genoType = "PLINK"
    
    if(is.null(AlleleOrder)) AlleleOrder = "alt-first"
    
    if(is.null(GenoFileIndex)){  
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bim' and 'fam' files
      GenoFileIndex = c(gsub("bed$", "bim", GenoFile),
                        gsub("bed$", "fam", GenoFile))
    }
    
    if(length(GenoFileIndex) != 2)
      stop("If Plink format is used, argument 'GenoFileIndex' should be 'NULL' or a character vector of c(bimFile, famFile).")
    
    bimFile = GenoFileIndex[1]
    famFile = GenoFileIndex[2]
    bedFile = GenoFile
    
    if(!file.exists(bimFile)) stop(paste("Cannot find bim file of", bimFile))
    markerInfo = data.table::fread(bimFile, header = F)
    markerInfo = as.data.frame(markerInfo)
    
    if(AlleleOrder == "alt-first")
      markerInfo = markerInfo[,c(1,4,2,6,5)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    if(AlleleOrder == "ref-first")
      markerInfo = markerInfo[,c(1,4,2,5,6)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    
    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT")
    markerInfo$genoIndex = 1:nrow(markerInfo) - 1  # -1 is to convert 'R' to 'C++' 
    
    if(!file.exists(famFile)) stop(paste("Cannot find fam file of", famFile))
    sampleInfo = data.table::fread(famFile)
    samplesInGeno = sampleInfo$V2
    SampleIDs = updateSampleIDs(SampleIDs, samplesInGeno)
      
    setPLINKobjInCPP(bimFile, famFile, bedFile, SampleIDs, AlleleOrder)
  }
  
  ########## ----------  BGEN format ---------- ##########
  
  if(GenoFileExt == "bgen"){
    
    genoType = "BGEN"
    bgenFile = GenoFile
    
    if(is.null(AlleleOrder)) AlleleOrder = "ref-first"
    
    if(is.null(GenoFileIndex)){  
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bgen.bgi' file
      GenoFileIndex = c(gsub("bgen$", "bgen.bgi", GenoFile),
                        gsub("bgen$", "bgen.samples", GenoFile))
    }
    
    if(length(GenoFileIndex) != 1 & length(GenoFileIndex) != 2)
      stop("For genotype input of BGEN format , 'GenoFileIndex' should be of length 1 or 2. Check 'Details' section in '?GRAB.ReadGeno' for more information.")
    
    if(length(GenoFileIndex) == 1){
      samplesInGeno = getSampleIDsFromBGEN(bgenFile)
    }
      
    if(length(GenoFileIndex) == 2){
      sampleFile = GenoFileIndex[2]
      if(!file.exists(sampleFile)){
        if(!checkIfSampleIDsExist(bgenFile)){
          stop(paste("Cannot find bgen.samples file of", sampleFile))
        }else{
          samplesInGeno = getSampleIDsFromBGEN(bgenFile)
        }
      }else{
        sampleData = data.table::fread(sampleFile, header=T, colClasses = c("character"))
        if(toupper(colnames(sampleData)[1]) != "GRAB_BGEN_SAMPLE")
          stop("The header of the first column in bgen.samples file should be 'GRAB_BGEN_SAMPLE'.")
        samplesInGeno = as.character(sampleData[,1])
      }
    }
    
    bgiFile = GenoFileIndex[1]
    
    if(!file.exists(bgiFile)) stop(paste("Cannot find bgi file of", bgiFile))
    
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
  
  ########## ----------  More format such as VCF ---------- ##########

  Files = c("IDsToIncludeFile", "IDsToExcludeFile", "RangesToIncludeFile", "RangesToExcludeFile")
  
  anyInclude = FALSE
  anyExclude = FALSE
  
  markersInclude = c()
  markersExclude = c()
  
  if(!is.null(control$IDsToIncludeFile)){
    IDsToInclude = data.table::fread(control$IDsToIncludeFile, 
                                     header = F, colClasses = c("character"))
    if(ncol(IDsToInclude) != 1)
      stop("IDsToIncludeFile should include one column.")
    IDsToInclude = IDsToInclude[,1]
    
    posRows = which(markerInfo$ID %in% IDsToInclude)
    if(length(posRows) != 0) 
      markersInclude = c(markersInclude, markerInfo$ID[posRows])
    anyInclude = TRUE
  }
  
  if(!is.null(control$RangesToIncludeFile)){
    RangesToInclude = data.table::fread(control$RangesToIncludeFile, 
                                        header = F, colClasses = c("character", "numeric", "numeric"))
    if(ncol(RangesToInclude) != 3)
      stop("RangesToIncludeFile should include three columns.")
    
    colnames(RangesToInclude) = c("CHROM","START","END")
    
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
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile'.")
    IDsToExclude = data.table::fread(control$IDsToExcludeFile, 
                                     header = F, colClasses = c("character"))
    if(ncol(IDsToExclude) != 1)
      stop("IDsToExcludeFile should include one column.")
    IDsToExclude = IDsToExclude[,1]
    
    posRows = which(markerInfo$ID %in% IDsToExclude)
    if(length(posRows) != 0) 
      markersExclude = c(markersExclude, markerInfo$ID[posRows])
    anyExclude = TRUE
  }
  
  if(!is.null(control$RangesToExcludeFile)){
    if(anyInclude) 
      stop("We currently do not support both 'IncludeFile' and 'ExcludeFile'.")
    
    RangesToExclude = data.table::fread(control$RangesToExcludeFile, 
                                        header = F, colClasses = c("character", "numeric", "numeric"))
    if(ncol(RangesToExclude) != 3)
      stop("RangesToExcludeFile should include three columns.")
    
    colnames(RangesToExclude) = c("CHROM","START","END")
    
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
  print(paste("Based on the 'GenoFile' and 'GenoFileIndex',", genoType, "format is used for genotype data."))
  
  if(anyInclude)
    markerInfo = subset(markerInfo, ID %in% markersInclude)
  
  if(anyExclude)
    markerInfo = subset(markerInfo, !ID %in% markersExclude)
  
  anyQueue = anyInclude | anyExclude
  
  genoList = list(genoType = genoType, markerInfo = markerInfo, SampleIDs = SampleIDs, AlleleOrder = AlleleOrder, GenoFile = GenoFile, GenoFileIndex = GenoFileIndex, anyQueue = anyQueue)
  
  return(genoList)
}

updateSampleIDs = function(SampleIDs, samplesInGeno)
{
  if(is.null(SampleIDs)){
    print("Since 'SampleIDs' not specified, we use all samples in 'GenoFile'.")
    SampleIDs = samplesInGeno;
  }
  
  if(any(!SampleIDs %in% samplesInGeno))
    stop("At least one sample from 'SampleIDs' are not in 'GenoFile' and 'GenoFileIndex'.")
  
  return(SampleIDs)
}

# https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
#' Extract sample identifiers from BGEN file (for BGEN v1.2)
#' 
#' Extract sample identifiers from BGEN file (for BGEN v1.2)
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
    stop("The BGEN file does not include sample identifiers. Please refer to ?checkIfSampleIDsExist and ?GRAB.ReadGeno for more details")
  con = file(bgenFile, "rb")
  seek(con, 4)
  LH = readBin(con, n = 1, what = "integer", size = 4)
  seek(con, 4 + LH + 4)
  N = readBin(con, n = 1, what = "integer", size = 4)  # number of samples
  samplesInGeno = rep(0, N)
  
  # cycle for all samples to extract IDs
  for(i in 1:N){
    LS = readBin(con, n = 1, what = "integer", size = 2)
    sample = readChar(con, nchars = LS)
    samplesInGeno[i] = sample
  }
  
  # close connection
  close(con)
  
  return(samplesInGeno)
}

#' Check if sample identifiers are stored in a BGEN file (for BGEN v1.2)
#' 
#' Check if sample identifiers are stored in a BGEN file (for BGEN v1.2)
#' 
#' @param bgenFile a character of BGEN file. Sometimes, BGEN file does not include sample IDs. This information can be extracted from BGEN file. Please refer to https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html for more details. 
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


