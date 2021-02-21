
## R codes to read in genotype data
## The package supports the following genotype format
## 1. plink files: bedFile as GenoFile, c(bimFile, famFile) as GenoFileIndex

#' Read genotype data into R as a matrix
#' 
#' GRAB package provides functions to read in genotype data from multiple format (including Plink, BGEN) into R
#' 
#' @param GenoFile a character of genotype file. See \code{Details} section for more information.
#' @param GenoFileIndex additional index file(s) corresponding to the \code{GenoFile}. See \code{Details} section for more information.
#' @param SampleIDs a character vector of sample IDs to extract. The default is NULL, that is, to use all samples in GenoFile.
#' @param MarkerIDs a character vector of marker IDs to extract. The default is NULL, the first 10 markers will be extracted.
#' @return An R list include an R genotype matrix (each row is for one sample and each column is for one marker) and an R SNP information matrix.
#' @details
#' We support three genotype format including Plink, BGEN, and VCF.
#' The program will check the format based on the filename extension.  
#' If \code{GenoFileIndex} is NULL (default), then it uses the same prefix as \code{GenoFile}.
#' \describe{
#'   \item{Plink}{\code{GenoFile}: "prefix.bed"; \code{GenoFileIndex}: c("prefix.bim", "prefix.fam")}
#'   \item{BGEN}{\code{GenoFile}: "prefix.bgen"; \code{GenoFileIndex}: "prefix.bgen.bgi"}
#'   \item{VCF}{Not available now. \code{GenoFile}: "prefix.vcf"; \code{GenoFileIndex}: "prefix.vcf.tbi"}
#' }
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
#' ## The below is for BGEN format
#' BGENFile = system.file("extdata", "example_bgen_1.2_8bits.bgen", package = "GRAB")
#' GenoList = GRAB.ReadGeno(BGENFile)
#' GenoMat = GenoList$GenoMat
#' markerInfo = GenoList$markerInfo
#' head(GenoMat)
#' markerInfo
#' 
#' @export
#' @import data.table, tidyr, dbplyr, RSQLite
GRAB.ReadGeno = function(GenoFile,
                         GenoFileIndex = NULL,
                         SampleIDs = NULL,
                         MarkerIDs = NULL)
{
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs)
  
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  SampleIDs = objGeno$SampleIDs
  
  if(is.null(MarkerIDs)){
    print("Since 'MarkerIDs' not specified, we use the first 10 markers in 'GenoFile'.")
    markerInfo = markerInfo[1:min(10,nrow(markerInfo)),]
    MarkerIDs = markerInfo$ID
  }
  
  posMarker = match(MarkerIDs, markerInfo$ID, 0)
  if(any(posMarker == 0))
    stop("At least one marker from 'MarkerIDs' are not in 'GenoFile' and 'GenoFileIndex'.")
  
  markerInfo = markerInfo[posMarker, ,drop=F]
  
  if(genoType == "PLINK")
  {
    GenoMat = getGenoInCPP(genoType, 
                           list(MarkerReqstd = MarkerIDs,
                                n = length(SampleIDs),
                                q = length(posMarker)))  # for more details about getGenoInCPP, please check Main.cpp
  }
  if(genoType == "BGEN")
  {
    GenoMat = getGenoInCPP(genoType, 
                           list(fileStartPosVec = markerInfo$StartPositionInBGEN,
                                n = length(SampleIDs),
                                q = length(posMarker)))  # for more details about getGenoInCPP, please check Main.cpp
  }
  
  colnames(GenoMat) = MarkerIDs;
  rownames(GenoMat) = SampleIDs;
  
  return(list(GenoMat = GenoMat,
              markerInfo = markerInfo))
}

## add something in setGenoInput(.) to select specific markers requested by users
setGenoInput = function(GenoFile, 
                        GenoFileIndex = NULL, 
                        SampleIDs = NULL)
{
  if(missing(GenoFile))
    stop("Argument 'GenoFile' is required.")
  
  if(!file.exists(GenoFile))
    stop(paste("Cannot find genotype file of", GenoFile))
  
  GenoFileExt = tools::file_ext(GenoFile);
  
  # support more genotype format later
  if(GenoFileExt != "bed" & GenoFileExt != "bgen")
    stop("Current version only supports genotype input of Plink (filename extension of '.bed') or BGEN (filename extension of '.bgen').")  
  
  ########## ----------  Plink format ---------- ##########
  
  if(GenoFileExt == "bed"){
    
    genoType = "PLINK"
    
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
    markerInfo = data.table::fread(bimFile)
    markerInfo = as.data.frame(markerInfo)
    markerInfo = markerInfo[,c(1,4,2,6,5)]  # https://www.cog-genomics.org/plink/2.0/formats#bim
    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT")
    
    if(!file.exists(famFile)) stop(paste("Cannot find fam file of", famFile))
    sampleInfo = data.table::fread(famFile)
    samplesInGeno = sampleInfo$V2
    SampleIDs = updateSampleIDs(SampleIDs, samplesInGeno)
      
    setPLINKobjInCPP(bimFile, famFile, bedFile, SampleIDs)
  }
  
  ########## ----------  BGEN format ---------- ##########
  
  if(GenoFileExt == "bgen"){
    
    genoType = "BGEN"
    
    if(is.null(GenoFileIndex)){  
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bgen.bgi' file
      GenoFileIndex = gsub("bgen$", "bgen.bgi", GenoFile)
    }
    
    if(length(GenoFileIndex) != 1)
      stop("If BGEN format is used, argument 'GenoFileIndex' should be 'NULL' or a character of bgiFile.")
    
    bgiFile = GenoFileIndex[1]
    bgenFile = GenoFile
    
    if(!file.exists(bgiFile)) stop(paste("Cannot find bgi file of", bgiFile))
    
    db_con <- RSQLite::dbConnect(RSQLite::SQLite(), bgiFile)
    on.exit(RSQLite::dbDisconnect(db_con), add = TRUE)
    bgiData = dplyr::tbl(db_con, "Variant")
    bgiData = as.data.frame(bgiData)
    
    markerInfo = bgiData[,c(1,2,3,6,5,7)]  # https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT","StartPositionInBGEN")
    
    samplesInGeno = getSampleIDsFromBGEN(bgenFile)
    SampleIDs = updateSampleIDs(SampleIDs, samplesInGeno)
    
    setBGENobjInCPP(bgenFile, bgiFile, samplesInGeno, SampleIDs, F, F)
  }
  
  ########## ----------  More format such as BGEN and VCF ---------- ##########
  
  
  
  # return genotype
  print(paste("Based on the 'GenoFile' and 'GenoFileIndex',", genoType, "format is used for genotype data."))
  
  genoList = list(genoType = genoType, markerInfo = markerInfo, SampleIDs = SampleIDs)
  
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

## split 'markers' into multiple chunks each of which includes no more than 'nMarkersEachChunk' markers
# Examples to check:
# splitMarker(1:10, 2)
# splitMarker(1:10, 3)
# splitMarker(1:10, 10)
# splitMarker(1:10, 11)
splitMarker = function(markers, nMarkersEachChunk)
{
  M = length(markers)
  
  idxStart = seq(1, M, nMarkersEachChunk)
  idxEnd = idxStart + nMarkersEachChunk - 1
  
  nChunks = length(idxStart)
  idxEnd[nChunks] = M
  
  markerList = list()
  for(i in 1:nChunks){
    idxMarker = idxStart[i]:idxEnd[i]
    markerList[[i]] = markers[idxMarker]
  }
  
  return(markerList)
}

# https://www.well.ox.ac.uk/~gav/bgen_format/spec/v1.2.html
getSampleIDsFromBGEN = function(bgenFile)
{
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
  
  return(samplesInGeno)
}




