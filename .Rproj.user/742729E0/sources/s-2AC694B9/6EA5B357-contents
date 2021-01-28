
## R codes to read in genotype data
## The package supports the following genotype format
## 1. plink files: bedFile as GenoFile, c(bimFile, famFile) as GenoFileIndex

#' Read genotype data into R
#' 
#' GRAB package provides functions to read genotype data from multiple formats (including plink) into R
#' 
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". The default is NULL, that is, to share the same prefix as GenoFile. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param SampleIDs a character vector of sample IDs to extract. The default is NULL, that is, to use the all samples in GenoFile.
#' @param MarkerIDs a character vector of marker IDs to extract. The default is NULL, the function will extract the first 10 markers.
#' @return An R matrix, each row is for one sample and each column is for one marker.
#' @examples
#' GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
#' GenoMat = readGeno(GenoFile)
#' head(GenoMat)
#'      
#' @export
#' @import data.table
readGeno = function(GenoFile,
                    GenoFileIndex = NULL,
                    SampleIDs = NULL,
                    MarkerIDs = NULL)
{
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs)
  
  genoType = objGeno$genoType
  markers = objGeno$markers
  samples = objGeno$samples
  
  if(is.null(SampleIDs)){
    print("Since 'SampleIDs' not specified, we use all samples in 'GenoFile'.")
    SampleIDs = samples 
  }
  
  if(is.null(MarkerIDs)){
    print("Since 'MarkerIDs' not specified, we use the first 10 markers in 'GenoFile'.")
    MarkerIDs = markers[1:min(10,length(markers))]
  }
  
  if(any(!SampleIDs %in% samples))
    stop("At least one sample from 'SampleIDs' are not in 'GenoFile' and 'GenoFileIndex'.")
  
  if(any(!MarkerIDs %in% markers))
    stop("At least one marker from 'MarkerIDs' are not in 'GenoFile' and 'GenoFileIndex'.")
  
  GenoMat = getGenoInCPP(MarkerIDs)  # for more details about getGenoInCPP, please check Main.cpp
  colnames(GenoMat) = MarkerIDs;
  rownames(GenoMat) = SampleIDs;
  return(GenoMat)
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
  
  if(GenoFileExt != "bed")
    stop("Current version only supports plink file input (file extension should be .bed).")  # support more genotype input later
  
  ########## ----------  Plink format ---------- ##########
  
  if(GenoFileExt == "bed"){
    genoType = "PLINK"
    if(is.null(GenoFileIndex)){  
      # If GenoFileIndex is not given, we use the same prefix for bim and fam files
      GenoFileIndex = c(gsub("bed$","bim",GenoFile),
                        gsub("bed$","fam",GenoFile))
    }
    if(length(GenoFileIndex) != 2)
      stop("If plink file is used, argument 'GenoFileIndex' should be 'NULL' or a character vector of c(bimFile, famFile).")
    
    bimFile = GenoFileIndex[1]
    famFile = GenoFileIndex[2]
    
    if(!file.exists(bimFile)) stop(paste("Cannot find bim file of", bimFile))
    markerInfo = data.table::fread(bimFile)
    markers = markerInfo$V2
    
    if(!file.exists(famFile)) stop(paste("Cannot find fam file of", famFile))
    sampleInfo = data.table::fread(famFile)
    samples = sampleInfo$V2
    
    if(is.null(SampleIDs))
      SampleIDs = samples;
    
    setPLINKobjInCPP(bimFile, famFile, GenoFile, SampleIDs)
  }
  
  ########## ----------  More formats such as BGEN and VCF ---------- ##########
  
  # return genotype
  print(paste("Based on the 'GenoFile' and 'GenoFileIndex',", genoType, "format is used for genotype data."))
  
  genoList = list(genoType = genoType, markers = markers, samples = samples)
  
  return(genoList)
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


