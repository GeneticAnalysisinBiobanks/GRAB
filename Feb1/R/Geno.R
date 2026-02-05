
#' Read in genotype data
#' 
#' @description
#' Reads genotype data from PLINK or BGEN files into an R matrix.
#' This function supports filtering markers by ID or range and handling missing data imputation.
#' 
#' @param GenoFile Character. Path to the genotype file. Supported formats: PLINK (.bed) and BGEN (.bgen).
#' @param GenoFileIndex Character vector. Paths to index files. 
#' \itemize{
#'   \item For PLINK: c(bimFile, famFile). If NULL, defaults to replacing .bed with .bim/.fam.
#'   \item For BGEN: c(bgiFile, sampleFile). If NULL, defaults to replacing .bgen with .bgen.bgi/.sample.
#' }
#' @param SampleIDs Character vector. subset of sample IDs to read. If NULL, all samples are read.
#' @param control List. Control parameters for reading genotype data.
#' \itemize{
#'   \item \code{AlleleOrder}: "alt-first" (default for PLINK) or "ref-first" (default for BGEN).
#'   \item \code{ImputeMethod}: Method for missing genotype imputation: "none", "bestguess", or "mean" (default "none").
#'   \item \code{IDsToIncludeFile}: File containing marker IDs to include.
#'   \item \code{IDsToExcludeFile}: File containing marker IDs to exclude.
#'   \item \code{RangesToIncludeFile}: File containing genomic ranges to include (columns: CHROM, START, END).
#'   \item \code{RangesToExcludeFile}: File containing genomic ranges to exclude.
#'   \item \code{AllMarkers}: Logical. If TRUE, includes all markers when no filters are specified (default FALSE).
#' }
#' @param sparse Logical. If TRUE, returns a sparse matrix (default FALSE).
#' 
#' @return A list containing:
#' \item{GenoMat}{Genotype matrix (Rows=Samples, Cols=Markers).}
#' \item{markerInfo}{Data frame with marker information (CHROM, POS, ID, REF, ALT).}
#' 
#' @examples
#' \dontrun{
#' # Read PLINK data
#' # Geno = SPAmixPlus.ReadGeno(GenoFile = "test.bed", 
#' #                      control = list(ImputeMethod = "mean"))
#' }
#' 
#' @export
SPAmixPlus.ReadGeno = function(GenoFile,
                         GenoFileIndex = NULL,
                         SampleIDs = NULL,
                         control = NULL,
                         sparse = FALSE)
{
  control = checkControl.ReadGeno(control)  # check 'Framework.R'
  
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

#' Check if sample identifiers exist in BGEN files
#' 
#' @description
#' Verifies if the BGEN file contains embedded sample identifiers.
#' 
#' @param bgenFile Character. Path to the BGEN file.
#' @return Logical. TRUE if sample IDs exist, FALSE otherwise.
#' @examples
#' \dontrun{
#' # BGENFile = system.file("extdata", "simuBGEN.bgen", package = "SPAmixPlus")
#' # checkIfSampleIDsExist(BGENFile)
#' }
#' @export
checkIfSampleIDsExist = function(bgenFile)
{
  if(missing(bgenFile))
    stop("Argument 'bgenFile' is required.")
  
  if(!file.exists(bgenFile))
    stop("Cannot find bgen file of ", bgenFile, ".")
  
  checkIfSampleIDsExistInBGEN(bgenFile)  # this function is in BGEN.cpp
}


updateSampleIDs = function(SampleIDs, samplesInGeno)
{
  if(!is.null(SampleIDs))
  {
    SampleIDs = as.character(unique(SampleIDs))
    pos = which(!SampleIDs %in% samplesInGeno)
    if(length(pos) != 0){
      message = paste("Among the required", length(SampleIDs), "SampleIDs,", length(pos), "samples are not in the GenoFile.")
      warning(message)
      # message = "The first 5 samples that are not in the GenoFile:"
      # print(head(SampleIDs[pos], 5))
    }
    
    SampleIDs = intersect(SampleIDs, samplesInGeno)
    if(length(SampleIDs) == 0) stop("No samples are extracted from the GenoFile.")
  }else{
    SampleIDs = samplesInGeno
  }
  
  return(SampleIDs)
}

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
    
    # if(getVersionFromBGEN(bgenFile) != "v1.2")
    #   stop("Package SPAmixPlus currently only supports version 1.2 of BGEN format.")
    
    if(is.null(AlleleOrder)) AlleleOrder = "ref-first"
    
    if(is.null(GenoFileIndex)){  
      # If 'GenoFileIndex' is not given, we use the same prefix for 'bgen.bgi' file
      GenoFileIndex = c(gsub("bgen$", "bgen.bgi", GenoFile),
                        gsub("bgen$", "sample", GenoFile))
    }
    
    if(length(GenoFileIndex) != 1 & length(GenoFileIndex) != 2)
      stop("For genotype input of BGEN format, 'GenoFileIndex' should be of length 1 or 2. Check 'Details' section for more details.")
    
    if(length(GenoFileIndex) == 1){
      samplesInGeno = getSampleIDsFromBGEN(bgenFile)
      if(samplesInGeno[[1]] == "0")   # check BGEN.cpp
        stop("The BGEN file does not include sample identifiers. Please provide a .sample file.")
    }
    
    if(length(GenoFileIndex) == 2){
      # https://www.cog-genomics.org/plink/2.0/formats#sample
      sampleFile = GenoFileIndex[2]
      if(!file.exists(sampleFile)) stop(paste("Cannot find sample file of", sampleFile))
      
      cat("Reading sample file:\t", sampleFile, "\n")
      sampleInfo = data.table::fread(sampleFile, header = T, sep = " ")
      if(ncol(sampleInfo) == 1)
        sampleInfo = data.table::fread(sampleFile, header = T, sep = "\t")
      
      # The first two lines of the file are header lines. 
      # The first line has column names (ID_1 ID_2 missing sex) 
      # The second line is 0 0 0 D (type of the field)
      
      samplesInGeno = sampleInfo$ID_2[-1]
    }
    
    SampleIDs = updateSampleIDs(SampleIDs, samplesInGeno)
    
    bgiFile = GenoFileIndex[1]
    if(!file.exists(bgiFile)) stop(paste("Cannot find bgi file of", bgiFile))
    
    cat("Reading bgi file:\t", bgiFile, "\n")
    # https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md
    
    # The bgi file is a sqlite3 database
    con = RSQLite::dbConnect(RSQLite::SQLite(), bgiFile)
    markerInfo = RSQLite::dbGetQuery(con, 'select * from Variant')
    # chromosome	position	rsid	number_of_alleles	allele1	allele2	file_start_position	size_in_bytes
    RSQLite::dbDisconnect(con)
    
    if(AlleleOrder == "alt-first")
      markerInfo = markerInfo[,c("chromosome", "position", "rsid", "allele2", "allele1")]
    if(AlleleOrder == "ref-first")
      markerInfo = markerInfo[,c("chromosome", "position", "rsid", "allele1", "allele2")]
    
    colnames(markerInfo) = c("CHROM", "POS", "ID", "REF", "ALT")
    markerInfo$genoIndex = 1:nrow(markerInfo) - 1  # -1 is to convert 'R' to 'C++' 
    
    cat("setting BGEN object in CPP....\n")
    setBGENobjInCPP(bgenFile, bgiFile, samplesInGeno, SampleIDs, FALSE, FALSE, AlleleOrder)
  }
  
  # make an interaction 
  IDsToIncludeFile = control$IDsToIncludeFile
  IDsToExcludeFile = control$IDsToExcludeFile
  RangesToIncludeFile = control$RangesToIncludeFile
  RangesToExcludeFile = control$RangesToExcludeFile
  
  anyQueue = FALSE
  if(!is.null(IDsToIncludeFile)){
    anyQueue = TRUE
    IDsToInclude = data.table::fread(IDsToIncludeFile, header = F)
    markerInfo = markerInfo[markerInfo$ID %in% IDsToInclude$V1,]
  }
  
  if(!is.null(IDsToExcludeFile)){
    anyQueue = TRUE
    IDsToExclude = data.table::fread(IDsToExcludeFile, header = F)
    markerInfo = markerInfo[!markerInfo$ID %in% IDsToExclude$V1,]
  }
  
  if(!is.null(RangesToIncludeFile)){
    anyQueue = TRUE
    RangesToInclude = data.table::fread(RangesToIncludeFile, header = F)
    # 3 columns: CHROM, START, END
    posSeries = 1:nrow(markerInfo)
    posKeep = c()
    for(i in 1:nrow(RangesToInclude)){
      chrom = RangesToInclude$V1[i]
      start = RangesToInclude$V2[i]
      end = RangesToInclude$V3[i]
      pos = posSeries[markerInfo$CHROM == chrom & markerInfo$POS >= start & markerInfo$POS <= end]
      posKeep = c(posKeep, pos)
    }
    posKeep = sort(unique(posKeep))
    markerInfo = markerInfo[posKeep,]
  }

  if(!is.null(RangesToExcludeFile)){
    anyQueue = TRUE
    RangesToExclude = data.table::fread(RangesToExcludeFile, header = F)
    # 3 columns: CHROM, START, END
    posSeries = 1:nrow(markerInfo)
    posRemove = c()
    for(i in 1:nrow(RangesToExclude)){
      chrom = RangesToExclude$V1[i]
      start = RangesToExclude$V2[i]
      end = RangesToExclude$V3[i]
      pos = posSeries[markerInfo$CHROM == chrom & markerInfo$POS >= start & markerInfo$POS <= end]
      posRemove = c(posRemove, pos)
    }
    posRemove = sort(unique(posRemove))
    if(length(posRemove) > 0)
      markerInfo = markerInfo[-posRemove,]
  }
  
  return(list(genoType = genoType,
              markerInfo = markerInfo,
              SampleIDs = SampleIDs,
              anyQueue = anyQueue))
}


splitMarker = function(genoIndex, nMarkersEachChunk, CHROM)
{
  # genoIndex is a vector of logic indices of markers to include in analysis
  # nMarkersEachChunk is an integer, number of markers in each chunk
  
  # split markers into chunks
  nMarkers = length(genoIndex)
  nChunks = ceiling(nMarkers / nMarkersEachChunk)
  
  # Group markers by chromosome
  # Markers in the same chromosome should be in the same chunk to avoid random access to genotype file
  
  if(length(unique(CHROM)) == 1){
    genoIndexList = split(genoIndex, ceiling(seq_along(genoIndex)/nMarkersEachChunk))
    
    # re-organize the list structure
    genoIndexList = lapply(genoIndexList, function(x){
      list(genoIndex = x, chrom = unique(CHROM))
    })
  }else{
    # if markers are from multiple chromosomes, we process them chromosome by chromosome
    
    # split marker indices by chromosome
    genoIndexList0 = split(genoIndex, CHROM)
    
    genoIndexList = list()
    for(i in 1:length(genoIndexList0)){
      x = genoIndexList0[[i]]
      chrom = unique(CHROM[genoIndex %in% x])
      
      # further split by chunk size
      subList = split(x, ceiling(seq_along(x)/nMarkersEachChunk))
      
      for(j in 1:length(subList)){
        genoIndexList[[length(genoIndexList)+1]] = list(genoIndex = subList[[j]], chrom = chrom)
      }
    }
  }
  
  return(genoIndexList)
}
