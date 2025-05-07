
#' Make temporary files to be passed to function \code{\link{getSparseGRM}}. 
#' 
#' Make temporary files to be passed to function \code{\link{getSparseGRM}}. We strongly suggest using parallel computing for different \code{partParallel}.
#' @param PlinkFile a path to PLINK files (without file extensions of bed/bim/fam). Note that the current version (gcta_1.93.1beta) of gcta software does not support different prefix names for bim, bed, and fam files. 
#' @param nPartsGRM a numeric value (e.g. 250): \code{GCTA} software can split subjects to multiple parts. For UK Biobank data analysis, it is recommended to set \code{nPartsGRM=250}. 
#' @param partParallel a numeric value (from 1 to \code{nPartsGRM}) to split all jobs for parallel computation.
#' @param tempDir a path to store temp files to be passed to \code{\link{getSparseGRM}}. This should be consistent to the input of \code{\link{getSparseGRM}}. Default is system.file("SparseGRM", "temp", package = "GRAB").
#' @param subjData a character vector to specify subject IDs to retain (i.e. IID). Default is \code{NULL}, i.e. all subjects are retained in sparse GRM. If the number of subjects is less than 1,000, the GRM estimation might not be accurate.
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from PLINK files) to make sparse GRM. *(default=0.01)*
#' @param maxMissingGRM Maximal value of missing rate to select markers (from PLINK files) to make sparse GRM. *(default=0.1)*
#' @param threadNum Number of threads (CPUs) to use.
#' @details 
#' \describe{
#' \itemize{
#'   \item \code{Step 1}: Run \code{getTempFilesFullGRM} to get temporary files.
#'   \item \code{Step 2}: Run \code{\link{getSparseGRM}} to combine the temporary files to make a \code{SparseGRMFile} to be passed to \code{\link{GRAB.NullModel}}.
#'   }
#' }
#' @examples 
#' ## Please check help(getSparseGRM) for an example.
#' @export
 
getTempFilesFullGRM = function(PlinkFile,
                               nPartsGRM,
                               partParallel,
                               tempDir = NULL,
                               subjData = NULL,
                               minMafGRM = 0.01,
                               maxMissingGRM = 0.1,
                               threadNum = 8)
{
  bimFile = paste0(PlinkFile,".bim")
  bedFile = paste0(PlinkFile,".bed")
  famFile = paste0(PlinkFile,".fam")
  
  if(!file.exists(bimFile)) stop("Could not find bimFile or paste0(PlinkFile,'.bim')")
  if(!file.exists(bedFile)) stop("Could not find bedFile or paste0(PlinkFile,'.bed')")
  if(!file.exists(famFile)) stop("Could not find famFile or paste0(PlinkFile,'.fam')")
  
  OS = Sys.info()['sysname']
  if(OS != "Linux") 
    stop(paste0("We generate GRM using gcta software which only support Linux. Current OS is ", OS))
  
  gcta64File = system.file("gcta_1.93.1beta", "gcta64", package = "GRAB");
  system(paste("chmod +x",gcta64File))
  
  if(is.null(tempDir)) 
    tempDir = system.file("SparseGRM", "temp", package = "GRAB")
  
  PlinkName = basename(PlinkFile)
  tempFile = paste0(tempDir, "/Plink-",PlinkName,"-autosome-minMaf-",minMafGRM,"-maxMissing-",maxMissingGRM)
  
  cmd = paste(gcta64File, 
              "--bfile", PlinkFile,
              "--out", tempFile,
              "--autosome",
              "--make-grm-part", nPartsGRM, partParallel,
              "--maf", minMafGRM,
              "--geno", maxMissingGRM,
              "--thread-num",threadNum)
  
  ## only retain parts of subjects
  if(!is.null(subjData)){
    if(length(subjData) < 1000) 
      stop("length(subjData) < 1000, the MAF estimate might be inaccurate.")
    
    subjFile = paste0(tempDir,"/subjData.txt")
    famData = read.table(famFile)
    posSubj = match(subjData, famData$V2, 0)
    if(any(posSubj==0)) 
      stop("All subjects in 'subjData' should be in IID column of 'famFile'.")
    
    write.table(famData[posSubj,c(1,2)], 
                subjFile, row.names=F, col.names=F, quote=F)
    cmd = paste(cmd, 
                "--keep", subjFile)
  }
  
  system(cmd)
  
  message = paste0("Temp files of Full GRM have been saved to ", tempFile)
  return(message)
}

#' Make a \code{SparseGRMFile} for \code{\link{GRAB.NullModel}}. 
#' 
#' If the sample size in analysis is greater than 100,000, we recommend using sparse GRM (instead of dense GRM) to adjust for sample relatedness. 
#' This function is to use \code{GCTA} software (gcta_1.93.1beta, [link](https://cnsgenomics.com/software/gcta/#Overview)) to make a \code{SparseGRMFile} to be passed to function \code{\link{GRAB.NullModel}}.
#' This function can only support \code{Linux} and \code{PLINK} files as required by \code{GCTA} software. To make a \code{SparseGRMFile}, two steps are needed. Please check \code{Details} section for more details.
#' @param PlinkFile a path to PLINK binary files (without file extension). Note that the current version (gcta_1.93.1beta) of \code{GCTA} software does not support different prefix names for BIM, BED, and FAM files. 
#' @param nPartsGRM a numeric value (e.g. 250): \code{GCTA} software can split subjects to multiple parts. For UK Biobank data analysis, it is recommended to set \code{nPartsGRM=250}. 
#' @param SparseGRMFile a path to file of output to be passed to \code{\link{GRAB.NullModel}}.
#' @param tempDir a path to store temp files from \code{\link{getTempFilesFullGRM}}. This should be consistent to the input of \code{\link{getTempFilesFullGRM}}. Default is \code{system.file("SparseGRM", "temp", package = "GRAB")}.
#' @param relatednessCutoff a cutoff for sparse GRM, only kinship coefficient greater than this cutoff will be retained in sparse GRM. *(default=0.05)*
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from PLINK files) to make sparse GRM. *(default=0.01)*
#' @param maxMissingGRM Maximal value of missing rate to select markers (from PLINK files) to make sparse GRM. *(default=0.1)*
#' @param rm.tempFiles a logical value indicating if the temp files generated in \code{\link{getTempFilesFullGRM}} will be deleted. *(default=FALSE)*
#' @details 
#' \describe{
#' \itemize{
#'   \item \code{Step 1}: Run \code{\link{getTempFilesFullGRM}} to save temporary files to \code{tempDir}.
#'   \item \code{Step 2}: Run \code{getSparseGRM} to combine the temporary files to make a \code{SparseGRMFile} to be passed to function \code{\link{GRAB.NullModel}}.
#'   }
#' }
#' Users can customize parameters including \code{(minMafGRM, maxMissingGRM, nPartsGRM)}, but functions \code{\link{getTempFilesFullGRM}} and \code{getSparseGRM} should use the same ones. 
#' Otherwise, package \code{GRAB} cannot accurately identify temporary files.
#' @examples 
#' ## Input data:
#' GenoFile = system.file("extdata", "example.bed", package = "GRAB")
#' PlinkFile = tools::file_path_sans_ext(GenoFile)   # remove file extension
#' nPartsGRM = 2;   # we recommend setting nPartsGRM = 250 for UK Biobank data analysis with 500K samples.
#' 
#' # Step 1:
#' # We strongly recommend parallel computing in high performance clusters (HPC). 
#' for(partParallel in 1:nPartsGRM){
#'   getTempFilesFullGRM(PlinkFile, nPartsGRM, partParallel)  # this function only supports Linux
#' }
#' 
#' # After step 1, the temporary files are in tempDir (default: system.file("SparseGRM", "temp", package = "GRAB")), which might needs a large amount of space.
#' 
#' # Step 2:
#' # Combine files in Step 1 to make a SparseGRMFile,
#' tempDir = system.file("SparseGRM", "temp", package = "GRAB")
#' SparseGRMFile = gsub("temp", "SparseGRM.txt", tempDir)
#' getSparseGRM(PlinkFile, nPartsGRM, SparseGRMFile)
#' 
#' @export

getSparseGRM = function(PlinkFile,
                        nPartsGRM,
                        SparseGRMFile,
                        tempDir = NULL,
                        relatednessCutoff = 0.05,
                        minMafGRM = 0.01,
                        maxMissingGRM = 0.1,
                        rm.tempFiles = FALSE)
{
  PlinkName = basename(PlinkFile)
  nDigits = floor(log10(nPartsGRM)) + 1 # 1-9 -> 1; 10-99 -> 2; 100:999 -> 3.
  
  AllIDs = c()
  n0 = 0
  SparseGRM = c()
  
  if(is.null(tempDir)) 
    tempDir = system.file("SparseGRM", "temp", package = "GRAB")
  
  ## cycle for nPartsGRM
  for(i in 1:nPartsGRM){
    cat("Analyzing part", i, "of total", nPartsGRM, "parts.\n")
    tempList = list()
    
    ##
    tempFile = paste0(tempDir, "/Plink-",PlinkName,"-autosome-minMaf-",minMafGRM,"-maxMissing-",maxMissingGRM)
    
    ## Three files generated by GCTA
    IDFile = paste0(tempFile,".part_",nPartsGRM,"_",formatC(i, width = nDigits, flag = '0'),".grm.id")
    BinFile = paste0(tempFile,".part_",nPartsGRM,"_",formatC(i, width = nDigits, flag = '0'),".grm.bin")
    NFile = paste0(tempFile,".part_",nPartsGRM,"_",formatC(i, width = nDigits, flag = '0'),".grm.N.bin")
    
    ## read in the three files
    IDs = read.table(IDFile, stringsAsFactors=F)
    ID = IDs$V2
    AllIDs = c(AllIDs, ID)
    n1 = n0 + length(ID)
    nData = (n1 - n0) * (n0 + n1 + 1) / 2
    
    grm = readBin(BinFile, n = nData, what = numeric(0), size = 4)
    # nMarkers = readBin(NFile, n = nData, what = numeric(0), size = 4)
    
    pos = which(grm > relatednessCutoff)
    value = grm[pos]
    pairs = getPairs(pos, value, AllIDs, n0, n1)
    SparseGRM = rbind(SparseGRM, pairs)
    
    n0 = n1 
  }
  
  # class(SparseGRM) = "SparseGRM"
  colnames(SparseGRM) = c("ID1", "ID2", "Value")
  write.table(SparseGRM, SparseGRMFile, row.names = F, quote=F, sep="\t", col.names = T)
  
  message = paste("The SparseGRM has been stored in", SparseGRMFile)
  
  return(message)
}


getPairs = function(pos, value, ID, n0, n1){
  start = c()
  end = c()
  temp = 0
  for(i in (n0+1):n1){
    start = c(start, temp)
    temp = temp + i;
    end = c(end, temp);
  }
  
  pos.row = sapply(pos, FUN=function(x){min(which(x<=end))})
  pos.col = pos - start[pos.row]
  
  ID1 = ID[n0 + pos.row];
  ID2 = ID[pos.col];
  pairs = cbind.data.frame(ID1, ID2, value);
  return(pairs)
}

# Suppose that subjData is only a subset of the subjects in SparseGRM
# this function is to extract subjects from SparseGRM
updateSparseGRM = function(SparseGRM, subjData)
{
  # later add another column to specify the relationship degree
  if(any(toupper(colnames(SparseGRM)) != c("ID1", "ID2", "VALUE")))
    stop("The header in 'SparseGRMFile' should be c('ID1','ID2','Value')")
  
  colnames(SparseGRM) = toupper(colnames(SparseGRM))
  
  tempGRM1 = SparseGRM;
  tempGRM2 = data.frame(ID1=tempGRM1$ID2,
                        ID2=tempGRM1$ID1,
                        VALUE=tempGRM1$VALUE)
  
  tempGRM = rbind(tempGRM1, tempGRM2)
  tempGRM = tempGRM[-1*which(duplicated(tempGRM)),]
  
  ID1 = tempGRM$ID1;
  ID2 = tempGRM$ID2;
  value = tempGRM$VALUE;
  
  if(any(!is.element(subjData, ID1)))
    stop("At least one of subjects is not in SparseGRM.")
  
  location1 = match(ID1, subjData);
  location2 = match(ID2, subjData);
  pos = which(!is.na(location1) & !is.na(location2))
  locations = rbind(location1[pos]-1,  # -1 is to convert R to C++
                    location2[pos]-1)
  
  value = value[pos];
  nSubj = length(subjData);
  KinMatListR = list(locations = locations,
                     values = value,
                     nSubj = nSubj)
  
  return(KinMatListR)
}

## set up Dense GRM and Sparse GRM
setGRM = function(GenoFile, GenoFileIndex, SparseGRMFile, subjData)
{
  # NOTE on 2022-01-27: for approaches using GRM, variance ratio is required, which needs 'GenoFile' and 'GenoFileIndex'.
  if(is.null(GenoFile))
    stop("Argument of 'GenoFile' is required to estimate variance ratio.")
  
  genoList = setGenoInput(GenoFile, GenoFileIndex, subjData)   # check Geno.R for more details
  
  if(!is.null(SparseGRMFile))
  {
    cat("Sparse GRM is used when fitting a null model.")
    SparseGRM = data.table::fread(SparseGRMFile)
    SparseGRM = as.data.frame(SparseGRM)
    
    KinMatListR = updateSparseGRM(SparseGRM, subjData)
    
    # the following function is in Main.cpp
    setSparseGRMInCPP(KinMatListR)
    optionGRM = "SparseGRM"
  }else{
    cat("Dense GRM is used when fitting a null model.")
    subjGeno = genoList$SampleIDs      # subjGeno should be the same as subjData
    if(genoList$genoType != "PLINK")
      stop("If DenseGRM is used when fitting a null model, then only PLINK format is supported.")
    
    memoryChunk = 2 # (GB)
    minMafGRM = 0.01
    maxMissingGRM = 0.1
    
    # the following function is in Main.cpp
    setDenseGRMInCPP(memoryChunk, minMafGRM, maxMissingGRM)
    optionGRM = "DenseGRM"
  }
  
  return(list(optionGRM = optionGRM,
              genoType = genoList$genoType,
              markerInfo = genoList$markerInfo))
}

setSparseGRMInStep2 = function(SparseGRMFile, objNull)
{
  SparseGRM = data.table::fread(SparseGRMFile)
  SparseGRM = as.data.frame(SparseGRM)
  KinMatListR = updateSparseGRM(SparseGRM, objNull$subjData)
  setSparseGRMInCPP(KinMatListR)    # check Main.cpp
}
