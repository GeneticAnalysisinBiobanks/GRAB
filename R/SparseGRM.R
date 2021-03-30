
#' Please run this function before running getSparseGRM(). 
#' 
#' Please run this function before running getSparseGRM(). We strongly suggest using parallel computing for different pairs of (chrParallel, partParallel).
#' @param PlinkFile a path to Plink files. The current version (gcta_1.93.1beta) of gcta software does not support different prefix names for bim, bed and fam files. 
#' @param partParallel a numeric value (from 1 to nPartsGRM)
#' @param nPartsGRM a numeric value (e.g. 250): GCTA software can split subjects to multiple parts. For UK-Biobank analysis, it is recommanded to use 250 parts. 
#' @param tempDir a path to store temp files from getSparseGRMParallel(). This should be consistent to the input of getSparseGRM(). Default is system.file("SparseGRM", "temp", package = "GRAB").
#' @param subjData a character vector to specify subject IDs to retain (i.e. IID). Default is NULL, i.e. all subjects are retained in sparse GRM. If the number of subjects is less than 1,000, the GRM estimation might not be accurate.
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from Plink files) to construct GRM.
#' @param maxMissingGRM Maximal value of missing rate to select markers (from Plink files) to construct GRM.
#' @param threadNum Number of threads (CPUs) to use.
#' @examples 
#' ## Please check help(getSparseGRM) for example codes.
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
    stop(paste0("We generate GRM using gcta software which is only supported in Linux. Current OS is ", OS))
  
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


#' Get an object of 'SparseGRM' for GRAB_Null_Model(). Only valid in Linux since GCTA software can only be used in Linux.
#' 
#' If the sample size is greater than 100,000, we recommend using sparse GRM to adjust for sample relatedness. 
#' This function is to use GCTA software (gcta_1.93.1beta, https://cnsgenomics.com/software/gcta/#Overview) to get an object of 'SparseGRM' to be passed to function GRAB_Null_Model(). 
#' \cr\cr
#' Step 1: Run getTempFilesFullGRM(); please check help(getSparseGRMParallel) for more details. 
#' \cr
#' Step 2: Run getSparseGRM().
#' @param PlinkFile a path to Plink files. The current version (gcta_1.93.1beta) of gcta software does not support different prefix names for bim, bed and fam files. 
#' @param nPartsGRM a numeric value (e.g. 250): GCTA software can split subjects to multiple parts. For UK-Biobank analysis, it is recommanded to use 250 parts. 
#' @param SparseGRMFile a path to file of output to record SparseGRM
#' @param tempDir a path to store temp files from getSparseGRMParallel(). This should be consistent to the input of getSparseGRMParallel(). Default is system.file("SparseGRM", "temp", package = "GRAB").
#' @param relatednessCutoff a cutoff for sparse GRM, only kinship coefficient greater than this cutoff will be retained in sparse GRM.
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from Plink files) to construct GRM.
#' @param maxMissingGRM Maximal value of missing rate to select markers (from Plink files) to construct GRM.
#' @param rm.tempFiles a logical value indicating if the temp files generated in getSparseGRMParallel will be deleted
#' @examples 
#' ## Input data:
#' GenoFile = system.file("extdata", "example.bed", package = "GRAB")
#' PlinkFile = tools::file_path_sans_ext(GenoFile)   # fam/bim/bed files should have the same prefix
#' nPartsGRM = 2;   # nPartsGRM = 250 for UK Biobank data analysis
#' 
#' ## Step 1:
#' ## We strongly suggest parallel computing in high performance clusters (HPC) for different pairs of (chrParallel, partParallel). 
#' for(partParallel in 1:nPartsGRM){
#'   getTempFilesFullGRM(PlinkFile, nPartsGRM, partParallel)
#' }
#' 
#' ## After step 1, the results are stored in "tempDir (default: system.file("SparseGRM", "temp", package = "GRAB"))", which might needs large amount of space.
#' 
#' ## Step 2:
#' ## Combine results in step 1 to calculate an object with class of SparseGRM for GRAB_Null_Model(),
#' tempDir = system.file("SparseGRM", "temp", package = "GRAB")
#' SparseGRMFile = gsub("temp", "SparseGRM.txt", tempDir)
#' SparseGRM = getSparseGRM(PlinkFile, nPartsGRM, SparseGRMFile)
#' 
#' ## NOTE: You can change some options such as (minMafGRM, maxMissingGRM, nPartsGRM), but keep in mind that functions getSparseGRMParallel() and getSparseGRM() should use the same options.
#' @export

getSparseGRM = function(PlinkFile,
                        nPartsGRM,
                        SparseGRMFile,
                        tempDir = NULL,
                        relatednessCutoff = 0.05,
                        minMafGRM = 0.01,
                        maxMissingGRM = 0.1,
                        rm.tempFiles = F)
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
    print(paste("Analyzing part",i,"of total",nPartsGRM,"parts."))
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
  write.table(SparseGRM, SparseGRMFile, row.names=F, quote=F)
  
  message = paste0("The SparseGRM has been stored in ", SparseGRMFile)
  
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
updateSparseGRM = function(SparseGRM, subjData){
  names = names(SparseGRM)
  KinMatListR = list();
  
  print("names(SparseGRM) is")
  print(names)
  for(excludeChr in names){
    
    print(paste0("Updating chromosome ", excludeChr, " in 'SparseGRM' based on 'subjData'."))
    
    tempGRM1 = SparseGRM[[excludeChr]]
    # updated on 05-14-2020
    tempGRM2 = data.frame(ID1=tempGRM1$ID2,
                          ID2=tempGRM1$ID1,
                          value=tempGRM1$value)
    
    tempGRM = rbind(tempGRM1, tempGRM2)
    tempGRM = tempGRM[-1*which(duplicated(tempGRM)),]
    
    ID1 = tempGRM$ID1;
    ID2 = tempGRM$ID2;
    value = tempGRM$value;
    
    if(any(!is.element(subjData, ID1)))
      stop("At least one of subjects is not in SparseGRM.")
    
    location1 = match(ID1, subjData);
    location2 = match(ID2, subjData);
    pos = which(!is.na(location1) & !is.na(location2))
    locations = rbind(location1[pos]-1,  # -1 is to convert R to C++
                      location2[pos]-1)
    
    value = value[pos];
    nSubj = length(subjData);
    KinMatListR[[excludeChr]] = list(locations = locations,
                                     values = value,
                                     nSubj = nSubj)
  }
  
  return(KinMatListR)
}


