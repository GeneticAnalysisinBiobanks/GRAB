
#' Conduct marker-level genetic association testing
#' 
#' Test for association between phenotype of interest and genetic marker.
#' 
#' @param objNull the output object of function \code{\link{GRAB.NullModel}}. 
#' @param GenoFile a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param GenoFileIndex additional index files corresponding to the \code{GenoFile}. If \code{NULL} (default), the prefix is the same as GenoFile. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param OutputFile a character of output file to save the analysis results. 
#' @param OutputFileIndex a character of output index file to record the end point. If the program ends unexpectedly, the end point can help \code{GRAB} package understand where to restart the analysis. If \code{NULL} (default), \code{OutputFileIndex = paste0(OutputFile, ".index")}. 
#' @param control a list of parameters for controlling function \code{GRAB.Marker}, more details can be seen in \code{Details} section. 
#' @details 
#' \code{GRAB} package supports \code{SAIGE}, \code{POLMM}, and \code{SPACox} methods. 
#' Detailed information about the analysis methods is given in the \code{Details} section of \code{\link{GRAB.NullModel}}. 
#' Users do not need to specify them since functions \code{GRAB.Marker} and \code{\link{GRAB.Region}} will check the \code{class(objNull)}.
#' 
#' ## The following details are about argument \code{control}
#' The below is to let users customize markers to include in analysis. 
#' If these parameters are not specified, \code{GRAB} package will include all markers in analysis. 
#' For PLINK files, the default \code{control$AlleleOrder = "alt-first"}; 
#' for BGEN files, the default \code{control$AlleleOrder = "ref-first"}.
#'   \itemize{
#'   \item \code{IDsToIncludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{IDsToExcludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{RangesToIncludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{RangesToExcludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{AlleleOrder}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   }
#' The below is to customize the quality-control (QC) process.
#'   \itemize{
#'   \item \code{omp_num_threads}: (To be added later) a numeric value (default: value from data.table::getDTthreads()) to specify the number of threads in OpenMP for parallel computation.
#'   \item \code{ImputeMethod}: a character, "mean" (default), "bestguess", or "drop" (to be added later). Please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{MissingRateCutoff}: a numeric value *(default=0.15)*. Markers with missing rate > this value will be excluded from analysis.  
#'   \item \code{MinMAFCutoff}: a numeric value *(default=0.001)*. Markers with MAF < this value will be excluded from analysis.  
#'   \item \code{MinMACCutoff}: a numeric value *(default=20)*. Markers with MAC < this value will be excluded from analysis.  
#'   \item \code{nMarkersEachChunk}: number of markers *(default=10000)* in one chunk to output.
#'   }
#'  The below is to customize the columns in the \code{OutputFile}. 
#'  Columns of \code{Marker}, \code{Info}, \code{AltFreq}, \code{AltCounts}, \code{MissingRate}, \code{Pvalue} are included for all methods.
#'  \itemize{
#'  \item \code{outputColumns}: For example, for POLMM method, users can set \code{control$outputColumns = c("beta", "seBeta", "AltFreqInGroup")}. 
#'     \itemize{
#'     \item \code{POLMM}: Default: \code{beta}, \code{seBeta}; Optional: \code{zScore}, \code{AltFreqInGroup}, \code{nSamplesInGroup}, \code{AltCountsInGroup}
#'     \item \code{SPACox}: Optional: \code{zScore}
#'     }
#'  }
#' @return The analysis results are written in a file of \code{OutputFile}, which includes the following columns.
#' \item{Marker}{Marker IDs extracted from \code{GenoFile} and \code{GenoFileIndex}.}
#' \item{Info}{Marker Information of "CHR:POS:REF:ALT". The order of REF/ALT depends on \code{control$AlleleOrder}: "ref-first" or "alt-first".}
#' \item{AltFreq}{Alternative allele frequency (before genotype imputation, might be > 0.5). If the \code{AltFreq} of most markers are > 0.5, you should consider resetting \code{control$AlleleOrder}.}
#' \item{AltCounts}{Alternative allele counts (before genotype imputation).}
#' \item{MissingRate}{Missing rate for each marker}
#' \item{Pvalue}{Association test p-value}
#' The following columns can be customized using \code{control$outputColumns}. Check \code{\link{makeGroup}} for details about phenotype grouping which are used for
#' \code{nSamplesInGroup}, \code{AltCountsInGroup}, and \code{AltFreqInGroup}.
#' \item{beta}{Estimated effect size of the ALT allele.}
#' \item{seBeta}{Estimated standard error (se) of the effect size.}
#' \item{zScore}{z score, standardized score statistics, usually follows a standard normal distribution.}
#' \item{nSamplesInGroup}{Number of samples in different phenotype groups. This can be slightly different from the original distribution due to the genotype missing.}
#' \item{AltCountsInGroup}{Alternative allele counts (before genotype imputation) in different phenotype groups.}
#' \item{AltFreqInGroup}{Alternative allele frequency (before genotype imputation) in different phenotype groups.}
#' @examples
#' objNullFile = system.file("results", "objNull.RData", package = "GRAB")
#' load(objNullFile)
#' class(objNull)    # "POLMM_NULL_Model", that indicates an object from POLMM method.
#' 
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/simuOUTPUT.txt")
#' GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' 
#' ## make sure the output files does not exist at first
#' file.remove(OutputFile)
#' file.remove(paste0(OutputFile, ".index"))
#' 
#' GRAB.Marker(objNull,
#'             GenoFile = GenoFile,
#'             OutputFile = OutputFile)
#'             
#' data.table::fread(OutputFile)
#' 
#' ## additional columns of "zScore", "nSamplesInGroup", "AltCountsInGroup", "AltFreqInGroup"
#' ## We do not recommend adding too many columns for all markers
#' 
#' file.remove(OutputFile)
#' file.remove(paste0(OutputFile, ".index"))
#' GRAB.Marker(objNull,
#'             GenoFile = GenoFile,
#'             OutputFile = OutputFile,
#'             control = list(outputColumns = c("beta", "seBeta", "zScore","nSamplesInGroup","AltCountsInGroup","AltFreqInGroup")))
#' data.table::fread(OutputFile)
#'             
#' @export
#' @import data.table
GRAB.Marker = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile,
                       OutputFileIndex = NULL,   
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);         # this function is in 'Util.R'
  
  if(is.null(OutputFileIndex)) 
    OutputFileIndex = paste0(OutputFile, ".index")
  
  # check the setting of control, if not specified, the default setting will be used
  # The following functions are in 'control.R'
  checkControl.ReadGeno(control)
  control = checkControl.Marker(control, NullModelClass)  
  nMarkersEachChunk = control$nMarkersEachChunk;
  outList = checkOutputFile(OutputFile, OutputFileIndex, "Marker", format(nMarkersEachChunk, scientific=F))    # this function is in 'Util.R'
  
  # added by XH-2023-05-09
  if(NullModelClass == "SPAGRM_NULL_Model")
  {
    if(length(objNull$MAF_interval) > 1)
    {
      if(control$min_maf_marker <= min(objNull$MAF_interval))
        stop("min_maf_marker is out of MAF_interval, Please reset min_maf_marker or check MAF_interval.")
    }
  }
  
  if(NullModelClass == "SAGELD_NULL_Model")
  {
    if(length(objNull$MAF_interval) > 1)
    {
      if(control$min_maf_marker <= min(objNull$MAF_interval))
        stop("min_maf_marker is out of MAF_interval, Please reset min_maf_marker or check MAF_interval.")
    }
  }
  
  
  
  # added by yuzhuoma
  if(NullModelClass == "SPAyuzhuoma_NULL_Model")
  {
    if(length(objNull$MAF_interval) > 1)
    {
      if(control$min_maf_marker <= min(objNull$MAF_interval))
        stop("min_maf_marker is out of MAF_interval, Please reset min_maf_marker or check MAF_interval.")
    }
  }
  
  
  
  
  
  indexChunk = outList$indexChunk
  Start = outList$Start
  End = outList$End
  
  if(End)
  {
    message = paste0("The analysis has been completed in earlier analysis. Results have been saved in '", OutputFile, "'. ",
                     "If you want to change parameters and restart the analysis, please use another 'OutputFile'.")
    return(message)
  }
  
  if(!Start){
    message = paste0("We detected that parts of analysis have been conducted from file:\t",
                     OutputFileIndex,"\n",
                     "We restart the analysis from chunk:\t",indexChunk+1,"\n");
    cat(message)
  }
  
  subjData = as.character(objNull$subjData);
  
  Group = makeGroup(objNull$yVec)  # this function is in Util.R: categorize subjects into multiple groups and check AltFreq/AltCoutns in each group.
  ifOutGroup = any(c("AltFreqInGroup", "AltCountsInGroup") %in% control$outputColumns)
  
  ## set up an object for genotype
  objGeno = setGenoInput(GenoFile, GenoFileIndex, subjData, control)  # this function is in 'Geno.R'
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  CHROM = markerInfo$CHROM
  genoIndex = markerInfo$genoIndex
  
  # all markers were split into multiple chunks, 
  # 1. SNPs in the same CHROM will be grouped into the chunk
  # 2. the chunks will be ordered based on CHROM
  genoIndexList = splitMarker(genoIndex, nMarkersEachChunk, CHROM);  
  
  nChunks = length(genoIndexList)
  
  cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
  cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
  cat("Number of chunks for all markers:\t", nChunks, "\n")
  
  chrom = "InitialChunk"
  for(i in (indexChunk+1):nChunks)
  {
    tempList = genoIndexList[[i]]
    genoIndex = tempList$genoIndex
    tempChrom = tempList$chrom
    
    # set up objects that do not change for different variants
    # print(c(tempChrom, chrom))
    
    if(tempChrom != chrom){
      # print("test1")
      obj.setMarker = setMarker(NullModelClass, objNull, control, chrom, Group, ifOutGroup)
      # print("test2")
      chrom = tempChrom
    }
    
    cat(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- \n"))
    
    # main function to calculate summary statistics for markers in one chunk
    resMarker = mainMarker(NullModelClass, genoType, genoIndex, control$outputColumns, objNull, obj.setMarker)
    
    writeOutputFile(Output = list(resMarker), 
                    OutputFile = list(OutputFile), 
                    OutputFileIndex = OutputFileIndex,
                    AnalysisType = "Marker",
                    nEachChunk = format(nMarkersEachChunk, scientific=F),
                    indexChunk = i,
                    Start = (i==1),
                    End = (i==nChunks))
  }
  
  # information to users
  output = paste0("Analysis done! The results have been saved to '", OutputFile,"'.")

  return(output)
}

setMarker = function(NullModelClass, objNull, control, chrom, Group, ifOutGroup)
{
  # Check Main.cpp
  nGroup = length(unique(Group))
  setMarker_GlobalVarsInCPP(control$impute_method,
                            control$missing_cutoff,
                            control$min_maf_marker,
                            control$min_mac_marker,
                            control$omp_num_threads,
                            Group, ifOutGroup, nGroup)
  
  # Check POLMM.R
  if(NullModelClass == "POLMM_NULL_Model")
    obj.setMarker = setMarker.POLMM(objNull, control, chrom)
  
  # Check SAIGE.R
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.setMarker = setMarker.SAIGE(objNull, control)
  
  # Check SPACox.R
  if(NullModelClass == "SPACox_NULL_Model")
    obj.setMarker = setMarker.SPACox(objNull, control)
  
  # Check SPAmix.R
  if(NullModelClass == "SPAmix_NULL_Model")
    obj.setMarker = setMarker.SPAmix(objNull, control)
  
  # Check SPAGRM.R
  if(NullModelClass == "SPAGRM_NULL_Model")
    obj.setMarker = setMarker.SPAGRM(objNull, control)
  
  # Check SPAyuzhuoma.R
  if(NullModelClass == "SPAyuzhuoma_NULL_Model")
    obj.setMarker = setMarker.SPAyuzhuoma(objNull, control)
  
  # Check SAGELD.R
  if(NullModelClass == "SAGELD_NULL_Model")
    obj.setMarker = setMarker.SAGELD(objNull, control)
  
  # Check WtSPAG.R
  if(NullModelClass == "WtSPAG_NULL_Model")
    obj.setMarker = setMarker.WtSPAG(objNull, control)
    
  return(obj.setMarker)
}

mainMarker = function(NullModelClass, genoType, genoIndex, outputColumns, objNull, obj.setMarker)
{
  # Check 'POLMM.R'
  if(NullModelClass == "POLMM_NULL_Model"){
    
    # Check 'Main.cpp'
    OutList = mainMarkerInCPP("POLMM", genoType, genoIndex);  
    
    obj.mainMarker = data.frame(Marker = OutList$markerVec,           # marker IDs
                                Info = OutList$infoVec,               # marker information: CHR:POS:REF:ALT
                                AltFreq = OutList$altFreqVec,         # alternative allele frequencies
                                AltCounts = OutList$altCountsVec,     # alternative allele counts
                                MissingRate = OutList$missingRateVec, # alternative allele counts
                                Pvalue = OutList$pvalVec)             # marker-level p-values
    
    optionalColumns = c("beta", "seBeta", "zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
    additionalColumns = intersect(optionalColumns, outputColumns)
    
    if(length(additionalColumns) > 0)
      obj.mainMarker = cbind.data.frame(obj.mainMarker, 
                                        as.data.frame(OutList[additionalColumns]))
  }

  # Check 'SAIGE.R'
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.mainMarker = mainMarker.SAIGE(genoType, genoIndex, outputColumns)
   
  # Check 'SPACox.R'
  if(NullModelClass == "SPACox_NULL_Model")
    obj.mainMarker = mainMarker.SPACox(genoType, genoIndex, outputColumns)
  
  # Check 'SPAmix.R'
  if(NullModelClass == "SPAmix_NULL_Model")
    obj.mainMarker = mainMarker.SPAmix(genoType, genoIndex, outputColumns, objNull)
  
  # Check 'SPAGRM.R'
  if(NullModelClass == "SPAGRM_NULL_Model")
    obj.mainMarker = mainMarker.SPAGRM(genoType, genoIndex, outputColumns)
  
  # Check 'SPAyuzhuoma.R'
  if(NullModelClass == "SPAyuzhuoma_NULL_Model")
    obj.mainMarker = mainMarker.SPAyuzhuoma(genoType, genoIndex, outputColumns)
  
  # Check 'SAGELD.R'
  if(NullModelClass == "SAGELD_NULL_Model")
    obj.mainMarker = mainMarker.SAGELD(genoType, genoIndex, outputColumns, objNull)
  
  # Check 'WtSPAG.R'
  if(NullModelClass == "WtSPAG_NULL_Model")
    obj.mainMarker = mainMarker.WtSPAG(genoType, genoIndex, outputColumns, obj.setMarker)
  
  return(obj.mainMarker)
}


## split 'markerInfo' into multiple chunks, each of which includes no more than 'nMarkersEachChunk' markers
splitMarker = function(genoIndex, nMarkersEachChunk, CHROM)
{
  genoIndexList = list()
  iTot = 1;
  
  uCHROM = unique(CHROM)
  for(chrom in uCHROM){
    pos = which(CHROM == chrom)
    gIdx = genoIndex[pos]
    M = length(gIdx)
    
    idxStart = seq(1, M, nMarkersEachChunk)
    idxEnd = idxStart + nMarkersEachChunk - 1
    
    nChunks = length(idxStart)
    idxEnd[nChunks] = M
    
    for(i in 1:nChunks){
      idxMarker = idxStart[i]:idxEnd[i]
      genoIndexList[[iTot]] = list(chrom = chrom,
                                   genoIndex = gIdx[idxMarker])
      iTot = iTot + 1;
    }
  }
  
  return(genoIndexList)
}