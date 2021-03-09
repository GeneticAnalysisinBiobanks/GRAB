
#' GRAB: Marker-level genetic analysis
#' 
#' Marker-level genetic analysis: test for association between a phenotype of interest and genome-wide genetic markers
#' 
#' @param objNull output object of the GRAB.NullModel() function. 
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". If Null (default), the same prefix as GenoFile is used. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param OutputFile a character of output file to store the analysis results. 
#' @param control a list of parameters for controlling the GRAB.Marker(). 
#' @details 
#' List of 'control' is different for different methods
#' \itemize{
#' \item{SAIGE: to analyze binary phenotype, check ?GRAB.SAIGE for more details.}
#' \item{POLMM: to analyze ordinal categorical phenotype, check ?GRAB.POLMM for more details.}
#' \item{SPACox: to analyze time-to-event phenotype, check ?GRAB.SPACox for more details.}
#' \item{SPAGE: to analyze gene-environment interaction effect for binary phenotype, check ?GRAB.SPAGE for more details.}
#' }
#' @return The results will be written in a file (if OutputFile != NULL) or be saved to an R data.frame (default, if OutputFile == NULL). The results include the following columns.
#' \item{Marker}{Marker IDs extracted from "GenoFile" and "GenoFileIndex".}
#' \item{Info}{Marker Infomation of "CHR:POS:REF:ALT". This information is from "GenoFile" or "GenoFileIndex" and does not change even if the REF/ALT alleles are flipped in analysis.}
#' \item{Freq}{Minor allele frequency (always < 0.5) in analysis.}
#' \item{Flip}{Logical value indicating if the REF/ALT alleles were switched in analysis. This information is useful to estimate the effect direction.}
#' \item{Beta}{Estimated effect size. The sign (positive or negative) depends on "Info" and "Flip". For example, if 'Flip' is false and 'Beta' is positive, then ALT allele is to increase continuous trait (trait.type == "quantitative"), to increase the risk of being cases (trait.type == "binary"), or to xxxx}
#' \item{seBeta}{Estimated standard error (se) of the effect size}
#' \item{Pval}{p-value from normal distribution approximation or saddlepoint approximation.}
#' @examples
#' # We put examples to the specific help pages for different methods. 
#' # For example, if you want to use "SPACox" method, please check ?GRAB.SPACox for more details. 
#' @export
#' @import data.table

GRAB.Marker = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile,
                       OutputFileIndex = NULL,   ## check it later: record the end point, avoid letting the program starts from the very beginning
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);         # this function is in 'Util.R'
  
  if(is.null(OutputFileIndex)) OutputFileIndex = paste0(OutputFile, ".index")
  outIndex = checkOutputFile(OutputFile, OutputFileIndex)    # this function is in 'Util.R'
  
  # check the setting of control, if not specified, the default setting will be used
  control = checkControl.Marker(control, NullModelClass)
  
  SampleIDs = as.character(objNull$SampleIDs);
  
  ## set up an object for genotype
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs)  # this function is in 'Geno.R'
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  
  # different genotype format corresponds to different index
  if(genoType == "PLINK")
    genoIndex = markerInfo$PositionInPLINK
  if(genoType == "BGEN")
    genoIndex = markerInfo$StartPositionInBGEN
  
  # all markers were split into multiple chunks, 
  nMarkersEachChunk = control$nMarkersEachChunk;
  genoIndexList = splitMarker(genoIndex, nMarkersEachChunk);
  nChunks = length(genoIndexList)
  
  cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
  cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
  cat("Number of chunks for all markers:\t", nChunks, "\n")
  
  # set up objects that do not change for different variants
  setMarker(NullModelClass, objNull, control)
  
  for(i in outIndex:nChunks)
  {
    print(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, " ---- "))
    genoIndex = genoIndexList[[i]]
    
    # main function to calculate summary statistics for markers in one chunk
    resMarker = mainMarker(NullModelClass, genoType, genoIndex)
    
    # write summary statistics to output file
    if(i == 1){
      data.table::fwrite(resMarker, OutputFile, quote = F, sep = "\t", append = F, col.names = T)
      write.table(matrix(c("GRAB.outIndex", "Please do not modify this file.", 1), ncol = 1), 
                  OutputFileIndex, col.names = F, row.names = F, quote = F, append = F)
    }else{
      data.table::fwrite(resMarker, OutputFile, quote = F, sep = "\t", append = T, col.names = F)
      write.table(matrix(i, ncol = 1), 
                  OutputFileIndex, col.names = F, row.names = F, quote = F, append = T)
    }
  }
  
  # information to users
  output = paste0("Done! The results have been saved to '", OutputFile,"'.")
  
  return(output)
}

checkControl.Marker = function(NullModelClass, control)
{
  # check if control is an R list
  if(!is.null(control))
    if(class(control) != "list")
      stop("If specified, the argument of 'control' should be an R 'list'.")
  
  # uniform default control setting for marker-level analysis
  default.marker.control = list(impute_method = "fixed",  
                                missing_cutoff = 0.15,
                                min_maf_marker = 0.001,
                                min_mac_marker = 20,
                                nMarkersEachChunk = 10000)
  
  control = updateControl(control, default.marker.control)
  
  # specific default control setting for different approaches
  if(NullModelClass == "POLMM_NULL_Model")
    control = checkControl.Marker.POLMM(control)    # This function is in 'POLMM.R'
  
  if(NullModelClass == "SPACox_NULL_Model")
    control = checkControl.Marker.SPACox(control)
  
  
  # print control list 
  print("The below is the list of control parameters used in marker-level genetic association analysis.")
  print(control)
  
  return(control)
}

setMarker = function(NullModelClass, objNull, control)
{
  
  setMarker_GlobalVarsInCPP(control$impute_method,
                            control$missing_cutoff,
                            control$min_maf_marker,
                            control$min_mac_marker)
  
  if(NullModelClass == "POLMM_NULL_Model")
    obj.setMarker = setMarker.POLMM(objNull, control)
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.setMarker = setMarker.SAIGE(objNull, control)
  
  if(NullModelClass == "SPACox_NULL_Model")
    obj.setMarker = setMarker.SPACox(objNull, control)
    
  return(obj.setMarker)
}

mainMarker = function(NullModelClass, genoType, genoIndex)
{
  if(NullModelClass == "POLMM_NULL_Model")
    obj.mainMarker = mainMarker.POLMM(genoType, genoIndex)
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.mainMarker = mainMarker.SAIGE(genoType, genoIndex)
  
  if(NullModelClass == "SPACox_NULL_Model")
    obj.mainMarker = mainMarker.SPACox(genoType, genoIndex)
  
  return(obj.mainMarker)
}


## split 'markerInfo' into multiple chunks, each of which includes no more than 'nMarkersEachChunk' markers
splitMarker = function(genoIndex, nMarkersEachChunk)
{
  M = length(genoIndex)
  
  idxStart = seq(1, M, nMarkersEachChunk)
  idxEnd = idxStart + nMarkersEachChunk - 1
  
  nChunks = length(idxStart)
  idxEnd[nChunks] = M
  
  genoIndexList = list()
  for(i in 1:nChunks){
    idxMarker = idxStart[i]:idxEnd[i]
    genoIndexList[[i]] = genoIndex[idxMarker]
  }
  
  return(genoIndexList)
}