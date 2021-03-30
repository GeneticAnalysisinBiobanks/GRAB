
#' GRAB: Marker-level score test
#' 
#' GRAB package uses score test for GWAS: in step 1, we fit a null model (check \code{?GRAB.NullModel}) including response variable, covariates, and GRM (if needed). In step 2, we perform score test for marker-level analysis (check \code{?GRAB.Marker}) and region-level analysis (check \code{?GRAB.Region}).
#' 
#' @param objNull output object of the GRAB.NullModel() function. 
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". If Null (default), the same prefix as GenoFile is used. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param OutputFile a character of output file to store the analysis results. 
#' @param OutputFileIndex a character of output index file to record the end point in case that program ends unexpectedly
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
#' # If you want to use "SPACox" method, please check ?GRAB.SPACox for more details.
#' # If you want to use "POLMM" method, please check ?GRAB.POLMM for more details.
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
  outIndex = checkOutputFile(OutputFile, OutputFileIndex, "Marker", format(nMarkersEachChunk, scientific=F))    # this function is in 'Util.R'
  
  subjData = as.character(objNull$subjData);
  
  ## set up an object for genotype
  objGeno = setGenoInput(GenoFile, GenoFileIndex, subjData, control)  # this function is in 'Geno.R'
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  CHROM = markerInfo$CHROM
  genoIndex = markerInfo$genoIndex
  
  # all markers were split into multiple chunks, 
  genoIndexList = splitMarker(genoIndex, nMarkersEachChunk, CHROM);
  nChunks = length(genoIndexList)
  
  cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
  cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
  cat("Number of chunks for all markers:\t", nChunks, "\n")
  if(outIndex != 1)
    cat("Restart the analysis from chunk:\t", outIndex, "\n")
  
  chrom = "InitialChunk"
  for(i in outIndex:nChunks)
  {
    tempList = genoIndexList[[i]]
    genoIndex = tempList$genoIndex
    tempChrom = tempList$chrom
    
    # set up objects that do not change for different variants
    if(tempChrom != chrom){
      setMarker(NullModelClass, objNull, control, chrom)
      chrom = tempChrom
    }
    
    print(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- "))
    
    # main function to calculate summary statistics for markers in one chunk
    resMarker = mainMarker(NullModelClass, genoType, genoIndex)
    
    # write summary statistics to output file
    if(i == 1){
      data.table::fwrite(resMarker, OutputFile, quote = F, sep = "\t", append = F, col.names = T, na="NA")
      write.table(matrix(c("GRAB.outIndex", "Please_do_not_modify_this_file.", "Marker", 
                           format(nMarkersEachChunk, scientific=F), 1), 
                         ncol = 1), 
                  OutputFileIndex, col.names = F, row.names = F, quote = F, append = F)
    }else{
      data.table::fwrite(resMarker, OutputFile, quote = F, sep = "\t", append = T, col.names = F, na="NA")
      write.table(matrix(i, ncol = 1), 
                  OutputFileIndex, col.names = F, row.names = F, quote = F, append = T)
    }
  }
  
  # information to users
  output = paste0("Done! The results have been saved to '", OutputFile,"'.")
  write.table(matrix(-1, ncol = 1), 
              OutputFileIndex, col.names = F, row.names = F, quote = F, append = T)
  
  return(output)
}

setMarker = function(NullModelClass, objNull, control, chrom)
{
  # The following function is in Main.cpp
  setMarker_GlobalVarsInCPP(control$impute_method,
                            control$missing_cutoff,
                            control$min_maf_marker,
                            control$min_mac_marker,
                            control$omp_num_threads)
  
  # The following function is in POLMM.R
  if(NullModelClass == "POLMM_NULL_Model")
    obj.setMarker = setMarker.POLMM(objNull, control, chrom)
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.setMarker = setMarker.SAIGE(objNull, control)
  
  # The following function is in SPACox.R
  if(NullModelClass == "SPACox_NULL_Model")
    obj.setMarker = setMarker.SPACox(objNull, control)
    
  return(obj.setMarker)
}

mainMarker = function(NullModelClass, genoType, genoIndex)
{
  # The following function is in 'POLMM.R'
  if(NullModelClass == "POLMM_NULL_Model")
    obj.mainMarker = mainMarker.POLMM(genoType, genoIndex)
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.mainMarker = mainMarker.SAIGE(genoType, genoIndex)
  
  # The following function is in SPACox.R
  if(NullModelClass == "SPACox_NULL_Model")
    obj.mainMarker = mainMarker.SPACox(genoType, genoIndex)
  
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