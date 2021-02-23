
#' GRAB: Marker-level genetic analysis
#' 
#' Marker-level genetic analysis: test for association between a phenotype of interest and genome-wide genetic markers
#' 
#' @param objNull output object of the GRAB.NullModel() function. 
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". If Null (default), the same prefix as GenoFile is used. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param OutputFile a character of output file to store the analysis results. If Null (default), the results are outputted to R. 
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
                       OutputFile = NULL, 
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);
  
  if(!is.null(OutputFile)){
    if(file.exists(OutputFile))
      stop(paste0("'OutputFile' of '", OutputFile, "' has existed. Please use another 'OutputFile' or remove the existing one."))
  }
    
  # check the setting of control, if not specified, the default setting will be used
  control = checkControl.Marker(control, NullModelClass)
  
  subjData = as.character(objNull$subjData);
  
  ## add something in setGenoInput(.) to select specific markers requested by users
  objGeno = setGenoInput(GenoFile, GenoFileIndex, subjData)
  genoType = objGeno$genoType
  markers = objGeno$markers
  nMarkersEachChunk = control$nMarkersEachChunk;
  markersList = splitMarker(markers, nMarkersEachChunk);
  
  cat("Number of chunks for all markers:\t", length(markersList), "\n")
  cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
  
  setMarker(NullModelClass, objNull, control)
  
  nChunks = length(markersList)
  output = c()
  for(i in 1:nChunks)
  {
    print(paste0("Analyzing Chunk ", i, "/", nChunks, " ......"))
    markers = markersList[[i]]
    resMarker = mainMarker(NullModelClass, objNull, control, markers, genoType)
    
    quote = F
    if(i == 1){
      append = F; col.names = T
    }else{
      append = T; col.names = F
    }
      
    if(is.null(OutputFile)){
      output = rbind(output, resMarker)
    }else{
      data.table::fwrite(resMarker, OutputFile, quote = quote, sep = "\t", append = append, col.names = col.names)
    }
  }
  
  if(!is.null(OutputFile)){
    output = paste0("The analysis results have been saved to '", OutputFile,"'.")
  }
  
  return(output)
}

setMarker = function(NullModelClass, objNull, control)
{
  if(NullModelClass == "POLMM_NULL_Model")
    obj.setMarker = setMarker.POLMM(objNull, control)
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.setMarker = setMarker.SAIGE(objNull, control)
  
  if(NullModelClass == "SPACox_NULL_Model")
    obj.setMarker = setMarker.SPACox(objNull, control)
    
  return(obj.setMarker)
}

mainMarker = function(NullModelClass, objNull, control, markers, genoType)
{
  if(NullModelClass == "POLMM_NULL_Model")
    obj.mainMarker = mainMarker.POLMM(objNull, control, markers, genoType)
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.mainMarker = mainMarker.SAIGE(objNull, control, markers, genoType)
  
  if(NullModelClass == "SPACox_NULL_Model")
    obj.mainMarker = mainMarker.SPACox(objNull, control, markers, genoType)
  
  return(obj.mainMarker)
}

checkControl.Marker = function(control, NullModelClass)
{
  # check if control is an R list
  if(!is.null(control))
    if(class(control) != "list")
      stop("If specified, the argument of 'control' should be an R 'list'.")
  
  default.marker.control = list(impute_method = "fixed",  
                                missing_cutoff = 0.15,
                                min_maf_marker = 0.001,
                                min_mac_marker = 20,
                                nMarkersEachChunk = 10000)
  
  control = updateControl(control, default.marker.control)
  
  ######### SPACox method
  if(NullModelClass == "SPACox_NULL_Model")
    control = checkControl.Marker.SPACox(control)
  
  
  ######### the below is for other methods
  
  print("The below is the list of control parameters used in marker-level genetic association analysis.")
  print(control)
  
  return(control)
}