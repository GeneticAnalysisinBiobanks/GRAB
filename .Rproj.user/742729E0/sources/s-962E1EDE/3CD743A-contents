
#' GRAB: Single-marker analysis
#' 
#' Single-marker analysis: Test for association between genetic variants and a phenotype of interest
#' 
#' @param objNull an output object of the GRAB.NullModel() function with a class of "POLMM_NULL_Model". 
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". The default is NULL, that is, to share the same prefix as GenoFile. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param OutputFile a character of output file to store the analysis results
#' @param chrom a character to specify chromosome of the markers in analysis. Must be specified unless LOCO = F when fitting the null model.
#' @param POLMM.control a list of parameters for controlling the POLMM.Marker(). The default is NULL, that is, to use the default parameters. Check 'Details' and 'Examples' for more details.
#' @details 
#' More information about the list of 'control'
#' \itemize{
#' \item{impute_method: imputation method when genotype data is missing [default="fixed"].}
#' \item{missing_cutoff: exclude markers with missing rate greater than this cutoff [default=0.15].}
#' \item{min_maf_marker: exclude markers with MAF less than this cutoff [default=0.001].}
#' \item{min_mac_marker: exclude markers with MAC less than this cutoff [default=20].}
#' \item{nMarkers_output: output results after analyzing each chunk of fixed number of markers [default=10000].}
#' \item{SPA_cutoff: a cutoff to determine "normal distribution approximation" or "saddlepoint approximation" [default=2].}
#' }
#' @return The analysis results will be written to OutputFile with the following columns
#' \item{Marker}{Marker IDs extracted from "GenoFile" or "GenoFileIndex".}
#' \item{Info}{Marker Infomation of "CHR:POS:REF:ALT". This information is from "GenoFile" or "GenoFileIndex" and does not change even if the REF/ALT alleles are flipped in analysis}
#' \item{Freq}{Minor allele frequency (always < 0.5) in analysis.}
#' \item{Flip}{a logical value indicating if the REF/ALT alleles were switched in analysis. This information is useful to estimate the effect direction.}
#' \item{Stat}{Score statistics (changed to beta later.) The direction depends on both "Info" and "Flip".}
#' \item{Var}{Estimated variance of the score statistics (changed to se.beta later)}
#' \item{Pvalue}{p-value from normal distribution approximation or saddlepoint approximation.}
#' @examples
#' famFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.fam", package = "POLMM")
#' GenoFile = gsub("-ext.fam", "-ext.bed", famFile)
#' AnnoFile = system.file("extdata", "AnnoFile.txt", package = "POLMM")
#' OutputFile = gsub("AnnoFile","OutputFile",AnnoFile)
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.RData", package = "POLMM")
#' load(SparseGRMFile)
#' objNullFile = system.file("objNull.RData", package = "POLMM")
#' load(objNullFile)
#' chrom = 1
#' 
#' OUTPUT = POLMM.Marker(objNull, AnnoFile, GenoFile, GenoFileIndex = NULL, OutputFile, 
#'                       SparseGRM, chrom, POLMM.control = list(max_maf_region = 0.5))
#'      
#' @export
#' @import data.table

GRAB.Marker = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile,
                       chrom,
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);
  checkOutputFile(OutputFile);
  
  # check the setting of control, if not specified, the default setting will be used
  control = checkControl(control)
  print("The below is the list of control parameters used in analysis.")
  print(control)
  
  SampleInModel = as.character(objNull$subjIDs);
  
  ## add something in setGenoInput(.) to select specific markers requested by users
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleInModel)  
  markers = objGeno$markers
  nMarkersEachChunk = control$nMarkersEachChunk;
  markersList = splitMarker(markers, nMarkersEachChunk);
  print(paste("Markers are splitted into", length(markersList), "chunks, each of which includes <=", nMarkersEachChunk, "markers."))
  
  setMarker(NullModelClass, objNull, control)
  
  nChunks = length(markersList)
  for(i in 1:nChunks)
  {
    print(paste0("Analyzing Chunk ", i, "/", nChunks, " ......"))
    markers = markersList[[i]]
    resMarker = mainMarker(markers, NullModelClass, objNull, control)
    
    quote = F
    if(i == 1){
      append = F; col.names = T
    }else{
      append = T; col.names = F
    }
      
    data.table::fwrite(resMarker, OutputFile, quote = quote, sep = "\t", append = append, col.names = col.names)
  }
  
  message = paste0("The analysis results have been saved to '", OutputFile,"'.")
  return(message)
}

setMarker = function(NullModelClass, objNull, control)
{
  if(NullModelClass == "POLMM")
    obj.setMarker = setMarker.POLMM(objNull, control)
  
  if(NullModelClass == "SAIGE")
    obj.setMarker = setMarker.SAIGE(objNull, control)
  
  if(NullModelClass == "SPACox")
    obj.setMarker = setMarker.SPACox(objNull, control)
    
  return(obj.setMarker)
}

mainMarker = function(NullModelClass, objNull, control)
{
  if(NullModelClass == "POLMM")
    obj.mainMarker = mainMarker.POLMM(objNull, control)
  
  if(NullModelClass == "SAIGE")
    obj.mainMarker = mainMarker.SAIGE(objNull, control)
  
  if(NullModelClass == "SPACox")
    obj.mainMarker = mainMarker.SPACox(objNull, control)
  
  return(obj.mainMarker)
}

