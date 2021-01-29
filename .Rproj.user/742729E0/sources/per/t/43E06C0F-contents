
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
#' # Simulation phenotype and genotype
#' GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
#' N = 100
#' Pheno = data.frame(ID = paste0("f",1:N,"_1"),
#'                    event=rbinom(N,1,0.5),
#'                    time=runif(N),
#'                    Cov1=rnorm(N),
#'                    Cov2=rbinom(N,1,0.5))
#' obj.SPACox = GRAB.NullModel(survival::Surv(time,event)~Cov1+Cov2, 
#'                             data=Pheno, subjData = Pheno$ID, method = "SPACox", GenoFile = GenoFile)
#' 
#' GRAB.Marker(obj.SPACox, GenoFile, chrom=1)
#' 
#' @export
#' @import data.table

GRAB.Marker = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile = NULL,
                       chrom,             
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);
  
  if(!is.null(OutputFile)){
    if(file.exists(OutputFile))
      stop(paste0("'OutputFile' of '",OutputFile,"' has existed. Please use another 'OutputFile' or remove the existing one."))
  }
    
  
  # check the setting of control, if not specified, the default setting will be used
  control = checkControl.Marker(control, NullModelClass)
  
  SampleInModel = as.character(objNull$subjIDs);
  
  ## add something in setGenoInput(.) to select specific markers requested by users
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleInModel)
  genoType = objGeno$genoType
  markers = objGeno$markers
  nMarkersEachChunk = control$nMarkersEachChunk;
  markersList = splitMarker(markers, nMarkersEachChunk);
  
  cat("Number of chunks for all markers:\t", length(markersList), "\n")
  cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
  
  setMarker(NullModelClass, objNull, control, chrom)
  
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

setMarker = function(NullModelClass, objNull, control, chrom)
{
  if(NullModelClass == "POLMM_NULL_Model")
    obj.setMarker = setMarker.POLMM(objNull, control, chrom)
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.setMarker = setMarker.SAIGE(objNull, control, chrom)
  
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