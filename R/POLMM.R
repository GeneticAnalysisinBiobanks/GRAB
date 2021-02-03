
#' How to run POLMM method in GRAB package
#' 
#' POLMM method is to analysis ordinal categorical phenotype for large-scale biobank data
#' 
#' @details 
#' List of 'control' in GRAB.NullModel function
#' \itemize{
#' \item{range: a two-element numeric vector [default=c(-100,100)] to specify the domain of the empirical CGF.}
#' \item{length.out: a positive integer [default=10000] for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.}
#' }
#' List of 'control' in GRAB.Marker function
#' \itemize{
#' \item{impute_method: a character string [default="fixed"] to specify the method to impute missing genotypes. "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).}
#' \item{missing_cutoff: a numeric value [default=0.15] to specify the cutoff of the missing rates. Any variant with missing rate higher than this cutoff will be excluded from the analysis.}
#' \item{min_maf_marker: a numeric value [default=0.001] to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.}
#' \item{min_mac_marker: a numeric value [default=20] to specify the cutoff of the minimal MAC. Any SNP with MAC < cutoff will be excluded from the analysis.}
#' \item{nMarkersEachChunk: a numeric value [default=10000] to specify the number of markers in each chunk.}
#' \item{pVal_covaAdj_Cutoff: a numeric value [default=5e-5]. If the p-value is less than this cutoff, then we would use an additional technic to adjust for covariates.}
#' \item{SPA_cutoff: a numeric value [default=2] to specify the standard deviation cutoff to be used. If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation, otherwise, its p value is calculated based on a saddlepoint approximation.".}
#' }
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
#' res.SPACox = GRAB.Marker(obj.SPACox, GenoFile)
#' head(res.SPACox)
#' @export
GRAB.POLMM = function(){
  print("Check ?GRAB.POLMM for more details about 'POLMM' method.")
}

# check the control list in null model fitting for POLMM method
checkControl.NullModel.POLMM = function(control)
{
  default.control = list(range = c(-100, 100),
                         length.out = 10000)
  
  control = updateControl(control, default.control)
  
  # check the parameters
  range = control$range
  length.out = control$length.out
  
  if(range[1] >= -50 | range[2] <= 50 | length.out <= 1000)
    stop("We suggest setting argument 'control$range=c(-100,100)' and 'control$length.out=10000'.")
  
  if(range[2]!=-1*range[1])
    stop("range[2] should be -1*range[1]")
  
  return(control)
}


# fit null model using POLMM method
fitNullModel.POLMM = function(formula, data, subjData, subjGeno, control, ...)
{
  ### extract information from control
  range = control$range
  length.out = control$length.out
  print(range)
  print(length.out)
}

setMarker.POLMM = function(objNull, control, chrom)
{
  if(objNull$controlList$LOCO){
    if(!chrom %in% names(objNull$LOCOList))
      stop("'chrom' should be in names(objNull$LOCOList).")
    objCHR = objNull$LOCOList[[chrom]]
  }else{
    # to be continued
  }
  
  # single marker analysis does not require sparse GRM any more 
  # Note: it might be not so accurate if min_mac_marker is very low
  SPmatR.CHR = list(locations = c(0,0), values = 1)
  
  setPOLMMobjInCPP(objCHR$muMat,
                   objCHR$iRMat,
                   objNull$Cova,
                   objNull$yVec,          # 1 to J
                   SPmatR.CHR,
                   objNull$tau,
                   control$printPCGInfo,
                   control$tolPCG,
                   control$maxiterPCG)
  
  print(paste0("The current POLMM.control$nMarkers_output is ", nMarkers_output,"."))
}

mainMarker.POLMM = function(objNull, control, markers, genoType)
{
  OutList = mainMarkerInCPP("POLMM",
                            genoType,
                            markers,
                            control$SPA_cutoff,
                            control$missing_cutoff,
                            control$min_maf_region,
                            control$min_mac_region)  
  
  markerVec = OutList$markerVec   # marker IDs
  infoVec = OutList$infoVec       # marker infomation: CHR:POS:REF:ALT
  flipVec = OutList$flipVec       # 
  freqVec = OutList$freqVec       # minor allele frequencies (freq of ALT if flip=F, freq of REF if flip=T)
  BetaVec = OutList$BetaVec       # beta for ALT if flip=F, beta for REF if flip=T
  seBetaVec = OutList$seBetaVec   # sebeta
  PvalVec = OutList$pvalVec;      # marker-level p-values
  
  obj.mainMarker = data.frame(Marker = markerVec,
                              Info = infoVec,
                              Flip = flipVec,
                              Freq = freqVec,
                              Beta = BetaVec,
                              seBeta = seBetaVec,
                              Pval = PvalVec)
  return(obj.mainMarker)
}