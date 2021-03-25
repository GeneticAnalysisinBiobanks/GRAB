
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
#' GenoFile = system.file("extdata", "example.bed", package = "GRAB")
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

################### This file includes the following functions

# ------------ used in 'GRAB_Marker.R' -----------
# 1. checkControl.Marker.POLMM(control)
# 2. setMarker.POLMM(objNull, control)
# 3. mainMarker.POLMM()

# check the control list in marker-level testing
checkControl.Marker.POLMM = function(control)
{
  default.control = list();
  control = updateControl(control, default.control)  # This file is in 'Util.R'
  
  # check the parameter
  
  return(control)
}

# check the control list in region-level testing
checkControl.Region.POLMM = function(control)
{
  default.control = list();
  control = updateControl(control, default.control)  # This file is in 'Util.R'
  
  # check the parameter
  
  return(control)
}

# set up an object in C++
setMarker.POLMM = function(objNull, control)
{
  muMat = objNull$muMat;
  iRMat = objNull$iRMat;
  Cova = objNull$Cova;
  yVec = objNull$yVec;
  SPmatR = objNull$SPmatR;
  tau = objNull$tau;
  printPCGInfo = control$printPCGInfo;
  tolPCG = control$tolPCG;
  maxiterPCG = control$maxiterPCG;
  varRatio = objNull$varRatio; 
  StdStat_cutoff = control$StaStat_cutoff;
  
  setPOLMMobjInCPP(muMat,
                   iRMat,
                   Cova,
                   yVec,
                   SPmatR,
                   tau,
                   printPCGInfo,
                   tolPCG,
                   maxiterPCG,
                   varRatio,
                   StdStat_cutoff)
}

# main function to calculae summary statistics
mainMarker.POLMM = function(genoType, genoIndex)
{
  OutList = mainMarkerInCPP("POLMM",
                            genoType,
                            genoIndex);  
  
  markerVec = OutList$markerVec   # marker IDs
  infoVec = OutList$infoVec       # marker infomation: CHR:POS:REF:ALT
  altFreqVec = OutList$altFreqVec       # minor allele frequencies (freq of ALT if flip=F, freq of REF if flip=T)
  BetaVec = OutList$BetaVec       # beta for ALT if flip=F, beta for REF if flip=T
  seBetaVec = OutList$seBetaVec   # sebeta
  pvalVec = OutList$pvalVec;      # marker-level p-values
  
  obj.mainMarker = data.frame(Marker = markerVec,
                              Info = infoVec,
                              AltFreq = altFreqVec,
                              Beta = BetaVec,
                              seBeta = seBetaVec,
                              Pval = pvalVec)
  return(obj.mainMarker)
}

# check the control list in null model fitting for POLMM method
checkControl.NullModel.POLMM = function(control)
{
  # default setting of control for POLMM method
  default.control = list(memoryChunk = 2,
                         seed = 12345678,
                         tracenrun = 30,
                         maxiter = 100,
                         tolBeta = 0.001,
                         tolTau = 0.002,
                         tau = 0.2,
                         maxiterPCG = 100,
                         tolPCG = 1e-6,
                         maxiterEps = 100,
                         tolEps = 1e-10,
                         minMafVarRatio = 0.1,
                         maxMissingVarRatio = 0.1, 
                         nSNPsVarRatio = 20,
                         CVcutoff = 0.0025,
                         LOCO = T,
                         numThreads = "auto",
                         stackSize = "auto",
                         grainSize = 1,
                         minMafGRM = 0.01,
                         maxMissingGRM = 0.1,
                         showInfo = T,
                         onlyCheckTime = F)
  
  control = updateControl(control, default.control)
  
  # check the parameters
  # range = control$range
  # length.out = control$length.out
  # 
  # if(range[1] >= -50 | range[2] <= 50 | length.out <= 1000)
  #   stop("We suggest setting argument 'control$range=c(-100,100)' and 'control$length.out=10000'.")
  # 
  # if(range[2]!=-1*range[1])
  #   stop("range[2] should be -1*range[1]")
  
  return(control)
}


# fit null model using POLMM method
fitNullModel.POLMM = function(response, designMat, subjData, control)
{
  ######## -------------- fit the null POLMM --------------  ###########
  
  bVec = rep(0, n)  # initiate random effect of 0
  
  objNull = fitPOLMMcpp(t_flagSparseGRM = flagSparseGRM,       # if 1, then use SparseGRM, otherwise, use DenseGRM
                        t_flagGMatRatio = flagGMatRatio,       # if 1, then use GMatRatio, otherwise, extract from Plink files
                        t_bimfile = bimFile,
                        t_famfile = famFile,
                        t_bedfile = bedFile,
                        t_posSampleInPlink = posSampleInPlink,
                        t_Cova = Cova,
                        t_yVec = yVec,                         # should be from 1 to J
                        t_beta = beta,
                        t_bVec = bVec,
                        t_eps = eps,
                        t_tau = control$tau,
                        t_GMatRatio = GMat,                    # only used if m_LOCO = FALSE
                        t_SparseGRM = SparseGRM,
                        t_controlList = control)
  return(objNull)
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


