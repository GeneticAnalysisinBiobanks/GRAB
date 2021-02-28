
#' How to run SPACox method in GRAB package
#' 
#' SPACox method is to analysis time-to-event phenotype for large-scale biobank data
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
GRAB.SPACox = function(){
  print("Check ?GRAB.SPACox for more details about 'SPACox' method.")
}

# check the control list in null model fitting for SPACox method
checkControl.NullModel.SPACox = function(control)
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

# fit null model using SPACox method
fitNullModel.SPACox = function(formula, data, subset, subjData, subjGeno, control, ...)
{
  ### extract information from control
  range = control$range
  length.out = control$length.out
  
  ### Fit a Cox model using survival package
  obj.coxph = survival::coxph(formula, data=data, subset=subset, x=T, na.action="na.omit", ...)
  
  ### The below is commented since it has been checked in fitNullModel() 
  ### By Wenjian Bi on 01-31-2021
  
  # if(!is.null(obj.coxph$na.action)){
  #   posNA = c(obj.coxph$na.action)
  #   if(any(posNA > length(subjData)))
  #     stop("Please check the consistency between 'formula' and 'subjData'.")
  #   
  #   print(paste0("Due to missing data in response/indicators, ",length(posNA)," entries are removed from analysis."))
  #   print("If concerned about the power loss, users can impute data first and then use SPACox package.")
  #   
  #   subjData = subjData[-1*posNA]  # remove IDs with missing data
  # }
  
  mresid = obj.coxph$residuals
  
  if(length(mresid)!=length(subjData))
    stop("Please check the consistency between 'formula' and 'subjData'.")
  
  ### Get the covariate matrix to adjust for genotype
  Cova = obj.coxph$x
  
  X = cbind(1, Cova)
  X.invXX = X %*% solve(t(X)%*%X)
  tX = t(X)
  
  ### calculate empirical CGF for martingale residuals
  idx0 = qcauchy(1:length.out/(length.out+1))
  idx1 = idx0 * max(range) / max(idx0)
  
  cumul = NULL
  print("Start calculating empirical CGF for martingale residuals...")
  c = 0
  for(i in idx1){
    c = c+1
    t = i
    e_resid = exp(mresid*t)
    M0 = mean(e_resid)
    M1 = mean(mresid*e_resid)
    M2 = mean(mresid^2*e_resid)
    K0 = log(M0)
    K1 = M1/M0
    K2 = (M0*M2-M1^2)/M0^2
    cumul = rbind(cumul, c(t, K0, K1, K2))
    if(c %% 1000 == 0) print(paste0("Complete ",c,"/",length.out,"."))
  }
  
  re = list(mresid = mresid,
            cumul = cumul,
            tX = tX,
            X.invXX = X.invXX,
            subjData = subjData, 
            subjGeno = subjGeno)
  
  class(re) = "SPACox_NULL_Model"
  return(re)
}


checkControl.Marker.SPACox = function(control)
{
  default.control = list(pVal_covaAdj_Cutoff = 5e-05,
                         SPA_Cutoff = 2);
  
  control = updateControl(control, default.control)
  
  return(control)
}

setMarker.SPACox = function(objNull, control)
{
  cumul = objNull$cumul
  mresid = objNull$mresid
  XinvXX = objNull$X.invXX
  tX = objNull$tX
  N = length(mresid)
  pVal_covaAdj_Cutoff = control$pVal_covaAdj_Cutoff
  SPA_Cutoff = control$SPA_Cutoff
  
  setSPACoxobjInCPP(cumul,
                    mresid,
                    XinvXX,
                    tX,
                    N,
                    pVal_covaAdj_Cutoff,
                    SPA_Cutoff)
}


mainMarker.SPACox = function(objNull, control, markers, genoType)
{
  OutList = mainMarkerInCPP("SPACox",
                            genoType,
                            markers,
                            control$missing_cutoff,
                            control$min_maf_marker,
                            control$min_mac_marker)  
  
  markerVec = OutList$markerVec   # marker IDs
  infoVec = OutList$infoVec       # marker infomation: CHR:POS:REF:ALT
  flipVec = OutList$flipVec       # 
  freqVec = OutList$freqVec       # minor allele frequencies (freq of ALT if flip=F, freq of REF if flip=T)
  # BetaVec = OutList$BetaVec     # beta for ALT if flip=F, beta for REF if flip=T
  # seBetaVec = OutList$seBetaVec # sebeta
  PvalVec = OutList$pvalVec;      # marker-level p-values
  
  obj.mainMarker = data.frame(Marker = markerVec,
                              Info = infoVec,
                              Flip = flipVec,
                              Freq = freqVec,
                              # Beta = BetaVec,
                              # seBeta = seBetaVec,
                              Pval = PvalVec)
  
}


check_input1 = function(obj.null, Geno.mtx, par.list)
{
  if(class(obj.null)!="SPACox_NULL_Model")
    stop("obj.null should be a returned outcome from SPACox_Null_Model()")

  if(any(obj.null$gIDs != rownames(Geno.mtx))) stop("gIDs should be the same as rownames(Geno.mtx).")
  if(is.null(rownames(Geno.mtx))) stop("Row names of 'Geno.mtx' should be given.")
  if(is.null(colnames(Geno.mtx))) stop("Column names of 'Geno.mtx' should be given.")
  if(!is.numeric(Geno.mtx)|!is.matrix(Geno.mtx)) stop("Input 'Geno.mtx' should be a numeric matrix.")

  if(!is.numeric(par.list$min.maf)|par.list$min.maf<0|par.list$min.maf>0.5) stop("Argument 'min.maf' should be a numeric value >= 0 and <= 0.5.")
  if(!is.numeric(par.list$Cutoff)|par.list$Cutoff<0) stop("Argument 'Cutoff' should be a numeric value >= 0.")
  # if(!is.element(par.list$impute.method,c("none","bestguess","random","fixed"))) stop("Argument 'impute.method' should be 'none', 'bestguess', 'random' or 'fixed'.")
  if(!is.element(par.list$impute.method,c("fixed"))) stop("Argument 'impute.method' should be 'fixed'.")
  if(!is.numeric(par.list$missing.cutoff)|par.list$missing.cutoff<0|par.list$missing.cutoff>1) stop("Argument 'missing.cutoff' should be a numeric value between 0 and 1.")
  if(!is.element(par.list$G.model,c("Add","Dom","Rec"))) stop("Argument 'G.model' should be 'Add', 'Dom' or 'Rec'.")
}

