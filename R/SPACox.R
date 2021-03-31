
#' SPACox method in GRAB package
#' 
#' SPACox method is to analysis time-to-event phenotype for unrelated samples in a large-scale biobank. 
#' 
#' @details 
#' Please check \code{?GRAB.control} for the generic list of \code{control} in \code{GRAB.NullModel()} and \code{GRAB.Marker()}.
#' 
#' Additional list of \code{control} in \code{GRAB.NullModel()} function
#' \itemize{
#' \item{\code{range}: a two-element numeric vector [default=c(-100,100)] to specify the domain of the empirical CGF.}
#' \item{\code{length.out}: a positive integer [default=10000] for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.}
#' }
#' Additional list of \code{control} in \code{GRAB.Marker()} function
#' \itemize{
#' \item{\code{pVal_covaAdj_Cutoff}: a numeric value [default=5e-5]. If the p-value is less than this cutoff, then we would use an additional technic to adjust for covariates.}
#' \item{\code{SPA_cutoff}: a numeric value [default=2] to specify the standard deviation cutoff to be used. If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation, otherwise, its p value is calculated based on a saddlepoint approximation.}
#' }
#' @examples 
#' # Step 1: fit a null model
#' PhenoData = read.table(system.file("extdata", "example.pheno", package = "GRAB"), header = T)
#' obj.SPACox = GRAB.NullModel(survival::Surv(time,event) ~ Cova1 + Cova2, 
#'                             data = PhenoData, subjData = PhenoData$IID, 
#'                             method = "SPACox", traitType = "time-to-event")
#' 
#' # Step 2: perform score test
#' GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/SPACoxMarkers.txt")
#' GRAB.Marker(obj.SPACox, GenoFile = GenoFile,
#'             OutputFile = OutputFile)
#' head(read.table(OutputFile, header=T))
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
# fitNullModel.SPACox = function(formula, data, subset, subjData, subjGeno, control, ...)
fitNullModel.SPACox = function(response, designMat, subjData, control, ...)
{
  if(class(response) != "Surv") stop("For SPACox, the response variable should be of class 'Surv'.")
  
  formula = response ~ designMat
  
  ### extract information from control
  range = control$range
  length.out = control$length.out
  
  ### Fit a Cox model using survival package
  # obj.coxph = survival::coxph(formula, data=data, subset=subset, x=T, na.action="na.omit", ...)
  obj.coxph = survival::coxph(formula, x=T, ...)
  
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
            subjData = subjData) 
            # subjGeno = subjGeno)
  
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
  
  # The following function is in Main.cpp
  setSPACoxobjInCPP(cumul,
                    mresid,
                    XinvXX,
                    tX,
                    N,
                    pVal_covaAdj_Cutoff,
                    SPA_Cutoff)
}


# mainMarker.SPACox = function(objNull, control, markers, genoType)
mainMarker.SPACox = function(genoType, genoIndex)
{
  # The following function is in Main.cpp
  OutList = mainMarkerInCPP("SPACox",
                            genoType,
                            genoIndex)
                            # control$missing_cutoff,
                            # control$min_maf_marker,
                            # control$min_mac_marker)  
  
  # Rcpp::List OutList = Rcpp::List::create(Rcpp::Named("markerVec") = markerVec,
  #                                         Rcpp::Named("infoVec") = infoVec,
  #                                         Rcpp::Named("altFreqVec") = altFreqVec,
  #                                         Rcpp::Named("missingRateVec") = missingRateVec,
  #                                         Rcpp::Named("BetaVec") = BetaVec,
  #                                         Rcpp::Named("seBetaVec") = seBetaVec,
  #                                         Rcpp::Named("pvalVec") = pvalVec);
  
  markerVec = OutList$markerVec   # marker IDs
  infoVec = OutList$infoVec       # marker infomation: CHR:POS:REF:ALT
  altFreqVec = OutList$altFreqVec       # 
  altCountsVec = OutList$altCountsVec       # 
  missingRateVec = OutList$missingRateVec       # minor allele frequencies (freq of ALT if flip=F, freq of REF if flip=T)
  # BetaVec = OutList$BetaVec     # beta for ALT if flip=F, beta for REF if flip=T
  # seBetaVec = OutList$seBetaVec # sebeta
  pvalVec = OutList$pvalVec;      # marker-level p-values
  zScoreVec = OutList$zScoreVec;
  
  obj.mainMarker = data.frame(Marker = markerVec,
                              Info = infoVec,
                              AltFreq = altFreqVec,
                              AltCounts = altCountsVec,
                              MissingRate = missingRateVec,
                              Pval = pvalVec,
                              zScore = zScoreVec)
  
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

