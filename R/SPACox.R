
#' SPACox method in GRAB package
#' 
#' SPACox method is an empirical approach to analyzing complex traits (including but not limited to time-to-event trait) for unrelated samples in a large-scale biobank. 
#' 
#' @details 
#' Additional list of \code{control} in \code{GRAB.NullModel()} function.
#' 
#' Additional list of \code{control} in \code{GRAB.Marker()} function.
#' 
#' @examples 
#' # Step 1: fit a null model
#' PhenoFile = system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData = data.table::fread(PhenoFile, header = T)
#' obj.SPACox = GRAB.NullModel(Surv(SurvTime, SurvEvent)~AGE+GENDER, 
#'                             data = PhenoData, 
#'                             subjData = IID, 
#'                             method = "SPACox", 
#'                             traitType = "time-to-event")
#' 
#' # Using model residuals performs exactly the same as the above. Note that confounding factors are still required in the right of the formula.
#' obj.coxph = coxph(Surv(SurvTime, SurvEvent)~AGE+GENDER, data = PhenoData, x=T)
#' obj.SPACox = GRAB.NullModel(obj.coxph$residuals~AGE+GENDER, 
#'                             data = PhenoData, 
#'                             subjData = IID, 
#'                             method = "SPACox", 
#'                             traitType = "Residual")
#' 
#' # Step 2: conduct score test
#' GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/Results_SPACox.txt")
#' GRAB.Marker(obj.SPACox, GenoFile = GenoFile, OutputFile = OutputFile, control = list(outputColumns = "zScore"))
#' data.table::fread(OutputFile)
#' @export
GRAB.SPACox = function(){
  print("Check ?GRAB.SPACox for more details about 'SPACox' method.")
}

# check the control list in null model fitting for SPACox method
checkControl.NullModel.SPACox = function(control, traitType, ...)
{
  if(!traitType %in% c("time-to-event", "Residual"))
    stop("For 'SPACox' method, only traitType of 'time-to-event' or 'Residual' is supported.")
  
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
  if(!class(response) %in% c("Surv", "Residual")) 
    stop("For SPAcox, the response variable should be of class 'Surv' or 'Residual'.")
  
  if(class(response) == "Surv")
  {
    formula = response ~ designMat
    obj.coxph = survival::coxph(formula, x=T, ...)
    
    y = obj.coxph$y
    yVec = y[,ncol(y)]
    
    mresid = obj.coxph$residuals
    # Cova = obj.coxph$x
    Cova = designMat
  }
  
  if(class(response) == "Residual")
  {
    yVec = mresid = response
    Cova = designMat
  }
  
  ### extract information from control
  range = control$range
  length.out = control$length.out
  
  ### Fit a Cox model using survival package
  # obj.coxph = survival::coxph(formula, data=data, subset=subset, x=T, na.action="na.omit", ...)
  
  
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
  
  if(length(mresid)!=length(subjData))
    stop("Please check the consistency between 'formula' and 'subjData'.")
  
  ### Get the covariate matrix to adjust for genotype
  
  
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
  
  re = list(N = length(mresid),
            mresid = mresid,
            cumul = cumul,
            tX = tX,
            yVec = yVec, 
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

setRegion.SPACox = function(objNull, control)
{
  setMarker.SPACox(objNull, control)
}


# mainMarker.SPACox = function(objNull, control, markers, genoType)
mainMarker.SPACox = function(genoType, genoIndex, outputColumns)
{
  # The following function is in Main.cpp
  OutList = mainMarkerInCPP("SPACox", genoType, genoIndex)
  obj.mainMarker = data.frame(Marker = OutList$markerVec,           # marker IDs
                              Info = OutList$infoVec,               # marker information: CHR:POS:REF:ALT
                              AltFreq = OutList$altFreqVec,         # alternative allele frequencies
                              AltCounts = OutList$altCountsVec,     # alternative allele counts
                              MissingRate = OutList$missingRateVec, # alternative allele counts
                              Pvalue = OutList$pvalVec)             # marker-level p-values
  
  # optionalColumns = c("beta", "seBeta", "zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  optionalColumns = c("zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  additionalColumns = intersect(optionalColumns, outputColumns)
  
  if(length(additionalColumns) > 0)
    obj.mainMarker = cbind.data.frame(obj.mainMarker, 
                                      as.data.frame(OutList[additionalColumns]))
  
  return(obj.mainMarker)
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

