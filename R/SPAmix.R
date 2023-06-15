#' SPAmix method in GRAB package
#' 
#' SPAmix method is an empirical approach to analyzing complex traits (including but not limited to time-to-event trait) for unrelated samples in a large-scale biobank. SPAmix extend SPACox to support an admixture population or multiple populations. 
#' 
#' @details 
#' For ```SPAmix```, the confounding factors of SNP-derived PCs are required and should be specified in ```control```.
#' 
#' @examples 
#' # Step 1: fit a null model
#' PhenoFile = system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData = data.table::fread(PhenoFile, header = T)
#' N = nrow(PhenoData)
#' PhenoData = PhenoData %>% mutate(PC1 = rnorm(N), PC2 = rnorm(N))  # add two PCs
#' obj.SPAmix = GRAB.NullModel(Surv(SurvTime, SurvEvent)~AGE+GENDER+PC1+PC2, 
#'                             data = PhenoData, 
#'                             subjData = IID, 
#'                             method = "SPAmix", 
#'                             traitType = "time-to-event",
#'                             control = list(PC_columns = "PC1,PC2"))
#' 
#' # Using model residuals performs exactly the same as the above. Note that confounding factors are still required in the right of the formula.
#' obj.coxph = coxph(Surv(SurvTime, SurvEvent)~AGE+GENDER+PC1+PC2, data = PhenoData, x=T)
#' obj.SPAmix = GRAB.NullModel(obj.coxph$residuals~AGE+GENDER+PC1+PC2, 
#'                             data = PhenoData, 
#'                             subjData = IID, 
#'                             method = "SPAmix", 
#'                             traitType = "Residual",
#'                             control = list(PC_columns = "PC1,PC2"))
#' 
#' # Step 2: conduct score test
#' GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/Results_SPAmix.txt")
#' GRAB.Marker(obj.SPAmix, GenoFile = GenoFile, OutputFile = OutputFile, control = list(outputColumns = "zScore"))
#' data.table::fread(OutputFile)
#' @export
GRAB.SPAmix = function(){
  print("Check ?GRAB.SPAmix for more details about 'SPAmix' method.")
}

################### This file includes the following functions

# ------------ used in 'GRAB_Marker.R' -----------
# 1. checkControl.Marker.SPAmix(control)
# 2. setMarker.SPAmix(objNull, control)
# 3. mainMarker.SPAmix()

# check the control list in marker-level testing
# unified control list (such as nMarkersEachChunk) can be found in checkControl.Marker() in contro.R
checkControl.Marker.SPAmix = function(control)
{
  default.control = list(SPA_Cutoff = 2,
                         dosage_option = "rounding_first")  # "rounding_first" or "rounding_last"
                    # list(SPA_Cutoff = 2,
                    #      outputColumns = c("beta", "seBeta"));
  
  control = updateControl(control, default.control)  # This file is in 'Util.R'
  
  # check the parameter
  if(!control$dosage_option %in% c("rounding_first", "rounding_last"))
    stop("control$dosage_option should be 'rounding_first' or 'rounding_last'.")
  
  return(control)
}

# SPAmixClass::SPAmixClass(arma::vec t_resid,
#                          arma::mat t_XinvXX,
#                          arma::mat t_tX,
#                          arma::mat t_PCs,
#                          int t_N,
#                          double t_SPA_Cutoff)

setMarker.SPAmix = function(objNull, control)
{
  # the below function is in 'Main.cpp'
  setSPAmixobjInCPP(objNull$resid,
                    # objNull$X.invXX,
                    # objNull$tX,
                    objNull$PCs,
                    objNull$N,
                    control$SPA_Cutoff,
                    objNull$posOutlier - 1,
                    setdiff(1:objNull$N, objNull$posOutlier) - 1)  # "-1" is to convert R style (index starting from 1) to C++ style (index starting from 0)
  
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,"."))
}

# mainMarker.SPAmix(genoType, genoIndex, outputColumns)
mainMarker.SPAmix = function(genoType, genoIndex, outputColumns)
{
  OutList = mainMarkerInCPP("SPAmix", genoType, genoIndex);  
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


# fit null model using SPAmix method
fitNullModel.SPAmix = function(response, designMat, subjData, control, ...)
{
  if(!class(response) %in% c("Surv", "Residual")) 
    stop("For SPAmix, the response variable should be of class 'Surv' or 'Residual'.")
  
  if(class(response) == "Surv")
  {
    formula = response ~ designMat
    
    obj.coxph = coxph(formula, x=T, ...)
    
    ### Check input arguments
    # p2g = check_input(pIDs, gIDs, obj.coxph, range)
    
    y = obj.coxph$y
    yVec = y[,ncol(y)]
    
    mresid = obj.coxph$residuals
    Cova = obj.coxph$x
  }
  
  if(class(response) == "Residual")
  {
    yVec = mresid = response
    Cova = designMat
  }
  
  if(length(mresid)!=length(subjData))
    stop("Please check the consistency between 'formula' and 'subjData'.")
  
  PC_columns = control$PC_columns
  
  # Remove the below if checked later
  cat("colnames(designMat):\n")
  print(colnames(designMat))
  cat("PC columns specified in 'control':\n")
  print(PC_columns)
  cat("dimension of 'designMat' and 'Cova':\n")
  print(dim(designMat))
  print(dim(Cova))
  
  if(any(!PC_columns %in% colnames(designMat)))
    stop("PC columns specified in 'control$PC_columns' should be in 'formula'.")
  
  pos_col = match(PC_columns, colnames(designMat))
  
  PCs = Cova[,pos_col,drop=F]  
  # X = cbind(1, PCs)
  # X.invXX = X %*% solve(t(X)%*%X)
  # tX = t(X)
  
  var.resid = var(mresid)
  
  ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
  q25 = quantile(mresid, 0.25)
  q75 = quantile(mresid, 0.75)
  IQR = q75 - q25
  r.outlier = 1.5    # put this to the control argument later
  cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
  posOutlier = which(mresid < cutoff[1] | mresid > cutoff[2])
  
  cat("The outlier of residuals will be passed to SPA analysis.\n")
  cat("Cutoffs to define residuals:\t", signif(cutoff,2),"\n")
  cat("Totally, ", length(posOutlier),"/", length(mresid), " are defined as outliers.\n")
    
  objNull = list(resid = mresid,
                 var.resid = var.resid,
                 # tX = tX,
                 # X.invXX = X.invXX,
                 N = nrow(Cova),
                 yVec = yVec,          # event variable: 0 or 1
                 PCs = PCs,
                 posOutlier = posOutlier)
  
  class(objNull) = "SPAmix_NULL_Model"
  return(objNull)
}

# check the control list in null model fitting for SPACox method
checkControl.NullModel.SPAmix = function(control, traitType, ...)
{
  if(!traitType %in% c("time-to-event", "Residual"))
    stop("For 'SPAmix' method, only traitType of 'time-to-event' or 'Residual' is supported.")

  if(is.null(control$PC_columns))
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') is required for 'SPAmix' method.")
  
  if(length(control$PC_columns) != 1)
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') should be a character, not a character vector.")
  
  control$PC_columns = unlist(strsplit(control$PC_columns, split=","))
  if(length(control$PC_columns) == 1)
    warning("We detected that only one PC column exsits, is that what you want? Note that control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') should be a character splitted using ','.")
  
  # default.control = list(range = c(-100, 100),
  #                        length.out = 10000)
  
  # control = updateControl(control, default.control)
  
  # check the parameters
  return(control)
}
