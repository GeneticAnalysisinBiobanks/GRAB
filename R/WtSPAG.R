
#' WtSPAG method in GRAB package
#' 
#' WtSPAG method is an empirical approach to analyzing complex traits (including but not limited to time-to-event trait) for unrelated samples in a large-scale biobank. 
#' 
#' @details 
#' Additional list of \code{control} in \code{GRAB.NullModel()} function.
#' 
#' Additional list of \code{control} in \code{GRAB.Marker()} function.
#' 
#' @import dplyr, data.table
#' @examples
#' setwd(system.file("WtSPAG", package = "GRAB"))
#' PhenoData = read.table(system.file("WtSPAG", "simuPHENO_WtSPAG.txt", package = "GRAB"), header = T)
#' RefPrevalence = 0.1
#' 
#' qcResult1 = QCforBatchEffect(GenoFile = "simuBGEN1.bgen", 
#'                              GenoFileIndex = c("simuBGEN1.bgen.bgi", 
#'                                                 "simuBGEN1.sample"),
#'                              OutputFile = "qcBGEN1.txt",
#'                              control=list(AlleleOrder = "ref-first", 
#'                                           AllMarkers = T,
#'                                           IndicatorColumn = "SurvEvent", SampleIDColumn = "IID"),
#'                              PhenoData=PhenoData,
#'                              RefAfFile = "RefMAFs.txt",
#'                              RefPrevalence = RefPrevalence,
#'                              SNPnum=1e4)
#' 
#' qcResult2 = QCforBatchEffect(GenoFile = "GenoMat2.bed",
#'                              OutputFile = "qcBGEN2.txt",
#'                              control=list(AlleleOrder = "ref-first", 
#'                                           AllMarkers = T,
#'                                           IndicatorColumn = "SurvEvent", SampleIDColumn = "IID"),
#'                              PhenoData=PhenoData,
#'                              RefAfFile = "RefMAFs.txt",
#'                              RefPrevalence = RefPrevalence,
#'                              SNPnum=1e4)
#' 
#' allQC = mergeQCresults(qcResult1, qcResult2)
#'  
#' obj.WtSPAG = GRAB.NullModel(Surv(SurvTime, SurvEvent) ~ Cov1 + Cov2,
#'                            data = PhenoData, 
#'                            subjData = IID, 
#'                            method = "WtSPAG", 
#'                            traitType = "time-to-event",
#'                            control = list(RefPrevalence = RefPrevalence))
#' 
#' names(obj.WtSPAG)
#' 
#' GRAB.Marker(obj.WtSPAG, 
#'             GenoFile = "simuBGEN1.bgen",
#'             GenoFileIndex = c("simuBGEN1.bgen.bgi", "simuBGEN1.sample"),
#'             OutputFile = "simuBGEN1.txt", 
#'             OutputFileIndex = NULL, 
#'             control = list(QCFile = "qcBGEN1.txt",
#'                            QCCutoff = qcResult1$cutoff,
#'                            AlleleOrder = "ref-first", 
#'                            AllMarkers = T))
#' 
#' @export
GRAB.WtSPAG = function(){
  print("Check ?GRAB.WtSPAG for more details about 'WtSPAG' method.")
}

# check the control list in null model fitting for SPACox method
checkControl.NullModel.WtSPAG = function(control, traitType, ...)
{
  if(!traitType %in% c("time-to-event", "binary", "Residual"))
    stop("For 'WtSPAG' method, only traitType of 'time-to-event', 'binary' or 'Residual' is supported.")
  
  default.control = list(RefPrevalence = 0)
  
  control = updateControl(control, default.control)
  
  if(control$RefPrevalence <= 0 | control$RefPrevalence >= 0.5)
    stop("control$Refprevalence is required and should be between (0, 0.5).")
  
  return(control)
}

getWeight.WtSPAG = function(Indicator, RefPrevalence)
{
  # BWJ (2023-08-08): Check NA later. 
  # if(any(!unique(Indicator) %in% c(0, 1, NA)))
  #   stop("The value of Indicator should be 0, 1, or NA")
  
  if(any(!unique(Indicator) %in% c(0, 1)))
    stop("The value of Indicator should be 0 or 1.")
  
  sumOnes = sum(Indicator)
  sumZeros = sum(1-Indicator)
  ratio = sumOnes / sumZeros
  
  weight = ifelse(Indicator == 1, 1,
                  (1-RefPrevalence) / RefPrevalence * ratio)
  return(weight)
}

# fit null model using WtSPAG method
fitNullModel.WtSPAG = function(response, designMat, subjData,
                               control=list(OutlierRatio=1.5), ...)
{
  if(!class(response) %in% c("Surv", "Residual"))   # later support binary trait and Residual
    stop("For WtSPAG, the response variable should be of class 'Surv' or 'Residual'.")
  
  if(class(response) == "Surv")
  {
    formula = response ~ designMat
    
    Indicator = as.matrix(response)[,"status"]
    RefPrevalence = control$RefPrevalence
    
    weight = getWeight.WtSPAG(Indicator, RefPrevalence)
    
    obj.coxph = survival::coxph(formula, x=T, weight = weight, robust = T, ...)
    
    y = obj.coxph$y
    yVec = y[,ncol(y)]  # status
    
    mresid = obj.coxph$residuals
    # Cova = obj.coxph$x
    Cova = designMat
  }else{
    stop("We only support 'time-to-event' trait for WtSPAG by 2023-08-08.")
  }
  
  # if(class(response) == "Residual")
  # {
  #   yVec = mresid = response
  #   Cova = designMat
  # }
  
  # var.resid = var(mresid, na.rm = T)
  ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
  q25 = quantile(mresid, 0.25, na.rm = T)
  q75 = quantile(mresid, 0.75, na.rm = T)
  IQR = q75 - q25
  
  # r.outlier = 1.5    # put this to the control argument later
  r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)  # put this to the control argument later
  cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
  posOutlier = which(mresid < cutoff[1] | mresid > cutoff[2])
  cat("The r.outlier is:",r.outlier,"\n")
  while(length(posOutlier)==0){
    r.outlier = r.outlier*0.8
    cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier = which(mresid < cutoff[1] | mresid > cutoff[2])
    cat("The curent outlier ratio is:",r.outlier,"\n")
    cat("The number of outlier is:",length(posOutlier),"\n")

  }
  
  # The original code is from SPAmix in which multiple residuals were analysis simultaneously 
  # posValue = which(!is.na(mresid))
  # posNonOutlier = setdiff(posValue, posOutlier)
  posValue = 1:length(mresid)
  posNonOutlier = setdiff(posValue, posOutlier)
  
  cat("The outlier of residuals will be passed to SPA analysis.\n")
  cat("Cutoffs to define residuals:\t", signif(cutoff,2),"\n")
  cat("Totally, ", length(posOutlier),"/", length(posValue), " are defined as outliers.\n")
  
  if(length(posOutlier) == 0)
    stop("No outlier is observed. SPA is not required in this case.")
  
  # "-1" is to convert R style (index starting from 1) to C++ style (index starting from 0)
  outLierList = list(posOutlier = posOutlier - 1,
                     posNonOutlier = posNonOutlier - 1,
                     resid = mresid,
                     resid2 = mresid^2,
                     residOutlier = mresid[posOutlier],
                     residNonOutlier = mresid[posNonOutlier],
                     resid2NonOutlier = mresid[posNonOutlier]^2)
  
  re = list(mresid = mresid, Cova = Cova, yVec = yVec, weight = weight, 
            RefPrevalence = RefPrevalence, 
            N = length(mresid),
            outLierList = outLierList)
  
  class(re) = "WtSPAG_NULL_Model"
  return(re)
}


checkControl.Marker.WtSPAG = function(control)
{
  default.control = list(SPA_Cutoff = 2);
  
  if(!"QCFile" %in% names(control))
    stop("control$QCFile is required for WtSPAG method. Please run function QCforBatchEffect for this purpose.")
  
  if(!"QCCutoff" %in% names(control))
    stop("control$QCCutoff is required for WtSPAG method. Please run function QCforBatchEffect for this purpose.")
  
  control$QCCutoff = as.numeric(control$QCCutoff)
  
  if(control$QCCutoff > 0.2)
    stop("control$QCCutoff > 0.2, the batch effect seems more than expected. The reference population might be different from the study cohort and WtSPAG is not valid in this case.")
  
  control = updateControl(control, default.control)
  
  return(control)
}

setMarker.WtSPAG = function(objNull, control)
{
  QCFile = control$QCFile
  QCCutoff = control$QCCutoff
  
  QCData = data.table::fread(QCFile)
  QCData = QCData %>% select(genoIndex, AF_ref, AN_ref, pvalue_bat)
  
  mresid = objNull$mresid
  N = objNull$N
  weight = objNull$weight
  
  SPA_Cutoff = control$SPA_Cutoff
  
  obj.setMarker = list(QCData = QCData,
                       QCCutoff = QCCutoff)
  
  # The following function is in Main.cpp
  setWtSPAGobjInCPP(mresid,
                    # weight,
                    N,
                    SPA_Cutoff,
                    objNull$outLierList)
  
  return(obj.setMarker)
}

# mainMarker.SPACox = function(objNull, control, markers, genoType)
mainMarker.WtSPAG = function(genoType, genoIndex, outputColumns, obj.setMarker)
{
  QCData = obj.setMarker$QCData
  QCCutoff = obj.setMarker$QCCutoff
  
  genoIndexSet = genoIndex
  QCData = QCData %>% filter(genoIndex %in% genoIndexSet)
  
  # QCData %>% head() %>% print()
  genoIndex = QCData$genoIndex
  AF_ref = QCData$AF_ref
  AN_ref = QCData$AN_ref
  pvalue_bat = QCData$pvalue_bat
  
  # The below function is in Main.cpp
  updateQCInCPP(AF_ref, AN_ref, pvalue_bat, QCCutoff)
  
  # The below function is in Main.cpp
  OutList = mainMarkerInCPP("WtSPAG", genoType, genoIndex)
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

#' Get model residuals using a weighted Cox regression model.
#'
#' Get model residuals using a weighted Cox regression model.
#' 
#' @param formula an refression formula: a symbolic description of the null model to be fitted
#' @param PhenoData a dataframe including at least 3 columns \code{SampleID}, \code{Indicator} and \code{time}
#' @param RefPrevalence a numeric of the event rate in the population, which is consisitent with the input in QCforBatchEffect.
#' @return an R object with a class of "QCforBatchEffect".
#' \itemize{
#'   \item{Residual}: A dataframe with 2 columns, \code{SampleID} and \code{Residuals}, and \code{Residuals} are calculated by fitting null model.
#' }
getResid.wCox = function(formula
                         ,PhenoData # an R data frame including at least 3 columns c("SampleID","Indicator","time")
                         ,RefPrevalence # consisitent with the RefPrevalence used in QCforBatchEffect
){
  suppressPackageStartupMessages(library("survival",quietly = T))
  
  # check 
  if(!"Indicator" %in% colnames(PhenoData))
    stop("Indicator is missing in PhenoData!")
  
  if(any(!unique(PhenoData$Indicator) %in% c(0,1,NA)))
    stop("The value of Indicator should be 0,1 or NA")
  
  if(!"SampleID" %in% colnames(PhenoData))
    stop("SampleID is missing in PhenoData!")
  
  if(!"time" %in% colnames(PhenoData))
    stop("time is missing in PhenoData!")
  
  #calculate weight
  PhenoData = PhenoData %>% mutate()
  #fit null model
  obj.null = coxph(formula, PhenoData, weight = weight, robust = T)
  Residual = data.frame(Resid = obj.null$residuals,SampleID = PhenoData$SampleID )
  return(Residual)
}