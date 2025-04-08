#' SPAmixPlusV4 method in GRAB package
#' 
#' SPAmixPlusV4 method is an empirical approach to analyzing complex traits (including but not limited to time-to-event trait) for unrelated samples in a large-scale biobank. SPAmixPlusV4 extend SPACox to support an admixture population or multiple populations. 
#' 
#' @details 
#' For ```SPAmixPlusV4```, the confounding factors of SNP-derived PCs are required and should be specified in ```control```.
#' 
#' @examples 
#' # Step 1: fit a null model
#' PhenoFile = system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData = data.table::fread(PhenoFile, header = T)
#' N = nrow(PhenoData)
#' PhenoData = PhenoData %>% mutate(PC1 = rnorm(N), PC2 = rnorm(N))  # add two PCs, which are required for SPAmixPlusV4
#' 
#' # Users can directly specify a time-to-event trait to analyze
#' obj.SPAmixPlusV4 = GRAB.NullModel(Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2, 
#'                             data = PhenoData, 
#'                             subjData = IID, 
#'                             method = "SPAmixPlusV4", 
#'                             traitType = "time-to-event",
#'                             control = list(PC_columns = "PC1,PC2"))
#' 
#' # Using model residuals performs exactly the same as the above. Note that confounding factors are still required in the right of the formula.
#' obj.coxph = coxph(Surv(SurvTime, SurvEvent)~AGE+GENDER+PC1+PC2, data = PhenoData)
#' obj.SPAmixPlusV4 = GRAB.NullModel(obj.coxph$residuals ~ AGE + GENDER + PC1 + PC2, 
#'                             data = PhenoData, 
#'                             subjData = IID, 
#'                             method = "SPAmixPlusV4", 
#'                             traitType = "Residual",
#'                             control = list(PC_columns = "PC1,PC2"))
#'                             
#' # SPAmixPlusV4 also supports multiple residuals as below                           
#' obj.coxph = coxph(Surv(SurvTime, SurvEvent)~AGE+GENDER+PC1+PC2, data = PhenoData)
#' obj.lm = lm(QuantPheno ~ AGE+GENDER+PC1+PC2, data = PhenoData)
#' obj.SPAmixPlusV4 = GRAB.NullModel(obj.coxph$residuals + obj.lm$residuals ~ AGE + GENDER + PC1 + PC2, 
#'                             data = PhenoData, 
#'                             subjData = IID, 
#'                             method = "SPAmixPlusV4", 
#'                             traitType = "Residual",
#'                             control = list(PC_columns = "PC1,PC2"))
#'                              
#' # Step 2: conduct score test
#' GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/Results_SPAmixPlusV4.txt")
#' GRAB.Marker(obj.SPAmixPlusV4, GenoFile = GenoFile, OutputFile = OutputFile, control = list(outputColumns = "zScore"))
#' data.table::fread(OutputFile)
#' @export
GRAB.SPAmixPlusV4 = function(){
  print("Check ?GRAB.SPAmixPlusV4 for more details about 'SPAmixPlusV4' method.")
}

################### This file includes the following functions

# ------------ used in 'GRAB_Marker.R' -----------
# 1. checkControl.Marker.SPAmixPlusV4(control)
# 2. setMarker.SPAmixPlusV4(objNull, control)
# 3. mainMarker.SPAmixPlusV4()

# check the control list in marker-level testing
# unified control list (such as nMarkersEachChunk) can be found in checkControl.Marker() in control.R
checkControl.Marker.SPAmixPlusV4 = function(control)
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

# SPAmixPlusV4Class::SPAmixPlusV4Class(arma::vec t_resid,
#                          arma::mat t_XinvXX,
#                          arma::mat t_tX,
#                          arma::mat t_PCs,
#                          int t_N,
#                          double t_SPA_Cutoff)

setMarker.SPAmixPlusV4 = function(objNull, control)
{
  # the below function is in 'Main.cpp'
  setSPAmixPlusV4objInCPP(objNull$resid,
                          # objNull$X.invXX,
                          # objNull$tX,
                          objNull$PCs,
                          objNull$N,
                          control$SPA_Cutoff,
                          objNull$outLierList)  
  
  # outLierList[[i]] = list(posValue = posValue - 1,
  #                         posOutlier = posOutlier - 1,
  #                         posNonOutlier = posNonOutlier - 1)
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,"."))
}

# mainMarker.SPAmixPlusV4(genoType, genoIndex, outputColumns)
mainMarker.SPAmixPlusV4 = function(genoType, genoIndex, outputColumns, objNull)
{
  OutList = mainMarkerInCPP("SPAmixPlusV4", genoType, genoIndex);
  
  nPheno = objNull$nPheno;
  obj.mainMarker = data.frame(Pheno = paste0("pheno_", 1:nPheno),
                              Marker = rep(OutList$markerVec, each = nPheno),           # marker IDs
                              Info = rep(OutList$infoVec, each = nPheno),               # marker information: CHR:POS:REF:ALT
                              AltFreq = rep(OutList$altFreqVec, each = nPheno),         # alternative allele frequencies
                              AltCounts = rep(OutList$altCountsVec, each = nPheno),     # alternative allele counts
                              MissingRate = rep(OutList$missingRateVec, each = nPheno), # alternative allele counts
                              Pvalue = OutList$pvalVec)             # marker-level p-values
  
  # optionalColumns = c("beta", "seBeta", "zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  optionalColumns = c("zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  additionalColumns = intersect(optionalColumns, outputColumns)
  
  if(length(additionalColumns) > 0)
    obj.mainMarker = cbind.data.frame(obj.mainMarker, 
                                      as.data.frame(OutList[additionalColumns]))
  
  return(obj.mainMarker)
}


# # fit null model using SPAmixPlusV4 method
# fitNullModel.SPAmixPlusV4 = function(response, designMat, subjData, control=list(OutlierRatio=1.5), ...)
# {
#   if(!class(response) %in% c("Surv", "Residual"))
#     stop("For SPAmixPlusV4, the response variable should be of class 'Surv' or 'Residual'.")
# 
#   if(class(response) == "Surv")
#   {
#     formula = response ~ designMat
# 
#     obj.coxph = coxph(formula, x=T, ...)
# 
#     ### Check input arguments
#     # p2g = check_input(pIDs, gIDs, obj.coxph, range)
# 
#     y = obj.coxph$y
#     yVec = y[,ncol(y)]
# 
#     mresid = obj.coxph$residuals
#     Cova = obj.coxph$x
# 
#     if(length(mresid) != length(subjData))
#       stop("Please check the consistency between 'formula' and 'subjData'.")
# 
#     mresid = matrix(mresid, ncol=1)
#   }
# 
#   if(class(response) == "Residual")
#   {
#     yVec = mresid = response
#     Cova = designMat
# 
#     print(head(mresid))
#     if(nrow(mresid) != length(subjData))
#       stop("Please check the consistency between 'formula' and 'subjData'.")
#   }
# 
#   PC_columns = control$PC_columns
# 
#   # Remove the below if checked later
#   cat("colnames(designMat):\n")
#   print(colnames(designMat))
#   cat("PC columns specified in 'control':\n")
#   print(PC_columns)
#   cat("dimension of 'designMat' and 'Cova':\n")
#   print(dim(designMat))
#   print(dim(Cova))
# 
#   if(any(!PC_columns %in% colnames(designMat)))
#     stop("PC columns specified in 'control$PC_columns' should be in 'formula'.")
# 
#   pos_col = match(PC_columns, colnames(designMat))
# 
#   PCs = Cova[,pos_col,drop=F]
#   # X = cbind(1, PCs)
#   # X.invXX = X %*% solve(t(X)%*%X)
#   # tX = t(X)
# 
#   outLierList = list()
#   nPheno = ncol(mresid)
#   for(i in 1:nPheno)
#   {
#     mresid.temp = mresid[,i]
# 
#     # var.resid = var(mresid.temp, na.rm = T)
#     ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
#     q25 = quantile(mresid.temp, 0.25, na.rm = T)
#     q75 = quantile(mresid.temp, 0.75, na.rm = T)
#     IQR = q75 - q25
#     r.outlier =   r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)  # put this to the control argument later
#     # put this to the control argument later
#     cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#     posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
# 
#     while(length(posOutlier)==0){
#       r.outlier = r.outlier*0.8
#       cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#       posOutlier = which(mresid < cutoff[1] | mresid > cutoff[2])
#       cat("The current outlier ratio is:",r.outlier,"\n")
#       cat("The number of outlier is:",length(posOutlier),"\n")
# 
#     }
# 
# 
#     posValue = which(!is.na(mresid.temp))
#     posNonOutlier = setdiff(posValue, posOutlier)
# 
#     cat("The outlier of residuals will be passed to SPA analysis.\n")
#     cat("Cutoffs to define residuals:\t", signif(cutoff,2),"\n")
#     cat("Totally, ", length(posOutlier),"/", length(posValue), " are defined as outliers.\n")
# 
#     if(length(posOutlier) == 0)
#       stop("No outlier is observed. SPA is not required in this case.")
# 
#     # "-1" is to convert R style (index starting from 1) to C++ style (index starting from 0)
#     outLierList[[i]] = list(posValue = posValue - 1,
#                             posOutlier = posOutlier - 1,
#                             posNonOutlier = posNonOutlier - 1,
#                             resid = mresid.temp[posValue],
#                             resid2 = mresid.temp[posValue]^2,
#                             residOutlier = mresid.temp[posOutlier],
#                             residNonOutlier = mresid.temp[posNonOutlier],
#                             resid2NonOutlier = mresid.temp[posNonOutlier]^2)
#   }
# 
#   objNull = list(resid = mresid,
#                  # var.resid = var.resid,
#                  # tX = tX,
#                  # X.invXX = X.invXX,
#                  N = nrow(Cova),
#                  yVec = yVec,          # event variable: 0 or 1
#                  PCs = PCs,
#                  # posOutlier = posOutlier,
#                  nPheno = nPheno,
#                  outLierList = outLierList)
# 
#   class(objNull) = "SPAmixPlusV4_NULL_Model"
#   return(objNull)
# }





# # ok 2025-03-28
# # fit null model using SPAmixPlusV4 method
# fitNullModel.SPAmixPlusV4 = function(response, designMat, subjData,
#                                      # sparseGRM = NULL, # update on 2025-03-27
#                                      # sparseGRMFile_SPAmixPlus = NULL, # update on 2025-03-27
#                                      # ResidMat = NULL,  # update on 2025-03-27
#                                      control=list(OutlierRatio=1.5),                                
#                                      sparseGRM_SPAmixPlus = NULL, # update on 2025-03-27
#                                      sparseGRMFile_SPAmixPlus = NULL, # update on 2025-03-27
#                                      ...)
# {
#   # #### # update on 2024-09-11 ###############################################################
#   
#   cat(paste0("sparseGRMFile is :", sparseGRMFile_SPAmixPlus, "\n"))  
#   
#   sparseGRM = data.table::fread(sparseGRMFile_SPAmixPlus)
#   
#   cat("sparseGRM is\n")
#   
#   print(head(sparseGRM))
#   
#   sparseGRM$ID1 = as.character(sparseGRM$ID1); sparseGRM$ID2 = as.character(sparseGRM$ID2)
#   
#   SubjID.In.GRM = unique(c(sparseGRM$ID1, sparseGRM$ID2))
#   
#   
#   cat("SubjID.In.GRM is\n")
#   
#   print(head(SubjID.In.GRM))
#   
#   
#   ############# old SPAmix code #############################################################
#   
#   if(!class(response) %in% c("Surv", "Residual")) 
#     stop("For SPAmixPlusV4, the response variable should be of class 'Surv' or 'Residual'.")
#   
#   if(class(response) == "Surv")
#   {
#     formula = response ~ designMat
#     
#     obj.coxph = coxph(formula, x=T, ...)
#     
#     ### Check input arguments
#     # p2g = check_input(pIDs, gIDs, obj.coxph, range)
#     
#     y = obj.coxph$y
#     yVec = y[,ncol(y)]
#     
#     mresid = obj.coxph$residuals
#     Cova = obj.coxph$x
#     
#     if(length(mresid) != length(subjData))
#       stop("Please check the consistency between 'formula' and 'subjData'.")
#     
#     mresid = matrix(mresid, ncol=1)
#   }
#   
#   if(class(response) == "Residual")
#   {
#     yVec = mresid = response
#     # yVec = mresid = as.matrix(response)    
#     Cova = designMat
#     
#     # 调试输出 mresid 的类和维度
#     cat("Class of mresid:", class(mresid), "\n")
#     cat("Dimension of mresid:", dim(mresid), "\n")
#     
#     
#     print(head(mresid))
#     if(nrow(mresid) != length(subjData))
#       stop("Please check the consistency between 'formula' and 'subjData'.")
#   }
#   
#   PC_columns = control$PC_columns
#   
#   # Remove the below if checked later
#   cat("colnames(designMat):\n")
#   print(colnames(designMat))
#   cat("PC columns specified in 'control':\n")
#   print(PC_columns)
#   cat("dimension of 'designMat' and 'Cova':\n")
#   print(dim(designMat))
#   print(dim(Cova))
#   
#   if(any(!PC_columns %in% colnames(designMat)))
#     stop("PC columns specified in 'control$PC_columns' should be in 'formula'.")
#   
#   pos_col = match(PC_columns, colnames(designMat))
#   
#   PCs = Cova[,pos_col,drop=F]  
#   # X = cbind(1, PCs)
#   # X.invXX = X %*% solve(t(X)%*%X)
#   # tX = t(X)
#   
#   outLierList = list()
#   nPheno = ncol(mresid)
#   for(i in 1:nPheno)
#   {
#     mresid.temp = mresid[,i]
#     
#     # var.resid = var(mresid.temp, na.rm = T)
#     ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
#     q25 = quantile(mresid.temp, 0.25, na.rm = T)
#     q75 = quantile(mresid.temp, 0.75, na.rm = T)
#     IQR = q75 - q25
#     r.outlier =   r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)  # put this to the control argument later
#     # put this to the control argument later
#     cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#     posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
#     
#     while(length(posOutlier)==0){
#       r.outlier = r.outlier*0.8
#       cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#       posOutlier = which(mresid < cutoff[1] | mresid > cutoff[2])
#       cat("The current outlier ratio is:",r.outlier,"\n")
#       cat("The number of outlier is:",length(posOutlier),"\n")
#       
#     }
#     
#     
#     posValue = which(!is.na(mresid.temp))
#     posNonOutlier = setdiff(posValue, posOutlier)
#     
#     cat("The outlier of residuals will be passed to SPA analysis.\n")
#     cat("Cutoffs to define residuals:\t", signif(cutoff,2),"\n")
#     cat("Totally, ", length(posOutlier),"/", length(posValue), " are defined as outliers.\n")
#     
#     if(length(posOutlier) == 0)
#       stop("No outlier is observed. SPA is not required in this case.")
#     
#     # "-1" is to convert R style (index starting from 1) to C++ style (index starting from 0)
#     outLierList[[i]] = list(posValue = posValue - 1,
#                             posOutlier = posOutlier - 1,
#                             posNonOutlier = posNonOutlier - 1,
#                             resid = mresid.temp[posValue],
#                             resid2 = mresid.temp[posValue]^2,
#                             residOutlier = mresid.temp[posOutlier],
#                             residNonOutlier = mresid.temp[posNonOutlier],
#                             resid2NonOutlier = mresid.temp[posNonOutlier]^2)
#   }
#   
#   
#   
#   
#   # ResidMat = data.table::data.table(SubjID = subjData, Resid = mresid) # update on 2025-03-27
#   # ResidMat = data.table::data.table(SubjID = subjData, mresid) # update on 2025-03-27
#   
#   # 假设 mresid 是矩阵或数据框，且每列对应一个表型（nPheno列）
#   # # 需要明确命名残差列（例如 Resid_1, Resid_2, ..., Resid_nPheno）
#   # ResidMat = data.table::data.table(
#   #   SubjID = subjData,
#   #   mresid %>% as.data.table() %>% setnames(paste0("Resid_", 1:nPheno))
#   # )
#   
#   # 修正后的代码（方案1：直接命名列）
#   ResidMat = data.table::data.table(
#     SubjID = subjData,
#     # mresid
#     # as.data.frame(as.matrix(mresid))  # 关键修复点
#     as.data.frame(unclass(mresid))  # 使用 unclass() 彻底剥离类属性
#   )
#   colnames(ResidMat)[2:(nPheno+1)] = paste0("Resid_", 1:nPheno)  # 假设 mresid 有 nPheno 列
#   
#   
#   cat("ResidMat is\n")
#   
#   print(head(ResidMat))
#   
#   # #### # update on 2024-09-11 ###############################################################
#   # 
#   # print(head(sparseGRM))
#   # 
#   # sparseGRM$ID1 = as.character(sparseGRM$ID1); sparseGRM$ID2 = as.character(sparseGRM$ID2)
#   # 
#   SubjID.In.Resid = ResidMat$SubjID
#   
#   print(head(SubjID.In.Resid))
#   
#   # SubjID.In.GRM = unique(c(sparseGRM$ID1, sparseGRM$ID2))
#   
#   # print(head(SubjID.In.GRM))
#   
#   #
#   if(any(!SubjID.In.Resid %in% SubjID.In.GRM))
#     stop("At least one subject in residual matrix does not have GRM information.")
#   #
#   # # SubjID = SubjID.In.Resid
#   sparseGRM_new = sparseGRM %>% filter(ID1 %in% SubjID.In.Resid & ID2 %in% SubjID.In.Resid)
#   
#   #####
#   
#   objNull = list(resid = mresid,
#                  ResidMat = ResidMat,   # update on 2025-03-27
#                  sparseGRM_new = sparseGRM_new, # update on 2025-03-27
#                  subjData = subjData,
#                  # ResidMat = data.table::data.table(SubjID = subjData, Resid = mresid), # update on 2025-03-27
#                  # var.resid = var.resid,
#                  # tX = tX,
#                  # X.invXX = X.invXX,
#                  N = nrow(Cova),
#                  yVec = yVec,          # event variable: 0 or 1
#                  PCs = PCs,
#                  # posOutlier = posOutlier,
#                  nPheno = nPheno,
#                  outLierList = outLierList)
#   
#   class(objNull) = "SPAmixPlusV4_NULL_Model"
#   return(objNull)
# }








#### 2025-03-28 map ID #####################################################################################

# # fit null model using SPAmixPlusV4 method
# fitNullModel.SPAmixPlusV4 = function(response, designMat, subjData,
#                                      # sparseGRM = NULL, # update on 2025-03-27
#                                      # sparseGRMFile_SPAmixPlus = NULL, # update on 2025-03-27
#                                      # ResidMat = NULL,  # update on 2025-03-27
#                                      control=list(OutlierRatio=1.5),
#                                      sparseGRM_SPAmixPlus = NULL, # update on 2025-03-27
#                                      sparseGRMFile_SPAmixPlus = NULL, # update on 2025-03-27
#                                      ...)
# {
#   # #### # update on 2024-09-11 ###############################################################
#
#   cat(paste0("sparseGRMFile is :", sparseGRMFile_SPAmixPlus, "\n"))
#
#   sparseGRM = data.table::fread(sparseGRMFile_SPAmixPlus)
#
#   cat("sparseGRM is\n")
#
#   print(head(sparseGRM))
#
#   sparseGRM$ID1 = as.character(sparseGRM$ID1); sparseGRM$ID2 = as.character(sparseGRM$ID2)
#
#   SubjID.In.GRM = unique(c(sparseGRM$ID1, sparseGRM$ID2))
#
#
#   cat("SubjID.In.GRM is\n")
#
#   print(head(SubjID.In.GRM))
#
#
#   ############# old SPAmix code #############################################################
#
#   if(!class(response) %in% c("Surv", "Residual"))
#     stop("For SPAmixPlusV4, the response variable should be of class 'Surv' or 'Residual'.")
#
#   if(class(response) == "Surv")
#   {
#     formula = response ~ designMat
#
#     obj.coxph = coxph(formula, x=T, ...)
#
#     ### Check input arguments
#     # p2g = check_input(pIDs, gIDs, obj.coxph, range)
#
#     y = obj.coxph$y
#     yVec = y[,ncol(y)]
#
#     mresid = obj.coxph$residuals
#     Cova = obj.coxph$x
#
#     if(length(mresid) != length(subjData))
#       stop("Please check the consistency between 'formula' and 'subjData'.")
#
#     mresid = matrix(mresid, ncol=1)
#   }
#
#   if(class(response) == "Residual")
#   {
#     yVec = mresid = response
#     # yVec = mresid = as.matrix(response)
#     Cova = designMat
#
#     # 调试输出 mresid 的类和维度
#     cat("Class of mresid:", class(mresid), "\n")
#     cat("Dimension of mresid:", dim(mresid), "\n")
#
#
#     print(head(mresid))
#     if(nrow(mresid) != length(subjData))
#       stop("Please check the consistency between 'formula' and 'subjData'.")
#   }
#
#   PC_columns = control$PC_columns
#
#   # Remove the below if checked later
#   cat("colnames(designMat):\n")
#   print(colnames(designMat))
#   cat("PC columns specified in 'control':\n")
#   print(PC_columns)
#   cat("dimension of 'designMat' and 'Cova':\n")
#   print(dim(designMat))
#   print(dim(Cova))
#
#   if(any(!PC_columns %in% colnames(designMat)))
#     stop("PC columns specified in 'control$PC_columns' should be in 'formula'.")
#
#   pos_col = match(PC_columns, colnames(designMat))
#
#   PCs = Cova[,pos_col,drop=F]
#   # X = cbind(1, PCs)
#   # X.invXX = X %*% solve(t(X)%*%X)
#   # tX = t(X)
#
#   outLierList = list()
#   nPheno = ncol(mresid)
#   for(i in 1:nPheno)
#   {
#     mresid.temp = mresid[,i]
#
#     # var.resid = var(mresid.temp, na.rm = T)
#     ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
#     q25 = quantile(mresid.temp, 0.25, na.rm = T)
#     q75 = quantile(mresid.temp, 0.75, na.rm = T)
#     IQR = q75 - q25
#     r.outlier =   r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)  # put this to the control argument later
#     # put this to the control argument later
#     cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#     posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
#
#     while(length(posOutlier)==0){
#       r.outlier = r.outlier*0.8
#       cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#       posOutlier = which(mresid < cutoff[1] | mresid > cutoff[2])
#       cat("The current outlier ratio is:",r.outlier,"\n")
#       cat("The number of outlier is:",length(posOutlier),"\n")
#
#     }
#
#
#     posValue = which(!is.na(mresid.temp))
#     posNonOutlier = setdiff(posValue, posOutlier)
#
#     cat("The outlier of residuals will be passed to SPA analysis.\n")
#     cat("Cutoffs to define residuals:\t", signif(cutoff,2),"\n")
#     cat("Totally, ", length(posOutlier),"/", length(posValue), " are defined as outliers.\n")
#
#     if(length(posOutlier) == 0)
#       stop("No outlier is observed. SPA is not required in this case.")
#
#     # "-1" is to convert R style (index starting from 1) to C++ style (index starting from 0)
#     outLierList[[i]] = list(posValue = posValue - 1,
#                             posOutlier = posOutlier - 1,
#                             posNonOutlier = posNonOutlier - 1,
#                             resid = mresid.temp[posValue],
#                             resid2 = mresid.temp[posValue]^2,
#                             residOutlier = mresid.temp[posOutlier],
#                             residNonOutlier = mresid.temp[posNonOutlier],
#                             resid2NonOutlier = mresid.temp[posNonOutlier]^2)
#   }
#
#
#
#
#   # ResidMat = data.table::data.table(SubjID = subjData, Resid = mresid) # update on 2025-03-27
#   # ResidMat = data.table::data.table(SubjID = subjData, mresid) # update on 2025-03-27
#
#   # 假设 mresid 是矩阵或数据框，且每列对应一个表型（nPheno列）
#   # # 需要明确命名残差列（例如 Resid_1, Resid_2, ..., Resid_nPheno）
#   # ResidMat = data.table::data.table(
#   #   SubjID = subjData,
#   #   mresid %>% as.data.table() %>% setnames(paste0("Resid_", 1:nPheno))
#   # )
#
#   # 修正后的代码（方案1：直接命名列）
#   ResidMat = data.table::data.table(
#     SubjID = subjData,
#     # mresid
#     # as.data.frame(as.matrix(mresid))  # 关键修复点
#     as.data.frame(unclass(mresid))  # 使用 unclass() 彻底剥离类属性
#   )
#   colnames(ResidMat)[2:(nPheno+1)] = paste0("Resid_", 1:nPheno)  # 假设 mresid 有 nPheno 列
#
#
#   cat("ResidMat is\n")
#
#   print(head(ResidMat))
#
#   # #### # update on 2024-09-11 ###############################################################
#   #
#   # print(head(sparseGRM))
#   #
#   # sparseGRM$ID1 = as.character(sparseGRM$ID1); sparseGRM$ID2 = as.character(sparseGRM$ID2)
#   #
#   SubjID.In.Resid = ResidMat$SubjID
#
#   print(head(SubjID.In.Resid))
#
#   # SubjID.In.GRM = unique(c(sparseGRM$ID1, sparseGRM$ID2))
#
#   # print(head(SubjID.In.GRM))
#
#   #
#   if(any(!SubjID.In.Resid %in% SubjID.In.GRM))
#     stop("At least one subject in residual matrix does not have GRM information.")
#   #
#   # # SubjID = SubjID.In.Resid
#   sparseGRM_new = sparseGRM %>% filter(ID1 %in% SubjID.In.Resid & ID2 %in% SubjID.In.Resid)
#
#   #####
#
#   objNull = list(resid = mresid,
#                  ResidMat = ResidMat,   # update on 2025-03-27
#                  sparseGRM_new = sparseGRM_new, # update on 2025-03-27
#                  subjData = subjData,
#                  # ResidMat = data.table::data.table(SubjID = subjData, Resid = mresid), # update on 2025-03-27
#                  # var.resid = var.resid,
#                  # tX = tX,
#                  # X.invXX = X.invXX,
#                  N = nrow(Cova),
#                  yVec = yVec,          # event variable: 0 or 1
#                  PCs = PCs,
#                  # posOutlier = posOutlier,
#                  nPheno = nPheno,
#                  outLierList = outLierList)
#
#   class(objNull) = "SPAmixPlusV4_NULL_Model"
#   return(objNull)
# }


#### 20250407 map ID only new ID ------------------------------------------------------------------

# library(data.table)  # 确保:=操作符可用
# 
# fitNullModel.SPAmixPlusV4 = function(response, designMat, subjData,
#                                      control=list(OutlierRatio=1.5),
#                                      sparseGRM_SPAmixPlus = NULL,
#                                      sparseGRMFile_SPAmixPlus = NULL,
#                                      ...)
# {
#   # ---- 1. 读取稀疏GRM文件 ----
#   cat(paste0("sparseGRMFile is :", sparseGRMFile_SPAmixPlus, "\n"))
#   sparseGRM = data.table::fread(sparseGRMFile_SPAmixPlus)
#   data.table::setDT(sparseGRM)  # 新增：确保为data.table
#   cat("Initial sparseGRM:\n")
#   print(head(sparseGRM))
# 
# 
#   cat("Part1:\n")
# 
# 
#   # ---- 2. 处理响应变量（保持原有逻辑）----
#   if(!class(response) %in% c("Surv", "Residual"))
#     stop("For SPAmixPlusV4, the response variable should be of class 'Surv' or 'Residual'.")
# 
#   if(class(response) == "Surv") {
#     formula = response ~ designMat
#     obj.coxph = survival::coxph(formula, x=TRUE, ...)
#     y = obj.coxph$y
#     yVec = y[,ncol(y)]
#     mresid = residuals(obj.coxph)
#     Cova = obj.coxph$x
#     if(length(mresid) != length(subjData))
#       stop("Please check the consistency between 'formula' and 'subjData'.")
#     mresid = matrix(mresid, ncol=1)
#     nPheno = 1
#   } else if (class(response) == "Residual") {
#     yVec = mresid = response
#     Cova = designMat
#     if(nrow(mresid) != length(subjData))
#       stop("Please check the consistency between 'formula' and 'subjData'.")
#     nPheno = ncol(mresid)
#   }
# 
#   cat("Part2:\n")
# 
# 
#   # ---- 3. 处理主成分（PC）列 ----
#   PC_columns = control$PC_columns
#   if(any(!PC_columns %in% colnames(designMat)))
#     stop("PC columns specified in 'control$PC_columns' should be in 'formula'.")
#   pos_col = match(PC_columns, colnames(designMat))
#   PCs = Cova[,pos_col,drop=FALSE]
# 
#   cat("Part3:\n")
# 
#   # ---- 4. 检测异常值（完整保留原有逻辑）----
#   outLierList = list()
#   for(i in 1:nPheno) {
#     mresid.temp = mresid[,i]
#     q25 = quantile(mresid.temp, 0.25, na.rm = TRUE)
#     q75 = quantile(mresid.temp, 0.75, na.rm = TRUE)
#     IQR = q75 - q25
#     r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)
#     cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#     posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
# 
#     # 动态调整异常值阈值
#     while(length(posOutlier) == 0) {
#       r.outlier = r.outlier * 0.8
#       cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#       posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
#       cat("Adjusted outlier ratio:", r.outlier, "| Outliers found:", length(posOutlier), "\n")
#     }
# 
#     posValue = which(!is.na(mresid.temp))
#     posNonOutlier = setdiff(posValue, posOutlier)
# 
#     outLierList[[i]] = list(
#       posValue = posValue - 1,  # R索引转C++索引（从0开始）
#       posOutlier = posOutlier - 1,
#       posNonOutlier = posNonOutlier - 1,
#       resid = mresid.temp[posValue],
#       resid2 = (mresid.temp[posValue])^2,
#       residOutlier = mresid.temp[posOutlier],
#       residNonOutlier = mresid.temp[posNonOutlier],
#       resid2NonOutlier = (mresid.temp[posNonOutlier])^2
#     )
#   }
# 
#   cat("Part4:\n")
# 
# 
#   # ---- 5. 创建整数ID映射表 ----
#   data.table::set(sparseGRM, j = "ID1", value = as.character(sparseGRM$ID1))
#   data.table::set(sparseGRM, j = "ID2", value = as.character(sparseGRM$ID2))
# 
#   # 检查列名和类型
#   cat("\nsparseGRM列名:", names(sparseGRM))
#   cat("\nID1类型:", class(sparseGRM$ID1))  # 应输出 "character"
#   cat("\nID2类型:", class(sparseGRM$ID2))  # 应输出 "character"
# 
#   # 生成唯一ID列表
#   subjData = as.character(subjData)
#   all_ids = unique(c(subjData, sparseGRM$ID1, sparseGRM$ID2))
# 
#   # 显式定义id_map列名
#   id_map = data.table::data.table(
#     OriginalID = all_ids,
#     Index = seq_along(all_ids) - 1
#   )
#   data.table::setDT(id_map)  # 强制转换
#   data.table::setkey(id_map, "OriginalID")
# 
#   # 调试：输出id_map结构
#   cat("\nid_map列名:", names(id_map))
#   cat("\nid_map前3行:\n")
#   print(head(id_map, 3))
# 
# 
#   cat("Part5:\n")
# 
# 
#   # ---- 6. 转换ResidMat的SubjID为整数 ----
#   ResidMat = data.table::data.table(
#     SubjID = subjData,
#     as.data.frame(mresid)
#   )
#   data.table::setDT(ResidMat)
#   colnames(ResidMat)[2:(nPheno+1)] = paste0("Resid_", 1:nPheno)
# 
#   # 合并操作
#   ResidMat = data.table::merge.data.table(
#     ResidMat,
#     id_map[, c("OriginalID", "Index")],
#     by.x = "SubjID",
#     by.y = "OriginalID",
#     all.x = TRUE
#   )
# 
#   # 显式删除列
#   data.table::set(ResidMat, j = "SubjID", value = NULL)
#   data.table::setnames(ResidMat, "Index", "SubjID")
#   if (anyNA(ResidMat$SubjID)) stop("Missing SubjID after conversion.")
# 
#   cat("Part6:\n")
# 
#   # ---- 7. 转换稀疏GRM的ID为整数 ----
#   cat("\n===== sparseGRM验证 =====\n")
#   cat("sparseGRM列名:", names(sparseGRM), "\n")
#   cat("Value列是否存在:", "Value" %in% names(sparseGRM), "\n")
# 
#   # 显式提取列并验证类型
#   id1_values <- sparseGRM$ID1
#   id2_values <- sparseGRM$ID2
#   value_values <- sparseGRM$Value
# 
#   # 筛选符合条件的行
#   filter_condition <- id1_values %in% id_map$OriginalID & id2_values %in% id_map$OriginalID
# 
#   # 转换ID为索引
#   id1_indices <- id_map$Index[match(id1_values[filter_condition], id_map$OriginalID)]
#   id2_indices <- id_map$Index[match(id2_values[filter_condition], id_map$OriginalID)]
# 
#   # 构建新数据表
#   sparseGRM_new <- data.table::data.table(
#     ID1 = id1_indices,
#     ID2 = id2_indices,
#     Value = value_values[filter_condition]
#   )
# 
#   # 检查转换结果
#   cat("\nsparseGRM_new列名:", names(sparseGRM_new))
#   cat("\nsparseGRM_new前6行:\n")
#   print(head(sparseGRM_new))
# 
#   cat("Part7:\n")
# 
# 
#   # ---- 8. 验证数据一致性 ----
#   if(nrow(sparseGRM_new) == 0) stop("No valid GRM entries after ID conversion!")
#   if(anyNA(ResidMat$SubjID)) stop("Missing SubjID in ResidMat!")
# 
# 
#   cat("Part8:\n")
# 
#   # ---- 9. 构建完整Null模型对象 ----
#   objNull = list(
#     resid = mresid,
#     ResidMat = ResidMat,
#     sparseGRM = sparseGRM_new,
#     id_map = id_map,
#     subjData = subjData,
#     N = nrow(Cova),
#     yVec = yVec,
#     PCs = PCs,
#     nPheno = nPheno,
#     outLierList = outLierList,
#     control = control
#   )
#   class(objNull) = "SPAmixPlusV4_NULL_Model"
# 
#   cat("Part9:\n")
# 
# 
#   # ---- 10. 调试输出 ----
#   cat("\n===== Final Object Structure =====")
#   cat("\nResidMat dimensions:", dim(ResidMat))
#   cat("\nsparseGRM entries:", nrow(sparseGRM_new))
#   cat("\nID mapping table size:", nrow(id_map), "\n")
# 
#   return(objNull)
# }












# #### 20250407 map ID new ID and old ID ------------------------------------------------------------------
# library(data.table)
# library(survival)
# 
# fitNullModel.SPAmixPlusV4 = function(response, designMat, subjData,
#                                      control=list(OutlierRatio=1.5),
#                                      sparseGRM_SPAmixPlus = NULL,
#                                      sparseGRMFile_SPAmixPlus = NULL,
#                                      ...)
# {
#   # ---- 1. 读取稀疏GRM文件 ----
#   cat(paste0("sparseGRMFile is :", sparseGRMFile_SPAmixPlus, "\n"))
#   sparseGRM = data.table::fread(sparseGRMFile_SPAmixPlus)
#   data.table::setDT(sparseGRM)
#   cat("Initial sparseGRM:\n")
#   print(head(sparseGRM))
# 
#   cat("Part1:\n")
# 
#   # ---- 2. 严格类型检查 ----
#   if(!inherits(response, c("Surv", "Residual")))
#     stop("Response must be either a Surv object or Residual object")
# 
#   # ---- 3. 处理生存分析 ----
#   if(inherits(response, "Surv")) {
#     formula = response ~ designMat
#     obj.coxph = survival::coxph(formula, x=TRUE, ...)
#     y = obj.coxph$y
#     yVec = y[,ncol(y)]
#     mresid = residuals(obj.coxph)
#     Cova = obj.coxph$x
#     if(length(mresid) != length(subjData))
#       stop("CoxPH residuals length must match subjData length")
#     mresid = matrix(mresid, ncol=1)
#     nPheno = 1
#   }
# 
#   # ---- 4. 处理残差对象 ----
#   else if(inherits(response, "Residual")) {
#     # 强制转换并验证维度
#     if(!is.matrix(response)) mresid = as.matrix(response)
#     else mresid = response
# 
#     if(nrow(mresid) != length(subjData))
#       stop(paste(nrow(mresid), "residual rows vs", length(subjData), "subjects"))
# 
#     yVec = mresid
#     Cova = designMat
#     nPheno = ncol(mresid)
#   }
# 
#   cat("Part2:\n")
# 
#   # ---- 5. 处理主成分（PC）列 ----
#   PC_columns = control$PC_columns
#   if(any(!PC_columns %in% colnames(designMat)))
#     stop("PC columns specified in 'control$PC_columns' should be in 'formula'.")
#   pos_col = match(PC_columns, colnames(designMat))
#   PCs = Cova[,pos_col,drop=FALSE]
# 
#   cat("Part3:\n")
# 
#   # ---- 6. 检测异常值 ----
#   outLierList = list()
#   for(i in 1:nPheno) {
#     mresid.temp = mresid[,i]
#     q25 = quantile(mresid.temp, 0.25, na.rm = TRUE)
#     q75 = quantile(mresid.temp, 0.75, na.rm = TRUE)
#     IQR = q75 - q25
#     r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)
#     cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#     posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
# 
#     # 动态调整阈值
#     while(length(posOutlier) == 0) {
#       r.outlier = r.outlier * 0.8
#       cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
#       posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
#       cat("Adjusted outlier ratio:", r.outlier, "| Outliers found:", length(posOutlier), "\n")
#     }
# 
#     posValue = which(!is.na(mresid.temp))
#     posNonOutlier = setdiff(posValue, posOutlier)
# 
#     outLierList[[i]] = list(
#       posValue = posValue - 1,
#       posOutlier = posOutlier - 1,
#       posNonOutlier = posNonOutlier - 1,
#       resid = mresid.temp[posValue],
#       resid2 = (mresid.temp[posValue])^2,
#       residOutlier = mresid.temp[posOutlier],
#       residNonOutlier = mresid.temp[posNonOutlier],
#       resid2NonOutlier = (mresid.temp[posNonOutlier])^2
#     )
#   }
# 
#   cat("Part4:\n")
# 
#   # ---- 7. 创建ID映射表 ----
#   data.table::set(sparseGRM, j = "ID1", value = as.character(sparseGRM$ID1))
#   data.table::set(sparseGRM, j = "ID2", value = as.character(sparseGRM$ID2))
# 
#   subjData = as.character(subjData)
#   all_ids = unique(c(subjData, sparseGRM$ID1, sparseGRM$ID2))
# 
#   id_map = data.table::data.table(
#     OriginalID = all_ids,
#     Index = seq_along(all_ids) - 1
#   )
#   data.table::setDT(id_map)
#   data.table::setkey(id_map, "OriginalID")
# 
#   cat("Part5:\n")
# 
#   # ---- 8. 构建ResidMat（关键部分）----
#   # 创建基础结构
#   ResidMat = data.table::data.table(
#     SubjID = subjData,
#     SubjID_Index = id_map$Index[match(subjData, id_map$OriginalID)]
#   )
# 
#   # 添加残差列（确保列顺序）
#   resid_cols = paste0("Resid_", 1:nPheno)
#   for(i in 1:nPheno) {
#     data.table::set(ResidMat, j = resid_cols[i], value = mresid[,i])
#   }
# 
#   # 验证顺序
#   if(!identical(ResidMat$SubjID, subjData))
#     stop("ResidMat ID顺序异常!")
#   if(anyNA(ResidMat$SubjID_Index))
#     stop("存在未映射的SubjID")
# 
#   cat("Part6:\n")
# 
#   # ---- 9. 处理稀疏GRM ----
#   id1_values <- sparseGRM$ID1
#   id2_values <- sparseGRM$ID2
#   value_values <- sparseGRM$Value
# 
#   # 筛选有效行
#   valid_pairs <- id1_values %in% id_map$OriginalID & id2_values %in% id_map$OriginalID
# 
#   sparseGRM_new <- data.table::data.table(
#     ID1 = id1_values[valid_pairs],
#     ID2 = id2_values[valid_pairs],
#     ID1_Index = id_map$Index[match(id1_values[valid_pairs], id_map$OriginalID)],
#     ID2_Index = id_map$Index[match(id2_values[valid_pairs], id_map$OriginalID)],
#     Value = value_values[valid_pairs]
#   )
# 
#   cat("Part7:\n")
# 
#   # ---- 10. 最终验证 ----
#   if(nrow(sparseGRM_new) == 0) stop("转换后GRM为空!")
#   if(anyNA(ResidMat$SubjID_Index)) stop("存在无效的SubjID索引")
# 
#   cat("Part8:\n")
# 
#   # ---- 11. 构建结果对象 ----
#   objNull = list(
#     resid = mresid,
#     ResidMat = ResidMat,
#     sparseGRM = sparseGRM_new,
#     id_map = id_map,
#     subjData = subjData,
#     N = nrow(Cova),
#     yVec = yVec,
#     PCs = PCs,
#     nPheno = nPheno,
#     outLierList = outLierList,
#     control = control
#   )
#   class(objNull) = "SPAmixPlusV4_NULL_Model"
# 
#   # ---- 12. 调试输出 ----
#   cat("\n===== 最终对象结构 =====\n")
#   cat("ResidMat列:", names(ResidMat), "\n")
#   cat("sparseGRM列:", names(sparseGRM_new), "\n")
#   cat("ID映射表记录数:", nrow(id_map), "\n")
# 
#   return(objNull)
# }


##### 20250407 map ID new ID and old ID v2 ------------------------------------------------------------------

library(data.table)
library(survival)

fitNullModel.SPAmixPlusV4 = function(response, designMat, subjData,
                                     control=list(OutlierRatio=1.5),
                                     sparseGRM_SPAmixPlus = NULL,
                                     sparseGRMFile_SPAmixPlus = NULL,
                                     ...) 
{
  # ---- 1. 读取稀疏GRM文件（修复列名问题）----
  cat(paste0("sparseGRMFile is :", sparseGRMFile_SPAmixPlus, "\n"))
  sparseGRM = data.table::fread(
    file = sparseGRMFile_SPAmixPlus,
    col.names = c("ID1", "ID2", "Value")  # 关键修复：明确指定列名
  )
  data.table::setDT(sparseGRM)
  
  ########################### ID一致性过滤 ###########################
  # 获取所有GRM中的ID
  grm_ids = unique(c(sparseGRM$ID1, sparseGRM$ID2))
  # 过滤subjData只保留与GRM的交集
  keep_idx = subjData %in% grm_ids
  subjData = subjData[keep_idx]
  designMat = designMat[keep_idx, , drop = FALSE]
  # 过滤GRM只保留与subjData的交集
  sparseGRM = sparseGRM[ID1 %in% subjData & ID2 %in% subjData]
  
  cat("Initial sparseGRM:\n")
  print(head(sparseGRM))
  
  cat("Part1:\n")
  
  # ---- 2. 严格类型检查 ----
  if(!inherits(response, c("Surv", "Residual")))
    stop("Response must be either a Surv object or Residual object")
  
  # ---- 3. 处理生存分析 ----
  if(inherits(response, "Surv")) {
    formula = response ~ designMat
    obj.coxph = survival::coxph(formula, x=TRUE, ...)
    y = obj.coxph$y
    yVec = y[,ncol(y)]
    mresid = residuals(obj.coxph)
    Cova = obj.coxph$x
    if(length(mresid) != length(subjData))
      stop("CoxPH residuals length must match subjData length")
    mresid = matrix(mresid, ncol=1)
    nPheno = 1
  }
  
  # ---- 4. 处理残差对象 ----
  else if(inherits(response, "Residual")) {
    if(!is.matrix(response)) mresid = as.matrix(response)
    else mresid = response
    
    if(nrow(mresid) != length(subjData))
      stop(paste(nrow(mresid), "residual rows vs", length(subjData), "subjects"))
    
    yVec = mresid
    Cova = designMat
    nPheno = ncol(mresid)
  }
  
  cat("Part2:\n")
  
  # ---- 5. 处理主成分（PC）列 ----
  PC_columns = control$PC_columns
  if(any(!PC_columns %in% colnames(designMat)))
    stop("PC columns specified in 'control$PC_columns' should be in 'formula'.")
  pos_col = match(PC_columns, colnames(designMat))
  PCs = Cova[,pos_col,drop=FALSE]
  
  cat("Part3:\n")
  
  # ---- 6. 检测异常值 ----
  outLierList = list()
  for(i in 1:nPheno) {
    mresid.temp = mresid[,i]
    q25 = quantile(mresid.temp, 0.25, na.rm = TRUE)
    q75 = quantile(mresid.temp, 0.75, na.rm = TRUE)
    IQR = q75 - q25
    r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)
    cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
    
    while(length(posOutlier) == 0) {
      r.outlier = r.outlier * 0.8
      cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
      posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
      cat("Adjusted outlier ratio:", r.outlier, "| Outliers found:", length(posOutlier), "\n")
    }
    
    posValue = which(!is.na(mresid.temp))
    posNonOutlier = setdiff(posValue, posOutlier)
    
    outLierList[[i]] = list(
      posValue = posValue - 1,
      posOutlier = posOutlier - 1,
      posNonOutlier = posNonOutlier - 1,
      resid = mresid.temp[posValue],
      resid2 = (mresid.temp[posValue])^2,
      residOutlier = mresid.temp[posOutlier],
      residNonOutlier = mresid.temp[posNonOutlier],
      resid2NonOutlier = (mresid.temp[posNonOutlier])^2
    )
  }
  
  cat("Part4:\n")
  
  # ---- 7. 创建ID映射表 ----
  data.table::set(sparseGRM, j = "ID1", value = as.character(sparseGRM$ID1))
  data.table::set(sparseGRM, j = "ID2", value = as.character(sparseGRM$ID2))
  
  subjData = as.character(subjData)
  all_ids = unique(c(subjData, sparseGRM$ID1, sparseGRM$ID2))
  
  id_map = data.table::data.table(
    OriginalID = all_ids,
    Index = seq_along(all_ids) - 1
  )
  data.table::setDT(id_map)
  data.table::setkey(id_map, "OriginalID")
  
  cat("Part5:\n")
  
  # ---- 8. 构建ResidMat（关键部分）----
  ResidMat = data.table::data.table(
    SubjID = subjData,
    SubjID_Index = id_map$Index[match(subjData, id_map$OriginalID)]
  )
  
  resid_cols = paste0("Resid_", 1:nPheno)
  for(i in 1:nPheno) {
    data.table::set(ResidMat, j = resid_cols[i], value = mresid[,i])
  }
  
  if(!identical(ResidMat$SubjID, subjData))
    stop("ResidMat ID顺序异常!")
  if(anyNA(ResidMat$SubjID_Index))
    stop("存在未映射的SubjID")
  
  cat("Part6:\n")
  
  # ---- 9. 处理稀疏GRM ----
  id1_values <- sparseGRM$ID1
  id2_values <- sparseGRM$ID2
  value_values <- sparseGRM$Value
  
  valid_pairs <- id1_values %in% id_map$OriginalID & id2_values %in% id_map$OriginalID
  
  sparseGRM_new <- data.table::data.table(
    ID1 = id1_values[valid_pairs],
    ID2 = id2_values[valid_pairs],
    ID1_Index = id_map$Index[match(id1_values[valid_pairs], id_map$OriginalID)],
    ID2_Index = id_map$Index[match(id2_values[valid_pairs], id_map$OriginalID)],
    Value = value_values[valid_pairs]
  )
  
  cat("Part7:\n")
  
  # ---- 10. 最终验证 ----
  if(nrow(sparseGRM_new) == 0) stop("转换后GRM为空!")
  if(anyNA(ResidMat$SubjID_Index)) stop("存在无效的SubjID索引")
  
  cat("Part8:\n")
  
  # ---- 11. 构建结果对象 ----
  objNull = list(
    resid = mresid,
    ResidMat = ResidMat,
    sparseGRM = sparseGRM_new,
    id_map = id_map,
    subjData = subjData,
    N = nrow(Cova),
    yVec = yVec,
    PCs = PCs,
    nPheno = nPheno,
    outLierList = outLierList,
    control = control
  )
  class(objNull) = "SPAmixPlusV4_NULL_Model"
  
  # ---- 12. 调试输出 ----
  cat("\n===== 最终对象结构 =====\n")
  cat("ResidMat列:", names(ResidMat), "\n")
  cat("sparseGRM列:", names(sparseGRM_new), "\n")
  cat("ID映射表记录数:", nrow(id_map), "\n")
  
  return(objNull)
}






#########################################################################################################

# check the control list in null model fitting for SPACox method
checkControl.NullModel.SPAmixPlusV4 = function(control, traitType, ...)
{
  if(!traitType %in% c("time-to-event", "Residual"))
    stop("For 'SPAmixPlusV4' method, only traitType of 'time-to-event' or 'Residual' is supported.")
  
  if(is.null(control$PC_columns))
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') is required for 'SPAmixPlusV4' method.")
  
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
