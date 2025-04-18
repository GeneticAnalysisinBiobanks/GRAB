#' SPAGxEmixPlus method in GRAB package
#' 
#' SPAGxEmixPlus method is an empirical approach to analyzing complex traits (including but not limited to time-to-event trait) for unrelated samples in a large-scale biobank. SPAGxEmixPlus extend SPACox to support an admixture population or multiple populations. 
#' 
#' @details 
#' For ```SPAGxEmixPlus```, the confounding factors of SNP-derived PCs are required and should be specified in ```control```.
#' 
#' @examples 
#' # Step 1: fit a null model
#' PhenoFile = system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData = data.table::fread(PhenoFile, header = T)
#' N = nrow(PhenoData)
#' PhenoData = PhenoData %>% mutate(PC1 = rnorm(N), PC2 = rnorm(N))  # add two PCs, which are required for SPAGxEmixPlus
#' 
#' # Users can directly specify a time-to-event trait to analyze
#' obj.SPAGxEmixPlus = GRAB.NullModel(Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2, 
#'                             data = PhenoData, 
#'                             subjData = IID, 
#'                             method = "SPAGxEmixPlus", 
#'                             traitType = "time-to-event",
#'                             control = list(PC_columns = "PC1,PC2"))
#' 
#' @export
GRAB.SPAGxEmixPlus = function(){
  print("Check ?GRAB.SPAGxEmixPlus for more details about 'SPAGxEmixPlus' method.")
}

################### This file includes the following functions

# ------------ used in 'GRAB_Marker.R' -----------
# 1. checkControl.Marker.SPAGxEmixPlus(control)
# 2. setMarker.SPAGxEmixPlus(objNull, control)
# 3. mainMarker.SPAGxEmixPlus()

# check the control list in marker-level testing
# unified control list (such as nMarkersEachChunk) can be found in checkControl.Marker() in control.R
checkControl.Marker.SPAGxEmixPlus = function(control)
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

# SPAGxEmixPlusClass::SPAGxEmixPlusClass(arma::vec t_resid,
#                          arma::mat t_XinvXX,
#                          arma::mat t_tX,
#                          arma::mat t_PCs,
#                          int t_N,
#                          double t_SPA_Cutoff)

setMarker.SPAGxEmixPlus = function(objNull, control)
{
  # the below function is in 'Main.cpp'
  setSPAGxEmixPlusobjInCPP(objNull$resid,
                           # objNull$X.invXX,
                           # objNull$tX,
                           objNull$PCs,
                           objNull$N,
                           control$SPA_Cutoff,
                           objNull$outLierList,
                           
                           objNull$sparseGRM,   # update by Yuzhuo Ma
                           objNull$ResidMat,    # update by Yuzhuo Ma
                           objNull$E            # 新增环境因子参数
  )  
  
  # outLierList[[i]] = list(posValue = posValue - 1,
  #                         posOutlier = posOutlier - 1,
  #                         posNonOutlier = posNonOutlier - 1)
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,"."))
}

# mainMarker.SPAGxEmixPlus(genoType, genoIndex, outputColumns)
mainMarker.SPAGxEmixPlus = function(genoType, genoIndex, outputColumns, objNull)
{
  OutList = mainMarkerInCPP("SPAGxEmixPlus", genoType, genoIndex);
  
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




##### 20250409 map ID new ID and old ID v3 ------------------------------------------------------------------
library(data.table)
library(survival)

fitNullModel.SPAGxEmixPlus = function(response, designMat, subjData,
                                      control = list(OutlierRatio = 1.5),
                                      sparseGRM_SPAmixPlus = NULL,
                                      sparseGRMFile_SPAmixPlus = NULL,
                                      EnviColName,            # update by Yuzhuo Ma : column name of Environmental factor in designMat or PhenoMat
                                      # EnviData,            # update by Yuzhuo Ma : the first column is SubjID, and the second column is environmental factor
                                      ...) 
{
  
  cat("[DEBUG] EnviColName received:", EnviColName, "\n")
  cat("[DEBUG] designMat columns:", paste(colnames(designMat), collapse = ", "), "\n")
  
  # ---- 1. 加载必要包 ----
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("请安装data.table包：install.packages('data.table')")
  }
  data.table::setDTthreads(threads = 0)  # 启用多线程加速
  
  # ---- 2. 读取稀疏GRM文件 ----
  cat(paste0("sparseGRMFile is :", sparseGRMFile_SPAmixPlus, "\n"))
  sparseGRM = data.table::fread(sparseGRMFile_SPAmixPlus)
  data.table::setDT(sparseGRM)
  
  # ---- 3. 校验GRM格式 ----
  if (ncol(sparseGRM) != 3) {
    stop("GRM文件必须是三列格式：ID1, ID2, Value")
  }
  data.table::setnames(sparseGRM, c("ID1", "ID2", "Value"))
  
  cat("Initial sparseGRM:\n")
  print(head(sparseGRM))
  
  # ---- 4. 处理响应变量（生存分析/残差） ----
  if (!inherits(response, c("Surv", "Residual"))) {
    stop("响应变量必须是Surv或Residual类型")
  }
  
  if (inherits(response, "Surv")) {
    # ---- 生存分析逻辑 ----
    formula = response ~ designMat
    obj.coxph = survival::coxph(formula, x = TRUE, ...)
    y = obj.coxph$y
    yVec = y[, ncol(y)]
    mresid = residuals(obj.coxph)
    Cova = obj.coxph$x
    
    cat(paste0("Class of is designMat is :", class(designMat), "\n"))
    head(designMat)
    
    # mresid.by.E = mresid * designMat$EnviColName # update by Yuzhuo Ma
    mresid.by.E = mresid * designMat[, EnviColName]  # 关键修正：矩阵列名索引
    # 在函数内添加检查
    if (!EnviColName %in% colnames(designMat)) {
      stop(paste("环境因子列", EnviColName, "不存在于 designMat 中"))
    }
    
    if (length(mresid) != length(subjData)) {
      stop("CoxPH残差长度必须与subjData一致")
    }
    
    mresid = matrix(mresid, ncol = 1)
    nPheno = 1
  } else if (inherits(response, "Residual")) {
    
    cat(paste0("Class of is designMat is :", class(designMat), "\n"))
    head(designMat)
    
    # ---- 残差对象逻辑 ----
    if (!is.matrix(response)) {
      mresid = as.matrix(response)
      mresid.by.E = mresid * designMat[, EnviColName] # update by Yuzhuo Ma
    } else {
      mresid = response
      # mresid.by.E = mresid * designMat$EnviColName # update by Yuzhuo Ma
      mresid.by.E = mresid * designMat[, EnviColName]  # 关键修正：矩阵列名索引
      # 在函数内添加检查
      if (!EnviColName %in% colnames(designMat)) {
        stop(paste("环境因子列", EnviColName, "不存在于 designMat 中"))
      }
      
    }
    
    if (nrow(mresid) != length(subjData)) {
      stop(paste("残差行数（", nrow(mresid), "）与样本数（", length(subjData), "）不匹配"))
    }
    
    yVec = mresid
    Cova = designMat
    nPheno = ncol(mresid)
  }
  
  # ---- 5. 处理主成分（PC）列（保持原样） ----
  PC_columns = control$PC_columns
  if (any(!PC_columns %in% colnames(designMat))) {
    stop("control$PC_columns中指定的PC列必须在设计矩阵中存在")
  }
  pos_col = match(PC_columns, colnames(designMat))
  PCs = Cova[, pos_col, drop = FALSE]
  
  # ---- 6. 检测异常值（保持原样） ----
  outLierList = list()
  for (i in 1:nPheno) {
    mresid.temp = mresid[, i]
    q25 = quantile(mresid.temp, 0.25, na.rm = TRUE)
    q75 = quantile(mresid.temp, 0.75, na.rm = TRUE)
    IQR = q75 - q25
    r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)
    cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
    
    # 动态调整阈值
    while (length(posOutlier) == 0) {
      r.outlier = r.outlier * 0.8
      cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
      posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
      cat("调整后异常值阈值:", r.outlier, "| 发现异常值:", length(posOutlier), "\n")
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
  
  
  # ---- 6. 检测异常值（Resid x Envi） ---- update by Yuzhuo Ma
  outLierList_ResidByE = list()
  for (i in 1:nPheno) {
    mresid.by.E.temp = mresid.by.E[, i]
    q25 = quantile(mresid.by.E.temp, 0.25, na.rm = TRUE)
    q75 = quantile(mresid.by.E.temp, 0.75, na.rm = TRUE)
    IQR = q75 - q25
    r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)
    cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier = which(mresid.by.E.temp < cutoff[1] | mresid.by.E.temp > cutoff[2])
    
    # 动态调整阈值
    while (length(posOutlier) == 0) {
      r.outlier = r.outlier * 0.8
      cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
      posOutlier = which(mresid.by.E.temp < cutoff[1] | mresid.by.E.temp > cutoff[2])
      cat("调整后异常值阈值:", r.outlier, "| 发现异常值:", length(posOutlier), "\n")
    }
    
    posValue = which(!is.na(mresid.by.E.temp))
    posNonOutlier = setdiff(posValue, posOutlier)
    
    outLierList_ResidByE[[i]] = list(
      posValue = posValue - 1,
      posOutlier = posOutlier - 1,
      posNonOutlier = posNonOutlier - 1,
      resid.by.Envi = mresid.by.E.temp[posValue],
      resid.by.Envi2 = (mresid.by.E.temp[posValue])^2,
      resid.by.EnviOutlier = mresid.by.E.temp[posOutlier],
      resid.by.EnviNonOutlier = mresid.by.E.temp[posNonOutlier],
      resid.by.Envi2NonOutlier = (mresid.by.E.temp[posNonOutlier])^2
    )
  }
  
  
  # ---- 7. 创建ID映射表（保持原样） ----
  data.table::set(sparseGRM, j = "ID1", value = as.character(sparseGRM$ID1))
  data.table::set(sparseGRM, j = "ID2", value = as.character(sparseGRM$ID2))
  subjData = as.character(subjData)
  
  all_ids = unique(c(subjData, sparseGRM$ID1, sparseGRM$ID2))
  
  id_map = data.table::data.table(
    OriginalID = all_ids,
    Index = as.integer(seq_along(all_ids) - 1)  # 强制为整型
  )
  data.table::setkey(id_map, "OriginalID")
  
  # ---- 8. 构建ResidMat（修复data.table作用域问题） ----
  ResidMat = data.table::data.table(
    SubjID = subjData,
    SubjID_Index = as.integer(id_map$Index[match(subjData, id_map$OriginalID)])
  )
  
  # 动态添加Resid_*列（使用data.table::set避免作用域问题）
  resid_cols = paste0("Resid_", 1:nPheno)
  for (i in 1:nPheno) {
    data.table::set(ResidMat, j = resid_cols[i], value = as.numeric(mresid[, i]))
  }
  
  # 验证ResidMat数据类型
  if (!all(sapply(ResidMat[, .SD, .SDcols = patterns("^Resid_")], is.numeric))) {
    stop("Resid_*列必须为双精度（numeric）")
  }
  if (!is.integer(ResidMat$SubjID_Index)) {
    stop("SubjID_Index必须为整型（integer）")
  }
  
  
  
  # ---- 8. 构建ResidByEnviMat ---- update by Yuzhuo Ma
  ResidByEnviMat = data.table::data.table(
    SubjID = subjData,
    SubjID_Index = as.integer(id_map$Index[match(subjData, id_map$OriginalID)])
  )
  
  # 动态添加ResidByEnvi_*列（使用data.table::set避免作用域问题）
  ResidByEnvi_cols = paste0("ResidByEnvi_", 1:nPheno)
  for (i in 1:nPheno) {
    data.table::set(ResidByEnviMat, j = ResidByEnvi_cols[i], value = as.numeric(mresid.by.E[, i]))
  }
  
  # 验证ResidByEnviMat数据类型
  if (!all(sapply(ResidByEnviMat[, .SD, .SDcols = patterns("^ResidByEnvi_")], is.numeric))) {
    stop("ResidByEnvi_*列必须为双精度（numeric）")
  }
  if (!is.integer(ResidByEnviMat$SubjID_Index)) {
    stop("SubjID_Index必须为整型（integer）")
  }
  
  
  
  
  
  
  
  # ---- 9. 处理稀疏GRM（保持原样） ----
  sparseGRM_new = data.table::data.table(
    ID1 = sparseGRM$ID1,
    ID2 = sparseGRM$ID2,
    ID1_Index = as.integer(id_map$Index[match(sparseGRM$ID1, id_map$OriginalID)]),
    ID2_Index = as.integer(id_map$Index[match(sparseGRM$ID2, id_map$OriginalID)]),
    Value = as.numeric(sparseGRM$Value)
  )
  
  # 移除无效行并验证
  sparseGRM_new = na.omit(sparseGRM_new)
  if (nrow(sparseGRM_new) == 0) {
    stop("转换后的稀疏GRM为空，请检查ID映射")
  }
  if (!is.integer(sparseGRM_new$ID1_Index) || !is.integer(sparseGRM_new$ID2_Index)) {
    stop("GRM索引列必须为整型")
  }
  
  # ---- 10. 构建最终对象 ----
  objNull = list(
    resid = mresid,
    resid.by.E = as.numeric(mresid.by.E),                             # update by Yuzhuo Ma
    ResidMat = as.data.frame(ResidMat),
    E = as.numeric(designMat[, EnviColName]),
    ResidByEnviMat = as.data.frame(ResidByEnviMat),       # update by Yuzhuo Ma
    sparseGRM = as.data.frame(sparseGRM_new),
    id_map = id_map,
    subjData = subjData,
    N = nrow(Cova),
    yVec = yVec,
    PCs = PCs,
    nPheno = nPheno,
    outLierList = outLierList,
    outLierList_ResidByE = outLierList_ResidByE,          # update by Yuzhuo Ma
    control = control
  )
  class(objNull) = "SPAGxEmixPlus_NULL_Model"
  
  # ---- 11. 调试输出 ----
  cat("\n===== 最终对象结构 =====\n")
  cat("ResidMat列类型:\n")
  print(sapply(ResidMat, class))
  cat("ResidByEnviMat列类型:\n")
  print(sapply(ResidByEnviMat, class))
  cat("\nsparseGRM列类型:\n")
  print(sapply(sparseGRM_new, class))
  
  return(objNull)
}


#########################################################################################################

# check the control list in null model fitting for SPACox method
checkControl.NullModel.SPAGxEmixPlus = function(control, traitType, ...) # update by Yuzhuo Ma 
{
  if(!traitType %in% c("time-to-event", "Residual"))
    stop("For 'SPAGxEmixPlus' method, only traitType of 'time-to-event' or 'Residual' is supported.")
  
  if(is.null(control$PC_columns))
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') is required for 'SPAGxEmixPlus' method.")
  
  if(length(control$PC_columns) != 1)
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') should be a character, not a character vector.")
  
  control$PC_columns = unlist(strsplit(control$PC_columns, split=","))
  if(length(control$PC_columns) == 1)
    warning("We detected that only one PC column exists, is that what you want? Note that control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') should be a character splitted using ','.")
  
  
  # # update by Yuzhuo Ma 
  # 
  # SubjID.In.Resid = ResidMat$SubjID
  # SubjID.In.GRM = unique(c(SparseGRM$ID1, SparseGRM$ID2))
  # 
  # if(any(!SubjID.In.Resid %in% SubjID.In.GRM))
  #   stop("At least one subject in residual matrix does not have GRM information.")
  
  # default.control = list(range = c(-100, 100),
  #                        length.out = 10000)
  
  # control = updateControl(control, default.control)
  
  # check the parameters
  return(control)
}
