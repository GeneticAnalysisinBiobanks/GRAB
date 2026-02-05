
#' SPAmixPlus method in SPAmixPlus package
#' 
#' SPAmixPlus method is an empirical approach to analyzing complex traits (including but not limited to time-to-event trait) for admixed populations in a large-scale biobank.
#' 
#' @param formula an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model.
#' @param subset an optional vector specifying a subset of observations to be used in the fitting process.
#' @param subjData a character vector specifing the subject IDs.
#' @param method a character of "SPAmixPlus" (default).
#' @param traitType a character of "time-to-event" (default) or "One of the supported trait types".
#' @param control a list of control parameters. 
#' @param sparseGRMFile a character of sparse GRM file.
#' @param sparseGRM_SPAmixPlus an R object for sparse GRM.
#' @param sparseGRMFile_SPAmixPlus a character of sparse GRM file for SPAmixPlus.
#' @param ... other arguments.
#' 
#' @details
#' The \code{SPAmixPlus.NullModel} function fits a null model for SPAmixPlus analysis.
#' It handles phenotype adjustment and prepares the object for subsequent marker-level or region-level analysis.
#' 
#' @return An object of class \code{SPAmixPlus_NULL_Model} containing the fitted null model information.
#' 
#' @examples
#' \dontrun{
#' # objNull = SPAmixPlus.NullModel(resid ~ Cov1 + Cov2 + PC1 + PC2,
#' #                                data = Pheno.mtx,
#' #                                subjData = IID_Vec, 
#' #                                sparseGRMFile_SPAmixPlus = sparse_GRM_file)
#' }
#' 
#' @export
SPAmixPlus.NullModel = function(formula,
                                data = NULL,
                                subset = NULL,
                                subjData,
                                method = "SPAmixPlus",
                                traitType = "time-to-event",
                                control = NULL,
                                sparseGRMFile = NULL,
                                sparseGRM_SPAmixPlus = NULL,     
                                sparseGRMFile_SPAmixPlus = NULL, 
                                ...)
{
  if(missing(subjData))
    stop("Argument 'subjData' is required to specify the subjects IDs in 'formula' and/or 'data'.")

  # Handle sparseGRMFile argument aliases
  if(is.null(sparseGRMFile_SPAmixPlus)) sparseGRMFile_SPAmixPlus = sparseGRMFile
  
  Call = match.call()
  
  # --- Formula Parsing (Copied from GRAB.NullModel) ---
  mf <- match.call(expand.dots = FALSE)
  
  LeftInFormula = deparse(formula[[2]])
  LeftIncludesAdd = grepl("\\+", LeftInFormula)

  if(LeftIncludesAdd){
    if(traitType != "Residual")
      stop("Only 'SPAmixPlus' method with traitType of 'Residual' supports multiple response variables in 'formula'.")
    
    nInLeft = length(strsplit(LeftInFormula, "\\+")[[1]])
    cat("SPAmixPlus method supports multiple response variables of model residuals.\n")
    RightInFormula = deparse(formula[[3]])
    NewLeftInFormla = paste0("paste(", gsub("\\+", ",", LeftInFormula), ")")
    NewRightInFormula = paste0(RightInFormula, collapse = " ")   
    mf$formula = as.formula(paste(NewLeftInFormla, "~", NewRightInFormula))
  }
  
  m <- match(x = c("formula", "data", "subset", "subjData"), 
             table = names(mf), nomatch = 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  
  mt = attr(x = mf, which = "terms")
  response = model.response(mf)
  designMat = model.matrix(object = mt, data = mf)
  
  # Ensure designMat has column names
  if(is.null(colnames(designMat)))
      colnames(designMat) = paste0("Cov", 1:ncol(designMat))

  subjData = model.extract(mf, "subjData")
  
  if(traitType == "Residual"){
    if(LeftIncludesAdd){
      noValueInAnyPheno = paste(rep(NA, nInLeft), collapse = " ")
      posNoValue = which(response == noValueInAnyPheno)
      response.temp = response
      
      if(length(posNoValue) > 0){
        cat("We remove", length(posNoValue), "individuals without any phenotyeps in analysis.\n")
        response.temp = response[-1*posNoValue]
        designMat = designMat[-1*posNoValue,,drop=F]
        subjData = subjData[-1*posNoValue]
      }
      
      nRes = length(response.temp)
      response = matrix(NA, nRes, nInLeft)
      for(i in 1:nRes)
        response[i,] = as.numeric(unlist(strsplit(response.temp[i], split = " ")))
      
    }else{
      response = matrix(response, ncol=1)
    }
    class(response) = "Residual" 
  }
  
  if(colnames(designMat)[1] == "(Intercept)")
    designMat = designMat[,-1,drop=F]
    
  if(any(duplicated(subjData))) 
    stop("Duplicated subject IDs in 'subjData' is not supported!")

  # --- Check Control ---
  control = checkControl.NullModel.SPAmixPlus(control, traitType)
  
  # --- Call Fitting Function ---
  objNull = fitNullModel.SPAmixPlus(response, designMat, subjData, control, sparseGRMFile_SPAmixPlus = sparseGRMFile_SPAmixPlus, ...)
  
  # objNull$subjData = subjData # fitNullModel.SPAmixPlus includes it
  objNull$Call = Call
  objNull$sessionInfo = sessionInfo()
  objNull$time = paste0("Complete Time: ",Sys.time())
  # objNull$control = control # fitNullModel.SPAmixPlus includes it
  objNull$method = "SPAmixPlus" 
  
  cat("Complete the null model fitting in package SPAmixPlus:\t", objNull$time,"\n")
  return(objNull)
}

# unified control list (such as nMarkersEachChunk) can be found in checkControl.Marker() in control.R
checkControl.Marker.SPAmixPlus = function(control)
{
  default.control = list(SPA_Cutoff = 2,
                         dosage_option = "rounding_first")  # "rounding_first" or "rounding_last"
  
  control = updateControl(control, default.control)  # This function is in 'Framework.R'
  
  # check the parameter
  if(!control$dosage_option %in% c("rounding_first", "rounding_last"))
    stop("control$dosage_option should be 'rounding_first' or 'rounding_last'.")
  
  return(control)
}

setMarker.SPAmixPlus = function(objNull, control)
{
  # the below function is in 'Main.cpp'
  setSPAmixPlusobjInCPP(objNull$resid,
                          objNull$PCs,
                          objNull$N,
                          control$SPA_Cutoff,
                          objNull$outLierList,
                          objNull$sparseGRM,  # update by Yuzhuo Ma
                          objNull$ResidMat    # update by Yuzhuo Ma
                          )  
  
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,"."))
}

# mainMarker.SPAmixPlus(genoType, genoIndex, outputColumns)
mainMarker.SPAmixPlus = function(genoType, genoIndex, outputColumns, objNull, control = NULL)
{
  # Handle step0 calling (optional)
  afFile = ""
  if(!is.null(control) && !is.null(control$afFile)) {
     afFile = control$afFile
     if(!file.exists(afFile)) warning(paste("AF File not found:", afFile))
  }
  
  OutList = mainMarkerInCPP("SPAmixPlus", genoType, genoIndex, character(), afFile);
  
  nPheno = objNull$nPheno;

  obj.mainMarker = data.frame(Pheno = paste0("pheno_", 1:nPheno),
                              Marker = rep(OutList$markerVec, each = nPheno),           # marker IDs
                              Info = rep(OutList$infoVec, each = nPheno),               # marker information: CHR:POS:REF:ALT
                              AltFreq = rep(OutList$altFreqVec, each = nPheno),         # alternative allele frequencies
                              AltCounts = rep(OutList$altCountsVec, each = nPheno),     # alternative allele counts
                              MissingRate = rep(OutList$missingRateVec, each = nPheno), # alternative allele counts
                              S = OutList$Stat,                                         # Score
                              S_mean = OutList$StatMean,                                # Expectation of Score
                              VarS = OutList$StatVar,                                   # Variance of Score
                              zScore = OutList$zScore,                                  # zScore
                              Pvalue = OutList$pvalVec,                                 # marker-level p-values
                              BetaG = OutList$BetaG)                                    # genetic effect size

  # optionalColumns = c("beta", "seBeta", "zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  optionalColumns = c("PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  additionalColumns = intersect(optionalColumns, outputColumns)
  
  if(length(additionalColumns) > 0)
    obj.mainMarker = cbind.data.frame(obj.mainMarker, 
                                      as.data.frame(OutList[additionalColumns]))
  
  return(obj.mainMarker)
}


#' Step 0 of SPAmixPlus: Pre-calculate individual-specific allele frequency coefficients
#' @export
SPAmixPlus.AF = function(GenoFile,
                         objNull,
                         outputFile,
                         GenoFileIndex = NULL,
                         control = list()){

    # Validate objNull
    if(!inherits(objNull, "SPAmixPlus_NULL_Model")) stop("objNull must be of class SPAmixPlus_NULL_Model")
    
    # Initialize control defaults if necessary
    if(is.null(control)) control = list()
    if(is.null(control$impute_method)) control$impute_method = "mean"
    if(is.null(control$missing_cutoff)) control$missing_cutoff = 0.15
    if(is.null(control$min_maf_marker)) control$min_maf_marker = 0.0001
    if(is.null(control$min_mac_marker)) control$min_mac_marker = 1
    if(is.null(control$omp_num_threads)) control$omp_num_threads = 1
    if(is.null(control$SPA_Cutoff)) control$SPA_Cutoff = 2
    
    # 1. Use internal setGenoInput to get ALL marker information
    # We pass 'AllMarkers=TRUE' to ensure no filtering happens at this stage unless specified
    control$AllMarkers = TRUE
    
    # Note: setGenoInput is internal in Geno.R. We can access it if inside package, or copy definition.
    # Since we are modifying the package, we can just call it.
    # However, setGenoInput is not exported. It is available in namespace.
    
    objGeno = setGenoInput(GenoFile, 
                           GenoFileIndex = GenoFileIndex, 
                           SampleIDs = as.character(objNull$subjData), 
                           control = control)
                           
    genoType = objGeno$genoType
    markerInfo = objGeno$markerInfo
    genoIndex = markerInfo$genoIndex
    
    if(length(genoIndex) == 0) stop("No available markers to process for Step 0.")
    
    cat("Step 0: Pre-calculating AF models for", length(genoIndex), "markers...\n")
    
    # 2. Initialize Global Vars in C++
    Group = makeGroup(objNull$yVec)
    
    setMarker_GlobalVarsInCPP(control$impute_method,
                              control$missing_cutoff,
                              control$min_maf_marker,
                              control$min_mac_marker,
                              control$omp_num_threads,
                              Group,
                              FALSE, # ifOutGroup
                              length(unique(Group)))

    # 3. Initialize SPAmixPlus Object using Null Model info
    setMarker.SPAmixPlus(objNull, control)
    
    # 4. Run Export
    exportAFModelInCPP("SPAmixPlus", genoType, genoIndex, outputFile);
    
    cat("Step 0 Completed. Model file saved to:", outputFile, "\n")
}

##### 20250409 map ID new ID and old ID v3 ------------------------------------------------------------------
#' @import data.table
#' @import survival
fitNullModel.SPAmixPlus = function(response, designMat, subjData,
                                     control = list(OutlierRatio = 1.5),
                                     sparseGRM_SPAmixPlus = NULL,
                                     sparseGRMFile_SPAmixPlus = NULL,
                                     ...) 
{
  # ---- 1. Load necessary packages ----
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install data.table package: install.packages('data.table')")
  }
  data.table::setDTthreads(threads = 0)  # Enable multi-threading
  
  # ---- 2. Read sparse GRM file ----
  if (is.null(sparseGRMFile_SPAmixPlus) || !file.exists(sparseGRMFile_SPAmixPlus)) {
      stop("sparseGRMFile_SPAmixPlus is required and must exist.")
  }
  cat(paste0("sparseGRMFile is :", sparseGRMFile_SPAmixPlus, "\n"))
  sparseGRM = data.table::fread(sparseGRMFile_SPAmixPlus)
  data.table::setDT(sparseGRM)
  
  # ---- 3. Validate GRM format ----
  if (ncol(sparseGRM) != 3) {
    stop("GRM file must have 3 columns: ID1, ID2, Value")
  }
  data.table::setnames(sparseGRM, c("ID1", "ID2", "Value"))
  
  cat("Initial sparseGRM:\n")
  print(head(sparseGRM))
  
  # ---- 4. Process response variable (Survival/Residual) ----
  if (!inherits(response, c("Surv", "Residual"))) {
    stop("Response must be of type Surv or Residual")
  }
  
  if (inherits(response, "Surv")) {
    # ---- Survival analysis logic ----
    formula = response ~ designMat
    obj.coxph = survival::coxph(formula, x = TRUE, ...)
    y = obj.coxph$y
    if(ncol(y) == 2){
        yVec = y[, 2] # Event
    } else {
        yVec = y[, ncol(y)] 
    }
    
    mresid = residuals(obj.coxph)
    Cova = obj.coxph$x
    
    if (length(mresid) != length(subjData)) {
      stop("CoxPH residuals length must match subjData length")
    }
    
    mresid = matrix(mresid, ncol = 1)
    nPheno = 1
  } else if (inherits(response, "Residual")) {
    # ---- Residual object logic ----
    if (!is.matrix(response)) {
      mresid = as.matrix(response)
    } else {
      mresid = response
    }
    
    if (nrow(mresid) != length(subjData)) {
      stop(paste("Residual rows (", nrow(mresid), ") do not match sample size (", length(subjData), ")"))
    }
    
    yVec = mresid[,1] # Just for shape
    Cova = designMat
    nPheno = ncol(mresid)
  }
  
  # ---- 5. Process Principal Components (PC) columns ----
  PC_columns = control$PC_columns
  # Remove extra spaces if any
  if(!is.null(PC_columns)) PC_columns = trimws(PC_columns)
  
  if (is.null(PC_columns) || any(!PC_columns %in% colnames(designMat))) {
    missing_pcs = PC_columns[!PC_columns %in% colnames(designMat)]
    stop(paste("PC columns specified in control$PC_columns must exist in the design matrix. Missing:", paste(missing_pcs, collapse=", ")))
  }
  pos_col = match(PC_columns, colnames(designMat))
  PCs = Cova[, pos_col, drop = FALSE]
  
  # ---- 6. Detect outliers ----
  outLierList = list()
  for (i in 1:nPheno) {
    mresid.temp = mresid[, i]
    q25 = quantile(mresid.temp, 0.25, na.rm = TRUE)
    q75 = quantile(mresid.temp, 0.75, na.rm = TRUE)
    IQR = q75 - q25
    r.outlier = ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio)
    cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
    
    # Dynamically adjust threshold
    while (length(posOutlier) == 0 && r.outlier > 0.1) {
      r.outlier = r.outlier * 0.8
      cutoff = c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
      posOutlier = which(mresid.temp < cutoff[1] | mresid.temp > cutoff[2])
      cat("Adjusted outlier threshold:", r.outlier, "| Outliers found:", length(posOutlier), "\n")
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
  
  # ---- 7. Create ID mapping table ----
  data.table::set(sparseGRM, j = "ID1", value = as.character(sparseGRM$ID1))
  data.table::set(sparseGRM, j = "ID2", value = as.character(sparseGRM$ID2))
  subjData = as.character(subjData)
  
  all_ids = unique(c(subjData, sparseGRM$ID1, sparseGRM$ID2))
  
  id_map = data.table::data.table(
    OriginalID = all_ids,
    Index = as.integer(seq_along(all_ids) - 1)  # Force integer
  )
  data.table::setkey(id_map, "OriginalID")
  
  # ---- 8. Build ResidMat ----
  ResidMat = data.table::data.table(
    SubjID = subjData,
    SubjID_Index = as.integer(id_map$Index[match(subjData, id_map$OriginalID)])
  )
  
  # Dynamically add Resid_* columns
  resid_cols = paste0("Resid_", 1:nPheno)
  for (i in 1:nPheno) {
    data.table::set(ResidMat, j = resid_cols[i], value = as.numeric(mresid[, i]))
  }
  
  # Validate ResidMat data types
  if (!all(sapply(ResidMat[, .SD, .SDcols = patterns("^Resid_")], is.numeric))) {
    stop("Resid_* columns must be numeric (double)")
  }
  if (!is.integer(ResidMat$SubjID_Index)) {
    stop("SubjID_Index must be integer")
  }
  
  # ---- 9. Process Sparse GRM ----
  sparseGRM_new = data.table::data.table(
    ID1 = sparseGRM$ID1,
    ID2 = sparseGRM$ID2,
    ID1_Index = as.integer(id_map$Index[match(sparseGRM$ID1, id_map$OriginalID)]),
    ID2_Index = as.integer(id_map$Index[match(sparseGRM$ID2, id_map$OriginalID)]),
    Value = as.numeric(sparseGRM$Value)
  )
  
  # Remove invalid rows and validate
  sparseGRM_new = na.omit(sparseGRM_new)
  if (nrow(sparseGRM_new) == 0) {
    stop("Converted sparse GRM is empty, please check ID mapping")
  }
  if (!is.integer(sparseGRM_new$ID1_Index) || !is.integer(sparseGRM_new$ID2_Index)) {
    stop("GRM index columns must be integer")
  }
  
  # ---- 10. Build final object ----
  objNull = list(
    resid = mresid,
    ResidMat = as.data.frame(ResidMat),
    sparseGRM = as.data.frame(sparseGRM_new),
    id_map = id_map,
    subjData = subjData,
    N = nrow(Cova),
    yVec = yVec,
    PCs = PCs,
    nPheno = nPheno,
    outLierList = outLierList,
    control = control
  )
  class(objNull) = "SPAmixPlus_NULL_Model"
  
  # ---- 11. Debug output ----
  cat("\n===== Final Object Structure =====\n")
  cat("ResidMat column types:\n")
  print(sapply(ResidMat, class))
  cat("\nsparseGRM column types:\n")
  print(sapply(sparseGRM_new, class))
  
  return(objNull)
}

# check the control list in null model fitting for SPACox method
checkControl.NullModel.SPAmixPlus = function(control, traitType, ...)
{
  if(!traitType %in% c("time-to-event", "Residual"))
    stop("For 'SPAmixPlus' method, only traitType of 'time-to-event' or 'Residual' is supported.")
  
  if(is.null(control$PC_columns))
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') is required for 'SPAmixPlus' method.")
  
  if(length(control$PC_columns) != 1)
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') should be a character, not a character vector.")
  
  control$PC_columns = unlist(strsplit(control$PC_columns, split=","))
  if(length(control$PC_columns) == 1 && !grepl(",", control$PC_columns))
    warning("We detected that only one PC column exsit. Note that control$PC_columns should be a character splitted using ','. e.g. 'PC1,PC2'")
  
  return(control)
}
