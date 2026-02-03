GRAB.SPAmixPlus <- function() {
  .message("?GRAB.SPAmixPlus for instructions")
}


checkControl.NullModel.SPAmixPlus <- function(traitType, GenoFile, SparseGRMFile, control) {

  if (!traitType %in% c("time-to-event", "Residual")) {
    stop("For 'SPAmixPlus' method, only traitType of 'time-to-event' or 'Residual' is supported.")
  }

  if (!is.null(GenoFile)) {
    warning("Argument 'GenoFile' is ignored for method 'SPAmixPlus'.")
  }

  if (!is.null(SparseGRMFile)) {
    warning("Argument 'SparseGRMFile' is ignored for method 'SPAmixPlus'.")
  }

  default.control <- list(
    OutlierRatio = 1.5
  )
  control <- updateControl(control, default.control)

  if (is.null(control$PC_columns)) {
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') is required for 'SPAmixPlus' method.")
  }

  if (length(control$PC_columns) != 1) {
    stop("control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') should be a character, not a character vector.")
  }

  control$PC_columns <- unlist(strsplit(control$PC_columns, split = ","))
  if (length(control$PC_columns) == 1) {
    warning("We detected that only one PC column exists, is that what you want? ",
            "Note that control$PC_columns (e.g. 'PC1,PC2,PC3,PC4') should be ",
            "a character string split using ','.")
  }

  return(list(control = control, optionGRM = NULL))
}



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

makeGroup = function(yVec)
{
  # yVec is categorical data
  m1 = length(unique(yVec))
  if(m1 <= 10)
    Group = as.numeric(as.factor(yVec))-1  # from 0 to (m1-1)
  
  # yVec is quantitative data
  if(length(unique(yVec)) > 10)
    Group = floor((rank(yVec, ties.method = "max")-1) / length(yVec) * 10)  # from 0 to 9
  
  return(Group)
}

SPAmixPlus.Marker = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile,
                       OutputFileIndex = NULL,   
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);         # this function is in 'Framework.R'
  
  cat(paste0("NullModelClass is ", NullModelClass, "\n")) 
  
  if(is.null(OutputFileIndex)) 
    OutputFileIndex = paste0(OutputFile, ".index")
  
  # check the setting of control
  checkControl.ReadGeno(control)
  control = checkControl.Marker(control, NullModelClass)  
  nMarkersEachChunk = control$nMarkersEachChunk;
  outList = checkOutputFile(OutputFile, OutputFileIndex, "Marker", format(nMarkersEachChunk, scientific=F))    # this function is in 'Framework.R'
  
  # method = "SPAmixPlus" derived from "SPAmixPlus_NULL_Model"
  method = gsub("_NULL_Model", "", NullModelClass)

  indexChunk = outList$indexChunk
  Start = outList$Start
  End = outList$End
  
  if(End)
  {
    message = paste0("The analysis has been completed in earlier analysis. Results have been saved in '", OutputFile, "'. ",
                     "If you want to change parameters and restart the analysis, please use another 'OutputFile'.")
    return(message)
  }
  
  if(!Start){
    message = paste0("We detected that parts of analysis have been conducted from file:\t",
                     OutputFileIndex,"\n",
                     "We restart the analysis from chunk:\t",indexChunk+1,"\n");
    cat(message)
  }
  
  subjData = as.character(objNull$subjData);
  
  Group = makeGroup(objNull$yVec)  # this function is in Framework.R
  ifOutGroup = any(c("AltFreqInGroup", "AltCountsInGroup") %in% control$outputColumns)
  
  # Set global variables in C++ for filtering (MAF, MAC, etc.) and grouping
  # Implemented based on GRAB package logic to ensure consistency
  setMarker_GlobalVarsInCPP(control$impute_method,
                            control$missing_cutoff,
                            control$min_maf_marker,
                            control$min_mac_marker,
                            control$omp_num_threads,
                            Group,
                            ifOutGroup,
                            length(unique(Group)))
  
  ## set up an object for genotype
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs = subjData, control = control)  # this function is in 'Geno.R'
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  CHROM = markerInfo$CHROM
  genoIndex = markerInfo$genoIndex
  
  # all markers were split into multiple chunks
  genoIndexList = splitMarker(genoIndex, nMarkersEachChunk, CHROM);  
  
  nChunks = length(genoIndexList)
  
  cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
  cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
  cat("Number of chunks for all markers:\t", nChunks, "\n")
  
  chrom = "InitialChunk"
  for(i in (indexChunk+1):nChunks)
  {
    tempList = genoIndexList[[i]]
    genoIndex = tempList$genoIndex
    tempChrom = tempList$chrom
    
    # Update Null Model Object (C++) if necessary. 
    # For SPAmixPlus, we call setMarker.SPAmixPlus every time, or at least when we start.
    # The original logic calls it every chunk if chrom changes, or maybe every chunk?
    # In GRAB.Marker:
    # if(chrom != tempChrom){ ... setMarker... }
    # Let's just call it once at the beginning if chrom changes, or always?
    # SPAmixPlus doesn't depend on Chromosome. But calling it is cheap (just setting pointer).
    # But wait, setSPAmixPlusobjInCPP involves passing large vectors? 
    # No, it passes them. Rcpp might copy them.
    # So we should only call it once at the start.
    
    if(i == (indexChunk+1)){
        setMarker.SPAmixPlus(objNull, control)
    }
    
    cat(paste0("Analysis of chunk ", i, "/", nChunks, " starts (Chrom: ", tempChrom, ")...\n"))
    
    # mainMarker function
    # call mainMarker.SPAmixPlus
    obj.mainMarker = mainMarker.SPAmixPlus(genoType, genoIndex, control$outputColumns, objNull, control)
    
    # if OutputFile includes other columns customized by users
    if(ifOutGroup){
      # To be implemented if needed. For now assume basic columns.
      # The original GRAB.Marker implemented this.
      # But mainMarker.SPAmixPlus returns what it returns.
    }
    
    # write the output
    writeOutputFile(list(obj.mainMarker), OutputFile, OutputFileIndex, "Marker", nMarkersEachChunk, i, i==1, i==nChunks)
    
    cat(paste0("Analysis of chunk ", i, "/", nChunks, " completed.\n"))
  }
  
  closeGenoInputInCPP(genoType)
  
  message = paste0("The analysis has been completed. Results have been saved in '", OutputFile, "'.")
  return(message)
}


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

