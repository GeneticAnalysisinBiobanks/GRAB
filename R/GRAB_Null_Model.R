#' GRAB: NULL model fitting
#'
#' GRAB package uses score test for GWAS: in step 1, we fit a null model (check \code{?GRAB.NullModel}) including response variable, covariates, and GRM (if needed). In step 2, we perform score test for marker-level analysis (check \code{?GRAB.Marker}) and region-level analysis (check \code{?GRAB.Region}).
#' 
#' @param formula a formula object, with the response on the left of a ~ operator and the covariates on the right. Do not add a column of intercept (e.g. a vector of ones) on the right. Missing values should be denoted by NA and the corresponding samples will be removed from analysis.
#' @param data a data.frame in which to interpret the variables named in the formula, or in the subset argument. Check \code{?model.frame} for more details.
#' @param subset a specification of the rows to be used: defaults to all rows. This can be any valid indexing vector for the rows of data or if that is not supplied, a data frame made up of the variables used in formula. Check \code{?model.frame} for more details.
#' @param subjData a character vector of subject IDs. Its order should be the same as the subject order in the formula and data. 
#' @param method a character: "SPACox", "SPAGE", "SAIGE", "POLMM", or "GATE"
#' @param traitType a character: "binary", "ordinal", "quantitative", or "time-to-event"
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK, BGEN, and VCF. More details are in \cdoe{?GRAB.ReadGeno}.
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". If Null (default), the same prefix as GenoFile is used. More details are in \cdoe{?GRAB.ReadGeno}.
#' @param SparseGRMFile a character of sparseGRM file. An example is \code{system.file("SparseGRM","SparseGRM.txt",package="GRAB")}
#' @param control a list of parameters for controlling the \code{GRAB.NullModel()}. For more details, please check \code{?GRAB.control}. 
#' @param ... other arguments passed to or from other methods. 
#' @return an R object with a class of "XXXXX_NULL_Model" in which XXXXX is the 'method' used in analysis. The object will be used in \code{GRAB.Marker()} and \code{GRAB.Region()}. Functions of \code{save()} and \code{load()} can be used to save and load the object for future usage.
#' @examples
#' # Examples have been put in the specific help pages for specific methods. 
#' # If you want to use "SPACox" method, please check ?GRAB.SPACox for more details.
#' # If you want to use "POLMM" method, please check ?GRAB.POLMM for more details.
#' @export
#' @import survival, data.table
GRAB.NullModel = function(formula,
                          data = NULL,
                          subset = NULL,
                          subjData,
                          method = "SPACox",
                          traitType = "time-to-event",  # "binary", "ordinal", "quantitative", "time-to-event"
                          GenoFile,
                          GenoFileIndex = NULL,
                          SparseGRMFile,
                          control = NULL,
                          ...)
{
  if(missing(subjData))
    stop("Argument 'subjData' is required to specify the subjects IDs in 'formula' and/or 'data'.")
  
  Call = match.call()
  
  # The following function is in 'control.R'
  control = checkControl.NullModel(control, method, traitType)
  
  #### START: formula.R
  #### input: formula, data, subset, and subjData
  #### output: response, designMat, subjData
  
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(x = c("formula", "data", "subset", "subjData"), 
             table = names(mf), nomatch = 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  
  mt = attr(x = mf, which = "terms")
  
  response = model.response(mf)
  designMat = model.matrix(object = mt, data = mf)
  subjData = model.extract(mf, "subjData")
  
  if(colnames(designMat)[1] == "(Intercept)")
    designMat = designMat[,-1,drop=F]
  
  nData = length(subjData)
  cat("Number of subjects in 'formula':\t", nData,"\n")
  
  if(any(duplicated(subjData))) 
    stop("Duplicated subject IDs in 'formula' and 'data', i.e., 'subjData', is not supported!")
  
  #### END: formula.R
  
  optionGRM = handleGRM(GenoFile, GenoFileIndex, SparseGRMFile, subjData)
  
  if(method == "POLMM"){
    # The following function is in 'POLMM.R'
    objNull = fitNullModel.POLMM(response, designMat, subjData, control, optionGRM)
  }
  
  if(method == "SPACox")
    objNull = fitNullModel.SPACox(response, designMat, subjData, control, ...)
  
  objNull$subjData = subjData
  
  objNull$Call = Call;
  objNull$sessionInfo = sessionInfo()
  objNull$time = Sys.time()
  objNull$control = control
  
  print(paste0("Complete null model fitting: ", objNull$time))
  
  return(objNull)
}


## to be updated later
handleGRM = function(GenoFile, GenoFileIndex, SparseGRMFile, subjData)
{
  if(!missing(SparseGRMFile)){
    print("Sparse GRM is used when fitting a null model.")
    SparseGRM = data.table::fread(SparseGRMFile)
    SparseGRM = as.data.frame(SparseGRM)
    
    KinMatListR = updateSparseGRM(SparseGRM, subjData)
    
    # the following function is in Main.cpp
    setSparseGRMInCPP(KinMatListR)
    optionGRM = "SparseGRM"
  }else{
    print("Dense GRM is used when fitting a null model.")
    genoList = setGenoInput(GenoFile, GenoFileIndex, subjData)   # check Geno.R for more details
    subjGeno = genoList$SampleIDs      # subjGeno should be the same as subjData
    if(genoList$genoType != "PLINK")
      stop("If DenseGRM is used when fitting a null model, then only Plink file is supported.")
    
    memoryChunk = 2 # (GB)
    minMafGRM = 0.01
    maxMissingGRM = 0.1
    
    # the following function is in Main.cpp
    setDenseGRMInCPP(memoryChunk, minMafGRM, maxMissingGRM)
    optionGRM = "DenseGRM"
  }
  
  return(optionGRM)
}







