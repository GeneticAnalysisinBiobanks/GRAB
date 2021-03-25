#' Fits a NULL model
#'
#' Fits a null model
#' @param formula a formula object, with the response on the left of a ~ operator and the covariates on the right. Do not add a column of intercept (e.g. a vector of ones) on the right. Missing values should be denoted by NA and the corresponding samples will be removed from analysis.
#' @param data a data.frame in which to interpret the variables named in the formula, or in the subset argument. Check ?model.frame for more details.
#' @param subset a specification of the rows to be used: defaults to all rows. This can be any valid indexing vector for the rows of data or if that is not supplied, a data frame made up of the variables used in formula. Check ?model.frame for more details.
#' @param subjData a character vector of subject IDs. Its order should be the same as the subjects order in the formula and data. 
#' @param method an R character, "SPACox", "SPAGE", "SAIGE", "POLMM", or "GATE"
#' @param traitType an R character, "binary", "ordinal", "quantitative", or "time-to-event"
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". If Null (default), the same prefix as GenoFile is used. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param control a list of parameters for controlling the GRAB.NullModel(). 
#' @param ... Other arguments passed to or from other methods(). 
#' @return an R object with a class of "XXXXX_NULL_Model" in which XXXXX is the 'method' used in analysis.
#' @examples
#' ## Example using POLMM to analyze ordinal categorical data
#' POLMM_data = read.csv(system.file("extdata", "POLMM_data.csv", package = "GRAB"))
#' formula  = outcome ~ Cova1 + Cova2
#' data = POLMM_data
#' GRAB.NullModel(formula, data)
#' # For example, if you want to use "SPACox" method, please check ?GRAB.SPACox for more details. 
#' @export
#' @import survival
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
  
  handleGRM(GenoFile, GenoFileIndex, SparseGRMFile, subjData)
  
  if(method == "POLMM"){
    objNull = fitNullModel.POLMM(response, designMat, subjData, control)
  }
  
  if(method == "SPACox")
    objNull = fitNullModel.SPACox(response, designMat, subjData, control, ...)
  
  objNull$subjData = subjData
  
  objNull$Call = Call;
  objNull$sessionInfo = sessionInfo()
  objNull$time = Sys.time()
  
  print(paste0("Complete null model fitting: ", objNull$time))
  
  return(objNull)
}


## to be updated later
handleGRM = function(GenoFile, GenoFileIndex, SparseGRMFile, subjData)
{
  if(!missing(SparseGRMFile)){
    OptionGRM = "Sparse"
  }else{
    OptionGRM = "Dense"
  }
  
  if(!missing(GenoFile)){
    genoList = setGenoInput(GenoFile, GenoFileIndex, subjData)   # check Geno.R for more details
    subjGeno = genoList$SampleIDs      # subjGeno should be the same as subjData
    if(genoList$genoType != "PLINK" & OptionGRM == "Dense")
      stop("If DenseGRM is used when fitting a null model, then only Plink file is supported.")
  }
}







