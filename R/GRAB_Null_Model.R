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
#' @param ... Other arguments passed to function XXXX(). 
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
  
  # check missing data in formula, data, and subjData
  objFormula = handleFormula(formula, data, subset, subjData)  # check formula.R for more details
  response = objFormula$response
  designMat = objFormula$designMat   # note that the intercept column is not included
  subjData = as.character(objFormula$subjData)
  
  nData = length(subjData)
  cat("Number of subjects in 'formula':\t", nData,"\n")
  
  if(any(duplicated(subjData))) 
    stop("Duplicated subject IDs in 'formula' and 'data', i.e., 'subjData', is not supported!")
  
  OptionGRM = "none"
  
  if(!missing(GenoFile) & !missing(SparseGRMFile))
    stop("If 'DenseGRM' is used, please specify 'GenoFile', if 'SparseGRM' is used, please specify 'SparseGRMFile'. Cannot specify both files.")
  
  if(!missing(GenoFile))
    OptionGRM = "Dense"
  
  if(!missing(SparseGRMFile))
    OptionGRM = "Sparse"
 
  IfDenseGRM = IfSparseGRM = F 
  
    
    {
    # DenseGRM (FullGRM) if GenoFile is given
    genoList = setGenoInput(GenoFile, GenoFileIndex, subjData)   # check Geno.R for more details
    subjGeno = genoList$SampleIDs      # subjGeno should be the same as subjData
    if(genoList$genoType != "PLINK")
      stop("If DenseGRM is used when fitting a null model, then only Plink file is supported.")
  }else{
    # SparseGRM if missing(GenoFile)
    
  }
  
  if(method == "POLMM"){
    objNull = fitNullModel.POLMM(response, designMat, subjData, control)
  }
  
  # # SPACox method
  # if(method == "SPACox")
  #   objNull = fitNullModel.SPACox(formula, data, subset, subjData, subjGeno, control, ...);
  # 
  # # POLMM method
  # if(method == "POLMM")
  #   objNull = fitNullModel.POLMM(formula, data, subjData, subjGeno, control);
  
  objNull = list()
  objNull$subjGeno = subjGeno
  objNull$subjData = subjData
  
  objNull$Call = Call;
  
  return(objNull)
}









