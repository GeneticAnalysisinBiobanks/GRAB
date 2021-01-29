#' Fits a NULL model
#'
#' Fits a null model
#' @param formula a formula object, with the response on the left of a ~ operator and the covariates on the right. 
#' @param data a data.frame in which to interpret the variables named in the formula, or in the sugset and the weights argument (to be added later).
#' @param subjData a character vector of subject IDs. Its order should be the same as the subjects order in the formula and data.
#' @param method an R character, "SPACox", "SPAGE", "SAIGE", "POLMM", or "GATE"
#' @param trait.type an R character, "binary", "categorical", "quantitative", or "time-to-event"
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". If Null (default), the same prefix as GenoFile is used. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param control a list of parameters for controlling the GRAB.NullModel(). 
#' @param ... Other arguments passed to function XXXX(). 
#' @return an R object with a class of "XXXXX_NULL_Model" in which XXXXX is the 'method' used in analysis.
#' @examples
#' We put examples to the specific help pages for different methods. 
#' For example, if you want to use "SPACox" method, please check ?GRAB.SPACox for more details. 
#' @export
#' @import survival
GRAB.NullModel = function(formula,
                          data = NULL,
                          subjData,
                          method = "SPACox",
                          trait.type = "time-to-event",  # "binary", "categorical", "quantitative", "time-to-event"
                          GenoFile,
                          GenoFileIndex = NULL,
                          control = NULL,
                          ...)
{
  if(missing(subjData))
    stop("Argument 'subjData' is required to specify the subjects IDs in 'formula' and/or 'data'.")
  
  Call = match.call()

  control = checkControl.NullModel(control, method)
  
  objGeno = setGenoInput(GenoFile, GenoFileIndex)
  subjGeno = objGeno$samples      # subject IDs in genotype files
  
  objNull = fitNullModel(formula, data, subjData, subjGeno, method, trait.type, control, ...)
  objNull$Call = Call;
  
  return(objNull)
}

checkControl.NullModel = function(control, method)
{
  # check if control is an R list
  if(!is.null(control))
    if(class(control) != "list")
      stop("If specified, the argument of 'control' should be an R 'list'.")
  
  ######### SPACox method
  if(method == "SPACox"){
    control = checkControl.NullModel.SPACox(control)
  }
    
  
  ######### the below is for other methods
  
  print("The below is the list of control parameters used in null model fitting.")
  print(control)
  
  return(control)
}

fitNullModel = function(formula, data, subjData, subjGeno, method, trait.type, control, ...)
{
  # check the overlap of subjects in phenotype and genotype
  subjData = as.character(subjData)
  if(any(duplicated(subjData))) 
    stop("Duplicated subject IDs in 'formula' and 'data', i.e., 'subjData', is not supported!")
  if(any(duplicated(subjGeno))) 
    stop("Duplicated subject IDs in 'GenoFile' and 'GenoFileIndex' is not supported!")
  
  nData = length(subjData)
  nGeno = length(subjGeno)
  nBoth = length(intersect(subjData, subjGeno))
  
  # Update it to remove the below limitation
  if(nBoth != nData) 
    stop("For the current version, all subjects in 'formula' and 'data' should be also in 'GenoFile'.")
  
  cat("Number of subjects in 'formula':\t", nData,"\n")
  cat("Number of subjects in 'GenoFile':\t", nGeno,"\n")
  cat("Number of subjects in both 'formula' and 'GenoFile':\t", nBoth,"\n")
  
  ######### SPACox method
  if(method == "SPACox")
    obj.NullModel = fitNullModel.SPACox(formula, data, subjData, subjGeno, control, ...);
  
  ########
  
  return(obj.NullModel)
}




