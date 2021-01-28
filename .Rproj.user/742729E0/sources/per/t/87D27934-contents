#' Fits a NULL model for SPACox
#'
#' Fits a null Cox proportional hazards model and then calculates the empirical cumulant generation function (CGF) of the martingale residuals
#' @param formula a formula to be passed to function coxph(). For more details, please refer to package survival.
#' @param data a data.frame in which to interpret the variables named in the formula
#' @param subjData a character vector of subject IDs. NOTE: its order should be the same as the subjects order in the formula.
#' @param method an R character, "SPACox", "SPAGE", "SAIGE", "POLMM", "GATE"
#' @param trait.type an R character, "binary", "categorical", "quantitative", "time-to-event"
#' @param GenoFile xxxxxxxxxxx
#' @param GenoFileIndex xxxxxxxxxxx
#' @param control xxxxxxxxxxx
#' @param range a two-element numeric vector (default: c(-100,100)) to specify the domain of the empirical CGF.
#' @param length.out a positive integer (default: 9999) for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.
#' @param ... Other arguments passed to function coxph(). For more details, please refer to package survival.
#' @return an object with a class of "SPACox_NULL_Model".
#' @examples
#' # Simulation phenotype and genotype
#' GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
#' N = 100
#' Pheno = data.frame(ID = paste0("f",1:N,"_1"),
#'                    event=rbinom(N,1,0.5),
#'                    time=runif(N),
#'                    Cov1=rnorm(N),
#'                    Cov2=rbinom(N,1,0.5))
#' obj.SPACox = GRAB.NullModel(survival::Surv(time,event)~Cov1+Cov2, 
#'                             data=Pheno, subjData = Pheno$ID, method = "SPACox", GenoFile = GenoFile)
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
  
  obj.NullModel = fitNullModel(formula, data, subjData, subjGeno, method, trait.type, control, ...)
  # obj.NullModel$Call = Call;
  
  return(obj.NullModel)
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




