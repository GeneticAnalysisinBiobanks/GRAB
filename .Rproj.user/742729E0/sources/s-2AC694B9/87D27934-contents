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
#' Please check help(SPACox) for a simulated example.
#' @export
#' @import survival
GRAB.NullModel = function(formula,
                          data = NULL,
                          subjData,
                          method = "SPACox",
                          trait.type = "time-to-event",  # "binary", "categorical", "quantitative", "time-to-event"
                          GenoFile,
                          GeniFileIndex = NULL,
                          control = NULL,
                          range = c(-100,100),
                          length.out = 10000,
                          ...)
{
  if(missing(subjData))
    stop("Argument 'subjData' is required to specify the subjects IDs in 'formula' and/or 'data'.")
  
  Call = match.call()

  control = checkControl.NullModel(control, method)
  
  objGeno = setGenoInput(GenoFile, GenoFileIndex)
  samples = objGeno$samples      # subject IDs in genotype files
  
  obj.NullModel = fitNullModel(formula, data, subjData, samples, method, trait.type, control)
  
  return(obj.NullModel)
}

checkControl.NullModel = function(control, method)
{
  if(method == "SPACox")
    control = checkControl.NullModel.SPACox(control)
  
  #########
  
  return(control)
}

fitNullModel(formula, data, subjData, samples, method, trait.type, control)
{
  if(method == "SPACox")
    obj.NullModel = fitNullModel.SPACox(formula, data, subjData, samples, control);
  
  ########
  
  return(obj.NullModel)
}




