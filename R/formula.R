
#' handle a formula (used in GRAB.NullModel function)
#'
#' handle a formula (used in GRAB.NullModel function), this function can help users better understand the input of GRAB.NullModel() function
#' @param formula a formula object, with the response on the left of a ~ operator and the covariates on the right. Do not add a column of intercept (e.g. a vector of ones) on the right. Missing values should be denoted by NA and the corresponding samples will be removed from analysis.
#' @param data a data.frame in which to interpret the variables named in the formula, or in the subset argument. Check ?model.frame for more details.
#' @param subset a specification of the rows to be used: defaults to all rows. This can be any valid indexing vector for the rows of data or if that is not supplied, a data frame made up of the variables used in formula. Check ?model.frame for more details.
#' @param subjData a character vector of subject IDs. Its order should be the same as the subjects order in the formula and data. 
#' @return an R list with elements of 'response', 'designMat', and 'subjData'.
#' @examples
#' n = 20
#' subjData = paste0("ID-",1:n)
#' pheno = rbinom(n, 1, 0.5)
#' x1 = rnorm(n)
#' x2 = rnorm(n)
#' x3 = rbinom(n, 2, 0.5)
#' objFormula = handleFormula(pheno ~ x1+x2*x3, subset = x2>0, subjData = subjData)
#' objFormula
#' @export
handleFormula = function(formula, data, subset, subjData)  {
  
  cl <- match.call()
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
  
  return(list(response = response,
              designMat = designMat,          
              subjData = subjData))
}


# 
# handleFormula(event~x1*x2+as.factor(x3), subset = x2>=(-1))
# handleFormula(time~x1*x2+as.factor(x3), subset = x2>=(-1))
# out = handleFormula(survival::Surv(time,event)~x1*x2+as.factor(x3), subset = x2>=(-1))
# 
# survival::coxph(out$response ~ out$designMat)







