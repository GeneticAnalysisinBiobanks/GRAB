
######## 2021-03-23: this file is not used in GRAB package

handleFormula = function(formula, data, subset, subjData)  {
  
  print(subset)
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

# n = 20
# subjData = paste0("abc",1:(n+1))
# time = runif(n)
# event = rbinom(n, 1, 0.5)
# x1 = rnorm(n)
# x2 = rnorm(n)
# x3 = rbinom(n, 2, 0.5)
# 
# a = handleFormula(event ~ x1+x2*x3, subset = x2>0, subjData = subjData)

# 
# handleFormula(event~x1*x2+as.factor(x3), subset = x2>=(-1))
# handleFormula(time~x1*x2+as.factor(x3), subset = x2>=(-1))
# out = handleFormula(survival::Surv(time,event)~x1*x2+as.factor(x3), subset = x2>=(-1))
# 
# survival::coxph(out$response ~ out$designMat)







