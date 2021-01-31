
handleFormula = function (formula, data, subset)  {
  
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data", "subset"), 
             table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  
  mt = attr(x = mf, which = "terms")
  x <- model.matrix(object = mt, data = mf)
  
  response = model.response(mf)
  pos = as.numeric(names(response))
  designMat = model.matrix(object = mt, data = mf)
  
  if(colnames(designMat)[1] == "(Intercept)")
    designMat = designMat[,-1,drop=F]
  
  return(list(response = response,
              pos = pos,
              designMat = designMat))
}

# n = 20
# time = runif(n)
# event = rbinom(n, 1, 0.5)
# x1 = rnorm(n)
# x2 = rnorm(n)
# x3 = rbinom(n, 2, 0.5)
# 
# handleFormula(event~x1*x2+as.factor(x3), subset = x2>=(-1))
# handleFormula(time~x1*x2+as.factor(x3), subset = x2>=(-1))
# out = handleFormula(survival::Surv(time,event)~x1*x2+as.factor(x3), subset = x2>=(-1))
# 
# survival::coxph(out$response ~ out$designMat)







