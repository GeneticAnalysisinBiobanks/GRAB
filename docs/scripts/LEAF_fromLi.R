args=commandArgs(TRUE)
print(args)
print(sessionInfo())
n.cpu=as.numeric(args)
##########################################

library(data.table)
library(dplyr)
source("Summix.R")


# function: batch effect testing########################################################
Batcheffect.Test = function(n0,                      # number of controls
                            n1,                      # number of cases
                            n.ext,                   # number of external dataset
                            maf0,                    # estimated MAF in controls
                            maf1,                    # estimated MAF in cases
                            maf.ext,                 # estimated MAF in external dataset
                            pop.prev,
                            var.ratio=1)
{
  er = n1/(n1+n0)
  w0 = (1-pop.prev)/pop.prev/((1-er)/er)
  w1 = 1
  
  weight.maf = sum(maf0*w0*n0 + maf1*w1*n1) / sum(w0*n0 + w1*n1)                           ## weighted mean of genotypes
  est.maf = sum(maf0*w0*n0 + maf1*w1*n1 + maf.ext*n.ext*w0) / sum(n1*w1+n0*w0+ n.ext*w0)   ## MAF estimates
  
  v =( (n1*w1^2+n0*w0^2)/(2*(n1*w1+n0*w0)^2)  + 1/(2*n.ext) ) * est.maf * (1-est.maf) ## variance of test statistics
  z = (weight.maf - maf.ext) / sqrt(v)  ## standardized statistics
  z.adj = z/sqrt(var.ratio)             ## adjusted statistics by variance ratio
  p = 2*pnorm(-abs(z.adj), lower.tail=TRUE)
  
  return(p)
}

#estimate TPR and sigma2################################################################
# Reparameterized optimization:
#   TPR    ∈ (0, 1)   ←  σ(η_T) = 1 / (1 + exp(−η_T))
#   sigma2 ∈ (0, ∞)   ←  exp(η_S)
# Unbounded Nelder-Mead on (η_T, η_S) ∈ R²; no post-hoc clamp and no σ²
# upper bound.  Boundaries are reached only asymptotically as η → ±∞.
fun.est.param = function(vec_p_bat,
                         vec_var_Sbat,
                         vec_cutoff=seq(0.01,0.4,0.1)
){
  
  ########  step1: the proportion of p_bat>cutoff
  vec_p_deno=lapply(vec_cutoff,function(p_cut){
    p_deno=mean(na.omit(vec_p_bat>p_cut))
  })%>%unlist()
  
  ######## optimization function — argument is η = (η_T, η_S) ∈ R²
  opti_fun_eta = function( eta, var_Sbat, vec_p_deno ){
    TPR    = 1 / (1 + exp(-eta[1]))
    sigma2 = exp(eta[2])
    diff=lapply(1: length(vec_cutoff),function(j){
      p_cut=vec_cutoff[j]
      lb = -qnorm(1-p_cut/2) * sqrt(var_Sbat)
      ub = qnorm(1-p_cut/2) *sqrt(var_Sbat)
      p_deno=vec_p_deno[j]
      
      c = pnorm(ub, 0, sqrt(var_Sbat+sigma2), log.p = T)
      d = pnorm(lb, 0, sqrt(var_Sbat+sigma2), log.p = T)
      
      pro.cut = TPR*( exp(d) * (exp(c-d) - 1) )+(1-TPR)*(1-p_cut)
      t = ((p_deno - pro.cut)/p_deno)^2
    })%>%do.call("sum",.)
    
    return(diff)
    
  }
  
  #######estimate TPR and sigma2 for each SNP
  ## Nelder-Mead on (η_T, η_S).  Reason for choosing NM over BFGS:
  ## when TPR is near its upper limit, σ² becomes asymptotically
  ## unidentifiable — the objective is essentially flat along the
  ## η_S axis.  BFGS's inverse-Hessian estimate degenerates in such
  ## flat regions and the line search can drift η_S far from the
  ## starting point (we observed η_S ≈ +20, i.e., σ² ≈ 10⁹).  NM
  ## constrains motion to a bounded simplex and remains numerically
  ## stable in degenerate regions.
  var.diff = lapply(1:length(vec_var_Sbat), function(i){
    
    if(i%%100==0)cat(i,"\n")
    
    obj = optim(par     = c(0, log(0.01)),
                fn      = opti_fun_eta,
                vec_p_deno = vec_p_deno,
                var_Sbat   = vec_var_Sbat[i],
                method  = "Nelder-Mead",
                control = list(reltol = 1e-10, maxit = 500))
    TPR    = 1 / (1 + exp(-obj$par[1]))
    sigma2 = exp(obj$par[2])
    return(cbind(TPR, sigma2))
    
  })%>%
    do.call("rbind",.)%>%as_tibble()
  
  return(var.diff)
}


#optimal weight for mu_ext ##########################################################
fun.optimalWeight = function(par, pop.prev, R, y, mu1, w , mu, N, n.ext, sigma2, TPR){
  b=par[1]
  
  p.fun= function(b,pop.prev, R, y, mu1, mu, w, N, n.ext, sigma2, TPR){
    
    meanR= mean(R)
    sumR = sum(R)
    
    mu0 = mu
    mu.pop = mu1*pop.prev+mu0*(1-pop.prev)
    
    mu.i = ifelse(y==1, 2*mu1, 2*mu0)
    
    S = sum((R-(1-b)*meanR)*mu.i)-sumR*2*b*mu.pop
    
    w1 = w/(2*sum(w))
    mu = mean(mu.i)/2
    
    var_mu_ext = mu*(1-mu)/(2*n.ext)
    var_Sbat = sum(w1^2)*2*mu*(1-mu) + var_mu_ext
    
    p_cut=0.1
    lb = -qnorm(1-p_cut/2) * sqrt(var_Sbat)
    ub = qnorm(1-p_cut/2) *sqrt(var_Sbat)
    c = pnorm(ub, 0, sqrt(var_Sbat+sigma2),log.p = T)
    d = pnorm(lb, 0, sqrt(var_Sbat+sigma2),log.p = T)
    p_deno = TPR*( exp(d) * (exp(c-d) - 1) )+(1-TPR)*(1-p_cut)
    
    ##sigma2=0
    var.int = sum((R-(1-b)*meanR)^2)*2*mu*(1-mu)
    var_S = var.int +  4*b^2 *sumR^2* var_mu_ext
    cov_Sbat_S = sum(w1*(R-(1-b)*meanR))*2*mu*(1-mu)+2*b*sumR*var_mu_ext
    VAR = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat),nrow=2)
    p0 = max(0,pmvnorm(lower=c( -Inf, lb),upper=c( -abs(S), ub), mean=c(0,0), sigma=VAR))
    
    ##sigma2!=0
    var_S1 = var.int +  4*b^2 *sumR^2* (var_mu_ext+sigma2)
    cov_Sbat_S1 = sum(w1*(R-(1-b)*meanR))*2*mu*(1-mu)+2*b*sumR*(var_mu_ext+sigma2)
    var_Sbat1 = var_Sbat+sigma2
    VAR1 = matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1),nrow=2)
    p1=max(0,pmvnorm(lower=c( -Inf, lb),upper=c( -abs(S), ub), mean=c(0,0), sigma=VAR1))
    
    p.con = 2*(TPR*p1+(1-TPR)*p0)/p_deno
    #diff = -log10(p.con)+log10(5e-8)
    diff = -log10(p.con/5e-8)
    
    return(diff)
    
  }
  
  mu1=uniroot(p.fun, lower=mu,upper=1
              
              , b=b, pop.prev=pop.prev,mu=mu,
              R=R, y=y, w=w,N=N, n.ext=n.ext,sigma2=sigma2,TPR=TPR)$root
  
  return(mu1)
  
}
# SPA method ###################################################################
WtCoxG.test = function(g,
                       R,
                       w,
                       p_bat,
                       TPR=NA,
                       sigma2=NA,
                       b=0,
                       var.ratio.int = 1,          ### variance ratio of S.int
                       var.ratio.w0 = 1,           ### variance ratio of S.bat when bathceffect = 0
                       var.ratio.w1 = 1,           ### variance ratio of S.bat when bathceffect = 1
                       var.ratio0 = 1,             ### variance ratio of Score when bathceffect = 0
                       var.ratio1 = 1,             ### variance ratio of Score when bathceffect = 1
                       mu.ext=NA,
                       n.ext=NA,
                       p_cut=0.1                   ### batcheffect p value cut off
){
  
  ##imputation missing SNP
  missing.rate = mean(is.na(g))
  pos.na = which(is.na(g))
  if(missing.rate != 0){
    g[pos.na] = mean(na.omit(g))
  }
  
  
  ### if external MAF is unavailable,  internal only
  if(is.na(mu.ext)){
    
    mu.int = mean(g)/2
    p.con = SPA_G.one.SNP_homo(g=g, R=R, mu.ext=NA, n.ext= 0,sigma2=0, var.ratio = var.ratio.int)[1]
    p_deno =NA
    S = sum(R*(g - mean(g)))
    return(cbind( p.con, S))
    
  }
  
  if(is.na(p_bat)|sum(g)<10|sum(2-g)<10 ){
    
    p.con = S =NA
    return(cbind( p.con, S))
    
  }else if(p_bat<p_cut){  # batch effect p<0.1, internal only
    R_tilde = R - mean(R)
    mu=mu.int = mean(g)/2
    S = sum(R *(g - 2*mu.int))
    w1 = w/(2*sum(w))
    var_mu_ext = mu*(1-mu)/(2*n.ext)
    var_Sbat = sum(w1^2)*2*mu*(1-mu) + var_mu_ext
    
    lb = -qnorm(1-p_cut/2) * sqrt(var_Sbat)*sqrt(var.ratio.w0)
    ub=-lb
    c = pnorm(lb/sqrt(var.ratio.w1), 0, sqrt(var_Sbat+sigma2), log.p = T)
    p_deno = TPR * 2 * exp(c) + (1-TPR) * p_cut
    
    ##sigma2=0
    p_spa_s0 = SPA_G.one.SNP_homo(g=g, R=R, var.ratio = var.ratio0)[1]
    var_S = (S/sqrt(var.ratio0))^2/qchisq(p_spa_s0, 1, ncp = 0, lower.tail = F)
    
    var.int = sum(R_tilde^2) * 2 * mu * (1-mu)
    # var_S = var.int +  4*b^2 *sumR^2* var_mu_ext
    cov_Sbat_S = sum(w1*R_tilde)*2*mu*(1-mu)
    cov_Sbat_S = cov_Sbat_S*sqrt(var_S/var.int)
    VAR = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat),nrow=2)
    p0 =max(0,min(1, pmvnorm(lower=c( -Inf, -Inf), upper=c( -abs(S/sqrt(var.ratio0)), lb/sqrt(var.ratio.w0)),
                             mean=c(0,0), sigma=VAR)))+
      max(0,min(1, pmvnorm(lower=c( -Inf, ub/sqrt(var.ratio.w0)), upper=c( -abs(S/sqrt(var.ratio0)), Inf), mean=c(0,0), sigma=VAR)))
    ##sigma2!=0
    var_S1 = var_S
    cov_Sbat_S1 = cov_Sbat_S
    var_Sbat1 = var_Sbat + sigma2
    VAR1 = matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1),nrow=2)
    p1 = max(0,min(1, pmvnorm(lower=c( -Inf, -Inf), upper=c( -abs(S/sqrt(var.ratio1)), lb/sqrt(var.ratio.w1)), 
                              mean=c(0,0), sigma = VAR1)))+
      max(0,min(1, pmvnorm(lower = c( -Inf, ub/sqrt(var.ratio.w1)), upper = c( -abs(S/sqrt(var.ratio1)), Inf),
                           mean = c(0,0), sigma = VAR1)))
    
    p.con = 2*(TPR*p1+(1-TPR)*p0)/p_deno
    return(cbind( p.con, S) )
    
  }else{ # batch effect p>0.1
    
    meanR= mean(R)
    sumR = sum(R)
    mu.int = mean(g)/2
    N = length(g)
    
    
    mu = (1-b) * mu.int + b * mu.ext
    S= sum(R*(g-2*mu))
    # S=sum((R-(1-b)*meanR)*g)-sumR*2*b*mu.ext
    
    w1 = w/(2*sum(w))
    
    var_mu_ext = mu * (1 - mu)/(2 * n.ext)
    var_Sbat = sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext
    
    
    lb = -qnorm(1-p_cut/2) * sqrt(var_Sbat) * sqrt(var.ratio.w0)
    ub = qnorm(1-p_cut/2) *sqrt(var_Sbat) * sqrt(var.ratio.w0)
    c = pnorm(ub/sqrt(var.ratio.w1), 0, sqrt(var_Sbat+sigma2), log.p = T)
    d = pnorm(lb/sqrt(var.ratio.w1), 0, sqrt(var_Sbat+sigma2), log.p = T)
    p_deno = TPR*(  exp(d) * (exp(c-d) - 1) )+(1-TPR)*(1-p_cut)
    
    ##sigma2=0
    
    p_spa_s0 = SPA_G.one.SNP_homo(g=g, R=R, b=b ,mu.ext=mu.ext, n.ext= n.ext, sigma2=0,
                                  var.ratio = var.ratio0)[1]
    var_S = S^2 / var.ratio0 / qchisq(p_spa_s0, 1, ncp = 0, lower.tail = F)
    
    
    var.int = sum((R-(1-b) * meanR)^2) * 2 * mu * (1-mu)
    # var_S = var.int +  4*b^2 *sumR^2* var_mu_ext
    cov_Sbat_S = sum(w1 * (R-(1-b) * meanR)) * 2*mu * (1-mu) + 2 * b * sumR * var_mu_ext
    cov_Sbat_S = cov_Sbat_S*sqrt(var_S/(var.int +  4 * b^2 * sumR^2 * var_mu_ext))
    VAR = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat),nrow=2)
    p0 = max(0,min(1, pmvnorm(lower=c( -Inf, lb/sqrt(var.ratio.w0)),
                              upper=c( -abs(S/sqrt(var.ratio0)), ub/sqrt(var.ratio.w0)), mean=c(0,0), sigma=VAR)))
    
    ##sigma2!=0
    p_spa_s1 = SPA_G.one.SNP_homo(g=g, R=R, b=b, mu.ext=mu.ext, n.ext= n.ext, sigma2=sigma2,
                                  var.ratio = var.ratio1)[1]
    
    var_S1 = S^2/var.ratio1/qchisq(p_spa_s1, 1, ncp = 0, lower.tail = F)
    
    #var_S1 = var.int +  4*b^2 *sumR^2* (var_mu_ext+sigma2)
    cov_Sbat_S1 = sum(w1*(R-(1-b)*meanR))*2*mu*(1-mu)+2*b*sumR*(var_mu_ext+sigma2)
    cov_Sbat_S1 = cov_Sbat_S1*sqrt(var_S1/(var.int +  4*b^2 *sumR^2* (var_mu_ext+sigma2)))
    var_Sbat1 = var_Sbat+sigma2
    VAR1 = matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1),nrow=2)
    p1= max(0,min(1,pmvnorm(lower=c( -Inf, lb/sqrt(var.ratio.w1)),
                            upper=c( -abs(S/sqrt(var.ratio1)), ub/sqrt(var.ratio.w1)), mean=c(0,0), sigma=VAR1)))
    
    p.con = 2*(TPR*p1+(1-TPR)*p0)/p_deno
    
    return(cbind( p.con, S))
    
  }
  
  
}
#############################################################################
# hidden function----------------------------------------------------------------------
# MGF function --------------------------------------------------
library(mvtnorm)
M_G0 = function(t, MAF){
  re = (1 - MAF + MAF * exp(t))^2
  return(re)
}
# The first derivative of the MGF of G (genotype)
M_G1 = function(t, MAF){
  re = 2*(MAF * exp(t))*(1 - MAF + MAF * exp(t))
  return(re)
}

# The second derivative of the MGF of G (genotype)
M_G2 = function(t, MAF){
  re = 2*(MAF * exp(t))^2 + 2*(MAF * exp(t))*(1 - MAF + MAF * exp(t))
  return(re)
}

# The CGF of G (genotype)
K_G0 = function(t, MAF){
  re = log(M_G0(t, MAF))
  return(re)
}

K_G1 = function(t, MAF){
  re = M_G1(t, MAF)/M_G0(t, MAF)
  return(re)
}

# The second derivative of the CGF of G (genotype)
K_G2 = function(t, MAF){
  re = M_G0(t, MAF)/M_G0(t, MAF) * M_G2(t, MAF)/M_G0(t, MAF) -(M_G1(t, MAF)/M_G0(t, MAF))^2
  return(re)
}

# The CGF of score test statistic
H_org = function(t, R, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est, meanR, b){
  
  n.t = length(t)
  out = rep(0,n.t)
  R_hat  = sumR/N.all
  
  mu.adj = -2 * b * sumR * MAF
  var.adj = 4 * b^2 * sumR^2 * var_mu_ext
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(K_G0(t1 * (R - (1 - b) * meanR), MAF)) + mu.adj * t1 + var.adj/2 * t1^2
    
  }
  
  return(out)
}

# The first derivative of the CGF of score test statistic
H1_adj = function(t, R, s, MAF, n.ext, N.all, sumR,var_mu_ext, g.var.est,meanR,b)
{
  n.t = length(t)
  out = rep(0,n.t)
  R_hat  = sumR/N.all
  
  mu.adj = -2*b*sumR*MAF
  var.adj = 4*b^2 *sumR^2* var_mu_ext
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(( R - (1 - b) * meanR) *K_G1(t1 * ((R - (1 - b) * meanR)), MAF)) + mu.adj + var.adj * t1 - s
  }
  return(out)
}

# The second derivative of the CGF of score test statistic
H2 = function(t, R, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est, meanR, b)
{
  n.t = length(t)
  out = rep(0,n.t)
  R_hat  = sumR/N.all
  
  mu.adj = -n.ext*R_hat*2*MAF
  var.adj = n.ext*R_hat^2*2*MAF*(1-MAF)
  
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum((R - (1 - b) * meanR)^2 * K_G2(t1 * (R - (1 - b) * meanR) , MAF)) + var.adj
  }
  return(out)
}

GetProb_SPA_G = function(MAF, R, s, n.ext, N.all, sumR,var_mu_ext, g.var.est,meanR,b, lower.tail){
  
  out = uniroot(H1_adj, c(-1,1), extendInt = "yes",
                R=R, s=s, MAF=MAF, n.ext=n.ext, N.all = N.all, sumR = sumR,
                var_mu_ext = var_mu_ext, g.var.est=g.var.est, meanR =meanR, b=b)
  zeta = out$root
  
  k1 = H_org(zeta, R=R, MAF=MAF, n.ext=n.ext, N.all = N.all, sumR = sumR,
             var_mu_ext = var_mu_ext, g.var.est=g.var.est, meanR =meanR,b=b)
  k2 = H2(zeta, R=R, MAF=MAF, n.ext=n.ext, N.all = N.all, sumR = sumR,
          var_mu_ext = var_mu_ext, g.var.est=g.var.est, meanR =meanR, b=b)
  
  temp1 = zeta * s - k1
  
  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  
  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  # pval.norm = pnorm(q2, lower.tail = lower.tail)
  re = pval
  return(re)
}


#saddlepoint approximation (SPA) to calicrate the p value in the case that batcheffect = 0 or 1
SPA_G.one.SNP_homo = function(g,                ### genotype vector
                              R,                ### null model residual vector
                              mu.ext = NA,      ### external MAF
                              n.ext = NA,       ### external sample size
                              b = 0,            ### weight of external MAF
                              sigma2 = NA,      ###
                              var.ratio = 1,
                              Cutoff = 2,
                              impute.method = "fixed",
                              missing.cutoff = 0.15,
                              min.mac = 10,          # update on 2022-08-16 : replace 0.0001 by 0.000001
                              G.model = "Add")
{
  ## calculate MAF and update genotype vector
  
  ##imputation missing SNP
  missing.rate = mean(is.na(g))
  pos.na = which(is.na(g))
  if(missing.rate != 0){
    g[pos.na] = mean(na.omit(g))
  }
  
  
  if(is.na(mu.ext)){
    mu.ext =0
    n.ext=0
    
  }
  
  
  if(sum(g)<min.mac|sum(2-g)<min.mac|missing.rate>missing.cutoff){
    
    MAF= mean(na.omit(g))/2
    
    return(c(NA, NA))
  }
  
  
  ######################################################################
  
  ## Score statistic
  N=length(g)
  mu.int = mean(g)/2
  MAF = (1-b) * mu.int + b * mu.ext
  sumR = sum(R)
  N.all = N + n.ext
  S = sum(R *(g-2*MAF))
  S = S / var.ratio
  
  ## estimated variance without adjusting for covariates
  g.var.est = 2 * MAF * (1 - MAF)
  var_mu_ext = ifelse(n.ext==0,0, MAF*(1-MAF)/(2*n.ext)+sigma2)
  
  
  #  S.var = sum(R^2 * g.var.est + sum(R)^2 * g.var.est/n.ext)
  meanR  = mean(R)
  S.var = sum((R-(1-b)*meanR)^2)*g.var.est +  4*b^2 *sumR^2* var_mu_ext
  
  z = S/sqrt(S.var)
  
  if(abs(z) < Cutoff){
    pval.norm = pnorm(abs(z), lower.tail = FALSE)*2
    return(c(pval.norm, pval.norm))  # update on 2022-10-05 : MAF.est.negative.num
  }else{
    pval1 = GetProb_SPA_G(MAF, R = R, abs(S), n.ext=n.ext, N.all=N.all,
                          var_mu_ext = var_mu_ext, g.var.est=g.var.est,meanR =meanR,
                          sumR=sumR, b=b, lower.tail = FALSE) # EmpSPA-G p value
    pval2 = GetProb_SPA_G(MAF, R = R, -abs(S), n.ext=n.ext, N.all=N.all,
                          var_mu_ext = var_mu_ext, g.var.est=g.var.est,meanR =meanR,
                          sumR=sumR, b=b, lower.tail = TRUE) # EmpSPA-G p value
    
    pval3 = pnorm(abs(z), lower.tail = FALSE) # Normal
    pval4 = pnorm(-abs(z), lower.tail = TRUE) # Normal
    
    pval.spa.G = pval1 + pval2
    pval.norm = pval3 + pval4
    
    # if(abs(z) < Cutoff){
    #   pval.spa.G = pval.norm
    # }
    
    return(c(pval.spa.G, pval.norm))
    
    
  }
  
}
################################################################################
prev= 0.05
ncluster=3
## read in pheno
pheno = fread("Benchmark_n42640.pheno")  %>%
  mutate( weight = ifelse(Event==1 , 1, (1-prev)/prev))

## fit null model
obj.wglm = glm(Event ~ SEX +AGE + PC1+PC2 +PC3+PC4, data = pheno,
               weight = weight, family = "binomial" )

pheno = pheno %>% mutate(R = obj.wglm$y - obj.wglm$fitted.values,
                         weight1  =  weight/(2*sum(pheno$weight)))

## Summix: cluster the who cohort into 3 subcohort
set.seed(1)
km <- kmeans(pheno%>%select(PC1, PC2, PC3, PC4  ), centers = ncluster, nstart = 25)
pheno$cluster = km$cluster

fwrite(pheno, file="pheno_cluster")

## read in geno
cat("read in geno \n")
G_all = fread("Benchmark_m127037xn42640.AltCounts.gz")
GenoInfo = G_all[,1:4]%>%mutate(newID = paste0(SNP,REF,ALT))
#Geno = t(G_all[,-(1:4)])
Geno = G_all%>%select(-CHR, -SNP, -REF,-ALT)%>%as.matrix()%>%t()
colnames(Geno) = GenoInfo$SNP
rm(G_all)

## read in external AFs -----------------------------------------------------------------------------

WB = fread("Benchmark_m127037_WB.afreq") %>%
  mutate(newID = paste0(ID,REF,ALT))%>% 
  merge(.,GenoInfo, by = c("newID"), all.y=T,sort=F)%>%
  mutate(ALT_FREQS = ifelse(REF.x==REF.y, ALT_FREQS, 1-ALT_FREQS ))%>%
  select(-REF.x,-ALT.x,-ID,-`#CHROM`) %>% rename(REF=REF.y, ALT = ALT.y )
Car = fread("Benchmark_m127037_Car.afreq")%>%
  mutate(newID = paste0(ID,REF,ALT))%>% 
  merge(.,GenoInfo, by = c("newID"), all.y=T,sort=F)%>%
  mutate(ALT_FREQS = ifelse(REF.x==REF.y, ALT_FREQS, 1-ALT_FREQS ))%>%
  select(-REF.x,-ALT.x,-ID,-`#CHROM`) %>% rename(REF=REF.y, ALT = ALT.y )
Ind = fread("Benchmark_m127037_Ind.afreq")%>%
  mutate(newID = paste0(ID,REF,ALT))%>% 
  merge(.,GenoInfo, by = c("newID"), all.y=T,sort=F)%>%
  mutate(ALT_FREQS = ifelse(REF.x==REF.y, ALT_FREQS, 1-ALT_FREQS ))%>%
  select(-REF.x,-ALT.x,-ID,-`#CHROM`) %>% rename(REF=REF.y, ALT = ALT.y )

## For each subcohort, perform WtCoxG

##read in sparseGRM
sparseGRM = fread("Benchmark_m127037xn42640.SparseGRM")

###For each subcohort, perform WtCoxG
gwas_list = lapply(1:ncluster, function(cl){
  
  cat("cluster", cl, "\n")
  
  # geno and pheno in each subcohort
  PhenoData = pheno%>%filter(cluster==cl)
  G = Geno[pheno$`#IID`[which(pheno$cluster==cl)],]
  
  ## allele frequency in cases,controls,external
  mergeGenoInfo = GenoInfo %>% mutate(mu0 = apply(G[PhenoData$Event==0,] , 2, mean)/2 ,
                                      mu1 = apply(G[PhenoData$Event==1,] , 2, mean)/2 ,
                                      mu.target = 0.5 * mu0 + 0.5 * mu1, ##Summix
                                      mu.int = ifelse(mu.target>0.5, 1-mu.target, mu.target), # for group AFs
                                      mu.w = prev * mu1 + (1-prev) * mu0 ,
                                      AF_WB = WB$ALT_FREQS,
                                      AN_WB = WB$OBS_CT,
                                      AF_Ind = Ind$ALT_FREQS,
                                      AN_Ind = Ind$OBS_CT,
                                      AF_Car = Car$ALT_FREQS,
                                      AN_Car = Car$OBS_CT
  )
  
  ##summix  --------------------------------------------------
  anc_prop= summix( data = na.omit(mergeGenoInfo),
                    reference=c("AF_WB","AF_Ind","AF_Car"),
                    observed="mu.target",
                    pi.start = c(0.3,0.3,0.4) )
  anc_prop = unlist(anc_prop)[c("AF_WB","AF_Ind","AF_Car")]
  anc_prop=ifelse(is.na(anc_prop),0, anc_prop)%>%setNames(c("AF_WB","AF_Ind","AF_Car"))
  
  mergeGenoInfo =mergeGenoInfo %>%
    mutate(AF_ref = anc_prop["AF_WB"] *AF_WB + anc_prop["AF_Ind"] * AF_Ind + anc_prop["AF_Car"] * AF_Car,
           AN_ref = 1/(sum(anc_prop["AF_WB"]^2 /AN_WB[1]  + anc_prop["AF_Ind"]^2 /AN_Ind[1] + anc_prop["AF_Car"]^2 /AN_Car[1])))
  
  ## batch effect testing  -----------------------------------------------------------
  ## variance ratio
  sparseGRMsub = sparseGRM %>% filter(ID1 %in% PhenoData$`#IID`,
                                      ID2 %in% PhenoData$`#IID` )
  
  cat("variance ratio ---------\n")
  if(!is.null(sparseGRMsub)){
    
    w1 = PhenoData$weight1
    names(w1) = PhenoData$`#IID`
    R_tilde = PhenoData$R - mean( PhenoData$R)
    names(R_tilde) = PhenoData$`#IID`
    meanR = mean(PhenoData$R)
    
    sparseGRMsub = sparseGRMsub %>% mutate(cov = Value * w1[as.character(ID1)] * w1[as.character(ID2)],
                                           cov_R = Value * R_tilde[as.character(ID1)] * R_tilde[as.character(ID2)]
    )%>%mutate(cov = ifelse(ID1==ID2, cov, 2*cov),
               cov_R = ifelse(ID1==ID2, cov_R, 2*cov_R))
    
    
    var.ratio.w0 = (sum(sparseGRMsub$cov) + 1/(2*mergeGenoInfo$AN_ref))/(sum(w1^2) + 1/(2*mergeGenoInfo$AN_ref))
    var.ratio.int = sum(sparseGRMsub$cov_R)/sum(R_tilde^2)
    
    
  }else{var.ratio.w0 = var.ratio.int= 1}
  event_sub = mean(PhenoData$Event)
  
  mergeGenoInfo = mergeGenoInfo %>%
    mutate(var.ratio.w0 = var.ratio.w0,
           var.ratio.int = var.ratio.int
    )
  
  mergeGenoInfo$pvalue_bat = lapply(1:nrow(mergeGenoInfo), function(i){
    if(i%%1000==0)cat(i,"\n")
    n1 = sum(PhenoData$Event)
    n0 = sum(1 - PhenoData$Event)
    p_bat = Batcheffect.Test(n0 = n0,
                             n1 = n1,
                             n.ext = mergeGenoInfo$AN_ref[i]/2, 
                             maf0 = mergeGenoInfo$mu0[i], 
                             maf1 = mergeGenoInfo$mu1[i], 
                             maf.ext = mergeGenoInfo$AF_ref[i] , 
                             pop.prev=prev,
                             var.ratio = mergeGenoInfo$var.ratio.w0[i]
    )
    
    return(p_bat)
  })%>%unlist()
  mergeGenoInfo = mergeGenoInfo %>% mutate(index=1:n())
  
  ### estimate batcheffect parameters  ----------------------------------------------
  cat("Estimate TPR and sigma2--------------\n")
  maf.group = c(seq(-0.00001, 0.4, 0.05),max(mergeGenoInfo$mu.int))
  
  mergeGenoInfo=lapply(1:(length(maf.group)-1), function(i){
    cat(i,"\n")
    
    ##assume that genotypes with MAF in [ maf.group[i] , maf.group[i+1]] have the same mixture distribution
    data = mergeGenoInfo %>% filter(mu.int > maf.group[i] & mu.int <= maf.group[i+1] )
    
    ##using batcheffect p-values with MAF in [maf.group[i]-0.1 , maf.group[i+1]+0.1] to estimate parameters
    data.ref = mergeGenoInfo %>%
      filter(mu.int >= max(maf.group[i]-0.1,0) & mu.int < min(1,maf.group[i+1]+0.1) )
    
    mu = (maf.group[i]+maf.group[i+1])/2
    
    n.ext = na.omit(data$AN_ref)[1]/2
    var_mu_ext = mu*(1-mu)/(2*n.ext)
    
    var_Sbat = ifelse(is.null(sparseGRM), sum(w1^2)*2*mu*(1-mu) + var_mu_ext,
                      na.omit(mergeGenoInfo$var.ratio.w0)[1] * (sum(w1^2)*2*mu*(1-mu) + var_mu_ext) )
    
    obj = fun.est.param(vec_p_bat=data.ref$pvalue_bat ,
                        vec_var_Sbat=var_Sbat)
    TPR = obj[1]
    sigma2 = obj[2]
    
    w.ext = optim(par = 0.5, method = "L-BFGS-B", lower = 0, upper = 1
                  , fn = fun.optimalWeight
                  , pop.prev = prev
                  , y = PhenoData$Event
                  , R = PhenoData$R
                  , w = PhenoData$weight
                  , mu = mu
                  , N = nrow(PhenoData)
                  , n.ext = n.ext
                  , sigma2 = obj$sigma2
                  , TPR = obj$TPR
    )$par[1]
    
    if(is.null(sparseGRMsub)){
      var.ratio.ext = 1
    } else{
      R_tilde_w = PhenoData$R - mean(PhenoData$R) * w.ext
      names(R_tilde_w) = PhenoData$`#IID`
      sparseGRMsub = sparseGRMsub %>% mutate(cov_Rext = Value * R_tilde_w[as.character(ID1)] * R_tilde_w[as.character(ID2)])%>%
        mutate(cov_Rext = ifelse(ID1==ID2, cov_Rext, 2*cov_Rext))
      var.ratio.ext = (sum(sparseGRMsub$cov_Rext)  + w.ext^2 * sum(PhenoData$R)^2/n.ext)/(sum(R_tilde_w^2) + w.ext^2 * sum(PhenoData$R)^2/n.ext)
      
    }
    
    data=data%>%cbind(.,TPR, sigma2, w.ext, var.ratio.ext)
    
  })%>%
    do.call("rbind",.) %>%
    as_tibble() %>%
    arrange(index) %>%
    select(-index)
  
  
  ### GWAS -----------------------------------
  GWAS = lapply(1:ncol(G),function(i){
    if(i%%1000==0)cat("Complete ",i,"/",ncol(G),"\n")
    
    g = G[,i]
    R = PhenoData$R
    w = PhenoData$weight
    
    mu.ext = mergeGenoInfo$AF_ref[i]
    n.ext = mergeGenoInfo$AN_ref[i]/2
    TPR = mergeGenoInfo$TPR[i]
    sigma2 = mergeGenoInfo$sigma2[i]
    p_bat = mergeGenoInfo$pvalue_bat[i]
    w.ext = mergeGenoInfo$w.ext[i]
    var.ratio.w0 = mergeGenoInfo$var.ratio.w0[i]
    var.ratio.int = mergeGenoInfo$var.ratio.int[i]
    var.ratio0 = mergeGenoInfo$var.ratio.ext[i]
    
    WtCoxG.ext = WtCoxG.test(g = g,
                             R = R,
                             w = w,
                             TPR=TPR,
                             sigma2 = sigma2,
                             b = w.ext,
                             var.ratio.w0 =var.ratio.w0,
                             var.ratio.w1 = var.ratio.w0,
                             var.ratio0 = var.ratio0,
                             var.ratio1 = var.ratio0,
                             mu.ext = mu.ext,
                             n.ext = n.ext,
                             p_bat = p_bat) # pvalue and Score statistics
    
    return( WtCoxG.ext )
  }) %>%
    do.call("rbind",.) %>%
    cbind(.,mergeGenoInfo)  
  
  fwrite(GWAS, file = paste0("gwas_cluster", cl, ".csv"))
  return(GWAS)
  
})



#### meta-analysis --------------------------------------------------------------
all_S = do.call(cbind, lapply(gwas_list, function(d) d$S))

# Extract all P-values ('WtCoxG.ext') into a single matrix
all_P = do.call(cbind, lapply(gwas_list, function(d) d$p.con))

# ---------------------------------------------------------
# 3. Variance Recovery & Calculation
# ---------------------------------------------------------
# Logic: To meta-analyze, we need the Variance of S.
# Since the input only provides S and P, we back-calculate Variance.
# Formula: ChiSq = S^2 / Var  =>  Var = S^2 / ChiSq

# Calculate Chi-Square statistics from P-values (1 d.f.)
all_ChiSq = qchisq(all_P, df = 1, lower.tail = FALSE)

# Calculate Variance
all_Var = (all_S^2) / all_ChiSq

# SAFETY: If S=0 and P=1, the result might be NaN (0/0). Set Variance to 0 or NA.
all_Var[is.na(all_Var)] = 0 

# ---------------------------------------------------------
# 4. Meta-Analysis (Score Test)
# ---------------------------------------------------------
message("Running Meta-analysis (Vectorized)...")

# Sum scores and variances across all studies (row-wise sum)
# rowSums is extremely fast compared to loops.
sum_S = rowSums(all_S, na.rm = TRUE)
sum_Var = rowSums(all_Var, na.rm = TRUE)

# Calculate Meta Z-score
# Formula: Z_meta = Sum(S) / Sqrt(Sum(Var))
Z_meta = sum_S / sqrt(sum_Var)

# Calculate two-tailed P-value from Z-score
P_meta = 2 * pnorm(-abs(Z_meta))

Meta_output = data.table(
  SNP = gwas_list[[1]]$SNP,
  Meta_S = sum_S,
  Meta_Var = sum_Var,
  Meta_Z = Z_meta,
  Meta_P = P_meta
)

fwrite(Meta_output, file = "metaGWAS.csv" )
