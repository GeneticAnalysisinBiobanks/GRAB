##### SPAmix functions with input of local ancestry information

# Gets NULL model residuals for SPA-G
SPA_G_Get_Resid = function(traits="survival/binary",
                           formula=NULL,
                           data=NULL,
                           pIDs=NULL,
                           gIDs=NULL,
                           range=c(-100,100),
                           length.out = 10000,
                           ...)
  
{
  if(traits=="survival"){
    Call = match.call()
    ### Fit a Cox model
    obj.coxph = coxph(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.coxph, range)
    
    ### Get the covariate matrix to adjust for genotype
    resid = obj.coxph$residuals
    
    re = resid
  }
  else if(traits=="binary"){
    Call = match.call()
    ### Fit a logistic model
    obj.logistic = glm(formula, data=data, x=T, ...)
    ### Check input arguments
    p2g = check_input(pIDs, gIDs, obj.logistic, range)
    
    ### Get the covariate matrix to adjust for genotype
    mu = obj.logistic$fitted.values
    resid = obj.logistic$y - mu
    re = resid
  }
  return(re)
}

# MAF is the minor allele frequency of the tested marker in the index ancestry
# MAFVec is the vector (length = sample size) of minor allele frequency of the tested marker in the index ancestry
# haplo_numVec is the vector (length = sample size) of the number haplotype of the tested marker in the index ancestry

# The MGF of G (genotype from the index ancestry)
M_G0 = function(t, MAF, haplo_num){
  re = (1 - MAF + MAF * exp(t))^haplo_num
  return(re)
}

# The first derivative of the MGF of G (genotype from the index ancestry)
M_G1 = function(t, MAF, haplo_num){
  re = haplo_num*(MAF * exp(t))*(1 - MAF + MAF * exp(t))^(haplo_num - 1)
  return(re)                           
}

# The second derivative of the MGF of G (genotype from the index ancestry)
M_G2 = function(t, MAF, haplo_num){
  re = haplo_num * (haplo_num-1) * (MAF * exp(t))^2 * (1 - MAF + MAF * exp(t))^(haplo_num - 2) +
    haplo_num * (MAF * exp(t)) * (1 - MAF + MAF * exp(t)) ^ (haplo_num - 1)
  return(re)
}


# The CGF of G (genotype from the index ancestry)
K_G0 = function(t, MAF, haplo_num){
  re = log(M_G0(t, MAF, haplo_num))
  return(re)
}

# The first derivative of the CGF of G (genotype from the index ancestry)
K_G1 = function(t, MAF, haplo_num){
  re = M_G1(t, MAF, haplo_num)/M_G0(t, MAF, haplo_num)
  return(re)
}

# The second derivative of the CGF of G (genotype from the index ancestry)
K_G2 = function(t, MAF, haplo_num){
  re = (M_G0(t, MAF, haplo_num)*M_G2(t, MAF, haplo_num)-M_G1(t, MAF, haplo_num)^2)/M_G0(t, MAF, haplo_num)^2
  return(re)
}

# The CGF of score test statistic for the index ancestry
H_org = function(t, R, MAFVec, haplo_numVec){
  
  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(K_G0(t1*R, MAFVec, haplo_numVec))
  }
  return(out)
}

# The first derivative of the CGF of score test statistic for the index ancestry
H1_adj = function(t, R, s, MAFVec, haplo_numVec)
{
  n.t = length(t)
  out = rep(0,n.t)
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(R*K_G1(t1*R, MAFVec, haplo_numVec)) - s
  }
  return(out)
}


# The second derivative of the CGF of score test statistic for the index ancestry
H2 = function(t, R, MAFVec, haplo_numVec)
{
  n.t = length(t)
  out = rep(0,n.t)
  
  for(i in 1:n.t){
    t1 = t[i]
    out[i] = sum(R^2*K_G2(t1*R, MAFVec, haplo_numVec))
  }
  return(out)
}


GetProb_SPA_G = function(MAFVec, R, haplo_numVec, s, lower.tail){
  
  out = uniroot(H1_adj, c(-20,20), extendInt = "upX",
                R=R, s=s, MAFVec=MAFVec, haplo_numVec=haplo_numVec)
  zeta = out$root
  
  k1 = H_org(t=zeta, R=R, MAFVec=MAFVec, haplo_numVec=haplo_numVec)
  k2 = H2(t=zeta, R=R, MAFVec=MAFVec, haplo_numVec=haplo_numVec)
  
  temp1 = zeta * s - k1
  
  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  
  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  re = pval
  return(re)
}




#### SPAmix #########################################################################

SPAmix_localance_one_SNP = function(g,
                                    R,
                                    haplo_numVec,
                                    # obj.null,
                                    Cutoff = 2,
                                    impute.method = "fixed",
                                    missing.cutoff = 0.15,
                                    min.maf = 0.00001,          # update on 2022-08-16 : replace 0.0001 by 0.000001
                                    G.model = "Add")
{
  ## calculate MAF and update genotype vector
  # MAF = mean(g, na.rm=T)/2     # MAF in the index ancestry
  MAF = sum(g)/sum(haplo_numVec)   # MAF in the index ancestry
  
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N
  
  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }
  ####################################################### AF, not MAF  
  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }
  
  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)
  
  if(MAF < min.maf)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA, NA, NA))
  
  ## Score statistic
  # R = obj.null$resid
  S = sum(g * R)
  
  ## estimated variance without adjusting for covariates
  # N1set = 1:N
  # N0 = 0
  
  N1set = which(g!=0)  # position of non-zero genotypes
  N0 = N-length(N1set)
  
  MAC = sum(g)    
  MAF.est.Vec = c(rep(MAF, N))
   
  S.mean = sum(MAF.est.Vec * haplo_numVec * R)
  
  g.var.est.Vec = haplo_numVec * MAF.est.Vec * (1 - MAF.est.Vec)
  S.var = sum(R^2 * g.var.est.Vec)
  
  z = (S - S.mean)/sqrt(S.var)
  
  if(abs(z) < Cutoff){
    pval.norm = pnorm(abs(z), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.mean, S.var, z, MAC))  
  }
  
  pval1 = GetProb_SPA_G(MAFVec=MAF.est.Vec, R=R, haplo_numVec=haplo_numVec, s=max(S, (2*S.mean-S)), lower.tail = FALSE) # SPA-G p value 
  pval2 = GetProb_SPA_G(MAFVec=MAF.est.Vec, R=R, haplo_numVec=haplo_numVec, s=min(S, (2*S.mean-S)), lower.tail = TRUE) # SPA-G p value 
  
  pval3 = pnorm(abs(z), lower.tail = FALSE) # Normal
  pval4 = pnorm(-abs(z), lower.tail = TRUE) # Normal
  
  pval.spa.G = pval1 + pval2
  pval.norm = pval3 + pval4
  
  pval = c(pval.spa.G, pval.norm) # 2 elements: element 1 is from SPA, element 2 is from Normal
  
  return(c(MAF, missing.rate, pval, S, S.mean, S.var, z, MAC))  
}


SPAmix_localance = function(Geno.mtx,
                            R,
                            haplo.mtx,
                            # obj.null,
                            Cutoff = 2,
                            impute.method = "fixed",
                            missing.cutoff = 0.15,
                            min.maf = 0.00001,          
                            G.model = "Add")
  
  
{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  G.model=G.model)
  
  # check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))
  
  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 9) 
  colnames(output) = c("MAF","missing.rate","Pvalue.SPAmix.index.ance","Pvalue.norm.index.ance",
                       "Stat","Mean","Var","z", "MAC") # update on 2022-10-05 : MAF.est.negative.num 
  rownames(output) = colnames(Geno.mtx)
  
  ## Start analysis
  print("Start Analyzing...")
  print(Sys.time())
  
  # Cycle for genotype matrix
  for(i in 1:n.Geno){
    # print(i)
    g = Geno.mtx[,i]
    haplo_numVec = haplo.mtx[,i]
    output.one.SNP = SPAmix_localance_one_SNP(g,
                                              R,
                                              haplo_numVec,
                                              Cutoff,
                                              impute.method,
                                              missing.cutoff,
                                              min.maf,         
                                              G.model)

    output[i,] = output.one.SNP
  }
  
  output = as.data.frame(cbind(rsID = colnames(Geno.mtx), output))
  
  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}
















check_input = function(pIDs, gIDs, obj, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = obj$residuals
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")
  
  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)
  
  return(p2g)
}

check_input_Resid = function(pIDs, gIDs, R, range)
{
  if(is.null(pIDs) & is.null(gIDs)) stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")
  if(any(sort(unique(pIDs))!=sort(unique(gIDs)))) stop("unique(pIDs) should be the same as unique(gIDs).")
  if(anyDuplicated(gIDs)!=0) stop("Argument 'gIDs' should not have a duplicated element.")
  if(range[2]!=-1*range[1]) stop("range[2] should be -1*range[1]")
  mresid = R
  if(length(mresid)!=length(pIDs)) stop("Argument 'pIDs' should be of the same length as input data.")
  
  if(all(pIDs == gIDs)) p2g = NULL
  else p2g = match(pIDs, gIDs)
  
  return(p2g)
}










