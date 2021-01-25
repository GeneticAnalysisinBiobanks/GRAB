
setMarker.SPACox = function(objNull, control)
{
  
}


mainMarker.SPACox = function(objNull, control)
{
  SPACoxinCPP();
  
  ## Score statistic
  S = sum(g * obj.null$resid)
  
  ## estimated variance without adjusting for covariates
  G1 = g - 2*MAF   # centered genotype (such that mean=0)
  S.var1 = obj.null$var.resid * sum(G1^2)
  z1 = S/sqrt(S.var1)
  
  if(abs(z1) < Cutoff){
    pval.norm = pnorm(abs(z1), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.var1, z1))
  }
  
  N1set = which(g!=0)  # position of non-zero genotypes
  N0 = N-length(N1set)
  
  G1norm = G1/sqrt(S.var1)  # normalized genotype (such that sd=1)
  
  G1N1 = G1norm[N1set]
  G1N0 = -2*MAF/sqrt(S.var1)   # all subjects with g=0 share the same normlized genotype, this is to reduce computation time
  
  pval1 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, abs(z1), lower.tail = FALSE)
  pval2 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, -abs(z1), lower.tail = TRUE)
  pval = pval1 + pval2
  
  if(pval[1] > CovAdj.cutoff)
    return(c(MAF, missing.rate, pval, S, S.var1, z1))
  
  ## estimated variance after adjusting for covariates
  
  G2 = g - obj.null$X.invXX %*% (obj.null$tX[,N1set,drop=F] %*% g[N1set])
  S.var2 = obj.null$var.resid * sum(G2^2)
  z2 = S/sqrt(S.var2)
  
  G2norm = G2/sqrt(S.var2)
  
  N1set = 1:N
  N0 = 0
  G2N1 = G2norm
  G2N0 = 0   # since N0=0, this value actually does not matter
  
  pval1 = GetProb_SPA(obj.null, G2N1, G2N0, N1set, N0, abs(z2), lower.tail = FALSE)
  pval2 = GetProb_SPA(obj.null, G2N1, G2N0, N1set, N0, -abs(z2), lower.tail = TRUE)
  pval = pval1 + pval2
  
}


############## The below are subfunctions in SPACox.R #################

GetProb_SPA = function(obj.null, G2NB, G2NA, NBset, N0, q2, lower.tail){
  
  out = uniroot(K1_adj, c(-20,20), extendInt = "upX",
                G2NB=G2NB, G2NA=G2NA, NBset=NBset,
                N0=N0, q2=q2, obj.null=obj.null)
  zeta = out$root
  
  k1 = K_org(zeta,  G2NB=G2NB, G2NA=G2NA, NBset=NBset, N0=N0, obj.null=obj.null)
  k2 = K2(zeta,  G2NB=G2NB, G2NA=G2NA, NBset=NBset, N0=N0, obj.null=obj.null)
  
  temp1 = zeta * q2 - k1
  
  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}
  
  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  pval.norm = pnorm(q2, lower.tail = lower.tail)
  
  re = c(pval, pval.norm)
  return(re)
}


K_org = function(t, G2NB, G2NA, NBset, N0, obj.null){
  
  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*obj.null$K_org_emp(t2NA) + sum(obj.null$K_org_emp(t2NB))
  }
  return(out)
}

K1_adj = function(t, G2NB, G2NA, NBset, N0, q2, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)
  
  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*G2NA*obj.null$K_1_emp(t2NA) + sum(G2NB*obj.null$K_1_emp(t2NB)) - q2
  }
  return(out)
}

K2 = function(t, G2NB, G2NA, NBset, N0, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)
  
  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*G2NA^2*obj.null$K_2_emp(t2NA) + sum(G2NB^2*obj.null$K_2_emp(t2NB))
  }
  return(out)
}



#' SaddlePoint Approximation implementation of Cox regression surival analysis (One-SNP-version)
#'
#' One-SNP-version SPACox function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param others the same as function SPACox. NOTE that we do not check subject order in this one-snp-version !!!
#' @return the same as function SPACox.
#' @export
SPACox.one.SNP = function(g,
                          obj.null,
                          Cutoff = 2,
                          impute.method = "fixed",
                          missing.cutoff = 0.15,
                          min.maf = 0.0001,
                          CovAdj.cutoff = 5e-5,
                          G.model = "Add")
{
  g[g==-9]=NA  # since we add plink input
  ## calculate MAF and update genotype vector
  MAF = mean(g, na.rm=T)/2
  N = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/N

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf || missing.rate > missing.cutoff)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA))

  if(!is.null(obj.null$p2g))
    g = g[obj.null$p2g]

  

  return(c(MAF, missing.rate, pval, S, S.var2, z2))
}





check_input = function(pIDs, gIDs, obj.coxph, range)
{
  if(is.null(pIDs) & is.null(gIDs))
    stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")

  pIDs = as.character(pIDs)
  gIDs = as.character(gIDs)
  if(!is.null(obj.coxph$na.action)){
    posNA = c(obj.coxph$na.action)
    if(any(posNA > length(pIDs)))
      stop("Number of input data is larger than length(pIDs).")
    pIDsNA = pIDs[posNA]

    print(paste0("Due to missing data in response/indicators, ",length(posNA)," entries are removed from analysis."))
    print("If concerned about the power loss, we suggest users impute data first and then use SPACox package.")
    print(head(cbind(posNA=posNA, pIDsNA=pIDsNA)))

    pIDs = pIDs[-1*posNA]  # remove IDs with missing data
  }

  if(any(!is.element(pIDs, gIDs)))
    stop("All elements in pIDs should be also in gIDs.")

  if(anyDuplicated(gIDs)!=0)
    stop("Argument 'gIDs' should not have a duplicated element.")

  if(range[2]!=-1*range[1])
    stop("range[2] should be -1*range[1]")

  mresid = obj.coxph$residuals

  if(length(mresid)!=length(pIDs))
    stop("length(mresid)!=length(pIDs) where mresid is the martingale residuals from coxph() in survival package.")

  p2g = NULL
  if(length(pIDs)!=length(gIDs)){
    p2g = match(pIDs, gIDs)
  }else{
    if(any(pIDs != gIDs))
      p2g = match(pIDs, gIDs)
  }

  return(list(p2g=p2g,pIDs=pIDs))
}

check_input1 = function(obj.null, Geno.mtx, par.list)
{
  if(class(obj.null)!="SPACox_NULL_Model")
    stop("obj.null should be a returned outcome from SPACox_Null_Model()")

  if(any(obj.null$gIDs != rownames(Geno.mtx))) stop("gIDs should be the same as rownames(Geno.mtx).")
  if(is.null(rownames(Geno.mtx))) stop("Row names of 'Geno.mtx' should be given.")
  if(is.null(colnames(Geno.mtx))) stop("Column names of 'Geno.mtx' should be given.")
  if(!is.numeric(Geno.mtx)|!is.matrix(Geno.mtx)) stop("Input 'Geno.mtx' should be a numeric matrix.")

  if(!is.numeric(par.list$min.maf)|par.list$min.maf<0|par.list$min.maf>0.5) stop("Argument 'min.maf' should be a numeric value >= 0 and <= 0.5.")
  if(!is.numeric(par.list$Cutoff)|par.list$Cutoff<0) stop("Argument 'Cutoff' should be a numeric value >= 0.")
  # if(!is.element(par.list$impute.method,c("none","bestguess","random","fixed"))) stop("Argument 'impute.method' should be 'none', 'bestguess', 'random' or 'fixed'.")
  if(!is.element(par.list$impute.method,c("fixed"))) stop("Argument 'impute.method' should be 'fixed'.")
  if(!is.numeric(par.list$missing.cutoff)|par.list$missing.cutoff<0|par.list$missing.cutoff>1) stop("Argument 'missing.cutoff' should be a numeric value between 0 and 1.")
  if(!is.element(par.list$G.model,c("Add","Dom","Rec"))) stop("Argument 'G.model' should be 'Add', 'Dom' or 'Rec'.")
}

