
#' GRAB: simulate random effect (i.e. bVec) based on family structure
#' 
#' Simulate random effect (i.e. bVec) based on family structure
#' 
#' @param n.fam number of families in simulation
#' @param fam.kin a matrix of GRM (kinship matrix)
#' @param tau variance component
#' @return a random effect following a multivariate normal distribution
#' @examples 
#' fam.kin = read.table(system.file("extdata", "example_10members.kin.txt", package = "GRAB"))
#' fam.kin = as.matrix(fam.kin)
#' n.fam = 100
#' tau = 2
#' bVec = GRAB.SimubVec(n.fam, fam.kin, tau)
#'      
#' @export
GRAB.SimubVec = function(n.fam,
                         fam.kin,
                         tau)
{
  fam.kin = as.matrix(fam.kin)
  n = n.fam * nrow(fam.kin)
  out.eigen = eigen(fam.kin)
  factor = t(out.eigen$vectors) * sqrt(out.eigen$values)
  kin.chol = diag(n.fam) %x% factor
  b.true = t(kin.chol) %*% rnorm(n) * sqrt(tau) 
  return(b.true)
}


#' GRAB: simulate genotype matrix based on family structure
#'
#' Simulate genotype matrix based on family structure using Haplotype stored in SKAT package
#'
#' @param n.fam number of families in simulation
#' @param n.subj number of unrelated subjects in simulation
#' @param fam.kin.mode "4-members" or "10-members" to simulate two family structures. "10-members" corresponds to system.file("extdata", "example_10members.kin.txt", package = "GRAB")
#' @param min.MAF minimal MAF in simulation
#' @param max.MAF maximal MAF in simulation
#' @param seedNum seed number for random simulation
#' @param region.length length of region to simulate. Default value is 2000, i.e., 2KB. 
#' @param n.region number of regions to simulate
#' @return a matrix of genotype data including multiple Genes
#' @examples
#' fam.kin.10.members = read.table(system.file("extdata", "example_10members.kin.txt", package = "GRAB"))
#' fam.kin.10.members = as.matrix(fam.kin.10.members)
#' n.fam = 500
#' n.subj = 5000
#' fam.kin.mode = "10-members"  # "10-members" corresponds to fam.kin.10.members
#' min.MAF = 0
#' max.MAF = 0.3
#' seedNum = 1234
#' region.length = 2000
#' n.region = 2
#' GMat = GRAB.SimuGMatSKAT(n.fam, n.subj, fam.kin.mode, min.MAF, max.MAF, seedNum, region.length, n.region)
#'
#' @export
#' @import SKAT
GRAB.SimuGMatSKAT = function(n.fam,
                             n.subj, 
                             fam.kin.mode,   # "4-members" or "10-members" 
                             min.MAF, 
                             max.MAF,
                             seedNum, 
                             region.length = 2000, 
                             n.region,
                             prob = 1)
{
  data(SKAT.haplotypes)
  attach(SKAT.haplotypes)
  SNPInfo = SKAT.haplotypes$SNPInfo
  Haplotype = SKAT.haplotypes$Haplotype
  
  set.seed(seedNum)
  min.pos = 79
  max.pos = 199956
  start.posVec = floor(runif(n.region, min.pos, max.pos-region.length))
  end.posVec = start.posVec + region.length
  
  G_FAM_Set = c()
  for(i in 1:n.region){
    G_FAM = Get_One_Set(start.posVec[i], end.posVec[i], Haplotype, SNPInfo, n.fam, n.subj, prob, fam.kin.mode)
    
    MAF = colMeans(G_FAM) / 2
    MAF = ifelse(MAF > 0.5, 1 - MAF, MAF)
    idx.MAF = which(MAF > min.MAF & MAF < max.MAF)
    G_FAM = G_FAM[,idx.MAF]
    colnames(G_FAM) = paste0("GENE-",i,"-",colnames(G_FAM))
    G_FAM_Set = cbind(G_FAM_Set, G_FAM)
  }
  
  return(G_FAM_Set)
}

Get_One_Set = function(start.pos, 
                       end.pos, 
                       Haplotype, 
                       SNPInfo, 
                       n.fam,
                       n.subj,
                       prob,
                       fam.kin.mode)
{
  if(fam.kin.mode == "10-members")
    G_FAM = Get_One_Set_10_members(start.pos, 
                                   end.pos, 
                                   Haplotype, 
                                   SNPInfo, 
                                   n.fam,
                                   n.subj,
                                   prob)
  
  return(G_FAM)
}

## 1+2->5+6; 3+5->7+8; 4+6->9+10 
Get_One_Set_10_members = function(start.pos, 
                                  end.pos, 
                                  Haplotype, 
                                  SNPInfo, 
                                  n.fam,
                                  n.subj,
                                  prob)
{
  idx = which(SNPInfo$CHROM_POS >= start.pos & SNPInfo$CHROM_POS <= end.pos);
  idx = idx[which(rbinom(length(idx), 1, prob)==1)]
  n = n.fam * 10 + n.subj
  m = length(idx)
  
  if(m <= 1) 
    stop("double check Get_One_Set_10_members()")
  
  ID = list()
  
  G_FAM = matrix(0, n, m)
  idxSubj = 1
  for(j in 1:n.fam){
    if(j %% 100 == 0) 
      sprintf("Finished simulation for %d/%d families.", j, n.fam)
    
    # random select 2 haplotypes from 10,000 haplotypes
    ID[[1]] = sample.int(10000, 2)
    ID[[2]] = sample.int(10000, 2)
    ID[[3]] = sample.int(10000, 2)
    ID[[4]] = sample.int(10000, 2)
    ID[[5]] = c(sample(ID[[1]],1), sample(ID[[2]], 1))
    ID[[6]] = c(sample(ID[[1]],1), sample(ID[[2]], 1))
    ID[[7]] = c(sample(ID[[3]],1), sample(ID[[5]], 1))
    ID[[8]] = c(sample(ID[[3]],1), sample(ID[[5]], 1))
    ID[[9]] = c(sample(ID[[4]],1), sample(ID[[6]], 1))
    ID[[10]] = c(sample(ID[[4]],1), sample(ID[[6]], 1))	
    
    for(i in 1:10){
      G = colSums(Haplotype[ID[[i]], idx])
      G_FAM[idxSubj,] = G
      idxSubj = idxSubj + 1
    }
  }
  
  if(n.subj > 0){
    for(j in 1:n.subj){
      if(j %% 1000 == 0) sprintf("Finished simulation for %d/%d unrelated subjects.", j, n.subj)
      
      ID = sample.int(10000, 2)
      G = colSums(Haplotype[ID, idx])
      G_FAM[idxSubj,] = G
      idxSubj = idxSubj + 1
    }
  }
  
  if(n.subj > 0){
    rownames(G_FAM) = c(paste0("f",rep(1:n.fam, each=10), "_", rep(1:10, n.fam)),
                        paste0("fid",1:n.subj,"_","iid",1:n.subj))
  }else{
    rownames(G_FAM) = paste0("f",rep(1:n.fam, each=10), "_", rep(1:10, n.fam))
  }
  
  colnames(G_FAM) = paste0("SNP_",SNPInfo$CHROM_POS[idx])
  
  return(G_FAM)
}


#' GRAB: make plink files using Geno.mat
#' 
#' Make plink files using Geno.mat, rownames(Geno.mat) should be subject IDs and colnames(Geno.mat) should be SNP IDs
#' 
#' @param Geno.mat Genotype matrix with row names of subject IDs and column names of SNP IDs
#' @param work.dir working directory to store the plink files
#' @param out.prefix prefix of the output Plink files
#' @param mCHRs if TRUE, then SNPs were randomly assigned to chromosomes 1-22; if FALSE, then all SNPs are in chromosome 1
#' @param BP if NULL, then Base Position (BP) is 1:length(SNPs), otherwise, the length(BP) should be the same as ncol(Geno.mat)
#' @return Plink files are stored in 'work.dir' with prefix of 'out.prefix'.
#' @examples 
#' n = 1000
#' m = 20
#' MAF = 0.3
#' Geno.mat = matrix(rbinom(n*m, 2, MAF), n, m)
#' rownames(Geno.mat) = paste0("subj-",1:n)
#' colnames(Geno.mat) = paste0("SNP-",1:m)
#' work.dir = system.file("results", package = "GRAB")
#' out.prefix = "test"
#' mCHRs = FALSE
#' BP = NULL
#' plink.make(Geno.mat, work.dir, out.prefix, mCHRs, BP)
#'      
#' # cmd = paste0("plink --file ", OUT.file,  " --make-bed --out ", OUT.file)
#' # system(cmd)
#' # system(paste("rm", MAP.file, PED.file))
#'      
#' @export
plink.make = function(Geno.mat,
                      work.dir,
                      out.prefix,
                      mCHRs = T,         # multiple CHRs
                      BP = NULL)
{
  SNPs = colnames(Geno.mat)
  IID = rownames(Geno.mat)
  if(is.null(SNPs) | is.null(IID))
    stop("rownames and colnames of Geno.mat should be given.")
  
  # FID = sapply(IID, FUN=function(x){unlist(strsplit(x,split = "_"))[1]})
  FID = IID
  
  nSNPs = length(SNPs)
  if(mCHRs){
    CHRs = sample(1:22, nSNPs, replace=T)
  }else{
    CHRs = 1
  }
  
  if(is.null(BP)){
    MAP = cbind(CHR=CHRs,SNP=SNPs,GeneDist=0,BP=1:nSNPs);
  }else{
    MAP = cbind(CHR=CHRs,SNP=SNPs,GeneDist=0,BP=BP);
  }
  
  PED = cbind(FID=FID, IID=IID, PID=0, MID=0, Sex=1, Phen=-9);
  m = length(SNPs)
  n = length(IID)
  Geno.ped = matrix(nrow = n, ncol = 2*m)
  Geno.ped[,seq(1,2*m,2)] = ifelse(Geno.mat>=1,"A","G");
  Geno.ped[,seq(2,2*m,2)] = ifelse(Geno.mat>=2,"A","G");
  
  PED = cbind(PED,Geno.ped)
  
  MAP.file = paste0(work.dir, "/", out.prefix,".map")
  PED.file = paste0(work.dir, "/", out.prefix,".ped")
  OUT.file = paste0(work.dir, "/", out.prefix)
  
  write.table(MAP, MAP.file, quote = F, col.names = F, row.names = F)
  write.table(PED, PED.file, quote = F, col.names = F, row.names = F)
  
  message = paste0("Check directory of ", work.dir, " for the PLINK files.")
  return(message)
  
  # cmd = paste0("plink --file ", OUT.file,  " --make-bed --out ", OUT.file)
  # system(cmd)
  # system(paste("rm", MAP.file, PED.file))
}

