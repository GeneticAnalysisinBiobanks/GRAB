
#' Simulate an R matrix of genotype data
#' 
#' \code{GRAB} package provides functions to simulate genotype data. We support simulations based on unrelated subjects and related subjects.
#' 
#' @param nSub the number of unrelated subjects in simulations, if \code{nSub = 0}, then all subjects are related to at least one of the others.
#' @param nFam the number of families in simulation, if \code{nFam = 0}, then all subjects are unrelated to each other.
#' @param FamMode \code{"4-members"}, \code{"10-members"}, or \code{"20-members"}. Check \code{Details} section for more details.
#' @param nSNP number of markers to simulate
#' @param MaxMAF a numeric value *(default=0.5)*, haplotype is simulated with allele frequency <= this value.
#' @param MinMAF a numeric value *(default=0.05)*, haplotype is simulated with allele frequency >= this value.
#' @param MAF a numeric vector with a length of *nSNP*. If this argument is given, then arguments of *MaxMAF* and *MinMAF* would be ignored.
#' @return an R list including genotype matrix and marker information
#' \itemize{
#'   \item \code{GenoMat} a numeric matrix of genotype: each row is for one subject and each column is for one SNP
#'   \item \code{markerInfo} a data frame with the following 2 columns: SNP ID and minor allele frequency
#' } 
#' @seealso \code{\link{GRAB.makePlink}} can make \code{PLINK} files using the genotype matrix.
#' @details 
#' Currently, function \code{GRAB.SimuGMat} supports both unrelated and related subjects. 
#' Genotype data is simulated following Hardy-Weinberg Equilibrium with allele frequency ~ \code{runif(MinMAF, MaxMAF)}.
#' 
#' ## If \code{FamMode = "4-members"}
#' Total number of subjects is \code{nSub + 4 * nFam}. Each family includes 4 members with the family structure as below: 1+2->3+4.
#' 
#' ## If \code{FamMode = "10-members"}
#' Total number of subjects is \code{nSub + 10 * nFam}. Each family includes 10 members with the family structure as below: 1+2->5+6, 3+5->7+8, 4+6->9+10.
#' 
#' ## If \code{FamMode = "20-members"}
#' Total number of subjects is \code{nSub + 20 * nFam}. Each family includes 20 members with the family structure as below: 1+2->9+10, 3+9->11+12, 4+10->13+14, 5+11->15+16, 6+12->17, 7+13->18, 8+14->19+20.
#' 
#' @examples 
#' nSub = 100
#' nFam = 10
#' FamMode = "10-members"
#' nSNP = 10000
#' OutList = GRAB.SimuGMat(nSub, nFam, FamMode, nSNP)      
#' GenoMat = OutList$GenoMat
#' markerInfo = OutList$markerInfo
#' GenoMat[1:10,1:10]
#' head(markerInfo)
#' 
#' ## The following is to calculate GRM
#' MAF = apply(GenoMat, 2, mean)/2
#' GenoMatSD = t((t(GenoMat) - 2*MAF)/sqrt(2*MAF*(1-MAF)))
#' GRM = GenoMatSD %*% t(GenoMatSD) / ncol(GenoMat)
#' GRM1 = GRM[1:10, 1:10];
#' GRM2 = GRM[100+1:10, 100+1:10];
#' GRM1
#' GRM2
#' @export
GRAB.SimuGMat = function(nSub,
                         nFam,
                         FamMode,
                         nSNP,
                         MaxMAF = 0.5,
                         MinMAF = 0.05,
                         MAF = NULL)
{
  if(missing(FamMode) & missing(nFam)){
    cat("Since both 'FamMode' and 'nFam' are not specified, we only simulate genotype for unrelated subjects.\n")
    nFam = 0;
    FamMode = "Unrelated";
  }
  
  if(missing(nSub)){
    cat("Since 'nSub' is not specified, we only simulate genotype for family members.\n")
    nSub = 0;
  }
    
  if(!is.element(FamMode, c("Unrelated", "4-members", "10-members", "20-members")))
    stop("FamMode should be one of 'Unrelated', '4-members', '10-members', and '20-members'. Other input is not supported.")
  
  nSubInEachFam = 0;
  
  if(FamMode == "4-members"){
    nSubInEachFam = 4
    nHaploInEachFam = 4
    fam.mat = example.fam.4.members(nFam)
  }
    
  if(FamMode == "10-members"){
    nSubInEachFam = 10
    nHaploInEachFam = 8
    fam.mat = example.fam.10.members(nFam)
  }
  
  if(FamMode == "20-members"){
    nSubInEachFam = 20
    nHaploInEachFam = 16
    fam.mat = example.fam.20.members(nFam)
  }
  
  n = nSub + nFam * nSubInEachFam
  nHaplo = nFam * nHaploInEachFam
  
  if(n == 0){
    stop("Please give at least one of 'nSub' and 'nFam'.")
  }
  
  cat("Number of unrelated subjects:\t", nSub, "\n")
  cat("Number of families:\t", nFam, "\n")
  cat("Number of subjects in each family:\t", nSubInEachFam, "\n")
  cat("Number of all subjects:\t", n, "\n")
  
  if(is.null(MAF)){
    MAF = runif(nSNP, MinMAF, MaxMAF)
  }else{
    if(length(MAF) != nSNP)
      stop("length(MAF) should equal to nSNP.")
    cat("Since argument 'MAF' is given, arguments of 'MaxMAF' and 'MinMAF' are ignored.\n")
  }
  
  SNP.info = make.SNP.info(nSNP, MAF)
  
  GenoMat1 = GenoMat2 = NULL
    
  if(nHaplo != 0){
    cat("Simulating haplotype data for related subjects....\n")
    haplo.mat = haplo.simu(nHaplo, SNP.info) 
    
    cat("Simulating genotype data for related subjects....\n")
    GenoMat1 = from.haplo.to.geno(haplo.mat, fam.mat)    # output of example.fam(): n x 5 where n is sample size
  }
  
  if(nSub != 0){
    cat("Simulationg Genotype data for unrelated subjects....\n")
    GenoMat2 = geno.simu(nSub, SNP.info)
  }

  GenoMat = rbind(GenoMat1, GenoMat2)
  
  return(list(GenoMat = GenoMat,
              markerInfo = SNP.info))
}

### an example of family structure including 20 memebers in each family
## 1+2->9+10; 3+9->11+12; 4+10->13+14; 5+11->15+16; 6+12->17; 7+13->18; 8+14->19+20 
example.fam.20.members = function(n.fam)           # family numbers
{
  m = 20  # family members in each family
  fam.mat = c()
  c.h = 0 # count of haplotype
  for(i in 1:n.fam){
    FID = paste0("f",i)
    IID = paste0(FID,"_",1:m)
    Role = c(rep("Founder",8),
             rep("Offspring",12))
    Source1 = c(paste0("haplo-",c.h+1:8),
                IID[c(1,1,3,3,4,4,5,5,6,7,8,8)])
    Source2 = c(paste0("haplo-",c.h+9:16),
                IID[c(2,2,9,9,10,10,11,11,12,13,14,14)])
    c.h = c.h+16
    fam.mat = rbind(fam.mat,
                    cbind(FID,IID,Role,Source1,Source2))
  }
  fam.mat = data.frame(fam.mat,stringsAsFactors = F)
  return(fam.mat)  # five columns of FID, IID, Role, Source1, and Source2
}

### an example of family structure including 10 memebers in each family
## 1+2->5+6; 3+5->7+8; 4+6->9+10 
example.fam.10.members = function(n.fam)           # family numbers
{
  m = 10  # family members in each family
  fam.mat=c()
  c.h = 0 # count of haplotype
  for(i in 1:n.fam){
    FID = paste0("f",i)
    IID = paste0(FID,"_",1:m)
    Role = c(rep("Founder",4),
             rep("Offspring",6))
    Source1 = c(paste0("haplo-",c.h+1:4),
                IID[c(1,1,3,3,4,4)])
    Source2 = c(paste0("haplo-",c.h+5:8),
                IID[c(2,2,5,5,6,6)])
    c.h = c.h+8
    fam.mat = rbind(fam.mat,
                    cbind(FID,IID,Role,Source1,Source2))
  }
  fam.mat = data.frame(fam.mat,stringsAsFactors = F)
  return(fam.mat)  # five columns of FID, IID, Role, Source1, and Source2
}

### an example of family structure including 4 memebers in each family
## 1+2->3+4
example.fam.4.members = function(n.fam)           # family numbers
{
  m = 4  # family members in each family
  fam.mat=c()
  c.h = 0 # count of haplotype
  for(i in 1:n.fam){
    FID = paste0("f",i)
    IID = paste0(FID,"_",1:m)
    Role = c(rep("Founder",2),
             rep("Offspring",2))
    Source1 = c(paste0("haplo-",c.h+1:2),
                IID[c(1,1)])
    Source2 = c(paste0("haplo-",c.h+3:4),
                IID[c(2,2)])
    c.h = c.h+4
    fam.mat = rbind(fam.mat,
                    cbind(FID,IID,Role,Source1,Source2))
  }
  fam.mat = data.frame(fam.mat,stringsAsFactors = F)
  return(fam.mat)  # five columns of FID, IID, Role, Source1, and Source2
}

## Make a data.frame with each row for one SNPs. Column 1: SNP id; Column 2: beta.g; Column 3: maf 
## Example: SNP.info = make.SNP.info.example(1000,1000,0.2,0.3)
make.SNP.info = function(nSNP,  # number of null SNPs
                         MAF)   # fixed value of MAF for all SNPs
{
  SNP.info = data.table::data.table(SNP=paste0("SNP_",1:nSNP),
                                    MAF=MAF,
                                    stringsAsFactors = F)
  return(SNP.info)
}

geno.simu = function(nSub, SNP.info)
{
  nSNPs = nrow(SNP.info)
  MAFs = SNP.info$MAF
  GenoMat = sapply(MAFs, FUN = function(x){rbinom(nSub, 2, x)})
  if(nSub == 1)
    GenoMat = matrix(GenoMat, 1, nSNPs)
  colnames(GenoMat) = SNP.info$SNP
  rownames(GenoMat) = paste0("Subj-",1:nSub)
  return(GenoMat)  # matrix of m x n, where m is number of SNPs, n is number of subjects
}

## haplotype simulation for all founders
haplo.simu = function(n.haplo,   # number of haplotypes
                      SNP.info)  # a number or a vector with length of n.SNPs
{
  # check input
  n.SNPs = nrow(SNP.info)
  MAFs = SNP.info$MAF
  
  # generate haplotype matrix with each row for one haplotype and each column for one SNP
  haplo.mat = sapply(MAFs, 
                     FUN = function(x){
                       rbinom(n.haplo,1,x)
                     })
  colnames(haplo.mat) = SNP.info$SNP
  rownames(haplo.mat) = paste0("haplo-",1:n.haplo)
  return(haplo.mat)   # matrix of m x p, where m is number of SNPs, p is number of founders
}

### genotype simulation based on haplotype
from.haplo.to.geno = function(haplo.mat,  # output of haplo.simu():  m x p where m is number of confounder alleles, and p is number of SNPs
                              fam.mat)    # output of example.fam(): n x 5 where n is sample size
{
  n = nrow(fam.mat)     # number of subjects
  m = ncol(haplo.mat)   # number of SNPs
  Haplo1.mat = matrix(nrow = n, ncol = m)
  Haplo2.mat = matrix(nrow = n, ncol = m)
  rownames(Haplo1.mat) = rownames(Haplo2.mat) = fam.mat$IID
  for(i in 1:n){  # cycle for all subjects
    if(i %% 1000 == 0) print(paste0("Complete Genotype Simulation for ",i," Subjects."))
    Role = fam.mat$Role[i]
    S1 = fam.mat$Source1[i]
    S2 = fam.mat$Source2[i]
    if(Role == "Founder"){  # directly extract from haplotype data
      Haplo1.mat[i,] = haplo.mat[S1,]
      Haplo2.mat[i,] = haplo.mat[S2,]
    }
    if(Role == "Offspring"){ # pass from Founders
      ## Haplotype 1 is from Founder S1, randomly selected from two haplotypes of Founder S1.
      S1.pass.H1=rbinom(m,1,0.5) 
      S1.pass.H2=1-S1.pass.H1
      Haplo1.mat[i,] = Haplo1.mat[S1,]*S1.pass.H1+Haplo2.mat[S1,]*S1.pass.H2
      ## Haplotype 2 is from Founder S2, randomly selected from two haplotypes of Founder S2.
      S2.pass.H1=rbinom(m,1,0.5) 
      S2.pass.H2=1-S2.pass.H1
      Haplo2.mat[i,] = Haplo1.mat[S2,]*S2.pass.H1+Haplo2.mat[S2,]*S2.pass.H2
    }
  }
  Geno.mat=Haplo1.mat+Haplo2.mat
  colnames(Geno.mat)=colnames(haplo.mat)
  # Geno.mat=data.frame(Geno.mat, stringsAsFactors = F)
  return(Geno.mat)
}

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


#' Make PLINK files using a numeric R matrix
#' 
#' Make PLINK files using a numeric matrix \code{GenoMat} (0,1,2,-9), \code{rownames(GenoMat)} are subject IDs and \code{colnames(GenoMat)} are marker IDs
#' 
#' @param GenoMat a numeric \code{n*m} genotype matrix (0,1,2,-9). Each row is for one subject and each column is for one marker. Row names of subject IDs and column names of marker IDs are required.
#' @param OutputPrefix a character, prefix of the PLINK files to output (including path).
#' @param A1 a character to specify allele 1 (*default="G"*), usually minor (ALT).
#' @param A2 a character to specify allele 2 (*default="A"*), usually major (REF).
#' @param CHR a character vector of the chromosome numbers for all markers. *Default=NULL*, that is, \code{CHR=rep(1, m)}.
#' @param BP a numeric vector of the base positions for all markers. *Default=NULL*, that is, \code{BP=1:m)}.
#' @param Pheno a character vector of the phenotypes for all subjects. *Default=NULL*, that is, \code{Pheno=rep(-9, n)}.
#' @param Sex a numeric vector of the sex for all subjects. *Default=NULL*, that is, \code{Sex=rep(1, n))}.
#' @return \code{PLINK} text files (PED and MAP) are stored in 'OutputPrefix'. Suppose A1 is "G" and A2 is "A", then genotype of 0,1,2,-9 will be coded as "GG", "AG", "AA", "00". If PLINK binary files (BED, BIM, and FAM) are required, please download PLINK software and use option of "--make-bed".
#' Please check \code{Details} section for the downstream process.
#' @details 
#' Check [link](https://www.cog-genomics.org/plink/2.0/) for detailed information of \code{PLINK} 2.00 alpha. Check [link](https://enkre.net/cgi-bin/code/bgen/doc/trunk/doc/wiki/bgenix.md) for detailed information of \code{bgenix} tool.
#' ## Convert PLINK text files to binary files
#' Run \code{plink --file simuPLINK --make-bed --out simuPLINK} to convert PLINK text files (MAP and PED) to binary files (BED, BIM, and FAM).
#' ## Convert PLINK binary files to raw files
#' Run \code{plink --bfile simuPLINK --recode A --out simuRAW} to convert PLINK binary files (BED, BIM, and FAM) to raw files (raw).
#' ## Convert PLINK binary files to bgen files
#' RUN \code{plink2 --bfile simuPLINK --export bgen-1.2 bits=8 ref-first --out simuBGEN} to convert PLINK binary files (BED, BIM, and FAM) to BGEN binary files (BGEN).
#' ## Make bgi file using \code{bgenix} tool
#' RUN \code{bgenix -g simuBGEN.bgen --index}
#' @examples 
#' ### Step 1: simulate a numeric genotype matrix
#' n = 1000
#' m = 20
#' MAF = 0.3
#' set.seed(123)
#' GenoMat = matrix(rbinom(n*m, 2, MAF), n, m)
#' rownames(GenoMat) = paste0("Subj-",1:n)
#' colnames(GenoMat) = paste0("SNP-",1:m)
#' outputDir = system.file("results", package = "GRAB")
#' outputPrefix = paste0(outputDir, "/simuPLINK")
#' 
#' ### Step 2(a): make PLINK files without missing genotype
#' GRAB.makePlink(GenoMat, outputPrefix)
#' 
#' ### Step 2(b): make PLINK files with genotype missing rate of 0.1
#' indexMissing = sample(n*m, 0.1*n*m)
#' GenoMat[indexMissing] = -9
#' GRAB.makePlink(GenoMat, outputPrefix)
#'      
#' ## The following are in shell environment
#' # plink --file simuPLINK --make-bed --out simuPLINK
#' # plink --bfile simuPLINK --recode A --out simuRAW 
#' # plink2 --bfile simuPLINK --export bgen-1.2 bits=8 ref-first --out simuBGEN  # UK Biobank use 'ref-first'
#' # bgenix -g simuBGEN.bgen --index
#'      
#' @export
GRAB.makePlink = function(GenoMat,
                          OutputPrefix,
                          A1 = "G",
                          A2 = "A",
                          CHR = NULL,      
                          BP = NULL,
                          Pheno = NULL,
                          Sex = NULL)
{
  if(!is.numeric(GenoMat))
    stop("'GenoMat' should be a numeric matrix.")
  
  if(any(!unique(as.numeric(GenoMat)) %in% c(0, 1, 2, -9)))
    stop("'GenoMat' should only include elements of 0, 1, 2, -9.")
  
  SNP = colnames(GenoMat)
  FID = IID = rownames(GenoMat)
  
  if(is.null(SNP) | is.null(IID))
    stop("rownames and colnames of GenoMat should be specified.")
  
  m = length(SNP)
  n = length(IID)
  
  cat("number of markers:\t", m, "\n")
  cat("number of samples:\t", n, "\n")
  
  if(is.null(Pheno)){
    Pheno = rep(-9, n)
  }else{
    if(length(Pheno != n))
      stop("length(Pheno) should be the same as nrow(GenoMat).")  
  }
  
  if(is.null(Sex)){
    Sex = rep(1, n)
  }else{
    if(length(Sex) != n)
      stop("length(Sex) should be the same as nrow(GenoMat).") 
  }
  
  if(is.null(BP)){
    BP = 1:m
  }else{
    if(length(BP) != m)
      stop("length(BP) should be the same as ncol(GenoMat).")
  }

  if(is.null(CHR)){
    CHR = rep(1, m)
  }else{
    if(length(CHR) != m)
      stop("length(CHR) should be the same as ncol(GenoMat).")
  }
  
  PED = cbind(FID=FID, 
              IID=IID, 
              PID=0, 
              MID=0, 
              Sex = Sex, 
              Phen = Pheno);
  
  MAP = cbind(CHR = CHR, 
              SNP = SNP, 
              GeneDist = 0,
              BP = BP);
  
  Geno.ped1 = Geno.ped2 = ifelse(GenoMat == -9, "0", A1)
  Geno.ped1 = ifelse(GenoMat>=1, A2, Geno.ped1)
  Geno.ped2 = ifelse(GenoMat>=2, A2, Geno.ped2)
  
  Geno.ped = matrix(nrow = n, ncol = 2*m)
  Geno.ped[,seq(1,2*m,2)] = Geno.ped1;
  Geno.ped[,seq(2,2*m,2)] = Geno.ped2;
  
  PED = cbind(PED, Geno.ped)
  
  MAP.file = paste0(OutputPrefix, ".map")
  PED.file = paste0(OutputPrefix, ".ped")
  # OUT.file = paste0(work.dir, "/", out.prefix)
  
  data.table::fwrite(MAP, MAP.file, quote = F, col.names = F, row.names = F, sep = " ")
  data.table::fwrite(PED, PED.file, quote = F, col.names = F, row.names = F, sep = " ")
  
  cat("Working directory:\t", getwd(), "\n")
  cat("PED file:\t", PED.file, "\n")
  cat("MAP file:\t", MAP.file, "\n")
  
  message = paste0("PLINK files have been saved to ", OutputPrefix, ".")
  return(message)
}

#' Simulate phenotype using linear predictor \code{eta}
#' 
#' \code{GRAB} package can help simulate a wide variaty of phenotypes
#' 
#' @param eta linear predictors, usually covar x beta.covar + genotype x beta.genotype 
#' @param traitType "quantitative", "binary", "ordinal", or "time-to-event"
#' @param control a list of parameters for controlling the simulation process
#' @return a numeric vector of phenotype
#' @export
GRAB.SimuPheno = function(eta, 
                          traitType = "binary", 
                          control = list(pCase = 0.1,
                                         sdError = 1,
                                         pEachGroup = c(1,1,1)))
{
  if(!traitType %in% c("quantitative", "binary", "ordinal", "time-to-event"))
    stop('traitType is limited to "quantitative", "binary", "ordinal", and "time-to-event".')
  
  if(traitType == "binary")
    if(!"pCase" %in% names(control))
      stop("For binary phenotype, argument 'control' should include 'pCase' which is the proportion of cases.")
  
  if(traitType == "quantitative")
    if(!"sdError" %in% names(control))
      cat("For quantitative phenotype, argument 'control' should include 'sdError' which is the stardard derivation of the error term.")
  
  if(traitType == "ordinal")
    if(!"pEachGroup" %in% names(control))
      cat("For quantitative phenotype, argument 'control' should include 'pEachGroup' which is ratio of sample size in each group.")
  
  eta = eta - mean(eta)
  n = length(eta)
  
  seed = sample(1e9,1)
  cat("Random number seed:\t", seed, "\n")
  
  ### quantitative trait
  if(traitType == "quantitative"){
    sdError = control$sdError
    set.seed(seed)
    error = rnorm(n, sd = sdError)
    pheno = eta + error
    return(pheno)
  }
  
  ### binary trait
  if(traitType == "binary"){
    pCase = control$pCase
    eta0 = uniroot(f.binary, c(-100,100), eta = eta, pCase = pCase,seed = seed)
    eta0 = eta0$root
    set.seed(seed)
    eta.new = eta0 + eta
    mu = exp(eta.new) / (1 + exp(eta.new))   # The probability being a case given the covariates, genotypes, and addition effect
    pheno = rbinom(n, 1, mu)                     # Case-control status
    return(pheno)
  }
  
  ### quantitative trait
  if(traitType == "ordinal"){
    pEachGroup = control$pEachGroup
    Eps = getEps(pEachGroup, eta, seed)
    set.seed(seed)
    pheno.latent = runif(n)
    pheno = rep(0, n)
    for(g in 1:length(Eps)){
      mu = exp(Eps[g]-eta)/(1+exp(Eps[g]-eta))
      pheno[pheno.latent > mu] = g
    }
    return(pheno)
  }
}


#### lower function to estimate eta0 given a prevalence. Will be used in data.simu.binary().
f.binary = function(eta,               # Sample size
                    pCase,             # Prevalence
                    eta0,              # Intercept
                    seed)     
{
  set.seed(seed)
  n = length(eta)
  eta.new = eta0 + eta
  mu = exp(eta.new) / (1 + exp(eta.new))   # The probability being a case given the covariates, genotypes, and addition effect
  Y = rbinom(n, 1, mu)                     # Case-control status
  re = mean(Y) - pCase
  return(re)
}

#### lower function to estimate epsilons given a ratio(phenotypic distribution). Will be used in data.simu.categorical().
getProb = function(eps,
                   eta.true,
                   prob,
                   seed)
{
  set.seed(seed)
  n = length(eta.true)
  mu = exp(eps-eta.true) / (1+exp(eps-eta.true))
  Y.latent = runif(n)
  diffprob = mean(Y.latent < mu) - prob
  return(diffprob)
}

getEps = function(ratios,
                  eta.true,
                  seed)
{
  sumR = sum(ratios)
  cumR = 0
  J = length(ratios)
  Eps = c()
  for(i in 1:(J-1)){
    cumR = cumR + ratios[i]
    eps = uniroot(getProb, c(-100,100), eta.true = eta.true, prob = cumR/sumR, seed = seed)
    Eps = c(Eps, eps$root)
  }
  return(Eps)
}