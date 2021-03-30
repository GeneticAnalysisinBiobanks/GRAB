
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


#' #' GRAB: simulate genotype matrix based on family structure
#' #' 
#' #' Simulate genotype matrix based on family structure
#' #' 
#' #' @param n.fam number of families in simulation
#' #' @param n.SNPs number of SNPs in simulation
#' #' @param fam.kin a matrix of GRM (kinship matrix)
#' #' @param min.MAF minimal MAF in simulation
#' #' @param max.MAF maximal MAF in simulation
#' #' @return a matrix of genotype
#' #' @examples 
#' #' fam.kin = read.table(system.file("extdata", "example_10members.kin.txt", package = "GRAB"))
#' #' fam.kin = as.matrix(fam.kin)
#' #' n.fam = 100
#' #' n.SNPs = 1000
#' #' min.MAF = 0.05
#' #' max.MAF = 0.3
#' #' GMat = GRAB.SimuGMat(n.fam, n.SNPs, fam.kin, min.MAF, max.MAF)
#' #'      
#' #' @export
#' GRAB.SimuGMat = function(n.fam, n.SNPs, fam.kin, min.MAF, max.MAF)
#' {
#'   ##
#' }




