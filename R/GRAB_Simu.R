

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
#' @import SKAT, data.table

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
