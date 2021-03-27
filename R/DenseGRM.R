
#' Set up a dense GRM (only for developers)
#'
#' Set up a dense GRM (only for developers), other users can ignore this function
#' @param GenoFile a character of genotype file. Three types of genotype files are supported: PLINK ("prefix.bed"), BGEN ("prefix.bgen"), and VCF ("prefix.vcf" or "prefix.vcf.gz"). 
#' @param GenoFileIndex additional index files corresponding to the "GenoFile". If Null (default), the same prefix as GenoFile is used. PLINK: c("prefix.bim", "prefix.fam"), BGEN: c("prefix.bgi"), and VCF: c("prefix.vcf.tbi") or c("prefix.vcf.gz.tbi").
#' @param subjData a character vector of subject IDs. Its order should be the same as the subjects order in the formula and data. 
#' @return no result is returned
#' @examples
#' # Check ?getDenseGRM() for an example.
#' @export

setDenseGRM = function(GenoFile, GenoFileIndex = NULL, subjData = NULL)
{
  genoList = setGenoInput(GenoFile, GenoFileIndex, subjData)   # check Geno.R for more details
  
  if(genoList$genoType != "PLINK")
    stop("Only PLINK format is supported")
  
  memoryChunk = 2 # (GB)
  minMafGRM = 0.01
  maxMissingGRM = 0.1
  
  setDenseGRMInCPP(memoryChunk, minMafGRM, maxMissingGRM)
  return(genoList)
}

#' Suppose that a dense GRM is Phi and input is bVec, return Phi * bVec (only for developers)
#'
#' Suppose that a dense GRM is Phi and input is bVec, return Phi * bVec (only for developers), users can simply ignore this function
#' @param bVec a numeric vector with the same length as in subjData (check the input of \code{setDenseGRM})
#' @return a numeric vector of Phi * bVec
#' @examples
#' # set up the dense GRM in C++
#' GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
#' famData = read.table(gsub("bed","fam", GenoFile))
#' subjData = famData$V2
#' genoList = setDenseGRM(GenoFile, subjData = subjData)
#' 
#' set.seed(1)
#' bVec = rnorm(1000)
#' KinbVec = getDenseGRM(bVec)
#' 
#' # The following is based on the definition of GRM to validate the DenseGRM object
#' PlinkFile = system.file("extdata", "example.bed", package = "GRAB")
#' IDsToIncludeFile = system.file("extdata", "example.IDsToIncludeFile.txt", package = "GRAB")
#' GenoList = GRAB.ReadGeno(PlinkFile, control = list(IDsToExcludeFile = IDsToIncludeFile))
#' GenoMat = GenoList$GenoMat
#' markerInfo = GenoList$markerInfo
#' pos = which(markerInfo$CHROM != 1)
#' GenoMat = GenoMat[,pos]
#' MAF = apply(GenoMat, 2, mean)/2
#' stdGenoMat = (t(GenoMat) - 2*MAF) / sqrt(2*MAF*(1-MAF)) / sqrt(ncol(GenoMat))
#' KinMat = t(stdGenoMat) %*% stdGenoMat
#' KinbVec1 = KinMat %*% bVec
#' plot(KinbVec, KinbVec1)
#' head(cbind(KinbVec, KinbVec1))
#' 
#' @export

getDenseGRM = function(bVec)
{
  excludeChr = "1"
  grainSize = 1
  
  KinbVec = getDenseGRMInCPP(bVec, excludeChr, grainSize)
  return(KinbVec)
}


