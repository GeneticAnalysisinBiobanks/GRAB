calRegrWeight <- function(Indicator, RefPrevalence) {
  sample_ratio <- sum(Indicator) / sum(1 - Indicator)
  population_ratio <- RefPrevalence / (1 - RefPrevalence)
  weight <- rep(1, length(Indicator))
  weight[Indicator == 0] <- sample_ratio / population_ratio
  return(weight)
}

Ncluster = 3
RefPrevalence <- 0.1
GenoFile = "inst/extdata/simuPLINK.bed"
SparseGRMFile <- "inst/extdata/SparseGRM.txt"
RefAfFile2 <- "inst/extdata/simuRefAf_2pop.txt"

PhenoData <- data.table::fread("inst/extdata/simuPHENO.txt")
subjData <- PhenoData$IID
indicator <- as.integer(PhenoData$SurvEvent)
PCmatrix <- as.matrix(PhenoData[, c("PC1", "PC2", "PC3", "PC4")])
designMat <- as.matrix(PhenoData[, c("AGE", "GENDER", "PC1", "PC2", "PC3", "PC4")])
clusterIdx <- kmeans(PCmatrix, centers = Ncluster, nstart = 25)$cluster

resid_file_lst <- character(Ncluster)
for (i in seq_len(Ncluster)) {
  idx_i <- which(clusterIdx == i)
  Indicator_i <- as.integer(indicator[idx_i])
  weight_i <- calRegrWeight(Indicator_i, RefPrevalence)
  
  residual_i <- glm.fit(
    designMat[idx_i, , drop = FALSE],
    Indicator_i,
    weights = weight_i,
    family = binomial(),
    intercept = TRUE
  )$residuals
  
  dt <- data.table(
    subject   = PhenoData$IID[idx_i],
    indicator = Indicator_i,
    residual  = residual_i,
    weight    = weight_i
  )
  resid_file_lst[i] <- tempfile()
  fwrite(dt, resid_file_lst[i], sep = "\t", col.names = FALSE)
}


To refactor the code to run LEAF as a single call cppLEAF from residuals (write new cpp functions):
1. Calculate a1 allele frequencies of each cluster (AF_cluster) from GenoFile
2. Parse a1 allele frequencies from RefAfFile2
3. RefAfFile2 is moving to new format, with the first 6 columns being the same as .bim file, 
   and then AF AN paires for each reference population. 
   No header is needed. Lines start with # are comments and will be ignored.
4. Match sample and reference AF by chr, bp, a1, a2; only exactly match markers will be used.
   Do not flip alleles in either sample or reference; if there is no match, the marker will be dropped.
5. Make summix input matrix with columns AF_cluster, AF_ref1, AF_ref2, ... AF_refN
6. run summix for each cluster by calling R function nloptr::slsqp from cpp
7. New cpp functions to test batch effects
8. Perform marker-level tests with current workflow.


OutputFile = tempfile()
cppLEAF(
  ResidFile = resid_file_lst,
  bfile = "inst/extdata/simuPLINK",
  RefAfFile = "inst/extdata/simuRefAf_2pop.txt",
  SparseGRMFile = "inst/extdata/SparseGRM.txt",
  RefPrevalence = 0.1,
  OutputFile = OutputFile,
  nthread = 3
)
data.table::fread(OutputFile)
