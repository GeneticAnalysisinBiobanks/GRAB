```R
calRegrWeight <- function(Indicator, RefPrevalence) {
  sample_ratio <- sum(Indicator) / sum(1 - Indicator)
  population_ratio <- RefPrevalence / (1 - RefPrevalence)
  weight <- rep(1, length(Indicator))
  weight[Indicator == 0] <- sample_ratio / population_ratio
  return(weight)
}

Ncluster = 3
RefPrevalence <- 0.1
GenoFile = "examples/simuPLINK.bed"
SparseGRMFile <- "examples/SparseGRM.txt"
RefAfFile2 <- "examples/simuRefAf_2pop.txt"

PhenoData <- data.table::fread("examples/simuPHENO.txt")
subjData <- PhenoData$IID
indicator <- as.integer(PhenoData$SurvEvent)
PCmatrix <- as.matrix(PhenoData[, c("PC1", "PC2", "PC3", "PC4")])
designMat <- as.matrix(PhenoData[, c("AGE", "GENDER", "PC1", "PC2", "PC3", "PC4")])
clusterIdx <- kmeans(PCmatrix, centers = Ncluster, nstart = 25)$cluster

resid_file_lst <- c("examples/simuResid1.txt","examples/simuResid2.txt","examples/simuResid3.txt")
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
  
  dt <- data.table::data.table(
    subject   = PhenoData$IID[idx_i],
    residual  = residual_i,
    weight    = weight_i,
    indicator = Indicator_i
  )
  data.table::fwrite(dt, resid_file_lst[i], sep = "\t", col.names = FALSE)
}
```

```sh
build/grab \
  --method LEAF \
  --null-resid examples/simuResid1.txt,examples/simuResid2.txt,examples/simuResid3.txt \
  --bfile examples/simuPLINK \
  --ref-af examples/simuRefAf_2pop.txt \
  --sparse-grm examples/SparseGRM.txt \
  --prevalence 0.1 \
  --out tmp/LEAF_output.txt \
  --chunk-size 256 \
  --threads 4
```
