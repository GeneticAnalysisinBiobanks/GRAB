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
    indicator = Indicator_i,
    residual  = residual_i,
    weight    = weight_i
  )
  data.table::fwrite(dt, resid_file_lst[i], sep = "\t", col.names = FALSE)
}


LEAF is a method related to WtCoxG. It cluaters the sample into several clusters based on PCs, and then runs WtCoxG in each cluster.

Commond line for LEAF is as follows:

build/grab \
  --method LEAF \
  --resid-file examples/simuResid1.txt,examples/simuResid2.txt,examples/simuResid3.txt \
  --bfile examples/simuPLINK \
  --ref-af-file examples/simuRefAf_2pop.txt \
  --sparse-grm-file examples/SparseGRM.txt \
  --ref-prevalence 0.1 \
  --output-file tmp/LEAF_output.txt \
  --nsnp-per-chunk 256 \
  --nthread 4

--resid-file same as wtcoxg, but with N files for N clusters.
--bfile is the prefix of PLINK binary fileset. The .bed, .bim, and .fam files should be in the same directory.
--ref-af-file has 6+2N columns (is diff from old format): with the first 6 columns being the same as .bim file,
then allele frequency of A1 and allele sample size of the reference population for each of the N reference populations.
No header is needed. Lines start with # are comments and will be ignored.

The work flow is as follows:
1. Parse bfiles and calculate a1 allele frequencies of each cluster
2. Parse RefAfFile2 and get a1 allele frequencies of each reference population
3. Match sample and reference AF by chr, bp, a1, a2; only exactly match markers will be used.
   Do not flip alleles in either sample or reference; if there is no match, the marker will be dropped.
4. Run summix for each cluster to estimate ancestry proportions.
5. Test batch effects
6. Perform marker-level tests by WtCoxG
7. CCT to combine p values across clusters

Except for step 6, the code needs to be freshly written in C++. Please read R/LEAF.R and src/LEAF.h in develop branch to understand the details.
I know nloptr::slsqp is complex to implement in C++. We can use a simple method first. Ask me if you have any questions. Also tell me more details you plan to work.


1. It seems you parse --ref-af-file by R logic. The file format is different now.
--ref-af-file has 6+2N columns (examples/simuRefAf_2pop.txt): with the first 6 columns being the same as .bim file,
then, pairs of allele frequency and allele sample size for each of the N reference populations.
No header is needed. Lines start with # are comments and will be ignored. 
So no need to parse *_AF/*_AN.
No need to do canonical allele matching. No need to flip alleles. The allele frequency is always for A1.
If the first 6 columns does not matche .bim, drop the marker.

2. Remove "Phase 1,2,3:" from loggings.

3. please add --nsnp-per-chunk option to the command line, which control the number of markers processed in each chunk in engine/marker.cpp.
Default to 8192, min to the number of markers read by plink curser once (is 256?).

4. how do you parse bed file now? are you parse just once?

5. For the summix, constrain nclusters <=6, so enumerate all active sets (K <= 6). For each set, solve equality-constrained least squares 
using a small KKT system (Lagrange multiplier). is this good?
