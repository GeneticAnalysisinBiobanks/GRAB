calRegrWeight <- function(Indicator, RefPrevalence) {
  sample_ratio <- sum(Indicator) / sum(1 - Indicator)
  population_ratio <- RefPrevalence / (1 - RefPrevalence)
  weight <- rep(1, length(Indicator))
  weight[Indicator == 0] <- sample_ratio / population_ratio
  return(weight)
}

# get residuals from weighted cox regression and save to a file
PhenoData <- data.table::fread("examples/simuPHENO.txt")
weight <- calRegrWeight(PhenoData$SurvEvent, 0.1)
obj.fit <- survival::coxph(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData, weight = weight, robust = TRUE
)
dt <- data.table::data.table(
  subject = PhenoData$IID,
  indicator = PhenoData$SurvEvent,
  residual = obj.fit$residuals,
  weight = weight
)
resid_file = tempfile()
data.table::fwrite(dt, "examples/simuResid.txt", sep = "\t", col.names = FALSE)


build/grab \
  --method WtCoxG \
  --resid-file examples/simuResid.txt \
  --bfile examples/simuPLINK \
  --ref-af-file examples/simuRefAf.txt \
  --sparse-grm-file examples/SparseGRM.txt \
  --ref-prevalence 0.1 \
  --output-file tmp/WtCoxG_output.txt \
  --nthread 3

--resid-file has 4 columns: subject, indicator, residual, weight. No header is needed. Lines start with '#' will be ignored.
--bfile is the prefix of PLINK binary fileset. The .bed, .bim, and .fam files should be in the same directory.
--ref-af-file has 8 columns (is diff from old format): with the first 6 columns being the same as .bim file,
then allele frequency of A1 and allele sample size of the reference population.
No header is needed. Lines start with # are comments and will be ignored.

The work flow is as follows:
1. (The code needs freshly written in C++). Parse resid-file, ref-af-file and bfile. Match sample (bfile) and reference (ref-af-file) markers by chr, bp, a1, a2; only exactly matched markers are kept.
Do not flip alleles in either sample or reference; if there is no match, the marker will be dropped (is diff from old workflow).
2. (The code needs freshly written in C++). Test batch effects. Please refer to R/WtCoxG.R in develop branch.
3. Perform marker-level tests. This part can be refactered from Reference/src/mtWtCoxG.cpp.

Please read the related files and fully understand the algorithm. Ask me if you have any questions. Also tell me more details you plan to work.
