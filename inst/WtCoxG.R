calRegrWeight <- function(Indicator, RefPrevalence) {
  sample_ratio <- sum(Indicator) / sum(1 - Indicator)
  population_ratio <- RefPrevalence / (1 - RefPrevalence)
  weight <- rep(1, length(Indicator))
  weight[Indicator == 0] <- sample_ratio / population_ratio
  return(weight)
}

# get residuals from weighted cox regression and save to a file
PhenoData <- data.table::fread("inst/extdata/simuPHENO.txt")
weight <- calRegrWeight(PhenoData$SurvEvent, RefPrevalence)
obj.fit <- survival::coxph(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData, weight = weight, robust = TRUE
)
dt <- data.table(
  subject = PhenoData$IID,
  indicator = PhenoData$SurvEvent,
  residual = obj.fit$residuals,
  weight = weight
)
resid_file = tempfile()
data.table::fwrite(dt, resid_file, sep = "\t", col.names = FALSE)


To refactor the code to run WtCoxG as a single call cppWtCoxG from residuals:
1. Write cpp functions to test batch effects
2. RefAfFile is moving to new format (\s+ deliminated), with the first 6 columns being the same as .bim file,
then allele frequency of A1 and allele sample size of the reference population.
Match sample and reference markers by chr, bp, a1, a2; only exactly match markers will be used.
Do not flip alleles in either sample or reference; if there is no match, the marker will be dropped.
No header is needed. Lines start with # are comments and will be ignored.
3. The ResidFile (\s+ deliminated) will have 4 columns: subject, indicator, residual, weight.
No header is needed. Lines start with # are comments and will be ignored.

OutputFile = tempfile()
cppWtCoxG(
  ResidFile = resid_file,
  bfile = "inst/extdata/simuPLINK",
  RefAfFile = "inst/extdata/simuRefAf.txt",
  SparseGRMFile = "inst/extdata/SparseGRM.txt",
  RefPrevalence = 0.1,
  OutputFile = OutputFile,
  nthread = 3
)
data.table::fread(OutputFile)
