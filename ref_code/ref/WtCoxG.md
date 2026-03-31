```R
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
  residual = obj.fit$residuals,
  weight = weight,
  indicator = PhenoData$SurvEvent
)
resid_file = tempfile()
data.table::fwrite(dt, "examples/simuResid.txt", sep = "\t", col.names = FALSE)
```

```sh
build/grab \
  --method WtCoxG \
  --null-resid examples/simuResid.txt \
  --bfile examples/simuPLINK \
  --ref-af examples/simuRefAf_1pop.txt \
  --sparse-grm examples/SparseGRM.txt \
  --prevalence 0.1 \
  --out tmp/WtCoxG_output.txt \
  --threads 3
```
[ERROR] examples/simuRefAf_1pop.txt line 14: missing N field