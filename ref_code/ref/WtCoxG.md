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
  --resid-file examples/simuResid.txt \
  --bfile examples/simuPLINK \
  --ref-af-file examples/simuRefAf.txt \
  --sparse-grm-file examples/SparseGRM.txt \
  --ref-prevalence 0.1 \
  --output-file tmp/WtCoxG_output.txt \
  --nthread 3
```
