The current WtCoxG workflow is that, users calculate residuals and save to file.

```R
calRegrWeight <- function(Indicator, RefPrevalence) {
  sample_ratio <- sum(Indicator) / sum(1 - Indicator)
  population_ratio <- RefPrevalence / (1 - RefPrevalence)
  weight <- rep(1, length(Indicator))
  weight[Indicator == 0] <- sample_ratio / population_ratio
  return(weight)
}
RefPrevalence <- 0.1
PhenoData <- data.table::fread("examples/simuPHENO.txt")
indicator <- as.integer(PhenoData$SurvEvent)
weight <- calRegrWeight(indicator, RefPrevalence)


## WtCoxG
resid1 <- survival::coxph(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData, weight = weight
)$residuals

fit2 <- glm(
  indicator ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData, weight = weight,
  family = binomial()
)
resid2 <- residuals(fit2, type="response")

dt <- data.table::data.table(
  subject = PhenoData$IID,
  residual = resid1,
  weight = weight,
  indicator = indicator
)
data.table::fwrite(dt, "examples/simuResid.txt", sep = "\t", col.names = FALSE)
```

Then grab read residuals and perform tests.

```sh
build/grab \
  --method WtCoxG \
  --null-resid examples/simuResid.txt \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --prevalence 0.1 \
  --out tmp/WtCoxG_output.txt \
  --threads 3
```

Refactor grab to do all the above steps in one go. Keep the option to input residuals.

```sh
build/grab \
  --method WtCoxG \
  --pheno examples/simuPHENO.txt \
  --pheno-binary BinaryPheno \
  --covar-name AGE,GENDER,PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --prevalence 0.1 \
  --out tmp/WtCoxG_binary.txt \
  --threads 3

build/grab \
  --method WtCoxG \
  --pheno examples/simuPHENO.txt \
  --pheno-survival SurvTime:SurvEvent \
  --covar-name AGE,GENDER,PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --prevalence 0.1 \
  --out tmp/WtCoxG_survival.txt \
  --threads 3
```
