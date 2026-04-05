# SAGELD

SAGELD (Scalable and Accurate algorithm for Gene-Environment interaction analysis
using Longitudinal Data) extends SPAGRM to support G×E interaction testing for
related samples in large-scale biobanks.

## Residual file format

```
#IID  R_G  R_<E1>  R_Gx<E1>  [R_<E2>  R_Gx<E2>  ...]
```

- `R_G` — main-effect residual: sum of outcome model residuals per subject
- `R_<E>` — environment residual: sum of E-model residuals per subject (e.g. `R_AGE`)
- `R_Gx<E>` — G×E residual: weighted sum `sum(r * E)` per subject (e.g. `R_GxAGE`)

Column names must follow this pattern when a header is provided; otherwise columns
are assigned by position and environments are named E1, E2, ...

## Generating the residual file (R)

```r
library(data.table)

LongPheno <- fread("examples/simuLongPHENO.txt")

# Fit null models
fit          <- lme4::lmer(Y ~ AGE + GENDER + (AGE | IID), data = LongPheno)
fit_E_age    <- lme4::lmer(AGE    ~ (1 | IID), data = LongPheno)
fit_E_gender <- lme4::lmer(GENDER ~ (1 | IID), data = LongPheno)

df <- data.table(
  IID        = LongPheno$IID,
  r          = residuals(fit),
  r_E_age    = residuals(fit_E_age),
  r_E_gender = residuals(fit_E_gender),
  AGE        = LongPheno$AGE,
  GENDER     = LongPheno$GENDER
)

agg <- df[, .(
  R_G       = sum(r),
  R_AGE     = sum(r_E_age),
  R_GxAGE   = sum(r * AGE),
  R_GENDER  = sum(r_E_gender),
  R_GxGENDER = sum(r * GENDER)
), by = .(`#IID` = IID)]

# Single environment (AGE only)
fwrite(
  agg[, .(`#IID`, R_G, R_AGE, R_GxAGE)],
  "examples/simuResid_SAGELD1.txt", sep = "\t"
)

# Two environments (AGE + GENDER)
fwrite(
  agg[, .(`#IID`, R_G, R_AGE, R_GxAGE, R_GENDER, R_GxGENDER)],
  "examples/simuResid_SAGELD2.txt", sep = "\t"
)
```

## Running SAGELD

```sh
# Single environment
build/grab \
  --method SAGELD \
  --null-resid examples/simuResid_SAGELD1.txt \
  --sp-grm-grab examples/SparseGRM.txt \
  --pairwise-ibd examples/PairwiseIBD.txt \
  --bfile examples/simuPLINK \
  --out tmp/SAGELD_output1.txt

# Two environments (single output file)
build/grab \
  --method SAGELD \
  --null-resid examples/simuResid_SAGELD2.txt \
  --sp-grm-grab examples/SparseGRM.txt \
  --pairwise-ibd examples/PairwiseIBD.txt \
  --bfile examples/simuPLINK \
  --out tmp/SAGELD_output2.txt
```

## Output columns

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
P_G  P_Gx<E1>  [P_Gx<E2>  ...]  Z_G  Z_Gx<E1>  [Z_Gx<E2>  ...]
```

- `P_G` — main genetic effect p-value (normal approximation)
- `P_Gx<E>` — G×E interaction p-value (saddlepoint approximation)
- `Z_G`, `Z_Gx<E>` — corresponding z-scores
