
```R
PhenoData <- data.table::fread("examples/simuPHENO.txt")
PhenoData$OrdinalPheno <- factor(PhenoData$OrdinalPheno, levels = c(0, 1, 2))

obj.POLMM <- GRAB.NullModel(
 OrdinalPheno ~ AGE + GENDER + PC1 + PC2,
 data = PhenoData,
 subjIDcol = "IID",
 method = "POLMM",
 traitType = "ordinal",
 GenoFile = "examples/simuPLINK.bed",
 SparseGRMFile = "examples/SparseGRM.txt"
)

GRAB.Marker(obj.POLMM, GenoFile, "tmp/POLMM_output.txt")
```

The above code performs POLMM by the R code in branch develop.
I plan to refactor to std cpp/eigen/bh in the current project.
In current branch, the cmd will be

```sh
build/grab \
  --method POLMM \
  --pheno examples/simuPHENO.txt \
  --pheno-ordinal OrdinalPheno \
  --covar-name AGE,GENDER,PC1,PC2 \
  --bfile examples/simuPLINK \
  --sp-grm-grab examples/SparseGRM.txt \
  --out tmp/POLMM_output.txt \
  --threads 3
```

Please read the R code and get the full picture.
Only consider marker-level analysis using sparse grm; no need loco; check if the default options are optimized.

Don't just translate to cpp. Learn POLMM logic and redesign the framework to be compatible with current cpp framework. Remove useless steps and optimize for speed if possible.
