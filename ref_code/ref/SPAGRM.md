# Refactor SPAGRM from R/Rcpp to pure C++

## Old workflow: R and Rcpp

```R
source('ref_code/src/SPAGRM.R')
source('ref_code/src/mtMarker.R')

getPairwiseIBD(
  PlinkPrefix = "examples/simuPLINK",
  SparseGRMFile = "examples/SparseGRM.txt",
  PairwiseIBDFile = "examples/PairwiseIBD.txt"
)

obj.SPAGRM <- SPAGRM.NullModel(
  ResidMatFile = "examples/simuResid_2cols.txt",
  SparseGRMFile = "examples/SparseGRM.txt",
  PairwiseIBDFile = "examples/PairwiseIBD.txt"
)

GRAB.mtMarker(
  obj.SPAGRM,
  "examples/simuPLINK.bed",
  "tmp/SPAGRM_output.txt"
)
```

## New wrokflow: pure cpp program do followings

```sh
build/grab \
  --cal-pairwise-ibd \
  --sparse-grm examples/SparseGRM.txt \
  --bfile examples/simuPLINK \
  --out examples/PairwiseIBD.txt

build/grab \
  --method SPAGRM \
  --null-resid examples/simuResid_2cols.txt \
  --sparse-grm examples/SparseGRM.txt \
  --pairwise-ibd examples/PairwiseIBD.txt \
  --bfile examples/simuPLINK \
  --out tmp/SPAGRM_output.txt
```
