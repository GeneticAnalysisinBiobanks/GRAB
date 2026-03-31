# Examples

## WtCoxG

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

## LEAF

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

## SPAsqr

```sh
build/grab \
  --method SPAsqr \
  --null-resid examples/simuResidMat.txt \
  --sparse-grm examples/SparseGRM.txt \
  --bfile examples/simuPLINK \
  --out tmp/SPAsqr_output.txt
```

## SPACox

```sh
build/grab \
  --method SPACox \
  --null-resid examples/simuResid_2cols.txt \
  --design examples/simuDesign.txt \
  --bfile examples/simuPLINK \
  --out tmp/SPACox_output.txt
```

## individual AF coefficient

```sh
build/grab \
  --cal-af-coef \
  --eigenvec examples/simuPCs.txt \
  --bfile examples/simuPLINK \
  --out tmp/SPAmixAfCoef.txt.gz
```

## SPAmix

```sh
build/grab \
  --method SPAmix \
  --null-resid examples/simuResid_2cols.txt \
  --eigenvec examples/simuPCs.txt \
  --bfile examples/simuPLINK \
  --out tmp/SPAmix_output.txt
```

## SPAmixPlus

```sh
build/grab \
  --method SPAmixPlus \
  --null-resid examples/simuResid_2cols.txt \
  --eigenvec examples/simuPCs.txt \
  --bfile examples/simuPLINK \
  --af-coef tmp/SPAmixAfCoef.txt.gz \
  --sparse-grm examples/SparseGRM.txt \
  --out tmp/SPAmixPlus_output.txt
```

## plink pairwise IBD and sparse GRM

```sh
plink --genome
plink2 --make-king 
plink2 --make-grm-sparse
```
