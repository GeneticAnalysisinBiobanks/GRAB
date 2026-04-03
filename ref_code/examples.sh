## WtCoxG

build/grab \
  --method WtCoxG \
  --null-resid examples/simuResid.txt \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --prevalence 0.1 \
  --out tmp/WtCoxG_output.txt \
  --threads 3


## LEAF

build/grab \
  --method LEAF \
  --null-resid examples/simuResid1.txt,examples/simuResid2.txt,examples/simuResid3.txt \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq,examples/simuPLINK.afreq \
  --sp-grm-grab examples/SparseGRM.txt \
  --prevalence 0.1 \
  --out tmp/LEAF_output.txt \
  --chunk-size 256 \
  --threads 4


## SPAGRM

build/grab \
  --cal-pairwise-ibd \
  --sp-grm-grab examples/SparseGRM.txt \
  --bfile examples/simuPLINK \
  --out examples/PairwiseIBD.txt

build/grab \
  --method SPAGRM \
  --null-resid examples/simuResid_2cols.txt \
  --sp-grm-grab examples/SparseGRM.txt \
  --pairwise-ibd examples/PairwiseIBD.txt \
  --bfile examples/simuPLINK \
  --out tmp/SPAGRM_output.txt


## SPAsqr

build/grab \
  --method SPAsqr \
  --null-resid examples/simuResidMat.txt \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --bfile examples/simuPLINK \
  --out tmp/SPAsqr_output.txt


## SPACox

build/grab \
  --method SPACox \
  --null-resid examples/simuResid_2resids.txt \
  --covar examples/simuDesign.txt \
  --bfile examples/simuPLINK \
  --out tmp/SPACox_output.txt


## individual AF coefficient

build/grab \
  --cal-ind-af-coef \
  --covar examples/simuPheno.txt \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --out tmp/SPAmixAfCoef.txt.gz


## SPAmix

build/grab \
  --method SPAmix \
  --null-resid examples/simuResid_1col.txt \
  --covar examples/simuPheno.txt \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --out tmp/SPAmix_output.txt


## SPAmixPlus

build/grab \
  --method SPAmixPlus \
  --null-resid examples/simuResid_2cols.txt \
  --covar examples/simuPheno.txt \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --ind-af-coef tmp/SPAmixAfCoef.txt.gz \
  --sp-grm-grab examples/SparseGRM.txt \
  --out tmp/SPAmixPlus_output.txt

## SPAsqr

```sh
build/grab \
  --method SPAsqr \
  --pheno examples/simuPHENO.txt \
  --pheno-quantitative QuantPheno \
  --spasqr-taus 0.1,0.3,0.5,0.7,0.9 \
  --covar-name AGE,GENDER,PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --out tmp/SPAsqr_output.txt \
  --threads 3
```