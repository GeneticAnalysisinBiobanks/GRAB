## SPACox

build/grab \
  --method SPACox \
  --null-resid examples/simuResid_2resids.txt \
  --covar examples/simuDesign.txt \
  --bfile examples/simuPLINK \
  --out tmp/SPACox_output.txt

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

## SPAmixPlus

build/grab \
  --cal-ind-af-coef \
  --covar examples/simuPheno.txt \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --out tmp/SPAmixAfCoef.txt.gz

build/grab \
  --method SPAmixPlus \
  --null-resid examples/simuResid_2cols.txt \
  --covar examples/simuPheno.txt \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --ind-af-coef tmp/SPAmixAfCoef.txt.gz \
  --sp-grm-grab examples/SparseGRM.txt \
  --out tmp/SPAmixPlus_output.txt
