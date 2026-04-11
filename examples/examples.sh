## POLMM

build/grab \
  --method POLMM \
  --pheno examples/simuPHENO.txt \
  --pheno-ordinal OrdinalPheno \
  --covar-name AGE,GENDER,PC1,PC2 \
  --bfile examples/simuPLINK \
  --sp-grm-grab examples/SparseGRM.txt \
  --out tmp/POLMM_output \
  --threads 3

## WtCoxG

build/grab \
  --method WtCoxG \
  --null-resid examples/simuResid.txt \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --prevalence 0.1 \
  --out tmp/WtCoxG_output \
  --threads 3

build/grab \
  --method WtCoxG \
  --pheno examples/simuPHENO.txt \
  --pheno-binary BinaryPheno \
  --covar-name AGE,GENDER,PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --prevalence 0.1 \
  --out tmp/WtCoxG_binary \
  --threads 3

build/grab \
  --method WtCoxG \
  --pheno examples/simuPHENO.txt \
  --pheno-surv SurvTime:SurvEvent \
  --covar-name AGE,GENDER,PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --prevalence 0.1 \
  --out tmp/WtCoxG_survival \
  --threads 3

## LEAF

build/grab \
  --method LEAF \
  --null-resid examples/simuResid1.txt,examples/simuResid2.txt,examples/simuResid3.txt \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq,examples/simuPLINK.afreq \
  --sp-grm-grab examples/SparseGRM.txt \
  --prevalence 0.1 \
  --out tmp/LEAF_output \
  --chunk-size 256 \
  --threads 3

build/grab \
  --method LEAF \
  --leaf-nclusters 3 \
  --pheno examples/simuPHENO.txt \
  --pheno-binary BinaryPheno \
  --covar-name AGE,GENDER,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq,examples/simuPLINK.afreq \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --prevalence 0.1 \
  --out tmp/LEAF_binary \
  --threads 3

build/grab \
  --method LEAF \
  --leaf-nclusters 3 \
  --pheno examples/simuPHENO.txt \
  --pheno-surv SurvTime:SurvEvent \
  --covar-name AGE,GENDER,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --ref-af examples/simuPLINK.afreq,examples/simuPLINK.afreq \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --prevalence 0.1 \
  --out tmp/LEAF_survival \
  --threads 3

## SPAGRM

build/grab \
  --cal-pairwise-ibd \
  --sp-grm-grab examples/SparseGRM.txt \
  --bfile examples/simuPLINK \
  --out examples/PairwiseIBD

build/grab \
  --method SPAGRM \
  --null-resid examples/simuResid_2cols.txt \
  --sp-grm-grab examples/SparseGRM.txt \
  --pairwise-ibd examples/PairwiseIBD.ibd \
  --bfile examples/simuPLINK \
  --out tmp/SPAGRM_output

build/grab \
  --method SPAGRM \
  --null-resid examples/simuResid_2resids.txt \
  --sp-grm-grab examples/SparseGRM.txt \
  --pairwise-ibd examples/PairwiseIBD.ibd \
  --bfile examples/simuPLINK \
  --out tmp/SPAGRM_output

## SPAsqr

build/grab \
  --method SPAsqr \
  --null-resid examples/simuResidMat.txt \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --bfile examples/simuPLINK \
  --out tmp/SPAsqr_output

build/grab \
  --method SPAsqr \
  --pheno examples/simuPHENO.txt \
  --pheno-quant QuantPheno \
  --spasqr-taus 0.1,0.3,0.5,0.7,0.9 \
  --covar-name AGE,GENDER,PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --out tmp/SPAsqr_output \
  --threads 3

## SPACox

build/grab \
  --method SPACox \
  --null-resid examples/simuResid_2resids.txt \
  --covar examples/simuDesign.txt \
  --bfile examples/simuPLINK \
  --out tmp/SPACox_output \
  --compression gz \
  --compression-level 6

## individual AF coefficient

build/grab \
  --cal-af-coef \
  --covar examples/simuPHENO.txt \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --out tmp/SPAmixAfCoef

## SPAmix

build/grab \
  --method SPAmix \
  --null-resid examples/simuResid_1col.txt \
  --covar examples/simuPHENO.txt \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --out tmp/SPAmix_output

## SPAmixPlus

build/grab \
  --method SPAmixPlus \
  --null-resid examples/simuResid_2cols.txt \
  --covar examples/simuPHENO.txt \
  --pc-cols PC1,PC2,PC3,PC4 \
  --bfile examples/simuPLINK \
  --ind-af-coef tmp/SPAmixAfCoef.afc \
  --sp-grm-grab examples/SparseGRM.txt \
  --out tmp/SPAmixPlus_output

## SPAmixLocalPlus

build/grab --make-abed \
  --vcf examples/spamixlocalp/admix_bfile/head10.bcf \
  --rfmix-msp examples/spamixlocalp/admix_bfile/head10.msp.gz \
  --out examples/spamixlocalp/admix_bfile/head10 
  
build/grab --make-abed \
  --admix-text-prefix examples/spamixlocalp/simuAncestry \
  --out examples/spamixlocalp/simuAncestry

build/grab \
  --cal-phi \
  --admix-bfile examples/spamixlocalp/simuAncestry \
  --sp-grm-plink2 examples/simuPLINK.grm.sp \
  --out examples/spamixlocalp/simuAncestry

build/grab \
  --method SPAmixLocalPlus \
  --null-resid examples/simuResid_2cols.txt \
  --admix-bfile examples/spamixlocalp/simuAncestry \
  --admix-phi examples/spamixlocalp/simuAncestry.phi \
  --out tmp/SPAmixLocalPlus_output

# SAGELD

build/grab \
  --method SAGELD \
  --null-resid examples/simuResid_SAGELD2.txt \
  --sp-grm-grab examples/SparseGRM.txt \
  --pairwise-ibd examples/PairwiseIBD.ibd \
  --bfile examples/simuPLINK \
  --out tmp/SAGELD_output2
  