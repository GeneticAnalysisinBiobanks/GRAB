set -e

## SPACox

build/grab --threads 2 \
  --method SPACox \
  --seed 2026 \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --compression zst \
  --out examples_output/1kg

build/grab --threads 2 \
  --method SPACox \
  --pheno examples/1kg.pheno \
  --resid-name Resid1,Resid2 \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --compression zst \
  --out examples_output/1kg

## SPAmix

build/grab \
  --cal-af-coef \
  --pheno examples/1kg.pheno \
  --pc-cols PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --compression zst \
  --out examples_output/1kg

build/grab --threads 2 \
  --method SPAmix \
  --seed 2026 \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --ind-af-coef examples_output/1kg.afc.zst \
  --compression zst \
  --out examples_output/1kg

build/grab --threads 2 \
  --method SPAmix \
  --pheno examples/1kg.pheno \
  --resid-name Resid1,Resid2 \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --ind-af-coef examples_output/1kg.afc.zst \
  --compression zst \
  --out examples_output/1kg

## SPAGRM

build/grab \
  --cal-pairwise-ibd \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pfile examples/1kg \
  --compression zst \
  --out examples_output/1kg

build/grab --threads 2 \
  --method SPAGRM \
  --seed 2026 \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd examples_output/1kg.ibd.zst \
  --pfile examples/1kg \
  --compression zst \
  --out examples_output/1kg

build/grab --threads 2 \
  --method SPAGRM \
  --pheno examples/1kg.pheno \
  --resid-name Resid1,Resid2 \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd examples_output/1kg.ibd.zst \
  --pfile examples/1kg \
  --compression zst \
  --out examples_output/1kg

## SAGELD

build/grab --threads 2 \
  --method SAGELD \
  --pheno examples/long_pheno \
  --covar-name MALE,TIME,PC1,PC2 \
  --pheno-name Long1,Long2 \
  --sageld-x TIME \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd examples_output/1kg.ibd.zst \
  --pfile examples/1kg \
  --compression zst \
  --out examples_output/1kg

build/grab --threads 2 \
  --method SAGELD \
  --pheno examples/long_pheno_resid \
  --resid-name R_G,R_TIME,R_GxTIME \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd examples_output/1kg.ibd.zst \
  --pfile examples/1kg \
  --compression zst \
  --out examples_output/1kg

## SPAsqr

build/grab --threads 2 \
  --method SPAsqr \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pred-list examples/loco_prs.list \
  --pfile examples/1kg \
  --compression zst \
  --out examples_output/1kg

## WtCoxG

build/grab --threads 2 \
  --method WtCoxG \
  --pheno examples/1kg.pheno \
  --pheno-name Time:Event,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --ref-af examples/ref_pop1.afreq \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --prevalence 0.1 \
  --compression zst \
  --out examples_output/1kg

## LEAF

build/grab --threads 2 \
  --method LEAF \
  --seed 2026 \
  --leaf-nclusters 3 \
  --remove examples/leaf_kmeans_remove \
  --pheno examples/1kg.pheno \
  --pheno-name Time:Event,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --ref-af examples/ref_pop1.afreq,examples/ref_pop2.afreq \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --prevalence 0.1 \
  --compression zst \
  --out examples_output/1kg
