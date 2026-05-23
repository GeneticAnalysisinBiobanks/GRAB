#!/usr/bin/env bash
# examples/run.sh — canonical regression script for the GRAB repository.
#
# This is the single end-to-end script that exercises every utility mode
# (cal-af-coef, cal-pairwise-ibd, int-pheno) and every analysis method
# (SPACox, SPAmix, SPAGRM, SAGELD, SPAsqr, WtCoxG, LEAF) against the
# bundled examples/1kg.* fixtures.  It serves two purposes:
#
#   1. Documentary.  Every command spells out every command-line flag the
#      method accepts.  Mandatory flags appear first; optional flags
#      follow a "# Optional flags below:" comment with each numeric or
#      categorical knob set to its built-in default value.  A first-time
#      user can copy the lines above the comment to get a working
#      invocation, and refer to the lines below for every available
#      tuning knob.
#
#   2. Regression baseline.  After any refactor (shared engine, SIMD
#      kernels, null-model fitting, genotype readers, output formatting,
#      per-method code), re-run this script and confirm that the
#      resulting examples_output/* artifacts are byte-identical (or
#      numerically identical up to documented tolerance) to the
#      pre-refactor baseline.  A passing build is not sufficient
#      evidence that a refactor preserved behavior; output equivalence
#      is.

set -e

OUT_DIR=examples_output
mkdir -p ${OUT_DIR}
OUT=${OUT_DIR}/fit          # output prefix for fit-mode runs
RESID_OUT=${OUT_DIR}/resid  # output prefix for residual-mode runs (SPACox, SAGELD)

## ── Utility: cal-af-coef ──────────────────────────────────────────────
# Produces ${OUT}.afc.zst, consumed by SPAmix via --ind-af-coef.

build/grab2 \
  --cal-af-coef \
  --pheno examples/1kg.pheno \
  --pc-cols PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --compression zst \
  --compression-level 3

## ── Utility: cal-pairwise-ibd ─────────────────────────────────────────
# Produces ${OUT}.ibd.zst, consumed by SPAGRM / SAGELD via --pairwise-ibd.

build/grab2 \
  --cal-pairwise-ibd \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --min-maf-ibd 0.01 \
  --threads 2 \
  --compression zst \
  --compression-level 3

## ── SPACox (fit mode, --pheno-name) ───────────────────────────────────

build/grab2 \
  --method SPACox \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --save-resid \
  --chr 1-2,3 \
  --covar-p-threshold 5e-5 \
  --spa-z-threshold 2.0 \
  --seed 2026 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0

## ── SPACox (residual mode, --resid-name) ──────────────────────────────
# Consumes the combined residual file produced by the SPACox fit-mode

build/grab2 \
  --method SPACox \
  --pheno ${OUT}.null.resid \
  --resid-name Quantitative,Time_Event,Binary,Ordinal \
  --covar examples/1kg.pheno \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --out ${RESID_OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --covar-p-threshold 5e-5 \
  --spa-z-threshold 2.0 \
  --seed 2026 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0

## ── SPAmix (fit mode, --pheno-name) ───────────────────────────────────

build/grab2 \
  --method SPAmix \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --ind-af-coef ${OUT}.afc.zst \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --chr 1-2,3 \
  --outlier-iqr-multiplier 1.5 \
  --spa-z-threshold 2.0 \
  --seed 2026 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --compression zst \
  --compression-level 3

## ── SPAGRM (fit mode, --pheno-name) ───────────────────────────────────

build/grab2 \
  --method SPAGRM \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT}.ibd.zst \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --seed 2026 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --compression zst \
  --compression-level 3

## ── SAGELD (fit mode: --pheno-name + --sageld-x) ────────────────────

build/grab2 \
  --method SAGELD \
  --pheno examples/long_pheno \
  --pheno-name Long1,Long2 \
  --covar-name MALE,TIME,PC1,PC2 \
  --sageld-x TIME \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT}.ibd.zst \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --save-resid \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --compression zst \
  --compression-level 3

## ── SAGELD (residual mode, --resid-name) ──────────────────────────────
# Consumes the residual file produced by the SAGELD fit mode

build/grab2 \
  --method SAGELD \
  --pheno ${OUT}.Long1.SAGELD.resid \
  --resid-name R_G,R_TIME,R_GxTIME \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT}.ibd.zst \
  --pfile examples/1kg \
  --out ${RESID_OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --compression zst \
  --compression-level 3

## ── Utility: int-pheno ────────────────────────────────────────────────
# Produces ${OUT}.int.txt, a phenotype file containing the
# INT-transformed Quantitative and Time columns only.

build/grab2 \
  --int-pheno \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time \
  --out ${OUT_DIR}/int_pheno

## ── SPAsqr (score mode, fit path) ─────────────────────────────────────
# Consumes the INT-transformed phenotype file produced above.

build/grab2 \
  --method SPAsqr \
  --pheno ${OUT_DIR}/int_pheno.txt \
  --pheno-name Quantitative,Time \
  --covar examples/1kg.pheno \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pred-list examples/loco_prs.list \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --spasqr-taus 0.1,0.3,0.5,0.7,0.9 \
  --spasqr-tol 1e-7 \
  --spasqr-h-scale 3 \
  --spasqr-solver qmme \
  --spasqr-mode score \
  --pheno-transform int \
  --outlier-iqr-multiplier 1.5 \
  --spasqr-outlier-abs-bound 0.55 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --compression zst \
  --compression-level 3

## ── WtCoxG ────────────────────────────────────────────────────────────

build/grab2 \
  --method WtCoxG \
  --pheno examples/1kg.pheno \
  --pheno-name Time:Event,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --ref-af examples/ref_pop1.afreq \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --prevalence 0.1 \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --chr 1-2,3 \
  --batch-effect-p-threshold 0.05 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --compression gz \
  --compression-level 6

## ── LEAF ──────────────────────────────────────────────────────────────

build/grab2 \
  --method LEAF \
  --pheno examples/1kg.pheno \
  --pheno-name Time:Event,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --ref-af examples/ref_pop1.afreq,examples/ref_pop2.afreq \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --prevalence 0.1 \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --remove examples/leaf_kmeans_remove \
  --leaf-nclusters 3 \
  --leaf-kmeans-nstart 25 \
  --seed 2026 \
  --chr 1-2,3 \
  --batch-effect-p-threshold 0.05 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-size 8192 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --compression gz \
  --compression-level 6
