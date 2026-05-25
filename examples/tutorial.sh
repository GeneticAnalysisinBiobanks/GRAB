#!/usr/bin/env bash
# examples/tutorial.sh — minimal walkthrough of each analysis method.
# One simple command per method, listing only the mandatory options.

set -e

OUT_DIR=examples_output/tutorial
mkdir -p ${OUT_DIR}

## ── Smoke test: report the version ────────────────────────────────────
build/grab2 --version

## ── Utility: cal-pairwise-ibd (input to SPAGRM / SAGELD) ──────────────
build/grab2 \
  --cal-pairwise-ibd \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pfile examples/1kg \
  --out ${OUT_DIR}/ibd

## ── SPACox ────────────────────────────────────────────────────────────
build/grab2 \
  --method SPACox \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spacox

## ── SPAGRM ────────────────────────────────────────────────────────────
build/grab2 \
  --method SPAGRM \
  --pheno examples/1kg.pheno \
  --pheno-name Binary,Time:Event \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT_DIR}/ibd.ibd \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spagrm

## ── SPAmix ────────────────────────────────────────────────────────────
# Ordinal regression is randomly seeded.
build/grab2 \
  --seed 2026 \
  --method SPAmix \
  --pheno examples/1kg.pheno \
  --pheno-name Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spamix

## ── SAGELD: longitudinal score test (GxE on TIME) ─────────────────────
build/grab2 \
  --method SAGELD \
  --pheno examples/long_pheno \
  --pheno-name Long1,Long2 \
  --covar-name MALE,TIME,PC1,PC2 \
  --sageld-x TIME \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT_DIR}/ibd.ibd \
  --pfile examples/1kg \
  --out ${OUT_DIR}/sageld

## ── Utility: int-pheno (inverse-normal transform; input to SPAsqr) ────
build/grab2 \
  --int-pheno \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time \
  --out ${OUT_DIR}/int_pheno

## ── SPAsqr (score mode): smoothed quantile-regression score test ──────
build/grab2 \
  --method SPAsqr \
  --pheno ${OUT_DIR}/int_pheno.txt \
  --pheno-name Quantitative,Time \
  --covar examples/1kg.pheno \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pred-list examples/loco_prs.list \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spasqr

## ── SPAsqr (wald mode): per-marker effect size + SE ───────────────────
build/grab2 \
  --threads 4 \
  --method SPAsqr \
  --spasqr-mode wald \
  --pheno ${OUT_DIR}/int_pheno.txt \
  --pheno-name Quantitative,Time \
  --covar examples/1kg.pheno \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pred-list examples/loco_prs.list \
  --extract examples/spasqr_wald_extract \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spasqr_wald

## ── WtCoxG: weighted Cox / logistic with external AF reference ────────
build/grab2 \
  --method WtCoxG \
  --pheno examples/1kg.pheno \
  --pheno-name Time:Event,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --ref-af examples/ref_pop1.afreq \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --prevalence 0.1 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/wtcoxg

## ── LEAF: Multi-ancestry extension of WtCoxG ──────────────────────────
# Kmeans clustering is randomly seeded.
build/grab2 \
  --seed 2026 \
  --method LEAF \
  --pheno examples/1kg.pheno \
  --pheno-name Time:Event,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --ref-af examples/ref_pop1.afreq,examples/ref_pop2.afreq \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --prevalence 0.1 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/leaf

echo
echo "Tutorial finished.  Output written to ${OUT_DIR}/."
