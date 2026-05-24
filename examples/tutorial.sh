#!/usr/bin/env bash
# examples/tutorial.sh — minimal walkthrough of each GRAB analysis method.
#
# One simple command per method, exercising two phenotypes.  Intended
# for new users who want to see each method run end-to-end with the
# fewest possible flags; every command lists only the mandatory inputs
# plus a phenotype list and an output prefix.  All other knobs fall
# back to grab2's built-in defaults.
#
# For exhaustive flag coverage and cross-format / cross-codec regression
# checks, see examples/baseline.sh.

set -e

OUT_DIR=examples_output/tutorial
mkdir -p ${OUT_DIR}

## ── Utility: cal-pairwise-ibd (input to SPAGRM / SAGELD) ──────────────
# Pairwise IBD probabilities derived from a sparse GRM; SPAGRM and
# SAGELD condition on these to retain power for related subjects.

build/grab2 \
  --cal-pairwise-ibd \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pfile examples/1kg \
  --out ${OUT_DIR}/ibd

## ── SPACox: saddlepoint-corrected score test ─────────────────────────
# Two phenotypes covering quantitative (linear regression) and survival
# (Cox proportional hazards).  --pheno-name selects columns from the
# phenotype file; "Time:Event" denotes a Cox spec.

build/grab2 \
  --method SPACox \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spacox

## ── SPAGRM: GRM-conditioned score test for related subjects ──────────
# Consumes the sparse GRM and the pairwise-IBD file produced above.
# Binary + Ordinal exercise the discrete-trait code paths; --seed
# fixes the SPA Monte-Carlo so the output is bit-reproducible.

build/grab2 \
  --method SPAGRM \
  --pheno examples/1kg.pheno \
  --pheno-name Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT_DIR}/ibd.ibd \
  --seed 2026 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spagrm

## ── SPAmix: score test for admixed populations ───────────────────────
# Without --ind-af-coef SPAmix fits the per-subject allele-frequency
# model on the fly from --pc-cols.  --seed fixes the SPA Monte-Carlo.

build/grab2 \
  --method SPAmix \
  --pheno examples/1kg.pheno \
  --pheno-name Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --seed 2026 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spamix

## ── SAGELD: longitudinal score test (GxE on TIME) ────────────────────
# Uses the longitudinal phenotype fixture (one row per subject × visit).

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

## ── Utility: int-pheno (inverse-normal transform; input to SPAsqr) ───
# Emits a phenotype file with INT-transformed Quantitative and Time
# columns.  SPAsqr below consumes this file via --pheno and pulls the
# original covariates via --covar.

build/grab2 \
  --int-pheno \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time \
  --out ${OUT_DIR}/int_pheno

## ── SPAsqr (score mode): smoothed quantile-regression score test ─────
# Consumes the INT-transformed phenotype file produced above; the
# --pheno-transform default is still `int` and is idempotent on already-
# INT'd inputs.

build/grab2 \
  --method SPAsqr \
  --pheno ${OUT_DIR}/int_pheno.txt \
  --pheno-name Quantitative,Time \
  --covar examples/1kg.pheno \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spasqr

## ── SPAsqr (wald mode): per-marker effect size + SE ──────────────────
# Wald refits the joint smoothed-QR model with [X | G] for every
# (marker, τ) and is appreciably slower than score mode; restrict to
# the 10 variant IDs in examples/spasqr_wald_extract for a quick demo.

build/grab2 \
  --method SPAsqr \
  --spasqr-mode wald \
  --pheno ${OUT_DIR}/int_pheno.txt \
  --pheno-name Quantitative,Time \
  --covar examples/1kg.pheno \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --extract examples/spasqr_wald_extract \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spasqr_wald

## ── WtCoxG: weighted Cox / logistic with external AF reference ───────
# Survival phenotype Time:Event paired with binary Binary; both consume
# the same ref-AF file and disease prevalence.

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

## ── LEAF: latent-cluster admixed Cox / logistic ──────────────────────
# Multi-ancestry extension of WtCoxG; consumes one --ref-af file per
# reference population.

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
  --out ${OUT_DIR}/leaf

echo
echo "Tutorial finished.  Output written to ${OUT_DIR}/."
