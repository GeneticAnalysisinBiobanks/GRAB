#!/usr/bin/env bash
# examples/baseline.sh — exhaustive regression script for the GRAB repository.
#
# This script exercises every utility mode (cal-af-coef, cal-pairwise-ibd,
# int-pheno) and every analysis method (SPACox, SPAmix, SPAGRM, SAGELD,
# SPAsqr score + wald, WtCoxG, LEAF, SPAGxE base + sparse-GRM, SPAGxEmix)
# against the bundled examples/1kg.* fixtures.  It serves two purposes:
#
#   1. Documentary.  Every command spells out every command-line flag the
#      method accepts.  Mandatory flags appear first; optional flags
#      follow a "# Optional flags below:" comment with each numeric or
#      categorical knob set to its built-in default value.  Every CLI
#      flag should appear somewhere in this script; when a new flag is
#      added, the corresponding block must gain a line listing it at
#      its default value.
#
#   2. Regression baseline.  After any refactor (shared engine, SIMD
#      kernels, null-model fitting, genotype readers, output formatting,
#      per-method code), re-run this script and confirm that the
#      resulting examples_output/* artifacts are byte-identical (or
#      numerically identical up to documented tolerance) to the
#      pre-refactor baseline.  A passing build is not sufficient
#      evidence that a refactor preserved behavior; output equivalence
#      is.
#
# For a minimal, user-facing walkthrough of each method with only the
# mandatory flags (no defaults spelled out, no regression cross-checks),
# see examples/tutorial.sh.

set -e

OUT_DIR=examples_output/baseline
mkdir -p ${OUT_DIR}
OUT=${OUT_DIR}/fit          # output prefix for fit-mode runs
RESID_OUT=${OUT_DIR}/resid  # output prefix for residual-mode runs (SPACox, SAGELD)
LOCO=${OUT_DIR}/loco        # output prefix for LOCO runs (SPACox/SPAmix/SPAGRM/WtCoxG)

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
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
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
#
# Output columns (spa_unify Stage 3): the nine meta columns, then
#   P  LOG10P  Z  Z_Norm  BETA  SE  SPA_STATUS
# LOG10P is -log10(P) computed in the log domain, so it stays meaningful
# past the point where the linear-scale normal tail underflows.  SPA_STATUS
# is the spa::Status enumerator of the saddlepoint that produced P, as an
# integer: 0 OK, 1 MAXITER, 2 GUARD_TEMP, 3 GUARD_CURV, 4 GUARD_W,
# 5 NONFINITE, 6 NORMAL (|Z| <= --spa-z-threshold, saddlepoint not
# attempted).  P and LOG10P are NA for every status other than 0 and 6, so
# a saddlepoint failure is reported rather than silently substituted.

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
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1

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
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1

## ── SPAmix (fit mode, --pheno-name) ───────────────────────────────────
#
# Output columns (spa_unify Stage 5): the nine meta columns, then
#   P  LOG10P  Z  Z_Norm  BETA  SE  SPA_STATUS
# with the same meaning and the same integer SPA_STATUS encoding as the
# SPACox and SPAGRM blocks above (0 OK, 1 MAXITER, 2 GUARD_TEMP,
# 3 GUARD_CURV, 4 GUARD_W, 5 NONFINITE, 6 NORMAL).  Stage 5 also stopped
# SPAmix returning P = 1 for a marker with a non-positive score variance:
# three markers of this fixture, all with ALT_FREQ > 0.99, have every
# per-individual q_i saturated at 1 by the AF model and therefore have no
# statistic at all.  They now report NA with SPA_STATUS 5 rather than a
# fabricated perfectly-null p-value alongside NA in Z, BETA and SE.

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
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── SPAmix on-the-fly AF (omit --ind-af-coef) ─────────────────────────
# Same SPAmix invocation as above, but with --ind-af-coef removed: the
# per-individual AF model is refit per marker from --pc-cols at run time
# instead of being read from the pre-computed ${OUT}.afc.zst file.  In a
# dataset with no missing genotypes the two paths fit the same logistic
# AF model on the same N subjects and therefore must produce identical
# AFVec, identical scores, and byte-identical output tables.  The
# cross-check below asserts md5 equivalence per phenotype.

build/grab2 \
  --method SPAmix \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary,Ordinal \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spamix_otf \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --chr 1-2,3 \
  --outlier-iqr-multiplier 1.5 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── SPAGRM (fit mode, --pheno-name) ───────────────────────────────────
#
# Output columns (spa_unify Stage 4): the nine meta columns, then
#   P  LOG10P  Z  Z_Norm  BETA  SE  SPA_STATUS
# with the same meaning and the same integer SPA_STATUS encoding as the
# SPACox block above (0 OK, 1 MAXITER, 2 GUARD_TEMP, 3 GUARD_CURV,
# 4 GUARD_W, 5 NONFINITE, 6 NORMAL).  Stage 4 also gave SPAGRM the shared
# guard set, so a saddlepoint that used to emit a bare NaN now reports NA in
# P with the reason named.
#
# GUARD_TEMP no longer occurs anywhere in this script, and that is a property
# worth watching rather than a coincidence.  zeta*Score - K(zeta) is the
# Legendre transform of a cumulant generating function and is non-negative by
# construction, so GUARD_TEMP cannot arise from a CGF that is one — it is a
# statement about the CGF's inputs, not about the solver.  Fourteen markers of
# the cross-format Binary fixture below used to report it, because the
# Chow-Liu family tables were built from IBD triples that did not sum to one
# and the class-3 block therefore carried K(0) != 0.  All fourteen now
# converge.  If GUARD_TEMP reappears, look at the pairwise-IBD input first.

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
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── SAGELD (fit mode: --pheno-name + --envir-name) ────────────────────
# Output columns (log10p_unify Stage 7): the nine meta columns, a six-wide
# G main-effect block  P_G LOG10P_G Z_G BETA_G SE_G SPA_STATUS_G, then a
# seven-wide block per environment
#   P_Gx<E>  LOG10P_Gx<E>  Z_Gx<E>  Z_Norm_Gx<E>  BETA_Gx<E>  SE_Gx<E>
#   SPA_STATUS_Gx<E>
# LOG10P and SPA_STATUS are new in Stage 7.  Both were already produced by
# SPAGRMClass::getMarkerPvalFromScore -- the routine SAGELD shares with
# --method spagrm -- and discarded at the call site, which CLAUDE.md listed as
# the first of the known gaps in the output-column contract.  SPA_STATUS_Gx is
# therefore the same saddlepoint outcome SPAGRM reports on the same cohort;
# SPA_STATUS_G is the constant 1 (NORMAL), the G main effect being a plain
# two-sided normal test that never attempts a saddlepoint.

build/grab2 \
  --method SAGELD \
  --pheno examples/long_pheno \
  --pheno-name Long1,Long2 \
  --covar-name MALE,TIME,PC1,PC2 \
  --envir-name TIME \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT}.ibd.zst \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --save-resid \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── SAGELD (residual mode, --resid-name; Long1 and Long2) ─────────────
# Consumes the per-phenotype residual files produced by the SAGELD fit
# mode (one --resid output per phenotype, with R_G/R_TIME/R_GxTIME).
# Output prefix is suffixed by the phenotype name so that fit and resid
# artefacts for each of Long1 and Long2 can be compared independently
# in the regression cross-checks block at the bottom of this script.

for sageld_pheno in Long1 Long2; do
build/grab2 \
  --method SAGELD \
  `# Fit mode wrote this with --compression zst, so the residual carries a` \
  `# .zst suffix; residual mode's TextReader auto-detects and decompresses it.` \
  --pheno ${OUT}.${sageld_pheno}.SAGELD.resid.zst \
  --resid-name R_G,R_TIME,R_GxTIME \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT}.ibd.zst \
  --pfile examples/1kg \
  --out ${RESID_OUT}_${sageld_pheno} \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3
done

## ── SAGELD GALLOP variant (--sageld-method gallop) ───────────────────
# The exact Sikorska et al. (2013) two-step Wald test (develop-R
# SAGELD.NullModel UsedMethod="GALLOP").  Same random-slope null model
# Y ~ X + (TIME | IID) as the SAGELD fit block above, but the per-marker
# test re-fits the genetic main effect (G) and the G x TIME interaction
# exactly, emitting BETA / SE / Pvalue per effect.  GALLOP captures
# relatedness through the random intercept/slope, so it takes NEITHER
# --sp-grm NOR --pairwise-ibd; supplying them only triggers a warning.
# Output: ${OUT}.<pheno>.GALLOP.gz (gzip path through multiPhenoEngine).
# Columns are two six-wide Wald blocks, P LOG10P Z BETA SE SPA_STATUS for the
# G main effect and for the G x TIME interaction.  SPA_STATUS is 1 (NORMAL) on
# both wherever the 2x2 solve produced a positive SE -- GALLOP runs no
# saddlepoint at all -- and 8 (NA_NO_TEST) where it did not.

build/grab2 \
  --method SAGELD \
  --sageld-method gallop \
  --pheno examples/long_pheno \
  --pheno-name Long1,Long2 \
  --covar-name MALE,TIME,PC1,PC2 \
  --envir-name TIME \
  --pfile examples/1kg \
  --out ${OUT} \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression gz \
  --compression-level 3

## ── Longitudinal input (--longitudinal): SPACox / SPAmix / SPAGRM ─────
# These three single-residual methods accept a long-format (repeated-
# measures) phenotype directly.  For each --pheno-name outcome a random-
# intercept null model  Y ~ X + (1 | IID)  is fit (X = intercept +
# --covar-name covariates), and the per-IID aggregated residual
# R_G = sum_j r_ij is fed to the marker test (marginal MAIN genetic
# effect; no --envir-name and no G x E).  Reuses the SAGELD long-format
# fixture examples/long_pheno (#IID MALE PC1 PC2 TIME Long1 Long2).

## ── SPACox longitudinal (plain-text output) ───────────────────────────

build/grab2 \
  --method SPACox \
  --pheno examples/long_pheno \
  --pheno-name Long1,Long2 \
  --longitudinal \
  --covar-name MALE,PC1,PC2 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/long_spacox \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --covar-p-threshold 5e-5 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1

## ── SPAmix longitudinal (zstd; on-the-fly AF) ─────────────────────────
# --pc-cols must be a subset of --covar-name in longitudinal mode, because
# the per-individual AF model sources its PCs from the long-format design.

build/grab2 \
  --method SPAmix \
  --pheno examples/long_pheno \
  --pheno-name Long1,Long2 \
  --longitudinal \
  --covar-name MALE,PC1,PC2 \
  --pc-cols PC1,PC2 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/long_spamix \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --outlier-iqr-multiplier 1.5 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── SPAGRM longitudinal (zstd) ────────────────────────────────────────

build/grab2 \
  --method SPAGRM \
  --pheno examples/long_pheno \
  --pheno-name Long1,Long2 \
  --longitudinal \
  --covar-name MALE,PC1,PC2 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT}.ibd.zst \
  --pfile examples/1kg \
  --out ${OUT_DIR}/long_spagrm \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
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
#
# Output columns: the nine meta columns, then LOG10P_CCT and five per-tau
# groups —
#   LOG10P_CCT  P_tau*  LOG10P_tau*  Z_tau*  Z_Norm_tau*  SPA_STATUS_tau*
# The saddlepoint is per tau, so both new columns are per tau too: a marker
# can converge at one quantile level and fail at another.  SPA_STATUS_tau*
# uses the same integer encoding as the SPACox and SPAGRM blocks.  The Cauchy
# combination is taken over the LOG10P_tau group (log10p_unify Stage 5): its
# statistic's terms are 1/(pi*p), so a linear combination overflowed and
# returned exactly 0 as soon as one tau reached p ~ 1e-308.  A tau whose test
# produced no answer is NA and drops out rather than poisoning the rest.

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
  --spasqr-tol 1e-6 \
  --spasqr-h-scale 3 \
  --spasqr-mode score \
  --pheno-transform int \
  --outlier-iqr-multiplier 1.5 \
  --spasqr-outlier-abs-bound 0.55 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── SPACox LOCO (fit path + per-chromosome LOCO PGS) ──────────────────
# For each chromosome the Cox null model is refit with that chromosome's
# LOCO PGS appended as an estimated covariate column; the augmented design
# also feeds the stage-2 covariate projection.  The pred.list is keyed on
# the canonical spec name (Cox: "Time_Event").  Output is plain text.

build/grab2 \
  --method SPACox \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pred-list examples/loco_prs.list2 \
  --pfile examples/1kg \
  --out ${LOCO} \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --chr 1-2,3 \
  --covar-p-threshold 5e-5 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1

## ── SPAmix LOCO (fit path + per-chromosome LOCO PGS) ──────────────────
# On-the-fly individual-AF model (no --ind-af-coef).  Only the residuals
# change per chromosome; the PC/OLS pools are built once.

build/grab2 \
  --method SPAmix \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --pc-cols PC1,PC2,PC3,PC4 \
  --pred-list examples/loco_prs.list2 \
  --pfile examples/1kg \
  --out ${LOCO} \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --chr 1-2,3 \
  --outlier-iqr-multiplier 1.5 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── SPAGRM LOCO (fit path + per-chromosome LOCO PGS) ──────────────────
# The retrospective GRM variance (R^T Phi R) is compatible with LOCO: the
# GRM topology / IBD are built once, and buildSPAGRMNullModel is re-run per
# chromosome from the refreshed residuals.

build/grab2 \
  --method SPAGRM \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT}.ibd.zst \
  --pred-list examples/loco_prs.list2 \
  --pfile examples/1kg \
  --out ${LOCO} \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── WtCoxG LOCO (fit path + per-chromosome LOCO PGS) ──────────────────
# Cox phenotype only (Binary has no bundled LOCO PGS).  Per chromosome the
# residuals and the batch-effect map (Phase D) are recomputed; the matched-
# marker scan (Phase C) is built once.

build/grab2 \
  --method WtCoxG \
  --pheno examples/1kg.pheno \
  --pheno-name Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --ref-af examples/ref_pop1.afreq \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --prevalence 0.1 \
  --pred-list examples/loco_prs.list2 \
  --pfile examples/1kg \
  --out ${LOCO} \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --chr 1-2,3 \
  `# NOT the built-in default, which is 0.1 (src/cli/cli.hpp).  0.05 is pinned` \
  `# here because it decides which batch-effect branch each marker takes, so` \
  `# changing it re-baselines every WtCoxG and LEAF output at once.` \
  --batch-effect-p-threshold 0.05 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression gz \
  --compression-level 6

## ── SPAsqr (wald mode, follow-up effect-size estimation) ──────────────
# Wald mode refits the joint smoothed-QR model with [X | G] per marker,
# emitting β̂_G + SE.  Per-marker QR refit is appreciably slower than
# score mode, so this block restricts to the 100 variant IDs in
# examples/spasqr_wald_extract.  Per-marker work runs on the shared
# marker-engine thread pool; chunk size auto-shrinks when --chunk-ksnp
# is left at its 8-ksnp default (8192 SNPs) so the worker pool stays fed even on
# small --extract subsets.  Output is plink2-style one-marker-per-line
# wide format (LOG10P_CCT + P_tau* + LOG10P_tau* + Z_tau* + BETA_tau* +
# SE_tau* + SPA_STATUS_tau* columns), written through TextWriter honoring
# --compression.  SPA_STATUS_tau is 1 (NORMAL) wherever the tau produced a
# test -- the Wald leg is a plain normal-reference z and never attempts a
# saddlepoint -- and 8 (NA_NO_TEST) where the sandwich variance is unusable
# (log10p_unify Stage 7).  A distinct --out
# prefix keeps the per-phenotype .SPAsqr.zst files from colliding with
# the score-mode artifacts produced above.

build/grab2 \
  --method SPAsqr \
  --pheno ${OUT_DIR}/int_pheno.txt \
  --pheno-name Quantitative,Time \
  --covar examples/1kg.pheno \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pred-list examples/loco_prs.list \
  --extract examples/spasqr_wald_extract \
  --pfile examples/1kg \
  --out ${OUT_DIR}/wald \
  `# Optional flags below (set to built-in defaults):` \
  --chr 1-2,3 \
  --spasqr-taus 0.1,0.3,0.5,0.7,0.9 \
  --spasqr-tol 1e-6 \
  --spasqr-h-scale 5 \
  --spasqr-mode wald \
  --pheno-transform int \
  --outlier-iqr-multiplier 1.5 \
  --spasqr-outlier-abs-bound 0.55 \
  --spa-z-threshold 2.0 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── WtCoxG ────────────────────────────────────────────────────────────
#
# Output columns (spa_unify Stage 6): the nine meta columns, then
#   P_EXT  LOG10P_EXT  P_NOEXT  LOG10P_NOEXT
#   Z_EXT  Z_NOEXT  Z_Norm_EXT  Z_Norm_NOEXT
#   P_BAT  PI_BAT  VAR_BAT  SPA_STATUS_EXT  SPA_STATUS_NOEXT
# LOG10P_* is -log10 of the p-value in the preceding column, computed in the
# log domain so it stays meaningful past the linear-scale underflow.
# SPA_STATUS_* is the spa::Status enumerator of the saddlepoint that produced
# that p-value, as an integer: 0 OK, 1 MAXITER, 2 GUARD_TEMP, 3 GUARD_CURV,
# 4 GUARD_W, 5 NONFINITE, 6 NORMAL (|Z| <= --spa-z-threshold, saddlepoint not
# attempted).  P and LOG10P are NA for every status other than 0 and 6.
#
# On the two batch-effect branches the reported P is a conditional probability
# assembled from a bivariate-normal integral rather than a saddlepoint tail,
# and the saddlepoint enters only through the variance recovered by inverting
# it; SPA_STATUS then describes that saddlepoint, and status 5 additionally
# covers a marker whose conditional integral does not exist because the
# assembled 2x2 covariance is not positive semi-definite (177 markers of
# Time_Event on this fixture, all in Branch B with sigma^2 >> var_Sbat).

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
  `# NOT the built-in default, which is 0.1 (src/cli/cli.hpp).  0.05 is pinned` \
  `# here because it decides which batch-effect branch each marker takes, so` \
  `# changing it re-baselines every WtCoxG and LEAF output at once.` \
  --batch-effect-p-threshold 0.05 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression gz \
  --compression-level 6

## ── LEAF ──────────────────────────────────────────────────────────────
#
# Output columns (spa_unify Stage 6): the nine meta columns, then
#   meta_P_EXT  meta_LOG10P_EXT  meta_P_NOEXT  meta_LOG10P_NOEXT
#   meta_SPA_STATUS_EXT  meta_SPA_STATUS_NOEXT
# and, per k-means cluster N,
#   clN_MAC  clN_P_EXT  clN_LOG10P_EXT  clN_P_NOEXT  clN_LOG10P_NOEXT
#   clN_P_BAT  clN_PI_BAT  clN_VAR_BAT  clN_SPA_STATUS_EXT  clN_SPA_STATUS_NOEXT
# The per-cluster columns are the WtCoxG block above, evaluated within that
# cluster; the encoding of SPA_STATUS is the same.  The META status describes
# the POOLING rather than the clusters: it is the worst status among the
# clusters that contributed to the fixed-effects pool, hence 0 or 6 whenever
# meta_P is a number and 5 when no cluster contributed and meta_P is NA.  A
# cluster with no informative subjects for a marker has no test, reports NA
# with status 5, and is simply left out of the pool -- that is the common case
# rather than the exception (2036 of 3000 markers for cluster 1 here).

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
  `# NOT the built-in default, which is 0.1 (src/cli/cli.hpp).  0.05 is pinned` \
  `# here because it decides which batch-effect branch each marker takes, so` \
  `# changing it re-baselines every WtCoxG and LEAF output at once.` \
  --batch-effect-p-threshold 0.05 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression gz \
  --compression-level 6

## ── SPAGxE (base, gene x environment interaction; env = MALE) ──────────
# Binary environment MALE (the paper's genetic-sex analyses) x three trait
# types (quantitative / survival / binary).  Every --envir-name column must
# also appear in --covar-name (it enters the genotype-independent null model
# trait ~ X + E).  Plain-text output exercises the uncompressed writer path.
#
# Output columns (log10p_unify Stage 7): the nine meta columns, the marginal
# block P_G Z_G BETA_G SE_G SPA_STATUS_G, then an eight-wide block per
# environment
#   P_Gx<E>  LOG10P_Gx<E>  LOG10P_Wald_Gx<E>  Z_Gx<E>  Z_Norm_Gx<E>
#   BETA_Gx<E>  SE_Gx<E>  SPA_STATUS_Gx<E>
# SPA_STATUS_Gx uses the same integer encoding as the SPACox and SPAGRM
# blocks above and always describes the saddlepoint leg, so a Branch-B
# marker whose SPA failed but whose Wald refit succeeded still shows a
# finite P_Gx together with the status that says the SPA dropped out of the
# Cauchy combination.  LOG10P_Gx is log-domain on every path: where the
# reported p is the saddlepoint p it comes from the tail assembly, and where
# Branch B adds a Wald leg the combination is taken over the two magnitudes
# (Stage 5), both of which are produced as magnitudes at the source (Stage 7).
# The marginal block is always the normal approximation, so it carries no
# LOG10P_G, and SPA_STATUS_G is the constant 1 (NORMAL) wherever Var(S_G) > 0
# and 8 (NA_NO_TEST) where it is not.
#
# A distinct --out prefix keeps the fitted residual and per-marker tables from
# colliding with the SPACox/SAGELD ${OUT}.* artifacts.

build/grab2 \
  --method SPAGxE \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --envir-name MALE \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spagxe \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --save-resid \
  --spagxe-marginal-cutoff 0.001 \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1

## ── SPAGxE+ (sparse-GRM relatedness correction; --sp-grm-plink2) ───────
# Passing an optional sparse GRM engages the SPAGxE+ variance correction (a
# retrospective GRM quadratic form; no --pairwise-ibd needed).  The GRM path
# keeps no Branch-B Wald leg (paper), so LOG10P_Wald_Gx<E> is NA throughout.  gzip
# output exercises the gz writer path.

build/grab2 \
  --method SPAGxE \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --envir-name MALE \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spagxe_plus \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --spagxe-marginal-cutoff 0.001 \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression gz \
  --compression-level 6

## ── SPAGxEmix (per-individual allele frequency; --pc-cols) ─────────────
# SPAGxEmix estimates a per-individual ALT frequency q_i from the --pc-cols
# principal components (the SPAmix cascade), so the retrospective genotype law
# is Binomial(2, q_i).  The --pc-cols columns must also appear in --covar-name
# (they adjust the null model), as must every --envir-name column.  A sparse
# GRM is not accepted (SPAGxEmix+ is out of scope).  zstd output exercises the
# zst writer path, completing the plain / gz / zst rotation across the three
# G x E blocks.  Output columns are those of the SPAGxE block above.
#
# Three of the 3000 markers here (rs3131973, rs300788, rs9310509; ALT_FREQ
# 0.992-0.996) are NA from P_G onward, and have been since before the
# saddlepoint rework.  The cause is upstream of the saddlepoint: their AF
# model takes the logistic branch, the binarised genotype is one for
# essentially every subject, IRLS runs to its cap under complete separation,
# and sigmoid of the resulting linear predictor rounds to exactly 1, so
# every q_i is 1, every 2q_i(1-q_i) is 0 and Var(S_G) vanishes.  The values
# are preserved and the row states why: SPA_STATUS_GxMALE = 8 (NA_NO_TEST)
# since the Stage 2 re-partition, and SPA_STATUS_G = 8 as well since Stage 7
# gave the marginal block a status of its own.  The same three markers are ordinary under the uniform-q
# SPAGxE block above, which does not use the AF model.

build/grab2 \
  --method SPAGxEmix \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time:Event,Binary \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --envir-name MALE \
  --pc-cols PC1,PC2,PC3,PC4 \
  --pfile examples/1kg \
  --out ${OUT_DIR}/spagxemix \
  `# Optional flags below (set to built-in defaults):` \
  --regression-model auto \
  --spagxe-marginal-cutoff 0.001 \
  --chr 1-2,3 \
  --spa-z-threshold 2.0 \
  --outlier-iqr-multiplier 1.5 \
  --threads 2 \
  --chunk-ksnp 8 \
  --geno 0.1 \
  --maf 1e-5 \
  --mac 10 \
  --hwe 0 \
  --hard-call-threshold 0.1 \
  --compression zst \
  --compression-level 3

## ── Regression cross-checks (md5/diff over compressed outputs) ────────
# The blocks below are not new analyses: they re-run a method on
# equivalent inputs (or in resid mode after the fit mode has saved its
# null-model residuals) and assert that the per-marker output table is
# bit-identical to the reference run.  Each comparison strips compression
# (zstdcat / zcat) before hashing so that codec framing does not mask
# content equivalence.  Failures are reported as "[FAIL]" but the script
# continues so that all checks run in one pass.

# Helper: md5sum every member of the list (raw bytes, no decompression)
# and report whether every member shares the same hash with the first.
# This relies on grab2's compressed writers (zstd, gzip) being byte-
# deterministic for identical inputs at a fixed codec level, which has
# been verified for the zstd path.  Always returns 0 so `set -e` does
# not abort the script when a check fails.
md5_equiv() {
  local label="$1"; shift
  local ref="" status="PASS" md
  for f in "$@"; do
    md=$(md5sum "${f}" | awk '{print $1}')
    printf "    %s  %s\n" "${md}" "${f}"
    if [ -z "${ref}" ]; then
      ref="${md}"
    elif [ "${md}" != "${ref}" ]; then
      status="FAIL"
    fi
  done
  echo "  [${status}] ${label}"
  echo
  return 0
}

# Helper: numeric_equiv TOL LABEL FILE_A FILE_B — compare two decompressed
# tables cell by cell, requiring exact equality of every non-numeric cell and
# agreement to relative tolerance TOL for every numeric one.
#
# Used ONLY for the SAGELD fit-vs-resid pair, and only because md5_equiv is the
# wrong instrument there.  Measured with a 17-significant-digit diagnostic build
# at spa_unify Stage 4 (see the block comment at that check): the two modes
# already disagreed in the last bits of P_GxTIME for 1423 of 3000 Long1 markers
# and 617 of 3000 Long2 markers BEFORE Stage 4, by up to 4.4e-15 relative.  Byte
# identity was therefore never a property of the computation, only of the
# 6-significant-figure print rounding both values to the same string.  Always
# returns 0 so `set -e` does not abort the script when a check fails.
numeric_equiv() {
  local tol="$1" label="$2" a="$3" b="$4"
  local status
  status=$(zstdcat "$a" > "${a}.__ne_a" 2>/dev/null || cat "$a" > "${a}.__ne_a"
           zstdcat "$b" > "${a}.__ne_b" 2>/dev/null || cat "$b" > "${a}.__ne_b"
    awk -v TOL="$tol" '
      function isnum(x) { return (x ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][+-]?[0-9]+)?$/) }
      NR==FNR { for (i=1;i<=NF;i++) A[FNR,i]=$i; NC[FNR]=NF; NRA=FNR; next }
      { if (NF != NC[FNR]) { bad++; next }
        for (i=1;i<=NF;i++) {
          x=A[FNR,i]; y=$i
          if (x==y) continue
          if (isnum(x) && isnum(y)) {
            d = x-y; if (d<0) d=-d
            m = (x<0?-x:x); if ((y<0?-y:y) > m) m = (y<0?-y:y)
            if (m>0 && d/m <= TOL+0) continue
            if (d <= TOL+0) continue
          }
          bad++; if (worst=="") worst=sprintf("row %d col %d: %s vs %s", FNR, i, x, y)
        } }
      END { if (FNR != NRA) { print "FAIL row-count"; exit }
            if (bad) printf "FAIL %d cell(s), first %s\n", bad, worst
            else print "PASS" }
    ' "${a}.__ne_a" "${a}.__ne_b")
  rm -f "${a}.__ne_a" "${a}.__ne_b"
  printf "    rtol %s  %s\n    rtol %s  %s\n" "$tol" "$a" "$tol" "$b"
  echo "  [${status}] ${label}"
  echo
  return 0
}

echo
echo "════════════════════════════════════════════════════════════════════"
echo "Regression cross-checks"
echo "════════════════════════════════════════════════════════════════════"

## ── Cross-format SPAGRM equivalence: pgen / bed / bcf / bgen ──────────
# Convert the bundled pgen fixture to BED, BCF, and BGEN with plink2,
# then run SPAGRM on each input format with a distinct --out prefix and
# verify that the per-phenotype output tables agree byte-for-byte across
# all four genotype readers in src/geno_factory/.
#
# This block also exercises --extract / --exclude / --keep / --remove
# uniformly across all four readers.  Because every run consumes the
# same four filter lists, the md5 cross-check at the bottom of this
# block verifies that variant-ID filtering and sample-IID filtering
# produce identical output across pgen, bed, bcf, and bgen.

CONV_DIR=${OUT_DIR}/converted
mkdir -p ${CONV_DIR}

plink2 --pfile examples/1kg --make-bed                       --out ${CONV_DIR}/1kg
plink2 --pfile examples/1kg --export bcf                     --out ${CONV_DIR}/1kg
plink2 --pfile examples/1kg --export bgen-1.2 bits=8         --out ${CONV_DIR}/1kg

# ── Filter lists for --extract / --exclude / --keep / --remove ──
# Variant lists: every 3rd .pvar variant goes into extract; the first
# ten IDs of the .pvar are also in exclude, so extract∩exclude is a
# proper subset and tests "exclude wins over extract".
awk 'NR>1 && $1 !~ /^#/ && (NR-1)%3==0 {print $3}' examples/1kg.pvar > ${OUT_DIR}/1kg.extract.txt
awk 'NR>1 && $1 !~ /^#/ && (NR-1)<=10        {print $3}' examples/1kg.pvar > ${OUT_DIR}/1kg.exclude.txt

# Subject lists: keep the first 2000 .psam IIDs (PLINK2-compatible
# FID+IID two-column format) and then remove the last 100 of those.
awk 'NR>1 && $1 !~ /^#/ {print $1"\t"$1}' examples/1kg.psam | head -2000 > ${OUT_DIR}/1kg.keep.txt
awk 'NR>1 && $1 !~ /^#/ {print $1"\t"$1}' examples/1kg.psam | sed -n '1901,2000p' > ${OUT_DIR}/1kg.remove.txt

# Shared SPAGRM invocation; only --pfile/--bfile/--bcf/--bgen and --out
# vary between the four runs below.
SPAGRM_COMMON=(
  --method SPAGRM
  --pheno examples/1kg.pheno
  --pheno-name Quantitative,Time:Event,Binary,Ordinal
  --covar-name MALE,PC1,PC2,PC3,PC4
  --sp-grm-plink2 examples/1kg.grm.sp
  --pairwise-ibd ${OUT}.ibd.zst
  --extract ${OUT_DIR}/1kg.extract.txt
  --exclude ${OUT_DIR}/1kg.exclude.txt
  --keep    ${OUT_DIR}/1kg.keep.txt
  --remove  ${OUT_DIR}/1kg.remove.txt
  --regression-model auto
  --chr 1-2,3
  --spa-z-threshold 2.0
  --outlier-iqr-multiplier 1.5
  --threads 2
  --chunk-ksnp 8
  --geno 0.1
  --maf 1e-5
  --mac 10
  --hwe 0 \
  --hard-call-threshold 0.1
  --compression zst
  --compression-level 3
)

build/grab2 "${SPAGRM_COMMON[@]}" --pfile examples/1kg                  --out ${OUT_DIR}/spagrm_pgen
build/grab2 "${SPAGRM_COMMON[@]}" --bfile ${CONV_DIR}/1kg               --out ${OUT_DIR}/spagrm_bed
build/grab2 "${SPAGRM_COMMON[@]}" --bcf   ${CONV_DIR}/1kg.bcf           --out ${OUT_DIR}/spagrm_bcf
build/grab2 "${SPAGRM_COMMON[@]}" --bgen  ${CONV_DIR}/1kg.bgen ref-last --out ${OUT_DIR}/spagrm_bgen

for phen in Quantitative Time_Event Binary Ordinal; do
  md5_equiv "SPAGRM cross-format ${phen}" \
    ${OUT_DIR}/spagrm_pgen.${phen}.SPAGRM.zst \
    ${OUT_DIR}/spagrm_bed.${phen}.SPAGRM.zst \
    ${OUT_DIR}/spagrm_bcf.${phen}.SPAGRM.zst \
    ${OUT_DIR}/spagrm_bgen.${phen}.SPAGRM.zst
done

## ── SPAmix: pre-computed AF vs on-the-fly AF ──────────────────────────
# The pre-computed-AF run (--ind-af-coef ${OUT}.afc.zst) and the
# on-the-fly run (no --ind-af-coef) refit the same logistic AF model on
# the same N subjects when the dataset has no missing genotypes, so
# their per-phenotype output tables must be byte-identical.

for phen in Quantitative Time_Event Binary Ordinal; do
  md5_equiv "SPAmix precomputed-vs-onthefly ${phen}" \
    ${OUT}.${phen}.SPAmix.zst \
    ${OUT_DIR}/spamix_otf.${phen}.SPAmix.zst
done

## ── SPACox: fit mode vs residual mode ─────────────────────────────────
# The fit-mode run above wrote ${OUT}.null.resid; the residual-mode run
# replayed the same scan against the same .pgen using those residuals.
# Per-phenotype tables must match exactly.

for phen in Quantitative Time_Event Binary Ordinal; do
  md5_equiv "SPACox fit-vs-resid ${phen}" \
    ${OUT}.${phen}.SPACox \
    ${RESID_OUT}.${phen}.SPACox
done

## ── SAGELD: fit mode vs residual mode (Long1 and Long2) ───────────────
# Each SAGELD residual-mode run consumes ${OUT}.${pheno}.SAGELD.resid.zst
# (fit mode compressed it because --compression zst was set) and emits
# ${RESID_OUT}_${pheno}.SAGELD.zst; the corresponding fit-mode result is
# ${OUT}.${pheno}.SAGELD.zst.
#
# This is the ONE cross-check in this script that is a tolerance comparison
# rather than md5_equiv, and the reason is measured, not assumed.  Residual mode
# does not merely replay fit mode: it rebuilds the GALLOP projection cache (Q and
# its factorization, Si, StS, XTs, AtS) from the ##sageld-cache-* header block
# and the Si_00.. columns of the residual file, and re-derives the quantities the
# fit had in registers.  A diagnostic build printing 17 significant digits shows
# that at the PRE-Stage-4 revision the two modes already disagreed in P_GxTIME
# for 1423 of 3000 Long1 markers and 617 of 3000 Long2 markers, by up to
# 4.416e-15 relative (Long1) and 2.373e-15 (Long2).  The md5 check passed only
# because 6-significant-figure printing rounded both sides to the same string.
#
# spa_unify Stage 4 replaced SPAGRM's Family-B Newton iteration with the shared
# bracketed safeguarded solver, whose accepted iterate is determined only to
# within the configured tolerance (--spagrm tol, 1e-6, relative on the residual).
# Two solves whose inputs differ in the last bits may therefore accept different
# iterates inside that tolerance, which lifts the pre-existing fit-vs-resid
# disagreement from ~4e-15 to ~5e-7 relative — still an order of magnitude
# tighter than the base-to-Stage-4 change of 8.6e-6, and far below any
# statistical threshold, but now large enough for 5 markers to straddle a print
# boundary.  Byte identity between the two modes cannot be restored without
# making the residual-mode cache reconstruction bit-exact, which is a SAGELD
# question and not a saddlepoint one.
#
# rtol 1e-5 is chosen to sit an order of magnitude above the measured 5e-7 and an
# order of magnitude below the 1e-4 at which a difference would start to matter
# for a reported p-value.  Every non-numeric cell must still match exactly, and
# every other cross-check in this script remains md5_equiv.

for sageld_pheno in Long1 Long2; do
  numeric_equiv 1e-5 "SAGELD fit-vs-resid ${sageld_pheno}" \
    ${OUT}.${sageld_pheno}.SAGELD.zst \
    ${RESID_OUT}_${sageld_pheno}.SAGELD.zst
done
