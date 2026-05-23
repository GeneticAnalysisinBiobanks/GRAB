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
#
# Notes on defaults shown below:
#
#   --chr 1-2,3              The bundled 1kg fixture contains markers only
#                            on chromosomes 1, 2, 3; the syntax "1-2,3"
#                            combines a range and a singleton and resolves
#                            to the set {1, 2, 3}.
#   --compression {gz,zst}   Output compression codec.  Methods are split
#                            across all three settings so the script
#                            exercises plain text (SPACox), gzip (LEAF,
#                            WtCoxG), and zstd (everything else).  For
#                            SPACox the --compression / --compression-level
#                            flags are deliberately omitted from the
#                            optional block to leave the output as plain
#                            text.
#   --compression-level INT  zstd default is 3; gzip default is 6.  The
#                            script passes each level explicitly so that
#                            the "Options in effect" log matches the
#                            value actually applied.
#   --save-resid             Parameterless switch that writes fitted
#                            residuals.  Per-method suffix avoids
#                            collisions across invocations and the
#                            residual-mode block of the same method
#                            reloads the file, exercising the fit →
#                            save → reload round trip:
#                              SPACox  → ${OUT}.null.resid
#                                        (one column per phenotype;
#                                        re-read by the SPACox
#                                        residual-mode block via
#                                        --pheno + --resid-name.)
#                              SAGELD  → ${OUT}.<pheno>.SAGELD.resid
#                                        (per-phenotype triplet R_G /
#                                        R_<E> / R_Gx<E>; the
#                                        residual-mode SAGELD block
#                                        consumes the Long1 file.)
#                            SPAmix / SPAGRM / SPAmixPlus /
#                            SPAmixLocalPlus also share the SPACox
#                            null.resid suffix; including --save-resid
#                            in their fit-mode blocks would overwrite
#                            the SPACox file, so the flag is omitted
#                            there.
#
# Two categories of flags are intentionally omitted because they have no
# default value and require an external companion file:
#
#   --keep FILE         --remove FILE
#   --extract FILE      --exclude FILE
#
# Add them by hand when a real subject or SNP subset restriction is
# required.
#
# Mutually exclusive flag groups are resolved by choosing one canonical
# default:
#   - Genotype input             --pfile
#   - Sparse GRM input           --sp-grm-plink2
#   - --pheno-name vs --resid-name: only SPACox and SAGELD are run in
#     both modes (the residual-mode block reloads the .resid file
#     produced by the fit-mode block above it).  SPAmix and SPAGRM
#     run only in fit-mode; their residual-mode invocations would
#     duplicate coverage already provided by SPACox.
#   - SPAsqr bandwidth           --spasqr-h-scale 3  (the score-mode default;
#                                --spasqr-h is unset)

set -e

OUT=examples_output/fit          # output prefix for fit-mode blocks
RESID_OUT=examples_output/resid  # output prefix for residual-mode blocks (SPACox, SAGELD)

mkdir -p examples_output

## ── Utility: cal-af-coef ──────────────────────────────────────────────
# Produces ${OUT}.afc.zst, consumed by SPAmix / SPAmixPlus via
# --ind-af-coef.

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
# Produces ${OUT}.ibd.zst, consumed by SPAGRM / SAGELD via
# --pairwise-ibd.

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
# Default --regression-model is 'auto': per-token inference from values
# and from the TIME:EVENT colon syntax.  Output left uncompressed
# (--compression omitted) to exercise the plain-text writer path.  This
# is the one block that keeps --save-resid: it writes ${OUT}.null.resid
# and no later block overwrites it.

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
# block above (--save-resid → ${OUT}.null.resid).  The file's data
# columns share the phenotype names from --pheno-name; covariates for
# the SPA adjustment are loaded separately from the original phenotype
# file via --covar.  Output is written under ${RESID_OUT} so the
# reloaded-residual outputs can be diffed against the fit-mode
# outputs under ${OUT}.

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
# Add --spagrm-control-outlier (parameterless) to enable iterative
# IQR-ratio adjustment that keeps the outlier share in (0, 5%].

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

## ── SAGELD (pheno mode: --pheno-name + --sageld-x) ────────────────────
# SAGELD does not accept --regression-model (auto is implicit).

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
# Consumes the per-phenotype residual file produced by the SAGELD
# pheno-mode block above (--save-resid → ${OUT}.Long1.SAGELD.resid);
# the column layout (R_G, R_TIME, R_GxTIME) is exactly what
# --resid-name expects in residual-input mode.  Output is written
# under ${RESID_OUT} so the reloaded-residual SAGELD output can be
# diffed against the pheno-mode Long1.SAGELD.zst output under ${OUT}.

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
# INT-transformed Quantitative and Time columns only.  The downstream
# SPAsqr block consumes this file via --pheno and loads the remaining
# covariates from the original phenotype file via --covar.

build/grab2 \
  --int-pheno \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time \
  --out ${OUT}.int

## ── SPAsqr (score mode, fit path) ─────────────────────────────────────
# Bandwidth: pick --spasqr-h-scale (the IQR divisor, default 3 in score
# mode).  --spasqr-h is left unset; the two flags are mutually
# exclusive.  SPAsqr does not consume --save-resid.  The phenotype file
# is the int-pheno output (Quantitative, Time only); MALE / PC1..PC4 are
# loaded from the original phenotype file via --covar.

build/grab2 \
  --method SPAsqr \
  --pheno examples_output/1kg.int.txt \
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
# WtCoxG does not consume --save-resid.  Output compressed with gzip to
# exercise the zlib writer path.

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
# --leaf-cluster-file is mutually exclusive with K-means restarts.
# This block uses internal K-means; supply --leaf-cluster-file FILE to
# skip K-means and read pre-computed cluster labels instead.  LEAF does
# not consume --save-resid.  Output compressed with gzip.

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
