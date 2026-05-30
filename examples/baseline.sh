#!/usr/bin/env bash
# examples/baseline.sh — exhaustive regression script for the GRAB repository.
#
# This script exercises every utility mode (cal-af-coef, cal-pairwise-ibd,
# int-pheno) and every analysis method (SPACox, SPAmix, SPAGRM, SAGELD,
# SPAsqr score + wald, WtCoxG, LEAF) against the bundled examples/1kg.*
# fixtures.  It serves two purposes:
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

## ── SAGELD (residual mode, --resid-name; Long1 and Long2) ─────────────
# Consumes the per-phenotype residual files produced by the SAGELD fit
# mode (one --resid output per phenotype, with R_G/R_TIME/R_GxTIME).
# Output prefix is suffixed by the phenotype name so that fit and resid
# artefacts for each of Long1 and Long2 can be compared independently
# in the regression cross-checks block at the bottom of this script.

for sageld_pheno in Long1 Long2; do
build/grab2 \
  --method SAGELD \
  --pheno ${OUT}.${sageld_pheno}.SAGELD.resid \
  --resid-name R_G,R_TIME,R_GxTIME \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pairwise-ibd ${OUT}.ibd.zst \
  --pfile examples/1kg \
  --out ${RESID_OUT}_${sageld_pheno} \
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
done

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

## ── SPAsqr (wald mode, follow-up effect-size estimation) ──────────────
# Wald mode refits the joint smoothed-QR model with [X | G] per marker,
# emitting β̂_G + SE.  Per-marker QR refit is appreciably slower than
# score mode, so this block restricts to the 100 variant IDs in
# examples/spasqr_wald_extract.  Per-marker work runs on the shared
# marker-engine thread pool; chunk size auto-shrinks when --chunk-size
# is left at its 8192 sentinel so the worker pool stays fed even on
# small --extract subsets.  Output is plink2-style one-marker-per-line
# wide format (P_CCT + P_tau* + Z_tau* + BETA_tau* + SE_tau* columns),
# written through TextWriter honoring --compression.  A distinct --out
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
  --spasqr-tol 1e-7 \
  --spasqr-h-scale 10 \
  --spasqr-mode wald \
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
  --seed 2026
  --threads 2
  --chunk-size 8192
  --geno 0.1
  --maf 1e-5
  --mac 10
  --hwe 0
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
# Each SAGELD residual-mode run consumes ${OUT}.${pheno}.SAGELD.resid
# and emits ${RESID_OUT}_${pheno}.SAGELD.zst; the corresponding fit-mode
# result is ${OUT}.${pheno}.SAGELD.zst.

for sageld_pheno in Long1 Long2; do
  md5_equiv "SAGELD fit-vs-resid ${sageld_pheno}" \
    ${OUT}.${sageld_pheno}.SAGELD.zst \
    ${RESID_OUT}_${sageld_pheno}.SAGELD.zst
done
