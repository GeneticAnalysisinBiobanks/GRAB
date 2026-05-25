# Performance and Workflow Reference

This page describes the internal workflows of every GRAB analysis mode
at function level, with notes on parallelization, threading, and the
optimization choices baked into the codebase.  The companion document
[engines.md](engines.md) describes the shared marker-streaming
infrastructure that all GWAS methods reuse; this page focuses on the
per-method specifics.

---

## Complexity and parallelization overview

### Notation

| Symbol  | Meaning                                | Typical range     |
| ------- | -------------------------------------- | ----------------- |
| N       | Subjects                               | 100K – 1M         |
| M       | Markers                                | 1M – 10M          |
| E       | GRM non-zero entries                   | 2–5 × N           |
| C       | Covariates / PCs                       | 5 – 20            |
| P       | Phenotypes                             | 1 – 20            |
| K_τ     | Quantile levels (SPAsqr)               | 5 – 9             |
| K_env   | Environments (SAGELD)                  | 1 – 5             |
| K_anc   | Ancestries (SPAmixLocalPlus)           | 2 – 3             |
| K_clust | Clusters (LEAF)                        | 3 – 5             |
| F       | Families; K_fam = max family size      | K_fam ≤ 10        |
| T       | Worker threads                         | 8 – 64            |
| M_ref   | Matched reference markers              | ~ M               |

### Per-method summary

All per-marker costs are stated **per marker**; total marker-phase wall
time ≈ (cost × M) / T.  Pre-marker stages run once (or once per
phenotype) before marker streaming begins.

| Method               | Phase        | Dominant stage                                | CPU cost                        | Parallel                                | Notes                                                                |
| -------------------- | ------------ | --------------------------------------------- | ------------------------------- | --------------------------------------- | -------------------------------------------------------------------- |
| **SPACox**           | Pre          | CGF interpolation grid                        | O(N)                            | Serial                                  | CGF table per phenotype (residual-dependent)                         |
|                      | Per marker   | Score + SPA tail                              | O(N)                            | Chunk ∥                                | SPA tail evaluated only when `\|Z\| > spaCutoff`                     |
| **SPAGRM**           | Pre          | Chow–Liu tree + MAF-bin tables                | O(E + F · 3^K_fam)              | Outlier families serial; bins ∥        | K_fam ≤ 10 bounds the exponent                                       |
|                      | Per marker   | Score + GRM variance + SPA                    | O(N + E)                        | Chunk ∥                                | O(1) MAF-bin interpolation                                           |
| **SAGELD**           | Pre          | K_env × SPAGRM null model                     | O(K_env · (E + F · 3^K_fam))    | Envs `min(T, K_env)` ∥; inner T/K_env ∥| Outer / inner thread split rebalanced from `--threads`               |
|                      | Per marker   | K_env × SPAGRM test                           | O(K_env · (N + E))              | Chunk ∥                                |                                                                      |
| **SPAmixPlus**       | Pre          | OLS matrices                                  | O(N · C²)                       | Serial                                  | `(X'X)⁻¹` deduped across phenotypes sharing the same NA pattern       |
|                      | Per marker   | Score + GRM variance                          | O(N + E)                        | Chunk ∥                                | O(N) without GRM                                                     |
| **SPAmixLocalPlus**  | Pre          | `MultiPhenoRprodSoA` build                    | O(E · K_anc)                    | Serial                                  | One pass per ancestry; multi-phenotype packing baked in              |
|                      | Per marker   | Batch score + phi variance                    | O(N · K_anc + E)                | Chunk ∥ + runtime SIMD                 | `PHI_BATCH = 8`                                                      |
| **WtCoxG**           | Phase A      | Shared loading + ref-AF + GRM                 | O(M_ref + N · C²)               | Serial                                  | `runWtCoxG` multi-pheno entry                                         |
|                      | Phase B      | P null models                                 | O(P · N · C²)                   | `min(T, P)` ∥                          | Parallel regression                                                  |
|                      | Phase C      | Single geno scan with P indicators            | O(M_ref · N · P)                | Serial                                  | Replaces P × `computeMarkerStats`                                    |
|                      | Phase D      | P batch-effect tests                          | O(P · M_ref · N)                | `min(T, P)` ∥                          | GRM shared                                                           |
|                      | Per marker   | P × score + bivariate SPA                     | O(P · N)                        | Chunk ∥                                | Single `multiPhenoEngine`                                            |
| **LEAF**             | Phase A      | Shared loading + K-means + GRM × K            | O(M_ref + 300 · N · K · C + K · GRM) | `min(T, 25)` ∥                    | K-means restarts in parallel                                         |
|                      | Phase B      | P × K null models                             | O(P · K · N · C²)               | `min(T, P · K)` ∥                       | Per (p, c) atomic work-stealing                                      |
|                      | Phase C      | Single geno scan + K summix                   | O(M_ref · N · P · K)            | `min(T, K)` ∥ summix                    | `intAF` shared; AF_ref synthesis shared                              |
|                      | Phase D      | P × K batch-effect tests                      | O(P · K · M_ref · N)            | `min(T, P · K)` ∥                       | GRMs preloaded                                                       |
|                      | Per marker   | P × (K tests + CCT)                           | O(P · K · N)                    | Chunk ∥                                | Single `multiPhenoEngine`                                            |
| **SPAsqr** (score)   | Pre (1)      | P × K_τ conquer / QMME QR                     | O(P · K_τ · N · C²)             | `min(T, P · K_τ)` ∥                     | Atomic work-stealing; GRM and genotype shared                        |
|                      | Pre (2)      | P × K_τ × GRM variance + outlier detect       | O(P · K_τ · E)                  | Serial                                  | No family work; negligible vs Pre (1)                                |
|                      | Per marker   | K_τ × SPAGRM test + CCT                       | O(K_τ · (N + E))                | Chunk ∥                                | Per-phenotype via `multiPhenoEngine`                                 |
| **SPAsqr** (LOCO)    | Pre          | same as score                                 |                                 |                                         | `locoEngine` rebuilds K_τ × P null models per chromosome             |
| **SPAsqr** (Wald)    | Per marker   | K_τ × full-model QMME refit + sandwich SE     | O(K_τ · N · C³)                 | Chunk ∥                                | No GRM; long-format output                                           |
| **`--cal-af-coef`**  | Pre          | Design matrix X                               | O(N · C²)                       | Serial                                  | One-time                                                             |
|                      | Per marker   | OLS / logistic AF model                       | O(N · C)                        | Chunk ∥                                | I/O bound                                                            |
| **`--cal-pairwise-ibd`** | Pre      | GRM pair index                                | O(E)                            | Serial                                  | One-time                                                             |
|                      | Per marker   | Pair accumulation                             | O(E_pairs)                      | Chunk ∥ + runtime SIMD                 | 4–8 pairs / cycle                                                    |
| **`--cal-phi`**      | Pre          | GRM + AdmixData load                          | O(E + M · K_anc)                | Serial                                  |                                                                      |
|                      | Per marker   | Per-ancestry phi accumulation                 | O(E_pairs · K_anc)              | Chunk ∥                                |                                                                      |

### Optimizations already applied

| Class      | Mechanism                                                                                                                                                       | Affected methods                                                                 | Source                                                  |
| ---------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- | ------------------------------------------------------- |
| Threading  | Chunk-level work-stealing via `atomic<size_t> nextChunk`; ordered writer thread                                                                                  | All GWAS methods, all utility modes                                              | [src/engine/marker.cpp](../../src/engine/marker.cpp)    |
| Threading  | LOCO outer loop wraps the multi-phenotype inner loop; one chromosome's chunks active at a time                                                                   | SPAsqr LOCO                                                                      | [src/engine/loco.cpp](../../src/engine/loco.cpp)        |
| Threading  | Shared genotype decode across phenotypes through the union mask (one decode replaces K decodes)                                                                  | SPACox, SPAGRM, SPAmixPlus, WtCoxG, LEAF, SPAsqr, SAGELD (pheno mode)            | `multiPhenoEngine`                                      |
| Threading  | Atomic work-stealing across the P × K_τ conquer / QMME fits (Pre-stage 1)                                                                                        | SPAsqr                                                                           | `runSPAsqr`                                             |
| Threading  | Atomic work-stealing across the P × K cluster null models and batch-effect tests                                                                                 | LEAF                                                                             | `runLEAF`                                               |
| Threading  | Parallel regression across P null models with `min(T, P)` threads                                                                                                | WtCoxG                                                                           | `runWtCoxG`                                             |
| Threading  | Parallel K-means restarts (`min(T, 25)`) and parallel summix (`min(T, K)`)                                                                                       | LEAF                                                                             | `runLEAF`                                               |
| Threading  | Atomic work-stealing across the MAF-bin probability tables of each outlier family; per-bin work is independent (each writes its own column of CLT and owns its `entropyMat` scratch) | SPAGRM, SAGELD, SPAsqr                                                           | `buildChowLiuTree` in [grm_null.cpp](../../src/spagrm/grm_null.cpp) |
| Threading  | Atomic work-stealing across K_env null-model constructions; inner thread budget is rebalanced to `max(1, T / K_env)` per builder                                  | SAGELD                                                                           | `runSAGELD` in [sageld.cpp](../../src/spagrm/sageld.cpp) |
| Threading  | Parallel null-model fits across phenotypes via `nullmodel::fitAll`; respects `--threads`                                                                          | SPACox / SPAGRM / SPAmix / SPAmixPlus / SPAmixLocalPlus / WtCoxG / LEAF (fit path)| [src/util/null_model.cpp](../../src/util/null_model.cpp) |
| Threading  | Per-ancestry per-batch per-phenotype score computed as three fused Eigen GEMMs (`S = bDosBig^T·R_mat`, `HR = bHapBig^T·R_mat`, `HR2 = bHapBig^T·R2_mat`)         | SPAmixLocalPlus                                                                  | `runUnifiedGWAS` in [spamixlocalp.cpp](../../src/localplus/spamixlocalp.cpp) |
| Numerical  | Pre-computed CGF interpolation grid (10,000 points) eliminates repeated `exp` calls inside the SPA root-finding loop                                              | SPACox                                                                           | `SPACoxMethod`                                          |
| Numerical  | MAF-bin interpolation (11 hard-coded bins) reduces per-marker family probability work to O(1)                                                                    | SPAGRM, SAGELD, SPAsqr                                                           | `SPAGRMClass`                                           |
| Numerical  | Chow–Liu tree caps family joint distribution at `3^K_fam` with `K_fam ≤ 10`                                                                                       | SPAGRM, SAGELD, SPAsqr                                                           | `SPAGRMClass`                                           |
| Numerical  | SPA tail computation gated by `\|Z\| > spaCutoff`; normal approximation is used otherwise                                                                        | Every SPA method                                                                 | All SPA paths                                           |
| Numerical  | Covariate-adjusted recomputation gated by `p < pvalCovAdjCut`                                                                                                    | SPACox                                                                           | `SPACoxMethod`                                          |
| Numerical  | Design-matrix inverse `(X'X)⁻¹`, PC matrices, and GRM sub-matrices deduplicated across phenotypes sharing the same non-missingness pattern                       | SPACox, SPAmixPlus                                                               | `runSPACox`, `runSPAmixPlus`                            |
| Numerical  | QMME caches the Cholesky decomposition of the Hessian upper bound once per phenotype × bandwidth and reuses it across all τ levels                                | SPAsqr (`--spasqr-solver qmme`)                                                  | `runSPAsqr`                                             |
| Memory     | Fused union-level GEMM (`AugResid^T × G`) computes scores for every fuseable phenotype in one batched matrix product; eliminates per-phenotype extraction        | SPAGRM, SPAmixPlus, WtCoxG, LEAF, SPAsqr, SAGELD (pheno mode)                    | `multiPhenoEngineRange` in [src/engine/marker.cpp](../../src/engine/marker.cpp) |
| Memory     | `MissBatch` groups non-fuseable phenotypes by identical missingness patterns; one extraction per pattern replaces one extraction per phenotype                    | SPACox                                                                           | `multiPhenoEngineRange`                                 |
| Memory     | Structure-of-Arrays phi layout (`MultiPhenoRprodSoA`) for cache-friendly variance scans, with `R[i]·R[j]·phi` pre-multiplied across phenotypes                    | SPAmixLocalPlus                                                                  | [src/localplus/spamixlocalp.cpp](../../src/localplus/spamixlocalp.cpp) |
| Memory     | Mini-batch of `PHI_BATCH = 8` markers scans phi entries once for the whole batch, amortizing L3 misses ≈ 8×                                                       | SPAmixLocalPlus                                                                  | `runUnifiedGWAS`                                        |
| Memory     | Pre-computed phi tables reduce per-marker variance from O(E · K_anc) to O(E)                                                                                     | SPAmixLocalPlus                                                                  | [src/localplus/spamixlocalp.cpp](../../src/localplus/spamixlocalp.cpp) |
| SIMD       | Runtime dispatch (AVX-512 / AVX2 / scalar) resolved once at process startup via `__builtin_cpu_supports`                                                          | SPAmixLocalPlus variance, pairwise-IBD accumulation, ABED decode/encode          | [src/util/simd_dispatch.hpp](../../src/util/simd_dispatch.hpp) |
| SIMD       | `popcnt` / `ctz` via GCC built-ins (part of the x86-64-v2 baseline)                                                                                              | LEAF cluster bitmask operations                                                  | [src/wtcoxg/leaf.cpp](../../src/wtcoxg/leaf.cpp)        |
| SIMD       | Eigen's blocked GEMM kernel drives all per-method matrix products; SSE/AVX/NEON intrinsics are bundled inside Eigen and require no external runtime               | Every method that uses Eigen `matrix.transpose() * matrix` patterns              | Eigen (vendored, header-only)                           |
| I/O        | Marker-list filters (`--extract`, `--exclude`, `--chr`) applied in the genotype-factory constructor before any decode                                             | All methods                                                                      | [src/geno_factory/](../../src/geno_factory/)            |
| I/O        | Shared single-pass genotype scan over `M_ref` markers for all P indicators (WtCoxG) or P × K indicators (LEAF)                                                    | WtCoxG, LEAF                                                                     | `runWtCoxG`, `runLEAF`                                  |
| I/O        | Per-thread `GenoCursor` enables independent streaming without lock contention                                                                                    | All methods                                                                      | [src/geno_factory/](../../src/geno_factory/)            |

### User-tunable knobs

| Knob                            | Default                              | Effect                                                                                                                                                            | Guidance                                                                                                                                                                |
| ------------------------------- | ------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--threads`                     | 1                                    | Number of *worker* threads in the chunk-level pool; marker-phase wall time scales as O(... / T).  The main thread is *not* counted in this number — it dispatches chunks, drains output, and briefly occupies one additional CPU during load / finalize / synchronization steps. | Set to the number of physical cores; budget `N + 1` logical cores when sizing the job because the main thread runs alongside the N workers.  Pre-marker work-stealing pools (e.g. `min(T, P·K_τ)` in SPAsqr) reuse the same budget and saturate when T exceeds the unit count. |
| `--chunk-size`                  | 8192 (minimum 256)                   | Markers per work-stealing unit; controls granularity of the writer pipeline and the size of in-flight output buffers.                                              | Reduce on memory-constrained hosts (peak output buffer ≈ C × K × chunkSize × 256 B).  Increase on very large genotype files to amortize per-chunk overhead.             |
| `--spa-z-threshold`             | 2.0                                  | Threshold on `\|Z\|` above which the SPA tail (O(N) or O(E)) replaces the normal approximation.                                                                    | Lowering increases tail accuracy at the cost of more saddlepoint solves; raising trades tail accuracy for speed when only screening is required.                       |
| `--outlier-iqr-multiplier`      | 1.5                                  | IQR multiplier used for residual outlier detection in SPAGRM / SAGELD / SPAmix / SPAmixPlus / SPAmixLocalPlus / WtCoxG / LEAF / SPAsqr.                           | Lower values flag more outliers and slow the per-family SPA tail; defaults are tuned for inverse-normal-transformed residuals.                                          |
| `--spasqr-outlier-abs-bound`    | 0.55                                 | Absolute residual cutoff used only by SPAsqr.                                                                                                                      | Adjust if quantile-residual scales differ markedly from the inverse-normal regime.                                                                                       |
| `--spagrm-control-outlier`      | off                                  | When on, iteratively shrinks / expands the IQR ratio so that the outlier share stays in (0, 5%].                                                                   | Enable if the residual distribution is highly skewed; leave off otherwise.                                                                                              |
| `--spasqr-solver`               | `qmme`                               | Selects the SQR null-model solver.                                                                                                                                 | Use `qmme` (default) for ill-conditioned designs; `conquer` is provided for backward comparison.                                                                        |
| `--spasqr-mode`                 | `score`                              | Selects the SPA score test (default) or per-marker × per-τ Wald refit.                                                                                              | Use `wald` only on a small `--extract` list when you need per-τ effect-size estimates.                                                                                  |
| `--pheno-transform`             | `int`                                | Inverse-normal / standardize / raw Y transform.                                                                                                                    | When `--pred-list` is supplied the transform must match the scale the LOCO PRS was trained on.                                                                          |
| `--compression {gz, zst}`       | plain text                           | Output compression format.                                                                                                                                          | `zst` is faster at comparable ratios; combine with `--compression-level` to trade CPU vs disk.                                                                          |
| `--compression-level`           | 0 (library default)                  | Compression aggressiveness (gz: 1–9, zst: 1–22).                                                                                                                    | Higher levels are I/O-friendly but CPU-heavy; on fast disks the default level is usually optimal.                                                                       |
| `GRAB_MARCH` (build)            | `-march=native`                      | Target ISA for GRAB sources; third-party retains its own SIMD flags.                                                                                               | Override to `-march=x86-64-v2` to build a portable binary that still benefits from runtime AVX-512 / AVX2 dispatch in the kernels that opt in.                          |

---

## Shared infrastructure

The shared marker-streaming infrastructure is described in
[engines.md](engines.md).  In brief:

- All GWAS methods plug into `markerEngine` (single phenotype),
  `multiPhenoEngine` (K phenotypes, single genotype pass with fused
  GEMM + MissBatch), or `locoEngine` (per-chromosome multi-phenotype).
  SPAmixLocalPlus is the sole exception: it ships its own
  `runUnifiedGWAS` driver.
- Each method derives from `MethodBase` and overrides
  `getResultVec` (for `markerEngine`) and / or `getResultBatch` plus
  the fused-GEMM hooks (`supportsFusedGemm`, `fillUnionResiduals`,
  `fillResidualSums`, `processScoreBatch`) for the multi-phenotype
  path.
- The genotype factory ([src/geno_factory/](../../src/geno_factory/))
  exposes the `GenoMeta` / `GenoCursor` interfaces over PLINK 1 BED,
  PLINK 2 PGEN, VCF/BCF, BGEN, and the admixed `.abed` format.
- The subject pipeline (`SubjectSet` + `SubjectData` in
  [src/io/](../../src/io/)) implements the intersection
  `(geno ∩ GRM) ∩ keep \ remove ∩ pheno-presence ∩ residual-validity`
  and provides the `unionToLocal` maps consumed by `multiPhenoEngine`.

---

## GWAS methods

### SPACox

**Source:** [src/spacox/spacox.cpp](../../src/spacox/spacox.cpp) — `runSPACox`

| Property        | Value                                                  |
| --------------- | ------------------------------------------------------ |
| Input           | `--pheno`, `--resid-name` or `--pheno-name + --regression-model`; `--covar` optional |
| GRM             | Not required                                           |
| Engine          | `multiPhenoEngine` (non-fuseable; uses `MissBatch`)    |
| `MethodBase`    | `SPACoxMethod`, `resultSize() = 4` (`P`, `Z`, `BETA`, `SE`) |

**Pre-marker setup:**

1. Build design matrix `X = [1 | covariates]`.
2. Compute `(X'X)⁻¹ X'` for covariate projection.
3. Build the empirical CGF interpolation table (10,000 grid points over
   `[-100, +100]`).  This enables O(1) cumulant lookups inside the SPA
   root-finding loop.

**Per-marker test** (two stages):

1. Unadjusted score test → normal approximation, or SPA tail when
   `|Z| > spaCutoff`.
2. If `p < pvalCovAdjCut`, recompute with a covariate-adjusted genotype
   via `G − X (X'X)⁻¹ X' G`.

**Performance:**

- Pre-marker: O(N) per phenotype — CGF grid dominates (10K × N
  summations).  Serial; CGF table is rebuilt per phenotype.  Design
  matrices are deduplicated across phenotypes sharing the same
  non-missingness pattern.
- Per-marker: O(N) score.  SPA tail (also O(N)) fires only when
  `|Z| > spaCutoff`.  Covariate adjustment adds O(N · P) for markers
  passing `pvalCovAdjCut`.
- Marker-phase wall time: O(M · N / T).

---

### SPAGRM

**Source:** [src/spagrm/spagrm.cpp](../../src/spagrm/spagrm.cpp) — `runSPAGRM`

| Property        | Value                                                            |
| --------------- | ---------------------------------------------------------------- |
| Input           | `--pheno`, `--resid-name` or `--pheno-name + --regression-model`        |
| GRM             | Required (`--sp-grm-grab` or `--sp-grm-plink2`) + `--pairwise-ibd`|
| Engine          | `multiPhenoEngine` (fused GEMM)                                  |
| `MethodBase`    | `SPAGRMMethod`, `resultSize() = 4` (`P`, `Z`, `BETA`, `SE`)      |

**Pre-marker setup:**

1. Load sparse GRM and pairwise-IBD probabilities.
2. Build `SPAGRMClass` null model:
   - Family classification (singletons, 2-subject, 3+-subject families).
   - IQR-based outlier detection (ratio = `--outlier-iqr-multiplier`).
     `--spagrm-control-outlier` iteratively rebalances the IQR ratio so
     the outlier share remains in (0, 5%].
   - Chow–Liu tree: entropy-based MST over family joint distributions
     with complexity O(F × 3^K_fam).
   - MAF-bin interpolation: 11 hard-coded bins from 0.0001 to 0.5.
     Pre-compute the `3^K` probability matrices per bin.  Within each
     outlier family the bins are computed in parallel via atomic
     work-stealing.

**Per-marker test:** score + family-aware variance from GRM quadratic
forms + interpolation between the two closest MAF bins; SPA tail via
`fastGetRoot` when `|Z| > spaCutoff`.

**Performance:**

- Pre-marker: O(E + F · 3^K_fam) — GRM load O(E), Chow–Liu tree
  O(F · K_fam² · 3^K_fam), MAF-bin tables O(11 · F · 3^K_fam).
  Exponential in K_fam but K_fam ≤ 10 (`3^10 ≈ 59K`).
- Per-marker: O(N) score + O(E) GRM variance + O(1) MAF-bin lookup.
  SPA tail O(E) only when `|Z| > spaCutoff`.
- Marker-phase wall time: O(M · (N + E) / T).

---

### SAGELD

**Source:** [src/spagrm/sageld.cpp](../../src/spagrm/sageld.cpp) — `runSAGELD`

| Property        | Value                                                                                              |
| --------------- | -------------------------------------------------------------------------------------------------- |
| Input modes     | Residual mode (`--resid-name R_G, R_<E1>, R_Gx<E1>, ...`) or pheno mode (`--pheno-name + --covar-name + --sageld-x`) |
| GRM             | Required                                                                                           |
| Engine          | `markerEngine` (residual mode) or `multiPhenoEngine` (pheno mode)                                  |
| `MethodBase`    | `SAGELDMethod`, `resultSize() = 4 × (1 + nEnv)`; columns `P_G Z_G BETA_G SE_G P_Gx<E1> Z_Gx<E1> BETA_Gx<E1> SE_Gx<E1> ...` |

**Residual mode:**

Per environment pair `(R_Ek, R_GxEk)`:

1. Compute `λ = (R_GxE · R_G) / (R_G · R_G)`.
2. Combined residual `R_comb = R_GxE − λ · R_G`.
3. Build a `SPAGRMClass` null model for each combined residual.  This
   replicates SPAGRM's per-marker cost K_env times.

**Pheno mode:**

For each `(phenotype, env)` pair, fit `Y ~ X + (E | IID)` by EM-ML and
aggregate BLUP residuals into per-IID `(R_G, R_<E>, R_Gx<E>)` triples
before the marker-level tests run.

**Performance:**

- Pre-marker (residual mode): O(K_env · (E + F · 3^K_fam)) — K_env
  independent `SPAGRMClass` constructions.  Environments are built in
  parallel via atomic work-stealing with `min(T, K_env)` outer workers
  and `max(1, T / K_env)` inner threads per builder.
- Per-marker: O(K_env · (N + E)).
- Marker-phase wall time: O(K_env · M · (N + E) / T).

---

### SPAmix / SPAmixPlus

**Source:** [src/spamix/spamixplus.cpp](../../src/spamix/spamixplus.cpp) — `runSPAmixPlus`

`SPAmix` and `SPAmixPlus` share one runner: the only difference is
whether `--sp-grm-*` is supplied (Plus = with GRM, base = without).

| Property        | Value                                                            |
| --------------- | ---------------------------------------------------------------- |
| Input           | `--pheno + --resid-name + --pc-cols`, or fit path via `--pheno-name + --regression-model` |
| GRM             | Optional (SPAmix) / Required (SPAmixPlus)                        |
| Engine          | `multiPhenoEngine` (fused GEMM)                                  |
| `MethodBase`    | `SPAmixPlusMethod`, `resultSize() = 4` (`P`, `Z`, `BETA`, `SE`)  |

**Pre-marker setup:**

1. Build design matrix `X = [1 | PCs]` for the per-individual AF model.
2. Per phenotype: compute OLS matrices `(X'X)⁻¹`, `(X'X)⁻¹ X'`, standard
   errors for fast per-marker AF fitting.
3. If `--ind-af-coef` is supplied, load the pre-computed AF model
   instead.
4. IQR-based outlier detection per residual column.

**Variance:**

- With GRM (SPAmixPlus): variance from the GRM-based covariance.
- Without GRM (SPAmix): variance =
  `diag(Σ resid² × 2 × AF × (1 − AF))`.

**Performance:**

- Pre-marker: O(N · P²) per unique non-missingness pattern for OLS
  matrices.  Design-matrix inverse, PC matrices, and GRM sub-matrices
  are deduplicated across phenotypes sharing the same valid-subject set.
- Per-marker (with GRM): O(N) score + O(E) GRM variance.
- Per-marker (without GRM): O(N) score + O(N) diagonal variance.
- Marker-phase wall time: O(M · (N + E) / T) with GRM, O(M · N / T)
  without.

---

### SPAmixLocalPlus

**Source:** [src/localplus/spamixlocalp.cpp](../../src/localplus/spamixlocalp.cpp) — `runSPAmixLocalPlus`

| Property        | Value                                                                              |
| --------------- | ---------------------------------------------------------------------------------- |
| Input           | `--admix-bfile + --admix-phi + --pheno + --resid-name` (or fit path)               |
| GRM             | Not directly used (variance comes from the pre-computed phi matrices)              |
| Engine          | Custom `runUnifiedGWAS` (not the standard marker engine)                           |
| Output          | Per-ancestry columns × per-residual output files (one file per residual column)    |

**Pre-marker setup** (per residual column):

1. IQR-based outlier detection.
2. Build `MultiPhenoRprodSoA` (one Structure-of-Arrays layout per
   ancestry) merging phi entries into flat arrays and pre-computing
   `R[i] × R[j] × phi × multiplier` for every phenotype.  The layout is
   cache-friendly for AVX2 / AVX-512 scanning with `PHI_BATCH = 8`.

**Custom engine — `runUnifiedGWAS`:**

Mini-batched processing of 8 markers at a time.  For each marker × each
ancestry the engine computes the score via three fused Eigen GEMMs
(`S_all`, `HR_all`, `HR2_all`), the off-diagonal variance via
`computeVarOffMultiPhenoBatch` (runtime AVX-512 / AVX2 / scalar
dispatch), and the SPA tail when `|Z| > spaCutoff`.

**Output columns** (per residual column,
`PREFIX.PHENO.LocalP[.gz|.zst]`):

```
CHROM  POS  ID  REF  ALT
anc0_MISS_RATE  anc0_ALT_FREQ  anc0_MAC  anc0_P  anc0_BETA  anc0_SE
anc1_MISS_RATE  anc1_ALT_FREQ  anc1_MAC  anc1_P  anc1_BETA  anc1_SE
...
```

**Performance:**

- Pre-marker: O(E · K_anc) per residual column for the `MultiPhenoRprodSoA`
  build (serial).
- Per-marker (amortized over `PHI_BATCH = 8`): O(N · K_anc) score +
  O(E) batch phi variance.
- Marker-phase wall time: O(M · (N · K_anc + E) / T).

---

### WtCoxG

**Source:** [src/wtcoxg/wtcoxg.cpp](../../src/wtcoxg/wtcoxg.cpp) — `runWtCoxG` / `runWtCoxGPheno`

| Property        | Value                                                                          |
| --------------- | ------------------------------------------------------------------------------ |
| Input           | `--pheno + --pheno-name` (binary or `TIME:EVENT` survival) + `--ref-af + --prevalence` |
| GRM             | Optional                                                                       |
| Engine          | `multiPhenoEngine` (fused GEMM)                                                |
| Multi-phenotype | P ≥ 1 via `runWtCoxG`                                                          |
| `MethodBase`    | `WtCoxGMethod`, `resultSize() = 5` (`p_ext`, `p_noext`, `z_ext`, `z_noext`, `p_batch`) |

**Per-phenotype setup** (inside `runWtCoxGPheno`):

1. Load phenotype + covariates; fit null model (logistic or Cox);
   compute residuals / weights / indicator.  Load reference AF file.
2. For every matched marker: per-marker case/ctrl AF statistics, fit
   `P(batch = 1 | G, AF_ref)` logistic, run batch-effect test, estimate
   TPR / σ² (extra variance), and run a 1-D Brent minimization for the
   optimal external weight.

**Multi-phenotype work sharing** (`runWtCoxG`):

| Step                                | Shared | Per-phenotype                  |
| ----------------------------------- | ------ | ------------------------------ |
| Load reference AF, match markers    | yes    | —                              |
| Build `GenoData`                    | yes    | —                              |
| Load sparse GRM                     | yes    | —                              |
| Fit null model (regression)         | —      | `min(T, P)` ∥                  |
| Single geno scan for all indicators | yes    | —                              |
| `testBatchEffects`                  | —      | `min(T, P)` ∥                  |
| Per-marker test                     | yes    | one `multiPhenoEngine` call     |

**Performance:**

- Pre-marker Phase A: O(M_ref + N · C²) — file I/O and null model are
  fast.
- Pre-marker Phase C: O(M_ref · N) shared genotype scan across P
  phenotypes; Brent minimization adds O(M_ref · 20).
- Per-marker: O(N) score + SPA + O(1) reference lookup.
- Marker-phase wall time: O(M · N / T).

---

### LEAF

**Source:** [src/wtcoxg/leaf.cpp](../../src/wtcoxg/leaf.cpp) — `runLEAF` / `runLEAFPheno`

| Property        | Value                                                                          |
| --------------- | ------------------------------------------------------------------------------ |
| Input           | `--pheno + --pheno-name + --pc-cols + --ref-af + --prevalence`                 |
| GRM             | Optional                                                                       |
| Engine          | `multiPhenoEngine` (fused GEMM)                                                |
| Multi-phenotype | P ≥ 1 via `runLEAF`                                                            |
| `MethodBase`    | `LEAFMethod`, `resultSize() = 2 + 3 · K_clust`; columns `meta.p_ext meta.p_noext cl1.p_ext cl1.p_noext cl1.p_batch cl2....` |

**Per-phenotype setup** (inside `runLEAFPheno`):

1. Load phenotype + covariates; fit null model; compute residuals /
   weights / indicator.
2. K-means clustering on PC columns for ancestry assignment.
3. Per cluster: load reference AF, estimate ancestry proportions via
   `summix`, compute ancestry-matched AF, match markers, compute
   per-cluster case/ctrl AF statistics, run the batch-effect test.

**Multi-phenotype work sharing** (`runLEAF`):

| Step                                       | Shared | Per-phenotype                          |
| ------------------------------------------ | ------ | -------------------------------------- |
| K-means clustering (25 restarts)           | yes    | `min(T, 25)` ∥                         |
| Load reference AF files, match markers     | yes    | —                                      |
| Summix estimation (K clusters)             | yes    | `min(T, K)` ∥                          |
| Build `GenoData`                           | yes    | —                                      |
| Load sparse GRM (× K clusters)             | yes    | loaded once per cluster, shared across P|
| Fit null model (× K clusters)              | —      | `min(T, P · K)` ∥ atomic work-stealing |
| Single geno scan for all P × K indicators  | yes    | —                                      |
| Per-cluster `testBatchEffects`             | —      | `min(T, P · K)` ∥                      |
| Per-marker test                            | yes    | one `multiPhenoEngine` call             |

**Per-marker test:** for each cluster, extract the cluster-specific
genotype subvector, run `WtCoxGMethod::computeDual` → (p_ext, p_noext,
z values), then perform CCT (Cauchy-combination) meta-analysis across
clusters.

**Performance:**

- Pre-marker (K-means): O(300 · N · K_clust · P) with deterministic
  per-restart RNG.
- Pre-marker (per-cluster batch effect): O(K_clust · M_ref · N) total.
- Per-marker: O(N) across K_clust clusters (each processes ~N / K_clust
  subjects) + O(K_clust) CCT combination.
- Marker-phase wall time: O(M · N / T).

---

### SPAsqr

**Source:** [src/spasqr/spasqr.cpp](../../src/spasqr/spasqr.cpp) — `runSPAsqr` / `runSPAsqrLoco` / `runSPAsqrWald`

| Property        | Value                                                                          |
| --------------- | ------------------------------------------------------------------------------ |
| Input           | `--pheno + --pheno-name + --spasqr-taus`                                       |
| GRM             | Optional (falls back to identity GRM with a warning when absent)               |
| Engines         | `multiPhenoEngine` (score, non-LOCO) / `locoEngine` (score + `--pred-list`) / standalone Wald loop |
| `MethodBase`    | `SPAsqrMethod`, `resultSize() = 2 · K_τ + 1`; columns `P_CCT P_tau{val}... Z_tau{val}...` (score mode) |

**Per-phenotype setup** (inside `runSPAsqrPheno`):

1. Fit smoothed quantile regression with the selected solver
   (`--spasqr-solver qmme` or `conquer`) at every τ in
   `--spasqr-taus`.  QMME shares one Cholesky decomposition of the
   Hessian upper bound across all τ levels; conquer refits from
   scratch per τ.
2. Per τ: outlier detection (IQR + absolute bound), build
   `SPAGRMClass` null model, pre-compute variance terms.

**Per-marker test (score mode):**

- For each τ, delegate to the corresponding `SPAGRMClass`.
- Collect Z and p from all τ.
- Combine p-values across τ via CCT.
- Output columns: `P_CCT  P_tau{val}...  Z_tau{val}...`.

**Wald mode (`--spasqr-mode wald`):**

For every `(marker, τ)`, refit the joint smoothed QR model with
`[X | G]` by QMME, then derive `β̂_G` and SE from the `(γ, γ)` entry of
the M-estimation sandwich `V = A⁻¹ B A⁻¹ / n`.  No GRM is used.
Output is long-format:

```
CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P TAU BETA SE Z P
```

with one row per `(marker, τ)`.

**LOCO mode (`--pred-list`):**

`runSPAsqrLoco` drives `locoEngine`: each chromosome rebuilds K × K_τ
`SPAGRMClass` instances against the LOCO-adjusted Y.  The conquer /
QMME fits themselves are typically chromosome-independent and are
cloned per chromosome rather than refit.

**Performance:**

- Pre-marker (conquer / QMME QR): O(K_τ · N · P²); parallelized across
  P · K_τ fits with `min(T, P · K_τ)` workers via atomic
  work-stealing; GRM and genotype data shared.
- Pre-marker (null models): O(P · K_τ · E); serial but negligible vs the
  QR fits.
- Per-marker: O(K_τ · (N + E)) + O(K_τ) CCT.
- Marker-phase wall time: O(K_τ · M · (N + E) / T) in score mode.
- LOCO marker-phase wall time: same as above plus per-chromosome
  null-model rebuild (22 × K × K_τ `SPAGRMClass` constructions for an
  autosomal scan).

---

## Utility modes

Utility modes do not use `MethodBase` or the shared marker engines.
Each one has its own streaming loop.

### `--cal-af-coef`

**Source:** [src/spamix/indiv_af.cpp](../../src/spamix/indiv_af.cpp) — `runSPAmixAF`

| Property        | Value                                                                          |
| --------------- | ------------------------------------------------------------------------------ |
| Input           | `--bfile` + `--pc-cols` (provided via `--pheno` or `--covar`)                  |
| GRM             | Not required                                                                   |
| Engine          | `computeAFModelsInMemory` (chunk-level work-stealing)                          |
| Output          | Binary `PREFIX.afc` (one AFModel record per marker)                            |

**Workflow:**

1. Load phenotype / covariate file, extract PC columns.
2. Build design matrix `X = [1 | PCs]`.
3. Multi-threaded per-marker fitting via three-cascade AF model:
   - MAC ≤ 20 → status = 0 (uniform AF, no model).
   - MAC > 20 + OLS good → status = 1 (full OLS coefficients).
   - OLS bad + no significant PCs → status = 0.
   - Otherwise → status = 2 (logistic on significant PCs only).
4. Write binary records: `int8_t status + double betas[1 + nPC]`.

**Performance:**

- Pre-marker: O(N · P²) for design matrix.
- Per-marker: O(N · P) for OLS or logistic fit; I/O bound for large
  genotype files.
- Marker-phase wall time: O(M · N · P / T).

---

### `--cal-pairwise-ibd`

**Source:** [src/spagrm/ibd.cpp](../../src/spagrm/ibd.cpp) — `runPairwiseIBD`

| Property        | Value                                                                          |
| --------------- | ------------------------------------------------------------------------------ |
| Input           | `--bfile/--pfile/--vcf/--bgen + --sp-grm-*`                                    |
| GRM             | Required (off-diagonal pairs)                                                  |
| Engine          | Own chunk-level work-stealing loop                                             |
| Output          | `PREFIX.ibd[.gz|.zst]` — TSV columns `ID1 ID2 pa pb pc`                        |

**Workflow:**

1. Build the filtered subject set via `SubjectSet` (no phenotype loading).
2. Load sparse GRM off-diagonal pairs.
3. Multi-threaded marker streaming:
   - Per chunk: decode genotypes for all subjects.
   - Per pair `(i, j)`: accumulate `X_ij += G_i × G_j × weight` where
     the weight is allele-frequency-dependent.
   - Runtime SIMD dispatch: AVX-512 (8 pairs/cycle), AVX2 (4 pairs/cycle),
     scalar fallback.
4. Post-compute: convert accumulated statistics to IBD sharing
   probabilities `(pa, pb, pc)`.

**Performance:**

- Pre-marker: O(E) for GRM pair index construction.
- Per-marker: O(E_pairs).
- Total wall time: O(M · E_pairs / T).

---

### `--cal-phi`

**Source:** [src/localplus/spamixlocalp.cpp](../../src/localplus/spamixlocalp.cpp) — `runPhiEstimation`

| Property        | Value                                                                          |
| --------------- | ------------------------------------------------------------------------------ |
| Input           | `--admix-bfile + --sp-grm-*`                                                   |
| GRM             | Required (off-diagonal pairs)                                                  |
| Engine          | Own chunk-level work-stealing loop                                             |
| Output          | `PREFIX.phi[.gz|.zst]` — wide TSV (phi matrices per ancestry)                  |

**Workflow:**

1. Build the filtered subject set via `SubjectSet`.
2. Load sparse GRM off-diagonal pairs.
3. Load `AdmixData` (`.abed`: per-ancestry dosage + hapcount).
4. Per ancestry, per off-diagonal pair `(i, j)`:
   - For haplotype scenarios `(h_i ∈ {1, 2}, h_j ∈ {1, 2})`, compute
     allele-frequency-weighted phi values.
   - Store 4 values (A, B, C, D) per pair per ancestry.
5. Write the wide-format TSV.

**Performance:**

- Pre-marker: O(E + M · K_anc).
- Per-marker: O(E_pairs · K_anc).
- Total wall time: O(M · E_pairs · K_anc / T).

---

### `--make-abed`

**Source:** [src/localplus/abed_convert_msp.cpp](../../src/localplus/abed_convert_msp.cpp), [src/localplus/abed_convert_txt.cpp](../../src/localplus/abed_convert_txt.cpp)

Converts either a phased VCF/BCF plus rfmix2 MSP local-ancestry file, or
the extract_tracts text dosage/hapcount output, into the
`.abed + .bim + .fam` triplet consumed by SPAmixLocalPlus and
`--cal-phi`.  Multi-threaded; threading parameter is `--threads`.

---

### `--int-pheno`

**Source:** [src/util/int_pheno.cpp](../../src/util/int_pheno.cpp) — `runIntPheno`

Reads a phenotype file `FID  IID  Y1  Y2  ...` and writes
`PREFIX.txt` with the same columns and row order; each `Y*` column is
independently inverse-normal-transformed (Blom plotting position,
average-rank ties) on its non-missing scope.  Missing entries
(`NA`, `.`, blank) remain missing in the output.

---

## SIMD and runtime dispatch

GRAB uses runtime SIMD dispatch via `__attribute__((target(...)))` and
`__builtin_cpu_supports`.  The resolver runs once at process startup and
caches a `SimdLevel` enum (`Scalar`, `AVX2`, `AVX512`) used by all
dispatch sites.

### Build configuration

GRAB's own sources compile with `GRAB_MARCH = -march=native` by default;
override with `GRAB_MARCH=-march=x86-64-v2` for portable distribution
binaries that still benefit from the runtime AVX-512 / AVX2 paths.
Third-party code (pgenlib, htslib, bgen) retains its own
`-mavx2 -mbmi -mbmi2 -mlzcnt -mfma` flags when the compiler supports
them.

### Dispatch header

[src/util/simd_dispatch.hpp](../../src/util/simd_dispatch.hpp) provides:

```cpp
enum class SimdLevel : int { Scalar = 0, AVX2 = 1, AVX512 = 2 };
inline SimdLevel simdLevel();  // cached, thread-safe
```

### SIMD usage sites

| Component                       | File                                                      | Functions                                  | Tiers                       | Width                |
| ------------------------------- | --------------------------------------------------------- | ------------------------------------------ | --------------------------- | -------------------- |
| SPAmixLocalPlus phi variance    | [localplus/spamixlocalp.cpp](../../src/localplus/spamixlocalp.cpp) | `computeVarOffSoA`, `computeVarOffSoABatch`, `computeVarOffMultiPhenoBatch` | AVX-512 / AVX2 / scalar     | 8 / 4 / 1 entries     |
| Pairwise IBD                    | [spagrm/ibd.cpp](../../src/spagrm/ibd.cpp)                | `ibdAccPairs`                              | AVX-512 / AVX2 / scalar     | 8 / 4 / 1 pairs       |
| ABED decode / encode            | [localplus/abed_io.cpp](../../src/localplus/abed_io.cpp)  | `decodeNoMissAVX2`, `encodeSamplesAVX2`    | AVX2 / scalar               | 32 / 1 samples        |
| LEAF cluster bitmask            | [wtcoxg/leaf.cpp](../../src/wtcoxg/leaf.cpp)              | bitmask operations                         | — (GCC built-ins)           | `popcnt` / `ctz`      |

**Notes:**

- `leaf.cpp` uses `__builtin_popcountll` / `__builtin_ctzll`, which
  GCC and Clang emit unconditionally on x86-64 (`popcnt` is part of the
  `x86-64-v2` baseline).
- The PLINK BED decoder ([src/geno_factory/plink.cpp](../../src/geno_factory/plink.cpp))
  uses plink2's SIMD primitives (`GenoarrCountFreqsUnsafe`,
  `GenoarrLookup16x8bx2`, `GenoarrCountSubsetFreqs2`) from pgenlib.
  These are compiled with the third-party SIMD flags and are not
  governed by GRAB's runtime dispatch.

---

## Method summary

| Method            | Input type                       | GRM      | Engine                | `MethodBase`           | Result columns                                                  | SIMD                       |
| ----------------- | -------------------------------- | -------- | --------------------- | ---------------------- | --------------------------------------------------------------- | -------------------------- |
| SPACox            | 1+ residual columns               | No       | `multiPhenoEngine`    | `SPACoxMethod`         | `P Z BETA SE`                                                   | —                          |
| SPAGRM            | 1+ residual columns               | Yes      | `multiPhenoEngine`    | `SPAGRMMethod`         | `P Z BETA SE`                                                   | —                          |
| SAGELD            | multi-col G×E residual or pheno  | Yes      | `markerEngine` / `multiPhenoEngine` | `SAGELDMethod` | `P_G Z_G BETA_G SE_G P_Gx<E> Z_Gx<E> BETA_Gx<E> SE_Gx<E> ...`   | —                          |
| SPAmix / SPAmixPlus | multi-col residual + PCs        | Opt/Req  | `multiPhenoEngine`    | `SPAmixPlusMethod`     | `P Z BETA SE`                                                   | —                          |
| SPAmixLocalPlus   | multi-col residual + `.abed` + `.phi` | (phi) | `runUnifiedGWAS`     | custom                 | `anc{k}_MISS_RATE anc{k}_ALT_FREQ anc{k}_MAC anc{k}_P anc{k}_BETA anc{k}_SE` | AVX-512 / AVX2 / scalar |
| WtCoxG            | binary / survival + ref AF       | Optional | `multiPhenoEngine`    | `WtCoxGMethod`         | `p_ext p_noext z_ext z_noext p_batch`                           | —                          |
| LEAF              | binary / survival + ref AF       | Optional | `multiPhenoEngine`    | `LEAFMethod`           | `meta.p_ext meta.p_noext cl{k}.p_ext cl{k}.p_noext cl{k}.p_batch` | popcnt / ctz built-ins    |
| SPAsqr (score)    | quantitative + τ list             | Optional | `multiPhenoEngine` or `locoEngine` | `SPAsqrMethod` | `P_CCT P_tau{val}... Z_tau{val}...`                             | —                          |
| SPAsqr (Wald)     | quantitative + τ list             | No       | standalone per-marker | —                      | `TAU BETA SE Z P` (long format)                                 | —                          |
| `--cal-af-coef`   | PCs + geno                        | No       | own                   | —                      | binary `.afc`                                                   | —                          |
| `--cal-pairwise-ibd` | geno + GRM                    | Yes      | own                   | —                      | `ID1 ID2 pa pb pc`                                              | AVX-512 / AVX2 / scalar    |
| `--cal-phi`       | `.abed` + GRM                     | Yes      | own                   | —                      | wide phi TSV                                                    | —                          |
| `--make-abed`     | VCF + MSP or extract_tracts text  | —        | own                   | —                      | `.abed + .bim + .fam`                                           | AVX2 / scalar (decode/encode)|
| `--int-pheno`     | phenotype file                    | —        | serial                | —                      | `PREFIX.txt`                                                    | —                          |
