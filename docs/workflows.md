# Analysis Workflow Reference

This document describes the internal workflows of every GRAB analysis mode
at function level, with notes on performance, threading, and optimization
opportunities.

---

## Complexity & Parallelization Overview

### Notation

| Symbol | Meaning | Typical Range |
|--------|---------|---------------|
| N | Subjects | 100K – 1M |
| M | Markers | 1M – 10M |
| E | GRM non-zero entries | 2–5 × N |
| C | Covariates / PCs | 5 – 20 |
| P | Phenotypes | 1 – 20 |
| J | Ordinal levels (POLMM) | 2 – 10 |
| K_τ | Quantile levels (SPAsqr) | 5 – 9 |
| K_env | Environments (SAGELD) | 1 – 5 |
| K_anc | Ancestries (SPAmixLocalPlus) | 2 – 3 |
| K_clust | Clusters (LEAF) | 3 – 5 |
| F | Families; K_fam = max family size | K_fam ≤ 10 |
| T | Worker threads | 8 – 64 |
| M_ref | Matched reference markers | ~ M |

### Summary Table

All per-marker costs are **per marker**; total marker-phase wall time ≈
Cost × M / T.  Pre-marker stages run once (or once per phenotype) before
marker streaming begins.

| Method | Phase | Dominant Stage | CPU Cost | Parallel | Notes |
|--------|-------|----------------|----------|----------|-------|
| **SPACox** | Pre | CGF interpolation grid | O(N) | Serial | CGF per phenotype (residual-dependent) |
| | Per marker | Score + SPA tail | O(N) | Chunk ∥ | SPA tail for extreme Z only |
| **SPAGRM** | Pre | Chow-Liu tree + MAF-bin tables | O(E + F·3^K_fam) | Serial | K_fam ≤ 10 bounds the exponent |
| | Per marker | Score + GRM variance + SPA | O(N + E) | Chunk ∥ | O(1) MAF-bin interpolation |
| **SAGELD** | Pre | K_env × SPAGRM null model | O(K_env · (E + F·3^K_fam)) | Serial | Envs sequential ᴮ |
| | Per marker | K_env × SPAGRM test | O(K_env · (N + E)) | Chunk ∥ | |
| **SPAmixPlus** | Pre | OLS matrices | O(N·C²) | Serial | (X′X)⁻¹ deduped by NA pattern |
| | Per marker | Score + GRM variance | O(N + E) | Chunk ∥ | O(N) without GRM |
| **SPAmixLocalP** | Pre | RprodSoA build | O(E·K_anc) | Serial | Per residual column |
| | Per marker | Batch score + phi variance | O(N·K_anc + E) | Chunk ∥ + AVX2 | 8-marker mini-batch |
| **POLMM** | Pre (1) | PQL null model (PCG) | O(5000·E·J) | min(T,P) ∥ + T/P probes | Atomic work-stealing; block-diagonal Σ preconditioner |
| | Pre (2) | Variance ratio estimation | O(480·E) | Serial | 120 SNPs × 4 MAC bins |
| | Per marker | Score + ordinal SPA | O(N·C + N·J) | Chunk ∥ | O(1) bin lookup |
| **WtCoxG** | Phase A | Shared loading + ref-AF + GRM | O(M_ref + N·C²) | Serial | `runWtCoxG` multi-pheno entry |
| | Phase B | P null models | O(P·N·C²) | min(T,P) ∥ | Parallel regression |
| | Phase C | 1× geno scan (P indicators) | O(M_ref·N·P) | Serial | Replaces P×`computeMarkerStats` |
| | Phase D | P batch-effect tests | O(P·M_ref·N) | min(T,P) ∥ | GRM shared |
| | Per marker | P × score + bivariate SPA | O(P·N) | Chunk ∥ | Single `multiPhenoEngine` |
| **LEAF** | Phase A | Shared loading + K-means + GRM×K | O(M_ref + 300·N·K·C + K·GRM) | min(T,25) ∥ | K-means restarts parallel |
| | Phase B | P × K null models | O(P·K·N·C²) | min(T,P·K) ∥ | Per (p,c) atomic work-stealing |
| | Phase C | 1× geno scan + K summix | O(M_ref·N·P·K) | min(T,K) ∥ summix | intAF shared; AF_ref synthesis shared |
| | Phase D | P × K batch-effect tests | O(P·K·M_ref·N) | min(T,P·K) ∥ | GRMs preloaded |
| | Per marker | P × (K tests + CCT) | O(P·K·N) | Chunk ∥ | Single `multiPhenoEngine` |
| **SPAsqr** | Pre (1) | P × K_τ × conquer QR | O(P·K_τ·N·C²) | min(T,P·K_τ) ∥ | Atomic work-stealing; GRM/geno shared |
| | Pre (2) | P × K_τ × GRM variance + outlier detect | O(P·K_τ·E) | Serial | No family work; negligible vs Pre (1) |
| | Per marker | K_τ × SPAGRM test + CCT | O(K_τ·(N + E)) | Chunk ∥ | Per-phenotype via multiPhenoEngine |
| **cal-af-coef** | Pre | Design matrix X | O(N·C²) | Serial | One-time |
| | Per marker | OLS / logistic AF model | O(N·C) | Chunk ∥ | I/O bound |
| **cal-pairwise-ibd** | Pre | GRM pair index | O(E) | Serial | One-time |
| | Per marker | Pair accumulation | O(E_pairs) | Chunk ∥ + AVX2 | 4 pairs/cycle |
| **cal-phi** | Pre | GRM + AdmixData load | O(E + M·K_anc) | Serial | |
| | Per marker | Per-ancestry phi accum. | O(E_pairs·K_anc) | Chunk ∥ | |

### Remaining Optimization Candidates

| ID | Issue | Affected Methods | Suggested Fix |
|----|-------|-----------------|---------------|
| ᴮ | SPAGRM null models built sequentially per environment | SAGELD | Parallelize across K_env (independent models) |

---

## Table of Contents

0. [Complexity & Parallelization Overview](#complexity--parallelization-overview)
1. [Shared Infrastructure](#shared-infrastructure)
   - [Marker Engines](#marker-engines)
   - [MethodBase Interface](#methodbase-interface)
   - [GenoCursor Interface](#genocursor-interface)
   - [Genotype Backends](#genotype-backends)
   - [Subject Pipeline](#subject-pipeline)
2. [GWAS Methods](#gwas-methods)
   - [SPACox](#spacox)
   - [SPAGRM](#spagrm)
   - [SAGELD](#sageld)
   - [SPAmixPlus](#spamixplus)
   - [SPAmixLocalPlus](#spamixlocalplus)
   - [POLMM](#polmm)
   - [WtCoxG](#wtcoxg)
   - [LEAF](#leaf)
   - [SPAsqr](#spasqr)
3. [Utility Modes](#utility-modes)
   - [cal-af-coef](#cal-af-coef)
   - [cal-pairwise-ibd](#cal-pairwise-ibd)
   - [cal-phi](#cal-phi)
4. [SIMD and Runtime Dispatch](#simd-and-runtime-dispatch)
5. [Summary Table](#summary-table)

---

## Shared Infrastructure

### Marker Engines

GRAB provides two marker-streaming engines in `src/engine/marker.cpp`.
Both use chunk-level **work-stealing** with ordered output.

#### `markerEngine` — single-phenotype engine

Used by methods that produce a single output file (SAGELD).

Thread model:

| Role           | Count | Description                                       |
| -------------- | ----- | ------------------------------------------------- |
| Worker threads | N     | Pull chunks via `atomic<size_t> nextChunk`         |
| Writer thread  | 1     | Drains completed chunks in strict index order      |

Per-thread resources (allocated once, reused across all chunks):

- Cloned `MethodBase` instance (`method->clone()`)
- Independent `GenoCursor` (genotype decoder)
- `Eigen::VectorXd GVec` (genotype vector, size = nUsed)
- `std::vector<uint32_t> indexForMissing`
- `std::vector<double> rv` (result buffer, size = `resultSize()`)
- `char fmtBuf[32]` (printf formatting buffer)

Synchronization:

- `PaddedFlag` array, 64-byte aligned (`alignas(64)`) to prevent false
  sharing between writer and workers.
- `std::condition_variable writeCv` wakes the writer when a chunk is ready.

Per-marker flow inside a worker:

1. Decode genotypes: `cursor->getGenotypes(gIndex, GVec, altFreq, ...)`
2. QC gate — skip if `missingRate > geno` or `maf < maf` or `mac < mac`
   or `(hwe > 0 && hweP < hwe)`.  Failed markers get `NA` result columns.
3. Impute missing genotypes: `G[j] = 2 * altFreq` for each missing index.
4. Call `method->getResultVec(GVec, altFreq, idx, rv)`.
5. Format 9 meta columns + result columns as a TSV line.

#### `multiPhenoEngine` — multi-phenotype engine

Used by methods that produce one output file **per phenotype**
(SPACox, SPAGRM, SPAmixPlus, POLMM, WtCoxG, LEAF, SPAsqr).

Differences from `markerEngine`:

- Accepts a `std::vector<PhenoTask>` (one per phenotype).
- Decodes genotypes once via the **union** subject mask, then extracts
  per-phenotype genotype subvectors.
- Applies per-phenotype QC independently (allele stats recomputed from
  per-phenotype non-missing sets).
- Opens K output files (one per phenotype); the writer thread interleaves
  them in chunk order.

Performance note: genotype I/O is shared across phenotypes, so multi-pheno
runs are substantially cheaper than K separate single-pheno runs.

### MethodBase Interface

```
class MethodBase
    clone()            → unique_ptr<MethodBase>     deep copy for per-thread use
    resultSize()       → int                        number of output columns
    getHeaderColumns() → string                     tab-separated column names
    prepareChunk(...)                               optional per-chunk hook
    getResultVec(GVec, altFreq, idx, result)         per-marker test
```

Each GWAS method implements a `MethodBase` subclass.  The marker engine
calls `clone()` N times (once per worker thread) so threads never share
mutable state.

### GenoCursor Interface

```
class GenoCursor
    beginSequentialBlock(firstMarker)
    getGenotypes(gIndex, out, altFreq, altCounts, missingRate,
                 hweP, maf, mac, indexForMissing)
```

Per-thread cursors are independent.  Each backend implements this
interface:

### Genotype Backends

| Backend | Flag | Format | Decoder | Notes |
| ------- | ---- | ------ | ------- | ----- |
| PLINK1  | `--bfile` | `.bed` 2-bit | plink2 SIMD primitives | Per-thread `ifstream` seek |
| PGEN    | `--pfile` | `.pgen` | pgenlib (plink2 C++) | Per-thread `PgenReader` |
| VCF/BCF | `--vcf`   | `.vcf.gz` / `.bcf` | htslib | Per-thread `htsFile` handle |
| BGEN    | `--bgen`  | `.bgen v1.2` | genfile library | Probabilistic dosage |

All backends produce a dense `VectorXd` of ALT dosages (0–2) for the
filtered subject set.  Missing genotypes are flagged via
`indexForMissing` and later imputed with `2 * altFreq`.

### Subject Pipeline

See [subject_marker_pipeline.md](subject_marker_pipeline.md) for the
full filtering order.  In brief:

```
final = (genotype ∩ GRM) ∩ (keep \ remove) ∩ (pheno/resid dropna)
```

Methods that need a GRM add a `setGrmSubjects()` step.
Utility modes use `SubjectSet` directly (no phenotype loading).

---

## GWAS Methods

### SPACox

**Source:** `src/spacox/spacox.cpp` — `runSPACox()`

| Property | Value |
| -------- | ----- |
| Input | Single-column `--pheno --resid-name`, optional `--covar` |
| GRM | Not required |
| Engine | `multiPhenoEngine` |
| MethodBase | `SPACoxMethod`, `resultSize() = 2` (P, Z) |

**Pre-marker setup:**

1. Build design matrix `X = [1 | covariates]`.
2. Compute `(X'X)⁻¹X'` for covariate projection.
3. Build **empirical CGF interpolation table**: 10,000 grid points over
   `[-100, +100]`, computing `K₀(t)`, `K₁(t)`, `K₂(t)` plus linear slopes.
   Enables O(1) SPA cumulant lookups.

**Per-marker test (two stages):**

1. Unadjusted score test → normal approximation or SPA tail (when
   `|Z| > spaCutoff`).
2. If p < `pvalCovAdjCut`, recompute with covariate-adjusted genotype
   via `X (X'X)⁻¹X' G` subtraction.

**Performance:**
- Pre-marker: O(N) per phenotype — CGF grid dominates (10K × N summations).
  Currently serial; CGF table is built per phenotype (residual-dependent).
  Design matrices are deduped across phenotypes sharing the same
  non-missingness pattern.
- Per-marker: O(N) score.  SPA tail (also O(N)) fires only when
  `|Z| > spaCutoff`.  Covariate adjustment adds O(N·P) for markers passing
  `pvalCovAdjCut`.
- Marker phase wall time: O(M·N / T).
- CGF table avoids repeated `exp()` calls inside SPA root-finding.

---

### SPAGRM

**Source:** `src/spagrm/spagrm.cpp` — `runSPAGRM()`

| Property | Value |
| -------- | ----- |
| Input | Single-column `--pheno --resid-name` |
| GRM | Required (`--sp-grm-*` + `--pairwise-ibd`) |
| Engine | `multiPhenoEngine` |
| MethodBase | `SPAGRMMethod`, `resultSize() = 2` (P, Z) |

**Pre-marker setup:**

1. Load sparse GRM + pairwise IBD.
2. Build `SPAGRMClass` null model:
   - **Family classification**: singletons, 2-subject, 3+-subject families.
   - **Outlier detection**: IQR-based (ratio = 1.5×).
   - **Chow-Liu tree**: entropy-based MST over family joint distributions.
     Complexity: O(F × 3^K) where F = number of families, K = family size.
   - **MAF-bin interpolation**: 11 hardcoded bins from 0.0001 to 0.5.
     Precompute 3^K probability matrices for each bin.

**Per-marker test:**
- Score + family-aware variance from GRM quadratic forms.
- Interpolate between two closest MAF bins for the observed AF.
- SPA tail via `fastGetRoot()` when `|Z| > spaCutoff`.

**Performance:**
- Pre-marker: O(E + F·3^K_fam) — GRM load O(E), Chow-Liu tree
  O(F·K_fam²·3^K_fam), MAF-bin tables O(11·F·3^K_fam).  Exponential in
  K_fam but K_fam ≤ 10 (3^10 ≈ 59K), so the constant is manageable.
  Entire pre-marker is serial.
- Per-marker: O(N) score + O(E) GRM variance + O(1) MAF-bin lookup.
  SPA tail O(E) only when `|Z| > cutoff`.
- Marker phase wall time: O(M·(N+E) / T).

**Optimization opportunity:** The MAF-bin precomputation is serial; it could
be parallelized across the 11 bins (independent computations).

---

### SAGELD

**Source:** `src/spagrm/sageld.cpp` — `runSAGELD()`

| Property | Value |
| -------- | ----- |
| Input | Multi-column `--pheno --resid-name` (R_G, R_E1, R_GxE1, ...) |
| GRM | Required |
| Engine | `markerEngine` |
| MethodBase | `SAGELDMethod`, `resultSize() = 2 × nEnv` |

**Residual format:**
```
#IID    R_G    R_E1    R_GxE1    R_E2    R_GxE2    ...
```

**Pre-marker setup:**
Per environment pair (R_Ek, R_GxEk):

1. Compute lambda: `λ = (R_GxE · R_G) / (R_G · R_G)`.
2. Combined residual: `R_comb = R_GxE − λ × R_G`.
3. Build a `SPAGRMClass` null model for each combined residual (same cost
   as SPAGRM × number of environments).

**Per-marker test:**
For each environment, compute score via corresponding `SPAGRMClass`.
Output: `P_G, P_GxE1, P_GxE2, ..., Z_G, Z_GxE1, Z_GxE2, ...`.

**Performance:**
- Pre-marker: O(K_env · (E + F·3^K_fam)) — builds K_env independent
  `SPAGRMClass` null models, one per combined residual.  Currently serial
  across environments; could parallelize (independent models).
- Per-marker: O(K_env · (N + E)) — K_env score tests, each with GRM
  variance computation.
- Marker phase wall time: O(K_env · M · (N+E) / T).

---

### SPAmixPlus

**Source:** `src/spamix/spamixplus.cpp` — `runSPAmixPlus()`

| Property | Value |
| -------- | ----- |
| Input | Multi-column `--pheno --resid-name` + `--pc-cols` |
| GRM | Optional |
| Engine | `multiPhenoEngine` |
| MethodBase | `SPAmixPlusMethod`, `resultSize() = 3` (P, Z, BETA) |

**Pre-marker setup:**

1. Build design matrix `X = [1 | PCs]`.
2. Per phenotype:
   - Compute OLS matrices: `(X'X)⁻¹`, `(X'X)⁻¹X'`, standard errors.
   - Purpose: fast per-marker AF model fitting.
3. If `--ind-af-coef` provided: load pre-computed AF models.
4. Outlier detection: IQR-based per residual column.

**Two modes:**
- **With GRM** (SPAmixPlus): variance from GRM-based covariance.
- **Without GRM** (SPAmix): variance = `diag(Σ resid² × 2 × AF × (1−AF))`.

**Performance:**
- Pre-marker: O(N·P²) per unique non-missingness pattern for OLS matrices.
  Design-matrix inverse (X′X)⁻¹, PC matrices, and GRM sub-matrices are
  deduped across phenotypes sharing the same valid-subject set.
- Per-marker (with GRM): O(N) score + O(E) GRM variance.
- Per-marker (without GRM): O(N) score + O(N) diagonal variance.
- Marker phase wall time: O(M·(N+E) / T) with GRM, O(M·N / T) without.

---

### SPAmixLocalPlus

**Source:** `src/localplus/spamixlocalp.cpp` — `runSPAmixLocalPlus()`

| Property | Value |
| -------- | ----- |
| Input | Multi-column `--pheno --resid-name` + `.abed` + `.phi` |
| GRM | Not directly (uses pre-computed phi matrices) |
| Engine | Custom `runUnifiedGWAS` (NOT standard marker engine) |
| Output | Per-ancestry columns × per-residual output files |

**Pre-marker setup:**
Per residual column:

1. Outlier detection (IQR-based).
2. Build `RprodSoA` (Structure-of-Arrays layout):
   - Merge phi entries into flat arrays.
   - Pre-compute: `R[i] × R[j] × phi × multiplier`.
   - Layout is cache-friendly for AVX2 scanning with `PHI_BATCH = 8`.

**Custom engine — `runUnifiedGWAS`:**

Instead of the standard marker engine, this method uses a mini-batched
approach:

1. Decode `PHI_BATCH` (8) markers at once from the `.abed` file.
2. Scan phi entries **once** for all 8 markers (amortizes L3 cache misses
   from random-access phi reads).
3. Per marker × ancestry: compute score via `RprodSoA`.

**SIMD usage:** `computeVarOffSoA()` and `computeVarOffSoABatch()` use
runtime SIMD dispatch with three tiers:
- **AVX-512** (512-bit): 8 phi entries per iteration
  (`computeVarOffSoA`), all 8 batch items in one `__m512d` register
  (`computeVarOffSoABatch`).  Uses `_mm256_cmpeq_epi32_mask` (AVX-512VL)
  for direct mask generation.
- **AVX2** (256-bit): 4 phi entries per iteration, two `__m256d`
  accumulators for the 8-marker batch.
- **Scalar**: portable fallback for pre-Haswell hardware.

Resolution is via `__builtin_cpu_supports()` at process startup.

**Performance:**
- The mini-batch approach reduces random memory access by ~8× compared
  to marker-by-marker processing.
- The SoA phi layout is critical: row-major variant processing would
  thrash cache with 4K × 4K phi matrices.
- Phi pre-computation makes per-marker variance O(E) instead of
  O(E × K_anc).

**Complexity:**
- Pre-marker: O(E·K_anc) per residual column for RprodSoA build (serial).
- Per-marker (amortized over PHI_BATCH=8): O(N·K_anc) score +
  O(E) batch phi variance.
- Marker phase wall time: O(M·(N·K_anc + E) / T).

**SIMD tiers:** The phi variance kernel already dispatches to AVX-512
(8 entries/iteration), AVX2 (4 entries/iteration), or scalar at runtime.
The per-ancestry score computation is not yet vectorized.

---

### POLMM

**Source:** `src/polmm/polmm.cpp` — `runPOLMM()`

| Property | Value |
| -------- | ----- |
| Input | Ordinal `--pheno` (0..J−1) + `--covar` |
| GRM | Required |
| Engine | `multiPhenoEngine` |
| Multi-phenotype | **K ≥ 1** (one PhenoTask per `--pheno-name` column) |
| MethodBase | `POLMMMethod`, `resultSize() = 4` (P, Z, BETA, SE) |

**Multi-phenotype work sharing:**

| Step | Shared / Per-phenotype | Detail |
| ---- | ---------------------- | ------ |
| Load phenotype file | Shared | `loadPhenoFile()` + `finalize()` on union of all subjects |
| Build `GenoData` | Shared | One union mask for all traits |
| Load sparse GRM | Shared | `unionGrm` loaded once in union ordering |
| Build per-phenotype mapping | Per-phenotype | `unionToLocal` / `localToUnion` excluding per-phenotype NaN subjects |
| Re-index GRM | Per-phenotype | `reindexGrm()` extracts sub-GRM via `unionToLocal` |
| Fit null model | Per-phenotype | `fitNullModel()` — PQL iteration (expensive) |
| Compute test matrices | Per-phenotype | `computeTestMatrices()` — RPsiR, RymuVec, projections |
| Estimate variance ratios | Per-phenotype | ~120 SNPs sampled per MAC bin; reads union genotypes, extracts per-phenotype subsets |
| Marker-level genotype I/O | Shared | Decoded once per chunk, extracted K times |

**Pre-marker setup — most expensive of all methods:**

1. **Null model fitting** (`fitNullModel`):
   - Cumulative logit initialization (β, ε cutpoints).
   - Penalized Quasi-Likelihood (PQL) iteration:
     - Inner: update β, random effects b, cutpoints ε (Newton-Raphson).
     - Outer: estimate τ (variance ratio) via REML.
   - **Hutchinson trace estimator**: stochastic trace of (V⁻¹ − ...) using
     PCG solves with random probe vectors.  Dominant cost.

2. **Test matrix precomputation**:
   - Per-subject weights `RPsiR_i`, response residuals `RymuVec_i`.
   - `XXR_Psi_RX`, `XR_Psi_R` matrices.

3. **Variance ratio estimation** (`estimateVarianceRatiosFromUnion`):
   - 4 MAC bins: [20,50), [50,100), [100,500), [500,∞).
   - Sample ~120 SNPs per bin from union genotypes; extract per-phenotype
     subsets via `localToUnion` mapping.
   - Involves repeated PQL-like calculations per sampled SNP.

**Per-marker test:**
- Covariate adjustment: `adjG = G − X(X'WX)⁻¹X'WG`.
- Score: `adjG · RymuVec`.
- MAC-binned variance ratio lookup.
- Custom ordinal SPA (supports J > 2 outcome levels).

**Performance:**
- Pre-marker (null model): O(5000·E·J) per phenotype — ~50 PQL outer
  iterations × ~(30 Hutchinson probes + 1) PCG solves per iteration ×
  ~100 CG iterations × O(E·J) per CG step.  Parallelized across P
  phenotypes with min(T,P) threads; intra-phenotype Hutchinson probes
  also parallelized (see below).
- Pre-marker (variance ratios): O(480·(N·P + E)) — 120 sampled SNPs
  × 4 MAC bins, each requiring PQL-like calculation.  Serial.
- Per-marker: O(N·P) covariate adjustment + O(N·J) ordinal SPA +
  O(1) MAC-bin variance ratio lookup.
- Marker phase wall time: O(M·N·J / T).

**Parallelization (pre-marker):**
- Null model fits run in parallel with `min(T, K)` threads via
  `std::thread` + `std::atomic<int>` work-stealing.
- Block-diagonal Σ preconditioner: Σ = (Ψ⁻¹ + τ·K_ii·I)⁻¹ via
  Sherman-Morrison, reducing CG iterations ~2–3× at O(J) cost per
  subject.
- Intra-phenotype parallelism: when K < T, each worker gets
  `nt = floor(T/K)` threads.  The 30 Hutchinson trace probes (each an
  independent PCG solve — the dominant cost) are dispatched across `nt`
  threads via atomic work-stealing, giving near-linear speedup on ~97%
  of null-model wall time.  Probe vectors are pre-generated sequentially
  for deterministic results regardless of thread count.

---

### WtCoxG

**Source:** `src/wtcoxg/wtcoxg.cpp` — `runWtCoxG()` / `runWtCoxGPheno()`

| Property | Value |
| -------- | ----- |
| Input | Binary or survival `--pheno` + `--ref-af` |
| GRM | Optional |
| Engine | `multiPhenoEngine` |
| Multi-phenotype | **P ≥ 1** via `runWtCoxG` |
| MethodBase | `WtCoxGMethod`, `resultSize() = 5` (p_ext, p_noext, z_ext, z_noext, p_batch) |

**Per-phenotype setup (inside `runWtCoxGPheno`):**

1. **Phase 1**: Load phenotype + covariates, fit null model
   (logistic or Cox regression), compute residuals / weights / indicator.
   Load reference AF file.  Match markers by (CHROM, ID) with allele
   flip detection.

2. **Phase 2** (expensive): For every matched marker:
   - Compute per-marker case/ctrl AF statistics.
   - Fit logistic: `P(batch=1 | G, AF_ref)`.
   - Test batch effect.
   - Estimate TPR, σ² (extra variance).
   - **1-D Brent minimization** for optimal external weight.
   - Nested root-finding for parameter search.
   - Store result in `WtCoxGRefInfo` map.

3. **Phase 3**: refInfoMap ready for lookup during marker streaming.

**Multi-phenotype work sharing (`runWtCoxG`):**

| Step | Shared | Per-phenotype |
| ---- | ------ | ------------- |
| Load ref AF, match markers | Yes | — |
| Build `GenoData` | Yes | — |
| Load sparse GRM | Yes | — |
| Fit null model (regression) | — | min(T,P) ∥ |
| 1× geno scan (all P indicators) | Yes | — |
| `testBatchEffects()` | — | min(T,P) ∥ |
| Per-marker test | Single `multiPhenoEngine` | Genotype decoded once |

**Per-marker test:**
- `prepareChunk()`: populate chunk-level ref info from map.
- With external ref: `wtCoxGTest(..., w_ext, var_ratio_ext, ...)`.
- Without external ref: `wtCoxGTest(..., 0.0, var_ratio_int, ...)`.
- Both paths: score + variance via bivariate normal + SPA.

**Performance:**
- Pre-marker Phase 1: O(M_ref) file I/O + O(N·P²) null model — fast.
- Pre-marker Phase 2: O(M_ref·N) — genotype scan for case/ctrl AF
  statistics dominates; Brent minimization adds O(M_ref·20) (fast
  convergence).  Shared scan across P phenotypes; batch-effect testing
  parallelized with min(T,P) threads.
- Per-marker: O(N) score + SPA + O(1) ref lookup.
- Marker phase wall time: O(M·N / T).

---

### LEAF

**Source:** `src/wtcoxg/leaf.cpp` — `runLEAF()` / `runLEAFPheno()`

| Property | Value |
| -------- | ----- |
| Input | Binary or survival `--pheno` + K `--ref-af` files |
| GRM | Optional |
| Engine | `multiPhenoEngine` |
| Multi-phenotype | **P ≥ 1** via `runLEAF` |
| MethodBase | `LEAFMethod` |

**Per-phenotype setup (inside `runLEAFPheno`):**

1. Load phenotype + covariates; fit null model; compute residuals /
   weights / indicator.
2. K-means clustering on PC columns for ancestry assignment.
3. Per cluster: load ref AF, estimate ancestry proportions via `summix()`,
   compute ancestry-matched AF, match markers, compute per-cluster
   case/ctrl AF stats, run batch-effect testing.

**Multi-phenotype work sharing (`runLEAF`):**

| Step | Shared | Per-phenotype |
| ---- | ------ | ------------- |
| K-means clustering (25 restarts) | Yes — min(T,25) ∥ | — |
| Load ref AF files, match markers | Yes | — |
| Summix estimation (K clusters) | Yes — min(T,K) ∥ | — |
| Build `GenoData` | Yes | — |
| Load sparse GRM (×K clusters) | Yes — loaded once, shared across P | — |
| Fit null model (×K clusters) | — | min(T,P·K) ∥ — atomic work-stealing |
| 1× geno scan (all P×K indicators) | Yes | — |
| Per-cluster `testBatchEffects` | — | min(T,P·K) ∥ |
| Per-marker test | Single `multiPhenoEngine` | Genotype decoded once |

**Per-marker test:**
- For each cluster: extract cluster-specific genotype subvector, run
  `WtCoxGMethod.computeDual()` → p_ext, p_noext, z values.
- **CCT meta-analysis** (Cauchy Combination Test) across clusters.

**Performance:**
- Pre-marker (K-means): O(300·N·K_clust·P) — parallelized across 25
  restarts with min(T,25) workers; deterministic per-restart RNG.
- Pre-marker (per-cluster batch effect): O(K_clust · M_ref · N) total.
  GRM loaded K times (once per cluster), shared across P phenotypes.
  Batch-effect testing parallelized with min(T,P·K) threads.
- Per-marker: O(N) total across K_clust clusters (each processes
  ~N/K_clust subjects) + O(K_clust) CCT combination.
- Marker phase wall time: O(M·N / T).

**Note:** Uses `__builtin_popcountll` and `__builtin_ctzll` for bitmask
operations (GCC/Clang built-ins, unconditional on x86-64).

---

### SPAsqr

**Source:** `src/spasqr/spasqr.cpp` — `runSPAsqr()` / `runSPAsqrPheno()`

| Property | Value |
| -------- | ----- |
| Input | `--pheno` + quantile levels (`--spasqr-taus`) |
| GRM | Required |
| Engine | `multiPhenoEngine` |
| Multi-phenotype | **P ≥ 1** via `runSPAsqr` |
| MethodBase | `SPAsqrMethod` |

**Per-phenotype setup (inside `runSPAsqrPheno`):**

Fits conquer quantile regression to produce residuals (one per tau).
Then per tau: outlier detection, build `SPAGRMClass` null model,
precompute variance terms.

**Multi-phenotype work sharing (`runSPAsqr`):**

| Step | Shared | Per-phenotype |
| ---- | ------ | ------------- |
| Load sparse GRM | Yes | — |
| Build `GenoData` | Yes | — |
| Conquer QR (×K_τ) | — | min(T,P·K_τ) ∥ — atomic work-stealing |
| Build SPAGRMClass (×K_τ) | — | Serial per tau |
| Per-marker test | Single `multiPhenoEngine` | Genotype decoded once |

**Per-marker test:**
- For each tau: delegate to corresponding `SPAGRMClass`.
- Collect Z/P from all taus.
- **CCT meta-analysis**: combine p-values across taus.
- Output: `CCT_P, combined_Z, tau1_Z, tau2_Z, ...`.

**Performance:**
- Pre-marker (conquer QR): O(K_τ · N · P²) — one quantile regression
  per τ level.  Parallelized across P·K_τ fits with min(T,P·K_τ)
  threads via atomic work-stealing; GRM and genotype data shared.
- Pre-marker (null models): O(P·K_τ·E) — builds P × K_τ `SPAGRMClass`
  instances with GRM variance terms only (no family/Chow-Liu work).
  Serial but negligible vs Pre (1).
- Per-marker: O(K_τ · (N + E)) — K_τ SPAGRM tests + O(K_τ) CCT.
- Marker phase wall time: O(K_τ · M · (N+E) / T).

---

## Utility Modes

Utility modes do not use `MethodBase` or the standard marker engines.
They have their own streaming loops.

### cal-af-coef

**Source:** `src/spamix/indiv_af.cpp` — `runSPAmixAF()`

| Property | Value |
| -------- | ----- |
| Input | `--bfile` + `--covar` + `--pc-cols` |
| GRM | Not required |
| Engine | Own engine via `computeAFModelsInMemory()` |
| Output | Binary `.afc` (one AFModel record per marker) |

**Workflow:**

1. Load phenotype/covariate files; extract PC columns.
2. Build design matrix `X = [1 | PCs]`.
3. `computeAFModelsInMemory()` — chunk-level work-stealing, multi-threaded:
   - Three-cascade AF model fitting per marker:
     - **MAC ≤ 20**: status=0 (uniform AF, no model).
     - **MAC > 20 + OLS good**: status=1 (full OLS coefficients).
     - **OLS bad + no significant PCs**: status=0.
     - Otherwise: status=2 (logistic on significant PCs only).
   - Write binary records: `int8_t status + double betas[1+nPC]`.

**Performance:**
- Pre-marker: O(N·P²) for design matrix — one-time.
- Per-marker: O(N·P) for OLS or logistic fit.  I/O bound for large
  genotype files (model fitting is cheap relative to disk reads).
- Marker phase wall time: O(M·N·P / T).  Threading scales well:
  each chunk is independent, no cross-marker state.

---

### cal-pairwise-ibd

**Source:** `src/spagrm/ibd.cpp` — `runPairwiseIBD()`

| Property | Value |
| -------- | ----- |
| Input | `--bfile` + `--sp-grm-*` |
| GRM | Required (off-diagonal pairs) |
| Engine | Own chunk-level work-stealing engine |
| Output | `.ibd` TSV (ID1, ID2, pa, pb, pc) |

**Workflow:**

1. Build filtered subject set via `SubjectSet` (not `SubjectData`).
2. Load sparse GRM off-diagonal pairs.
3. Multi-threaded marker streaming:
   - Per chunk: decode genotypes for all subjects.
   - Per pair (i,j): accumulate `X_ij += G_i × G_j × weight` where
     weight is allele-frequency-dependent.
   - Runtime SIMD dispatch: AVX-512 (8 pairs), AVX2 (4 pairs), or
     scalar per iteration.
4. Post-compute: convert accumulated statistics to IBD sharing
   probabilities (pa, pb, pc).

**Performance:**
- Pre-marker: O(E) for GRM pair index construction — one-time.
- Per-marker: O(E_pairs) pair accumulation, where E_pairs = number of
  off-diagonal GRM pairs.
- Total wall time: O(M·E_pairs / T).  AVX-512 gather processes 8 pairs
  per cycle (~8× scalar), AVX2 processes 4 pairs (~4× scalar).
- Memory: O(E_pairs) accumulators per thread, reduced at end.

---

### cal-phi

**Source:** `src/localplus/spamixlocalp.cpp` — `runPhiEstimation()`

| Property | Value |
| -------- | ----- |
| Input | `--admix-bfile` + `--sp-grm-*` |
| GRM | Required (off-diagonal pairs) |
| Engine | Own chunk-level work-stealing engine |
| Output | `.phi` wide TSV (phi matrices per ancestry) |

**Workflow:**

1. Build filtered subject set via `SubjectSet`.
2. Load sparse GRM off-diagonal pairs.
3. Load AdmixData (.abed: per-ancestry dosage + hapcount).
4. Per ancestry, per off-diagonal pair (i,j):
   - For haplotype scenarios (h_i ∈ {1,2}, h_j ∈ {1,2}):
     compute allele-frequency-weighted phi values.
   - Store 4 values (A, B, C, D) per pair per ancestry.
5. Write wide-format TSV.

**Performance:**
- Pre-marker: O(E + M·K_anc) for GRM pair load + AdmixData.
- Per-marker: O(E_pairs · K_anc) — per-ancestry phi accumulation for
  each off-diagonal pair.
- Total wall time: O(M · E_pairs · K_anc / T).  Threading via
  chunk-level work-stealing over marker dimension.

---

## SIMD and Runtime Dispatch

GRAB uses **runtime SIMD dispatch** via `__attribute__((target(...)))`
and `__builtin_cpu_supports()` for its hand-written SIMD kernels.  The
resolver runs once at process startup and caches a `SimdLevel` enum
(`Scalar`, `AVX2`, `AVX512`) used by all dispatch sites.

### Build Configuration

GRAB's own sources compile with `-march=x86-64-v2` (SSE4.2 + POPCNT
baseline).  Third-party code (pgenlib, htslib, bgen) retains global
`-mavx2 -mbmi -mbmi2 -mlzcnt -mfma` flags when the compiler supports
them.

```makefile
# GRAB sources: runtime dispatch, no compile-time SIMD
GRAB_CXXFLAGS := ... -march=x86-64-v2 ...
# Third-party: compile-time SIMD for best codegen
CXXFLAGS      := ... $(SIMD_FLAGS) ...
```

### Dispatch Header

`src/util/simd_dispatch.hpp` provides:

```cpp
enum class SimdLevel : int { Scalar = 0, AVX2 = 1, AVX512 = 2 };
inline SimdLevel simdLevel();  // cached, thread-safe
```

### SIMD Usage Sites

| Component | File | Functions | Tiers | Width |
| --------- | ---- | --------- | ----- | ----- |
| SPAmixLocalP variance | `localplus/spamixlocalp.cpp` | `computeVarOffSoA`, `computeVarOffSoABatch` | AVX-512 / AVX2 / Scalar | 8 / 4 / 1 entries |
| Pairwise IBD | `spagrm/ibd.cpp` | `ibdAccPairs` | AVX-512 / AVX2 / Scalar | 8 / 4 / 1 pairs |
| ABED decode | `localplus/abed_io.cpp` | `decodeNoMissAVX2`, `encodeSamplesAVX2` | AVX2 / Scalar | 32 / 1 samples |
| LEAF bitmask | `wtcoxg/leaf.cpp` | Cluster bitmask ops | — | GCC built-ins |

**Notes:**
- `leaf.cpp` uses `__builtin_popcountll` / `__builtin_ctzll` which are
  unconditional GCC/Clang built-ins on x86-64 (`popcnt` instruction,
  part of x86-64-v2 baseline).
- The PLINK BED decoder (`geno_factory/plink.cpp`) uses plink2's SIMD
  primitives (`GenoarrCountFreqsUnsafe`, `GenoarrLookup16x8bx2`,
  `GenoarrCountSubsetFreqs2`) from pgenlib — compiled with third-party
  flags, not covered by runtime dispatch.

---

## Summary Table

| Method | Input Type | GRM | Engine | MethodBase | Result Cols | SIMD |
| ------ | ---------- | --- | ------ | ---------- | ----------- | ---- |
| SPACox | 1-col resid | No | `multiPhenoEngine` | SPACoxMethod | 2 | — |
| SPAGRM | 1-col resid | Yes | `multiPhenoEngine` | SPAGRMMethod | 2 | — |
| SAGELD | Multi-col G×E resid | Yes | `markerEngine` | SAGELDMethod | 2×nEnv | — |
| SPAmixPlus | Multi-col resid + PCs | Optional | `multiPhenoEngine` | SPAmixPlusMethod | 3 | — |
| SPAmixLocalPlus | Multi-col resid + .abed + .phi | Via phi | `runUnifiedGWAS` | Custom | 6×K | 512/2/— |
| POLMM | Ordinal pheno | Yes | `multiPhenoEngine` | POLMMMethod | 4 | — |
| WtCoxG | Binary/surv pheno + ref AF | Optional | `multiPhenoEngine` | WtCoxGMethod | 5 | — |
| LEAF | Binary/surv pheno + ref AF | Optional | `multiPhenoEngine` | LEAFMethod | 3+K | Built-ins |
| SPAsqr | Quant pheno + taus | Yes | `multiPhenoEngine` | SPAsqrMethod | 3+ | — |
| cal-af-coef | PCs + geno | No | Own | — | Binary | — |
| cal-pairwise-ibd | geno + GRM | Yes | Own | — | 3 (a,b,c) | 512/2/— |
| cal-phi | .abed + GRM | Yes | Own | — | 4×K | — |
