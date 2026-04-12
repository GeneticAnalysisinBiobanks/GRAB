# Analysis Workflow Reference

This document describes the internal workflows of every GRAB analysis mode
at function level, with notes on performance, threading, and optimization
opportunities.

---

## Table of Contents

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
4. [AVX2 and SIMD](#avx2-and-simd)
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
- Pre-marker is fast (O(n) for the grid, one-time).
- Per-marker is O(n) score + O(n) SPA tail (only for significant markers).
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
- Pre-marker is expensive: Chow-Liu tree + per-bin precomputation scales
  with family size (exponential in K, but K is typically small ≤ 10).
- Per-marker is efficient: O(n) score, O(pairs) variance, O(1) bin lookup.

**Optimization opportunity:** The MAF-bin precomputation is serial; it could
be parallelized across bins.

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
- Pre-marker cost scales linearly with number of environments.
- Per-marker cost = nEnv × SPAGRM per-marker cost.

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
- Pre-marker OLS is cheap: O(nPC² × n) per phenotype.
- Per-marker: O(n) score + O(pairs) GRM variance (if GRM present).

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

**AVX2 usage:** `computeVarOffSoA()` and `computeVarOffSoABatch()` use
SSE/AVX2 gather intrinsics (`_mm_i32gather_epi32`, mask operations) for
phi scanning, with a scalar fallback when AVX2 is unavailable.

**Performance:**
- The mini-batch approach reduces random memory access by ~8× compared
  to marker-by-marker processing.
- The SoA phi layout is critical: row-major variant processing would
  thrash cache with 4K × 4K phi matrices.
- Phi pre-computation makes per-marker variance O(pairs) instead of
  O(pairs × K_ancestry).

**Optimization opportunity:** The per-ancestry score computation inside
the inner loop could benefit from wider SIMD (AVX-512) on supported
hardware.

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
- Null model fitting can take minutes for large samples with many
  ordinal levels (J > 5).  The PCG solves dominate.
- Variance ratio estimation is also expensive (~120 SNPs × 4 bins ×
  PQL solve each).
- Per-marker is relatively cheap: O(n) score, O(1) bin lookup.

**Optimization opportunities:**
- The K null model fits are independent — the serial loop
  could be parallelized with `min(nthreads, K)` threads.
- PCG convergence could be accelerated with a better preconditioner
  (currently diagonal).

---

### WtCoxG

**Source:** `src/wtcoxg/wtcoxg.cpp` — `runWtCoxGPheno()`

| Property | Value |
| -------- | ----- |
| Input | Binary or survival `--pheno` + `--ref-af` |
| GRM | Optional |
| Engine | `multiPhenoEngine` |
| Multi-phenotype | **K = 1** |
| MethodBase | `WtCoxGMethod`, `resultSize() = 5` (p_ext, p_noext, z_ext, z_noext, p_batch) |

**Pre-marker setup — 3 phases:**

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

**K > 1 extension (not yet implemented):**

| Step | Would be shared | Would be per-phenotype |
| ---- | --------------- | ---------------------- |
| Load ref AF, match markers | Yes | — |
| Build `GenoData` | Yes | — |
| Load sparse GRM | Yes | — |
| Fit null model (regression) | — | Yes |
| `computeMarkerStats()` | — | Yes (depends on indicator per phenotype) |
| `testBatchEffects()` | — | Yes (depends on residuals / weights) |

**Per-marker test:**
- `prepareChunk()`: populate chunk-level ref info from map.
- With external ref: `wtCoxGTest(..., w_ext, var_ratio_ext, ...)`.
- Without external ref: `wtCoxGTest(..., 0.0, var_ratio_int, ...)`.
- Both paths: score + variance via bivariate normal + SPA.

**Performance:**
- Phase 2 cost is proportional to number of matched markers.  Brent
  minimization converges fast (typically < 20 iterations).
- Per-marker is O(n) + O(1) ref lookup.

---

### LEAF

**Source:** `src/wtcoxg/leaf.cpp` — `runLEAFPheno()`

| Property | Value |
| -------- | ----- |
| Input | Binary or survival `--pheno` + K `--ref-af` files |
| GRM | Optional |
| Engine | `multiPhenoEngine` |
| Multi-phenotype | **K = 1** |
| MethodBase | `LEAFMethod` |

**Pre-marker setup:**

1. Load phenotype + covariates; fit null model; compute residuals /
   weights / indicator.
2. K-means clustering (Lloyd's + K-means++ init) on PC columns for
   ancestry assignment.
3. Per cluster:
   - Load reference AF file per reference population.
   - Ancestry proportion estimation via `summix()`.
   - Ancestry-matched AF: `AF_matched = Σ p_k × AF_pop_k`.
   - Match markers against reference AF.
   - Compute per-cluster case/ctrl AF stats.
   - Per-cluster batch-effect testing via `WtCoxGMethod`.

**K > 1 extension (not yet implemented):**

| Step | Would be shared | Would be per-phenotype |
| ---- | --------------- | ---------------------- |
| K-means clustering | Yes (PC-based, trait-agnostic) | — |
| Load ref AF files, match markers | Yes | — |
| Summix estimation | Yes (genotype AF, trait-agnostic) | — |
| Build `GenoData` | Yes | — |
| Fit null model | — | Yes |
| Per-cluster `computeMarkerStats` | — | Yes (depends on indicator) |
| Per-cluster `testBatchEffects` | — | Yes |
| Load sparse GRM | Currently loaded per cluster (redundant) | Could load once, re-index per cluster |

**Per-marker test:**
- For each cluster: extract cluster-specific genotype subvector, run
  `WtCoxGMethod.computeDual()` → p_ext, p_noext, z values.
- **CCT meta-analysis** (Cauchy Combination Test) across clusters.

**Performance:**
- Pre-marker K-means is cheap.  Batch-effect pre-computation per cluster
  has the same cost as a single WtCoxG run per cluster.
- Per-marker: K_clusters × WtCoxG per-marker cost + CCT combination.

**Optimization opportunity:** The sparse GRM is currently loaded
once per cluster.  Could load once and re-index per cluster.

**Note:** Uses `__builtin_popcountll` and `__builtin_ctzll` for bitmask
operations (GCC/Clang built-ins, unconditional on x86-64).

---

### SPAsqr

**Source:** `src/spasqr/spasqr.cpp` — `runSPAsqrPheno()`

| Property | Value |
| -------- | ----- |
| Input | `--pheno` + quantile levels (`--spasqr-taus`) |
| GRM | Required |
| Engine | `multiPhenoEngine` |
| Multi-phenotype | **K = 1** (multiple taus are within one SPAsqrMethod) |
| MethodBase | `SPAsqrMethod` |

**Pre-marker setup:**
`runSPAsqrPheno` fits conquer quantile regression internally to produce
residuals (one per tau).  Then per tau:

1. Outlier detection: IQR-based + absolute bound.
2. Load sparse GRM.
3. Build `SPAGRMClass` null model (same as SPAGRM, per tau).
4. Precompute variance terms.

**K > 1 extension (not yet implemented):**

| Step | Would be shared | Would be per-phenotype |
| ---- | --------------- | ---------------------- |
| Load sparse GRM | Yes (currently inside `buildSPAsqrMethod`) | — |
| Build `GenoData` | Yes | — |
| Conquer QR (per tau) | — | Yes (τ × phenotype fits) |
| Build SPAGRMClass (per tau) | — | Yes |

**Per-marker test:**
- For each tau: delegate to corresponding `SPAGRMClass`.
- Collect Z/P from all taus.
- **CCT meta-analysis**: combine p-values across taus.
- Output: `CCT_P, combined_Z, tau1_Z, tau2_Z, ...`.

**Performance:**
- Pre-marker cost = nTau × SPAGRM null model cost.
- Per-marker cost = nTau × SPAGRM per-marker cost + CCT.
- The conquer quantile regression adds a one-time
  O(n × p²) cost per tau.

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
- I/O bound for large genotype files.  Model fitting is cheap per marker
  (OLS or logistic with few PCs).
- Threading scales well: each chunk is independent.

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
   - AVX2: process 4 pairs/iteration via `_mm256_i32gather_pd` + FMA.
4. Post-compute: convert accumulated statistics to IBD sharing
   probabilities (pa, pb, pc).

**Performance:**
- Hot loop is O(nMarkers × nPairs).  AVX2 gather amortizes the pair
  inner loop by 4×.
- Memory: O(nPairs) accumulators per thread, reduced at end.

**Note:** No scalar fallback for the AVX2 gather path.  The code
compiles with `#ifdef __AVX2__` but has no `#else` branch, so this
mode requires AVX2 hardware.

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
- O(nMarkers × nPairs × K_ancestry) — scales with number of ancestries.
- Threading via chunk-level work-stealing over marker dimension.

---

## AVX2 and SIMD

GRAB optionally compiles with `-mavx2 -mbmi -mbmi2 -mlzcnt -mfma`.
The Makefile probes for AVX2 support at compile time:

```makefile
AVX2_OK = $(shell $(CXX) ... -mavx2 ... && echo 1 || echo 0)
```

If AVX2 is available, these flags are applied **globally** to all
translation units.  There is no runtime dispatch.

### AVX2 Usage Sites

| Component | File | Functions | Fallback |
| --------- | ---- | --------- | -------- |
| Pairwise IBD | `spagrm/ibd.cpp` | Worker inner loop | **Scalar** (`#else`) |
| SPAmixLocalP variance | `localplus/spamixlocalp.cpp` | `computeVarOffSoA`, `computeVarOffSoABatch` | **Scalar** (`#else`) |
| LEAF bitmask | `wtcoxg/leaf.cpp` | Cluster bitmask ops | **None** (GCC built-ins) |

**Fallback status:**
- `ibd.cpp` and `spamixlocalp.cpp` both have proper `#ifdef __AVX2__` /
  `#else` scalar fallbacks.
- `leaf.cpp` uses `__builtin_popcountll` / `__builtin_ctzll` which are
  unconditional GCC/Clang built-ins on x86-64 (no AVX2 needed, uses
  `popcnt` instruction).
- The PLINK BED decoder (`geno_factory/plink.cpp`) uses plink2's SIMD
  primitives (`GenoarrCountFreqsUnsafe`, `GenoarrLookup16x8bx2`,
  `GenoarrCountSubsetFreqs2`) from pgenlib for vectorized counting and
  decoding.  The bitmask scatter path remains scalar.

---

## Summary Table

| Method | Input Type | GRM | Engine | MethodBase | Result Cols | AVX2 |
| ------ | ---------- | --- | ------ | ---------- | ----------- | ---- |
| SPACox | 1-col resid | No | `multiPhenoEngine` | SPACoxMethod | 2 | — |
| SPAGRM | 1-col resid | Yes | `multiPhenoEngine` | SPAGRMMethod | 2 | — |
| SAGELD | Multi-col G×E resid | Yes | `markerEngine` | SAGELDMethod | 2×nEnv | — |
| SPAmixPlus | Multi-col resid + PCs | Optional | `multiPhenoEngine` | SPAmixPlusMethod | 3 | — |
| SPAmixLocalPlus | Multi-col resid + .abed + .phi | Via phi | `runUnifiedGWAS` | Custom | 6×K | Yes |
| POLMM | Ordinal pheno | Yes | `multiPhenoEngine` | POLMMMethod | 4 | — |
| WtCoxG | Binary/surv pheno + ref AF | Optional | `multiPhenoEngine` | WtCoxGMethod | 5 | — |
| LEAF | Binary/surv pheno + ref AF | Optional | `multiPhenoEngine` | LEAFMethod | 3+K | Built-ins |
| SPAsqr | Quant pheno + taus | Yes | `multiPhenoEngine` | SPAsqrMethod | 3+ | — |
| cal-af-coef | PCs + geno | No | Own | — | Binary | — |
| cal-pairwise-ibd | geno + GRM | Yes | Own | — | 3 (a,b,c) | Yes |
| cal-phi | .abed + GRM | Yes | Own | — | 4×K | — |
