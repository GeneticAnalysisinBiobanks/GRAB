# Engine Architecture

GRAB ships three reusable marker-scan engines plus one method-specific
engine for SPAmixLocalPlus.  All shared engines use the same
chunk-level, work-stealing thread pool defined in
[src/engine/](../../src/engine/).  They differ only in (i) how many
phenotypes they process per marker, (ii) whether they iterate over
chromosomes, and (iii) how the method-side score statistic is
computed.

| Property              | `markerEngine`         | `multiPhenoEngine` (+ `multiPhenoEngineRange`) | `locoEngine`              | `runUnifiedGWAS`           |
| --------------------- | ---------------------- | ---------------------------------------------- | ------------------------- | -------------------------- |
| Source                | [marker.cpp](../../src/engine/marker.cpp) | [marker.cpp](../../src/engine/marker.cpp) | [loco.cpp](../../src/engine/loco.cpp) | [spamixlocalp.cpp](../../src/localplus/spamixlocalp.cpp) |
| Phenotypes per call   | 1                      | K                                              | K                         | K (custom)                 |
| Output files          | 1                      | K                                              | K                         | K                          |
| Geno decode per marker| 1                      | 1 (union)                                      | 1 (union)                 | per marker (mini-batch)    |
| Score computation     | per-marker `getResultVec` | fused GEMM + MissBatch (see below)          | same as multiPhenoEngine  | fused Eigen GEMMs          |
| Task rebuild          | never                  | never                                          | per chromosome            | never                      |
| Chunk-buffer scope    | all chunks             | all chunks                                     | one chromosome at a time  | all chunks                 |
| Writer architecture   | dedicated writer thread| dedicated writer thread                        | main thread per chromosome| dedicated writer thread    |
| Used by               | SAGELD (residual mode) | SPACox, SPAGRM, SPAmixPlus, WtCoxG, LEAF, SPAsqr (non-LOCO), SAGELD (pheno mode) | SPAsqr LOCO | SPAmixLocalPlus            |

The `--pheno-missing` flag and the dedicated `imputeMultiPhenoEngine`
that previous releases shipped have been removed.  Missing residuals are
now handled by the unified fused-GEMM path described below; no toggle
is required and there is no equivalent flag.

---

## Shared infrastructure (`marker_impl.hpp`)

Every engine reuses the helpers defined in
[src/engine/marker_impl.hpp](../../src/engine/marker_impl.hpp):

- **`PaddedFlag`** — 64-byte cache-line-aligned completion flag, used to
  signal per-chunk readiness without false sharing between workers and
  the writer thread.
- **`statsFromGVec`** / **`statsFromUnionVec`** — per-marker QC
  statistics (`altFreq`, `MAC`, `missRate`, `sumSq`, `hweP`) computed
  from a dense genotype vector with or without an optional 0/1 mask
  column.
- **`extractPhenoGVec`** — scatters a union-ordered genotype vector
  into per-phenotype order via `unionToLocal`.  Used by the
  `MissBatch` path in `multiPhenoEngine` and by `locoEngine`.
- **Formatting helpers** — `formatLine`, `formatLineNA`, `appendMeta`.
  These append directly into a pre-reserved `std::string`, using a
  stack-local `char[32]` scratch buffer; no heap allocations occur in
  the hot path.
- **`META_HEADER`** / **`buildHeader`** — shared TSV column header:
  `CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P`.
- **`makeNaSuffix`** — produces the QC-fail row suffix consisting of
  one `\tNA` per method-specific result column.

The common data structure consumed by `multiPhenoEngine`,
`multiPhenoEngineRange`, and `locoEngine` is `PhenoTask`:

```cpp
struct PhenoTask {
    std::string phenoName;
    std::unique_ptr<MethodBase> method;   // single-phenotype method
    std::vector<uint32_t> unionToLocal;   // UINT32_MAX = absent
    uint32_t nUsed;
};
```

---

## 1. `markerEngine` — single phenotype

```cpp
void markerEngine(
    const GenoMeta &genoData,
    const MethodBase &method,
    const std::string &outputFile,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff);
```

### Function-level flow

```
markerEngine
  ├── allocate chunkOutput[C], chunkReady[C], nextChunk = 0
  ├── launch writerThread
  │     ├── open outputFile
  │     ├── write header
  │     └── for ci = 0..C-1:
  │           wait(chunkReady[ci])
  │           write chunkOutput[ci]            ← sequential ordering
  ├── launch N workerThreads (or inline when N ≤ 1)
  │     each worker:
  │       ├── ThreadContext(method.clone(), gd.makeCursor())
  │       ├── allocate GVec(nUsed), indexForMissing, rv, fmtBuf
  │       └── loop: cidx = nextChunk++          ← work stealing
  │             ├── cursor.beginSequentialBlock
  │             ├── method.prepareChunk
  │             └── for each marker in chunk:
  │                   ├── cursor.getGenotypes   (1 decode)
  │                   ├── QC gate → formatLineNA on fail
  │                   ├── impute missing with 2·altFreq
  │                   ├── method.getResultVec
  │                   └── formatLine → append to out
  │             publish: chunkOutput[cidx] = move(out)
  ├── join workers
  ├── signal stopWriter, join writer
  └── rethrow if a worker threw
```

### Per-thread allocations (one-time)

| Object              | Size                                         |
| ------------------- | -------------------------------------------- |
| `MethodBase` clone  | method-dependent                             |
| `GenoCursor`        | genotype-format-dependent                    |
| `GVec`              | `nUsed × 8` bytes                            |
| `indexForMissing`   | reserved `nUsed / 10`                        |
| `rv`                | reserved 16 doubles                          |
| `fmtBuf`            | 32-byte stack array                          |

### Shared allocations

| Object              | Size                                                | Lifetime              |
| ------------------- | --------------------------------------------------- | --------------------- |
| `chunkOutput[C]`    | C strings (≈ `chunkSize × 256` bytes each)          | freed as writer drains |
| `chunkReady[C]`     | `C × 64` bytes                                      | entire run            |

---

## 2. `multiPhenoEngine` — K phenotypes, single genotype pass

```cpp
void multiPhenoEngine(
    const GenoMeta &genoData,
    std::vector<PhenoTask> &tasks,
    const std::string &outPrefix,
    const std::string &methodName,
    const std::string &compression,
    int compressionLevel,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff);
```

`multiPhenoEngine` opens K writers (one per phenotype), writes their
headers, and delegates the chunk-level loop to
`multiPhenoEngineRange(start, end)`.  `locoEngine` reuses
`multiPhenoEngineRange` directly with a per-chromosome chunk range.

Internally `multiPhenoEngineRange` partitions the K phenotypes into two
groups based on `MethodBase::supportsFusedGemm()`:

- **Fuseable phenotypes** — those whose method overrides
  `supportsFusedGemm` to return `true`.  Currently SPAGRM, SAGELD
  (pheno-input mode), SPAmixPlus, WtCoxG, LEAF, and SPAsqr.  Their
  residuals are gathered into a single shared `AugResid` matrix and the
  score statistic for every fuseable phenotype is computed by **one
  fused GEMM** per marker window.
- **Non-fuseable phenotypes** — currently only SPACox.  Their genotype
  vectors are extracted per unique missingness pattern via the
  `MissBatch` path.

The two groups coexist within one engine call; for mixed runs both
phases execute per window.

### Fused-GEMM path

Each fuseable method contributes its residual columns to an
`AugResid` matrix of shape `N_union × (totalFusedCols + nFuseable)`:

```
AugResid = [ residCols_0 | residCols_1 | ... | residCols_{K-1} | mask_0 | mask_1 | ... | mask_{K-1} ]
```

- The residual columns are filled by `fillUnionResiduals`, which
  scatters per-phenotype residuals into union order and zero-pads rows
  that are absent for that phenotype.
- The trailing `nFuseable` mask columns are 0/1 indicators of presence;
  the engine builds these directly from `unionToLocal`.

For each window of `B = max(preferredBatchSize)` markers, the engine
decodes the genotypes into a `GBatch_union` matrix and runs

```
allScoresAndSums = AugResid^T × GBatch_union
```

in a single Eigen GEMM call.  The result simultaneously yields:

- per-phenotype raw scores `Sₚ = residᵀ × G` (residual columns), and
- per-phenotype genotype sums `gSumₚ = maskᵀ × G` (mask columns), used
  for per-phenotype `altFreq` and `MAC`.

The engine then groups fuseable phenotypes by identical subject sets
(`FusedStatsGroup`) and computes per-marker QC once per group.  Each
method consumes its score slice through `processScoreBatch`, which
receives the pre-computed scores, gSums, gSumSqs, and altFreqs.

### `MissBatch` (non-fuseable) path

Non-fuseable phenotypes are grouped by identical `(nUsed,
unionToLocal)` tuples.  Each unique pattern triggers one extraction
from the union genotype matrix and one QC computation, regardless of
how many phenotypes share that pattern.  Phenotypes within a
`MissBatch` then iterate `getResultBatch` (or `getResultVec` for
methods that have not overridden the batched path) on the extracted
per-phenotype matrix.

### Per-thread allocations

| Object                      | Size                                    |
| --------------------------- | --------------------------------------- |
| K `MethodBase` clones       | method-dependent × K                    |
| `GenoCursor`                | genotype-format-dependent               |
| `GBatch_union`              | `N_union × B × 8` bytes                 |
| `allScoresAndSums`          | `(totalFusedCols + nFuseable) × B`      |
| `passScoresBuf`             | `max(fusedGemmColumns) × B`             |
| `GBatch_pheno_nf`           | `maxNonFusedN × B` (when MissBatch used)|
| `nfBatchMissings[ni]`       | reserved `nUsed / 10` each              |

### Shared allocations

| Object              | Size                                                | Lifetime              |
| ------------------- | --------------------------------------------------- | --------------------- |
| `AugResid`          | `N_union × (totalFusedCols + nFuseable)`            | entire engine call    |
| `chunkOutput[C][K]` | C · K strings                                       | freed as writer drains |
| `chunkReady[C]`     | `C × 64` bytes                                      | entire run            |

### Performance characteristics vs `markerEngine`

- **Saved**: only one genotype decode per marker instead of K, and the
  per-phenotype score computation for fuseable methods is replaced by
  one large GEMM that benefits from Eigen's blocked SIMD kernels.
- **Extra per-thread memory**: K method clones, the `GBatch_union`
  window buffer, and (for non-fuseable phenotypes) a `GBatch_pheno_nf`
  buffer sized by the maximum non-fuseable `nUsed`.
- **Extra per-marker work**: one fused GEMM (cheaper than K separate
  per-phenotype dot products for fuseable methods) plus the
  per-MissBatch extractions for non-fuseable methods (one extraction
  per unique missingness pattern, not per phenotype).

---

## 3. `locoEngine` — per-chromosome LOCO with K phenotypes

```cpp
void locoEngine(
    const GenoMeta &genoData,
    const std::unordered_set<std::string> &locoChroms,
    const std::vector<std::string> &phenoNames,
    LocoTaskBuilder buildTasks,
    const std::string &outPrefix,
    const std::string &methodName,
    const std::string &compression,
    int compressionLevel,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff);
```

### `LocoTaskBuilder` callback

```cpp
using LocoTaskBuilder = std::function<void(
    const std::string &chr,
    std::vector<PhenoTask> &tasks)>;
```

The callback is invoked once per chromosome and is responsible for
rebuilding the K `PhenoTask`s with LOCO-adjusted null models for that
chromosome.

### Function-level flow

```
locoEngine
  ├── extract chrOrder from genoData.markerInfo()
  ├── intersect with locoChroms → activeChroms
  ├── partition chunks by chromosome → chrChunks[chr] = {start, end}
  ├── buildTasks(activeChroms[0], tasks)   ← peek to get headers
  ├── open K persistent TextWriters, write K headers
  └── for each chr in activeChroms:        ← SERIAL per-chromosome
        ├── look up ChrRange; skip if empty (warning logged)
        ├── buildTasks(chr, tasks)         ← rebuild null models
        ├── rebuild per-task naSuffixes
        └── multiPhenoEngineRange(genoData, tasks, naSuffixes,
                                  chrStart, chrEnd, writers, ...)
  close K writers
```

`multiPhenoEngineRange` provides the worker pool, writer thread, fused
GEMM, and MissBatch infrastructure already described.  `locoEngine`'s
role is to drive the per-chromosome rebuild of `tasks` and to keep one
chromosome's results draining before the next chromosome begins.

### Why per-chromosome rebuild?

For SPAsqr LOCO the per-phenotype residual is
`Y_transformed − loco_chr`, which differs per chromosome.  The
`buildTasks` callback refits the null model and rebuilds the
`SPAGRMClass` instances accordingly.  The conquer / QMME quantile-
regression fits themselves are typically chromosome-independent and are
cloned per chromosome rather than refit.

### Overhead vs `multiPhenoEngine`

- **Per-chromosome task rebuild** — `buildTasks(chr, tasks)` runs the
  per-chromosome null-model preparation.  For SPAsqr-LOCO this includes
  one `SPAGRMClass` construction per `(phenotype, τ)` per chromosome,
  which is the dominant additional cost vs the non-LOCO path.
- **Worker re-launch** — threads are created and joined per
  chromosome.  Thread creation cost is microseconds and negligible
  compared with the null-model refit.
- **Writer architecture** — outputs are appended to persistent writers
  in genotype-encounter chromosome order.  The writer thread is shared
  across chromosomes within each `multiPhenoEngineRange` call.
- **Peak memory** — only one chromosome's chunk buffers are live at any
  time, in contrast to `multiPhenoEngine` which keeps all chunks alive
  until the writer drains them.

---

## 4. `runUnifiedGWAS` — SPAmixLocalPlus mini-batched engine

`runUnifiedGWAS` is the method-specific engine for SPAmixLocalPlus
([src/localplus/spamixlocalp.cpp](../../src/localplus/spamixlocalp.cpp)).
It does not use `MethodBase`; instead it implements its own work-stealing
worker pool and writer thread with a structure that resembles
`multiPhenoEngine` but with admixed-ancestry-specific decoding and a
mini-batch of `PHI_BATCH = 8` markers per phi-table scan.

Per marker × per ancestry, the engine:

1. Decodes the dosage and hapcount tracks from `.abed`.
2. Computes per-ancestry QC: `missing rate = #missing / N_used`,
   `allele freq = dosage_sum / hapcount_sum`,
   `MAC = min(dosage_sum, hapcount_sum − dosage_sum)`.  Markers that fail
   `--geno` / `--maf` / `--mac` produce `NA` for that ancestry's
   columns.
3. Computes the per-ancestry per-phenotype score via three fused Eigen
   GEMMs (`S_all = bDosBig^T × R_mat`, `HR_all = bHapBig^T × R_mat`,
   `HR2_all = bHapBig^T × R2_mat`).
4. Scans the phi table once per mini-batch via
   `computeVarOffMultiPhenoBatch`, dispatching at runtime to AVX-512,
   AVX2, or a scalar fallback.
5. Applies the SPA tail with outlier correction when `|Z| > spaCutoff`.

Output (one file per residual column, `PREFIX.PHENO.LocalP[.gz|.zst]`):

```
CHROM  POS  ID  REF  ALT
anc0_MISS_RATE  anc0_ALT_FREQ  anc0_MAC  anc0_P  anc0_BETA  anc0_SE
anc1_MISS_RATE  anc1_ALT_FREQ  anc1_MAC  anc1_P  anc1_BETA  anc1_SE
...
```

---

## Engine selection cheat sheet

| Caller                                         | Engine                          |
| ---------------------------------------------- | ------------------------------- |
| SPACox                                         | `multiPhenoEngine` (MissBatch path; SPACox is the only non-fuseable method today) |
| SPAGRM                                         | `multiPhenoEngine` (fused)      |
| SPAmix / SPAmixPlus                            | `multiPhenoEngine` (fused)      |
| SAGELD residual mode                           | `markerEngine`                  |
| SAGELD pheno-input mode                        | `multiPhenoEngine` (fused)      |
| WtCoxG                                         | `multiPhenoEngine` (fused)      |
| LEAF                                           | `multiPhenoEngine` (fused)      |
| SPAsqr (score mode, no `--pred-list`)          | `multiPhenoEngine` (fused)      |
| SPAsqr (score mode, with `--pred-list`)        | `locoEngine` (fused)            |
| SPAsqr (`--spasqr-mode wald`)                  | per-marker per-τ refit, no SPA, no engine helper |
| SPAmixLocalPlus                                | `runUnifiedGWAS`                |
| `--cal-af-coef` / `--cal-pairwise-ibd` / `--cal-phi` | utility-specific work-stealing loops (see [performance.md](performance.md)) |

---

## Overhead summary

| Cost                | `markerEngine` | `multiPhenoEngine`           | `locoEngine`                    |
| ------------------- | -------------- | ---------------------------- | ------------------------------- |
| Geno decode         | M              | M                            | M                               |
| Method clone        | T              | T × K                        | `nChr × T × K`                  |
| Cursor create       | T              | T                            | `nChr × T`                      |
| QC compute          | M              | M × #FusedStatsGroups + #MissBatch | same as multiPhenoEngine  |
| Fused GEMM          | —              | `M × O(N_union × totalFusedCols × B) / B` | same                |
| MissBatch extract   | —              | `M × #MissBatch × O(N_union)`| same                            |
| Thread create/join  | `T + 1`        | `T + 1`                      | `nChr × T`                      |
| Task rebuild        | —              | —                            | `nChr × buildTasks`             |
| Output buffer peak  | C × 1          | C × K                        | `max(C_chr) × K`                |
| Output files        | 1              | K                            | K                               |

`M` = total markers, `T` = threads, `K` = phenotypes,
`C` = total chunks, `C_chr` = chunks in one chromosome,
`B` = window size (`max(preferredBatchSize)`, currently 16 for fused
methods).

For genome-wide scans the dominant cost in every engine is genotype
decoding plus the per-method work done inside `processScoreBatch` /
`getResultBatch` (the SPA test itself).  Engine overhead (thread
management, MissBatch extraction, QC recomputation) is typically below
1% of wall time.  In `locoEngine`, the `buildTasks` callback cost
becomes significant: for SPAsqr-LOCO with K phenotypes and `K_τ`
quantiles, each chromosome requires K × K_τ `SPAGRMClass`
constructions, making the total null-model preparation cost
`22 × K × K_τ` per autosomal scan versus `K × K_τ` for the non-LOCO
path.
