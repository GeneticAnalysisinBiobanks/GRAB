# Engine Comparison

GRAB has four marker-scan engines.  All share the work-stealing
chunk-parallel design in `src/engine/` but differ in how they handle
phenotype multiplexing, chromosome iteration, and residual missingness.

| Property | `markerEngine` | `multiPhenoEngine` | `imputeMultiPhenoEngine` | `locoEngine` |
|---|---|---|---|---|
| Entry point | `marker.cpp` | `marker.cpp` | `marker.cpp` | `loco.cpp` |
| Phenotypes | 1 | K | K | K |
| Output files | 1 | K | K | K |
| Geno decode per marker | 1 | 1 (union) | 1 (union) | 1 (union) |
| Per-phenotype GVec buffer | N/A | 1 per thread | **none** | 1 per thread / none (impute) |
| Per-phenotype stats | same as union | independent | **union (shared)** | independent / union (impute) |
| Per-phenotype QC | same as union | independent | **union (shared)** | independent / union (impute) |
| Residual imputation | N/A | N/A | **0-pad missing** | N/A / 0-pad (impute) |
| Task rebuild | never | never | never | per chromosome |
| Chunk buffer scope | all chunks | all chunks | all chunks | one chromosome |
| Writer architecture | dedicated thread | dedicated thread | dedicated thread | main thread (after join) |
| Used by | SAGELD | SPACox, SPAGRM (drop), SPAmixPlus, POLMM, SPAsqr (drop), WtCoxG, LEAF | SPAGRM (impute, K>1), SPAsqr (impute, K>1) | SPAsqr-LOCO |

### `--pheno-missing` flag

The `--pheno-missing` flag controls how missing residuals are handled in
multi-phenotype runs:

- **`impute`**: Residuals for missing subjects are zero-padded
  to the union sample space.  All K phenotypes share **one union GVec**
  — no per-phenotype genotype extraction or per-phenotype QC.  This
  eliminates K per-phenotype GVec buffers per thread, reducing memory
  from O(T × K × N) to O(T × N) and removing K extraction loops per
  marker.

- **`drop`**: Original behaviour.  Each phenotype uses only its
  non-missing subjects (per-phenotype GVec extraction, per-phenotype
  QC stats).

When only one phenotype/residual is in the analysis, the flag has no
effect and is silently ignored: the single-phenotype path is always used.

### Which methods support `--pheno-missing`?

| Method | Impute engine support | Why not? |
|---|---|---|
| SPAGRM | ✅ Yes | `SPAGRMClass::padToUnionSpace()` expands residuals to union size |
| SPAsqr | ✅ Yes | Pads each inner `SPAGRMClass` and rebuilds `m_residMat` |
| SPAsqr-LOCO | ✅ Yes | Via `locoEngine(…, imputeMode=true)` |
| SPACox | ❌ | Loop bound `m_N`, const ref residuals, pre-computed `m_varResid` |
| SPAmixPlus | ❌ | Outlier indices are per-phenotype, per-subject AF vector sized by `m_N`, GRM variance call uses `m_N` |
| POLMM | ❌ | Internal `m_N` loops, ordinal model state not resizable |
| WtCoxG | ❌ | Multi-phenotype dispatch is per-phenotype entry, internal size assumptions |
| LEAF | ❌ | Same as WtCoxG (shared codebase) |
| SAGELD | N/A | Single-phenotype only (`markerEngine`) |

Passing `--pheno-missing` with an unsupported method is a hard error.
Only SPAGRM and SPAsqr accept this flag.

---

## Shared infrastructure (`marker_impl.hpp`)

All three engines share:

- **`PaddedFlag`** — 64-byte cache-line-aligned completion flag to avoid
  false sharing between workers.
- **`statsFromGVec`** — per-marker genotype stats (altFreq, MAC, missRate,
  HWE p-value).
- **`extractPhenoGVec`** — scatters union-ordered genotypes into
  per-phenotype order via `unionToLocal` mapping.  Used by
  `multiPhenoEngine` and `locoEngine`.
- **Formatting helpers** — `formatLine`, `formatLineNA`, `appendMeta`.
  Zero-allocation in the hot path (stack `char[32]`, append to
  pre-reserved `std::string`).
- **`META_HEADER`** / `buildHeader` — shared TSV column header.

Data structure consumed by all engines:

```
struct PhenoTask {
    string              phenoName;
    unique_ptr<MethodBase> method;
    vector<uint32_t>    unionToLocal;   // UINT32_MAX = absent
    uint32_t            nUsed;
};
```

---

## 1. `markerEngine` — single phenotype

```
void markerEngine(
    const GenoMeta &genoData,
    const MethodBase &method,
    const string &outputFile,
    int nthreads,
    double missingCutoff, double minMafCutoff,
    double minMacCutoff,  double hweCutoff);
```

### Function-level flow

```
markerEngine
  ├── allocate chunkOutput[C], chunkReady[C], nextChunk=0
  ├── launch writerThread
  │     ├── open outputFile
  │     ├── write header
  │     └── for i = 0..C-1:
  │           wait(chunkReady[i])
  │           write chunkOutput[i]       ← sequential ordering
  ├── launch N workerThreads (or inline if N≤1)
  │     each worker:
  │       ├── ThreadContext(method.clone(), gd.makeCursor())
  │       ├── allocate GVec(nSubj), indexForMissing, rv, fmtBuf
  │       └── loop: cidx = nextChunk++   ← work stealing
  │             ├── cursor.beginSequentialBlock
  │             ├── method.prepareChunk
  │             └── for each marker in chunk:
  │                   ├── cursor.getGenotypes   ← 1 decode
  │                   ├── QC check → formatLineNA if fail
  │                   ├── impute missing with 2·altFreq
  │                   ├── method.getResultVec
  │                   └── formatLine → append to out
  │             publish: chunkOutput[cidx] = move(out)
  ├── join workers
  ├── signal stopWriter, join writer
  └── rethrow if error
```

### Per-thread allocation (one-time)

| Object | Size |
|---|---|
| `MethodBase` clone | method-dependent |
| `GenoCursor` | genotype-format-dependent |
| `GVec` | N × 8 bytes |
| `indexForMissing` | reserved N/10 |
| `rv` | reserved 16 doubles |
| `fmtBuf` | 32 bytes (stack) |

### Shared allocation

| Object | Size | Lifetime |
|---|---|---|
| `chunkOutput[C]` | C strings (each ≈ chunkSize × 256 bytes) | freed as writer drains |
| `chunkReady[C]` | C × 64 bytes | entire run |

---

## 2. `multiPhenoEngine` — K phenotypes, one geno pass

```
void multiPhenoEngine(
    const GenoMeta &genoData,
    vector<PhenoTask> &tasks,
    const string &outPrefix, const string &methodName,
    const string &compression, int compressionLevel,
    int nthreads,
    double missingCutoff, double minMafCutoff,
    double minMacCutoff,  double hweCutoff);
```

### Function-level flow

```
multiPhenoEngine
  ├── K = tasks.size()
  ├── build K headers, K NA-suffixes, K output paths
  ├── allocate chunkOutput[C][K], chunkReady[C], nextChunk=0
  ├── launch writerThread
  │     ├── open K writers
  │     ├── write K headers
  │     └── for i = 0..C-1:
  │           wait(chunkReady[i])
  │           write chunkOutput[i][0..K-1]    ← sequential
  ├── launch N workerThreads
  │     each worker:
  │       ├── cursor = makeCursor()
  │       ├── methods[K] = tasks[0..K-1].method.clone()
  │       ├── allocate GVec_union(nUnion), GVec_pheno(maxPhenoN)
  │       ├── allocate unionMissing, phenoMissing, rv, fmtBuf
  │       └── loop: cidx = nextChunk++
  │             ├── cursor.beginSequentialBlock
  │             ├── methods[0..K-1].prepareChunk
  │             └── for each marker:
  │                   ├── cursor.getGenotypes(GVec_union)  ← 1 decode
  │                   └── for p = 0..K-1:
  │                         ├── extractPhenoGVec → GVec_pheno
  │                         ├── statsFromGVec (per-pheno stats)
  │                         ├── per-pheno QC → NA if fail
  │                         ├── impute missing
  │                         ├── methods[p].getResultVec
  │                         └── formatLine → phenoOut[p]
  │             publish: chunkOutput[cidx][0..K-1] = move
  ├── join workers
  ├── signal stopWriter, join writer
  └── rethrow if error
```

### Per-thread allocation (one-time)

| Object | Size |
|---|---|
| K `MethodBase` clones | method-dependent × K |
| `GenoCursor` | genotype-format-dependent |
| `GVec_union` | nUnion × 8 bytes |
| `GVec_pheno` | max(nUsed_p) × 8 bytes |
| `unionMissing`, `phenoMissing` | reserved nUnion/10, maxN/10 |

### Shared allocation

| Object | Size | Lifetime |
|---|---|---|
| `chunkOutput[C][K]` | C × K strings | freed as writer drains |
| `chunkReady[C]` | C × 64 bytes | entire run |

### Overhead vs markerEngine

- **Extra per-thread memory**: K method clones instead of 1;
  `GVec_pheno(maxN)` buffer.
- **Extra per-marker work**: K × `extractPhenoGVec` (O(nUnion) each) +
  K × `statsFromGVec` (O(nPheno) each).
- **Saved**: only 1 genotype decode per marker instead of K.
  For PLINK .bed this saves K−1 bit-unpack + disk-seek passes.

---

## 3. `imputeMultiPhenoEngine` — K phenotypes, shared union GVec (impute mode)

```
void imputeMultiPhenoEngine(
    const GenoMeta &genoData,
    vector<PhenoTask> &tasks,
    const string &outPrefix, const string &methodName,
    const string &compression, int compressionLevel,
    int nthreads,
    double missingCutoff, double minMafCutoff,
    double minMacCutoff,  double hweCutoff);
```

Selected when `--pheno-missing impute` and K > 1 for methods that
support it (SPAGRM, SPAsqr).  Wraps `imputeMultiPhenoEngineRange`.

### Key difference from `multiPhenoEngine`

Before entering the engine loop, each task's method is padded to union
space via `padToUnionSpace(unionToLocal, nUnion)`.  This zero-pads
residuals for missing subjects, so the score statistic naturally cancels
absent subjects.

### Function-level flow

```
imputeMultiPhenoEngine
  ├── for p = 0..K-1:
  │     tasks[p].method.padToUnionSpace(utl, nUnion)   ← expand residuals
  │     tasks[p].nUsed = nUnion
  ├── build K headers, K NA-suffixes, K output paths
  ├── allocate chunkOutput[C][K], chunkReady[C], nextChunk=0
  ├── launch writerThread
  ├── launch N workerThreads
  │     each worker:
  │       ├── cursor = makeCursor()
  │       ├── methods[K] = tasks[0..K-1].method.clone()
  │       ├── allocate GVec_union(nUnion), indexForMissing ← NO GVec_pheno
  │       └── loop: cidx = nextChunk++
  │             ├── cursor.beginSequentialBlock
  │             ├── methods[0..K-1].prepareChunk
  │             └── for each marker:
  │                   ├── cursor.getGenotypesSimple(GVec_union) ← 1 decode
  │                   ├── statsFromGVec(GVec_union)        ← 1 stats call
  │                   ├── QC check (shared for all K)
  │                   ├── impute missing: GVec[j] = 2·altFreq
  │                   └── for p = 0..K-1:
  │                         ├── methods[p].getResultVec(GVec_union)
  │                         └── formatLine → phenoOut[p]
  │             publish: chunkOutput[cidx][0..K-1] = move
  ├── join workers, signal stopWriter, join writer
  └── rethrow if error
```

### Per-thread allocation (one-time)

| Object | Size |
|---|---|
| K `MethodBase` clones | method-dependent × K |
| `GenoCursor` | genotype-format-dependent |
| `GVec_union` | nUnion × 8 bytes |
| `indexForMissing` | reserved nUnion/10 |

No `GVec_pheno` buffer.  No `unionMissing` / `phenoMissing` arrays.

### Memory savings vs `multiPhenoEngine`

| Component | `multiPhenoEngine` | `imputeMultiPhenoEngine` |
|---|---|---|
| GVec buffers per thread | 1 union + 1 pheno | 1 union |
| Stats calls per marker | K | 1 |
| Extract calls per marker | K | 0 |
| QC checks per marker | K | 1 |

For K=100 phenotypes, N=400k, T=16 threads:
- `multiPhenoEngine`: T × (nUnion + maxN) × 8 ≈ 100 MB for GVec alone,
  plus T × K × method clones.
- `imputeMultiPhenoEngine`: T × nUnion × 8 ≈ 50 MB for GVec, plus
  T × K × method clones.  K-fold reduction in stats/extract overhead.

### Residual imputation correctness

The score statistic is `Score = resid · GVec`.  For a subject j that is
missing from phenotype p:

- `padToUnionSpace` sets `resid_padded[j] = 0`
- After imputation: `GVec[j] = 2 · altFreq`
- Contribution: `0 × 2·altFreq = 0` — no effect on the score

The score variance (`R_GRM_R`, `m_varResid`, cumulant tables) is
computed from the original per-phenotype residuals before padding, so
the SPA adjustment uses the correct per-phenotype null distribution.

The union-level QC stats (altFreq, MAC, missRate, HWE) are computed
from all union subjects, not per-phenotype.  This uses a slightly
larger denominator than per-phenotype stats but is acceptable for QC
filtering purposes.

---

## 4. `locoEngine` — per-chromosome LOCO with K phenotypes

```
void locoEngine(
    const GenoMeta &genoData,
    const unordered_set<string> &locoChroms,
    const vector<string> &phenoNames,
    LocoTaskBuilder buildTasks,
    const string &outPrefix, const string &methodName,
    const string &compression, int compressionLevel,
    int nthreads,
    double missingCutoff, double minMafCutoff,
    double minMacCutoff,  double hweCutoff);
```

### LocoTaskBuilder callback

```
using LocoTaskBuilder = function<void(
    const string &chr,
    vector<PhenoTask> &tasks)>;
```

Called once per chromosome.  The callback rebuilds K `PhenoTask`s with
LOCO-adjusted null models.

### Function-level flow

```
locoEngine
  ├── extract chrOrder from genoData.markerInfo()
  ├── intersect with locoChroms → activeChroms
  ├── partition chunks by chromosome → chrChunks[chr] = {start, end}
  ├── buildTasks(activeChroms[0], tasks)     ← get headers
  ├── open K persistent TextWriters, write K headers
  └── for each chr in activeChroms:           ← SERIAL per-chromosome
        ├── look up ChrRange; skip if empty
        ├── buildTasks(chr, tasks)            ← rebuild null models
        ├── rebuild K NA-suffixes
        ├── allocate chunkOutput[nChrChunks][K], chunkReady, nextChunk=0
        ├── launch N workerThreads
        │     each worker:                    ← same body as multiPhenoEngine
        │       ├── cursor, K method clones, GVec_union, GVec_pheno
        │       └── work-stealing loop over chromosome's chunks
        ├── join workers, rethrow if error
        └── drain: for ci = 0..nChrChunks-1:
              write chunkOutput[ci][0..K-1] to persistent writers
  close K writers
```

### Per-chromosome allocation

| Object | Size | Lifetime |
|---|---|---|
| `chunkOutput[nChrChunks][K]` | varies per chr | freed at end of chr loop iteration |
| `chunkReady[nChrChunks]` | nChrChunks × 64 bytes | freed at end of chr loop iteration |
| K `PhenoTask`s | method-dependent × K | rebuilt per chr |
| N × (cursor + K method clones) | per-thread | freed when workers join |

### Overhead vs multiPhenoEngine

- **Per-chromosome task rebuild**: `buildTasks(chr, tasks)` runs
  per-chromosome null-model fitting (e.g., K × ntaus conquer fits for
  SPAsqr-LOCO).  This is the dominant cost difference.
- **Worker re-launch**: threads are created and joined per chromosome
  (22× for autosomal).  Thread creation cost is ~µs per thread,
  negligible vs conquer fitting.
- **No dedicated writer thread**: main thread drains output after workers
  join.  Eliminates pipeline overlap between writing and computing the
  next chromosome, but writing is I/O-bound and fast relative to
  conquer.
- **Lower peak memory**: only one chromosome's output buffers are live at
  a time vs all chunks for multiPhenoEngine.

---

## Overhead summary

| Cost | `markerEngine` | `multiPhenoEngine` | `locoEngine` |
|---|---|---|---|
| Geno decode | M | M | M |
| Method clone | T | T × K | 22 × T × K |
| Cursor create | T | T | 22 × T |
| Stats compute | M | M × K | M × K |
| Extract scatter | — | M × K × O(nUnion) | M × K × O(nUnion) |
| Thread create/join | T + 1 | T + 1 | 22 × T |
| Task rebuild | — | — | 22 × buildTasks |
| Output buffer peak | C × 1 | C × K | max(C_chr) × K |
| Output files | 1 | K | K |

M = total markers, T = threads, K = phenotypes, C = total chunks,
C_chr = chunks in one chromosome.

The dominant cost for all engines is genotype decoding + `getResultVec`
(the SPA test itself).  Engine overhead (thread management, scatter,
stats recomputation) is typically <1% of wall time for genome-wide scans.

For `locoEngine`, the `buildTasks` callback cost (null-model refitting per
chromosome) is significant and method-dependent.  For SPAsqr-LOCO with K
phenotypes and ntaus quantiles, each chromosome requires K × ntaus conquer
fits — making the total null-model cost 22 × K × ntaus fits vs K × ntaus
for the non-LOCO path.
