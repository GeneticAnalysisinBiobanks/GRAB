# Data Loading and Filtering Pipeline

This page describes how GRAB loads input files, filters subjects and
markers, and passes the resulting data to the analysis engines.  The
authoritative implementations live in [src/io/](../../src/io/) and
[src/engine/](../../src/engine/); this page documents user-visible
behavior.

---

## Overview

GRAB executes two sequential pipelines before entering the marker
engine:

1. **Subject pipeline** — determines which subjects participate in the
   analysis.  Runs once at startup.
2. **Marker pipeline** — determines which markers pass QC.  Runs
   per-marker inside the parallel marker engine.

---

## Subject pipeline

The subject pipeline is split across two cooperating classes:

- [SubjectSet](../../src/io/subject_set.cpp) — the lightweight
  genotype → GRM → keep → remove intersection.  Utility modes that have
  no phenotype loading (`--cal-pairwise-ibd`, `--cal-phi`) use this
  class directly.
- [SubjectData](../../src/io/subject_data.cpp) — wraps `SubjectSet` and
  layers the phenotype / residual / covariate intersection on top.
  Every GWAS method uses this class.

### Phase 1: genotype → GRM → keep → remove (`SubjectSet::finalize`)

Starting from the genotype file's full subject list, four filters are
applied in a single pass over `.fam` order:

```
for each subject f in .fam order:
    if GRM provided and f ∉ GRM subjects   → skip
    if --keep provided and f ∉ keep set     → skip
    if --remove provided and f ∈ remove set → skip
    mark f as used
```

The result is a 64-bit bitmask (`usedMask`), where bit `f % 64` in word
`f / 64` is set when `.fam` index `f` passed, together with a dense
index array (`usedFamIndices`).

### Phase 2: phenotype / residual intersection (`SubjectData::finalize`)

After Phase 1, `SubjectData::finalize` narrows the bitmask further:

- For each subject still set in the mask, look up its IID in the
  phenotype / residual file.
- Single-column residual mode: drop subjects with `NaN` in the
  residual column.
- Multi-column residual mode (one residual file with multiple
  columns selected by `--resid-name`): keep subjects with at least one
  non-`NaN` residual.  Per-phenotype masks are built later via
  `buildPerColumnMasks`, which records each subject's status
  independently per column.
- Phenotype-input mode (`--pheno-name`): drop subjects absent from the
  phenotype file.  Per-column NaN handling is delegated to
  `dropNaInColumns` (see below) once the analyzed columns are known.

Dense Eigen arrays for residuals, covariates, and phenotype data are
then built in `.fam` order for the surviving subjects.

For covariates: subjects that pass Phase 1 but are absent from
`--covar` receive per-column mean values.  Covariates therefore do not
participate in the subject intersection.

### Resulting union

```
union = (genotype ∩ GRM) ∩ keep \ remove ∩ pheno-presence ∩ residual-validity
```

In multi-residual mode the union covers all subjects non-`NaN` in at
least one residual column; each individual phenotype's effective subject
set (built by `buildPerColumnMasks`) is a subset of the union.

### GRM subject pre-parsing

Methods that consume a sparse GRM pre-parse subject IDs before loading
the full matrix:

- **GRAB format** (`--sp-grm-grab`): scans the `ID1` / `ID2` columns.
- **plink2 format** (`--sp-grm-plink2`): reads the companion `.grm.id`
  file.  When `.grm.id` is absent, all genotype subjects are assumed to
  belong to the GRM and a warning is logged.

Methods that do not require a GRM (e.g., SPACox) skip the GRM
intersection entirely.

### Post-finalize: `dropNaInColumns`

For methods that consume named phenotype or covariate columns (the
`--pheno-name` fit path; WtCoxG; LEAF; SPAsqr), `dropNaInColumns` is
invoked after `finalize` to drop subjects that have `NaN` in any of the
selected columns.  This rebuilds the bitmask, `nUsed`, and every dense
array.

### Pipeline log

GRAB prints a summary of the subject pipeline:

```
  5000 subjects in --bfile my_data
  4900 subjects in --sp-grm-grab grm.txt
  (no --keep provided)
  (no --remove provided)
  4900 subjects used in analysis

  Union subjects: 4800
  Phenotype 'Y1': 4800 subjects
  Phenotype 'Y2': 4700 subjects
```

---

## Marker pipeline

Markers traverse two filtering stages: list-level filtering at load time
and per-marker QC inside the engine.

### Stage 1: extract / exclude / chr

When loading the genotype file, the genotype factory constructors
(`PlinkData`, `PgenData`, `VcfData`, `BgenData` in
[src/geno_factory/](../../src/geno_factory/)) apply three list filters
before any genotype is decoded:

- `--extract` keeps only the listed SNP IDs.
- `--exclude` drops the listed SNP IDs.
- `--chr` keeps only the listed chromosomes (single values or ranges).

SNP IDs are matched against `.bim` column 2 (PLINK), `.pvar` ID column
(PGEN), the `ID` field (VCF), or the BGEN variant metadata.

### Stage 2: per-marker QC

Inside the marker engine ([src/engine/marker.cpp](../../src/engine/marker.cpp)),
each marker is tested against QC thresholds computed from the active
subject set:

```cpp
bool passQC = !(
    missingRate > missingCutoff   ||   // --geno (default 0.1)
    maf         < minMafCutoff    ||   // --maf  (default 1e-5)
    mac         < minMacCutoff    ||   // --mac  (default 10)
    (hweCutoff > 0 && hweP < hweCutoff) // --hwe (default 0, disabled)
);
```

**QC statistic definitions:**

| Statistic     | Formula                                                |
| ------------- | ------------------------------------------------------ |
| missing rate  | `#missing / N_used`                                    |
| alt freq      | `alt_allele_count / (2 × N_non_missing)`               |
| MAF           | `min(alt_freq, 1 − alt_freq)`                          |
| MAC           | `MAF × 2 × N_non_missing`                              |
| HWE p-value   | exact test (SNPHWE2, Wigginton et al. 2005)            |

- `--hwe 0` (the default) disables HWE filtering.
- Markers that fail QC are emitted with `NA` for every method-specific
  result column; the QC columns themselves (`MISS_RATE`, `ALT_FREQ`,
  `MAC`, `HWE_P`) are still printed.

### Imputation of missing genotypes

Markers that pass QC have their missing entries imputed in place with
the mean dosage `2 × altFreq` before being passed to the method.  The
imputation is performed once per marker by the engine, not by the
method.

### Per-phenotype QC for multi-phenotype runs

`multiPhenoEngine` performs per-phenotype QC.  The implementation
distinguishes two cases (see [engines.md](engines.md) for the full
architecture):

- **Fuseable phenotypes** (those whose method overrides
  `supportsFusedGemm`, currently SPACox, SPAGRM, SPAmix, and SPAsqr):
  QC statistics are computed once per *group* of fuseable phenotypes
  that share the same subject set (`FusedStatsGroup`).  In particular,
  when all fuseable phenotypes share one subject set (the common case
  for residual-input mode), QC is computed once per marker and reused
  across the group.
- **Non-fuseable phenotypes** (WtCoxG, LEAF, SPAmixPlus, SAGELD's
  multi-env G×E): QC is computed per unique missingness pattern via the
  `MissBatch` path.  Phenotypes that share an identical
  `unionToLocal` mapping share one extraction and one QC computation.

In both cases, per-phenotype QC may produce different `MAC` /
`ALT_FREQ` / `MISS_RATE` values from the union-level statistics because
the denominators differ.

---

## Genotype processing workflow

### Standard pipeline (all methods except SPAmixLocalPlus)

```
┌──────────────────────────────────────────────────────────────────┐
│  1. Parse subject IDs                                            │
│     --bfile → .fam   |  --pfile → .psam                          │
│     --vcf → header   |  --bgen → sample block / .sample          │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  2. Subject pipeline                                             │
│     Phase 1 (SubjectSet::finalize):                              │
│       genotype ∩ GRM ∩ keep \ remove → usedMask                  │
│     Phase 2 (SubjectData::finalize):                             │
│       ∩ phenotype presence ∩ residual validity                   │
│       → dense Eigen arrays in .fam order                         │
│     Optional post-finalize: dropNaInColumns                      │
│     Covariates: mean-fill for missing subjects                   │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  3. Optional null-model fit (--pheno-name + --regression-model)  │
│     nullmodel::fitAll runs one regression per PhenoSpec in       │
│     parallel (min(--threads, K_pheno) workers).  --save-resid    │
│     writes PREFIX.null.resid for round-trip reload.              │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  4. Load marker list                                             │
│     Read .bim / .pvar / VCF records / BGEN variant metadata      │
│     Apply --extract / --exclude / --chr (marker list filters)    │
│     Partition surviving markers into chunks                      │
│     (--chunk-size, default 8192, minimum 256)                    │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  5. Marker engine (--threads worker threads + 1 writer thread)   │
│                                                                  │
│     Worker threads use an atomic counter for lock-free chunk     │
│     stealing.  Each worker owns a cloned MethodBase and an       │
│     independent GenoCursor.                                      │
│                                                                  │
│     For each marker in a chunk:                                  │
│                                                                  │
│       a. Decode genotypes for the used subjects only             │
│          (2-bit PLINK / dosage from PGEN / BCF / BGEN)           │
│                                                                  │
│       b. Compute QC statistics                                   │
│          (missing rate, alt freq, MAC, HWE p-value)              │
│                                                                  │
│       c. Apply QC filter; failed markers → output NA suffix      │
│                                                                  │
│       d. Impute missing genotypes with 2 × altFreq               │
│                                                                  │
│       e. Pass the imputed genotype vector to the method          │
│          (SPACox, SPAGRM, SPAmix, ..., SPAsqr)                   │
│                                                                  │
│       f. Format the output line                                  │
│                                                                  │
│     Writer thread writes chunks in input order using             │
│     condition-variable synchronisation with cache-line-padded    │
│     ready flags.                                                 │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
                    Output file(s) written
```

**Output column layout** (all standard methods; `META_HEADER` is
defined in [src/engine/marker_impl.hpp](../../src/engine/marker_impl.hpp)):

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  [method-specific columns]
```

### Multi-phenotype extension

When `--resid-name` selects multiple columns, or when WtCoxG / LEAF /
SPAsqr is invoked with multiple `--pheno-name` columns:

1. **Union mask** — subjects present in *any* selected phenotype
   column.
2. Genotypes are decoded once per marker using the union mask.
3. The fused-GEMM path computes
   `Scores = AugResid^T × G` in a single matrix-matrix product, where
   `AugResid` is the union-level residual matrix concatenated across
   every fuseable phenotype.  Each fuseable phenotype reads its score
   slice via `processScoreBatch`.
4. Non-fuseable phenotypes fall back to per-phenotype extraction via
   the `MissBatch` path (one extraction per unique missingness pattern,
   not per phenotype).
5. Each phenotype writes to its own output file
   `PREFIX.PHENO.METHOD[.gz|.zst]`.

The fused-GEMM and `MissBatch` paths are always active for
multi-phenotype runs; there are no user-facing flags to toggle them.

### ABED pipeline (SPAmixLocalPlus)

```
┌──────────────────────────────────────────────────────────────────┐
│  1. Load subjects from .fam, intersect with --pheno              │
│     Apply --keep / --remove                                      │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  2. Load admix data (.abed + .bim)                               │
│     Apply --extract / --exclude / --chr                          │
│     Detect K ancestries from header (nAnc & 0x7F)                │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  3. Load phi matrices (--admix-phi)                              │
│     Four kinship scenarios per ancestry (A/B/C/D)                │
│     Build the multi-phenotype Structure-of-Arrays rprod table    │
│     once per ancestry; entries pre-multiply phi by Rᵢ·Rⱼ across  │
│     all residual columns.                                        │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  4. runUnifiedGWAS — per-marker mini-batch processing            │
│     (batch of PHI_BATCH = 8 markers, K_pheno phenotypes,         │
│      K ancestries)                                               │
│                                                                  │
│     For each marker, for each ancestry k:                        │
│                                                                  │
│       a. Decode dosage and hapcount tracks from .abed            │
│                                                                  │
│       b. Compute per-ancestry QC:                                │
│          missing rate = #missing / N_used                        │
│          allele freq  = dosage_sum / hapcount_sum                │
│          MAC          = min(dosage_sum, hapcount_sum − dosage_sum)│
│                                                                  │
│       c. Apply QC filter (per ancestry):                         │
│          PASS if missing_rate ≤ --geno                           │
│                 AND MAF ≥ --maf AND MAF ≤ (1 − --maf)            │
│                 AND MAC ≥ --mac                                  │
│          FAIL → output NA for this ancestry's result columns     │
│                                                                  │
│       d. Compute the per-ancestry score via three fused Eigen    │
│          GEMMs (S, HR, HR2) across PHI_BATCH × K_pheno           │
│                                                                  │
│       e. Compute the off-diagonal variance once for the batch    │
│          via computeVarOffMultiPhenoBatch (AVX-512 / AVX2 /      │
│          scalar dispatch)                                        │
│                                                                  │
│       f. SPA tail with outlier correction when |Z| > spaCutoff   │
│                                                                  │
│     Output per phenotype (one file per residual column):         │
│       CHROM  POS  ID  REF  ALT                                   │
│       anc0_MISS_RATE  anc0_ALT_FREQ  anc0_MAC                    │
│       anc0_P  anc0_BETA  anc0_SE                                 │
│       anc1_...  (repeated for each ancestry k)                   │
└──────────────────────────────────────────────────────────────────┘
```

---

## Dispatch call sequence

All GWAS methods follow the same dispatch pattern
([src/cli/dispatch.cpp](../../src/cli/dispatch.cpp)):

```
1. Parse CLI arguments (parseArgs).
2. Split comma-separated column-name lists (--pheno-name,
   --covar-name, --resid-name, --pc-cols, --sageld-x, --spasqr-taus).
3. Resolve genotype format → GenoSpec.

4. Construct SubjectData(famIIDs):
   a. loadResidOne / loadPhenoFile   (residual vs phenotype-input mode)
   b. loadCovar                       (when --covar is supplied)
   c. setKeepRemove(keepFile, removeFile)
   d. setGrmSubjects(grmIDs)          (when a sparse GRM is supplied)
   e. finalize()                       (Phase 1 + Phase 2)
   f. dropNaInColumns(...)             (for fit-path / WtCoxG / LEAF /
                                        SPAsqr)

5. Optional null-model fit: nullmodel::fitAll → setResidualsFromFit.
   When --save-resid is set, writeResidualsFile emits PREFIX.null.resid.

6. Load the full sparse GRM (re-indexed to the final subject order).

7. Construct GenoData (PlinkData / PgenData / VcfData / BgenData /
                       AdmixData) and apply --extract / --exclude /
                       --chr; build chunk indices.

8. Initialise the per-method state (CGF tables, Chow-Liu tree,
   conquer / QMME fits, K-means clusters, ...).

9. Invoke markerEngine / multiPhenoEngine / locoEngine
   (or runUnifiedGWAS for SPAmixLocalPlus) with QC thresholds and
   thread count.
```

---

## MAF interpolation grid (SPAGRM / SAGELD / SPAsqr)

The methods that use the Chow–Liu tree for family-based SPA build a
dynamic MAF grid for probability interpolation.  The grid starts at
`min(--maf, --mac / (2 × N_used))` and proceeds in half-decade steps
(×3, ×3.33, ...) up to 0.1, then in linear steps to 0.5.  This grid is
internal; no user-facing flag controls it.
