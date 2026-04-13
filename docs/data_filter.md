# Data Loading and Filtering Pipeline

This page describes how GRAB loads input files, filters subjects and markers,
and passes the resulting data to analysis engines.

---

## Overview

GRAB processes data in two sequential pipelines before entering the marker
engine:

1. **Subject pipeline** — determines which subjects are analysed.
2. **Marker pipeline** — determines which markers pass QC.

The subject pipeline runs once at startup.  The marker pipeline runs
per-marker inside the parallel marker engine.

---

## Subject pipeline

The subject pipeline is implemented across `SubjectSet` (`src/io/subject_set.cpp`)
and `SubjectData` (`src/io/subject_data.cpp`).  It proceeds in two phases.

### Phase 1: genotype → GRM → keep → remove (SubjectSet::finalize)

Starting from the genotype file's full subject list, four filters are applied
in a single pass over `.fam` order:

```
for each subject f in .fam order:
    if GRM provided and f ∉ GRM subjects  → skip
    if --keep provided and f ∉ keep set   → skip
    if --remove provided and f ∈ remove set → skip
    mark f as used
```

This produces a 64-bit bitmask (`usedMask`) where bit `f % 64` in word
`f / 64` is set if `.fam` index `f` passed, plus a dense index array
(`usedFamIndices`).

### Phase 2: phenotype / residual intersection (SubjectData::finalize)

After Phase 1, SubjectData narrows the mask further:

- For each used subject, look up its IID in the phenotype/residual file.
- If the subject is absent from the pheno/resid file → clear its bit.
- Build dense Eigen arrays (residuals, covariates, phenotype data) in
  `.fam` order for the surviving subjects.

For covariates: subjects that passed Phase 1 but are **absent** from `--covar`
receive per-column mean values.  Covariates do not participate in the
intersection.

### Formula

```
final = (genotype ∩ GRM) ∩ keep \ remove ∩ pheno/resid(dropna)
```

### GRM subject pre-parsing

Methods that use a sparse GRM pre-parse subject IDs before loading the full
matrix:

- **GRAB format** (`--sp-grm-grab`): scans the `ID1`/`ID2` columns.
- **plink2 format** (`--sp-grm-plink2`): reads the companion `.grm.id` file.
  If `.grm.id` is missing, all genotype subjects are assumed to be in the
  GRM (a warning is logged).

Methods that do not require a GRM (e.g., SPACox) skip the GRM intersection
step entirely.

### Post-finalize: dropNaInColumns

For methods using `--pheno-name`, `SubjectData::dropNaInColumns()` is called
after `finalize()` to drop subjects that have `NaN` in any of the selected
phenotype columns.  This rebuilds the bitmask and re-filters all dense
arrays.

### Pipeline log

GRAB prints a summary after subject filtering:

```
  5000 subjects in --bfile my_data
  4900 subjects after GRM intersection
  (no --keep provided)
  (no --remove provided)
  4900 subjects after filtering
  Residual1: 4800 valid, 4700 after filtering
  Residual2: 4900 valid, 4700 after filtering
```

---

## Marker pipeline

Markers pass through two stages: list filtering and per-marker QC.

### Stage 1: extract / exclude

When loading the genotype file, `--extract` and `--exclude` filters are
applied to restrict the marker list before any genotypes are decoded.  SNP
IDs are matched against `.bim` column 2 (or `.pvar` ID column).

These filters are applied in the genotype factory constructors
(`PlinkData`, `PgenData`, `VcfData`, `BgenData` in `src/geno_factory/`).

### Stage 2: per-marker QC

Inside the marker engine (`src/engine/marker.cpp`), each marker is tested
against QC thresholds computed from the used subjects:

```cpp
bool passQC = !(
    missingRate > missingCutoff ||    // --geno (default 0.1)
    maf < minMafCutoff ||             // --maf  (default 1e-5)
    mac < minMacCutoff ||             // --mac  (default 10)
    (hweCutoff > 0 && hweP < hweCutoff)  // --hwe (default 0, disabled)
);
```

**QC stat definitions:**

| Stat          | Formula                                                |
| ------------- | ------------------------------------------------------ |
| missing rate  | #missing / N_used                                      |
| alt freq      | alt_allele_count / (2 × N_non_missing)                 |
| MAF           | min(alt_freq, 1 − alt_freq)                            |
| MAC           | MAF × 2 × N_non_missing                               |
| HWE p-value   | exact test (SNPHWE2, Wigginton et al. 2005)            |

- `--hwe 0` (the default) disables HWE filtering.
- Markers that fail QC are written to output with `NA` for all result columns;
  the QC stats themselves are still printed.

---

## Genotype processing workflow

### Standard pipeline (all methods except SPAmixLocalPlus)

```
┌──────────────────────────────────────────────────────────────────┐
│  1. Parse subject IDs                                            │
│     --bfile → .fam   |  --pfile → .psam                         │
│     --vcf → header   |  --bgen → sample block / .sample         │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  2. Subject pipeline                                             │
│     Phase 1: genotype ∩ GRM ∩ keep \ remove → usedMask          │
│     Phase 2: ∩ pheno/resid → dense arrays                        │
│     Covariates: mean-fill for missing subjects                   │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  3. Load marker list                                             │
│     Read .bim / .pvar / VCF records / BGEN variant metadata      │
│     Apply --extract / --exclude marker filters                   │
│     Partition markers into chunks (--chunk-size, default 8192)   │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  4. Marker engine (N worker threads + 1 writer thread)           │
│                                                                  │
│     Worker threads use atomic counter for lock-free chunk         │
│     stealing.  Each worker owns a cloned method instance and     │
│     an independent GenoCursor.                                   │
│                                                                  │
│     For each marker in a chunk:                                  │
│                                                                  │
│       a. Decode genotypes for used subjects only                 │
│          (2-bit PLINK / dosage from PGEN / BCF / BGEN)           │
│                                                                  │
│       b. Compute QC stats (missing rate, alt freq, MAF, MAC,     │
│          HWE p-value)                                            │
│                                                                  │
│       c. Apply QC filter → FAIL: output NA, skip to next         │
│                                                                  │
│       d. Impute missing genotypes with mean (2 × alt_freq)       │
│                                                                  │
│       e. Pass imputed genotype vector to method                  │
│          (SPACox, SPAGRM, SPAmix, WtCoxG, POLMM, etc.)          │
│                                                                  │
│       f. Format output line                                      │
│                                                                  │
│     Writer thread writes chunks in input order using             │
│     condition-variable synchronisation with cache-line-padded    │
│     ready flags.                                                 │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
                    Output file written
```

**Output columns:**

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  [method-specific columns]
```

### Multi-phenotype extension

When `--resid-name` selects multiple columns, or multiple phenotypes are
analysed (POLMM, WtCoxG, LEAF, SPAsqr):

1. **Union mask**: subjects present in **any** phenotype column.
2. Genotypes decoded once per marker using the union mask.
3. For each phenotype column independently:
   - Extract per-phenotype genotype subset via index map.
   - Recompute QC stats (different N per phenotype).
   - Apply QC filters independently.
   - Run method, write to `PREFIX.PHENO.METHOD[.gz|.zst]`.

### ABED pipeline (SPAmixLocalPlus)

```
┌──────────────────────────────────────────────────────────────────┐
│  1. Load subjects from .fam, intersect with --pheno              │
│     Apply --keep / --remove                                      │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  2. Load admix data (.abed + .bim)                               │
│     Apply --extract / --exclude                                  │
│     Detect K ancestries from header                              │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  3. Load phi matrices (--admix-phi)                              │
│     Four kinship scenarios per ancestry (A/B/C/D)                │
│     Build SoA arrays for fast per-marker variance computation    │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  4. Per-marker mini-batch processing (batch of 8 markers)        │
│                                                                  │
│     For each marker, for each ancestry k:                        │
│                                                                  │
│       a. Decode dosage and hapcount tracks from .abed            │
│                                                                  │
│       b. Compute per-ancestry QC:                                │
│          missing rate = #missing / N_used                        │
│          allele freq  = dosage_sum / hapcount_sum                │
│          MAC = min(dosage_sum, hapcount_sum − dosage_sum)        │
│                                                                  │
│       c. Apply QC filter (per ancestry):                         │
│          PASS if missing_rate ≤ --geno                           │
│                 AND MAF ≥ --maf AND MAF ≤ (1 − --maf)           │
│                 AND MAC ≥ --mac                                  │
│          FAIL → output NA for this ancestry's result columns     │
│                                                                  │
│       d. Compute score statistic S = dosage · R                  │
│          Compute variance from phi-based off-diagonal terms      │
│                                                                  │
│       e. SPA p-value with outlier correction                     │
│                                                                  │
│     Output: CHROM POS ID REF ALT                                 │
│       anc0_MISS_RATE anc0_ALT_FREQ anc0_MAC anc0_P anc0_Z ...   │
└──────────────────────────────────────────────────────────────────┘
```

---

## Call sequence in dispatch

All GWAS methods follow the same dispatch pattern
(`src/cli/dispatch.cpp`):

```
1. Parse CLI arguments
2. Split comma-separated column names (--pheno-name, --covar-name, etc.)
3. Resolve genotype format → GenoSpec

4. Build SubjectData:
   a. SubjectData(famIIDs)       ← from genotype file
   b. loadPhenoFile() / loadResidOne()
   c. loadCovar()
   d. setKeepRemove(keepFile, removeFile)
   e. setGrmSubjects(grmIDs)    ← pre-parsed from GRM file
   f. finalize()                ← Phase 1 + Phase 2

5. Load full sparse GRM (reindexed to final subject order)
6. Construct GenoData (PlinkData / PgenData / VcfData / BgenData)
   - Apply --extract / --exclude marker filters
   - Build chunk indices

7. Initialise method (null model, etc.)

8. Call markerEngine(genoData, method, outputFile,
                     nthreads, missingCutoff, minMafCutoff,
                     minMacCutoff, hweCutoff)
```

---

## MAF interval grid (SPAGRM / SAGELD / SPAsqr)

Methods that use the Chow–Liu tree for family-based SPA build a dynamic MAF
grid for probability interpolation.  The grid starts at
`min(--maf, --mac / (2 × N_used))` and uses half-decade steps
(×3, ×3.33, ...) up to 0.1, then linear steps to 0.5.
