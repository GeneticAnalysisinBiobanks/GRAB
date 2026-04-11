# Genotype Input and Filtering

This page describes how GRAB loads, filters, and processes genotype data
across all supported formats and analysis methods.

---

## Supported formats

Exactly one of the following is required for GWAS methods and most utility
modes:

| Flag             | Format                        | Subject ID source                                   |
| ---------------- | ----------------------------- | --------------------------------------------------- |
| `--bfile PREFIX` | PLINK 1 (.bed/.bim/.fam)      | `.fam` column 2 (IID)                               |
| `--pfile PREFIX` | PLINK 2 (.pgen/.pvar/.psam)   | `.psam` IID column                                  |
| `--vcf FILE`     | VCF/BCF (.vcf, .vcf.gz, .bcf) | VCF header sample IDs                               |
| `--bgen FILE`    | BGEN v1.2 (.bgen, .sample)    | Embedded sample-ID block → companion `.sample` file |

For `--bgen`: if the BGEN file contains an embedded sample-identifier block,
those IDs are used. Otherwise GRAB looks for a companion `.sample` file (Oxford
format: header `ID_1 ID_2 missing`, type row `0 0 0`, then `FID IID ...` data
lines). If neither is found, numeric IDs (`"0"`, `"1"`, ...) are generated
with a warning.

### Admixed ancestry genotypes (.abed)

SPAmixLocalPlus and `--cal-phi` use `--admix-bfile PREFIX` instead of
standard genotype input. The `.abed` format stores per-ancestry dosage and
hapcount tracks in a BGZF-compressed binary. See
[file_formats.md](file_formats.md#admixed-ancestry-genotypes-abed----admix-bfile-prefix)
for the binary layout.

---

## QC filtering flags

These flags control per-marker quality filtering during GWAS and are shared
across all genotype formats (standard and ABED).

| Flag     | Default   | Description                    |
| -------- | --------- | ------------------------------ |
| `--geno` | 0.1       | Per-marker missing rate cutoff |
| `--maf`  | 1e-5      | Minimum minor allele frequency |
| `--mac`  | 10        | Minimum minor allele count     |

Markers failing **any** cutoff are reported in the output with `NA` for all
result columns (the QC stats themselves are still printed).

### Marker-level filtering flags

| Flag             | Description                                       |
| ---------------- | ------------------------------------------------- |
| `--extract FILE` | Restrict to SNP IDs listed in file (one per line) |
| `--exclude FILE` | Exclude SNP IDs listed in file                    |

These are applied when loading the genotype data, before any QC filtering.

### Subject-level filtering

Subject filtering (`--keep`, `--remove`) is described in
[file_formats.md](file_formats.md#subject-set-intersection).

---

## Genotype processing workflow

The diagram below shows how GRAB processes genotypes from input to
per-marker association results.

### Standard genotype pipeline (all methods except SPAmixLocalPlus)

```
┌──────────────────────────────────────────────────────────────────┐
│  1. Load subject IDs                                             │
│     --bfile → .fam   |  --pfile → .psam                         │
│     --vcf → header   |  --bgen → sample block / .sample         │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
┌──────────────────────────────────────────────────────────────────┐
│  2. Subject intersection                                         │
│     Intersect with --null-resid / --pheno / --covar subjects     │
│     Drop NaN residuals, apply --keep / --remove                  │
│     Build usedMask bitset for fast subject subsetting             │
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
│  4. Parallel marker engine (N worker threads + 1 writer thread)  │
│                                                                  │
│     For each chunk of markers:                                   │
│       For each marker in chunk:                                  │
│                                                                  │
│         a. Decode genotypes for used subjects                    │
│            (2-bit PLINK / dosage from PGEN / BCF / BGEN)         │
│                                                                  │
│         b. Compute QC stats:                                     │
│            • missing rate = #missing / N_used                    │
│            • alt freq = alt_allele_count / (2 × N_non_missing)   │
│            • MAF = min(alt_freq, 1 - alt_freq)                   │
│            • MAC = MAF × 2 × N_non_missing                      │
│            • HWE p-value (exact test)                            │
│                                                                  │
│         c. Apply QC filter:                                      │
│            PASS if missing_rate ≤ --geno                        │
│                 AND MAF ≥ --maf                                  │
│                 AND MAC ≥ --mac                                  │
│            FAIL → output NA for result columns, skip to next     │
│                                                                  │
│         d. Impute missing genotypes with mean (2 × alt_freq)     │
│                                                                  │
│         e. Pass imputed genotype vector to method                │
│            (SPACox, SPAGRM, SPAmix, WtCoxG, POLMM, etc.)        │
│                                                                  │
│         f. Format output: CHROM POS ID REF ALT                   │
│            MISS_RATE ALT_FREQ MAC HWE_P  [method results]        │
│                                                                  │
│     Writer thread writes chunks in order                         │
└────────────────────────────┬─────────────────────────────────────┘
                             ▼
                    Output file written
```

### Multi-phenotype extension

When `--null-resid` contains multiple columns (SPACox, SPAGRM, SPAmix,
SPAmixPlus, SPAmixLocalPlus):

1. Union mask: subjects present in **any** phenotype column
2. Genotypes decoded once per marker using the union mask
3. For each phenotype column:
   - Extract per-phenotype genotype subset via index map
   - Recompute QC stats per phenotype (different N)
   - Apply QC filters independently
   - Run method, write to `PREFIX.PHENO.METHOD[.gz|.zst]`

### ABED pipeline (SPAmixLocalPlus)

```
┌──────────────────────────────────────────────────────────────────┐
│  1. Load subjects from .fam, intersect with --null-resid         │
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
│  3. Load pre-computed phi matrices (--admix-phi)                 │
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
│          • missing rate = #missing / N_used                      │
│          • allele freq = dosage_sum / hapcount_sum                │
│          • MAC = min(dosage_sum, hapcount_sum - dosage_sum)       │
│                                                                  │
│       c. Apply QC filter (per ancestry):                         │
│          PASS if missing_rate ≤ --geno                        │
│                 AND MAF ≥ --maf AND MAF ≤ (1 - --maf)           │
│                 AND MAC ≥ --mac                                  │
│          FAIL → output NA for this ancestry's result columns     │
│                                                                  │
│       d. Compute score statistic S = dosage · R                  │
│          Compute variance using phi-based off-diagonal terms     │
│                                                                  │
│       e. SPA p-value with outlier correction                     │
│                                                                  │
│     Output: CHROM POS ID REF ALT                                 │
│       anc0_MISS_RATE anc0_ALT_FREQ anc0_MAC anc0_P anc0_Z ...   │
└──────────────────────────────────────────────────────────────────┘
```

---

## MAF interval grid (SPAGRM / SAGELD / SPAsqr)

Methods that use the Chow-Liu tree for family-based SPA (SPAGRM, SAGELD,
SPAsqr) build a dynamic MAF grid for probability interpolation. The grid
starts at `min(--maf, --mac / (2 × N_used))` and uses half-decade steps
(×3, ×3.33, ...) up to 0.1, then linear steps to 0.5.

With default `--maf 1e-5` and `--mac 10` for N=100,000 subjects, the
effective minimum is `min(1e-5, 10/200000) = 5e-5`, giving a grid like:

```
5e-5, 1.5e-4, 5e-4, 1.5e-3, 5e-3, 1.5e-2, 5e-2, 0.1, 0.2, 0.3, 0.4, 0.5
```

This provides uniform spacing on a log scale in the rare-variant region,
which improves interpolation accuracy compared to a fixed grid.
