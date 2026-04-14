# Data Input Reference

This page describes the input files accepted by GRAB and the CLI flags that
control them.

Input text files are **whitespace-delimited** (tabs and/or spaces).
Lines beginning with `##` are skipped as comments.

---

## Genotype files

Exactly one of the following flags is required for GWAS methods:

| Flag             | Format                         | Subject ID source                                     |
| ---------------- | ------------------------------ | ----------------------------------------------------- |
| `--bfile PREFIX` | PLINK 1 (.bed/.bim/.fam)       | `.fam` column 2 (IID)                                 |
| `--pfile PREFIX` | PLINK 2 (.pgen/.pvar/.psam)    | `.psam` IID column                                    |
| `--vcf FILE`     | VCF/BCF (.vcf, .vcf.gz, .bcf)  | VCF header sample IDs                                 |
| `--bgen FILE`    | BGEN v1.2 (.bgen, .sample)     | Embedded sample-ID block or companion `.sample` file  |

For `--bgen`: if the BGEN file contains an embedded sample-identifier block,
those IDs are used.  Otherwise GRAB looks for a companion `.sample` file
(Oxford format: header `ID_1 ID_2 missing`, type row `0 0 0`, then
`FID IID ...` data lines).  If neither is found, numeric IDs (`"0"`, `"1"`,
...) are generated with a warning.

### Admixed ancestry genotypes — `--admix-bfile PREFIX`

SPAmixLocalPlus and `--cal-phi` use `--admix-bfile` instead of standard
genotype input.  The `.abed` format stores 2K tracks (dosage + hapcount per
ancestry) in a BGZF-compressed binary sharing PLINK `.fam` and `.bim` files.

Header layout (8 bytes):

| Offset | Size | Field       | Value                                        |
| ------ | ---- | ----------- | -------------------------------------------- |
| 0      | 2    | magic       | `0xAD 0x4D`                                  |
| 2      | 1    | version     | `0x02`                                       |
| 3      | 1    | nAnc        | K in bits 0–6; bit 7 = NO_MISSING flag       |
| 4      | 4    | nSamples    | N (little-endian)                             |

Each marker stores 2K tracks of `ceil(N/4)` bytes using PLINK 2-bit encoding
(`00`→0, `10`→1, `11`→2, `01`→missing).

The `.abed` file is produced by `grab --make-abed` from one of:

- `--vcf FILE --rfmix-msp FILE` (rfmix2 MSP local-ancestry output)
- `--admix-text-prefix PREFIX` (extract_tracts text dosage/hapcount files)

---

## Phenotype file — `--pheno FILE`

Whitespace-delimited with a mandatory header line.  One row per subject.
All data column names must match `[0-9A-Za-z_\-.]+`.

The subject ID column is determined from the header:

| Header columns 1–2        | Subject ID   | Data columns start at |
| ------------------------- | ------------ | --------------------- |
| `#FID  IID  …` or `FID  IID  …` | column 2 (IID) | column 3            |
| anything else             | column 1     | column 2              |

When `#FID`/`FID` + `IID` is detected, the FID column is silently ignored.

| Flag                        | Purpose                                                    |
| --------------------------- | ---------------------------------------------------------- |
| `--pheno-name COL[,COL2]`  | Select phenotype column(s) from `--pheno`                  |
| `--resid-name COL1[,COL2]` | Select pre-computed residual columns from `--pheno`        |

`--pheno-name` selects columns by name.  Required for phenotype-based methods.
The interpretation depends on the method:

| Method          | Column(s)                                |
| --------------- | ---------------------------------------- |
| SPAsqr          | quantitative phenotype                   |
| POLMM           | ordinal phenotype                        |
| WtCoxG, LEAF    | binary *or* survival `TIME,EVENT`        |

For WtCoxG/LEAF survival analysis, pass `TIME,EVENT` as two comma-separated
column names.

`--resid-name` selects pre-computed residual columns for residual-based
methods (SPACox, SPAGRM, SAGELD, SPAmix, SPAmixPlus, SPAmixLocalPlus).

Example:

```
#FID  IID       AGE  BMI  BinaryPheno  SurvTime  SurvEvent  QuantPheno
FAM1  SAMPLE01  55   1.2  1            4.5       1          12.3
FAM2  SAMPLE02  43   0.8  0            NA        0          8.1
```

### Missing values

Any of the following tokens (case-sensitive where shown) are treated as
missing:

    .    NA    Na    na    nA    NaN    nan    Nan    NAN    N/A    n/a    NULL    null    Null

Missing values in phenotype and residual columns are preserved as `NaN`
internally; subjects with `NaN` in the active column are dropped during
subject filtering (see [data_filter.md](data_filter.md)).

---

## Covariate file — `--covar FILE`

Same format as `--pheno`: whitespace-delimited, mandatory header.
Subject ID detection follows the same `#FID`/`FID` + `IID` rule as `--pheno`.
All data column names must match `[0-9A-Za-z_\-.]+`.

| Flag                          | Purpose                                          |
| ----------------------------- | ------------------------------------------------ |
| `--covar-name COL1,COL2,...`  | Select named columns as covariates               |
| `--pc-cols COL1,COL2,...`     | PC columns (default: `PC1,PC2,PC3,PC4`)          |

### Covariate loading rules

| Combination                        | Behavior                                              |
| ---------------------------------- | ----------------------------------------------------- |
| `--covar FILE --covar-name COLS`   | Load named columns from the covariate file            |
| `--covar FILE` (no `--covar-name`) | Load **all** data columns from the covariate file     |
| `--covar-name COLS` (no `--covar`) | Load named columns from the `--pheno` file instead    |

Missing covariate values are filled with the per-column mean (computed from
the covariate file, excluding `NaN` entries).  If an entire column is missing,
an error is thrown.

**Covariates are not part of the subject intersection.**  A subject that
passes all filtering steps but is absent from `--covar` will have its
covariate values filled with the column means.  An intercept is added
automatically.

After loading, a covariate imputation summary is logged:

```
--covar AGE: filled 12 missing values with mean 54.312000
--covar PC1: filled 3 missing values with mean 0.002140
```

Example:

```
#FID  IID       AGE  SEX  PC1      PC2
FAM1  SAMPLE01  55   1    0.012   -0.003
FAM2  SAMPLE02  43   0    0.005    0.021
```

---

## Pre-computed residuals — `--pheno FILE --resid-name COLS`

Pre-computed residuals from an external null model (e.g., R).  The residual
columns are selected from the `--pheno` file using `--resid-name COL1,COL2,...`.
Both `--pheno` and `--resid-name` are required for residual-based methods.

Subject order need not match the genotype file; subjects are matched by IID.

### Column layout by method

| Method                                               | `--resid-name` columns                                    |
| ---------------------------------------------------- | --------------------------------------------------------- |
| SPACox, SPAGRM, SPAmix, SPAmixPlus, SPAmixLocalPlus | `RESID[,RESID2,...]` — one GWAS per column                |
| SAGELD                                               | `R_G,R_<E1>,R_Gx<E1>[,R_<E2>,R_Gx<E2>,...]`              |

For multi-column residuals, each column produces a separate output file:
`PREFIX.PHENO.METHOD[.gz|.zst]`.  Per-phenotype subject masks, allele counts,
and QC stats are computed independently.

Example:

```
#IID       Residual1   Residual2
SAMPLE01   0.523       -0.112
SAMPLE02  -0.311        0.847
```

---

## Sparse GRM

Exactly one of these is accepted (mutually exclusive):

### GRAB format — `--sp-grm-grab FILE`

Whitespace-delimited with a mandatory header.  The header must start with
one of: `ID1`, `#ID1`, `IID1`, `#IID1`.

```
#ID1      ID2       VALUE
SAMPLE01  SAMPLE01  1.000
SAMPLE01  SAMPLE02  0.125
SAMPLE02  SAMPLE02  1.000
```

Subjects not present in the genotype file are silently dropped.

### plink2 format — `--sp-grm-plink2 FILE`

The `.grm.sp` file uses 0-based indices (`idx1 idx2 value`).  A companion
`.grm.id` file (`FID IID` per line) is auto-detected from the file path.
If `.grm.id` is absent, indices are assumed to match `.fam` order (a warning
is logged).

---

## LOCO predictions — `--pred-list FILE`

Regenie step 1 `pred.list` file for LOCO (leave-one-chromosome-out) analysis.
Each line maps a phenotype name to its LOCO prediction file:

```
pheno1  /path/to/pheno1_loco.txt
pheno2  /path/to/pheno2_loco.txt
```

Each LOCO file is whitespace-delimited with a header row of subject IDs followed
by one row per chromosome.  The first column is the chromosome label and the
remaining columns are per-subject LOCO scores:

```
IID_1  IID_2  IID_3  ...
1      0.123  0.456  0.789  ...
2      0.234  0.567  0.890  ...
...
22     0.345  0.678  0.901  ...
```

Requirements:

- All 22 autosomal chromosomes (1–22) must be present; GRAB throws an error
  otherwise.
- The header subject IDs are matched to `--covar` subjects; every covariate
  subject must appear in the LOCO file.
- Non-autosomal rows (chr 23, X, Y, etc.) are silently skipped.

Used by `--method SPAsqr` when LOCO analysis is requested.

---

## Pairwise IBD — `--pairwise-ibd FILE`

Tab-separated with a mandatory header:

```
#ID1      ID2       pa      pb      pc
SAMPLE01  SAMPLE02  0.750   0.250   0.000
```

Produced by `grab --cal-pairwise-ibd`.

---

## Reference allele frequency — `--ref-af FILE`

Two formats are accepted:

**(A) plink2 `.afreq`:**

```
#CHROM  ID   REF  ALT  ALT_FREQS  OBS_CT
21      rs1  A    G    0.15       5000
```

Matched to `.bim` by `(CHROM, ID)` with allele flip detection.

**(B) Two-column numeric** (rows in `.bim` order):

```
0.15  5000
0.32  5000
```

For LEAF: comma-separated, one file per reference population.

---

## Individual AF model — `--ind-af-coef FILE`

Pre-computed model produced by `grab --cal-af-coef`.  Text or gzip format:

```
#STATUS  BETA0  BETA1  ...
```

Rows in filtered-marker order (positional matching with `.bim`).

---

## Admixed phi — `--admix-phi FILE`

Pre-computed phi kinship file produced by `grab --cal-phi`.  Tab-separated
with 0-based indices into `.fam` row order:

```
idx1  idx2  anc0_A  anc0_B  anc0_C  anc0_D  anc1_A  anc1_B  anc1_C  anc1_D
0     1     0.25    0.00    0.00    0.00    0.10    0.00    0.00    0.00
```

Four kinship scenarios (A/B/C/D) per ancestry.

---

## Subject filter files — `--keep FILE` / `--remove FILE`

Subject ID files for include/exclude filtering.  The format is auto-detected
from the first line:

| Format                                  | IID column |
| --------------------------------------- | ---------- |
| `#FID IID ...` (header with `#`)        | column 2   |
| `#IID ...` (header with `#`)            | column 1   |
| Two columns, no `#` header              | column 2   |
| Single column, no `#` header            | column 1   |

PLINK2-compatible.

---

## Marker filter files — `--extract FILE` / `--exclude FILE`

Single-column files of SNP IDs (one per line, matching `.bim` column 2).
Applied to both `--bfile` and `--admix-bfile` inputs.

---

## Chromosome filter — `--chr`

Restrict analysis to the specified chromosomes.  Accepts a comma-separated
list of chromosome numbers and/or ranges:

```
--chr 5
--chr 2,3
--chr 1-4,6-8,22
```

Chromosomes are matched as strings against the CHROM column (`.bim` column 1,
`.pvar` `#CHROM`, VCF `CHROM`, or BGEN chromosome field).  Applied to all
genotype formats (`--bfile`, `--pfile`, `--vcf`, `--bgen`).

---

## Marker QC flags

| Flag     | Default | Description                                        |
| -------- | ------- | -------------------------------------------------- |
| `--geno` | 0.1     | Per-marker missing rate cutoff                     |
| `--maf`  | 1e-5    | Minimum minor allele frequency                     |
| `--mac`  | 10      | Minimum minor allele count                         |
| `--hwe`  | 0       | HWE p-value threshold (0 = disabled)               |

Markers failing any cutoff are reported with `NA` for all result columns.
The QC stats (MISS_RATE, ALT_FREQ, MAC, HWE_P) are still printed.

---

## Output flags

| Flag                   | Default      | Description                                   |
| ---------------------- | ------------ | --------------------------------------------- |
| `--out PREFIX`         | —            | Output file prefix                            |
| `--compression gz|zst` | plain text   | Output compression format                     |
| `--compression-level`  | 0 (default)  | Compression level (gz: 1–9, zst: 1–22)        |

Output files are tab-delimited with one header line.  Plain text is the
default.

---

## Runtime flags

| Flag             | Default                | Description                                      |
| ---------------- | ---------------------- | ------------------------------------------------ |
| `--threads INT`  | 1                      | Number of worker threads                         |
| `--chunk-size`   | 8192 (min: 256)        | Markers per processing chunk                     |
| `--seed INT`     | 0 (random device)      | Random seed for reproducibility                  |

---

## Method-specific flags

| Flag                             | Default              | Used by                              |
| -------------------------------- | -------------------- | ------------------------------------ |
| `--prevalence FLOAT`             | —                    | —                                    |
| `--batch-effect-p-threshold`     | 0.05                 | —                                    |
| `--covar-p-threshold`            | 5e-5                 | —                                    |
| `--spa-z-threshold`              | 2.0                  | SPA methods                          |
| `--outlier-iqr-threshold`        | 1.5                  | —                                    |
| `--outlier-abs-bound`            | 0.55                 | —                                    |
| `--min-maf-ibd`                  | 0.01                 | `--cal-pairwise-ibd`                 |
| `--leaf-nclusters INT`           | from `--ref-af` count | LEAF                                |
| `--spasqr-taus LIST`             | 0.1,0.3,0.5,0.7,0.9 | SPAsqr                               |
| `--spasqr-tol FLOAT`            | 1e-7                 | SPAsqr                               |
| `--spasqr-h FLOAT`              | —                    | SPAsqr (exclusive with `--spasqr-h-scale`) |
| `--spasqr-h-scale FLOAT`        | 3                    | SPAsqr (exclusive with `--spasqr-h`)       |

---

## GWAS methods and their required inputs

| Method          | Genotype          | Phenotype          | GRM               | Other                         |
| --------------- | ----------------- | ------------------ | ------------------ | ----------------------------- |
| SPACox          | `--bfile/pfile/…` | `--resid-name`     | —                  | —                             |
| SPAGRM          | `--bfile/pfile/…` | `--resid-name`     | `--sp-grm-*`       | `--pairwise-ibd`              |
| SAGELD          | `--bfile/pfile/…` | `--resid-name`     | `--sp-grm-*`       | `--pairwise-ibd`              |
| SPAmix          | `--bfile/pfile/…` | `--resid-name`     | —                  | `--ind-af-coef`               |
| SPAmixPlus      | `--bfile/pfile/…` | `--resid-name`     | `--sp-grm-*`       | `--ind-af-coef`               |
| SPAmixLocalPlus | `--admix-bfile`   | `--resid-name`     | —                  | `--admix-phi`                 |
| POLMM           | `--bfile/pfile/…` | `--pheno-name`     | `--sp-grm-*`       | —                             |
| SPAsqr          | `--bfile/pfile/…` | `--pheno-name`     | `--sp-grm-*`       | `--pred-list` (LOCO)          |
| WtCoxG          | `--bfile/pfile/…` | `--pheno-name`     | `--sp-grm-*`       | `--ref-af`                    |
| LEAF            | `--bfile/pfile/…` | `--pheno-name`     | `--sp-grm-*`       | `--ref-af`                    |
