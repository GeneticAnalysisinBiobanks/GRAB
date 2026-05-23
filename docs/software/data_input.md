# Data Input Reference

This page enumerates the input files accepted by GRAB and the command-line
flags that govern them.  The authoritative definitions reside in
[src/cli/flags.hpp](../../src/cli/flags.hpp) and [src/cli/parse.cpp](../../src/cli/parse.cpp);
this page documents user-visible behavior.

All textual input files are **whitespace-delimited** (tabs and/or spaces).
Lines beginning with `##` are skipped as comments.

---

## Genotype files

Exactly one of the following flags is required for every GWAS method
except SPAmixLocalPlus (which consumes `--admix-bfile` instead):

| Flag             | Format                                       | Subject ID source                                            |
| ---------------- | -------------------------------------------- | ------------------------------------------------------------ |
| `--bfile PREFIX` | PLINK 1 (`.bed/.bim/.fam`)                   | `.fam` column 2 (IID)                                        |
| `--pfile PREFIX` | PLINK 2 (`.pgen/.pvar/.psam`)                | `.psam` IID column                                           |
| `--vcf FILE`     | VCF / BCF (`.vcf`, `.vcf.gz`, `.bcf`)        | VCF header sample IDs                                        |
| `--bgen FILE {ref-first\|ref-last\|ref-unknown}` | BGEN v1.1 / v1.2 / v1.3 (`.bgen`, `.sample`) | Embedded sample-identifier block or companion `.sample` file |

`--bfile`, `--pfile`, `--vcf`, and `--bgen` are mutually exclusive; the
dispatcher rejects any combination that specifies more than one.

For `--bgen`: BGEN does not encode which of the two listed alleles is
REF, so a REF/ALT mode token is mandatory and follows the plink2
`--bgen` syntax:

- `ref-first` — `alleles[0]` is REF, `alleles[1]` is ALT.  Use for
  IMPUTE / qctool / UK Biobank exports.
- `ref-last`  — `alleles[0]` is ALT, `alleles[1]` is REF.  Use for
  `.bgen` files produced by plink2 default `--export bgen-1.x`.
- `ref-unknown` — REF status unknown.  Accepted for plink2
  compatibility; GRAB treats it identically to `ref-last` and emits a
  "REF is provisional" warning because the output table has no PR
  marker column.

When the BGEN file contains an embedded sample-identifier block, those
IDs are used.  Otherwise GRAB looks for a companion `.sample` file
(Oxford format: header `ID_1 ID_2 missing`, type row `0 0 0`, then
`FID IID ...` data lines).  If neither source is available, numeric IDs
(`"0"`, `"1"`, ...) are generated and a warning is emitted.

### Admixed-ancestry genotypes — `--admix-bfile PREFIX`

SPAmixLocalPlus and `--cal-phi` consume the `.abed` format instead of
standard genotype input.  The `.abed` file stores 2K tracks
(dosage + hapcount per ancestry) in a BGZF-compressed binary that shares
the standard PLINK `.fam` and `.bim` files.

The 8-byte header layout (defined in [src/localplus/abed_io.hpp](../../src/localplus/abed_io.hpp)):

| Offset | Size | Field    | Value                                                |
| ------ | ---- | -------- | ---------------------------------------------------- |
| 0      | 2    | magic    | `0xAD 0x4D`                                          |
| 2      | 1    | version  | `0x02`                                               |
| 3      | 1    | nAnc     | K in bits 0–6; bit 7 = `NO_MISSING` flag             |
| 4      | 4    | nSamples | N (little-endian)                                    |

`nMarkers` is derived from the companion `.bim` file.  Each marker stores
2K tracks of `⌈N/4⌉` bytes using PLINK-compatible 2-bit encoding
(`00`→0, `10`→1, `11`→2, `01`→missing).  Track order per marker:
`[dosage_anc0][hapcount_anc0]...[dosage_anc(K-1)][hapcount_anc(K-1)]`.

The `.abed` file is produced by `grab --make-abed` from one of:

- `--vcf FILE --rfmix-msp FILE`  (phased VCF/BCF + rfmix2 MSP local-ancestry output)
- `--admix-text-prefix PREFIX`   (extract_tracts text dosage/hapcount files)

---

## Phenotype file — `--pheno FILE`

Whitespace-delimited with a mandatory header line and one row per subject.
Data-column names must match `[0-9A-Za-z_\-.]+`.

The subject-ID column is inferred from the header:

| Header columns 1–2                | Subject ID column | Data columns start at |
| --------------------------------- | ----------------- | --------------------- |
| `#FID  IID  …` or `FID  IID  …`   | column 2 (IID)    | column 3              |
| anything else                     | column 1          | column 2              |

When `#FID`/`FID` plus `IID` is detected, the FID column is silently
ignored; GRAB matches subjects by IID alone.

| Flag                          | Purpose                                                      |
| ----------------------------- | ------------------------------------------------------------ |
| `--pheno-name COL[,COL2]`     | Select phenotype column(s) from `--pheno`                    |
| `--resid-name COL1[,COL2]`    | Select pre-computed residual columns from `--pheno`          |
| `--regression-model MODEL`    | Declare the null-model regression family for `--pheno-name`  |
| `--save-resid`                | Write fitted residuals to `PREFIX.null.resid`                |

`--pheno-name` selects phenotype columns for in-process null-model
fitting (`--regression-model` controls the regression family) or for the
external-reference methods (WtCoxG, LEAF).  `--resid-name` selects
pre-computed residual columns from a previous fit.  The two flags are
**mutually exclusive**: pick one of the residual-input mode or the
phenotype-input (fit) mode.

### `--regression-model` and the fit path

`--regression-model` accepts one of:

| Value      | Regression family                                  | `--pheno-name` syntax                  |
| ---------- | -------------------------------------------------- | -------------------------------------- |
| `auto`     | per-token inference (default)                      | per spec                               |
| `linear`   | least-squares regression (intercept added automatically) | `Y1,Y2,...`                      |
| `logistic` | logistic regression (0/1 response; intercept added) | `Y1,Y2,...`                           |
| `cox`      | Cox proportional-hazards survival                  | `TIME1:EVENT1,TIME2:EVENT2,...`        |
| `ordinal`  | cumulative-logit (proportional odds) on `0..J−1`   | `Y1,Y2,...`                            |

Per-method support:

| Method                          | Accepted `--regression-model`                            |
| ------------------------------- | -------------------------------------------------------- |
| SPACox / SPAGRM / SPAmix /      | `auto`, `linear`, `logistic`, `cox`, `ordinal`           |
| SPAmixPlus / SPAmixLocalPlus    | (null model is fit in-process from `--pheno-name`)        |
| WtCoxG / LEAF                   | `auto`, `logistic`, `cox`                                |
| SPAsqr                          | `auto`, `linear` (purely declarative)                    |
| SAGELD                          | `auto` only (rejects any other value)                    |

`auto`-inference rules (see [src/util/null_model.cpp](../../src/util/null_model.cpp)):
constant or all-NaN columns are rejected; two distinct values → logistic
(recoded to 0/1); integer columns with ≤ 10 contiguous distinct values →
ordinal (shifted to start at 0); everything else → linear.  Cox survival
is signaled by the `TIME:EVENT` colon syntax in `--pheno-name`.

### `--save-resid`

When `--save-resid` is supplied with `--pheno-name`, GRAB writes the
fitted residuals to `PREFIX.null.resid` in the strict format consumed by
`--resid-name`.  The file can be reloaded in a later invocation via
`--pheno PREFIX.null.resid --resid-name ...`, avoiding repeated null-model
fits across follow-up analyses.

Example:

```
#FID  IID       AGE  BMI  BinaryPheno  SurvTime  SurvEvent  QuantPheno
FAM1  SAMPLE01  55   1.2  1            4.5       1          12.3
FAM2  SAMPLE02  43   0.8  0            NA        0          8.1
```

### Missing values

The following tokens are treated as missing:

```
.    NA    Na    na    nA    NaN    nan    Nan    NAN    N/A    n/a    NULL    null    Null
```

Missing entries in phenotype and residual columns are preserved as `NaN`
internally; subjects with `NaN` in the active column are dropped during
subject filtering (see [data_filter.md](data_filter.md)).

---

## Covariate file — `--covar FILE`

The format mirrors `--pheno`: whitespace-delimited with a mandatory
header, subject-ID detection following the same `#FID`/`FID` plus `IID`
rule, data-column names matching `[0-9A-Za-z_\-.]+`.

| Flag                          | Purpose                                          |
| ----------------------------- | ------------------------------------------------ |
| `--covar-name COL1,COL2,...`  | Select named columns as covariates               |
| `--pc-cols COL1,COL2,...`     | PC columns (default: `PC1,PC2,PC3,PC4`)          |

### Covariate-loading rules

| Combination                          | Behavior                                                  |
| ------------------------------------ | --------------------------------------------------------- |
| `--covar FILE --covar-name COLS`     | Load named columns from `--covar`                         |
| `--covar FILE` (no `--covar-name`)   | Load **all** data columns from `--covar`                  |
| `--covar-name COLS` (no `--covar`)   | Load named columns from `--pheno` instead                 |

Missing covariate values are imputed with the per-column mean (computed
from the covariate file, excluding `NaN` entries).  An entirely missing
column raises an error.

Covariates are **not** part of the subject intersection: a subject that
passes all upstream filters but is absent from `--covar` retains the
mean-imputed values.  An intercept column is appended automatically by the
regression routines that require one (linear, logistic, ordinal); Cox
does not include an intercept.

### `--pc-cols` semantics

`--pc-cols` selects columns from `--covar` (or `--pheno` if `--covar` is
absent) for use as principal components.  Its consumers are:

- LEAF (K-means clustering of subjects on PC space).
- SPAmix / SPAmixPlus (per-individual ALT-allele-frequency model).
- `--cal-af-coef` (precomputes the per-individual AF model).

`--pc-cols` does **not** enter the null-model design.  To adjust the
null model for PCs as well, list them explicitly in `--covar-name`.

Example:

```
#FID  IID       AGE  SEX  PC1      PC2
FAM1  SAMPLE01  55   1    0.012   -0.003
FAM2  SAMPLE02  43   0    0.005    0.021
```

---

## Pre-computed residuals — `--pheno FILE --resid-name COLS`

Pre-computed residuals from an external null model.  The residual columns
are selected from the `--pheno` file by name.  Both `--pheno` and
`--resid-name` are required for the residual-input mode of every
residual-based method.

Subject order in `--pheno` need not match the genotype file; subjects are
matched by IID.

### Column layout by method

| Method                                                | `--resid-name` columns                                       |
| ----------------------------------------------------- | ------------------------------------------------------------ |
| SPACox, SPAGRM, SPAmix, SPAmixPlus, SPAmixLocalPlus  | `RESID[,RESID2,...]` — one independent GWAS per column       |
| SAGELD                                                | `R_G, R_<E1>, R_Gx<E1>[, R_<E2>, R_Gx<E2>, ...]`             |

For SAGELD residual mode, the column layout is positional: `R_G` first,
followed by `(R_<Ee>, R_Gx<Ee>)` pairs in any number.  Environment names
are read from the header.

For multi-column residuals, each phenotype yields a separate output file
`PREFIX.PHENO.METHOD[.gz|.zst]`; per-phenotype subject masks, allele
counts, and QC statistics are computed independently.

Example:

```
#IID       Residual1   Residual2
SAMPLE01   0.523       -0.112
SAMPLE02  -0.311        0.847
```

---

## Sparse GRM

Exactly one of these is accepted; they are mutually exclusive.

### GRAB format — `--sp-grm-grab FILE`

Whitespace-delimited with a mandatory header.  The header must begin with
one of: `ID1`, `#ID1`, `IID1`, `#IID1`.

```
#ID1      ID2       VALUE
SAMPLE01  SAMPLE01  1.000
SAMPLE01  SAMPLE02  0.125
SAMPLE02  SAMPLE02  1.000
```

Subjects not present in the genotype file are silently dropped.

### plink2 format — `--sp-grm-plink2 FILE`

The `.grm.sp` file uses 0-based indices (`idx1 idx2 value`).  The
companion `.grm.id` file (`FID IID` per line) is auto-detected from the
file path.  If `.grm.id` is absent, the indices are assumed to match
`.fam` order and a warning is logged.

---

## SAGELD environment columns — `--sageld-x COL_IDS`

`--sageld-x` activates SAGELD's phenotype-input mode.  It accepts a
comma-separated list of environment column names; each name must
(a) match a numeric column of `--pheno`, and (b) be listed in
`--covar-name` so it enters the fixed-effect design.  For every
`(phenotype, env)` pair the null model `Y ~ X + (E | IID)` is fit by
EM-ML internally, and BLUP residuals are aggregated into per-IID
`(R_G, R_<E>, R_Gx<E>)` triples before the marker-level G and G×E
score tests run.  Multiple environments trigger one model per
environment.  `--sageld-x` is incompatible with `--resid-name`.

---

## LOCO predictions — `--pred-list FILE`

Step-1 polygenic predictions used by SPAsqr's LOCO mode.  The file is
space-separated with two columns: phenotype name and path to its `.loco`
file.

```
pheno1  /path/to/pheno1_loco.txt
pheno2  /path/to/pheno2_loco.txt
```

The `.loco` file format is auto-detected per phenotype.  Two formats are
supported:

- **Regenie** — header begins with `FID_IID`; chromosome-major rows
  (first column is the chromosome label, remaining columns are
  per-subject LOCO scores).
- **LDAK-KVIK** — header `FID IID Chr1 Chr2 ... Chr22`; subject-major
  rows (one row per subject, one column per chromosome).

Used by `--method SPAsqr` when `--pred-list` is supplied.  Every
phenotype listed in `--pheno-name` must appear in the pred.list,
otherwise dispatch fails.

---

## Pairwise IBD — `--pairwise-ibd FILE`

Tab-separated with a mandatory header:

```
#ID1      ID2       pa      pb      pc
SAMPLE01  SAMPLE02  0.750   0.250   0.000
```

Produced by `grab --cal-pairwise-ibd`.  Required by SPAGRM and SAGELD.

---

## Reference allele frequency — `--ref-af FILE`

Two formats are accepted:

**(A) plink2 `.afreq`:**

```
#CHROM  ID   REF  ALT  ALT_FREQS  OBS_CT
21      rs1  A    G    0.15       5000
```

Matched to the `.bim` rows by `(CHROM, ID)` with allele-flip detection.

**(B) Two-column numeric** (rows in `.bim` order):

```
0.15  5000
0.32  5000
```

For LEAF, `--ref-af` accepts a comma-separated list, one file per
reference population.

---

## Individual AF model — `--ind-af-coef FILE`

Pre-computed per-marker per-individual AF model produced by
`grab --cal-af-coef`.  The file may be plain text or gzip-compressed.

```
#STATUS  BETA0  BETA1  ...
```

Records are positionally matched to the filtered `.bim` rows.  Consumed
by SPAmix and SPAmixPlus.

---

## Admixed phi — `--admix-phi FILE`

Pre-computed phi-kinship file produced by `grab --cal-phi`.
Tab-separated with 0-based indices into `.fam` row order:

```
idx1  idx2  anc0_A  anc0_B  anc0_C  anc0_D  anc1_A  anc1_B  anc1_C  anc1_D
0     1     0.25    0.00    0.00    0.00    0.10    0.00    0.00    0.00
```

Four kinship scenarios (A/B/C/D) are stored per ancestry.  Consumed by
SPAmixLocalPlus.

---

## Subject filter files — `--keep FILE` / `--remove FILE`

Subject-ID files for include / exclude filtering.  The format is
auto-detected from the first line:

| Format                                  | IID column |
| --------------------------------------- | ---------- |
| `#FID IID ...` (header begins with `#`) | column 2   |
| `#IID ...` (header begins with `#`)     | column 1   |
| Two columns, no `#` header              | column 2   |
| Single column, no `#` header            | column 1   |

PLINK2-compatible.

---

## Marker filter files — `--extract FILE` / `--exclude FILE`

Single-column files of SNP IDs (one per line, matching `.bim` column 2
or `.pvar` ID column).  Applied to all genotype inputs (`--bfile`,
`--pfile`, `--vcf`, `--bgen`, `--admix-bfile`) in the genotype-factory
constructor, before any genotype is decoded.

---

## Chromosome filter — `--chr SPEC`

Restricts analysis to specific chromosomes.  `SPEC` is a comma-separated
list of chromosome numbers and/or ranges:

```
--chr 5
--chr 2,3
--chr 1-4,6-8,22
```

Chromosomes are matched as strings against the CHROM column (`.bim`
column 1, `.pvar` `#CHROM`, VCF `CHROM`, or the BGEN chromosome field).
Applied to every supported genotype format.

---

## Marker QC flags

| Flag     | Default | Description                                        |
| -------- | ------- | -------------------------------------------------- |
| `--geno` | 0.1     | Per-marker missing-rate cutoff                     |
| `--maf`  | 1e-5    | Minimum minor-allele frequency                     |
| `--mac`  | 10      | Minimum minor-allele count                         |
| `--hwe`  | 0       | HWE p-value threshold (0 = disabled)               |

Markers that fail any cutoff are reported with `NA` for every result
column; the QC statistics (`MISS_RATE`, `ALT_FREQ`, `MAC`, `HWE_P`) are
still printed.

---

## Output flags

| Flag                     | Default      | Description                                   |
| ------------------------ | ------------ | --------------------------------------------- |
| `--out PREFIX`           | —            | Output file prefix (required)                 |
| `--compression gz\|zst`  | plain text   | Output compression format                     |
| `--compression-level`    | 0 (default)  | Compression level (gz: 1–9, zst: 1–22)        |

Output files are tab-delimited with one header line.  Plain text is the
default.  Per-phenotype output paths follow
`PREFIX.PHENO.METHOD[.gz|.zst]` (see each method's documentation for the
method tag).

---

## Runtime flags

| Flag             | Default                | Description                                      |
| ---------------- | ---------------------- | ------------------------------------------------ |
| `--threads INT`  | 1                      | Number of worker threads (see note below)        |
| `--chunk-size`   | 8192 (min: 256)        | Markers per processing chunk                     |
| `--seed INT`     | 0 (random device)      | Random seed for reproducibility                  |

`--threads N` sizes the worker pool used by the chunk-level marker
engine and by every pre-marker work-stealing pool (null-model fits,
Chow-Liu MAF-bin tables, K-means restarts, batch-effect tests, etc.).
The **main thread is not counted in `N`**: it dispatches chunks, drains
output, and stays mostly idle, but briefly occupies one additional CPU
during loading, finalize, and synchronization steps.  When sizing a job
to physical cores, plan for `N + 1` logical cores; oversubscribing by
one is usually harmless because the main thread is rarely on-CPU at the
same time as all `N` workers.

---

## Method-specific flags

### SPA / outlier knobs

| Flag                             | Default | Used by                                       |
| -------------------------------- | ------- | --------------------------------------------- |
| `--spa-z-threshold FLOAT`        | 2.0     | All SPA methods                               |
| `--covar-p-threshold FLOAT`      | 5e-5    | SPACox (covariate-adjusted refit gate)        |
| `--outlier-iqr-multiplier FLOAT` | 1.5     | SPAGRM, SAGELD, SPAmix, SPAmixPlus,           |
|                                  |         | SPAmixLocalPlus, WtCoxG, LEAF, SPAsqr         |
| `--spagrm-control-outlier`       | off     | SPAGRM (iterative IQR-ratio adjustment)       |
| `--spasqr-outlier-abs-bound FLOAT` | 0.55  | SPAsqr (absolute outlier cutoff)              |

### Batch-effect / reference

| Flag                          | Default | Used by               |
| ----------------------------- | ------- | --------------------- |
| `--prevalence FLOAT`          | —       | WtCoxG, LEAF (required) |
| `--batch-effect-p-threshold`  | 0.05    | WtCoxG, LEAF          |

### IBD / clustering

| Flag                          | Default                     | Used by               |
| ----------------------------- | --------------------------- | --------------------- |
| `--min-maf-ibd FLOAT`         | 0.01                        | `--cal-pairwise-ibd`  |
| `--leaf-nclusters INT`        | from `--ref-af` count       | LEAF                  |

### SPAsqr

| Flag                          | Default                                  | Description                                 |
| ----------------------------- | ---------------------------------------- | ------------------------------------------- |
| `--spasqr-taus LIST`          | `0.1,0.3,0.5,0.7,0.9` (max 20)           | Quantile levels                             |
| `--spasqr-solver NAME`        | `qmme`                                   | `qmme` or `conquer`                         |
| `--spasqr-mode MODE`          | `score`                                  | `score` (default) or `wald`                 |
| `--spasqr-tol FLOAT`          | 1e-7 (QMME tightens to min(tol, 1e-9))   | SQR convergence tolerance                   |
| `--spasqr-h FLOAT`            | auto                                     | Explicit SQR bandwidth                      |
| `--spasqr-h-scale FLOAT`      | 3 (score) / 10 (wald)                    | IQR-divisor when `--spasqr-h` is not set    |
| `--pheno-transform MODE`      | `int`                                    | `raw`, `int`, or `standardize`              |
| `--pred-list FILE`            | —                                        | Enables LOCO mode                           |

`--spasqr-h` and `--spasqr-h-scale` are mutually exclusive.

`--spasqr-solver`:

- `qmme` (default) — Quadratic Majorization-Minimization with
  Extrapolation (Heng & Wang, 2025).  The Hessian upper bound is
  Cholesky-decomposed once per phenotype × bandwidth and reused across
  all τ levels; the solver is more robust on ill-conditioned designs.
- `conquer` — Convolution-type smoothed quantile regression (He et al.,
  2021): Huber initialization followed by Barzilai–Borwein gradient
  descent.  Each τ is refit from scratch.

`--spasqr-mode`:

- `score` (default) — score test on null-model residuals.  One null QR
  fit per phenotype × τ; markers are streamed and tested via
  `S = Σᵢ Rᵢ Gᵢ` with M-estimation sandwich variance.  Output
  columns: `P_CCT  P_tau{val}...  Z_tau{val}...`.
- `wald` — per-marker × per-τ full-model refit by QMME, with β̂_G and
  SE drawn from the (γ, γ) entry of the M-estimation sandwich
  `V = A⁻¹ B A⁻¹ / n`.  Output is long-format
  `CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P TAU BETA SE Z P`,
  with one row per `(marker, τ)`.  No GRM is used; this mode is intended
  for follow-up effect-size estimation on a small SNP list selected by
  `--extract`.

`--pheno-transform`:

- `raw` — no transform; SQR fits on Y as supplied.
- `int` (default) — inverse-normal transform (Blom plotting position,
  average-rank ties) applied per phenotype on its non-missing scope.
- `standardize` — Y is centered and scaled to unit variance per
  phenotype.

When `--pred-list` is supplied, the transformed Y has `loco_chr`
subtracted as an offset (β = 1, α = 0).  The LOCO PRS scale must match
the chosen transform; mixing scales produces inconsistent residuals and
the dispatcher warns when `--pheno-transform raw` or `standardize` is
combined with `--pred-list`.

### SAGELD

| Flag                  | Used by                                          |
| --------------------- | ------------------------------------------------ |
| `--sageld-x COL_IDS`  | Activates SAGELD pheno-input mode (see above)    |

---

## GWAS methods and their required inputs

| Method            | Genotype                | Phenotype                                      | GRM           | Other                       |
| ----------------- | ----------------------- | ---------------------------------------------- | ------------- | --------------------------- |
| SPACox            | `--bfile`/`--pfile`/`--vcf`/`--bgen` | `--resid-name` or `--pheno-name`     | —             | —                           |
| SPAGRM            | `--bfile`/`--pfile`/`--vcf`/`--bgen` | `--resid-name` or `--pheno-name`     | `--sp-grm-*`  | `--pairwise-ibd`            |
| SAGELD            | `--bfile`/`--pfile`/`--vcf`/`--bgen` | `--resid-name` or `--pheno-name + --sageld-x`| `--sp-grm-*`  | `--pairwise-ibd`            |
| SPAmix            | `--bfile`/`--pfile`/`--vcf`/`--bgen` | `--resid-name` or `--pheno-name + --pc-cols` | optional `--sp-grm-*` | `--ind-af-coef` (optional) |
| SPAmixPlus        | `--bfile`/`--pfile`/`--vcf`/`--bgen` | `--resid-name` or `--pheno-name + --pc-cols` | `--sp-grm-*`  | `--ind-af-coef` (optional) |
| SPAmixLocalPlus   | `--admix-bfile`         | `--resid-name` or `--pheno-name`               | (via `--admix-phi`) | `--admix-phi`         |
| SPAsqr            | `--bfile`/`--pfile`/`--vcf`/`--bgen` | `--pheno-name`                       | optional `--sp-grm-*` | `--pred-list` (optional, LOCO) |
| WtCoxG            | `--bfile`/`--pfile`/`--vcf`/`--bgen` | `--pheno-name`                       | optional `--sp-grm-*` | `--ref-af`, `--prevalence` |
| LEAF              | `--bfile`/`--pfile`/`--vcf`/`--bgen` | `--pheno-name + --pc-cols`           | optional `--sp-grm-*` | `--ref-af`, `--prevalence` |

---

## Utility modes

These modes do not run a GWAS; they pre-compute auxiliary files consumed
by the GWAS methods, or convert ancillary data into the formats GRAB
accepts.  Each utility flag is mutually exclusive with `--method`.

| Flag                  | Purpose                                                  |
| --------------------- | -------------------------------------------------------- |
| `--cal-af-coef`       | Pre-compute per-marker per-individual AF model           |
|                       | (`PREFIX.afc`; consumed by `--ind-af-coef`)              |
| `--cal-pairwise-ibd`  | Pre-compute pairwise IBD probabilities                   |
|                       | (`PREFIX.ibd`; consumed by `--pairwise-ibd`)             |
| `--cal-phi`           | Compute phi kinship matrices from admixed genotypes      |
|                       | (`PREFIX.phi`; consumed by `--admix-phi`)                |
| `--make-abed`         | Build `.abed` admixed binary from VCF+MSP or             |
|                       | extract_tracts text output                               |
| `--int-pheno`         | Apply inverse-normal transform to every column of a      |
|                       | phenotype file (writes `PREFIX.txt`)                     |
