# GRAB File Format Reference

This page describes the file formats accepted and produced by GRAB
(Genome-Wide Robust Analysis for Biobank Data).

Input text files are **whitespace-delimited**
(tabs and/or spaces). Lines beginning with `##` are skipped as comments.
Output files are tab-delimited with one header line.

Compressed output (`.gz` or `.zst`) is controlled by `--compression` and
`--compression-level`. Plain text is the default.

---

## Subject set intersection

GRAB builds the analysis subject set through a multi-step pipeline.
See [subject_marker_pipeline.md](subject_marker_pipeline.md) for the full
specification and log format.

Key points (filtering order):

1. **genotype** — all subjects in `.fam` / `.psam` / VCF / BGEN sample block.
2. **∩ GRM** — intersect with sparse GRM subjects (if a GRM is provided).
   GRM subject IDs are pre-parsed before the full GRM load.
3. **∩ keep** — apply `--keep` filter (if provided).
4. **\ remove** — apply `--remove` filter (if provided).
5. **∩ pheno/resid** — intersect with phenotype / residual (subjects with
   non-NaN values).  Applied last so multi-residual analyses can compute
   per-phenotype masks on the already-filtered set.

**Covariates are not part of the intersection.**  A subject that passes
all pipeline steps but is missing from `--covar` will have its covariate
values filled with the per-column mean.  This avoids discarding subjects
just because a covariate is unavailable.

### Example

Suppose the genotype file has subjects {A, B, C, D, E}, the GRM has
{B, C, D}, and the residual file has {A, B, C, F}.

* genotype = {A, B, C, D, E}
* ∩ GRM = {B, C, D}
* No `--keep` / `--remove`, so stays {B, C, D}
* ∩ pheno/resid = {B, C}   ← residual has {A, B, C, F}
* If `--keep` lists {A, B}: after keep = {B}, final ∩ resid = {B}.
* The covariate file is not part of the intersection — if it has {B, C, D},
  subject B gets its values from the file; if it were missing, B would get
  per-column means.

### Input file categories

Files are grouped by what they are keyed on. Header and IID requirements
differ across groups.

#### (a) Per-subject files

These files have one row per subject. Two lines are shown when both headered
and headerless modes are supported.

| Flag                              | Subject matching                                      | plink2-compatible   |
| --------------------------------- | ----------------------------------------------------- | ------------------- |
| `--pheno` (headered)              | match by `IID` column                                 | Yes                 |
| `--covar` (headered)              | match by `IID` column                                 | Yes                 |
| `--covar` (numeric matrix)        | assume identical subjects as genotype file            | Yes                 |
| `--null-resid` (headered)         | match by `IID` column                                 | —                   |
| `--null-resid` (numeric matrix)   | assume identical subjects as genotype file            | —                   |
| `--keep` / `--remove` (auto)      | —                                                     | Yes                 |

#### (b) Per-subject-pair files

These files have one row per pair of subjects.

| Flag                             | Subject matching                                                   | plink2-compatible                        |
| -------------------------------- | ------------------------------------------------------------------ | ---------------------------------------- |
| `--sp-grm-grab` (headered)       | match by `ID1` `ID2` columns                                       | —                                        |
| `--sp-grm-plink2` (headerless)   | match by `.grm.id` or assume identical subjects as genotype file   | produced by `plink2 --make-grm-sparse`   |
| `--pairwise-ibd` (headered)      | match by `ID1` `ID2` columns                                       | —                                        |
| `--admix-phi` (headerless)       | assume identical subjects as genotype file                         | —                                        |

---

## Genotype input

See [genotype.md](genotype.md) for supported genotype formats, QC filtering,
and the complete genotype processing workflow.

---

## Phenotype file — `--pheno FILE`  *(plink2-compatible)*

Whitespace-delimited with a mandatory header containing `IID` (or `#IID`).
One row per subject. See [plink2 phenotype file format](https://www.cog-genomics.org/plink/2.0/input#pheno)
for full specification.

Related flags:

| Flag                        | Purpose                                             |
| --------------------------- | --------------------------------------------------- |
| `--pheno-binary COL`        | Select a 0/1 case/control column (WtCoxG, LEAF)     |
| `--pheno-surv TIME:EVENT`   | Select survival time + event columns (WtCoxG, LEAF) |
| `--pheno-quant COL`         | Select a quantitative column (SPAsqr)               |
| `--pheno-ordinal COL`       | Select an ordinal column (POLMM)                    |

Example:

```
#FID  IID       AGE  BMI  BinaryPheno  SurvTime  SurvEvent  QuantPheno
FAM1  SAMPLE01  55   1.2  1            4.5       1          12.3
FAM2  SAMPLE02  43   0.8  0            NA        0          8.1
```

### Missing values

Any of the following tokens represent a missing value:

    .    NA    na    Na    nA    NaN    nan    Nan    NAN    N/A    n/a    NULL    null    Null

Missing values are preserved as `NaN` internally.

---

## Covariate file — `--covar FILE`  *(plink2-compatible)*

Two formats are supported:

**Format A: IID-keyed (recommended).** Same plink2-compatible format as
`--pheno`; subjects are matched by `IID`.

**Format B: Pure numeric matrix.** All columns numeric, no header. Row count
must equal the genotype subject count; rows in genotype order.

Related flags:

| Flag                            | Purpose                                             |
| ------------------------------- | --------------------------------------------------- |
| `--covar-name COL1,COL2,...`    | Select named columns as covariates                  |
| `--covar-col-nums 3,5,7-10`     | Select by 1-based column position                   |
| `--not-covar COL1,COL2,...`     | Exclude named columns                               |
| `--pc-cols COL1,COL2,...`       | PC columns (default: `PC1,PC2,PC3,PC4`)             |

When none of `--covar-name`, `--covar-col-nums`, `--not-covar` are given, all
columns are loaded **except** `FID`, `IID`, `SID`, `PAT`, `MAT`, `SEX`, and
columns starting with `PHENO`. Missing values are imputed with the column mean.

Example:

```
#FID  IID       AGE  SEX  PC1      PC2
FAM1  SAMPLE01  55   1    0.012   -0.003
FAM2  SAMPLE02  43   0    0.005    0.021
```

---

## Null-model residual file — `--null-resid FILE`

Pre-computed residuals from an external null model (e.g., R). Two formats are
supported: **IID-keyed** (with header) and **pure numeric matrix** (no header).

### Format A: IID-keyed (recommended)

A whitespace-delimited text file with a header containing `IID` (or `#IID`).
The remaining columns are residual values.

```
#IID       RESID1   RESID2
SAMPLE01   0.523    -0.112
SAMPLE02  -0.311     0.847
```

Subject order need not match the genotype file; subjects are matched by IID.

### Format B: Pure numeric matrix

All columns are numeric (no IID/FID header). Row count must equal the genotype
subject count, and rows are assumed to be in genotype order.

```
 0.523  -0.112
-0.311   0.847
```

### Column layout by method

| Method                                                | Columns                                                       |
| ----------------------------------------------------- | ------------------------------------------------------------- |
| SPACox, SPAGRM, SPAmix, SPAmixPlus, SPAmixLocalPlus   | `RESID [RESID2 ...]` — one GWAS per column                    |
| WtCoxG, LEAF                                          | `RESID  WEIGHT  INDICATOR` — 3 fixed columns                  |
| SPAsqr                                                | `R_tau1  R_tau2  ...  R_tauK` — one column per quantile level |
| SAGELD                                                | `R_G  R_<E1>  R_Gx<E1>  [R_<E2>  R_Gx<E2>  ...]`              |

For multi-column residual files (SPACox, SPAGRM, SPAmix, SPAmixPlus,
SPAmixLocalPlus), each column produces a separate output file:
`PREFIX.PHENO.METHOD[.gz|.zst]`. Per-phenotype subject masks, allele counts,
and QC stats are independent.

---

## Reference allele frequency — `--ref-af FILE`

Two formats accepted (see also [Marker set intersection](#marker-set-intersection)):

**Format A: PLINK 2 `.afreq`** *(plink2-compatible, produced by `plink2 --freq`)* —
matched to genotype variants by `(CHROM, ID)` with automatic allele flip
detection.

```
#CHROM  ID  REF  ALT  ALT_FREQS  OBS_CT
1       rs1 A    G    0.35       1000
```

**Format B: Two-column numeric** — no header; rows in genotype variant order.

```
0.35  1000
```

For LEAF: multiple files are comma-separated on the command line, one per
reference population.

---

## Sparse GRM

Provide exactly one of the following.

### `--sp-grm-grab FILE`

GRAB 3-column format: whitespace-delimited.

**Header required:** Yes. The first non-comment line must be a header starting
with `ID1` (or `#ID1`). The three columns are `ID1`, `ID2`, and `VALUE`.
Lines beginning with `#` (after the header) are skipped as comments.

```
ID1       ID2       VALUE
SAMPLE01  SAMPLE01  1.002
SAMPLE01  SAMPLE02  0.498
```

Subjects not in the genotype file are silently dropped.

### `--sp-grm-plink2 FILE`  *(plink2-compatible)*

PLINK 2 sparse GRM (`.grm.sp`) format. A companion `.grm.id` file (headerless,
`FID  IID` per line — as produced by `plink2 --make-grm-sparse`) is
auto-detected by replacing `.grm.sp` with `.grm.id` in the file path. If
`.grm.id` is absent, 0-based indices matching genotype order are assumed.

---

## Subject filtering

### `--keep FILE` / `--remove FILE`  *(plink2-compatible)*

Restrict or exclude analysis subjects. Four formats are auto-detected:

| Format                  | Example                  |
| ----------------------- | ------------------------ |
| PLINK 2 with `#FID`     | `#FID  IID\nFAM1  ID1`   |
| PLINK 2 with `#IID`     | `#IID\nID1`              |
| Two-column (FID IID)    | `FAM1  ID1`              |
| Single-column (IID)     | `ID1`                    |

### `--extract FILE` / `--exclude FILE`  *(plink2-compatible)*

Single-column file of SNP IDs (one per line, matching `.bim` column 2).

---

## Marker set intersection

Markers (variants) follow a similar intersection logic as subjects. The
genotype file's variant list (`.bim`, `.pvar`, VCF records, or BGEN variant
metadata) is the starting set.

1. **Apply `--extract` / `--exclude`.** Restrict or remove markers by SNP ID.
   Both filter by exact ID string match against the genotype variant list.
2. **Match per-marker input files.** Files like `--ref-af` (`.afreq` format)
   are matched to genotypes by `(CHROM, ID)` with allele-flip detection.
   Positional per-marker files (`--ref-af` 2-column numeric, `--ind-af-coef`)
   assume rows in the same order as the (filtered) genotype variant list.
3. **Intersect `--bfile` with `--admix-bfile`.** When both are provided
   (SPAmixLocalPlus), `--extract` / `--exclude` are applied to both, and only
   markers present in both `.bim` files are analysed. Markers are matched by
   SNP ID.

### Per-marker input files

| Flag                           | Header                            | Key columns                        | plink2-compatible                        |
| ------------------------------ | --------------------------------- | ---------------------------------- | ---------------------------------------- |
| `--ref-af` (`.afreq`)          | Mandatory                         | `#CHROM  ID`                       | produced by `plink2 --freq`              |
| `--ref-af` (2-col numeric)     | —                                 | —                                  | —                                        |
| `--ind-af-coef`                | `#STATUS` (skipped)               | — (positional)                     | —                                        |
| `--extract` / `--exclude`      | No                                | Single SNP ID column               | Yes                                      |

---

## SNP-level individual AF model — `--ind-af-coef FILE`

Produced by `grab --cal-af-coef --out PREFIX`. Output is `PREFIX.afc[.gz|.zst]`.

```
#STATUS  BETA0  BETA1  ...
0        0.35   0.02   ...
1        0.35   0.02   0.001  ...
2        0.35   0.02  -0.003  ...
```

One row per filtered marker, in the same order as the genotype variant list
(positional matching — no CHROM or ID columns). The header line `#STATUS ...`
is skipped on read.

`STATUS` is an integer: `0` = uniform (MAC ≤ 20), `1` = OLS, `2` = logistic.
`BETA0..BETAK` are the regression coefficients (intercept + one per PC).
Row count must equal the number of filtered markers.

---

## Pairwise IBD file — `--pairwise-ibd FILE`

Produced by `grab --cal-pairwise-ibd --out PREFIX`. Output is `PREFIX.ibd[.gz|.zst]`.

**Header required:** Yes. The first non-comment line must be a header starting
with `ID1` (or `#ID1`). The five columns are `ID1`, `ID2`, `pa`, `pb`, `pc`.

Tab-separated with header:

```
#ID1      ID2       pa          pb          pc
SAMPLE01  SAMPLE02  0.00123     0.24500     0.75377
```

`pa`, `pb`, `pc` are probabilities of sharing 2, 1, 0 alleles IBD respectively.

---

## Phi kinship file — `--admix-phi FILE`

Produced by `grab --cal-phi --out PREFIX`. Output is `PREFIX.phi[.gz|.zst]`.

Wide tab-separated format with one row per related pair:

```
idx1  idx2  anc0_A     anc0_B     anc0_C     anc0_D     anc1_A     anc1_B     ...
0     5     0.12345    0.00000    0.00000    0.12345    0.23456    0.00000    ...
```

Indices are 0-based into `.fam` row order. Four phi matrix values (A, B, C, D)
per ancestry.

---

## Admixed ancestry genotypes (.abed) — `--admix-bfile PREFIX`

A GRAB-specific binary format for ancestry-specific dosage and hapcount data.
Produced by `grab --make-abed --out PREFIX`.

| File            | Description                                                   |
| --------------- | ------------------------------------------------------------- |
| `PREFIX.abed`   | BGZF-compressed binary ancestry dosage/hapcount matrix        |
| `PREFIX.bim`    | Standard PLINK `.bim` variant file (shared with PLINK)        |
| `PREFIX.fam`    | Standard PLINK `.fam` subject file (shared with PLINK)        |

### Binary layout

The `.abed` file is BGZF-compressed (indexed via `.abed.gzi`).

**Header (8 bytes):**

| Offset   | Size   | Field         | Description                                                      |
| -------- | ------ | ------------- | ---------------------------------------------------------------- |
| 0        | 2      | `magic`       | `{0xAD, 0x4D}`                                                   |
| 2        | 1      | `version`     | `0x02`                                                           |
| 3        | 1      | `nAnc`        | Bits 0–6: number of ancestries K (1–127). Bit 7: NO_MISSING flag |
| 4        | 4      | `nSamples`    | Subject count N (little-endian uint32)                           |

Marker count M is derived from the companion `.bim` file.

**Body (SNP-major):**

For each of the M markers, 2K tracks are stored sequentially:

```
[dosage_anc0][hapcount_anc0][dosage_anc1][hapcount_anc1]...[dosage_anc(K-1)][hapcount_anc(K-1)]
```

Each track is `ceil(N / 4)` bytes using **PLINK-compatible 2-bit encoding**
(4 subjects per byte, LSB first):

| 2-bit code   | Value                    |
| ------------ | ------------------------ |
| `00`         | 0 (reference homozygote) |
| `10`         | 1 (heterozygote)         |
| `11`         | 2 (alternate homozygote) |
| `01`         | missing                  |

When the NO_MISSING flag is set (bit 7 of `nAnc`), code `01` is treated as 0.

**Total file size:** `8 + M × 2K × ceil(N/4)` bytes (before BGZF compression).

### Production

Two input modes (mutually exclusive):

- `grab --make-abed --vcf FILE --rfmix-msp FILE --out PREFIX` — phased VCF/BCF
  + rfmix2 MSP → `.abed`
- `grab --make-abed --admix-text-prefix PREFIX --out PREFIX` — extract_tracts
  text output → `.abed`

### Usage

Pass `PREFIX` as `--admix-bfile PREFIX` to SPAmixLocalPlus or `--cal-phi`.

---

## MSP local-ancestry file — `--rfmix-msp FILE`

rfmix2 output format accepted by `--make-abed`:

```
#Subpopulation order/codes: 0=POP0\t1=POP1\t...
#chm  spos  epos  sgpos  egpos  n snps  IID0.0  IID0.1  ...
1     0     50000 0.0    0.5    100      0       1       ...
```

Line 1 defines ancestry codes (K inferred from count).
Line 2 is the column header.
Data lines: chromosome, start/end positions (0-based, exclusive end),
genetic positions, SNP count, then per-haplotype ancestry calls.

---

## Extract-tracts text format — `--admix-text-prefix PREFIX`

Text files produced by extract_tracts:

```
{PREFIX}.anc{k}.dosage[.gz]
{PREFIX}.anc{k}.hapcount[.gz]
```

For k = 0, 1, ... (K auto-detected from file presence). Each file has:

```
CHROM  POS  ID  REF  ALT  SAMPLE1  SAMPLE2  ...
1      100  rs1 A    G    0        2        ...
```

Values are integers 0–2.

---

## Output files

All GWAS methods write to `--out PREFIX` with a suffix determined by the method.

### GWAS output naming

| Method               | Suffix                | Output path                                  |        |
| -------------------- | --------------------- | -------------------------------------------- | ------ |
| SPACox               | `.PHENO.SPACox`       | `PREFIX.PHENO.SPACox[.gz\                    | .zst]` |
| SPAGRM               | `.PHENO.SPAGRM`       | `PREFIX.PHENO.SPAGRM[.gz\                    | .zst]` |
| SPAmix               | `.PHENO.SPAmix`       | `PREFIX.PHENO.SPAmix[.gz\                    | .zst]` |
| SPAmixPlus           | `.PHENO.SPAmixP`      | `PREFIX.PHENO.SPAmixP[.gz\                   | .zst]` |
| SPAmixLocalPlus      | `.PHENO.LocalP`       | `PREFIX.PHENO.LocalP[.gz\                    | .zst]` |
| SAGELD               | `.SAGELD`             | `PREFIX.SAGELD[.gz\                          | .zst]` |
| POLMM                | `.POLMM`              | `PREFIX.POLMM[.gz\                           | .zst]` |
| SPAsqr               | `.SPAsqr`             | `PREFIX.SPAsqr[.gz\                          | .zst]` |
| WtCoxG               | `.WtCoxG`             | `PREFIX.WtCoxG[.gz\                          | .zst]` |
| LEAF                 | `.LEAF`               | `PREFIX.LEAF[.gz\                            | .zst]` |

For multi-residual methods (SPACox, SPAGRM, SPAmix, SPAmixPlus,
SPAmixLocalPlus), `PHENO` is the column name from `--null-resid`.

### Utility output naming

| Mode                      | Suffix                    | Output path                                |        |
| ------------------------- | ------------------------- | ------------------------------------------ | ------ |
| `--cal-af-coef`            | `.afc`                    | `PREFIX.afc[.gz\                           | .zst]` |
| `--cal-pairwise-ibd`      | `.ibd`                    | `PREFIX.ibd[.gz\                           | .zst]` |
| `--cal-phi`               | `.phi`                    | `PREFIX.phi[.gz\                           | .zst]` |
| `--make-abed`             | `.abed` / `.bim` / `.fam` | `PREFIX.abed`, `PREFIX.bim`, `PREFIX.fam`  |        |

### Common GWAS output columns

All GWAS output files are tab-separated with one header line. The first columns
are always marker metadata:

```
CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
```

| Column       | Description                                        |
| ------------ | -------------------------------------------------- |
| `CHROM`      | Chromosome                                         |
| `POS`        | Base-pair position                                 |
| `ID`         | Variant ID (from `.bim`)                           |
| `REF`        | Reference allele                                   |
| `ALT`        | Alternate allele                                   |
| `MISS_RATE`  | Per-marker missing rate                            |
| `ALT_FREQ`   | Alternate allele frequency (after QC filtering)    |
| `MAC`        | Minor allele count                                 |
| `HWE_P`      | Hardy-Weinberg equilibrium p-value                 |

Method-specific columns follow the metadata columns (see `grab --help METHOD`
for each method's output column details).

---

## Compression

Output compression is controlled by two flags:

| Flag                   | Values                              | Default               |
| ---------------------- | ----------------------------------- | --------------------- |
| `--compression`        | `gz` (gzip) or `zst` (Zstandard)    | plain text            |
| `--compression-level`  | gz: 1–9; zst: 1–22                  | 0 (library default)   |

Input files with `.gz` or `.zst` extensions are auto-detected and decompressed
transparently.
