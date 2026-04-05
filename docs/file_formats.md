# GRAB Input File Formats

Reference for every input file format accepted by `grab`.  
Each section describes the on-disk layout, column semantics, and
subject-matching rules.

---

## PLINK Binary Genotypes (`--bfile PREFIX`)

Three companion files must share the same PREFIX.

### PREFIX.bed — Binary Genotype Matrix

| Property | Value |
|---|---|
| Encoding | Binary, 2-bit per sample, SNP-major |
| Magic bytes | `0x6C 0x1B 0x01` |

Standard PLINK 1 binary BED format.  Each marker is stored as a
packed array of 2-bit codes:

| Code | Meaning |
|------|---------|
| `00` | Homozygous A1 |
| `01` | Missing |
| `10` | Heterozygous |
| `11` | Homozygous A2 |

Rows are padded to full bytes.

### PREFIX.bim — Variant Information

Whitespace-delimited, one line per variant.

| Column | Name | Description |
|--------|------|-------------|
| 1 | CHROM | Chromosome code |
| 2 | ID | Variant identifier |
| 3 | CM | Genetic distance (centiMorgans), ignored |
| 4 | POS | Base-pair position |
| 5 | A1 | Allele 1 (usually minor/alt) |
| 6 | A2 | Allele 2 (usually major/ref) |

### PREFIX.fam — Sample Information

Whitespace-delimited, one line per sample.

| Column | Name | Description |
|--------|------|-------------|
| 1 | FID | Family ID |
| 2 | IID | Individual ID (primary key for subject matching) |
| 3 | PAT | Paternal ID (ignored) |
| 4 | MAT | Maternal ID (ignored) |
| 5 | SEX | Sex code (ignored) |
| 6 | PHENO | Phenotype (ignored) |

The `.fam` file defines the **canonical subject order**.  All other
per-subject files are re-aligned to this order during loading.

---

## Null Model Residuals (`--null-resid FILE`)

Accepted formats, in detection order:

### Format A — Header with IID column

```
#IID    Residual
Subj-1  -0.006972
Subj-2   0.034234
```

- Lines starting with `##` are skipped (plink2-style comments).
- The first non-`##` line starting with `#FID`, `FID`, `#IID`, or `IID`
  is treated as a header.
- The `IID` (or `#IID`) column identifies subjects.
- `#FID` / `FID` columns are skipped; all other columns are value columns.
- **Header names (apart from IID/FID) are recorded but do not affect
  parsing** — the parser uses column order, not column names, to
  determine which values are residuals/weights/indicators.

### Format B — Pure numeric matrix (no IID, no header)

```
-0.006972    0.890109    0
 0.034234    0.890109    0
```

- First token of the first data line is numeric → pure numeric mode.
- All columns are value columns; no IID column exists.
- **Row count must equal the number of subjects in `.fam`.**
  Rows are assumed to be in `.fam` order.
- If no header is found and the first token is **not** numeric, an error
  is raised.  Add a header line with an `IID` column to fix this.

### Column layout by method

| Method | Columns after IID |
|--------|------------------|
| SPACox, SPAGRM, SPAmix, SPAmixPlus | `Residual` (1 or more; multi-column → multi-GWAS) |
| WtCoxG, LEAF | `Residual Weight Indicator` (exactly 3) |
| SPAsqr | `R_tau1 R_tau2 ... R_tauK` (K quantile residual columns) |

When multiple residual columns are present (SPACox/SPAGRM/SPAmix/SPAmixPlus),
each column is run as a separate GWAS.  Single-column output goes to `--out`
directly; multi-column output produces `PREFIX.1.gz`, `PREFIX.2.gz`, etc.

### Subject matching

- Format A: subjects are intersected with `.fam` by IID.
  Only subjects present in **all** loaded files are retained.
- Format B: rows are assumed to be in `.fam` order; no intersection
  is performed.

---

## Covariate File (`--covar FILE`)

Accepted formats:

### Format A — plink2 .cov header

```
#IID    AGE     GENDER
Subj-1  59.61   0
Subj-2  61.32   1
```

- A header line is required (first non-blank line).
- The `IID` (or `#IID`) column is identified for subject matching.
- Automatically skipped columns: `#FID`, `FID`, `SID`, `PAT`, `MAT`, `SEX`,
  `PHENO*` (any column whose name starts with "PHENO").
- **Column names of remaining columns do not matter** — all non-skipped,
  non-IID columns are used as covariates, identified by position.
- An intercept column is added automatically by the method code.

### Format B — Pure numeric matrix (no IID, no header)

```
59.61   0
61.32   1
```

- When no `IID` or `#IID` column is found in the first line, and all
  tokens are numeric, the file is treated as a pure numeric matrix.
- **Row count must equal the number of subjects in `.fam`.**
  Rows are assumed to be in `.fam` order.

---

## Eigenvector File (`--eigenvec FILE`)

Accepted formats:

### Format A — plink2 .eigenvec header

```
#IID    PC1     PC2
Subj-1  -0.8293 -0.2644
Subj-2   0.5602  0.2659
```

- A header line is required (first non-blank line).
- The `IID` (or `#IID`) column is identified for subject matching.
- Skipped columns: `#FID`, `FID`, `SID`.
- **PC columns are identified by column name**: any header matching
  `^[Pp][Cc]` (e.g., `PC1`, `pc2`, `PC10`) is treated as a PC column.
  Unlike `--null-resid` and `--covar`, non-matching column names are
  silently ignored, so column names *do* matter here.

Produced by `plink2 --pca`

### Format B — Pure numeric matrix (no IID, no header)

```
-0.8293  -0.2644
 0.5602   0.2659
```

- When no `IID` or `#IID` column is found, and all tokens are numeric,
  the file is treated as a pure numeric matrix.
- **Row count must equal the number of subjects in `.fam`.**
  Rows are assumed to be in `.fam` order.
- All columns are treated as PCs.

---

## Sparse GRM (`--sp-grm-grab FILE` or `--sp-grm-plink2 PREFIX`)

Exactly one of the two formats must be provided.  They are mutually
exclusive.

### `--sp-grm-grab FILE` — GRAB format

Three-column whitespace-delimited text file.

```
#IID1   IID2    COEF
f1_1    f1_1    0.955062
f1_2    f1_2    1.027229
f1_1    f1_2    0.482510
```

| Column | Name | Description |
|--------|------|-------------|
| 1 | ID1 | Subject ID (matches `.fam` IID) |
| 2 | ID2 | Subject ID (matches `.fam` IID) |
| 3 | VALUE | GRM coefficient |

- Lines starting with `#` are skipped.
- Subjects not present in `.fam` are silently dropped.
- Diagonal entries (`ID1 == ID2`) and off-diagonal entries are both accepted.
- For off-diagonal entries, the symmetric counterpart `(ID2, ID1)` is
  added automatically when needed.

This format is produced by the GRAB R package
(`GRAB.SparseGRM()`).

### `--sp-grm-plink2 PREFIX` — plink2 / GCTA format

Two companion files:

#### PREFIX.grm.id — Subject Index (optional)

```
f1_1    f1_1
f1_2    f1_2
f1_3    f1_3
```

| Column | Name | Description |
|--------|------|-------------|
| 1 | FID | Family ID (ignored) |
| 2 | IID | Individual ID (maps 0-based line index to IID) |

If `PREFIX.grm.id` **does not exist**, the 0-based indices in `.grm.sp`
are assumed to correspond to `.fam` subject order.  This requires that the
maximum index in `.grm.sp` is less than the number of subjects in `.fam`.

#### PREFIX.grm.sp — Sparse GRM Entries

```
0       0       0.92288343
1       1       1.04717390
0       1       0.48251000
```

| Column | Name | Description |
|--------|------|-------------|
| 1 | idx1 | 0-based row index |
| 2 | idx2 | 0-based column index |
| 3 | value | GRM coefficient |

- Indices reference the line order in `.grm.id` (or `.fam` when
  `.grm.id` is absent).
- Lower-triangular + diagonal entries are expected.

Produced by `plink2 --make-grm-sparse`.

---

## Reference Allele Frequency (`--ref-af FILE`)

### Format A — plink2 `.afreq` header (default)

Produced by `plink2 --freq`.

```
#CHROM  ID      REF     ALT     ALT_FREQS       OBS_CT
1       SNP_1   A       G       0.251           1000
```

| Column | Name | Description |
|--------|------|-------------|
| 1 | #CHROM | Chromosome code |
| 2 | ID | Variant identifier |
| 3 | REF | Reference allele |
| 4 | ALT | Alternate allele |
| 5 | ALT_FREQS | Alternate allele frequency in the reference population |
| 6 | OBS_CT | Total allele number in the reference population (sample size = OBS_CT / 2) |

### Variant matching (Format A)

Markers are matched to `.bim` by **(CHROM, ID)** — both must agree.
After a CHROM+ID match the alleles are checked:

- **Same orientation**: `.afreq` ALT == `.bim` col 5, `.afreq` REF == `.bim` col 6
  → `AF_ref = ALT_FREQS`.
- **Flipped orientation**: `.afreq` REF == `.bim` col 5, `.afreq` ALT == `.bim` col 6
  → `AF_ref = 1 − ALT_FREQS`.
- **No allele match** → marker is dropped.

Alleles are compared case-insensitively.

### Format B — two-column numeric fallback

A headerless file with exactly two numeric columns per line.
Rows are assumed to be in the same order as the `.bim` file, and the
row count must equal the number of `.bim` markers.  No variant-ID or
allele matching is performed — frequencies are taken as the
`.bim` col 5 allele frequency directly.

```
0.251   1000
0.490   998
```

| Column | Description |
|--------|-------------|
| 1 | `ALT_FREQS` — frequency of the `.bim` col 5 allele in the reference population |
| 2 | `OBS_CT` — total allele number in the reference population (sample size = OBS_CT / 2) |

### LEAF multi-population

For **LEAF**, provide one `--ref-af` file per **reference population**
and one `--null-resid` file per **ancestry cluster** as comma-separated
lists.  The number of reference populations and clusters are independent:

```
--null-resid cluster1_resid.txt,cluster2_resid.txt,cluster3_resid.txt
--ref-af     pop1.afreq,pop2.afreq
```

Summix uses the `nPop` reference AFs to estimate ancestry proportions
for each cluster; per-cluster analysis then proceeds as WtCoxG.
Either ref-af format (A or B) may be used.

---

## Pairwise IBD (`--pairwise-ibd FILE`)

Tab-separated output from `grab --cal-pairwise-ibd`.

```
#ID1    ID2     pa      pb      pc
f1_1    f1_2    0.002   0.480   0.518
```

| Column | Name | Description |
|--------|------|-------------|
| 1 | ID1 | Subject ID |
| 2 | ID2 | Subject ID |
| 3 | pa | P(IBD = 0) |
| 4 | pb | P(IBD = 1) |
| 5 | pc | P(IBD = 2) |

- The header line starts with `#` and is skipped by the reader.
- Lines starting with `#` are ignored, so headerless files from
  older versions are also accepted.

Required by SPAGRM.

---

## Individual AF Coefficients (`--ind-af-coef FILE`)

Pre-computed allele-frequency model produced by
`grab --cal-ind-af-coef`.

Markers are **not** matched by ID.  Both text and binary formats assume
the same marker order as the `.bim` file used when the model was
computed.  The marker count must match exactly.

### Text/gzip format (`.txt` / `.gz`)

Tab-separated with header:

```
#STATUS BETA0   BETA1   BETA2
0       0.25    0       0
1       0.25    -0.03   0.01
2       -1.10   0       0.08
```

| Column | Name | Description |
|--------|------|-------------|
| 1 | STATUS | Model fit status (integer, see below) |
| 2+ | BETA0 ... BETAp | Regression coefficients (intercept + PC betas) |

STATUS values:

| Value | Meaning |
|-------|---------|
| 0 | Uniform — MAC too low, use population allele frequency |
| 1 | OLS — linear regression fit used |
| 2 | Logistic — logistic regression on significant PCs |

- Lines starting with `#` are treated as header/comments and skipped.
- One data row per filtered marker, in `.bim` order.

### Binary format (`.bin`)

Random-access binary file.  No header.  Records are indexed by raw
`.bim` line number (genoIndex), so the file is pre-allocated to
`nBimMarkers * recordSize` bytes.

**Record layout** (one per `.bim` marker):

| Offset | Type | Description |
|--------|------|-------------|
| 0 | `int8_t` | STATUS (0, 1, or 2) |
| 1 | `double[1+nPC]` | BETA0, BETA1, ..., BETAp |

`recordSize = sizeof(int8_t) + (1 + nPC) * sizeof(double)`

- Unprocessed markers (e.g. filtered out) are all-zero (status 0, betas 0).
- The number of PCs (`nPC`) must be known at read time (same value used
  when the file was written).
- Byte order is machine-native (not portable across architectures).

---

## extract_tracts Text Files (`--admix-text-prefix PREFIX`)

Input for `grab --make-abed --admix-text-prefix PREFIX`, produced by
`extract_tracts_fast_pgzip.py`.

### File naming

For each ancestry `k` (0-based, auto-detected by scanning from `k=0`):

| File | Description |
|------|-------------|
| `{PREFIX}.anc{k}.dosage[.gz]` | Dosage values (0/1/2 per sample) for ancestry k |
| `{PREFIX}.anc{k}.hapcount[.gz]` | Haplotype count values (0/1/2 per sample) for ancestry k |

Both plain text and gzip-compressed files are accepted transparently.
K is determined by scanning from `k=0` until no dosage file is found.

### File format

Tab-delimited text with a header row:

```
CHROM   POS     ID      REF     ALT     SAMPLE1 SAMPLE2 ...
chr1    100000  rs1234  A       G       0       1       ...
chr1    200000  rs5678  C       T       2       0       ...
```

| Column | Name | Description |
|--------|------|-------------|
| 1 | CHROM | Chromosome |
| 2 | POS | Base-pair position |
| 3 | ID | Variant identifier |
| 4 | REF | Reference allele |
| 5 | ALT | Alternate allele |
| 6+ | SAMPLE_IDs | Integer value per sample (0, 1, or 2) |

- The header of `{PREFIX}.anc0.dosage[.gz]` defines the `.bim` and `.fam` of the output.
- All ancestry files must share the same marker and sample order.
- Missing values (`NA`, `.`, `-9`) are encoded as missing in the `.abed` file.

### Output

`grab --make-abed --admix-text-prefix PREFIX --out-prefix OUT` writes:

| File | Contents |
|------|----------|
| `OUT.abed` | Binary ancestry track matrix (2K tracks per marker: dosage + hapcount per ancestry) |
| `OUT.bim` | Variant information (CHROM ID 0 POS REF ALT) from the first dosage file |
| `OUT.fam` | Sample list from the first dosage file header (FID=IID) |

Pass `OUT` as `--admix-bfile` to `SPAmixLocalPlus` or `--cal-admix-phi`.

---

## Summary of Subject Matching Rules

| File | IID source | Fallback | Parser | Source file |
|------|-----------|----------|--------|-------------|
| `.fam` | Column 2 (IID) | — | `parseFamIIDs()` | `io/subject_data.cpp` |
| `--null-resid` | Header `IID` / `#IID` column | Pure numeric: `.fam` order (row count must match) | `SubjectData::parseIIDFile()` | `io/subject_data.cpp` |
| `--covar` | Header `IID` / `#IID` column | Pure numeric: `.fam` order (row count must match) | `SubjectData::loadCovar()` | `io/subject_data.cpp` |
| `--eigenvec` | Header `IID` / `#IID` column | Pure numeric: `.fam` order (row count must match) | `SubjectData::loadEigenVecs()` | `io/subject_data.cpp` |
| `--sp-grm-grab` | ID1 / ID2 columns | — | `SparseGRM()` constructor | `io/sparse_grm.cpp` |
| `--sp-grm-plink2` | `.grm.id` file | `.fam` order when `.grm.id` absent | `SparseGRM::fromGCTA()` | `io/sparse_grm.cpp` |
| `--pairwise-ibd` | ID1 / ID2 columns | — | `loadIBD()` | `spagrm/geno_prob.cpp` |
| `--ref-af` | CHROM + ID matching, allele flip | Two-column numeric: `.bim` order (row count must match) | `loadRefAfFile()` + `matchMarkers()` / `matchMarkersNumeric()` | `wtcoxg/wtcoxg.cpp` |
| `--ind-af-coef` | Same marker order as `.bim` | — | `loadAFModels()` / `IndivAFReader` | `spamix/indiv_af.cpp` |

When subject IIDs are provided, the analysis uses the **intersection**
of all loaded files and `.fam`.  When a file uses `.fam`-order fallback
(pure numeric matrix or missing `.grm.id`), all `.fam` subjects are
implicitly included from that file.
