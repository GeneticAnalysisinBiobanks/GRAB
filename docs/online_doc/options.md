# GRAB — Complete command-line option reference

This document enumerates every command-line option accepted by the `grab2`
binary, including the analysis methods it dispatches and the utility modes
shared across the pipeline.  It is intended as an authoritative companion
to the `grab2 --help` output: each option entry lists the metavar, default
value (when one exists), the set of methods that consume the option, and
the validation rules that the dispatcher applies.

The conventions used below are:

- `PREFIX` — a file path prefix; companion files are derived from it.
- `FILE` — a single concrete file path.
- `INT`, `FLOAT`, `LIST`, `COL_IDS`, `NUMS`, `NAME`, `MODE` — typed values.
- A flag without a metavar is a parameterless boolean switch.
- Defaults given as `0 = library default` indicate that the value `0` is a
  sentinel forwarded to the underlying library to request its built-in
  default.

The dispatcher rejects any flag that appears more than once on a single
command line.  An option that is consumed by a method but supplied without
a value triggers a "missing value" error.  Numeric values that are not
parseable as the expected type trigger a "requires a numeric value" error
with the offending token echoed.

## 1. Analysis modes

`grab2` operates in one of two top-level dispatch modes:

1. **Method mode** — `--method NAME` selects a GWAS method.  The
   case of `NAME` is normalized to the canonical spelling (`spacox`,
   `SPACox`, and `SPACOX` are all accepted).  Supported names:

   | Method            | Purpose                                                            |
   | ----------------- | ------------------------------------------------------------------ |
   | `SPACox`          | Saddlepoint Approximation for Cox proportional hazards             |
   | `SPAGRM`          | SPA with sparse GRM relatedness correction                         |
   | `SAGELD`          | G x E interaction for longitudinal data with GRM correction        |
   | `SPAmix`          | SPA with individual ancestry-based allele frequencies              |
   | `SPAmixPlus`      | SPAmix with additional sparse GRM relatedness correction           |
   | `SPAmixLocalPlus` | Local-ancestry-specific GWAS with SPA + phi-based variance         |
   | `SPAsqr`          | SPA for quantile regression residuals (multiple tau levels)        |
   | `WtCoxG`          | Weighted Cox regression for time-to-event GWAS                     |
   | `LEAF`            | Local Ethnicity-Aware GWAS — WtCoxG per ancestry cluster           |

2. **Utility mode** — one of the parameterless switches below replaces
   `--method` and selects a pipeline pre-/post-processing step:

   | Switch                | Output                                                            |
   | --------------------- | ----------------------------------------------------------------- |
   | `--cal-af-coef`       | Per-marker individual-ancestry AF model (`PREFIX.afc[.gz\|.zst]`) |
   | `--cal-pairwise-ibd`  | Pairwise IBD probabilities (`PREFIX.ibd[.gz\|.zst]`)              |
   | `--cal-phi`           | Phi kinship matrices (`PREFIX.phi[.gz\|.zst]`)                    |
   | `--make-abed`         | Admixed ancestry binary (`PREFIX.abed`, `PREFIX.bim`, `PREFIX.fam`) |
   | `--int-pheno`         | Inverse-normal-transformed phenotype file (`PREFIX.txt`)          |

   Utility switches are mutually exclusive with `--method` and with each
   other.  Where a utility mode shares an option with the GWAS methods,
   the semantics are identical.

The default option `grab2 --help` prints a short usage summary together
with the list of method and utility topics.  Detailed per-method or
per-flag help is available via `grab2 --help <topic>`, where `<topic>` is
a method name, a utility-mode name without the `--` prefix, the literal
`options`, or one of the file-flag topics
(`pheno`, `covar`, `ref-af`, `sp-grm`, `pairwise-ibd`, `ind-af-coef`,
`admix-phi`).

## 2. Genotype input

Every analysis mode that reads genotypes accepts exactly one of the
following:

```
--bfile PREFIX
--pfile PREFIX
--vcf FILE
--bcf FILE
--bgen FILE <REF/ALT mode>
```

| Flag                           | Metavar                | Description                                                                                                                                                                       |
| ------------------------------ | ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--bfile`                      | `PREFIX`               | PLINK binary genotype prefix (reads `PREFIX.bed`, `PREFIX.bim`, `PREFIX.fam`).                                                                                                    |
| `--pfile`                      | `PREFIX`               | PLINK2 genotype prefix (reads `PREFIX.pgen`, `PREFIX.pvar`, `PREFIX.psam`).                                                                                                       |
| `--vcf`                        | `FILE`                 | VCF (Variant Call Format) input, plain text or BGZF-compressed (`.vcf.gz`).  Rejects BCF2 content with a redirect to `--bcf`.                                                     |
| `--bcf`                        | `FILE`                 | BCF2 (binary VCF).  Rejects VCF text with a redirect to `--vcf`.                                                                                                                  |
| `--bgen`                       | `FILE <REF/ALT mode>`  | BGEN genotype file.  The trailing mode token is mandatory and must be one of `ref-first`, `ref-last`, `ref-unknown` (matches `plink2 --bgen`).                                    |
| `--admix-bfile`                | `PREFIX`               | Admixed ancestry binary genotype prefix (reads `PREFIX.abed`, `PREFIX.bim`, `PREFIX.fam`).  Used by `SPAmixLocalPlus` and by `--cal-phi`.                                          |

The `--bgen` REF/ALT mode semantics are:

- `ref-first`  — the first allele of each variant is REF (raw UK Biobank /
  IMPUTE / qctool convention).
- `ref-last`   — the last allele of each variant is REF (default for
  plink2-exported `.bgen`; `alleles[0]` is ALT).
- `ref-unknown` — REF is provisional; treated identically to `ref-last`
  with a `[WARN]` noting the provisional REF status.

`--bfile`, `--pfile`, `--vcf`, `--bcf`, and `--bgen` are mutually
exclusive.  `--admix-bfile` is consumed only by `SPAmixLocalPlus`,
`--cal-phi`, and `--make-abed`.

## 3. Subject input files

| Flag                 | Metavar    | Description                                                                                                                       |
| -------------------- | ---------- | --------------------------------------------------------------------------------------------------------------------------------- |
| `--pheno`            | `FILE`     | Phenotype file.  Strict format: mandatory header, first column = subject ID.  Column names must match `[0-9A-Za-z_\-.]+`.         |
| `--covar`            | `FILE`     | Covariate file.  Same strict format as `--pheno`.  When `--covar` is absent, `--covar-name` selects columns from `--pheno` instead. |
| `--covar-name`       | `COL_IDS`  | Comma-separated list of covariate column names.  An intercept is added automatically.                                              |
| `--pheno-name`       | `COL_IDS`  | Comma-separated list of phenotype column names selected for analysis.                                                              |
| `--resid-name`       | `COL_IDS`  | Comma-separated list of residual column names (pre-computed residuals).                                                            |
| `--pc-cols`          | `COL_IDS`  | Comma-separated list of principal-component column names.  Default `PC1,PC2,PC3,PC4`.  Does **not** enter the null-model design.   |
| `--keep`             | `FILE`     | Restrict analysis to subjects listed in file (one `FID IID` pair per line).                                                       |
| `--remove`           | `FILE`     | Exclude subjects listed in file (one `FID IID` pair per line).  Mutually inclusive with `--keep` is allowed; the final set is `keep \ remove`. |

`--pheno-name` and `--resid-name` are mutually exclusive on every method
that supports both.  When `--pheno-name` is supplied, the dispatcher
fits a null model internally and feeds the resulting residuals to the
score test; when `--resid-name` is supplied, the named columns of
`--pheno` are used verbatim as residuals.

## 4. Reference & relatedness files

| Flag                  | Metavar  | Description                                                                                                                                                                                |
| --------------------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--ref-af`            | `FILE`   | Reference allele frequency file.  Two accepted formats: (A) plink2 `.afreq` (`#CHROM  ID  REF  ALT  ALT_FREQS  OBS_CT`); (B) two-column numeric (`ALT_FREQS  OBS_CT`, rows in `.bim` order). For `LEAF`: comma-separated, one file per reference population. |
| `--sp-grm-grab`       | `FILE`   | Sparse GRM in GRAB format (whitespace-delimited; header `ID1  ID2  VALUE`).                                                                                                                |
| `--sp-grm-plink2`     | `FILE`   | Sparse GRM in plink2 `.grm.sp` format; companion `.grm.id` (`FID  IID`) auto-detected.                                                                                                     |
| `--ind-af-coef`       | `FILE`   | Pre-computed individual-ancestry AF model produced by `--cal-af-coef`.  Text or gz; format `#STATUS  BETA0  BETA1 ...`.                                                                    |
| `--pairwise-ibd`      | `FILE`   | Pairwise IBD probabilities produced by `--cal-pairwise-ibd`.  Tab-separated `ID1  ID2  pa  pb  pc`.                                                                                        |
| `--admix-phi`         | `FILE`   | Pre-computed phi file from `--cal-phi`.  Wide tab-separated `idx1  idx2  anc0_A  anc0_B  anc0_C  anc0_D  [anc1_A ...]`.                                                                    |
| `--rfmix-msp`         | `FILE`   | MSP local-ancestry file from rfmix2.  Used only with `--make-abed`.                                                                                                                        |
| `--admix-text-prefix` | `PREFIX` | extract_tracts text output prefix.  Reads `{PREFIX}.anc{k}.dosage[.gz]` and `{PREFIX}.anc{k}.hapcount[.gz]`.  Used only with `--make-abed`.                                                |

`--sp-grm-grab` and `--sp-grm-plink2` are mutually exclusive.

## 5. Filtering files

| Flag        | Metavar | Description                                                              |
| ----------- | ------- | ------------------------------------------------------------------------ |
| `--extract` | `FILE`  | Restrict analysis to SNPs listed in file (one ID per line; matches `.bim` column 2).      |
| `--exclude` | `FILE`  | Exclude SNPs listed in file (one ID per line).                                              |
| `--chr`     | `NUMS`  | Restrict analysis to specified chromosomes.  Comma-separated values or ranges (`1-4,6,22`). |

These three flags are propagated through the `GenoSpec` descriptor and
apply to every method that opens a genotype file.

## 6. Null-model fitting

The following flags drive the internal null-model engine that is invoked
whenever a method's residuals are computed from a phenotype column
(`--pheno-name`):

| Flag                  | Metavar  | Default  | Description                                                                                                                                                                                                                                                                |
| --------------------- | -------- | -------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--regression-model`  | `MODEL`  | `auto`   | Null-model regression family.  Accepted: `auto`, `linear`, `logistic`, `cox`, `ordinal`.  In `auto` mode, the family is inferred per phenotype token from the column values and from colon syntax (`TIME:EVENT`).                                                            |
| `--save-resid`        | (none)   | off      | Write fitted residuals to `PREFIX.null.resid` (re-loadable via `--pheno PREFIX.null.resid --resid-name ...` in a later invocation).  Requires `--pheno-name` and `--out`.                                                                                                   |
| `--sageld-x`          | `COL_IDS`| —        | Comma-separated environment column name(s) for SAGELD pheno mode.  Each name must also appear in `--covar-name`.                                                                                                                                                            |

Per-method whitelist of `--regression-model` values:

| Method        | Accepted values                  |
| ------------- | -------------------------------- |
| `SPACox`      | `auto`, `linear`, `logistic`, `cox`, `ordinal` |
| `SPAGRM`      | `auto`, `linear`, `logistic`, `cox`, `ordinal` |
| `SPAmix`      | `auto`, `linear`, `logistic`, `cox`, `ordinal` |
| `SPAmixPlus`  | `auto`, `linear`, `logistic`, `cox`, `ordinal` |
| `SPAmixLocalPlus` | `auto`, `linear`, `logistic`, `cox`, `ordinal` |
| `SAGELD`      | `auto` only (rejects any other value)          |
| `SPAsqr`      | `auto`, `linear` only                          |
| `WtCoxG`      | `auto`, `logistic`, `cox` only                 |
| `LEAF`        | `auto`, `logistic`, `cox` only                 |

For `cox`, `--pheno-name` must use the `TIME:EVENT` syntax
(e.g. `SURVTIME1:STATUS1,SURVTIME2:STATUS2`).  For `ordinal`, phenotype
columns must be integer-coded levels `0, 1, ..., J-1`.

## 7. Output

| Flag                  | Metavar    | Default | Description                                                                                                       |
| --------------------- | ---------- | ------- | ----------------------------------------------------------------------------------------------------------------- |
| `--out`               | `PREFIX`   | —       | Output prefix; method-specific suffixes are appended automatically (see per-method sections below).               |
| `--compression`       | `gz\|zst`  | none    | Compress output files with gzip or Zstandard.  Default produces plain text.                                       |
| `--compression-level` | `INT`      | `0`     | Compression level.  Gzip accepts `1..9`; Zstandard accepts `1..22`.  `0` requests the library default.            |

## 8. Marker quality-control thresholds

| Flag                          | Metavar | Default | Description                                                                                       |
| ----------------------------- | ------- | ------- | ------------------------------------------------------------------------------------------------- |
| `--geno`                      | `FLOAT` | `0.1`   | Per-marker missing-rate cutoff; markers with `MISS_RATE > geno` are excluded.                     |
| `--maf`                       | `FLOAT` | `1e-5`  | Minimum minor allele frequency.                                                                   |
| `--mac`                       | `FLOAT` | `10`    | Minimum minor allele count.                                                                       |
| `--hwe`                       | `FLOAT` | `0`     | Exclude markers with HWE p-value below the threshold.  `0` disables the filter.                   |
| `--min-maf-ibd`               | `FLOAT` | `0.01`  | Minimum MAF for the IBD-calculation marker set.  Used only by `--cal-pairwise-ibd`.               |

## 9. Statistical thresholds and tuning

| Flag                              | Metavar   | Default | Used by                                                                  | Description                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| --------------------------------- | --------- | ------- | ------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--prevalence`                    | `FLOAT`   | —       | `WtCoxG`, `LEAF`                                                         | Disease prevalence in the reference population.  Required; must be positive.                                                                                                                                                                                                                                                                                                                                                                |
| `--batch-effect-p-threshold`      | `FLOAT`   | `0.05`  | `WtCoxG`, `LEAF`                                                         | Batch-effect p-value cutoff.                                                                                                                                                                                                                                                                                                                                                                                                                |
| `--covar-p-threshold`             | `FLOAT`   | `5e-5`  | `SPACox`                                                                  | Covariate p-value cutoff used to admit covariate effects into the null model.                                                                                                                                                                                                                                                                                                                                                                |
| `--spa-z-threshold`               | `FLOAT`   | `2.0`   | all SPA methods                                                          | Magnitude of the score-test z-statistic above which the saddlepoint approximation is invoked instead of the Gaussian tail.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| `--outlier-iqr-multiplier`        | `FLOAT`   | `1.5`   | `SPAmix`, `SPAmixPlus`, `SPAmixLocalPlus`, `WtCoxG`, `LEAF`, `SPAsqr`, `SPAGRM` | Tukey IQR multiplier for residual outlier detection.  A residual is flagged as an outlier iff `r < q25 − k·IQR` or `r > q75 + k·IQR`.  Larger `k` → narrower outlier set; `k ≤ −0.5` forces every residual through the exact CGF.                                                                                                                                                                                                                                                                                                                                                          |
| `--spasqr-outlier-abs-bound`      | `FLOAT`   | `0.55`  | `SPAsqr`                                                                  | SPAsqr-only absolute outlier cutoff.                                                                                                                                                                                                                                                                                                                                                                                                       |
| `--spagrm-control-outlier`        | (none)    | off     | `SPAGRM`                                                                  | Enable iterative IQR-ratio adjustment so the outlier share stays in `(0, 5%]`.  When enabled, the IQR ratio is shrunk by `0.8` until at least one outlier appears, and increased by `0.5` until the share is at most `5%`.                                                                                                                                                                                                                                                                                                                                                       |

## 10. SPAsqr-specific knobs

| Flag                  | Metavar  | Default                     | Description                                                                                                                                                                                                                                                                                                                                                                                                                                  |
| --------------------- | -------- | --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--spasqr-taus`       | `LIST`   | `0.1,0.3,0.5,0.7,0.9`       | Comma-separated tau levels for SPAsqr.  At most 20 entries; each value must lie in `(0, 1)`.                                                                                                                                                                                                                                                                                                                                                |
| `--spasqr-tol`        | `FLOAT`  | `1e-7`                      | Convergence tolerance for the SPAsqr null-model SQR solver.  QMME tightens this internally to `min(--spasqr-tol, 1e-9)`.                                                                                                                                                                                                                                                                                                                    |
| `--spasqr-h`          | `FLOAT`  | —                           | Explicit bandwidth for the SQR null model.  Mutually exclusive with `--spasqr-h-scale`.                                                                                                                                                                                                                                                                                                                                                      |
| `--spasqr-h-scale`    | `FLOAT`  | `3` (score) / `10` (wald)   | Divisor for IQR-based bandwidth: `h = IQR(Y) / SCALE`.  Mutually exclusive with `--spasqr-h`.                                                                                                                                                                                                                                                                                                                                                  |
| `--spasqr-solver`     | `NAME`   | `qmme`                      | Smoothed quantile regression solver.  `qmme` — Quadratic Majorization-Minimization with Extrapolation; caches the Cholesky of the Hessian upper bound once per `(phenotype, bandwidth)` and reuses it across all τ levels.  `conquer` — Convolution-type smoothed QR; Huber initialization followed by BB gradient descent, refit per τ.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| `--spasqr-mode`       | `MODE`   | `score`                     | Inference mode.  `score` — score test on null-model residuals; output columns `... P_CCT  P_tau{val}...  Z_tau{val}...`.  `wald` — per-marker full-model Wald test; output columns `... TAU  BETA  SE  Z  P` (one row per marker, τ).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |
| `--pred-list`         | `FILE`   | —                           | `pred.list` file for LOCO analysis.  Two whitespace-separated columns: phenotype name and the path to the corresponding `.loco` file.  Format of the `.loco` file is auto-detected per phenotype (Regenie chromosome-major or LDAK-KVIK subject-major).                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                |
| `--pheno-transform`   | `MODE`   | `int`                       | Phenotype pre-transform for SPAsqr.  `raw` — no transform.  `int` — inverse normal transform (Blom plotting position, average-rank ties) applied per phenotype on its non-missing scope.  `standardize` — phenotype centered and scaled to unit variance per phenotype.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          |

## 11. LEAF clustering knobs

| Flag                      | Metavar | Default                       | Description                                                                                                                            |
| ------------------------- | ------- | ----------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `--leaf-nclusters`        | `INT`   | from `--ref-af` file count    | Number of K-means clusters for LEAF.                                                                                                  |
| `--leaf-cluster-file`     | `FILE`  | —                             | Pre-computed cluster labels.  Two-column file: `IID` (or `#IID`) and `cluster` (integers `1..K`).  Disables internal K-means.          |
| `--leaf-kmeans-nstart`    | `INT`   | `25`                          | Number of K-means++ random restarts.  Ignored when `--leaf-cluster-file` is provided.                                                  |

## 12. Engine and runtime controls

| Flag                | Metavar | Default | Description                                                                                                                                                                                                                                                                                                                                                              |
| ------------------- | ------- | ------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--threads`         | `INT`   | `1`     | Number of worker threads for the chunk-level marker engine and for the pre-marker work-stealing pools (null-model fits, Chow-Liu MAF bins, K-means restarts, etc.).  The main thread is not counted in N.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `--chunk-size`      | `INT`   | `8192`  | Markers per chunk.  Minimum permitted value is `256`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |
| `--seed`            | `INT`   | `0`     | Random seed for reproducibility.  `0` uses `std::random_device`.  Affects the K-means initialization in LEAF and the ordinal null-model fit in `SPACox`, `SPAGRM`, `SPAmix`, `SPAmixPlus`, and `LEAF`.                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       |

## 13. Per-method option summary

The following tables enumerate every option accepted by each analysis
mode, classified as required or optional.  An option marked as required
must be supplied; an option marked as optional has a sensible default
(see sections 6–12) and may be omitted.  Combined "exactly-one" entries
indicate that any one — but only one — of the listed flags must be
present.

### 13.1 `SPACox`

**Required**: genotype input (one of `--bfile`, `--pfile`, `--vcf`,
`--bcf`, `--bgen`), `--pheno`, `--out`.

**Optional**: `--covar`, `--covar-name`, `--resid-name`, `--pheno-name`,
`--regression-model`, `--save-resid`, `--covar-p-threshold`,
`--spa-z-threshold`, `--seed`, `--threads`, `--chunk-size`,
`--compression`, `--compression-level`, `--keep`, `--remove`,
`--extract`, `--exclude`, `--geno`, `--maf`, `--mac`, `--hwe`, `--chr`.

Either `--pheno-name` (fit a null model internally) or `--resid-name`
(supply pre-computed residuals) must be present.

Output: `PREFIX.PHENO.SPACox[.gz|.zst]`
with columns
`CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P  Z`.

### 13.2 `SPAGRM`

**Required**: genotype input, `--pheno`, `--out`, sparse GRM (one of
`--sp-grm-grab`, `--sp-grm-plink2`), `--pairwise-ibd`.

**Optional**: `--resid-name`, `--pheno-name`, `--regression-model`,
`--save-resid`, `--covar`, `--covar-name`, `--spa-z-threshold`,
`--outlier-iqr-multiplier`, `--spagrm-control-outlier`, `--seed`,
`--threads`, `--chunk-size`, `--compression`, `--compression-level`,
`--keep`, `--remove`, `--extract`, `--exclude`, `--geno`, `--maf`,
`--mac`, `--hwe`, `--chr`.

The pairwise IBD file is produced once with
`grab2 --cal-pairwise-ibd`.

Output: `PREFIX.PHENO.SPAGRM[.gz|.zst]`
with columns
`CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P  Z`.

### 13.3 `SAGELD`

**Required**: genotype input, `--pheno`, `--out`, sparse GRM,
`--pairwise-ibd`.

**Optional**: `--resid-name`, `--pheno-name`, `--covar-name`,
`--sageld-x`, `--save-resid`, `--spa-z-threshold`, `--threads`,
`--chunk-size`, `--compression`, `--compression-level`, `--keep`,
`--remove`, `--extract`, `--exclude`, `--geno`, `--maf`, `--mac`,
`--hwe`, `--chr`.

Two input modes are supported:

- **Residual mode** — supply `lme4::lmer()` residuals via
  `--resid-name R_G,R_<E1>,R_Gx<E1>[,...]`.
- **Pheno mode** — supply long-format Y, X, E via `--pheno-name` and
  `--covar-name`, and name the environment columns via `--sageld-x`.
  `--covar-name` must include every variable in `--sageld-x`.

Output: `PREFIX.SAGELD` (residual mode, single file) or
`PREFIX.<phenoName>.SAGELD` (pheno mode, one file per phenotype)
with columns
`CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
 P_G  P_Gx<E1>  [...]  Z_G  Z_Gx<E1>  [...]`.

### 13.4 `SPAmix`

**Required**: genotype input, `--pheno`, `--out`, `--pc-cols`.

**Optional**: `--covar`, `--covar-name`, `--resid-name`, `--pheno-name`,
`--regression-model`, `--save-resid`, `--ind-af-coef`,
sparse GRM (one of `--sp-grm-grab`, `--sp-grm-plink2`),
`--outlier-iqr-multiplier`, `--spa-z-threshold`, `--seed`,
`--threads`, `--chunk-size`, `--compression`, `--compression-level`,
`--keep`, `--remove`, `--extract`, `--exclude`, `--geno`, `--maf`,
`--mac`, `--hwe`, `--chr`.

`--pc-cols` is **not** automatically added to the null-model design;
list any PCs you want adjusted in `--covar-name` explicitly.  The
optional sparse GRM enables GRM-based variance (same machinery as
`SPAmixPlus`).  A per-marker AF model can be pre-computed for speed
with `grab2 --cal-af-coef`.

Output: `PREFIX.PHENO.SPAmix[.gz|.zst]`
with columns
`CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P  Z  BETA  SE`.

### 13.5 `SPAmixPlus`

**Required**: genotype input, `--pheno`, `--out`, `--pc-cols`,
sparse GRM.

**Optional**: `--covar`, `--covar-name`, `--resid-name`, `--pheno-name`,
`--regression-model`, `--save-resid`, `--ind-af-coef`,
`--outlier-iqr-multiplier`, `--spa-z-threshold`, `--seed`,
`--threads`, `--chunk-size`, `--compression`, `--compression-level`,
`--keep`, `--remove`, `--extract`, `--exclude`, `--geno`, `--maf`,
`--mac`, `--hwe`, `--chr`.

Output: `PREFIX.PHENO.SPAmixP[.gz|.zst]`
with columns
`CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P  P  Z  BETA  SE`.

### 13.6 `SPAmixLocalPlus`

**Required**: `--admix-bfile`, `--pheno`, `--admix-phi`, `--out`.

**Optional**: `--covar`, `--covar-name`, `--resid-name`, `--pheno-name`,
`--regression-model`, `--save-resid`, `--keep`, `--remove`,
`--extract`, `--exclude`, `--outlier-iqr-multiplier`,
`--spa-z-threshold`, `--threads`, `--chunk-size`, `--compression`,
`--compression-level`, `--geno`, `--maf`, `--mac`, `--hwe`, `--chr`.

Two-phase workflow:

1. `grab2 --cal-phi --admix-bfile PREFIX --sp-grm-plink2 FILE --out OUTPUT_PREFIX`
2. `grab2 --method SPAmixLocalPlus --admix-bfile PREFIX --admix-phi OUTPUT_PREFIX.phi --pheno FILE --out PREFIX`

Output: `PREFIX.PHENO.LocalP[.gz|.zst]`
with columns
`CHROM  POS  ID  REF  ALT  P_CCT
 anc0_AltFreq  anc0_MissingRate  anc0_P  anc0_Pnorm  anc0_Stat  anc0_Var  anc0_zScore  anc0_AltCounts  anc0_BetaG
 anc1_AltFreq  ...  (repeated for each ancestry k)`.

### 13.7 `SPAsqr`

**Required**: genotype input, `--out`, sparse GRM.

**Optional**: `--pheno`, `--covar`, `--covar-name`, `--resid-name`,
`--pheno-name`, `--spasqr-taus`, `--spasqr-tol`, `--spasqr-h`,
`--spasqr-h-scale`, `--outlier-iqr-multiplier`,
`--spasqr-outlier-abs-bound`, `--spa-z-threshold`, `--threads`,
`--chunk-size`, `--compression`, `--compression-level`, `--keep`,
`--remove`, `--extract`, `--exclude`, `--geno`, `--maf`, `--mac`,
`--hwe`, `--chr`, `--pred-list`, `--pheno-transform`,
`--spasqr-solver`, `--spasqr-mode`.

Output columns:

- `--spasqr-mode score` (default):
  `CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
   P_CCT  P_tau{val}...  Z_tau{val}...`.
- `--spasqr-mode wald`:
  `CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
   TAU  BETA  SE  Z  P` — one row per `(marker, τ)`.

### 13.8 `WtCoxG`

**Required**: genotype input, `--out`, `--ref-af`, `--prevalence`.

**Optional**: `--pheno`, `--covar`, `--covar-name`, `--resid-name`,
`--pheno-name`, `--regression-model`, sparse GRM,
`--batch-effect-p-threshold`, `--spa-z-threshold`,
`--outlier-iqr-multiplier`, `--threads`, `--chunk-size`,
`--compression`, `--compression-level`, `--keep`, `--remove`,
`--extract`, `--exclude`, `--geno`, `--maf`, `--mac`, `--hwe`, `--chr`.

`--pheno-name` accepts either a binary column or a `TIME:EVENT` pair.
`--regression-model` accepts `auto`, `logistic`, `cox`.

Output:
`CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
 p_ext  p_noext  z_ext  z_noext  p_batch`.

### 13.9 `LEAF`

**Required**: genotype input, `--out`, `--ref-af`, `--prevalence`.

**Optional**: `--pheno`, `--covar`, `--covar-name`, `--resid-name`,
`--pheno-name`, `--regression-model`, `--pc-cols`, `--leaf-nclusters`,
`--leaf-cluster-file`, `--leaf-kmeans-nstart`, `--seed`, sparse GRM,
`--batch-effect-p-threshold`, `--spa-z-threshold`,
`--outlier-iqr-multiplier`, `--threads`, `--chunk-size`,
`--compression`, `--compression-level`, `--keep`, `--remove`,
`--extract`, `--exclude`, `--geno`, `--maf`, `--mac`, `--hwe`, `--chr`.

`--ref-af` accepts a comma-separated list, one file per reference
population.  The number of clusters and the number of reference
populations may differ; Summix estimates per-cluster ancestry
proportions from the reference set.

Output:
`CHROM  POS  ID  REF  ALT  MISS_RATE  ALT_FREQ  MAC  HWE_P
 meta.p_ext  meta.p_noext  cl1.p_ext  cl1.p_noext  cl1.p_batch  [cl2 ...]`.

## 14. Utility modes

### 14.1 `--cal-af-coef`

Compute the per-marker individual-ancestry AF model.

**Required**: genotype input, `--pc-cols`, `--out`.

**Optional**: `--pheno`, `--covar`, `--keep`, `--remove`, `--extract`,
`--exclude`, `--compression`, `--compression-level`, `--threads`,
`--chunk-size`, `--geno`, `--maf`, `--mac`, `--hwe`, `--chr`.

At least one of `--pheno` or `--covar` must be provided; the PC columns
named in `--pc-cols` are read from whichever file supplies them.

Output: `PREFIX.afc[.gz|.zst]`.  Pass to `--ind-af-coef` for
`SPAmix` / `SPAmixPlus`.

### 14.2 `--cal-pairwise-ibd`

Compute pairwise IBD probabilities from a sparse GRM plus genotypes.

**Required**: genotype input, `--out`, sparse GRM.

**Optional**: `--keep`, `--remove`, `--extract`, `--exclude`, `--chr`,
`--min-maf-ibd`, `--threads`, `--compression`, `--compression-level`.

Output: `PREFIX.ibd[.gz|.zst]` — tab-separated
`ID1  ID2  pa  pb  pc`.  Pass to `--pairwise-ibd` for `SPAGRM` /
`SAGELD`.

### 14.3 `--cal-phi`

Compute phi kinship matrices from admixed genotypes and a sparse GRM.

**Required**: `--admix-bfile`, sparse GRM, `--out`.

**Optional**: `--keep`, `--remove`, `--extract`, `--exclude`,
`--threads`, `--compression`, `--compression-level`.

Output: `PREFIX.phi[.gz|.zst]` — wide tab-separated
`idx1  idx2  anc0_A  anc0_B  anc0_C  anc0_D  [anc1_A ...]`.  Pass to
`--admix-phi` for `SPAmixLocalPlus`.

### 14.4 `--make-abed`

Build a `.abed` admixed-ancestry binary from a VCF/BCF + MSP pair or
from `extract_tracts` text output.

**Required**: `--out`.

**Optional**: `--vcf`, `--bcf`, `--rfmix-msp`, `--admix-text-prefix`,
`--keep`, `--remove`, `--threads`.

Modes (mutually exclusive):

- `--vcf FILE --rfmix-msp FILE` — phased VCF + rfmix2 MSP.
- `--bcf FILE --rfmix-msp FILE` — phased BCF2 + rfmix2 MSP.
- `--admix-text-prefix PREFIX`  — `extract_tracts` text output.

Output: `PREFIX.abed`, `PREFIX.bim`, `PREFIX.fam`.  Pass `PREFIX`
to `--admix-bfile` for `SPAmixLocalPlus` or `--cal-phi`.

### 14.5 `--int-pheno`

Apply the inverse normal transform to every trait column of a
phenotype file.

**Required**: `--pheno`, `--out`.

Output: `PREFIX.txt` — header `FID  IID  Y1  Y2  ...`, same row order
as the input.  Each `Y*` column is independently INT-transformed using
the Blom plotting position with average-rank ties on its non-missing
scope; missing entries (`NA`, `"."`, blank) remain missing in the
output.

## 15. Help-topic flags

The following help topics are accepted by `grab2 --help <topic>`:

- Method names: `SPACox`, `SPAGRM`, `SAGELD`, `SPAmix`, `SPAmixPlus`,
  `SPAmixLocalPlus`, `SPAsqr`, `WtCoxG`, `LEAF`.
- Utility names (without the `--` prefix): `cal-af-coef`,
  `cal-pairwise-ibd`, `cal-phi`, `make-abed`, `int-pheno`.
- The literal `options`, which prints every option in a single grouped
  listing.
- File-flag topics that expand the format documentation: `pheno`,
  `covar`, `ref-af`, `sp-grm`, `pairwise-ibd`, `ind-af-coef`,
  `admix-phi`.
- Any single-flag topic that names a documented flag (e.g.
  `--help outlier-iqr-multiplier`).
