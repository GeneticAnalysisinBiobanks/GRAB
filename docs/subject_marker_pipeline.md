# Subject & Marker Filtering Pipeline

GRAB uses a structured, multi-step pipeline to determine which subjects
and markers enter the marker-level association analysis.

---

## Subject Pipeline

Each analysis method applies successive filtering steps in this order:

| Step              | Source                  | Description                                                   |
| ----------------- | ----------------------- | ------------------------------------------------------------- |
| genotype          | Genotype file           | All subjects in `.fam` / `.psam` / VCF / BGEN sample block    |
| ∩ GRM             | Sparse GRM              | Subjects present in the sparse GRM (diagonal entries)         |
| ∩ keep            | `--keep` filter         | Restrict to keep-list (no-op when `--keep` is not provided)   |
| \ remove          | `--remove` filter       | Exclude remove-list (no-op when `--remove` is not provided)   |
| ∩ pheno/resid     | Phenotype / Residual    | Subjects with non-`NaN` values in `--null-resid` or `--pheno` |

The final step produces the working sample for that analysis.

> **Design note:** phenotype/residual intersection is applied last so
> that multi-residual analyses (e.g., SPACox with multiple columns)
> can apply per-phenotype subject masks on the already-filtered set.

### Formula

```
final = (genotype ∩ GRM) ∩ (keep \ remove) ∩ (pheno/resid dropna)
```

### GRM Intersection

Methods that use a sparse GRM (`--sp-grm-grab` or `--sp-grm-plink2`) require
subjects to exist in both the genotype file and the GRM file.  The GRM
subject IDs are pre-parsed (lightweight scan) before the full GRM load:

- **GRAB format** (`--sp-grm-grab`): scans the `ID1`/`ID2` columns.
- **GCTA format** (`--sp-grm-plink2`): reads the `.grm.id` sidecar.
  If `.grm.id` is missing, all genotype subjects are assumed to be in
  the GRM (a warning is logged).

Methods without a GRM (e.g., SPACox) skip the GRM step.

### Covariates Are Not in the Intersection

Covariates (`--covar`) are **not** part of the subject intersection.  A
subject that passes all pipeline steps but is missing from the covariate file
will have its covariate values filled with the per-column mean computed
from the covariate file.  This avoids discarding otherwise valid subjects
just because a covariate is unavailable.

After finalize, a covariate imputation table is printed:

```
── Covariate imputation (mean-fill) ─────────────────
  Column         N_filled / N_total
  AGE                  12 / 4700
  SEX                   0 / 4700
  PC1                   3 / 4700
─────────────────────────────────────────────────────
```

### Pipeline Log

GRAB prints a summary table after subject intersection:

```
── Subject pipeline ──────────────────────────────────
  Step              Count  Description
  genotype           5000  .fam / .psam / VCF / BGEN subjects
  ∩ GRM              4900  sparse GRM subjects
  ∩ keep             4900  (no --keep provided)
  \ remove           4900  (no --remove provided)
  ∩ pheno/resid      4700  phenotype / residual (non-NaN, final)
──────────────────────────────────────────────────────
```

---

## Marker Filtering

Markers pass through a QC gate before being tested.  A marker is
**excluded** (result columns set to `NA`) if any of the following hold:

| Flag         | Default   | Condition                     |
| ------------ | --------- | ----------------------------- |
| `--geno`     | 0.1       | missing rate > threshold      |
| `--maf`      | 1e-5      | MAF < threshold               |
| `--mac`      | 0.5       | MAC < threshold               |
| `--hwe`      | 0         | HWE p-value < threshold       |

- `--hwe 0` (default) disables HWE filtering.  Positive values use the
  plink2 convention: markers whose HWE p-value falls **below** the
  threshold are excluded.
- HWE p-values are computed using an exact test when any genotype count
  is < 5, otherwise a chi-squared approximation is used.

### Marker Include / Exclude

- `--extract FILE` — restrict analysis to markers listed in the file.
- `--exclude FILE` — exclude markers listed in the file.

These filters are applied before QC.

### Multi-Phenotype QC

When multiple residual columns are loaded (`--null-resid` with > 1
value column), each phenotype receives independent marker QC.  The union
bitmask covers subjects non-`NaN` in **any** column; per-phenotype
allele stats are computed from the per-phenotype non-missing set.

---

## Cross-References

- Genotype file format details: [genotype.md](genotype.md)
- CLI flag reference: `grab --help <method>`
