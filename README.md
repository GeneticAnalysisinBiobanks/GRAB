# GRAB 2.0 alpha

**GRAB** (Genome-wide Robust Analysis for Biobank data) is a free,
open-source toolkit of GWAS methods designed for biobank-scale data.  Up
to version 0.2.4, GRAB was distributed as an R package
([CRAN](https://CRAN.R-project.org/package=GRAB)).  The code has since
been re-implemented from scratch in pure C++17; the new line is
numbered from **2.0** onwards and is referred to as **GRAB2** when
disambiguation is needed. For detailed instructions,
see the [GRAB 2.0 manual page](https://wenjianbi.github.io/grab.github.io/).

## Methods provided

Selected with `--method <name>`.

| Method | Supported trait types | Notes                                                                     |
|--------|-----------------------|---------------------------------------------------------------------------|
| SPACox | any [1]               | Baseline residual-based score test using saddlepoint approximation (SPA)  |
| SPAmix | any [1]               | Extends SPACox to admixed population                                      |
| SPAGRM | any [1]               | Extends SPACox to account for sample relatedness                          |
| SAGELD | longitudinal          | Dedicated for testing longitudinal gene–environment interactions          |
| SPAsqr | quantitative          | Smoothed quantile regression with LOCO PRS to improve statistical power   |
| WtCoxG | time-to-event, binary | Leverages allele frequencies from an external reference population to improve statistical power |
| LEAF   | time-to-event, binary | Extends WtCoxG to multiple reference panels and heterogeneous cohorts     |

**Notes:**

- **[1]** For SPACox, SPAmix, and SPAGRM, the phenotype may be supplied
  in either of two ways: as a **raw phenotype column** via
  `--pheno-name`, in which case GRAB's built-in null-model fitter
  accepts quantitative, binary, time-to-event, and ordinal traits; or
  as **pre-computed residuals** via `--resid-name`, in which case any
  trait type is supported (the user is responsible for producing the
  residuals from an appropriate model upstream).

## Build and install

![Linux](https://img.shields.io/badge/Linux-000?logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/macOS-000?logo=apple&logoColor=white)
![Windows](https://img.shields.io/badge/Windows-0078D6?logo=windows&logoColor=white)

GRAB is a self-contained C++17 application.  All third-party libraries
are bundled in the source tree, so the only things you need on your
machine are a recent C++ compiler and GNU `make` — no extra installs.

```bash
git clone --depth=1 https://github.com/GeneticAnalysisinBiobanks/GRAB.git
cd GRAB
make -j
```

This produces a single binary at `build/grab`, tuned for the CPU you built on.
Copy it anywhere on `PATH` and you are done.

To produce a binary that can be shared across any AVX2-capable machine, which is useful when the build host differs from the run host and the run host is older:

```bash
make -j GRAB_MARCH=-march=x86-64-v2
```

## Quick start

GRAB consumes genotype files in the same formats as PLINK 2:

- PLINK 2 `--pfile <prefix>` triples (`*.pgen` + `*.pvar` + `*.psam`).
- PLINK 1 `--bfile <prefix>` triples (`*.bed` + `*.bim` + `*.fam`).
- BGEN 1.1 / 1.2 / 1.3 (`--bgen <filename> <REF/ALT mode>`).
- VCF / BCF (`--vcf <filename>`).

The example below runs the **SPAsqr** method on two quantitative
phenotypes (`Quantitative` and `Time`) jointly.  The options `--pheno`,
`--pheno-name`, `--covar-name`, and `--out` are compatible
with the PLINK 2 phenotype and covariate input conventions.
`--sp-grm-plink2` consumes the output of `plink2 --make-grm-sparse`.

`--pred-list` points at a text file with one row per phenotype, each
row giving a phenotype name and the path to an LOCO PRS file for that
phenotype.  Two layouts are auto-detected: **Regenie** (chromosome-major,
produced by `regenie --step 1`) and **LDAK-KVIK** (subject-major,
produced by `ldak6 --kvik-step1`).

```bash
build/grab \
  --method SPAsqr \
  --pheno examples/1kg.pheno \
  --pheno-name Quantitative,Time \
  --covar-name MALE,PC1,PC2,PC3,PC4 \
  --sp-grm-plink2 examples/1kg.grm.sp \
  --pred-list examples/loco_prs.list \
  --pfile examples/1kg \
  --out examples_output/1kg
```

The companion script `examples/run.sh` exercises every method end-to-end using the files under `examples/`.

## Documentation

Per-method documentation, including the statistical model, the input file
specification, and the output column reference, lives under
[`docs/methods/`](docs/methods/) and [`docs/software/`](docs/software/).

## Licence

![Licence](https://img.shields.io/badge/licence-GPL--3.0+-blue)

GRAB 2.0 is released under the **GNU General Public License, version 3
or later**. The distribution comprises both the full source tree and precompiled
binaries built from it; the same licence applies to both.

The third-party libraries vendored under
[`third_party/`](third_party/) are shipped as source and statically
linked into the precompiled binaries.  Each retains its own upstream
licence:

| Library            | Purpose                        | Licence                                          |
|--------------------|--------------------------------|--------------------------------------------------|
| **pgenlib**        | PLINK 2 file-format reader     | LGPL-3-or-later                                  |
| **Eigen**          | dense linear algebra           | MPL-2.0                                          |
| **htslib**         | VCF / BCF reader               | MIT                                              |
| **libdeflate**     | fast BGZF (pgenlib, htslib)    | MIT                                              |
| **zlib**           | gzip I/O (GRAB, bgen, htslib)  | zlib License                                     |
| **zstd**           | Zstandard compression          | BSD-3-Clause (BSD branch of zstd's dual licence) |
| **Boost** (subset) | Boost.Math distributions       | Boost Software License 1.0                       |
| **bgen**           | BGEN file-format reader        | Boost Software License 1.0                       |

Verbatim copies of each upstream licence are preserved under the
corresponding subdirectory of `third_party/`.  Any redistribution of a
GRAB binary — whether one shipped from this repository or one rebuilt
from this source tree — must carry all of these licences alongside it.

## Citation

If GRAB contributes to a publication, please cite the method-specific
paper(s) listed below:

- **SPACox** — Bi *et al*. (2020).  Fast and accurate method for
  genome-wide time-to-event data analysis and its application to UK
  Biobank.  *Am. J. Hum. Genet.*
  [doi:10.1016/j.ajhg.2020.06.003](https://doi.org/10.1016/j.ajhg.2020.06.003)
- **SPAmix** — Ma *et al*. (2025).  Sparse estimation of high-dimensional
  genetic correlation and its application to global biobank
  meta-analysis.  *Genome Biol.*
  [doi:10.1186/s13059-025-03827-9](https://doi.org/10.1186/s13059-025-03827-9)
- **SPAGRM** — Xu *et al*. (2025).  Scalable and accurate variance
  component analysis with large sample relatedness.
  *Nat. Commun.*
  [doi:10.1038/s41467-025-56669-1](https://doi.org/10.1038/s41467-025-56669-1)
- **WtCoxG** — Li *et al*. (2025).  High-powered, robust, and versatile
  survival analysis via weighted Cox regression.
  *Nat. Comput. Sci.*
  [doi:10.1038/s43588-025-00864-z](https://doi.org/10.1038/s43588-025-00864-z)
- **SAGELD** — Xu *et al*. (in preparation)
- **SPAsqr** — Heng *et al*. (in preparation)
- **LEAF** — Li *et al*. (in preparation)
