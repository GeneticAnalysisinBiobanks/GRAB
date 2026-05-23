# GRAB 2.0 alpha

**GRAB** (Genome-wide Robust Analysis for Biobank data) is a free,
open-source toolkit of GWAS methods designed for biobank-scale data.  Up
to version 0.2.4, GRAB was distributed as an R package
([CRAN](https://CRAN.R-project.org/package=GRAB)). The code has since
been re-implemented from scratch in pure C++17; the new line is
numbered from **2.0** onwards and is referred to as **GRAB2** when
disambiguation is needed.

Note that GRAB2 does not include POLMM or POLMM-GENE. Users who wish to use these methods should continue using the original R package. For general marker-level analyses (including ordinal traits), we recommend SPACox, SPAMIX, or SPAGRM.

## Highlights

- **Single-thread performance.** GRAB2 is reimplemented in pure C++17 and runs approximately 10× to 1000× faster than the original R package on a single thread, with the precise factor depending on the method.
- **Near-linear parallel speedup.** Marker-level analysis is parallelized via SNP-chunk-based work stealing and scales nearly linearly with the number of threads.
- **Joint multi-phenotype analysis.** Multiple phenotypes may be specified in a single invocation and analyzed jointly, amortizing genotype I/O and null-model fitting across phenotypes.
- **Automatic LOCO handling.** Leave-one-chromosome-out (LOCO) predictions are routed to the corresponding chromosome automatically; no per-chromosome scripting is required.
- **Reduced memory footprint.** The engine streams genotype data in chunks and applies SIMD-vectorized kernels in place, keeping peak resident-set size substantially below that of the R package.

## Methods provided

| Method     | Supported trait types | Notes                                                                        |
|------------|-----------------------|------------------------------------------------------------------------------|
| **SPACox** | any [1]               | Baseline residual-based score test using saddlepoint approximation (SPA)     |
| **SPAmix** | any [1]               | Extends SPACox to admixed population                                         |
| **SPAGRM** | any [1]               | Extends SPACox to account for sample relatedness                             |
| **SAGELD** | longitudinal          | Dedicated for testing longitudinal gene–environment interactions             |
| **SPAsqr** | quantitative          | Smoothed quantile regression with LOCO PRS to improve statistical power      |
| **WtCoxG** | time-to-event, binary | Leverages external reference allele frequencies to improve statistical power |
| **LEAF**   | time-to-event, binary | Extends WtCoxG to multiple reference panels and heterogeneous cohorts        |

**Notes:**

- [1] For SPACox, SPAmix, and SPAGRM, the phenotype may be supplied
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

GRAB is self-contained: all third-party libraries are bundled in the
source tree. To build it, you only need a standard compiler toolchain
with C++17 support (a recent `gcc`/`g++` on Linux or MSYS2/MinGW, or
the Xcode Command Line Tools on macOS) and GNU `make`.

```bash
git clone --depth=1 https://github.com/GeneticAnalysisinBiobanks/GRAB.git
cd GRAB
make -j
```

This produces a single binary `build/grab2`, tuned for the CPU you built on.
Copy it anywhere on `PATH` and you are done.

To produce a binary that can be shared across any AVX2-capable machine, which is useful when the run host is older than the build host:

```bash
make -j GRAB_MARCH=-march=x86-64-v2
```

## Quick start

GRAB consumes genotype files in the same formats as PLINK 2:

- PLINK 2 `--pfile <prefix>` triples (`*.pgen` + `*.pvar` + `*.psam`).
- PLINK 1 `--bfile <prefix>` triples (`*.bed` + `*.bim` + `*.fam`).
- BGEN 1.1 / 1.2 / 1.3 (`--bgen <filename> <REF/ALT mode>`).
- VCF, optionally BGZF-compressed (`--vcf <filename>`).
- BCF2 binary (`--bcf <filename>`).

The example below runs the **SPAsqr** method on two quantitative
phenotypes (`Quantitative` and `Time`) jointly.  The options `--pheno`,
`--pheno-name`, `--covar-name`, and `--out` are compatible
with the PLINK 2 input and output conventions.
`--sp-grm-plink2` consumes the output of `plink2 --make-grm-sparse`, for example:

```bash
plink2 --pfile examples/1kg --make-grm-sparse 0.125 --out examples_output/1kg
```

`--pred-list` points at a text file with one row per phenotype, each
row giving a phenotype name and the path to an LOCO PRS file for that
phenotype.  Two layouts are auto-detected: **Regenie** (chromosome-major,
produced by `regenie --step 1`) and **LDAK-KVIK** (subject-major,
produced by `ldak6 --kvik-step1`).

```bash
build/grab2 \
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

## Licence

![Licence](https://img.shields.io/badge/licence-GPL--3.0+-blue)

GRAB 2.0 is released under the **GNU General Public License, version 3
or later**.

Vendored libraries under `third_party/` retain their upstream licences:

| Library                       | Licence                    |
|-------------------------------|----------------------------|
| **pgenlib**                   | LGPL-3-or-later            |
| **Eigen**                     | MPL-2.0                    |
| **htslib**, **libdeflate**    | MIT                        |
| **zlib**                      | zlib License               |
| **zstd**                      | BSD-3-Clause               |
| **Boost** (subset), **bgen**  | Boost Software License 1.0 |

## Citation

If GRAB contributes to a publication, please cite the method-specific
paper(s) listed below:

- **SPACox** — Bi *et al* (2020). A Fast and Accurate Method for Genome-Wide Time-to-Event Data Analysis and Its Application to UK Biobank. *Am. J. Hum. Genet.*
  [doi:10.1016/j.ajhg.2020.06.003](https://doi.org/10.1016/j.ajhg.2020.06.003)
- **SPAmix** — Ma *et al* (2025). SPAmix: a scalable, accurate, and universal analysis framework for large-scale genetic association studies in admixed populations. *Genome Biol.*
  [doi:10.1186/s13059-025-03827-9](https://doi.org/10.1186/s13059-025-03827-9)
- **SPAGRM** — Xu *et al* (2025). SPA(GRM): effectively controlling for sample relatedness in large-scale genome-wide association studies of longitudinal traits *Nat. Commun.*
  [doi:10.1038/s41467-025-56669-1](https://doi.org/10.1038/s41467-025-56669-1)
- **WtCoxG** — Li *et al* (2025). Applying weighted Cox regression to genome-wide association studies of time-to-event phenotypes *Nat. Comput. Sci.*
  [doi:10.1038/s43588-025-00864-z](https://doi.org/10.1038/s43588-025-00864-z)
- **SAGELD** — Xu *et al* (in prep). Leveraging longitudinal data to boost statistical power for gene-environment interaction analysis.
- **SPAsqr** — Heng *et al* (in prep). Discovering and Dissecting Heterogeneous Genetic Associations with Genome-wide Smoothed Quantile Regression
- **LEAF** — Li *et al* (in prep). Leveraging external allele frequency to boost powers of genome-wide association studies.
