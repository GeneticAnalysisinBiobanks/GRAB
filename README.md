# GRAB

[![GitHub source](https://img.shields.io/badge/GitHub-source-black?logo=github)](https://github.com/GeneticAnalysisinBiobanks/GRAB)
[![License: GPL v2+](https://img.shields.io/badge/License-GPL%20v2%2B-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)

**GRAB** is an R package that provides a comprehensive suite of GWAS methods for biobank-scale data. For detailed instructions, see the [GRAB manual page](https://wenjianbi.github.io/grab.github.io/).

Version 0.1.2 (the last version before v0.2.0 and prior to June 2025) is archived in branch [release/v0.1.2](https://github.com/GeneticAnalysisinBiobanks/GRAB/tree/release/v0.1.2).

## Installation

GRAB is an R package, with part of its code written in C++ for improved performance. GRAB can be installed on Linux, Windows, or macOS via CRAN, Conda, or from source code.

### Install via CRAN

[![CRAN Status](https://www.r-pkg.org/badges/version/GRAB)](https://CRAN.R-project.org/package=GRAB)
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/GRAB)](https://CRAN.R-project.org/package=GRAB)

Install GRAB from CRAN in your R console:

```r
install.packages("GRAB")
```

### Install via Conda

[![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/r-grab.svg)](https://anaconda.org/conda-forge/r-grab)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/r-grab/badges/downloads.svg)](https://anaconda.org/conda-forge/r-grab)

Install GRAB in a new Conda environment named `grab_env` from the `conda-forge` channel:

```sh
conda create -n grab_env -c conda-forge r-grab
```

## Quick tutorial

Here is a quick tutorial for GWAS of a time-to-event trait using SPAmix.

### Step 1: fit a null model

```r
library(GRAB)
PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)

obj.SPAmix <- GRAB.NullModel(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData,
  subjData = IID,
  method = "SPAmix",
  traitType = "time-to-event",
  control = list(PC_columns = "PC1,PC2")
)
```

### Step 2: conduct score test

```r
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
OutputFile <- file.path(tempdir(), "Results_SPACox.txt")

GRAB.Marker(
  objNull = obj.SPACox,
  GenoFile = GenoFile,
  OutputFile = OutputFile,
  control = list(outputColumns = "zScore")
)
data.table::fread(OutputFile)
```
