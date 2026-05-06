# GRAB

**GRAB** is an R package that provides a comprehensive suite of GWAS methods for biobank-scale data. For detailed instructions, see the [GRAB manual page](https://wenjianbi.github.io/grab.github.io/).

Version 0.1.2 (the last version before v0.2.0 and prior to June 2025) is archived in branch [release/v0.1.2](https://github.com/GeneticAnalysisinBiobanks/GRAB/tree/release/v0.1.2).

## Installation

### [![CRAN Status](https://www.r-pkg.org/badges/version/GRAB)](https://CRAN.R-project.org/package=GRAB) CRAN

![Linux](https://img.shields.io/badge/Linux-000?logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/macOS-000?logo=apple&logoColor=white)
![Windows](https://img.shields.io/badge/Windows-0078D6?logo=windows&logoColor=white)

Install GRAB from CRAN in your R console:

```r
install.packages("GRAB", dependencies = TRUE)
```

### [![Conda-Forge](https://img.shields.io/conda/vn/conda-forge/r-grab.svg)](https://anaconda.org/conda-forge/r-grab) Conda

![Linux](https://img.shields.io/badge/Linux-000?logo=linux&logoColor=white)
![macOS](https://img.shields.io/badge/macOS-000?logo=apple&logoColor=white)

Install GRAB in a new Conda environment named `grab_env` from the `conda-forge` channel:

```sh
conda create -n grab_env -c conda-forge r-grab r-skat r-dbplyr r-tidyr
```

### [![Docker Image Version](https://img.shields.io/docker/v/geneticanalysisinbiobanks/grab?sort=semver&label=Docker%20latest)](https://hub.docker.com/r/geneticanalysisinbiobanks/grab) Docker

![Linux](https://img.shields.io/badge/Linux-000?logo=linux&logoColor=white)

Pull the latest GRAB Docker image from Docker Hub:

```sh
docker pull geneticanalysisinbiobanks/grab:latest
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
  subjIDcol = "IID",
  method = "SPAmix",
  traitType = "time-to-event",
  control = list(PC_columns = "PC1,PC2")
)
```

### Step 2: conduct score test

```r
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
OutputFile <- file.path(tempdir(), "Results_SPAmix.txt")

GRAB.Marker(obj.SPAmix, GenoFile = GenoFile, OutputFile = OutputFile)
data.table::fread(OutputFile)
```
