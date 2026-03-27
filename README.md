# GRABmtMarker

The feat/mt branch introduces a multithreaded implementation of GRAB.Marker, provided as a new function GRAB.mtMarker. The function interface is largely unchanged, with the following differences:

- Added nthreads to the control list to specify the number of threads.
- Replaced OutputFileIndex with an overwrite = TRUE/FALSE argument.
- Added support for gzipped output when the output filename ends with the .gz suffix.

## Installation

Install directly from the feat/mt branch:

```r
remotes::install_github(
  "GeneticAnalysisinBiobanks/GRAB",
  ref = "feat/mt",
  dependencies = FALSE
)
```

## Examples

```r
# SPAsqr
extdir <- system.file("extdata", package = "GRABmtMarker")
GenoFile <- paste0(extdir, "/simuPLINK.bed")
OutputFile <- tempfile(fileext = ".tsv.gz")

control <- list(nMarkersEachChunk = 256, nthreads = 3)
objNull <- readRDS(paste0(extdir, "/obj.SPAsqr.rds"))

GRABmtMarker::GRAB.mtMarker(objNull, GenoFile, OutputFile, control = control)
head(read.table(OutputFile, header = TRUE))

# SPAmixPlus
objNull <- readRDS(paste0(extdir, "/obj.SPAmixPlus.rds"))
control$afFilePath = paste0(extdir, "/afModels.bin")

GRABmtMarker::GRAB.mtMarker(objNull, GenoFile, OutputFile, control = control, overwrite = TRUE)
head(read.table(OutputFile, header = TRUE))
```
