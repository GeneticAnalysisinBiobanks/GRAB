devtools::load_all(quiet = TRUE)

PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)
SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")

run_wtcoxg <- function() {
  cat("===== WtCoxG =====\n")
  RefAfFile <- system.file("extdata", "simuRefAf.txt", package = "GRAB")
  out <- file.path(tempdir(), "resultWtCoxG_marker4.txt")

  obj <- GRAB.NullModel(
    survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
    data = PhenoData,
    subjIDcol = "IID",
    method = "WtCoxG",
    traitType = "time-to-event",
    GenoFile = GenoFile,
    SparseGRMFile = SparseGRMFile,
    RefAfFile = RefAfFile,
    RefPrevalence = 0.1
  )

  out_base <- file.path(tempdir(), "resultWtCoxG_marker.txt")
  GRAB.Marker(obj, GenoFile, out_base, control = list(nMarkersEachChunk = 1000))
  a <- data.table::fread(out_base)

  GRAB.Marker4(obj, GenoFile, out, nthreads = 1, overwrite = TRUE, control = list(nMarkersEachChunk = 1000))
  b <- data.table::fread(out)
  cat("WtCoxG Marker vs Marker4(1):", isTRUE(all.equal(a, b, tolerance = 1e-10)), "\n")

  GRAB.Marker4(obj, GenoFile, out, nthreads = 4, overwrite = TRUE, control = list(nMarkersEachChunk = 1000))
  c <- data.table::fread(out)
  cat("WtCoxG Marker vs Marker4(4):", isTRUE(all.equal(a, c, tolerance = 1e-10)), "\n")
  cat("WtCoxG Marker4(1) vs Marker4(4):", isTRUE(all.equal(b, c, tolerance = 1e-10)), "\n")
}

run_leaf <- function() {
  cat("===== LEAF =====\n")
  RefAfFile <- system.file("extdata", "simuRefAf_2pop.txt", package = "GRAB")
  out <- file.path(tempdir(), "resultLEAF_marker4.txt")

  obj <- GRAB.NullModel(
    BinaryPheno ~ AGE + GENDER,
    data = PhenoData,
    subjIDcol = "IID",
    method = "LEAF",
    traitType = "binary",
    GenoFile = GenoFile,
    SparseGRMFile = SparseGRMFile,
    RefAfFile = RefAfFile,
    RefPrevalence = 0.1,
    Ncluster = 2,
    PCmatrix = PhenoData[, c("PC1", "PC2", "PC3", "PC4")]
  )

  out_base <- file.path(tempdir(), "resultLEAF_marker.txt")
  GRAB.Marker(obj, GenoFile, out_base, control = list(nMarkersEachChunk = 1000))
  a <- data.table::fread(out_base)

  GRAB.Marker4(obj, GenoFile, out, nthreads = 1, overwrite = TRUE, control = list(nMarkersEachChunk = 1000))
  b <- data.table::fread(out)
  cat("LEAF Marker vs Marker4(1):", isTRUE(all.equal(a, b, tolerance = 1e-10)), "\n")

  GRAB.Marker4(obj, GenoFile, out, nthreads = 4, overwrite = TRUE, control = list(nMarkersEachChunk = 1000))
  c <- data.table::fread(out)
  cat("LEAF Marker vs Marker4(4):", isTRUE(all.equal(a, c, tolerance = 1e-10)), "\n")
  cat("LEAF Marker4(1) vs Marker4(4):", isTRUE(all.equal(b, c, tolerance = 1e-10)), "\n")
}

run_wtcoxg()
run_leaf()
