devtools::load_all(quiet = TRUE)
PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)
SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
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
GRAB.Marker4(obj, GenoFile, out, nthreads = 1, overwrite = TRUE, control = list(nMarkersEachChunk = 1000))
cat("WtCoxG Marker4(1) finished\n")
GRAB.Marker4(obj, GenoFile, out, nthreads = 4, overwrite = TRUE, control = list(nMarkersEachChunk = 1000))
cat("WtCoxG Marker4(4) finished\n")
