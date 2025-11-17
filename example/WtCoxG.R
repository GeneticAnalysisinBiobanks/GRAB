# Step 1: fit null model and test batch effect
PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)
SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
RefAfFile <- system.file("extdata", "simuRefAf.txt", package = "GRAB")
OutputFile <- file.path(tempdir(), "resultWtCoxG.txt")

obj.WtCoxG <- GRAB.NullModel(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
  data = PhenoData,
  subjIDcol = "IID",
  method = "WtCoxG",
  traitType = "time-to-event",
  GenoFile = GenoFile,
  SparseGRMFile = SparseGRMFile,
  control = list(RefPrevalence = 0.01, SNPnum = 1e3),
  RefAfFile = RefAfFile,
  SurvTimeColumn = "SurvTime",
  IndicatorColumn = "SurvEvent"
)

GRAB.Marker(obj.WtCoxG, GenoFile, OutputFile,
  control = list(nMarkersEachChunk = 5000)
)

head(data.table::fread(OutputFile))