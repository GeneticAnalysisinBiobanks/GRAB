PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
OutputFile <- file.path(tempdir(), "resultSPACox.txt")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)

# Step 1 option 1
obj.SPACox <- GRAB.NullModel(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
  data = PhenoData,
  subjIDcol = "IID",
  method = "SPACox",
  traitType = "time-to-event"
)

# Step 1 option 2
residuals <- survival::coxph(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
  data = PhenoData,
  x = TRUE
)$residuals

obj.SPACox <- GRAB.NullModel(
  residuals ~ AGE + GENDER,
  data = PhenoData,
  subjIDcol = "IID",
  method = "SPACox",
  traitType = "Residual"
)

# Step 2
GRAB.Marker(obj.SPACox, GenoFile, OutputFile)

head(data.table::fread(OutputFile))