PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
OutputFile <- file.path(tempdir(), "resultSPAmix.txt")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)

# Step 1 option 1
obj.SPAmix <- GRAB.NullModel(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData,
  subjIDcol = "IID",
  method = "SPAmix",
  traitType = "time-to-event",
  control = list(PC_columns = "PC1,PC2")
)

# Step 1 option 2
residuals <- survival::coxph(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData
)$residuals

obj.SPAmix <- GRAB.NullModel(
  residuals ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData,
  subjIDcol = "IID",
  method = "SPAmix",
  traitType = "Residual",
  control = list(PC_columns = "PC1,PC2")
)

# Step 1 option 2: analyze multiple traits at once
res_cox <- survival::coxph(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData
)$residuals

res_lm <- lm(QuantPheno ~ AGE + GENDER + PC1 + PC2, data = PhenoData)$residuals

obj.SPAmix <- GRAB.NullModel(
  res_cox + res_lm ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData,
  subjIDcol = "IID",
  method = "SPAmix",
  traitType = "Residual",
  control = list(PC_columns = "PC1,PC2")
)

# Step 2
GRAB.Marker(obj.SPAmix, GenoFile, OutputFile)

head(data.table::fread(OutputFile))