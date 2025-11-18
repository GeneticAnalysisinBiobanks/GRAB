PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
GenoFileStep1 <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
GenoFileStep2 <- system.file("extdata", "simuPLINK_RV.bed", package = "GRAB")
GroupFile <- system.file("extdata", "simuPLINK_RV.group", package = "GRAB")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)
PhenoData$OrdinalPheno <- factor(PhenoData$OrdinalPheno, levels = c(0, 1, 2))

# Step 1
obj.POLMM <- GRAB.NullModel(
  factor(PhenoData$OrdinalPheno) ~ AGE + GENDER,
  data = PhenoData,
  subjIDcol = "IID",
  method = "POLMM",
  traitType = "ordinal",
  GenoFile = GenoFile,
  SparseGRMFile = SparseGRMFile,
  control = list(tolTau = 0.2, tolBeta = 0.1)
)

# Step 2
OutputFile <- file.path(tempdir(), "resultPOLMMmarker.txt")
GRAB.Marker(obj.POLMM, GenoFile, OutputFile,
  control = list(ifOutGroup = TRUE))

head(data.table::fread(OutputFile))


# Step 1
obj.POLMM <- GRAB.NullModel(
  PhenoData$OrdinalPheno ~ AGE + GENDER,
  data = data.table::fread(PhenoFile, header = TRUE),
  subjIDcol = "IID",
  method = "POLMM",
  traitType = "ordinal",
  GenoFile = GenoFileStep1,
  control = list(tolTau = 0.2, tolBeta = 0.1)
)

# Step 2
OutputFile <- file.path(tempdir(), "resultPOLMMregion2.txt")
GRAB.Region(obj.POLMM, GenoFileStep2, OutputFile,
  GroupFile = GroupFile,
  SparseGRMFile = SparseGRMFile,
  MaxMAFVec = "0.01,0.005"
)

head(data.table::fread(OutputFile))
head(data.table::fread(paste0(OutputFile, ".markerInfo")))
head(data.table::fread(paste0(OutputFile, ".otherMarkerInfo")))
head(data.table::fread(paste0(OutputFile, ".infoBurdenNoWeight")))