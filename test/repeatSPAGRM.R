devtools::load_all()

PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)
SparseGRMFile <- file.path("test/SparseGRM.txt")
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
OutputFile <- file.path("test/repeatSPAGRM.output.txt")

obj.SPAGRM_list <- GRAB.NullModel(
  QuantPheno ~ AGE + GENDER,
  data = PhenoData,
  subjIDcol = "IID",
  method = "SPAsqr",
  traitType = "quantitative",
  GenoFile = GenoFile,
  SparseGRMFile = SparseGRMFile,
  control = list(
    taus = c(0.05, 0.2, 0.5, 0.8, 0.95),
    h = 0,
    sqr_tol = 1e-9,
    ControlOutlier = FALSE
  )
)

saveRDS(obj.SPAGRM_list, file = "test/repeatSPAGRM.obj_list.rds")

cat("\n####################\n\n")
GRAB.Marker(obj.SPAGRM_list, GenoFile, OutputFile = OutputFile, control = list(tol = 1e-9))
