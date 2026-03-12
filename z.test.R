devtools::load_all()

SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
PairwiseIBDFile <- system.file("extdata", "PairwiseIBD.txt", package = "GRAB")

PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)

GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
OutputFile <- file.path(tempdir(), "resultSPAsqr.txt")
#'
obj.SPAsqr <- GRAB.NullModel(
  QuantPheno ~ AGE + GENDER,
  data = PhenoData,
  subjIDcol = "IID",
  method = "SPAsqr",
  traitType = "quantitative",
  SparseGRMFile = SparseGRMFile,
  control = list(taus = c(0.2, 0.5, 0.8)),
  PairwiseIBDFile = PairwiseIBDFile
)

print("==================================")
print("Marker1 ...")
GRAB.Marker(obj.SPAsqr, GenoFile, OutputFile)
a <- data.table::fread(OutputFile)

print("==================================")
print("Marker4 one thread ...")
GRAB.Marker4(
    obj.SPAsqr, GenoFile, OutputFile,
    nthreads=1, overwrite=T, 
    control=list(nMarkersEachChunk = 100)
)
b <- data.table::fread(OutputFile)
print(all.equal(a, b))

print("==================================")
GRAB.Marker4(
  obj.SPAsqr, GenoFile, OutputFile, 
  nthreads=10, overwrite=T, 
  control=list(nMarkersEachChunk = 100)
)
c <- data.table::fread(OutputFile)
print(all.equal(a, c))
