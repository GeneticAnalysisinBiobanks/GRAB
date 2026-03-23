library(data.table)
devtools::load_all()

extdata_dir <- system.file("extdata", package = "GRAB")
SparseGRMFile <- file.path(extdata_dir, "SparseGRM.txt")
GenoFile <- file.path(extdata_dir, "simuPLINK.bed")

ResidMatFile <- file.path(extdata_dir, "ResidMat.txt")
PairwiseIBDFile <- file.path(extdata_dir, "PairwiseIBD.txt")
RefAfFile <- file.path(extdata_dir, "simuRefAf.txt")
RefAfFile2 <- file.path(extdata_dir, "simuRefAf_2pop.txt")

PhenoFile <- file.path(extdata_dir, "simuPHENO.txt")
PhenoData <- fread(PhenoFile, header = TRUE)
PCmatrix <- as.matrix(PhenoData[, c("PC1", "PC2", "PC3", "PC4")])
PhenoData$OrdinalPheno <- factor(PhenoData$OrdinalPheno, levels = c(0, 1, 2))

LongDataFile <- file.path(extdata_dir, "simuLongPHENO.txt")
LongPheno <- fread(LongDataFile)
afFileOutput <- ".test/afModels.bin"
SPAmixPlus.AF(GenoFile, PCmatrix, PhenoData$IID, afFileOutput)

nullobj_dir <- ".test/nullobj"
nullobj_spagrm <- file.path(nullobj_dir, "obj.SPAGRM.rds")
nullobj_spasqr <- file.path(nullobj_dir, "obj.SPAsqr.rds")
nullobj_wtcoxg <- file.path(nullobj_dir, "obj.WtCoxG.rds")
nullobj_leaf <- file.path(nullobj_dir, "obj.LEAF.rds")

nullobj_spacox <- file.path(nullobj_dir, "obj.SPACox.rds")
nullobj_spamix <- file.path(nullobj_dir, "obj.SPAmix.rds")
nullobj_spamixplus <- file.path(nullobj_dir, "obj.SPAmixPlus.rds")
nullobj_polmm <- file.path(nullobj_dir, "obj.POLMM.rds")
nullobj_sageld <- file.path(nullobj_dir, "obj.SAGELD.rds")

result_dir <- ".test/dev_results"
devel_spagrm <- file.path(result_dir, "resultSPAGRM.txt")
devel_spasqr <- file.path(result_dir, "resultSPAsqr.txt")
devel_wtcoxg <- file.path(result_dir, "resultWtCoxG.txt")
devel_leaf <- file.path(result_dir, "resultLEAF.txt")

devel_spacox <- file.path(result_dir, "resultSPACox.txt")
devel_spamix <- file.path(result_dir, "resultSPAmix.txt")
devel_spamixplus <- file.path(result_dir, "resultSPAmixPlus1.txt")
devel_polmm <- file.path(result_dir, "resultPOLMM.txt")
devel_sageld <- file.path(result_dir, "resultSAGELD.txt")

mt_dir <- ".test/mt_results"
mt_spagrm <- file.path(mt_dir, "resultSPAGRM.txt.gz")
mt_spasqr <- file.path(mt_dir, "resultSPAsqr.txt.gz")
mt_wtcoxg <- file.path(mt_dir, "resultWtCoxG.txt.gz")
mt_leaf <- file.path(mt_dir, "resultLEAF.txt.gz")

mt_spacox <- file.path(mt_dir, "resultSPACox.txt.gz")
mt_spamix <- file.path(mt_dir, "resultSPAmix.txt.gz")
mt_spamixplus <- file.path(mt_dir, "resultSPAmixPlus1.txt.gz")
mt_polmm <- file.path(mt_dir, "resultPOLMM.txt.gz")
mt_sageld <- file.path(mt_dir, "resultSAGELD.txt.gz")


## SPAGRM SPAsqr
obj.SPAGRM <- SPAGRM.NullModel(
  ResidMatFile = ResidMatFile,
  SparseGRMFile = SparseGRMFile,
  PairwiseIBDFile = PairwiseIBDFile,
  control = list(ControlOutlier = FALSE)
)

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

## WtCoxG LEAF
obj.WtCoxG <- GRAB.NullModel(
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

obj.LEAF <- GRAB.NullModel(
  BinaryPheno ~ AGE + GENDER,
  data = PhenoData,
  subjIDcol = "IID",
  method = "LEAF",
  traitType = "binary",
  GenoFile = GenoFile,
  SparseGRMFile = SparseGRMFile,
  RefAfFile = RefAfFile2,
  RefPrevalence = 0.1,
  Ncluster = 2,
  PCmatrix = PCmatrix
)

## SPACox SPAmix SPAmixPlus
residuals <- survival::coxph(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData
)$residuals

obj.SPACox <- GRAB.NullModel(
  residuals ~ AGE + GENDER,
  data = PhenoData,
  subjIDcol = "IID",
  method = "SPACox",
  traitType = "Residual"
)

obj.SPAmix <- GRAB.NullModel(
  residuals ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData,
  subjIDcol = "IID",
  method = "SPAmix",
  traitType = "Residual",
  control = list(PC_columns = "PC1,PC2")
)

obj.SPAmixPlus <- GRAB.NullModel(
  residuals ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData,
  subjIDcol = "IID",
  SparseGRMFile = SparseGRMFile,
  method = "SPAmixPlus",
  traitType = "Residual",
  control = list(PC_columns = "PC1,PC2")
)

## POLMM
obj.POLMM <- GRAB.NullModel(
 OrdinalPheno ~ AGE + GENDER,
 data = PhenoData,
 subjIDcol = "IID",
 method = "POLMM",
 traitType = "ordinal",
 GenoFile = GenoFile,
 SparseGRMFile = SparseGRMFile
)

## SAGELD
lme4null <- lme4::lmer(LongPheno ~ AGE + GENDER + (AGE|IID), data = LongPheno)
obj.SAGELD <- SAGELD.NullModel(
  NullModel = lme4null,
  UsedMethod = "SAGELD",
  PlinkFile = GenoFile,
  SparseGRMFile = SparseGRMFile,
  PairwiseIBDFile = PairwiseIBDFile,
)

## save null models
saveRDS(obj.SPAGRM, file = nullobj_spagrm)
saveRDS(obj.SPAsqr, file = nullobj_spasqr)
saveRDS(obj.WtCoxG, file = nullobj_wtcoxg)
saveRDS(obj.LEAF, file = nullobj_leaf)
saveRDS(obj.SPACox, file = nullobj_spacox)
saveRDS(obj.SPAmix, file = nullobj_spamix)
saveRDS(obj.SPAmixPlus, file = nullobj_spamixplus)
saveRDS(obj.POLMM, file = nullobj_polmm)
saveRDS(obj.SAGELD, file = nullobj_sageld)

## run Marker in develop branch
GRAB.Marker(obj.SPAGRM, GenoFile, devel_spagrm)
GRAB.Marker(obj.SPAsqr, GenoFile, devel_spasqr)
GRAB.Marker(obj.WtCoxG, GenoFile, devel_wtcoxg)
GRAB.Marker(obj.LEAF, GenoFile, devel_leaf)
GRAB.Marker(obj.SPACox, GenoFile, devel_spacox)
GRAB.Marker(obj.SPAmix, GenoFile, devel_spamix)
GRAB.Marker(obj.SPAmixPlus, GenoFile, devel_spamixplus, control = list(afFilePath = afFileOutput))
GRAB.Marker(obj.POLMM, GenoFile, devel_polmm)
GRAB.Marker(obj.SAGELD, GenoFile, devel_sageld)

########## load results
obj.SPAGRM <- readRDS(nullobj_spagrm)
obj.SPAsqr <- readRDS(nullobj_spasqr)
obj.WtCoxG <- readRDS(nullobj_wtcoxg)
obj.LEAF <- readRDS(nullobj_leaf)
obj.SPACox <- readRDS(nullobj_spacox)
obj.SPAmix <- readRDS(nullobj_spamix)
obj.SPAmixPlus <- readRDS(nullobj_spamixplus)
obj.POLMM <- readRDS(nullobj_polmm)
obj.SAGELD <- readRDS(nullobj_sageld)

dev_res_spagrm <- fread(devel_spagrm)
dev_res_spasqr <- fread(devel_spasqr)
dev_res_wtcoxg <- fread(devel_wtcoxg)
dev_res_leaf <- fread(devel_leaf)
dev_res_leaf$meta.p_ext[dev_res_leaf$meta.p_ext == 1] <- NA
dev_res_spacox <- fread(devel_spacox)
dev_res_spamix <- fread(devel_spamix)
dev_res_spamixplus <- fread(devel_spamixplus)
dev_res_polmm <- fread(devel_polmm)
dev_res_sageld <- fread(devel_sageld)

plink_acount <- fread(".test/simuPLINK.acount")
plink_hardy <- fread(".test/simuPLINK.hardy")

compare_df <- function(df1, df2, name = "") {
  common_cols <- intersect(names(df1), names(df2))
  cat("====", name, "====\n")

  for (col in common_cols) {
    res <- all.equal(df1[[col]], df2[[col]], tolerance = 1e-5)
    cat("Column:", col, res, "\n")
  }
  cat("\n")
}

########## run Marker6
control <- list(nthreads = 4, nMarkersEachChunk = 256, afFilePath = afFileOutput)

GRAB.Marker6(obj.SPAGRM, GenoFile, mt_spagrm, control = control, overwrite = TRUE)
GRAB.Marker6(obj.SPAsqr, GenoFile, mt_spasqr, control = control, overwrite = TRUE)
GRAB.Marker6(obj.WtCoxG, GenoFile, mt_wtcoxg, control = control, overwrite = TRUE)
GRAB.Marker6(obj.LEAF, GenoFile, mt_leaf, control = control, overwrite = TRUE)
GRAB.Marker6(obj.SPACox, GenoFile, mt_spacox, control = control, overwrite = TRUE)
GRAB.Marker6(obj.SPAmix, GenoFile, mt_spamix, control = control, overwrite = TRUE)
GRAB.Marker6(obj.SPAmixPlus, GenoFile, mt_spamixplus, control = control, overwrite = TRUE)
GRAB.Marker6(obj.POLMM, GenoFile, mt_polmm, control = control, overwrite = TRUE)
GRAB.Marker6(obj.SAGELD, GenoFile, mt_sageld, control = control, overwrite = TRUE)

mt_res_spagrm <- fread(mt_spagrm)
mt_res_spasqr <- fread(mt_spasqr)
mt_res_wtcoxg <- fread(mt_wtcoxg)
mt_res_leaf <- fread(mt_leaf)
mt_res_spacox <- fread(mt_spacox)
mt_res_spamix <- fread(mt_spamix)
mt_res_spamixplus <- fread(mt_spamixplus)
mt_res_polmm <- fread(mt_polmm)
mt_res_sageld <- fread(mt_sageld)

compare_df(dev_res_spagrm, mt_res_spagrm, "spagrm")
compare_df(dev_res_spasqr, mt_res_spasqr, "spasqr")
compare_df(dev_res_wtcoxg, mt_res_wtcoxg, "wtcoxg")
compare_df(dev_res_leaf, mt_res_leaf, "leaf")
compare_df(dev_res_polmm, mt_res_polmm, "polmm")
compare_df(dev_res_sageld, mt_res_sageld, "sageld")
compare_df(dev_res_spacox, mt_res_spacox, "spacox")
compare_df(dev_res_spamix, mt_res_spamix, "spamix")
compare_df(dev_res_spamixplus, mt_res_spamixplus, "spamixplus")
