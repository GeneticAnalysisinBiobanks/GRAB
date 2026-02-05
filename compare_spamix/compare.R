devtools::load_all()

SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
GRM <- data.table::fread(SparseGRMFile)

PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
PhenoData <- data.table::fread(PhenoFile)
subjData <- PhenoData$IID

extdata_dir <- system.file("extdata", package = "GRAB")
dosagePrefix <- paste0(extdata_dir, "/simuAncestry")
phiOutputPrefix <- "./compare_spamix/phi_"
outPrefix <- phiOutputPrefix

# Estimate phi for ancestries 1 and 2
SPAmixLocalPlus.EstimatePhi(
  GRM = GRM,
  dosagePrefix = dosagePrefix,
  haploPrefix = NULL,  # Uses dosagePrefix if NULL
  ancIdx = c(1, 2),
  SampleIDs = subjData,
  phiOutputPrefix = phiOutputPrefix,
  MAF_cutoff = 0.01
)

# Fit null model to get residuals
residuals <- survival::coxph(
  survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
  data = PhenoData
)$residuals

# Run SPAmixPlus tests for ancestries 1 and 2
result_lst <- SPAmixLocalPlus.Marker(
  resid = residuals,
  subjData = subjData,
  dosagePrefix = dosagePrefix,
  haploPrefix = NULL,  # Uses dosagePrefix if NULL
  phiPrefix = phiOutputPrefix,
  ancIdx = c(1, 2),
  outPrefix = outPrefix
)

# View results
OutputFile1 <- paste0(phiOutputPrefix, "Ancestry1.txt")
head(data.table::fread(OutputFile1))
