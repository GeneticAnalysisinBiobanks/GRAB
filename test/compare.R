devtools::load_all()

PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
PhenoData <- data.table::fread(PhenoFile, header = TRUE)
SparseGRMFile <- file.path("test/SparseGRM.txt")
GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")

obj2 <- GRAB.NullModel(
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
    ControlOutlier = FALSE #, OutlierRatio = 50 # 1.5 for compare with obj1_lst
  )
)

# Step 1: compare obj1_lst and obj2
obj1_lst <- readRDS("test/repeatSPAGRM.obj_list.rds")
sapply(1:5, function(i) identical(obj1_lst[[i]]$Resid, obj2$Resid_mat[, i]))
sapply(1:5, function(i) identical(obj1_lst[[i]]$Resid.unrelated.outliers, obj2$Resid.unrelated.outliers[[i]]))
sapply(1:5, function(i) identical(obj1_lst[[i]]$R_GRM_R, obj2$R_GRM_R_vec[[i]]))
sapply(1:5, function(i) identical(obj1_lst[[i]]$R_GRM_R_TwoSubjOutlier, obj2$R_GRM_R_TwoSubjOutlier_vec[[i]]))
sapply(1:5, function(i) identical(obj1_lst[[i]]$sum_R_nonOutlier, obj2$sum_R_nonOutlier_vec[[i]]))
sapply(1:5, function(i) identical(obj1_lst[[i]]$R_GRM_R_nonOutlier, obj2$R_GRM_R_nonOutlier_vec[[i]]))
sapply(1:5, function(i) identical(obj1_lst[[i]]$TwoSubj_list, obj2$TwoSubj_list_lst[[i]]))
sapply(1:5, function(i) identical(obj1_lst[[i]]$ThreeSubj_list, obj2$ThreeSubj_list_lst[[i]]))


# Step 2: compare P-values
OutputFile <- file.path("test/refactorSPAsqr.output.txt")
GRAB.Marker(obj2, GenoFile, OutputFile = OutputFile, control = list(tol = 1e-9))

res1df <- data.table::fread("test/repeatSPAGRM.output.txt")
res2df <- data.table::fread(OutputFile)
  
all.equal(res1df$Pvalue_tau0.05, res2df$P_tau0.05)
all.equal(res1df$Pvalue_tau0.2, res2df$P_tau0.2)
all.equal(res1df$Pvalue_tau0.5, res2df$P_tau0.5)
all.equal(res1df$Pvalue_tau0.8, res2df$P_tau0.8)
all.equal(res1df$Pvalue_tau0.95, res2df$P_tau0.95)
all.equal(res1df$Pvalue_CCT, res2df$P_CCT)
max(abs(res1df$Pvalue_CCT - res2df$P_CCT) / pmin(abs(res1df$Pvalue_CCT), abs(res2df$P_CCT)))
