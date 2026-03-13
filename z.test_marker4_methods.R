devtools::load_all(quiet = TRUE)

library(data.table)

pheno_file <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
pheno_data <- fread(pheno_file, header = TRUE)
pheno_data$OrdinalPheno <- factor(pheno_data$OrdinalPheno, levels = c(0, 1, 2))

long_file <- system.file("extdata", "simuLongPHENO.txt", package = "GRAB")
long_data <- fread(long_file, header = TRUE)

geno_file <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
sparse_grm_file <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
pairwise_ibd_file <- system.file("extdata", "PairwiseIBD.txt", package = "GRAB")
resid_mat_file <- system.file("extdata", "ResidMat.txt", package = "GRAB")
ref_af_file <- system.file("extdata", "simuRefAf.txt", package = "GRAB")

marker4_control <- list(nMarkersEachChunk = 100)

results <- data.frame(
  method = character(),
  status = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

record_result <- function(method_name, status, detail = "") {
  results <<- rbind(
    results,
    data.frame(method = method_name, status = status, detail = detail, stringsAsFactors = FALSE)
  )
}

compare_marker_paths <- function(method_name, obj_null, marker_control = NULL) {
  marker_out <- file.path(tempdir(), paste0("marker_", method_name, ".txt"))
  marker4_out <- file.path(tempdir(), paste0("marker4_", method_name, ".txt"))

  GRAB.Marker(
    objNull = obj_null,
    GenoFile = geno_file,
    OutputFile = marker_out,
    control = marker_control
  )

  GRAB.Marker4(
    objNull = obj_null,
    GenoFile = geno_file,
    OutputFile = marker4_out,
    control = utils::modifyList(marker4_control, marker_control %||% list()),
    nthreads = 4,
    overwrite = TRUE
  )

  marker_res <- fread(marker_out)
  marker4_res <- fread(marker4_out)
  same <- isTRUE(all.equal(marker_res, marker4_res, tolerance = 1e-10, check.attributes = FALSE))

  cat(sprintf("[%s] equal: %s\n", method_name, same))
  if (!same) {
    detail <- paste(all.equal(marker_res, marker4_res, tolerance = 1e-10, check.attributes = FALSE), collapse = " | ")
    record_result(method_name, "mismatch", detail)
  } else {
    record_result(method_name, "equal")
  }

  invisible(list(marker = marker_res, marker4 = marker4_res))
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

run_case <- function(method_name, expr) {
  cat(sprintf("===== %s =====\n", method_name))
  tryCatch(
    force(expr),
    error = function(e) {
      msg <- conditionMessage(e)
      cat(sprintf("[%s] error: %s\n", method_name, msg))
      record_result(method_name, "error", msg)
      NULL
    }
  )
}

run_case("POLMM", {
  obj_polmm <- GRAB.NullModel(
    OrdinalPheno ~ AGE + GENDER,
    data = pheno_data,
    subjIDcol = "IID",
    method = "POLMM",
    traitType = "ordinal",
    GenoFile = geno_file,
    SparseGRMFile = sparse_grm_file
  )
  compare_marker_paths("POLMM", obj_polmm)
})

run_case("SPACox", {
  obj_spacox <- GRAB.NullModel(
    survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
    data = pheno_data,
    subjIDcol = "IID",
    method = "SPACox",
    traitType = "time-to-event"
  )
  compare_marker_paths("SPACox", obj_spacox)
})

run_case("SPAmix", {
  obj_spamix <- GRAB.NullModel(
    survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
    data = pheno_data,
    subjIDcol = "IID",
    method = "SPAmix",
    traitType = "time-to-event",
    control = list(PC_columns = "PC1,PC2")
  )
  compare_marker_paths("SPAmix", obj_spamix)
})

run_case("SPAmixPlus", {
  af_file_output <- file.path(tempdir(), "afModels_marker4.bin")
  SPAmixPlus.AF(
    GenoFile = geno_file,
    PCs = as.matrix(pheno_data[, c("PC1", "PC2", "PC3")]),
    subjData = pheno_data$IID,
    afFileOutput = af_file_output,
    control = list(afFilePrecision = "double")
  )

  residuals_spamixplus <- survival::coxph(
    survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER + PC1 + PC2,
    data = pheno_data
  )$residuals

  obj_spamixplus <- GRAB.NullModel(
    residuals_spamixplus ~ AGE + GENDER + PC1 + PC2,
    data = pheno_data,
    subjIDcol = "IID",
    SparseGRMFile = sparse_grm_file,
    method = "SPAmixPlus",
    traitType = "Residual",
    control = list(PC_columns = "PC1,PC2")
  )
  compare_marker_paths(
    "SPAmixPlus",
    obj_spamixplus,
    marker_control = list(afFilePath = af_file_output, afFilePrecision = "double")
  )
})

run_case("SPAGRM", {
  obj_spagrm <- SPAGRM.NullModel(
    ResidMatFile = resid_mat_file,
    SparseGRMFile = sparse_grm_file,
    PairwiseIBDFile = pairwise_ibd_file,
    control = list(ControlOutlier = FALSE)
  )
  compare_marker_paths("SPAGRM", obj_spagrm)
})

run_case("SAGELD", {
  nullmodel_sageld <- lme4::lmer(LongPheno ~ AGE + GENDER + (AGE | IID), data = long_data)
  obj_sageld <- SAGELD.NullModel(
    NullModel = nullmodel_sageld,
    UsedMethod = "SAGELD",
    PlinkFile = geno_file,
    SparseGRMFile = sparse_grm_file,
    PairwiseIBDFile = pairwise_ibd_file
  )
  compare_marker_paths("SAGELD", obj_sageld)
})

run_case("SPAsqr", {
  obj_spasqr <- GRAB.NullModel(
    QuantPheno ~ AGE + GENDER,
    data = pheno_data,
    subjIDcol = "IID",
    method = "SPAsqr",
    traitType = "quantitative",
    SparseGRMFile = sparse_grm_file,
    PairwiseIBDFile = pairwise_ibd_file,
    control = list(taus = c(0.2, 0.5, 0.8))
  )
  compare_marker_paths("SPAsqr", obj_spasqr)
})

print(results)