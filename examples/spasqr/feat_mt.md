#

## The R/Rcpp implementation, mtMarker() is in the feat/mt branch

```R
PhenoData <- data.table::fread("examples/simuPHENO.txt")

y = PhenoData$QuantPheno,
X = as.matrix(PhenoData[, c("AGE", "GENDER", "PC1", "PC2")]),
subjData = PhenoData$IID,
SparseGRMFile = "examples/SparseGRM.txt",
taus  = c(0.1, 0.3, 0.5, 0.7, 0.9),
h     = 0,
sqr_tol = 1e-7

# ========== Fit null model ==========
y[is.na(y)] <- median(y, na.rm = TRUE)
ntaus    <- length(taus)
ResidMat <- matrix(0, nrow = length(y), ncol = ntaus)

if (h == 0) h <- IQR(y) / 3
for (i in seq_along(taus)) {
fit <- conquer::conquer(X, y, tau = taus[i], kernel = "Gaussian", h = h, tol = sqr_tol)
resid <- as.numeric(y - fit$coeff[1] - X %*% fit$coeff[2:(ncol(X) + 1)])
ResidMat[, i] <- taus[i] - pnorm((-resid) / h)
}

# ========== Identify outliers ==========
Quant       <- apply(ResidMat, 2, quantile, probs = c(0.25, 0.75), na.rm = TRUE)
Range       <- Quant[2, ] - Quant[1, ]
cutoffLower <- Quant[1, ] - 1.5 * Range
cutoffUpper <- Quant[2, ] + 1.5 * Range
cutoffLower <- ifelse(cutoffLower < -0.55, -0.55, cutoffLower)
cutoffUpper <- ifelse(cutoffUpper >  0.55,  0.55, cutoffUpper)
tooSmall <- sweep(ResidMat, 2, cutoffLower, "<")
tooLarge <- sweep(ResidMat, 2, cutoffUpper, ">")
Outlier  <- tooSmall | tooLarge

outlier_prop <- colMeans(Outlier, na.rm = TRUE)
cat("Outlier ratio for each column in the residual matrix:\n")
print(round(outlier_prop, 2))

# ========== Load GRM (3 columns: ID1, ID2, Value; always includes diagonal) ==========
SparseGRM <- data.table::fread(SparseGRMFile)
# Assuming: always 3 columns (id1, id2, coef), id1/id2 are character, diagonal always present
id1 <- as.character(SparseGRM[[1]])
id2 <- as.character(SparseGRM[[2]])
val <- SparseGRM[[3]]

pos1   <- match(id1, subjData)
pos2   <- match(id2, subjData)
# off-diagonal entries stored once; factor=2 accounts for symmetric matrix
factor <- ifelse(id1 == id2, 1, 2)

# ========== Accumulate GRM-based variance terms ==========
R_GRM_R_vec              <- numeric(ntaus)
R_GRM_R_nonOutlier_vec   <- numeric(ntaus)
sum_R_nonOutlier_vec     <- numeric(ntaus)
Resid.unrelated.outliers_lst <- lapply(seq_len(ntaus), function(x) numeric(0))

for (i in seq_along(taus)) {
R_col   <- ResidMat[, i]
Out_col <- Outlier[, i]

contrib        <- factor * val * R_col[pos1] * R_col[pos2]
R_GRM_R_vec[i] <- sum(contrib)

both_non_out           <- !Out_col[pos1] & !Out_col[pos2]
R_GRM_R_nonOutlier_vec[i] <- sum(contrib * both_non_out)

sum_R_nonOutlier_vec[i]           <- sum(R_col[!Out_col])
Resid.unrelated.outliers_lst[[i]] <- R_col[Out_col]
}

obj <- list(
taus                         = taus,
Resid_mat                    = ResidMat,
subjData                     = subjData,
N                            = length(subjData),
R_GRM_R_vec                  = R_GRM_R_vec,
R_GRM_R_TwoSubjOutlier_vec   = numeric(ntaus),
sum_R_nonOutlier_vec         = sum_R_nonOutlier_vec,
R_GRM_R_nonOutlier_vec       = R_GRM_R_nonOutlier_vec,
Resid.unrelated.outliers_lst = Resid.unrelated.outliers_lst,
TwoSubj_list_lst             = lapply(seq_len(ntaus), function(x) list()),
CLT_union_lst                = list(),
ThreeSubj_family_idx_lst     = lapply(seq_len(ntaus), function(x) integer(0)),
ThreeSubj_stand_S_lst        = lapply(seq_len(ntaus), function(x) list()),
MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
)
class(obj) <- "SPAsqr_NULL_Model"

GRAB.mtMarker(obj, "examples/simuPLINK", "examples/spasqr/r_out.tsv")

colnames(ResidMat) <- paste0("tau", taus)
df <- data.frame(IID = subjData, ResidMat)

write.table(
    df,
    "examples/spasqr/ResidMat.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
```

```sh
build/grab \
  --method SPAsqr \
  --null-resid examples/spasqr/ResidMat.tsv \
  --sp-grm-grab examples/SparseGRM.txt \
  --bfile examples/simuPLINK \
  --out examples/spasqr/cpp_out1.tsv


build/grab \
  --method SPAsqr \
  --pheno examples/simuPHENO.txt \
  --pheno-quant QuantPheno \
  --spasqr-taus 0.1,0.3,0.5,0.7,0.9 \
  --covar-name AGE,GENDER,PC1,PC2 \
  --sp-grm-grab examples/SparseGRM.txt \
  --bfile examples/simuPLINK \
  --out tmp/cpp_out2.tsv
```
