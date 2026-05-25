#!/usr/bin/env Rscript
# LEAF.R — Stand-alone reference R implementation of LEAF for the GRAB2
# `examples/` 1kg fixture.  Used to validate the C++ port in src/wtcoxg/leaf.cpp.
#
# Reads (all under EXAMPLES_DIR):
#   1kg.{pgen,pvar,psam}     genotypes (plink2 native)
#   1kg.pheno                Binary phenotype + covariates (MALE, PC1..4)
#   ref_pop{1,2}.afreq       reference allele frequencies (two populations)
#   1kg.grm.{sp,id}          plink2 sparse GRM and the IID lookup table
#   leaf_kmeans_remove       subjects excluded by GRAB2's --remove
#
# Writes (under OUT_DIR):
#   1kg.Binary.LEAF.tsv      per-marker GWAS output (GRAB2-matched layout,
#                            without the MISS_RATE and HWE_P columns)
#   1kg.KmeansCluster.tsv    subject-level K-means cluster assignments
#                            (same #IID\tcluster format as GRAB2 emits)
#   1kg.traw                 intermediate genotype matrix (plink2 --export Av
#                            with --remove applied at the export step)
#
# Dependencies:
#   - R packages: data.table (fread/fwrite + columnar mutation), mvtnorm,
#     nloptr.  Workflow logic is otherwise pure base R; no dplyr / tibble.
#   - plink2 in PATH (used once to convert 1kg.pgen → 1kg.traw)
#
# Assumptions:
#   - Script is run from the GRAB project root so the relative paths below
#     ("examples/...", "tmp/...") resolve correctly.
#   - The 1kg fixture has no missing genotypes, so no per-genotype NA
#     imputation is performed.
#
# Usage:  Rscript LEAF.R [n.cpu]
# n.cpu is consumed only by Summix (default 1).

# ── IO configuration (edit paths here) ─────────────────────────────────────
# Paths are resolved relative to the current working directory; run this
# script from the GRAB project root (the directory containing examples/).
EXAMPLES_DIR = "examples"
OUT_DIR      = "tmp"
PLINK2       = "plink2"

pheno_path   = file.path(EXAMPLES_DIR, "1kg.pheno")
pvar_path    = file.path(EXAMPLES_DIR, "1kg.pvar")
pgen_prefix  = file.path(EXAMPLES_DIR, "1kg")
ref_paths    = c(pop1 = file.path(EXAMPLES_DIR, "ref_pop1.afreq"),
                 pop2 = file.path(EXAMPLES_DIR, "ref_pop2.afreq"))
grm_path     = file.path(EXAMPLES_DIR, "1kg.grm.sp")
grm_id_path  = file.path(EXAMPLES_DIR, "1kg.grm.id")
remove_path  = file.path(EXAMPLES_DIR, "leaf_kmeans_remove")

pheno_col    = "Binary"
covar_rhs    = "MALE + PC1 + PC2 + PC3 + PC4"
pc_cols      = c("PC1", "PC2", "PC3", "PC4")
prev         = 0.1
ncluster     = 3
seed         = 1        # chosen so R's stats::kmeans(nstart=25) produces the
                        # same cluster labels (no permutation) as GRAB2's
                        # K-means++ / Lloyd implementation seeded with the
                        # same value via --seed 2026 in examples/run.sh.
pop_names    = names(ref_paths)
af_cols      = paste0("AF_", pop_names)
an_cols      = paste0("AN_", pop_names)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
out_main     = file.path(OUT_DIR, paste0("1kg.", pheno_col, ".LEAF.tsv"))
out_cluster  = file.path(OUT_DIR, "1kg.KmeansCluster.tsv")
traw_path    = file.path(OUT_DIR, "1kg.traw")

args  = commandArgs(TRUE)
n.cpu = ifelse(length(args) >= 1, as.numeric(args[1]), 1)

library(data.table)
library(mvtnorm)

# ── plink2: 1kg.pgen → 1kg.traw, applying --remove at the export step ──────
# All --remove filtering happens here; the R code below operates on the
# subject set actually written to the .traw, and never re-applies the list.
if (!file.exists(traw_path)) {
    cmd = paste(PLINK2, "--pfile", pgen_prefix,
                "--remove", remove_path,
                "--export Av --out", file.path(OUT_DIR, "1kg"))
    cat("[plink2] ", cmd, "\n")
    rc = system(cmd)
    if (rc != 0L)
        stop("plink2 failed (rc=", rc, "); confirm PLINK2 binary in PATH.")
}

# ─────────────────────────────────────────────────────────────────────────────
# Section A — Vendored Summix (Hendricks lab; subset that LEAF actually
# calls: summix(), summix_calc(), calc_scaledObj()).  Inlined verbatim
# from Summix.R to remove the source() dependency.
# ─────────────────────────────────────────────────────────────────────────────

summix <- function(data,
                   reference,
                   observed,
                   pi.start = NA,
                   goodness.of.fit = TRUE,
                   override_removeSmallRef = FALSE,
                   network = FALSE,
                   N_reference = NA,
                   reference_colors = NA) {
  start_time = Sys.time()
  if (network) {
    if (length(N_reference) != length(reference))
      stop("ERROR: Please make sure N_reference is the same length as reference")
  }
  if (length(reference_colors) != 1) {
    if (length(reference_colors) != length(reference))
      stop("ERROR: Please make sure reference_colors is the same length as reference")
  }
  if (!override_removeSmallRef) {
    globTest <- invisible(summix_calc(data = data, reference = reference,
                                      observed = observed, pi.start = pi.start))
    if (sum(globTest[1, 5:ncol(globTest)] < 0.01) > 0) {
      toremove <- which(globTest[1, 5:ncol(globTest)] < 0.01)
      reference <- reference[-toremove]
    }
    pi.start = NA
    sum_res <- summix_calc(data = data, reference = reference,
                           observed = observed)
  }
  sum_res <- summix_calc(data = data, reference = reference,
                         observed = observed, pi.start = pi.start)
  if (goodness.of.fit) {
    new_obj <- calc_scaledObj(data = data, observed = observed,
                              reference = reference, pi.start = pi.start)
    sum_res[1] <- new_obj
  }
  end_time = Sys.time()
  ttime = end_time - start_time
  sum_res[3] <- ttime
  if (goodness.of.fit) {
    if (sum_res[1] >= 0.5 & sum_res[1] < 1.5) {
      print(paste0("CAUTION: Objective/1000SNPs = ", round(sum_res[1], 4),
                   " which is within the moderate fit range"))
    } else if (sum_res[1] >= 1.5) {
      print(paste0("WARNING: Objective/1000SNPs = ", round(sum_res[1], 4),
                   " which is above the poor fit threshold"))
    }
  } else {
    if (sum_res[1] / nrow(data) / 1000 >= 0.5 &
        sum_res[1] / nrow(data) / 1000 < 1.5) {
      print(paste0("CAUTION: Objective/1000SNPs = ",
                   round(sum_res[1] / nrow(data) / 1000, 4),
                   " which is within the moderate fit range"))
    } else if (sum_res[1] / nrow(data) / 1000 >= 1.5) {
      print(paste0("WARNING: Objective/1000SNPs = ",
                   round(sum_res[1] / nrow(data) / 1000, 4),
                   " which is above the poor fit threshold"))
    }
  }
  return(sum_res)
}

summix_calc = function(data, reference, observed, pi.start = NA) {
  if (!is(object = data, class2 = "data.frame"))
    stop("ERROR: data must be a data.frame as described in the vignette")
  if (typeof(observed) != "character")
    stop("ERROR: 'observed' must be a character string for the column name of the observed group in data")
  if (!(observed %in% names(data)))
    stop("ERROR: 'observed' must be the column name of the observed group in data")
  if (typeof(reference) != "character")
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  if (all(reference %in% names(data)) == FALSE)
    stop("ERROR: 'reference' must be a vector of column names from data to be used in the reference")
  data <- as.data.frame(data)
  filteredNA <- length(which(is.na(data[[observed]] == TRUE)))
  keep_rows  <- which(!is.na(data[[observed]]))
  observed.b <- data.frame(data[keep_rows, observed, drop = FALSE])
  refmatrix  <- data.frame(data[keep_rows, reference, drop = FALSE])
  if (!is.na(pi.start)[1]) {
    if (is.numeric(pi.start) == FALSE)
      stop("ERROR: Please make sure pi.start is a numeric vector")
    if (length(pi.start) != length(reference))
      stop("ERROR: Please make sure pi.start is the same length as reference")
    if (all(pi.start > 0) == FALSE)
      stop("ERROR: Please make sure pi.start is a positive numeric vector")
    if (sum(pi.start) != 1)
      stop("ERROR: Please make sure pi.start sums to one")
    starting = pi.start
  } else {
    starting = rep(1 / ncol(refmatrix), ncol(refmatrix))
  }
  fn.refmix = function(x) {
    expected = x %*% t(as.matrix(refmatrix))
    minfunc = sum((expected - observed.b)**2)
    return(minfunc)
  }
  gr.refmix <- function(x) {
    gradfunc = x %*% t(as.matrix(refmatrix)) - observed.b
    gradvec <- apply(2 * refmatrix * t(gradfunc), 2, sum)
    return(gradvec)
  }
  heq.refmix = function(x) sum(x) - 1
  hin.refmix <- function(x) { h = numeric(ncol(refmatrix)); h = x; h }
  start_time = Sys.time()
  S = suppressMessages(nloptr::slsqp(starting, fn = fn.refmix, gr = gr.refmix,
                                     hin = hin.refmix, heq = heq.refmix))
  end_time = Sys.time()
  ttime = end_time - start_time
  d <- data.frame(matrix(ncol = length(reference) + 4, nrow = 1))
  colnames(d) <- c("goodness.of.fit", "iterations", "time", "filtered",
                   colnames(refmatrix))
  d[1] <- S$value
  d[2] <- S$iter
  d[3] <- ttime
  d[4] <- filteredNA
  d[5:(length(reference) + 4)] <- round(S$par[1:length(reference)], 6)
  return(d)
}

calc_scaledObj <- function(data, reference, observed, pi.start) {
  start_time <- Sys.time()
  multiplier <- c(5, 1.5, 1)
  data$obs <- data[[observed]]
  data$maf_obs <- ifelse(data$obs > 0.5, 1 - data$obs, data$obs)
  breaks <- c(0, 0.1, 0.3, 0.5)
  data$bin <- cut(data$maf_obs, breaks = breaks, include.lowest = TRUE,
                  right = FALSE, labels = c(1, 2, 3))
  for (b in 1:3) {
    subdata <- data[data$bin == b, ]
    subdata <- subdata[complete.cases(subdata), ]
    if (nrow(subdata) > 0) {
      res <- summix_calc(subdata, reference = reference, observed = observed,
                         pi.start)
      res$obj_adj <- res$goodness.of.fit / (nrow(subdata) / 1000)
      res$nSNPs <- nrow(subdata)
      res$bin <- b
      if (b == 1) sum_res <- res else sum_res <- rbind(sum_res, res)
    } else {
      res <- summix_calc(data, reference = reference, observed = observed,
                         pi.start)
      res$obj_adj <- NA
      res$nSNPs <- 0
      res$bin <- b
      if (b == 1) sum_res <- res else sum_res <- rbind(sum_res, res)
    }
  }
  objective_scaled <- 0
  nonNA <- 0
  for (b in 1:3) {
    if (!is.na(sum_res[b, ]$obj_adj)) {
      nonNA <- nonNA + 1
      objective_scaled <- objective_scaled + (sum_res[b, ]$obj_adj * multiplier[b])
    }
  }
  objective_scaled <- objective_scaled / nonNA
  end_time <- Sys.time()
  return(objective_scaled)
}

# ─────────────────────────────────────────────────────────────────────────────
# Section B — LEAF helpers (batch-effect test, parameter estimation,
# optimal external-weight search, WtCoxG SPA + MGF/CGF kernels).
# ─────────────────────────────────────────────────────────────────────────────

# Batch-effect testing using internal vs external MAF (case-control weighted)
Batcheffect.Test = function(n0, n1, n.ext, maf0, maf1, maf.ext,
                            pop.prev, var.ratio = 1) {
  er = n1 / (n1 + n0)
  w0 = (1 - pop.prev) / pop.prev / ((1 - er) / er)
  w1 = 1
  weight.maf = sum(maf0 * w0 * n0 + maf1 * w1 * n1) / sum(w0 * n0 + w1 * n1)
  est.maf    = sum(maf0 * w0 * n0 + maf1 * w1 * n1 + maf.ext * n.ext * w0) /
               sum(n1 * w1 + n0 * w0 + n.ext * w0)
  v =((n1 * w1^2 + n0 * w0^2) / (2 * (n1 * w1 + n0 * w0)^2) +
      1 / (2 * n.ext)) * est.maf * (1 - est.maf)
  z = (weight.maf - maf.ext) / sqrt(v)
  z.adj = z / sqrt(var.ratio)
  p = 2 * pnorm(-abs(z.adj), lower.tail = TRUE)
  return(p)
}

# Estimate (TPR, sigma2) for each MAF bin.  Reparameterized Nelder-Mead in
# unbounded η-space:  TPR = σ(η_T) ∈ (0,1);  σ² = exp(η_S) ∈ (0,∞).
fun.est.param = function(vec_p_bat, vec_var_Sbat,
                         vec_cutoff = seq(0.01, 0.4, 0.1)) {
  vec_p_deno = vapply(vec_cutoff, function(p_cut)
                      mean(na.omit(vec_p_bat > p_cut)),
                      numeric(1))
  opti_fun_eta = function(eta, var_Sbat, vec_p_deno) {
    TPR    = 1 / (1 + exp(-eta[1]))
    sigma2 = exp(eta[2])
    diff = 0
    for (j in seq_along(vec_cutoff)) {
      p_cut = vec_cutoff[j]
      lb = -qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
      ub =  qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
      p_deno = vec_p_deno[j]
      c = pnorm(ub, 0, sqrt(var_Sbat + sigma2), log.p = TRUE)
      d = pnorm(lb, 0, sqrt(var_Sbat + sigma2), log.p = TRUE)
      pro.cut = TPR * (exp(d) * (exp(c - d) - 1)) + (1 - TPR) * (1 - p_cut)
      diff = diff + ((p_deno - pro.cut) / p_deno)^2
    }
    diff
  }
  rows = lapply(seq_along(vec_var_Sbat), function(i) {
    if (i %% 100 == 0) cat(i, "\n")
    obj = optim(par = c(0, log(0.01)),
                fn = opti_fun_eta,
                vec_p_deno = vec_p_deno,
                var_Sbat   = vec_var_Sbat[i],
                method = "Nelder-Mead",
                control = list(reltol = 1e-10, maxit = 500))
    TPR    = 1 / (1 + exp(-obj$par[1]))
    sigma2 = exp(obj$par[2])
    c(TPR = TPR, sigma2 = sigma2)
  })
  as.data.frame(do.call(rbind, rows))
}

# Optimal external-data weight b for combining mu.ext into the score.
fun.optimalWeight = function(par, pop.prev, R, y, mu1, w, mu, N, n.ext,
                             sigma2, TPR) {
  b = par[1]
  p.fun = function(b, pop.prev, R, y, mu1, mu, w, N, n.ext, sigma2, TPR) {
    meanR = mean(R); sumR = sum(R)
    mu0    = mu
    mu.pop = mu1 * pop.prev + mu0 * (1 - pop.prev)
    mu.i   = ifelse(y == 1, 2 * mu1, 2 * mu0)
    S = sum((R - (1 - b) * meanR) * mu.i) - sumR * 2 * b * mu.pop
    w1 = w / (2 * sum(w))
    mu = mean(mu.i) / 2
    var_mu_ext = mu * (1 - mu) / (2 * n.ext)
    var_Sbat   = sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext
    p_cut = 0.1
    lb = -qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
    ub =  qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
    c = pnorm(ub, 0, sqrt(var_Sbat + sigma2), log.p = TRUE)
    d = pnorm(lb, 0, sqrt(var_Sbat + sigma2), log.p = TRUE)
    p_deno = TPR * (exp(d) * (exp(c - d) - 1)) + (1 - TPR) * (1 - p_cut)
    var.int = sum((R - (1 - b) * meanR)^2) * 2 * mu * (1 - mu)
    var_S = var.int + 4 * b^2 * sumR^2 * var_mu_ext
    cov_Sbat_S = sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) +
                 2 * b * sumR * var_mu_ext
    VAR = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat), nrow = 2)
    p0 = max(0, pmvnorm(lower = c(-Inf, lb), upper = c(-abs(S), ub),
                        mean = c(0, 0), sigma = VAR))
    var_S1 = var.int + 4 * b^2 * sumR^2 * (var_mu_ext + sigma2)
    cov_Sbat_S1 = sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) +
                  2 * b * sumR * (var_mu_ext + sigma2)
    var_Sbat1 = var_Sbat + sigma2
    VAR1 = matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1), nrow = 2)
    p1 = max(0, pmvnorm(lower = c(-Inf, lb), upper = c(-abs(S), ub),
                        mean = c(0, 0), sigma = VAR1))
    p.con = 2 * (TPR * p1 + (1 - TPR) * p0) / p_deno
    diff = -log10(p.con / 5e-8)
    return(diff)
  }
  mu1 = uniroot(p.fun, lower = mu, upper = 1,
                b = b, pop.prev = pop.prev, mu = mu,
                R = R, y = y, w = w, N = N, n.ext = n.ext,
                sigma2 = sigma2, TPR = TPR)$root
  return(mu1)
}

# ── MGF/CGF kernels for the SPA p-value (homogeneous-genotype SPA) ─────────
M_G0 = function(t, MAF) (1 - MAF + MAF * exp(t))^2
M_G1 = function(t, MAF) 2 * (MAF * exp(t)) * (1 - MAF + MAF * exp(t))
M_G2 = function(t, MAF) 2 * (MAF * exp(t))^2 + 2 * (MAF * exp(t)) *
                         (1 - MAF + MAF * exp(t))
K_G0 = function(t, MAF) log(M_G0(t, MAF))
K_G1 = function(t, MAF) M_G1(t, MAF) / M_G0(t, MAF)
K_G2 = function(t, MAF) M_G2(t, MAF) / M_G0(t, MAF) -
                        (M_G1(t, MAF) / M_G0(t, MAF))^2

H_org = function(t, R, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est,
                 meanR, b) {
  n.t = length(t); out = rep(0, n.t)
  mu.adj  = -2 * b * sumR * MAF
  var.adj =  4 * b^2 * sumR^2 * var_mu_ext
  for (i in 1:n.t) {
    t1 = t[i]
    out[i] = sum(K_G0(t1 * (R - (1 - b) * meanR), MAF)) +
             mu.adj * t1 + var.adj / 2 * t1^2
  }
  return(out)
}

H1_adj = function(t, R, s, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est,
                  meanR, b) {
  n.t = length(t); out = rep(0, n.t)
  mu.adj  = -2 * b * sumR * MAF
  var.adj =  4 * b^2 * sumR^2 * var_mu_ext
  for (i in 1:n.t) {
    t1 = t[i]
    out[i] = sum((R - (1 - b) * meanR) *
                 K_G1(t1 * (R - (1 - b) * meanR), MAF)) +
             mu.adj + var.adj * t1 - s
  }
  return(out)
}

H2 = function(t, R, MAF, n.ext, N.all, sumR, var_mu_ext, g.var.est,
              meanR, b) {
  n.t = length(t); out = rep(0, n.t)
  R_hat = sumR / N.all
  var.adj = n.ext * R_hat^2 * 2 * MAF * (1 - MAF)
  for (i in 1:n.t) {
    t1 = t[i]
    out[i] = sum((R - (1 - b) * meanR)^2 *
                 K_G2(t1 * (R - (1 - b) * meanR), MAF)) + var.adj
  }
  return(out)
}

GetProb_SPA_G = function(MAF, R, s, n.ext, N.all, sumR, var_mu_ext,
                         g.var.est, meanR, b, lower.tail) {
  out  = uniroot(H1_adj, c(-1, 1), extendInt = "yes",
                 R = R, s = s, MAF = MAF, n.ext = n.ext, N.all = N.all,
                 sumR = sumR, var_mu_ext = var_mu_ext,
                 g.var.est = g.var.est, meanR = meanR, b = b)
  zeta = out$root
  k1 = H_org(zeta, R = R, MAF = MAF, n.ext = n.ext, N.all = N.all,
             sumR = sumR, var_mu_ext = var_mu_ext, g.var.est = g.var.est,
             meanR = meanR, b = b)
  k2 = H2(zeta, R = R, MAF = MAF, n.ext = n.ext, N.all = N.all,
          sumR = sumR, var_mu_ext = var_mu_ext, g.var.est = g.var.est,
          meanR = meanR, b = b)
  temp1 = zeta * s - k1
  w = sign(zeta) * (2 * temp1)^{1 / 2}
  v = zeta * (k2)^{1 / 2}
  pval = pnorm(w + 1 / w * log(v / w), lower.tail = lower.tail)
  return(pval)
}

# Saddle-point approximation for homogeneous (b = 0 or σ² = 0) baseline p-values
SPA_G.one.SNP_homo = function(g, R, mu.ext = NA, n.ext = NA, b = 0,
                              sigma2 = NA, var.ratio = 1, Cutoff = 2,
                              min.mac = 10, G.model = "Add") {
  if (is.na(mu.ext)) { mu.ext = 0; n.ext = 0 }
  if (sum(g) < min.mac | sum(2 - g) < min.mac) return(c(NA, NA))
  N = length(g)
  mu.int = mean(g) / 2
  MAF    = (1 - b) * mu.int + b * mu.ext
  sumR   = sum(R)
  N.all  = N + n.ext
  S = sum(R * (g - 2 * MAF))
  S = S / var.ratio
  g.var.est  = 2 * MAF * (1 - MAF)
  var_mu_ext = ifelse(n.ext == 0, 0, MAF * (1 - MAF) / (2 * n.ext) + sigma2)
  meanR = mean(R)
  S.var = sum((R - (1 - b) * meanR)^2) * g.var.est +
          4 * b^2 * sumR^2 * var_mu_ext
  z = S / sqrt(S.var)
  if (abs(z) < Cutoff) {
    pval.norm = pnorm(abs(z), lower.tail = FALSE) * 2
    return(c(pval.norm, pval.norm))
  }
  pval1 = GetProb_SPA_G(MAF, R = R,  abs(S), n.ext = n.ext, N.all = N.all,
                        var_mu_ext = var_mu_ext, g.var.est = g.var.est,
                        meanR = meanR, sumR = sumR, b = b, lower.tail = FALSE)
  pval2 = GetProb_SPA_G(MAF, R = R, -abs(S), n.ext = n.ext, N.all = N.all,
                        var_mu_ext = var_mu_ext, g.var.est = g.var.est,
                        meanR = meanR, sumR = sumR, b = b, lower.tail = TRUE)
  pval3 = pnorm( abs(z), lower.tail = FALSE)
  pval4 = pnorm(-abs(z), lower.tail = TRUE)
  return(c(pval1 + pval2, pval3 + pval4))
}

# Main WtCoxG p-value for one SNP within one cluster.  Branches on whether
# a usable external MAF is available and whether the batch-effect p-value
# clears the cutoff.
WtCoxG.test = function(g, R, w, p_bat,
                       TPR = NA, sigma2 = NA, b = 0,
                       var.ratio.int = 1, var.ratio.w0 = 1, var.ratio.w1 = 1,
                       var.ratio0 = 1, var.ratio1 = 1,
                       mu.ext = NA, n.ext = NA, p_cut = 0.1) {
  if (is.na(mu.ext)) {
    mu.int = mean(g) / 2
    p.con  = SPA_G.one.SNP_homo(g = g, R = R, mu.ext = NA, n.ext = 0,
                                sigma2 = 0, var.ratio = var.ratio.int)[1]
    S = sum(R * (g - mean(g)))
    return(cbind(p.con, S))
  }
  if (is.na(p_bat) | sum(g) < 10 | sum(2 - g) < 10) {
    return(cbind(p.con = NA, S = NA))
  } else if (p_bat < p_cut) {
    R_tilde = R - mean(R)
    mu = mu.int = mean(g) / 2
    S = sum(R * (g - 2 * mu.int))
    w1 = w / (2 * sum(w))
    var_mu_ext = mu * (1 - mu) / (2 * n.ext)
    var_Sbat   = sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext
    lb = -qnorm(1 - p_cut / 2) * sqrt(var_Sbat) * sqrt(var.ratio.w0)
    ub = -lb
    c = pnorm(lb / sqrt(var.ratio.w1), 0, sqrt(var_Sbat + sigma2),
              log.p = TRUE)
    p_deno = TPR * 2 * exp(c) + (1 - TPR) * p_cut
    p_spa_s0 = SPA_G.one.SNP_homo(g = g, R = R, var.ratio = var.ratio0)[1]
    var_S = S^2 / var.ratio0 / qchisq(p_spa_s0, 1, ncp = 0, lower.tail = FALSE)
    var.int = sum(R_tilde^2) * 2 * mu * (1 - mu)
    cov_Sbat_S = sum(w1 * R_tilde) * 2 * mu * (1 - mu)
    cov_Sbat_S = cov_Sbat_S * sqrt(var_S / var.int)
    VAR  = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat), nrow = 2)
    p0 = max(0, min(1, pmvnorm(lower = c(-Inf, -Inf),
                               upper = c(-abs(S / sqrt(var.ratio0)),
                                         lb / sqrt(var.ratio.w0)),
                               mean = c(0, 0), sigma = VAR))) +
         max(0, min(1, pmvnorm(lower = c(-Inf, ub / sqrt(var.ratio.w0)),
                               upper = c(-abs(S / sqrt(var.ratio0)), Inf),
                               mean = c(0, 0), sigma = VAR)))
    var_Sbat1 = var_Sbat + sigma2
    VAR1 = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat1), nrow = 2)
    p1 = max(0, min(1, pmvnorm(lower = c(-Inf, -Inf),
                               upper = c(-abs(S / sqrt(var.ratio1)),
                                         lb / sqrt(var.ratio.w1)),
                               mean = c(0, 0), sigma = VAR1))) +
         max(0, min(1, pmvnorm(lower = c(-Inf, ub / sqrt(var.ratio.w1)),
                               upper = c(-abs(S / sqrt(var.ratio1)), Inf),
                               mean = c(0, 0), sigma = VAR1)))
    p.con = 2 * (TPR * p1 + (1 - TPR) * p0) / p_deno
    return(cbind(p.con, S))
  } else {
    meanR = mean(R); sumR = sum(R)
    mu.int = mean(g) / 2
    N = length(g)
    mu = (1 - b) * mu.int + b * mu.ext
    S  = sum(R * (g - 2 * mu))
    w1 = w / (2 * sum(w))
    var_mu_ext = mu * (1 - mu) / (2 * n.ext)
    var_Sbat   = sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext
    lb = -qnorm(1 - p_cut / 2) * sqrt(var_Sbat) * sqrt(var.ratio.w0)
    ub =  qnorm(1 - p_cut / 2) * sqrt(var_Sbat) * sqrt(var.ratio.w0)
    c = pnorm(ub / sqrt(var.ratio.w1), 0, sqrt(var_Sbat + sigma2),
              log.p = TRUE)
    d = pnorm(lb / sqrt(var.ratio.w1), 0, sqrt(var_Sbat + sigma2),
              log.p = TRUE)
    p_deno = TPR * (exp(d) * (exp(c - d) - 1)) + (1 - TPR) * (1 - p_cut)
    p_spa_s0 = SPA_G.one.SNP_homo(g = g, R = R, b = b, mu.ext = mu.ext,
                                  n.ext = n.ext, sigma2 = 0,
                                  var.ratio = var.ratio0)[1]
    var_S = S^2 / var.ratio0 / qchisq(p_spa_s0, 1, ncp = 0, lower.tail = FALSE)
    var.int = sum((R - (1 - b) * meanR)^2) * 2 * mu * (1 - mu)
    cov_Sbat_S = sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) +
                 2 * b * sumR * var_mu_ext
    cov_Sbat_S = cov_Sbat_S *
                 sqrt(var_S / (var.int + 4 * b^2 * sumR^2 * var_mu_ext))
    VAR = matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat), nrow = 2)
    p0 = max(0, min(1, pmvnorm(lower = c(-Inf, lb / sqrt(var.ratio.w0)),
                               upper = c(-abs(S / sqrt(var.ratio0)),
                                         ub / sqrt(var.ratio.w0)),
                               mean = c(0, 0), sigma = VAR)))
    p_spa_s1 = SPA_G.one.SNP_homo(g = g, R = R, b = b, mu.ext = mu.ext,
                                  n.ext = n.ext, sigma2 = sigma2,
                                  var.ratio = var.ratio1)[1]
    var_S1 = S^2 / var.ratio1 / qchisq(p_spa_s1, 1, ncp = 0, lower.tail = FALSE)
    cov_Sbat_S1 = sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) +
                  2 * b * sumR * (var_mu_ext + sigma2)
    cov_Sbat_S1 = cov_Sbat_S1 *
                  sqrt(var_S1 / (var.int +
                                 4 * b^2 * sumR^2 * (var_mu_ext + sigma2)))
    var_Sbat1 = var_Sbat + sigma2
    VAR1 = matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1), nrow = 2)
    p1 = max(0, min(1, pmvnorm(lower = c(-Inf, lb / sqrt(var.ratio.w1)),
                               upper = c(-abs(S / sqrt(var.ratio1)),
                                         ub / sqrt(var.ratio.w1)),
                               mean = c(0, 0), sigma = VAR1)))
    p.con = 2 * (TPR * p1 + (1 - TPR) * p0) / p_deno
    return(cbind(p.con, S))
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Section C — Workflow (phenotype, K-means, genotype, per-cluster LEAF,
# meta-analysis, GRAB2-shaped TSV).
# ─────────────────────────────────────────────────────────────────────────────

cat(sprintf("LEAF.R [1kg]  prev=%g  ncluster=%d  seed=%d  n.cpu=%g\n",
            prev, ncluster, seed, n.cpu))

# Phenotype + null model + K-means.  Subjects removed by --remove were
# dropped at the plink2 --export step, so the .traw genotype matrix already
# excludes them; restrict the phenotype table to the .traw subject set.
pheno = fread(pheno_path)
traw_header = scan(traw_path, what = character(), nlines = 1, sep = "\t",
                   quiet = TRUE)
keep_iids = sub("^[^_]+_", "", traw_header[-(1:6)])
pheno = pheno[as.character(`#IID`) %in% keep_iids]
y = pheno[[pheno_col]]
if (!all(y %in% c(0L, 1L)))
  stop("pheno column '", pheno_col,
       "' is not binary (expected 0/1 for the logistic null model).")
pheno$weight = ifelse(y == 1, 1, (1 - prev) / prev)
fmla = as.formula(paste(pheno_col, "~", covar_rhs))
obj.wglm = glm(fmla, data = pheno, weight = weight, family = "binomial")
pheno$R       = obj.wglm$y - obj.wglm$fitted.values
pheno$weight1 = pheno$weight / (2 * sum(pheno$weight))
set.seed(seed)
km = kmeans(as.matrix(pheno[, ..pc_cols]),
            centers = ncluster, nstart = 25)
pheno$cluster = km$cluster
fwrite(pheno[, .(`#IID`, cluster)], file = out_cluster, sep = "\t")

# Genotype matrix (ALT-dosage, NA where missing) — .traw orientation via .pvar
cat("read in geno (plink2 .traw + .pvar orientation alignment)\n")
traw = fread(traw_path)
pvar = fread(pvar_path)
setnames(pvar, c("#CHROM", "ID"), c("CHR", "SNP"))
GenoInfo = merge(traw[, .(CHR, SNP, COUNTED, ALT)],
                 pvar[, .(SNP, REF_pvar = REF, ALT_pvar = ALT)],
                 by = "SNP", sort = FALSE)
GenoInfo[, `:=`(REF   = REF_pvar,
                ALT   = ALT_pvar,
                newID = paste0(SNP, REF_pvar, ALT_pvar))]
GenoInfo = GenoInfo[match(traw$SNP, GenoInfo$SNP), ]
flip = GenoInfo$COUNTED == GenoInfo$REF
# .traw stores 6 meta columns (CHR, SNP, (C)M, POS, COUNTED, ALT) before
# the per-subject dosages; drop them as a single column slice.
Gmat = as.matrix(traw[, -(1:6)])
Gmat[flip, ] = 2 - Gmat[flip, ]
Geno = t(Gmat)
rownames(Geno) = sub("^[^_]+_", "", colnames(traw)[-(1:6)])
colnames(Geno) = GenoInfo$SNP
GenoInfo = GenoInfo[, .(CHR, SNP, REF, ALT, newID)]
# Carry POS through from pvar (needed for the GRAB2-shaped TSV).
GenoInfo$POS = pvar$POS[match(GenoInfo$SNP, pvar$SNP)]
rm(traw, Gmat)

# Reference allele frequencies (one table per population)
pop_tables = lapply(pop_names, function(nm) {
  d = fread(ref_paths[[nm]])
  d$newID = paste0(d$ID, d$REF, d$ALT)
  d = merge(d, GenoInfo, by = "newID", all.y = TRUE, sort = FALSE)
  d$ALT_FREQS = ifelse(d$REF.x == d$REF.y, d$ALT_FREQS, 1 - d$ALT_FREQS)
  d[, c("REF.x", "ALT.x", "ID", "#CHROM") := NULL]
  setnames(d, c("REF.y", "ALT.y"), c("REF", "ALT"))
  d
})
names(pop_tables) = pop_names

# Sparse GRM (.grm.sp 0-based indices + .grm.id IID table → (ID1, ID2, Value))
grm_id  = fread(grm_id_path, header = FALSE, col.names = c("FID", "IID"))
grm_raw = fread(grm_path,    header = FALSE,
                col.names = c("idx1", "idx2", "Value"))
sparseGRM = data.table(
  ID1   = grm_id$IID[grm_raw$idx1 + 1L],
  ID2   = grm_id$IID[grm_raw$idx2 + 1L],
  Value = grm_raw$Value
)
rm(grm_id, grm_raw)

# Per-marker overall QC (ALT_FREQ, MAC).  The 1kg fixture has no missing
# genotypes, so the counts are computed directly from sum(g) and length(g).
cat("compute per-marker overall QC (ALT_FREQ, MAC)\n")
Geno_kept    = Geno[as.character(pheno$`#IID`), , drop = FALSE]
nKept        = nrow(Geno_kept)
altCount_vec = colSums(Geno_kept)
altFreq_vec  = altCount_vec / (2 * nKept)
mac_vec      = pmin(altCount_vec, 2 * nKept - altCount_vec)

# Per-cluster WtCoxG run (unchanged logic from the original LEAF.R)
gwas_list = lapply(1:ncluster, function(cl) {

  cat("cluster", cl, "\n")
  PhenoData = pheno[pheno$cluster == cl, ]
  G    = Geno[pheno$`#IID`[which(pheno$cluster == cl)], ]
  y_cl = PhenoData[[pheno_col]]

  # mergeGenoInfo carries one row per marker; start from GenoInfo (data.table)
  # and tack on the per-cluster MAF columns.  `:=` updates by reference.
  mergeGenoInfo = copy(GenoInfo)
  mergeGenoInfo[, `:=`(
    mu0       = colMeans(G[y_cl == 0, ], na.rm = TRUE) / 2,
    mu1       = colMeans(G[y_cl == 1, ], na.rm = TRUE) / 2
  )]
  mergeGenoInfo[, `:=`(
    mu.target = 0.5 * mu0 + 0.5 * mu1,
    mu.w      = prev * mu1 + (1 - prev) * mu0
  )]
  mergeGenoInfo[, mu.int := pmin(mu.target, 1 - mu.target)]
  for (nm in pop_names) {
    mergeGenoInfo[[paste0("AF_", nm)]] = pop_tables[[nm]]$ALT_FREQS
    mergeGenoInfo[[paste0("AN_", nm)]] = pop_tables[[nm]]$OBS_CT
  }

  # Summix: solve mixture proportions, then collapse populations → AF_ref/AN_ref
  K_pop    = length(pop_names)
  anc_prop = summix(data      = na.omit(mergeGenoInfo),
                    reference = af_cols,
                    observed  = "mu.target",
                    pi.start  = rep(1 / K_pop, K_pop))
  anc_prop = unlist(anc_prop)[af_cols]
  anc_prop = ifelse(is.na(anc_prop), 0, anc_prop)
  names(anc_prop) = af_cols
  af_mat   = as.matrix(mergeGenoInfo[, ..af_cols])
  an_vec   = vapply(an_cols, function(nm) mergeGenoInfo[[nm]][1], numeric(1))
  mergeGenoInfo$AF_ref = as.numeric(af_mat %*% anc_prop)
  mergeGenoInfo$AN_ref = 1 / sum(anc_prop^2 / an_vec)

  # Variance ratios via the cluster-restricted sparse GRM
  sparseGRMsub = sparseGRM[ID1 %in% PhenoData$`#IID` &
                           ID2 %in% PhenoData$`#IID`, ]
  cat("variance ratio ---------\n")
  if (nrow(sparseGRMsub) > 0L) {
    w1 = PhenoData$weight1
    names(w1) = PhenoData$`#IID`
    R_tilde = PhenoData$R - mean(PhenoData$R)
    names(R_tilde) = PhenoData$`#IID`
    sym = ifelse(sparseGRMsub$ID1 == sparseGRMsub$ID2, 1, 2)
    sparseGRMsub[, `:=`(
      cov   = sym * Value * w1[as.character(ID1)] * w1[as.character(ID2)],
      cov_R = sym * Value * R_tilde[as.character(ID1)] * R_tilde[as.character(ID2)]
    )]
    var.ratio.w0  = (sum(sparseGRMsub$cov)   + 1 / (2 * mergeGenoInfo$AN_ref)) /
                    (sum(w1^2)               + 1 / (2 * mergeGenoInfo$AN_ref))
    var.ratio.int =  sum(sparseGRMsub$cov_R) / sum(R_tilde^2)
  } else {
    var.ratio.w0  = 1
    var.ratio.int = 1
  }
  mergeGenoInfo$var.ratio.w0  = var.ratio.w0
  mergeGenoInfo$var.ratio.int = var.ratio.int

  # Batch-effect p-value per SNP
  mergeGenoInfo$pvalue_bat = vapply(seq_len(nrow(mergeGenoInfo)), function(i) {
    if (i %% 1000 == 0) cat(i, "\n")
    n1 = sum(y_cl); n0 = sum(1 - y_cl)
    Batcheffect.Test(n0 = n0, n1 = n1,
                     n.ext = mergeGenoInfo$AN_ref[i] / 2,
                     maf0 = mergeGenoInfo$mu0[i],
                     maf1 = mergeGenoInfo$mu1[i],
                     maf.ext = mergeGenoInfo$AF_ref[i],
                     pop.prev = prev,
                     var.ratio = mergeGenoInfo$var.ratio.w0[i])
  }, numeric(1))
  mergeGenoInfo$index = seq_len(nrow(mergeGenoInfo))

  # Per (MAF bin) TPR/sigma2 estimation
  cat("Estimate TPR and sigma2--------------\n")
  maf.group = c(seq(-0.00001, 0.4, 0.05), max(mergeGenoInfo$mu.int))
  bin_rows = lapply(seq_len(length(maf.group) - 1), function(i) {
    cat(i, "\n")
    data     = mergeGenoInfo[mergeGenoInfo$mu.int >  maf.group[i] &
                             mergeGenoInfo$mu.int <= maf.group[i + 1], ]
    data.ref = mergeGenoInfo[mergeGenoInfo$mu.int >= max(maf.group[i] - 0.1, 0) &
                             mergeGenoInfo$mu.int <  min(1, maf.group[i + 1] + 0.1), ]
    mu = (maf.group[i] + maf.group[i + 1]) / 2
    n.ext = na.omit(data$AN_ref)[1] / 2
    var_mu_ext = mu * (1 - mu) / (2 * n.ext)
    var_Sbat = if (is.null(sparseGRM)) {
                 sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext
               } else {
                 na.omit(mergeGenoInfo$var.ratio.w0)[1] *
                   (sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext)
               }
    obj = fun.est.param(vec_p_bat = data.ref$pvalue_bat,
                        vec_var_Sbat = var_Sbat)
    TPR = obj[1]; sigma2 = obj[2]
    w.ext = optim(par = 0.5, method = "L-BFGS-B", lower = 0, upper = 1,
                  fn = fun.optimalWeight,
                  pop.prev = prev, y = y_cl,
                  R = PhenoData$R, w = PhenoData$weight,
                  mu = mu, N = nrow(PhenoData),
                  n.ext = n.ext,
                  sigma2 = obj$sigma2, TPR = obj$TPR)$par[1]
    if (nrow(sparseGRMsub) == 0L) {
      var.ratio.ext = 1
    } else {
      R_tilde_w = PhenoData$R - mean(PhenoData$R) * w.ext
      names(R_tilde_w) = PhenoData$`#IID`
      sym_ext = ifelse(sparseGRMsub$ID1 == sparseGRMsub$ID2, 1, 2)
      sparseGRMsub[, cov_Rext := sym_ext * Value *
                                 R_tilde_w[as.character(ID1)] *
                                 R_tilde_w[as.character(ID2)]]
      var.ratio.ext = (sum(sparseGRMsub$cov_Rext) +
                       w.ext^2 * sum(PhenoData$R)^2 / n.ext) /
                      (sum(R_tilde_w^2) +
                       w.ext^2 * sum(PhenoData$R)^2 / n.ext)
    }
    cbind(data, TPR = TPR, sigma2 = sigma2,
                w.ext = w.ext, var.ratio.ext = var.ratio.ext)
  })
  mergeGenoInfo = as.data.frame(do.call(rbind, bin_rows))
  mergeGenoInfo = mergeGenoInfo[order(mergeGenoInfo$index), ]
  mergeGenoInfo$index = NULL

  # Per-SNP WtCoxG (ext + noext)
  gwas_rows = lapply(seq_len(ncol(G)), function(i) {
    if (i %% 1000 == 0) cat("Complete ", i, "/", ncol(G), "\n")
    g = G[, i]; R = PhenoData$R; w = PhenoData$weight
    mu.ext = mergeGenoInfo$AF_ref[i]; n.ext = mergeGenoInfo$AN_ref[i] / 2
    TPR    = mergeGenoInfo$TPR[i];     sigma2 = mergeGenoInfo$sigma2[i]
    p_bat  = mergeGenoInfo$pvalue_bat[i]; w.ext = mergeGenoInfo$w.ext[i]
    var.ratio.w0  = mergeGenoInfo$var.ratio.w0[i]
    var.ratio.int = mergeGenoInfo$var.ratio.int[i]
    var.ratio0    = mergeGenoInfo$var.ratio.ext[i]

    WtCoxG.ext = WtCoxG.test(g = g, R = R, w = w,
                             TPR = TPR, sigma2 = sigma2, b = w.ext,
                             var.ratio.w0 = var.ratio.w0,
                             var.ratio.w1 = var.ratio.w0,
                             var.ratio0   = var.ratio0,
                             var.ratio1   = var.ratio0,
                             mu.ext = mu.ext, n.ext = n.ext, p_bat = p_bat)
    noext_out = WtCoxG.test(g = g, R = R, w = w,
                            TPR = TPR, sigma2 = sigma2, b = w.ext,
                            var.ratio.int = var.ratio.int,
                            var.ratio.w0  = var.ratio.w0,
                            var.ratio.w1  = var.ratio.w0,
                            var.ratio0    = var.ratio0,
                            var.ratio1    = var.ratio0,
                            mu.ext = NA, n.ext = NA, p_bat = p_bat)
    p.noext = noext_out[1, "p.con"]
    S.noext = noext_out[1, "S"]

    # Per-cluster MAC = min(sum(g), 2N − sum(g)); no missing genotypes here.
    gSum_clu = sum(g); N_clu = length(g)
    mac_clu  = min(gSum_clu, 2 * N_clu - gSum_clu)

    cbind(WtCoxG.ext, p.noext = p.noext, S.noext = S.noext,
          mac_clu = mac_clu)
  })
  cbind(as.data.frame(do.call(rbind, gwas_rows)), mergeGenoInfo)
})

# Meta-analysis (recover Var from S and p, then Z = ΣS / √ΣVar)
all_S_ext = do.call(cbind, lapply(gwas_list, function(d) d$S))
all_P_ext = do.call(cbind, lapply(gwas_list, function(d) d$p.con))
all_S_no  = do.call(cbind, lapply(gwas_list, function(d) d$S.noext))
all_P_no  = do.call(cbind, lapply(gwas_list, function(d) d$p.noext))

meta_p_from_S_and_P = function(S_mat, P_mat) {
  # Mirror src/util/meta_pvalue.hpp: symmetric p-value clamp before
  # back-recovering Var_c = S_c^2 / qchisq(p_c, 1, lower.tail = FALSE).
  P_FLOOR = 1e-300
  P_CEIL  = 1 - 1e-15
  na_mask = is.na(S_mat) | is.na(P_mat)
  P_safe  = pmin(pmax(P_mat, P_FLOOR), P_CEIL)
  ChiSq   = qchisq(P_safe, df = 1, lower.tail = FALSE)
  Var     = (S_mat^2) / ChiSq
  # Treat NaN-score / NaN-p clusters as absent (S = 0, Var = 0 contribute
  # nothing to either sum) — the same semantics as the C++ `continue`.
  S_eff = S_mat; S_eff[na_mask] = 0
  V_eff = Var;   V_eff[na_mask] = 0
  sum_S = rowSums(S_eff)
  sum_V = rowSums(V_eff)
  Z = ifelse(sum_V > 0, sum_S / sqrt(sum_V), NA_real_)
  2 * pnorm(-abs(Z))
}
meta_P_EXT   = meta_p_from_S_and_P(all_S_ext, all_P_ext)
meta_P_NOEXT = meta_p_from_S_and_P(all_S_no,  all_P_no)

# ── Assemble GRAB2-shaped TSV (MISS_RATE / HWE_P dropped) ──────────────────
out = data.table(
  CHROM        = GenoInfo$CHR,
  POS          = GenoInfo$POS,
  ID           = GenoInfo$SNP,
  REF          = GenoInfo$REF,
  ALT          = GenoInfo$ALT,
  ALT_FREQ     = altFreq_vec,
  MAC          = mac_vec,
  meta_P_EXT   = meta_P_EXT,
  meta_P_NOEXT = meta_P_NOEXT
)
for (cl in 1:ncluster) {
  d = gwas_list[[cl]]
  out[[paste0("cl", cl, "_MAC")]]     = d$mac_clu
  out[[paste0("cl", cl, "_P_EXT")]]   = d$p.con
  out[[paste0("cl", cl, "_P_NOEXT")]] = d$p.noext
  out[[paste0("cl", cl, "_P_BAT")]]   = d$pvalue_bat
  out[[paste0("cl", cl, "_PI_BAT")]]  = d$TPR
  out[[paste0("cl", cl, "_VAR_BAT")]] = d$sigma2
}
fwrite(out, file = out_main, sep = "\t")
cat(sprintf("wrote %s  (%d markers × %d columns)\n",
            out_main, nrow(out), ncol(out)))
