#' WtCoxG method in GRAB package
#'
#' WtCoxG is an accurate, powerful, and computationally efficient Cox-based approach to perform genome-wide time-to-event data analyses in study cohorts with case ascertainment.
#'
#' @details
#' Additional arguments in \code{GRAB.NullModel()}:
#' \itemize{
#'   \item \code{RefAfFile}: A character string specifying a reference allele frequency file, which is a csv file (with a header) and includes columns of \code{CHROM}, \code{POS}, \code{ID}, \code{REF}, \code{ALT}, \code{AF_ref}, and \code{AN_ref}.
#'   \item \code{OutputFile}: A character string specifying the output file name.
#'   \item \code{SampleIDColumn}: A character string specifying the column name in the input data that contains sample IDs.
#'   \item \code{SurvTimeColumn}: A character string specifying the column name in the input data that contains survival time information.
#'   \item \code{IndicatorColumn}: A character string specifying the column name in the input data that indicates case-control status (should be 0 for controls and 1 for cases).
#' }
#'
#' Additional arguments in list \code{control} in \code{GRAB.NullModel()}:
#' \itemize{
#'   \item \code{RefPrevalence}: A numeric value specifying the population-level disease prevalence used for weighting in the analysis.
#'   \item \code{SNPnum}: Minimum number of SNPs. Default is 1e4.
#' }
#'
#' Additional arguments in list \code{control} in \code{GRAB.Marker()}:
#' \itemize{
#'   \item \code{cutoff}: A numeric value specifying the batch effect p-value cutoff for method selection of an assocition test. Default is 0.1.
#' }
#'
#' @examples
#' # Step0&1: fit a null model and estimate parameters according to batch effect p-values
#' PhenoFile <- system.file("extdata", "simuPHENO.txt", package = "GRAB")
#' PhenoData <- data.table::fread(PhenoFile, header = TRUE)
#' SparseGRMFile <- system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#'
#' GenoFile <- system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' RefAfFile <- system.file("extdata", "simuRefAf.txt", package = "GRAB")
#' RefPrevalence <- 0.1 # population-level disease prevalence
#'
#' OutputDir <- system.file("results", package = "GRAB")
#' OutputStep1 <- paste0(OutputDir, "/WtCoxG_step1_out.txt")
#' OutputStep2 <- paste0(OutputDir, "/WtCoxG_step2_out.txt")
#'
#' obj.WtCoxG <- GRAB.NullModel(
#'   formula = survival::Surv(SurvTime, SurvEvent) ~ AGE + GENDER,
#'   data = PhenoData,
#'   subjData = PhenoData$IID,
#'   method = "WtCoxG",
#'   traitType = "time-to-event",
#'   GenoFile = GenoFile,
#'   SparseGRMFile = SparseGRMFile,
#'   control = list(AlleleOrder = "ref-first", AllMarkers = TRUE, RefPrevalence = RefPrevalence,
#'                  SNPnum = 1000), # minimum number of SNPs for to call TestforBatchEffect
#'   RefAfFile = RefAfFile,
#'   OutputFile = OutputStep1,
#'   SampleIDColumn = "IID",
#'   SurvTimeColumn = "SurvTime",
#'   IndicatorColumn = "SurvEvent"
#' )
#'
#' resultStep1 <- data.table::fread(OutputStep1)
#' resultStep1[, c("CHROM", "POS", "pvalue_bat")]
#'
#' # Step2: conduct association testing
#' GRAB.Marker(
#'   objNull = obj.WtCoxG,
#'   GenoFile = GenoFile,
#'   OutputFile = OutputStep2,
#'   control = list(AlleleOrder = "ref-first", AllMarkers = TRUE,
#'                  cutoff = 0.1, nMarkersEachChunk = 5000)
#' )
#'
#' resultStep2 <- data.table::fread(OutputStep2)
#' resultStep2[, c("CHROM", "POS", "WtCoxG.noext", "WtCoxG.ext")]
#' @export
GRAB.WtCoxG <- function() {
  print("Check ?GRAB.WtCoxG for more details about 'WtCoxG' method.")
}


# check the control list in null model fitting for SPACox method
checkControl.NullModel.WtCoxG <- function(control, traitType, ...) {
  if (!traitType %in% c("time-to-event", "binary", "Residual")) {
    stop("For 'WtCoxG' method, only traitType of 'time-to-event', 'binary' or 'Residual' is supported.")
  }

  default.control <- list(RefPrevalence = 0)

  control <- updateControl(control, default.control)

  if (control$RefPrevalence <= 0 | control$RefPrevalence >= 0.5) {
    stop("control$Refprevalence is required and should be between (0, 0.5).")
  }

  return(control)
}


getWeight.WtCoxG <- function(Indicator, RefPrevalence) {
  # BWJ (2023-08-08): Check NA later.
  # if(any(!unique(Indicator) %in% c(0, 1, NA)))
  #   stop("The value of Indicator should be 0, 1, or NA")

  if (any(!unique(Indicator) %in% c(0, 1))) {
    stop("The value of Indicator should be 0 or 1.")
  }

  sumOnes <- sum(Indicator)
  sumZeros <- sum(1 - Indicator)
  ratio <- sumOnes / sumZeros

  weight <- ifelse(Indicator == 1, 1,
    (1 - RefPrevalence) / RefPrevalence * ratio
  )
  return(weight)
}


# fit null model using WtCoxG method
fitNullModel.WtCoxG <- function(response, designMat, subjData,
                                control = list(OutlierRatio = 1.5), ...) {
  if (!(inherits(response, "Surv") || inherits(response, "Residual"))) { # later support binary trait and Residual
    stop("For WtCoxG, the response variable should be of class 'Surv' or 'Residual'.")
  }

  if (inherits(response, "Surv")) {
    formula <- response ~ designMat

    Indicator <- as.matrix(response)[, "status"]
    RefPrevalence <- control$RefPrevalence

    weight <- getWeight.WtCoxG(Indicator, RefPrevalence)

    obj.coxph <- survival::coxph(formula, x = T, weight = weight, robust = T, ...)

    y <- obj.coxph$y
    yVec <- y[, ncol(y)] # status

    mresid <- obj.coxph$residuals
    # Cova = obj.coxph$x
    Cova <- designMat
  } else {
    stop("We only support 'time-to-event' trait for WtCoxG by 2023-08-08.")
  }

  # if(class(response) == "Residual")
  # {
  #   yVec = mresid = response
  #   Cova = designMat
  # }

  # var.resid = var(mresid, na.rm = T)
  ## outliers or not depending on the residuals (0.25%-1.5IQR, 0.75%+1.5IQR)
  q25 <- quantile(mresid, 0.25, na.rm = T)
  q75 <- quantile(mresid, 0.75, na.rm = T)
  IQR <- q75 - q25

  # r.outlier = 1.5    # put this to the control argument later
  r.outlier <- ifelse(is.null(control$OutlierRatio), 1.5, control$OutlierRatio) # put this to the control argument later
  cutoff <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
  posOutlier <- which(mresid < cutoff[1] | mresid > cutoff[2])
  cat("The r.outlier is:", r.outlier, "\n")
  while (length(posOutlier) == 0) {
    r.outlier <- r.outlier * 0.8
    cutoff <- c(q25 - r.outlier * IQR, q75 + r.outlier * IQR)
    posOutlier <- which(mresid < cutoff[1] | mresid > cutoff[2])
    cat("The curent outlier ratio is:", r.outlier, "\n")
    cat("The number of outlier is:", length(posOutlier), "\n")
  }

  # The original code is from SPAmix in which multiple residuals were analysis simultaneously
  # posValue = which(!is.na(mresid))
  # posNonOutlier = setdiff(posValue, posOutlier)
  posValue <- 1:length(mresid)
  posNonOutlier <- setdiff(posValue, posOutlier)

  cat("The outlier of residuals will be passed to SPA analysis.\n")
  cat("Cutoffs to define residuals:\t", signif(cutoff, 2), "\n")
  cat("Totally, ", length(posOutlier), "/", length(posValue), " are defined as outliers.\n")

  if (length(posOutlier) == 0) {
    stop("No outlier is observed. SPA is not required in this case.")
  }

  # "-1" is to convert R style (index starting from 1) to C++ style (index starting from 0)
  outLierList <- list(
    posOutlier = posOutlier - 1,
    posNonOutlier = posNonOutlier - 1,
    resid = mresid,
    resid2 = mresid^2,
    residOutlier = mresid[posOutlier],
    residNonOutlier = mresid[posNonOutlier],
    resid2NonOutlier = mresid[posNonOutlier]^2
  )

  re <- list(
    mresid = mresid, Cova = Cova, yVec = yVec, weight = weight,
    RefPrevalence = RefPrevalence,
    N = length(mresid),
    outLierList = outLierList
  )

  class(re) <- "WtCoxG_NULL_Model"
  return(re)
}


#' Quality control to check batch effect between study cohort and reference population.
#'
#' This function performs quality control to test for the batch effect between a study cohort and a reference population. And fit a weighted null model.
#'
#' @param objNull a \code{WtCoxG_NULL_Model} object, which is the output of \code{\link{GRAB.NullModel}}.
#' @param data a data.frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data.frame), containing the variables in formula. Neither a matrix nor an array will be accepted.
#' @param GenoFile A character string of the genotype file. See Details section for more details.
#' @param GenoFileIndex Additional index file(s) corresponding to GenoFile. See Details section for more details.
#' @param Geno.mtx A matrix of genotype data. If provided, it will be used instead of GenoFile. The matrix should have samples in rows and markers in columns.
#' @param SparseGRMFile a path to file of output to be passed to \code{\link{GRAB.NullModel}}.
#' @param RefAfFile A character string of the reference file. The reference file must be a \code{txt} file (header required) including at least 7 columns: \code{CHROM}, \code{POS}, \code{ID}, \code{REF}, \code{ALT}, \code{AF_ref}, \code{AN_ref}.
#' @param OutputFile A character string of the output file name. The output file will be a \code{txt} file.
#' @param IndicatorColumn A character string of the column name in \code{data} that indicates the case-control status. The value should be 0 for controls and 1 for cases.
#' @param SurvTimeColumn A character string of the column name in \code{data} that indicates the survival time.
#' @param SampleIDColumn A character string of the column name in \code{data} that indicates the sample ID.
#' @return A dataframe of marker info and reference MAF.
#' @export
TestforBatchEffect <- function(
    objNull,
    data,
    GenoFile = NULL, # a character of file names of genotype files
    GenoFileIndex = NULL, # additional index file(s) corresponding to GenoFile
    Geno.mtx = NULL, # genotype matrix, if provided, will be used instead of GenoFile
    SparseGRMFile = NULL, # sparse genotype relatedness matrix
    RefAfFile, # header should include c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref","AN_ref")
    OutputFile,
    IndicatorColumn,
    SurvTimeColumn,
    SampleIDColumn) {
  if (!is.null(SparseGRMFile)) {
    sparseGRM <- data.table::fread(SparseGRMFile)
  } else {
    sparseGRM <- NULL
  }

  control <- objNull$control
  RefPrevalence <- control$RefPrevalence # refernce population prevalence, the proportion of indicator == 1.

  if (is.null(control$SNPnum)) {
    SNPnum <- 1e4 # default least number of SNPs is 1e4
  } else {
    SNPnum <- control$SNPnum
  }

  posCol <- which(colnames(data) == IndicatorColumn)
  colnames(data)[posCol] <- "Indicator"

  posCol <- which(colnames(data) == SurvTimeColumn)
  colnames(data)[posCol] <- "SurvTime"

  posCol <- which(colnames(data) == SampleIDColumn)
  colnames(data)[posCol] <- "SampleID"


  # step1: quality control--------------------------------------------------------
  ## reference genoInfo----------------------------------------------------------
  refGenoInfo <- data.table::fread(RefAfFile) %>% as_tibble()

  # check if there are 7 columns in RefAfFile
  for (colname in c("CHROM", "POS", "ID", "REF", "ALT", "AF_ref", "AN_ref")) {
    if (!colname %in% colnames(refGenoInfo)) {
      stop(paste0(colname, " is missing in RefAfFile!"))
    }
  }

  ## merge sample genoInfo and ref genoInfo--------------------------------------
  if(is.null(Geno.mtx)){
    GenoInfo.ctrl <- GRAB.getGenoInfo(
      GenoFile = GenoFile,
      GenoFileIndex = GenoFileIndex,
      SampleIDs = with(data, SampleID[Indicator == 0]), # MAF in cases
      control = control
    ) %>%
      rename(mu0 = altFreq, mr0 = missingRate) %>%
      select(mu0, mr0)

    if (nrow(GenoInfo.ctrl) < SNPnum) {
      stop("The number of genetic variants < ", SNPnum)
    }

    GenoInfo <- GRAB.getGenoInfo(
      GenoFile = GenoFile,
      GenoFileIndex = GenoFileIndex,
      SampleIDs = with(data, SampleID[Indicator == 1]), # MAF in controls
      control = control
    ) %>%
      rename(mu1 = altFreq, mr1 = missingRate) %>%
      cbind(., GenoInfo.ctrl) %>%
      as_tibble() %>%
      mutate(RA = paste0(pmin(REF, ALT), pmax(REF, ALT))) %>%
      mutate(index = 1:n())

    mergeGenoInfo <- refGenoInfo %>%
      mutate(RA = paste0(pmin(REF, ALT), pmax(REF, ALT))) %>%
      merge(., GenoInfo, by = c("CHROM", "POS", "RA"), all.y = T, sort = F) %>%
      rename(REF = REF.y, ALT = ALT.y, ID = ID.y) %>%
      mutate(AF_ref = ifelse(REF == REF.x, AF_ref, 1 - AF_ref)) %>%
      select(-REF.x, -ALT.x, -ID.x, -RA) %>%
      filter(!duplicated(index)) %>%
      mutate(
        n1 = sum(data$Indicator) * (1 - mr1),
        n0 = sum(1 - data$Indicator) * (1 - mr0),
        mu.int = 0.5 * mu1 + 0.5 * mu0,
        mu.int = ifelse(mu.int > 0.5, 1 - mu.int, mu.int)
      )
  }else{
    GenoInfo <- data.frame(
      ID = colnames(Geno.mtx),
      mu0 = apply(Geno.mtx[data$Indicator==0,], 2 , function(x){mean(na.omit(x)/2)}),
      mu1 = apply(Geno.mtx[data$Indicator==1,], 2 ,  function(x){mean(na.omit(x)/2)}),
      n0 = apply(Geno.mtx[data$Indicator==0,], 2 , function(x){sum(!is.na(x))}),
      n1 = apply(Geno.mtx[data$Indicator==1,], 2 , function(x){sum(!is.na(x))})) %>%
      mutate(mu.int = 0.5*mu1 + 0.5*mu0,
             mu.int = ifelse(mu.int>0.5, 1-mu.int, mu.int),
             index= 1:n()
      )
    mergeGenoInfo <- merge(GenoInfo, refGenoInfo, by = "ID", all.x=T, sort=F ) %>%
      filter( !duplicated(index) )
  }

  #### calculate batch effect p-value for each genetic variant------------------------------------
  w1 <- objNull$weight / (2 * sum(objNull$weight))
  names(w1) <- data$SampleID
  meanR <- mean(objNull$mresid)
  R_tilde <- objNull$mresid - meanR
  names(R_tilde) <- data$SampleID

  if (!is.null(sparseGRM)) {
    sparseGRM <- sparseGRM %>% mutate(
      cov = Value * w1[as.character(ID1)] * w1[as.character(ID2)],
      cov_R = Value * R_tilde[as.character(ID1)] * R_tilde[as.character(ID2)]
    )
    var.ratio.w0 <- (sum(sparseGRM$cov) + 1 / (2 * mergeGenoInfo$AN_ref)) / (sum(w1^2) + 1 / (2 * mergeGenoInfo$AN_ref))
    var.ratio.int <- sum(sparseGRM$cov_R) / sum(R_tilde^2)
  } else {
    var.ratio.w0 <- var.ratio.int <- 1
  }

  mergeGenoInfo <- mergeGenoInfo %>%
    mutate(
      var.ratio.w0 = var.ratio.w0,
      var.ratio.int = var.ratio.int
    )

  pvalue_bat <- lapply(1:nrow(mergeGenoInfo), function(ind) {
    p.test <- Batcheffect.Test(
      n0 = mergeGenoInfo$n0[ind],
      n1 = mergeGenoInfo$n1[ind],
      n.ext = mergeGenoInfo$AN_ref[ind] / 2,
      maf0 = mergeGenoInfo$mu0[ind],
      maf1 = mergeGenoInfo$mu1[ind],
      maf.ext = mergeGenoInfo$AF_ref[ind],
      pop.prev = RefPrevalence,
      var.ratio = mergeGenoInfo$var.ratio.w0[ind]
    )
  }) %>% unlist()

  mergeGenoInfo <- mergeGenoInfo %>% mutate(pvalue_bat)
  rm(pvalue_bat)

  #### estimate unknown parameters according to batch effect p-values---------------------------------
  cat("Estimate TPR and sigma2--------------\n")
  maf.group <- c(seq(-1e-4, 0.4, 0.05), max(mergeGenoInfo$mu.int))
  mergeGenoInfo <- lapply(1:(length(maf.group) - 1), function(i) {
    cat(i, "\n")

    ## assume that genotypes with MAF in [ maf.group[i] , maf.group[i+1]] have the same mixture distribution
    mergeGenoInfo_1 <- mergeGenoInfo %>% filter(mu.int > maf.group[i] & mu.int <= maf.group[i + 1])

    ## using batcheffect p-values with MAF in [maf.group[i]-0.1 , maf.group[i+1]+0.1] to estimate parameters
    mergeGenoInfo_2 <- mergeGenoInfo %>%
      filter(mu.int >= max(maf.group[i] - 0.1, 0) & mu.int < min(1, maf.group[i + 1] + 0.1))

    mu <- (maf.group[i] + maf.group[i + 1]) / 2

    n.ext <- mean(na.omit(mergeGenoInfo_1$AN_ref)[1]) / 2
    var_mu_ext <- mu * (1 - mu) / (2 * n.ext)

    var_Sbat <- ifelse(is.null(sparseGRM), sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext,
      na.omit(mergeGenoInfo$var.ratio.w0)[1] * (sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext)
    )

    obj <- fun.est.param(
      vec_p_bat = mergeGenoInfo_2$pvalue_bat,
      vec_var_Sbat = var_Sbat
    )
    TPR <- obj[1]
    sigma2 <- obj[2]

    w.ext <- optim(
      par = 0.5, method = "L-BFGS-B", lower = 0, upper = 1,
      fn = fun.optimalWeight,
      pop.prev = RefPrevalence,
      y = data$Indicator,
      R = objNull$mresid,
      w = objNull$weight,
      mu = mu,
      N = nrow(data),
      n.ext = n.ext,
      sigma2 = obj$sigma2,
      TPR = obj$TPR
    )$par[1]

    if (is.null(sparseGRM)) {
      var.ratio.ext <- 1
    } else {
      R_tilde_w <- objNull$mresid - mean(objNull$mresid) * w.ext
      names(R_tilde_w) <- data$SampleID
      sparseGRM <- sparseGRM %>% mutate(cov_Rext = Value * R_tilde_w[as.character(ID1)] * R_tilde_w[as.character(ID2)])
      var.ratio.ext <- (sum(sparseGRM$cov_Rext) + w.ext^2 * sum(objNull$mresid)^2 / n.ext) / (sum(R_tilde_w^2) + w.ext^2 * sum(objNull$mresid)^2 / n.ext)
    }

    mergeGenoInfo_1 <- mergeGenoInfo_1 %>% cbind(., TPR, sigma2, w.ext, var.ratio.ext)
  }) %>%
    do.call("rbind", .) %>%
    as_tibble() %>%
    arrange(index) %>%
    select(-index)

  #### output-----------------------------------------------------------
  data.table::fwrite(
    mergeGenoInfo,
    OutputFile,
    row.names = F, col.names = T, quote = F, sep = "\t"
  )
  return(mergeGenoInfo)
}


## function to test for batch effect--------------------------------------------
#' Test for batch effect
#'
#' This function test for the allele frequency difference between the study cohort and the external datasets.
#'
#' @param n0 A numeric. The sample size of cases in the study cohort
#' @param n1 A numeric. The sample size of controls in the study cohort
#' @param n.ext A numeric. The sample size of external datasets
#' @param maf0 A numeric. The MAF of the cases.
#' @param maf1 A numeric. The MAF of the controls
#' @param maf.ext A numeric. The MAF of the external datasets.
#' @param pop.prev A numeric. The population prevalence of the disease.
#' @param var.ratio A numeric. The variance ratio calculated by sparseGRM.
#' @return A numeric of batch effect p-value
Batcheffect.Test <- function(n0, # number of controls
                             n1, # number of cases
                             n.ext, # number of external dataset
                             maf0, # estimated MAF in controls
                             maf1, # estimated MAF in cases
                             maf.ext, # estimated MAF in external dataset
                             pop.prev,
                             var.ratio = 1) {
  er <- n1 / (n1 + n0)
  w0 <- (1 - pop.prev) / pop.prev / ((1 - er) / er)
  w1 <- 1

  weight.maf <- sum(maf0 * w0 * n0 + maf1 * w1 * n1) / sum(w0 * n0 + w1 * n1) ## weighted mean of genotypes
  est.maf <- sum(maf0 * w0 * n0 + maf1 * w1 * n1 + maf.ext * n.ext * w0) / sum(n1 * w1 + n0 * w0 + n.ext * w0) ## MAF estimates

  v <- ((n1 * w1^2 + n0 * w0^2) / (2 * (n1 * w1 + n0 * w0)^2) + 1 / (2 * n.ext)) * est.maf * (1 - est.maf) ## variance of test statistics
  z <- (weight.maf - maf.ext) / sqrt(v) ## standardized statistics
  z.adj <- z / sqrt(var.ratio) ## adjusted statistics by variance ratio
  p <- 2 * pnorm(-abs(z.adj), lower.tail = TRUE)

  return(p)
}

# estimate TPR and sigma2----------------------------------------------------
fun.est.param <- function(vec_p_bat,
                          vec_var_Sbat,
                          vec_cutoff = seq(0.01, 0.4, 0.1)) {
  ########  step1: the proportion of p_bat>cutoff
  vec_p_deno <- lapply(vec_cutoff, function(p_cut) {
    p_deno <- mean(na.omit(vec_p_bat > p_cut))
  }) %>% unlist()

  ######## optimization function
  opti_fun <- function(var_Sbat, vec_p_deno, par) {
    diff <- lapply(1:length(vec_cutoff), function(j) {
      p_cut <- vec_cutoff[j]
      lb <- -qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
      ub <- qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
      p_deno <- vec_p_deno[j]

      c <- pnorm(ub, 0, sqrt(var_Sbat + par[2]), log.p = T)
      d <- pnorm(lb, 0, sqrt(var_Sbat + par[2]), log.p = T)

      pro.cut <- par[1] * (exp(d) * (exp(c - d) - 1)) + (1 - par[1]) * (1 - p_cut)
      t <- ((p_deno - pro.cut) / p_deno)^2
    }) %>% do.call("sum", .)

    return(diff)
  }

  ####### estimate TPR and sigma2 for each SNP
  var.diff <- lapply(1:length(vec_var_Sbat), function(i) {
    if (i %% 100 == 0) cat(i, "\n")

    obj <- optim(
      par = c(0.01, 0.01),
      # ,method = "SANN"
      # , lower = 0, upper = 1
      fn = opti_fun,
      vec_p_deno = vec_p_deno,
      var_Sbat = vec_var_Sbat[i]
    )
    TPR <- min(1, max(0, obj$par[1]))
    sigma2 <- min(1, max(0, obj$par[2]))
    return(cbind(TPR, sigma2))
  }) %>%
    do.call("rbind", .) %>%
    as_tibble()

  return(var.diff)
}

# optimal weight for mu_ext -----------------------------------------------
fun.optimalWeight <- function(par, pop.prev, R, y, mu1, w, mu, N, n.ext, sigma2, TPR) {
  b <- par[1]

  p.fun <- function(b, pop.prev, R, y, mu1, mu, w, N, n.ext, sigma2, TPR) {
    meanR <- mean(R)
    sumR <- sum(R)

    mu0 <- mu
    mu.pop <- mu1 * pop.prev + mu0 * (1 - pop.prev)

    mu.i <- ifelse(y == 1, 2 * mu1, 2 * mu0)

    S <- sum((R - (1 - b) * meanR) * mu.i) - sumR * 2 * b * mu.pop

    w1 <- w / (2 * sum(w))
    mu <- mean(mu.i) / 2

    var_mu_ext <- mu * (1 - mu) / (2 * n.ext)
    var_Sbat <- sum(w1^2) * 2 * mu * (1 - mu) + var_mu_ext

    p_cut <- 0.1
    lb <- -qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
    ub <- qnorm(1 - p_cut / 2) * sqrt(var_Sbat)
    c <- pnorm(ub, 0, sqrt(var_Sbat + sigma2), log.p = T)
    d <- pnorm(lb, 0, sqrt(var_Sbat + sigma2), log.p = T)
    p_deno <- TPR * (exp(d) * (exp(c - d) - 1)) + (1 - TPR) * (1 - p_cut)

    ## sigma2=0
    var.int <- sum((R - (1 - b) * meanR)^2) * 2 * mu * (1 - mu)
    var_S <- var.int + 4 * b^2 * sumR^2 * var_mu_ext
    cov_Sbat_S <- sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) + 2 * b * sumR * var_mu_ext
    VAR <- matrix(c(var_S, cov_Sbat_S, cov_Sbat_S, var_Sbat), nrow = 2)
    p0 <- max(0, mvtnorm::pmvnorm(lower = c(-Inf, lb), upper = c(-abs(S), ub), mean = c(0, 0), sigma = VAR))

    ## sigma2!=0
    var_S1 <- var.int + 4 * b^2 * sumR^2 * (var_mu_ext + sigma2)
    cov_Sbat_S1 <- sum(w1 * (R - (1 - b) * meanR)) * 2 * mu * (1 - mu) + 2 * b * sumR * (var_mu_ext + sigma2)
    var_Sbat1 <- var_Sbat + sigma2
    VAR1 <- matrix(c(var_S1, cov_Sbat_S1, cov_Sbat_S1, var_Sbat1), nrow = 2)
    p1 <- max(0, mvtnorm::pmvnorm(lower = c(-Inf, lb), upper = c(-abs(S), ub), mean = c(0, 0), sigma = VAR1))

    p.con <- 2 * (TPR * p1 + (1 - TPR) * p0) / p_deno
    # diff = -log10(p.con)+log10(5e-8)
    diff <- -log10(p.con / 5e-8)

    return(diff)
  }

  mu1 <- uniroot(p.fun,
    lower = mu, upper = 1,
    b = b, pop.prev = pop.prev, mu = mu,
    R = R, y = y, w = w, N = N, n.ext = n.ext, sigma2 = sigma2, TPR = TPR
  )$root

  return(mu1)
}


checkControl.Marker.WtCoxG <- function(control) {
  return(control)
}

setMarker.WtCoxG <- function(objNull, control) {
  ImputeMethod <- if (is.null(control$ImputeMethod)) "none" else control$ImputeMethod
  cutoff <- if (is.null(control$cutoff)) 0.1 else control$cutoff
  setWtCoxGobjInCPP(
    objNull$mresid,
    objNull$weight,
    ImputeMethod,
    cutoff
  )
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk, "."))
}

mainMarker.WtCoxG <- function(genoType, genoIndex, control, objNull) {
  mergeGenoInfo <- objNull$mergeGenoInfo
  mergeGenoInfo_subset <- mergeGenoInfo[mergeGenoInfo$genoIndex %in% genoIndex, ]

  # Use match to reorder, but check for missing values
  match_indices <- match(genoIndex, mergeGenoInfo_subset$genoIndex)
  if (any(is.na(match_indices))) {
    missing_indices <- genoIndex[is.na(match_indices)]
    stop(paste0(
      "Missing marker info for genoIndex values: ",
      paste(missing_indices[1:min(10, length(missing_indices))], collapse = ", "),
      if (length(missing_indices) > 10) " ..." else ""
    ))
  }

  mergeGenoInfo_subset <- mergeGenoInfo_subset[match_indices, ]

  # Safety check: ensure sizes match
  if (nrow(mergeGenoInfo_subset) != length(genoIndex)) {
    stop(paste0(
      "Size mismatch: genoIndex has ", length(genoIndex),
      " elements, but mergeGenoInfo_subset has ", nrow(mergeGenoInfo_subset),
      " rows. Some markers in genoIndex may not be found in mergeGenoInfo."
    ))
  }

  # Update WtCoxG object with marker information for current chunk
  updateWtCoxGChunkInCPP(mergeGenoInfo_subset)

  OutList <- mainMarkerInCPP("WtCoxG", genoType, genoIndex)
  pvals <- data.frame(matrix(OutList$pvalVec, ncol = 2, byrow = TRUE))
  colnames(pvals) <- c("WtCoxG.ext", "WtCoxG.noext")

  obj.mainMarker <- cbind(pvals, mergeGenoInfo_subset)
  return(obj.mainMarker)
}
