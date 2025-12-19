## ------------------------------------------------------------------------------
## PairwiseIBD.R
## Utilities for constructing sparse GRMs and pairwise IBD estimates via GCTA,
## and combining intermediate outputs into analysis-ready matrices.
##
## Functions:
##   getTempFilesFullGRM : Generate per-part GRM files using GCTA (parallelizable).
##   getSparseGRM        : Merge per-part outputs into a sparse GRM file.
##   getPairwiseIBD      : Compute/assemble pairwise IBD metrics (pa, pb, pc).
## ------------------------------------------------------------------------------

#' Make temporary files to be passed to function \code{\link{getSparseGRM}}.
#'
#' Make temporary files to be passed to function \code{\link{getSparseGRM}}.
#' We strongly suggest using parallel computing for different \code{partParallel}.
#' @param PlinkPrefix a path to PLINK files (without file extensions of bed/bim/fam).
#'   Note that the current version (gcta_1.93.1beta) of gcta software does not support
#'   different prefix names for bim, bed, and fam files.
#' @param nPartsGRM a numeric value (e.g. 250): \code{GCTA} software can split subjects
#'   to multiple parts. For UK Biobank data analysis, it is recommended to set \code{nPartsGRM=250}.
#' @param partParallel a numeric value (from 1 to \code{nPartsGRM}) to split all jobs
#'   for parallel computation.
#' @param gcta64File a path to \code{GCTA} program. GCTA can be downloaded from
#'   [link](https://yanglab.westlake.edu.cn/software/gcta/#Download).
#' @param tempDir a path to store temp files to be passed to \code{\link{getSparseGRM}}.
#'   This should be consistent to the input of \code{\link{getSparseGRM}}. Default is \code{tempdir()}.
#' @param subjData a character vector to specify subject IDs to retain (i.e. IID).
#'   Default is \code{NULL}, i.e. all subjects are retained in sparse GRM. If the number
#'   of subjects is less than 1,000, the GRM estimation might not be accurate.
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from PLINK files)
#'   to make sparse GRM. *(default=0.01)*
#' @param maxMissingGRM Maximal value of missing rate to select markers (from PLINK files)
#'   to make sparse GRM. *(default=0.1)*
#' @param threadNum Number of threads (CPUs) to use.
#' @details
#' \itemize{
#'   \item \code{Step 1}: Run \code{getTempFilesFullGRM} to get temporary files.
#'   \item \code{Step 2}: Run \code{\link{getSparseGRM}} to combine the temporary files
#'     to make a \code{SparseGRMFile} to be passed to \code{\link{GRAB.NullModel}}.
#' }
#'
#' @return A character string message indicating the completion status and location
#'   of the temporary files.
#' @examples
#' ## Please check help(getSparseGRM) for an example.
#'
getTempFilesFullGRM <- function(
  PlinkPrefix,
  nPartsGRM,
  partParallel,
  gcta64File,
  tempDir = NULL,
  subjData = NULL,
  minMafGRM = 0.01,
  maxMissingGRM = 0.1,
  threadNum = 8
) {
  bimFile <- paste0(PlinkPrefix, ".bim")
  bedFile <- paste0(PlinkPrefix, ".bed")
  famFile <- paste0(PlinkPrefix, ".fam")

  if (!file.exists(bimFile)) stop("Could not find bimFile or paste0(PlinkPrefix,'.bim')")
  if (!file.exists(bedFile)) stop("Could not find bedFile or paste0(PlinkPrefix,'.bed')")
  if (!file.exists(famFile)) stop("Could not find famFile or paste0(PlinkPrefix,'.fam')")

  if (is.null(tempDir)) {
    tempDir <- tempdir()
  }

  PlinkName <- basename(PlinkPrefix)
  tempFile <- paste0(
    tempDir, "/Plink-", PlinkName, "-autosome-minMaf-",
    minMafGRM, "-maxMissing-", maxMissingGRM
  )

  cmd <- paste(
    gcta64File,
    "--bfile", PlinkPrefix,
    "--out", tempFile,
    "--autosome",
    "--make-grm-part", nPartsGRM, partParallel,
    "--maf", minMafGRM,
    "--geno", maxMissingGRM,
    "--thread-num", threadNum
  )

  ## only retain parts of subjects
  if (!is.null(subjData)) {
    if (length(subjData) < 1000) {
      stop("length(subjData) < 1000, the MAF estimate might be inaccurate.")
    }

    subjFile <- paste0(tempDir, "/subjData.txt")
    famData <- read.table(famFile)
    posSubj <- match(subjData, famData$V2, 0)
    if (any(posSubj == 0)) {
      stop("All subjects in 'subjData' should be in IID column of 'famFile'.")
    }

    write.table(famData[posSubj, c(1, 2)],
      subjFile,
      row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    cmd <- paste(
      cmd,
      "--keep", subjFile
    )
  }

  system(cmd)

  message <- paste0("Temp files of Full GRM have been saved to ", tempFile)
  return(message)
}


#' Make a \code{SparseGRMFile} for \code{\link{GRAB.NullModel}}.
#'
#' If the sample size in analysis is greater than 100,000, we recommend using sparse GRM
#' (instead of dense GRM) to adjust for sample relatedness.
#' This function is to use \code{GCTA} ([link](https://cnsgenomics.com/software/gcta/#Overview))
#' to make a \code{SparseGRMFile} to be passed to function \code{\link{GRAB.NullModel}}.
#' This function can only support \code{Linux} and \code{PLINK} files as required by
#' \code{GCTA} software. To make a \code{SparseGRMFile}, two steps are needed.
#' Please check \code{Details} section for more details.
#' @param PlinkPrefix a path to PLINK binary files (without file extension).
#'   Note that the current version (gcta_1.93.1beta) of \code{GCTA} software does not support
#'   different prefix names for BIM, BED, and FAM files.
#' @param nPartsGRM a numeric value (e.g. 250): \code{GCTA} software can split subjects
#'   to multiple parts. For UK Biobank data analysis, it is recommended to set \code{nPartsGRM=250}.
#' @param SparseGRMFile a path to file of output to be passed to \code{\link{GRAB.NullModel}}.
#' @param tempDir a path to store temp files from \code{\link{getTempFilesFullGRM}}.
#'   This should be consistent to the input of \code{\link{getTempFilesFullGRM}}.
#'   Default is \code{tempdir()}.
#' @param relatednessCutoff a cutoff for sparse GRM, only kinship coefficient greater
#'   than this cutoff will be retained in sparse GRM. *(default=0.05)*
#' @param minMafGRM Minimal value of MAF cutoff to select markers (from PLINK files)
#'   to make sparse GRM. *(default=0.01)*
#' @param maxMissingGRM Maximal value of missing rate to select markers (from PLINK files)
#'   to make sparse GRM. *(default=0.1)*
#' @param rm.tempFiles a logical value indicating if the temp files generated in
#'   \code{\link{getTempFilesFullGRM}} will be deleted. *(default=FALSE)*
#' @details
#' \itemize{
#'   \item \code{Step 1}: Run \code{\link{getTempFilesFullGRM}} to save temporary files
#'     to \code{tempDir}.
#'   \item \code{Step 2}: Run \code{getSparseGRM} to combine the temporary files to make
#'     a \code{SparseGRMFile} to be passed to function \code{\link{GRAB.NullModel}}.
#' }
#' Users can customize parameters including \code{(minMafGRM, maxMissingGRM, nPartsGRM)},
#' but functions \code{\link{getTempFilesFullGRM}} and \code{getSparseGRM} should use the same ones.
#' Otherwise, package \code{GRAB} cannot accurately identify temporary files.
#'
#' # The following shows a typical workflow for creating a sparse GRM:
#'
#' \code{# Input data (We recommend setting nPartsGRM=250 for UKBB with N=500K):}
#'
#' \code{GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")}
#'
#' \code{PlinkPrefix = tools::file_path_sans_ext(GenoFile)}
#'
#' \code{nPartsGRM = 2}
#'
#' # Step 1: We strongly recommend parallel computing in high performance clusters (HPC).
#'
#' \code{# For Linux, get the file path of gcta64 by which command:}
#'
#' \code{gcta64File <- system("which gcta64", intern = TRUE)}
#'
#' \code{# For Windows, set the file path directly:}
#'
#' \code{gcta64File <- "C:\\\\path\\\\to\\\\gcta64.exe"}
#'
#' \code{# The temp outputs (may be large) will be in tempdir() by default:}
#'
#' \code{for(partParallel in 1:nPartsGRM) getTempFilesFullGRM(PlinkPrefix, nPartsGRM, }
#' \code{partParallel, gcta64File)}
#'
#' # Step 2: Combine files in Step 1 to make a SparseGRMFile
#'
#' \code{tempDir = tempdir()}
#'
#' \code{SparseGRMFile = file.path(tempDir, "SparseGRM.txt")}
#'
#' \code{getSparseGRM(PlinkPrefix, nPartsGRM, SparseGRMFile)}
#'
#' @return
#' A character string containing a message with the path to the output file
#' where the sparse Genetic Relationship Matrix (SparseGRM) has been stored.
#'
getSparseGRM <- function(
  PlinkPrefix,
  nPartsGRM,
  SparseGRMFile,
  tempDir = NULL,
  relatednessCutoff = 0.05,
  minMafGRM = 0.01,
  maxMissingGRM = 0.1,
  rm.tempFiles = FALSE
) {
  PlinkName <- basename(PlinkPrefix)
  # Number of digits needed for formatting part numbers (1-9: 1 digit, 10-99: 2 digits, etc.)
  nDigits <- floor(log10(nPartsGRM)) + 1

  AllIDs <- c()
  n0 <- 0
  SparseGRM <- c()

  if (is.null(tempDir)) {
    tempDir <- tempdir()
  }

  ## cycle for nPartsGRM
  for (i in 1:nPartsGRM) {
    .message("Processing GRM part %d/%d", i, nPartsGRM)

    tempFile <- paste0(
      tempDir, "/Plink-", PlinkName, "-autosome-minMaf-",
      minMafGRM, "-maxMissing-", maxMissingGRM
    )

    ## Three files generated by GCTA
    part_suffix <- paste0(".part_", nPartsGRM, "_", formatC(i, width = nDigits, flag = "0"))
    IDFile <- paste0(tempFile, part_suffix, ".grm.id")
    BinFile <- paste0(tempFile, part_suffix, ".grm.bin")

    ## read in the three files
    IDs <- read.table(IDFile, stringsAsFactors = FALSE)
    ID <- IDs$V2
    AllIDs <- c(AllIDs, ID)
    n1 <- n0 + length(ID)
    nData <- (n1 - n0) * (n0 + n1 + 1) / 2

    grm <- readBin(BinFile, n = nData, what = numeric(0), size = 4)

    pos <- which(grm > relatednessCutoff)
    value <- grm[pos]
    # ---- BEGIN inlined: getPairs ----
    ID <- AllIDs
    start <- c()
    end <- c()
    temp <- 0
    for (ii in (n0 + 1):n1) {
      start <- c(start, temp)
      temp <- temp + ii
      end <- c(end, temp)
    }

    pos.row <- sapply(pos, FUN = function(x) {
      min(which(x <= end))
    })
    pos.col <- pos - start[pos.row]

    ID1 <- ID[n0 + pos.row]
    ID2 <- ID[pos.col]
    pairs <- cbind.data.frame(ID1, ID2, value)
    # ---- END inlined: getPairs ----
    SparseGRM <- rbind(SparseGRM, pairs)

    n0 <- n1
  }

  colnames(SparseGRM) <- c("ID1", "ID2", "Value")
  write.table(SparseGRM, SparseGRMFile,
              row.names = FALSE, quote = FALSE, sep = "\t", col.names = TRUE)

  message <- paste("The SparseGRM has been stored in", SparseGRMFile)

  return(message)
}


#' Calculate Pairwise IBD (Identity By Descent)
#'
#' This function calculates pairwise IBD probabilities for related samples using
#' PLINK genotype data and sparse GRM. It follows the getSparseGRM() function workflow.
#'
#' @param PlinkPrefix Character. Path to PLINK file (without file extensions .bed/.bim/.fam).
#' @param SparseGRMFile Character. Path to sparse GRM file from getSparseGRM() function.
#' @param PairwiseIBDOutput Character. Output path to save pairwise IBD results.
#'   If NULL (default), return a data.frame instead of saving to file.
#' @param frqFile Character. Path to frequency file corresponding to PLINK file.
#'   If NULL (default), uses PlinkPrefix.frq.
#' @param tempDir Character. Directory to save temporary files. If NULL (default), uses tempdir().
#' @param maxSampleNums Integer. Maximum number of subjects' genotypes to read for analysis (default: 2500).
#' @param minMafIBD Numeric. Minimum MAF cutoff to select markers (default: 0.01).
#' @param rm.tempFile Logical. Whether to delete temporary files (default: FALSE).
#'
#' @return If PairwiseIBDOutput is NULL, returns a data.frame with columns ID1, ID2, pa, pb, pc.
#'   Otherwise, returns a character message indicating where the pairwise IBD results have been stored.
#'
#' @examples
#' PlinkPrefix <- file.path(system.file(package = "GRAB"), "extdata", "simuPLINK")
#' SparseGRMFile <- system.file("extdata", "SparseGRM.txt", package = "GRAB")
#' PairwiseIBDOutput <- file.path(tempdir(), "PairwiseIBD.txt")
#' getPairwiseIBD(PlinkPrefix, SparseGRMFile, PairwiseIBDOutput)
#'
#' @export
getPairwiseIBD <- function(
  PlinkPrefix,
  SparseGRMFile,
  PairwiseIBDOutput = NULL,
  frqFile = NULL,
  tempDir = NULL,
  maxSampleNums = 2500,
  minMafIBD = 0.01,
  rm.tempFile = FALSE
) {

  bedFile <- paste0(PlinkPrefix, ".bed")
  bimFile <- paste0(PlinkPrefix, ".bim")
  famFile <- paste0(PlinkPrefix, ".fam")

  if (is.null(frqFile)) {
    frqFile <- paste0(PlinkPrefix, ".frq")
  } else if (substr(frqFile, nchar(frqFile) - 3, nchar(frqFile)) != ".frq") {
    frqFile <- paste0(frqFile, ".frq")
  }
  if (!file.exists(frqFile)) {
    stop("Please check frqFile name ", frqFile, " is correct or this file exists!")
  }

  if (is.null(tempDir)) {
    tempDir <- tempdir()
  }

  # read all genotype and pass to QC.
  GenoInfoMat <- data.table::fread(frqFile)

  # read in the Sparse GRM.
  if (is.data.frame(SparseGRMFile)) {
    SparseGRMData <- SparseGRMFile
  } else {
    SparseGRMData <- data.table::fread(SparseGRMFile)
  }

  SparseGRMData$ID1 <- as.character(SparseGRMData$ID1)
  SparseGRMData$ID2 <- as.character(SparseGRMData$ID2)

  if (any(colnames(SparseGRMData) != c("ID1", "ID2", "Value"))) {
    stop("The column names of SparseGRMFile should be ['ID1', 'ID2', 'Value']!")
  }

  SubjID_related <- SparseGRMData %>%
    filter(ID1 != ID2) %>%
    select(ID1, ID2) %>%
    unlist(use.names = FALSE) %>%
    unique()

  if (length(SubjID_related) == 0) {
    PairwiseIBD <- as.data.frame(matrix(NA, 0, 5))

    colnames(PairwiseIBD) <- c("ID1", "ID2", "pa", "pb", "pc")

    if (!is.null(PairwiseIBDOutput)) {
      data.table::fwrite(PairwiseIBD, PairwiseIBDOutput,
        row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
      )
    }
  } else {
    # read in the bim and fam data.
    bim <- data.table::fread(bimFile)
    fam <- data.table::fread(famFile)

    # check SNPs of frqFile and subjects of SparseGRMFile correspond to PlinkFile.
    if (any(GenoInfoMat$SNP != bim$V2)) {
      stop("Please check if SNPs in FrqFile match with PlinkFile")
    }

    if (!all(SubjID_related %in% fam$V2)) {
      stop("Please check if related subjects in SparseGRMFile but not in PlinkFile")
    }

    GenoInfoMat <- GenoInfoMat %>% filter(MAF > minMafIBD)

    .message("IBD analysis: %d subjects, %d markers", length(SubjID_related), nrow(GenoInfoMat))

    if (nrow(GenoInfoMat) < 1e4) {
      warning("Number of Markers is a bit small, we recommend nSNPs > 10,000.\n")
    }

    if (nrow(GenoInfoMat) > 1e5) {
      warning("Number of Markers is a bit large, we recommend maxSampleNums < 2,500.\n")
    }

    # metrics used in pairwise IBD calculation.
    altFreq <- GenoInfoMat$MAF
    pro_var <- 2 * (altFreq * (1 - altFreq))^2 # 2*pi^2*(1-pi)^2, where pi comes from Binom(2, pi).
    wi <- sqrt(pro_var / (1 - pro_var)) # weights of each SNP

    # write the passed SNPIDs into IDsToInclude.
    IDsToIncludeFile <- paste0(tempDir, "/IDsToInclude.txt")
    data.table::fwrite(data.table::data.table(GenoInfoMat$SNP), IDsToIncludeFile,
      row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t"
    )

    edges <- t(SparseGRMData[, c("ID1", "ID2")])
    graph_GRM <- igraph::make_graph(edges, directed = FALSE)
    graph_list_all <- graph_GRM %>% igraph::decompose()
    graph_length <- lapply(graph_list_all, length)

    graph_list <- graph_list_all[graph_length > 1]
    graph_length <- lapply(graph_list, length) %>% unlist()

    # initialize parameters
    PairwiseIBD <- c()
    tSampleNums <- 0
    tSampleIDs <- c()
    nParts <- 1

    # cycle for calculating pairwise IBD
    for (i in seq_along(graph_list)) {
      tSampleNums <- tSampleNums + graph_length[i]
      tSampleIDs <- c(tSampleIDs, igraph::V(graph_list[[i]])$name)

      if (tSampleNums >= maxSampleNums || i == length(graph_list)) {
        .message("Processing block %d with %d subjects", nParts, tSampleNums)

        GenoList <- GRAB.ReadGeno(bedFile,
          SampleIDs = tSampleIDs,
          control = list(
            IDsToIncludeFile = IDsToIncludeFile,
            imputeMethod = "mean"
          )
        )

        tempGRM <- SparseGRMData %>%
          filter(ID1 %in% tSampleIDs & ID2 %in% tSampleIDs) %>%
          filter(ID1 != ID2) %>%
          mutate(idxID1 = match(ID1, rownames(GenoList$GenoMat))) %>%
          mutate(idxID2 = match(ID2, rownames(GenoList$GenoMat)))

        for (j in seq_len(nrow(tempGRM))) {
          tempmetrics <- GenoList$GenoMat[tempGRM$idxID1[j], ] - GenoList$GenoMat[tempGRM$idxID2[j], ]

          pc <- 0.5 * weighted.mean(((abs(tempmetrics - 1) + abs(tempmetrics + 1) - 2) / pro_var), wi, na.rm = TRUE)

          # Apply constraints to pc value
          upper_bound <- (1 - tempGRM$Value[j])^2 - 1e-10
          lower_bound <- 1 - 2 * tempGRM$Value[j]

          pc <- ifelse(pc > (1 - tempGRM$Value[j])^2, upper_bound, ifelse(pc < lower_bound, lower_bound, pc))
          pb <- 2 - 2 * pc - 2 * tempGRM$Value[j]
          pa <- 2 * tempGRM$Value[j] + pc - 1

          if (pb < 0) {
            pa <- pa + 0.5 * pb
            pb <- 0
            pc <- 0
          }

          PairwiseIBD <- rbind(
            PairwiseIBD,
            c(ID1 = tempGRM$ID1[j], ID2 = tempGRM$ID2[j], pa = pa, pb = pb, pc = pc)
          )
        }
        .message("Completed block %d", nParts)

        tSampleNums <- 0
        tSampleIDs <- c()
        nParts <- nParts + 1
      }
    }
  }

  PairwiseIBD <- data.table::as.data.table(PairwiseIBD)
  if (rm.tempFile) file.remove(IDsToIncludeFile)

  if (is.null(PairwiseIBDOutput)) {
    data.table::set(PairwiseIBD, j = "ID1", value = as.character(PairwiseIBD$ID1))
    data.table::set(PairwiseIBD, j = "ID2", value = as.character(PairwiseIBD$ID2))

    data.table::set(PairwiseIBD, j = "pa", value = as.numeric(PairwiseIBD$pa))
    data.table::set(PairwiseIBD, j = "pb", value = as.numeric(PairwiseIBD$pb))
    data.table::set(PairwiseIBD, j = "pc", value = as.numeric(PairwiseIBD$pc))
    
    return(PairwiseIBD)
  } else {
    data.table::fwrite(PairwiseIBD, PairwiseIBDOutput,
      row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
    )
    .message("The PairwiseIBD has been stored in %s", PairwiseIBDOutput)
    return(invisible(NULL))
  }
}
