#' Simulate genotype data matrix for related and unrelated subjects
#'
#' Generates genotype data for association studies, supporting both unrelated subjects
#' and family-based designs with various pedigree structures.
#'
#' @param nSub Number of unrelated subjects. If 0, all subjects are related.
#' @param nFam Number of families. If 0, all subjects are unrelated.
#' @param FamMode Family structure: "4-members", "10-members", or "20-members".
#'   See \code{Details} for pedigree structures.
#' @param nSNP Number of genetic markers to simulate.
#' @param MaxMAF Maximum minor allele frequency for simulation (default: 0.5).
#' @param MinMAF Minimum minor allele frequency for simulation (default: 0.05).
#' @param MAF Optional vector of specific MAF values for each marker. If provided,
#'   \code{MaxMAF} and \code{MinMAF} are ignored.
#' @return List containing:
#' \describe{
#'   \item{GenoMat}{Numeric genotype matrix (subjects × markers) with values 0, 1, 2.}
#'   \item{markerInfo}{Data frame with marker IDs and MAF values.}
#' }
#' @details
#' Genotypes are simulated under Hardy-Weinberg equilibrium with MAF ~ Uniform(MinMAF, MaxMAF).
#'
#' **Family Structures:**
#' \itemize{
#'   \item \strong{4-members}: 1+2→3+4 (parents 1,2 → offspring 3,4)
#'   \item \strong{10-members}: 1+2→5+6, 3+5→7+8, 4+6→9+10
#'   \item \strong{20-members}: Complex multi-generational pedigree with 20 members
#' }
#'
#' Total subjects: \code{nSub + nFam × family_size}
#'
#' @seealso \code{\link{GRAB.makePlink}} for converting to PLINK format.
#'
#' @examples
#' nSub <- 100
#' nFam <- 10
#' FamMode <- "10-members"
#' nSNP <- 10000
#' OutList <- GRAB.SimuGMat(nSub, nFam, FamMode, nSNP)
#' GenoMat <- OutList$GenoMat
#' markerInfo <- OutList$markerInfo
#' GenoMat[1:10, 1:10]
#' head(markerInfo)
#'
#' ## The following is to calculate GRM
#' MAF <- apply(GenoMat, 2, mean) / 2
#' GenoMatSD <- t((t(GenoMat) - 2 * MAF) / sqrt(2 * MAF * (1 - MAF)))
#' GRM <- GenoMatSD %*% t(GenoMatSD) / ncol(GenoMat)
#' GRM1 <- GRM[1:10, 1:10]
#' GRM2 <- GRM[100 + 1:10, 100 + 1:10]
#' GRM1
#' GRM2
#'
GRAB.SimuGMat <- function(
  nSub,
  nFam,
  FamMode,
  nSNP,
  MaxMAF = 0.5,
  MinMAF = 0.05,
  MAF = NULL
) {
  # Validate input parameters and family structure
  inputList <- checkInput(nSub, nFam, FamMode)

  nSubInEachFam <- inputList$nSubInEachFam
  nHaploInEachFam <- inputList$nHaploInEachFam
  fam.mat <- inputList$fam.mat
  nSub <- inputList$nSub
  nFam <- inputList$nFam

  # Calculate total number of subjects and haplotypes needed
  n <- nSub + nFam * nSubInEachFam
  nHaplo <- nFam * nHaploInEachFam

  if (n == 0) {
    stop("Please give at least one of 'nSub' and 'nFam'.")
  }

  .message("Simulation design: %d unrelated subjects, %d families (%d per family)",
           nSub, nFam, nSubInEachFam)
  .message("Total subjects: %d", n)

  # Generate or validate minor allele frequencies
  if (is.null(MAF)) {
    MAF <- runif(nSNP, MinMAF, MaxMAF)
  } else {
    if (length(MAF) != nSNP) {
      stop("length(MAF) should equal to nSNP.")
    }
    .message("Using provided MAF values (ignoring MaxMAF and MinMAF)")
  }

  # Create marker information table with SNP IDs and MAF values
  SNP.info <- make.SNP.info(nSNP, MAF)

  GenoMat1 <- GenoMat2 <- NULL

  # Simulate genotypes for related subjects (family members)
  if (nHaplo != 0) {
    .message("Simulating haplotype data for %d related subjects", nHaplo)
    haplo.mat <- haplo.simu(nHaplo, SNP.info)

    .message("Converting haplotypes to genotypes")
    GenoMat1 <- from.haplo.to.geno(haplo.mat, fam.mat)
  }

  # Simulate genotypes for unrelated subjects
  if (nSub != 0) {
    .message("Simulating genotype data for %d unrelated subjects", nSub)
    GenoMat2 <- geno.simu(nSub, SNP.info)
  }

  # Combine genotype matrices from related and unrelated subjects
  GenoMat <- rbind(GenoMat1, GenoMat2)

  return(list(
    GenoMat = GenoMat,
    markerInfo = SNP.info
  ))
}


# Validate input parameters and set default family structure
checkInput <- function(
  nSub,
  nFam,
  FamMode
) {
  if (missing(FamMode) && missing(nFam)) {
    .message("Simulating unrelated subjects only (no family structure specified)")
    nFam <- 0
    FamMode <- "Unrelated"
  }

  if (missing(nSub)) {
    .message("Simulating family members only (no unrelated subjects)")
    nSub <- 0
  }

  # Validate family mode against supported options
  if (!is.element(FamMode, c("Unrelated", "4-members", "10-members", "20-members"))) {
    stop(
      "FamMode should be one of 'Unrelated', '4-members', '10-members', and '20-members'. ",
      "Other input is not supported."
    )
  }

  # Initialize family structure parameters
  nSubInEachFam <- 0
  nHaploInEachFam <- 0
  fam.mat <- NULL

  # Set family-specific parameters based on pedigree structure
  if (FamMode == "4-members") {
    nSubInEachFam <- 4
    nHaploInEachFam <- 4
    fam.mat <- example.fam.4.members(nFam)
  }

  if (FamMode == "10-members") {
    nSubInEachFam <- 10
    nHaploInEachFam <- 8
    fam.mat <- example.fam.10.members(nFam)
  }

  if (FamMode == "20-members") {
    nSubInEachFam <- 20
    nHaploInEachFam <- 16
    fam.mat <- example.fam.20.members(nFam)
  }

  # Return structured list of family parameters
  inputList <- list(
    nSubInEachFam = nSubInEachFam,
    nHaploInEachFam = nHaploInEachFam,
    fam.mat = fam.mat,
    nSub = nSub,
    nFam = nFam,
    FamMode = FamMode
  )
  return(inputList)
}


# Generate 20-member family structure with complex pedigree
# Pedigree: 1+2→9+10; 3+9→11+12; 4+10→13+14; 5+11→15+16; 6+12→17; 7+13→18; 8+14→19+20
example.fam.20.members <- function(n.fam) { # Number of families to generate
  m <- 20 # Number of family members per family
  fam.mat <- c()
  c.h <- 0 # Haplotype counter across all families

  for (i in 1:n.fam) {
    FID <- paste0("f", i)
    IID <- paste0(FID, "_", 1:m)

    # Define roles: first 8 are founders (independent), rest are offspring
    Role <- c(
      rep("Founder", 8),
      rep("Offspring", 12)
    )

    # Source 1 (first parent/haplotype source)
    Source1 <- c(
      paste0("haplo-", c.h + 1:8),
      IID[c(1, 1, 3, 3, 4, 4, 5, 5, 6, 7, 8, 8)]
    )

    # Source 2 (second parent/haplotype source)
    Source2 <- c(
      paste0("haplo-", c.h + 9:16),
      IID[c(2, 2, 9, 9, 10, 10, 11, 11, 12, 13, 14, 14)]
    )

    c.h <- c.h + 16  # Update haplotype counter for next family

    # Combine family information into matrix
    fam.mat <- rbind(
      fam.mat,
      cbind(FID, IID, Role, Source1, Source2)
    )
  }

  fam.mat <- data.frame(fam.mat, stringsAsFactors = FALSE)
  return(fam.mat) # Data frame with columns: FID, IID, Role, Source1, Source2
}


# Generate 10-member family structure with two-generation pedigree
# Pedigree: 1+2→5+6; 3+5→7+8; 4+6→9+10
example.fam.10.members <- function(n.fam) { # Number of families to generate
  m <- 10 # Number of family members per family
  fam.mat <- c()
  c.h <- 0 # Haplotype counter across all families

  for (i in 1:n.fam) {
    FID <- paste0("f", i)
    IID <- paste0(FID, "_", 1:m)

    # Define roles: first 4 are founders, rest are offspring
    Role <- c(
      rep("Founder", 4),
      rep("Offspring", 6)
    )

    # Source 1 (first parent/haplotype source)
    Source1 <- c(
      paste0("haplo-", c.h + 1:4),
      IID[c(1, 1, 3, 3, 4, 4)]
    )

    # Source 2 (second parent/haplotype source)
    Source2 <- c(
      paste0("haplo-", c.h + 5:8),
      IID[c(2, 2, 5, 5, 6, 6)]
    )

    c.h <- c.h + 8  # Update haplotype counter for next family

    # Combine family information into matrix
    fam.mat <- rbind(
      fam.mat,
      cbind(FID, IID, Role, Source1, Source2)
    )
  }

  fam.mat <- data.frame(fam.mat, stringsAsFactors = FALSE)
  return(fam.mat) # Data frame with columns: FID, IID, Role, Source1, Source2
}


# Generate 4-member family structure with simple nuclear family
# Pedigree: 1+2→3+4 (parents 1,2 produce offspring 3,4)
example.fam.4.members <- function(n.fam) { # Number of families to generate
  m <- 4 # Number of family members per family
  fam.mat <- c()
  c.h <- 0 # Haplotype counter across all families

  for (i in 1:n.fam) {
    FID <- paste0("f", i)
    IID <- paste0(FID, "_", 1:m)

    # Define roles: first 2 are founders (parents), last 2 are offspring
    Role <- c(
      rep("Founder", 2),
      rep("Offspring", 2)
    )

    # Source 1 (first parent/haplotype source)
    Source1 <- c(
      paste0("haplo-", c.h + 1:2),
      IID[c(1, 1)]
    )

    # Source 2 (second parent/haplotype source)
    Source2 <- c(
      paste0("haplo-", c.h + 3:4),
      IID[c(2, 2)]
    )

    c.h <- c.h + 4  # Update haplotype counter for next family

    # Combine family information into matrix
    fam.mat <- rbind(
      fam.mat,
      cbind(FID, IID, Role, Source1, Source2)
    )
  }

  fam.mat <- data.frame(fam.mat, stringsAsFactors = FALSE)
  return(fam.mat) # Data frame with columns: FID, IID, Role, Source1, Source2
}


# Create marker information table with SNP IDs and MAF values
make.SNP.info <- function(
  nSNP, # Number of SNPs to create
  MAF   # Vector of minor allele frequencies
) {
  SNP.info <- data.table::data.table(
    SNP = paste0("SNP_", 1:nSNP),
    MAF = MAF,
    stringsAsFactors = FALSE
  )
  return(SNP.info)
}


# Simulate genotype matrix for unrelated subjects under Hardy-Weinberg equilibrium
geno.simu <- function(
  nSub,
  SNP.info
) {
  nSNPs <- nrow(SNP.info)
  MAFs <- SNP.info$MAF

  # Generate genotypes (0, 1, 2) using binomial distribution
  GenoMat <- sapply(MAFs, FUN = function(x) {
    rbinom(nSub, 2, x)
  })

  # Handle case where nSub = 1 (sapply returns vector instead of matrix)
  if (nSub == 1) {
    GenoMat <- matrix(GenoMat, 1, nSNPs)
  }

  # Add row and column names for identification
  colnames(GenoMat) <- SNP.info$SNP
  rownames(GenoMat) <- paste0("Subj-", 1:nSub)

  return(GenoMat) # Matrix: subjects × SNPs
}


# Simulate haplotype matrix for founder individuals
haplo.simu <- function(
  n.haplo, # Number of haplotypes to simulate
  SNP.info # Marker information with MAF values
) {
  # Extract marker information
  MAFs <- SNP.info$MAF

  # Generate haplotype matrix: each row = one haplotype, each column = one SNP
  # Alleles drawn independently from Bernoulli distribution with probability = MAF
  haplo.mat <- sapply(MAFs,
    FUN = function(x) {
      rbinom(n.haplo, 1, x)
    }
  )

  # Add identifiers for haplotypes and SNPs
  colnames(haplo.mat) <- SNP.info$SNP
  rownames(haplo.mat) <- paste0("haplo-", 1:n.haplo)

  return(haplo.mat) # Matrix: haplotypes × SNPs
}


# Convert genotype matrix to haplotype matrix (for phasing)
from.geno.to.haplo <- function(geno.mat) { # Matrix: subjects × markers
  # Initialize haplotype matrices by splitting genotypes
  haplo.mat1 <- haplo.mat2 <- geno.mat / 2 # 2→1; 1→0.5; 0→0; -9→-4.5

  # Handle missing values
  posMissing <- which(geno.mat == -9)
  haplo.mat1[posMissing] <- haplo.mat2[posMissing] <- -9

  # Randomly phase heterozygous genotypes (genotype = 1)
  posHetero <- which(geno.mat == 1)
  haplo.mat1[posHetero] <- rbinom(length(posHetero), 1, 0.5)
  haplo.mat2[posHetero] <- 1 - haplo.mat1[posHetero]

  # Combine both haplotype matrices (each subject contributes 2 haplotypes)
  haplo.mat <- rbind(haplo.mat1, haplo.mat2)

  return(haplo.mat) # Matrix: (2 × subjects) × markers
}


# Convert haplotype data to genotype matrix based on family structure
from.haplo.to.geno <- function(
  haplo.mat, # Matrix: haplotypes × SNPs
  fam.mat    # Family structure: subjects × family_info
) {
  n <- nrow(fam.mat) # Number of subjects to generate
  m <- ncol(haplo.mat) # Number of SNPs

  # Initialize haplotype matrices for each subject
  Haplo1.mat <- matrix(nrow = n, ncol = m)
  Haplo2.mat <- matrix(nrow = n, ncol = m)
  rownames(Haplo1.mat) <- rownames(Haplo2.mat) <- fam.mat$IID

  # Generate genotypes for each subject based on family relationships
  for (i in 1:n) {
    if (i %% 1000 == 0) .message("Genotype simulation progress: %d subjects", i)

    Role <- fam.mat$Role[i]
    S1 <- fam.mat$Source1[i]  # First parent/haplotype source
    S2 <- fam.mat$Source2[i]  # Second parent/haplotype source

    if (Role == "Founder") {
      # Founders: directly use haplotype data
      Haplo1.mat[i, ] <- haplo.mat[S1, ]
      Haplo2.mat[i, ] <- haplo.mat[S2, ]
    }

    if (Role == "Offspring") {
      # Offspring: inherit haplotypes from parents through Mendelian transmission

      # Randomly select one haplotype from parent S1
      S1.pass.H1 <- rbinom(m, 1, 0.5)  # Random transmission: 50% chance each haplotype
      S1.pass.H2 <- 1 - S1.pass.H1
      Haplo1.mat[i, ] <- Haplo1.mat[S1, ] * S1.pass.H1 + Haplo2.mat[S1, ] * S1.pass.H2

      # Randomly select one haplotype from parent S2
      S2.pass.H1 <- rbinom(m, 1, 0.5)
      S2.pass.H2 <- 1 - S2.pass.H1
      Haplo2.mat[i, ] <- Haplo1.mat[S2, ] * S2.pass.H1 + Haplo2.mat[S2, ] * S2.pass.H2
    }
  }

  # Combine two haplotypes to form genotypes (0, 1, or 2 copies of minor allele)
  Geno.mat <- Haplo1.mat + Haplo2.mat
  colnames(Geno.mat) <- colnames(haplo.mat)

  return(Geno.mat)
}


#' Simulate random effects based on family structure
#'
#' Generates random effect vectors (bVec) that account for family relationships
#' through kinship-based correlation structures.
#'
#' @param nSub Number of unrelated subjects. If 0, all subjects are related.
#' @param nFam Number of families. If 0, all subjects are unrelated.
#' @param FamMode Family structure: "4-members", "10-members", or "20-members".
#'   See \code{\link{GRAB.SimuGMat}} for pedigree details.
#' @param tau Variance component for random effects.
#' @return Data frame with columns:
#' \describe{
#'   \item{IID}{Subject identifiers.}
#'   \item{bVec}{Random effect values following appropriate correlation structure.}
#' }
#' @details
#' For related subjects, random effects follow a multivariate normal distribution
#' with covariance proportional to kinship coefficients. For unrelated subjects,
#' effects are independent normal random variables.
#'
GRAB.SimubVec <- function(
  nSub,
  nFam,
  FamMode,
  tau
) {
  # Validate input parameters and family structure
  inputList <- checkInput(nSub, nFam, FamMode)

  nSubInEachFam <- inputList$nSubInEachFam
  nSub <- inputList$nSub
  nFam <- inputList$nFam
  FamMode <- inputList$FamMode
  fam.mat <- inputList$fam.mat

  n <- nSub + nFam * nSubInEachFam

  if (n == 0) {
    stop("Please give at least one of 'nSub' and 'nFam'.")
  }

  .message("bVec simulation: %d unrelated subjects, %d families (%d per family)",
           nSub, nFam, nSubInEachFam)
  .message("Total subjects: %d", n)

  # Generate random effects for related subjects using kinship-based covariance
  if (FamMode == "Unrelated") {
    bVec.Related <- data.table::data.table()
  } else {
    # Load appropriate kinship matrix for the family structure
    if (FamMode == "4-members") {
      fam.kin.file <- system.file("extdata", "example_4-members.kin.txt", package = "GRAB")
    }

    if (FamMode == "10-members") {
      fam.kin.file <- system.file("extdata", "example_10-members.kin.txt", package = "GRAB")
    }

    if (FamMode == "20-members") {
      fam.kin.file <- system.file("extdata", "example_20-members.kin.txt", package = "GRAB")
    }

    # Read kinship matrix and generate correlated random effects
    fam.kin <- data.table::fread(fam.kin.file)
    fam.kin <- as.matrix(fam.kin)

    n <- nFam * nrow(fam.kin)
    out.eigen <- eigen(fam.kin)
    factor <- t(out.eigen$vectors) * sqrt(out.eigen$values)
    kin.chol <- diag(nFam) %x% factor  # Kronecker product for multiple families
    b.true <- t(kin.chol) %*% rnorm(n) * sqrt(tau)
    bVec.Related <- data.table::data.table(
      IID = fam.mat$IID,
      bVec = as.numeric(b.true)
    )
  }

  # Generate independent random effects for unrelated subjects
  if (nSub != 0) {
    bVec.Unrelated <- data.table::data.table(
      IID = paste0("Subj-", 1:nSub),
      bVec = rnorm(nSub, sd = tau)
    )
  }

  # Combine random effects from related and unrelated subjects
  bVec <- rbind(bVec.Related, bVec.Unrelated)

  return(bVec)
}


#' Simulate genotype matrix from external genotype file
#'
#' Generates genotype matrices for families and unrelated subjects using
#' haplotype data from existing genotype files. Primarily designed for
#' rare variant analysis simulations.
#'
#' @param nFam Number of families to simulate.
#' @param nSub Number of unrelated subjects to simulate.
#' @param FamMode Family structure: "4-members", "10-members", or "20-members".
#'   See Details for pedigree structures.
#' @param GenoFile Path to genotype file (passed to \code{\link{GRAB.ReadGeno}}).
#' @param GenoFileIndex Index file for genotype data (optional, for BGEN files).
#' @param SampleIDs Vector of sample IDs to include (optional).
#' @param control List of control parameters passed to \code{\link{GRAB.ReadGeno}}.
#' @return List containing:
#' \describe{
#'   \item{GenoMat}{Simulated genotype matrix (subjects × variants).}
#'   \item{SubjIDs}{Subject identifiers for simulated samples.}
#'   \item{markerInfo}{Variant information from input file.}
#' }
#' @details
#' This function supports both unrelated and related subjects. Founder genotypes
#' are sampled from the input genotype file, and offspring genotypes are generated
#' through Mendelian inheritance.
#'
#' **Note**: When simulating related subjects, alleles are randomly assigned to
#' haplotypes during the phasing process.
#'
#' ## Family Structures:
#' - **4-members**: Total subjects = nSub + 4×nFam. Structure: 1+2→3+4
#' - **10-members**: Total subjects = nSub + 10×nFam.
#'   Structure: 1+2→5+6, 3+5→7+8, 4+6→9+10
#' - **20-members**: Total subjects = nSub + 20×nFam.
#'   Structure: 1+2→9+10, 3+9→11+12, 4+10→13+14, 5+11→15+16,
#'   6+12→17, 7+13→18, 8+14→19+20
#'
#' @examples
#' nFam <- 50
#' nSub <- 500
#' FamMode <- "10-members"
#'
#' # PLINK data format. Currently, this function does not support BGEN data format.
#' PLINKFile <- system.file("extdata", "example_n1000_m236.bed", package = "GRAB")
#' IDsToIncludeFile <- system.file("extdata", "example_n1000_m236.IDsToInclude", package = "GRAB")
#'
#' GenoList <- GRAB.SimuGMatFromGenoFile(nFam, nSub, FamMode, PLINKFile,
#'   control = list(IDsToIncludeFile = IDsToIncludeFile)
#' )
#'
#'
GRAB.SimuGMatFromGenoFile <- function(
  nFam,
  nSub,
  FamMode,
  GenoFile,
  GenoFileIndex = NULL,
  SampleIDs = NULL,
  control = NULL
) {
  # Validate input parameters and prepare family structure
  inputList <- checkInput(nSub, nFam, FamMode)

  nSubInEachFam <- inputList$nSubInEachFam
  nHaploInEachFam <- inputList$nHaploInEachFam
  fam.mat <- inputList$fam.mat
  nSub <- inputList$nSub
  nFam <- inputList$nFam

  n <- nSub + nFam * nSubInEachFam  # Total number of subjects
  nHaplo <- nFam * nHaploInEachFam  # Total number of haplotypes needed

  if (n == 0) {
    stop("Please give at least one of 'nSub' and 'nFam'.")
  }

  .message("Genotype matrix simulation: %d unrelated subjects, %d families (%d per family)",
           nSub, nFam, nSubInEachFam)
  .message("Total subjects: %d", n)

  # Check file format compatibility for family simulations
  GenoFileExt <- tools::file_ext(GenoFile)
  if (GenoFileExt != "bed" && nFam != 0) {
    stop("Current version of 'GRAB.SimuGMatFromGenoFile()' only supports ",
         "PLINK files when simulating related subjects.")
  }

  # Read genotype data from input file
  GenoList <- GRAB.ReadGeno(GenoFile, GenoFileIndex, SampleIDs, control)
  GenoMat <- GenoList$GenoMat
  markerInfo <- GenoList$markerInfo
  nSubInGeno <- nrow(GenoMat)

  # Ensure sufficient samples are available in the input file
  if (nSubInGeno < nSub + nHaplo / 2) {
    stop("Number of subjects in Genotype File (", nSubInGeno,
         ") < Number of subjects requested (", nSub + nHaplo / 2, ").")
  }

  # Initialize genotype matrices for related and unrelated subjects
  GenoMat1 <- GenoMat2 <- NULL

  # Randomly sample subjects from input file for simulation
  randomRow <- sample(nSubInGeno)
  rowForHaplo <- randomRow[1:(nHaplo / 2)]  # Samples for haplotype extraction
  rowForSub <- randomRow[1:nSub + nHaplo / 2]  # Samples for unrelated subjects

  # Process related subjects through haplotype-based inheritance
  if (nHaplo != 0) {
    .message("Extracting haplotype data for %d related subjects", nHaplo)
    GenoMatTemp <- GenoMat[rowForHaplo, ]

    # Convert genotypes to haplotypes for Mendelian transmission
    haplo.mat <- from.geno.to.haplo(GenoMatTemp)
    rownames(haplo.mat) <- paste0("haplo-", 1:nHaplo)

    .message("Converting haplotypes to genotypes for related subjects")
    # Generate family genotypes through inheritance patterns
    GenoMat1 <- from.haplo.to.geno(haplo.mat, fam.mat)
  }

  # Process unrelated subjects by direct sampling
  if (nSub != 0) {
    .message("Extracting genotype data for %d unrelated subjects", nSub)
    GenoMat2 <- GenoMat[rowForSub, ]
    rownames(GenoMat2) <- paste0("Subj-", 1:nSub)
  }

  # Combine genotype data from related and unrelated subjects
  GenoMat <- rbind(GenoMat1, GenoMat2)

  return(list(
    GenoMat = GenoMat,
    markerInfo = markerInfo
  ))
}


#' Convert genotype matrix to PLINK format files
#'
#' Converts a numeric genotype matrix to PLINK text files (PED and MAP format)
#' for use with PLINK software and other genetic analysis tools.
#'
#' @param GenoMat Numeric genotype matrix (n×m) with values 0, 1, 2, or -9.
#'   Rows = subjects, columns = markers. Row and column names are required.
#' @param OutputPrefix Output file prefix including path (without extension).
#' @param A1 Allele 1 character, usually minor/ALT allele (default: "G").
#' @param A2 Allele 2 character, usually major/REF allele (default: "A").
#' @param CHR Chromosome numbers for markers (default: all chromosome 1).
#' @param BP Base positions for markers (default: 1:m).
#' @param Pheno Phenotype values for subjects (default: all missing as -9).
#' @param Sex Sex codes for subjects (default: all coded as 1).
#' @return Character message confirming file creation location.
#' @details
#' **Genotype Encoding:**
#' - 0, 1, 2 → copies of minor allele
#' - -9 → missing genotype (coded as "00" in PED)
#' - A1="G", A2="A": 0→"GG", 1→"AG", 2→"AA", -9→"00"
#'
#' **Output Files:**
#' - `.ped`: Pedigree file with genotype data
#' - `.map`: Marker map file with positions
#'
#' **Downstream Processing:**
#' ```bash
#' # Convert to binary format
#' plink --file prefix --make-bed --out prefix
#'
#' # Convert to raw format
#' plink --bfile prefix --recode A --out prefix_raw
#'
#' # Convert to BGEN format
#' plink2 --bfile prefix --export bgen-1.2 bits=8 ref-first --out prefix_bgen
#'
#' # Create BGEN index
#' bgenix -g prefix_bgen.bgen --index
#' ```
#'
#' @examples
#' ### Step 1: simulate a numeric genotype matrix
#' n <- 1000
#' m <- 20
#' MAF <- 0.3
#' set.seed(123)
#' GenoMat <- matrix(rbinom(n * m, 2, MAF), n, m)
#' rownames(GenoMat) <- paste0("Subj-", 1:n)
#' colnames(GenoMat) <- paste0("SNP-", 1:m)
#' OutputDir <- tempdir()
#' outputPrefix <- file.path(OutputDir, "simuPLINK")
#'
#' ### Step 2(a): make PLINK files without missing genotype
#' GRAB.makePlink(GenoMat, outputPrefix)
#'
#' ### Step 2(b): make PLINK files with genotype missing rate of 0.1
#' indexMissing <- sample(n * m, 0.1 * n * m)
#' GenoMat[indexMissing] <- -9
#' GRAB.makePlink(GenoMat, outputPrefix)
#'
#' ## The following are in shell environment
#' # plink --file simuPLINK --make-bed --out simuPLINK
#' # plink --bfile simuPLINK --recode A --out simuRAW
#' # plink2 --bfile simuPLINK --export bgen-1.2 bits=8 ref-first --out simuBGEN
#' # UK Biobank use 'ref-first'
#' # bgenix -g simuBGEN.bgen --index
#'
#'
GRAB.makePlink <- function(
  GenoMat,
  OutputPrefix,
  A1 = "G",
  A2 = "A",
  CHR = NULL,
  BP = NULL,
  Pheno = NULL,
  Sex = NULL
) {
  # Validate input genotype matrix format
  if (!is.numeric(GenoMat)) {
    stop("'GenoMat' should be a numeric matrix.")
  }

  if (any(!unique(as.numeric(GenoMat)) %in% c(0, 1, 2, -9))) {
    stop("'GenoMat' should only include elements of 0, 1, 2, -9.")
  }

  # Validate allele specification
  if (length(A1) != 1 || length(A2) != 1) {
    stop("Argument A1 and A2 should be a character, not a character vector.")
  }

  # Extract marker and subject identifiers
  SNP <- colnames(GenoMat)
  FID <- IID <- rownames(GenoMat)

  if (is.null(SNP) || is.null(IID)) {
    stop("rownames and colnames of GenoMat should be specified.")
  }

  m <- length(SNP)  # Number of markers
  n <- length(IID)  # Number of subjects

  .message("Creating PLINK files: %d markers, %d samples", m, n)

  # Set default phenotype values (missing) if not provided
  if (is.null(Pheno)) {
    Pheno <- rep(-9, n)
  } else {
    if (length(Pheno) != n) {
      stop("length(Pheno) should be the same as nrow(GenoMat).")
    }
  }

  # Set default sex codes if not provided
  if (is.null(Sex)) {
    Sex <- rep(1, n)
  } else {
    if (length(Sex) != n) {
      stop("length(Sex) should be the same as nrow(GenoMat).")
    }
  }

  # Set default base positions if not provided
  if (is.null(BP)) {
    BP <- 1:m
  } else {
    if (length(BP) != m) {
      stop("length(BP) should be the same as ncol(GenoMat).")
    }
  }

  # Set default chromosome numbers if not provided
  if (is.null(CHR)) {
    CHR <- rep(1, m)
  } else {
    if (length(CHR) != m) {
      stop("length(CHR) should be the same as ncol(GenoMat).")
    }
  }

  # Create PED file header with pedigree information
  PED <- cbind(
    FID = FID,      # Family ID
    IID = IID,      # Individual ID
    PID = 0,        # Paternal ID (0 = missing)
    MID = 0,        # Maternal ID (0 = missing)
    Sex = Sex,      # Sex code
    Phen = Pheno    # Phenotype
  )

  # Create MAP file with marker information
  MAP <- cbind(
    CHR = CHR,        # Chromosome
    SNP = SNP,        # SNP identifier
    GeneDist = 0,     # Genetic distance (cM, set to 0)
    BP = BP           # Base position
  )

  # Convert numeric genotypes to allele pairs for PED format
  # Initialize both alleles with A1 (or "0" for missing)
  Geno.ped1 <- Geno.ped2 <- ifelse(GenoMat == -9, "0", A1)

  # Set first allele: 1 or 2 copies → A2
  Geno.ped1 <- ifelse(GenoMat >= 1, A2, Geno.ped1)

  # Set second allele: 2 copies → A2
  Geno.ped2 <- ifelse(GenoMat >= 2, A2, Geno.ped2)

  # Interleave alleles for PED format (allele1, allele2, allele1, allele2, ...)
  Geno.ped <- matrix(nrow = n, ncol = 2 * m)
  Geno.ped[, seq(1, 2 * m, 2)] <- Geno.ped1  # Odd columns: first alleles
  Geno.ped[, seq(2, 2 * m, 2)] <- Geno.ped2  # Even columns: second alleles

  # Combine pedigree info with genotype data
  PED <- cbind(PED, Geno.ped)

  # Define output file paths
  MAP.file <- paste0(OutputPrefix, ".map")
  PED.file <- paste0(OutputPrefix, ".ped")

  # Write PLINK files
  data.table::fwrite(MAP, MAP.file, quote = FALSE, col.names = FALSE,
                     row.names = FALSE, sep = " ")
  data.table::fwrite(PED, PED.file, quote = FALSE, col.names = FALSE,
                     row.names = FALSE, sep = " ")

  .message("PLINK files saved to: %s", OutputPrefix)
  .message("Working directory: %s", getwd())

  message <- paste0("PLINK files have been saved to ", OutputPrefix, ".")
  return(message)
}


#' Simulate phenotypes from linear predictors
#'
#' Generates various types of phenotypes (quantitative, binary, ordinal,
#' time-to-event) from linear predictors using appropriate link functions.
#'
#' @param eta Vector of linear predictors (typically covariates×beta + genotypes×beta).
#' @param traitType Type of phenotype: "quantitative", "binary", "ordinal", or "time-to-event".
#' @param control List of simulation parameters specific to each trait type:
#' \describe{
#'   \item{pCase}{Proportion of cases (binary traits).}
#'   \item{sdError}{Error term standard deviation (quantitative traits).}
#'   \item{pEachGroup}{Group proportions (ordinal traits).}
#'   \item{eventRate}{Event rate (time-to-event traits).}
#' }
#' @param seed Random seed for reproducibility (optional).
#' @return Simulated phenotype vector or data frame:
#' \describe{
#'   \item{quantitative}{Numeric vector of continuous values.}
#'   \item{binary}{Numeric vector of 0/1 values.}
#'   \item{ordinal}{Numeric vector of categorical levels.}
#'   \item{time-to-event}{Data frame with SurvTime and SurvEvent columns.}
#' }
#' @details
#' **Trait Type Details:**
#' - **Quantitative**: Y = eta + error, where error ~ N(0, sdError²)
#' - **Binary**: Logistic model with specified case proportion
#' - **Ordinal**: Proportional odds model with specified group proportions
#' - **Time-to-event**: Weibull survival model with specified event rate
#'
#' For more details, see: https://wenjianbi.github.io/grab.github.io/docs/simulation_phenotype.html
#'
#'
GRAB.SimuPheno <- function(
  eta,
  traitType = "binary",
  control = list(
    pCase = 0.1,
    sdError = 1,
    pEachGroup = c(1, 1, 1),
    eventRate = 0.1
  ),
  seed = NULL
) {
  # Validate trait type
  if (!traitType %in% c("quantitative", "binary", "ordinal", "time-to-event")) {
    stop('"traitType" is limited to "quantitative", "binary", "ordinal", and "time-to-event".')
  }

  # Validate required control parameters for each trait type
  if (traitType == "binary") {
    if (!"pCase" %in% names(control)) {
      stop("For binary phenotype, argument 'control' should include 'pCase' ",
           "which is the proportion of cases.")
    }
  }

  if (traitType == "quantitative") {
    if (!"sdError" %in% names(control)) {
      .message("Note: For quantitative phenotype, 'control$sdError' should ",
               "specify error term standard deviation")
    }
  }

  if (traitType == "ordinal") {
    if (!"pEachGroup" %in% names(control)) {
      .message("Note: For ordinal categorical phenotype, 'control$pEachGroup' ",
               "should specify group proportions")
    }
  }

  if (traitType == "time-to-event") {
    if (!"eventRate" %in% names(control)) {
      .message("Note: For time-to-event phenotype, 'control$eventRate' ",
               "should specify event rate")
    }
  }

  # Center linear predictors and set random seed
  eta <- eta - mean(eta)
  n <- length(eta)

  if (!is.null(seed)) set.seed(seed)

  # Generate quantitative phenotypes: Y = eta + error
  if (traitType == "quantitative") {
    sdError <- control$sdError
    error <- rnorm(n, sd = sdError)
    pheno <- eta + error
    return(pheno)
  }

  # Generate binary phenotypes using logistic model
  if (traitType == "binary") {
    pCase <- control$pCase
    # Find intercept that achieves desired case proportion
    eta0 <- uniroot(f.binary, c(-100, 100), eta = eta, pCase = pCase, seed = seed)
    eta0 <- eta0$root

    # Apply logistic transformation
    eta.new <- eta0 + eta
    mu <- exp(eta.new) / (1 + exp(eta.new))  # Case probability
    pheno <- rbinom(n, 1, mu)  # Generate case-control status
    return(pheno)
  }

  # Generate ordinal phenotypes using proportional odds model
  if (traitType == "ordinal") {
    pEachGroup <- control$pEachGroup
    # Calculate threshold values for each ordinal category
    Eps <- getEps(pEachGroup, eta, seed)

    # Generate latent continuous variable
    pheno.latent <- runif(n)
    pheno <- rep(0, n)

    # Assign ordinal categories based on thresholds
    for (g in seq_along(Eps)) {
      mu <- exp(Eps[g] - eta) / (1 + exp(Eps[g] - eta))
      pheno[pheno.latent > mu] <- g
    }
    return(pheno)
  }

  # Generate time-to-event phenotypes using Weibull survival model
  if (traitType == "time-to-event") {
    eventRate <- control$eventRate

    # Fixed survival model parameters
    shape0 <- 2      # Weibull shape parameter
    cens.shape <- 1  # Censoring distribution shape
    cens.scale <- 0.15  # Censoring distribution scale

    # Find scale parameter that achieves desired event rate
    scale0 <- uniroot(f.surv, c(-100, 100),
                      eta.true = eta, event.rate = eventRate, seed = seed,
                      shape0 = shape0, cens.shape = cens.shape,
                      cens.scale = cens.scale)
    scale0 <- scale0$root

    # Generate censoring times and survival times
    cens <- rweibull(length(eta), shape = cens.shape, scale = cens.scale)
    eps <- runif(length(eta), 0, 1)
    # Generate survival times with Weibull distribution
    time <- (-log(eps) * exp(-1 * eta))^(1 / shape0) * scale0
    surv.time <- pmin(time, cens)  # Apply censoring
    event <- ifelse(time < cens, 1, 0)  # Event indicator (1=event, 0=censored)

    # Return survival data frame
    pheno <- data.frame(SurvTime = surv.time, SurvEvent = event)
    return(pheno)
  }
}


# Helper function: Find intercept for binary traits to achieve target prevalence
f.binary <- function(
  eta,     # Linear predictors
  pCase,   # Target case proportion
  eta0,    # Intercept to optimize
  seed     # Random seed
) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(eta)
  eta.new <- eta0 + eta
  mu <- exp(eta.new) / (1 + exp(eta.new))  # Case probabilities
  Y <- rbinom(n, 1, mu)  # Simulate case-control status
  re <- mean(Y) - pCase  # Difference from target proportion
  return(re)
}


# Helper function: Calculate probability difference for ordinal traits
getProb <- function(
  eps,       # Threshold parameter
  eta.true,  # Linear predictors
  prob,      # Target probability
  seed       # Random seed
) {
  if (!is.null(seed)) set.seed(seed)
  n <- length(eta.true)
  mu <- exp(eps - eta.true) / (1 + exp(eps - eta.true))
  Y.latent <- runif(n)
  diffprob <- mean(Y.latent < mu) - prob
  return(diffprob)
}


# Helper function: Calculate threshold parameters for ordinal traits
getEps <- function(
  ratios,    # Group proportion ratios
  eta.true,  # Linear predictors
  seed       # Random seed
) {
  sumR <- sum(ratios)
  cumR <- 0
  J <- length(ratios)
  Eps <- c()

  # Find threshold for each ordinal category (except the last)
  for (i in 1:(J - 1)) {
    cumR <- cumR + ratios[i]
    eps <- uniroot(getProb, c(-100, 100),
                   eta.true = eta.true, prob = cumR / sumR, seed = seed)
    Eps <- c(Eps, eps$root)
  }
  return(Eps)
}


# Helper function: Find scale parameter for survival traits to achieve target event rate
f.surv <- function(
  scale0,      # Scale parameter to optimize
  eta.true,    # Linear predictors
  event.rate,  # Target event rate
  seed,        # Random seed
  shape0 = 2,  # Weibull shape parameter
  cens.shape = 1,    # Censoring shape
  cens.scale = 0.15  # Censoring scale
) {
  if (!is.null(seed)) set.seed(seed)

  # Generate censoring times and survival times
  cens <- rweibull(length(eta.true), shape = cens.shape, scale = cens.scale)
  eps <- runif(length(eta.true), 0, 1)
  time <- (-log(eps) * exp(-1 * eta.true))^(1 / shape0) * scale0
  event <- ifelse(time < cens, 1, 0)

  # Return difference from target event rate
  re <- mean(event) - event.rate
  return(re)
}
