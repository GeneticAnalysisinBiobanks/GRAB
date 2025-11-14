## ------------------------------------------------------------------------------
## Util.R
## Common utilities shared across GRAB analyses: logging, control merging,
## GRM helpers, output checkpointing, and p-value aggregation.
##
## Functions:
##   .message         : Structured log messages with timestamps
##   .printParameters : Print formatted analysis parameters with timestamp
##   updateControl    : Merge user control options into defaults
##   updateSparseGRM  : Filter/reshape sparse GRM for a subject subset
##   checkOutputFile  : Validate/resume output with an index file
##   writeOutputFile  : Append outputs and progress safely
##   CCT              : Cauchy Combination Test (combine p-values)
## ------------------------------------------------------------------------------

# Helper function for structured logging
.message <- function(msg, ...) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- paste0("[INFO] ", timestamp, " ", sprintf(msg, ...))
  message(formatted_msg)
}


# Helper function to print formatted analysis parameters
.printParameters <- function(title, params, control) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  
  # Print title with timestamp
  message(sprintf("[INFO] %s %s ...", timestamp, title))
  
  # Print each parameter with indentation
  for (name in names(params)) {
    value <- params[[name]]
    if (!is.null(value)) {
      message(sprintf("    %s: %s", name, paste(value, collapse = ", ")))
    }
  }
  
  # Print control parameters with indentation
  if (!is.null(control) && length(control) > 0) {
    message("    Control parameters:")
    for (name in names(control)) {
      value <- control[[name]]
      # Format value based on type
      if (is.numeric(value)) {
        message(sprintf("      %s: %g", name, value))
      } else if (is.logical(value)) {
        message(sprintf("      %s: %s", name, as.character(value)))
      } else if (is.character(value)) {
        message(sprintf("      %s: %s", name, value))
      } else {
        message(sprintf("      %s: %s", name, paste(value, collapse = ", ")))
      }
    }
  }
}


# Function to update control parameters from defaults
updateControl <- function(control, default.control) {
  if (is.null(default.control)) {
    return(control)
  }

  # use the default setting or update it
  if (!is.null(control)) {
    ctrl.nm <- names(control)
    for (nm in ctrl.nm) {
      default.control[[nm]] <- control[[nm]]
    }
  }

  control <- default.control
  return(control)
}


# Suppose that subjData is only a subset of the subjects in SparseGRM.
# This function is to extract subjects from SparseGRM.
updateSparseGRM <- function(SparseGRM, subjData) {
  # later add another column to specify the relationship degree
  if (any(toupper(colnames(SparseGRM)) != c("ID1", "ID2", "VALUE"))) {
    stop("The header in 'SparseGRMFile' should be c('ID1','ID2','Value')")
  }

  colnames(SparseGRM) <- toupper(colnames(SparseGRM))

  tempGRM1 <- SparseGRM
  tempGRM2 <- data.frame(
    ID1 = tempGRM1$ID2,
    ID2 = tempGRM1$ID1,
    VALUE = tempGRM1$VALUE
  )

  tempGRM <- rbind(tempGRM1, tempGRM2)
  tempGRM <- tempGRM[-1 * which(duplicated(tempGRM)), ]

  ID1 <- tempGRM$ID1
  ID2 <- tempGRM$ID2
  value <- tempGRM$VALUE

  if (any(!is.element(subjData, ID1))) {
    stop("At least one of subjects is not in SparseGRM.")
  }

  location1 <- match(ID1, subjData)
  location2 <- match(ID2, subjData)
  pos <- which(!is.na(location1) & !is.na(location2))
  locations <- rbind(
    location1[pos] - 1, # -1 is to convert R to C++
    location2[pos] - 1
  )

  value <- value[pos]
  nSubj <- length(subjData)
  KinMatListR <- list(
    locations = locations,
    values = value,
    nSubj = nSubj
  )

  return(KinMatListR)
}


# Check and validate output file settings for analysis restart capability
checkOutputFile <- function(
  OutputFile, # Character string specifying the output file path.
  OutputFileIndex, # Character string specifying the index file path for tracking progress.
  AnalysisType, # Character string indicating analysis type ("Marker" or "Region").
  nEachChunk # Integer specifying the number of items per chunk.
) {
  ## The following messages are for 'OutputFileIndex'
  message1 <- paste("This is the output index file for GRAB package to record",
                    "the end point in case users want to restart the analysis.",
                    "Please do not modify this file.")
  message2 <- paste("This is a", AnalysisType, "level analysis.")
  message3 <- paste("nEachChunk =", nEachChunk)
  message5 <- "Have completed the analyses of all chunks."

  if (file.exists(OutputFile)) {
    if (!file.exists(OutputFileIndex)) {
      stop(paste0("'OutputFile' of '", OutputFile, "' has existed. ",
                  "Please use another 'OutputFile' or remove the existing one."))
    } else {
      outIndexData <- read.table(
        OutputFileIndex,
        header = FALSE,
        sep = "\t"
      )

      # Validate index file structure
      if (outIndexData[1, 1] != message1 ||
            outIndexData[2, 1] != message2 ||
            outIndexData[3, 1] != message3) {
        stop(paste0("'OutputFileIndex' of '", OutputFileIndex, "' is not as expected. ",
                    "Probably, it has been modified by user, which is not permitted. ",
                    "Please remove the existing files of 'OutputFile' and 'OutputFileIndex'."))
      }

      lastMessage <- outIndexData[nrow(outIndexData), 1]
      if (lastMessage == message5) {
        # Analysis already completed
        indexChunk <- outIndexData[nrow(outIndexData) - 1, 1]
        indexChunk <- as.numeric(gsub("Have completed the analysis of chunk ", "", indexChunk))
        
        stop(
          "Analysis completed in an earlier run. Results saved in '",
          OutputFile,
          "'. Use a different 'OutputFile' to restart analysis."
        )
      } else {
        # Partial analysis - restart from next chunk
        indexChunk <- lastMessage
        indexChunk <- as.numeric(gsub("Have completed the analysis of chunk ", "", indexChunk))
        
        .message(
          "Part of analysis completed and saved in %s. Restarting from chunk %d",
          OutputFileIndex, indexChunk + 1
        )
      }
    }
  } else {
    # New analysis - start from beginning
    indexChunk <- 0
  }

  return(indexChunk)
}


# Write analysis output to files with progress tracking
writeOutputFile <- function(
  Output, # List of output data to write to files.
  OutputFile, # Character vector of output file paths.
  OutputFileIndex, # Character string specifying the index file path.
  AnalysisType, # Character string indicating analysis type.
  nEachChunk, # Integer specifying the number of items per chunk.
  indexChunk, # Integer indicating the current chunk index.
  Start, # Logical indicating if this is the first output to save.
  End # Logical indicating if this is the last output to save.
) {
  ## The following messages are for 'OutputFileIndex'
  message1 <- paste("This is the output index file for GRAB package to record",
                    "the end point in case users want to restart the analysis.",
                    "Please do not modify this file.")
  message2 <- paste("This is a", AnalysisType, "level analysis.")
  message3 <- paste("nEachChunk =", nEachChunk)
  message4 <- paste("Have completed the analysis of chunk", indexChunk)
  message5 <- "Have completed the analyses of all chunks."

  n1 <- length(Output)
  n2 <- length(OutputFile)

  if (n1 != n2) {
    stop("length(Output) != length(OutputFile)")
  }

  if (n1 != 0) {
    for (i in 1:n1) {
      if (Start) {
        write.table(
          Output[[i]], OutputFile[[i]],
          quote = FALSE, sep = "\t", append = FALSE,
          col.names = TRUE, row.names = FALSE
        )
      } else {
        write.table(
          Output[[i]], OutputFile[[i]],
          quote = FALSE, sep = "\t", append = TRUE,
          col.names = FALSE, row.names = FALSE
        )
      }
    }
  }

  if (Start) {
    write.table(
      c(message1, message2, message3), OutputFileIndex,
      quote = FALSE, sep = "\t", append = FALSE,
      col.names = FALSE, row.names = FALSE
    )
  }

  write.table(
    message4, OutputFileIndex,
    quote = FALSE, sep = "\t", append = TRUE,
    col.names = FALSE, row.names = FALSE
  )

  if (End) {
    write.table(
      message5, OutputFileIndex,
      quote = FALSE, sep = "\t", append = TRUE,
      col.names = FALSE, row.names = FALSE
    )
  }
}


#' Cauchy Combination Test for p-value aggregation
#'
#' Combines multiple p-values using the Cauchy distribution method, which provides
#' analytical p-value calculation under arbitrary dependency structures.
#'
#' @param pvals Numeric vector of p-values to combine (each between 0 and 1).
#'   P-values equal to 1 are automatically adjusted to 0.999. P-values equal
#'   to 0 will cause an error.
#' @param weights Numeric vector of non-negative weights for each p-value.
#'   If NULL, equal weights are used. Must have same length as pvals.
#' @return Single aggregated p-value combining all input p-values.
#' @details
#' The Cauchy Combination Test (CCT) transforms p-values using the inverse Cauchy
#' distribution and combines them with specified weights. This method is particularly
#' powerful because it:
#' \itemize{
#'   \item Works under arbitrary dependency structures
#'   \item Provides exact analytical p-values (no simulation needed)
#'   \item Maintains good power properties across different scenarios
#' }
#'
#' **Special Cases:**
#' - If any p-value equals 0, returns 0 immediately
#' - P-values equal to 1 are adjusted to 0.999 with a warning
#' - Very small p-values (< 1e-16) receive special numerical treatment
#'
#' @examples
#' # Basic usage with equal weights
#' pvalues <- c(0.02, 0.0004, 0.2, 0.1, 0.8)
#' CCT(pvals = pvalues)
#'
#' # Usage with custom weights
#' weights <- c(2, 3, 1, 1, 1)
#' CCT(pvals = pvalues, weights = weights)
#'
#' @references
#' Liu, Y., & Xie, J. (2020). Cauchy combination test: a powerful test
#' with analytic p-value calculation under arbitrary dependency structures.
#' \emph{Journal of the American Statistical Association}, 115(529), 393-402.
#' \doi{10.1080/01621459.2018.1554485}
#'
CCT <- function(
  pvals,
  weights = NULL
) {
  # Validate p-values for missing values
  if (sum(is.na(pvals)) > 0) {
    stop("Cannot have NAs in the p-values!")
  }

  # Ensure all p-values are in valid range [0, 1]
  if ((sum(pvals < 0) + sum(pvals > 1)) > 0) {
    stop("All p-values must be between 0 and 1!")
  }

  # Handle edge cases: p-values exactly 0 or 1
  is.zero <- (sum(pvals == 0) >= 1)
  which.one <- which(pvals == 1)
  is.one <- (length(which.one) >= 1)

  if (is.zero && is.one) {
    stop("Cannot have both 0 and 1 p-values!")
  }
  if (is.zero) {
    return(0)  # If any p-value is 0, combined p-value is 0
  }
  if (is.one) {
    warning("There are p-values that are exactly 1!")
    pvals[which.one] <- 0.999  # Adjust p-values of 1 to avoid numerical issues
  }

  # Validate and standardize weights
  if (is.null(weights)) {
    weights <- rep(1 / length(pvals), length(pvals))  # Equal weights
  } else if (length(weights) != length(pvals)) {
    stop("The length of weights should be the same as that of the p-values!")
  } else if (sum(weights < 0) > 0) {
    stop("All the weights must be positive!")
  } else {
    weights <- weights / sum(weights)  # Normalize weights to sum to 1
  }

  # Calculate CCT test statistic with special handling for very small p-values
  is.small <- (pvals < 1e-16)
  if (sum(is.small) == 0) {
    # Standard calculation when no extremely small p-values
    cct.stat <- sum(weights * tan((0.5 - pvals) * pi))
  } else {
    # Special handling for very small p-values to avoid numerical issues
    cct.stat <- sum((weights[is.small] / pvals[is.small]) / pi)
    cct.stat <- cct.stat + sum(weights[!is.small] * tan((0.5 - pvals[!is.small]) * pi))
  }

  # Calculate final p-value from CCT statistic
  if (cct.stat > 1e+15) {
    # Asymptotic approximation for very large test statistics
    pval <- (1 / cct.stat) / pi
  } else {
    # Standard Cauchy distribution calculation
    pval <- 1 - pcauchy(cct.stat)
  }

  return(pval)
}
