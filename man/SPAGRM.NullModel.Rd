% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPAGRM.R
\name{SPAGRM.NullModel}
\alias{SPAGRM.NullModel}
\title{Fit a SPAGRM Null Model}
\usage{
SPAGRM.NullModel(
  ResidMatFile,
  SparseGRMFile,
  PairwiseIBDFile,
  control = list(MaxQuantile = 0.75, MinQuantile = 0.25, OutlierRatio = 1.5,
    ControlOutlier = TRUE, MaxNuminFam = 5, MAF_interval = c(1e-04, 5e-04, 0.001, 0.005,
    0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5))
)
}
\arguments{
\item{ResidMatFile}{A file path (character) or data.frame containing residuals. If a file path, it should point to a tab-delimited file with two columns: 'SubjID' (subject IDs) and 'Resid' (residual values). If a data.frame, it should have the same structure with columns named 'SubjID' and 'Resid'.}

\item{SparseGRMFile}{A file path (character) to a sparse genetic relationship matrix (GRM) file. This file should be generated using the \code{getSparseGRM()} function and contain three columns: 'ID1', 'ID2', and 'Value' representing the genetic relationships between pairs of individuals.}

\item{PairwiseIBDFile}{A file path (character) to a pairwise identity-by-descent (IBD) file. This file should be generated using the \code{getPairwiseIBD()} function and contain five columns: 'ID1', 'ID2', 'pa', 'pb', and 'pc' representing IBD probabilities between pairs of individuals.}

\item{control}{A list of control parameters for the null model fitting process. Available options include:}
}
\value{
A SPAGRM null model object
}
\description{
Fit a SPAGRM Null Model
}
