% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GRAB_Null_Model.R
\name{GRAB.NullModel}
\alias{GRAB.NullModel}
\title{GRAB: NULL model fitting}
\usage{
GRAB.NullModel(
  formula,
  data = NULL,
  subset = NULL,
  subjData,
  method = "SPACox",
  traitType = "time-to-event",
  GenoFile,
  GenoFileIndex = NULL,
  SparseGRMFile,
  control = NULL,
  ...
)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a ~ operator and the covariates on the right. Do not add a column of intercept (e.g. a vector of ones) on the right. Missing values should be denoted by NA and the corresponding samples will be removed from analysis.}

\item{data}{a data.frame in which to interpret the variables named in the formula, or in the subset argument. Check \code{?model.frame} for more details.}

\item{subset}{a specification of the rows to be used: defaults to all rows. This can be any valid indexing vector for the rows of data or if that is not supplied, a data frame made up of the variables used in formula. Check \code{?model.frame} for more details.}

\item{subjData}{a character vector of subject IDs. Its order should be the same as the subject order in the formula and data.}

\item{method}{a character: "SPACox", "SPAGE", "SAIGE", "POLMM", or "GATE"}

\item{traitType}{a character: "binary", "ordinal", "quantitative", or "time-to-event"}

\item{GenoFile}{a character of genotype file. Three types of genotype files are supported: PLINK, BGEN, and VCF. More details are in \code{?GRAB.ReadGeno}.}

\item{GenoFileIndex}{additional index files corresponding to the "GenoFile". If Null (default), the same prefix as GenoFile is used. More details are in \code{?GRAB.ReadGeno}.}

\item{SparseGRMFile}{a character of sparseGRM file. An example is \code{system.file("SparseGRM","SparseGRM.txt",package="GRAB")}}

\item{control}{a list of parameters for controlling the \code{GRAB.NullModel()}. For more details, please check \code{?GRAB.control}.}

\item{...}{other arguments passed to or from other methods.}
}
\value{
an R object with a class of "XXXXX_NULL_Model" in which XXXXX is the 'method' used in analysis. The object will be used in \code{GRAB.Marker()} and \code{GRAB.Region()}. Functions of \code{save()} and \code{load()} can be used to save and load the object for future usage.
}
\description{
GRAB package uses score test for GWAS: in step 1, we fit a null model (check \code{?GRAB.NullModel}) including response variable, covariates, and GRM (if needed). In step 2, we perform score test for marker-level analysis (check \code{?GRAB.Marker}) and region-level analysis (check \code{?GRAB.Region}).
}
\examples{
# Examples have been put in the specific help pages for specific methods. 
# If you want to use "SPACox" method, please check ?GRAB.SPACox for more details.
# If you want to use "POLMM" method, please check ?GRAB.POLMM for more details.
}
