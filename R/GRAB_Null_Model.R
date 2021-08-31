#' Fit a null model to estimate parameters and residuals
#'
#' We fit a null model including response variable, covariates, and Genetic Relationship Matrix (GRM, if needed) to estimate parameters and residuals. 
#' 
#' @param formula a formula object, with the response on the left of a ~ operator and the covariates on the right. Do not add a column of intercept (i.e. a vector of ones) on the right. Missing values should be denoted by NA and the corresponding samples will be removed from analysis. Other values (e.g. -9, -999) will be treated as ordinary numeric values in analysis.
#' @param data a data.frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data.frame), containing the variables in formula. Neither a matrix nor an array will be accepted.
#' @param subset a specification of the rows to be used: defaults to all rows. This can be any valid indexing vector for the rows of data or if that is not supplied, a data frame made up of the variables used in formula.
#' @param subjData a character vector of subject IDs. Its order should be the same as the subject order in the formula and data (before subset process). 
#' @param method a character: "SPACox" (check \code{\link{GRAB.SPACox}}), "POLMM" (check \code{\link{GRAB.POLMM}}), "SPAGE" (will be supported later), "SAIGE" (will be supported later), or "GATE" (will be supported later).
#' @param traitType a character: "binary", "ordinal" (check \code{\link{GRAB.POLMM}}), "quantitative", or "time-to-event" (check \code{\link{GRAB.SPACox}}).
#' @param GenoFile a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param GenoFileIndex additional index files corresponding to the \code{GenoFile}. If Null (default), the same prefix as GenoFile is used. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param SparseGRMFile a character of sparseGRM file. An example is \code{system.file("SparseGRM","SparseGRM.txt",package="GRAB")}
#' @param control a list of parameters for controlling the model fitting process. For more details, please check \code{?GRAB.control}. 
#' @param ... other arguments passed to or from other methods. 
#' @return an R object with a class of "XXXXX_NULL_Model" in which XXXXX is the 'method' used in analysis. 
#' @details 
#' \code{GRAB} package uses score tests which consists of two steps. In Step 1, function \code{GRAB.NullModel} fits a null model including response variable, covariates, and Genetic Relationship Matrix (GRM, if needed).
#' In Step 2, functions \code{\link{GRAB.Marker}} and \code{\link{GRAB.Region}} perform genome-wide marker-level analysis and region-level analysis, respectively. 
#' Step 1 fits a null model to get an R object, which can be passed to Step 2 for association testing. Functions of \code{\link{save}} and \code{\link{load}} can save and load the object.
#' 
#' \describe{
#' \code{GRAB} package includes multiple methods for large-scale genome-wide analysis of a wide variety of phenotypes as follows.
#' \itemize{
#'   \item \code{SAIGE}: Support traitType = "quantitative" or "binary". Will be added later.
#'   \item \code{POLMM}: Support traitType = "ordinal". Check \code{\link{GRAB.POLMM}} for more details.
#'   \item \code{SPACox}: Support traitType = "time-to-event". Check \code{\link{GRAB.SPACox}} for more details.
#' }
#' }
#' 
#' \code{GRAB} package supports both dense and sparse GRM to adjust for family relatedness. 
#' If Dense GRM is used, then \code{GenoFile} is required to get genome-wide marker information for GRM.
#' If Sparse GRM is used, then users should first use \code{\link{getSparseGRM}} to get a \code{SparseGRMFile}, which is required in analysis.
#' 
#' \describe{
#' Argument \code{control} includes a list of parameters for controlling the null model fitting process.  
#'   \itemize{
#'     \item \code{memoryChunk}: Size (Gb) for each memory chunk when reading in PLINK file [default=2].
#'     \item \code{seed}: An integer as a random seed [default=12345678].
#'     \item \code{tracenrun}: Number of runs for trace estimator [default=30].
#'     \item \code{maxiter}: Maximum number of iterations used to fit the null POLMM [default=100].
#'     \item \code{tolBeta}: Positive tolerance: the iterations converge when |beta - beta_old| / (|beta| + |beta_old| + tolBeta) < tolBeta [default=0.001].
#'     \item \code{tolTau}: Positive tolerance: the iterations converge when |tau - tau_old| / (|tau| + |tau_old| + tolTau) < tolTau [default=0.002].
#'     \item \code{tau}: Initial value of the variance component (tau) [default=0.2].
#'     \item \code{maxiterPCG}: Maximum number of iterations for PCG to converge [default=100].
#'     \item \code{tolEps}: Positive tolerance for PCG to converge [default=1e-6].
#'     \item \code{minMafVarRatio}: Minimal value of MAF cutoff to select markers (from Plink file) to estimate variance ratio [default=0.1].
#'     \item \code{maxMissingVarRatio}: Maximal value of missing rate cutoff to select markers (from Plink file) to estimate variance ratio [default=0.1].
#'     \item \code{nSNPsVarRatio}: Initial number of the selected markers to estimate variance ratio [default=20], the number will be automatically added by 10 until the coefficient of variantion (CV) of the variance ratio estimate is below CVcutoff.
#'     \item \code{CVcutoff}: Minimal cutoff of coefficient of variantion (CV) for variance ratio estimation [default=0.0025].
#'     \item \code{LOCO}: Whether to apply the leave-one-chromosome-out (LOCO) approach [default=TRUE].
#'     \item \code{numThreads}: Number of threads (CPUs) to use. Only valid if dense GRM is used, default is "auto", that is, RcppParallel::defaultNumThreads() [default="auto"].
#'     \item \code{stackSize}: Stack size (in bytes) to use for worker threads. For more details, check help(RcppParallel::setThreadOptions) [default="auto"].
#'     \item \code{grainSize}: Grain size of a parallel algorithm sets a minimum chunk size for parallelization. In other words, at what point to stop processing input on separate threads [default=1].
#'     \item \code{minMafGRM}: Minimal value of MAF cutoff to select markers (from Plink files) to construct dense GRM [default=0.01].
#'     \item \code{maxMissingGRM}: Maximal value of missing rate to select markers (from Plink files) to construct dense GRM [default=0.1].
#'     \item \code{showInfo }: Whether to show more detailed information for trouble shooting [default=TRUE].
#'     \item \code{onlyCheckTime}: Not fit the null model, only check the computation time of reading Plink files and running 30 KinbVec() functions [default=FALSE].
#'   }
#' }
#' @examples
#' # For POLMM method (ordinal categorical data analysis while adjusting for sample relatedness)
#' # Step 1(a): fit a null model using a dense (full) GRM
#' PhenoData = read.table(system.file("extdata", "example.pheno", package = "GRAB"), header = T)
#' GenoFile = system.file("extdata", "example.bed", package = "GRAB")
#' obj.POLMM = GRAB.NullModel(factor(Ordinal) ~ Cova1 + Cova2,
#'                            data = PhenoData, subjData = PhenoData$IID, method = "POLMM", traitType = "ordinal",
#'                            GenoFile = GenoFile,
#'                            control = list(showInfo = FALSE, LOCO = FALSE, tolTau = 0.2, tolBeta = 0.1))
#' 
#' names(obj.POLMM)
#' obj.POLMM$tau    # 1.820102
#'
#' # Step 1(b): fit a null model using a sparse GRM
#' # First use getSparseGRM() function to get a sparse GRM file
#' PhenoData = read.table(system.file("extdata", "example.pheno", package = "GRAB"), header = T)
#' GenoFile = system.file("extdata", "example.bed", package = "GRAB")
#' SparseGRMFile =  system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#' obj.POLMM = GRAB.NullModel(factor(Ordinal) ~ Cova1 + Cova2,
#'                            data = PhenoData, subjData = PhenoData$IID, method = "POLMM", traitType = "ordinal",
#'                            GenoFile = GenoFile,
#'                            SparseGRMFile = SparseGRMFile,
#'                            control = list(showInfo = FALSE, LOCO = FALSE, tolTau = 0.2, tolBeta = 0.1))
#' 
#' names(obj.POLMM)
#' obj.POLMM$tau    # 1.870175
#' 
#' # save(obj.POLMM, "obj.POLMM.RData")  # save the object for analysis in step 2
#' 
#' # For SPACox method (time-to-event data analysis)
#' 
#' @export
#' @import survival, data.table
GRAB.NullModel = function(formula,
                          data = NULL,
                          subset = NULL,
                          subjData,
                          method = "SPACox",
                          traitType = "time-to-event",  # "binary", "ordinal", "quantitative", "time-to-event"
                          GenoFile = NULL,
                          GenoFileIndex = NULL,
                          SparseGRMFile = NULL,
                          control = NULL,
                          ...)
{
  if(missing(subjData))
    stop("Argument 'subjData' is required to specify the subjects IDs in 'formula' and/or 'data'.")
  
  Call = match.call()
  
  # check 'control.R'
  control = checkControl.NullModel(control, method, traitType)
  
  #### START: formula.R
    
  #### input: formula, data, subset, subjData
  #### output: response, designMat, subjData
  
  mf <- match.call(expand.dots = FALSE)
  
  m <- match(x = c("formula", "data", "subset", "subjData"), 
             table = names(mf), nomatch = 0L)
  
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  
  mt = attr(x = mf, which = "terms")
  
  response = model.response(mf)
  designMat = model.matrix(object = mt, data = mf)
  subjData = model.extract(mf, "subjData")
  
  if(colnames(designMat)[1] == "(Intercept)")
    designMat = designMat[,-1,drop=F]
  
  nData = length(subjData)
  cat("Number of subjects in 'formula':\t", nData,"\n")
  
  if(any(duplicated(subjData))) 
    stop("Duplicated subject IDs in 'subjData' is not supported!")
  
  #### END: formula.R
  
  if(method %in% c("POLMM", "SAIGE", "GATE"))
    optionGRM = handleGRM(GenoFile, GenoFileIndex, SparseGRMFile, subjData)  # Check 'SparseGRM.R'
  
  if(method == "POLMM"){
    # Check 'POLMM.R'
    objNull = fitNullModel.POLMM(response, designMat, subjData, control, optionGRM)
  }
  
  if(method == "SPACox"){
    # Check 'SPACox.R'
    objNull = fitNullModel.SPACox(response, designMat, subjData, control, ...)
  }
  
  objNull$subjData = subjData
  
  objNull$Call = Call;
  objNull$sessionInfo = sessionInfo()
  objNull$time = Sys.time()
  objNull$control = control
  
  print(paste0("Complete the null model fitting in package GRAB: ", objNull$time))
  
  return(objNull)
}






