
#' Information about the control in GRAB package
#' 
#' Information about the control in GRAB package
#' 
#' @details
#' In GRAB package, we use \code{control} to specify parameters in multiple functions. Here, we only list generic parameters. 
#' For the parameters only used in specific method, please refer to its own help page. 
#' For example, check \code{?GRAB.SPACox} for more information about SPACox method.
#' \describe{
#'   \item{GRAB.ReadGeno}{
#'   We support two include files of (IDsToIncludeFile, RangesToIncludeFile) and two exclude files of (IDsToExcludeFile, RangesToExcludeFile), but do not support both include and exclude files.
#'   \itemize{
#'   \item \code{IDsToIncludeFile}: a file of marker IDs to include, one column (no header). Check \code{system.file("extdata", "IDsToInclude.txt", package = "GRAB")} for an example. 
#'   \item \code{IDsToExcludeFile}: a file of marker IDs to exclude, one column (no header). 
#'   \item \code{RangesToIncludeFile}: a file of ranges to include, three columns (no headers): chromosome, start position, end position. Check \code{system.file("extdata", "RangesToInclude.txt", package = "GRAB")} for an example.
#'   \item \code{RangesToExcludeFile}: a file of ranges to exclude, three columns (no headers): chromosome, start position, end position.
#'   \item \code{AlleleOrder}: a character, "ref-first" or "alt-first", to determine whether the REF/major allele should appear first or second. Default is "alt-first" for PLINK and "ref-first" for BGEN. If the analysis results show the ALT allele frequencies of most markers are > 0.5, you might should change this option. NOTE, if you use plink2 to convert Plink file to BGEN file, you probably need to set it as 'ref-first'.
#'   }
#'   }
#'   \item{GRAB.NullModel}{
#'   The following parameters are to handle dense GRM for mixed model approaches (\code{SAIGE}, \code{GATE}, and \code{POLMM})
#'     \itemize{
#'     \item \code{memoryChunk}: Size (Gb) for each memory chunk when reading in Plink file [default=2].
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
#'     }
#'   }
#'   \item{GRAB.Marker}{
#'   \itemize{
#'   The following parameters are to specify Markers to analyze.
#'   \item \code{IDsToIncludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{IDsToExcludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{RangesToIncludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{RangesToExcludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{AlleleOrder}: please refer to section \code{GRAB.ReadGeno}.
#'   }
#'   \itemize{
#'   The following parameter are for score test.
#'   \item \code{omp_num_threads}: a numeric value (default: value from data.table::getDTthreads()) to specify the number of threads in OpenMP for parallel computation.
#'   \item \code{impute_method}: a character, "mean", "minor", or "drop". If "mean", impute genotype using 2 * AlleleFreq; if "minor", impute genotype using minor alleles; if "drop", drop the subject whose genotype is missing.
#'   \item \code{missing_cutoff}: a numeric value (default: 0.15). Any variant with missing rate > this value will be excluded from analysis.  
#'   \item \code{min_maf_marker}: a numeric value (default: 0.001). Any variants with MAF < this value will be excluded from analysis.  
#'   \item \code{min_mac_marker}: a numeric value (default: 20). Any variants with MAC < this value will be excluded from analysis.  
#'   \item \code{nMarkersEachChunk}: number of markers (default: 10000) in one chunk to be outputted
#'   }
#'   }
#'   \item{GRAB.Region}{
#'   \itemize{
#'   \item \code{xxx}: xxx
#'   }
#'   }
#' }
#' @export
GRAB.control = function(){
  return("Check ?GRAB.control for more information about the 'control' in package GRAB.")
}

# update 'control' or use 'default.control' if the corresponding elements are not specified 
updateControl = function(control, default.control)
{
  if(is.null(default.control))
    return(control)
  
  # use the default setting or update it
  if(!is.null(control)){
    ctrl.nm = names(control)
    for(nm in ctrl.nm){
      default.control[[nm]] = control[[nm]]
    }
  }
  
  control = default.control
  return(control)
}

checkControl.ReadGeno = function(control)
{
  # check if control is an R list
  if(!is.null(control))
    if(class(control) != "list")
      stop("If specified, the argument of 'control' should be an R 'list'.")
  
  if(!is.null(control$allele.order))
    if(control$allele.order != "ref-first" & control$allele.order != "alt-first")
      stop("control$allele.order should be 'ref-first' or 'alt-first'.")
  
  FileType = c("IDsToIncludeFile", "IDsToExcludeFile", "RangesToIncludeFile", "RangesToExcludeFile")
  
  # check if the files specified exist
  for(ft in FileType){
    if(ft %in% names(control)){
      file = control[[ft]]
      if(!file.exists(file))
        stop(paste0("Cannot find the file of ",file,"..."))
    }
  }
}

checkControl.Marker = function(control, NullModelClass)
{
  # check if control is an R list
  if(!is.null(control))
    if(class(control) != "list")
      stop("If specified, the argument of 'control' should be an R 'list'.")
  
  # uniform default control setting for marker-level analysis
  default.marker.control = list(impute_method = "mean",  
                                missing_cutoff = 0.15,
                                min_maf_marker = 0.001,
                                min_mac_marker = 20,
                                nMarkersEachChunk = 10000,
                                omp_num_threads = data.table::getDTthreads())    # if 0, value is from omp_get_num_threads()
  
  control = updateControl(control, default.marker.control)
  
  # check if 'control' is reasonable
  if(!control$impute_method %in% c("mean", "minor", "drop"))
    stop("control$impute_method should be 'mean', 'minor', or 'drop'.")
  
  if(!is.numeric(control$missing_cutoff) | control$missing_cutoff < 0 | control$missing_cutoff > 0.5)
    stop("control$missing_cutoff should be a numeric value ranging from 0 to 0.5.")
  
  if(!is.numeric(control$min_maf_marker) | control$min_maf_marker < 0 | control$min_maf_marker > 0.1)
    stop("control$min_maf_marker should be a numeric value ranging from 0 to 0.1.")
  
  if(!is.numeric(control$min_mac_marker) | control$min_mac_marker < 0 | control$min_mac_marker > 100)
    stop("control$min_mac_marker should be a numeric value ranging from 0 to 100.")
  
  if(!is.numeric(control$nMarkersEachChunk) | control$nMarkersEachChunk < 1e3 | control$nMarkersEachChunk > 1e5)
    stop("control$nMarkersEachChunk should be a numeric value ranging from 1e3 to 1e5.")
  
  if(control$omp_num_threads < 0)
    stop("control$omp_num_threads should be a positive integral value.")
  
  # specific default control setting for different approaches
  if(NullModelClass == "POLMM_NULL_Model")
    control = checkControl.Marker.POLMM(control)    # This function is in 'POLMM.R'
  
  if(NullModelClass == "SPACox_NULL_Model")
    control = checkControl.Marker.SPACox(control)   # This function is in 'SPACox.R'
  
  
  # print control list 
  print("The below is the list of control parameters used in marker-level genetic association analysis.")
  print(control)
  
  return(control)
}


checkControl.Region = function(control, NullModelClass)
{
  # check if control is an R list
  if(!is.null(control))
    if(class(control) != "list")
      stop("If specified, the argument of 'control' should be an R 'list'.")
  
  # uniform default control setting for region-level analysis
  default.region.control = list(impute_method = "minor",  
                                missing_cutoff = 0.15,
                                max_maf_region = 0.01,
                                min_mac_region = 10,
                                max_mem_region = 4,
                                r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),
                                weights.beta = c(1, 25),
                                omp_num_threads = data.table::getDTthreads())
  
  control = updateControl(control, default.region.control)
  
  # check if 'control' is reasonable
  if(!control$impute_method %in% c("mean", "minor", "drop"))
    stop("control$impute_method should be 'mean', 'minor', or 'drop'.")
  
  if(!is.numeric(control$missing_cutoff) | control$missing_cutoff < 0 | control$missing_cutoff > 0.5)
    stop("control$missing_cutoff should be a numeric value ranging from 0 to 0.5.")
  
  if(!is.numeric(control$max_maf_region) | control$max_maf_region < 0 | control$max_maf_region > 0.3)
    stop("control$max_maf_region should be a numeric value ranging from 0 to 0.3.")
  
  if(!is.numeric(control$min_mac_region) | control$min_mac_region < 0)
    stop("control$min_mac_region should be a numeric value >= 0.")
  
  if(!is.numeric(control$max_mem_region) | control$max_mem_region <= 0)
    stop("control$min_mac_marker should be a numeric value > 0.")
  
  if(!is.numeric(control$r.corr) | min(control$r.corr) < 0 | max(control$r.corr) > 1)
    stop("control$r.corr should be a numeric vector whose elements are between 0 and 1.")
  
  if(!is.numeric(control$weights.beta) | length(control$weights.beta) != 2 | min(control$weights.beta) < 0)
    stop("control$weights.beta should be a numeric vector with two non-negative elements.")
  
  # specific default control setting for different approaches
  if(NullModelClass == "POLMM_NULL_Model")
    control = checkControl.Region.POLMM(control)    # This function is in 'POLMM.R'
  
  if(NullModelClass == "SPACox_NULL_Model")
    control = checkControl.Region.SPACox(control)
  
  # print control list 
  print("The below is the list of control parameters used in region-level genetic association analysis.")
  print(control)
  
  return(control)
}

# check control list in null model fitting
checkControl.NullModel = function(control, method, traitType)
{
  # check if control is an R list
  if(!is.null(control))
    if(class(control) != "list")
      stop("If specified, the argument of 'control' should be an R 'list'.")
  
  # SPACox method
  if(method == "SPACox"){
    if(traitType != "time-to-event")
      stop("For method of 'SPACox', only traitType of 'time-to-event' is supported.")
    control = checkControl.NullModel.SPACox(control)
  }
  
  # POLMM method
  if(method == "POLMM"){
    if(traitType != "ordinal")
      stop("For method of 'POLMM', only traitType of 'ordinal' is supported.")
    # The following function is in 'POLMM.R'
    control = checkControl.NullModel.POLMM(control)
  }
  
  # to be updated for other methods
  #
  # ------------------------------
  #
  # to be updated for other methods
  
  print("The below are the list of control parameters used in null model fitting.")
  print(control)
  
  return(control)
}

