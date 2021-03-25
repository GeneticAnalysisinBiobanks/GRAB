
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
#'   \item \code{IDsToIncludeFile}: a file of marker IDs to include, one column (no header) 
#'   \item \code{IDsToExcludeFile}: a file of marker IDs to exclude, one column (no header) 
#'   \item \code{RangesToIncludeFile}: a file of ranges to include, three columns (no header): chromosome, start position, end position
#'   \item \code{RangesToExcludeFile}: a file of ranges to exclude, three columns (no header): chromosome, start position, end position
#'   \item \code{AlleleOrder}: a character, "ref-first" or "alt-first", to determine whether the REF/major allele should appear first or second. Default is "alt-first" for PLINK and "ref-first" for BGEN. If the analysis results show the REF allele frequencies of most markers are > 0.5, you might should change this option. NOTE, if you use plink2 to convert Plink file to BGEN file, you probably need to set it as 'ref-first'.
#'   }
#'   }
#'   \item{GRAB.NullModel}{
#'     \itemize{
#'     \item \code{GenoFile}: "prefix.bgen"; 
#'     \item \code{GenoFileIndex}: "prefix.bgen.bgi" or c("prefix.bgen.bgi", "prefix.bgen.samples").
#'     \item If only one element is given for \code{GenoFileIndex}, then we assume it should be "prefix.bgen.bgi". 
#'     \item Sometimes, BGEN file does not include sample identifiers, and thus file of "prefix.bgen.samples" is required.
#'     \item NOTE that "prefix.bgen.samples" should be of only one column with the column name of "GRAB_BGEN_SAMPLE" (case insensitive). One example can be found in \code{system.file("extdata", "example_bgen_1.2_8bits.bgen.samples", package = "GRAB")}.
#'     \item If you are not sure if sample identifiers are in BGEN file, you can try function \code{?checkIfSampleIDsExist}.
#'     }
#'   }
#'   \item{GRAB.Marker}{
#'   \itemize{
#'   \item \code{IDsToIncludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{IDsToExcludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{RangesToIncludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{RangesToExcludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{AlleleOrder}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{impute_method}: a character, "mean", "minor", or "drop". If "mean", impute genotype using 2 * AlleleFreq; if "minor", impute genotype using minor alleles; if "drop", drop the subject whose genotype is missing.
#'   \item \code{missing_cutoff}: a numeric value (default: 0.15). Any variant with missing rate > this value will be excluded from analysis.  
#'   \item \code{min_maf_marker}: a numeric value (default: 0.001). Any variants with MAF < this value will be excluded from analysis.  
#'   \item \code{min_mac_marker}: a numeric value (default: 20). Any variants with MAC < this value will be excluded from analysis.  
#'   \item \code{nMarkersEachChunk}: number of markers (default: 10000) in one chunk to be outputted
#'   }
#'   }
#'   \item{GRAB.Region}{
#'   \itemize{
#'   \item \code{IDsToIncludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{IDsToExcludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{RangesToIncludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{RangesToExcludeFile}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{AlleleOrder}: please refer to section \code{GRAB.ReadGeno}.
#'   \item \code{impute_method}: a character, "mean", "minor", or "drop". If "mean", impute genotype using 2 * AlleleFreq; if "minor", impute genotype using minor alleles; if "drop", drop the subject whose genotype is missing.
#'   \item \code{missing_cutoff}: a numeric value (default: 0.15). Any variant with missing rate > this value will be excluded from analysis.  
#'   \item \code{max_maf_region}: a numeric value (default: 0.001). Any variants with MAF < this value will be excluded from analysis.  
#'   \item \code{nRegionsEachChunk}: number of markers (default: 10000) in one chunk to be outputted
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
                                nMarkersEachChunk = 10000)
  
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
  
  # specific default control setting for different approaches
  if(NullModelClass == "POLMM_NULL_Model")
    control = checkControl.Marker.POLMM(control)    # This function is in 'POLMM.R'
  
  if(NullModelClass == "SPACox_NULL_Model")
    control = checkControl.Marker.SPACox(control)
  
  
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
  
  # uniform default control setting for marker-level analysis
  default.region.control = list(impute_method = "fixed",  
                                missing_cutoff = 0.15,
                                max_maf_region = 0.01,
                                max_mem_region = 4)
  
  control = updateControl(control, default.region.control)
  
  # check if 'control' is reasonable
  if(!control$impute_method %in% c("mean", "minor", "drop"))
    stop("control$impute_method should be 'mean', 'minor', or 'drop'.")
  
  if(!is.numeric(control$missing_cutoff) | control$missing_cutoff < 0 | control$missing_cutoff > 0.5)
    stop("control$missing_cutoff should be a numeric value ranging from 0 to 0.5.")
  
  if(!is.numeric(control$max_maf_region) | control$max_maf_region < 0 | control$max_maf_region > 0.3)
    stop("control$min_maf_marker should be a numeric value ranging from 0 to 0.3.")
  
  if(!is.numeric(control$max_mem_region) | control$max_mem_region <= 0)
    stop("control$min_mac_marker should be a numeric value > 0.")
  
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

