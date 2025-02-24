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
  
  if(!is.null(control$AlleleOrder))
    if(control$AlleleOrder != "ref-first" & control$AlleleOrder != "alt-first")
      stop("control$AlleleOrder should be 'ref-first' or 'alt-first'.")
  
  if(is.null(control$ImputeMethod)){
    control$ImputeMethod = "none"
  }
  
  if(!control$ImputeMethod %in% c("none", "bestguess", "mean"))
    stop("control$ImputeMethod should be 'none', 'bestguess', or 'mean'.")
  
  if(is.null(control$AllMarkers))
    control$AllMarkers = FALSE
  
  FileType = c("IDsToIncludeFile", "IDsToExcludeFile", "RangesToIncludeFile", "RangesToExcludeFile")
  
  # check if the files specified exist
  for(ft in FileType){
    if(ft %in% names(control)){
      file = control[[ft]]
      if(!file.exists(file))
        stop(paste0("Cannot find the file of ",file,"..."))
    }
  }
  
  return(control)
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
                                omp_num_threads = data.table::getDTthreads())    # if 0, value is from omp_get_num_threads(). Not supported on 2022-02-07
  
  control = updateControl(control, default.marker.control)
  
  # check if argument of 'control' is reasonable
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
    
  method = gsub("_NULL_Model", "", NullModelClass)  # updated on 2022-04-26: "POLMM_NULL_Model" -> "POLMM"
  
  textToParse = paste0("control = checkControl.Marker.", method, "(control)")
  eval(parse(text = textToParse))
  
  # specific default control setting for different approaches
  # if(NullModelClass == "POLMM_NULL_Model")
  #   control = checkControl.Marker.POLMM(control)    # This function is in 'POLMM.R'
  # 
  # if(NullModelClass == "SPACox_NULL_Model")
  #   control = checkControl.Marker.SPACox(control)   # This function is in 'SPACox.R'
  
  
  # print control list 
  print("The below is the list of control parameters used in marker-level genetic association analysis.")
  print(control)
  
  return(control)
}


checkControl.Region = function(control)
{
  # check if control is an R list
  if(!is.null(control))
    if(class(control) != "list")
      stop("If specified, the argument of 'control' should be an R 'list'.")
  
  # uniform default control setting for region-level analysis
  default.region.control = list(impute_method = "minor",  
                                missing_cutoff = 0.15,
                                min_mac_region = 5,
                                max_markers_region = 100,
                                r.corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.5, 1),
                                weights.beta = c(1, 25),
                                omp_num_threads = data.table::getDTthreads(),
                                min_nMarker = 3)
  
  control = updateControl(control, default.region.control)
  
  # check if 'control' is reasonable
  if(!control$impute_method %in% c("mean", "minor", "drop"))
    stop("control$impute_method should be 'mean', 'minor', or 'drop'.")
  
  if(!is.numeric(control$missing_cutoff) | control$missing_cutoff < 0 | control$missing_cutoff > 0.5)
    stop("control$missing_cutoff should be a numeric value ranging from 0 to 0.5.")
  
  if(!is.numeric(control$min_mac_region) | control$min_mac_region < 0)
    stop("control$min_mac_region should be a numeric value >= 0.")
  
  if(!is.numeric(control$max_markers_region) | control$max_markers_region < 50)
    stop("control$max_markers_region should be a integer >= 50.")
  
  if(!is.numeric(control$r.corr) | min(control$r.corr) < 0 | max(control$r.corr) > 1)
    stop("control$r.corr should be a numeric vector whose elements are between 0 and 1.")
  
  if(!is.numeric(control$weights.beta) | length(control$weights.beta) != 2 | min(control$weights.beta) < 0)
    stop("control$weights.beta should be a numeric vector with two non-negative elements.")
  
  if(!is.numeric(control$min_nMarker) | control$min_nMarker <= 0)
    stop("control$min_nMarker should be a positive integer.")
  
  return(control)
}

# check control list in null model fitting
checkControl.NullModel = function(control, method, traitType, optionGRM)
{
  # check if control is an R list
  if(!is.null(control))
    if(class(control) != "list")
      stop("If specified, the argument of 'control' should be an R 'list'.")
  
  # if(method == "POLMM"){control = checkControl.NullModel.POLMM(control, traitType, optionGRM)}
  textToParse = paste0("control = checkControl.NullModel.", method, "(control, traitType, optionGRM)")
  eval(parse(text = textToParse))
  
  # to be updated for other methods
  #
  # ------------------------------
  #
  # to be updated for other methods
  
  cat("The below are the list of control parameters used in null model fitting.\n")
  print(control)
  
  return(control)
}

