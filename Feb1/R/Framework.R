
# update 'control' or use 'default.control' if the corresponding elements are not specified 
#' @export
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

#' @export
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

#' @export
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
                                omp_num_threads = data.table::getDTthreads())    
  
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
    
  # Specific check for SPAmixPlus
  control = checkControl.Marker.SPAmixPlus(control)
  
  return(control)
}


#' A lower function to make groups based on phenotype
#' @export
makeGroup = function(yVec)
{
  # yVec is categorical data
  m1 = length(unique(yVec))
  if(m1 <= 10)
    Group = as.numeric(as.factor(yVec))-1  # from 0 to (m1-1)
  
  # yVec is quantitative data
  if(length(unique(yVec)) > 10)
    Group = floor((rank(yVec, ties.method = "max")-1) / length(yVec) * 10)  # from 0 to 9
  
  return(Group)
}

#' @export
checkObjNull = function(objNull)
{
  NullModelClass = class(objNull)
  nm = names(objNull)
  
  if(!NullModelClass %in% c("SPAmixPlus_NULL_Model"))      
  {
    stop('class(objNull) should be "SPAmixPlus_NULL_Model"')
  }
    
  if(any(!c("subjData", "N") %in% nm))
  {
    stop("c('subjData', 'N') should be in names(objNull).")
  }
  
  return(NullModelClass)
}

#' @export
checkOutputFile = function(OutputFile, 
                           OutputFileIndex, 
                           AnalysisType,      ## "Marker" or "Region"
                           nEachChunk)
{
  ## The following messages are for 'OutputFileIndex'
  message1 = "This is the output index file for SPAmixPlus package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)
  message5 = "Have completed the analyses of all chunks."
  
  ## an R list of output
  if(missing(OutputFile))
    stop("Argument of 'OutputFile' is required.")
  
  if(file.exists(OutputFile)){
    if(!file.exists(OutputFileIndex)){
      stop(paste0("'OutputFile' of '", OutputFile,"' has existed. 
                  Please use another 'OutputFile' or remove the existing one."))
    }
    else{
      outIndexData = read.table(OutputFileIndex, header = F, sep="\t")
      
      if(outIndexData[1,1] != message1 | outIndexData[2,1] != message2 | outIndexData[3,1] != message3)
        stop(paste0("'OutputFileIndex' of '", OutputFileIndex, "' is not as expected. 
                    Probably, it has been modified by user, which is not permitted. 
                    Please remove the existing files of 'OutputFile' and 'OutputFileIndex'."))
      
      lastMessage = outIndexData[nrow(outIndexData), 1]
      if(lastMessage == message5){
        End = TRUE
        indexChunk = outIndexData[nrow(outIndexData)-1, 1];
        indexChunk = as.numeric(gsub("Have completed the analysis of chunk ", "", indexChunk))
        cat("Based on 'OutputFile' and 'OutputFileIndex', the analysis has been completed for the toal", indexChunk, "chunks.\n")
      }else{
        End = FALSE;
        indexChunk = lastMessage;
        indexChunk = as.numeric(gsub("Have completed the analysis of chunk ", "", indexChunk))
        cat("Based on 'OutputFile' and 'OutputFileIndex', we restart the analysis from the", indexChunk+1, "chunk.\n")
      }
    }
    Start = FALSE
  }else{
    Start = TRUE;
    End = FALSE;
    indexChunk = 0;
  }
  
  returnList = list(Start = Start, End = End, indexChunk = indexChunk)
  return(returnList)
}

#' @export
writeOutputFile = function(Output,
                           OutputFile,
                           OutputFileIndex,
                           AnalysisType,
                           nEachChunk,
                           indexChunk,
                           Start,   # TRUE or FALSE, to indicate is the 'Output' is the first one to save into 'OutputFile'
                           End)     # TRUE or FALSE, to indicate is the 'Output' is the last one to save into 'OutputFile'
{
  ## The following messages are for 'OutputFileIndex'
  message1 = "This is the output index file for SPAmixPlus package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)
  message4 = paste("Have completed the analysis of chunk", indexChunk)
  message5 = "Have completed the analyses of all chunks."
  
  n1 = length(Output)
  n2 = length(OutputFile)
  
  if(n1 != n2)
    stop("length(Output) != length(OutputFile)")
  
  if(n1 != 0){
    for(i in 1:n1){
      if(Start){
        write.table(Output[[i]], OutputFile[[i]], quote = F, sep = "\t", append = F, col.names = T, row.names = F)
      }else{
        write.table(Output[[i]], OutputFile[[i]], quote = F, sep = "\t", append = T, col.names = F, row.names = F)
      }
    }
  }

  if(Start)
    write.table(c(message1, message2, message3), OutputFileIndex, 
                quote = F, sep = "\t", append = F, col.names = F, row.names = F)
  
  write.table(message4, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
  
  if(End)
    write.table(message5, OutputFileIndex, quote = F, sep = "\t", append = T, col.names = F, row.names = F)
}
