
#' Conduct marker-level genetic association testing
#' 
#' Test for association between phenotype of interest and genetic marker.
#' 
#' @param objNull the output object of function \code{\link{GRAB.NullModel}}. 
#' @param GenoFile a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param GenoFileIndex additional index files corresponding to the \code{GenoFile}. If \code{NULL} (default), the prefix is the same as GenoFile. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param OutputFile a character of output file to save the analysis results. 
#' @param OutputFileIndex a character of output index file to record the end point. If the program ends unexpectedly, the end point can help \code{GRAB} package understand where to restart the analysis. If \code{NULL} (default), \code{OutputFileIndex = paste0(OutputFile, ".index")}. 
#' @param control a list of parameters for controlling function \code{GRAB.Marker}, more details can be seen in \code{Details} section. 
#' @details 
#' \code{GRAB} package supports \code{SAIGE}, \code{POLMM}, and \code{SPACox} methods. 
#' Detailed information about the analysis methods is given in the \code{Details} section of \code{\link{GRAB.NullModel}}. 
#' Users do not need to specify them since functions \code{GRAB.Marker} and \code{\link{GRAB.Region}} will check the \code{class(objNull)}.
#' 
#' ## The following details are about argument \code{control}
#' The below is to let users customize markers to include in analysis. If these parameters are not specified, \code{GRAB} package will include all markers in analysis. 
#' For PLINK files, the default \code{control$AlleleOrder = "alt-first"}; for BGEN files, the default \code{control$AlleleOrder = "alt-first"}.
#'   \itemize{
#'   \item \code{IDsToIncludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{IDsToExcludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{RangesToIncludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{RangesToExcludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{AlleleOrder}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   }
#' The below is to customize the quality-control (QC) process.
#'   \itemize{
#'   \item \code{omp_num_threads}: (To be added later) a numeric value (default: value from data.table::getDTthreads()) to specify the number of threads in OpenMP for parallel computation.
#'   \item \code{ImputeMethod}: a character, "mean", "bestguess", or "drop" (to be added later). please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{MissingRateCutoff}: a numeric value *(default=0.15)*. Markers with missing rate > this value will be excluded from analysis.  
#'   \item \code{MinMAFCutoff}: a numeric value *(default=0.001)*. Markers with MAF < this value will be excluded from analysis.  
#'   \item \code{MinMACCutoff}: a numeric value *(default=20)*. Markers with MAC < this value will be excluded from analysis.  
#'   \item \code{nMarkersEachChunk}: number of markers *(default=10000)* in one chunk to output.
#'   }
#'  The below is to customize the columns in the \code{OutputFile}. 
#'  Columns of \code{Marker}, \code{Info}, \code{AltFreq}, \code{AltCounts}, \code{MissingRate}, \code{Pvalue} are included for methods.
#'  \itemize{
#'  \item \code{outputColumns}: For example, for POLMM method, users can set \code{control$outputColumns = c("beta", "seBeta", "AltFreqInGroup")}. 
#'     \itemize{
#'     \item \code{POLMM}: Default: \code{beta}, \code{seBeta}; Optional: \code{PvalueNorm}, \code{zScore}, \code{AltFreqInGroup}
#'     \item \code{SPACox}: Optional: \code{PvalueNorm}, \code{zScore}
#'     }
#'  }
#' @return The results is written in a file of \code{OutputFile}, which includes the following columns.
#' \item{Marker}{Marker IDs extracted from \code{GenoFile} and \code{GenoFileIndex}.}
#' \item{Info}{Marker Information of "CHR:POS:REF:ALT". The order of REF/ALT depends on \code{control$AlleleOrder}: "ref-first" or "alt-first".}
#' \item{AltFreq}{Alternative allele frequency (before genotype imputation, might be > 0.5). If the \code{AltFreq} of most markers are > 0.5, you should consider resetting \code{control$AlleleOrder}.}
#' \item{AltCounts}{Alternative allele counts (before genotype imputation).}
#' \item{MissingRate}{Missing rate for each marker}
#' \item{Pvalue}{Association test p-value}
#' The following columns can be customized in \code{control$outputColumns}.
#' \item{Beta}{Estimated effect size of the ALT allele.}
#' \item{seBeta}{Estimated standard error (se) of the effect size.}
#' \item{zScore}{z score, standardized score statistics, usually follows a standard normal distribution.}
#' \item{altFreqInGroup}{Alternative frequency (before genotype imputation) in different phenotype groups.}
#' @examples
#' # We put examples to the specific help pages for different methods. 
#' # If you want to use "SPACox" method, please check ?GRAB.SPACox for more details.
#' # If you want to use "POLMM" method, please check ?GRAB.POLMM for more details.
#' @export
#' @import data.table

GRAB.Marker = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile,
                       OutputFileIndex = NULL,   
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);         # this function is in 'Util.R'
  
  if(is.null(OutputFileIndex)) 
    OutputFileIndex = paste0(OutputFile, ".index")
  
  # check the setting of control, if not specified, the default setting will be used
  # The following functions are in 'control.R'
  checkControl.ReadGeno(control)
  control = checkControl.Marker(control, NullModelClass)
  nMarkersEachChunk = control$nMarkersEachChunk;
  outIndex = checkOutputFile(OutputFile, OutputFileIndex, "Marker", format(nMarkersEachChunk, scientific=F))    # this function is in 'Util.R'
  
  subjData = as.character(objNull$subjData);
  
  ## set up an object for genotype
  objGeno = setGenoInput(GenoFile, GenoFileIndex, subjData, control)  # this function is in 'Geno.R'
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  CHROM = markerInfo$CHROM
  genoIndex = markerInfo$genoIndex
  
  # all markers were split into multiple chunks, 
  genoIndexList = splitMarker(genoIndex, nMarkersEachChunk, CHROM);
  nChunks = length(genoIndexList)
  
  cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
  cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
  cat("Number of chunks for all markers:\t", nChunks, "\n")
  if(outIndex != 1)
    cat("Restart the analysis from chunk:\t", outIndex, "\n")
  
  chrom = "InitialChunk"
  for(i in outIndex:nChunks)
  {
    tempList = genoIndexList[[i]]
    genoIndex = tempList$genoIndex
    tempChrom = tempList$chrom
    
    # set up objects that do not change for different variants
    if(tempChrom != chrom){
      setMarker(NullModelClass, objNull, control, chrom)
      chrom = tempChrom
    }
    
    print(paste0("(",Sys.time(),") ---- Analyzing Chunk ", i, "/", nChunks, ": chrom ", chrom," ---- "))
    
    # main function to calculate summary statistics for markers in one chunk
    resMarker = mainMarker(NullModelClass, genoType, genoIndex)
    
    # write summary statistics to output file
    if(i == 1){
      data.table::fwrite(resMarker, OutputFile, quote = F, sep = "\t", append = F, col.names = T, na="NA")
      write.table(matrix(c("GRAB.outIndex", "Please_do_not_modify_this_file.", "Marker", 
                           format(nMarkersEachChunk, scientific=F), 1), 
                         ncol = 1), 
                  OutputFileIndex, col.names = F, row.names = F, quote = F, append = F)
    }else{
      data.table::fwrite(resMarker, OutputFile, quote = F, sep = "\t", append = T, col.names = F, na="NA")
      write.table(matrix(i, ncol = 1), 
                  OutputFileIndex, col.names = F, row.names = F, quote = F, append = T)
    }
  }
  
  # information to users
  output = paste0("Done! The results have been saved to '", OutputFile,"'.")
  write.table(matrix(-1, ncol = 1), 
              OutputFileIndex, col.names = F, row.names = F, quote = F, append = T)
  
  return(output)
}

setMarker = function(NullModelClass, objNull, control, chrom)
{
  # The following function is in Main.cpp
  setMarker_GlobalVarsInCPP(control$impute_method,
                            control$missing_cutoff,
                            control$min_maf_marker,
                            control$min_mac_marker,
                            control$omp_num_threads)
  
  # The following function is in POLMM.R
  if(NullModelClass == "POLMM_NULL_Model")
    obj.setMarker = setMarker.POLMM(objNull, control, chrom)
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.setMarker = setMarker.SAIGE(objNull, control)
  
  # The following function is in SPACox.R
  if(NullModelClass == "SPACox_NULL_Model")
    obj.setMarker = setMarker.SPACox(objNull, control)
    
  return(obj.setMarker)
}

mainMarker = function(NullModelClass, genoType, genoIndex)
{
  # The following function is in 'POLMM.R'
  if(NullModelClass == "POLMM_NULL_Model")
    obj.mainMarker = mainMarker.POLMM(genoType, genoIndex)
  
  if(NullModelClass == "SAIGE_NULL_Model")
    obj.mainMarker = mainMarker.SAIGE(genoType, genoIndex)
  
  # The following function is in SPACox.R
  if(NullModelClass == "SPACox_NULL_Model")
    obj.mainMarker = mainMarker.SPACox(genoType, genoIndex)
  
  return(obj.mainMarker)
}


## split 'markerInfo' into multiple chunks, each of which includes no more than 'nMarkersEachChunk' markers
splitMarker = function(genoIndex, nMarkersEachChunk, CHROM)
{
  genoIndexList = list()
  iTot = 1;
  
  uCHROM = unique(CHROM)
  for(chrom in uCHROM){
    pos = which(CHROM == chrom)
    gIdx = genoIndex[pos]
    M = length(gIdx)
    
    idxStart = seq(1, M, nMarkersEachChunk)
    idxEnd = idxStart + nMarkersEachChunk - 1
    
    nChunks = length(idxStart)
    idxEnd[nChunks] = M
    
    for(i in 1:nChunks){
      idxMarker = idxStart[i]:idxEnd[i]
      genoIndexList[[iTot]] = list(chrom = chrom,
                                   genoIndex = gIdx[idxMarker])
      iTot = iTot + 1;
    }
  }
  
  return(genoIndexList)
}