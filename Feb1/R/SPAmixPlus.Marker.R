
#' Conduct marker-level genetic association testing using SPAmixPlus
#' 
#' Test for association between phenotype of interest and genetic marker using the SPAmixPlus method.
#' 
#' @param objNull the output object of function \code{\link{SPAmixPlus.NullModel}}. 
#' @param GenoFile a character of genotype file. Currently, two types of genotype formats are supported: PLINK and BGEN. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param GenoFileIndex additional index files corresponding to the \code{GenoFile}. If \code{NULL} (default), the prefix is the same as GenoFile. Check \code{\link{GRAB.ReadGeno}} for more details.
#' @param OutputFile a character of output file to save the analysis results. 
#' @param OutputFileIndex a character of output index file to record the end point. If the program ends unexpectedly, the end point can help \code{SPAmixPlus} package understand where to restart the analysis. If \code{NULL} (default), \code{OutputFileIndex = paste0(OutputFile, ".index")}. 
#' @param control a list of parameters for controlling function \code{SPAmixPlus.Marker}, more details can be seen in \code{Details} section. 
#' @details 
#' \code{SPAmixPlus} package supports analysis of complex triats in admixed populations.
#' Detailed information about the analysis methods is given in the \code{Details} section of \code{\link{SPAmixPlus.NullModel}}. 
#' 
#' ## The following details are about argument \code{control}
#' The below is to let users customize markers to include in analysis. 
#' If these parameters are not specified, \code{SPAmixPlus} package will include all markers in analysis. 
#' For PLINK files, the default \code{control$AlleleOrder = "alt-first"}; 
#' for BGEN files, the default \code{control$AlleleOrder = "ref-first"}.
#'   \itemize{
#'   \item \code{IDsToIncludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{IDsToExcludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{RangesToIncludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{RangesToExcludeFile}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{AlleleOrder}: please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   }
#' The below is to customize the quality-control (QC) process.
#'   \itemize{
#'   \item \code{omp_num_threads}: a numeric value (default: value from data.table::getDTthreads()) to specify the number of threads in OpenMP for parallel computation.
#'   \item \code{ImputeMethod}: a character, "mean" (default), "bestguess", or "drop". Please refer to the \code{Details} section of \code{\link{GRAB.ReadGeno}}.
#'   \item \code{MissingRateCutoff}: a numeric value *(default=0.15)*. Markers with missing rate > this value will be excluded from analysis.  
#'   \item \code{MinMAFCutoff}: a numeric value *(default=0.001)*. Markers with MAF < this value will be excluded from analysis.  
#'   \item \code{MinMACCutoff}: a numeric value *(default=20)*. Markers with MAC < this value will be excluded from analysis.  
#'   \item \code{nMarkersEachChunk}: number of markers *(default=10000)* in one chunk to output.
#'   }
#'  The below is to customize the columns in the \code{OutputFile}. 
#'  Columns of \code{Marker}, \code{Info}, \code{AltFreq}, \code{AltCounts}, \code{MissingRate}, \code{Pvalue} are included for all methods.
#'  \itemize{
#'  \item \code{outputColumns}: Users can set \code{control$outputColumns} to include additional columns.
#'     \itemize{
#'     \item Default columns: \code{Marker}, \code{Info}, \code{AltFreq}, \code{AltCounts}, \code{MissingRate}, \code{Pvalue}
#'     \item Optional: \code{beta}, \code{seBeta}, \code{zScore}, \code{PvalueNorm}, \code{AltFreqInGroup}, \code{nSamplesInGroup}, \code{AltCountsInGroup}
#'     }
#'  }
#' @return The analysis results are written in a file of \code{OutputFile}, which includes the following columns.
#' \item{Marker}{Marker IDs extracted from \code{GenoFile} and \code{GenoFileIndex}.}
#' \item{Info}{Marker Information of "CHR:POS:REF:ALT". The order of REF/ALT depends on \code{control$AlleleOrder}: "ref-first" or "alt-first".}
#' \item{AltFreq}{Alternative allele frequency (before genotype imputation, might be > 0.5). If the \code{AltFreq} of most markers are > 0.5, you should consider resetting \code{control$AlleleOrder}.}
#' \item{AltCounts}{Alternative allele counts (before genotype imputation).}
#' \item{MissingRate}{Missing rate for each marker}
#' \item{Pvalue}{Association test p-value}
#' The following columns can be customized using \code{control$outputColumns}. Check \code{\link{makeGroup}} for details about phenotype grouping which are used for
#' \code{nSamplesInGroup}, \code{AltCountsInGroup}, and \code{AltFreqInGroup}.
#' \item{beta}{Estimated effect size of the ALT allele.}
#' \item{seBeta}{Estimated standard error (se) of the effect size.}
#' \item{zScore}{z score, standardized score statistics, usually follows a standard normal distribution.}
#' \item{nSamplesInGroup}{Number of samples in different phenotype groups. This can be slightly different from the original distribution due to the genotype missing.}
#' \item{AltCountsInGroup}{Alternative allele counts (before genotype imputation) in different phenotype groups.}
#' \item{AltFreqInGroup}{Alternative allele frequency (before genotype imputation) in different phenotype groups.}
#' @examples
#' \dontrun{
#' # Assuming you have run SPAmixPlus.NullModel and have an objNull
#' # objNullFile = system.file("results", "objNull.RData", package = "SPAmixPlus")
#' # load(objNullFile)
#' 
#' # OutputDir = system.file("results", package = "SPAmixPlus")
#' # OutputFile = paste0(OutputDir, "/simuOUTPUT.txt")
#' # GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' 
#' # SPAmixPlus.Marker(objNull,
#' #             GenoFile = GenoFile,
#' #             OutputFile = OutputFile)
#'             
#' # data.table::fread(OutputFile)
#' }
#' 
#' @export
#' @import data.table
SPAmixPlus.Marker = function(objNull,
                       GenoFile,
                       GenoFileIndex = NULL,
                       OutputFile,
                       OutputFileIndex = NULL,   
                       control = NULL)
{
  NullModelClass = checkObjNull(objNull);         # this function is in 'Framework.R'
  
  cat(paste0("NullModelClass is ", NullModelClass, "\n")) 
  
  if(is.null(OutputFileIndex)) 
    OutputFileIndex = paste0(OutputFile, ".index")
  
  # check the setting of control
  checkControl.ReadGeno(control)
  control = checkControl.Marker(control, NullModelClass)  
  nMarkersEachChunk = control$nMarkersEachChunk;
  outList = checkOutputFile(OutputFile, OutputFileIndex, "Marker", format(nMarkersEachChunk, scientific=F))    # this function is in 'Framework.R'
  
  # method = "SPAmixPlus" derived from "SPAmixPlus_NULL_Model"
  method = gsub("_NULL_Model", "", NullModelClass)

  indexChunk = outList$indexChunk
  Start = outList$Start
  End = outList$End
  
  if(End)
  {
    message = paste0("The analysis has been completed in earlier analysis. Results have been saved in '", OutputFile, "'. ",
                     "If you want to change parameters and restart the analysis, please use another 'OutputFile'.")
    return(message)
  }
  
  if(!Start){
    message = paste0("We detected that parts of analysis have been conducted from file:\t",
                     OutputFileIndex,"\n",
                     "We restart the analysis from chunk:\t",indexChunk+1,"\n");
    cat(message)
  }
  
  subjData = as.character(objNull$subjData);
  
  Group = makeGroup(objNull$yVec)  # this function is in Framework.R
  ifOutGroup = any(c("AltFreqInGroup", "AltCountsInGroup") %in% control$outputColumns)
  
  # Set global variables in C++ for filtering (MAF, MAC, etc.) and grouping
  # Implemented based on GRAB package logic to ensure consistency
  setMarker_GlobalVarsInCPP(control$impute_method,
                            control$missing_cutoff,
                            control$min_maf_marker,
                            control$min_mac_marker,
                            control$omp_num_threads,
                            Group,
                            ifOutGroup,
                            length(unique(Group)))
  
  ## set up an object for genotype
  objGeno = setGenoInput(GenoFile, GenoFileIndex, SampleIDs = subjData, control = control)  # this function is in 'Geno.R'
  genoType = objGeno$genoType
  markerInfo = objGeno$markerInfo
  CHROM = markerInfo$CHROM
  genoIndex = markerInfo$genoIndex
  
  # all markers were split into multiple chunks
  genoIndexList = splitMarker(genoIndex, nMarkersEachChunk, CHROM);  
  
  nChunks = length(genoIndexList)
  
  cat("Number of all markers to test:\t", nrow(markerInfo), "\n")
  cat("Number of markers in each chunk:\t", nMarkersEachChunk, "\n")
  cat("Number of chunks for all markers:\t", nChunks, "\n")
  
  chrom = "InitialChunk"
  for(i in (indexChunk+1):nChunks)
  {
    tempList = genoIndexList[[i]]
    genoIndex = tempList$genoIndex
    tempChrom = tempList$chrom
    
    # Update Null Model Object (C++) if necessary. 
    # For SPAmixPlus, we call setMarker.SPAmixPlus every time, or at least when we start.
    # The original logic calls it every chunk if chrom changes, or maybe every chunk?
    # In GRAB.Marker:
    # if(chrom != tempChrom){ ... setMarker... }
    # Let's just call it once at the beginning if chrom changes, or always?
    # SPAmixPlus doesn't depend on Chromosome. But calling it is cheap (just setting pointer).
    # But wait, setSPAmixPlusobjInCPP involves passing large vectors? 
    # No, it passes them. Rcpp might copy them.
    # So we should only call it once at the start.
    
    if(i == (indexChunk+1)){
        setMarker.SPAmixPlus(objNull, control)
    }
    
    cat(paste0("Analysis of chunk ", i, "/", nChunks, " starts (Chrom: ", tempChrom, ")...\n"))
    
    # mainMarker function
    # call mainMarker.SPAmixPlus
    obj.mainMarker = mainMarker.SPAmixPlus(genoType, genoIndex, control$outputColumns, objNull, control)
    
    # if OutputFile includes other columns customized by users
    if(ifOutGroup){
      # To be implemented if needed. For now assume basic columns.
      # The original GRAB.Marker implemented this.
      # But mainMarker.SPAmixPlus returns what it returns.
    }
    
    # write the output
    writeOutputFile(list(obj.mainMarker), OutputFile, OutputFileIndex, "Marker", nMarkersEachChunk, i, i==1, i==nChunks)
    
    cat(paste0("Analysis of chunk ", i, "/", nChunks, " completed.\n"))
  }
  
  closeGenoInputInCPP(genoType)
  
  message = paste0("The analysis has been completed. Results have been saved in '", OutputFile, "'.")
  return(message)
}
