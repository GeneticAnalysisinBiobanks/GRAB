#' A lower function to make groups based on phenotype
#' 
#' In functions \code{\link{GRAB.Marker}} and \code{\link{GRAB.Region}}, users can get detailed information for each markers in different groups.
#' 
#' @param yVec the phenotype recorded in \code{objNull$yVec}, the output object of function \code{\link{GRAB.NullModel}}. 
#' @return a numeric vector (\code{Group}, starting from 0) for group information.
#' @details 
#' If \code{yVec} is categorical with groups <= 10, then \code{Group} is the same as \code{yVec}. Otherwise, \code{Group} is calcualted based on the rank of \code{yVec}.
 
makeGroup = function(yVec)
{
  # yVec is categorical data
  m1 = length(unique(yVec))
  if(m1 <= 10)
    Group = as.numeric(as.factor(yVec))-1  # from 0 to (m1-1)
  
  # yVec is quantitative data
  if(length(unique(yVec)) > 10)
    Group = floor((rank(yVec, ties.method = "max")-1) / length(yVec) * 10)  # from 0 to 9
  
  # yVec is a time-to-event data (tbc)
  
  return(Group)
}

checkObjNull = function(objNull)
{
  NullModelClass = class(objNull)
  nm = names(objNull)
  
  ## check objNull
  if(!NullModelClass %in% c("SPAGE_NULL_Model",       # SPAGE: GxE analysis 
                            "SPACox_NULL_Model",      # SPACox: Survival analysis for unrelated subjects
                            "POLMM_NULL_Model",       # POLMM: categorical data analysis
                            "SPAmix_NULL_Model",      # SPAmix: mixture population analysis
                            "SPAGRM_NULL_Model",      # SPAGRM: related subjects
                            
                            "SPAyuzhuoma_NULL_Model", # SPAyuzhuoma: related subjects
                            
                            "SAGELD_NULL_Model",      # SAGELD: GxE for longitudinal data
                            "WtSPAG_NULL_Model"))      
  {
    stop('class(objNull) should be one of 
         c("SPAGE_NULL_Model", "SPACox_NULL_Model", "POLMM_NULL_Model", "SPAmix_NULL_Model", "SPAGRM_NULL_Model", "SPAyuzhuoma_NULL_Model", "SAGELD_NULL_Model", "WtSPAG_NULL_Model")')
  }
    
  if(any(!c("subjData", "N") %in% nm))
  {
    stop("c('subjData', 'N') should be in names(objNull).")
  }
  
  ## More detailed information can be checked
  ##
  ##
  ##
    
  return(NullModelClass)
}

checkOutputFile = function(OutputFile, 
                           OutputFileIndex, 
                           AnalysisType,      ## "Marker" or "Region"
                           nEachChunk)
{
  ## The following messages are for 'OutputFileIndex'
  message1 = "This is the output index file for GRAB package to record the end point in case users want to restart the analysis. Please do not modify this file."
  message2 = paste("This is a", AnalysisType, "level analysis.")
  message3 = paste("nEachChunk =", nEachChunk)
  # message4 = paste("Have completed the analysis of chunk", indexChunk)
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
  message1 = "This is the output index file for GRAB package to record the end point in case users want to restart the analysis. Please do not modify this file."
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

# This is the output index file for GRAB package. This file is to record the end point in case users want to restart the analysis. Please do not modify this file.
# This is a region/marker level analysis.
# The analysis has been completed.

checkControl = function(control = NULL)  
{
  default.control = list(impute_method = "fixed",  # the below are shared parameters for both single marker testing and region-based testing
                         missing_cutoff = 0.15,
                         min_maf_marker = 0.001,
                         min_mac_marker = 20,
                         max_maf_region = 0.01,
                         nMarkersEachChunk = 10000,
                         memory_chunk = 4,
                         kernel = "linear.weighted",   # the below are parameters for region-based testing
                         method_region = "SKAT-O",
                         weights_beta = c(1, 25),
                         r_corr = NULL,
                         printPCGInfo = FALSE,
                         tolPCG = 1e-5,
                         maxiterPCG = 100,
                         SPA_cutoff = 2)
  
  # use the default setting or update it
  if(!is.null(control)){
    if(class(control) != "list"){
      stop("Argument of 'control' should be 'list'.")
    }
    ctrl.nm = names(control)
    for(nm in ctrl.nm){
      default.control[[nm]] = control[[nm]]
    }
  }
      
  control = default.control
  
  if(! control$impute_method %in% c("fixed", "bestguess", "random"))
    stop("'impute.method' should be 'fixed','bestguess', or 'random'. Check 'Details' for more details.")
  
  if(control$missing_cutoff > 1 | control$missing_cutoff < 0)
    stop("'missing_cutoff' should be between 0 and 1. The default setting is 0.15.")
  
  if(control$SPA_cutoff < 0)
    stop("'SPA_cutoff' should be greater than or equal to 0. The default setting is 2. Check 'Details' for more details.")
  
  # only used in single-marker analysis
  if(control$min_maf_marker > 1 | control$min_maf_marker <= 0)
    stop("'min_maf_marker' should be between 0 and 1. The default setting is 0.001.")
  
  if(control$min_mac_marker <= 5)
    stop("'min_mac_marker' should be greater than 5. The default setting is 20.")
  
  if(control$nMarkers_output <= 999)
    stop("nMarkers_output should be greater than or equal to 1000. The default setting is 10000.")
  
  # only used in region-based analysis
  if(control$max_maf_region > 1 | control$max_maf_region <= 0)
    stop("'max_maf_region' should be between 0 and 1. The default setting is 0.01.")
  
  if(! control$kernel %in% c("linear", "linear.weighted"))
    stop("'kernel' should be 'linear' or 'linear.weighted'. Check 'Details' for more details.")
  
  if(length(control$weights_beta) != 2 | any(control$weights_beta < 0))
    stop("length of 'weights_beta' should be 2. The two elements in 'weights_beta' should be non-negative. Check 'Details' for more details.")
  
  method_region = control$method_region
  
  if(! method_region %in% c("SKAT", "Burden", "SKAT-O"))
    stop("'method' should be 'SKAT', 'Burden', or 'SKAT-O'. Check 'Details' for more details.")
  
  if(is.null(control$r_corr)){
    if(method_region == "SKAT") control$r_corr = 0;
    if(method_region == "Burden") control$r_corr = 1;
    if(method_region == "SKAT-O") control$r_corr = c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);  # r_corr = 0 is SKAT, r_corr = 1 is Burden Test
  }else{
    message("Since 'r_corr' is specified, the 'method_region' is ignored.")
  }
  
  if(any(control$r_corr < 0 | control$r_corr > 1))
    stop("'r_corr' should be a numeric vector in which each element is between 0 and 1. Check 'Details' for more details.")
  
  return(control)
}


getAnnoList = function(AnnoFile, AnnoHeader)
{
  if(!file.exists(AnnoFile))
    stop(paste("Cannot find AnnoFile in", AnnoFile))
  
  AnnoData = data.table::fread(AnnoFile, header = T, stringsAsFactors = F);
  AnnoData = as.data.frame(AnnoData)
  HeaderInAnnoFile = colnames(AnnoData)
  
  if(any(HeaderInAnnoFile[1:2] != c("GENE","SNP")))
    stop("The first two elements in the header of AnnoFile should be c('GENE', 'SNP').")
  
  if(!is.null(AnnoHeader)){
    if(any(!AnnoHeader %in% HeaderInAnnoFile))
      stop("At least one element in AnnoHeader is not in the header of AnnoFile.")
    posAnno = which(HeaderInAnnoFile %in% AnnoHeader)
  }else{
    posAnno = NULL
  }
  
  AnnoList = list()
  uGENE = unique(AnnoData$GENE)
  for(g in uGENE){
    pos = which(AnnoData$GENE == g)
    SNP = AnnoData$SNP[pos]
    AnnoMat = cbind(All=1, AnnoData[pos, posAnno, drop=F])
    rownames(AnnoMat) = SNP
    if(any(duplicated(SNP)))
      stop(paste0("Please check AnnoFile: in gene ",g,", duplicated SNPs exist."))
    AnnoList[[g]] = list(SNP = SNP,
                         AnnoMat = AnnoMat)
  }
  
  return(AnnoList)
}

# make a matrix that can be passed to arma::sp_mat
makeSPmatR = function(SparseGRM,   # three columns of ID1, ID2, and value
                      subjData)
{
  SparseGRM = subset(SparseGRM, ID1 %in% subjData & ID2 %in% subjData)
  
  row.diag = which(SparseGRM$ID1 == SparseGRM$ID2)
  SparseGRM.diag = SparseGRM[row.diag,]
  SparseGRM.off.d1 = SparseGRM[-1*row.diag,]
  SparseGRM.off.d2 = data.frame(ID1 = SparseGRM.off.d1$ID2,
                                ID2 = SparseGRM.off.d1$ID1,
                                value = SparseGRM.off.d1$value)
  
  SparseGRM = rbind(SparseGRM.diag, SparseGRM.off.d1, SparseGRM.off.d2)
  
  ID1 = SparseGRM$ID1;
  ID2 = SparseGRM$ID2;
  value = SparseGRM$value;
  
  if(any(!is.element(subjData, ID1)))
    stop("At least one of subjects in `subjData` is not in `SparseGRM`.")
  
  location1 = match(ID1, subjData);
  location2 = match(ID2, subjData);
  
  locations = rbind(location1 - 1,  # -1 is to convert R to C++
                    location2 - 1)
  
  SPmatR = list(locations = locations,
                values = value)
  
  return(SPmatR)
}

getMaxMarkers = function(memory_chunk,
                         n, J, p)
{
  # THE BELOW include some large matrix that takes most of the memory usage
  
  # VarSMat = mat(indexPassingQC, indexPassingQC)
  # adjGMat = fmat(n, maxMarkers): 4*n*maxMarkers
  # ZPZ_adjGMat = fmat(n, maxMarkers): 4*n*maxMarkers
  # CovaMat = mat(n(J-1) x p): 8*n*(J-1)*p
  # arma::mat m_XXR_Psi_RX;  // XXR_Psi_RX ( n x p )
  # arma::mat m_XR_Psi_R;    // XR_Psi_R ( p x n ), sum up XR_Psi_R ( p x n(J-1) ) for each subject 
  # arma::vec m_RymuVec;     // n x 1: row sum of the n x (J-1) matrix R %*% (yMat - muMat)
  # arma::mat m_iSigmaX_XSigmaX; // n(J-1) x p
  # arma::mat m_CovaMat;     // n(J-1) x p
  
  fixed.memory = 8*(3*(n*(J-1)*p)+2*(n*p)) / 1e9
  if(memory_chunk < fixed.memory)
    stop(paste0("Please give control$memory_chunk greater than ", fixed.memory,"."))
  
  maxMarkers = (memory_chunk - fixed.memory) * 1e9 / (8*n)
  maxMarkers = maxMarkers / 2;
  maxMarkers = floor(maxMarkers)
  
  return(maxMarkers)
}

####### ---------- Get Weights from MAF ---------- #######

Get_Weights = function(kernel, freqVec, weights_beta)
{
  if(kernel == "linear"){
    weights = rep(1, length(freqVec)) 
  }
  
  if(kernel == "linear.weighted"){
    weights = dbeta(freqVec, weights_beta[1], weights_beta[2])
  }
  
  return(weights)
}
