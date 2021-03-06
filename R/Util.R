

checkObjNull = function(objNull)
{
  NullModelClass = class(objNull)
  nm = names(objNull)
  
  ## check objNull
  if(!NullModelClass %in% c("SPAGE_NULL_Model",       # SPAGE: GxE analysis 
                            "SPACox_NULL_Model",      # SPACox: Survival analysis for unrelated subjects
                            "POLMM_NULL_Model"))      # POLMM: categorical data analysis
  {
    stop('class(objNull) should be from c("SPAGE_NULL_Model", "SPACox_NULL_Model", "POLMM_NULL_Model")')
  }
    
  if(any(nm %in% c("SampleIDs","n")))
  {
    stop("c('SampleIDs', 'n') should be in names(objNull).")
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
  if(missing(OutputFile))
    stop("Argument of 'OutputFile' is required.")
  
  if(file.exists(OutputFile)){
    if(!file.exists(OutputFileIndex)){
      stop(paste0("'OutputFile' of '", OutputFile, "' has existed. Please use another 'OutputFile' or remove the existing one."))
    }
    else{
      outIndexData = read.table(OutputFileIndex)
      
      if(outIndexData[1,1] != "GRAB.outIndex" | outIndexData[2,1] != "Please_do_not_modify_this_file." | outIndexData[3,1] != AnalysisType | outIndexData[4,1] != nEachChunk)
        stop("'OutputFileIndex' is not expected. Please use another 'OutputFileIndex' or remove the existing one.")
      
      outIndex = as.numeric(outIndexData[nrow(outIndexData),1]) + 1
      if(outIndex == 0)
        stop("Based on 'OutputFile' and 'OutputFileIndex', the analysis has been completed.")
      
      warning(paste0("Based on 'OutputFile' and 'OutputFileIndex', we restart the analysis from the ", outIndex, "-th chunks."))
    }
  }else{
    outIndex = 1;
  }
  
  return(outIndex)
}


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
