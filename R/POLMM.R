
#' POLMM method in GRAB package
#' 
#' POLMM method is to analyze ordinal categorical data for related samples in a large-scale biobank.
#' 
#' @details 
#' Please check \code{?GRAB.control} for the generic list of \code{control} in \code{GRAB.NullModel()} and \code{GRAB.Marker()}.
#' 
#' Additional list of \code{control} in \code{GRAB.NullModel()} function
#' \itemize{
#' \item{range: a two-element numeric vector [default=c(-100,100)] to specify the domain of the empirical CGF.}
#' \item{length.out: a positive integer [default=10000] for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.}
#' }
#' Additional list of \code{control} in \code{GRAB.Marker()} function
#' \itemize{
#' \item{impute_method: a character string [default="fixed"] to specify the method to impute missing genotypes. "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).}
#' }
#' Additional list of \code{control} in \code{GRAB.Region()} function
#' \itemize{
#' \item{xxx: xxx}
#' }
#' 
#' @examples
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
#' # First use getSparseGRM() function and plink file to get a sparse GRM file
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
#' # Step 2(a): perform marker-level score test
#' GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/POLMMMarkers.txt")
#' GRAB.Marker(obj.POLMM, GenoFile = GenoFile,
#'             OutputFile = OutputFile)
#'             
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#' 
#' # Step 2(b): perform region-level score test
#' objNull = obj.POLMM
#' GenoFile = system.file("extdata", "nSNPs-10000-nsubj-1000-ext.bed", package = "GRAB")
#' GenoFileIndex = NULL
#' 
#' RegionFile = system.file("extdata", "example.RegionFile.txt", package = "GRAB")
#' RegionAnnoHeader = c("ANNO1")
#' 
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/POLMM_Regions.txt")
#' OutputFileIndex = NULL
#' 
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#' control = list(max_maf_region = 0.3)
#' 
#' GRAB.Region(objNull,
#' GenoFile,
#' GenoFileIndex,
#' OutputFile,
#' OutputFileIndex,
#' RegionFile,              # column 1: marker Set ID, column 2: SNP ID, columns 3-n: Annotations similar as in STAAR
#' RegionAnnoHeader,
#' SparseGRMFile,
#' control)
#'            
#' @export
#' @import ordinal
GRAB.POLMM = function(){
  print("Check ?GRAB.POLMM for more details about 'POLMM' method.")
}



################### This file includes the following functions

# ------------ used in 'GRAB_Marker.R' -----------
# 1. checkControl.Marker.POLMM(control)
# 2. setMarker.POLMM(objNull, control)
# 3. mainMarker.POLMM()

# check the control list in marker-level testing
checkControl.Marker.POLMM = function(control)
{
  default.control = list(SPA_Cutoff = 2,
                         outputColumns = c("beta", "seBeta"));
  
  control = updateControl(control, default.control)  # This file is in 'Util.R'
  
  # check the parameter
  
  return(control)
}

# check the control list in region-level testing
checkControl.Region.POLMM = function(control)
{
  default.control = list(SPA_Cutoff = 2,
                         outputColumns = c("beta", "seBeta"));
  
  control = updateControl(control, default.control)  # This file is in 'Util.R'
  
  # check the parameter
  if(!is.numeric(control$SPA_Cutoff) | control$SPA_Cutoff <=0)
    stop("control$SPA_Cutoff should be a numeric value > 0")
  
  return(control)
}

# main function to calculate summary statistics
mainMarker.POLMM = function(genoType, genoIndex, outputColumns)
{
  # Check 'Main.cpp'
  OutList = mainMarkerInCPP("POLMM",
                            genoType,
                            genoIndex);  
  
  obj.mainMarker = data.frame(Marker = OutList$markerVec,           # marker IDs
                              Info = OutList$infoVec,               # marker infomation: CHR:POS:REF:ALT
                              AltFreq = OutList$altFreqVec,         # alternative allele frequencies
                              AltCounts = OutList$altCountsVec,     # alternative allele counts
                              MissingRate = OutList$missingRateVec, # alternative allele counts
                              Pvalue = OutList$pvalVec)             # marker-level p-values
  
  ## optional ouput columns include c("beta", "seBeta", "zScore", "PvalueNorm", "altFreqInGroup", "AltCountsInGroup")
  
  optionalColumns = c("beta", "seBeta", "zScore", "PvalueNorm", "altFreqInGroup", "AltCountsInGroup")
  additionalColumns = intersect(optionalColumns, outputColumns)
  
  if(length(additionalColumns) > 0)
    obj.mainMarker = cbind.data.frame(obj.mainMarker, 
                                      as.data.frame(OutList[additionalColumns]))
  
  return(obj.mainMarker)
}

# main function to calculate summary statistics for region-based test
mainRegion.POLMM = function(genoType, genoIndex)
{
  # add something here
  
}

# check the control list in null model fitting for POLMM method
checkControl.NullModel.POLMM = function(control, optionGRM)
{
  # default setting of control for POLMM method
  default.control = list(memoryChunk = 2,
                         seed = 12345678,
                         tracenrun = 30,
                         maxiter = 100,
                         tolBeta = 0.001,
                         tolTau = 0.002,
                         tau = 0.2,
                         maxiterPCG = 100,
                         tolPCG = 1e-6,
                         maxiterEps = 100,
                         tolEps = 1e-10,
                         minMafVarRatio = 0.1,
                         maxMissingVarRatio = 0.1, 
                         nSNPsVarRatio = 20,
                         CVcutoff = 0.0025,
                         # LOCO = TRUE,
                         numThreads = "auto",
                         stackSize = "auto",
                         grainSize = 1,
                         minMafGRM = 0.01,
                         maxMissingGRM = 0.1,
                         showInfo = T,
                         onlyCheckTime = F)
  
  if(optionGRM == "DenseGRM")
    default.control$LOCO = TRUE;
  
  if(optionGRM == "SparseGRM")
    default.control$LOCO = FALSE;
  
  control = updateControl(control, default.control)
  
  # check the parameters
  # range = control$range
  # length.out = control$length.out
  # 
  # if(range[1] >= -50 | range[2] <= 50 | length.out <= 1000)
  #   stop("We suggest setting argument 'control$range=c(-100,100)' and 'control$length.out=10000'.")
  # 
  # if(range[2]!=-1*range[1])
  #   stop("range[2] should be -1*range[1]")
  
  return(control)
}


# fit null model using POLMM method
fitNullModel.POLMM = function(response, designMat, subjData, control, optionGRM)
{
  ######## -------------- first set up the object in C++ -------- ########
  
  if(class(response) != "factor")
    stop("response in POLMM method should be a factor.")
  
  obj.clm = summary(ordinal::clm(response ~ designMat))
  beta = c(-1 * obj.clm$alpha[1], obj.clm$beta)
  eps = c(0, obj.clm$alpha[-1] - obj.clm$alpha[1])
  bVec = rep(0, length(response))  # initiate random effect of 0
  
  yVec = as.numeric(response) - 1; # "-1" means change from R to C++
  Cova = cbind(1, designMat)
  tau = control$tau
  
  # This value is not used any more, remove it later. 03/31/2021
  SPmatR = list(locations = matrix(c(0,0),2,1),
                values = rep(0,1))
  
  # default.control = list(memoryChunk = 2,
  #                        seed = 12345678,
  #                        tracenrun = 30,
  #                        maxiter = 100,
  #                        tolBeta = 0.001,
  #                        tolTau = 0.002,
  #                        tau = 0.2,
  #                        maxiterPCG = 100,
  #                        tolPCG = 1e-6,
  #                        maxiterEps = 100,
  #                        tolEps = 1e-10,
  #                        minMafVarRatio = 0.1,
  #                        maxMissingVarRatio = 0.1, 
  #                        nSNPsVarRatio = 20,
  #                        CVcutoff = 0.0025,
  #                        LOCO = T,
  #                        numThreads = "auto",
  #                        stackSize = "auto",
  #                        grainSize = 1,
  #                        minMafGRM = 0.01,
  #                        maxMissingGRM = 0.1,
  #                        showInfo = T,
  #                        onlyCheckTime = F)
  
  controlList = control
  flagSparseGRM = ifelse(optionGRM == "SparseGRM", TRUE, FALSE)
  
  # The following function is in 'Main.cpp'
  objNull = setPOLMMobjInCPP_NULL(flagSparseGRM,
                                  Cova,
                                  yVec,
                                  beta,
                                  bVec,
                                  eps,
                                  tau,
                                  SPmatR,
                                  controlList)
  
  class(objNull) = "POLMM_NULL_Model"
  return(objNull)
  
  # void setPOLMMobjInCPP_NULL(bool t_flagSparseGRM,       // if 1, then use SparseGRM, otherwise, use DenseGRM
  #                            arma::mat t_Cova,
  #                            arma::uvec t_yVec,     // should be from 0 to J-1
  #                            arma::vec t_beta,
  #                            arma::vec t_bVec,
  #                            arma::vec t_eps,           // 
  #                              double t_tau,
  #                            arma::mat t_GMatRatio,     // only used if m_LOCO = FALSE
  #                            Rcpp::List t_SPmatR,    // output of makeSPmatR()
  #                            Rcpp::List t_controlList)
  
  ######## -------------- fit the null POLMM --------------  ###########
  
  # bVec = rep(0, n)  # initiate random effect of 0
  
  # objNull = fitPOLMMcpp(t_flagSparseGRM = flagSparseGRM,       # if 1, then use SparseGRM, otherwise, use DenseGRM
  #                       t_flagGMatRatio = flagGMatRatio,       # if 1, then use GMatRatio, otherwise, extract from Plink files
  #                       t_bimfile = bimFile,
  #                       t_famfile = famFile,
  #                       t_bedfile = bedFile,
  #                       t_posSampleInPlink = posSampleInPlink,
  #                       t_Cova = Cova,
  #                       t_yVec = yVec,                         # should be from 1 to J
  #                       t_beta = beta,
  #                       t_bVec = bVec,
  #                       t_eps = eps,
  #                       t_tau = control$tau,
  #                       t_GMatRatio = GMat,                    # only used if m_LOCO = FALSE
  #                       t_SparseGRM = SparseGRM,
  #                       t_controlList = control)
  # return(objNull)
}

setMarker.POLMM = function(objNull, control, chrom)
{
  if(objNull$control$LOCO){
    if(!chrom %in% names(objNull$LOCOList))
      stop("If control$LOCO == TRUE, then 'chrom' should be in names(objNull$LOCOList).")
    objCHR = objNull$LOCOList[[chrom]]
  }else{
    objCHR = objNull$LOCOList[["LOCO=F"]]
  }
  
  # marker-level analysis does not require the following parameters 
  # Note: it might be not so accurate if min_mac_marker is very low
  flagSparseGRM = FALSE;
  # SPmatR.CHR = list(locations = matrix(c(0,0), 2, 1), values = 1)
  printPCGInfo = FALSE
  tolPCG = 0.001
  maxiterPCG = 100;
  
  # Check 'Main.cpp'
  setPOLMMobjInCPP(objCHR$muMat,
                   objCHR$iRMat,
                   objNull$Cova,
                   objNull$yVec,          # 0 to J-1
                   # SPmatR.CHR,
                   objNull$tau,
                   printPCGInfo,
                   tolPCG,
                   maxiterPCG,
                   objCHR$VarRatio, 
                   control$SPA_Cutoff,
                   flagSparseGRM)
  
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,"."))
}

# Used in setRegion() function in GRAB_Region.R
setRegion.POLMM = function(objNull, control, chrom, SparseGRMFile)
{
  if(objNull$control$LOCO){
    if(!chrom %in% names(objNull$LOCOList))
      stop("If control$LOCO == TRUE, then 'chrom' should be in names(objNull$LOCOList).")
    objCHR = objNull$LOCOList[[chrom]]
  }else{
    objCHR = objNull$LOCOList[["LOCO=F"]]
  }
  
  # Since region-level analysis mainly focuses on rare variants, we use sparse GRM for all markers
  
  print("Sparse GRM is used for POLMM-GENE method.")
  SparseGRM = data.table::fread(SparseGRMFile)
  SparseGRM = as.data.frame(SparseGRM)
  KinMatListR = updateSparseGRM(SparseGRM, objNull$subjData)
  setSparseGRMInCPP(KinMatListR)    # check Main.cpp
  
  # The following parameters are not used any more
  flagSparseGRM = TRUE;
  printPCGInfo = FALSE
  tolPCG = 0.001
  maxiterPCG = 100;
  VarRatio = 1
  
  # Check 'Main.cpp'
  setPOLMMobjInCPP(objCHR$muMat,
                   objCHR$iRMat,
                   objNull$Cova,
                   objNull$yVec,          # 0 to J-1
                   objNull$tau,
                   printPCGInfo,
                   tolPCG,
                   maxiterPCG,
                   VarRatio, 
                   control$SPA_Cutoff,
                   flagSparseGRM)

  # print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,"."))
}

