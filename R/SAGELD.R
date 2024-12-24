#' SAGELD method in GRAB package
#' 
#' SAGELD method is Scalable and Accurate algorithm for Gene-Environment interaction analysis using Longitudinal Data for related samples in a large-scale biobank. SAGELD extended SPA<sub>GRM</sub> to support gene-environment interaction analysis.
#' 
#' @details 
#' Additional list of \code{control} in \code{SPAGRM.NullModel()} function.
#' 
#' Additional list of \code{control} in \code{GRAB.Marker()} function.
#' 
#' @examples 
#' # Step 2a: process model residuals
#' ResidMatFile = system.file("extdata", "ResidMat.txt", package = "GRAB")
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#' PairwiseIBDFile = system.file("PairwiseIBD", "PairwiseIBD.txt", package = "GRAB")
#' obj.SPAGRM = SPAGRM.NullModel(ResidMatFile = ResidMatFile, 
#'                               SparseGRMFile = SparseGRMFile, 
#'                               PairwiseIBDFile = PairwiseIBDFile,
#'                               control = list(ControlOutlier = FALSE))
#' 
#' # Step 2b: perform score test
#' GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/SPAGRMMarkers.txt")
#' GRAB.Marker(objNull = obj.SPAGRM,
#'             GenoFile = GenoFile,
#'             OutputFile = OutputFile)
#' head(read.table(OutputFile, header=T))
#' @export
GRAB.SAGELD = function(){
  print("Check ?GRAB.SAGELD for more details about 'SAGELD' method.")
}

################### This file includes the following functions

# ------------ used in 'GRAB_Marker.R' -----------
# 1. checkControl.Marker.SAGELD(control)
# 2. setMarker.SAGELD(objNull, control)
# 3. mainMarker.SAGELD()

# check the control list in marker-level testing
checkControl.Marker.SAGELD = function(control)
{
  default.control = list(SPA_Cutoff = 2,
                         zeta = 0,
                         tol = 1e-4)
  
  control = updateControl(control, default.control)  # This file is in 'control.R'
  
  return(control)
}

checkControl.SAGELD.NullModel = function(control,
                                         ResidMat,
                                         SparseGRM,
                                         PairwiseIBD)
{
  default.control = list(MaxQuantile = 0.75,
                         MinQuantile = 0.25,
                         OutlierRatio = 1.5,
                         ControlOutlier = TRUE,
                         MaxNuminFam = 5,
                         MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5))
  
  control = updateControl(control, default.control)  # This file is in 'control.R'
  
  if(control$MaxQuantile < control$MinQuantile)
    stop("MaxQuantile(default is 0.75) should be larger than MinQuantile(default is 0.25).")
  
  if(control$OutlierRatio < 0)
    stop("OutlierRatio should be larger than or equal 0 (default is 1.5).")
  
  if(any(colnames(SparseGRM) != c("ID1", "ID2", "Value")))
    stop("The column names of SparseGRM should be ['ID1', 'ID2', 'Value'].")
  
  if(any(colnames(PairwiseIBD) != c("ID1", "ID2", "pa", "pb", "pc")))
    stop("The column names of PairwiseIBD should be ['ID1', 'ID2', 'pa', 'pb', 'pc'].")
  
  SubjID.In.Resid = ResidMat$SubjID
  SubjID.In.GRM = unique(c(SparseGRM$ID1, SparseGRM$ID2))
  SubjID.In.IBD = unique(c(PairwiseIBD$ID1, PairwiseIBD$ID2))
  
  if(any(!SubjID.In.Resid %in% SubjID.In.GRM))
    stop("At least one subject in residual matrix does not have GRM information.")
  
  if(any(!SubjID.In.IBD %in% SubjID.In.GRM))
    stop("At least one subject has IBD information but does not have GRM information.")
  
  return(control)
}

setMarker.SAGELD = function(objNull, control)
{
  # the below function is in 'Main.cpp'
  setSAGELDobjInCPP(objNull$Method, 
                    objNull$XTs, 
                    objNull$SS, 
                    objNull$AtS, 
                    objNull$Q, 
                    objNull$A21, 
                    objNull$TTs, 
                    objNull$Tys, 
                    objNull$sol, 
                    objNull$blups, 
                    objNull$sig,
                    objNull$Resid,
                    objNull$Resid_G,
                    objNull$Resid_GxE,
                    objNull$Resid_E,
                    objNull$Resid.unrelated.outliers,
                    objNull$Resid.unrelated.outliers_G,
                    objNull$Resid.unrelated.outliers_GxE,
                    objNull$sum_R_nonOutlier,
                    objNull$sum_R_nonOutlier_G,
                    objNull$sum_R_nonOutlier_GxE,
                    objNull$R_GRM_R,
                    objNull$R_GRM_R_G,
                    objNull$R_GRM_R_GxE,
                    objNull$R_GRM_R_G_GxE,
                    objNull$R_GRM_R_E,
                    objNull$R_GRM_R_nonOutlier,
                    objNull$R_GRM_R_nonOutlier_G,
                    objNull$R_GRM_R_nonOutlier_GxE,
                    objNull$R_GRM_R_nonOutlier_G_GxE,
                    objNull$R_GRM_R_TwoSubjOutlier,
                    objNull$R_GRM_R_TwoSubjOutlier_G,
                    objNull$R_GRM_R_TwoSubjOutlier_GxE,
                    objNull$R_GRM_R_TwoSubjOutlier_G_GxE,
                    objNull$TwoSubj_list,
                    objNull$ThreeSubj_list,
                    objNull$MAF_interval,
                    objNull$zScoreE_cutoff,
                    control$SPA_Cutoff,
                    control$zeta,
                    control$tol)
  
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,".")) # This file is in 'control.R'
}

mainMarker.SAGELD = function(genoType, genoIndex, outputColumns, objNull)
{
  OutList = mainMarkerInCPP("SAGELD", genoType, genoIndex);
  
  Method = objNull$Method
  
  n_marker = length(OutList$markerVec)
  
  if(Method == "SAGELD")
  {
    obj.mainMarker = data.frame(Marker = OutList$markerVec,                     # marker IDs
                                Info = OutList$infoVec,                         # marker information: CHR:POS:REF:ALT
                                AltFreq = OutList$altFreqVec,                   # alternative allele frequencies
                                AltCounts = OutList$altCountsVec,               # alternative allele counts
                                MissingRate = OutList$missingRateVec,           # alternative allele counts
                                Method = Method,                                # method
                                zScore_G = OutList$zScore[2*1:n_marker-1], # standardized score statistics for G
                                zScore_GxE = OutList$zScore[2*1:n_marker],   # standardized score statistics for GxE
                                Pvalue_G = OutList$pvalVec[2*1:n_marker-1],   # marker-level p-value for G
                                Pvalue_GxE = OutList$pvalVec[2*1:n_marker],   # marker-level p-value for GxE
                                hwepval = OutList$hwepvalVec)                   # marker-level HWE pvalue
  }else if(Method == "GALLOP")
  {
    obj.mainMarker = data.frame(Marker = OutList$markerVec,                     # marker IDs
                                Info = OutList$infoVec,                         # marker information: CHR:POS:REF:ALT
                                AltFreq = OutList$altFreqVec,                   # alternative allele frequencies
                                AltCounts = OutList$altCountsVec,               # alternative allele counts
                                MissingRate = OutList$missingRateVec,           # alternative allele counts
                                Method = Method,                                # method
                                Beta_G = OutList$beta[2*1:n_marker-1],     # beta estimate for G
                                Beta_GxE = OutList$beta[2*1:n_marker],       # beta estimate for GxE
                                SE_G = OutList$seBeta[2*1:n_marker-1],       # SE estimate for G
                                SE_GxE = OutList$seBeta[2*1:n_marker],         # SE estimate for GxE
                                Pvalue_G = OutList$pvalVec[2*1:n_marker-1],   # marker-level p-value for G
                                Pvalue_GxE = OutList$pvalVec[2*1:n_marker],   # marker-level p-value for GxE
                                hwepval = OutList$hwepvalVec)                   # marker-level HWE pvalue
  }
                              
  return(obj.mainMarker)
}

SAGELD.NullModel = function(NullModel,             # a fitted null model from lme4 or glmmTMB.
                            UsedMethod = "SAGELD", # default running "SAGELD", user can also run "GALLOP" using unrelated samples.
                            PlinkFile,             # a PLINK file path to read in some genotypes (without file suffix like ".bim", "bed" or "fam").
                            SparseGRMFile,         # a path of SparseGRMFile get from getSparseGRM() function.
                            PairwiseIBDFile,       # a path of PairwiseIBDFile get from getPairwiseIBD() function.
                            PvalueCutoff = 0.001,  # a p value cutoff for marginal genetic effect on environmental variable.
                            control = list())      # control command used in 'SPAGRM.NullModel', see also 'SPAGRM.NullModel'.
{
  cat("Process the null model product...\n")
  
  # Extract variance components and compute penalty matrix (P)
  if(inherits(NullModel, "merMod"))
  {
    Pheno_data = NullModel@frame
    
    SubjIDColname = colnames(Pheno_data)[ncol(Pheno_data)]
    
    varcor = lme4::VarCorr(NullModel)
    cmd = paste0("varcor$", SubjIDColname); G = eval(parse(text = cmd))
    sig = attr(varcor, "sc")
    P = solve(G / sig ^ 2)
  }else if(inherits(NullModel, "glmmTMB"))
  {
    cat("Warning: when you use glmmTMB to fit the null model, you should library(glmmTMB), before you library(GRAB)!\n")
    
    Pheno_data = NullModel$frame
    
    SubjIDColname = colnames(Pheno_data)[ncol(Pheno_data)]
    
    # varcor = glmmTMB::VarCorr(NullModel)$cond
    cmd="varcor = glmmTMB::VarCorr(NullModel)$cond"; eval(parse(text = cmd))
    cmd = paste0("varcor$", SubjIDColname); G = eval(parse(text = cmd))
    sig = attr(varcor, "sc")
    P = solve(G / sig ^ 2)
  }else
  {
    stop("Currently we only support fitted models fitted by 'LME4' and 'glmmTMB'.")
  }
  
  PhenoColname = colnames(Pheno_data)[1]
  Envcolname = colnames(G)[2]
  CovaColname = colnames(Pheno_data)[-c(1, ncol(Pheno_data))]; CovaColname = setdiff(CovaColname, Envcolname)
  
  # extract residuals of G-P; GxE-P; and G-E:
  model_formula = as.formula(paste(Envcolname, "~(1|", SubjIDColname, ")"))
  null_model = lme4::lmer(model_formula, data = Pheno_data)
  
  Resid_data = Pheno_data %>% 
    mutate(Resid_G = residuals(NullModel)) %>%
    mutate(Resid_GxE = Resid_G * !!sym(Envcolname)) %>% 
    mutate(Resid_E = residuals(null_model)) %>%
    group_by(!!sym(SubjIDColname)) %>% 
    summarize(Resid_G = sum(Resid_G), Resid_GxE = sum(Resid_GxE), Resid_E = sum(Resid_E)) %>%
    rename(SubjID = !!sym(SubjIDColname)) %>% ungroup()
  
  # Put data in convenient arrays and vector
  cmd = paste0("Pheno_data = Pheno_data %>% arrange(", SubjIDColname, ")"); eval(parse(text = cmd))
  cmd = paste0("unique(Pheno_data$", SubjIDColname, ")"); SubjID = eval(parse(text = cmd)); SubjID = as.character(SubjID)
  cmd = paste0("table(Pheno_data$", SubjIDColname, ")"); k = eval(parse(text = cmd))
  n = length(SubjID)
  k = k[SubjID]
  XY = bind_cols(intercept = 1, Pheno_data %>% select(all_of(c(Envcolname, CovaColname, PhenoColname))))
  XY = as.matrix(XY)
  TT = XY[, 1:2]
  nxy = ncol(XY)
  X = XY[, -nxy]
  nx = ncol(X)
  y = XY[, nxy]
  ncov = nx - ncol(TT)
  
  # Compute components of block-diagonal system with covariates  (without SNP)
  # Additionally compute and store object SS = Rot %*% Si, which is used later on
  A21 = matrix(NA, 2*n, 2 + ncov)
  q2 = rep(NA, 2*n)
  SS = matrix(NA, 2*n, 2)
  
  uk = 0
  for(i in 1:n)
  {
    ki = k[i]
    uk = max(uk) + 1:ki
    u2 = (i - 1)*2 + 1:2
    Ti = TT[uk,]
    if(ki == 1) Ti = t(as.matrix(Ti))
    Si = crossprod(Ti, Ti)
    sv = svd(Si + P)
    Rot = sqrt(1/sv$d) * sv$u
    Q = Rot %*% t(Ti)
    SS[u2,] = Rot %*% Si
    A21[u2,] = Q %*% X[uk, ]
    q2[u2] = Q %*% y[uk]
  }
  
  q1 = crossprod(X, y)
  A11 = crossprod(X)
  Q = A11 - crossprod(A21)
  q = q1 - crossprod(A21, q2)
  sol = solve(Q, q)
  blups = q2 - A21 %*% sol
  
  # Compute sums of products per subject involved in the crossprod(X, G), crossprod(G) and crossprod(G, y).
  ex = matrix(1, 1, ncov + 2)
  et = matrix(1, 1, 2)
  XTk = kronecker(et, X) * kronecker(TT, ex)
  TTk = kronecker(et, TT) * kronecker(TT, et)
  Tyk = y * TT
  XTs = matrix(0, n, ncol(XTk))
  TTs = matrix(0, n, ncol(TTk))
  Tys = matrix(0, n, 2)
  AtS = matrix(0, n, 2 * nx)
  
  uk = 0
  for (i in 1:n) 
  {
    ki = k[i]
    uk = max(uk) + 1:ki
    if(ki == 1)
    {
      XTs[i, ] = XTk[uk, ]
      TTs[i, ] = TTk[uk, ]
      Tys[i, ] = Tyk[uk, ]
    }else
    {
      XTs[i, ] = apply(XTk[uk, ], 2, sum)
      TTs[i, ] = apply(TTk[uk, ], 2, sum)
      Tys[i, ] = apply(Tyk[uk, ], 2, sum)
    }
    u2 = (i - 1) * 2 + (1 : 2)
    AtS[i, ] = c(crossprod(A21[u2, ], SS[u2, ]))
  }
  
  if(UsedMethod == "GALLOP")
  {
    obj = list(subjData = SubjID, N = length(SubjID), Method = "GALLOP",
               XTs = XTs, SS = SS, AtS = AtS, Q = Q, A21 = A21, TTs = TTs, Tys = Tys, sol = sol, blups = blups, sig = sig,
               Resid = numeric(0), Resid_G = numeric(0), Resid_GxE = numeric(0), Resid_E = numeric(0), Resid.unrelated.outliers = numeric(0),
               Resid.unrelated.outliers_G = numeric(0), Resid.unrelated.outliers_GxE = numeric(0),
               R_GRM_R = 0, R_GRM_R_G = 0, R_GRM_R_GxE = 0, R_GRM_R_G_GxE = 0, R_GRM_R_E = 0,
               R_GRM_R_TwoSubjOutlier = 0, R_GRM_R_TwoSubjOutlier_G = 0,
               R_GRM_R_TwoSubjOutlier_GxE = 0, R_GRM_R_TwoSubjOutlier_G_GxE = 0,
               sum_R_nonOutlier = 0, sum_R_nonOutlier_G = 0, sum_R_nonOutlier_GxE = 0,
               R_GRM_R_nonOutlier = 0, R_GRM_R_nonOutlier_G = 0,
               R_GRM_R_nonOutlier_GxE = 0, R_GRM_R_nonOutlier_G_GxE = 0,
               TwoSubj_list = list(), ThreeSubj_list = list(), MAF_interval = numeric(0), zScoreE_cutoff = 0)
    
    class(obj) = "SAGELD_NULL_Model"
    
  }else if(UsedMethod == "SAGELD")
  {
    SparseGRM = data.table::fread(SparseGRMFile)
    PairwiseIBD = data.table::fread(PairwiseIBDFile)
    
    Resid_data$SubjID = as.character(Resid_data$SubjID)
    SparseGRM$ID1 = as.character(SparseGRM$ID1); SparseGRM$ID2 = as.character(SparseGRM$ID2)
    PairwiseIBD$ID1 = as.character(PairwiseIBD$ID1); PairwiseIBD$ID2 = as.character(PairwiseIBD$ID2)
    
    control = checkControl.SAGELD.NullModel(control, Resid_data, SparseGRM, PairwiseIBD)
    
    MaxQuantile = control$MaxQuantile;
    MinQuantile = control$MinQuantile;
    OutlierRatio = control$OutlierRatio;
    ControlOutlier = control$ControlOutlier;
    MaxNuminFam = control$MaxNuminFam;
    MAF_interval = control$MAF_interval;
    
    Resid_G = Resid_data$Resid_G;
    Resid_GxE = Resid_data$Resid_GxE;
    Resid_E = Resid_data$Resid_E;
    
    SparseGRM = SparseGRM %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
    PairwiseIBD = PairwiseIBD %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
    
    # Decompose the subjects based on family structure and use a greedy algorithm to reduce family size if needed
    SparseGRM1 = SparseGRM
    SparseGRM1$G_pos1 = Resid_data$Resid_G[match(SparseGRM$ID1, Resid_data$SubjID)]
    SparseGRM1$G_pos2 = Resid_data$Resid_G[match(SparseGRM$ID2, Resid_data$SubjID)]
    SparseGRM1$GxE_pos1 = Resid_data$Resid_GxE[match(SparseGRM$ID1, Resid_data$SubjID)]
    SparseGRM1$GxE_pos2 = Resid_data$Resid_GxE[match(SparseGRM$ID2, Resid_data$SubjID)]
    SparseGRM1$E_pos1 = Resid_data$Resid_E[match(SparseGRM$ID1, Resid_data$SubjID)]
    SparseGRM1$E_pos2 = Resid_data$Resid_E[match(SparseGRM$ID2, Resid_data$SubjID)]
    SparseGRM1 = SparseGRM1 %>% mutate(G_Cov = Value * G_pos1 * G_pos2, 
                                       GxE_Cov = Value * GxE_pos1 * GxE_pos2,
                                       G_GxE_Cov1 = Value * GxE_pos1 * G_pos2,
                                       G_GxE_Cov2 = Value * GxE_pos2 * G_pos1,
                                       E_Cov = Value * E_pos1 * E_pos2)
    
    SparseGRM2 = rbind(SparseGRM1 %>% filter(ID1 == ID2),
                       SparseGRM1 %>% filter(ID1 != ID2),
                       SparseGRM1 %>% filter(ID1 != ID2))
    
    R_GRM_R_G = sum(SparseGRM2$G_Cov)
    R_GRM_R_GxE = sum(SparseGRM2$GxE_Cov)
    R_GRM_R_G_GxE = sum(SparseGRM2$G_GxE_Cov1 + SparseGRM2$G_GxE_Cov2)
    R_GRM_R_E = sum(SparseGRM2$E_Cov)
    # R_GRM_R = R_GRM_R_GxE + lambda^2 * R_GRM_R_G - lambda * R_GRM_R_G_GxE
    
    # threshold for G to E
    if(PvalueCutoff >= 0 & PvalueCutoff <= 1)
    {
      zScoreE_cutoff = qnorm(PvalueCutoff/2,lower.tail=FALSE)
    }else
    {
      stop("Please check ", PvalueCutoff, ". It should be a p value.")
    }
    
    # read in the Plink file to random select SNPs to calculate mean lambda.
    bedfile = paste0(PlinkFile, ".bed")
    bimfile = paste0(PlinkFile, ".bim")
    
    totalSNPs = data.table::fread(bimfile, header = FALSE); totalSNPs = totalSNPs$V2
    
    if(length(totalSNPs) > 2e3)
    {
      cat("Too many SNPs in the PLINK file! We randomly select 2000 SNPs.\n")
      
      random_SNPs = sample(totalSNPs, 2000)
      
      SNPIDfile = system.file("PairwiseIBD", "temp", package = "GRAB")
      SNPIDfile = paste0(SNPIDfile, "/tempSNPID", sample(1:1e9, 1), ".txt")
      
      data.table::fwrite(data.frame(random_SNPs), file = SNPIDfile, 
                         row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
      
      GenoMatInfo = GRAB.ReadGeno(GenoFile = bedfile,
                                  SampleIDs = SubjID,
                                  control = list(IDsToIncludeFile = SNPIDfile,
                                                 ImputeMethod = "mean"))
      file.remove(SNPIDfile)
      
    }else
    {
      GenoMatInfo = GRAB.ReadGeno(GenoFile = bedfile,
                                  SampleIDs = SubjID,
                                  control = list(AllMarkers = TRUE,
                                                 ImputeMethod = "mean"))
    }
    
    lambdaObs = c()
    for (i in 1:nrow(GenoMatInfo$markerInfo))
    {
      si = GenoMatInfo$GenoMat[, i]; mu = mean(si)/2
      
      zScoreE = sum(Resid_E * si) / sqrt(2*mu*(1-mu) * R_GRM_R_E)
      
      if(mu > 0.05 & mu < 0.95 & abs(zScoreE) < zScoreE_cutoff)
      {
        snp2 = rep(si, each = 2)
        H1 = matrix(crossprod(si, XTs), nx, 2)
        H2 = snp2 * SS
        AtH = matrix(crossprod(si, AtS), nx, 2)
        R = H1 - AtH
        Cfix = solve(Q, R)
        Cran = H2 - A21 %*% Cfix
        GtG = matrix(crossprod(si ^ 2, TTs), 2, 2)
        V = GtG - crossprod(H1, Cfix) - crossprod(H2, Cran)
        
        lambda = V[1,2]/V[1,1]
        
        lambdaObs = c(lambdaObs, lambda)
      }
    }
    
    if(length(lambdaObs) > 1e2)
    {
      lambda = mean(lambdaObs, na.rm = TRUE)
      cat("Mean lambda:\t", lambda, "\n")
    }else
    {
      stop("Less than 100 common SNPs (MAF > 0.05 & Pvalue(G-E) > ", PvalueCutoff, ") in the PLINK file!\n")
    }
    
    Resid_data = Resid_data %>% mutate(Resid = Resid_GxE - lambda * Resid_G)
    Resid = Resid_data$Resid
    
    SparseGRM1$pos1 = Resid_data$Resid[match(SparseGRM$ID1, Resid_data$SubjID)]
    SparseGRM1$pos2 = Resid_data$Resid[match(SparseGRM$ID2, Resid_data$SubjID)]
    SparseGRM1 = SparseGRM1 %>% mutate(Cov = abs(Value * pos1 * pos2))
    
    # Use residual information to define outliers / non-outliers
    Quant = quantile(Resid, probs = c(MinQuantile, MaxQuantile))
    Range = max(Quant) - min(Quant)
    cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
    
    cat("cutoffVec:\t",cutoffVec,"\n")
    Resid_data$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                                TRUE, FALSE)
    
    if(ControlOutlier)
    {
      cat("ControlOutlier = TRUE (default) to keep the outliers < 5%;\nSet ControlOutlier = FALSE for higher accuracy.\n")
      
      while(sum(Resid_data$Outlier) == 0)
      {
        OutlierRatio = OutlierRatio*0.8
        cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
        cat("cutoffVec:\t",cutoffVec,"\n")
        Resid_data$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                                    TRUE, FALSE)
        cat("The number of outlier is:", sum(Resid_data$Outlier),"\n")
      }
      
      while(sum(Resid_data$Outlier)/nrow(Resid_data) > 0.05)
      {
        OutlierRatio = OutlierRatio + 0.5
        cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
        cat("cutoffVec:\t",cutoffVec,"\n")
        Resid_data$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                                    TRUE, FALSE)
        cat("The number of outlier is:", sum(Resid_data$Outlier),"\n")
      }
    }
    
    cat("Outliers information is as below\n")
    print(Resid_data %>% filter(Outlier == TRUE) %>% dplyr::select(SubjID, Resid, Outlier) %>% arrange(Resid))
    
    # group samples into families
    edges = t(SparseGRM1[, c("ID1", "ID2")])
    graph_GRM = make_graph(edges, directed = FALSE)
    graph_list_all = graph_GRM %>% decompose()
    graph_length = lapply(graph_list_all, length)
    
    graph_list_1 = graph_list_all[graph_length == 1]
    SubjID.unrelated = lapply(graph_list_1, get.vertex.attribute) %>% unlist(use.names = FALSE)
    Resid_data.unrelated = Resid_data %>% filter(SubjID %in% SubjID.unrelated)
    SubjID.unrelated.nonOutlier = Resid_data.unrelated %>% filter(Outlier == FALSE) %>% select(SubjID) %>% unlist(use.names = FALSE)
    
    # Values used in association analysys
    R_GRM_R = R_GRM_R_GxE + lambda^2 * R_GRM_R_G - lambda * R_GRM_R_G_GxE
    
    sum_R_nonOutlier_G = Resid_data.unrelated %>% filter(Outlier == FALSE) %>% select(Resid_G) %>% sum
    sum_R_nonOutlier_GxE = Resid_data.unrelated %>% filter(Outlier == FALSE) %>% select(Resid_GxE) %>% sum
    sum_R_nonOutlier = sum_R_nonOutlier_GxE - lambda * sum_R_nonOutlier_G
    
    R_GRM_R_nonOutlier_G = SparseGRM1 %>% filter(ID1 %in% SubjID.unrelated.nonOutlier) %>% select(G_Cov) %>% sum
    R_GRM_R_nonOutlier_GxE = SparseGRM1 %>% filter(ID1 %in% SubjID.unrelated.nonOutlier) %>% select(GxE_Cov) %>% sum
    R_GRM_R_nonOutlier_G_GxE = SparseGRM1 %>% filter(ID1 %in% SubjID.unrelated.nonOutlier) %>% select(G_GxE_Cov1, G_GxE_Cov2) %>% sum
    R_GRM_R_nonOutlier = R_GRM_R_nonOutlier_GxE + lambda^2 * R_GRM_R_nonOutlier_G - lambda * R_GRM_R_nonOutlier_G_GxE
    
    Resid.unrelated.outliers_G = Resid_data.unrelated %>% filter(Outlier == TRUE) %>% select(Resid_G) %>% unlist(use.names = FALSE)
    Resid.unrelated.outliers_GxE = Resid_data.unrelated %>% filter(Outlier == TRUE) %>% select(Resid_GxE) %>% unlist(use.names = FALSE)
    Resid.unrelated.outliers = Resid.unrelated.outliers_GxE - lambda * Resid.unrelated.outliers_G
    
    R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier_G = R_GRM_R_TwoSubjOutlier_GxE = R_GRM_R_TwoSubjOutlier_G_GxE = 0; 
    
    TwoSubj_list = ThreeSubj_list = list();
    
    # initialize parameters
    graph_list_updated = list()
    graph_list = graph_list_all[graph_length > 1]
    nGraph = length(graph_list)
    index.outlier = 1
    
    if(nGraph != 0)
    {
      cat("Start process the related residual information.\n")
      
      for(i in 1:nGraph)
      {
        if(i %% 1000 == 0)
          cat("Processing the related residual information:\t", i,"/",nGraph,"\n")
        
        comp1 = graph_list[[i]]
        comp3 = V(comp1)$name
        
        # Step 0: calculate variance for the family
        pos1 = match(comp3, SubjID)
        outlierInFam = any(Resid_data$Outlier[pos1])
        
        if(!outlierInFam)
        {
          block_GRM = make.block.GRM(comp1, SparseGRM)
          
          R_GRM_R_G.temp = as.numeric(t(Resid_G[pos1]) %*% block_GRM %*% Resid_G[pos1])
          R_GRM_R_GxE.temp = as.numeric(t(Resid_GxE[pos1]) %*% block_GRM %*% Resid_GxE[pos1])
          R_GRM_R_G_GxE.temp = as.numeric(t(Resid_G[pos1]) %*% block_GRM %*% Resid_GxE[pos1]) * 2
          R_GRM_R.temp = R_GRM_R_GxE.temp + lambda^2 * R_GRM_R_G.temp - lambda * R_GRM_R_G_GxE.temp
          
          sum_R_nonOutlier_G = sum_R_nonOutlier_G + sum(Resid_G[pos1])
          sum_R_nonOutlier_GxE = sum_R_nonOutlier_GxE + sum(Resid_GxE[pos1])
          sum_R_nonOutlier = sum_R_nonOutlier + sum(Resid[pos1])
          
          R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
          R_GRM_R_nonOutlier_G = R_GRM_R_nonOutlier_G + R_GRM_R_G.temp
          R_GRM_R_nonOutlier_GxE = R_GRM_R_nonOutlier_GxE + R_GRM_R_GxE.temp
          R_GRM_R_nonOutlier_G_GxE = R_GRM_R_nonOutlier_G_GxE + R_GRM_R_G_GxE.temp
          next
        }
        
        # cat("Family ", i, " (with outliers) includes ", length(comp3), " subjects:", comp3, "\n")
        
        vcount = vcount(comp1)   # number of vertices 
        
        if(vcount <= MaxNuminFam)
        {
          graph_list_updated[[index.outlier]] = comp1
          index.outlier = index.outlier + 1
          next
        }
        
        # Step 1: remove the edges until the largest family size is <= MaxNuminFam, default is 5.
        
        comp1.temp = comp1
        tempGRM1 = SparseGRM1 %>% filter(ID1 %in% comp3 | ID2 %in% comp3) %>% arrange(Cov)
        for(j in 1:nrow(tempGRM1))
        {
          # cat("j:\t",j,"\n")
          edgesToRemove = paste0(tempGRM1$ID1[j],"|",tempGRM1$ID2[j])
          comp1.temp = delete.edges(comp1.temp, edgesToRemove)
          vcount = decompose(comp1.temp) %>% sapply(vcount)  # vertices count for the new graph after edge removal
          # cat("vcount:\t",vcount,"\n")
          if(max(vcount) <= MaxNuminFam)
            break;
        }
        
        # cat("Edge removal complete. Counts of vertices:\t", vcount,"\n")
        
        # Step 2: add the (removed) edges while keeping the largest family size <= MaxNuminFam, default is 5.
        
        tempGRM1 = tempGRM1[1:j,] %>% arrange(desc(Cov))
        comp1 = comp1.temp
        for(k in 1:nrow(tempGRM1))
        {
          # cat("k:\t",k,"\n")
          edgesToAdd = c(tempGRM1$ID1[k], tempGRM1$ID2[k])
          comp1.temp = add.edges(comp1, edgesToAdd)
          
          vcount = decompose(comp1.temp) %>% sapply(vcount)  # vertices count for the new graph after edge removal
          # cat("vcount:\t",vcount,"\n")
          
          if(max(vcount) <= MaxNuminFam)
            comp1 = comp1.temp
        }
        
        comp1 = decompose(comp1)
        
        # cat("Edge add complete. Counts of vertices:\t", comp1 %>% sapply(vcount),"\n")
        
        for(k in 1:length(comp1))
        {
          comp11 = comp1[[k]]
          comp13 = V(comp11)$name
          
          pos2 = match(comp13, SubjID)
          outlierInFam = any(Resid_data$Outlier[pos2])
          
          if(!outlierInFam){
            block_GRM = make.block.GRM(comp11, SparseGRM)
            
            R_GRM_R_G.temp = as.numeric(t(Resid_G[pos2]) %*% block_GRM %*% Resid_G[pos2])
            R_GRM_R_GxE.temp = as.numeric(t(Resid_GxE[pos2]) %*% block_GRM %*% Resid_GxE[pos2])
            R_GRM_R_G_GxE.temp = as.numeric(t(Resid_G[pos2]) %*% block_GRM %*% Resid_GxE[pos2]) * 2
            R_GRM_R.temp = R_GRM_R_GxE.temp + lambda^2 * R_GRM_R_G.temp - lambda * R_GRM_R_G_GxE.temp
            
            sum_R_nonOutlier_G = sum_R_nonOutlier_G + sum(Resid_G[pos2])
            sum_R_nonOutlier_GxE = sum_R_nonOutlier_GxE + sum(Resid_GxE[pos2])
            sum_R_nonOutlier = sum_R_nonOutlier + sum(Resid[pos2])
            
            R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
            R_GRM_R_nonOutlier_G = R_GRM_R_nonOutlier_G + R_GRM_R_G.temp
            R_GRM_R_nonOutlier_GxE = R_GRM_R_nonOutlier_GxE + R_GRM_R_GxE.temp
            R_GRM_R_nonOutlier_G_GxE = R_GRM_R_nonOutlier_G_GxE + R_GRM_R_G_GxE.temp
          }else{
            graph_list_updated[[index.outlier]] = comp11
            index.outlier = index.outlier + 1;
          }
        }
      }
      
      cat("Start process the Chow-Liu tree.\n")
      
      # Make a list of array index.
      arr.index = list()
      for(n in 1:MaxNuminFam)
      {
        temp = c()
        for(i in 1:n)
        {
          indexString = rep("c(1, 1, 1)", n)
          indexString[i] = "0:2"
          indexString = paste0(indexString, collapse = "%o%")
          cmd = paste0("temp = c(temp, list(arr.index", i, "=", indexString, "))")
          eval(parse(text = cmd))
        }
        arr.index[[n]] = temp
      }
      
      # build chou-liu-tree.
      n.outliers = length(graph_list_updated)
      if(n.outliers != 0)
      {
        ## The below values are only used in chou.liu.tree
        TwofamID.index = ThreefamID.index = 0
        for(index.outlier in 1:n.outliers)
        {
          if(index.outlier %% 1000 == 0)
            cat("Processing the CLT for families with outliers:\t", TwofamID.index, ", ", ThreefamID.index, "/", nGraph, "\n")
          
          comp1 = graph_list_updated[[index.outlier]]
          comp3 = V(comp1)$name
          n1 = length(comp3)
          pos3 = match(comp3, SubjID)
          
          Resid.temp = Resid[pos3]; Resid_G.temp = Resid_G[pos3]; Resid_GxE.temp = Resid_GxE[pos3]
          
          if(n1 == 1)
          {
            Resid.unrelated.outliers = c(Resid.unrelated.outliers, Resid.temp)
            Resid.unrelated.outliers_G = c(Resid.unrelated.outliers_G, Resid_G.temp)
            Resid.unrelated.outliers_GxE = c(Resid.unrelated.outliers_GxE, Resid_GxE.temp)
            next;
          }
          
          block_GRM = make.block.GRM(comp1, SparseGRM)
          
          tempIBD = PairwiseIBD %>% filter(ID1 %in% comp3 & ID2 %in% comp3)
          
          if(n1 == 2)
          {
            TwofamID.index = TwofamID.index + 1;
            
            R_GRM_R_TwoSubjOutlier_G.temp = as.numeric(t(Resid_G.temp) %*% block_GRM %*% Resid_G.temp)
            R_GRM_R_TwoSubjOutlier_GxE.temp = as.numeric(t(Resid_GxE.temp) %*% block_GRM %*% Resid_GxE.temp)
            R_GRM_R_TwoSubjOutlier_G_GxE.temp = as.numeric(t(Resid_G.temp) %*% block_GRM %*% Resid_GxE.temp) * 2
            R_GRM_R_TwoSubjOutlier.temp = R_GRM_R_TwoSubjOutlier_GxE.temp + lambda^2 * R_GRM_R_TwoSubjOutlier_G.temp - lambda * R_GRM_R_TwoSubjOutlier_G_GxE.temp
            
            R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier + R_GRM_R_TwoSubjOutlier.temp
            R_GRM_R_TwoSubjOutlier_G = R_GRM_R_TwoSubjOutlier_G + R_GRM_R_TwoSubjOutlier_G.temp
            R_GRM_R_TwoSubjOutlier_GxE = R_GRM_R_TwoSubjOutlier_GxE + R_GRM_R_TwoSubjOutlier_GxE.temp
            R_GRM_R_TwoSubjOutlier_G_GxE = R_GRM_R_TwoSubjOutlier_G_GxE + R_GRM_R_TwoSubjOutlier_G_GxE.temp
            
            Rho.temp = tempIBD$pa + 0.5*tempIBD$pb
            midterm = sqrt(Rho.temp^2 - tempIBD$pa)
            
            TwoSubj_list[[TwofamID.index]] = list(Resid = Resid.temp, Resid_G = Resid_G.temp, Resid_GxE = Resid_GxE.temp,
                                                  Rho = c(Rho.temp + midterm, Rho.temp - midterm))
            next;
          }
          
          ThreefamID.index = ThreefamID.index + 1;
          
          CLT = chow.liu.tree(N = n1,
                              IBD = tempIBD,
                              IDs = comp3,
                              MAF_interval = MAF_interval)
          
          stand.S.temp = array(rowSums(mapply(function(x, y) x*y, arr.index[[n1]], Resid.temp)), rep(3, n1))
          stand.S_G.temp = array(rowSums(mapply(function(x, y) x*y, arr.index[[n1]], Resid_G.temp)), rep(3, n1))
          stand.S_GxE.temp = array(rowSums(mapply(function(x, y) x*y, arr.index[[n1]], Resid_GxE.temp)), rep(3, n1))
          
          ThreeSubj_list[[ThreefamID.index]] = list(CLT = CLT,
                                                    stand.S = c(stand.S.temp),
                                                    stand.S_G = c(stand.S_G.temp),
                                                    stand.S_GxE = c(stand.S_GxE.temp))
        }
        cat("Completed processing the CLT for families with outliers:\t", TwofamID.index, ", ", ThreefamID.index, "/", nGraph, "\n")
      }
    }
    
    obj = list(subjData = SubjID, N = length(SubjID), Method = "SAGELD",
               XTs = XTs, SS = SS, AtS = AtS, Q = Q, A21 = A21, TTs = TTs, Tys = Tys, sol = sol, blups = blups, sig = sig,
               Resid = Resid, Resid_G = Resid_G, Resid_GxE = Resid_GxE, Resid_E = Resid_E, Resid.unrelated.outliers = Resid.unrelated.outliers,
               Resid.unrelated.outliers_G = Resid.unrelated.outliers_G, Resid.unrelated.outliers_GxE = Resid.unrelated.outliers_GxE,
               R_GRM_R = R_GRM_R, R_GRM_R_G = R_GRM_R_G, R_GRM_R_GxE = R_GRM_R_GxE, R_GRM_R_G_GxE = R_GRM_R_G_GxE, R_GRM_R_E = R_GRM_R_E,
               R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier, R_GRM_R_TwoSubjOutlier_G = R_GRM_R_TwoSubjOutlier_G,
               R_GRM_R_TwoSubjOutlier_GxE = R_GRM_R_TwoSubjOutlier_GxE, R_GRM_R_TwoSubjOutlier_G_GxE = R_GRM_R_TwoSubjOutlier_G_GxE,
               sum_R_nonOutlier = sum_R_nonOutlier, sum_R_nonOutlier_G = sum_R_nonOutlier_G, sum_R_nonOutlier_GxE = sum_R_nonOutlier_GxE,
               R_GRM_R_nonOutlier = R_GRM_R_nonOutlier, R_GRM_R_nonOutlier_G = R_GRM_R_nonOutlier_G,
               R_GRM_R_nonOutlier_GxE = R_GRM_R_nonOutlier_GxE, R_GRM_R_nonOutlier_G_GxE = R_GRM_R_nonOutlier_G_GxE,
               TwoSubj_list = TwoSubj_list, ThreeSubj_list = ThreeSubj_list, MAF_interval = MAF_interval, zScoreE_cutoff = zScoreE_cutoff)
    
    class(obj) = "SAGELD_NULL_Model"
    
  }else
  {
    stop("We currently only support 'SAGELD' and 'GALLOP', don't support '", UsedMethod, "'.")
  }
  
  return(obj)
}
