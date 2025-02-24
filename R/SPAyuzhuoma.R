#' SPAyuzhuoma method in GRAB package
#' 
#' SPAyuzhuoma method is an empirical approach to analyzing complex traits for related samples in a large-scale biobank. 
#' 
#' @details 
#' Additional list of \code{control} in \code{SPAyuzhuoma.NullModel()} function.
#' 
#' Additional list of \code{control} in \code{GRAB.Marker()} function.
#' 
#' @examples 
#' # Step 2a: process model residuals
#' ResidMatFile = system.file("extdata", "ResidMat.txt", package = "GRAB")
#' SparseGRMFile = system.file("SparseGRM", "SparseGRM.txt", package = "GRAB")
#' PairwiseIBDFile = system.file("PairwiseIBD", "PairwiseIBD.txt", package = "GRAB")
#' obj.SPAyuzhuoma = SPAyuzhuoma.NullModel(ResidMatFile = ResidMatFile, 
#'                               SparseGRMFile = SparseGRMFile, 
#'                               PairwiseIBDFile = PairwiseIBDFile,
#'                               control = list(ControlOutlier = FALSE))
#' 
#' # Step 2b: perform score test
#' GenoFile = system.file("extdata", "simuPLINK.bed", package = "GRAB")
#' OutputDir = system.file("results", package = "GRAB")
#' OutputFile = paste0(OutputDir, "/SPAyuzhuomaMarkers.txt")
#' GRAB.Marker(objNull = obj.SPAyuzhuoma,
#'             GenoFile = GenoFile,
#'             OutputFile = OutputFile)
#' head(read.table(OutputFile, header=T))
#' @export
GRAB.SPAyuzhuoma = function(){
  print("Check ?GRAB.SPAyuzhuoma for more details about 'SPAyuzhuoma' method.")
}

################### This file includes the following functions

# ------------ used in 'GRAB_Marker.R' -----------
# 1. checkControl.Marker.SPAyuzhuoma(control)
# 2. setMarker.SPAyuzhuoma(objNull, control)
# 3. mainMarker.SPAyuzhuoma()

# check the control list in marker-level testing
checkControl.Marker.SPAyuzhuoma = function(control)
{
  default.control = list(SPA_Cutoff = 2,
                         zeta = 0,
                         tol = 1e-5)
  
  control = updateControl(control, default.control)  # This file is in 'control.R'
  
  return(control)
}

checkControl.SPAyuzhuoma.NullModel = function(control,
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
  
  if(any(colnames(ResidMat) != c("SubjID", "Resid")))
    stop("The column names of ResidMat should be ['SubjID', 'Resid'].")
  
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

setMarker.SPAyuzhuoma = function(objNull, control)
{
  # the below function is in 'Main.cpp'
  setSPAyuzhuomaobjInCPP(objNull$Resid,
                    objNull$Resid.unrelated.outliers,
                    objNull$sum_R_nonOutlier,
                    objNull$R_GRM_R_nonOutlier,
                    objNull$R_GRM_R_TwoSubjOutlier,
                    objNull$R_GRM_R,
                    objNull$MAF_interval,
                    objNull$TwoSubj_list,
                    objNull$ThreeSubj_list,
                    control$SPA_Cutoff,
                    control$zeta,
                    control$tol)
  
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,".")) # This file is in 'control.R'
}

mainMarker.SPAyuzhuoma = function(genoType, genoIndex, outputColumns)
{
  OutList = mainMarkerInCPP("SPAyuzhuoma", genoType, genoIndex);
  
  obj.mainMarker = data.frame(Marker = OutList$markerVec,           # marker IDs
                              Info = OutList$infoVec,               # marker information: CHR:POS:REF:ALT
                              AltFreq = OutList$altFreqVec,         # alternative allele frequencies
                              AltCounts = OutList$altCountsVec,     # alternative allele counts
                              MissingRate = OutList$missingRateVec, # alternative allele counts
                              zScore = OutList$zScore,              # standardized score statistics
                              Pvalue = OutList$pvalVec,             # marker-level p-value
                              hwepval = OutList$hwepvalVec)
  
  return(obj.mainMarker)
}

SPAyuzhuoma.NullModel = function(ResidMatFile,    # two columns: column 1 is subjID, column 2 is Resid
                            SparseGRMFile,   # a path of SparseGRMFile get from getSparseGRM() function.
                            PairwiseIBDFile, # a path of PairwiseIBDFile get from getPairwiseIBD() function.
                            control = list(MaxQuantile = 0.75,
                                           MinQuantile = 0.25,
                                           OutlierRatio = 1.5,
                                           ControlOutlier = TRUE,
                                           MaxNuminFam = 5,
                                           MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)))
{
  if(is.data.frame(ResidMatFile))
  {
    ResidMat = ResidMatFile
  }else
  {
    ResidMat = data.table::fread(ResidMatFile)
  }
  SparseGRM = data.table::fread(SparseGRMFile)
  PairwiseIBD = data.table::fread(PairwiseIBDFile)
  
  ResidMat$SubjID = as.character(ResidMat$SubjID)
  SparseGRM$ID1 = as.character(SparseGRM$ID1); SparseGRM$ID2 = as.character(SparseGRM$ID2)
  PairwiseIBD$ID1 = as.character(PairwiseIBD$ID1); PairwiseIBD$ID2 = as.character(PairwiseIBD$ID2)
  
  control = checkControl.SPAyuzhuoma.NullModel(control, ResidMat, SparseGRM, PairwiseIBD)
  
  MaxQuantile = control$MaxQuantile;
  MinQuantile = control$MinQuantile;
  OutlierRatio = control$OutlierRatio;
  ControlOutlier = control$ControlOutlier;
  MaxNuminFam = control$MaxNuminFam;
  MAF_interval = control$MAF_interval;
  
  SubjID = ResidMat$SubjID
  SparseGRM = SparseGRM %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
  PairwiseIBD = PairwiseIBD %>% filter(ID1 %in% SubjID & ID2 %in% SubjID)
  
  # Use residual information to define outliers / non-outliers
  Resid = ResidMat$Resid
  Quant = quantile(Resid, probs = c(MinQuantile, MaxQuantile))
  Range = max(Quant) - min(Quant)
  cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
  
  cat("cutoffVec:\t",cutoffVec,"\n")
  ResidMat$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                            TRUE, FALSE)
  
  if(ControlOutlier)
  {
    cat("ControlOutlier = TRUE (default) to keep the outliers < 5%;\nSet ControlOutlier = FALSE for higher accuracy.\n")
    
    while(sum(ResidMat$Outlier) == 0)
    {
      OutlierRatio = OutlierRatio*0.8
      cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
      cat("cutoffVec:\t",cutoffVec,"\n")
      ResidMat$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                                TRUE, FALSE)
      cat("The number of outlier is:", sum(ResidMat$Outlier),"\n")
    }
    
    while(sum(ResidMat$Outlier)/nrow(ResidMat) > 0.05)
    {
      OutlierRatio = OutlierRatio + 0.5
      cutoffVec = c(min(Quant) - OutlierRatio * Range, max(Quant) + OutlierRatio * Range)
      cat("cutoffVec:\t",cutoffVec,"\n")
      ResidMat$Outlier = ifelse(Resid < cutoffVec[1] | Resid > cutoffVec[2],
                                TRUE, FALSE)
      cat("The number of outlier is:", sum(ResidMat$Outlier),"\n")
    }
  }
  
  cat("Outliers information is as below\n")
  print(ResidMat %>% filter(Outlier == TRUE) %>% dplyr::select(SubjID, Resid, Outlier) %>% arrange(Resid))
  
  # Decompose the subjects based on family structure and use a greedy algorithm to reduce family size if needed
  SparseGRM1 = SparseGRM
  SparseGRM1$pos1 = ResidMat$Resid[match(SparseGRM$ID1, ResidMat$SubjID)]
  SparseGRM1$pos2 = ResidMat$Resid[match(SparseGRM$ID2, ResidMat$SubjID)]
  SparseGRM1 = SparseGRM1 %>% mutate(Cov = abs(Value * pos1 * pos2))
  
  edges = t(SparseGRM1[, c("ID1", "ID2")])
  graph_GRM = make_graph(edges, directed = F)
  graph_list_all = graph_GRM %>% decompose()
  graph_length = lapply(graph_list_all, length)
  
  graph_list_1 = graph_list_all[graph_length == 1]
  SubjID.unrelated = lapply(graph_list_1, get.vertex.attribute) %>% unlist(use.names = FALSE)
  ResidMat.unrelated = ResidMat %>% filter(SubjID %in% SubjID.unrelated)
  SubjID.unrelated.nonOutlier = ResidMat.unrelated %>% filter(Outlier == FALSE) %>% select(SubjID) %>% unlist(use.names = F)
  
  # Values used in association analysys
  R_GRM_R = SparseGRM1 %>% filter(ID1 %in% SubjID.unrelated) %>% select(Cov) %>% sum
  sum_R_nonOutlier = ResidMat.unrelated %>% filter(Outlier == FALSE) %>% select(Resid) %>% sum
  R_GRM_R_nonOutlier = SparseGRM1 %>% filter(ID1 %in% SubjID.unrelated.nonOutlier) %>% select(Cov) %>% sum
  Resid.unrelated.outliers = ResidMat.unrelated %>% filter(Outlier == TRUE) %>% select(Resid) %>% unlist(use.names = F)
  R_GRM_R_TwoSubjOutlier = 0; TwoSubj_list = ThreeSubj_list = list();
  
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
      outlierInFam = any(ResidMat$Outlier[pos1])
      
      block_GRM = make.block.GRM(comp1, SparseGRM)
      
      R_GRM_R.temp = as.numeric(t(Resid[pos1]) %*% block_GRM %*% Resid[pos1])
      R_GRM_R = R_GRM_R + R_GRM_R.temp
      
      if(!outlierInFam)
      {
        sum_R_nonOutlier = sum_R_nonOutlier + sum(ResidMat$Resid[pos1])
        R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
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
        outlierInFam = any(ResidMat$Outlier[pos2])
        
        block_GRM = make.block.GRM(comp11, SparseGRM)
        
        R_GRM_R.temp = as.numeric(t(Resid[pos2]) %*% block_GRM %*% Resid[pos2])
        
        if(!outlierInFam){
          sum_R_nonOutlier = sum_R_nonOutlier + sum(ResidMat$Resid[pos2])
          R_GRM_R_nonOutlier = R_GRM_R_nonOutlier + R_GRM_R.temp
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
        
        Resid.temp = ResidMat$Resid[pos3]
        
        if(n1 == 1)
        {
          Resid.unrelated.outliers = c(Resid.unrelated.outliers, Resid.temp)
          next;
        }
        
        block_GRM = make.block.GRM(comp1, SparseGRM)
        
        tempIBD = PairwiseIBD %>% filter(ID1 %in% comp3 & ID2 %in% comp3)
        
        if(n1 == 2)
        {
          TwofamID.index = TwofamID.index + 1;
          
          R_GRM_R_TwoSubjOutlier.temp = as.numeric(t(Resid.temp) %*% block_GRM %*% Resid.temp)
          R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier + R_GRM_R_TwoSubjOutlier.temp
          
          Rho.temp = tempIBD$pa + 0.5*tempIBD$pb
          midterm = sqrt(Rho.temp^2 - tempIBD$pa)
          
          TwoSubj_list[[TwofamID.index]] = list(Resid = Resid.temp,
                                                Rho = c(Rho.temp + midterm, Rho.temp - midterm))
          next;
        }
        
        ThreefamID.index = ThreefamID.index + 1;
        
        CLT = chow.liu.tree(N = n1,
                            IBD = tempIBD,
                            IDs = comp3,
                            MAF_interval = MAF_interval)
        
        stand.S.temp = array(rowSums(mapply(function(x, y) x*y, arr.index[[n1]], Resid.temp)), rep(3, n1))
        
        ThreeSubj_list[[ThreefamID.index]] = list(CLT = CLT,
                                                  stand.S = c(stand.S.temp))
      }
      cat("Completed processing the CLT for families with outliers:\t", TwofamID.index, ", ", ThreefamID.index, "/", nGraph, "\n")
    }
  }
  
  obj = list(Resid = Resid, subjData = SubjID, N = length(SubjID), Resid.unrelated.outliers = Resid.unrelated.outliers,
             R_GRM_R = R_GRM_R, R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier,
             sum_R_nonOutlier = sum_R_nonOutlier, R_GRM_R_nonOutlier = R_GRM_R_nonOutlier,
             TwoSubj_list = TwoSubj_list, ThreeSubj_list = ThreeSubj_list, 
             MAF_interval = MAF_interval)
  
  class(obj) = "SPAyuzhuoma_NULL_Model"
  
  return(obj)
}

SPAyuzhuomaGE.NullModel = function(NullModel = NULL,   # a fitted null model from lme4.
                              PhenoFile,          # a file path to read in the phenotype.
                              SubjIDColname,      # a character to specifie the column name of the subject ID.
                              PhenoColname,       # a character to specifie the column name of the phenotype.
                              CovaColname,        # a character (vector) to specifie the column name of the covariates (except for Envcolname).
                              Envcolname,         # a character to specifie the column name of the environment variable.
                              PlinkFile,          # a PLINK file path to read in some genotypes (without file suffix like ".bim", "bed" or "fam").
                              SparseGRMFile,      # a path of SparseGRMFile get from getSparseGRM() function.
                              PairwiseIBDFile,    # a path of PairwiseIBDFile get from getPairwiseIBD() function.
                              control = list())   # control command used in 'SPAyuzhuoma.NullModel', see also 'SPAyuzhuoma.NullModel'.
{
  if(is.data.frame(PhenoFile))
  {
    Pheno_data = PhenoFile
  }else
  {
    Pheno_data = data.table::fread(PhenoFile)
  }
  
  cmd = paste0("Pheno_data = Pheno_data %>% arrange(", SubjIDColname, ")"); eval(parse(text = cmd))
  
  if(is.null(NullModel))
  {
    cat("Model formula is...\n")
    model_formula = as.formula(paste(PhenoColname, "~", paste(CovaColname, collapse = "+"), "+", Envcolname, "+ (", Envcolname, "|", SubjIDColname, ")"))
    print(model_formula)
    
    cat("Fit the null model...\n")
    null_model = lme4::lmer(model_formula, data = Pheno_data)
  }else
  {
    null_model = NullModel
  }
  
  cat("Process the null model product...\n")
  
  # Extract variance components and compute penalty matrix (P)
  if(inherits(null_model, "merMod"))
  {
    varcor = VarCorr(null_model)
    cmd = paste0("varcor$", SubjIDColname); G = eval(parse(text = cmd))
    sig = attr(varcor, "sc")
    P = solve(G / sig ^ 2)
  }else if(inherits(null_model, "glmmTMB"))
  {
    varcor = VarCorr(null_model)$cond
    cmd = paste0("varcor$", SubjIDColname); G = eval(parse(text = cmd))
    sig = attr(varcor, "sc")
    P = solve(G / sig ^ 2)
  }else
  {
    stop("Currently we only support fitted models fitted by 'LME4' and 'glmmTMB'.")
  }
  
  # Put data in convenient arrays and vector
  cmd = paste0("unique(Pheno_data$", SubjIDColname, ")"); SubjID = eval(parse(text = cmd)); SubjID = as.character(SubjID)
  cmd = paste0("table(Pheno_data$", SubjIDColname, ")"); k = eval(parse(text = cmd))
  n = length(SubjID)
  k = k[SubjID]
  XY = bind_cols(intercept = 1, Pheno_data %>% select(all_of(Envcolname)), 
                 Pheno_data %>% select(all_of(CovaColname)), Pheno_data %>% select(all_of(PhenoColname)))
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
    
    if(mu > 0.05)
    {
      snp2 = rep(si, each = 2)
      H1 = matrix(crossprod(si, XTs), nx, 2)
      H2 = snp2 * SS
      AtH = matrix(crossprod(si, AtS), nx, 2)
      R = H1 - AtH
      Cfix = solve(Q, R)
      Cran = H2 - A21 %*% Cfix
      GtG = matrix(crossprod(si ^ 2, TTs), 2, 2)
      Gty = matrix(crossprod(si, Tys), 2, 1)
      V = GtG - crossprod(H1, Cfix) - crossprod(H2, Cran)
      v = Gty - crossprod(H1, sol) - crossprod(H2, blups)
      
      lambda = V[1,2]/V[1,1]
      
      lambdaObs = c(lambdaObs, lambda)
    }
  }
  
  if(length(lambdaObs) > 1e2)
  {
    lambda = mean(lambdaObs, na.rm = TRUE)
    cat("lambda:\t", lambda, "\n")
  }else
  {
    stop("Less than 100 common SNPs (MAF > 0.05) in the PLINK file!\n")
  }
  
  cat("Calculate model residuals...\n")
  
  if(inherits(null_model, "merMod"))
  {
    coeffs = summary(null_model)$coefficients[,1]
  }else if(inherits(null_model, "glmmTMB"))
  {
    coeffs = summary(null_model)$coefficients$cond[,1]
  }else
  {
    stop("Currently we only support fitted models fitted by 'LME4' and 'glmmTMB'.")
  }
  coeffs = c(coeffs[1], coeffs[Envcolname], coeffs[CovaColname], -1)
  update_residuals = - as.numeric(XY %*% coeffs)
  
  uk = 0
  for(i in 1:n)
  {
    if(i %% 10000 == 0) cat(i, "\n")
    
    ki = k[i]
    uk = max(uk) + 1:ki
    
    if(ki > 1)
    {
      tempmatrix = Matrix::Diagonal(k[i]) * sig^2 + TT[uk,] %*% G %*% t(TT[uk,])
      
      update_residuals[uk] = as.numeric(solve(tempmatrix, update_residuals[uk]))
    }else
    {
      tempmatrix = sig^2 + t(TT[uk,]) %*% G %*% TT[uk,]
      
      update_residuals[uk] = as.numeric(solve(tempmatrix, update_residuals[uk]))
    }
  }
  
  Resid_data = Pheno_data %>% select(all_of(SubjIDColname)) %>% mutate(Resid = update_residuals * (TT[, 2] - lambda))
  colnames(Resid_data) = c("SubjID", "Resid")
  Resid_data = Resid_data %>% group_by(SubjID) %>% summarize(Resid = sum(Resid)) %>% ungroup()
  
  output = SPAyuzhuoma.NullModel(ResidMatFile = Resid_data, SparseGRMFile = SparseGRMFile, PairwiseIBDFile = PairwiseIBDFile, control = control)
  
  return(output)
}

make.block.GRM = function(graph, 
                          GRM)    # three columns: "ID1", "ID2", and "Value"
{
  comp2 = get.data.frame(graph)
  
  # igraph gives an unexpected additional loop, which may change the block GRM
  # the below is to remove the additional loop
  comp2 = comp2[!duplicated(comp2),]
  
  comp3 = V(graph)$name
  
  colnames(GRM) = c("to", "from", "Value")
  
  n1 = nrow(comp2)
  comp2 = merge(comp2, GRM)
  n2 = nrow(comp2)
  
  if(n1 != n2)
    stop("Ask Wenjian Bi (wenjianb@pku.edu.cn) to check why 'n1 != n2'.")
  
  block_GRM = sparseMatrix(i = match(comp2$from, comp3),
                           j = match(comp2$to, comp3),
                           x = comp2$Value,
                           symmetric = T)
  return(block_GRM)
}

chow.liu.tree = function(N,
                         IBD,
                         IDs,
                         MAF_interval)
{
  CLT = c()
  
  # MAF_interval = c(0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
  for(index in 1:length(MAF_interval))
  {
    mu = MAF_interval[index]
    
    # p = c(G0, G1, G2)
    p0 = c((1-mu)^2, 2*mu*(1-mu), mu^2)
    
    # p = c(G00, G10, G20, G01, G11, G21, G02, G12, G22)
    pa.allele2 = c((1-mu)^2, 0, 0, 0, 2*mu*(1-mu), 0, 0, 0, mu^2)
    
    pb.allele1 = c((1-mu)^3, mu*(1-mu)^2, 0, mu*(1-mu)^2, mu*(1-mu), mu^2*(1-mu), 0, mu^2*(1-mu), mu^3)
    
    pc.allele0 = c((1-mu)^4, 2*mu*(1-mu)^3, mu^2*(1-mu)^2, 2*mu*(1-mu)^3, 4*mu^2*(1-mu)^2, 
                   2*mu^3*(1-mu), mu^2*(1-mu)^2, 2*mu^3*(1-mu), mu^4)
    
    # calculate entropy I(Gi, Gj). Noting that entropy of unrelated pairs is zero.
    for(j in 1:nrow(IBD))
    {
      pro = IBD$pa[j] * pa.allele2 + IBD$pb[j] * pb.allele1 + IBD$pc[j] * pc.allele0
      
      entropy = sum(pro * log(pro/pc.allele0), na.rm = T)
      IBD$entropy[j] = entropy
    }
    
    # use the "prim" lgorithm to bulid a maximum spanning tree.
    Max_span_tree = IBD %>% graph_from_data_frame(directed = T) %>% 
      mst(weights = - IBD$entropy, algorithm = "prim") %>% get.edgelist() %>% 
      data.table::as.data.table() %>% rename(ID1 = V1, ID2 = V2)
    
    mst.IBD = merge(Max_span_tree, IBD, all.x = T) %>%
      mutate(idxID1 = match(ID1, IDs), idxID2 = match(ID2, IDs))
    
    arr.prob = array(1, dim = rep(3, N))
    for(i in 1:N)
      dimnames(arr.prob)[[i]] = paste0("ID",i,":",0:2) 
    
    vec = c(mst.IBD$idxID1, mst.IBD$idxID2); vec = vec[duplicated(vec)]
    
    for(k in 1:(N - 1))
    {
      pro = mst.IBD$pa[k] * pa.allele2 + mst.IBD$pb[k] * pb.allele1 + mst.IBD$pc[k] * pc.allele0
      
      matrix.prob = matrix(pro, 3, 3)
      matrix.index1 = mst.IBD$idxID1[k]; matrix.index2 = mst.IBD$idxID2[k]
      for(i in 1:3){
        for(j in 1:3){
          indexString = rep("", N)
          indexString[matrix.index1] = i
          indexString[matrix.index2] = j
          indexString = paste0(indexString, collapse = ",")
          cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] * matrix.prob[",i,",",j,"]")
          # "arr.prob[1,1,] = arr.prob[1,1,] * matrix.prob[1,1]"
          eval(parse(text = cmd))
        }
      }
    }
    
    for(k in 1:(N - 2))
    {
      vector.prob = p0
      vector.index = vec[k]
      for(i in 1:3){
        indexString = rep("", N)
        indexString[vector.index] = i
        indexString = paste0(indexString, collapse = ",")
        cmd = paste0("arr.prob[",indexString,"] = arr.prob[", indexString, "] / vector.prob[",i,"]")
        # arr.prob[,,1] = arr.prob[,,1] / vector.prob[1]"
        eval(parse(text = cmd))
      }
    }
    
    CLT = cbind(CLT, c(arr.prob))
  }
  
  return(CLT)
}

