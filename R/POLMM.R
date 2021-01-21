
setMarker.POLMM = function(objNull, control)
{
  
  if(NullModelClass == "POLMM"){
    if(objNull$controlList$LOCO){
      if(!chrom %in% names(objNull$LOCOList))
        stop("'chrom' should be in names(objNull$LOCOList).")
      obj.CHR = objNull$LOCOList[[chrom]]
    }else{
      # to be continued
    }
    
    # single marker analysis does not require sparse GRM any more 
    # Note: it might be not so accurate if min_mac_marker is very low
    SPmatR.CHR = list(locations = c(0,0), values = 1)
    
    setPOLMMobjInR(obj.CHR$muMat,
                   obj.CHR$iRMat,
                   objNull$Cova,
                   objNull$yVec,          # 1 to J
                   SPmatR.CHR,
                   objNull$tau,
                   POLMM.control$printPCGInfo,
                   POLMM.control$tolPCG,
                   POLMM.control$maxiterPCG)
    
    print(paste0("The current POLMM.control$nMarkers_output is ", nMarkers_output,"."))
  }
}

mainMarker.POLMM = function(objNull, control)
{
  OutList = MAIN_Marker(markers,
                        POLMM.control$SPA_cutoff,
                        POLMM.control$missing_cutoff,
                        POLMM.control$min_maf_region,
                        POLMM.control$min_mac_region,
                        obj.CHR$varRatio)
  
  StatVec = OutList$StatVec
  VarSVec = OutList$VarSVec
  adjPVec = OutList$pvalVec;
  markerVec = OutList$markerVec
  freqVec = OutList$freqVec
  flipVec = OutList$flipVec
  infoVec = OutList$infoVec       # marker infomation: CHR:POS:REF:ALT
  
  OUT.Marker = data.frame(Marker = markerVec,
                          Info = infoVec,
                          Freq = freqVec,
                          Flip = flipVec,
                          Stat = StatVec,
                          Var = VarSVec,
                          Pvalue = adjPVec)
}
