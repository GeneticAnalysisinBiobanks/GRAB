
setMarker.SPAGRM = function(objNull, control)
{
  # the below function is in 'Main.cpp'
  setSPAGRMobjInCPP(objNull$Cova,
                    objNull$yVec,          
                    control$SPA_Cutoff)
  
  print(paste0("The current control$nMarkersEachChunk is ", control$nMarkersEachChunk,"."))
}


mainMarker.SPAGRM = function(objNull, control)
{
  OutList = mainMarkerInCPP("SPAGRM", genoType, genoIndex);  
  
  obj.mainMarker = data.frame(Marker = OutList$markerVec,           # marker IDs
                              Info = OutList$infoVec,               # marker information: CHR:POS:REF:ALT
                              AltFreq = OutList$altFreqVec,         # alternative allele frequencies
                              AltCounts = OutList$altCountsVec,     # alternative allele counts
                              MissingRate = OutList$missingRateVec, # alternative allele counts
                              Pvalue = OutList$pvalVec)             # marker-level p-values
  
  optionalColumns = c("beta", "seBeta", "zScore", "PvalueNorm", "AltFreqInGroup", "AltCountsInGroup", "nSamplesInGroup")
  additionalColumns = intersect(optionalColumns, outputColumns)
  
  if(length(additionalColumns) > 0)
    obj.mainMarker = cbind.data.frame(obj.mainMarker, 
                                      as.data.frame(OutList[additionalColumns]))
}