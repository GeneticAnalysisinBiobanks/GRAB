// mtMain.h -- Pure-C++ chunk-parallel marker execution engine interface
//
// Functions in this file:
//
//   [Engine]
//     mainMarkerChunksCore — thread-pool dispatcher: reads genotypes, runs statistics,
//                            and writes results for all marker chunks

#ifndef MAIN_H
#define MAIN_H

#include "mtPLINK.h"
#include "mtPOLMM.h"
#include "mtWtCoxG.h"
#include "mtLEAF.h"
#include "mtSPAGRM.h"
#include "mtSAGELD.h"
#include "mtSPAsqr.h"
#include "mtSPACox.h"
#include "mtSPAmix.h"
#include "mtSPAmixPlus.h"


namespace mtMain {

// Global method pointers — defined in mtMain.cpp.
extern POLMM::POLMMClass*            ptr_gPOLMMobj;
extern WtCoxG::WtCoxGClass*          ptr_gWtCoxGobj;
extern LEAF::LEAFClass*              ptr_gLEAFobj;
extern SPACox::SPACoxClass*          ptr_gSPACoxobj;
extern SPAmix::SPAmixClass*          ptr_gSPAmixobj;
extern SPAGRM::SPAGRMClass*          ptr_gSPAGRMobj;
extern SAGELD::SAGELDClass*          ptr_gSAGELDobj;
extern SPAsqr::SPAsqrClass*          ptr_gSPAsqrobj;
extern SPAmixPlus::SPAmixPlusClass*  ptr_gSPAmixPlusobj;

// ---- Engine ----
void mainMarkerChunksCore(
  const std::string&                        method,
  const std::vector<std::vector<uint64_t>>& chunkMarkers,
  const std::string&                        outputFile,
  unsigned int                              nthreads,
  const std::string&                        impute_method,
  double                                    missing_cutoff,
  double                                    min_maf_marker,
  double                                    min_mac_marker,
  const mtPLINK::PlinkData&                 plinkData
);

} // namespace mtMain

#endif // MAIN_H
