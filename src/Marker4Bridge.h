

#ifndef MARKER4BRIDGE_H
#define MARKER4BRIDGE_H

// Marker4Bridge.h -- R/Rcpp boundary: global method pointers set before engine dispatch

#include <RcppArmadillo.h>
#include "Marker4Engine.h"


extern POLMM::POLMMClass*        ptr_gPOLMMobj;
extern SPACox::SPACoxClass*       ptr_gSPACoxobj;
extern SPAmix::SPAmixClass*       ptr_gSPAmixobj;
extern SPAGRM::SPAGRMClass*      ptr_gSPAGRMobj;
extern SAGELD::SAGELDClass*       ptr_gSAGELDobj;
extern WtCoxG::WtCoxGClass*      ptr_gWtCoxGobj;
extern SPAsqr::SPAsqrClass*       ptr_gSPAsqrobj;
extern LEAF::LEAFClass*           ptr_gLEAFobj;
extern SPAmixPlus::SPAmixPlusClass* ptr_gSPAmixPlusobj;

#endif
