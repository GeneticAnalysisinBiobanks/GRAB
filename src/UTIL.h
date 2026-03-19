#ifndef UTIL_H
#define UTIL_H

// UTIL.h -- Shared utility functions

#include <RcppArmadillo.h>

double hwe_exact(int obs_hets, int obs_hom1, int obs_hom2);
void gethwepval(arma::vec GVec, double& hwepval, double hwepvalCutoff);

#endif
