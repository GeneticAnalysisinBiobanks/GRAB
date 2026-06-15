// impute_kernel.hpp — SADGE stage-1 parental-genotype imputation.
//
// Port of precompute_func / dsg2prb / generate_Gpar_1par2sib / calc_*_func
// from example_sadge_zouyu/impute_func.R.  All computations are phenotype-
// independent: given a marker's dosages for a family's genotyped members and
// the per-marker coefficients, return the expected parental genotype G_par for
// one sibling.  Scalar double arithmetic mirrors the R reference for numerical
// reproducibility (no exp/log → no SIMD-dispatch kernel needed).
#pragma once

#include "sadge/sample_info.hpp"

namespace sadge {

// Per-marker coefficients (precompute_func).  Only A,B,D,E,F enter the
// imputation kernels; H,I,J are the parental-genotype variance-inflation
// factors used by the score variance / rho (see precond.hpp).
struct PerMarkerCoeffs {
    double mu;
    double A, B, C, D, E, F, G, H, I, J;
};

void precomputeCoeffs(double mu, PerMarkerCoeffs &c);

// Convert a dosage d ∈ [0,2] to a 3-genotype probability vector (dsg2prb).
inline void dsg2prb(double d, double &p0, double &p1, double &p2) {
    if (d <= 1.0) {
        p0 = 1.0 - d;
        p1 = d;
        p2 = 0.0;
    } else {
        p0 = 0.0;
        p1 = 2.0 - d;
        p2 = d - 1.0;
    }
}

// Expected parental genotype G_par for one sibling under family structure fs
// (the genotyped/stage-1 structure).  gSelf is this sibling's dosage; gCoSib
// the genotyped co-sibling's (used by 0par2sib / 1par2sib); gPar1/gPar2 the
// genotyped parents' dosages (the lone parent of a 1-parent family is gPar1).
// Singletons (s0par1sib) return 0 and are handled as G_observed + 2*mu by the
// caller.
double imputeGpar(
    FamStruct fs,
    double gSelf,
    double gCoSib,
    double gPar1,
    double gPar2,
    const PerMarkerCoeffs &c);

} // namespace sadge
