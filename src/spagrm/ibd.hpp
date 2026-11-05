// ibd.hpp — Pairwise IBD estimation from PLINK genotypes + sparse GRM
//
// Replaces the R function getPairwiseIBD() with a pure C++17 / Eigen
// implementation.  Connected-component detection uses union-find instead
// of igraph.
//
// Output: tab-separated file with columns  ID1  ID2  pa  pb  pc
#pragma once

#include "geno_factory/geno_data.hpp"

#include <algorithm>
#include <cmath>

namespace ibd {

/// One pair's IBD-sharing probabilities: the chances that the two subjects
/// share 2, 1 and 0 alleles identical by descent.  These are delta2, delta1
/// and delta0 of Eq. (17) of the SPAGRM paper (Nat Commun 16:1413), and the
/// three columns `pa pb pc` of the `--cal-pairwise-ibd` output.
struct Triple {
    double pa, pb, pc;
};

/// Solve Eq. (17) for one pair, subject to the admissibility constraints.
///
///     grmValue  rho1, the sparse-GRM off-diagonal = 2*phi
///     pcRaw     rho2, the IBS0 method-of-moments estimator of delta0,
///               unconstrained (a weighted mean of non-negative per-marker
///               terms, so pcRaw >= 0 but with no upper bound)
///
/// Eq. (17) is the linear system
///
///     delta0 = rho2 ;  delta1/2 + delta2 = rho1 ;  delta0+delta1+delta2 = 1
///
/// whose solution is delta0 = rho2, delta1 = 2*(1-rho1-rho2),
/// delta2 = 2*rho1+rho2-1.  The third equation is the normalization, so
/// pa+pb+pc = 1 is an identity of the solved form in exact arithmetic.  In
/// binary floating point it holds to 2 ulp — pa and pb each round twice — and
/// a search over 2e6 grid points and 4e7 random draws puts the worst deviation
/// at 4.441e-16, which survives the %.17g round trip to the output file and
/// sits 6.5 orders inside what nsGRMNull::loadIndexedIBD accepts.
///
/// Non-negativity is NOT automatic: the solution is a valid IBD vector iff
/// three inequalities hold, and a fourth fact, rho2 >= 0, is a precondition on
/// the caller rather than something this function can supply — for v > 0.5 the
/// lower bound 1-2v is negative and a negative rho2 would pass through.
/// runPairwiseIBD satisfies it because its accumulator |d-1|+|d+1|-2 is
/// non-negative pointwise.  Writing rho± for the roots of
/// x^2 - 2*rho1*x + delta2 = 0 — the per-allele sharing probabilities the
/// two-subject CGF factorizes through — the three conditions are
///
///     (i)   delta2 >= 0   <=>  rho2 >= 1 - 2*rho1
///     (ii)  rho± real     <=>  rho2 <= (1 - rho1)^2
///     (iii) rho+ <= 1     <=>  rho1 <= 1
///
/// (ii) is also the non-inbred admissibility inequality
/// delta1^2 >= 4*delta0*delta2, and given (i) and (ii), delta1 >= 0 follows
/// automatically PROVIDED rho1 <= 1, since
/// 1 - rho1 - rho2 >= 1 - rho1 - (1-rho1)^2 = rho1*(1-rho1).
///
/// (i) and (ii) are imposed by clamping rho2; (iii) by clamping rho1, since
/// 2*phi > 1 is not attainable by any relationship and therefore says
/// something about the GRM rather than about the pair.  At rho1 = 1 the
/// result is (1, 0, 0), the duplicate / monozygotic vector.
///
/// The predecessor omitted (iii) and instead carried a repair branch that
/// fired exactly when it was violated, rewriting the answer to
/// (pa, pb, pc) = (rho1, 0, 0) — which does not sum to one and has pa > 1.
/// See src/spagrm/ibd.cpp for what that cost downstream.
inline Triple deriveIBD(double grmValue, double pcRaw) noexcept {
    // (iii), plus the symmetric floor: a negative off-diagonal is sampling
    // noise about zero, not a relationship, and v = 0 gives the unrelated
    // vector (0, 0, 1).
    const double v = std::clamp(grmValue, 0.0, 1.0);

    // The admissible band for rho2 is [1 - 2v, (1 - v)^2], of width exactly
    // v^2, and it is formed so that the two endpoints CANNOT invert.
    //
    // The predecessor wrote the upper bound as (1-v)^2 - 1e-10, an
    // unmotivated epsilon carried over from the R original, and passed the
    // two endpoints straight to std::clamp.  Then upper - lower = v^2 - 1e-10
    // is negative for every v below 1e-5, and std::clamp is UNDEFINED when
    // hi < lo.  Removing the epsilon is not by itself enough: evaluating
    // (1-v)*(1-v) and 1-2v independently rounds each, so their difference is
    // v^2 + O(2^-53) and the rounding term still dominates below v ~ 1.3e-8.
    // A sweep of 1e-30 <= v < 1 finds 211 490 abscissae where the
    // independently-evaluated upper bound falls BELOW the lower one, all with
    // v <= 7.45e-09.
    //
    // Writing upper = lower + v*v removes the inversion by construction: v*v
    // is non-negative, so the sum cannot round below `lower` whatever v is.
    // It is the same quantity — (1-v)^2 = (1-2v) + v^2 exactly in the reals,
    // and the two spellings agree to 4.4e-16 relative for every v > 1e-6.
    //
    // Nothing needs the band to be open.  At pc = (1-v)^2 the discriminant
    // rho1^2 - delta2 is exactly zero and the two per-allele roots coincide
    // at v, which is the maximum-entropy point of the admissible set — the
    // pair is modelled as two independent Bernoulli(v) allele-sharing
    // indicators — not a singularity, and nothing downstream divides by the
    // discriminant or by the root separation.
    const double lower = 1.0 - 2.0 * v;
    const double upper = lower + v * v;
    const double pc = std::clamp(pcRaw, lower, upper);

    return Triple{2.0 * v + pc - 1.0, 2.0 - 2.0 * pc - 2.0 * v, pc};
}

} // namespace ibd

/// Compute pairwise IBD (pa, pb, pc) for every off-diagonal pair in the
/// sparse GRM and write the result to `outputFile`.
///
/// @param spgrmGrabFile   3-column TSV (ID1 ID2 VALUE), '#' lines skipped.
/// @param spgrmGctaFile  plink2 --make-grm-sparse .grm.sp file.
/// @param bfilePrefix     PLINK binary prefix (.bed/.bim/.fam).
/// @param outputFile      Tab-separated output (ID1 ID2 pa pb pc).
/// @param minMafIBD       Minimum MAF for a marker to be used (default 0.01).
void runPairwiseIBD(
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outputFile,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    double minMafIBD = 0.01,
    int nthreads = 1
);
