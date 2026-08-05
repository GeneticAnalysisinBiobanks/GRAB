// wtcoxg_cgf.hpp — the saddlepoint CGF of WtCoxG and LEAF (tier 3)
//
// Stage 6 of the saddlepoint rework.  WtCoxG and LEAF share one kernel: LEAF
// is K independent WtCoxG analyses, one per k-means cluster, whose per-cluster
// score statistics are pooled afterwards by math::metaPvalueScorePool.  There
// is exactly one saddlepoint here, and this header is it.
//
// The statistic is S = sum_i R_i (G_i - 2*MAF) with G_i ~ Binomial(2, MAF) for
// a single MAF shared by every subject, so the tier-2 kernel is
// `spa_cgf::binomUniform`.  As in SPAmix and SPAsqr the subjects are split by
// the 1.5*IQR rule: the outliers keep the exact binomial CGF and the
// non-outliers contribute a second-order (Gaussian) approximation in closed
// form.  What is specific to WtCoxG is that the Gaussian block carries two
// DIFFERENT variances — see D4 below.
//
// ══════════════════════════════════════════════════════════════════════
// What this replaces
// ══════════════════════════════════════════════════════════════════════
//
// `outlierCgf_{scalar,avx2,avx512}` + `fastGetRootK1_wt` + `getProbSpaG_wt`
// in wtcoxg.cpp.  Line for line:
//
//   * The kernel already subtracted per element (`K2 += mgf2/mgf0 -
//     ratio*ratio`, wtcoxg.cpp:81 before this change), which 01_findings.md D1
//     names as the least bad of the four binomial kernels — but it is still a
//     difference of two O(4 r^2) quantities whose gap closes as t*r grows.
//     `spa_cgf` forms K'' as 2 r^2 e u / alpha^2, in which every operation is a
//     product or a quotient of positive numbers and no cancellation is
//     possible anywhere in the domain.
//
//   * The Newton loop computed `K0 += std::log(mgf0)` on every iteration and
//     read only K1 and K2 (01_findings.md P2).  The k12 / kFull split removes
//     one logarithm per outlier per iteration; the logarithm now runs once per
//     tail, in kFull.
//
//   * The finder returned `{root, K2, converge}` and the caller read only
//     `.root` (D5), then re-evaluated the whole CGF at that root anyway
//     (P3).  The unified solver returns cumulants evaluated AT the root
//     together with a Status that the caller cannot drop, because it is part
//     of the value rather than an out-parameter.
//
//     P3 is worth stating precisely, because the finding reads as though the
//     terminal pass could be deleted.  It cannot.  Once P2 removes the
//     logarithm from the loop, the loop no longer HAS a K0 to reuse, so the
//     terminal evaluation stops being redundant and becomes the only place K
//     is computed.  What P3 actually buys here, and what is realized, is that
//     the terminal cumulants are evaluated at the abscissa the solver returns
//     rather than at whatever abscissa the loop happened to stop on, and that
//     the loop pass costs one exp and one reciprocal per outlier, vectorized,
//     instead of one exp, one log and two divides, scalar.  The measured
//     evaluation count is 2 tails * (bracket + Newton) k12 calls plus exactly
//     2 kFull calls per marker.  `spa::detail::CachedFull` (the single-callable
//     overload that skips the terminal evaluation when the loop's last
//     abscissa equals the root) is deliberately NOT used, for the same reason
//     it is not used in spamix_cgf.hpp: with split callables the loop's last
//     evaluation carries no K0.
//
//   * The tail had the strongest guard set in the codebase and the weakest
//     reporting of it: five conditions, each returning a bare NaN, and a
//     caller that wrote `std::min(1.0, pval1 + pval2)` — which returns 1.0
//     when either tail is NaN, so every failure surfaced as a perfectly null
//     marker (01_findings.md D2).  `spa::bnTailLog` keeps every one of those
//     guards and names which one fired; `spa::combineTailsLog` propagates the
//     failure as NaN plus a status instead of masking it.
//
// ══════════════════════════════════════════════════════════════════════
// D4 — the (K, K'') pair is deliberately inconsistent.  PRESERVE IT.
// ══════════════════════════════════════════════════════════════════════
//
// The Gaussian block of WtCoxG's CGF is not one Gaussian.  wtcoxg.cpp builds
//
//     var_n_resid    = 2*MAF*(1-MAF) * sum_nonoutlier (R_i - bm)^2
//     var_adj_batch  = 4*b^2 * sumR^2 * var_mu_ext            (batch effect)
//     var_adj_finite = obs_ct * (sumR/N_all)^2 * MAF*(1-MAF)  (finite panel)
//
//     var_n_K01 = var_n_resid + var_adj_batch     -> enters K and K'
//     var_n_K2  = var_n_resid + var_adj_finite    -> enters K''
//
// so K'' is not the second derivative of K whenever var_adj_batch differs from
// var_adj_finite, which happens exactly when obs_ct > 0 — Branch B, the branch
// that retains the external reference MAF.  The same asymmetry sits between
// the Newton residual and its Jacobian.  It is documented as intentional at
// the head of the old `fastGetRootK1_wt`, and 01_findings.md D4 and
// 02_design.md both rule that it is a statistical-modelling question for the
// maintainer, not a coding defect to be repaired inside a migration.
//
// It is preserved here exactly, in three places, and each is load-bearing:
//
//   1. `k12` returns K'' = binomial + varK2 while its residual uses varK01, so
//      the Newton step divides by the same non-Jacobian the old loop divided
//      by (wtcoxg.cpp:464-467 before this change).
//   2. `kFull` returns K built from varK01 and K'' built from varK2, so the
//      tail receives the same mismatched pair as before.
//   3. `spa::bnTailLog` takes K0 and K2 as opaque scalars and its documentation
//      says why: a kernel that assumed K2 = K0'' and "corrected" one of them
//      would compute a different statistic from the one WtCoxG defines.
//
// Do not "simplify" Context by collapsing varK01 and varK2 into one field.
//
// ══════════════════════════════════════════════════════════════════════
// The sign constraint, and the initial abscissa
// ══════════════════════════════════════════════════════════════════════
//
// `spa::SolveOpts::scoreSign` restricts the search to sgn(zeta) == sgn(s).
// Algebraically it would be admissible here — K'(0) is identically zero for
// this CGF, since the outlier block contributes 2*MAF*(sumR - N*bm) =
// 2*MAF*b*sumR and the non-outlier mean is mean_n = 2*MAF*shifted_sum -
// 2*b*sumR*MAF, which cancel — but "identically zero" holds in exact
// arithmetic only.  K'(0) is a difference of two O(sumR) quantities and its
// computed value is a rounding residue, so for a marker whose |S| is close to
// that residue the true root can sit on the far side of the origin and the
// constraint would exclude it.  It is therefore left off.  Nothing is lost:
// K' is strictly increasing, the root is unique, and the solver brackets it.
//
// The initial abscissa is the first-order estimate |dev| / var, capped at 1.2
// exactly as SPAsqr, SPAGRM and SPAmix cap theirs.  The predecessor always
// started at t = 0.  The cap is not cosmetic: an uncapped guess on a marker
// with a very small variance lands where spa_cgf's t*r clamp has saturated
// K'' to ~1e-307, and the solver's first bracket probe from there overflows.

#pragma once

#include <cmath>
#include <cstdint>

#include "util/spa.hpp"
#include "util/spa_cgf.hpp"

namespace wtcoxg_cgf {

// ══════════════════════════════════════════════════════════════════════
// Context
// ══════════════════════════════════════════════════════════════════════
//
// `resid` points at the outlier residuals already centred by bm = (1-b)*meanR.
// Centring depends only on the cluster and on b, never on the marker, so the
// array is a per-(cluster, b) precompute owned by WtCoxGShared rather than an
// Eigen temporary rebuilt per marker per EXT/NOEXT variant (01_findings.md P7).
struct Context {
    const double *resid = nullptr;   // centred outlier residuals, nOutlier entries
    int nOutlier = 0;
    double maf = 0.0;                // the single allele frequency, in (0, 1)

    double mean = 0.0;               // mean_n      — enters K and K'
    double varK01 = 0.0;             // var_n_K01   — enters K and K'
    double varK2 = 0.0;              // var_n_K2    — enters K''  (D4)
};

// Newton-loop payload: {K'(t) - s, K''(t)}.  No logarithm (P2).
inline spa::K12 k12(double t, const Context &c, double s) noexcept {
    const spa_cgf::Cgf12 d =
        spa_cgf::binomUniformK12(t, c.resid, c.nOutlier, c.maf);
    return spa::K12{d.K1 + c.mean + c.varK01 * t - s, d.K2 + c.varK2};
}

// Terminal payload: {K(t), K'(t) - s, K''(t)}.  Runs once per tail.
//
// The Gaussian block contributes mean*t + varK01*t^2/2 to K, written
// t * (mean + 0.5 * varK01 * t) so that a large |t| multiplies once rather
// than squaring an already-large quantity.
inline spa::K012 kFull(double t, const Context &c, double s) noexcept {
    const spa_cgf::Cgf012 d =
        spa_cgf::binomUniformKFull(t, c.resid, c.nOutlier, c.maf);
    return spa::K012{d.K0 + t * (c.mean + 0.5 * c.varK01 * t),
                     d.K1 + c.mean + c.varK01 * t - s,
                     d.K2 + c.varK2};
}

// ══════════════════════════════════════════════════════════════════════
// Two-sided saddlepoint p-value
// ══════════════════════════════════════════════════════════════════════
//
// The tails are placed by reflection about zero, upper = +|dev| and
// lower = -|dev|, which is what the predecessor did (`getProbSpaG_wt(+|S|,
// lower_tail=false)` and `(-|S|, lower_tail=true)`) and is correct because
// K'(0) = 0 for this CGF.
//
// `indepVar` sizes the initial abscissa only; pass 0 to start at the origin.
//
// Failure semantics are spa::assemble's: a failure in EITHER tail discards the
// saddlepoint for both, and the row reports the two-sided normal tail at
// `zNorm` under a status naming the substitution (decision D5).  Never a
// half-sized or a masked p-value.
//
// `zNorm` is the caller's raw score z = S / sqrt(Var S).  It is a parameter
// rather than absDev/sqrt(indepVar) because the caller has already divided S
// by the variance ratio, so the two are not the same number.
inline spa::Result twoSidedSpa(
    const Context &c,
    double absDev,
    double indepVar,
    double zNorm,
    double rtol = 1e-6
) noexcept {
    double zeta0 = 0.0;
    if (indepVar > 0.0) {
        const double e = absDev / indepVar;
        if (std::isfinite(e) && e > 0.0) zeta0 = (e > 1.2) ? 1.2 : e;
    }

    double logp[2];
    spa::Status st[2];

    for (int k = 0; k < 2; ++k) {
        const bool lowerTail = (k == 1);
        const double s = lowerTail ? -absDev : absDev;

        spa::SolveOpts opt;
        opt.init = lowerTail ? -zeta0 : zeta0;
        opt.rtol = rtol;
        // opt.scoreSign deliberately left at 0; see the header note.

        const spa::Saddle sd = spa::solveSaddlepoint(
            s,
            [&](double t) { return k12(t, c, s); },
            [&](double t) { return kFull(t, c, s); },
            opt);

        spa::Status stLog = spa::Status::SpaOk;
        logp[k] = spa::bnTailLog(sd.zeta, s, sd.K0, sd.K2, lowerTail, stLog);
        st[k] = spa::worseStatus(sd.status, stLog);
    }

    // The two-sided assembly and the D5 fallback when either tail failed —
    // one call, in spa.hpp.  One tail evaluation per side: the linear sibling
    // of `bnTailLog` that used to run beside it is gone with the P column's
    // parallel assembly (log10p_unify Stage 3).
    return spa::assemble(logp[0], logp[1], st[0], st[1], zNorm);
}

// The SPA_STATUS output column.  MethodBase hands the engine a
// std::vector<double> and marker_impl.hpp formats every cell through
// `numToChars`, so the column is the numeric enumerator rather than the token
// `spa::statusName` returns: 0 SPA_OK, 1 NORMAL, 2 SPA_W_SINGULAR,
// 3..6 the FALLBACK_* codes, 7 NA_POST_FAIL, 8 NA_NO_TEST.  The order is a
// contract: <= 2 means LOG10P is trustworthy, 3..6 that it is the substituted
// normal tail, >= 7 that it is NA.  See the Status note in util/spa.hpp.
inline double statusCode(spa::Status s) noexcept {
    return static_cast<double>(static_cast<uint8_t>(s));
}

}  // namespace wtcoxg_cgf
