// spamix_cgf.hpp — the independence CGF of SPAmix / SPAGxE / SPAGxEmix (tier 3)
//
// Stage 5 of the saddlepoint rework.  These three methods share ONE saddlepoint
// kernel: the score S = sum_i r_i G_i under an independent binomial genotype
// law, split by the IQR outlier rule into
//
//   * an exact binomial block over the outliers, and
//   * a Gaussian block over the non-outliers, contributing
//     mean*t + var*t^2/2 to K, mean + var*t to K', and var to K''.
//
// The binomial block is the shared tier-2 kernel, in the variant the method's
// genotype model calls for:
//
//   SPAGxE (base and +)   G_i ~ Binomial(2, q)     one q for all subjects
//                         -> spa_cgf::binomUniform
//   SPAmix / SPAmixPlus   G_i ~ Binomial(2, q_i)   per-individual AF
//   SPAGxEmix             G_i ~ Binomial(2, q_i)   -> spa_cgf::binomIndiv
//
// 02_design.md lists all three of these methods under binomIndiv.  That is right
// for two of them and wrong for SPAGxE base, whose q is a scalar: the previous
// code reached the per-individual kernel only because it materialized
// `std::vector<double> mafOut(nOut, q)` per marker per environment and passed
// the constant vector in.  Selecting binomUniform there removes that allocation
// from the hot path, hoists the p endpoints and the `1 - p` out of the loop, and
// drops the per-lane alpha == 0 guard the per-individual variant must carry.
//
// ══════════════════════════════════════════════════════════════════════
// What this replaces, and the residual-factor convention
// ══════════════════════════════════════════════════════════════════════
//
// The predecessor was `spa::getProbSpaG` in spamix/common.cpp: two scalar loops
// over the outliers calling `math::kG012`, a private Newton iteration, and a
// private and only partly guarded copy of the Barndorff-Nielsen tail.  None of
// it was vectorized — this family had no SIMD at all (01_findings.md P6).
//
// `math::kG012(x, p, K0, K1, K2)` and `spa_cgf::subjectK12(t, r, p, h)` do NOT
// have the same shape, and the difference is easy to miss because both are
// "the diploid binomial CGF":
//
//     math::kG012      arguments differentiated with respect to x = t*r, so the
//                      caller must supply the chain-rule factors itself:
//                          K'  += r   * k1
//                          K'' += r*r * k2
//                      and the value returned in K2 is (m0*m2 - m1*m1)/(m0*m0),
//                      a difference of products (01_findings.md D1).
//
//     spa_cgf          differentiated with respect to t, with the residual
//                      weight r a kernel argument.  K' and K'' already carry
//                      their r and r^2.
//
// Porting D1's cancellation-free rewrite into `kG012` while leaving the call
// sites alone would therefore have multiplied K'' by a spurious r^2.  The
// convention adopted here is spa_cgf's — the residual belongs to the kernel,
// because that is what lets the SIMD tiers load `resid[i]` as a vector lane
// instead of applying a scalar factor after the reduction.  The external
// `resid[i] * k1` / `resid[i]*resid[i] * k2` factors are deleted at the call
// sites, not moved.
//
// ══════════════════════════════════════════════════════════════════════
// P2 — the discarded logarithm
// ══════════════════════════════════════════════════════════════════════
//
// `evalOutlierK1adjK2` called `math::kG012`, which computes K0 = log(m0)
// unconditionally, and then used only k1 and k2 (01_findings.md P2,
// common.cpp:44-46).  One `std::log` per outlier per Newton iteration, thrown
// away every time.  The k12 / kFull split removes it: the loop calls
// `binom*K12`, which computes no logarithm at all, and the terminal evaluation
// calls `binom*KFull` once per tail.
//
// ══════════════════════════════════════════════════════════════════════
// P3 — the extra terminal pass
// ══════════════════════════════════════════════════════════════════════
//
// 01_findings.md P3 observes that the old finder broke on `|diffX| < tol`
// BEFORE applying the step (common.cpp:97-98), so the returned abscissa was
// exactly where the last CGF evaluation had been made, and that the caller then
// discarded that evaluation and ran a second full O(nOutlier) pass at the same
// point (common.cpp:125) to recover K0 — the one cumulant the Newton loop never
// accumulated.
//
// Only half of that is recoverable, and this is worth stating plainly because
// the finding reads as though the whole pass could be deleted.  Once P2 removes
// the logarithm from the loop, the loop no longer HAS a K0 to reuse: the second
// pass stops being redundant and becomes the only place K is computed.  What P3
// actually buys, and what is realized here, is
//
//   * accuracy — the terminal cumulants are evaluated at the root the solver
//     returns rather than at a stale abscissa, and the step is applied before
//     the convergence test rather than after it; and
//   * per-iteration cost — the loop pass is now one exp and one reciprocal per
//     subject, vectorized, instead of one exp, one log and a divide, scalar.
//
// The evaluation count is 2 tails * (bracketing + Newton) k12 calls plus
// exactly 2 kFull calls per marker, and the second of those two is genuinely
// required.  `spa::detail::CachedFull` — the single-callable solver overload
// that skips the terminal evaluation when the loop's last abscissa equals the
// returned root — is deliberately NOT used here, because with the split
// callables the loop's last evaluation carries no K0 to reuse.
//
// ══════════════════════════════════════════════════════════════════════
// The sign constraint, and why this family does not use it
// ══════════════════════════════════════════════════════════════════════
//
// `spa::SolveOpts::scoreSign` restricts the search to sgn(zeta) == sgn(Score).
// SPAGRM and SPAsqr set it because their statistic is centred, so the tail
// direction and the sign of zeta coincide.  Here the two tails are placed by
// reflection about the fitted mean, upper = sMean + |dev| and lower =
// sMean - |dev|, and sMean is NOT zero — so `s` itself carries no useful sign.
// The correct constraint would be sgn(zeta) = +/-1 by tail rather than sgn(s),
// and even that is not safe on the variance-ratio path: SPAmixPlus and SPAGxE+
// reflect about the RESCALED mean sMean*sqrt(indepVar/grmVar) while the CGF's
// own mean K'(0) is the unscaled sMean, so for a rescaling far from one the
// root of the upper-tail problem is not guaranteed to be positive.  The
// constraint is therefore left off.  Nothing is lost: K' is non-decreasing, so
// the root is unique, and the solver brackets it and cannot converge to
// anything else.  The constraint existed to rescue unbracketed Family-B
// iterations.
//
// ══════════════════════════════════════════════════════════════════════
// The initial abscissa
// ══════════════════════════════════════════════════════════════════════
//
// The predecessor always started at t = 0.  The first-order saddlepoint is
// zeta ~ (s - K'(0)) / K''(0) = |dev| / indepVar, which is exact for the
// Gaussian block and a good estimate whenever the outlier block is a small part
// of the variance, so it is used here, capped at 1.2 exactly as SPAsqr and
// SPAGRM cap theirs.
//
// The cap is not cosmetic.  Uncapped, a marker with a very small independence
// variance produces an enormous initial abscissa, at which spa_cgf's t*r clamp
// has saturated K'' to ~1e-307; the solver's first bracket probe is the Newton
// step |K'/K''| from there, which is then astronomically large, overflows to
// non-finite and abandons the expansion.  Capping keeps the initial evaluation
// in the well-conditioned part of the domain and costs at most a few extra
// bracket probes when the true root is larger.

#pragma once

#include <cmath>
#include <cstdint>

#include "util/spa.hpp"
#include "util/spa_cgf.hpp"

namespace spamix_cgf {

// ══════════════════════════════════════════════════════════════════════
// Context
// ══════════════════════════════════════════════════════════════════════
//
// `af == nullptr` selects the uniform-frequency kernel with the scalar `q`;
// otherwise `af` supplies one allele frequency per outlier, in the same order
// as `resid`, and `q` is unused.  The branch is resolved once per Newton
// iteration, not once per subject.
struct Context {
    const double *resid = nullptr;   // outlier residual weights, n entries
    const double *af = nullptr;      // per-outlier allele frequency, or null
    int nOutlier = 0;
    double q = 0.0;                  // uniform allele frequency, when af == null

    // Gaussian block over the non-outliers.
    double mean = 0.0;               // 2 * sum_non r_i q_i
    double var = 0.0;                // 2 * sum_non r_i^2 q_i (1 - q_i)
};

// Newton-loop payload: {K'(t) - s, K''(t)}.  No logarithm (P2).
inline spa::K12 k12(double t, const Context &c, double s) noexcept {
    const spa_cgf::Cgf12 d =
        c.af ? spa_cgf::binomIndivK12(t, c.resid, c.af, c.nOutlier)
             : spa_cgf::binomUniformK12(t, c.resid, c.nOutlier, c.q);
    return spa::K12{d.K1 + c.mean + c.var * t - s, d.K2 + c.var};
}

// Terminal payload: {K(t), K'(t) - s, K''(t)}.  Runs once per tail.
//
// The Gaussian block contributes mean*t + var*t^2/2 to K.  It is written
// t * (mean + 0.5 * var * t) so that a large |t| multiplies once rather than
// squaring an already-large quantity, and so that the two terms are combined
// before the (typically much larger) binomial K is added.
inline spa::K012 kFull(double t, const Context &c, double s) noexcept {
    const spa_cgf::Cgf012 d =
        c.af ? spa_cgf::binomIndivKFull(t, c.resid, c.af, c.nOutlier)
             : spa_cgf::binomUniformKFull(t, c.resid, c.nOutlier, c.q);
    return spa::K012{d.K0 + t * (c.mean + 0.5 * c.var * t),
                     d.K1 + c.mean + c.var * t - s,
                     d.K2 + c.var};
}

// ══════════════════════════════════════════════════════════════════════
// Two-sided saddlepoint p-value
// ══════════════════════════════════════════════════════════════════════
//
// The two tails are placed by reflection about the fitted mean, which is the
// convention all three callers already used:
//
//     upper = P(S >= sMean + |dev|)      lower = P(S <= sMean - |dev|)
//
// `indepVar` is the INDEPENDENCE variance K''(0) — the diagonal one, before any
// GRM quadratic form — and is used only to size the initial abscissa.  Pass 0
// to start at the origin.
//
// Failure semantics are spa::combineTails': a failure in either tail yields NaN
// in p and a named status, never a half-sized p-value.  The predecessor summed
// the two tails unconditionally, so a NaN in one of them propagated to the sum
// with no diagnostic; that is D5 and D6 as they appear in this family.
inline spa::TwoSided twoSidedSpa(
    const Context &c,
    double sMean,
    double absDev,
    double indepVar,
    double rtol = 1e-6
) noexcept {
    double zeta0 = 0.0;
    if (indepVar > 0.0) {
        const double e = absDev / indepVar;
        if (std::isfinite(e) && e > 0.0) zeta0 = (e > 1.2) ? 1.2 : e;
    }

    double p[2], logp[2];
    spa::Status st[2];

    for (int k = 0; k < 2; ++k) {
        const bool lowerTail = (k == 1);
        const double s = lowerTail ? (sMean - absDev) : (sMean + absDev);

        spa::SolveOpts opt;
        opt.init = lowerTail ? -zeta0 : zeta0;
        opt.rtol = rtol;
        // opt.scoreSign deliberately left at 0; see the header note.

        const spa::Saddle sd = spa::solveSaddlepoint(
            s,
            [&](double t) { return k12(t, c, s); },
            [&](double t) { return kFull(t, c, s); },
            opt);

        spa::Status stLin = spa::Status::Converged;
        spa::Status stLog = spa::Status::Converged;
        p[k]    = spa::bnTail(sd.zeta, s, sd.K0, sd.K2, lowerTail, stLin);
        logp[k] = spa::bnTailLog(sd.zeta, s, sd.K0, sd.K2, lowerTail, stLog);
        st[k] = spa::worseStatus(sd.status, spa::worseStatus(stLin, stLog));
    }

    // P from the linear-scale assembly (the sum of two positive tails, as
    // before); -log10(P) from the log-domain assembly, which stays meaningful
    // past the point where both linear tails underflow to zero.
    const spa::TwoSided lin = spa::combineTails(p[0], p[1], st[0], st[1]);
    const spa::TwoSided lg  = spa::combineTailsLog(logp[0], logp[1], st[0], st[1]);
    return spa::TwoSided{lin.p, lg.negLog10p, lin.status};
}

// The SPA_STATUS output column.  MethodBase hands the engine a
// std::vector<double> and marker_impl.hpp formats every cell through
// `numToChars`, so the column is the numeric enumerator rather than the token
// `spa::statusName` returns: 0 OK, 1 MAXITER, 2 GUARD_TEMP, 3 GUARD_CURV,
// 4 GUARD_W, 5 NONFINITE, 6 NORMAL.  This is the encoding Stages 3 and 4
// established for SPACox, SPAGRM and SPAsqr.
inline double statusCode(spa::Status s) noexcept {
    return static_cast<double>(static_cast<uint8_t>(s));
}

}  // namespace spamix_cgf
