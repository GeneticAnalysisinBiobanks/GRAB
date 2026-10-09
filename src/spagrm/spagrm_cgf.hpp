// spagrm_cgf.hpp — SPAGRM's cumulant generating function (tier 3)
//
// Stage 4 of the saddlepoint rework.  SPAGRM's CGF is a sum over THREE
// structurally different term classes, and 02_design.md puts it in tier 3 —
// per-method, deliberately NOT forced into the shared binomial variants of
// util/spa_cgf:
//
//   class 1  unrelated outliers          G_i ~ Binomial(2, MAF), independent
//   class 2  two-subject families        a four-point-mass correlated mixture
//   class 3  three-or-more-subject       moments of a general discrete law
//            families                    tabulated over a Chow-Liu table
//
// plus the analytic non-outlier block, a Gaussian approximation contributing
// `mean·t + var·t²/2` to K.
//
// Only class 1 is the shared binomial CGF, and it is delegated to
// spa_cgf::binomUniform — which supplies D1's exact cancellation-free K'', N1's
// log1p/expm1 K, and the AVX2 / AVX-512 tiers.  Classes 2 and 3 are handled
// here on their own terms.
//
// ══════════════════════════════════════════════════════════════════════
// Why D1's rewrite does not generalize, and what replaces it
// ══════════════════════════════════════════════════════════════════════
//
// 01_findings.md D1 prescribes K'' = h·r²·e·u/alpha², reached by substituting
// alpha - e = u = 1 - p.  That substitution exists only because the binomial
// MGF factorizes as alpha^h.  It is valid for class 1 and for nothing else
// here: the previous implementation nevertheless accumulated
//
//     CGF2 = sum_k MGF2_k/MGF0_k  -  sum_k (MGF1_k/MGF0_k)^2  +  var
//
// with the two sums taken GLOBALLY over all three classes and differenced once
// at the end, so the cancellation ran across the whole outlier set rather than
// per element.  Applying D1's substitution at that aggregation point would
// apply it to two classes it does not cover.
//
// For a general discrete law the per-element term is a variance under the
// exponentially tilted measure,
//
//     ptilde_k = w_k·e^{t·v_k} / sum_j w_j·e^{t·v_j}
//     K'_elem  = sum_k ptilde_k·v_k                    =: m
//     K''_elem = sum_k ptilde_k·(v_k - m)^2
//
// and the stable form is the MEAN-CENTRED weighted sum on the right, not a
// difference of raw moments.  Every summand is a non-negative weight times a
// square, so K''_elem cannot be driven negative or to noise by cancellation;
// the only subtraction is v_k - m, a difference of two O(|v|) quantities whose
// result is then squared.  This is the class-2 and class-3 treatment below, and
// it is also what makes the accumulation per element rather than global — the
// pattern 01_findings.md D1 identifies at wtcoxg.cpp:80 as the one to copy.
//
// ══════════════════════════════════════════════════════════════════════
// class 2 — the two-subject block, read as a distribution
// ══════════════════════════════════════════════════════════════════════
//
// The previous code built, per family i and per IBD-mixture component j,
//
//     tj  = (1 - rho_j)·MAF·(1-MAF)
//     g0  = e^{tR1}·tj + e^{tR2}·tj + e^{t(R1+R2)}·(MAF - tj) - tj + (1-MAF)
//
// That is the MGF of a four-point law, with the last two terms being the
// value-0 mass written as a difference:
//
//     value    R1        R2        R1+R2          0
//     weight   tj        tj        MAF - tj       1 - MAF - tj
//
// The four weights sum to one (tj + tj + MAF - tj + 1 - MAF - tj = 1), which is
// the check that the reading is right.  Grouping the constant as the single
// weight w0 = (1-MAF) - tj makes g0 a sum of four NON-NEGATIVE terms and
// removes the `- tj` subtraction, which was a cancelling difference in its own
// right at t = 0, where the three exponential terms sum to tj + MAF and the
// whole expression collapses to 1.
//
// Non-negativity of all four weights requires rho_j in [0, 1]: then
// tj <= MAF·(1-MAF), so MAF - tj >= MAF² >= 0 and 1 - MAF - tj >= (1-MAF)² > 0.
// rho comes from the `--cal-pairwise-ibd` file via grm_null.cpp, which does not
// validate it, so the bound is conditional on well-formed input rather than
// proven in code.  A rho outside [0, 1] makes two weights negative, K'' can
// then be non-positive, and the solver reports GUARD_CURV with P = NA instead
// of the silent NaN the unguarded predecessor emitted.
//
// K for this class uses the same log1p identity util/spa_cgf uses for class 1.
// Because the weights sum to one exactly (symbolically),
//
//     g0 = 1 + w1·expm1(a1) + w2·expm1(a2) + w3·expm1(a3)  =:  1 + delta
//
// with a1 = t·R1, a2 = t·R2, a3 = a1 + a2, so K_elem = log1p(delta) never
// materializes 1 + delta and keeps full relative precision at small t.  The
// usual danger of a log1p form — delta approaching -1 — cannot arise here:
// g0 >= w0 >= (1-MAF)² >= 1/4 for MAF <= 1/2, so delta >= -3/4.  The value-0
// support point is what pins the bound.
//
// ══════════════════════════════════════════════════════════════════════
// class 3 — the Chow-Liu block
// ══════════════════════════════════════════════════════════════════════
//
// Support `stand_S` (3^n standardized genotype sums for an n-member family,
// n <= MAX_NUM_IN_FAM = 5, so at most 243 points) with probabilities `arr_prob`
// interpolated per marker between two MAF columns of the family's Chow-Liu
// table.  There is no closed form of any shape here, so:
//
//   * K'' is the mean-centred weighted sum, as above.  It needs the tilted
//     weights twice, so they are written to a caller-supplied scratch buffer on
//     the first pass; the arrays are at most 243 doubles and stay in L1.
//   * K keeps log(g0) with g0 the same sum of non-negative terms the previous
//     code formed.  The log1p treatment is deliberately NOT applied: it would
//     require sum_j arr_prob_j == 1 to machine precision after the two-column
//     interpolation, which nothing in the code establishes, and the bound that
//     keeps delta away from -1 for class 2 has no analogue here.  g0 is
//     instead guarded: a non-positive or non-finite g0 propagates a non-finite
//     K, which spa::bnTailLog reports as FALLBACK_NONFINITE rather than handing
//     `std::log` a non-positive argument as the predecessor did.
//
// ══════════════════════════════════════════════════════════════════════
// Why classes 2 and 3 are NOT domain-clamped the way class 1 is
// ══════════════════════════════════════════════════════════════════════
//
// util/spa_cgf clamps t·r to +/-708 before `exp`, which keeps class 1 finite over
// the whole real line.  Classes 2 and 3 deliberately do not do this, and the
// asymmetry is a decision rather than an oversight.
//
// A clamp applied per support point of a MIXTURE changes the mixture.  If
// t·R1 = 800 and t·R2 = 750 the true weight ratio between those two support
// points is e^{50}; clamping both to 708 makes it 1.  The result would be finite,
// plausible and wrong, and nothing downstream could detect it.  Class 1 is exempt
// because each subject's binomial factor is a separate factor of the MGF, not a
// competing term in one sum, so saturating it costs a bounded amount in a
// quantity below 1e-307.
//
// What happens instead: past |t·(R1+R2)| ~ 709 the class-2 exponential overflows
// and K'/K'' become non-finite.  spa::solveSaddlepoint tests isfinite on every
// evaluation — it halves the probe step during bracket expansion, falls back to
// bisection inside the Newton loop, and retreats toward the origin from a bad
// initial abscissa — so a non-finite region far out is walked around, and if no
// usable abscissa exists at all the result is Status::FallbackNonFinite: the
// saddlepoint is discarded and the normal tail reported under that code.  A
// refusal is the correct behaviour there; a distorted finite value is not.
//
// The regime is in any case unreachable in a converged solve: the initial
// abscissa is capped at 1.2, K' is monotone, and the analytic non-outlier block
// contributes var·t, so the root for any representable score lies far inside
// |t·r| < 700.  Over the 37 512 marker tests of examples/baseline.sh no SPAGRM
// solve reports FALLBACK_MAXITER or FALLBACK_NONFINITE.
//
// ══════════════════════════════════════════════════════════════════════
// Cost structure
// ══════════════════════════════════════════════════════════════════════
//
// `k12` is the Newton-loop payload and computes no logarithm at all; `kFull`
// adds K and runs once per tail, at the converged root (util/spa.hpp § 4a).
// The predecessor walked the workspace three times per tail in getProbSpa and
// FOUR times per Newton iteration in fastGetRoot, materializing MGF0, MGF1,
// MGF2 and a temporary quotient vector — 8 x nOutlier doubles of scratch per
// clone.  Fusing to a single pass per class removes all four vectors (P6);
// the only scratch left is the class-3 tilted-weight buffer.

#pragma once

#include <array>
#include <cmath>
#include <limits>
#include <vector>

#include "util/spa.hpp"
#include "util/spa_cgf.hpp"

namespace spagrm_cgf {

// ══════════════════════════════════════════════════════════════════════
// Per-marker data for one three-or-more-subject family
// ══════════════════════════════════════════════════════════════════════
//
// `stand_S` is marker-invariant; `arr_prob` is refreshed per marker from the
// family's Chow-Liu table at the marker's MAF.  Formerly
// nsSPAGRM::UpdatedThreeSubj, which remains as an alias.
struct ThreeSubjTable {
    std::vector<double> stand_S;
    std::vector<double> arr_prob;
};

// ══════════════════════════════════════════════════════════════════════
// Everything the CGF reads, gathered once per marker
// ══════════════════════════════════════════════════════════════════════
//
// Pointers are non-owning and must outlive the solve.  `scratch` must have room
// for the largest `stand_S` across the class-3 families; it may be null when
// there are none.
struct Context {
    // class 1 — unrelated outliers
    const double *outlierResid = nullptr;
    int nOutlier = 0;

    // class 2 — two-subject families; twoResid[i] = {R1, R2}, twoRho[i] = rho_j
    const std::array<double, 2> *twoResid = nullptr;
    const std::vector<double> *twoRho = nullptr;
    int nTwo = 0;

    // class 3 — Chow-Liu families
    const ThreeSubjTable *three = nullptr;
    int nThree = 0;
    double *scratch = nullptr;

    // shared and analytic
    double MAF = 0.0;
    double mean = 0.0;   // 2·MAF·sum_R_nonOutlier
    double var = 0.0;    // 2·MAF·(1-MAF)·R_GRM_R_nonOutlier
};

// ══════════════════════════════════════════════════════════════════════
// Newton-loop payload:  {K'(t) - s, K''(t)}
// ══════════════════════════════════════════════════════════════════════
//
// The residual is assembled here, in the owning method's association order, as
// spa::K12 requires.

inline spa::K12 k12(double t, const Context &c, double s) noexcept {
    double K1 = 0.0;
    double K2 = 0.0;

    // ── class 1 — the shared binomial CGF, SIMD-dispatched ──────────
    if (c.nOutlier > 0) {
        const spa_cgf::Cgf12 d =
            spa_cgf::binomUniformK12(t, c.outlierResid, c.nOutlier, c.MAF);
        K1 += d.K1;
        K2 += d.K2;
    }

    // ── class 2 — four-point-mass mixture, mean-centred K'' ─────────
    if (c.nTwo > 0) {
        const double maf01 = c.MAF * (1.0 - c.MAF);
        const double oneMinusMAF = 1.0 - c.MAF;
        for (int i = 0; i < c.nTwo; ++i) {
            const double R1 = c.twoResid[i][0];
            const double R2 = c.twoResid[i][1];
            const double Rs = R1 + R2;
            const double e1 = std::exp(t * R1);
            const double e2 = std::exp(t * R2);
            const double e3 = e1 * e2;                 // e^{t·(R1+R2)}

            const std::vector<double> &rho = c.twoRho[i];
            const size_t ki = rho.size();
            for (size_t j = 0; j < ki; ++j) {
                const double tj = (1.0 - rho[j]) * maf01;
                const double w3 = c.MAF - tj;
                const double w0 = oneMinusMAF - tj;

                const double E1 = tj * e1;
                const double E2 = tj * e2;
                const double E3 = w3 * e3;
                // g0 as a sum of four non-negative terms — no `- tj`.
                const double g0 = (E1 + E2) + (E3 + w0);
                const double inv = 1.0 / g0;

                const double m = (R1 * E1 + R2 * E2 + Rs * E3) * inv;
                const double d1 = R1 - m;
                const double d2 = R2 - m;
                const double d3 = Rs - m;
                K1 += m;
                K2 += ((E1 * d1 * d1 + E2 * d2 * d2) +
                       (E3 * d3 * d3 + w0 * m * m)) * inv;
            }
        }
    }

    // ── class 3 — tabulated discrete law, mean-centred K'' ──────────
    for (int i = 0; i < c.nThree; ++i) {
        const double *ss = c.three[i].stand_S.data();
        const double *ap = c.three[i].arr_prob.data();
        const size_t ns = c.three[i].stand_S.size();
        double *buf = c.scratch;

        double g0 = 0.0, g1 = 0.0;
        for (size_t j = 0; j < ns; ++j) {
            const double e = std::exp(t * ss[j]) * ap[j];
            buf[j] = e;
            g0 += e;
            g1 += ss[j] * e;
        }
        const double inv = 1.0 / g0;
        const double m = g1 * inv;
        double v = 0.0;
        for (size_t j = 0; j < ns; ++j) {
            const double d = ss[j] - m;
            v += buf[j] * d * d;
        }
        K1 += m;
        K2 += v * inv;
    }

    return spa::K12{K1 + c.mean + c.var * t - s, K2 + c.var};
}

// ══════════════════════════════════════════════════════════════════════
// Terminal payload:  {K(t), K'(t) - s, K''(t)}
// ══════════════════════════════════════════════════════════════════════
//
// Invoked once per tail, at the converged root.  Adds the logarithms: class 1
// through spa_cgf's log1p/expm1 form, class 2 through the log1p identity
// derived in the header comment, class 3 through a guarded log.

inline spa::K012 kFull(double t, const Context &c, double s) noexcept {
    double K0 = 0.0;
    double K1 = 0.0;
    double K2 = 0.0;

    if (c.nOutlier > 0) {
        const spa_cgf::Cgf012 d =
            spa_cgf::binomUniformKFull(t, c.outlierResid, c.nOutlier, c.MAF);
        K0 += d.K0;
        K1 += d.K1;
        K2 += d.K2;
    }

    if (c.nTwo > 0) {
        const double maf01 = c.MAF * (1.0 - c.MAF);
        const double oneMinusMAF = 1.0 - c.MAF;
        for (int i = 0; i < c.nTwo; ++i) {
            const double R1 = c.twoResid[i][0];
            const double R2 = c.twoResid[i][1];
            const double Rs = R1 + R2;
            const double a1 = t * R1;
            const double a2 = t * R2;
            const double e1 = std::exp(a1);
            const double e2 = std::exp(a2);
            const double e3 = e1 * e2;
            // expm1 for K only, never to reconstruct the exponential itself
            // (util/spa_cgf.hpp, the N1 caveat on logAlpha).
            const double q1 = std::expm1(a1);
            const double q2 = std::expm1(a2);
            const double q3 = std::expm1(a1 + a2);

            const std::vector<double> &rho = c.twoRho[i];
            const size_t ki = rho.size();
            for (size_t j = 0; j < ki; ++j) {
                const double tj = (1.0 - rho[j]) * maf01;
                const double w3 = c.MAF - tj;
                const double w0 = oneMinusMAF - tj;

                const double E1 = tj * e1;
                const double E2 = tj * e2;
                const double E3 = w3 * e3;
                const double g0 = (E1 + E2) + (E3 + w0);
                const double inv = 1.0 / g0;

                const double m = (R1 * E1 + R2 * E2 + Rs * E3) * inv;
                const double d1 = R1 - m;
                const double d2 = R2 - m;
                const double d3 = Rs - m;
                K1 += m;
                K2 += ((E1 * d1 * d1 + E2 * d2 * d2) +
                       (E3 * d3 * d3 + w0 * m * m)) * inv;
                K0 += std::log1p(tj * q1 + tj * q2 + w3 * q3);
            }
        }
    }

    for (int i = 0; i < c.nThree; ++i) {
        const double *ss = c.three[i].stand_S.data();
        const double *ap = c.three[i].arr_prob.data();
        const size_t ns = c.three[i].stand_S.size();
        double *buf = c.scratch;

        double g0 = 0.0, g1 = 0.0;
        for (size_t j = 0; j < ns; ++j) {
            const double e = std::exp(t * ss[j]) * ap[j];
            buf[j] = e;
            g0 += e;
            g1 += ss[j] * e;
        }
        const double inv = 1.0 / g0;
        const double m = g1 * inv;
        double v = 0.0;
        for (size_t j = 0; j < ns; ++j) {
            const double d = ss[j] - m;
            v += buf[j] * d * d;
        }
        K1 += m;
        K2 += v * inv;
        // The guard 01_findings.md P6 asks for at spagrm.cpp:231.  A
        // non-positive or non-finite g0 yields a non-finite K, which
        // spa::bnTailLog turns into NaN plus Status::FallbackNonFinite; the predecessor
        // handed std::log a non-positive argument with no diagnostic.
        K0 += (g0 > 0.0) ? std::log(g0)
                         : -std::numeric_limits<double>::infinity();
    }

    return spa::K012{K0 + c.mean * t + 0.5 * c.var * t * t,
                     K1 + c.mean + c.var * t - s,
                     K2 + c.var};
}

// Largest class-3 table, i.e. the size `Context::scratch` must have.
inline size_t scratchSize(const ThreeSubjTable *three, int nThree) noexcept {
    size_t mx = 0;
    for (int i = 0; i < nThree; ++i)
        if (three[i].stand_S.size() > mx) mx = three[i].stand_S.size();
    return mx;
}

}  // namespace spagrm_cgf
