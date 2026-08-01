// math_helper.cpp — Out-of-line implementations for math_helper.hpp
#include "math_helper.hpp"

#include <boost/math/special_functions/beta.hpp>

namespace math {

namespace {

constexpr double kLn2      = 0.69314718055994530942;  // ln 2
constexpr double kLnPi     = 1.14472988584940017414;  // ln π
constexpr double kLn2OverPi = kLn2 - kLnPi;           // ln(2/π)

// 20-point Gauss-Legendre nodes / weights on [-1, 1] (positive half; each
// value is shared by ±x).  Source: standard tables, e.g. Abramowitz & Stegun
// 25.4.30.  One copy for the whole file: `pmvnorm2dHalfRect` and
// `pmvnorm2dHalfRectLog` must use the same rule, since the claim that the
// log-domain routine is the same quadrature carried through a log-sum-exp
// rests on it.  (`bvnCdf` above keeps its own table, whose entries are the
// negatives of these; it consumes them as `is * x20[i]` over is = ±1, so the
// set of abscissae is identical but the summation ORDER is not, and changing
// it would move that function's last bit for no benefit.)
constexpr double kGl20Node[10] = {
    0.0765265211334973, 0.2277858511416451, 0.3737060887154195,
    0.5108670019508271, 0.6360536807265150, 0.7463319064601508,
    0.8391169718222188, 0.9122344282513259, 0.9639719272779138,
    0.9931285991850949
};
constexpr double kGl20Weight[10] = {
    0.1527533871307258, 0.1491729864726037, 0.1420961093183820,
    0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
    0.0832767415767048, 0.0626720483341091, 0.0406014298003869,
    0.0176140071391521
};

} // namespace

// ──────────────────────────────────────────────────────────────────────
// § 1b  Log-domain p-value inversion and evaluation
// ──────────────────────────────────────────────────────────────────────

namespace {

// Below this L the linear argument 0.5·10^(−L) is an ordinary double and
// Boost's quantile is used; above it the Newton iteration on ln P runs.  The
// seam is continuous to a few ULP (measured: |Δz|/z ≤ 3e-16 over
// L ∈ [0.9, 1.1]), because both sides are accurate there — the same
// overlap-in-a-regime-where-both-are-right argument that fixes pnormLog's own
// branch at |t| = 37.
constexpr double kZNewtonMinL = 1.0;
constexpr int    kZMaxIter    = 64;

} // namespace

double zFromNegLog10P(
    double negLog10p,
    double zNormForSign
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    if (std::isnan(negLog10p) || std::isnan(zNormForSign)) return nan;
    const double sign = (zNormForSign >= 0.0) ? 1.0 : -1.0;

    // P ≥ 1 gives z = 0 and P = 0 gives z = ∞; both are exact and both lie
    // outside the iteration's domain.
    if (!(negLog10p > 0.0)) return sign * 0.0;
    if (std::isinf(negLog10p)) return sign * negLog10p;
    if (negLog10p < kZNewtonMinL)
        return sign * qnorm(0.5 * std::pow(10.0, -negLog10p), 0.0, 1.0,
                            /*lower_tail=*/false);

    const double target = -negLog10p * kLn10;          // ln P
    const double A      = 2.0 * negLog10p * kLn10;
    // Analytic inversion of the Mills asymptote ln P = ln 2 − z²/2 − ln z
    // − ln√(2π), with ln z ≈ ½ ln A on the right-hand side.
    double z = std::sqrt(std::fmax(A - std::log(A) + kLn2OverPi, 1.0));

    // Relative criterion, no absolute constant: see the header comment and
    // 04_validation.md §2 invariant 3.
    const double tol = 4.0 * std::numeric_limits<double>::epsilon() * std::fabs(target);
    for (int it = 0; it < kZMaxIter; ++it) {
        const double lnPhi = pnormLog(-z, 0.0, 1.0, /*lower_tail=*/true);
        const double f     = kLn2 + lnPhi - target;
        if (std::fabs(f) <= tol) break;
        // d/dz ln(2Φ(−z)) = −φ(z)/Φ(−z), formed in the log domain so that
        // neither φ nor Φ underflows on its own.
        const double deriv = -std::exp((-0.5 * z * z - kLogSqrt2Pi) - lnPhi);
        if (!(deriv < 0.0)) break;
        double zNext = z - f / deriv;
        if (!(zNext > 0.0)) zNext = 0.5 * z;
        if (!std::isfinite(zNext) || zNext == z) break;
        z = zNext;
    }
    return sign * z;
}

double ptLog(
    double t,
    double df
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    if (std::isnan(t) || std::isnan(df) || !(df > 0.0)) return nan;

    const double at = std::fabs(t);
    if (!(at > 0.0)) return 0.0;      // P(|T| > 0) = 1
    if (std::isinf(at)) return -inf;

    const double a     = 0.5 * df;
    const double b     = 0.5;
    const double denom = df + at * at;
    const double x     = df / denom;
    // ln x and ln(1−x) built from the closed forms rather than from x, so
    // that neither loses digits as x → 1 (small t, large df).
    const double lnx   = std::log(df) - std::log(denom);
    const double ln1mx = 2.0 * std::log(at) - std::log(denom);
    const double lnBeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);
    // ln of the leading factor x^a (1−x)^b / (a·B(a,b)) of DLMF 8.17.7.  The
    // ₂F₁ that multiplies it is ≥ 1, so this is a lower bound on ln I_x and
    // decides whether Boost's ibeta can represent the answer at all.
    const double lnLead = a * lnx + b * ln1mx - std::log(a) - lnBeta;

    if (lnLead >= -690.0) {
        const double ib = boost::math::ibeta(a, b, x);
        if (ib > 0.0) return std::log(ib);
    }

    // ₂F₁(a+b, 1; a+1; x) = Σ_{n≥0} [(a+b)_n/(a+1)_n]·xⁿ, all terms positive.
    // The ratio of successive terms tends to x, so the term count is data
    // dependent (about 300 at df = 10 000) and the stopping rule has to be
    // relative rather than a fixed truncation.
    constexpr long kMaxTerms = 100000;
    double term = 1.0, series = 1.0;
    for (long n = 0; n < kMaxTerms; ++n) {
        const double k = static_cast<double>(n);
        term *= x * (a + b + k) / (a + 1.0 + k);
        series += term;
        if (term < 1e-18 * series) break;
    }
    return lnLead + std::log(series);
}

double cauchyCombineLog10(
    const double *negLog10p,
    int n
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    if (n <= 0) return nan;

    // The p ≥ 1 → 0.999 clamp of `cauchyCombine`, expressed on the L scale.
    const double clampL = -std::log10(0.999);

    int nValid = 0;
    double maxL = -inf;
    for (int i = 0; i < n; ++i) {
        if (std::isnan(negLog10p[i])) continue;
        ++nValid;
        const double Li = std::fmax(negLog10p[i], clampL);
        if (Li > maxL) maxL = Li;
    }
    if (nValid == 0) return nan;
    if (std::isinf(maxL)) return inf;   // some p is 0; so is the combination

    // T = 10^maxL · A, never formed.  Terms with pᵢ ≤ 1e-15 use
    // cot(π p) = 1/(π p) − π p/3 − …, whose dropped correction is 1e-29
    // relative to the retained one.
    const double scale = std::pow(10.0, -maxL);
    double A = 0.0;
    for (int i = 0; i < n; ++i) {
        if (std::isnan(negLog10p[i])) continue;
        const double Li = std::fmax(negLog10p[i], clampL);
        if (Li >= 15.0) {
            A += std::pow(10.0, Li - maxL) / M_PI;
        } else {
            const double p = std::pow(10.0, -Li);
            A += (std::cos(M_PI * p) / std::sin(M_PI * p)) * scale;
        }
    }
    A /= static_cast<double>(nValid);

    // Upper tail of the standard Cauchy at T, returned as −log10.
    if (!(A > 0.0)) {
        // T ≤ 0, and |T| is bounded by cot(π·0.999) = 318.3 because a single
        // Lᵢ ≥ 15 would already have made A positive.  Written through log1p
        // on the tail's distance from 1 so that L keeps full RELATIVE accuracy
        // as p → 1 (the direct −log10(0.5 − atan T/π) loses it there).
        const double negT = -(A * std::pow(10.0, maxL));
        if (!(negT > 0.0)) return kLn2 / kLn10;                 // T = 0 → p = ½
        return -std::log1p(-std::atan(1.0 / negT) / M_PI) / kLn10;
    }

    const double lnT = maxL * kLn10 + std::log(A);
    if (lnT <= 1.0)
        return -std::log10(0.5 - std::atan(std::exp(lnT)) / M_PI);
    if (lnT < 350.0) {
        // p = atan(1/T)/π.  The direct form is accurate over the whole of
        // this interval; the truncated log(atan u)/u series an earlier draft
        // used is wrong by 3.5e-5 at lnT = 1 and only reaches double
        // precision at lnT ≥ 5 (01_numerics §2.2).
        return -(std::log(std::atan(std::exp(-lnT))) - kLnPi) / kLn10;
    }
    // exp(−lnT) underflows past lnT = 708; atan(u) = u to well within double
    // precision once lnT > 20, so the asymptote is exact here.
    return (lnT + kLnPi) / kLn10;
}

// § 2  Bivariate normal CDF  (Genz 2004)
//
// Internal helper.  Returns Φ₂(dh, dk; r) = P(X ≤ dh, Y ≤ dk) for
// (X, Y) ~ BVN(0, 0, 1, 1, r).  Used by pmvnorm2dHalfRect when one Y
// bound is infinite (the half-infinite case reduces to a single Φ₂
// evaluation with no inclusion-exclusion), and when |ρ| has been clamped
// to ±1 in the finite-rectangle case.
//
// ══════════════════════════════════════════════════════════════════════
// The |r| ≥ 0.925 branch was wrong.  What was wrong, and how it showed.
// ══════════════════════════════════════════════════════════════════════
//
// Genz's routine is written for the UPPER tail,
//     BVNU(H, K, r) = P(X > H, Y > K),
// and the CDF is recovered as Φ₂(dh, dk; r) = BVNU(−dh, −dk, r).  The
// |r| < 0.925 branch below implements Plackett's identity directly and
// is correct; the |r| ≥ 0.925 branch is Genz's asymptotic expansion about
// the degenerate line, and its previous transcription departed from the
// published algorithm in three places:
//
//   1. The quadrature term read
//          exp(−HK(1−r) / (2(1 + (±xᵢ+1)·(a/2)²/AS))) − 1 − C·XS·(1 + D·XS)
//      where Genz has
//          exp(−HK·XS / (2(1+RS)²)) / RS − (1 + C·XS·(1 + D·XS)),
//      with XS = (a/2·(±xᵢ+1))² and RS = √(1 − XS).  Both the exponent and
//      the missing 1/RS matter: the series then does not converge to the
//      expansion it is truncating and the result is unbounded.  Measured
//      against an independent Simpson reduction of
//      ∫_{−∞}^{h} φ(x)·Φ((k − r x)/√(1−r²)) dx over a 784-point grid, the
//      previous code was in error by up to 9.9e-01 in absolute probability
//      and returned a value outside [0, 1] at 56 of those 784 points.
//
//   1b. The Mills-ratio term read `PHI(+B/A)` where Genz has `PHI(−B/A)`.
//      B = √BS ≥ 0 and A > 0, so the transcribed factor is at least ½ where
//      the correct one is at most ½; at (h, k, r) = (−4, −0.392, 0.925) that
//      alone moved the answer by 1.0 in absolute probability.
//
//   2. The r > 0 terminal read `bvn += Φ(h) + Φ(k) − 1 + min(Φ(h), Φ(k))`.
//      Genz's line is `BVN = BVN + PHI(−max(H, K))`, which under the
//      substitution H = −dh, K = −dk is `bvn += min(Φ(dh), Φ(dk))`.  The
//      extra `Φ(h) + Φ(k) − 1` was spurious.
//
//   3. The r < 0 terminal negated the series twice: it computed
//      `bvn = −bvn` and then `bvn = Φ(h) − Φ(kk) − bvn`, i.e.
//      +series + (Φ(h) − Φ(kk)), where Genz has −series + max(0, ·).
//
//   4. The quadrature used the 6-point rule at every |r|.  Genz selects
//      6 / 12 / 20 points for |r| < 0.3 / < 0.75 / else; the |r| ≥ 0.925
//      branch is always in the last bin.
//
// Observable consequence in GRAB.  `pmvnorm2dHalfRect` clamps ρ into
// [−1, 1] and routes |ρ| ≥ 1 − 1e-12 through this function; WtCoxG's
// Branch B builds ρ from a variance recovered by inverting an SPA p-value
// through qchisq, which is not constrained to keep the 2×2 covariance
// positive semi-definite, so ρ > 1 does occur (max 1.073 on the bundled
// fixture).  With defect (2), Φ₂(h, b; 1) − Φ₂(h, a; 1) came back as
// 0.4557 where the exact degenerate answer is 0.1505 — larger than the
// S_bat marginal 0.3052 that bounds it — and the conditional p-value
// 2·(TPR·p₁ + (1−TPR)·p₀)/p_deno exceeded 1, reaching 2.9866.  That is the
// origin of the 166 markers with P_EXT > 1 in the pinned baseline.
//
// The |r| < 0.925 branch is left byte-for-byte as it was: it is correct
// (agreeing with the Simpson reference to ~1e-7, the accuracy of the
// 6-point rule Genz specifies for that regime), and every call in the
// bundled fixture that is not clamped to |ρ| = 1 lands there, so keeping it
// untouched confines the numeric change to the calls that were wrong.
//
// After the repair the same 784-point grid gives a worst absolute error of
// 6.5e-13 and no return outside [0, 1].  See tests/wtcoxg_cgf_test.cpp.
double bvnCdf(
    double dh,
    double dk,
    double r
) {
    // Independent case
    if (std::abs(r) < 1e-15) return 0.5 * std::erfc(-dh / std::sqrt(2.0)) * 0.5 * std::erfc(-dk / std::sqrt(2.0));

    static constexpr double w6[] = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910};
    static constexpr double x6[] = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969};

    const double h = dh, k = dk, hk = h * k;
    double bvn = 0.0;

    if (std::abs(r) < 0.925) {
        const double hs = (h * h + k * k) / 2.0;
        const double asr = std::asin(r);
        for (int i = 0; i < 3; ++i) {
            for (int is = -1; is <= 1; is += 2) {
                double sn = std::sin(asr * (is * x6[i] + 1.0) / 2.0);
                bvn += w6[i] * std::exp((sn * hk - hs) / (1.0 - sn * sn));
            }
        }
        bvn *= asr / (4.0 * M_PI);
        bvn += 0.5 * std::erfc(-h / std::sqrt(2.0)) * 0.5 * std::erfc(-k / std::sqrt(2.0));
    } else {
        // Genz's asymptotic branch.  It is written in terms of BS = (H − K)²
        // and HK = H·K only, and both are invariant under the joint sign flip
        // H = −dh, K = −dk that turns the upper tail into the CDF, so the
        // series may be evaluated with (dh, dk) directly.  Only the terminal
        // needs the substitution spelled out.
        //
        // Genz specifies 20 points (his NG = 3 table) for every |r| in this
        // branch; the same nodes and weights already serve the Gauss-Legendre
        // rectangle integration in pmvnorm2dHalfRect below.
        static constexpr double w20[10] = {
            0.1527533871307258, 0.1491729864726037, 0.1420961093183820,
            0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
            0.0832767415767048, 0.0626720483341091, 0.0406014298003869,
            0.0176140071391521
        };
        static constexpr double x20[10] = {
            -0.0765265211334973, -0.2277858511416451, -0.3737060887154195,
            -0.5108670019508271, -0.6360536807265150, -0.7463319064601508,
            -0.8391169718222188, -0.9122344282513259, -0.9639719272779138,
            -0.9931285991850949
        };

        double kk = k;
        double hkk = hk;
        if (r < 0.0) {
            kk = -kk;
            hkk = -hkk;
        }

        if (std::abs(r) < 1.0) {
            const double as_ = (1.0 - r) * (1.0 + r);
            double a = std::sqrt(as_);
            const double bs = (h - kk) * (h - kk);
            const double c = (4.0 - hkk) / 8.0;
            const double d = (12.0 - hkk) / 16.0;
            const double asr = -(bs / as_ + hkk) / 2.0;
            if (asr > -100.0)
                bvn = a * std::exp(asr) * (1.0 - c * (bs - as_) * (1.0 - d * bs / 5.0) / 3.0 + c * d * as_ * as_ / 5.0);
            if (-hkk < 100.0) {
                // Genz: BVN -= exp(-HK/2)*sqrt(2*pi)*PHI(-B/A)*B*(...).
                // The argument of PHI is NEGATIVE B/A.  The previous
                // transcription wrote 0.5*erfc(-b/(sqrt(2)*a)), which is
                // Phi(+b/a); since b = sqrt(bs) >= 0 and a > 0 that is >= 1/2
                // where the correct factor is <= 1/2, so the subtracted term
                // was too large by up to a factor of the Mills ratio.
                const double b = std::sqrt(bs);
                const double phiNegBA = 0.5 * std::erfc(b / (std::sqrt(2.0) * a));
                bvn -= std::exp(-hkk / 2.0) * std::sqrt(2.0 * M_PI) * phiNegBA * b *
                       (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
            }
            a /= 2.0;
            for (int i = 0; i < 10; ++i) {
                for (int is = -1; is <= 1; is += 2) {
                    const double t = a * (is * x20[i] + 1.0);
                    const double xs = t * t;
                    const double rs = std::sqrt(1.0 - xs);
                    const double asr2 = -(bs / xs + hkk) / 2.0;
                    if (asr2 > -100.0)
                        bvn += a * w20[i] * std::exp(asr2) *
                               (std::exp(-hkk * xs / (2.0 * (1.0 + rs) * (1.0 + rs))) / rs -
                                (1.0 + c * xs * (1.0 + d * xs)));
                }
            }
            bvn = -bvn / (2.0 * M_PI);
        }

        // Terminal, with H = −dh and (after the r < 0 flip) K = −kk:
        //   r > 0 :  BVN += PHI(−max(H, K)) = PHI(min(dh, dk))
        //                                   = min(Φ(dh), Φ(dk))
        //   r < 0 :  BVN  = −BVN + max(0, PHI(−H) − PHI(−K))
        //                 = −BVN + max(0, Φ(dh) − Φ(kk))
        // At |r| = 1 the series above is skipped and these are the exact
        // degenerate answers, min(Φ(h), Φ(k)) and max(0, Φ(h) + Φ(k) − 1).
        const double phih  = 0.5 * std::erfc(-h / std::sqrt(2.0));
        const double phikk = 0.5 * std::erfc(-kk / std::sqrt(2.0));
        if (r > 0.0) bvn += std::min(phih, phikk);
        else         bvn = -bvn + std::max(0.0, phih - phikk);

    }
    // A probability, on both branches.  Plackett's identity under a 6-point
    // rule (the |r| < 0.925 arm) is accurate to about 1e-7 and can return a
    // value a few times 1e-8 below zero where the true probability is zero;
    // every call site already wrapped the result in std::clamp for exactly
    // that reason, so hoisting the clamp here changes no caller's answer and
    // makes the function's own contract true.
    if (!(bvn > 0.0)) bvn = std::isnan(bvn) ? bvn : 0.0;
    if (bvn > 1.0) bvn = 1.0;
    return bvn;
}

double pmvnorm2dHalfRect(
    double s_hi,
    double sb_lo,
    double sb_hi,
    double var1,
    double cov12,
    double var2
) {
    if (var1 <= 0.0 || var2 <= 0.0) return 0.0;

    const double sd1 = std::sqrt(var1);
    const double sd2 = std::sqrt(var2);
    double rho = cov12 / (sd1 * sd2);

    // ── |ρ| > 1 means there is no such bivariate normal ──────────────
    //
    // The guard here used to be a bare clamp with the comment "cov12 may
    // marginally exceed sd1·sd2 due to round-off".  That is a correct thing to
    // absorb — the three arguments arrive from separate accumulations, so a
    // few ULP either way is expected — but it is not what was happening.
    //
    // Measured over every WtCoxG and LEAF invocation of examples/baseline.sh
    // (22 776 calls): |ρ| reached or exceeded 1 in 177 of them, and in all 177
    // it exceeded 1 by between 0.14 % and 172 % (median 7.3 %).  The largest
    // |ρ| anywhere below the threshold was 0.9956, so the two populations do
    // not even touch.  These are not rounding artefacts.
    //
    // They come from WtCoxG's Branch B, where the covariance and the second
    // variance are built with different σ² dependence:
    //     cov  = w1·(R−bm)·2μ(1−μ) + 2·b·sumR·(var_mu_ext + σ²)
    //     var2 = w1Sq·2μ(1−μ) + var_mu_ext + σ²
    // so as σ² grows the covariance grows like σ² while sd2 grows like σ, and
    // ρ grows without bound.  On the fixture ρ > 1 occurs exactly where
    // σ²/var_Sbat is large (median 24, against 8 elsewhere) — the degenerate,
    // weakly identified regime.  Note that ρ does not depend on var1 at all:
    // the sqrt(var_S/var_int) rescaling the caller applies to cov12 cancels
    // against sd1 exactly, so this is a property of WtCoxG's model, not of the
    // variance it recovers by inverting the saddlepoint p-value.
    //
    // Clamping |ρ| to 1 does not rescue such a marker, it replaces it: at
    // ρ = 1 the two coordinates are the same random variable, so the joint
    // event {X ≤ h} ∩ {a ≤ Y ≤ b} is either the whole rectangle or empty.
    // Three markers of the bundled fixture land on the empty side, and the
    // conditional p-value 2·(TPR·p₁+(1−TPR)·p₀)/p_deno then collapses to
    // 2·(1−TPR)·p₀/p_deno ≈ 8e-14 — a fabricated genome-wide-significant hit
    // manufactured out of TPR differing from 1 in the fourteenth digit.  The
    // previous code reported the same three markers as p ≈ 2, which is not a
    // probability either.  Neither number should be produced.
    //
    // So: round-off is still absorbed, and a genuinely indefinite covariance
    // is reported.  NaN propagates through the caller's assembly to an NA
    // p-value with SPA_STATUS = NONFINITE, which is the honest answer — the
    // integral being asked for does not exist.  Whether WtCoxG's Branch B
    // should build a positive semi-definite pair in the first place is a
    // modelling decision for whoever owns that formula, and is deliberately
    // NOT taken here (compare D4, which this file also declines to "fix").
    constexpr double kRhoRoundoff = 1e-8;
    if (!(std::abs(rho) <= 1.0 + kRhoRoundoff))
        return std::numeric_limits<double>::quiet_NaN();
    if (rho > 1.0) rho = 1.0;
    if (rho < -1.0) rho = -1.0;

    const double h     = s_hi / sd1;
    const bool lo_inf  = std::isinf(sb_lo) && sb_lo < 0.0;
    const bool hi_inf  = std::isinf(sb_hi) && sb_hi > 0.0;

    constexpr double inv_sqrt_2  = 0.7071067811865475;  // 1 / √2
    constexpr double inv_sqrt_2pi = 0.3989422804014327; // 1 / √(2π)
    const double phi_h = 0.5 * std::erfc(-h * inv_sqrt_2);

    if (lo_inf && hi_inf) return phi_h;
    if (lo_inf) {
        // P(X ≤ s_hi, Y ≤ sb_hi) = Φ₂(h, sb_hi/sd2, ρ)
        return std::clamp(bvnCdf(h, sb_hi / sd2, rho), 0.0, 1.0);
    }
    if (hi_inf) {
        // P(X ≤ s_hi, Y ≥ sb_lo).  Substitute Y' = −Y (ρ → −ρ):
        //   = P(X ≤ s_hi, Y' ≤ −sb_lo) = Φ₂(h, −sb_lo/sd2, −ρ)
        return std::clamp(bvnCdf(h, -sb_lo / sd2, -rho), 0.0, 1.0);
    }

    // Both Y bounds finite — direct 1-D integration of the conditional
    // tail probability over [a, b] = [sb_lo/sd2, sb_hi/sd2]:
    //
    //   P(X ≤ s_hi, sb_lo ≤ Y ≤ sb_hi)
    //     = ∫_a^b φ(u) · Φ((h − ρ u) / √(1 − ρ²)) du
    //
    // Evaluated by 20-point Gauss-Legendre quadrature on the finite
    // interval [a, b].  Subtraction-free.
    const double a = sb_lo / sd2;
    const double b = sb_hi / sd2;
    if (b <= a) return 0.0;

    // Edge case: |ρ| → 1.  The conditional CDF degenerates to a step
    // function and Gauss-Legendre loses accuracy.  Fall back to the
    // 2-term inclusion-exclusion (bvnCdf handles this regime via its
    // own asymptotic expansion).
    if (std::abs(rho) >= 1.0 - 1e-12) {
        const double p = bvnCdf(h, b, rho) - bvnCdf(h, a, rho);
        return std::clamp(p, 0.0, 1.0);
    }

    // 20-point Gauss-Legendre rule; the nodes and weights are the file-scope
    // kGl20Node / kGl20Weight so that this routine and pmvnorm2dHalfRectLog
    // demonstrably share one quadrature.
    const double *const x20 = kGl20Node;
    const double *const w20 = kGl20Weight;

    const double half_len      = 0.5 * (b - a);
    const double mid           = 0.5 * (a + b);
    const double inv_sqrt_1mr2 = 1.0 / std::sqrt((1.0 - rho) * (1.0 + rho));

    double sum = 0.0;
    for (int i = 0; i < 10; ++i) {
        const double dx   = half_len * x20[i];
        const double u_p  = mid + dx;
        const double u_m  = mid - dx;
        const double arg_p = (h - rho * u_p) * inv_sqrt_1mr2;
        const double arg_m = (h - rho * u_m) * inv_sqrt_1mr2;
        const double f_p   = inv_sqrt_2pi * std::exp(-0.5 * u_p * u_p)
                             * 0.5 * std::erfc(-arg_p * inv_sqrt_2);
        const double f_m   = inv_sqrt_2pi * std::exp(-0.5 * u_m * u_m)
                             * 0.5 * std::erfc(-arg_m * inv_sqrt_2);
        sum += w20[i] * (f_p + f_m);
    }
    const double p = half_len * sum;
    return std::clamp(p, 0.0, 1.0);
}

// ──────────────────────────────────────────────────────────────────────
// § 2b  The same probability in the log domain
// ──────────────────────────────────────────────────────────────────────

namespace {

// Truncation depth: the integration range is cut where the integrand has
// fallen this far below its peak in the LOG.  The integrand is log-concave
// with second derivative ≤ −1 everywhere (see below), so the omitted mass is
// at most e^(−kLnDropTrunc)/√(2·kLnDropTrunc) times the peak height, against
// a total that is at least a fixed fraction of it: at 64 that bound is 4e-28,
// twelve orders below the last bit of the answer.  01_numerics §3.4(b)
// suggests 750, which is the depth at which the LINEAR integrand would
// underflow; in the log domain that depth is unreachable long before it is
// useful, and paying for it would multiply the node count by six for no
// change in any digit.
constexpr double kLnDropTrunc = 64.0;

// Target log-variation of the integrand across one panel.  The 20-point rule
// integrates a factor of e^(kLnPanelVar/2) of variation over a half-panel to
// well below double precision.
constexpr double kLnPanelVar = 8.0;

// Ceiling on the panel width in units of the standardized variable, so that a
// flat integrand still gets several panels across its support.
constexpr double kPanelWidthMax = 2.0;

// Ceiling on how far the Φ argument may move within one panel while it is
// near zero, and on how far it may move toward zero from far away.  Without
// it a panel chosen from the local gradient and curvature can step clean over
// the transition of Φ, which for |ρ| → 1 is a cliff of width 1/(c|ρ|) that
// neither the gradient nor the curvature at the panel's near edge can see.
constexpr double kPhiArgStep = 4.0;

// Hard cap on panels per side.  Never approached in practice (the measured
// maximum over a 2400-point (h, ρ, bound) sweep is 41 panels for both sides
// together); it exists so that a pathological argument cannot spin.
constexpr int kMaxPanels = 4096;

// ln ∫_lo^hi φ(u)·Φ(c·(h − ρu)) du, with c = 1/√(1−ρ²).  `lo` may be −∞.
//
// The integrand is positive, so the quadrature sum is a sum of positive
// terms and its logarithm is a log-sum-exp: normalizing by the value at the
// peak keeps every exponential argument ≤ 0, so nothing overflows and the
// terms that underflow are exactly the ones that do not matter.
//
// Two structural facts are used and both are exact:
//   * ln of the integrand is strictly concave, so the peak is unique and a
//     bisection on its derivative locates it unconditionally;
//   * its second derivative is −1 + ρ²c²·H′(c(h−ρu)) where H = φ/Φ is the
//     reverse hazard and H′ < 0 everywhere, hence ≤ −1.  The truncation
//     bound above rests on that constant.
double lnCondTailIntegral(
    double h,
    double rho,
    double lo,
    double hi
) {
    const double inf = std::numeric_limits<double>::infinity();
    if (!(hi > lo)) return -inf;

    const double c   = 1.0 / std::sqrt((1.0 - rho) * (1.0 + rho));
    const double rc  = rho * c;
    const double rc2 = rc * rc;
    const double xRate = std::fabs(rc);          // |dx/du| for x = c(h − ρu)

    const auto lnF = [&](double u) {
        return -0.5 * u * u - kLogSqrt2Pi + pnormLog(c * (h - rho * u), 0.0, 1.0, true);
    };

    // Everything the panel walk needs at a point, sharing the one pnormLog the
    // three quantities have in common.  `curv` is |d²/du² ln F| =
    // 1 + (ρc)²·(xH + H²) with H = φ/Φ the reverse hazard, which is 1 where
    // Φ ≈ 1 and 1 + (ρc)² deep in Φ's tail; using the local value rather than
    // the global bound is what keeps the panel count finite as |ρ| → 1.
    struct Local { double ln, d1, curv; };
    const auto localAt = [&](double u) {
        const double x     = c * (h - rho * u);
        const double lnPhi = pnormLog(x, 0.0, 1.0, /*lower_tail=*/true);
        // H is formed in the log domain: for x → −∞ both φ and Φ underflow
        // separately while their ratio tends to −x.
        const double H     = std::exp((-0.5 * x * x - kLogSqrt2Pi) - lnPhi);
        double q = x * H + H * H;
        if (!(q > 0.0)) q = 0.0;
        return Local{-0.5 * u * u - kLogSqrt2Pi + lnPhi, -u - rc * H, 1.0 + rc2 * q};
    };

    // Unconstrained peak, by bracket expansion plus bisection on d/du ln F,
    // which is strictly decreasing.
    double uStar = 0.0;
    {
        const double g0 = localAt(0.0).d1;
        double a = 0.0, b = 0.0;
        if (g0 > 0.0) {
            a = 0.0; b = 1.0;
            for (int i = 0; i < 200 && localAt(b).d1 > 0.0; ++i) { a = b; b *= 2.0; }
        } else if (g0 < 0.0) {
            b = 0.0; a = -1.0;
            for (int i = 0; i < 200 && localAt(a).d1 < 0.0; ++i) { b = a; a *= 2.0; }
        }
        for (int i = 0; i < 100 && b > a; ++i) {
            const double mid = 0.5 * (a + b);
            if (mid == a || mid == b) break;
            if (localAt(mid).d1 > 0.0) a = mid; else b = mid;
        }
        uStar = 0.5 * (a + b);
    }
    const double uPeak  = std::fmin(std::fmax(uStar, lo), hi);
    const Local  atPeak = localAt(uPeak);
    const double lnPeak = atPeak.ln;
    if (!std::isfinite(lnPeak)) return -inf;

    const auto panelWidth = [&](double u, const Local &st) {
        const double g = std::fabs(st.d1);
        // Largest w with g·w + curv·w²/2 ≤ kLnPanelVar, written so that no
        // cancellation occurs when g dominates.
        double w = 2.0 * kLnPanelVar / (g + std::sqrt(g * g + 2.0 * st.curv * kLnPanelVar));
        if (w > kPanelWidthMax) w = kPanelWidthMax;
        if (xRate > 0.0) {
            const double x  = c * (h - rho * u);
            const double wx = std::fmax(kPhiArgStep, 0.5 * std::fabs(x)) / xRate;
            if (wx < w) w = wx;
        }
        return w;
    };

    double total = 0.0;
    for (int dir = -1; dir <= 1; dir += 2) {
        const double bound = (dir < 0) ? lo : hi;
        double u = uPeak;
        Local st = atPeak;
        for (int panel = 0; panel < kMaxPanels; ++panel) {
            if (dir < 0 ? !(u > bound) : !(u < bound)) break;
            const double w = panelWidth(u, st);
            double next = u + dir * w;
            if (dir < 0 ? next < bound : next > bound) next = bound;
            if (next == u) break;
            const double halfLen = 0.5 * std::fabs(next - u);
            const double mid     = 0.5 * (u + next);
            double sum = 0.0;
            for (int i = 0; i < 10; ++i) {
                const double d = halfLen * kGl20Node[i];
                sum += kGl20Weight[i] * (std::exp(lnF(mid + d) - lnPeak) +
                                         std::exp(lnF(mid - d) - lnPeak));
            }
            total += halfLen * sum;
            u  = next;
            st = localAt(u);
            if (st.ln < lnPeak - kLnDropTrunc) break;
        }
    }
    if (!(total > 0.0)) return -inf;
    return lnPeak + std::log(total);
}

// log of a linear probability, with an underflowed or empty probability
// reported as −∞ rather than as NaN or +∞.
double lnOfLinearProbability(double p) {
    if (std::isnan(p)) return p;
    if (!(p > 0.0)) return -std::numeric_limits<double>::infinity();
    return std::log(p > 1.0 ? 1.0 : p);
}

} // namespace

double pmvnorm2dHalfRectLog(
    double s_hi,
    double sb_lo,
    double sb_hi,
    double var1,
    double cov12,
    double var2
) {
    const double inf = std::numeric_limits<double>::infinity();
    if (var1 <= 0.0 || var2 <= 0.0) return -inf;

    const double sd1 = std::sqrt(var1);
    const double sd2 = std::sqrt(var2);
    double rho = cov12 / (sd1 * sd2);

    // Same contract as the linear routine: round-off in ρ is absorbed, an
    // indefinite covariance is reported as NaN rather than clamped into a
    // different problem.  See the long comment on pmvnorm2dHalfRect.
    constexpr double kRhoRoundoff = 1e-8;
    if (!(std::fabs(rho) <= 1.0 + kRhoRoundoff))
        return std::numeric_limits<double>::quiet_NaN();
    if (rho > 1.0) rho = 1.0;
    if (rho < -1.0) rho = -1.0;

    const double h      = s_hi / sd1;
    const bool   lo_inf = std::isinf(sb_lo) && sb_lo < 0.0;
    const bool   hi_inf = std::isinf(sb_hi) && sb_hi > 0.0;

    if (lo_inf && hi_inf) return pnormLog(h, 0.0, 1.0, /*lower_tail=*/true);

    const double a = lo_inf ? -inf : sb_lo / sd2;
    const double b = hi_inf ?  inf : sb_hi / sd2;
    if (!(b > a)) return -inf;

    // |ρ| → 1: the conditional CDF degenerates to a step and the quadrature
    // has nothing smooth to resolve.  This is the one branch that stays on
    // the linear scale, exactly as 01_numerics §3.4(c) leaves it; what it may
    // NOT do is report +∞ upward, so an underflow becomes −∞ here and NA at
    // the caller.
    if (std::fabs(rho) >= 1.0 - 1e-12) {
        if (lo_inf) return lnOfLinearProbability(bvnCdf(h, b, rho));
        if (hi_inf) return lnOfLinearProbability(bvnCdf(h, -a, -rho));
        return lnOfLinearProbability(bvnCdf(h, b, rho) - bvnCdf(h, a, rho));
    }

    // P(X ≤ s_hi, Y ≥ sb_lo) becomes P(X ≤ s_hi, Y′ ≤ −sb_lo) with Y′ = −Y
    // and correlation −ρ, so all three surviving cases are the one integral.
    if (hi_inf) return lnCondTailIntegral(h, -rho, -inf, -a);
    return lnCondTailIntegral(h, rho, a, b);
}

// § 5  Logistic regression (IRLS)

Eigen::VectorXd logisticRegressionBeta(
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &y,
    double tol,
    int maxIter
) {

    const Eigen::Index n = X.rows();
    const Eigen::Index p = X.cols();

    // Design matrix with intercept column prepended
    Eigen::MatrixXd Xd(n, p + 1);
    Xd.col(0).setOnes();
    Xd.rightCols(p) = X;

    Eigen::VectorXd beta = Eigen::VectorXd::Zero(p + 1);

    for (int iter = 0; iter < maxIter; ++iter) {
        // Linear predictor → predicted probabilities
        Eigen::VectorXd eta = Xd * beta;
        Eigen::ArrayXd mu = (1.0 / (1.0 + (-eta.array()).exp()));

        // Clamp to avoid division-by-zero in weights
        mu = mu.max(1e-10).min(1.0 - 1e-10);

        // IRLS weight W = mu * (1 - mu)
        Eigen::ArrayXd w = mu * (1.0 - mu);

        // Working response  z = eta + (y - mu) / w   (scalar Eigen-array ops)
        Eigen::VectorXd z = eta.array() + (y.array() - mu) / w;

        // Weighted least squares:  (X^T W X) β_new = X^T W z
        // Build X^T diag(w) efficiently using column-wise scaling.
        Eigen::MatrixXd XtW = Xd.transpose() * w.matrix().asDiagonal();
        Eigen::MatrixXd XtWX = XtW * Xd;
        Eigen::VectorXd XtWz = XtW * z;

        Eigen::VectorXd betaNew = XtWX.selfadjointView<Eigen::Lower>().llt().solve(XtWz);

        if ((betaNew - beta).squaredNorm() < tol * tol) {
            beta = betaNew;
            break;
        }
        beta = betaNew;
    }
    return beta;
}

Eigen::VectorXd logisticRegression(
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    const Eigen::Ref<const Eigen::VectorXd> &y
) {

    Eigen::VectorXd beta = logisticRegressionBeta(X, y);

    const Eigen::Index n = X.rows();
    const Eigen::Index p = X.cols();

    Eigen::MatrixXd Xd(n, p + 1);
    Xd.col(0).setOnes();
    Xd.rightCols(p) = X;

    Eigen::ArrayXd eta = (Xd * beta).array();
    Eigen::ArrayXd mu = 1.0 / (1.0 + (-eta).exp());

    return (1.0 - (1.0 - mu).sqrt()).matrix();
}

// § 6  Nelder-Mead simplex optimiser

OptimResult nelderMead(
    std::function<double(const std::vector<double> &)> f,
    const std::vector<double> &init,
    double tol,
    int maxIter
) {

    const int n = static_cast<int>(init.size());
    const int nv = n + 1; // simplex has n+1 vertices

    // Build initial simplex: init + unit perturbations
    std::vector<std::vector<double> > simplex(nv, init);
    std::vector<double> fvals(nv);
    for (int i = 0; i < n; ++i) {
        double delta = (init[i] == 0.0) ? 0.05 : 0.05 * std::abs(init[i]);
        simplex[i + 1][i] += delta;
    }
    for (int i = 0; i < nv; ++i)
        fvals[i] = f(simplex[i]);

    // Scratch vectors (allocated once)
    std::vector<double> centroid(n), xr(n), xe(n), xc(n);

    int iter = 0;
    for (; iter < maxIter; ++iter) {
        // Find indices of best (ilo), worst (ihi), second-worst (inhi)
        int ilo = 0, ihi = 0, inhi = 0;
        for (int i = 0; i < nv; ++i) {
            if (fvals[i] < fvals[ilo]) ilo = i;
            if (fvals[i] > fvals[ihi]) ihi = i;
        }
        inhi = ilo;
        for (int i = 0; i < nv; ++i)
            if (i != ihi && fvals[i] > fvals[inhi]) inhi = i;

        // Convergence: simplex diameter
        double diam = 0.0;
        for (int i = 0; i < nv; ++i)
            for (int j = 0; j < n; ++j) {
                double d = std::abs(simplex[i][j] - simplex[ilo][j]);
                if (d > diam) diam = d;
            }
        if (diam < tol) break;

        // Centroid (exclude worst)
        for (int j = 0; j < n; ++j) {
            double s = 0.0;
            for (int i = 0; i < nv; ++i)
                if (i != ihi) s += simplex[i][j];
            centroid[j] = s / n;
        }

        // Reflection
        for (int j = 0; j < n; ++j)
            xr[j] = centroid[j] + (centroid[j] - simplex[ihi][j]);
        double fr = f(xr);

        if (fr < fvals[ilo]) {
            // Expansion
            for (int j = 0; j < n; ++j)
                xe[j] = centroid[j] + 2.0 * (xr[j] - centroid[j]);
            double fe = f(xe);
            if (fe < fr) {
                simplex[ihi] = xe;
                fvals[ihi] = fe;
            } else {
                simplex[ihi] = xr;
                fvals[ihi] = fr;
            }
        } else if (fr < fvals[inhi]) {
            simplex[ihi] = xr;
            fvals[ihi] = fr;
        } else {
            // Contraction
            bool outside = fr < fvals[ihi];
            const auto &xref = outside ? xr : simplex[ihi];
            double fref = outside ? fr : fvals[ihi];
            for (int j = 0; j < n; ++j)
                xc[j] = centroid[j] + 0.5 * (xref[j] - centroid[j]);
            double fc = f(xc);
            if (fc <= fref) {
                simplex[ihi] = xc;
                fvals[ihi] = fc;
            } else {
                // Shrink towards best
                for (int i = 0; i < nv; ++i) {
                    if (i == ilo) continue;
                    for (int j = 0; j < n; ++j)
                        simplex[i][j] = simplex[ilo][j] + 0.5 * (simplex[i][j] - simplex[ilo][j]);
                    fvals[i] = f(simplex[i]);
                }
            }
        }
    }

    int ilo = 0;
    for (int i = 1; i < nv; ++i)
        if (fvals[i] < fvals[ilo]) ilo = i;

    return {simplex[ilo], fvals[ilo], iter};
}

// § 6.5  BFGS quasi-Newton optimiser with numerical gradient
//
// Inverse Hessian approximation updated by the BFGS rank-2 formula:
//
//   H_{k+1} = (I − ρ s yᵀ) H_k (I − ρ y sᵀ) + ρ s sᵀ
//
// where s = x_{k+1} − x_k,  y = ∇f(x_{k+1}) − ∇f(x_k),  ρ = 1 / yᵀs.
// Gradient computed by central differences with step h_i = ε·(|x_i| + 1).
// Step length found by Armijo backtracking with initial α = 1.

OptimResult bfgs(
    std::function<double(const std::vector<double> &)> f,
    const std::vector<double> &init,
    double reltol,
    int maxIter
) {
    const int n = static_cast<int>(init.size());

    Eigen::VectorXd x  = Eigen::Map<const Eigen::VectorXd>(init.data(), n);
    Eigen::MatrixXd H  = Eigen::MatrixXd::Identity(n, n);

    auto evalAt = [&](const Eigen::VectorXd &v) {
        std::vector<double> tmp(v.data(), v.data() + n);
        return f(tmp);
    };
    auto centralGrad = [&](const Eigen::VectorXd &v, double /*fv*/) {
        constexpr double eps = 1e-6;
        Eigen::VectorXd g(n);
        Eigen::VectorXd vt = v;
        for (int i = 0; i < n; ++i) {
            const double h  = eps * (std::abs(v[i]) + 1.0);
            const double xi = v[i];
            vt[i] = xi + h;
            const double fp = evalAt(vt);
            vt[i] = xi - h;
            const double fm = evalAt(vt);
            vt[i] = xi;
            g[i] = (fp - fm) / (2.0 * h);
        }
        return g;
    };

    double fx = evalAt(x);
    Eigen::VectorXd g = centralGrad(x, fx);

    int iter = 0;
    for (; iter < maxIter; ++iter) {
        Eigen::VectorXd d = -(H * g);
        double dg = d.dot(g);
        // If the search direction is not a descent direction, reset the
        // Hessian estimate to identity and take a steepest-descent step.
        if (!(dg < 0.0)) {
            H.setIdentity();
            d  = -g;
            dg = d.dot(g);
            if (!(dg < 0.0)) break; // gradient is zero (or NaN) → done
        }

        // Armijo backtracking line search.
        constexpr double c1 = 1e-4;
        double alpha = 1.0;
        Eigen::VectorXd x_new(n);
        double f_new = fx;
        bool ok = false;
        for (int ls = 0; ls < 30; ++ls) {
            x_new = x + alpha * d;
            f_new = evalAt(x_new);
            if (std::isfinite(f_new) && f_new <= fx + c1 * alpha * dg) {
                ok = true;
                break;
            }
            alpha *= 0.5;
        }
        if (!ok) break; // line search failed → return current best

        Eigen::VectorXd g_new = centralGrad(x_new, f_new);
        Eigen::VectorXd s = x_new - x;
        Eigen::VectorXd y = g_new - g;
        const double sy = s.dot(y);

        // BFGS inverse-Hessian update (rank-2).  Skip when sᵀy is not
        // positive (curvature condition violated).
        if (sy > 1e-14) {
            Eigen::VectorXd Hy = H * y;
            const double yHy = y.dot(Hy);
            H.noalias() += ((sy + yHy) / (sy * sy)) * (s * s.transpose())
                         - (1.0 / sy) * (Hy * s.transpose() + s * Hy.transpose());
        }

        const double df = std::abs(fx - f_new);
        x  = x_new;
        g  = g_new;
        fx = f_new;
        if (df <= reltol * (std::abs(fx) + reltol)) {
            ++iter;
            break;
        }
    }

    return { std::vector<double>(x.data(), x.data() + n), fx, iter };
}

// § 5  Inverse rank normal transform (Blom, average-rank ties)

Eigen::VectorXd inverseRankNormal(const Eigen::Ref<const Eigen::VectorXd> &y) {
    const Eigen::Index N = y.size();
    if (N == 0) throw std::runtime_error("inverseRankNormal: empty input");

    std::vector<Eigen::Index> idx(static_cast<size_t>(N));
    for (Eigen::Index i = 0; i < N; ++i) idx[static_cast<size_t>(i)] = i;
    std::stable_sort(idx.begin(), idx.end(),
                     [&](Eigen::Index a, Eigen::Index b) { return y[a] < y[b]; });

    Eigen::VectorXd out(N);
    const double denom = static_cast<double>(N) + 0.25;

    Eigen::Index i = 0;
    while (i < N) {
        Eigen::Index j = i + 1;
        while (j < N && y[idx[static_cast<size_t>(j)]] == y[idx[static_cast<size_t>(i)]]) ++j;
        // 1-based ranks i+1 .. j ; midpoint = (i+1 + j) / 2.
        const double avgRank = (static_cast<double>(i + 1) + static_cast<double>(j)) * 0.5;
        const double p = (avgRank - 0.375) / denom;
        const double z = qnorm(p);
        for (Eigen::Index k = i; k < j; ++k) out[idx[static_cast<size_t>(k)]] = z;
        i = j;
    }
    return out;
}

} // namespace math
