// math_helper.cpp — Out-of-line implementations for math_helper.hpp
#include "math_helper.hpp"

namespace math {

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

    // 20-point Gauss-Legendre nodes / weights on [-1, 1] (positive half;
    // each value is shared by ±x).  Source: standard tables, e.g.
    // Abramowitz & Stegun 25.4.30.
    static constexpr double x20[10] = {
        0.0765265211334973, 0.2277858511416451, 0.3737060887154195,
        0.5108670019508271, 0.6360536807265150, 0.7463319064601508,
        0.8391169718222188, 0.9122344282513259, 0.9639719272779138,
        0.9931285991850949
    };
    static constexpr double w20[10] = {
        0.1527533871307258, 0.1491729864726037, 0.1420961093183820,
        0.1316886384491766, 0.1181945319615184, 0.1019301198172404,
        0.0832767415767048, 0.0626720483341091, 0.0406014298003869,
        0.0176140071391521
    };

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
