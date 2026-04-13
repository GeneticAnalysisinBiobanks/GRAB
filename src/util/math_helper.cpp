// math_helper.cpp — Out-of-line implementations for math_helper.hpp
#include "math_helper.hpp"

namespace math {

// § 2  Bivariate normal CDF  (Genz 2004, 6-point quadrature)

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
        double kk = k;
        double hkk = hk;
        if (r < 0.0) {
            kk = -kk;
            hkk = -hkk;
        }

        if (std::abs(r) < 1.0) {
            const double as_ = (1.0 - r) * (1.0 + r);
            const double a = std::sqrt(as_);
            const double bs = (h - kk) * (h - kk);
            const double c = (4.0 - hkk) / 8.0;
            const double d = (12.0 - hkk) / 16.0;
            const double asr = -(bs / as_ + hkk) / 2.0;
            if (asr > -100.0)
                bvn = a * std::exp(asr) * (1.0 - c * (bs - as_) * (1.0 - d * bs / 5.0) / 3.0 + c * d * as_ * as_ / 5.0);
            if (-hkk < 100.0) {
                const double b = std::sqrt(bs);
                bvn -= std::exp(-hkk / 2.0) * std::sqrt(2.0 * M_PI) * 0.5 * std::erfc(-b / (std::sqrt(2.0) * a)) * b *
                       (1.0 - c * bs * (1.0 - d * bs / 5.0) / 3.0);
            }
            const double a2 = a / 2.0;
            for (int i = 0; i < 3; ++i) {
                for (int is = -1; is <= 1; is += 2) {
                    double xs = a2 * (is * x6[i] + 1.0);
                    xs *= xs;
                    const double asr2 = -(bs / xs + hkk) / 2.0;
                    if (asr2 > -100.0)
                        bvn += a2 * w6[i] * std::exp(asr2) *
                               (std::exp(-hkk * (1.0 - r) / (2.0 * (1.0 + (is * x6[i] + 1.0) * a2 * a2 / as_))) - 1.0 -
                                c * xs * (1.0 + d * xs));
                }
            }
            bvn = -bvn / (2.0 * M_PI);
        }

        if (r > 0.0) {
            double phih = 0.5 * std::erfc(-h / std::sqrt(2.0));
            double phik = 0.5 * std::erfc(-kk / std::sqrt(2.0));
            bvn += phih + phik - 1.0 + std::min(phih, phik);
            if (bvn < 0.0) bvn = 0.0;
        } else {
            bvn = -bvn;
            double phih = 0.5 * std::erfc(-h / std::sqrt(2.0));
            double phik = 0.5 * std::erfc(-kk / std::sqrt(2.0));
            if (phih - phik >= 0.0) {
                bvn = phih - phik - bvn;
                if (bvn < 0.0) bvn = 0.0;
            } else {
                bvn = 0.0;
            }
        }
    }
    return bvn;
}

double pmvnorm2d(
    double lo1,
    double hi1,
    double lo2,
    double hi2,
    double var1,
    double cov12,
    double var2
) {
    const double sd1 = std::sqrt(var1);
    const double sd2 = std::sqrt(var2);
    const double rho = cov12 / (sd1 * sd2);

    auto standardise = [](double v, double sd) -> double {
        if (std::isinf(v)) return v > 0 ? 1e15 : -1e15;
        return v / sd;
    };

    const double a1 = standardise(lo1, sd1), b1 = standardise(hi1, sd1);
    const double a2 = standardise(lo2, sd2), b2 = standardise(hi2, sd2);

    double p = bvnCdf(b1, b2, rho) - bvnCdf(a1, b2, rho) - bvnCdf(b1, a2, rho) + bvnCdf(a1, a2, rho);
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

} // namespace math
