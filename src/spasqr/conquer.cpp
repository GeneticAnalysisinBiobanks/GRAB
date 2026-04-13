// conquer.cpp — Smoothed quantile regression (Gaussian kernel, low-dimensional)
//
// Pure C++17 / Eigen port of the conquer R package (He et al., 2021).
// Only the Gaussian kernel variant (`smqrGauss`) is ported — the only
// kernel used by SPAsqr.
//
// Reference:
//   Xuming He, Xiaoou Pan, Kean Ming Tan, Wen-Xin Zhou (2021).
//   "conquer: Convolution-Type Smoothed Quantile Regression."
//   R package version 1.3.3, https://CRAN.R-project.org/package=conquer

#include "spasqr/conquer.hpp"
#include "util/math_helper.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace conquer {
namespace {

// ── Helpers ──────────────────────────────────────────────────────────

// Median of a vector (copies and partially sorts).
double eigMedian(const Eigen::VectorXd &v) {
    const Eigen::Index n = v.size();
    std::vector<double> buf(n);
    Eigen::VectorXd::Map(buf.data(), n) = v;
    if (n % 2 == 1) {
        std::nth_element(buf.begin(), buf.begin() + n / 2, buf.end());
        return buf[n / 2];
    }
    std::nth_element(buf.begin(), buf.begin() + n / 2 - 1, buf.end());
    double lo = buf[n / 2 - 1];
    std::nth_element(buf.begin() + n / 2, buf.begin() + n / 2, buf.end());
    double hi = buf[n / 2];
    return 0.5 * (lo + hi);
}

// MAD: 1.482602 * median(|x - median(x)|)
double eigMad(const Eigen::VectorXd &v) {
    double med = eigMedian(v);
    Eigen::VectorXd absdev = (v.array() - med).abs().matrix();
    return 1.482602 * eigMedian(absdev);
}

// R type-7 quantile
double eigQuantile(
    const Eigen::VectorXd &v,
    double prob
) {
    const Eigen::Index n = v.size();
    std::vector<double> buf(n);
    Eigen::VectorXd::Map(buf.data(), n) = v;
    std::sort(buf.begin(), buf.end());
    double idx = prob * (n - 1);
    Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
    Eigen::Index hi = std::min(lo + 1, n - 1);
    double frac = idx - lo;
    return buf[lo] * (1.0 - frac) + buf[hi] * frac;
}

// Element-wise Φ(x) (standard normal CDF) into an existing vector.
void vecPnorm(
    const Eigen::VectorXd &x,
    Eigen::VectorXd &out
) {
    out.resize(x.size());
    for (Eigen::Index i = 0; i < x.size(); ++i)
        out(i) = math::pnorm(x(i));
}

// ── Asymmetric Huber gradient ────────────────────────────────────────

void updateHuber(
    const Eigen::MatrixXd &Z,
    const Eigen::VectorXd &res,
    double tau,
    Eigen::VectorXd &der,
    Eigen::VectorXd &grad,
    int n,
    double rob,
    double n1
) {
    for (int i = 0; i < n; ++i) {
        double cur = res(i);
        if (cur > rob)
            der(i) = -2.0 * tau * rob;
        else if (cur > 0.0)
            der(i) = -2.0 * tau * cur;
        else if (cur > -rob)
            der(i) = 2.0 * (tau - 1.0) * cur;
        else
            der(i) = 2.0 * (1.0 - tau) * rob;
    }
    grad.noalias() = (n1)*Z.transpose() * der;
}

// ── Gaussian kernel gradient ─────────────────────────────────────────

void updateGauss(
    const Eigen::MatrixXd &Z,
    const Eigen::VectorXd &res,
    Eigen::VectorXd &der,
    Eigen::VectorXd &grad,
    double tau,
    double n1,
    double h1
) {
    Eigen::VectorXd arg = -res * h1;
    vecPnorm(arg, der);
    der.array() -= tau;
    grad.noalias() = n1 * Z.transpose() * der;
}

// ── Huber regression (BB gradient descent) ───────────────────────────

Eigen::VectorXd huberReg(
    const Eigen::MatrixXd &Z,
    const Eigen::VectorXd &Y,
    double tau,
    Eigen::VectorXd &der,
    Eigen::VectorXd &gradOld,
    Eigen::VectorXd &gradNew,
    int n,
    int p,
    double n1,
    double tol,
    double constTau,
    int iteMax,
    double stepMax
) {
    double rob = constTau * eigMad(Y);
    updateHuber(Z, Y, tau, der, gradOld, n, rob, n1);
    Eigen::VectorXd beta = -gradOld;
    Eigen::VectorXd betaDiff = -gradOld;
    Eigen::VectorXd res = Y - Z * beta;
    rob = constTau * eigMad(res);
    updateHuber(Z, res, tau, der, gradNew, n, rob, n1);
    Eigen::VectorXd gradDiff = gradNew - gradOld;
    int ite = 1;
    while (gradNew.lpNorm<Eigen::Infinity>() > tol && ite <= iteMax) {
        double alpha = 1.0;
        double cross = betaDiff.dot(gradDiff);
        if (cross > 0.0) {
            double a1 = cross / gradDiff.squaredNorm();
            double a2 = betaDiff.squaredNorm() / cross;
            alpha = std::min({a1, a2, stepMax});
        }
        gradOld = gradNew;
        betaDiff = -alpha * gradNew;
        beta += betaDiff;
        res -= Z * betaDiff;
        rob = constTau * eigMad(res);
        updateHuber(Z, res, tau, der, gradNew, n, rob, n1);
        gradDiff = gradNew - gradOld;
        ++ite;
    }
    return beta;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// smqrGauss — main entry point
// ══════════════════════════════════════════════════════════════════════

Eigen::VectorXd smqrGauss(
    const Eigen::MatrixXd &X,
    const Eigen::VectorXd &Y,
    double tau,
    double h,
    Eigen::VectorXd *residOut,
    double tol,
    int iteMax,
    double stepMax
) {
    const int n = static_cast<int>(X.rows());
    const int p = static_cast<int>(X.cols());
    const double constTau = 1.345;

    if (h <= 0.0) h = std::max(std::pow((std::log(n) + p) / static_cast<double>(n), 0.4), 0.05);

    const double n1 = 1.0 / n;
    const double h1 = 1.0 / h;

    // Standardize X → Z = [1 | (X - mean) / sd ]
    Eigen::RowVectorXd mx = X.colwise().mean();
    Eigen::VectorXd sx(p);
    for (int j = 0; j < p; ++j) {
        double s = (X.col(j).array() - mx(j)).matrix().norm() / std::sqrt(static_cast<double>(n - 1));
        sx(j) = (s > 0.0) ? 1.0 / s : 1.0;
    }

    Eigen::MatrixXd Z(n, p + 1);
    Z.col(0).setOnes();
    for (int j = 0; j < p; ++j)
        Z.col(j + 1) = (X.col(j).array() - mx(j)).matrix() * sx(j);

    // Center Y
    double my = Y.mean();
    Eigen::VectorXd Yc = Y.array() - my;

    // Scratch vectors
    Eigen::VectorXd der(n), gradOld(p + 1), gradNew(p + 1);

    // Phase 1: Huber initialization
    Eigen::VectorXd beta = huberReg(Z, Yc, tau, der, gradOld, gradNew, n, p, n1, tol, constTau, iteMax, stepMax);

    // Quantile intercept adjustment
    Eigen::VectorXd resNoIntercept = Yc - Z.rightCols(p) * beta.tail(p);
    beta(0) = eigQuantile(resNoIntercept, tau);

    // Phase 2: Gaussian kernel smoothed quantile regression (BB descent)
    Eigen::VectorXd res = Yc - Z * beta;
    updateGauss(Z, res, der, gradOld, tau, n1, h1);
    beta -= gradOld;
    Eigen::VectorXd betaDiff = -gradOld;
    res -= Z * betaDiff;
    updateGauss(Z, res, der, gradNew, tau, n1, h1);
    Eigen::VectorXd gradDiff = gradNew - gradOld;

    int ite = 1;
    while (gradNew.lpNorm<Eigen::Infinity>() > tol && ite <= iteMax) {
        double alpha = 1.0;
        double cross = betaDiff.dot(gradDiff);
        if (cross > 0.0) {
            double a1 = cross / gradDiff.squaredNorm();
            double a2 = betaDiff.squaredNorm() / cross;
            alpha = std::min({a1, a2, stepMax});
        }
        gradOld = gradNew;
        betaDiff = -alpha * gradNew;
        beta += betaDiff;
        res -= Z * betaDiff;
        updateGauss(Z, res, der, gradNew, tau, n1, h1);
        gradDiff = gradNew - gradOld;
        ++ite;
    }

    // Un-standardize coefficients
    beta.tail(p).array() *= sx.array();
    beta(0) += my - (mx.array() * beta.tail(p).transpose().array()).sum();

    if (residOut) *residOut = res;

    return beta;
}

} // namespace conquer
