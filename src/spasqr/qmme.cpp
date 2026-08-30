// qmme.cpp — QMME solver implementation for Gaussian smoothed QR.
//
// Iteration (paper Algorithm 1):
//   y^k = β^k + (l/(l+2)) (β^k − β^{k−1})
//   β^{k+1} = y^k − H^{-1} ∇f(y^k)
// with restart when f(β^{k+1}) > f(β^k) or l == P.
//
// Loss / gradient (averaged form):
//   ℓ_{h,τ}(u) = (h/√(2π)) e^{−u²/(2h²)} + (u/2)[1 − 2Φ(−u/h)] + (τ−½) u
//   ℓ'_{h,τ}(u) = τ − Φ(−u/h)
//   ∇f(β) = −(1/n) Z^T ψ = (1/n) Z^T [Φ(−r/h) − τ]   with r = Y − Zβ
//   H upper bound: (1/(n√(2π)h)) Z^T Z + δI
//
// Internal coordinate system: standardize X to z-scores per column,
// prepend an all-ones intercept column, center Y. Solve in (Z, Yc)
// space, un-standardize at the end.

#include "spasqr/qmme.hpp"
#include "util/math_helper.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace qmme {

namespace {

inline double smoothedQuantileLossAvg(
    const Eigen::VectorXd &r,
    double tau,
    double h
) {
    const int n = static_cast<int>(r.size());
    const double inv_sqrt2pi_h = 1.0 / std::sqrt(2.0 * M_PI);
    const double h_over_sqrt2pi = h * inv_sqrt2pi_h;
    const double inv_2h2 = 1.0 / (2.0 * h * h);
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
        const double u = r(i);
        const double phi_neg = math::pnorm(-u / h);
        const double term1 = h_over_sqrt2pi * std::exp(-u * u * inv_2h2);
        const double term2 = 0.5 * u * (1.0 - 2.0 * phi_neg);
        const double term3 = (tau - 0.5) * u;
        sum += term1 + term2 + term3;
    }
    return sum / static_cast<double>(n);
}

inline void smoothedQrGradAvg(
    const Eigen::MatrixXd &Z,
    const Eigen::VectorXd &r,
    double tau,
    double h,
    Eigen::VectorXd &derBuf,
    Eigen::VectorXd &gradOut
) {
    const int n = static_cast<int>(r.size());
    for (int i = 0; i < n; ++i)
        derBuf(i) = math::pnorm(-r(i) / h) - tau;
    gradOut.noalias() = (1.0 / static_cast<double>(n)) * Z.transpose() * derBuf;
}

inline double empiricalQuantile(
    const Eigen::VectorXd &v,
    double tau
) {
    const int n = static_cast<int>(v.size());
    std::vector<double> buf(n);
    Eigen::VectorXd::Map(buf.data(), n) = v;
    std::sort(buf.begin(), buf.end());
    const double idx = tau * (n - 1);
    const int lo = static_cast<int>(std::floor(idx));
    const int hi = std::min(lo + 1, n - 1);
    const double frac = idx - lo;
    return buf[lo] * (1.0 - frac) + buf[hi] * frac;
}

} // namespace

SqrSolver::SqrSolver(
    const Eigen::MatrixXd &X,
    double delta
)
    : m_n(static_cast<int>(X.rows())),
      m_p(static_cast<int>(X.cols())),
      m_delta(delta)
{
    // Standardize: m_Z = [1 | (X - mean) / sd]
    m_mx = X.colwise().mean();
    m_sx.resize(m_p);
    for (int j = 0; j < m_p; ++j) {
        const double s = (X.col(j).array() - m_mx(j)).matrix().norm() /
                         std::sqrt(static_cast<double>(m_n - 1));
        m_sx(j) = (s > 0.0) ? 1.0 / s : 1.0;
    }

    m_Z.resize(m_n, m_p + 1);
    m_Z.col(0).setOnes();
    for (int j = 0; j < m_p; ++j)
        m_Z.col(j + 1) = (X.col(j).array() - m_mx(j)).matrix() * m_sx(j);

    // Cache Z^T Z / n  (used for H assembly per bandwidth)
    m_ZtZ_over_n.noalias() = (1.0 / static_cast<double>(m_n)) *
                             m_Z.transpose() * m_Z;
}

void SqrSolver::prepareBandwidth(double h) {
    if (h == m_currentH) return;
    const double scale = 1.0 / (std::sqrt(2.0 * M_PI) * h);
    Eigen::MatrixXd H = scale * m_ZtZ_over_n;
    H.diagonal().array() += m_delta;
    m_chol.compute(H);
    m_currentH = h;
}

Eigen::VectorXd SqrSolver::solve(
    const Eigen::VectorXd &Y,
    double tau,
    Eigen::VectorXd *residOut,
    double tol,
    int maxIter,
    int restartPeriod,
    SolverStatus *statusOut,
    const Eigen::VectorXd *initBetaOrig
) {
    if (m_currentH <= 0.0)
        throw std::runtime_error("qmme::SqrSolver::solve: bandwidth not prepared");

    const int n = m_n;
    const int dim = m_p + 1;
    const double h = m_currentH;

    // Center Y
    const double my = Y.mean();
    Eigen::VectorXd Yc = Y.array() - my;

    // Initialize beta in standardized space.
    Eigen::VectorXd beta_curr(dim);
    if (initBetaOrig && initBetaOrig->size() == dim) {
        // Warm-start the SLOPES from β̂(τ_prev) — under the null this is
        // ≈ 0 (no help, no harm), and under real covariate effects this is
        // close to β̂(τ_new) for adjacent τ.
        //
        // The INTERCEPT is always reset to the empirical τ-quantile of
        // centered Y.  β̂_0(τ) tracks F_Y^{-1}(τ) — under H_0 with iid Y
        // this is ≈ Φ^{-1}(τ), which moves by ≈ 0.76 between τ=0.1 and
        // τ=0.3.  Carrying over β̂_0(τ_prev) would force QMME to spend
        // iterations re-finding the intercept; resetting it here is free
        // and avoids that.
        beta_curr.tail(m_p) = initBetaOrig->tail(m_p).array() / m_sx.array();
        beta_curr(0) = empiricalQuantile(Yc, tau);
    } else {
        // Cold start: intercept = empirical τ-quantile of Yc; slopes = 0.
        beta_curr.setZero();
        beta_curr(0) = empiricalQuantile(Yc, tau);
    }
    Eigen::VectorXd beta_prev = beta_curr;

    Eigen::VectorXd y(dim), beta_new(dim), grad(dim), step(dim), der(n);
    Eigen::VectorXd r_curr = Yc - m_Z * beta_curr;
    double f_curr = smoothedQuantileLossAvg(r_curr, tau, h);

    int l = 1;
    int iter = 0;
    double gradNorm = std::numeric_limits<double>::infinity();

    for (; iter < maxIter; ++iter) {
        const double extr = static_cast<double>(l) / static_cast<double>(l + 2);
        y = beta_curr + extr * (beta_curr - beta_prev);

        // Gradient at extrapolated point y
        Eigen::VectorXd r_y = Yc - m_Z * y;
        smoothedQrGradAvg(m_Z, r_y, tau, h, der, grad);
        gradNorm = grad.lpNorm<Eigen::Infinity>();
        if (gradNorm < tol) {
            // y is already stationary; accept.
            beta_prev = beta_curr;
            beta_curr = y;
            r_curr = r_y;
            f_curr = smoothedQuantileLossAvg(r_curr, tau, h);
            ++iter;
            break;
        }

        // QMME step: β_new = y − H^{-1} ∇f(y)
        step.noalias() = m_chol.solve(grad);
        beta_new = y - step;

        Eigen::VectorXd r_new = Yc - m_Z * beta_new;
        const double f_new = smoothedQuantileLossAvg(r_new, tau, h);

        // Restart: reset momentum on ascent or after P steps. On ascent we
        // additionally replace β^{k+1} with a plain MM step from β^k
        // (β^{k+1} = β^k − H^{-1}∇f(β^k)), which is a guaranteed-descent
        // quasi-Newton update. β^k itself is unchanged. This monotone variant
        // preserves the auxiliary-sequence descent that drives convergence.
        const bool ascent = (f_new > f_curr);
        if (ascent || l >= restartPeriod) {
            l = 1;
        } else {
            ++l;
        }
        if (ascent) {
            Eigen::VectorXd grad_curr(dim), der_curr(n);
            smoothedQrGradAvg(m_Z, r_curr, tau, h, der_curr, grad_curr);
            step.noalias() = m_chol.solve(grad_curr);
            beta_new = beta_curr - step;
            r_new = Yc - m_Z * beta_new;
        }

        beta_prev = beta_curr;
        beta_curr = beta_new;
        r_curr = r_new;
        f_curr = smoothedQuantileLossAvg(r_curr, tau, h);
    }

    // Final gradient at beta_curr (for honest convergence reporting)
    smoothedQrGradAvg(m_Z, r_curr, tau, h, der, grad);
    gradNorm = grad.lpNorm<Eigen::Infinity>();

    // Un-standardize: β_orig(0) = β_std(0) + my − Σ_j m_x(j) · sx(j) · β_std(j+1)
    //                 β_orig(j) = β_std(j+1) · sx(j-1)
    Eigen::VectorXd beta_orig(dim);
    beta_orig.tail(m_p) = beta_curr.tail(m_p).array() * m_sx.array();
    beta_orig(0) = beta_curr(0) + my -
                   (m_mx.array() * beta_orig.tail(m_p).transpose().array()).sum();

    if (residOut)
        *residOut = r_curr;
    if (statusOut) {
        statusOut->iter = iter;
        statusOut->converged = (gradNorm <= tol) && (iter <= maxIter);
        statusOut->finalGradNorm = gradNorm;
    }
    m_lastIters = iter;
    return beta_orig;
}

} // namespace qmme
