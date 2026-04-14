// polmm.cpp — POLMM: Proportional Odds Logistic Mixed Model (C++17/Eigen)
//
// Implements null model fitting (PQL with PCG), variance ratio estimation,
// and per-marker association testing with SPA for ordinal phenotypes.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "polmm/polmm.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <thread>
#include <vector>

#include <Eigen/Dense>

// ══════════════════════════════════════════════════════════════════════
// Constants
// ══════════════════════════════════════════════════════════════════════

namespace {

constexpr int kMaxIterPQL = 50;
constexpr double kTolBeta = 1e-4;
constexpr double kTolTau = 0.01;
constexpr int kMaxIterEps = 100;
constexpr double kTolEps = 1e-4;
constexpr int kMaxIterPCG = 100;
constexpr double kTolPCG = 1e-6;
constexpr int kTraceNRun = 30;   // Hutchinson trace estimator runs
constexpr int kNSnpsPerBin = 30; // SNPs per MAC bin for variance ratio
constexpr double kMinVarRatioCV = 0.0025;

// MAC bin boundaries for variance ratio estimation (Optimization 4)
// Bin i covers [kMacBins[i], kMacBins[i+1]), last bin covers [kMacBins.back(), ∞)
static const std::vector<double> kMacBins = {20.0, 50.0, 100.0, 500.0};

// ══════════════════════════════════════════════════════════════════════
// Logistic helpers
// ══════════════════════════════════════════════════════════════════════

inline double logistic(double x) {
    if (x >= 0.0) {
        double ez = std::exp(-x);
        return 1.0 / (1.0 + ez);
    }
    double ez = std::exp(x);
    return ez / (1.0 + ez);
}

inline double dlogistic(double x) {
    double p = logistic(x);
    return p * (1.0 - p);
}

// Clamp probability to [eps, 1-eps]
inline double clampProb(double p) {
    constexpr double lo = 1e-10;
    constexpr double hi = 1.0 - 1e-10;
    return std::max(lo, std::min(hi, p));
}

// ══════════════════════════════════════════════════════════════════════
// Batch SpMV for interleaved layout (Optimization 3).
// Vectors are interleaved: x[i * stride + j] for subject i, component j.
// Iterates over COO entries once instead of `stride` times.
// When stride==1, delegates to the common SparseGRM::multiply().
// ══════════════════════════════════════════════════════════════════════

void batchSpmvInterleaved(
    const SparseGRM &grm,
    const double *x,
    double *result,
    uint32_t n,
    int stride
) {
    if (stride == 1) {
        grm.multiply(x, result, n);
        return;
    }
    const int total = static_cast<int>(n) * stride;
    std::fill(result, result + total, 0.0);
    for (const auto &e : grm.entries()) {
        const double v = e.value;
        const int ri = static_cast<int>(e.row) * stride;
        const int ci = static_cast<int>(e.col) * stride;
        for (int j = 0; j < stride; ++j)
            result[ri + j] += v * x[ci + j];
        if (e.row != e.col) {
            for (int j = 0; j < stride; ++j)
                result[ci + j] += v * x[ri + j];
        }
    }
}

// Accumulate-mode batch SpMV: result += alpha * (K ⊗ I_stride) * x
// Zero-copy — no intermediate buffer allocation.
void batchSpmvInterleavedAdd(
    const SparseGRM &grm,
    const double *x,
    double *result,
    uint32_t n,
    int stride,
    double alpha
) {
    if (stride == 1) {
        for (const auto &e : grm.entries()) {
            result[e.row] += alpha * e.value * x[e.col];
            if (e.row != e.col) result[e.col] += alpha * e.value * x[e.row];
        }
        return;
    }
    for (const auto &e : grm.entries()) {
        const double av = alpha * e.value;
        const int ri = static_cast<int>(e.row) * stride;
        const int ci = static_cast<int>(e.col) * stride;
        for (int j = 0; j < stride; ++j)
            result[ri + j] += av * x[ci + j];
        if (e.row != e.col) {
            for (int j = 0; j < stride; ++j)
                result[ci + j] += av * x[ri + j];
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
// Working matrix updates: given (beta, bVec, eps) → muMat, iRMat, etc.
// ══════════════════════════════════════════════════════════════════════

// muMat(i,j) = P(y=j | eta_i) = F(eps_j - eta_i) - F(eps_{j-1} - eta_i)
//   where F = logistic, eps_0 = -inf, eps_J = +inf
// iRMat(i,j) = 1/sqrt(muMat(i,j) * cumMuBar(i,j))
//   where cumMuBar(i,j) = 1 - sum_{k<=j} muMat(i,k)

void updateMuMat(
    const Eigen::MatrixXd &X,
    const Eigen::VectorXd &beta,
    const Eigen::VectorXd &bVec,
    const Eigen::VectorXd &eps,
    int n,
    int J,
    Eigen::MatrixXd &muMat,
    Eigen::MatrixXd &iRMat
) {
    // eta_i = X[i,:]*beta + bVec[i]
    Eigen::VectorXd eta = X * beta + bVec;

    muMat.resize(n, J);
    iRMat.resize(n, J - 1);

    for (int i = 0; i < n; ++i) {
        double prevCdf = 0.0;
        for (int j = 0; j < J; ++j) {
            double cdf_j;
            if (j == J - 1)cdf_j = 1.0;
            else cdf_j = logistic(eps(j) - eta(i));
            double mu_ij = clampProb(cdf_j - prevCdf);
            muMat(i, j) = mu_ij;
            prevCdf = cdf_j;
        }
        // Compute iRMat: working correlation inverse scales
        double cumMu = 0.0;
        for (int j = 0; j < J - 1; ++j) {
            cumMu += muMat(i, j);
            double cumMuBar = 1.0 - cumMu; // = sum_{k > j} mu(i,k)
            cumMuBar = std::max(cumMuBar, 1e-10);
            double v = muMat(i, j) * cumMuBar;
            iRMat(i, j) = 1.0 / std::sqrt(std::max(v, 1e-20));
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
// Psi operations: Psi is the working covariance of the cumulative model.
//   Psi = diag(mu) - mu*mu' per subject (J-1 × J-1 blocks)
//
// For a vector w of length n*(J-1) (stacked J-1 blocks per subject):
//   (Psi * w)[i, j] = muMat(i,j)*(1 - cumP(i,j)) * w[i,j]
//                    - sum over k≠j of cross terms
// But a simpler formulation uses the fact that for cumulative logit:
//   Psi_i = D_i - c_i * c_i'
// where D_i = diag(g'(threshold_j - eta_i)), c_i = (g'(...))
// with g = logistic density.
//
// For efficiency, we implement Psi*v and Psi^{-1}*v per subject.
// ══════════════════════════════════════════════════════════════════════

// Compute Psi_i * w_i where w_i is a block of length (J-1) for subject i
// Psi_i[j,k] = mu_i[min(j,k)] * (1 - mu_i[0..max(j,k)] cumsum)  when represented
//            = delta(j,k)*f_j - f_j*f_k  where f_j = dlogistic(eps_j - eta_i)
// Actually from the original POLMM code:
//   Psi_i = WMat_i - mMat_i * mMat_i'   (J-1 × J-1)
// where WMat_i[j] = d logistic(eps_j - eta_i) = mu_(j) * (1 - cumP <= j)
//   and mMat_i[j,k] = ... complex
//
// For simplicity and correctness, use the concrete formula:
//   Psi_i[j,k] = mu_i[min(j,k)+1..J-1 cumulative] * mu_i[0..min(j,k) cumulative]
//              = CumMu_i(min(j,k)) * (1 - CumMu_i(max(j,k)))  [for j != k this is negative]
// Actually: Psi_i[j,k] = mu(j)*(1{j>k} - CumF(j)) for the cumulative probit/logit
// Let me use the standard result for cumulative link models:
//   Psi_i = A_i - c_i c_i'
//   where A_i = diag(h_0, h_1, ..., h_{J-2}), h_j = dF(eps_j - eta_i)
//   and c_i = (h_0, h_1, ..., h_{J-2})'
// This makes Psi_i = diag(h) - h * h'  which is a rank-1 update of a diagonal.

struct PsiBlock {
    // h_j = dlogistic(eps_j - eta_i) for j = 0..J-2
    Eigen::VectorXd h; // (J-1)

    void compute(
        double eta_i,
        const Eigen::VectorXd &eps,
        int Jm1
    ) {
        h.resize(Jm1);
        for (int j = 0; j < Jm1; ++j)
            h(j) = dlogistic(eps(j) - eta_i);
    }

    // y = Psi * x = diag(h)*x - h*(h'x)
    // Optimization 1: K=2 scalar fast path when Jm1==1
    void mulVec(
        const double *x,
        double *y,
        int Jm1
    ) const {
        if (Jm1 == 1) {
            y[0] = h(0) * (1.0 - h(0)) * x[0];
            return;
        }
        double hdotx = 0.0;
        for (int j = 0; j < Jm1; ++j)
            hdotx += h(j) * x[j];
        for (int j = 0; j < Jm1; ++j)
            y[j] = h(j) * x[j] - h(j) * hdotx;
    }

    // y = Psi^{-1} * x = diag(1/h)*x + 11' * x / (1 - sum(h))
    // By Sherman-Morrison: (diag(h) - h h')^{-1} = diag(1/h) + 11'/(1 - sum(h))
    // Optimization 1: K=2 scalar fast path when Jm1==1
    void solveVec(
        const double *x,
        double *y,
        int Jm1
    ) const {
        if (Jm1 == 1) {
            double psi = std::max(h(0) * (1.0 - h(0)), 1e-20);
            y[0] = x[0] / psi;
            return;
        }
        double sumH = 0.0;
        double sumXoverH = 0.0;
        for (int j = 0; j < Jm1; ++j) {
            sumH += h(j);
            double hi = std::max(h(j), 1e-20);
            y[j] = x[j] / hi;
            sumXoverH += y[j];
        }
        double denom = std::max(1.0 - sumH, 1e-20);
        double add = sumXoverH / denom;
        for (int j = 0; j < Jm1; ++j)
            y[j] += add;
    }

};

// ══════════════════════════════════════════════════════════════════════
// PCG solver: solves Sigma * x = rhs, where Sigma = Psi^{-1} + tau * K
// Returns x = Sigma^{-1} * rhs
// Vectors are of length n*(J-1), with blocks of (J-1) per subject.
// ══════════════════════════════════════════════════════════════════════

// Sigma * v = Psi^{-1} * v + tau * (K ⊗ I_{J-1}) * v
// where (K ⊗ I_{J-1}) * v means: for each dimension j in 0..J-2,
// extract the j-th element from each subject's block to form an n-vector,
// multiply by K, then scatter back.
void sigmaMultiply(
    const std::vector<PsiBlock> &psiBlocks,
    const SparseGRM &grm,
    double tau,
    int n,
    int Jm1,
    const double *v,
    double *result
) {
    // Part 1: Psi^{-1} * v (block diagonal, per subject)
    for (int i = 0; i < n; ++i)
        psiBlocks[i].solveVec(v + i * Jm1, result + i * Jm1, Jm1);

    // Part 2: result += tau * (K ⊗ I_{J-1}) * v
    // Zero-copy accumulate — no intermediate buffer.
    batchSpmvInterleavedAdd(grm, v, result, static_cast<uint32_t>(n), Jm1, tau);
}

// Preconditioner: block diagonal of Sigma.
// For each subject, the (J-1)×(J-1) block is M_i = Psi_i^{-1} + tau*K_ii*I.
// We solve M_i * z_i = r_i per block using Sherman-Morrison.
//
// Psi_i^{-1} = diag(1/h) + 11'/(1 - sum(h))  (from PsiBlock::solveVec)
// So M_i = diag(1/h_j + tau*d_i) + 11'/(1 - sum(h))
//        = diag(D_j) + uu'   where D_j = 1/h_j + tau*d_i,  u = 1/sqrt(1-sum(h))
//
// By Sherman-Morrison:
//   M_i^{-1} = diag(a) - a*a' / ((1-sum(h))^{-1} + sum(a))
//   where a_j = 1/D_j = 1/(1/h_j + tau*d_i) = h_j/(1 + tau*d_i*h_j)
void precondSolve(
    const std::vector<PsiBlock> &psiBlocks,
    const std::vector<double> &diagK,
    double tau,
    int n,
    int Jm1,
    const double *r,
    double *z
) {
    for (int i = 0; i < n; ++i) {
        const double *ri = r + i * Jm1;
        double *zi = z + i * Jm1;
        const auto &h = psiBlocks[i].h;
        const double td = tau * diagK[i];

        if (Jm1 == 1) {
            // J=2 fast path: M = 1/(h*(1-h)) + tau*d = (1 + tau*d*h*(1-h)) / (h*(1-h))
            double psi = std::max(h(0) * (1.0 - h(0)), 1e-20);
            zi[0] = ri[0] / (1.0 / psi + td);
            return;
        }

        // General case: a_j = h_j / (1 + td * h_j)
        double sumH = 0.0;
        double sumA = 0.0;
        double adotR = 0.0;
        for (int j = 0; j < Jm1; ++j) {
            double hj = std::max(h(j), 1e-20);
            double aj = hj / (1.0 + td * hj);
            sumH += hj;
            sumA += aj;
            zi[j] = aj * ri[j];   // diag(a) * r
            adotR += aj * ri[j];
        }
        // Subtract rank-1 correction: a*a'*r / ((1-sumH)^{-1} + sumA)
        double denom = 1.0 / std::max(1.0 - sumH, 1e-20) + sumA;
        double scale = adotR / denom;
        for (int j = 0; j < Jm1; ++j) {
            double hj = std::max(h(j), 1e-20);
            double aj = hj / (1.0 + td * hj);
            zi[j] -= aj * scale;
        }
    }
}

Eigen::VectorXd pcgSolve(
    const std::vector<PsiBlock> &psiBlocks,
    const SparseGRM &grm,
    const std::vector<double> &diagK,
    double tau,
    int n,
    int Jm1,
    const Eigen::VectorXd &rhs
) {
    const int dim = n * Jm1;
    Eigen::VectorXd x = Eigen::VectorXd::Zero(dim);
    Eigen::VectorXd r = rhs; // r = rhs - Sigma * x, but x=0 so r=rhs
    Eigen::VectorXd z(dim), p(dim), Ap(dim);

    // z = M^{-1} r
    precondSolve(psiBlocks, diagK, tau, n, Jm1, r.data(), z.data());
    p = z;
    double rz = r.dot(z);

    for (int iter = 0; iter < kMaxIterPCG; ++iter) {
        // Ap = Sigma * p
        sigmaMultiply(psiBlocks, grm, tau, n, Jm1, p.data(), Ap.data());
        double pAp = p.dot(Ap);
        if (std::abs(pAp) < 1e-30) break;
        double alpha = rz / pAp;
        x += alpha * p;
        r -= alpha * Ap;

        double rNorm = r.norm();
        if (rNorm < kTolPCG * rhs.norm()) break;

        precondSolve(psiBlocks, diagK, tau, n, Jm1, r.data(), z.data());
        double rz_new = r.dot(z);
        double beta = rz_new / std::max(rz, 1e-30);
        p = z + beta * p;
        rz = rz_new;
    }
    return x;
}

// ══════════════════════════════════════════════════════════════════════
// CLM initialization — cumulative logistic model via IRLS
// (replaces R's ordinal::clm)
// Fits: logit(P(y <= j)) = eps_j - X*beta
// ══════════════════════════════════════════════════════════════════════

void fitCLM(
    const Eigen::VectorXi &yVec,
    const Eigen::MatrixXd &X,
    int n,
    int J,
    int p,
    Eigen::VectorXd &beta,
    Eigen::VectorXd &eps
) {
    // Initialize cutpoints from marginal proportions
    eps.resize(J - 1);
    Eigen::VectorXd cumProp(J);
    cumProp.setZero();
    for (int i = 0; i < n; ++i)
        cumProp(yVec(i)) += 1.0;
    cumProp /= static_cast<double>(n);

    double cumSum = 0.0;
    for (int j = 0; j < J - 1; ++j) {
        cumSum += cumProp(j);
        cumSum = std::max(0.01, std::min(0.99, cumSum));
        eps(j) = std::log(cumSum / (1.0 - cumSum)); // logit
    }
    beta.setZero(p);

    // IRLS iterations
    const int maxIter = 50;
    const int nTheta = (J - 1) + p; // total parameters: eps + beta

    for (int iter = 0; iter < maxIter; ++iter) {
        // Build working response and weights
        // For cumulative logit: P(y <= j) = logistic(eps_j - X*beta)
        Eigen::VectorXd eta = X * beta;

        // Gradient and Hessian w.r.t. [eps; beta]
        Eigen::VectorXd grad = Eigen::VectorXd::Zero(nTheta);
        Eigen::MatrixXd hess = Eigen::MatrixXd::Zero(nTheta, nTheta);

        for (int i = 0; i < n; ++i) {
            int yi = yVec(i);
            // Compute P(y=yi) and derivatives
            // gamma_j = eps_j - eta_i
            // F_j = logistic(gamma_j), f_j = F_j*(1-F_j)
            // P(y=yi) = F_{yi} - F_{yi-1}  (F_{-1}=0, F_{J-1}=1)

            double F_yi = (yi < J - 1) ? logistic(eps(yi) - eta(i)) : 1.0;
            double F_yim1 = (yi > 0) ? logistic(eps(yi - 1) - eta(i)) : 0.0;
            double p_yi = clampProb(F_yi - F_yim1);

            double f_yi = (yi < J - 1) ? dlogistic(eps(yi) - eta(i)) : 0.0;
            double f_yim1 = (yi > 0) ? dlogistic(eps(yi - 1) - eta(i)) : 0.0;

            // Gradient w.r.t eps_j: dlog(p_yi)/deps_j
            //   = (f_yi * 1{j==yi} - f_yim1 * 1{j==yi-1}) / p_yi
            if (yi < J - 1) grad(yi) += f_yi / p_yi;
            if (yi > 0) grad(yi - 1) -= f_yim1 / p_yi;

            // Gradient w.r.t beta: dlog(p_yi)/dbeta = -(f_yi - f_yim1) / p_yi * X[i,:]
            double dLogPdEta = -(f_yi - f_yim1) / p_yi;
            grad.tail(p) += dLogPdEta * X.row(i).transpose();

            // Approximate Hessian (Fisher information, negative expected)
            // For simplicity, use the outer product of scores (BFGS-like)
            // or the exact Fisher info. Use the diagonal of working weights.
            double w = (f_yi - f_yim1) * (f_yi - f_yim1) / (p_yi * p_yi);
            w = std::max(w, 1e-10);

            // Hessian contributions for eps-eps, eps-beta, beta-beta blocks
            // Use score outer product as a simple approximation
            Eigen::VectorXd score_i(nTheta);
            score_i.setZero();
            if (yi < J - 1) score_i(yi) = f_yi / p_yi;
            if (yi > 0) score_i(yi - 1) = -f_yim1 / p_yi;
            score_i.tail(p) = dLogPdEta * X.row(i).transpose();
            hess += score_i * score_i.transpose();
        }

        // Newton step: theta_new = theta + H^{-1} * grad
        // (using observed info = outer product of scores)
        // Add small ridge for stability
        hess.diagonal().array() += 1e-6;
        Eigen::VectorXd delta = hess.ldlt().solve(grad);

        // Update
        for (int j = 0; j < J - 1; ++j)
            eps(j) += delta(j);
        beta += delta.tail(p);

        if (delta.norm() < 1e-6) break;
    }
    // Ensure cutpoints are ordered
    for (int j = 1; j < J - 1; ++j)
        eps(j) = std::max(eps(j), eps(j - 1) + 0.01);
}

// ══════════════════════════════════════════════════════════════════════
// PQL null model fitting
// ══════════════════════════════════════════════════════════════════════

// Y indicator matrix: yMat(i, j) = 1 if yVec(i) == j
Eigen::MatrixXd makeYMat(
    const Eigen::VectorXi &yVec,
    int n,
    int J
) {
    Eigen::MatrixXd yMat = Eigen::MatrixXd::Zero(n, J);
    for (int i = 0; i < n; ++i)
        yMat(i, yVec(i)) = 1.0;
    return yMat;
}

// Construct the n*(J-1) working response vector Y_w for PQL
// Y_w = Psi^{-1} * (yMat[:, 0..J-2] - muMat[:, 0..J-2])_stacked + eta_expanded
// Actually in POLMM: WVec = Psi^{-1}*(y_tilde - mu_tilde) where
// y_tilde_i = (1{y_i <= 0}, ..., 1{y_i <= J-2})' and mu_tilde_i = cumulative mu.
// But for our formulation, the working variate for PQL is:
//   Y_star_i = X*beta + b_i + Psi_i^{-1} * (y_tilde_i - cumMu_i)
// and we solve the mixed model equations on that.
//
// Instead, let's follow the POLMM reference implementation more closely.
// The PQL inner loop directly updates beta, bVec, and eps.

void fitNullModel(
    const Eigen::VectorXi &yVec,
    const Eigen::MatrixXd &X,
    const SparseGRM &grm,
    int n,
    int J,
    int p,
    POLMMNullModel &null,
    const char *tag = "",
    int nLocalThreads = 1
) {
    const int Jm1 = J - 1;

    // Use cached GRM diagonal (self-relatedness per subject)
    const auto &diagK = grm.diagonal();

    // 1. CLM initialization
    infoMsg("[%s] CLM initialization (J=%d levels, p=%d covariates)", tag, J, p);
    Eigen::VectorXd beta, eps;
    fitCLM(yVec, X, n, J, p, beta, eps);
    infoMsg("[%s] CLM done. eps[0]=%.4f, eps[J-2]=%.4f", tag, eps(0), eps(Jm1 - 1));

    Eigen::VectorXd bVec = Eigen::VectorXd::Zero(n);
    double tau = 0.1; // initial variance component

    Eigen::MatrixXd muMat, iRMat;
    std::vector<PsiBlock> psiBlocks(n);

    // Make indicator matrices
    Eigen::MatrixXd yMat = makeYMat(yVec, n, J);

    // 2. PQL iteration (outer loop over tau, inner loop over beta, b, eps)
    infoMsg("[%s] Starting PQL iteration", tag);
    for (int outerIter = 0; outerIter < kMaxIterPQL; ++outerIter) {
        double tau_old = tau;

        // Inner loop: update beta, bVec, eps with fixed tau
        for (int innerIter = 0; innerIter < 10; ++innerIter) {
            Eigen::VectorXd beta_old = beta;

            // Update working matrices
            updateMuMat(X, beta, bVec, eps, n, J, muMat, iRMat);

            // Compute PsiBlocks
            Eigen::VectorXd eta = X * beta + bVec;
            for (int i = 0; i < n; ++i)
                psiBlocks[i].compute(eta(i), eps, Jm1);

            // Working response: Y_tilde = cumulative indicators - cumulative mu
            // Y_tilde_i[j] = 1{y_i <= j} - F(eps_j - eta_i) for j = 0..J-2
            // Then: working variate W = Psi^{-1} * Y_tilde
            Eigen::VectorXd ytilde(n * Jm1);
            for (int i = 0; i < n; ++i) {
                double cumY = 0.0, cumMu = 0.0;
                for (int j = 0; j < Jm1; ++j) {
                    cumY += yMat(i, j);
                    cumMu += muMat(i, j);
                    ytilde(i * Jm1 + j) = cumY - cumMu;
                }
            }

            // Working response in the mixed model:
            // Y_star = X_exp * beta_exp + b_exp + Psi^{-1} * ytilde
            // where X_exp is the expanded design matrix (n*Jm1 × (Jm1+p))
            // and b_exp repeats b_i for each j block.
            //
            // The POLMM approach directly solves:
            //   Sigma^{-1} * (Y_w - X_exp * beta - b_exp) = 0
            // where Y_w = eta_exp + Psi^{-1} * ytilde
            //
            // For the beta update: solve normal equations
            //   (X_exp' Sigma^{-1} X_exp) beta = X_exp' Sigma^{-1} Y_w
            //
            // To simplify, we work with n-vectors by contracting over the J-1 dimension.
            // The score for beta is: sum_j X' * (something per j)
            // Following POLMM: S_beta = X' * R * Psi * R * (y_tilde - mu_tilde)
            // which reduces to an n-vector operation.

            // Compute R * ytilde: R_i[j] = iRMat(i,j), multiply element-wise
            // Then Psi * R * ytilde, then R * that → gives RPsiR * something

            // Simpler approach: update beta via the score equation
            // Score_beta = X' * W * (eta + W^{-1} * (y - mu) - X*beta - b)
            // where W is diagonal working weights (collapsed to n-vector)
            // This is the standard IRLS for the mean parameters.

            // Collapse to n-vector: summing over J-1 thresholds
            // Working weight per subject: w_i = sum_j h_i[j]^2 / p_ij - (sum h)^2
            // This is RPsiR from the POLMM paper.
            Eigen::VectorXd RPsiR(n);
            for (int i = 0; i < n; ++i) {
                double sumH = 0.0, sumH2 = 0.0;
                for (int j = 0; j < Jm1; ++j) {
                    double h = psiBlocks[i].h(j);
                    sumH += h;
                    sumH2 += h * h;
                }
                RPsiR(i) = std::max(sumH2 - sumH * sumH, 1e-10);
                // Note: this is actually trace of the Psi block = sum of eigenvalues
                // which equals sum(h_j) - sum(h_j)^2 for rank-1 update
                // Actually for Psi = diag(h) - h*h':
                //   sum diag = sum(h) - sum(h^2)   <- NO
                //   trace(Psi) = sum(h_j) - sum(h_j^2)  <- NO
                //   trace(Psi) = sum(h_j) - ||h||^2  <- that's sum(h) - sum(h^2)
                // For RPsiR where R = I (no iR scaling at this point):
                //   RPsiR_i = sum_j sum_k Psi_i[j,k] = sum_j (h_j - h_j * sum(h))
                //           = sum(h) - sum(h) * sum(h) = sum(h) * (1 - sum(h))
                RPsiR(i) = std::max(sumH * (1.0 - sumH), 1e-10);
            }

            // Score for beta: X' * RPsiR * (eta + RPsiR^{-1} * RymuVec - X*beta - bVec)
            // RymuVec_i = sum_j f_j * (1{y<=j} - F(eps_j-eta_i))
            Eigen::VectorXd RymuVec(n);
            for (int i = 0; i < n; ++i) {
                double val = 0.0;
                double cumY = 0.0, cumMu = 0.0;
                for (int j = 0; j < Jm1; ++j) {
                    cumY += yMat(i, j);
                    cumMu += muMat(i, j);
                    val += psiBlocks[i].h(j) * (cumY - cumMu);
                }
                RymuVec(i) = val;
            }

            // IRLS for beta: beta_new = (X' W X)^{-1} X' W z
            // where z_i = eta_i + RPsiR_i^{-1} * RymuVec_i - bVec_i
            //       W_i = RPsiR_i
            Eigen::VectorXd z(n);
            for (int i = 0; i < n; ++i)
                z(i) = eta(i) + RymuVec(i) / RPsiR(i) - bVec(i);

            // (X' diag(RPsiR) X) beta = X' diag(RPsiR) z
            Eigen::MatrixXd XtWX = X.transpose() * RPsiR.asDiagonal() * X;
            XtWX.diagonal().array() += 1e-8; // regularize
            Eigen::VectorXd XtWz = X.transpose() * (RPsiR.array() * z.array()).matrix();
            beta = XtWX.ldlt().solve(XtWz);

            // Update bVec via BLUP:
            //   b = tau * K * Sigma^{-1} * (Y_w - X*beta)
            // Using the collapsed formulation:
            //   b = tau * K * RPsiR * r / (RPsiR + tau * diagK) approximately
            // Or use PCG for the full n*(J-1) system, but that's expensive.
            //
            // Use the simpler diagonal approximation for the inner loop:
            Eigen::VectorXd resid = z - X * beta;
            Eigen::VectorXd rhs_b(n);
            for (int i = 0; i < n; ++i)
                rhs_b(i) = tau * RPsiR(i) * resid(i);
            // b ≈ tau * K * diag(RPsiR / (RPsiR + tau*diagK)) * resid
            // Direct GRM multiply
            std::vector<double> tmp_in(n), tmp_out(n);
            for (int i = 0; i < n; ++i)
                tmp_in[i] = RPsiR(i) / (RPsiR(i) + tau * diagK[i]) * resid(i);
            grm.multiply(tmp_in.data(), tmp_out.data(), static_cast<uint32_t>(n));
            for (int i = 0; i < n; ++i)
                bVec(i) = tau * tmp_out[i];

            // Update eps via Newton-Raphson
            eta = X * beta + bVec;
            updateMuMat(X, beta, bVec, eps, n, J, muMat, iRMat);

            for (int epsIter = 0; epsIter < kMaxIterEps; ++epsIter) {
                Eigen::VectorXd epsGrad(Jm1);
                Eigen::MatrixXd epsHess = Eigen::MatrixXd::Zero(Jm1, Jm1);
                epsGrad.setZero();

                for (int i = 0; i < n; ++i) {
                    int yi = yVec(i);
                    for (int j = 0; j < Jm1; ++j) {
                        double f_j = dlogistic(eps(j) - eta(i));

                        // P(y = j) contribution
                        double p_j = clampProb(muMat(i, j));
                        double p_jp1 = (j + 1 < J) ? clampProb(muMat(i, j + 1)) : 0.0;

                        // Gradient: (1{yi==j}/p_j - 1{yi==j+1}/p_{j+1}) * f_j
                        double g = 0.0;
                        if (yi == j) g += f_j / p_j;
                        if (yi == j + 1) g -= f_j / p_jp1;
                        epsGrad(j) += g;

                        // Hessian (diagonal approximation)
                        epsHess(j, j) -= f_j * f_j / (p_j * p_jp1);
                    }
                }

                // Regularize
                epsHess.diagonal().array() -= 1e-6;
                Eigen::VectorXd epsDelta = epsHess.ldlt().solve(epsGrad);
                // Negate because we used -Hessian
                eps -= epsDelta;

                // Ensure ordering
                for (int j = 1; j < Jm1; ++j)
                    eps(j) = std::max(eps(j), eps(j - 1) + 1e-4);

                updateMuMat(X, beta, bVec, eps, n, J, muMat, iRMat);

                if (epsDelta.norm() < kTolEps) break;
            }

            // Check beta convergence
            if ((beta - beta_old).cwiseAbs().maxCoeff() < kTolBeta) break;
        }

        // Outer loop: update tau
        // tau = b' * b / (n - trace(tau * K * Sigma^{-1}))
        // Use Hutchinson's trace estimator for trace(tau * K * Sigma^{-1})

        // Simple tau estimate: tau = b'Kb / tr(K * something)
        // Use REML-style: tau_new = (b' K^{-1}_reg b) / n_eff
        // Simplified: tau = sum(b^2) / (n - p)  (like REML for LMM)
        // A better estimate: use quadratic form
        double btb = bVec.squaredNorm();

        // Hutchinson trace: trace(tau * K * Sigma^{-1}) ≈ sum of random probes
        // For now, use a simpler estimate:
        //   tau_new = bKb / (n - p)
        // This is a rough REML estimate.
        // Actually: use the update from POLMM paper:
        //   tau_new = (b' * Sigma^{-1} * K * Sigma^{-1} * b) / trace(K * Sigma^{-2})
        // which is intractable. Use the simpler:
        //   tau_new = tau * bKb / btb  (ratio scaling)

        // Even simpler Brent-style: tau_new = bKb / n
        // The proper REML uses trace estimation. Let's implement Hutchinson:
        Eigen::VectorXd eta = X * beta + bVec;
        for (int i = 0; i < n; ++i)
            psiBlocks[i].compute(eta(i), eps, Jm1);

        double traceEst = 0.0;
        {
            // Pre-generate all probe vectors for deterministic results
            // regardless of thread count.
            std::mt19937 rng(42 + outerIter);
            std::normal_distribution<double> normal(0.0, 1.0);
            const int dim = n * Jm1;
            std::vector<Eigen::VectorXd> probeVecs(kTraceNRun);
            for (int run = 0; run < kTraceNRun; ++run) {
                probeVecs[run].resize(dim);
                for (int k = 0; k < dim; ++k)
                    probeVecs[run](k) = normal(rng);
            }

            const int probeThreads = std::min(nLocalThreads, kTraceNRun);
            if (probeThreads <= 1) {
                // Serial path
                for (int run = 0; run < kTraceNRun; ++run) {
                    Eigen::VectorXd SinvU = pcgSolve(psiBlocks, grm, diagK, tau, n, Jm1, probeVecs[run]);
                    Eigen::VectorXd KSinvU(dim);
                    batchSpmvInterleaved(grm, SinvU.data(), KSinvU.data(), static_cast<uint32_t>(n), Jm1);
                    traceEst += probeVecs[run].dot(KSinvU);
                }
            } else {
                // Parallel probes: each thread runs its share of PCG solves.
                // pcgSolve is thread-safe (reads shared psiBlocks/grm/diagK,
                // all mutable state is stack-local).
                std::vector<double> probeTraces(kTraceNRun, 0.0);
                std::atomic<int> nextProbe{0};

                auto probeWork = [&]() {
                    for (;;) {
                        int run = nextProbe.fetch_add(1, std::memory_order_relaxed);
                        if (run >= kTraceNRun) break;
                        Eigen::VectorXd SinvU = pcgSolve(psiBlocks, grm, diagK, tau, n, Jm1, probeVecs[run]);
                        Eigen::VectorXd KSinvU(dim);
                        batchSpmvInterleaved(grm, SinvU.data(), KSinvU.data(), static_cast<uint32_t>(n), Jm1);
                        probeTraces[run] = probeVecs[run].dot(KSinvU);
                    }
                };

                std::vector<std::thread> threads;
                threads.reserve(probeThreads - 1);
                for (int t = 0; t < probeThreads - 1; ++t)
                    threads.emplace_back(probeWork);
                probeWork(); // caller participates
                for (auto &th : threads)
                    th.join();

                for (int run = 0; run < kTraceNRun; ++run)
                    traceEst += probeTraces[run];
            }
            traceEst *= tau / kTraceNRun;
        }

        // REML update: tau_new = b'b / (n - trace(tau*K*Sigma^{-1}))
        double neff = std::max(static_cast<double>(n) - traceEst, 1.0);
        double tau_new = btb / neff;
        tau_new = std::max(tau_new, 1e-6);

        infoMsg(
            "[%s] PQL iter %d: tau=%.6f -> %.6f, |delta|=%.2e",
            tag,
            outerIter + 1,
            tau_old,
            tau_new,
            std::abs(tau_new - tau_old) / std::max(tau_old, 1e-10)
        );

        tau = tau_new;

        if (std::abs(tau - tau_old) / std::max(tau_old, 1e-10) < kTolTau) break;
    }

    infoMsg("[%s] Null model converged. tau=%.6f", tag, tau);

    // Store results
    null.n = n;
    null.J = J;
    null.p = p;
    null.beta = beta;
    null.eps = eps;
    null.bVec = bVec;
    null.tau = tau;
    null.muMat = muMat;
    null.iRMat = iRMat;
}

// ══════════════════════════════════════════════════════════════════════
// getobjP: precompute marker-test matrices
// ══════════════════════════════════════════════════════════════════════

void computeTestMatrices(
    const Eigen::VectorXi &yVec,
    const Eigen::MatrixXd &X,
    POLMMNullModel &null,
    const char *tag =
    ""
) {
    const int n = null.n;
    const int J = null.J;
    const int Jm1 = J - 1;

    Eigen::VectorXd eta = X * null.beta + null.bVec;

    // RPsiR(i) = sum_j sum_k R_ij * Psi_i[j,k] * R_ik
    //   where R_ij = 1 / iRMat(i,j)  (i.e., iRMat stores the inverse)
    // With our Psi formulation (Psi = diag(h) - h*h'):
    //   RPsiR_i = sum_j (R_ij^2 * h_j) - (sum_j R_ij * h_j)^2
    // But R_ij = 1/iRMat(i,j) = sqrt(mu(i,j) * cumMuBar(i,j))
    // So R_ij^2 * h_j = mu(i,j) * cumMuBar(i,j) * f(eps_j - eta_i)

    // Simpler: in POLMM, RPsiR is defined as the per-subject diagonal weight
    // of the collapsed working model. From getPsixMat and getRPsiR:
    //   RPsiR_i = sum_{j=0}^{J-2} f_j^2 / mu(i,j) + f_j^2 / mu(i,j+1)
    //           = sum_{j=0}^{J-2} f_j^2 * (1/mu(i,j) + 1/mu(i,j+1))
    // Wait, let me use the approach from the actual paper/code:
    //   getPsixMat does Psi * x for the n*(J-1) system
    //   getRPsiR gives n-vector: RPsiR_i = sum_{j=0}^{J-2} h_j - (sum h_j)^2
    //   where h_j = dF(eps_j - eta_i)
    // From our PsiBlock: trace(Psi_i) = sum(h) - sum(h)^2 = sum(h)*(1-sum(h))
    // But we need the quadratic form R' Psi R, not just trace.

    // Actually from the POLMM reference:
    //   RymuVec_i = sum_{j=0}^{J-2} (y_cum_j - F_j) * f_j       (simplified)
    //   RPsiR_i = sum_{j=0}^{J-2} f_j - (sum f_j)^2              (simplified)
    // Let me compute these from first principles.

    Eigen::VectorXd RPsiR(n);
    Eigen::VectorXd RymuVec(n);

    for (int i = 0; i < n; ++i) {
        double sumF = 0.0;
        double rym = 0.0;
        double cumY = 0.0;
        for (int j = 0; j < Jm1; ++j) {
            double f_j = dlogistic(null.eps(j) - eta(i));
            sumF += f_j;
            cumY += (yVec(i) == j) ? 1.0 : 0.0;
            double cumMu_j = logistic(null.eps(j) - eta(i));
            rym += f_j * (cumY - cumMu_j);
        }
        RPsiR(i) = std::max(sumF * (1.0 - sumF), 1e-10);
        RymuVec(i) = rym;
    }

    // XR_Psi_R_new (p × n): X' * diag(RPsiR)
    Eigen::MatrixXd XR = X.transpose() * RPsiR.asDiagonal(); // p × n

    // (X' diag(RPsiR) X)^{-1} * X' diag(RPsiR) = (XtWX)^{-1} * XR
    Eigen::MatrixXd XtWX = XR * X; // p × p
    XtWX.diagonal().array() += 1e-8;
    Eigen::MatrixXd XtWX_inv = XtWX.inverse();

    // XXR_Psi_RX_new (n × p) = X * (XtWX)^{-1}
    Eigen::MatrixXd XXR = X * XtWX_inv; // n × p

    // XR_Psi_R_new (p × n) = XtWX_inv * XR
    Eigen::MatrixXd XR_new = XtWX_inv * XR; // p × n

    null.RPsiR = RPsiR;
    null.RymuVec = RymuVec;
    null.XXR_Psi_RX_new = XXR;  // n × p
    null.XR_Psi_R_new = XR_new; // p × n

    infoMsg("[%s] Test matrices computed. mean(RPsiR)=%.6f", tag, RPsiR.mean());
}

// ══════════════════════════════════════════════════════════════════════
// Re-index a GRM for a per-phenotype subject subset.
// unionToLocal maps union-dense index → per-phenotype index (UINT32_MAX = excluded).
// ══════════════════════════════════════════════════════════════════════

SparseGRM reindexGrm(
    const SparseGRM &unionGrm,
    const std::vector<uint32_t> &unionToLocal,
    uint32_t nPheno
) {
    std::vector<SparseGRM::Entry> entries;
    entries.reserve(unionGrm.nnz());
    for (const auto &e : unionGrm.entries()) {
        uint32_t lr = unionToLocal[e.row];
        uint32_t lc = unionToLocal[e.col];
        if (lr == UINT32_MAX || lc == UINT32_MAX) continue;
        // Maintain lower-triangle invariant: row >= col
        if (lr < lc) std::swap(lr, lc);
        entries.push_back({lr, lc, e.value});
    }
    return SparseGRM::fromEntries(nPheno, std::move(entries));
}

// ══════════════════════════════════════════════════════════════════════
// estimateVarianceRatios with union genotype extraction.
// Reads union-sized genotypes, extracts per-phenotype dense vector via
// localToUnion mapping, then computes variance ratios.
// ══════════════════════════════════════════════════════════════════════

void estimateVarianceRatiosFromUnion(
    POLMMNullModel &null,
    const Eigen::MatrixXd &X,
    const Eigen::VectorXi &yVec,
    const SparseGRM &grm,
    GenoMeta &genoData,
    const std::vector<uint32_t> &localToUnion,
    const char *tag = ""
) {
    const int n = null.n;
    const double tau = null.tau;
    const int nBins = static_cast<int>(kMacBins.size());

    infoMsg("[%s] Estimating MAC-binned variance ratios (%d bins, %d SNPs/bin)", tag, nBins, kNSnpsPerBin);

    const uint32_t M = genoData.nMarkers();
    const uint32_t nUnion = genoData.nSubjUsed();
    const uint32_t nPheno = static_cast<uint32_t>(n);

    std::vector<uint32_t> candidates(M);
    std::iota(candidates.begin(), candidates.end(), 0);
    std::mt19937 rng(123);
    std::shuffle(candidates.begin(), candidates.end(), rng);

    auto cursor = genoData.makeCursor();
    Eigen::VectorXd unionGVec(nUnion);
    Eigen::VectorXd GVec(nPheno);

    std::vector<std::vector<double> > binRatios(nBins);
    std::vector<int> binTarget(nBins, kNSnpsPerBin);
    int totalNeeded = nBins * kNSnpsPerBin;

    for (uint32_t ci = 0; ci < candidates.size() && totalNeeded > 0; ++ci) {
        uint32_t midx = candidates[ci];

        double altFreq, altCounts, missingRate, hweP, maf, mac;
        std::vector<uint32_t> indexForMissing;
        cursor->beginSequentialBlock(midx);
        cursor->getGenotypes(midx, unionGVec, altFreq, altCounts, missingRate, hweP, maf, mac, indexForMissing);

        for (uint32_t mi : indexForMissing)
            unionGVec(mi) = 2.0 * altFreq;

        // Extract per-phenotype genotypes
        for (uint32_t i = 0; i < nPheno; ++i)
            GVec[i] = unionGVec[localToUnion[i]];

        // Compute per-phenotype MAC
        double phenoAltSum = GVec.sum();
        double phenoMac = std::min(phenoAltSum, 2.0 * nPheno - phenoAltSum);

        if (phenoMac < kMacBins[0]) continue;

        int bin = nBins - 1;
        for (int b = 0; b < nBins - 1; ++b) {
            if (phenoMac < kMacBins[b + 1]) {
                bin = b;
                break;
            }
        }

        if (static_cast<int>(binRatios[bin].size()) >= binTarget[bin]) continue;

        Eigen::VectorXd tmpP = null.XR_Psi_R_new * GVec;
        Eigen::VectorXd adjG = GVec - null.XXR_Psi_RX_new * tmpP;

        double VarW = (null.RPsiR.array() * adjG.array().square()).sum();
        if (VarW < 1e-10) continue;

        double GKG = grm.quadForm(adjG.data(), static_cast<uint32_t>(n));
        double VarP = VarW + tau * GKG;
        double ratio = VarP / VarW;

        binRatios[bin].push_back(ratio);
        if (static_cast<int>(binRatios[bin].size()) >= binTarget[bin]) --totalNeeded;
    }

    std::vector<double> varRatios(nBins, 0.0);
    for (int b = 0; b < nBins; ++b) {
        if (!binRatios[b].empty()) {
            varRatios[b] = std::accumulate(binRatios[b].begin(), binRatios[b].end(), 0.0) / binRatios[b].size();
            infoMsg(
                "[%s] VarRatio bin %d [MAC>=%.0f]: %.6f (from %zu SNPs)",
                tag,
                b,
                kMacBins[b],
                varRatios[b],
                binRatios[b].size()
            );
        }
    }

    for (int b = 0; b < nBins; ++b) {
        if (binRatios[b].empty()) {
            double best = 0.0;
            int bestDist = nBins + 1;
            for (int ob = 0; ob < nBins; ++ob) {
                if (!binRatios[ob].empty() && std::abs(ob - b) < bestDist) {
                    bestDist = std::abs(ob - b);
                    best = varRatios[ob];
                }
            }
            varRatios[b] = (best > 0.0) ? best : 1.0;
            infoMsg("[%s] VarRatio bin %d [MAC>=%.0f]: %.6f (borrowed from neighbor)", tag, b, kMacBins[b], varRatios[b]
            );
        }
    }

    null.macBounds = kMacBins;
    null.varRatios = varRatios;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// POLMMMethod implementation

POLMMMethod::POLMMMethod(
    const POLMMNullModel &null,
    double spaCutoff
)
    : m_n(null.n),
      m_J(null.J),
      m_XXR_Psi_RX_new(null.XXR_Psi_RX_new),
      m_XR_Psi_R_new(null.XR_Psi_R_new),
      m_RymuVec(null.RymuVec),
      m_RPsiR(null.RPsiR),
      m_spaCutoff(spaCutoff),
      m_macBounds(null.macBounds),
      m_varRatios(null.varRatios),
      m_muMat(null.muMat),
      m_iRMat(null.iRMat),
      m_adjG(null.n),
      m_tmpP(null.p)
{
}

double POLMMMethod::lookupVarRatio(double mac) const {
    const int nBins = static_cast<int>(m_macBounds.size());
    int bin = nBins - 1;
    for (int b = 0; b < nBins - 1; ++b) {
        if (mac < m_macBounds[b + 1]) {
            bin = b;
            break;
        }
    }
    return m_varRatios[bin];
}

std::unique_ptr<MethodBase> POLMMMethod::clone() const {
    return std::make_unique<POLMMMethod>(*this);
}

void POLMMMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int /*markerInChunkIdx*/,
    std::vector<double> &result
) {
    result.clear();
    result.resize(4, std::numeric_limits<double>::quiet_NaN());

    // Flip allele if altFreq > 0.5 for SPA numerical stability
    bool flipped = (altFreq > 0.5);
    if (flipped) {
        GVec.array() = 2.0 - GVec.array();
        altFreq = 1.0 - altFreq;
    }

    // Covariate adjustment: adjG = G - XXR * (XR_new * G)
    m_tmpP.noalias() = m_XR_Psi_R_new * GVec;
    m_adjG.noalias() = GVec - m_XXR_Psi_RX_new * m_tmpP;

    // Score statistic
    double Stat = m_adjG.dot(m_RymuVec);

    // Working variance
    double VarW = (m_RPsiR.array() * m_adjG.array().square()).sum();
    if (VarW < 1e-20) return;

    // Optimization 4: MAC-binned variance ratio
    double mac = std::min(altFreq, 1.0 - altFreq) * 2.0 * m_n;
    double varRatio = lookupVarRatio(mac);
    double VarP = VarW * varRatio;

    // Z-score and normal approximation
    double Z = Stat / std::sqrt(VarP);
    double absZ = std::abs(Z);

    double pval;
    if (absZ < m_spaCutoff) {
        pval = 2.0 * 0.5 * std::erfc(absZ / M_SQRT2);
    } else {
        // Optimization 2: use binary SPA fast path when J=2
        pval = (m_J == 2) ? spaTestBinary(Stat, VarW) : spaTest(Stat, VarW);
        if (std::isnan(pval)) {
            // Fallback to normal
            pval = 2.0 * 0.5 * std::erfc(absZ / M_SQRT2);
        }
    }

    // Effect size
    double BETA = Stat / VarP;
    double SE = 1.0 / std::sqrt(VarP);

    // Flip sign back if allele was flipped
    if (flipped) {
        Z = -Z;
        BETA = -BETA;
    }

    result[0] = pval;
    result[1] = Z;
    result[2] = BETA;
    result[3] = SE;
}

// ══════════════════════════════════════════════════════════════════════
// SPA for ordinal phenotype
// ══════════════════════════════════════════════════════════════════════

// CGF K(t) for the ordinal distribution with allele-frequency weighting
// Uses the "partial normal approximation": G==0 subjects contribute a
// normal component, G!=0 subjects contribute the exact ordinal CGF.

double POLMMMethod::K0(
    double t,
    const double *cMat,
    const double *muSub,
    int nSub,
    int Jm1,
    double Ratio0,
    double m1
)
const {
    // K(t) = 0.5 * Ratio0 * t^2  (normal part from G==0)
    //      + sum_i log(1 + sum_j mu_ij * (exp(c_ij * t) - 1))  (ordinal part)
    //      - m1 * t
    double K = 0.5 * Ratio0 * t * t - m1 * t;
    for (int i = 0; i < nSub; ++i) {
        double s = 0.0;
        for (int j = 0; j < Jm1; ++j) {
            double c = cMat[i * Jm1 + j];
            double mu = muSub[i * (Jm1 + 1) + j]; // muMat columns 0..J-2
            s += mu * (std::exp(c * t) - 1.0);
        }
        // Add the last category (J-1): mu_{J-1} has c=0, so contribution = 0
        // Actually last category doesn't appear in threshold-based cMat
        K += std::log(1.0 + s);
    }
    return K;
}

double POLMMMethod::K1(
    double t,
    const double *cMat,
    const double *muSub,
    int nSub,
    int Jm1,
    double Ratio0,
    double m1
)
const {
    // K'(t) = Ratio0 * t - m1
    //       + sum_i [sum_j mu_ij * c_ij * exp(c_ij*t)] / [1 + sum_j mu_ij*(exp(c_ij*t)-1)]
    double K1val = Ratio0 * t - m1;
    for (int i = 0; i < nSub; ++i) {
        double num = 0.0, den = 1.0;
        for (int j = 0; j < Jm1; ++j) {
            double c = cMat[i * Jm1 + j];
            double mu = muSub[i * (Jm1 + 1) + j];
            double ect = std::exp(c * t);
            num += mu * c * ect;
            den += mu * (ect - 1.0);
        }
        K1val += num / den;
    }
    return K1val;
}

double POLMMMethod::K2(
    double t,
    const double *cMat,
    const double *muSub,
    int nSub,
    int Jm1,
    double Ratio0
) const {
    // K''(t) = Ratio0
    //        + sum_i [(sum mu*c^2*exp(ct))*den - (sum mu*c*exp(ct))^2] / den^2
    double K2val = Ratio0;
    for (int i = 0; i < nSub; ++i) {
        double s1 = 0.0, s2 = 0.0, den = 1.0;
        for (int j = 0; j < Jm1; ++j) {
            double c = cMat[i * Jm1 + j];
            double mu = muSub[i * (Jm1 + 1) + j];
            double ect = std::exp(c * t);
            s1 += mu * c * ect;
            s2 += mu * c * c * ect;
            den += mu * (ect - 1.0);
        }
        K2val += (s2 * den - s1 * s1) / (den * den);
    }
    return K2val;
}

double POLMMMethod::spaTest(
    double Stat,
    double VarW
) const {
    const int n = m_n;
    const int J = m_J;
    const int Jm1 = J - 1;

    // Identify subjects with non-zero adjusted genotype
    std::vector<int> nzIdx;
    nzIdx.reserve(n);
    double VarW0 = 0.0;
    for (int i = 0; i < n; ++i) {
        if (std::abs(m_adjG(i)) > 1e-10) {
            nzIdx.push_back(i);
        } else {
            VarW0 += m_RPsiR(i) * m_adjG(i) * m_adjG(i);
        }
    }
    int nSub = static_cast<int>(nzIdx.size());
    if (nSub == 0) return std::numeric_limits<double>::quiet_NaN();

    double Ratio0 = VarW0 / VarW;
    double sqrtVarW = std::sqrt(VarW);

    // Build cMat and muSub for non-zero subjects
    // cMat[i * Jm1 + j] = adjG[nzIdx[i]] / iRMat(nzIdx[i], j) / sqrtVarW
    // muSub[i * (Jm1+1) + j] = muMat(nzIdx[i], j)  (J columns: 0..J-1)
    std::vector<double> cMat(nSub * Jm1);
    std::vector<double> muSub(nSub * (Jm1 + 1));

    double m1 = 0.0; // = sum_i sum_j muMat(i,j) * cMat(i,j)

    for (int si = 0; si < nSub; ++si) {
        int i = nzIdx[si];
        for (int j = 0; j < Jm1; ++j) {
            double c = m_adjG(i) / m_iRMat(i, j) / sqrtVarW;
            cMat[si * Jm1 + j] = c;
            double mu = m_muMat(i, j);
            muSub[si * (Jm1 + 1) + j] = mu;
            m1 += mu * c;
        }
        // Last category
        muSub[si * (Jm1 + 1) + Jm1] = m_muMat(i, J - 1);
    }

    double qNorm = Stat / sqrtVarW;

    // Solve for both tails
    auto solveTail = [&](double target) -> double {
        // Newton-Raphson: find zeta such that K'(zeta) = target
        double zeta = target;                  // initial guess
        for (int iter = 0; iter < 100; ++iter) {
            double k1 = K1(zeta, cMat.data(), muSub.data(), nSub, Jm1, Ratio0, m1);
            double k2 = K2(zeta, cMat.data(), muSub.data(), nSub, Jm1, Ratio0);
            if (std::abs(k2) < 1e-30) break;
            double delta = (k1 - target) / k2;
            zeta -= delta;
            if (std::abs(delta) < 1e-6) break;
        }
        return zeta;
    };

    auto tailProb = [&](double zeta, double target) -> double {
        // Lugannani-Rice formula
        double k0 = K0(zeta, cMat.data(), muSub.data(), nSub, Jm1, Ratio0, m1);
        double k2 = K2(zeta, cMat.data(), muSub.data(), nSub, Jm1, Ratio0);

        double w2 = 2.0 * (zeta * target - k0);
        if (w2 < 0.0) return std::numeric_limits<double>::quiet_NaN();
        double w = std::copysign(std::sqrt(w2), zeta);
        double v = zeta * std::sqrt(k2);

        if (std::abs(w) < 1e-10 || std::abs(v) < 1e-10) return std::numeric_limits<double>::quiet_NaN();

        // Lugannani-Rice: P ≈ Phi(-w) + phi(w)*(1/w - 1/v)
        // where Phi = normal CDF, phi = normal PDF
        double phi_w = std::exp(-0.5 * w * w) / std::sqrt(2.0 * M_PI);
        double Phi_w = 0.5 * std::erfc(w / M_SQRT2);                 // Phi(-w) = 0.5*erfc(w/sqrt2)
        // Wait: Phi_w should be the LOWER tail for the negative side
        // For p_upper: P(Z > w) = 1 - Phi(w) = 0.5 * erfc(w/sqrt2)
        // For p_lower: P(Z < -w) = 0.5 * erfc(w/sqrt2)  (same by symmetry)
        double p_tail = Phi_w + phi_w * (1.0 / w - 1.0 / v);

        return std::max(p_tail, 0.0);
    };

    // Upper tail: target = qNorm
    double zeta_upper = solveTail(qNorm);
    double p_upper = tailProb(zeta_upper, qNorm);

    // Lower tail: target = -qNorm
    double zeta_lower = solveTail(-qNorm);
    double p_lower = tailProb(zeta_lower, -qNorm);

    double pval = std::numeric_limits<double>::quiet_NaN();
    if (!std::isnan(p_upper) && !std::isnan(p_lower))pval = p_upper + p_lower;
    else if (!std::isnan(p_upper))pval = 2.0 * p_upper;
    else if (!std::isnan(p_lower))pval = 2.0 * p_lower;

    return std::min(pval, 1.0);
}

// ══════════════════════════════════════════════════════════════════════
// Optimization 2: K=2 binomial SPA fast path
// For binary phenotype (J=2), the ordinal CGF simplifies to the standard
// binomial CGF:  K(t) = sum_i log(1 - mu_i + mu_i * exp(g_i * t))
// where mu_i = P(y=0) and g_i = adjG(i) / iRMat(i,0) / sqrt(VarW)
// ══════════════════════════════════════════════════════════════════════

double POLMMMethod::spaTestBinary(
    double Stat,
    double VarW
) const {
    const int n = m_n;

    // Identify subjects with non-zero adjusted genotype
    std::vector<int> nzIdx;
    nzIdx.reserve(n);
    double VarW0 = 0.0;
    for (int i = 0; i < n; ++i) {
        if (std::abs(m_adjG(i)) > 1e-10) {
            nzIdx.push_back(i);
        } else {
            VarW0 += m_RPsiR(i) * m_adjG(i) * m_adjG(i);
        }
    }
    int nSub = static_cast<int>(nzIdx.size());
    if (nSub == 0) return std::numeric_limits<double>::quiet_NaN();

    double Ratio0 = VarW0 / VarW;
    double sqrtVarW = std::sqrt(VarW);

    // Precompute per-subject scalars: mu and g
    std::vector<double> muVec(nSub), gVec(nSub);
    double m1 = 0.0;
    for (int si = 0; si < nSub; ++si) {
        int i = nzIdx[si];
        muVec[si] = m_muMat(i, 0); // P(y = 0)
        gVec[si] = m_adjG(i) / m_iRMat(i, 0) / sqrtVarW;
        m1 += muVec[si] * gVec[si];
    }

    // Binomial CGF and derivatives (partial normal approximation)
    auto K0bin = [&](double t) -> double {
        double K = 0.5 * Ratio0 * t * t - m1 * t;
        for (int si = 0; si < nSub; ++si) {
            double egt = std::exp(gVec[si] * t);
            K += std::log(1.0 - muVec[si] + muVec[si] * egt);
        }
        return K;
    };

    auto K1bin = [&](double t) -> double {
        double K1 = Ratio0 * t - m1;
        for (int si = 0; si < nSub; ++si) {
            double egt = std::exp(gVec[si] * t);
            double den = 1.0 - muVec[si] + muVec[si] * egt;
            K1 += muVec[si] * gVec[si] * egt / den;
        }
        return K1;
    };

    auto K2bin = [&](double t) -> double {
        double K2 = Ratio0;
        for (int si = 0; si < nSub; ++si) {
            double egt = std::exp(gVec[si] * t);
            double den = 1.0 - muVec[si] + muVec[si] * egt;
            double g2 = gVec[si] * gVec[si];
            K2 += muVec[si] * (1.0 - muVec[si]) * g2 * egt / (den * den);
        }
        return K2;
    };

    double qNorm = Stat / sqrtVarW;

    // Newton-Raphson to find saddlepoint
    auto solveTail = [&](double target) -> double {
        double zeta = target;
        for (int iter = 0; iter < 100; ++iter) {
            double k1 = K1bin(zeta);
            double k2 = K2bin(zeta);
            if (std::abs(k2) < 1e-30) break;
            double delta = (k1 - target) / k2;
            zeta -= delta;
            if (std::abs(delta) < 1e-6) break;
        }
        return zeta;
    };

    // Lugannani-Rice tail probability
    auto tailProb = [&](double zeta, double target) -> double {
        double k0 = K0bin(zeta);
        double k2 = K2bin(zeta);

        double w2 = 2.0 * (zeta * target - k0);
        if (w2 < 0.0) return std::numeric_limits<double>::quiet_NaN();
        double w = std::copysign(std::sqrt(w2), zeta);
        double v = zeta * std::sqrt(k2);

        if (std::abs(w) < 1e-10 || std::abs(v) < 1e-10) return std::numeric_limits<double>::quiet_NaN();

        double phi_w = std::exp(-0.5 * w * w) / std::sqrt(2.0 * M_PI);
        double Phi_w = 0.5 * std::erfc(w / M_SQRT2);
        double p_tail = Phi_w + phi_w * (1.0 / w - 1.0 / v);
        return std::max(p_tail, 0.0);
    };

    double zeta_upper = solveTail(qNorm);
    double p_upper = tailProb(zeta_upper, qNorm);

    double zeta_lower = solveTail(-qNorm);
    double p_lower = tailProb(zeta_lower, -qNorm);

    double pval = std::numeric_limits<double>::quiet_NaN();
    if (!std::isnan(p_upper) && !std::isnan(p_lower))pval = p_upper + p_lower;
    else if (!std::isnan(p_upper))pval = 2.0 * p_upper;
    else if (!std::isnan(p_lower))pval = 2.0 * p_lower;

    return std::min(pval, 1.0);
}

// ══════════════════════════════════════════════════════════════════════
// runPOLMM — top-level orchestrator
// ══════════════════════════════════════════════════════════════════════

void runPOLMM(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const GenoSpec &geno,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const int K = static_cast<int>(phenoNames.size());

    // ── 1. Load phenotype/covariate data (union mask) ───────────────
    infoMsg("POLMM: Loading phenotype and covariate data (%d phenotypes)", K);
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    sd.loadPhenoFile(phenoFile);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();

    const uint32_t N_union = sd.nUsed();
    infoMsg("Subjects after intersection (union): %u", N_union);

    // ── 2. Extract phenotype values and build per-phenotype mappings ─
    Eigen::MatrixXd phenoMat(N_union, K);
    sd.fillColumnsInto(phenoNames, phenoMat);

    // Extract covariates once (union-sized)
    Eigen::MatrixXd unionCovar;
    if (!covarNames.empty()) {
        unionCovar.resize(N_union, static_cast<Eigen::Index>(covarNames.size()));
        sd.fillColumnsInto(covarNames, unionCovar);
    } else if (sd.hasCovar()) {
        unionCovar = sd.covar();
    }
    const int nCovar = static_cast<int>(unionCovar.cols());

    // Per-phenotype mappings
    struct PerPhenoMap {
        std::vector<uint32_t> unionToLocal; // UINT32_MAX = excluded
        std::vector<uint32_t> localToUnion;
        uint32_t nUsed;
    };

    std::vector<PerPhenoMap> phenoMap(K);
    for (int k = 0; k < K; ++k) {
        phenoMap[k].unionToLocal.resize(N_union, UINT32_MAX);
        uint32_t n = 0;
        for (uint32_t i = 0; i < N_union; ++i) {
            if (!std::isnan(phenoMat(i, k))) {
                phenoMap[k].unionToLocal[i] = n;
                phenoMap[k].localToUnion.push_back(i);
                ++n;
            }
        }
        phenoMap[k].nUsed = n;
        infoMsg("  Phenotype '%s': %u subjects (non-NA)", phenoNames[k].c_str(), n);
        if (n == 0) throw std::runtime_error("POLMM: phenotype '" + phenoNames[k] + "' has no non-missing subjects");
    }

    // ── 3. Build GenoData on union mask ─────────────────────────────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), N_union, nSnpPerChunk);

    // ── 4. Load sparse GRM in union ordering ────────────────────────
    auto usedIIDs = sd.usedIIDs();
    SparseGRM unionGrm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, usedIIDs, sd.famIIDs());
    infoMsg("Sparse GRM: %zu entries, %u subjects", unionGrm.nnz(), unionGrm.nSubjects());

    // ── 5. Per-phenotype: fit null model + variance ratios (parallel) ─
    std::vector<PhenoTask> tasks(K);
    const int nFitThreads = std::min(nthreads, K);
    // When K < T, each worker gets multiple threads for intra-phenotype
    // parallelism (Hutchinson probes).  Total ≤ T (no oversubscription).
    const int nLocalThreads = std::max(1, nthreads / std::max(K, 1));
    infoMsg("POLMM: Fitting %d null models with %d threads (%d threads/phenotype)",
            K, nFitThreads, nLocalThreads);

    std::atomic<int> nextK{0};
    std::vector<std::string> fitErrors(K);

    auto fitWorker = [&]() {
        for (;;) {
            int k = nextK.fetch_add(1, std::memory_order_relaxed);
            if (k >= K) break;

            try {
                const std::string &name = phenoNames[k];
                const uint32_t n_k = phenoMap[k].nUsed;
                const auto &localToUnion = phenoMap[k].localToUnion;
                infoMsg("[%s] Fitting null model (%u subjects)", name.c_str(), n_k);

                // Extract per-phenotype dense Y
                Eigen::VectorXd yRaw(n_k);
                for (uint32_t i = 0; i < n_k; ++i)
                    yRaw[i] = phenoMat(localToUnion[i], k);

                // Auto-detect ordinal levels: sort unique values, remap to 0..J-1
                std::set<int> levels;
                for (uint32_t i = 0; i < n_k; ++i)
                    levels.insert(static_cast<int>(std::round(yRaw[i])));
                std::vector<int> levelVec(levels.begin(), levels.end());

                const int J = static_cast<int>(levelVec.size());
                if (J < 2) throw std::runtime_error("POLMM: ordinal phenotype '" + name + "' must have >= 2 levels");

                for (int j = 0; j < J; ++j)
                    infoMsg("[%s] Level %d -> category %d", name.c_str(), levelVec[j], j);

                Eigen::VectorXi yVec(n_k);
                for (uint32_t i = 0; i < n_k; ++i) {
                    int val = static_cast<int>(std::round(yRaw[i]));
                    auto it = std::lower_bound(levelVec.begin(), levelVec.end(), val);
                    yVec[i] = static_cast<int>(it - levelVec.begin());
                }

                infoMsg("[%s] %u subjects, %d ordinal levels", name.c_str(), n_k, J);

                // Build per-phenotype covariate matrix (intercept + user covariates)
                const int p = 1 + nCovar;
                Eigen::MatrixXd X(n_k, p);
                X.col(0).setOnes();
                if (nCovar > 0) {
                    for (uint32_t i = 0; i < n_k; ++i)
                        X.row(i).tail(nCovar) = unionCovar.row(localToUnion[i]);
                }

                // Re-index GRM for per-phenotype ordering
                SparseGRM phenoGrm = reindexGrm(unionGrm, phenoMap[k].unionToLocal, n_k);
                infoMsg("[%s] Per-phenotype GRM: %zu entries", name.c_str(), phenoGrm.nnz());

                // Fit null model
                POLMMNullModel null;
                fitNullModel(yVec, X, phenoGrm, static_cast<int>(n_k), J, p, null, name.c_str(), nLocalThreads);
                computeTestMatrices(yVec, X, null, name.c_str());

                // Estimate MAC-binned variance ratios (using union genotypes + extraction)
                estimateVarianceRatiosFromUnion(null, X, yVec, phenoGrm, *genoData, localToUnion, name.c_str());

                // Build PhenoTask
                tasks[k].phenoName = name;
                tasks[k].method = std::make_unique<POLMMMethod>(null, spaCutoff);
                tasks[k].unionToLocal = phenoMap[k].unionToLocal;
                tasks[k].nUsed = n_k;
            } catch (const std::exception &ex) {
                fitErrors[k] = ex.what();
            }
        }
    };

    {
        std::vector<std::thread> threads;
        threads.reserve(nFitThreads);
        for (int t = 0; t < nFitThreads; ++t)
            threads.emplace_back(fitWorker);
        for (auto &th : threads)
            th.join();
    }

    // Check for errors
    for (int k = 0; k < K; ++k) {
        if (!fitErrors[k].empty())
            throw std::runtime_error("POLMM: phenotype '" + phenoNames[k] + "' failed: " + fitErrors[k]);
    }

    // ── 6. Run multi-phenotype engine ───────────────────────────────
    infoMsg("POLMM: Starting multi-phenotype association (%d phenotypes, %d threads)", K, nthreads);
    multiPhenoEngine(
        *genoData,
        tasks,
        outPrefix,
        "POLMM",
        compression,
        compressionLevel,
        nthreads,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        hweCutoff
    );
}
