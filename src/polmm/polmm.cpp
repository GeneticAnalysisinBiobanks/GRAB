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

// Recommended values: aligned with GRAB 0.2.4 reference defaults
// (see examples_1kg/vs_grab024/GRAB/R/POLMM.R::checkControl.NullModel.POLMM).
constexpr int kMaxIterPQL = 100;
constexpr double kTolBeta = 0.001;
constexpr double kTolTau = 0.002;
constexpr int kMaxIterEps = 100;
constexpr double kTolEps = 1e-10;
constexpr int kMaxIterPCG = 100;
constexpr double kTolPCG = 1e-6;
constexpr int kTraceNRun = 30;
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
// Reference-style POLMM null-model fit (faithful port of GRAB 0.2.4
// implementation in examples_1kg/vs_grab024/GRAB/src/POLMM.cpp).
//
// Convention (per-subject working covariance via multinomial Psi):
//   nuMat(i,j) = logistic(eps_j - eta_i)                  cumulative
//   muMat(i,j) = nuMat_j - nuMat_{j-1}                    category prob
//   mMat(i,j)  = nuMat_j + nuMat_{j-1} - 1
//   WMat(i,j)  = nuMat_j (1 - nuMat_j)
//   iRMat(i,j) = 1 / (mMat(i,j) - mMat(i, J-1))           signed
//   YMat(i,j)  = eta_i + iRMat(i,j) · iPsi[(yMat - muMat)_{0..J-2}](i,j)
//
// Sigma operator (acts on n × (J-1) matrices or n*(J-1) vectors with
// row-major flattening idx = i*(J-1) + j):
//   Sigma x = iRMat ⊙ Psi^{-1}(iRMat ⊙ x) + tau · tZ(K · Z x)
// where iPsi(z)(i,j) = sum_k z(i,k) / muMat(i,J-1) + z(i,j) / muMat(i,j);
// Z x produces an n-vector by per-subject row-sum;
// tZ v duplicates each n-entry across (J-1) slots.
//
// Tau update: AI-REML on the marginal log-likelihood.  Hutchinson
// trace estimator uses Rademacher probes precomputed once at setup
// (TraceRandMat and V_TRM = tZ(K · Z(TraceRandMat))).
//
// eps(0) is constrained to 0; the intercept of beta absorbs the shift.
// ──────────────────────────────────────────────────────────────────────
// Layout helpers: vector ↔ matrix conversion in row-major flattening.
//   vec[i*Jm1 + j] = mat(i, j)
// ──────────────────────────────────────────────────────────────────────

inline Eigen::VectorXd matToVec(
    const Eigen::MatrixXd &mat,
    int n,
    int Jm1
) {
    Eigen::VectorXd v(n * Jm1);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < Jm1; ++j)
            v(i * Jm1 + j) = mat(i, j);
    return v;
}

inline void vecToMat(
    const Eigen::VectorXd &v,
    int n,
    int Jm1,
    Eigen::MatrixXd &out
) {
    out.resize(n, Jm1);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < Jm1; ++j)
            out(i, j) = v(i * Jm1 + j);
}

// Expand n × p covariate matrix to n*(J-1) × p by row-replication.
inline Eigen::MatrixXd buildExpandedX(
    const Eigen::MatrixXd &X,
    int n,
    int p,
    int Jm1
) {
    Eigen::MatrixXd Xe(n * Jm1, p);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < Jm1; ++j)
            Xe.row(i * Jm1 + j) = X.row(i);
    return Xe;
}

// Per-subject row-sum: for v of length n*(J-1), w(i) = sum_j v(i*(J-1)+j).
inline Eigen::VectorXd ZsumVec(
    const Eigen::VectorXd &v,
    int n,
    int Jm1
) {
    Eigen::VectorXd w(n);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        for (int j = 0; j < Jm1; ++j) s += v(i * Jm1 + j);
        w(i) = s;
    }
    return w;
}

// y indicator matrix: y(i, yVec(i)) = 1.
inline Eigen::MatrixXd makeYMat(
    const Eigen::VectorXi &yVec,
    int n,
    int J
) {
    Eigen::MatrixXd Y = Eigen::MatrixXd::Zero(n, J);
    for (int i = 0; i < n; ++i) Y(i, yVec(i)) = 1.0;
    return Y;
}

// ──────────────────────────────────────────────────────────────────────
// fitCLM_initial: cumulative-logit fit, mirrors R's ordinal::clm output.
// Produces (beta_raw, eps_raw) with logit P(y<=j) = eps_raw_j - X*beta_raw.
// ──────────────────────────────────────────────────────────────────────

void fitCLM_initial(
    const Eigen::VectorXi &yVec,
    const Eigen::MatrixXd &X,
    int n,
    int J,
    int p,
    Eigen::VectorXd &beta,
    Eigen::VectorXd &eps
) {
    const int Jm1 = J - 1;
    eps.resize(Jm1);
    beta.setZero(p);

    Eigen::VectorXd cumProp = Eigen::VectorXd::Zero(J);
    for (int i = 0; i < n; ++i) cumProp(yVec(i)) += 1.0;
    cumProp /= static_cast<double>(n);
    double cs = 0.0;
    for (int j = 0; j < Jm1; ++j) {
        cs += cumProp(j);
        double q = std::min(0.999, std::max(0.001, cs));
        eps(j) = std::log(q / (1.0 - q));
    }

    const int nTheta = Jm1 + p;
    constexpr int maxIter = 50;
    for (int it = 0; it < maxIter; ++it) {
        Eigen::VectorXd eta = X * beta;
        Eigen::VectorXd grad = Eigen::VectorXd::Zero(nTheta);
        Eigen::MatrixXd H = Eigen::MatrixXd::Zero(nTheta, nTheta);

        for (int i = 0; i < n; ++i) {
            int yi = yVec(i);
            double F_yi = (yi < Jm1) ? logistic(eps(yi) - eta(i)) : 1.0;
            double F_yim1 = (yi > 0) ? logistic(eps(yi - 1) - eta(i)) : 0.0;
            double p_yi = clampProb(F_yi - F_yim1);
            double f_yi = (yi < Jm1) ? dlogistic(eps(yi) - eta(i)) : 0.0;
            double f_yim1 = (yi > 0) ? dlogistic(eps(yi - 1) - eta(i)) : 0.0;
            double dLogPdEta = -(f_yi - f_yim1) / p_yi;

            Eigen::VectorXd score = Eigen::VectorXd::Zero(nTheta);
            if (yi < Jm1) score(yi) = f_yi / p_yi;
            if (yi > 0) score(yi - 1) -= f_yim1 / p_yi;
            score.tail(p) = dLogPdEta * X.row(i).transpose();
            grad += score;
            H += score * score.transpose();
        }
        H.diagonal().array() += 1e-6;
        Eigen::VectorXd delta = H.ldlt().solve(grad);
        for (int j = 0; j < Jm1; ++j) eps(j) += delta(j);
        beta += delta.tail(p);
        if (delta.norm() < 1e-7) break;
    }
    for (int j = 1; j < Jm1; ++j)
        eps(j) = std::max(eps(j), eps(j - 1) + 0.01);
}

// ══════════════════════════════════════════════════════════════════════
// POLMMFitter — PQL / AI-REML null-model fit per reference algorithm.
// ══════════════════════════════════════════════════════════════════════

class POLMMFitter {
  public:
    POLMMFitter(
        const Eigen::VectorXi &yVec,
        const Eigen::MatrixXd &X,
        const SparseGRM &grm,
        int n,
        int J,
        int p,
        const char *tag
    )
        : m_yVec(yVec),
          m_X(X),
          m_grm(grm),
          m_n(n),
          m_J(J),
          m_p(p),
          m_Jm1(J - 1),
          m_tag(tag) {
        m_diagK = grm.diagonal();
        m_yMat = makeYMat(yVec, n, J);
        m_X_exp = buildExpandedX(X, n, p, m_Jm1);

        // CLM seed, then shift such that eps(0) = 0 (intercept absorbs).
        Eigen::VectorXd beta_raw, eps_raw;
        fitCLM_initial(yVec, X, n, J, p, beta_raw, eps_raw);
        const double shift = eps_raw(0);
        m_beta = beta_raw;
        m_beta(0) -= shift;
        m_eps = eps_raw;
        for (int j = 0; j < m_Jm1; ++j) m_eps(j) -= shift;
        m_eps(0) = 0.0;

        m_bVec = Eigen::VectorXd::Zero(n);
        m_tau = 0.2;

        // Rademacher trace probes and V_TRM = tZ(K · Z(probes)) — fixed for the fit.
        buildTraceRandMat();
    }

    void fit(POLMMNullModel &out) {
        infoMsg(
            "[%s] Fitting POLMM null model (n=%d, J=%d, p=%d, tau0=%.3f)",
            m_tag,
            m_n,
            m_J,
            m_p,
            m_tau
        );

        updateMats();
        for (m_iter = 0; m_iter < kMaxIterPQL; ++m_iter) {
            const double tau_old = m_tau;
            updateParaConv();
            updateTau();
            if (std::isnan(m_tau)) {
                infoMsg("[%s] tau became NaN; aborting fit", m_tag);
                break;
            }
            const double rel = std::abs(m_tau - tau_old)
                               / (std::abs(m_tau) + std::abs(tau_old) + kTolTau);
            infoMsg(
                "[%s] PQL iter %d: tau %.6f -> %.6f (rel=%.2e)",
                m_tag,
                m_iter + 1,
                tau_old,
                m_tau,
                rel
            );
            if (rel < kTolTau) {
                ++m_iter;
                break;
            }
        }
        infoMsg("[%s] Null model converged. tau=%.6f after %d iters", m_tag, m_tau, m_iter);

        out.n = m_n;
        out.J = m_J;
        out.p = m_p;
        out.beta = m_beta;
        out.eps = m_eps;
        out.bVec = m_bVec;
        out.tau = m_tau;
        out.muMat = m_muMat;
        out.iRMat = m_iRMat;
    }

    // Expose state for variance-ratio computation.
    const Eigen::MatrixXd &muMat() const {
        return m_muMat;
    }

    const Eigen::MatrixXd &iRMat() const {
        return m_iRMat;
    }

    const Eigen::MatrixXd &X_exp() const {
        return m_X_exp;
    }

    double tau() const {
        return m_tau;
    }

    int n() const {
        return m_n;
    }

    int J() const {
        return m_J;
    }

    int p() const {
        return m_p;
    }

    int Jm1() const {
        return m_Jm1;
    }

    // Build preconditioner at the current tau/iRMat.
    void buildPrecondAtCurrent() {
        buildPreconditioner();
    }

    // After fit(), recompute PCG-derived state at the converged tau so that
    // external pcgSolveVec / iSigmaX_XSigmaX / X_exp are consistent for the
    // variance-ratio step.
    void prepareForVarianceRatio() {
        buildPreconditioner();
        pcgSolveColumns(m_X_exp, m_iSigma_CovaMat);
        Eigen::MatrixXd XtiSX = m_X_exp.transpose() * m_iSigma_CovaMat;
        XtiSX.diagonal().array() += 1e-12;
        m_iSigmaX_XSigmaX = m_iSigma_CovaMat * XtiSX.inverse();
    }

    // PCG-solve Σ x = rhs (n × Jm1 form).  Requires preconditioner built.
    void pcgSolveMat(
        const Eigen::MatrixXd &rhs,
        Eigen::MatrixXd &x
    ) const {
        pcgSolveMatInternal(rhs, x);
    }

    // PCG-solve Σ x = rhs (n*(J-1) vector form).  Requires preconditioner built.
    Eigen::VectorXd pcgSolveVec(const Eigen::VectorXd &rhs_v) const {
        Eigen::MatrixXd rhs, x;
        vecToMat(rhs_v, m_n, m_Jm1, rhs);
        pcgSolveMatInternal(rhs, x);
        return matToVec(x, m_n, m_Jm1);
    }

    // Cached: (X' Σ^{-1} X)^{-1} and iSigmaX_XSigmaX = Σ^{-1}X · (X'Σ^{-1}X)^{-1}.
    // Available after a successful updatePara/updateTau call.
    const Eigen::MatrixXd &iSigmaX_XSigmaX() const {
        return m_iSigmaX_XSigmaX;
    }

    const Eigen::MatrixXd &iSigma_CovaMat() const {
        return m_iSigma_CovaMat;
    }

  private:
    // Inputs (references, lifetime owned by caller).
    const Eigen::VectorXi &m_yVec;
    const Eigen::MatrixXd &m_X;
    const SparseGRM &m_grm;
    int m_n, m_J, m_p, m_Jm1;
    const char *m_tag;

    std::vector<double> m_diagK;
    Eigen::MatrixXd m_yMat;  // n × J indicator
    Eigen::MatrixXd m_X_exp; // n*(J-1) × p

    // Random probes and their V-image (fixed for the fit).
    Eigen::MatrixXd m_TraceRandMat; // n*(J-1) × tracenrun
    Eigen::MatrixXd m_V_TRM;        // n*(J-1) × tracenrun

    // Parameters.
    Eigen::VectorXd m_beta; // p
    Eigen::VectorXd m_eps;  // J-1
    Eigen::VectorXd m_bVec; // n
    double m_tau = 0.0;

    // Working matrices.
    Eigen::VectorXd m_eta; // n
    Eigen::MatrixXd m_muMat, m_nuMat, m_mMat, m_WMat; // n × J
    Eigen::MatrixXd m_iRMat; // n × (J-1)
    Eigen::MatrixXd m_YMat;  // n × (J-1)

    // Preconditioner (per-subject (J-1)×(J-1) inverse blocks).
    std::vector<Eigen::MatrixXd> m_InvBlock;

    // PCG-derived cached values.
    Eigen::MatrixXd m_iSigma_CovaMat;  // n*(J-1) × p
    Eigen::MatrixXd m_iSigmaX_XSigmaX; // n*(J-1) × p

    int m_iter = 0;

    // ── core math ─────────────────────────────────────────────────

    void updateMats() {
        const int Jm1 = m_Jm1, J = m_J, n = m_n;
        m_eta = m_X * m_beta + m_bVec;
        m_muMat.resize(n, J);
        m_nuMat.resize(n, J);
        m_mMat.resize(n, J);
        m_WMat.resize(n, J);
        m_iRMat.resize(n, Jm1);

        for (int i = 0; i < n; ++i) {
            double nu0 = 0.0;
            for (int j = 0; j < Jm1; ++j) {
                const double nu1 = logistic(m_eps(j) - m_eta(i));
                m_muMat(i, j) = std::max(nu1 - nu0, 1e-20);
                m_WMat(i, j) = nu1 * (1.0 - nu1);
                m_mMat(i, j) = nu1 + nu0 - 1.0;
                m_nuMat(i, j) = nu1;
                nu0 = nu1;
            }
            const double nu1 = 1.0;
            m_muMat(i, Jm1) = std::max(nu1 - nu0, 1e-20);
            m_WMat(i, Jm1) = nu1 * (1.0 - nu1);
            m_mMat(i, Jm1) = nu1 + nu0 - 1.0;
            m_nuMat(i, Jm1) = nu1;
        }

        for (int i = 0; i < n; ++i) {
            const double mJ = m_mMat(i, Jm1);
            for (int j = 0; j < Jm1; ++j) {
                double denom = m_mMat(i, j) - mJ;
                if (std::abs(denom) < 1e-20) denom = std::copysign(1e-20, denom == 0.0 ? -1.0 : denom);
                m_iRMat(i, j) = 1.0 / denom;
            }
        }

        // YMat = eta + iRMat ⊙ iPsi[(yMat - muMat)_{0..J-2}]
        const Eigen::MatrixXd xMat = m_yMat.leftCols(Jm1) - m_muMat.leftCols(Jm1);
        Eigen::MatrixXd iPsi_xMat;
        applyIPsi(xMat, iPsi_xMat);
        m_YMat.resize(n, Jm1);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < Jm1; ++j)
                m_YMat(i, j) = m_eta(i) + m_iRMat(i, j) * iPsi_xMat(i, j);
    }

    void applyIPsi(
        const Eigen::MatrixXd &xMat,
        Eigen::MatrixXd &out
    ) const {
        const int n = m_n, Jm1 = m_Jm1;
        out.resize(n, Jm1);
        for (int i = 0; i < n; ++i) {
            const double inv_muJ = 1.0 / m_muMat(i, Jm1);
            const double base = xMat.row(i).sum() * inv_muJ;
            for (int j = 0; j < Jm1; ++j)
                out(i, j) = base + xMat(i, j) / m_muMat(i, j);
        }
    }

    void applySigma(
        const Eigen::MatrixXd &xMat,
        Eigen::MatrixXd &out
    ) const {
        // out = iRMat ⊙ iPsi(iRMat ⊙ x) + tau · broadcast(K · rowSum(x))
        const Eigen::MatrixXd iR_x = m_iRMat.array() * xMat.array();
        Eigen::MatrixXd iPsi_iR_x;
        applyIPsi(iR_x, iPsi_iR_x);
        out = m_iRMat.array() * iPsi_iR_x.array();
        if (m_tau != 0.0) {
            const Eigen::VectorXd tZ = xMat.rowwise().sum();
            Eigen::VectorXd V_tZ(m_n);
            m_grm.multiply(tZ.data(), V_tZ.data(), static_cast<uint32_t>(m_n));
            for (int j = 0; j < m_Jm1; ++j) out.col(j) += m_tau * V_tZ;
        }
    }

    void buildPreconditioner() {
        const int n = m_n, Jm1 = m_Jm1;
        m_InvBlock.assign(n, Eigen::MatrixXd());
        if (Jm1 == 1) {
            for (int i = 0; i < n; ++i) {
                const double iR = m_iRMat(i, 0);
                const double inv_muJ = 1.0 / m_muMat(i, 1);
                const double inv_mu0 = 1.0 / m_muMat(i, 0);
                const double M = iR * iR * inv_muJ + m_tau * m_diagK[i] + iR * iR * inv_mu0;
                m_InvBlock[i].resize(1, 1);
                m_InvBlock[i](0, 0) = 1.0 / (M > 1e-20 ? M : 1e-20);
            }
            return;
        }
        Eigen::MatrixXd M(Jm1, Jm1);
        for (int i = 0; i < n; ++i) {
            const double inv_muJ = 1.0 / m_muMat(i, Jm1);
            const double td = m_tau * m_diagK[i];
            for (int j2 = 0; j2 < Jm1; ++j2) {
                for (int j1 = 0; j1 < Jm1; ++j1) {
                    double v = m_iRMat(i, j2) * inv_muJ * m_iRMat(i, j1);
                    if (j1 == j2) {
                        v += m_iRMat(i, j2) * m_iRMat(i, j1) / m_muMat(i, j1);
                        v += td;
                    }
                    M(j2, j1) = v;
                }
            }
            m_InvBlock[i] = M.inverse();
        }
    }

    void applyPrecond(
        const Eigen::MatrixXd &xMat,
        Eigen::MatrixXd &out
    ) const {
        const int n = m_n, Jm1 = m_Jm1;
        out.resize(n, Jm1);
        if (Jm1 == 1) {
            for (int i = 0; i < n; ++i)
                out(i, 0) = xMat(i, 0) * m_InvBlock[i](0, 0);
            return;
        }
        for (int i = 0; i < n; ++i)
            out.row(i) = xMat.row(i) * m_InvBlock[i];
    }

    void pcgSolveMatInternal(
        const Eigen::MatrixXd &rhs,
        Eigen::MatrixXd &x
    ) const {
        x.setZero(m_n, m_Jm1);
        const double scale = static_cast<double>(m_n * m_Jm1);

        Eigen::MatrixXd r2;
        applySigma(x, r2);
        r2 = rhs - r2;
        double meanL2 = std::sqrt((r2.array() * r2.array()).sum() / scale);
        if (meanL2 <= kTolPCG) return;

        Eigen::MatrixXd z2;
        applyPrecond(r2, z2);
        Eigen::MatrixXd pMat = z2;
        Eigen::MatrixXd ApMat;
        applySigma(pMat, ApMat);

        auto ip = [](const Eigen::MatrixXd &A, const Eigen::MatrixXd &B) {
            return (A.array() * B.array()).sum();
        };

        double rz = ip(z2, r2);
        double pAp = ip(pMat, ApMat);
        if (std::abs(pAp) < 1e-30) return;
        double alpha = rz / pAp;
        x += alpha * pMat;
        Eigen::MatrixXd r1 = r2;
        Eigen::MatrixXd z1 = z2;
        r2 = r1 - alpha * ApMat;
        meanL2 = std::sqrt((r2.array() * r2.array()).sum() / scale);

        int iter = 1;
        while (meanL2 > kTolPCG && iter < kMaxIterPCG) {
            ++iter;
            applyPrecond(r2, z2);
            const double rz_new = ip(z2, r2);
            const double rz_old = ip(z1, r1);
            const double beta1 = rz_new / (std::abs(rz_old) > 1e-30 ? rz_old : 1e-30);
            pMat = z2 + beta1 * pMat;
            applySigma(pMat, ApMat);
            pAp = ip(pMat, ApMat);
            if (std::abs(pAp) < 1e-30) break;
            alpha = rz_new / pAp;
            x += alpha * pMat;
            r1 = r2;
            z1 = z2;
            r2 = r1 - alpha * ApMat;
            meanL2 = std::sqrt((r2.array() * r2.array()).sum() / scale);
        }
    }

    // PCG-solve every column of M individually.  Result has shape M.rows() × M.cols().
    void pcgSolveColumns(
        const Eigen::MatrixXd &M,
        Eigen::MatrixXd &out
    ) const {
        const int dim = m_n * m_Jm1;
        out.resize(dim, M.cols());
        Eigen::MatrixXd rhsMat, xMat;
        for (int c = 0; c < M.cols(); ++c) {
            Eigen::VectorXd v = M.col(c);
            vecToMat(v, m_n, m_Jm1, rhsMat);
            pcgSolveMatInternal(rhsMat, xMat);
            for (int i = 0; i < m_n; ++i)
                for (int j = 0; j < m_Jm1; ++j)
                    out(i * m_Jm1 + j, c) = xMat(i, j);
        }
    }

    void buildTraceRandMat() {
        const int dim = m_n * m_Jm1;
        m_TraceRandMat.resize(dim, kTraceNRun);
        std::mt19937 rng(12345);
        std::uniform_int_distribution<int> coin(0, 1);
        for (int t = 0; t < kTraceNRun; ++t)
            for (int k = 0; k < dim; ++k)
                m_TraceRandMat(k, t) = (coin(rng) == 0) ? -1.0 : 1.0;

        m_V_TRM.resize(dim, kTraceNRun);
        Eigen::VectorXd zCol(m_n), kzCol(m_n);
        for (int t = 0; t < kTraceNRun; ++t) {
            for (int i = 0; i < m_n; ++i) {
                double s = 0.0;
                for (int j = 0; j < m_Jm1; ++j) s += m_TraceRandMat(i * m_Jm1 + j, t);
                zCol(i) = s;
            }
            m_grm.multiply(zCol.data(), kzCol.data(), static_cast<uint32_t>(m_n));
            for (int i = 0; i < m_n; ++i) {
                const double v = kzCol(i);
                for (int j = 0; j < m_Jm1; ++j) m_V_TRM(i * m_Jm1 + j, t) = v;
            }
        }
    }

    // β, b update via the n*(J-1)-dimensional working LMM (reference updatePara).
    void updatePara() {
        buildPreconditioner();
        pcgSolveColumns(m_X_exp, m_iSigma_CovaMat);
        Eigen::VectorXd YVec = matToVec(m_YMat, m_n, m_Jm1);

        // single-vector PCG
        Eigen::MatrixXd yMat_rhs, yMat_sol;
        vecToMat(YVec, m_n, m_Jm1, yMat_rhs);
        pcgSolveMatInternal(yMat_rhs, yMat_sol);
        Eigen::VectorXd iSigma_Y = matToVec(yMat_sol, m_n, m_Jm1);

        Eigen::MatrixXd XtiSX = m_X_exp.transpose() * m_iSigma_CovaMat;
        XtiSX.diagonal().array() += 1e-12;
        const Eigen::MatrixXd XSX_inv = XtiSX.inverse();

        const Eigen::VectorXd Xt_iSigma_Y = m_X_exp.transpose() * iSigma_Y;
        m_beta = XSX_inv * Xt_iSigma_Y;

        m_iSigmaX_XSigmaX = m_iSigma_CovaMat * XSX_inv;

        const Eigen::VectorXd Z_iSY = ZsumVec(iSigma_Y, m_n, m_Jm1);
        const Eigen::VectorXd iSXbeta = m_iSigma_CovaMat * m_beta;
        const Eigen::VectorXd Z_iSXb = ZsumVec(iSXbeta, m_n, m_Jm1);
        const Eigen::VectorXd tmp = Z_iSY - Z_iSXb;
        Eigen::VectorXd Kt(m_n);
        m_grm.multiply(tmp.data(), Kt.data(), static_cast<uint32_t>(m_n));
        m_bVec = m_tau * Kt;
    }

    void updateEpsOneStep() {
        const int Jm1 = m_Jm1;
        const int Jm2 = Jm1 - 1;
        if (Jm2 < 1) return;

        Eigen::VectorXd d1 = Eigen::VectorXd::Zero(Jm2);
        Eigen::MatrixXd d2 = Eigen::MatrixXd::Zero(Jm2, Jm2);

        for (int k = 1; k < Jm1; ++k) {
            for (int i = 0; i < m_n; ++i) {
                const double mu_k = m_muMat(i, k);
                const double mu_kp1 = m_muMat(i, k + 1);
                if (mu_k < 1e-20 || mu_kp1 < 1e-20) continue;
                const double yk = m_yMat(i, k);
                const double ykp1 = m_yMat(i, k + 1);
                const double t1 = yk / mu_k - ykp1 / mu_kp1;
                const double t2 = -yk / (mu_k * mu_k) - ykp1 / (mu_kp1 * mu_kp1);
                const double W = m_WMat(i, k);
                d1(k - 1) += W * t1;
                d2(k - 1, k - 1) += W * (1.0 - 2.0 * m_nuMat(i, k)) * t1 + W * W * t2;
                if (k < Jm2) {
                    const double t3 = W * m_WMat(i, k + 1) * ykp1 / (mu_kp1 * mu_kp1);
                    d2(k - 1, k) += t3;
                    d2(k, k - 1) += t3;
                }
            }
        }
        // deps = -inv(d2) · d1
        const Eigen::VectorXd deps = -d2.fullPivLu().solve(d1);
        for (int k = 1; k < Jm1; ++k) m_eps(k) += deps(k - 1);
    }

    void updateEps() {
        for (int it = 0; it < kMaxIterEps; ++it) {
            const Eigen::VectorXd eps0 = m_eps;
            updateEpsOneStep();
            updateMats();
            double diff = 0.0;
            for (int j = 0; j < m_Jm1; ++j) {
                const double num = std::abs(m_eps(j) - eps0(j));
                const double den = std::abs(m_eps(j)) + std::abs(eps0(j)) + kTolEps;
                diff = std::max(diff, num / den);
            }
            if (diff < kTolEps) break;
        }
    }

    void updateParaConv() {
        for (int it = 0; it < kMaxIterPQL; ++it) {
            const Eigen::VectorXd beta0 = m_beta;
            updatePara();
            updateMats();
            updateEps();
            updateMats();
            double diff = 0.0;
            for (int k = 0; k < m_p; ++k) {
                const double num = std::abs(m_beta(k) - beta0(k));
                const double den = std::abs(m_beta(k)) + std::abs(beta0(k)) + kTolBeta;
                diff = std::max(diff, num / den);
            }
            if (diff < kTolBeta) break;
        }
    }

    void updateTau() {
        buildPreconditioner();
        pcgSolveColumns(m_X_exp, m_iSigma_CovaMat);

        Eigen::VectorXd YVec = matToVec(m_YMat, m_n, m_Jm1);
        Eigen::MatrixXd yMat_rhs, yMat_sol;
        vecToMat(YVec, m_n, m_Jm1, yMat_rhs);
        pcgSolveMatInternal(yMat_rhs, yMat_sol);
        Eigen::VectorXd iSigma_Y = matToVec(yMat_sol, m_n, m_Jm1);

        Eigen::MatrixXd XtiSX = m_X_exp.transpose() * m_iSigma_CovaMat;
        XtiSX.diagonal().array() += 1e-12;
        const Eigen::MatrixXd XSX_inv = XtiSX.inverse();
        m_iSigmaX_XSigmaX = m_iSigma_CovaMat * XSX_inv;

        const Eigen::VectorXd Xt_iSigma_Y = m_X_exp.transpose() * iSigma_Y;
        const Eigen::VectorXd PY = iSigma_Y - m_iSigmaX_XSigmaX * Xt_iSigma_Y;
        const Eigen::VectorXd ZPY = ZsumVec(PY, m_n, m_Jm1);
        Eigen::VectorXd K_ZPY(m_n);
        m_grm.multiply(ZPY.data(), K_ZPY.data(), static_cast<uint32_t>(m_n));
        Eigen::VectorXd VPY(m_n * m_Jm1);
        for (int i = 0; i < m_n; ++i)
            for (int j = 0; j < m_Jm1; ++j) VPY(i * m_Jm1 + j) = K_ZPY(i);

        Eigen::MatrixXd vpy_rhs, vpy_sol;
        vecToMat(VPY, m_n, m_Jm1, vpy_rhs);
        pcgSolveMatInternal(vpy_rhs, vpy_sol);
        Eigen::VectorXd iSigma_VPY = matToVec(vpy_sol, m_n, m_Jm1);

        const Eigen::VectorXd Xt_iSigma_VPY = m_X_exp.transpose() * iSigma_VPY;
        const Eigen::VectorXd PVPY = iSigma_VPY - m_iSigmaX_XSigmaX * Xt_iSigma_VPY;

        const double YPVPY = YVec.dot(PVPY);
        const double YPVPVPY = VPY.dot(PVPY);

        // tr(P V): iSigma · V_TRM, project, dot with TraceRandMat.
        Eigen::MatrixXd iSigma_VTRM;
        pcgSolveColumns(m_V_TRM, iSigma_VTRM);
        const Eigen::MatrixXd XtiSV = m_X_exp.transpose() * iSigma_VTRM;
        const Eigen::MatrixXd P_VTRM = iSigma_VTRM - m_iSigmaX_XSigmaX * XtiSV;

        double tracePV = 0.0;
        for (int t = 0; t < kTraceNRun; ++t)
            tracePV += m_TraceRandMat.col(t).dot(P_VTRM.col(t));
        tracePV /= static_cast<double>(kTraceNRun);

        const double score = 0.5 * YPVPY - 0.5 * tracePV;
        const double AI = 0.5 * YPVPVPY;
        if (!std::isfinite(score) || !std::isfinite(AI) || std::abs(AI) < 1e-12) {
            infoMsg(
                "[%s] AI-REML step skipped: score=%.3e, AI=%.3e",
                m_tag,
                score,
                AI
            );
            return;
        }
        double dtau = score / AI;
        const double tau0 = m_tau;
        double tau_new = tau0 + dtau;
        int halve_count = 0;
        while (tau_new < 0.0 && halve_count < 40) {
            dtau *= 0.5;
            tau_new = tau0 + dtau;
            ++halve_count;
        }
        if (tau_new < 1e-4) tau_new = 0.0;
        m_tau = tau_new;
    }
};

// ══════════════════════════════════════════════════════════════════════
// Driver: build a POLMMNullModel via POLMMFitter.
// ══════════════════════════════════════════════════════════════════════

// Forward declarations for the combined driver.
void computeTestMatrices(
    const Eigen::VectorXi &yVec,
    const Eigen::MatrixXd &X,
    POLMMNullModel &null,
    const char *tag
);
void fillRymuVec(const Eigen::VectorXi &yVec, POLMMNullModel &null);
void estimateVarianceRatiosFromUnion(
    POLMMFitter &fitter,
    POLMMNullModel &null,
    const Eigen::MatrixXd &X,
    const Eigen::VectorXi &yVec,
    const SparseGRM &grm,
    GenoMeta &genoData,
    const std::vector<uint32_t> &localToUnion,
    const char *tag
);

// Combined driver: fits null model and estimates variance ratios sharing
// one POLMMFitter instance (so PCG state from the converged fit is reused).
void fitNullModelAndVarianceRatios(
    const Eigen::VectorXi &yVec,
    const Eigen::MatrixXd &X,
    const SparseGRM &grm,
    int n,
    int J,
    int p,
    POLMMNullModel &null,
    GenoMeta &genoData,
    const std::vector<uint32_t> &localToUnion,
    const char *tag = ""
) {
    POLMMFitter fitter(yVec, X, grm, n, J, p, tag);
    fitter.fit(null);
    computeTestMatrices(yVec, X, null, tag);
    fillRymuVec(yVec, null);
    estimateVarianceRatiosFromUnion(fitter, null, X, yVec, grm, genoData, localToUnion, tag);
}

// ══════════════════════════════════════════════════════════════════════
// Test matrices for marker-level analysis (mirrors reference getobjP).
//
//   XR_Psi_R (p × n(J-1)):
//     for each covariate column X_k expanded to n*(J-1), compute
//       Psi · (X_exp_k / iRMat)  then divide by iRMat element-wise.
//   HessP = XR_Psi_R · X_exp (p × p), and its inverse.
//   XXR_Psi_RX_new (n × p) = X · HessP^{-1}
//   XR_Psi_R_new  (p × n)  = per-subject row-sum of XR_Psi_R.
//   RymuVec (n)     = per-subject row-sum of (y - mu)_{0..J-2} / iRMat.
//   RPsiR  (n)      = quadratic form  sum_{j1,j2} (1/iR)·Psi·(1/iR)
//                  = sum_j (mu_j/iR_j²) - (sum_j mu_j/iR_j)²
// ══════════════════════════════════════════════════════════════════════

void computeTestMatrices(
    const Eigen::VectorXi & /*yVec*/,
    const Eigen::MatrixXd &X,
    POLMMNullModel &null,
    const char *tag = ""
) {
    const int n = null.n, J = null.J, p = null.p;
    const int Jm1 = J - 1;
    const Eigen::MatrixXd &muMat = null.muMat;
    const Eigen::MatrixXd &iRMat = null.iRMat;

    // XR_Psi_R: p × n(J-1)
    Eigen::MatrixXd XR_Psi_R(p, n * Jm1);
    for (int k = 0; k < p; ++k) {
        // x_over_iR(i,j) = X(i,k) / iRMat(i,j)
        Eigen::MatrixXd x_over_iR(n, Jm1);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < Jm1; ++j)
                x_over_iR(i, j) = X(i, k) / iRMat(i, j);

        // Psi_x = mu ⊙ x - mu * (mu' x) per subject (multinomial-Psi)
        Eigen::MatrixXd Psi_x(n, Jm1);
        for (int i = 0; i < n; ++i) {
            double smux = 0.0;
            for (int j = 0; j < Jm1; ++j) smux += muMat(i, j) * x_over_iR(i, j);
            for (int j = 0; j < Jm1; ++j)
                Psi_x(i, j) = muMat(i, j) * x_over_iR(i, j) - muMat(i, j) * smux;
        }

        // out = Psi_x / iRMat, stored as XR_Psi_R(k, i*(J-1)+j)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < Jm1; ++j)
                XR_Psi_R(k, i * Jm1 + j) = Psi_x(i, j) / iRMat(i, j);
    }

    // HessP = XR_Psi_R · X_exp  (p × p)
    Eigen::MatrixXd X_exp = buildExpandedX(X, n, p, Jm1);
    Eigen::MatrixXd HessP = XR_Psi_R * X_exp;
    HessP.diagonal().array() += 1e-12;
    const Eigen::MatrixXd HessP_inv = HessP.inverse();

    // XXR_Psi_RX_new = X · HessP^{-1}  (n × p)
    null.XXR_Psi_RX_new = X * HessP_inv;

    // XR_Psi_R_new = per-subject row-sum  (p × n)
    Eigen::MatrixXd XR_new = Eigen::MatrixXd::Zero(p, n);
    for (int k = 0; k < p; ++k)
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < Jm1; ++j)
                XR_new(k, i) += XR_Psi_R(k, i * Jm1 + j);
    null.XR_Psi_R_new = XR_new;

    // RymuVec(i) = sum_j (y(i,j) - mu(i,j)) / iRMat(i,j)
    Eigen::MatrixXd yIndic = Eigen::MatrixXd::Zero(n, J);
    for (int i = 0; i < n; ++i) {
        // Reconstruct yIndic from fitted muMat via known yVec inside fitter.
        // We need yVec to rebuild yIndic — pass it in via the caller.
    }
    // The above placeholder requires yVec; we recompute it here from yMat is not
    // available, so we expect caller to supply it.  Instead, store RymuVec / RPsiR
    // directly using formulas that depend only on mu and iR:
    //   RPsiR(i) = sum_j mu_j/iR_j² - (sum_j mu_j/iR_j)²
    Eigen::VectorXd RPsiR(n);
    for (int i = 0; i < n; ++i) {
        double s_muR = 0.0;
        double s_diag = 0.0;
        for (int j = 0; j < Jm1; ++j) {
            const double muR = muMat(i, j) / iRMat(i, j);
            s_muR += muR;
            s_diag += muR / iRMat(i, j);
        }
        RPsiR(i) = s_diag - s_muR * s_muR;
    }
    null.RPsiR = RPsiR;
    null.RymuVec.resize(n);     // populated by caller via reconstructFromYVec below.

    (void)tag;
}

// Helper: populate RymuVec from yVec after computeTestMatrices.
void fillRymuVec(
    const Eigen::VectorXi &yVec,
    POLMMNullModel &null
) {
    const int n = null.n, J = null.J;
    const int Jm1 = J - 1;
    const Eigen::MatrixXd &muMat = null.muMat;
    const Eigen::MatrixXd &iRMat = null.iRMat;
    null.RymuVec.resize(n);
    for (int i = 0; i < n; ++i) {
        double s = 0.0;
        const int yi = yVec(i);
        for (int j = 0; j < Jm1; ++j) {
            const double y_ij = (yi == j) ? 1.0 : 0.0;
            s += (y_ij - muMat(i, j)) / iRMat(i, j);
        }
        null.RymuVec(i) = s;
    }
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
    POLMMFitter &fitter,
    POLMMNullModel &null,
    const Eigen::MatrixXd & /*X*/,
    const Eigen::VectorXi & /*yVec*/,
    const SparseGRM & /*grm*/,
    GenoMeta &genoData,
    const std::vector<uint32_t> &localToUnion,
    const char *tag = ""
) {
    const int n = null.n;
    const int J = null.J;
    const int Jm1 = J - 1;
    const int nBins = static_cast<int>(kMacBins.size());

    // Ensure PCG cache is consistent with the converged tau.
    fitter.prepareForVarianceRatio();

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

        // PCG-based VarP, matching reference getVarP():
        //   adjGLong = tZ(adjG)         duplicate adjG across (J-1) slots
        //   iSigmaG  = Σ^{-1} adjGLong
        //   PadjG    = iSigmaG - iSigmaX_XSigmaX · (X_exp' iSigmaG)
        //   VarP     = adjGLong' · PadjG
        Eigen::VectorXd adjGLong(n * Jm1);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < Jm1; ++j)
                adjGLong(i * Jm1 + j) = adjG(i);
        Eigen::VectorXd iSigmaG = fitter.pcgSolveVec(adjGLong);
        Eigen::VectorXd Xt_iSigmaG = fitter.X_exp().transpose() * iSigmaG;
        Eigen::VectorXd PadjG = iSigmaG - fitter.iSigmaX_XSigmaX() * Xt_iSigmaG;
        double VarP = adjGLong.dot(PadjG);
        if (!std::isfinite(VarP) || VarP < 1e-10) continue;
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

    // Reference convention (Main.cpp::mainMarkerInCPP): when the alt allele
    // had freq > 0.5, GVec has been flipped to test the minor allele.  The
    // output Beta is then negated so it refers to the original alt allele,
    // but the Z-score is reported as-is (referring to the tested/minor
    // allele direction).  See examples_1kg/vs_grab024/GRAB/src/Main.cpp:591.
    if (flipped) BETA = -BETA;

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

    // Newton-Raphson on K'(zeta) = target, with sign-flip damping (mirrors
    // reference fastgetroot_K1 in examples_1kg/vs_grab024/GRAB/src/POLMM.cpp).
    auto solveTail = [&](double target, double init) -> std::pair<double, bool> {
        constexpr double tol = 1e-4;
        constexpr int maxiter = 100;
        double x = init;
        double diffX = std::numeric_limits<double>::infinity();
        double oldK1 = 0.0;
        bool converged = true;
        for (int iter = 0; iter < maxiter; ++iter) {
            double oldX = x;
            double oldDiffX = diffX;
            double k1 = K1(x, cMat.data(), muSub.data(), nSub, Jm1, Ratio0, m1) - target;
            double k2 = K2(x, cMat.data(), muSub.data(), nSub, Jm1, Ratio0);
            if (!std::isfinite(k1)) {
                x = std::copysign(std::numeric_limits<double>::infinity(), target);
                break;
            }
            diffX = -k1 / k2;
            if (iter > 0 && ((k1 > 0.0) != (oldK1 > 0.0))) {
                while (std::abs(diffX) > std::abs(oldDiffX) - tol)
                    diffX *= 0.5;
            }
            oldK1 = k1;
            if (std::abs(diffX) < tol) {
                converged = true;
                break;
            }
            x = oldX + diffX;
            if (iter == maxiter - 1) converged = false;
        }
        return {x, converged};
    };

    // Barndorff-Nielsen saddlepoint tail probability (matches reference
    // fastGet_Saddle_Prob).
    //   lowerTail=false → P(S > target),  lowerTail=true → P(S < target).
    auto tailProb = [&](double zeta, double target, bool lowerTail) -> double {
        double k0 = K0(zeta, cMat.data(), muSub.data(), nSub, Jm1, Ratio0, m1);
        double k2 = K2(zeta, cMat.data(), muSub.data(), nSub, Jm1, Ratio0);
        if (!std::isfinite(k0) || !std::isfinite(k2)) return 0.0;
        double w2 = 2.0 * (zeta * target - k0);
        if (w2 < 0.0) return 0.0;
        double w = std::copysign(std::sqrt(w2), zeta);
        double v = zeta * std::sqrt(k2);
        if (std::abs(w) < 1e-14 || std::abs(v) < 1e-14) return 0.0;
        double Z_BN = w + (1.0 / w) * std::log(v / w);
        double sign = lowerTail ? 1.0 : -1.0;
        return 0.5 * std::erfc(-sign * Z_BN / M_SQRT2);
    };

    const double absStat = std::abs(qNorm);
    auto upper = solveTail(absStat, 3.0);
    auto lower = solveTail(-absStat, -3.0);

    if (!upper.second || !lower.second)
        return std::numeric_limits<double>::quiet_NaN();

    double p1 = tailProb(upper.first, absStat, false);
    double p2 = tailProb(lower.first, -absStat, true);
    double pval = p1 + p2;
    if (!std::isfinite(pval) || pval <= 0.0) return std::numeric_limits<double>::quiet_NaN();
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

    auto solveTail = [&](double target, double init) -> std::pair<double, bool> {
        constexpr double tol = 1e-4;
        constexpr int maxiter = 100;
        double x = init;
        double diffX = std::numeric_limits<double>::infinity();
        double oldK1 = 0.0;
        bool converged = true;
        for (int iter = 0; iter < maxiter; ++iter) {
            double oldX = x;
            double oldDiffX = diffX;
            double k1 = K1bin(x) - target;
            double k2 = K2bin(x);
            if (!std::isfinite(k1)) {
                x = std::copysign(std::numeric_limits<double>::infinity(), target);
                break;
            }
            diffX = -k1 / k2;
            if (iter > 0 && ((k1 > 0.0) != (oldK1 > 0.0))) {
                while (std::abs(diffX) > std::abs(oldDiffX) - tol)
                    diffX *= 0.5;
            }
            oldK1 = k1;
            if (std::abs(diffX) < tol) {
                converged = true;
                break;
            }
            x = oldX + diffX;
            if (iter == maxiter - 1) converged = false;
        }
        return {x, converged};
    };

    auto tailProb = [&](double zeta, double target, bool lowerTail) -> double {
        double k0 = K0bin(zeta);
        double k2 = K2bin(zeta);
        if (!std::isfinite(k0) || !std::isfinite(k2)) return 0.0;
        double w2 = 2.0 * (zeta * target - k0);
        if (w2 < 0.0) return 0.0;
        double w = std::copysign(std::sqrt(w2), zeta);
        double v = zeta * std::sqrt(k2);
        if (std::abs(w) < 1e-14 || std::abs(v) < 1e-14) return 0.0;
        double Z_BN = w + (1.0 / w) * std::log(v / w);
        double sign = lowerTail ? 1.0 : -1.0;
        return 0.5 * std::erfc(-sign * Z_BN / M_SQRT2);
    };

    const double absStat = std::abs(qNorm);
    auto upper = solveTail(absStat, 3.0);
    auto lower = solveTail(-absStat, -3.0);

    if (!upper.second || !lower.second)
        return std::numeric_limits<double>::quiet_NaN();

    double p_upper = tailProb(upper.first, absStat, false);
    double p_lower = tailProb(lower.first, -absStat, true);
    double pval = p_upper + p_lower;
    if (!std::isfinite(pval) || pval <= 0.0) return std::numeric_limits<double>::quiet_NaN();
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

                // Fit null model + variance ratios (one fitter instance, shared PCG state)
                POLMMNullModel null;
                fitNullModelAndVarianceRatios(
                    yVec, X, phenoGrm,
                    static_cast<int>(n_k), J, p,
                    null,
                    *genoData, localToUnion,
                    name.c_str()
                );

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
