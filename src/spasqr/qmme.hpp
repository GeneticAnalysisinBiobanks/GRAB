// qmme.hpp — Smoothed quantile regression via Quadratic Majorization
// Minimization with Extrapolation (Heng & Wang, 2025).
//
//   H = (1 / (n √(2π) h)) Z^T Z + δI
// where Z = [1 | standardized X]. Phenotype-level construction caches
// Z^T Z / n once. Per-bandwidth Cholesky is rebuilt by prepareBandwidth().
#pragma once

#include <Eigen/Dense>

namespace qmme {

// Per-fit convergence diagnostics for SqrSolver::solve.
struct SolverStatus {
    int    iter          = 0;     // QMME iterations executed
    bool   converged     = false; // true iff final ||grad||_∞ ≤ tol
    double finalGradNorm = 0.0;   // final ||grad||_∞
};

class SqrSolver {
  public:
    // Construct from raw X (n × p, no intercept). Computes and caches
    // standardization (column mean/sd) and Z^T Z / n where Z = [1|X_std].
    SqrSolver(
        const Eigen::MatrixXd &X,
        double delta = 1e-6
    );

    // Update Cholesky cache for new bandwidth h. O((p+1)^3); cheap for
    // small p. No-op if h already matches the cached one.
    void prepareBandwidth(double h);

    // Solve smoothed QR with the cached Cholesky.  Returns (p+1)
    // coefficients in the ORIGINAL (un-standardized) space. residOut
    // (if non-null) gets the final residual vector y - X·beta_orig.
    // statusOut (if non-null) is populated with iteration count,
    // convergence flag, and final gradient norm.
    //
    // initBetaOrig (optional): warm-start coefficients in ORIGINAL space
    // (size = p+1; layout = [intercept, β_X_1, ..., β_X_p]).  When non-null,
    // converted to standardized space and used as the QMME iterate start
    // instead of the cold-start (β = 0 with intercept = empirical τ-quantile
    // of centered Y).  Used by wald-mode τ-chaining to recycle β̂(τ_{k-1})
    // as init for τ_k.
    Eigen::VectorXd solve(
        const Eigen::VectorXd &Y,
        double tau,
        Eigen::VectorXd *residOut = nullptr,
        double tol = 1e-7,
        int maxIter = 5000,
        int restartPeriod = 50,
        SolverStatus *statusOut = nullptr,
        const Eigen::VectorXd *initBetaOrig = nullptr
    );

    // Diagnostic: number of QMME iterations executed in the most recent
    // solve(); -1 before any solve.
    int lastIters() const {
        return m_lastIters;
    }

  private:
    int m_n;                                // sample size
    int m_p;                                // raw covariate count (no intercept)
    double m_delta;
    Eigen::MatrixXd m_Z;                    // n × (p+1)
    Eigen::RowVectorXd m_mx;                // 1 × p
    Eigen::VectorXd m_sx;                   // p (1/sd per column)
    Eigen::MatrixXd m_ZtZ_over_n;           // (p+1) × (p+1)

    double m_currentH = -1.0;
    Eigen::LLT<Eigen::MatrixXd> m_chol;

    int m_lastIters = -1;
};

} // namespace qmme
