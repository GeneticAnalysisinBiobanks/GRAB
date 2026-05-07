// test_qmme_vs_conquer.cpp — compare QMME vs conquer on synthetic
// smoothed quantile regression, including large-scale and ill-conditioned
// scenarios. Reports β diff, gradient norms, iter counts, AND wall-clock
// timing.
//
// Compile:
//   g++ -std=c++17 -O3 -DNDEBUG -Isrc -Ithird_party/eigen-5.0.0 \
//       -Ithird_party/boost-1.90.0 \
//       test_qmme_vs_conquer.cpp \
//       src/spasqr/qmme.cpp src/spasqr/conquer.cpp \
//       src/util/math_helper.cpp \
//       -o test_qmme_vs_conquer

#include "spasqr/conquer.hpp"
#include "spasqr/qmme.hpp"
#include "util/math_helper.hpp"

#include <Eigen/Dense>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

static double gradInfNormOrig(
    const Eigen::MatrixXd &X,
    const Eigen::VectorXd &Y,
    double tau,
    double h,
    const Eigen::VectorXd &beta
) {
    const int n = static_cast<int>(X.rows());
    const int p = static_cast<int>(X.cols());
    Eigen::MatrixXd Xa(n, p + 1);
    Xa.col(0).setOnes();
    Xa.rightCols(p) = X;
    const Eigen::VectorXd r = Y - Xa * beta;
    Eigen::VectorXd der(n);
    for (int i = 0; i < n; ++i)
        der(i) = math::pnorm(-r(i) / h) - tau;
    Eigen::VectorXd g = (1.0 / n) * Xa.transpose() * der;
    return g.lpNorm<Eigen::Infinity>();
}

// AR(1)-style correlated design: column j+1 = ρ · column j + sqrt(1-ρ²)·noise.
// Larger ρ → higher condition number of X^T X.
static Eigen::MatrixXd makeAR1Design(
    int n,
    int p,
    double rho,
    std::mt19937_64 &rng
) {
    std::normal_distribution<double> norm(0.0, 1.0);
    Eigen::MatrixXd X(n, p);
    for (int i = 0; i < n; ++i) X(i, 0) = norm(rng);
    const double sqrt_one_minus_rho2 = std::sqrt(1.0 - rho * rho);
    for (int j = 1; j < p; ++j)
        for (int i = 0; i < n; ++i)
            X(i, j) = rho * X(i, j - 1) + sqrt_one_minus_rho2 * norm(rng);
    return X;
}

struct ScenarioResult {
    double timeC_ms = 0.0;
    double timeQ_ms = 0.0;
    int iterC_total = 0; // huber + gauss
    int iterQ = 0;
    double gnormC = 0.0;
    double gnormQ = 0.0;
    bool convC = false;
    bool convQ = false;
    double betaRelDiff = 0.0;
    double cond = 0.0; // approximate condition number of X^T X / n
};

static ScenarioResult runScenario(
    const char *name,
    int n,
    int p,
    double tau,
    double hScale,
    double rho,
    uint64_t seed
) {
    std::mt19937_64 rng(seed);
    Eigen::MatrixXd X = makeAR1Design(n, p, rho, rng);
    std::student_t_distribution<double> tdist(3.0);

    // True signal — use a random direction so collinearity actually matters
    Eigen::VectorXd beta_true(p);
    {
        std::normal_distribution<double> norm(0.0, 1.0);
        for (int j = 0; j < p; ++j) beta_true(j) = norm(rng);
        beta_true *= 1.0 / beta_true.norm();
    }

    Eigen::VectorXd Y(n);
    for (int i = 0; i < n; ++i)
        Y(i) = 1.0 + X.row(i).dot(beta_true) + tdist(rng);

    // Bandwidth: hScale × IQR(Y) / 3
    std::vector<double> Ysort(n);
    Eigen::VectorXd::Map(Ysort.data(), n) = Y;
    std::sort(Ysort.begin(), Ysort.end());
    auto q = [&](double prob) {
        const double idx = prob * (n - 1);
        const int lo = static_cast<int>(std::floor(idx));
        const int hi = std::min(lo + 1, n - 1);
        const double frac = idx - lo;
        return Ysort[lo] * (1.0 - frac) + Ysort[hi] * frac;
    };
    const double iqr = q(0.75) - q(0.25);
    const double h = hScale * iqr / 3.0;

    // Condition number of standardized X^T X / n (approx, p ≤ 64)
    double condNum = 0.0;
    {
        Eigen::RowVectorXd mx = X.colwise().mean();
        Eigen::MatrixXd Xs(n, p);
        for (int j = 0; j < p; ++j) {
            const double s = (X.col(j).array() - mx(j)).matrix().norm() /
                             std::sqrt(static_cast<double>(n - 1));
            const double inv = (s > 0.0) ? 1.0 / s : 1.0;
            Xs.col(j) = (X.col(j).array() - mx(j)).matrix() * inv;
        }
        Eigen::MatrixXd XtX = Xs.transpose() * Xs / n;
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(XtX);
        const auto &ev = es.eigenvalues();
        condNum = ev(p - 1) / std::max(ev(0), 1e-30);
    }

    const double tol = 1e-9;
    const int maxIter = 5000;
    ScenarioResult R;
    R.cond = condNum;

    // ── conquer (single fit timing) ──
    auto t0 = std::chrono::steady_clock::now();
    conquer::ConquerStatus statC;
    Eigen::VectorXd residC;
    Eigen::VectorXd betaC = conquer::smqrGauss(X, Y, tau, h, &residC, tol,
                                               maxIter, 100.0, &statC);
    auto t1 = std::chrono::steady_clock::now();
    R.timeC_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    R.iterC_total = statC.huberIter + statC.gaussIter;
    R.gnormC = gradInfNormOrig(X, Y, tau, h, betaC);
    R.convC = statC.converged;

    // ── QMME (construction + bandwidth + solve) ──
    t0 = std::chrono::steady_clock::now();
    qmme::SqrSolver solver(X, /*delta*/ 1e-6);
    solver.prepareBandwidth(h);
    conquer::ConquerStatus statQ;
    Eigen::VectorXd residQ;
    Eigen::VectorXd betaQ = solver.solve(Y, tau, &residQ, tol, maxIter,
                                         /*restartPeriod*/ 50, &statQ);
    t1 = std::chrono::steady_clock::now();
    R.timeQ_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    R.iterQ = statQ.gaussIter;
    R.gnormQ = gradInfNormOrig(X, Y, tau, h, betaQ);
    R.convQ = statQ.converged;

    R.betaRelDiff = (betaC - betaQ).norm() / std::max(betaC.norm(), 1e-12);

    std::printf(
        "  %-32s  cond=%.2e   "
        "C: %5dit %8.1fms ‖g‖=%.2e %s   "
        "Q: %4dit %8.1fms ‖g‖=%.2e %s   "
        "rel‖β‖=%.2e\n",
        name, condNum,
        R.iterC_total, R.timeC_ms, R.gnormC, R.convC ? "✓" : "✗",
        R.iterQ, R.timeQ_ms, R.gnormQ, R.convQ ? "✓" : "✗",
        R.betaRelDiff);
    return R;
}

// Multi-tau scenario: showcases QMME's amortized advantage when 9 taus
// share one Cholesky and one X^T X / n. conquer must run from scratch
// for each tau (no cross-tau cache).
static void runMultiTau(
    const char *name,
    int n,
    int p,
    double rho,
    double hScale,
    uint64_t seed
) {
    std::mt19937_64 rng(seed);
    Eigen::MatrixXd X = makeAR1Design(n, p, rho, rng);
    std::student_t_distribution<double> tdist(3.0);
    Eigen::VectorXd beta_true(p);
    {
        std::normal_distribution<double> norm(0.0, 1.0);
        for (int j = 0; j < p; ++j) beta_true(j) = norm(rng);
        beta_true *= 1.0 / beta_true.norm();
    }
    Eigen::VectorXd Y(n);
    for (int i = 0; i < n; ++i)
        Y(i) = 1.0 + X.row(i).dot(beta_true) + tdist(rng);

    std::vector<double> Ysort(n);
    Eigen::VectorXd::Map(Ysort.data(), n) = Y;
    std::sort(Ysort.begin(), Ysort.end());
    auto q = [&](double prob) {
        const double idx = prob * (n - 1);
        const int lo = static_cast<int>(std::floor(idx));
        return Ysort[lo] * (1.0 - (idx - lo)) +
               Ysort[std::min(lo + 1, n - 1)] * (idx - lo);
    };
    const double h = hScale * (q(0.75) - q(0.25)) / 3.0;
    const std::vector<double> taus = {0.1, 0.2, 0.3, 0.4, 0.5,
                                      0.6, 0.7, 0.8, 0.9};
    const double tol = 1e-9;

    // conquer: 9 independent fits, each rebuilds everything internally
    auto t0 = std::chrono::steady_clock::now();
    int totalIterC = 0;
    bool allConvC = true;
    for (double tau : taus) {
        conquer::ConquerStatus s;
        Eigen::VectorXd resid;
        conquer::smqrGauss(X, Y, tau, h, &resid, tol, 5000, 100.0, &s);
        totalIterC += s.huberIter + s.gaussIter;
        allConvC = allConvC && s.converged;
    }
    auto t1 = std::chrono::steady_clock::now();
    const double tC = std::chrono::duration<double, std::milli>(t1 - t0).count();

    // QMME: one solver instance, one prepareBandwidth(h), 9 solves
    t0 = std::chrono::steady_clock::now();
    qmme::SqrSolver solver(X, 1e-6);
    solver.prepareBandwidth(h);
    int totalIterQ = 0;
    bool allConvQ = true;
    for (double tau : taus) {
        conquer::ConquerStatus s;
        Eigen::VectorXd resid;
        solver.solve(Y, tau, &resid, tol, 5000, 50, &s);
        totalIterQ += s.gaussIter;
        allConvQ = allConvQ && s.converged;
    }
    t1 = std::chrono::steady_clock::now();
    const double tQ = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::printf("  %-32s  C: %5dit %8.1fms %s   Q: %4dit %8.1fms %s   speedup×%.2f\n",
                name, totalIterC, tC, allConvC ? "✓" : "✗",
                totalIterQ, tQ, allConvQ ? "✓" : "✗", tC / tQ);
}

int main() {
    std::printf("\nQMME vs conquer — large-scale & ill-conditioned\n");
    std::printf("(each row: median τ=0.5 unless noted; t-distributed noise df=3)\n\n");

    std::printf("Single-fit (τ=0.5):\n");
    runScenario("n=10k p=10 ρ=0.0",       10000, 10, 0.5, 1.0, 0.0, 1);
    runScenario("n=10k p=10 ρ=0.9",       10000, 10, 0.5, 1.0, 0.9, 1);
    runScenario("n=10k p=10 ρ=0.99",      10000, 10, 0.5, 1.0, 0.99, 1);
    runScenario("n=10k p=20 ρ=0.99",      10000, 20, 0.5, 1.0, 0.99, 1);
    runScenario("n=50k p=20 ρ=0.99",      50000, 20, 0.5, 1.0, 0.99, 2);
    runScenario("n=50k p=20 ρ=0.99 tight",50000, 20, 0.5, 0.3, 0.99, 2);
    runScenario("n=50k p=50 ρ=0.95",      50000, 50, 0.5, 1.0, 0.95, 3);
    runScenario("n=50k p=50 ρ=0.99",      50000, 50, 0.5, 1.0, 0.99, 3);
    runScenario("n=100k p=20 ρ=0.99",    100000, 20, 0.5, 1.0, 0.99, 4);
    runScenario("n=100k p=50 ρ=0.99",    100000, 50, 0.5, 1.0, 0.99, 4);

    std::printf("\nTail quantiles (n=50k p=20 ρ=0.99):\n");
    runScenario("τ=0.1",                   50000, 20, 0.1, 1.0, 0.99, 2);
    runScenario("τ=0.05",                  50000, 20, 0.05, 1.0, 0.99, 2);
    runScenario("τ=0.95",                  50000, 20, 0.95, 1.0, 0.99, 2);

    std::printf("\nMulti-τ (9 taus shared X, h):\n");
    runMultiTau("n=50k p=20 ρ=0.99",       50000, 20, 0.99, 1.0, 5);
    runMultiTau("n=50k p=50 ρ=0.99",       50000, 50, 0.99, 1.0, 6);
    runMultiTau("n=100k p=50 ρ=0.99",     100000, 50, 0.99, 1.0, 7);
    runMultiTau("n=100k p=20 ρ=0.99 tight",100000, 20, 0.99, 0.3, 8);

    return 0;
}
