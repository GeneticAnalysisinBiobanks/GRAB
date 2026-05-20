// sageld_fit.hpp — Long-format LMM fitting for SAGELD's direct-phenotype mode.
//
// Provides:
//   parseLongPheno()       — read a long-format phenotype/covariate file,
//                            group rows by IID, restrict to subjects also
//                            present in the .fam IID list.
//   fitRandomSlopeML()     — fit  Y ~ X + (E | IID)  via REML profile
//                            likelihood + 3-D Nelder-Mead simplex search;
//                            return per-row BLUP residuals and variance
//                            components.
//   fitRandomInterceptML() — fit  Y ~ 1 + (1 | IID)  via EM-ML; kept for
//                            source compatibility but no longer called by
//                            the SAGELD pheno-mode pipeline (R_E is not
//                            consumed by the closed-form lambda).
//   aggregatePerIID()      — Σ_j r_ij      (weight ≡ 1)
//   aggregateWeightedPerIID() — Σ_j w_ij r_ij
//
// fitRandomSlopeML profiles σ² out of the REML log-likelihood and uses
// math::nelderMead in 3 parameters of a lower-triangular Cholesky factor
// of τ = D/σ², so τ is always SPD without explicit constraint handling.
// Per evaluation cost ≈ one EM iteration; typical convergence in 50–200
// evaluations, ~25–100× faster than the previous EM-ML loop for our
// N = O(10⁴), n_subj = O(10³), q = 2 regime.  Matches lme4's REML default
// up to round-off; differs from ML by O(p/N) on σ̂² only.  See Bates et al.
// (2015), J Stat Softw 67(1), §3 for the parameterisation.

#pragma once

#include <Eigen/Dense>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

namespace nsSAGELDFit {

// Long-format phenotype representation in .fam order (after IID restriction).
//
// Layout invariants:
//   uniqueIIDs.size()              == nSubj            (subjects with ≥ 1 row)
//   subjStart.size()               == nSubj + 1        (CSR-style index)
//   X.rows() == E.size() == Y.rows() == N_total       (total rows = Σ_i n_i)
//   Rows for subject i occupy [subjStart[i], subjStart[i+1]).
//   uniqueIIDs is ordered by .fam IID order (subjects absent from the file
//   are excluded; subjects absent from .fam are dropped from the file).
struct LongPhenoData {
    std::vector<std::string> uniqueIIDs;
    std::vector<uint32_t> subjStart;       // size nSubj + 1
    std::vector<uint32_t> famSubjIdx;      // size nSubj; per-subject .fam index
    Eigen::MatrixXd X;                     // N_total × p (includes intercept column 0)
    Eigen::MatrixXd Y;                     // N_total × q (one column per phenotype)
    Eigen::VectorXd E;                     // N_total (the env variable for random slope)
    std::vector<std::string> phenoNames;   // q phenotype names
    std::vector<std::string> covarNames;   // p − 1 covariate names (X column 0 is intercept)
    std::string envName;                   // env column header
};

// Read long-format file, group rows by IID, restrict subjects to famIIDs.
// keptSubjects (when non-empty) further restricts to IIDs in the set —
// pass the intersection of GRM ID list, --keep, and --remove sets here.
LongPhenoData parseLongPheno(
    const std::string &filename,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::string &envName,
    const std::vector<std::string> &famIIDs,
    const std::unordered_set<std::string> &keptSubjects = {});

// One LMM fit's output.
struct LMMFit {
    Eigen::VectorXd beta;       // p (fixed-effects coefficients)
    Eigen::MatrixXd D;          // dRank × dRank random-effects covariance
    double sigma2;              // residual variance
    Eigen::VectorXd residPerRow; // N_total: Y - X β̂ - Z b̂
    int iterations;             // EM iterations actually taken
    double logLik;              // last evaluated ML log-likelihood
};

// Fit  y ~ X + (Z | IID)  with random Z = [1, E].
// X must contain an intercept column.
// y has length N_total.
LMMFit fitRandomSlopeML(
    const LongPhenoData &data,
    const Eigen::VectorXd &y,
    int maxIter = 10000,
    double tol = 1e-9);

// Fit  y ~ 1 + (1 | IID).
// y has length N_total.
LMMFit fitRandomInterceptML(
    const LongPhenoData &data,
    const Eigen::VectorXd &y,
    int maxIter = 10000,
    double tol = 1e-9);

// Σ_j r_ij — per-subject sum of per-row residuals.
// Output length = data.uniqueIIDs.size().
Eigen::VectorXd aggregatePerIID(
    const LongPhenoData &data,
    const Eigen::VectorXd &residPerRow);

// Σ_j w_ij r_ij — per-subject weighted sum of per-row residuals.
Eigen::VectorXd aggregateWeightedPerIID(
    const LongPhenoData &data,
    const Eigen::VectorXd &residPerRow,
    const Eigen::VectorXd &weightPerRow);

} // namespace nsSAGELDFit
