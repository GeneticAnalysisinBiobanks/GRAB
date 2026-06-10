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

// Fit  y ~ X + (1 | IID)  with the full fixed-effects design X = data.X
// (intercept in column 0) and a single subject random intercept.  Returns
// per-row BLUP residuals r_ij = y_ij − X_ij β̂ − b̂_i; aggregatePerIID then
// yields R_G_i = Σ_j r_ij.  Drives the --longitudinal main-effect input mode
// of SPACox / SPAmix / SPAGRM.  The p = 1 case reproduces
// fitRandomInterceptML.
LMMFit fitRandomInterceptWithCovar(
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

// GallopCache — per-marker variance-ratio λ estimation machinery.
//
// Replicates the precomputed objects in R/SAGELD.R:175-228:
//   • per-subject Si = TT_i' TT_i, StS_i = SS_i' SS_i  (TT = [1, E])
//   • XTs[i,] = vec(X_i' TT_i)                  (flattened to row, length 2·nx)
//   • AtS[i,] = vec(A21_i' SS_i)                (flattened to row, length 2·nx)
//   • Q = X'X − Σ_i A21_i' A21_i                (nx × nx, LDLT-factored)
// where Rot_i = sqrt(1/sv.d) sv.u for sv = svd(Si + P), SS_i = Rot_i Si,
// A21_i = Rot_i TT_i' X_i, and P = σ̂² · D̂⁻¹ (D̂ from fitRandomSlopeML).
//
// markerLambda(si, cache) → V[1,2] / V[1,1] for a single marker.
//
// The same cache also drives the exact GALLOP (Sikorska et al. 2013) per-marker
// Wald test via markerGallop().  V (the 2×2 reduced information computed for λ)
// equals σ²·I_γ, where I_γ = T_G' P T_G is the GALLOP information for the two
// genotype effects T_G = [G, G·E].  The Wald solve needs only the numerator
// ũ = (si' Tr0)ᵀ, where Tr0_i = TT_i' r0_i and r0 is the random-slope BLUP
// residual y − Xβ̂ − Zb̂ (LMMFit::residPerRow); then  β̂ = V⁻¹ũ  (the σ² cancels)
// and  SE_k = sig · √(V⁻¹)_kk  with sig = √σ².  See sageld_fit.cpp for the
// per-subject Woodbury derivation.
struct GallopCache {
    int nx = 0;        // ncol(X)
    int nSubj = 0;
    std::vector<Eigen::Matrix2d> Si;     // nSubj entries
    std::vector<Eigen::Matrix2d> StS;    // nSubj entries: SS_i' SS_i
    Eigen::MatrixXd XTs;                 // nSubj × (2·nx)
    Eigen::MatrixXd AtS;                 // nSubj × (2·nx)
    Eigen::MatrixXd Q;                   // nx × nx
    Eigen::LDLT<Eigen::MatrixXd> Qsolver;
    Eigen::MatrixXd Tr0;                 // nSubj × 2: per-subject TT_i' r0_i (GALLOP numerator)
    double sig = 0.0;                    // residual SD = √σ² (GALLOP SE scale)
};

GallopCache buildGallopCache(
    const LongPhenoData &data,
    const Eigen::MatrixXd &D,            // 2×2 random-effects covariance from fitMain
    double sigma2,                       // residual variance from fitMain
    const Eigen::VectorXd &residPerRow   // BLUP residual y − Xβ̂ − Zb̂ (LMMFit::residPerRow)
);

// Returns NaN if V[1,1] is non-positive (degenerate marker).
double markerLambda(
    const Eigen::Ref<const Eigen::VectorXd> &si,
    const GallopCache &cache
);

// Per-marker exact GALLOP Wald estimate for the two genotype effects
// [G (main), G×E (interaction)].  ok=false when V is not invertible (degenerate
// marker) — caller NaN-fills.
struct GallopMarker {
    Eigen::Vector2d beta;   // (β_G, β_GxE)
    Eigen::Vector2d se;     // (SE_G, SE_GxE)
    bool ok = false;
};

GallopMarker markerGallop(
    const Eigen::Ref<const Eigen::VectorXd> &si,
    const GallopCache &cache
);

} // namespace nsSAGELDFit
