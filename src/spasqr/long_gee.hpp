// long_gee.hpp — marginal GEE null models for LoQus
//
// A long-format phenotype has m_i records per subject.  The null model is a
// marginal GEE with a working within-subject correlation R_i, and each subject
// is reduced to ONE score weight
//
//     a_i = 1ᵀ R_i⁻¹ r_i ,
//
// which the SPAsqr score machinery then treats exactly like a per-subject
// residual (the genotype is constant within subject, so the SNP score is
// Σ_i a_i G_i).  Two record-level residuals r_it are supported:
//
//   quantile  r_it = ψ_it = τ − Φ(−e_it / h)   (smoothed quantile GEE)
//   linear    r_it = e_it                        (identity-link GEE)
//
// with e_it = y_it − x_itᵀβ.  Working correlation: independence (R = I) or
// exchangeable, R_i = (1−ρ) I + ρ 11ᵀ.  For exchangeable the products needed
// have a closed form, R_i⁻¹ v = (v − c_i 11ᵀ v)/(1−ρ), c_i = ρ/(1+(m_i−1)ρ),
// so a_i = Σ_t r_it / (1 + (m_i−1)ρ).
//
// Fitting alternates "solve the estimating equation for fixed ρ" with "re-
// estimate ρ by moments from the current residuals" until both settle.  ρ is
// always estimated from the residual the estimating equation uses (ψ for the
// quantile model, e for the linear one): the optimal subject weight depends on
// Cor(ψ), not on the correlation of the raw errors, and ψ is bounded and a
// continuous function of β, so the iteration has no limit cycles.
//
// The quantile estimating equation with R ≠ I is not the gradient of any loss
// (its Jacobian Σ Zᵢᵀ R_i⁻¹ D_i Z_i is not symmetric), so QMME does not apply;
// it is solved by Newton on S with step halving on ‖S‖.  The independence fit
// IS the smoothed QR, and is taken from the QMME solver.
//
// Failure is named: if the exchangeable fit cannot be completed — ρ̂ outside
// the admissible range (−1/(m_max−1), 1), a singular Jacobian, or no
// convergence — the phenotype falls back to the independence weights, which
// are always valid (a working correlation only redistributes weight; the
// retrospective test is valid for any fixed weights), and NullFit.fellBack /
// NullFit.reason say so.  ρ is never clamped.
#pragma once

#include "spasqr/qmme.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <string>
#include <vector>

namespace longgee {

enum class WorkingCorr { Independence, Exchangeable };

// One phenotype's records grouped by subject: subject i owns rows
// [start[i], start[i+1]).  Every subject has at least one record.
struct Records {
    Eigen::MatrixXd X;              // nRows × p covariates, no intercept column
    Eigen::VectorXd y;              // nRows
    std::vector<uint32_t> start;    // nSubj + 1

    int nSubj() const { return static_cast<int>(start.size()) - 1; }
    int maxRecords() const;
};

struct NullFit {
    Eigen::VectorXd weight;   // nSubj: a_i = 1ᵀ R_i⁻¹ r_i (overall scale is irrelevant)
    double rho = 0.0;         // working correlation in the final weights (0 = independence)
    int outer = 0;            // solve / update-ρ rounds (0 for independence)
    bool fellBack = false;    // exchangeable requested but independence weights returned
    std::string reason;       // why, when fellBack
    bool qmmeConverged = true;
};

// Rank check of the null-model design [1 | X] on this record set.  Throws
// std::runtime_error, prefixed with `context`, naming the first column of X
// (labelled by `names`) that is constant over the records or a linear
// combination of the intercept and the columns before it (1 − R² < 1e-10 on
// the standardized columns).  Both GEE fits assume a full-rank design; a
// deficient one is a user input error, reported once up front rather than as
// a singular solve in one model and a silent ridge fit in the other.
void checkDesign(
    const Records &rec,
    const std::vector<std::string> &names,
    const std::string &context
);

// Whether y has variation outside the intercept/covariate span at working
// precision. A perfectly explained outcome has no informative score; tiny
// solver residuals must not be normalized into a unit-variance test.
bool hasResidualVariation(const Records &rec);

// Smoothed quantile GEE.  `solver` must be built on rec.X with its bandwidth
// prepared at h; it provides the independence (smoothed QR) fit.
NullFit fitQuantile(
    const Records &rec,
    const qmme::SqrSolver &solver,
    double tau,
    double h,
    WorkingCorr corr,
    double qmmeTol
);

// Identity-link GEE (Gaussian working model), ρ by Pearson moments of e.
NullFit fitLinear(
    const Records &rec,
    WorkingCorr corr
);

} // namespace longgee
