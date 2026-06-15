// precond.hpp — SADGE stage-2 per-phenotype preconditioning.
//
// Port of var_geno_par_func / compute.rho.func / precond from
// example_sadge_zouyu/score_func_vec.R.  For one phenotype's residuals it
// computes the per-family-structure residual aggregates (sum_square, sum_prod),
// the rho(mu) interpolation grid, and the fs2 (phenotyped-structure) residual
// split that drives S_par.  All quantities are residual-side and per-phenotype;
// the per-marker mu / H / I / J are applied later in the score.
//
// fs1 vs fs2: the imputation kernels use the genotyped structure (fs1, fixed),
// but the variance and the singleton classification use the phenotyped
// structure (fs2 = genotyped-parents × this-phenotype-present-siblings), which
// is rebuilt here per phenotype.
#pragma once

#include "sadge/sample_info.hpp"

#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <vector>

namespace sadge {

// One family's union members, used to rebuild the phenotyped (stage-2)
// structure: genoNPar is fixed (genotyped parents), while the effective
// sibling count is the number of this-phenotype non-NaN members.
struct FamilyGroup {
    int genoNPar;          // 0/1/2
    int n;                 // number of union siblings in this family (1 or 2)
    std::array<int, 2> u;  // union indices of those siblings (u[1] valid iff n==2)
};

// var_geno_par_func(mu,H,I,J): per-structure parental-genotype variance,
// out[fs] = factor[fs] * mu*(1-mu) with factor = {2,H,4,J,I,4}.
void varGenoPar(double mu, double H, double I, double J,
                std::array<double, kNStruct> &out);

// compute.rho.func: rho(mu) given the per-structure sum_square weights, using
// the analytic (mu-derived) parental-genotype variance.
double computeRho(double mu, const std::array<double, kNStruct> &sumSquare);

// Read-only per-phenotype data shared across the engine's per-thread method
// clones (via shared_ptr).
struct PerPhenoData {
    Eigen::VectorXd residUnion;        // nUnion; 0 at subjects absent from this pheno (GEMM col 0)
    Eigen::VectorXd residSing;         // nUnion; residUnion at fs2-singletons, else 0 (GEMM col 1)
    Eigen::VectorXd residNonSingComp;  // nNonSing; residUnion at fs2-non-singletons (compact)
    double sumRSing = 0.0;             // sum of residSing
    double sumR = 0.0;                 // sum of residUnion (present subjects)
    double Rsq = 0.0;                 // sum over structures of sum_square
    double Rpr = 0.0;                 // sum over 2-sib structures of sum_prod
    std::array<double, kNStruct> sumSquare{};
    std::array<double, kNStruct> sumProd{};
    std::array<double, 1000> rhoGrid{}; // rho at mu = 0.5*k/999, k=0..999
};

// Build PerPhenoData from one phenotype's union residual vector (NaN at absent
// subjects), the family grouping over the union, and the compaction map for
// fs1-non-singleton siblings (unionToNonSing[u] = compact index, or -1).
PerPhenoData buildPerPhenoData(
    const Eigen::VectorXd &residUnionRaw,
    const std::vector<FamilyGroup> &families,
    uint32_t nUnion,
    const std::vector<int> &unionToNonSing,
    int nNonSing);

// Linear interpolation of the rho grid at mu (rule=2 clamping to [0,0.5]).
inline double interpRho(const std::array<double, 1000> &grid, double mu) {
    double x = mu;
    if (x < 0.0) x = 0.0;
    if (x > 0.5) x = 0.5;
    const double pos = x * (999.0 / 0.5); // grid spacing 0.5/999
    int i0 = static_cast<int>(pos);
    if (i0 >= 999) return grid[999];
    const double frac = pos - static_cast<double>(i0);
    return grid[i0] * (1.0 - frac) + grid[i0 + 1] * frac;
}

} // namespace sadge
