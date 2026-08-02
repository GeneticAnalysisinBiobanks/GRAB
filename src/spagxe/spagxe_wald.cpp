// spagxe_wald.cpp — Branch-B Wald leg of SPAGxE_CCT (per-marker G×E dispatcher).
//
// SPAGxE_CCT (Ma et al., Nat. Commun. 2025) combines, for the small fraction of
// variants with a significant marginal genetic effect (Branch B, p_marg ≤ ε), a
// retrospective saddlepoint G×E p-value with a prospective Wald p-value of the
// interaction coefficient in the FULL model  trait ~ [covar] + g + g:E , then
// reports  P = CCT(p_spa, p_wald)  (Cauchy combination; the log-domain
// math::cauchyCombineLog10, over the two magnitudes).
//
// This file assembles the full-interaction design  M = [covar | g | g∘E]  per
// marker (the interaction g∘E is the appended last column), drops incomplete-
// case rows, and calls the matching standard-model Wald fitter in
// src/util/wald.hpp to obtain the two-sided Wald p-value of the G:E coefficient.
// The fitters themselves (OLS / logistic IRLS / Breslow-Cox Newton /
// proportional-odds Fisher scoring) live in namespace wald so any method can
// share them without depending on SPAGxE.
//
// A singular information matrix, non-convergence, or a degenerate design returns
// NaN, which the caller folds into the Cauchy combination as a fall-back to the
// SPA p (never as p = 0 — GRAB2 convention; model doc §7 B4/B7).

#include "spagxe/spagxe_wald.hpp"

#include "util/wald.hpp" // wald::lastCoef*Pval, wald::OrdinalInfo

#include <cmath>
#include <limits>
#include <vector>

namespace spagxe_wald {

namespace {
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();
} // namespace

// ══════════════════════════════════════════════════════════════════════════
// Dispatcher — assemble  M = [covar | g | g∘E], drop NaN rows, fit per trait.
// ══════════════════════════════════════════════════════════════════════════

double waldInteractionPval(
    const WaldData &wd,
    const Eigen::Ref<const Eigen::VectorXd> &g,
    const Eigen::Ref<const Eigen::VectorXd> &E,
    wald::OrdinalInfo ordInfo
) {
    if (wd.trait == TraitType::None) return kNaN;
    const Eigen::Index n = wd.covar.rows();
    const Eigen::Index c = wd.covar.cols();
    if (g.size() != n || E.size() != n) return kNaN;

    // M = [covar | g | g∘E]; the interaction g∘E is the last column.
    Eigen::MatrixXd M(n, c + 2);
    if (c > 0) M.leftCols(c) = wd.covar;
    M.col(c) = g;
    M.col(c + 1) = g.array() * E.array();

    const bool isCox = (wd.trait == TraitType::Cox);

    // Complete-case rows: no NaN in M or in the trait's response(s).
    std::vector<Eigen::Index> keep;
    keep.reserve(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i) {
        bool ok = M.row(i).allFinite();
        if (ok) {
            if (isCox)
                ok = std::isfinite(wd.time[i]) && std::isfinite(wd.event[i]);
            else
                ok = std::isfinite(wd.y[i]);
        }
        if (ok) keep.push_back(i);
    }
    const Eigen::Index m = static_cast<Eigen::Index>(keep.size());
    if (m == 0) return kNaN;

    Eigen::MatrixXd Mf(m, M.cols());
    for (Eigen::Index i = 0; i < m; ++i) Mf.row(i) = M.row(keep[i]);

    switch (wd.trait) {
    case TraitType::Linear:
    case TraitType::Logistic: {
        // Prepend the intercept: Z = [1 | covar | g | g∘E].
        Eigen::MatrixXd Z(m, M.cols() + 1);
        Z.col(0).setOnes();
        Z.rightCols(M.cols()) = Mf;
        Eigen::VectorXd yf(m);
        for (Eigen::Index i = 0; i < m; ++i) yf[i] = wd.y[keep[i]];
        return (wd.trait == TraitType::Linear) ? wald::lastCoefLinearPval(Z, yf)
                                               : wald::lastCoefLogisticPval(Z, yf);
    }
    case TraitType::Cox: {
        Eigen::VectorXd tf(m), ef(m);
        for (Eigen::Index i = 0; i < m; ++i) {
            tf[i] = wd.time[keep[i]];
            ef[i] = wd.event[keep[i]];
        }
        return wald::lastCoefCoxPval(tf, ef, Mf);
    }
    case TraitType::Ordinal: {
        Eigen::VectorXi yf(m);
        for (Eigen::Index i = 0; i < m; ++i)
            yf[i] = static_cast<int>(std::lround(wd.y[keep[i]]));
        return wald::lastCoefOrdinalPval(yf, Mf, ordInfo);
    }
    default:
        return kNaN;
    }
}

} // namespace spagxe_wald
