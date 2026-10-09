// spagxe.cpp — SPAGxE_CCT G×E test implementation.
//
// Base (unrelated) and SPAGxE+ (sparse-GRM) share one code path: the GRM Φ
// enters only through the retrospective score variance as a quadratic form
// xᵀΦx = SparseGRM::quadForm(x), and the base path is the
// Φ = identity special case (phiQuad → Σx²).  Per phenotype a single
// genotype-independent residual R is fit once (trait ~ X + E).  Per variant:
//   1. marginal screen  S_G = Σ G_i R_i ,  Var = 2q(1−q)·RᵀΦR  (always normal;
//      routes Branch A vs Branch B on p_marg vs ε).
//   2a. Branch A (p_marg > ε): λ_GRM = RᵀΦR_E / RᵀΦR, w = (E − λ_GRM)R,
//       S_GxE = Σ G_i w_i, Var = 2q(1−q)·wᵀΦw, SPA / normal hybrid.
//   2b. Branch B (p_marg ≤ ε): R̃ = R − [1,G]([1,G]ᵀ[1,G])⁻¹[1,G]ᵀR (unweighted
//       even for +), Ŝ_GxE = Σ G_i E_i R̃_i, Var = 2q(1−q)·(E∘R̃)ᵀΦ(E∘R̃).
// The tail SPA rescales the statistic by √(indepVar/grmVar) (SAIGE variance
// ratio, the paper's non-mix convention) and applies the independence-CGF SPA
// via spamix_cgf::twoSidedSpa with the IQR outlier / Gaussian split.  A
// degenerate saddlepoint yields NaN, never 0 (GRAB2 convention).  Reference:
// tmp/SPAGxECCT/.../R  (SPAGxE_CCT_one_SNP, SPA_G_one_SNP_homo_new,
// SPAGxE_Plus_one_SNP, SPAGxE_Plus_Nullmodel).

#include "spagxe/spagxe.hpp"

#include "geno_factory/geno_data.hpp" // parseGenoIIDs, makeGenoData, GenoMeta
#include "io/sparse_grm.hpp"          // SparseGRM
#include "io/subject_data.hpp"        // SubjectData, extractPhenoVec/Mat, PerPhenoInfo
#include "spagxe/spagxe_wald.hpp"     // spagxe_wald::{WaldData, waldInteractionLog10P}
#include "spamix/spamix_cgf.hpp"      // spamix_cgf::{Context, twoSidedSpa, statusCode}
#include "spamix/indiv_af.hpp"        // AFContext, computeAFVec (SPAGxEmix)
#include "util/outlier.hpp"           // OutlierData, detectOutliers
#include "util/logging.hpp"           // infoMsg, warnMsg
#include "util/math_helper.hpp"       // math::pnorm, zFromNegLog10P, cauchyCombineLog10
#include "util/spa.hpp"               // spa::normalTwoSidedLog, spa::Status
#include "util/null_model.hpp"        // nullmodel::fitAll and friends

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Two-sided SPA p-value for the genotype score  S = Σ G_i w_i  under the
// homogeneous binomial law  G_i ~ Bin(2, q).  The GRM enters only through
// grmQuad = wᵀΦw (== wSq = Σw² in the base / unrelated case), which sets the
// score variance and — in the tail — a SAIGE-style variance-ratio rescaling of
// the statistic before the independence-CGF SPA (spamix_cgf::twoSidedSpa with
// the IQR outlier / Gaussian non-outlier split).  Mirrors SPA_G_one_SNP_homo_new /
// SPAGxE_Plus_one_SNP: normal approximation when |z| ≤ spaCutoff, otherwise the
// reflection-about-mean two-sided saddlepoint.
//   part     — IQR partition of w  (residOutlier = w[outlier],
//              resid2NonOutlier = w²[nonOutlier])
//   wSum     — Σ w      wSq — Σ w²      grmQuad — wᵀΦw   (all over every i)
//   rawScore — Σ G_i w_i
// Sets zScore (raw score z) and outVar = Var(S); returns the two-sided result
// (p, −log10 p, spa::Status).  p is NaN, with a status naming the reason,
// whenever the variance is non-positive or either tail's saddlepoint fails.
//
// spa_unify Stage 5.  The two `spa::getProbSpaG` calls whose sum was returned
// unguarded are replaced by spamix_cgf::twoSidedSpa — the shared bracketed
// solver, spa::bnTailLog's full guard set, and spa::combineTailsLog's single clamp and
// NaN policy.  The genotype law here has ONE q for every subject, so the
// binomial block is spa_cgf::binomUniform: 02_design.md lists SPAGxE under
// binomIndiv, but that is only because the predecessor materialized
// `std::vector<double> mafOut(nOut, q)` per marker per environment and passed
// the constant vector to the per-individual kernel.  That allocation is gone.
spa::Result spaScorePval(
    const OutlierData &part,
    double wSum,
    double wSq,
    double grmQuad,
    double rawScore,
    double q,
    double spaCutoff,
    double &zScore,
    double &outVar
) {
    const double f = 2.0 * q * (1.0 - q);
    const double scoreVar = f * grmQuad; // GRM variance (== f·wSq in the base case)
    if (!(scoreVar > 0.0)) {             // monomorphic / degenerate weight
        zScore = 0.0;
        outVar = 0.0;
        return spa::Result{kNaN, spa::Status::NaNoTest};
    }
    const double sMean = 2.0 * q * wSum; // retrospective mean E[S] = 2q·Σw
    zScore = (rawScore - sMean) / std::sqrt(scoreVar);
    outVar = scoreVar;
    if (!std::isfinite(zScore))
        return spa::Result{kNaN, spa::Status::NaNoTest};
    // spa::normalBranch reports -log10(2*Phi(-|z|)) in the log domain, so the
    // branch keeps its magnitude past |z| ~ 38.5 where the predecessor's
    // `2 * pnorm(-|z|)` was exactly zero (log10p_unify Stage 3).
    if (std::abs(zScore) <= spaCutoff) return spa::normalBranch(zScore);

    // ── Saddlepoint tail: variance-ratio rescale + outlier / non-outlier split.
    // sqrtRatio = √(indepVar / grmVar) = √(wSq / grmQuad) (the 2q(1−q) cancels);
    // == 1 in the base case, so this reduces to the plain independence SPA.
    const double sqrtRatio = (grmQuad > 0.0) ? std::sqrt(wSq / grmQuad) : 1.0;
    const double absDev = std::abs(rawScore - sMean) * sqrtRatio;

    // Gaussian non-outlier block by complement: Σ_non w = wSum − Σ_out w and
    // Σ_non w² = wSq − Σ_out w², both already available over the whole cohort.
    // The predecessor reduced over the non-outlier sub-vectors instead, which is
    // O(N) per marker per environment — and in Branch A those two reductions are
    // marker-independent, so it was O(N) work repeated unchanged 20000 times.
    // See the longer note in src/spamix/spamixplus.cpp::markerPvalFromAF.
    const double outSum = part.residOutlier.sum();
    const double outSq = part.residOutlier.squaredNorm();

    spamix_cgf::Context cgf;
    cgf.resid = part.residOutlier.data();
    cgf.af = nullptr;                                // uniform q: binomUniform
    cgf.nOutlier = static_cast<int>(part.posOutlier.size());
    cgf.q = q;
    cgf.mean = 2.0 * q * (wSum - outSum);
    cgf.var = f * (wSq - outSq);
    if (!(cgf.var > 0.0)) cgf.var = 0.0;

    // f·wSq is the independence K''(0) and only sizes the initial abscissa.
    return spamix_cgf::twoSidedSpa(cgf, sMean, absDev, f * wSq, zScore);
}

// Two-sided SPA p-value for the genotype score  S = Σ G_i w_i  under the
// per-individual binomial law  G_i ~ Bin(2, q̂_i) (SPAGxEmix).  This is the
// SPAmix / SPAmixPlus kernel (src/spamix/spamixplus.cpp::markerPvalFromAF, the
// no-GRM path) restated for an arbitrary weight vector w:
//   mean  = 2·Σ w_i q̂_i          var = Σ w_i²·2q̂_i(1−q̂_i)
// and, in the tail, the reflection-about-mean two-sided saddlepoint over the
// IQR outlier partition of w, with a Gaussian block for the non-outliers.
// There is NO GRM and hence no variance-ratio rescaling (mix stays diagonal).
//   part   — IQR partition of w  (residOutlier = w[outlier])
//   w      — the full weight vector (length N)
//   afVec  — per-individual q̂_i    wVec — 2 q̂_i(1−q̂_i)
//   s      — the observed score Σ G_i w_i
//   afOut  — caller-owned scratch, resized here to nOutlier and filled with
//            the outlier-position gather of afVec (the per-clone buffer that
//            replaces the predecessor's per-marker std::vector)
// Sets zScore (raw score z) and outVar = Var(S); returns NaN with a naming
// status on non-positive variance (GRAB2 convention: never 0).
//
// This is the per-individual genotype law, so the binomial block is
// spa_cgf::binomIndiv.
spa::Result spaScorePvalMix(
    const OutlierData &part,
    const Eigen::VectorXd &w,
    const Eigen::VectorXd &afVec,
    const Eigen::VectorXd &wVec,
    double s,
    double spaCutoff,
    std::vector<double> &afOut,
    double &zScore,
    double &outVar
) {
    const double sMean = 2.0 * w.dot(afVec);
    const double scoreVar = (w.array().square() * wVec.array()).sum();
    if (!(scoreVar > 0.0)) {
        zScore = 0.0;
        outVar = 0.0;
        return spa::Result{kNaN, spa::Status::NaNoTest};
    }
    zScore = (s - sMean) / std::sqrt(scoreVar);
    outVar = scoreVar;
    if (!std::isfinite(zScore))
        return spa::Result{kNaN, spa::Status::NaNoTest};
    if (std::abs(zScore) <= spaCutoff) return spa::normalBranch(zScore);

    // Gaussian non-outlier block with per-individual q̂_i, by complement:
    // sMean and scoreVar above are the whole-cohort sums, so the block is the
    // total minus an O(nOutlier) sum folded into the gather that has to run
    // anyway.  The outlier variance terms use wVec[pos] rather than a
    // recomputed 2*af*(1-af), so they are bit-identical to scoreVar's.  See the
    // longer note in src/spamix/spamixplus.cpp::markerPvalFromAF.
    const int nOut = static_cast<int>(part.posOutlier.size());
    if (static_cast<int>(afOut.size()) < nOut) afOut.resize(static_cast<size_t>(nOut));
    double outMeanSum = 0.0, outVarSum = 0.0;
    for (int i = 0; i < nOut; ++i) {
        const uint32_t pos = part.posOutlier[i];
        const double af = afVec[pos];
        const double wi = part.residOutlier[i];
        afOut[static_cast<size_t>(i)] = af;
        outMeanSum += wi * af;
        outVarSum += wi * wi * wVec[pos];
    }

    const double meanNon = sMean - 2.0 * outMeanSum;
    double varNon = scoreVar - outVarSum;
    if (!(varNon > 0.0)) varNon = 0.0;

    spamix_cgf::Context cgf;
    cgf.resid = part.residOutlier.data();
    cgf.af = afOut.data();
    cgf.nOutlier = nOut;
    cgf.mean = meanNon;
    cgf.var = varNon;

    // Reflect about the fitted mean: upper = sMean + |Δ|, lower = sMean − |Δ|.
    // There is no GRM here, so scoreVar IS the independence K''(0).
    return spamix_cgf::twoSidedSpa(cgf, sMean, std::abs(s - sMean), scoreVar, zScore);
}

// Resolve the trait family of one phenotype, matching nullmodel::fitAll's own
// dispatch: a Cox spec (TIME:EVENT) → Cox; otherwise the explicit
// --regression-model, or per-column inference in Auto mode.  Drives the
// Branch-B Wald refit; `None` disables it.
spagxe_wald::TraitType resolveTraitType(
    const SubjectData &sd,
    const nullmodel::PhenoSpec &spec,
    nullmodel::RegressionModel regModel
) {
    using TT = spagxe_wald::TraitType;
    using RM = nullmodel::RegressionModel;
    if (nullmodel::isCoxSpec(spec)) return TT::Cox;
    RM m = regModel;
    if (m == RM::Auto)
        m = nullmodel::inferModelFromColumn(sd.getColumn(spec.yColumn), spec.yColumn,
                                            sd.usedIIDs())
                .model;
    switch (m) {
    case RM::Linear:   return TT::Linear;
    case RM::Logistic: return TT::Logistic;
    case RM::Cox:      return TT::Cox;
    case RM::Ordinal:  return TT::Ordinal;
    default:           return TT::None;
    }
}

// Recode a raw phenotype column into the form the Wald fitter expects.  The two-
// sided interaction Wald p is invariant to the recode direction, so only the
// family matters: Logistic → {0,1} (below the midpoint → 0), Ordinal → {0..J−1}
// (shift by the minimum), Linear/Cox → unchanged (OLS is affine-invariant; the
// Cox response is time/event, handled separately).  NaN entries are preserved
// (the fitter drops those rows).
Eigen::VectorXd recodeResponse(const Eigen::VectorXd &col, spagxe_wald::TraitType t) {
    using TT = spagxe_wald::TraitType;
    if (t == TT::Logistic) {
        double lo = std::numeric_limits<double>::infinity();
        double hi = -std::numeric_limits<double>::infinity();
        for (Eigen::Index i = 0; i < col.size(); ++i)
            if (std::isfinite(col[i])) { lo = std::min(lo, col[i]); hi = std::max(hi, col[i]); }
        const double mid = 0.5 * (lo + hi);
        Eigen::VectorXd out(col.size());
        for (Eigen::Index i = 0; i < col.size(); ++i)
            out[i] = std::isnan(col[i]) ? col[i] : (col[i] > mid ? 1.0 : 0.0);
        return out;
    }
    if (t == TT::Ordinal) {
        double lo = std::numeric_limits<double>::infinity();
        for (Eigen::Index i = 0; i < col.size(); ++i)
            if (std::isfinite(col[i])) lo = std::min(lo, col[i]);
        Eigen::VectorXd out(col.size());
        for (Eigen::Index i = 0; i < col.size(); ++i)
            out[i] = std::isnan(col[i]) ? col[i] : std::round(col[i] - lo);
        return out;
    }
    return col; // Linear (affine-invariant); Cox uses time/event, not y
}

// Branch B's weight  u = E ∘ R̃, with the direction that contributes nothing to
// the statistic removed.  `VR0` is V·R̃ up to any positive scalar, with V the
// genotype covariance the caller's variance uses: Φ·R̃ under the uniform law
// (R̃ itself when Φ = I), and 2q̂(1−q̂) ∘ R̃ under SPAGxEmix's per-individual law.
//
// WHY THIS IS NOT A TUNING KNOB.  The Branch-B projection solves the [1, G]
// normal equations, so at the observed G
//
//     Σ_i R̃_i = 0        and        Σ_i G_i R̃_i = 0
//
// hold exactly, hence so does Σ_i (G_i − 2q) R̃_i = 0.  The score is therefore
// the same number for every a:
//
//     S = Σ_i (G_i − 2q)·u_i = Σ_i (G_i − 2q)·(u_i − a·R̃_i)
//
// while the plug-in variance u'V u is not.  Every member of that family is a
// representation of the SAME statistic, so choosing among them is not choosing
// among tests; it is choosing which representation the plug-in is applied to.
// The weight u = E ∘ R̃ generally has a component along R̃, and the plug-in
// charges the statistic for that component's variance even though it
// contributes exactly zero.  The minimiser over a,
//
//     a = (u' V R̃) / (R̃' V R̃)
//
// makes u ⊥_V R̃ and is therefore the tightest representation available.
//
// It is also the standard efficient-score correction for the nuisance
// parameters (α, β) the projection estimates: the correct variance is
// I_ββ − I_βη I_ηη⁻¹ I_ηβ and the naive plug-in is I_ββ, whose excess is a
// non-negative rank-1 term.  Omitting it is therefore always CONSERVATIVE and
// never anti-conservative, which is why the deficit had a consistent sign.
//
// HOW MUCH IT MOVES.  With E ⊥ R̃² the ratio of corrected to naive variance is
// Var(E)/E[E²] = 1/(1 + (Ē/SD_E)²).  An already-centred environment loses
// nothing; a 0/1 indicator of prevalence 0.5 has its Branch-B variance
// overstated by exactly 2.  `examples/baseline.sh` uses `MALE`, so its
// Branch-B columns sit at that end.  Shifting E by a constant leaves Y, R, R̃
// and S untouched and moves only this variance, which makes the attribution
// checkable: on chr22 of a 50 000-subject cohort λ_GC of the uncorrected
// Branch B ran 0.936 / 0.505 / 0.212 for E shifted by 0 / 1 / 2 against the
// predicted 1.000 / 0.500 / 0.200, while Branch A held at 0.945 throughout.
//
// On feat/sqrgxe this function is gxe_score::branchBWeight, shared with
// --method sqrgxe; here it is file-local because that tier does not exist yet
// on this branch.  Keep the two bodies identical.
void branchBWeight(
    const Eigen::Ref<const Eigen::VectorXd> &env,
    const Eigen::Ref<const Eigen::VectorXd> &R0,
    const Eigen::Ref<const Eigen::VectorXd> &VR0,
    Eigen::VectorXd &u
) {
    u = env.array() * R0.array();
    const double den = R0.dot(VR0);
    if (!(den > 0.0)) return; // no usable metric: leave the naive weight
    const double a = u.dot(VR0) / den;
    if (!std::isfinite(a)) return;
    u.noalias() -= a * R0;
}

} // namespace

// ══════════════════════════════════════════════════════════════════════
// SPAGxEMethod
// ══════════════════════════════════════════════════════════════════════

double SPAGxEMethod::phiQuad(const Eigen::VectorXd &x) const {
    if (m_grm)
        return m_grm->quadForm(x.data(), static_cast<uint32_t>(x.size()));
    return x.squaredNorm();
}

SPAGxEMethod::SPAGxEMethod(
    Eigen::VectorXd resid,
    std::vector<std::string> envNames,
    std::vector<Eigen::VectorXd> envVecs,
    std::shared_ptr<const SparseGRM> grm,
    std::shared_ptr<const spagxe_wald::WaldData> wald,
    double marginalCutoff,
    double spaCutoff,
    double outlierRatio,
    std::shared_ptr<const AFData> afData
)
    : m_resid(std::move(resid)),
      m_envNames(std::move(envNames)),
      m_grm(std::move(grm)),
      m_wald(std::move(wald)),
      m_afData(std::move(afData)),
      m_marginalCutoff(marginalCutoff),
      m_spaCutoff(spaCutoff),
      m_outlierRatio(outlierRatio)
{
    const Eigen::Index N = m_resid.size();
    m_resid2 = m_resid.array().square();
    m_residSum = m_resid.sum();
    m_R_phiQuad = phiQuad(m_resid); // RᵀΦR (= ΣR² in the base case)
    const int nEnv = static_cast<int>(m_envNames.size());
    const bool mix = static_cast<bool>(m_afData);

    if (mix) {
        m_afVec.resize(N);
        m_wVec.resize(N);
        m_wScratch.resize(N);
    }

    // m_W is the GEMM/GEMV operand: one product against a genotype batch yields
    // the marginal score (col 0) and every per-env score column at once.
    //   base:  m_W = [ R | w_1 | … | w_nEnv ],  w_e = (E_e − λ_e)∘R  (λ_e fixed)
    //   mix:   m_W = [ R | E_1∘R | … ],         λ_e is per-marker so only the
    //          naive weight E_e∘R is precomputable (S_GxE formed per marker).
    m_W.resize(N, 1 + nEnv);
    m_W.col(0) = m_resid;

    m_envs.reserve(nEnv);
    for (int e = 0; e < nEnv; ++e) {
        EnvData ed;
        ed.name = m_envNames[e];
        ed.E = std::move(envVecs[e]);

        if (mix) {
            // Per-individual AF: λ is variance-weighted and per-marker, so the
            // Branch-A weight cannot be precomputed.  Only the marker-independent
            // naive weight E∘R feeds the GEMV (rawScores[1+e] = Σ G_i E_i R_i).
            ed.lambda = 0.0;
            ed.wSum = ed.wSq = ed.wPhiQuad = 0.0;
            m_W.col(1 + e) = (ed.E.array() * m_resid.array()).matrix();
            m_envs.push_back(std::move(ed));
            continue;
        }

        // λ = RᵀΦR_E / RᵀΦR, R_E = R∘E.  RᵀΦR_E is the bilinear Φ-form obtained
        // by polarization (== Σ E R² in the base case; kept as the direct sum
        // there so the base numbers match the pure-Σ path exactly).
        double R_GRM_RE;
        if (m_grm) {
            const Eigen::VectorXd RE = (m_resid.array() * ed.E.array()).matrix();
            const Eigen::VectorXd RpRE = m_resid + RE;
            R_GRM_RE = 0.5 * (phiQuad(RpRE) - m_R_phiQuad - phiQuad(RE));
        } else {
            R_GRM_RE = (m_resid.array().square() * ed.E.array()).sum();
        }
        ed.lambda = (m_R_phiQuad > 0.0) ? R_GRM_RE / m_R_phiQuad : 0.0;
        ed.w = (ed.E.array() - ed.lambda) * m_resid.array();
        ed.wSum = ed.w.sum();
        ed.wSq = ed.w.squaredNorm();
        ed.wPhiQuad = phiQuad(ed.w); // wᵀΦw (= Σw² in the base case)
        ed.wOutlier = detectOutliers(ed.w, m_outlierRatio);
        m_W.col(1 + e) = ed.w;
        m_envs.push_back(std::move(ed));
    }
}

std::unique_ptr<MethodBase> SPAGxEMethod::clone() const {
    return std::make_unique<SPAGxEMethod>(*this);
}

std::string SPAGxEMethod::getHeaderColumns() const {
    // Leading marginal block LOG10P_G Z_G BETA_G SE_G SPA_STATUS_G
    // (normal-approx, so its Z is already p-consistent and needs no Z_Norm, and
    // its status is the constant 1 NORMAL wherever a statistic exists), then a
    // 7-wide G×E block per environment: the final p as the magnitude LOG10P_Gx
    // (CCT in Branch B, SPA otherwise); the Branch-B Wald magnitude
    // LOG10P_Wald_Gx (NaN when no Wald ran); the p-consistent Z and the
    // raw-score Z_Norm (which differ in the tails); BETA/SE; and the
    // saddlepoint outcome SPA_STATUS_Gx.  The linear P_G / P_Gx columns left in
    // log10p_unify Stage 8 (decision D1).
    std::string h = "\tLOG10P_G\tZ_G\tBETA_G\tSE_G\tSPA_STATUS_G";
    for (const auto &n : m_envNames)
        h += "\tLOG10P_Gx" + n + "\tLOG10P_Wald_Gx" + n + "\tZ_Gx" + n +
             "\tZ_Norm_Gx" + n + "\tBETA_Gx" + n + "\tSE_Gx" + n +
             "\tSPA_STATUS_Gx" + n;
    return h;
}

namespace {

// The seven cells of one environment's G×E block.
//
// `ts` is the saddlepoint (or normal-branch) result for the score test; `lWald`
// is the magnitude −log10 of the Branch-B prospective Wald p, NaN where no Wald
// leg ran.  When a Wald leg is present the reported p is the Cauchy combination
// of the two, and since log10p_unify Stage 5 that combination is taken in the
// log domain (math::cauchyCombineLog10) rather than on the linear scale.
// LOG10P_Gx is therefore the magnitude the CCT itself produces, NOT −log10 of a
// linear combination: the saddlepoint leg enters as its own L and is no longer
// truncated at the point where the linear p underflows.  Stage 8 removed the
// derived P_Gx column, so the magnitude is now the only thing reported.
//
// Since Stage 7 the Wald leg is a magnitude at the source as well
// (wald::lastCoef*Log10P through spagxe_wald::waldInteractionLog10P), so the
// last ceiling in this path is gone: a Wald tail below the linear underflow
// used to enter the combination as −log10(0) = +∞ and carry that infinity into
// LOG10P_Gx, which was a C1 violation waiting on an input extreme enough to
// trigger it.
//
// SPA_STATUS always describes the SADDLEPOINT, independently of what reached
// the p-value column — which matters precisely because the combination skips a
// NaN, so a failed saddlepoint with a successful Wald refit still yields a
// finite LOG10P_Gx and the status is the only record that the SPA leg dropped
// out.
void writeEnvBlock(
    std::vector<double> &out,
    int base,
    const spa::Result &ts,
    double lWald,
    bool waldEnabled,
    double z,
    double var,
    double score
) {
    out[base + 6] = spamix_cgf::statusCode(ts.status);
    if (!(var > 0.0)) return;

    double negLog = ts.negLog10p;
    if (waldEnabled) {
        // The CCT in the log domain (Stage 5), over two magnitudes (Stage 7).
        // The -0.0 normalization the linear form needed goes with it: the
        // combination of n copies of one p is that p, so with the
        // L >= -log10(0.999) clamp on every input the result is bounded below
        // by 4.34e-4 and cannot come back as a zero.
        const double Ls[2] = {ts.negLog10p, lWald};
        negLog = math::cauchyCombineLog10(Ls, 2);
    }
    out[base + 0] = negLog;                      // LOG10P_Gx<E> (final)
    out[base + 1] = lWald;                       // LOG10P_Wald_Gx<E> (NaN if none)
    // Z from LOG10P_Gx, not from P_Gx (Stage 4): the linear inversion
    // saturated at |Z| = 37.0470962993612 for every L >= 299.698970.  Since
    // Stage 5, `negLog` is the CCT's own magnitude, so the saturation is gone
    // from the combined branch as well rather than merely displaced onto it.
    out[base + 2] = math::zFromNegLog10P(negLog, z);  // Z_Gx<E> (p-consistent)
    out[base + 3] = z;                           // Z_Norm_Gx<E> (raw score z)
    out[base + 4] = score / var;                 // BETA_Gx<E>
    out[base + 5] = 1.0 / std::sqrt(var);        // SE_Gx<E>
}

}  // namespace

void SPAGxEMethod::evalMarker(
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    const Eigen::VectorXd &rawScores,
    std::vector<double> &out
) {
    const int nEnv = static_cast<int>(m_envNames.size());
    out.assign(static_cast<size_t>(5 + 7 * nEnv), kNaN);
    // Pre-set so that a marker leaving the marginal block empty still says why
    // (decision D4): no statistic exists there, it is not a failed one.
    out[4] = spamix_cgf::statusCode(spa::Status::NaNoTest);   // SPA_STATUS_G

    const double q = altFreq;

    // ── Marginal genetic block (always the normal approximation) ────────
    // S_G = Σ G_i R_i (uncentered, matching the reference S1; E[S_G]=2q·ΣR≈0),
    // variance 2q(1−q)·RᵀΦR.
    const double varSG = m_R_phiQuad * 2.0 * q * (1.0 - q);
    double pMarg = 1.0; // degenerate variance ⇒ take Branch A (no projection)
    if (varSG > 0.0) {
        const double sG = rawScores[0];
        const double sdG = std::sqrt(varSG);
        const double zG = sG / sdG;
        pMarg = 2.0 * math::pnorm(-std::abs(zG));
        // LOG10P_G in the log domain, not −log10(pMarg): the linear tail is
        // exactly zero past |zG| = 38.6 and its logarithm would be the +Inf
        // invariant C1 forbids.  pMarg itself stays linear because its only
        // other consumer is the Branch A / Branch B routing against
        // --spagxe-marginal-cutoff, an input threshold (decision D8).
        out[0] = -spa::normalTwoSidedLog(zG) / math::kLn10;   // LOG10P_G
        out[1] = zG;          // Z_G (normal ⇒ p-consistent)
        out[2] = sG / varSG;  // BETA_G
        out[3] = 1.0 / sdG;   // SE_G
        // 1 NORMAL: this block never attempts a saddlepoint, and decision D4
        // covers exactly that case.  Var(S_G) <= 0 leaves the pre-set
        // 8 NA_NO_TEST — there is no statistic to report, not a failed one.
        out[4] = spamix_cgf::statusCode(spa::Status::Normal);   // SPA_STATUS_G
    }

    const bool branchB = (pMarg <= m_marginalCutoff);

    // Branch B: genotype-adjusted residual R̃ = R − α − β·G, projecting the
    // marginal genetic effect out.  [α, β]ᵀ = ([1,G]ᵀ[1,G])⁻¹ [1,G]ᵀ R.
    // R̃ is env-independent, so it is formed once per marker.  (Unweighted even
    // for +; the GRM re-enters only through the variance — paper eq. 3.)
    Eigen::VectorXd R0;
    Eigen::VectorXd VR0; // V*R0, the metric branchBWeight orthogonalises in
    if (branchB) {
        const double n = static_cast<double>(m_resid.size());
        const double sg = GVec.sum();
        const double sgg = GVec.squaredNorm();
        const double det = n * sgg - sg * sg;
        if (det > 0.0) {
            const double sr = m_residSum;
            const double sgr = rawScores[0];
            const double alpha = (sgg * sr - sg * sgr) / det;
            const double beta = (n * sgr - sg * sr) / det;
            R0 = m_resid.array() - alpha - beta * GVec.array();
            // Φ·R̃, the metric branchBWeight projects in.  Env-independent like
            // R̃ itself, so it costs one GRM multiply per marker rather than
            // one per environment; Φ = I on the base path makes it a copy.
            if (m_grm) {
                VR0.resize(R0.size());
                m_grm->multiply(R0.data(), VR0.data(),
                                static_cast<uint32_t>(R0.size()));
            } else {
                VR0 = R0;
            }
        }
        // det ≤ 0 (monomorphic G): leave R0 empty ⇒ per-env NaN below.
    }

    // Branch B fires the prospective Wald leg only on the base (non-GRM) path
    // (SPAGxE+ keeps no Wald — paper) and only when raw phenotype data is
    // present (fit mode; residual mode leaves m_wald null → SPA-only).
    const bool waldEnabled = branchB && m_wald && !m_grm &&
                             m_wald->trait != spagxe_wald::TraitType::None;

    for (int e = 0; e < nEnv; ++e) {
        const int base = 5 + 7 * e;
        // A marker that never reaches the score test still gets a status: a
        // bare row of NA says nothing about why, and this family's NA are not
        // hypothetical (see the SPAGxEmix note in evalMarkerMix).
        out[base + 6] = spamix_cgf::statusCode(spa::Status::NaNoTest);

        double z = 0.0, var = 0.0, score;
        double lWald = kNaN;
        spa::Result ts{kNaN, spa::Status::NaNoTest};
        if (!branchB) {
            // Branch A: precomputed λ-orthogonalised weight and its Φ-form.
            const EnvData &ed = m_envs[e];
            score = rawScores[1 + e]; // Σ G_i w_{e,i}
            ts = spaScorePval(ed.wOutlier, ed.wSum, ed.wSq, ed.wPhiQuad, score, q,
                              m_spaCutoff, z, var);
        } else {
            if (R0.size() == 0) continue; // degenerate G ⇒ leave NaN
            // Branch B per-marker weight u = E ∘ R̃; Ŝ_GxE = Σ G_i u_i.
            Eigen::VectorXd u;
            branchBWeight(m_envs[e].E, R0, VR0, u);
            score = GVec.dot(u);
            const OutlierData ub = detectOutliers(u, m_outlierRatio);
            ts = spaScorePval(ub, u.sum(), u.squaredNorm(), phiQuad(u), score, q,
                              m_spaCutoff, z, var);
            // Prospective Wald p of the G:E coefficient in trait ~ X+E+g+g:E,
            // over the same environment used for the score weight.
            if (waldEnabled)
                lWald = spagxe_wald::waldInteractionLog10P(*m_wald, GVec, m_envs[e].E);
        }
        writeEnvBlock(out, base, ts, lWald, waldEnabled, z, var, score);
    }
}

// ══════════════════════════════════════════════════════════════════════
// evalMarkerMix — SPAGxEmix per-individual-AF marker evaluation
// ══════════════════════════════════════════════════════════════════════
//
// Mirrors evalMarker but with the per-individual genotype law G_i ~ Bin(2, q̂_i)
// (Ma et al. 2025, SPAGxEmix_CCT_one_SNP).  q̂_i is estimated per marker via the
// SPAmix AF cascade (computeAFVec); the marginal screen is re-centred by
// 2·Σ R_i q̂_i, λ is variance-weighted and recomputed per marker, and every SPA
// uses the per-individual afVec (spaScorePvalMix).  Branch B reuses the exact
// same genotype-adjusted-residual projection and the Phase-3 Wald + CCT leg.
void SPAGxEMethod::evalMarkerMix(
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    const Eigen::VectorXd &rawScores,
    std::vector<double> &out
) {
    const int nEnv = static_cast<int>(m_envNames.size());
    out.assign(static_cast<size_t>(5 + 7 * nEnv), kNaN);
    // Pre-set so that a marker leaving the marginal block empty still says why
    // (decision D4): no statistic exists there, it is not a failed one.
    out[4] = spamix_cgf::statusCode(spa::Status::NaNoTest);   // SPA_STATUS_G

    const double q = altFreq;

    // ── Per-individual ALT frequency q̂_i (SPAmix cascade) ──────────────
    // [1|PCs] design + OLS matrices from m_afData; identical to SPAmix.
    AFContext ctx{
        m_afData->onePlusPCs,
        m_afData->XtX_inv_Xt,
        m_afData->sqrt_XtX_inv_diag,
        m_afData->onePlusPCs.rightCols(m_afData->nPC),
        static_cast<int>(m_resid.size()),
        m_afData->nPC
    };
    computeAFVec(GVec, q, ctx, m_afVec);                       // q̂_i  → m_afVec
    m_wVec = 2.0 * m_afVec.array() * (1.0 - m_afVec.array());  // 2q̂(1−q̂) → m_wVec

    // ── Marginal genetic block (always the normal approximation) ────────
    // S_G = Σ G_i R_i, re-centred: E[S_G] = 2·Σ R_i q̂_i, Var = Σ R_i²·2q̂(1−q̂).
    // varSG = Σ R_i²·2q̂_i(1−q̂_i) can be exactly zero, and is so for three of
    // the 3000 markers in examples/1kg — the whole row is NA in the pinned
    // baseline, and has been since before this rework.  The cause is NOT the
    // [0,1] clamp in computeAFVec's OLS branch, which is where both the audit
    // and its re-verification place it; all three markers take the STATUS-2
    // logistic branch (indiv_af.cpp:189).  Their ALT frequency is ~0.995, so
    // the binarised response g0 = 1{G > 0.5} is one for essentially every
    // subject, the logistic MLE is at +∞, and IRLS stops at its iteration cap
    // having driven the intercept to ~101.2.  sigmoid(101.2) rounds to exactly
    // 1.0 in double, so AF = 1 − sqrt(1 − μ) = 1.0 exactly for all 2504
    // subjects, every q̂(1−q̂) is zero and varSG vanishes.  It is complete
    // separation in the AF model, not a clamp.
    //
    // That is a modelling defect in the AF cascade, not a saddlepoint one — the
    // saddlepoint is never reached on these markers — so this stage preserves
    // the values exactly and only makes the failure legible: SPA_STATUS_Gx is
    // 8 NA_NO_TEST rather than a bare unexplained NA.  Repairing it (a
    // separation-aware AF fit, or a Firth penalty) changes the statistic for
    // every high-frequency marker and belongs to whoever owns indiv_af.cpp.
    const double sG = rawScores[0];
    const double meanSG = 2.0 * m_resid.dot(m_afVec);
    const double varSG = m_resid2.dot(m_wVec);
    double pMarg = 1.0; // degenerate variance ⇒ take Branch A (no projection)
    if (varSG > 0.0) {
        const double sdG = std::sqrt(varSG);
        const double zG = (sG - meanSG) / sdG;
        pMarg = 2.0 * math::pnorm(-std::abs(zG));
        // The magnitude in the log domain; see the base path above.
        out[0] = -spa::normalTwoSidedLog(zG) / math::kLn10;   // LOG10P_G
        out[1] = zG;                    // Z_G (normal ⇒ p-consistent)
        out[2] = (sG - meanSG) / varSG; // BETA_G
        out[3] = 1.0 / sdG;             // SE_G
        out[4] = spamix_cgf::statusCode(spa::Status::Normal);   // SPA_STATUS_G
    }

    const bool branchB = (pMarg <= m_marginalCutoff);

    // Branch B: genotype-adjusted residual R̃ = R − α − β·G (env-independent).
    Eigen::VectorXd R0;
    Eigen::VectorXd VR0; // V*R0, the metric branchBWeight orthogonalises in
    if (branchB) {
        const double n = static_cast<double>(m_resid.size());
        const double sg = GVec.sum();
        const double sgg = GVec.squaredNorm();
        const double det = n * sgg - sg * sg;
        if (det > 0.0) {
            const double sr = m_residSum;
            const double sgr = sG;
            const double alpha = (sgg * sr - sg * sgr) / det;
            const double beta = (n * sgr - sg * sr) / det;
            R0 = m_resid.array() - alpha - beta * GVec.array();
            // V*R0 for the per-individual law, whose score variance is
            // sum u_i^2 * 2 qhat_i (1-qhat_i); the metric is diag(m_wVec), not
            // Phi (mix carries no GRM).
            VR0 = (m_wVec.array() * R0.array()).matrix();
        }
    }

    // Mix has no GRM (m_grm is null); the Wald leg fires on the same fit-mode
    // base condition as evalMarker.
    const bool waldEnabled = branchB && m_wald && !m_grm &&
                             m_wald->trait != spagxe_wald::TraitType::None;

    for (int e = 0; e < nEnv; ++e) {
        const int base = 5 + 7 * e;
        out[base + 6] = spamix_cgf::statusCode(spa::Status::NaNoTest);

        const EnvData &ed = m_envs[e];
        double z = 0.0, var = 0.0, score;
        double lWald = kNaN;
        spa::Result ts{kNaN, spa::Status::NaNoTest};
        if (!branchB) {
            // Branch A: variance-weighted λ (per marker) on the per-individual q̂.
            //   λ = Σ 2q̂(1−q̂) E R² / Σ 2q̂(1−q̂) R²
            const double denom = (m_wVec.array() * m_resid2.array()).sum();
            const double numer = (m_wVec.array() * ed.E.array() * m_resid2.array()).sum();
            const double lambda = (denom > 0.0) ? numer / denom : 0.0;
            // S_GxE = S2 − λ·S_G  (S2 = Σ G_i E_i R_i from the GEMV col 1+e).
            score = rawScores[1 + e] - lambda * sG;
            m_wScratch = (ed.E.array() - lambda) * m_resid.array(); // weight (E−λ)∘R
            const OutlierData part = detectOutliers(m_wScratch, m_outlierRatio);
            ts = spaScorePvalMix(part, m_wScratch, m_afVec, m_wVec, score,
                                 m_spaCutoff, m_afOutlier, z, var);
        } else {
            if (R0.size() == 0) continue; // degenerate G ⇒ leave NaN
            // Branch B per-marker weight u = E ∘ R̃; Ŝ_GxE = Σ G_i u_i.
            Eigen::VectorXd u;
            branchBWeight(ed.E, R0, VR0, u);
            score = GVec.dot(u);
            const OutlierData part = detectOutliers(u, m_outlierRatio);
            ts = spaScorePvalMix(part, u, m_afVec, m_wVec, score, m_spaCutoff,
                                 m_afOutlier, z, var);
            if (waldEnabled)
                lWald = spagxe_wald::waldInteractionLog10P(*m_wald, GVec, ed.E);
        }
        writeEnvBlock(out, base, ts, lWald, waldEnabled, z, var, score);
    }
}

void SPAGxEMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int /*markerInChunkIdx*/,
    std::vector<double> &result
) {
    const Eigen::VectorXd rawScores = m_W.transpose() * GVec;
    if (m_afData)
        evalMarkerMix(GVec, altFreq, rawScores, result);
    else
        evalMarker(GVec, altFreq, rawScores, result);
}

void SPAGxEMethod::getResultBatch(
    const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
    const std::vector<double> &altFreqs,
    const std::vector<int> & /*chunkIdxs*/,
    std::vector<std::vector<double> > &results
) {
    const int B = static_cast<int>(GBatch.cols());
    results.resize(B);
    // One GEMM supplies the marginal (col 0) and every per-env score column.
    //   base:  scores(b, 1+e) = Σ G_b·w_e         (Branch-A λ-orthogonalised)
    //   mix:   scores(b, 1+e) = Σ G_b·(E_e∘R)     (naive S2; λ applied per marker)
    const Eigen::MatrixXd scores = GBatch.transpose() * m_W; // B × (1+nEnv)
    const bool mix = static_cast<bool>(m_afData);
    for (int b = 0; b < B; ++b) {
        const Eigen::VectorXd rs = scores.row(b).transpose();
        if (mix)
            evalMarkerMix(GBatch.col(b), altFreqs[b], rs, results[b]);
        else
            evalMarker(GBatch.col(b), altFreqs[b], rs, results[b]);
    }
}

// ══════════════════════════════════════════════════════════════════════
// runSPAGxEImpl — full workflow (mirrors runSPAGRM); serves both
//   --method spagxe    (base / SPAGxE+; pcColNames empty)
//   --method spagxemix (per-individual AF; pcColNames = the --pc-cols columns)
// ══════════════════════════════════════════════════════════════════════

namespace {

void runSPAGxEImpl(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::vector<std::string> &pcColNames, // empty → base/+, non-empty → mix
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double marginalCutoff,
    double spaCutoff,
    double outlierIqrRatio,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    bool saveResid
) {
    const bool mix = !pcColNames.empty();
    const bool hasGrm = !spgrmGrabFile.empty() || !spgrmGctaFile.empty();
    // Display label (messages) and output-file suffix.  Base and SPAGxE+ share
    // the "SPAGxE" file suffix; the "+" is display only.  Mix writes ".SPAGxEmix".
    const std::string dispLabel = mix ? "SPAGxEmix" : (hasGrm ? "SPAGxE+" : "SPAGxE");
    const char *fileSuffix = mix ? "SPAGxEmix" : "SPAGxE";
    const int nPC = static_cast<int>(pcColNames.size());

    const bool fitPath = !phenoNameSpec.empty();
    nullmodel::RegressionModel regModel{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        regModel = nullmodel::parseRegressionModel(regressionModelStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
        infoMsg("%s: fitting %s null model for %zu phenotype(s)",
                dispLabel.c_str(), nullmodel::regressionModelName(regModel),
                phenoSpecs.size());
    }

    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (fitPath) {
        sd.loadPhenoFile(phenoFile, nullmodel::columnsNeeded(phenoSpecs));
    } else {
        sd.loadResidOne(phenoFile, residNames);
        if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile);
    }
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    // (+) variant: restrict to subjects present in the sparse GRM (the GRM
    // supplies the retrospective genotype covariance among relatives).
    if (hasGrm)
        sd.setGrmSubjects(
            SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();

    // Environment columns must enter the null model in fit mode (so R ⟂ E),
    // exactly as SAGELD requires --envir-name ⊆ --covar-name.
    if (fitPath) {
        for (const auto &en : envNames)
            if (std::find(covarNames.begin(), covarNames.end(), en) == covarNames.end())
                throw std::runtime_error(
                    dispLabel + ": environment '" + en +
                    "' must also appear in --covar-name (it enters the "
                    "genotype-independent null model).");
        // SPAGxEmix: the PC columns feeding the per-individual AF model must
        // also adjust the null model (matching SPAGxEmix_CCT's topPCs ⊆ Cova).
        if (mix)
            for (const auto &pc : pcColNames)
                if (std::find(covarNames.begin(), covarNames.end(), pc) == covarNames.end())
                    throw std::runtime_error(
                        dispLabel + ": PC column '" + pc +
                        "' must also appear in --covar-name (the per-individual "
                        "allele-frequency PCs adjust the null model).");
    }
    // Drop subjects with a missing environment value (a no-op in fit mode,
    // where E is a covariate; necessary in residual mode, where E is only
    // used to form λ and is not otherwise NA-filtered).
    sd.dropNaInColumns(envNames);
    if (mix) sd.dropNaInColumns(pcColNames); // AF design must be complete
    const uint32_t N = sd.nUsed();
    infoMsg("  %u subjects in union mask", N);

    // Covariate union (PCs + env main effects).  In fit mode it is both the
    // null-model design and — reused per marker — the Branch-B Wald design; the
    // per-phenotype trait family (below) selects the Wald fitter.
    Eigen::MatrixXd covarUnion;
    std::vector<spagxe_wald::TraitType> traitTypes;
    if (fitPath) {
        if (!covarNames.empty()) covarUnion = sd.getColumns(covarNames);
        else covarUnion.resize(N, 0);
        traitTypes.reserve(phenoSpecs.size());
        for (const auto &spec : phenoSpecs)
            traitTypes.push_back(resolveTraitType(sd, spec, regModel));
    }

    // SPAGxEmix: union-level [1 | PCs] design for the per-individual AF cascade.
    Eigen::MatrixXd unionOnePlusPCs;
    if (mix) {
        unionOnePlusPCs.resize(N, 1 + nPC);
        unionOnePlusPCs.col(0).setOnes();
        unionOnePlusPCs.rightCols(nPC) = sd.getColumns(pcColNames);
        infoMsg("SPAGxEmix: per-individual AF from %d PC(s)", nPC);
    }

    if (fitPath) {
        nullmodel::EngineOptions eo;
        eo.nthreads = nthreads;
        auto fits = nullmodel::fitAll(sd, phenoSpecs, regModel, covarUnion, eo);
        std::vector<Eigen::VectorXd> rs;
        std::vector<std::string> ns;
        rs.reserve(fits.size());
        ns.reserve(fits.size());
        for (auto &f : fits) {
            infoMsg("  Fitted '%s': %d subjects after NaN removal", f.name.c_str(),
                    f.nUsedRows);
            rs.push_back(std::move(f.residuals));
            ns.push_back(f.name);
        }
        if (saveResid) {
            std::vector<nullmodel::NullModelFit> dumpFits(rs.size());
            for (size_t i = 0; i < rs.size(); ++i) {
                dumpFits[i].name = ns[i];
                dumpFits[i].residuals = rs[i];
                dumpFits[i].nUsedRows = static_cast<int>(rs[i].size());
            }
            nullmodel::writeResidualsFile(outPrefix + ".null.resid", sd, dumpFits,
                                          compression, compressionLevel);
        }
        sd.setResidualsFromFit(std::move(rs), std::move(ns));
    }

    // Environment vectors in union-dense order (aligned to the residuals).
    std::vector<Eigen::VectorXd> envUnion;
    envUnion.reserve(envNames.size());
    for (const auto &en : envNames)
        envUnion.push_back(sd.getColumn(en));

    // (+) variant: load the union-level sparse GRM once; re-index per phenotype.
    std::shared_ptr<const SparseGRM> unionGrm;
    if (hasGrm) {
        auto subjIDs = sd.usedIIDs();
        unionGrm = std::make_shared<const SparseGRM>(
            SparseGRM::load(spgrmGrabFile, spgrmGctaFile, subjIDs, sd.famIIDs()));
        infoMsg("Sparse GRM: %zu entries (SPAGxE+ variance correction)",
                unionGrm->nnz());
    }

    auto genoData =
        makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);
    infoMsg("Genotype data: %u markers, %u subjects", genoData->nMarkers(),
            genoData->nSubjUsed());

    auto phenoInfos = sd.buildPerColumnMasks();
    const int K = sd.residOneCols();
    if (K > 1) infoMsg("Multi-column residual file: %d phenotypes", K);

    std::vector<PhenoTask> tasks(K);
    for (int rc = 0; rc < K; ++rc) {
        const auto &pi = phenoInfos[rc];
        Eigen::VectorXd phenoResid =
            (K > 1) ? extractPhenoVec(sd.residMatrix().col(rc), pi) : sd.residuals();
        std::vector<Eigen::VectorXd> envVecs;
        envVecs.reserve(envUnion.size());
        for (const auto &eu : envUnion)
            envVecs.push_back((K > 1) ? extractPhenoVec(eu, pi) : eu);

        // Per-phenotype GRM: for a single phenotype the union GRM is used
        // directly; for K > 1 the union entries are re-indexed to the
        // phenotype's dense subject order (as runSPAGRM does).
        std::shared_ptr<const SparseGRM> phenoGrm;
        if (hasGrm) {
            if (K == 1) {
                phenoGrm = unionGrm;
            } else {
                const auto &u2l = pi.unionToLocal;
                std::vector<SparseGRM::Entry> pEntries;
                for (const auto &en : unionGrm->entries()) {
                    const uint32_t li = u2l[en.row], lj = u2l[en.col];
                    if (li != UINT32_MAX && lj != UINT32_MAX)
                        pEntries.push_back({li, lj, en.value});
                }
                phenoGrm = std::make_shared<const SparseGRM>(
                    SparseGRM::fromEntries(pi.nUsed, std::move(pEntries)));
            }
        }

        // Per-phenotype Branch-B Wald data: only in fit mode on the base (non-
        // GRM) path — SPAGxE+ keeps no Wald, residual mode has no raw phenotype.
        // Extract the raw phenotype / covariate design into the same dense
        // subject order as phenoResid (identical (K>1)?extract:union pattern),
        // so the Wald design aligns with the engine's per-phenotype genotype.
        std::shared_ptr<const spagxe_wald::WaldData> phenoWald;
        if (fitPath && !hasGrm) {
            const spagxe_wald::TraitType tt = traitTypes[rc];
            if (tt != spagxe_wald::TraitType::None) {
                auto wd = std::make_shared<spagxe_wald::WaldData>();
                wd->trait = tt;
                wd->covar = (K > 1) ? extractPhenoMat(covarUnion, pi) : covarUnion;
                const auto &spec = phenoSpecs[rc];
                if (tt == spagxe_wald::TraitType::Cox) {
                    Eigen::VectorXd tU = sd.getColumn(spec.timeColumn);
                    Eigen::VectorXd eU = sd.getColumn(spec.eventColumn);
                    wd->time = (K > 1) ? extractPhenoVec(tU, pi) : tU;
                    wd->event = (K > 1) ? extractPhenoVec(eU, pi) : eU;
                } else {
                    Eigen::VectorXd yU = recodeResponse(sd.getColumn(spec.yColumn), tt);
                    wd->y = (K > 1) ? extractPhenoVec(yU, pi) : yU;
                }
                phenoWald = std::move(wd);
            }
        }

        // SPAGxEmix: per-phenotype [1|PCs] design + OLS matrices for computeAFVec.
        // Extracted into the phenotype's dense subject order (as covarUnion is),
        // so the AF cascade aligns with the residual / environment / genotype.
        std::shared_ptr<const AFData> phenoAF;
        if (mix) {
            auto af = std::make_shared<AFData>();
            af->nPC = nPC;
            af->onePlusPCs =
                (K > 1) ? extractPhenoMat(unionOnePlusPCs, pi) : unionOnePlusPCs;
            const Eigen::MatrixXd &X = af->onePlusPCs;
            const Eigen::MatrixXd XtX = X.transpose() * X;
            const Eigen::MatrixXd XtX_inv =
                XtX.ldlt().solve(Eigen::MatrixXd::Identity(1 + nPC, 1 + nPC));
            af->XtX_inv_Xt = XtX_inv * X.transpose();
            af->sqrt_XtX_inv_diag = XtX_inv.diagonal().cwiseSqrt();
            phenoAF = std::move(af);
        }

        tasks[rc].phenoName = pi.name;
        tasks[rc].method = std::make_unique<SPAGxEMethod>(
            std::move(phenoResid), envNames, std::move(envVecs), phenoGrm,
            std::move(phenoWald), marginalCutoff, spaCutoff, outlierIqrRatio,
            std::move(phenoAF));
        tasks[rc].unionToLocal = pi.unionToLocal;
        tasks[rc].nUsed = pi.nUsed;
        infoMsg("  Phenotype '%s': %u subjects", pi.name.c_str(), pi.nUsed);
    }

    infoMsg("Running %s marker tests (%d thread(s), %d phenotype(s))...",
            dispLabel.c_str(), nthreads, K);
    multiPhenoEngine(*genoData, tasks, outPrefix, fileSuffix, compression,
                     compressionLevel, nthreads, missingCutoff, minMafCutoff,
                     minMacCutoff, hweCutoff);
}

} // namespace

// ══════════════════════════════════════════════════════════════════════
// runSPAGxE / runSPAGxEmix — thin public entry points over runSPAGxEImpl
// ══════════════════════════════════════════════════════════════════════

void runSPAGxE(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double marginalCutoff,
    double spaCutoff,
    double outlierIqrRatio,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    bool saveResid
) {
    runSPAGxEImpl(phenoFile, residNames, envNames, /*pcColNames=*/{}, spgrmGrabFile,
                  spgrmGctaFile, geno, outPrefix, compression, compressionLevel,
                  marginalCutoff, spaCutoff, outlierIqrRatio, nthreads, nSnpPerChunk,
                  missingCutoff, minMafCutoff, minMacCutoff, hweCutoff, keepFile,
                  removeFile, regressionModelStr, phenoNameSpec, covarFile, covarNames,
                  saveResid);
}

void runSPAGxEmix(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::vector<std::string> &pcColNames,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double marginalCutoff,
    double spaCutoff,
    double outlierIqrRatio,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    bool saveResid
) {
    // SPAGxEmix has no sparse-GRM path (SPAGxEmix+ is out of scope).
    runSPAGxEImpl(phenoFile, residNames, envNames, pcColNames, /*spgrmGrabFile=*/{},
                  /*spgrmGctaFile=*/{}, geno, outPrefix, compression, compressionLevel,
                  marginalCutoff, spaCutoff, outlierIqrRatio, nthreads, nSnpPerChunk,
                  missingCutoff, minMafCutoff, minMacCutoff, hweCutoff, keepFile,
                  removeFile, regressionModelStr, phenoNameSpec, covarFile, covarNames,
                  saveResid);
}
