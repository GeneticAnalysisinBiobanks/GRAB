// spasqr.cpp — SPAsqr: SPA-squared multi-tau marker association (pure C++17 / Eigen / Boost)

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "spasqr/spasqr.hpp"
#include "engine/loco.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "spasqr/long_gee.hpp"
#include "spasqr/null_model.hpp"
#include "spasqr/qmme.hpp"
#include "spagrm/longitudinal_resid.hpp"
#include "spagrm/sageld_fit.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/spa.hpp"
#include "util/spa_cgf.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

// No <immintrin.h>, util/simd_dispatch.hpp or util/simd_math.hpp any more: the
// five hand-written CGF copies that needed them are gone, and the scalar / AVX2
// / AVX-512 triple now lives once in util/spa_cgf.cpp with its own dispatch
// site.  See the block comment below.
#include "util/empirical_cgf.hpp"

// ══════════════════════════════════════════════════════════════════════
// SPAsqr's saddlepoint — on the shared solver, tail kernel and CGF
// ══════════════════════════════════════════════════════════════════════
//
// spa_unify Stage 4.  What used to live here was, in order:
//
//   * `scalarOutlierCgf` plus `simdOutlierCgf_avx512` and
//     `simdOutlierCgf_avx2` with their masked / scalar tails and a
//     `pickOutlierCgfFn` dispatch — FIVE copies of the same binomial CGF, each
//     forming K'' as the cancelling difference MGF2/MGF0 - (MGF1/MGF0)^2 with
//     the two sums accumulated GLOBALLY and differenced once at the end
//     (01_findings.md D1, measured 36 % relative error at t*r = 35);
//   * `scalarFastGetRoot`, a private Family-B Newton iteration with no
//     convergence flag (D5) and no check on the K'' it divided by;
//   * `scalarGetProbSpa` and, inlined twice inside `fusedGetPvalSpa`, THREE
//     copies of the Barndorff-Nielsen modified signed root, all three with no
//     isfinite(zeta) check, no zeta*s - K >= 0 check, no K'' > 0 check, no
//     w != 0 check and no v/w > 0 check (D6);
//   * `scalarGetPvalFromScore`, a near-duplicate of `fusedGetPvalSpa`.
//
// All of it is deleted.  The CGF is spa_cgf::binomUniform — SPAsqr's outlier
// block IS the shared uniform-MAF binomial CGF with h = 2 — which supplies
// D1's exact cancellation-free K'' = 2 r^2 e u / alpha^2, N1's
// K = 2 log1p(p expm1(t r)), and the same scalar + AVX2 + AVX-512 triple this
// file used to carry by hand.  The root finder is spa::solveSaddlepoint and the
// tail is spa::bnTailLog.  Five CGF copies plus three tail copies
// plus two p-value copies collapse to ONE p-value routine, `tauPvalue`, which
// all three MethodBase entry points call.
//
// ── On deleting rather than migrating the "scalar tier" ───────────────
//
// `scalarGetProbSpa`, `scalarGetPvalFromScore` and the D6/D7 copies they
// carried were UNREACHABLE in any GRAB run: SPAsqrMethod::supportsFusedGemm()
// returns true, both entry points route through multiPhenoEngine or locoEngine
// (which delegates to multiPhenoEngineRange), and in marker.cpp a fuseable
// phenotype is dispatched exclusively through processScoreBatch while
// getResultBatch is called only for the nonFusedPhenos set.  A repair confined
// to those copies would therefore have changed no output and would have been
// undetectable by examples/baseline.sh.  Rather than migrate dead code, the
// duplicates are removed and the three MethodBase overrides — which the
// contract requires SPAsqr to implement whether or not the engine calls them —
// are pointed at the single live implementation.  CLAUDE.md's alpha-release
// policy is explicit that removing a dead surface beats maintaining it.
//
// ── D7 ───────────────────────────────────────────────────────────────
//
// The upper tail was called with a bare literal 1e-4 and the lower tail with
// the configured `spa.tol` (1e-6), so the two probabilities being summed were
// not computed to the same accuracy and --spasqr-tol did not do what its name
// implies.  Both tails now use `spa.tol`, and the criterion is relative
// (|K'(t) - s| <= tol*sqrt(K''(t))) rather than an absolute test on the pending
// step.  This is a real numeric change and is the point of the stage.
//
// ── D6, and why SPAsqr is exposed to it exactly as SPAGRM is ─────────
//
// The first-pass audit note claimed "EmpVar is non-negative by construction, so
// SPAsqr is not exposed".  At the time EmpVar was built from a sparse-GRM
// quadratic form restricted to the non-outlier block, accumulated as
// `e.factor * e.value * ResidMat(e.row, col) * ResidMat(e.col, col)` over a
// thresholded GRM read from a file, and such a matrix is not guaranteed
// positive semidefinite.  The CGF's non-outlier block is now the plain inner
// product R_non^T R_non (see SPAsqrPerTau), which is non-negative, so that
// route is closed; the GRM still enters Score_var = G_var * R^T Phi R, and a
// non-positive Score_var surfaces as a named spa::Status with P = NA.

using namespace spasqr_null;

namespace {

// ── Per-tau SPA data for SPAsqr (shared across clones via shared_ptr) ──

struct SPAsqrPerTau {
    Eigen::VectorXd outlierResid;        // residuals of outlier subjects
    double sum_unrelated_outliers2;      // outlierResid.squaredNorm()
    double sum_R_nonOutlier;
    // R_non^T R_non — plain residual inner product over non-outliers, per the
    // manuscript's Var(S_non) = 2mu(1-mu)*R_non^T R_non.  The sparse GRM
    // deliberately does NOT enter the CGF: the whole CGF is built under
    // independence and all relatedness is carried by the variance-ratio
    // denominator Var(S) = G_var * R_GRM_R.
    double R_R_nonOutlier;
    double R_GRM_R;                      // full variance term (all subjects), R^T Phi R

    // Empirical CGF of this (chromosome, tau)'s residuals, used instead of the
    // Binomial MGF below the MAC cutoff.  Empty when the branch is disabled.
    ecgf::EmpiricalCgf empCgf;
};

// `zeta` is gone: it was assigned 0.01 in buildSPAsqrMethod and never read,
// because both tails derived their initial abscissa from
// min(|Score_adj|/Score_var, 1.2) instead.  See 01_findings.md D7's companion
// note and the lower-tail convention comment in spagrm.cpp.

// Below this minor-allele count the Binomial-MGF CGF is replaced by the
// empirical residual CGF.  The Binomial branch approximates the whole
// non-outlier block by a Gaussian, which is what makes rare variants
// conservative; the empirical CGF is exact in the genotype configuration and
// only assumes the residuals are exchangeable.  Measured on i.i.d. synthetic
// null markers (50k subjects, 9 taus, obs/exp at p < 1e-3): the Binomial
// branch sits at 0.57 for MAC 10-15, 0.8-0.9 for 20-50, 0.93 for 50-100 and
// 0.86-0.9 for 100-140, and is calibrated from MAC ~150 on; the empirical
// branch is calibrated at every MAC tried (10-2000).  The Binomial deficit
// does not end at 200: at MAC 200-400 (5000 unrelated subjects, one record
// each, 200k i.i.d. markers x 200 null phenotypes -- normal, t3, chi2(2) --
// 9 taus) the Binomial branch was 0.93 at 1e-4 and 0.91 at 1e-5, the
// empirical branch 1.00 and 1.04, consistently in two independent batches of
// phenotypes.  Hence 400.  A heavy-tailed residual column (the linear-GEE
// column of LoQus under t3 errors) is the exception -- the
// empirical CGF is conservative there -- and was judged not to matter.
constexpr double kEmpiricalMacCutoff = 400.0;

struct SPAsqrSPAShared {
    std::vector<SPAsqrPerTau> perTau;    // one per tau
    double SPA_Cutoff;
    double tol;
};

// ── The two-sided saddlepoint, once ────────────────────────────────────
//
// Both tails share the MAF / variance setup; only the sign of s, the initial
// abscissa and the tail direction differ.  The CGF is the outlier block through
// spa_cgf::binomUniform plus SPAsqr's analytic non-outlier Gaussian block,
// whose mean and variance enter K as mean*t + var*t^2/2 — exactly as before,
// and in the same association order.

inline spa::Result spaTwoSided(
    const double *oResid,
    int nOutlier,
    double MAF,
    double mean,
    double var,
    double absAdj,
    double initZeta,
    double tol,
    double zNorm
) {
    double logp[2];
    spa::Status st[2];

    for (int k = 0; k < 2; ++k) {
        const bool lowerTail = (k == 1);
        const double s = lowerTail ? -absAdj : absAdj;

        spa::SolveOpts opt;
        opt.init = lowerTail ? -initZeta : initZeta;
        opt.scoreSign = s;   // only the sign is read
        opt.rtol = tol;      // D7: the configured tolerance reaches both tails

        const spa::Saddle sd = spa::solveSaddlepoint(
            s,
            [&](double t) {
                const spa_cgf::Cgf12 d =
                    spa_cgf::binomUniformK12(t, oResid, nOutlier, MAF);
                return spa::K12{d.K1 + mean + var * t - s, d.K2 + var};
            },
            [&](double t) {
                const spa_cgf::Cgf012 d =
                    spa_cgf::binomUniformKFull(t, oResid, nOutlier, MAF);
                return spa::K012{d.K0 + mean * t + 0.5 * var * t * t,
                                 d.K1 + mean + var * t - s,
                                 d.K2 + var};
            },
            opt);

        spa::Status stLog = spa::Status::SpaOk;
        logp[k] = spa::bnTailLog(sd.zeta, s, sd.K0, sd.K2, lowerTail, stLog);
        st[k] = spa::worseStatus(sd.status, stLog);
    }

    // The two-sided assembly and the D5 fallback when either tail failed —
    // one call, in spa.hpp.  One tail evaluation per side: the linear sibling
    // of `bnTailLog` is gone with the P column's parallel assembly.
    // `zNorm` is the raw score z, not absAdj/sqrt(Var): the statistic has
    // already been rescaled to the CGF's natural variance by this point.
    return spa::assemble(logp[0], logp[1], st[0], st[1], zNorm);
}

// ── Binomial-MGF branch ──────────────────────────────────────────────
//
// G_var = empirical Var(G) over post-impute genotypes, passed in by the caller,
// drives the outer Score_var = G_var * R_GRM_R.  EmpVar = K''(0) uses the
// model-based 2p(1-p) for both the non-outlier partial-Gaussian piece and the
// outlier Bernoulli factor; the rescaling
// Score_adj = Score * sqrt(K''(0) / Var(S)_observed) is the SAIGE-style
// variance ratio that brings the score onto the CGF's natural scale, dropping
// the implicit HWE assumption (G_var ~ 2p(1-p)) earlier code relied on.
// EmpVar is built entirely under independence; the rescaling carries all of
// the relatedness, matching the manuscript (Phi enters only via Var(S)).

inline spa::Result binomTwoSided(
    double Score,
    double MAF,
    double Score_var,
    double zNorm,
    const SPAsqrPerTau &tau,
    const SPAsqrSPAShared &spa_
) {
    const double EmpVar    = 2.0 * MAF * (1.0 - MAF) *
                              (tau.R_R_nonOutlier + tau.sum_unrelated_outliers2);
    const double Score_adj = Score * std::sqrt(EmpVar / Score_var);
    const double absAdj    = std::abs(Score_adj);

    // First-order saddlepoint estimate zeta ~ s/K''(0), capped; the lower tail
    // starts at its negation.
    const double initZeta = std::min(absAdj / Score_var, 1.2);

    return spaTwoSided(tau.outlierResid.data(),
                       static_cast<int>(tau.outlierResid.size()),
                       MAF,
                       2.0 * MAF * tau.sum_R_nonOutlier,
                       2.0 * MAF * (1.0 - MAF) * tau.R_R_nonOutlier,
                       absAdj, initZeta, spa_.tol, zNorm);
}

// ══════════════════════════════════════════════════════════════════════
// Empirical-CGF branch (MAC < kEmpiricalMacCutoff)
// ══════════════════════════════════════════════════════════════════════
//
// Same solver, same tail kernel, same two-sided assembly as the Binomial
// branch; all that changes is the source of K, K' and K'' — from the outlier
// Binomial loop to the carrier sum of util/empirical_cgf.
// The branch reports through the same spa::Status vocabulary, so a marker
// that could not be solved here is visible as such rather than laundered.

inline spa::Result empTwoSided(
    const ecgf::EmpiricalCgf &tab,
    const ecgf::CarrierSet &cs,
    double absAdj,
    double initZeta,
    double tol,
    double zNorm
) {
    double logp[2];
    spa::Status st[2];

    for (int k = 0; k < 2; ++k) {
        const bool lowerTail = (k == 1);
        const double s = lowerTail ? -absAdj : absAdj;

        spa::SolveOpts opt;
        opt.init = lowerTail ? -initZeta : initZeta;
        opt.scoreSign = s;
        opt.rtol = tol;

        const spa::Saddle sd = spa::solveSaddlepoint(
            s,
            [&](double t) {
                double k1, k2;
                // empK1K2 subtracts the target itself: k1 = K'(t) - s.
                ecgf::empK1K2(tab, cs, t, s, k1, k2);
                return spa::K12{k1, k2};
            },
            [&](double t) {
                double k1, k2;
                ecgf::empK1K2(tab, cs, t, s, k1, k2);
                return spa::K012{ecgf::empK0(tab, cs, t), k1, k2};
            },
            opt);

        spa::Status stLog = spa::Status::SpaOk;
        logp[k] = spa::bnTailLog(sd.zeta, s, sd.K0, sd.K2, lowerTail, stLog);
        st[k] = spa::worseStatus(sd.status, stLog);
    }

    return spa::assemble(logp[0], logp[1], st[0], st[1], zNorm);
}

// Two-sided p-value from the empirical CGF.  Returns Status::NaNoTest when the
// branch has no answer to give (degenerate variance, score outside the support
// of the residual-randomised statistic); the caller then takes the Binomial
// branch rather than dropping the marker.
inline spa::Result empiricalTwoSided(
    const ecgf::EmpiricalCgf &tab,
    const ecgf::CarrierSet &cs,
    double Score,
    double Score_var,
    double tol,
    double zNorm
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // K''(0) = Var(R) * sum_i Gtilde_i^2.  Score_adj puts the score on the CGF's own scale
    // — the same variance-ratio step the Binomial branch performs with
    // 2p(1-p)*||R||^2, so the two branches meet continuously at the cutoff.
    const double K2at0 = ecgf::empK2AtZero(tab, cs);
    if (!(K2at0 > 0.0) || !(Score_var > 0.0))
        return spa::Result{nan, spa::Status::NaNoTest};
    const double absAdj = std::abs(Score) * std::sqrt(K2at0 / Score_var);
    if (!std::isfinite(absAdj) || absAdj <= 0.0)
        return spa::Result{nan, spa::Status::NaNoTest};

    // K'(t) is bounded — the residual-randomised score has bounded support, so
    // outside it the saddlepoint equation has no solution at all.  Reject here
    // instead of letting the solver exhaust its budget on an unsolvable root.
    double supLo, supHi;
    ecgf::empSupport(tab, cs, supLo, supHi);
    if (!(absAdj < supHi) || !(-absAdj > supLo))
        return spa::Result{nan, spa::Status::NaNoTest};

    const double zeta0 = absAdj / K2at0;   // exact root of a quadratic CGF
    return empTwoSided(tab, cs, absAdj, zeta0, tol, zNorm);
}

// ══════════════════════════════════════════════════════════════════════
// SPAsqrMethod — MethodBase adapter with zero per-clone heap allocation
//
// All per-tau SPA data is shared via shared_ptr (read-only).  Cloning copies two
// shared_ptr's; the only per-clone storage is the four B x ntaus result buffers,
// which processScoreBatch resizes once and reuses.
// ══════════════════════════════════════════════════════════════════════

class SPAsqrMethod : public MethodBase {
  public:
    SPAsqrMethod(
        int ntaus,
        std::shared_ptr<const SPAsqrSPAShared> spaShared,
        Eigen::MatrixXd residMat,
        Eigen::VectorXd residSums,
        std::vector<std::string> tauLabels,
        int nCct
    )
        : m_ntaus(ntaus),
          m_nCct(nCct),
          m_spaShared(std::move(spaShared)),
          m_tauLabels(std::move(tauLabels))
    {
        auto sd = std::make_shared<SharedMethodData>();
        sd->residMat  = std::move(residMat);
        sd->residSums = std::move(residSums);
        m_methodShared = std::move(sd);
        // Only ask the engine for genotypes if a table was actually built.
        m_hasEmpCgf = !m_spaShared->perTau.empty() &&
                      !m_spaShared->perTau[0].empCgf.empty();
    }

    std::unique_ptr<MethodBase> clone() const override {
        // Both shared_ptr's are read-only — no per-clone scratch at all.
        return std::make_unique<SPAsqrMethod>(*this);
    }

    int resultSize() const override {
        return 4 * m_ntaus + (m_nCct > 0 ? 1 : 0);
    }

// LOG10P_CCT, then four per-tau groups: LOG10P, Z, Z_Norm, SPA_STATUS.
//
// The saddlepoint is per tau, so the status column is per tau as well: a marker
// can converge at one quantile level and fail at another, and a single
// aggregated status would hide which.  The linear P group left in log10p_unify
// Stage 8 (decision D1); LOG10P_tau is the sole p-value per quantile level and
// a consumer that needs the linear p takes 10^(-LOG10P_tau).
//
// The combination column is LOG10P_CCT, not P_CCT (log10p_unify Stage 5).  The
// Cauchy statistic's terms are cot(pi*p) -> 1/(pi*p), so the statistic itself
// overflows once any tau's p reaches 1e-308; the log-domain routine carries it
// as T = 10^M * A and never forms it (math::cauchyCombineLog10).  The
// combination is therefore taken over the LOG10P group rather than the P group,
// which is also why the P group could leave in Stage 8 without taking the
// combination with it.
//
// SPA_STATUS is static_cast<uint8_t>(spa::Status): 0 SPA_OK, 1 NORMAL,
// 2 SPA_W_SINGULAR, 3..6 the FALLBACK_* codes, 7 NA_POST_FAIL, 8 NA_NO_TEST.
// The order is a contract (log10p_unify D4): <= 2 means LOG10P is
// trustworthy, 3..6 that it is the substituted normal tail, >= 7 that it is
// NA.  Integer rather than the
// spa::statusName token because MethodBase hands the engine a
// std::vector<double> and marker_impl.hpp formats every cell through
// numToChars; a token column would need a new hook in the MethodBase contract,
// which dev-notes/methods/spa_unify/02_design.md places out of scope for the
// per-method migration stages.  This is the encoding Stage 3 established for
// SPACox, kept identical so all methods agree.
    std::string getHeaderColumns() const override {
        std::ostringstream oss;
        if (m_nCct > 0) oss << "\tLOG10P_CCT";
        static const char *const kGroups[] = {"LOG10P_", "Z_", "Z_Norm_",
                                              "SPA_STATUS_"};
        for (const char *g : kGroups) {
            for (int i = 0; i < m_ntaus; ++i) {
                oss << '\t' << g << m_tauLabels[i];
            }
        }
        return oss.str();
    }

    // SPAsqr is unconditionally fuseable (see supportsFusedGemm below), and the
    // engine routes purely on that flag, so every marker arrives through
    // processScoreBatch and this per-marker entry point is unreachable.  It
    // exists only because MethodBase declares it pure virtual.
    //
    // It is a throw rather than a second implementation on purpose.  A working
    // copy here would have to duplicate the score centering, the empirical
    // Var(G) and — since the low-MAC branch needs per-subject genotypes that
    // this signature does happen to carry — the whole empirical-CGF path.  The
    // version that used to live here did the first two and skipped the third,
    // so it silently disagreed with the fused path below MAC 100, and no test
    // could catch that because nothing calls it.  If SPAsqr ever needs a
    // non-fused path, implement it against processScoreBatch's logic rather
    // than reviving a parallel one.
    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> /*GVec*/,
        double /*altFreq*/,
        int /*markerInChunkIdx*/,
        std::vector<double> & /*result*/
    ) override {
        throw std::logic_error(
            "SPAsqr::getResultVec: SPAsqr only runs through the fused-GEMM path "
            "(processScoreBatch); reaching here means the engine's fuseable/"
            "non-fuseable split changed.");
    }

    int preferredBatchSize() const override {
        return std::min(std::max(4, 2 * m_ntaus), 16);
    }

    // ── Fused union-level GEMM interface ───────────────────────────────
    bool supportsFusedGemm() const override {
        return true;
    }

    int fusedGemmColumns() const override {
        return m_ntaus;
    }

    void fillUnionResiduals(
        Eigen::Ref<Eigen::MatrixXd> dest,
        const std::vector<uint32_t> &unionToLocal
    ) const override {
        // The low-MAC branch needs to pull this phenotype's genotypes out of
        // the union-space matrix later; the mapping is only handed to us here.
        m_unionToLocal = unionToLocal;
        const auto &residMat = m_methodShared->residMat;
        const uint32_t nUnion = static_cast<uint32_t>(unionToLocal.size());
        for (uint32_t i = 0; i < nUnion; ++i) {
            const uint32_t li = unionToLocal[i];
            if (li != UINT32_MAX)
                dest.row(i) = residMat.row(li);
        }
    }

    void fillResidualSums(double *dest) const override {
        const auto &rs = m_methodShared->residSums;
        std::copy(rs.data(), rs.data() + m_ntaus, dest);
    }

    void processScoreBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &scores,
        const double *gSums,
        const double *gSumSqs,
        uint32_t nUsed,
        const std::vector<double> &altFreqs,
        const std::vector<int> & /*chunkIdxs*/,
        const UnionGenotypes &geno,
        std::vector<std::vector<double> > &results
    ) override {
        const int B = static_cast<int>(scores.cols());
        results.resize(B);

        // ── Center scores ──────────────────────────────────────────
        m_centeredBuf = scores;
        const double nUsedD = static_cast<double>(nUsed);
        const double invN   = 1.0 / nUsedD;
        for (int b = 0; b < B; ++b) {
            const double gMean = gSums[b] * invN;
            for (int t = 0; t < m_ntaus; ++t)
                m_centeredBuf(t, b) -= gMean * m_methodShared->residSums[t];
        }

        // Layout for the per-(marker, tau) buffers: [b * ntaus + t].  Heap-backed
        // per-thread scratch (was a [20*16] stack array) so the fused GEMM
        // window B is not capped at 16 / B·ntaus ≤ 320.
        const size_t nCell = static_cast<size_t>(B) * m_ntaus;
        m_zBuf.resize(nCell);
        m_lpBuf.resize(nCell);
        m_stBuf.resize(nCell);
        m_slow.assign(nCell, 0);

        // ── Pass 1: z, the normal branch, and who needs the saddlepoint ────
        // Everything here is O(1) per cell.  Var(G) depends only on the marker,
        // so it is computed once rather than once per tau.
        m_gVar.resize(static_cast<size_t>(B));
        m_needCarriers.assign(static_cast<size_t>(B), 0);
        const double cutoff = m_spaShared->SPA_Cutoff;
        const double nan = std::numeric_limits<double>::quiet_NaN();

        auto put = [&](size_t c, double z, const spa::Result &r) {
            m_zBuf[c]  = z;
            m_lpBuf[c] = r.negLog10p;
            m_stBuf[c] = static_cast<double>(static_cast<uint8_t>(r.status));
        };

        for (int b = 0; b < B; ++b) {
            // Empirical Var(G) over post-impute genotypes:
            //   Var(g) = (Σg² − (Σg)²/n) / (n−1)
            const double gSum  = gSums[b];
            const double G_var = (nUsedD > 1.0)
                ? (gSumSqs[b] - gSum * gSum / nUsedD) / (nUsedD - 1.0)
                : 0.0;
            m_gVar[static_cast<size_t>(b)] = G_var;

            const double MAF = std::min(altFreqs[b], 1.0 - altFreqs[b]);
            for (int t = 0; t < m_ntaus; ++t) {
                const double Svar = G_var * m_spaShared->perTau[t].R_GRM_R;
                const size_t c = static_cast<size_t>(b) * m_ntaus + t;

                if (Svar <= 0.0 || MAF <= 0.0) {
                    put(c, 0.0, spa::Result{nan, spa::Status::NaNoTest});
                    continue;
                }
                const double z = m_centeredBuf(t, b) / std::sqrt(Svar);
                if (!std::isfinite(z)) {
                    put(c, z, spa::Result{nan, spa::Status::NaNoTest});
                } else if (std::abs(z) <= cutoff) {
                    // Normal approximation.  The cutoff is not only a speed knob
                    // — Lugannani-Rice is singular at the mean (w, v -> 0 as the
                    // saddlepoint approaches 0), so every marker below it must
                    // stay here.  Measured: routing |z| <= 2 through the SPA
                    // drives the empirical branch to obs/exp 2.5 at 1e-3 and 16
                    // at 1e-4.
                    put(c, z, spa::normalBranch(z));
                } else {
                    m_zBuf[c] = z;
                    m_slow[c] = 1;
                    m_needCarriers[static_cast<size_t>(b)] = 1;
                }
            }
        }

        // ── Pass 2: gather carriers, but only where they will be read ──────
        gatherCarriers(B, gSums, invN, geno);

        // ── Pass 3: saddlepoint, tau-first ─────────────────────────────────
        // Only the cells flagged above reach here — a few percent under the
        // null.  Tau-first so one tau's outlierResid and CGF table stay hot
        // across all the markers that need them.
        for (int t = 0; t < m_ntaus; ++t) {
            const SPAsqrPerTau &tau = m_spaShared->perTau[t];
            for (int b = 0; b < B; ++b) {
                const size_t c = static_cast<size_t>(b) * m_ntaus + t;
                if (!m_slow[c]) continue;

                const double Score = m_centeredBuf(t, b);
                const double G_var = m_gVar[static_cast<size_t>(b)];
                const double Svar  = G_var * tau.R_GRM_R;
                const double MAF   = std::min(altFreqs[b], 1.0 - altFreqs[b]);
                const double z     = m_zBuf[c];

                // Below the MAC cutoff the Binomial MGF is replaced by the
                // empirical residual CGF; if that branch cannot deliver a
                // usable saddlepoint (SPA_STATUS > 2) fall back to the
                // Binomial branch rather than substitute or lose the marker.
                spa::Result r{nan, spa::Status::NaNoTest};
                const ecgf::CarrierSet &cs = m_carriers[static_cast<size_t>(b)];
                if (!cs.empty() && !tau.empCgf.empty())
                    r = empiricalTwoSided(tau.empCgf, cs, Score, Svar,
                                          m_spaShared->tol, z);
                if (!spa::statusIsUsable(r.status))
                    r = binomTwoSided(Score, MAF, Svar, z, tau, *m_spaShared);
                put(c, z, r);
            }
        }

        for (int b = 0; b < B; ++b) {
            const size_t off = static_cast<size_t>(b) * m_ntaus;
            assemble(m_lpBuf.data() + off, m_zBuf.data() + off,
                     m_stBuf.data() + off, results[b]);
        }
    }

  private:
    // Gather the carrier sets of the markers pass 1 flagged, and only those.
    //
    // Reading a genotype column costs a pass over all nUnion subjects, so doing
    // it for every marker under the MAC cutoff would spend that pass on markers
    // no tau will ever ask about; under the null only a few percent of
    // (marker, tau) cells clear the SPA cutoff.
    void gatherCarriers(
        int B,
        const double *gSums,
        double invN,
        const UnionGenotypes &geno
    ) {
        m_carriers.resize(static_cast<size_t>(B));
        for (auto &cs : m_carriers) { cs.gt.clear(); cs.n0 = 0.0; }
        if (!m_hasEmpCgf || geno.empty()) return;

        const uint32_t nUnion = static_cast<uint32_t>(m_unionToLocal.size());
        if (static_cast<Eigen::Index>(nUnion) != geno.rows) return;

        const int nMarkers = std::min(B, static_cast<int>(geno.nCols));
        for (int b = 0; b < nMarkers; ++b) {
            if (!m_needCarriers[static_cast<size_t>(b)]) continue;
            if (!(geno.macs[b] < kEmpiricalMacCutoff)) continue;
            m_carriers[static_cast<size_t>(b)].gather(
                geno.column(static_cast<uint32_t>(b)),
                m_unionToLocal.data(), nUnion, gSums[b] * invN);
        }
    }

    int m_ntaus;
    // LOG10P_CCT combines the first m_nCct columns only (0 ⇒ no CCT column).
    // Plain SPAsqr combines all of them; LoQus appends a
    // linear-GEE column that is reported beside the quantile levels but is not
    // a quantile level, so it stays out of the combination.
    int m_nCct;
    std::shared_ptr<const SPAsqrSPAShared> m_spaShared;
    std::vector<std::string> m_tauLabels;

    struct SharedMethodData {
        Eigen::MatrixXd residMat;   // N × ntaus
        Eigen::VectorXd residSums;  // ntaus
    };

    std::shared_ptr<const SharedMethodData> m_methodShared;
    // Low-MAC branch state.  m_unionToLocal is filled by fillUnionResiduals
    // (const, hence mutable); the genotype block itself arrives as a
    // processScoreBatch argument and is not retained.
    bool m_hasEmpCgf = false;
    mutable std::vector<uint32_t> m_unionToLocal;
    std::vector<ecgf::CarrierSet> m_carriers;   // B, filled by gatherCarriers


    // Per-thread scratch, all reused across processScoreBatch calls.
    Eigen::MatrixXd m_centeredBuf;
    std::vector<double>  m_zBuf, m_lpBuf, m_stBuf;  // B × ntaus: Z_Norm, LOG10P, SPA_STATUS
    std::vector<uint8_t> m_slow;          // B × ntaus: cell needs the saddlepoint
    std::vector<double>  m_gVar;          // B
    std::vector<uint8_t> m_needCarriers;  // B: some tau of this marker is slow

    // Column assembly, shared by all three entry points so the header and the
    // row can only ever be built by the same code.
    void assemble(
        const double *lgs,
        const double *zScores,
        const double *stats,
        std::vector<double> &result
    ) const {
        result.clear();
        result.reserve(4 * m_ntaus + 1);

        // math::cauchyCombineLog10 skips NaN entries, so a tau whose
        // saddlepoint failed drops out of the combination rather than poisoning
        // it; that rule is the linear routine's, carried over unchanged.  The
        // inputs are the magnitudes, which is what makes the column meaningful
        // past L = 308: the tail is L_CCT ~ L_max - log10(ntaus), so the
        // combination is dominated by the SMALLEST p and inherits its
        // magnitude, where the linear statistic overflowed and returned
        // P_CCT = 0 (log10p_unify Stage 5, 01_numerics §2.4).
        if (m_nCct > 0) result.push_back(math::cauchyCombineLog10(lgs, m_nCct));

        for (int i = 0; i < m_ntaus; ++i) result.push_back(lgs[i]);
        for (int i = 0; i < m_ntaus; ++i) {
            const double zr = zScores[i];
            // D3: on the fast path (|z| ≤ SPA_Cutoff) the p-value is exactly
            // 2·Φ(−|z|), so the inversion returns the raw z; emit it directly
            // and skip the quantile round-trip (which also removes the spurious
            // round-trip ULP gap that made Z_τ differ from Z_Norm_τ on fast-path
            // entries).  A marker with no p-value (degenerate, or a saddlepoint
            // failure that D5 does not fall back from) has L = NaN and
            // propagates as NaN, exactly as the inversion itself would.
            // Stage 4: inverted from L rather than from the linear p, so Z_τ no
            // longer saturates at 37.0470962993612 once LOG10P_τ passes 299.7.
            if (std::isnan(lgs[i]))
                result.push_back(std::numeric_limits<double>::quiet_NaN());
            else if (std::abs(zr) <= m_spaShared->SPA_Cutoff)
                result.push_back(zr);                                 // Z_τ == Z_Norm_τ
            else
                result.push_back(math::zFromNegLog10P(lgs[i], zr));   // SPA-recalibrated
        }
        for (int i = 0; i < m_ntaus; ++i) result.push_back(zScores[i]);  // Z_Norm_τ
        for (int i = 0; i < m_ntaus; ++i) result.push_back(stats[i]);    // SPA_STATUS
    }

};

// ══════════════════════════════════════════════════════════════════════
// Outlier detection (IQR-based, per column)
// ══════════════════════════════════════════════════════════════════════

// Returns an N × K boolean matrix (as std::vector<std::vector<bool>>).
// outlierIqrRatio  = multiplier for IQR (default 1.5)
// outlierAbsBound  = absolute clamp for cutoffs (default 0.55)
struct OutlierInfo {
    // per-column: number of outlier subjects.  Only the count is ever wanted —
    // callers that need the values themselves walk isOutlier — so the indices
    // are not materialised.
    std::vector<uint32_t> outlierCount;
    // per-column: boolean mask (size N)
    std::vector<std::vector<bool> > isOutlier;
};

OutlierInfo detectOutliers(
    const Eigen::MatrixXd &ResidMat,
    double outlierIqrRatio,
    double outlierAbsBound
) {
    const Eigen::Index N = ResidMat.rows();
    const Eigen::Index K = ResidMat.cols();

    OutlierInfo info;
    info.outlierCount.assign(static_cast<size_t>(K), 0);
    info.isOutlier.resize(K);

    // Scratch for sorting
    std::vector<double> scratch(N);

    for (Eigen::Index col = 0; col < K; ++col) {
        // Copy column for sorting
        for (Eigen::Index i = 0; i < N; ++i)
            scratch[i] = ResidMat(i, col);

        std::sort(scratch.begin(), scratch.end());

        // Q1, Q3 via linear interpolation (same as R type=7)
        auto quantile = [&](double prob) -> double {
            const double idx = prob * (N - 1);
            const Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
            const Eigen::Index hi = std::min(lo + 1, N - 1);
            const double frac = idx - lo;
            return scratch[lo] * (1.0 - frac) + scratch[hi] * frac;
        };

        const double Q1 = quantile(0.25);
        const double Q3 = quantile(0.75);
        const double IQR = Q3 - Q1;

        double cutLo = Q1 - outlierIqrRatio * IQR;
        double cutHi = Q3 + outlierIqrRatio * IQR;
        cutLo = std::max(cutLo, -outlierAbsBound);
        cutHi = std::min(cutHi, outlierAbsBound);

        info.isOutlier[col].resize(N, false);
        for (Eigen::Index i = 0; i < N; ++i) {
            const double v = ResidMat(i, col);
            if (v < cutLo || v > cutHi) {
                info.isOutlier[col][i] = true;
                ++info.outlierCount[col];
            }
        }
    }
    return info;
}

// ══════════════════════════════════════════════════════════════════════
// Shared pipeline: outlier detection → GRM → method → marker engine
// ══════════════════════════════════════════════════════════════════════

// One sparse-GRM entry, pre-multiplied by its symmetry weight: the file stores
// only the lower triangle, so an off-diagonal pair contributes twice to
// R^T Phi R.  Folding that into the entry keeps the quadratic form a plain sum.
struct GRMEntry {
    uint32_t row, col;
    double value;
    double factor; // 1 for diagonal, 2 for off-diagonal
};

// Load GRM entries from disk and convert to flat GRMEntry vector.
std::vector<GRMEntry> loadGrmEntries(
    const std::vector<std::string> &subjOrder,
    const std::vector<std::string> &famIIDs,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile
) {
    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, subjOrder, famIIDs);
    infoMsg("Sparse GRM: %zu entries used", grm.nnz());
    std::vector<GRMEntry> entries;
    entries.reserve(grm.nnz());
    for (const auto &e : grm.entries())
        entries.push_back({e.row, e.col, e.value, (e.row == e.col) ? 1.0 : 2.0});
    return entries;
}

// Write one cross-tau correlation matrix as a small tab-separated text file:
// header = tau labels, then the ntaus x ntaus matrix.
static void writeOmegaFile(const std::string &path, const std::string &phenoName,
                           const std::vector<std::string> &tauLabels,
                           const Eigen::MatrixXd &omega)
{
    const int ntaus = static_cast<int>(tauLabels.size());
    std::ofstream ofs(path);
    if (!ofs)
        throw std::runtime_error("SPAsqr: cannot write " + path);
    for (int t = 0; t < ntaus; ++t)
        ofs << (t ? "\t" : "") << tauLabels[t];
    ofs << "\n" << std::setprecision(12);
    for (int a = 0; a < ntaus; ++a) {
        for (int b = 0; b < ntaus; ++b)
            ofs << (b ? "\t" : "") << omega(a, b);
        ofs << "\n";
    }
    infoMsg("[%s] Cross-tau residual correlation written to: %s",
            phenoName.c_str(), path.c_str());
}

// Build SPAsqrMethod from a pre-computed residual matrix and pre-loaded GRM entries.
std::unique_ptr<MethodBase> buildSPAsqrMethod(
    Eigen::MatrixXd &ResidMat,
    const std::vector<GRMEntry> &grmEntries,
    uint32_t nUsed,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    std::vector<std::string> tauLabels,
    std::vector<double> *outlierRatiosOut,
    Eigen::MatrixXd *omegaOut,
    int nCct
)
{
    const Eigen::Index N = ResidMat.rows();
    const Eigen::Index K = ResidMat.cols();
    const int ntaus = static_cast<int>(K);

    // ── 3. Outlier detection ───────────────────────────────────────────
    OutlierInfo outlierInfo = detectOutliers(ResidMat, outlierIqrRatio, outlierAbsBound);

    // Populate outlier ratios for caller to format
    if (outlierRatiosOut) {
        outlierRatiosOut->resize(K);
        for (Eigen::Index c = 0; c < K; ++c)
            (*outlierRatiosOut)[c] = static_cast<double>(outlierInfo.outlierCount[c]) / N;
    }

    // Cross-tau residual correlation, Omega_ab = R_a^T Phi R_b normalized to a
    // correlation: the quantile factor of the score-statistic covariance
    // Cov(Z) = LD (x) Omega, needed by multi-quantile summary-statistic
    // methods (ghost knockoffs, joint tests) alongside the LD matrix.  With
    // no GRM supplied, grmEntries is the identity and this is R^T R.  The
    // diagonal of the same quadratic form is R_GRM_R below.  Off-diagonal
    // GRM entries are stored once, so the a != b form is symmetrized.
    if (omegaOut) {
        Eigen::MatrixXd W = Eigen::MatrixXd::Zero(K, K);
        for (const auto &e : grmEntries) {
            for (Eigen::Index a = 0; a < K; ++a)
                for (Eigen::Index b = 0; b < K; ++b) {
                    if (e.row == e.col)
                        W(a, b) += e.value * ResidMat(e.row, a) * ResidMat(e.row, b);
                    else
                        W(a, b) += e.value * (ResidMat(e.row, a) * ResidMat(e.col, b) +
                                              ResidMat(e.col, a) * ResidMat(e.row, b));
                }
        }
        omegaOut->resize(K, K);
        for (Eigen::Index a = 0; a < K; ++a)
            for (Eigen::Index b = 0; b < K; ++b)
                (*omegaOut)(a, b) = W(a, b) / std::sqrt(W(a, a) * W(b, b));
    }

    // ── 4. Compute per-column variance terms + build SPAsqrSPAShared ──
    auto spaShared = std::make_shared<SPAsqrSPAShared>();
    spaShared->perTau.resize(ntaus);
    spaShared->SPA_Cutoff = spaCutoff;
    spaShared->tol        = 1e-6;

    for (int col = 0; col < ntaus; ++col) {
        const auto &isOut = outlierInfo.isOutlier[col];
        auto &pt = spaShared->perTau[col];

        // R_GRM_R = R^T Phi R over all subjects.  This is the ONLY place the
        // sparse GRM is used: it forms the normal-approximation variance
        // Var(S) = G_var * R^T Phi R, which is also the denominator of the
        // variance ratio.  The CGF itself never sees Phi.
        double rgrm_r = 0.0;
        for (const auto &e : grmEntries)
            rgrm_r += e.factor * e.value * ResidMat(e.row, col) * ResidMat(e.col, col);
        pt.R_GRM_R = rgrm_r;

        // sum_R_nonOutlier, R_non^T R_non, and outlier residual values
        double sumNO   = 0.0;
        double sumSqNO = 0.0;
        std::vector<double> outVals;
        outVals.reserve(outlierInfo.outlierCount[col]);
        for (uint32_t i = 0; i < nUsed; ++i) {
            const double r = ResidMat(i, col);
            if (!isOut[i]) {
                sumNO   += r;
                sumSqNO += r * r;
            } else {
                outVals.push_back(r);
            }
        }
        pt.sum_R_nonOutlier = sumNO;
        pt.R_R_nonOutlier   = sumSqNO;
        pt.outlierResid =
            Eigen::Map<Eigen::VectorXd>(outVals.data(), static_cast<Eigen::Index>(outVals.size()));
        pt.sum_unrelated_outliers2 = pt.outlierResid.squaredNorm();

        // Empirical residual CGF for the low-MAC branch.  buildSPAsqrMethod is
        // invoked once per chromosome by the LOCO driver, so this yields one
        // table per (chromosome, tau) — exactly what LOCO needs, the residuals
        // being chromosome-specific.
        pt.empCgf = ecgf::buildEmpiricalCgf(ResidMat.col(col).data(),
                                            static_cast<int>(nUsed));
    }

    // ── 5. Build residMat / residSums for fused dot products ──────────
    Eigen::VectorXd residSums(ntaus);
    for (int t = 0; t < ntaus; ++t)
        residSums[t] = ResidMat.col(t).sum();

    // ── 6. Build method ──────────────────────────────────────────────
    std::unique_ptr<SPAsqrMethod> method = std::make_unique<SPAsqrMethod>(
        ntaus,
        std::move(spaShared),
        Eigen::MatrixXd(ResidMat),   // copy — caller may reuse ResidMat
        std::move(residSums),
        std::move(tauLabels),
        nCct
    );

    return method;
}

// ══════════════════════════════════════════════════════════════════════
// Re-index GRM entries from union-dense space to pheno-dense space
// ══════════════════════════════════════════════════════════════════════

std::vector<GRMEntry> reindexGrm(
    const std::vector<GRMEntry> &unionGrm,
    const std::vector<uint32_t> &unionToLocal, // union index → pheno index (UINT32_MAX = absent)
    uint32_t nUnion
) {
    std::vector<GRMEntry> out;
    out.reserve(unionGrm.size());
    for (const auto &e : unionGrm) {
        if (e.row >= nUnion || e.col >= nUnion) continue;
        uint32_t lr = unionToLocal[e.row];
        uint32_t lc = unionToLocal[e.col];
        if (lr == UINT32_MAX || lc == UINT32_MAX) continue;
        out.push_back({lr, lc, e.value, e.factor});
    }
    return out;
}

// ══════════════════════════════════════════════════════════════════════
// Shared null-model setup for score mode and LOCO mode
// ══════════════════════════════════════════════════════════════════════

// The analysis set, resolved.  Both score entry points open by doing exactly
// this and then diverge: plain score mode fits the null model once here and
// now, LOCO mode defers it to the per-chromosome callback.
//
// SubjectData is held by pointer only because the caller keeps using it after
// the split (usedIIDs / usedMask / nFam all feed the genotype reader).
struct ScoreSetup {
    std::unique_ptr<SubjectData> sd;
    std::vector<PhenoWork> pw;   // one per phenotype, already transformed
    int nCov = 0;
    double hScale = 0.0;
};

// Union = subjects with genotype ∩ GRM ∩ keep/remove.  Per-phenotype NA
// filtering is deferred: each phenotype gets its own non-missing subset of the
// union, its own transformed Y and its own bandwidth.  `label` only prefixes
// the log lines and the error messages.
ScoreSetup buildScoreSetup(const SPAsqrConfig &cfg, const char *label) {
    ScoreSetup s;
    infoMsg("%s: pheno-transform = %s, solver = qmme", label, cfg.phenoTransform.c_str());
    infoMsg("%s: Loading phenotype and covariate data (%zu phenotypes, %zu taus)",
            label, cfg.phenoNames.size(), cfg.taus.size());

    s.sd = std::make_unique<SubjectData>(parseGenoIIDs(cfg.geno));
    SubjectData &sd = *s.sd;
    if (!cfg.phenoFile.empty()) sd.loadPhenoFile(cfg.phenoFile, cfg.phenoNames);
    if (!cfg.covarFile.empty()) sd.loadCovar(cfg.covarFile, cfg.covarNames);
    sd.setKeepRemove(cfg.keepFile, cfg.removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(cfg.spgrmGrabFile, cfg.spgrmGctaFile,
                                                 sd.famIIDs()));
    sd.setGenoLabel(cfg.geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(cfg.spgrmGrabFile, cfg.spgrmGctaFile));
    sd.finalize();

    const Eigen::Index N = static_cast<Eigen::Index>(sd.nUsed());
    const Eigen::MatrixXd unionX = cfg.covarNames.empty()
        ? (sd.hasCovar() ? Eigen::MatrixXd(sd.covar()) : Eigen::MatrixXd(N, 0))
        : sd.getColumns(cfg.covarNames);
    s.nCov = static_cast<int>(unionX.cols());

    s.pw = buildPhenoWorkspaces(sd, unionX, cfg.phenoNames, cfg.phenoTransform, label);
    s.hScale = (cfg.spasqrHScale >= 0.0) ? cfg.spasqrHScale : 3.0;
    for (auto &p : s.pw)
        p.h = (cfg.spasqrH >= 0.0) ? cfg.spasqrH
                                   : iqrBandwidth(p.Y, s.hScale, s.nCov);
    return s;
}

// Run `fit(k, t)` for every (phenotype, tau) pair across a work-stealing pool.
//
// The scaffolding is subtle enough to be worth having once: an exception must
// not escape a std::thread (that calls terminate), so each fit's failure is
// stashed by index, every thread is joined, and only then is the first failure
// rethrown with the phenotype and tau that produced it.  `label` prefixes both
// the progress line and that error.
template <typename Fit>
static void runQmmeFits(
    int K,
    int ntaus,
    int nthreads,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &colLabels,
    const std::string &label,
    Fit &&fit
) {
    const int totalFits = K * ntaus;
    const int nWorkers  = std::min(nthreads, totalFits);
    infoMsg("%s: Running %d null-model fits with %d threads", label.c_str(), totalFits, nWorkers);

    std::atomic<int> nextFit{0};
    std::vector<std::string> fitErrors(totalFits);

    auto worker = [&]() {
        for (;;) {
            const int idx = nextFit.fetch_add(1, std::memory_order_relaxed);
            if (idx >= totalFits) break;
            try {
                fit(idx / ntaus, idx % ntaus);
            } catch (const std::exception &ex) {
                fitErrors[idx] = ex.what();
            }
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(nWorkers - 1);
    for (int w = 0; w < nWorkers - 1; ++w)
        threads.emplace_back(worker);
    worker();
    for (auto &th : threads)
        th.join();

    for (int idx = 0; idx < totalFits; ++idx)
        if (!fitErrors[idx].empty())
            throw std::runtime_error(label + ": null-model fit failed for phenotype '" +
                                     phenoNames[idx / ntaus] + "' column " +
                                     colLabels[idx % ntaus] + ": " + fitErrors[idx]);
}

// ══════════════════════════════════════════════════════════════════════
// runSPAsqr — multi-phenotype entry point with parallel QMME fits
// ══════════════════════════════════════════════════════════════════════

} // anonymous namespace

void runSPAsqr(const SPAsqrConfig &cfg) {
    const auto &phenoNames = cfg.phenoNames;
    const auto &taus       = cfg.taus;
    const int K     = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());
    // Score mode fits the null model once and reuses it for every marker, so
    // --spasqr-tol (default 1e-6) is tight enough; apply it directly.
    const double qmmeTol = cfg.spasqrTol;

    // ── 1-3. Analysis set, covariates, per-phenotype split, bandwidth ──
    ScoreSetup setup = buildScoreSetup(cfg, "SPAsqr");
    SubjectData &sd = *setup.sd;
    std::vector<PhenoWork> &pw = setup.pw;
    const uint32_t nUnion = sd.nUsed();

    // Score mode fits one residual matrix per phenotype up front.
    std::vector<Eigen::MatrixXd> residMats(K);
    for (int k = 0; k < K; ++k)
        residMats[k].resize(static_cast<Eigen::Index>(pw[k].nk), ntaus);

    logPhenoTable(phenoNames, pw, "bandwidth");

    // ── 4. Parallel QMME fits: K × ntaus ────────────────────────────
    // Pre-construct one QMME solver per phenotype (X^T X / n cached) and
    // prepare Cholesky for that phenotype's bandwidth. Solver is then
    // shared read-only across worker threads — solve() only reads m_chol.
    std::vector<std::unique_ptr<qmme::SqrSolver> > qmmeSolvers(K);
    for (int k = 0; k < K; ++k) {
        qmmeSolvers[k] = std::make_unique<qmme::SqrSolver>(pw[k].X, /*delta*/ 1e-6);
        qmmeSolvers[k]->prepareBandwidth(pw[k].h);
    }

    runQmmeFits(K, ntaus, cfg.nthreads, phenoNames, makeTauLabels(taus), "SPAsqr", [&](int k, int t) {
        const double h1 = 1.0 / pw[k].h;
        Eigen::VectorXd resid;
        qmme::SolverStatus st;
        const Eigen::VectorXd beta = qmmeSolvers[k]->solve(
            pw[k].Y, taus[t], &resid, qmmeTol, /*maxIter*/ 50000,
            /*restartPeriod*/ 50, &st);
        infoMsg("[%s] tau=%.4f intercept=%.6f", phenoNames[k].c_str(), taus[t], beta(0));
        if (!st.converged)
            warnMsg("[%s] tau=%.4f qmme did not converge: %d iters (||g||=%.3e); tol=%.3e",
                    phenoNames[k].c_str(), taus[t], st.iter, st.finalGradNorm, qmmeTol);

        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        for (Eigen::Index i = 0; i < Nk; ++i)
            residMats[k](i, t) = taus[t] - math::pnorm(-resid(i) * h1);
    });

    // ── 5. Build tau labels ────────────────────────────────────────────
    const std::vector<std::string> tauLabels = makeTauLabels(taus);

    // ── 6. Load genotype data and GRM once (shared, union space) ────
    auto genoData = makeGenoData(cfg.geno, sd.usedMask(), sd.nFam(), sd.nUsed(), cfg.nSnpPerChunk);

    std::vector<GRMEntry> unionGrm = loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), cfg.spgrmGrabFile, cfg.spgrmGctaFile);

    // ── 7. Build per-phenotype SPAsqrMethod + PhenoTask ─────────────
    std::vector<PhenoTask> tasks(K);
    std::vector<std::vector<double> > allOutlierRatios(K);

    for (int k = 0; k < K; ++k) {
        infoMsg("[%s] Building SPAsqr method (%d taus, %u subjects)",
                phenoNames[k].c_str(), ntaus, pw[k].nk);

        // Re-index GRM to pheno-dense space
        auto phenoGrm = reindexGrm(unionGrm, pw[k].unionToLocal, nUnion);

        Eigen::MatrixXd omega;
        auto method = buildSPAsqrMethod(
            residMats[k],
            phenoGrm,
            pw[k].nk,
            cfg.spaCutoff,
            cfg.outlierIqrRatio,
            cfg.outlierAbsBound,
            tauLabels,
            &allOutlierRatios[k],
            cfg.writeOmega ? &omega : nullptr,
            ntaus
        );

        if (cfg.writeOmega)
            writeOmegaFile(cfg.outPrefix + "." + phenoNames[k] + ".SPAsqr.omega",
                           phenoNames[k], tauLabels, omega);

        tasks[k].phenoName = phenoNames[k];
        tasks[k].method = std::move(method);
        tasks[k].unionToLocal = pw[k].unionToLocal;
        tasks[k].nUsed = pw[k].nk;
    }

    logOutlierTable("Outlier ratios", phenoNames, tauLabels, allOutlierRatios,
                    cfg.outlierIqrRatio, cfg.outlierAbsBound);

    // Free per-phenotype work data
    pw.clear();

    // ── 8. Run multi-phenotype engine ───────────────────────────────
    infoMsg("SPAsqr: Starting multi-phenotype association (%d phenotypes, %d taus, %d threads)", K, ntaus, cfg.nthreads);
    multiPhenoEngine(
        *genoData, tasks, cfg.outPrefix, "SPAsqr", cfg.compression, cfg.compressionLevel,
        cfg.nthreads, cfg.missingCutoff, cfg.minMafCutoff, cfg.minMacCutoff, cfg.hweCutoff
    );
}

// ══════════════════════════════════════════════════════════════════════
// runSPAsqrLoco — LOCO entry point
//
// Per-chromosome workflow:
//   - Y is pre-transformed once via applyPhenoTransform (raw/int/standardize).
//   - For each chromosome, the buildTasks callback rebuilds the conquer fits
//     using y_adj = Y_transformed - loco_chr (β=1, α=0).
//   - All modes feed conquer with X = pw[k].X (original covariates only,
//     no LOCO augmentation).
//   - Per-chr bandwidth h_chr[k] = IQR(y_adj) / scale (chr-specific).
//   - locoEngine iterates chromosomes, testing each chromosome's markers via
//     the multi-phenotype work-stealing engine.
//   - The LOCO PRS scale must match the chosen transform — caller is
//     responsible for ensuring this (dispatch.cpp warns on raw + LOCO).
// ══════════════════════════════════════════════════════════════════════

void runSPAsqrLoco(const SPAsqrConfig &cfg) {
    const auto &phenoNames = cfg.phenoNames;
    const auto &taus       = cfg.taus;
    const int K     = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());
    const double qmmeTol = cfg.spasqrTol;

    // ── 1-3. Analysis set, covariates, per-phenotype split, bandwidth ──
    // The bandwidth set here is only the global reference value; buildTasks
    // recomputes it per chromosome from IQR(y_adj).
    ScoreSetup setup = buildScoreSetup(cfg, "SPAsqr-LOCO");
    SubjectData &sd = *setup.sd;
    std::vector<PhenoWork> &pw = setup.pw;
    const uint32_t nUnion = sd.nUsed();
    const int nCov        = setup.nCov;
    const double hScale   = setup.hScale;

    // QMME solver per phenotype: caches X^T X / n once; per-chr Cholesky
    // rebuilt inside buildTasks below. Read-only on solve() so worker
    // threads can share the same instance per phenotype.
    std::vector<std::unique_ptr<qmme::SqrSolver> > qmmeSolvers(K);
    for (int k = 0; k < K; ++k)
        qmmeSolvers[k] = std::make_unique<qmme::SqrSolver>(pw[k].X, /*delta*/ 1e-6);

    // Per-chr h_chr = IQR(y_adj)/scale is computed inside the chr loop; the value
    // shown here is the GLOBAL reference h derived from IQR(pw[k].Y) (post-transform).
    // Actual per-chr h is logged via "SPAsqr-LOCO chr%s [...]: per-chr h = ...".
    logPhenoTable(phenoNames, pw, "global_h(ref)");

    // ── 4. Load LOCO predictions (union space) ─────────────────────
    LocoData loco = LocoData::load(cfg.predListFile, phenoNames, sd.usedIIDs(), sd.famIIDs());
    auto locoChroms = loco.availableChromosomes();
    infoMsg("LOCO: %zu chromosomes available across all phenotypes", locoChroms.size());

    // ── 5. Build tau labels ─────────────────────────────────────────
    const std::vector<std::string> tauLabels = makeTauLabels(taus);

    // ── 6. Load genotype data and GRM once (union space) ────────────
    auto genoData = makeGenoData(cfg.geno, sd.usedMask(), sd.nFam(), sd.nUsed(), cfg.nSnpPerChunk);

    std::vector<GRMEntry> unionGrm = loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), cfg.spgrmGrabFile, cfg.spgrmGctaFile);

    // Pre-compute per-phenotype re-indexed GRM (shared across chromosomes)
    std::vector<std::vector<GRMEntry> > phenoGrms(K);
    for (int k = 0; k < K; ++k)
        phenoGrms[k] = reindexGrm(unionGrm, pw[k].unionToLocal, nUnion);

    // ── 7. Build LocoTaskBuilder callback ──────────────────────────
    // For each chromosome, build per-pheno y_adj (mode-dependent), run
    // K × ntaus QMME fits with X = baseX, and build K SPAsqrMethods.
    auto buildTasks = [&](const std::string &chr, std::vector<PhenoTask> &tasks) {
        tasks.resize(K);

        // y_adj = pw[k].Y(transformed) - loco_chr_pheno_dense  (β=1, α=0)
        // The transform was applied once in §3; only the per-chr LOCO subtraction
        // happens here. Per-chr bandwidth h_chr[k] = IQR(y_adj) / scale.
        std::vector<Eigen::VectorXd> y_adjs(K);
        std::vector<double> h_chr(K, 0.0);

        for (int k = 0; k < K; ++k) {
            const Eigen::VectorXd loco_dense = locoDense(
                loco.scores.at(phenoNames[k]).at(chr), pw[k], phenoNames[k], chr,
                "SPAsqr-LOCO");

            y_adjs[k] = pw[k].Y - loco_dense;

            // Diagnostic: corr(Y_transformed, loco_chr)^2 in pheno-dense space.
            double r2_loco = 0.0;
            {
                const double mean_y = pw[k].Y.mean();
                const double mean_l = loco_dense.mean();
                const Eigen::ArrayXd dev_y = pw[k].Y.array() - mean_y;
                const Eigen::ArrayXd dev_l = loco_dense.array() - mean_l;
                const double cov_yl  = (dev_y * dev_l).sum();
                const double var_y   = dev_y.square().sum();
                const double var_l   = dev_l.square().sum();
                if (var_y > 0.0 && var_l > 0.0)
                    r2_loco = (cov_yl * cov_yl) / (var_y * var_l);
            }

            // Per-chromosome bandwidth from IQR(y_adj).  --spasqr-h overrides auto-IQR.
            h_chr[k] = (cfg.spasqrH >= 0.0) ? cfg.spasqrH
                                        : iqrBandwidth(y_adjs[k], hScale, nCov);

            infoMsg("SPAsqr-LOCO chr%s [%s][transform=%s]: per-chr h = IQR(y_adj)/scale = %.6f (Nk=%u),"
                    " corr(Y, loco)^2 = %.4f",
                    chr.c_str(), phenoNames[k].c_str(), cfg.phenoTransform.c_str(),
                    h_chr[k], pw[k].nk, r2_loco);
        }

        // Rebuild Cholesky for each phenotype's per-chr bandwidth.
        // Single-threaded; subsequent solve() calls are read-only.
        for (int k = 0; k < K; ++k)
            qmmeSolvers[k]->prepareBandwidth(h_chr[k]);

        // Parallel QMME fits: K × ntaus
        std::vector<Eigen::MatrixXd> ResidMats(K);
        for (int k = 0; k < K; ++k)
            ResidMats[k].resize(static_cast<Eigen::Index>(pw[k].nk), ntaus);

        runQmmeFits(K, ntaus, cfg.nthreads, phenoNames, tauLabels,
                    "SPAsqr-LOCO chr" + chr, [&](int k, int t) {
            const double h1 = 1.0 / h_chr[k];
            Eigen::VectorXd resid;
            qmme::SolverStatus st;
            // SQR on (X, y_adjs[k]) with y_adjs[k] = Y_transformed - loco_chr.
            // Y_transformed comes from §3 applyPhenoTransform (raw/int/standardize).
            qmmeSolvers[k]->solve(y_adjs[k], taus[t], &resid, qmmeTol,
                                  /*maxIter*/ 50000, /*restartPeriod*/ 50, &st);
            if (!st.converged)
                warnMsg("[%s] chr%s tau=%.4f qmme did not converge: %d iters (||g||=%.3e); tol=%.3e",
                        phenoNames[k].c_str(), chr.c_str(), taus[t],
                        st.iter, st.finalGradNorm, qmmeTol);

            const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
            for (Eigen::Index i = 0; i < Nk; ++i)
                ResidMats[k](i, t) = taus[t] - math::pnorm(-resid(i) * h1);
        });

        // Build SPAsqrMethod for each phenotype (pheno-dense space)
        std::vector<std::vector<double> > allOutlierRatios(K);
        for (int k = 0; k < K; ++k) {
            Eigen::MatrixXd omega;
            auto method = buildSPAsqrMethod(
                ResidMats[k],
                phenoGrms[k],
                pw[k].nk,
                cfg.spaCutoff,
                cfg.outlierIqrRatio,
                cfg.outlierAbsBound,
                tauLabels,
                &allOutlierRatios[k],
                cfg.writeOmega ? &omega : nullptr,
                ntaus
            );

            if (cfg.writeOmega)
                writeOmegaFile(cfg.outPrefix + "." + phenoNames[k] + ".SPAsqr.chr" + chr + ".omega",
                               phenoNames[k], tauLabels, omega);

            tasks[k].phenoName = phenoNames[k];
            tasks[k].method = std::move(method);
            tasks[k].unionToLocal = pw[k].unionToLocal;
            tasks[k].nUsed = pw[k].nk;
        }

        logOutlierTable("chr" + chr + " outlier ratios", phenoNames, tauLabels,
                        allOutlierRatios, cfg.outlierIqrRatio, cfg.outlierAbsBound);
    };

    // ── 8. Run LOCO engine ─────────────────────────────────────────
    const int nChroms = static_cast<int>(locoChroms.size());
    infoMsg("SPAsqr-LOCO: Starting LOCO association (%d phenotypes, %d taus, %d chroms, %d threads)",
            K, ntaus, nChroms, cfg.nthreads);
    locoEngine(
        *genoData,
        locoChroms,
        phenoNames,
        buildTasks,
        cfg.outPrefix,
        "SPAsqr",
        cfg.compression,
        cfg.compressionLevel,
        cfg.nthreads,
        cfg.missingCutoff,
        cfg.minMafCutoff,
        cfg.minMacCutoff,
        cfg.hweCutoff
    );
}

// ══════════════════════════════════════════════════════════════════════
// runLoQus — LoQus (LOngitudinal QUantile Score test): long-format phenotype, marginal GEE null model
//
// Input protocol is the --longitudinal one shared with SPACox / SPAmix /
// SPAGRM: the kept set GRM ∩ keep − remove is built on the .fam IIDs, the
// long-format --pheno file (FID/IID + named columns) is read against it, and
// SubjectData then applies the standard union checks (genotype ∩ GRM ∩
// keep/remove) to the subjects that have records.
//
// Per phenotype, records with a missing value of that phenotype or of any
// covariate are dropped (so each phenotype has its own record set); the time
// column is required, orders each subject's records, and (IID, time) must be
// unique.  Whether a covariate varies over time is not declared — the
// estimating equation does not need it.
//
// Each phenotype is reduced to one weight per subject and column,
// a_i = 1ᵀ R_i⁻¹ r_i (see long_gee.hpp): one column per quantile level
// (smoothed quantile GEE) and/or one "linear" column (identity-link GEE).
// That matrix enters buildSPAsqrMethod exactly as a per-subject residual
// matrix, so the score, the GRM variance R^T Φ R, the Binomial / empirical
// CGF saddlepoint and the output contract are SPAsqr's own.  LOG10P_CCT
// combines the quantile columns only.
//
// With --pred-list the chromosome's LOCO PGS enters the null model as a
// covariate (its coefficient is estimated), not as an offset; the GEE is
// refitted per chromosome.
// ══════════════════════════════════════════════════════════════════════

namespace {

// Sort each subject's records by time and reject duplicated (IID, time).
void orderRecordsByTime(nsSAGELDFit::LongPhenoData &d, const std::string &timeName,
                        const std::string &file) {
    const int nSubj = static_cast<int>(d.uniqueIIDs.size());
    std::vector<uint32_t> perm;
    for (int i = 0; i < nSubj; ++i) {
        const uint32_t a = d.subjStart[i], b = d.subjStart[i + 1];
        perm.resize(b - a);
        std::iota(perm.begin(), perm.end(), a);
        std::stable_sort(perm.begin(), perm.end(),
                         [&](uint32_t u, uint32_t v) { return d.E(u) < d.E(v); });
        for (size_t t = 1; t < perm.size(); ++t)
            if (d.E(perm[t]) == d.E(perm[t - 1]))
                throw std::runtime_error(
                    file + ": IID " + d.uniqueIIDs[i] + " has two records with " + timeName +
                    " = " + std::to_string(d.E(perm[t])) +
                    "; (IID, time) must identify a record uniquely");
        const Eigen::MatrixXd X = d.X.middleRows(a, b - a);
        const Eigen::MatrixXd Y = d.Y.middleRows(a, b - a);
        const Eigen::VectorXd E = d.E.segment(a, b - a);
        for (size_t t = 0; t < perm.size(); ++t) {
            d.X.row(a + t) = X.row(perm[t] - a);
            d.Y.row(a + t) = Y.row(perm[t] - a);
            d.E(a + t) = E(perm[t] - a);
        }
    }
}

// One phenotype's longitudinal analysis set in union coordinates.
struct LongWork {
    PhenoWork index;          // unionToLocal / nk only (reused by locoDense)
    longgee::Records rec;     // covariate records (no PGS)
    double h = 0.0;
};

} // namespace

void runLoQus(const SPAsqrConfig &cfg) {
    const std::string label = "LoQus";
    const auto &phenoNames = cfg.phenoNames;
    const auto &taus = cfg.taus;
    const int K = static_cast<int>(phenoNames.size());
    const bool doQ = cfg.geeModel == "qr" || cfg.geeModel == "both";
    const bool doL = cfg.geeModel == "linear" || cfg.geeModel == "both";
    const int nQ = doQ ? static_cast<int>(taus.size()) : 0;
    const int nCol = nQ + (doL ? 1 : 0);
    const longgee::WorkingCorr corr = (cfg.workingCorr == "exchangeable")
        ? longgee::WorkingCorr::Exchangeable : longgee::WorkingCorr::Independence;

    std::vector<std::string> colLabels;
    if (doQ) colLabels = makeTauLabels(taus);
    if (doL) colLabels.push_back("linear");

    infoMsg("%s: gee-model = %s, working-corr = %s, time = %s, %d phenotype(s), %d column(s)",
            label.c_str(), cfg.geeModel.c_str(), cfg.workingCorr.c_str(), cfg.timeName.c_str(),
            K, nCol);

    // ── 1. Kept set and long-format records (one read per phenotype) ────
    std::vector<std::string> famIIDs = parseGenoIIDs(cfg.geno);
    std::unordered_set<std::string> grmIDs =
        SparseGRM::parseSubjectIDs(cfg.spgrmGrabFile, cfg.spgrmGctaFile, famIIDs);
    const auto kept = nsLongitudinal::buildKeptSet(cfg.keepFile, cfg.removeFile, famIIDs, grmIDs);

    std::vector<nsSAGELDFit::LongPhenoData> ld(K);
    std::unordered_set<std::string> withRecords;
    for (int k = 0; k < K; ++k) {
        infoMsg("[%s] reading long-format records", phenoNames[k].c_str());
        ld[k] = nsSAGELDFit::parseLongPheno(cfg.phenoFile, {phenoNames[k]}, cfg.covarNames,
                                            cfg.timeName, famIIDs, kept, /*envRequired*/ true,
                                            "--time-name");
        orderRecordsByTime(ld[k], cfg.timeName, cfg.phenoFile);
        withRecords.insert(ld[k].uniqueIIDs.begin(), ld[k].uniqueIIDs.end());
    }

    // ── 2. Union: standard SubjectData checks on the subjects with records ─
    SubjectData sd(famIIDs);
    sd.setKeepSubjects(std::move(withRecords));
    sd.setKeepRemove(cfg.keepFile, cfg.removeFile);
    sd.setGrmSubjects(std::move(grmIDs));
    sd.setGenoLabel(cfg.geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(cfg.spgrmGrabFile, cfg.spgrmGctaFile));
    sd.finalize();
    const uint32_t nUnion = sd.nUsed();
    const std::vector<std::string> usedIIDs = sd.usedIIDs();
    std::unordered_map<std::string, uint32_t> unionIdx;
    unionIdx.reserve(usedIIDs.size() * 2);
    for (uint32_t i = 0; i < nUnion; ++i) unionIdx.emplace(usedIIDs[i], i);

    // ── 3. Per-phenotype records in union order ─────────────────────────
    const int nCov = static_cast<int>(cfg.covarNames.size());
    const double hScale = (cfg.spasqrHScale >= 0.0) ? cfg.spasqrHScale : 5.0;
    std::vector<LongWork> lw(K);
    for (int k = 0; k < K; ++k) {
        const auto &d = ld[k];
        LongWork &w = lw[k];
        w.index.unionToLocal.assign(nUnion, UINT32_MAX);
        // Subjects are in .fam order in both d and the union, so local
        // indices assigned in d's order are increasing in union order.
        std::vector<int> keepSubj;
        uint32_t nRows = 0;
        for (int i = 0; i < static_cast<int>(d.uniqueIIDs.size()); ++i) {
            auto it = unionIdx.find(d.uniqueIIDs[i]);
            if (it == unionIdx.end()) continue;
            w.index.unionToLocal[it->second] = static_cast<uint32_t>(keepSubj.size());
            keepSubj.push_back(i);
            nRows += d.subjStart[i + 1] - d.subjStart[i];
        }
        w.index.nk = static_cast<uint32_t>(keepSubj.size());
        if (w.index.nk == 0)
            throw std::runtime_error(label + ": phenotype '" + phenoNames[k] +
                                     "' has no records for any analysed subject");
        w.rec.X.resize(nRows, nCov);
        w.rec.y.resize(nRows);
        w.rec.start.assign(1, 0);
        uint32_t r = 0;
        for (int i : keepSubj) {
            for (uint32_t t = d.subjStart[i]; t < d.subjStart[i + 1]; ++t, ++r) {
                if (nCov > 0) w.rec.X.row(r) = d.X.row(t).tail(nCov);
                w.rec.y(r) = d.Y(t, 0);
            }
            w.rec.start.push_back(r);
        }
        w.h = (cfg.spasqrH >= 0.0) ? cfg.spasqrH : iqrBandwidth(w.rec.y, hScale, nCov);
        longgee::checkDesign(w.rec, cfg.covarNames, label + ": phenotype '" + phenoNames[k] + "'");
        infoMsg("[%s] %u subjects, %u records (mean %.2f, max %d per subject), h = %.6g",
                phenoNames[k].c_str(), w.index.nk, nRows,
                static_cast<double>(nRows) / w.index.nk, w.rec.maxRecords(), w.h);
        ld[k] = nsSAGELDFit::LongPhenoData{};   // release the parsed copy
    }

    // ── 4. Genotypes and GRM (union space), per-phenotype GRM ───────────
    auto genoData = makeGenoData(cfg.geno, sd.usedMask(), sd.nFam(), sd.nUsed(), cfg.nSnpPerChunk);
    std::vector<GRMEntry> unionGrm =
        loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), cfg.spgrmGrabFile, cfg.spgrmGctaFile);
    std::vector<std::vector<GRMEntry> > phenoGrms(K);
    for (int k = 0; k < K; ++k) phenoGrms[k] = reindexGrm(unionGrm, lw[k].index.unionToLocal, nUnion);

    // Absolute outlier bound, in units of the weight SD (the weights are
    // centred and scaled to unit SD before outlier detection, below): |a_i| >
    // 1.2 is always an outlier, on top of the IQR rule.  SPAsqr's own bound,
    // 0.55 on the psi scale, cannot be reused as is — a_i is a sum over
    // records with no fixed range — but it is needed for the same reason.  A
    // subject with one record has a_i = psi, which is two-valued (τ and τ−1,
    // blurred by h); at τ = 0.3 the IQR fences fall inside the two clusters,
    // so the IQR rule alone flags nobody, the CGF is the Gaussian block only,
    // and the saddlepoint degenerates to the normal approximation.  Measured
    // with every m_i = 1 at MAC 200–400 (5000 subjects, 100 null phenotypes):
    // LOG10P exceeded the SPA value on every test, by +0.05 at 1e-4..1e-5 and
    // +0.2 below 1e-6 (τ 0.3/0.7; τ 0.5 unaffected).  1.2 SD is the psi-scale
    // 0.55 at τ = 0.3 or 0.7 (sd(psi) ≈ √(τ(1−τ)) = 0.46), so with m_i = 1
    // the outlier sets match plain SPAsqr's there.
    constexpr double kAbsBoundSd = 1.2;
    const double absBound = kAbsBoundSd;

    // ── 5. Null-model fits → per-phenotype SPAsqr methods ───────────────
    auto buildTasks = [&](const std::vector<longgee::Records> &recs, const std::string &tag,
                          std::vector<PhenoTask> &tasks) {
        tasks.resize(K);
        std::vector<Eigen::MatrixXd> weights(K);
        std::vector<std::vector<longgee::NullFit> > fits(K, std::vector<longgee::NullFit>(nCol));
        std::vector<std::unique_ptr<qmme::SqrSolver> > solvers(K);
        std::vector<bool> informative(K);
        for (int k = 0; k < K; ++k) {
            weights[k].resize(static_cast<Eigen::Index>(lw[k].index.nk), nCol);
            informative[k] = longgee::hasResidualVariation(recs[k]);
            if (!informative[k])
                warnMsg("[%s]%s outcome is fully explained by the null-model covariates;"
                        " reporting NA (SPA_STATUS=8) for all association tests",
                        phenoNames[k].c_str(), tag.c_str());
            if (nQ > 0 && informative[k]) {
                solvers[k] = std::make_unique<qmme::SqrSolver>(recs[k].X, /*delta*/ 1e-6);
                solvers[k]->prepareBandwidth(lw[k].h);
            }
        }
        runQmmeFits(K, nCol, cfg.nthreads, phenoNames, colLabels, label + tag, [&](int k, int c) {
            if (!informative[k]) {
                weights[k].col(c).setZero();
                return;
            }
            longgee::NullFit f = (c < nQ)
                ? longgee::fitQuantile(recs[k], *solvers[k], taus[c], lw[k].h, corr, cfg.spasqrTol)
                : longgee::fitLinear(recs[k], corr);
            weights[k].col(c) = f.weight;
            fits[k][c] = std::move(f);
        });

        // Centre every weight column and rescale it to unit SD.
        //
        // Centring: with an intercept in the null model Σ a_i is already ≈ 0
        // (the intercept row of the estimating equation is Σ_i 1ᵀR_i⁻¹r_i,
        // proportional to Σ a_i), to the fit's convergence tolerance
        // (measured |Σa|/‖a‖: ~1e-9 independence, ≤ 5e-7 exchangeable).  The
        // score is unaffected either way (the genotype is centred), but the
        // variance term aᵀΦa and the CGF's non-outlier mean are not, so the
        // residual Σ a_i is removed explicitly rather than left to tolerance.
        //
        // Scaling: the test is invariant to the overall scale of a_i, but the
        // empirical-CGF table
        // (util/empirical_cgf) tabulates K0 on a t-grid of fixed absolute
        // width, which is only unit-free for a dimensionless residual such as
        // psi.  The linear-GEE weight carries phenotype units; without this
        // step rescaling the phenotype moved low-MAC LOG10P_linear by up to
        // 0.35.  Dividing by the SD makes every column dimensionless; after
        // it, Y x 1024 moves no LOG10P column beyond output rounding (1e-5).
        for (int k = 0; k < K; ++k)
            for (int c = 0; c < nCol; ++c) {
                auto col = weights[k].col(c);
                const double mean = col.mean();
                col.array() -= mean;
                const double sd = std::sqrt(col.squaredNorm() / static_cast<double>(col.size()));
                if (sd > 0.0) col /= sd;
            }

        std::vector<std::vector<double> > outlierRatios(K);
        for (int k = 0; k < K; ++k) {
            std::ostringstream rhos;
            for (int c = 0; c < nCol; ++c) {
                const auto &f = fits[k][c];
                rhos << ' ' << colLabels[c] << '=' << std::setprecision(4) << f.rho;
                if (!f.qmmeConverged)
                    warnMsg("[%s]%s %s: QMME (independence start) did not converge",
                            phenoNames[k].c_str(), tag.c_str(), colLabels[c].c_str());
                if (f.fellBack)
                    warnMsg("[%s]%s %s: exchangeable fit failed (%s); using INDEPENDENCE weights"
                            " for this column", phenoNames[k].c_str(), tag.c_str(),
                            colLabels[c].c_str(), f.reason.c_str());
            }
            if (corr == longgee::WorkingCorr::Exchangeable)
                infoMsg("[%s]%s working rho:%s", phenoNames[k].c_str(), tag.c_str(), rhos.str().c_str());

            Eigen::MatrixXd omega;
            auto method = buildSPAsqrMethod(weights[k], phenoGrms[k], lw[k].index.nk,
                                            cfg.spaCutoff, cfg.outlierIqrRatio, absBound,
                                            colLabels, &outlierRatios[k],
                                            cfg.writeOmega ? &omega : nullptr, nQ);
            if (cfg.writeOmega)
                writeOmegaFile(cfg.outPrefix + "." + phenoNames[k] + ".LoQus" +
                                   (tag.empty() ? std::string() : "." + tag.substr(1)) + ".omega",
                               phenoNames[k], colLabels, omega);
            tasks[k].phenoName = phenoNames[k];
            tasks[k].method = std::move(method);
            tasks[k].unionToLocal = lw[k].index.unionToLocal;
            tasks[k].nUsed = lw[k].index.nk;
        }
        logOutlierTable((tag.empty() ? std::string("Outlier ratios") : tag.substr(1) + " outlier ratios"),
                        phenoNames, colLabels, outlierRatios, cfg.outlierIqrRatio, absBound);
    };

    if (cfg.predListFile.empty()) {
        // Without LOCO the records are not needed again: hand them over
        // rather than hold a second copy of every phenotype's design.
        std::vector<longgee::Records> recs(K);
        for (int k = 0; k < K; ++k) recs[k] = std::move(lw[k].rec);
        std::vector<PhenoTask> tasks;
        buildTasks(recs, "", tasks);
        infoMsg("%s: starting association (%d phenotypes, %d columns, %d threads)",
                label.c_str(), K, nCol, cfg.nthreads);
        multiPhenoEngine(*genoData, tasks, cfg.outPrefix, "LoQus", cfg.compression,
                         cfg.compressionLevel, cfg.nthreads, cfg.missingCutoff, cfg.minMafCutoff,
                         cfg.minMacCutoff, cfg.hweCutoff);
        return;
    }

    // ── LOCO: the chromosome's PGS is a covariate of the null model ─────
    LocoData loco = LocoData::load(cfg.predListFile, phenoNames, sd.usedIIDs(), sd.famIIDs());
    auto locoChroms = loco.availableChromosomes();
    infoMsg("LOCO: %zu chromosomes available across all phenotypes; PGS enters the null model"
            " as a covariate", locoChroms.size());
    auto locoRecords = [&](const std::string &chr) {
        std::vector<longgee::Records> recs(K);
        for (int k = 0; k < K; ++k) {
            const Eigen::VectorXd pgs = locoDense(loco.scores.at(phenoNames[k]).at(chr),
                                                  lw[k].index, phenoNames[k], chr, label.c_str());
            const longgee::Records &base = lw[k].rec;
            longgee::Records &rc = recs[k];
            rc.start = base.start;
            rc.y = base.y;
            rc.X.resize(base.X.rows(), base.X.cols() + 1);
            rc.X.leftCols(base.X.cols()) = base.X;
            for (int i = 0; i < base.nSubj(); ++i)
                rc.X.col(base.X.cols()).segment(base.start[i], base.start[i + 1] - base.start[i])
                    .setConstant(pgs(i));
        }
        return recs;
    };
    // The PGS is a covariate, so it must not make the design rank-deficient
    // (e.g. a PGS that is constant on some chromosome).  Checked before any
    // association output is written, for exactly the chromosomes locoEngine
    // will process: those with markers AND a PGS.  A pred file routinely
    // carries all 22 autosomes, and the ones without markers here (often all
    // zero, as in LDAK-KVIK output for unanalysed chromosomes) are never used.
    {
        std::vector<std::string> names = cfg.covarNames;
        names.push_back("LOCO PGS");
        std::vector<std::string> activeChroms;
        {
            std::unordered_set<std::string> seen;
            for (const auto &m : genoData->markerInfo())
                if (locoChroms.count(m.chrom) && seen.insert(m.chrom).second)
                    activeChroms.push_back(m.chrom);
        }
        for (const auto &chr : activeChroms) {
            const auto recs = locoRecords(chr);
            for (int k = 0; k < K; ++k)
                longgee::checkDesign(recs[k], names, label + ": phenotype '" + phenoNames[k] +
                                                         "', chromosome " + chr);
        }
    }
    auto locoBuild = [&](const std::string &chr, std::vector<PhenoTask> &tasks) {
        buildTasks(locoRecords(chr), " chr" + chr, tasks);
    };
    infoMsg("%s: starting LOCO association (%d phenotypes, %d columns, %zu chroms, %d threads)",
            label.c_str(), K, nCol, locoChroms.size(), cfg.nthreads);
    locoEngine(*genoData, locoChroms, phenoNames, locoBuild, cfg.outPrefix, "LoQus",
               cfg.compression, cfg.compressionLevel, cfg.nthreads, cfg.missingCutoff,
               cfg.minMafCutoff, cfg.minMacCutoff, cfg.hweCutoff);
}
