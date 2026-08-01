// spasqr.cpp — SPAsqr: SPA-squared multi-tau marker association (pure C++17 / Eigen / Boost)

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "spasqr/spasqr.hpp"
#include "engine/loco.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "spasqr/qmme.hpp"
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
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <Eigen/Dense>

// No <immintrin.h>, util/simd_dispatch.hpp or util/simd_math.hpp any more: the
// five hand-written CGF copies that needed them are gone, and the scalar / AVX2
// / AVX-512 triple now lives once in util/spa_cgf.cpp with its own dispatch
// site.  See the block comment below.

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
// SPAsqr is not exposed".  That is wrong.  EmpVar is built from
// tau.R_GRM_R_nonOutlier, accumulated as
// `e.factor * e.value * ResidMat(e.row, col) * ResidMat(e.col, col)` over a
// thresholded sparse GRM read from a file, and such a matrix is not guaranteed
// positive semidefinite.  A negative EmpVar makes Score_adj NaN, and a negative
// `var` term makes the CGF's K'' negative — the most reachable route to
// std::sqrt of a negative number in the old tail.  Both now surface as a named
// spa::Status with P = NA.

namespace {

// ── Per-tau SPA data for SPAsqr (shared across clones via shared_ptr) ──

struct SPAsqrPerTau {
    Eigen::VectorXd outlierResid;        // residuals of outlier subjects
    double sum_unrelated_outliers2;      // outlierResid.squaredNorm()
    double sum_R_nonOutlier;
    double R_GRM_R_nonOutlier;
    double R_GRM_R;                      // full variance term (all subjects)
};

// `zeta` is gone: it was assigned 0.01 in buildSPAsqrMethod and never read,
// because both tails derived their initial abscissa from
// min(|Score_adj|/Score_var, 1.2) instead.  See 01_findings.md D7's companion
// note and the lower-tail convention comment in spagrm.cpp.
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

// ── Per-tau p-value — the single live implementation ───────────────────
//
// G_var = empirical Var(G) over post-impute genotypes, passed in by the caller,
// drives the outer Score_var = G_var * R_GRM_R.  EmpVar = K''(0) uses the
// model-based 2p(1-p) for both the non-outlier partial-Gaussian piece and the
// outlier Bernoulli factor; the rescaling
// Score_adj = Score * sqrt(K''(0) / Var(S)_observed) is the SAIGE-style
// variance ratio that brings the score onto the CGF's natural scale, dropping
// the implicit HWE assumption (G_var ~ 2p(1-p)) earlier code relied on.

inline spa::Result tauPvalue(
    double Score,
    double altFreq,
    double G_var,
    double &zScore,
    const SPAsqrPerTau &tau,
    const SPAsqrSPAShared &spa_
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double MAF       = std::min(altFreq, 1.0 - altFreq);
    const double Score_var = G_var * tau.R_GRM_R;

    if (Score_var <= 0.0 || MAF <= 0.0) {
        zScore = 0.0;
        return spa::Result{nan, spa::Status::NaNoTest};
    }

    zScore = Score / std::sqrt(Score_var);

    if (!std::isfinite(zScore))
        return spa::Result{nan, spa::Status::NaNoTest};

    if (std::abs(zScore) <= spa_.SPA_Cutoff) return spa::normalBranch(zScore);

    const double EmpVar    = 2.0 * MAF * (1.0 - MAF) *
                              (tau.R_GRM_R_nonOutlier + tau.sum_unrelated_outliers2);
    const double Score_adj = Score * std::sqrt(EmpVar / Score_var);
    const double absAdj    = std::abs(Score_adj);

    // First-order saddlepoint estimate zeta ~ s/K''(0), capped; the lower tail
    // starts at its negation.
    const double initZeta = std::min(absAdj / Score_var, 1.2);

    return spaTwoSided(tau.outlierResid.data(),
                       static_cast<int>(tau.outlierResid.size()),
                       MAF,
                       2.0 * MAF * tau.sum_R_nonOutlier,
                       2.0 * MAF * (1.0 - MAF) * tau.R_GRM_R_nonOutlier,
                       absAdj, initZeta, spa_.tol, zScore);
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
        std::vector<std::string> tauLabels
    )
        : m_ntaus(ntaus),
          m_spaShared(std::move(spaShared)),
          m_tauLabels(std::move(tauLabels)),
          m_hasLabels(true)
    {
        auto sd = std::make_shared<SharedMethodData>();
        sd->residMat  = std::move(residMat);
        sd->residSums = std::move(residSums);
        m_methodShared = std::move(sd);
    }

    std::unique_ptr<MethodBase> clone() const override {
        // Both shared_ptr's are read-only — no per-clone scratch at all.
        return std::make_unique<SPAsqrMethod>(*this);
    }

    int resultSize() const override {
        return 5 * m_ntaus + 1;
    }

// P_CCT, then five per-tau groups: P, LOG10P, Z, Z_Norm, SPA_STATUS.
//
// The saddlepoint is per tau, so both new columns are per tau as well: a marker
// can converge at one quantile level and fail at another, and a single
// aggregated status would hide which.  The LOG10P group is inserted directly
// after the P group and the SPA_STATUS group appended last, so the relative
// order of the pre-existing columns is unchanged and tests/regress.py reports
// exactly one structural line per file.
//
// There is deliberately no LOG10P_CCT.  P_CCT comes from math::cauchyCombine,
// which has no log-domain form — its statistic is a mean of tan((0.5-p)*pi)
// evaluated on the linear scale — so a LOG10P_CCT column could only ever be
// -log10 of the linear P_CCT and would carry no information the P_CCT column
// does not already carry.  Giving CCT a genuine log-domain path is a separate
// change to math_helper.hpp, not part of this migration.
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
        oss << "\tP_CCT";
        static const char *const kGroups[] = {"P_", "LOG10P_", "Z_", "Z_Norm_",
                                              "SPA_STATUS_"};
        for (const char *g : kGroups) {
            for (int i = 0; i < m_ntaus; ++i) {
                oss << '\t' << g;
                if (m_hasLabels) oss << m_tauLabels[i];
                else             oss << "tau" << (i + 1);
            }
        }
        return oss.str();
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int /*markerInChunkIdx*/,
        std::vector<double> &result
    ) override {
        result.clear();

        const double n      = static_cast<double>(GVec.size());
        const double gSum   = GVec.sum();
        const double gMean  = gSum / n;
        // Empirical Var(G) over post-impute genotypes: (Σg² − (Σg)²/n) / (n−1).
        const double gSumSq = GVec.squaredNorm();
        const double G_var  = (n > 1.0)
            ? (gSumSq - gSum * gSum / n) / (n - 1.0)
            : 0.0;

        Eigen::VectorXd scores = m_methodShared->residMat.transpose() * GVec;
        for (int i = 0; i < m_ntaus; ++i)
            scores[i] -= gMean * m_methodShared->residSums[i];

        processOneMarker(scores.data(), altFreq, G_var, result);
    }

    // ── Batched analysis: B markers at once ────────────────────────────
    void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        const std::vector<int> & /*chunkIdxs*/,
        std::vector<std::vector<double> > &results
    ) override {
        const int B = static_cast<int>(GBatch.cols());
        results.resize(B);

        Eigen::MatrixXd scoreMatrix;
        scoreMatrix.noalias() = m_methodShared->residMat.transpose() * GBatch;

        const Eigen::VectorXd gMeans  = GBatch.colwise().mean();
        const Eigen::VectorXd gSumSqs = GBatch.colwise().squaredNorm();
        scoreMatrix.noalias() -= m_methodShared->residSums * gMeans.transpose();

        const double n = static_cast<double>(GBatch.rows());
        for (int b = 0; b < B; ++b) {
            // Empirical Var(G) over post-impute genotypes: (Σg² − n·ḡ²) / (n−1).
            const double G_var = (n > 1.0)
                ? (gSumSqs[b] - n * gMeans[b] * gMeans[b]) / (n - 1.0)
                : 0.0;
            processOneMarker(scoreMatrix.col(b).data(), altFreqs[b], G_var, results[b]);
        }
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

        // ── C2: Tau-first SPA computation ──────────────────────────
        // Process all B markers for one tau before moving to the next.
        // This keeps the tau's outlierResid data hot in cache.
        // Layout: buf[b * ntaus + t].
        // Heap-backed per-thread buffers (was a [20*16] stack array) so the
        // fused GEMM window B is not capped at 16/B·ntaus ≤ 320.
        const size_t nCell = static_cast<size_t>(B) * m_ntaus;
        m_zBuf.resize(nCell);
        m_pBuf.resize(nCell);
        m_lpBuf.resize(nCell);
        m_stBuf.resize(nCell);

        for (int t = 0; t < m_ntaus; ++t) {
            const SPAsqrPerTau &tau = m_spaShared->perTau[t];

            for (int b = 0; b < B; ++b) {
                // Empirical Var(G) over post-impute genotypes:
                //   Var(g) = (Σg² − (Σg)²/n) / (n−1)
                const double gSum  = gSums[b];
                const double G_var = (nUsedD > 1.0)
                    ? (gSumSqs[b] - gSum * gSum / nUsedD) / (nUsedD - 1.0)
                    : 0.0;

                double z;
                const spa::Result ts = tauPvalue(m_centeredBuf(t, b), altFreqs[b],
                                                   G_var, z, tau, *m_spaShared);
                const size_t k = static_cast<size_t>(b) * m_ntaus + t;
                m_zBuf[k]  = z;
                // P from LOG10P: the tier stops producing a linear tail in
                // log10p_unify Stage 3, and the column goes in Stage 8.
                m_pBuf[k]  = spa::pFromNegLog10P(ts.negLog10p);
                m_lpBuf[k] = ts.negLog10p;
                m_stBuf[k] = static_cast<double>(static_cast<uint8_t>(ts.status));
            }
        }

        for (int b = 0; b < B; ++b) {
            const size_t off = static_cast<size_t>(b) * m_ntaus;
            assemble(m_pBuf.data() + off, m_lpBuf.data() + off,
                     m_zBuf.data() + off, m_stBuf.data() + off, results[b]);
        }
    }

  private:
    int m_ntaus;
    std::shared_ptr<const SPAsqrSPAShared> m_spaShared;
    std::vector<std::string> m_tauLabels;
    bool m_hasLabels;

    struct SharedMethodData {
        Eigen::MatrixXd residMat;   // N × ntaus
        Eigen::VectorXd residSums;  // ntaus
    };

    std::shared_ptr<const SharedMethodData> m_methodShared;
    Eigen::MatrixXd m_centeredBuf;  // reused across processScoreBatch calls
    // B × ntaus each, reused across processScoreBatch calls.
    std::vector<double> m_zBuf, m_pBuf, m_lpBuf, m_stBuf;

    // Column assembly, shared by all three entry points so the header and the
    // row can only ever be built by the same code.
    void assemble(
        const double *pvals,
        const double *lgs,
        const double *zScores,
        const double *stats,
        std::vector<double> &result
    ) const {
        result.clear();
        result.reserve(5 * m_ntaus + 1);

        // math::cauchyCombine skips NaN entries, so a tau whose saddlepoint
        // failed drops out of the combination rather than poisoning it; that
        // behaviour is unchanged by this stage.
        result.push_back(math::cauchyCombine(pvals, m_ntaus));

        for (int i = 0; i < m_ntaus; ++i) result.push_back(pvals[i]);
        for (int i = 0; i < m_ntaus; ++i) result.push_back(lgs[i]);
        for (int i = 0; i < m_ntaus; ++i) {
            const double zr = zScores[i];
            // D3: on the fast path (|z| ≤ SPA_Cutoff) the p-value is exactly
            // 2·Φ(−|z|), so zFromPval inverts back to the raw z; emit it
            // directly and skip the long-double qnorm round-trip (which also
            // removes the spurious round-trip ULP gap that made Z_τ differ from
            // Z_Norm_τ on fast-path entries).  NaN p (degenerate marker or a
            // saddlepoint failure) propagates as NaN, matching zFromPval.
            if (std::isnan(pvals[i]))
                result.push_back(std::numeric_limits<double>::quiet_NaN());
            else if (std::abs(zr) <= m_spaShared->SPA_Cutoff)
                result.push_back(zr);                                 // Z_τ == Z_Norm_τ
            else
                result.push_back(math::zFromPval(pvals[i], zr));      // SPA-recalibrated
        }
        for (int i = 0; i < m_ntaus; ++i) result.push_back(zScores[i]);  // Z_Norm_τ
        for (int i = 0; i < m_ntaus; ++i) result.push_back(stats[i]);    // SPA_STATUS
    }

    void processOneMarker(
        const double *centeredScores,
        double altFreq,
        double G_var,
        std::vector<double> &result
    ) {
        // Stack arrays — ntaus is always small (≤20).
        double zScores[20], pvals[20], lgs[20], stats[20];

        for (int i = 0; i < m_ntaus; ++i) {
            double z;
            const spa::Result ts = tauPvalue(centeredScores[i], altFreq, G_var, z,
                                               m_spaShared->perTau[i], *m_spaShared);
            zScores[i] = z;
            pvals[i]   = spa::pFromNegLog10P(ts.negLog10p);
            lgs[i]     = ts.negLog10p;
            stats[i]   = static_cast<double>(static_cast<uint8_t>(ts.status));
        }

        assemble(pvals, lgs, zScores, stats, result);
    }

};

// ══════════════════════════════════════════════════════════════════════
// Outlier detection (IQR-based, per column)
// ══════════════════════════════════════════════════════════════════════

// Returns an N × K boolean matrix (as std::vector<std::vector<bool>>).
// outlierIqrRatio  = multiplier for IQR (default 1.5)
// outlierAbsBound  = absolute clamp for cutoffs (default 0.55)
struct OutlierInfo {
    // per-column: indices of outlier subjects
    std::vector<std::vector<int> > outlierIdx;
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
    info.outlierIdx.resize(K);
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
                info.outlierIdx[col].push_back(static_cast<int>(i));
            }
        }
    }
    return info;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// Shared pipeline: outlier detection → GRM → SPAGRM → marker engine
// ══════════════════════════════════════════════════════════════════════

// GRMEntry is now declared in spasqr.hpp

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

// Build SPAsqrMethod from a pre-computed residual matrix and pre-loaded GRM entries.
std::unique_ptr<MethodBase> buildSPAsqrMethod(
    Eigen::MatrixXd &ResidMat,
    const std::vector<GRMEntry> &grmEntries,
    uint32_t nUsed,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    double /*minMafCutoff*/,
    double /*minMacCutoff*/,
    std::vector<std::string> tauLabels,
    std::vector<double> *outlierRatiosOut
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
            (*outlierRatiosOut)[c] = static_cast<double>(outlierInfo.outlierIdx[c].size()) / N;
    }

    // ── 4. Compute per-column variance terms + build SPAsqrSPAShared ──
    auto spaShared = std::make_shared<SPAsqrSPAShared>();
    spaShared->perTau.resize(ntaus);
    spaShared->SPA_Cutoff = spaCutoff;
    spaShared->tol        = 1e-6;

    for (int col = 0; col < ntaus; ++col) {
        const auto &isOut = outlierInfo.isOutlier[col];
        auto &pt = spaShared->perTau[col];

        // R_GRM_R and R_GRM_R_nonOutlier
        double rgrm_r = 0.0;
        double rgrm_r_no = 0.0;
        for (const auto &e : grmEntries) {
            const double contrib = e.factor * e.value * ResidMat(e.row, col) * ResidMat(e.col, col);
            rgrm_r += contrib;
            if (!isOut[e.row] && !isOut[e.col]) rgrm_r_no += contrib;
        }
        pt.R_GRM_R            = rgrm_r;
        pt.R_GRM_R_nonOutlier = rgrm_r_no;

        // sum_R_nonOutlier and outlier residual values
        double sumNO = 0.0;
        std::vector<double> outVals;
        outVals.reserve(outlierInfo.outlierIdx[col].size());
        for (uint32_t i = 0; i < nUsed; ++i) {
            if (!isOut[i]) sumNO += ResidMat(i, col);
            else outVals.push_back(ResidMat(i, col));
        }
        pt.sum_R_nonOutlier = sumNO;
        pt.outlierResid =
            Eigen::Map<Eigen::VectorXd>(outVals.data(), static_cast<Eigen::Index>(outVals.size()));
        pt.sum_unrelated_outliers2 = pt.outlierResid.squaredNorm();
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
        std::move(tauLabels)
    );

    return method;
}

// ══════════════════════════════════════════════════════════════════════
// Re-index GRM entries from union-dense space to pheno-dense space
// ══════════════════════════════════════════════════════════════════════

static std::vector<GRMEntry> reindexGrm(
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

// ── Phenotype pre-transform helper ──────────────────────────────────
// Applied per-phenotype on the non-missing scope.
//   raw         — Y unchanged
//   int         — inverse normal transform (Blom, average-rank ties)
//   standardize — (Y − mean) / sd  (sample sd, n−1 denom)
// Caller is responsible for validating the mode string upstream
// (dispatch.cpp); this function asserts on unknown modes.
static void applyPhenoTransform(
    Eigen::VectorXd &Y,
    const std::string &mode
) {
    if (mode == "raw") return;
    if (mode == "int") {
        Y = math::inverseRankNormal(Y);
        return;
    }
    if (mode == "standardize") {
        const Eigen::Index n = Y.size();
        if (n <= 1) return;
        const double mean = Y.mean();
        const double ssq  = (Y.array() - mean).square().sum();
        const double sd   = std::sqrt(ssq / static_cast<double>(n - 1));
        if (sd > 0.0) Y = (Y.array() - mean) / sd;
        return;
    }
    throw std::runtime_error("applyPhenoTransform: unknown mode '" + mode + "'");
}

// ══════════════════════════════════════════════════════════════════════
// runSPAsqr — multi-phenotype entry point with parallel QMME fits
// ══════════════════════════════════════════════════════════════════════

void runSPAsqr(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<double> &taus,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    double spasqrTol,
    double spasqrH,
    double spasqrHScale,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &phenoTransform
) {
    const int K = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());

    // ── 1. Load phenotype/covariate data (union mask) ───────────────
    // Union = subjects with genotype ∩ GRM ∩ keep/remove.  Per-phenotype
    // NA filtering is deferred — each phenotype uses its own non-missing
    // subset of the union.
    infoMsg("SPAsqr: pheno-transform = %s, solver = qmme",
            phenoTransform.c_str());
    infoMsg("SPAsqr: Loading phenotype and covariate data (%d phenotypes, %d taus)", K, ntaus);
    // Score mode fits the null model once and reuses it for every marker, so
    // --spasqr-tol (default 1e-6) is tight enough; apply it directly.
    const double qmmeTol = spasqrTol;
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile, phenoNames);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();

    const uint32_t nUnion = sd.nUsed();
    const Eigen::Index N = static_cast<Eigen::Index>(nUnion);

    // ── 2. Extract union-space covariates ───────────────────────────
    Eigen::MatrixXd unionX = covarNames.empty()
        ? (sd.hasCovar() ? Eigen::MatrixXd(sd.covar()) : Eigen::MatrixXd(N, 0))
        : sd.getColumns(covarNames);
    const int nCov = static_cast<int>(unionX.cols());

    // ── 3. Per-phenotype: build non-missing mask, extract Y/X ───────
    struct PhenoWork {
        std::vector<uint32_t> unionToLocal; // size nUnion; UINT32_MAX = absent
        uint32_t nk;                        // non-missing count
        Eigen::VectorXd Y;                  // nk
        Eigen::MatrixXd X;                  // nk × nCov
        double h;                           // bandwidth
        Eigen::MatrixXd ResidMat;           // nk × ntaus (filled by QMME)
    };

    std::vector<PhenoWork> pw(K);

    for (int k = 0; k < K; ++k) {
        Eigen::VectorXd fullY = sd.getColumn(phenoNames[k]);
        pw[k].unionToLocal.resize(nUnion, UINT32_MAX);
        uint32_t localIdx = 0;
        for (uint32_t i = 0; i < nUnion; ++i) {
            if (!std::isnan(fullY[i])) {
                pw[k].unionToLocal[i] = localIdx++;
            }
        }
        pw[k].nk = localIdx;
        if (pw[k].nk == 0)
            throw std::runtime_error("SPAsqr: phenotype '" + phenoNames[k] + "' has no non-missing subjects");

        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        pw[k].Y.resize(Nk);
        pw[k].X.resize(Nk, nCov);
        for (uint32_t i = 0; i < nUnion; ++i) {
            uint32_t li = pw[k].unionToLocal[i];
            if (li == UINT32_MAX) continue;
            pw[k].Y[li] = fullY[i];
            if (nCov > 0) pw[k].X.row(li) = unionX.row(i);
        }

        // Apply pheno transform (raw / int / standardize) — per-phenotype non-missing scope.
        applyPhenoTransform(pw[k].Y, phenoTransform);

        // Per-phenotype bandwidth
        if (spasqrH >= 0.0) {
            pw[k].h = spasqrH;
        } else {
            std::vector<double> ysorted(Nk);
            Eigen::VectorXd::Map(ysorted.data(), Nk) = pw[k].Y;
            std::sort(ysorted.begin(), ysorted.end());
            auto quantile = [&](double prob) -> double {
                double idx = prob * (Nk - 1);
                Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
                Eigen::Index hi = std::min(lo + 1, Nk - 1);
                double frac = idx - lo;
                return ysorted[lo] * (1.0 - frac) + ysorted[hi] * frac;
            };
            double iqr = quantile(0.75) - quantile(0.25);
            double scale = (spasqrHScale >= 0.0) ? spasqrHScale : 3.0;
            pw[k].h = iqr / scale;
            if (pw[k].h <= 0.0)
                pw[k].h = std::max(std::pow((std::log(Nk) + nCov) / static_cast<double>(Nk), 0.4), 0.05);
        }

        pw[k].ResidMat.resize(Nk, ntaus);
    }

    // Log N/bandwidth table
    {
        infoMsg("Sample size and smooth bandwidth per phenotype:");
        size_t nameW = 4;
        for (int k = 0; k < K; ++k)
            nameW = std::max(nameW, phenoNames[k].size());
        nameW += 2;
        char row[256];
        std::snprintf(row, sizeof(row), "    %-*s %10s %12s\n", static_cast<int>(nameW), "", "N", "bandwidth");
        fprintf(stderr, "%s", row);
        for (int k = 0; k < K; ++k) {
            std::snprintf(row, sizeof(row), "    %-*s %10u %12.6f\n",
                          static_cast<int>(nameW), phenoNames[k].c_str(), pw[k].nk, pw[k].h);
            fprintf(stderr, "%s", row);
        }
    }

    // ── 4. Parallel QMME fits: K × ntaus ────────────────────────────
    // Pre-construct one QMME solver per phenotype (X^T X / n cached) and
    // prepare Cholesky for that phenotype's bandwidth. Solver is then
    // shared read-only across worker threads — solve() only reads m_chol.
    std::vector<std::unique_ptr<qmme::SqrSolver> > qmmeSolvers(K);
    for (int k = 0; k < K; ++k) {
        qmmeSolvers[k] = std::make_unique<qmme::SqrSolver>(pw[k].X, /*delta*/ 1e-6);
        qmmeSolvers[k]->prepareBandwidth(pw[k].h);
    }

    const int totalFits = K * ntaus;
    const int nWorkers = std::min(nthreads, totalFits);
    infoMsg("SPAsqr: Running %d QMME fits with %d threads", totalFits, nWorkers);

    std::atomic<int> nextFit{0};
    std::vector<std::string> fitErrors(totalFits);

    auto fitWorker = [&]() {
        for (;;) {
            int idx = nextFit.fetch_add(1, std::memory_order_relaxed);
            if (idx >= totalFits) break;
            int k = idx / ntaus;
            int t = idx % ntaus;
            const double h1 = 1.0 / pw[k].h;

            try {
                Eigen::VectorXd resid;
                qmme::SolverStatus st;
                Eigen::VectorXd beta = qmmeSolvers[k]->solve(
                    pw[k].Y, taus[t], &resid,
                    qmmeTol, /*maxIter*/ 50000,
                    /*restartPeriod*/ 50, &st);
                infoMsg("[%s] tau=%.4f intercept=%.6f", phenoNames[k].c_str(), taus[t], beta(0));
                if (!st.converged) {
                    warnMsg("[%s] tau=%.4f qmme did not converge: %d iters (||g||=%.3e); tol=%.3e",
                            phenoNames[k].c_str(), taus[t],
                            st.iter, st.finalGradNorm, qmmeTol);
                }

                const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
                for (Eigen::Index i = 0; i < Nk; ++i)
                    pw[k].ResidMat(i, t) = taus[t] - math::pnorm(-resid(i) * h1);
            } catch (const std::exception &ex) {
                fitErrors[idx] = ex.what();
            }
        }
    };

    {
        std::vector<std::thread> threads;
        threads.reserve(nWorkers - 1);
        for (int t = 0; t < nWorkers - 1; ++t)
            threads.emplace_back(fitWorker);
        fitWorker();
        for (auto &th : threads)
            th.join();
    }

    for (int idx = 0; idx < totalFits; ++idx) {
        if (!fitErrors[idx].empty()) {
            int k = idx / ntaus;
            int t = idx % ntaus;
            throw std::runtime_error("SPAsqr: QMME failed for phenotype '" +
                                     phenoNames[k] + "' tau=" + std::to_string(taus[t]) + ": " + fitErrors[idx]);
        }
    }

    // ── 5. Build tau labels ────────────────────────────────────────────
    std::vector<std::string> tauLabels;
    tauLabels.reserve(ntaus);
    for (double tau : taus) {
        std::ostringstream oss;
        oss << "tau" << tau;
        tauLabels.push_back(oss.str());
    }

    // ── 6. Load genotype data and GRM once (shared, union space) ────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    std::vector<GRMEntry> unionGrm = loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), spgrmGrabFile, spgrmGctaFile);

    // ── 7. Build per-phenotype SPAsqrMethod + PhenoTask ─────────────
    std::vector<PhenoTask> tasks(K);
    std::vector<std::vector<double> > allOutlierRatios(K);

    for (int k = 0; k < K; ++k) {
        infoMsg("[%s] Building SPAsqr method (%d taus, %u subjects)",
                phenoNames[k].c_str(), ntaus, pw[k].nk);

        // Re-index GRM to pheno-dense space
        auto phenoGrm = reindexGrm(unionGrm, pw[k].unionToLocal, nUnion);

        auto method = buildSPAsqrMethod(
            pw[k].ResidMat,
            phenoGrm,
            pw[k].nk,
            spaCutoff,
            outlierIqrRatio,
            outlierAbsBound,
            minMafCutoff,
            minMacCutoff,
            tauLabels,
            &allOutlierRatios[k]
        );

        tasks[k].phenoName = phenoNames[k];
        tasks[k].method = std::move(method);
        tasks[k].unionToLocal = pw[k].unionToLocal;
        tasks[k].nUsed = pw[k].nk;
    }

    // Print outlier ratio table line-by-line
    {
        size_t nameW = 4;
        for (int k = 0; k < K; ++k)
            nameW = std::max(nameW, phenoNames[k].size());
        nameW += 2;

        infoMsg("Outlier ratios (IQR=%.2f, bound=%.2f):", outlierIqrRatio, outlierAbsBound);

        std::ostringstream hdr;
        hdr << "    " << std::setw(static_cast<int>(nameW)) << std::left << "";
        for (const auto &tl : tauLabels)
            hdr << std::setw(10) << std::right << tl;
        fprintf(stderr, "%s\n", hdr.str().c_str());

        for (int k = 0; k < K; ++k) {
            std::ostringstream row;
            row << "    " << std::setw(static_cast<int>(nameW)) << std::left << phenoNames[k];
            for (double r : allOutlierRatios[k])
                row << std::setw(10) << std::right << std::fixed << std::setprecision(4) << r;
            fprintf(stderr, "%s\n", row.str().c_str());
        }
    }

    // Free per-phenotype work data
    pw.clear();

    // ── 8. Run multi-phenotype engine ───────────────────────────────
    infoMsg("SPAsqr: Starting multi-phenotype association (%d phenotypes, %d taus, %d threads)", K, ntaus, nthreads);
    multiPhenoEngine(
        *genoData, tasks, outPrefix, "SPAsqr", compression, compressionLevel,
        nthreads, missingCutoff, minMafCutoff, minMacCutoff, hweCutoff
    );
}

// ══════════════════════════════════════════════════════════════════════
// runSPAsqrLoco — LOCO entry point
//
// Per-chromosome workflow:
//   - Y is pre-transformed once via applyPhenoTransform (raw/int/standardize).
//   - For each chromosome, the buildTasks callback rebuilds the conquer fits
//     using y_adj = Y_transformed - loco_chr (β=1, α=0).
//   - All modes feed conquer with X = pw[k].baseX (original covariates only,
//     no LOCO augmentation).
//   - Per-chr bandwidth h_chr[k] = IQR(y_adj) / scale (chr-specific).
//   - locoEngine iterates chromosomes, testing each chromosome's markers via
//     the multi-phenotype work-stealing engine.
//   - The LOCO PRS scale must match the chosen transform — caller is
//     responsible for ensuring this (dispatch.cpp warns on raw + LOCO).
// ══════════════════════════════════════════════════════════════════════

void runSPAsqrLoco(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<double> &taus,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &predListFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    double spasqrTol,
    double spasqrH,
    double spasqrHScale,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &phenoTransform
) {
    const int K = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());

    // ── 1. Load phenotype/covariate data (union mask) ───────────────
    // Union = subjects with genotype ∩ GRM ∩ keep/remove.  Per-phenotype
    // NA filtering is deferred — each phenotype uses its own non-missing
    // subset of the union.
    infoMsg("SPAsqr-LOCO pheno-transform: %s, solver: qmme",
            phenoTransform.c_str());
    infoMsg("SPAsqr-LOCO: Loading phenotype and covariate data (%d phenotypes, %d taus)", K, ntaus);
    const double qmmeTol = spasqrTol;
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile, phenoNames);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();

    const uint32_t nUnion = sd.nUsed();
    const Eigen::Index N = static_cast<Eigen::Index>(nUnion);

    // ── 2. Extract union-space covariates ───────────────────────────
    Eigen::MatrixXd unionX = covarNames.empty()
        ? (sd.hasCovar() ? Eigen::MatrixXd(sd.covar()) : Eigen::MatrixXd(N, 0))
        : sd.getColumns(covarNames);
    const int nCov = static_cast<int>(unionX.cols());

    // ── 3. Per-phenotype: build non-missing mask, extract Y/baseX ──
    struct PhenoWork {
        std::vector<uint32_t> unionToLocal; // size nUnion; UINT32_MAX = absent
        uint32_t nk;                        // non-missing count
        Eigen::VectorXd Y;                  // nk
        Eigen::MatrixXd baseX;              // nk × nCov
        double h;                           // bandwidth
    };

    std::vector<PhenoWork> pw(K);

    for (int k = 0; k < K; ++k) {
        Eigen::VectorXd fullY = sd.getColumn(phenoNames[k]);
        pw[k].unionToLocal.resize(nUnion, UINT32_MAX);
        uint32_t localIdx = 0;
        for (uint32_t i = 0; i < nUnion; ++i) {
            if (!std::isnan(fullY[i])) {
                pw[k].unionToLocal[i] = localIdx++;
            }
        }
        pw[k].nk = localIdx;
        if (pw[k].nk == 0)
            throw std::runtime_error("SPAsqr-LOCO: phenotype '" + phenoNames[k] + "' has no non-missing subjects");

        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        pw[k].Y.resize(Nk);
        pw[k].baseX.resize(Nk, nCov);
        for (uint32_t i = 0; i < nUnion; ++i) {
            uint32_t li = pw[k].unionToLocal[i];
            if (li == UINT32_MAX) continue;
            pw[k].Y[li] = fullY[i];
            if (nCov > 0) pw[k].baseX.row(li) = unionX.row(i);
        }

        // Apply pheno transform (raw / int / standardize) — per-phenotype non-missing scope.
        // All downstream math (bandwidth, QMME fit, y_adj per chromosome) operates
        // on the transformed Y. The LOCO PRS scale must match the chosen transform.
        applyPhenoTransform(pw[k].Y, phenoTransform);

        // Per-phenotype bandwidth
        if (spasqrH >= 0.0) {
            pw[k].h = spasqrH;
        } else {
            std::vector<double> ysorted(Nk);
            Eigen::VectorXd::Map(ysorted.data(), Nk) = pw[k].Y;
            std::sort(ysorted.begin(), ysorted.end());
            auto quantile = [&](double prob) -> double {
                double idx = prob * (Nk - 1);
                Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
                Eigen::Index hi = std::min(lo + 1, Nk - 1);
                double frac = idx - lo;
                return ysorted[lo] * (1.0 - frac) + ysorted[hi] * frac;
            };
            double iqr = quantile(0.75) - quantile(0.25);
            double scale = (spasqrHScale >= 0.0) ? spasqrHScale : 3.0;
            pw[k].h = iqr / scale;
            if (pw[k].h <= 0.0)
                pw[k].h = std::max(std::pow((std::log(Nk) + nCov) / static_cast<double>(Nk), 0.4), 0.05);
        }
    }

    // QMME solver per phenotype: caches X^T X / n once; per-chr Cholesky
    // rebuilt inside buildTasks below. Read-only on solve() so worker
    // threads can share the same instance per phenotype.
    std::vector<std::unique_ptr<qmme::SqrSolver> > qmmeSolvers(K);
    for (int k = 0; k < K; ++k)
        qmmeSolvers[k] = std::make_unique<qmme::SqrSolver>(pw[k].baseX, /*delta*/ 1e-6);

    // Log N/bandwidth table
    // Per-chr h_chr = IQR(y_adj)/scale is computed inside the chr loop; the value
    // shown here is the GLOBAL reference h derived from IQR(pw[k].Y) (post-transform).
    // Actual per-chr h is logged via "SPAsqr-LOCO chr%s [...]: per-chr h = ...".
    {
        const char *bwLabel = "global_h(ref)";
        infoMsg("Sample size and smooth bandwidth per phenotype:");
        size_t nameW = 4;
        for (int k = 0; k < K; ++k)
            nameW = std::max(nameW, phenoNames[k].size());
        nameW += 2;
        char row[256];
        std::snprintf(row, sizeof(row), "    %-*s %10s %14s\n",
                      static_cast<int>(nameW), "", "N", bwLabel);
        fprintf(stderr, "%s", row);
        for (int k = 0; k < K; ++k) {
            std::snprintf(row, sizeof(row), "    %-*s %10u %14.6f\n",
                          static_cast<int>(nameW), phenoNames[k].c_str(), pw[k].nk, pw[k].h);
            fprintf(stderr, "%s", row);
        }
    }

    // ── 4. Load LOCO predictions (union space) ─────────────────────
    LocoData loco = LocoData::load(predListFile, phenoNames, sd.usedIIDs(), sd.famIIDs());
    auto locoChroms = loco.availableChromosomes();
    infoMsg("LOCO: %zu chromosomes available across all phenotypes", locoChroms.size());

    // ── 5. Build tau labels ─────────────────────────────────────────
    std::vector<std::string> tauLabels;
    tauLabels.reserve(ntaus);
    for (double tau : taus) {
        std::ostringstream oss;
        oss << "tau" << tau;
        tauLabels.push_back(oss.str());
    }

    // ── 6. Load genotype data and GRM once (union space) ────────────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    std::vector<GRMEntry> unionGrm = loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), spgrmGrabFile, spgrmGctaFile);

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
            const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);

            // Map LOCO scores from union → pheno-dense (Nk-length) space.
            const auto &locoVec = loco.scores.at(phenoNames[k]).at(chr);
            Eigen::VectorXd loco_dense(Nk);
            for (uint32_t i = 0; i < nUnion; ++i) {
                uint32_t li = pw[k].unionToLocal[i];
                if (li != UINT32_MAX)
                    loco_dense[li] = locoVec[i];
            }

            // LDAK Step 1 / Regenie Step 1 output should include a PGS for every
            // subject in the analysis set. If any non-missing-Y subject is absent
            // from the LOCO file, parseLdakLocoFile / parseRegenieLocoFile leaves
            // NaN at that position. Hard-fail instead of silently corrupting the
            // per-chr QMME fit (NaN propagates through the QMME solve).
            if (!loco_dense.allFinite()) {
                const Eigen::Index nBad = Nk - loco_dense.array().isFinite().count();
                throw std::runtime_error(
                    "SPAsqr-LOCO: LOCO file for phenotype '" + phenoNames[k] +
                    "' chr " + chr + " is missing " + std::to_string(nBad) +
                    " subject(s) that have non-missing Y. The LOCO PGS file must "
                    "contain every subject in the --pheno analysis set. Re-run "
                    "LDAK / Regenie Step 1 on the same sample set, or remove "
                    "those subjects from --pheno.");
            }

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
            if (spasqrH >= 0.0) {
                h_chr[k] = spasqrH;
            } else {
                std::vector<double> ysorted(static_cast<size_t>(Nk));
                Eigen::VectorXd::Map(ysorted.data(), Nk) = y_adjs[k];
                std::sort(ysorted.begin(), ysorted.end());
                auto quant = [&](double prob) -> double {
                    double idx = prob * (Nk - 1);
                    Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
                    Eigen::Index hi = std::min(lo + 1, Nk - 1);
                    double frac = idx - lo;
                    return ysorted[lo] * (1.0 - frac) + ysorted[hi] * frac;
                };
                double iqr = quant(0.75) - quant(0.25);
                double scale = (spasqrHScale >= 0.0) ? spasqrHScale : 3.0;
                h_chr[k] = iqr / scale;
                if (h_chr[k] <= 0.0)
                    h_chr[k] = std::max(std::pow((std::log(Nk) + nCov) / static_cast<double>(Nk), 0.4), 0.05);
            }

            infoMsg("SPAsqr-LOCO chr%s [%s][transform=%s]: per-chr h = IQR(y_adj)/scale = %.6f (Nk=%u),"
                    " corr(Y, loco)^2 = %.4f",
                    chr.c_str(), phenoNames[k].c_str(), phenoTransform.c_str(),
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

        const int totalFits = K * ntaus;
        const int nWorkers = std::min(nthreads, totalFits);
        infoMsg("SPAsqr-LOCO chr%s: Running %d QMME fits with %d threads",
                chr.c_str(), totalFits, nWorkers);

        std::atomic<int> nextFit{0};
        std::vector<std::string> fitErrors(totalFits);

        auto fitWorker = [&]() {
            for (;;) {
                int idx = nextFit.fetch_add(1, std::memory_order_relaxed);
                if (idx >= totalFits) break;
                int k = idx / ntaus;
                int t = idx % ntaus;
                const double h_use = h_chr[k];
                const double h1    = 1.0 / h_use;

                try {
                    Eigen::VectorXd resid;
                    qmme::SolverStatus st;
                    // SQR on (baseX, y_adjs[k]) with y_adjs[k] = Y_transformed - loco_chr.
                    // Y_transformed comes from §3 applyPhenoTransform (raw/int/standardize).
                    qmmeSolvers[k]->solve(y_adjs[k], taus[t], &resid,
                                          qmmeTol, /*maxIter*/ 50000,
                                          /*restartPeriod*/ 50, &st);
                    if (!st.converged) {
                        warnMsg("[%s] chr%s tau=%.4f qmme did not converge: %d iters (||g||=%.3e); tol=%.3e",
                                phenoNames[k].c_str(), chr.c_str(), taus[t],
                                st.iter, st.finalGradNorm, qmmeTol);
                    }

                    const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
                    for (Eigen::Index i = 0; i < Nk; ++i)
                        ResidMats[k](i, t) = taus[t] - math::pnorm(-resid(i) * h1);
                } catch (const std::exception &ex) {
                    fitErrors[idx] = ex.what();
                }
            }
        };

        {
            std::vector<std::thread> threads;
            threads.reserve(nWorkers - 1);
            for (int t = 0; t < nWorkers - 1; ++t)
                threads.emplace_back(fitWorker);
            fitWorker();
            for (auto &th : threads)
                th.join();
        }

        for (int idx = 0; idx < totalFits; ++idx) {
            if (!fitErrors[idx].empty()) {
                int k = idx / ntaus;
                int t = idx % ntaus;
                throw std::runtime_error("SPAsqr-LOCO chr" + chr + ": QMME failed for phenotype '" +
                                         phenoNames[k] + "' tau=" + std::to_string(taus[t]) + ": " + fitErrors[idx]);
            }
        }

        // Build SPAsqrMethod for each phenotype (pheno-dense space)
        std::vector<std::vector<double> > allOutlierRatios(K);
        for (int k = 0; k < K; ++k) {
            auto method = buildSPAsqrMethod(
                ResidMats[k],
                phenoGrms[k],
                pw[k].nk,
                spaCutoff,
                outlierIqrRatio,
                outlierAbsBound,
                minMafCutoff,
                minMacCutoff,
                tauLabels,
                &allOutlierRatios[k]
            );

            tasks[k].phenoName = phenoNames[k];
            tasks[k].method = std::move(method);
            tasks[k].unionToLocal = pw[k].unionToLocal;
            tasks[k].nUsed = pw[k].nk;
        }

        // Print outlier ratio table for this chromosome line-by-line
        {
            size_t nameW = 4;
            for (int k = 0; k < K; ++k)
                nameW = std::max(nameW, phenoNames[k].size());
            nameW += 2;

            infoMsg("chr%s outlier ratios (IQR=%.2f, bound=%.2f):",
                    chr.c_str(), outlierIqrRatio, outlierAbsBound);

            std::ostringstream hdr;
            hdr << "    " << std::setw(static_cast<int>(nameW)) << std::left << "";
            for (const auto &tl : tauLabels)
                hdr << std::setw(10) << std::right << tl;
            fprintf(stderr, "%s\n", hdr.str().c_str());

            for (int k = 0; k < K; ++k) {
                std::ostringstream row;
                row << "    " << std::setw(static_cast<int>(nameW)) << std::left << phenoNames[k];
                for (double r : allOutlierRatios[k])
                    row << std::setw(10) << std::right << std::fixed << std::setprecision(4) << r;
                fprintf(stderr, "%s\n", row.str().c_str());
            }
        }
    };

    // ── 8. Run LOCO engine ─────────────────────────────────────────
    const int nChroms = static_cast<int>(locoChroms.size());
    infoMsg("SPAsqr-LOCO: Starting LOCO association (%d phenotypes, %d taus, %d chroms, %d threads)",
            K, ntaus, nChroms, nthreads);
    locoEngine(
        *genoData,
        locoChroms,
        phenoNames,
        buildTasks,
        outPrefix,
        "SPAsqr",
        compression,
        compressionLevel,
        nthreads,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        hweCutoff
    );
}
