// spamixplus.cpp — Unified SPAmix / SPAmixPlus method implementation
//
// When a sparse GRM is provided, variance uses GRM-based covariance and
// the SPA tail is applied to variance-ratio-corrected S (SPAmixPlus).
// Without a GRM the diagonal variance Σ resid²·2·AF·(1−AF) is used (SPAmix).

#include "spamix/spamixplus.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <ctime>
#include <exception>
#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

#include <zlib.h>

#include "engine/loco.hpp"
#include "geno_factory/geno_data.hpp"
#include "spamix/indiv_af.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "spagrm/longitudinal_resid.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/null_model.hpp"

// ======================================================================
// SPAmixPlusMethod — construction / clone
// ======================================================================

// ── With GRM (SPAmixPlus) ──────────────────────────────────────────

// Pre-computed AF + GRM
SPAmixPlusMethod::SPAmixPlusMethod(
    const Eigen::VectorXd &residuals,
    const Eigen::VectorXd &resid2,
    const Eigen::MatrixXd &onePlusPCs,
    const OutlierData &outlier,
    double spaCutoff,
    const SparseGRM &grm,
    const std::vector<AFModel> &afModels,
    const std::vector<uint32_t> &genoToFlat,
    int maskIdx
)
    : m_resid(residuals),
      m_resid2(resid2),
      m_onePlusPCs(onePlusPCs),
      m_outlier(outlier),
      m_spaCutoff(spaCutoff),
      m_hasGRM(true),
      m_grm(&grm),
      m_N(static_cast<int>(residuals.size())),
      m_nPC(static_cast<int>(onePlusPCs.cols()) - 1),
      m_maskIdx(maskIdx),
      m_residSum(residuals.sum()),
      m_afModels(&afModels),
      m_genoToFlat(&genoToFlat),
      m_XtX_inv_Xt(nullptr),
      m_sqrt_XtX_inv_diag(nullptr),
      m_AFVec(m_N),
      m_WVec(m_N),
      m_R_new(m_N),
      m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
      m_chunkGenoIndices(nullptr)
{
}

// On-the-fly AF + GRM
SPAmixPlusMethod::SPAmixPlusMethod(
    const Eigen::VectorXd &residuals,
    const Eigen::VectorXd &resid2,
    const Eigen::MatrixXd &onePlusPCs,
    const OutlierData &outlier,
    double spaCutoff,
    const SparseGRM &grm,
    const Eigen::MatrixXd &XtX_inv_Xt,
    const Eigen::VectorXd &sqrt_XtX_inv_diag,
    int nPC,
    int maskIdx
)
    : m_resid(residuals),
      m_resid2(resid2),
      m_onePlusPCs(onePlusPCs),
      m_outlier(outlier),
      m_spaCutoff(spaCutoff),
      m_hasGRM(true),
      m_grm(&grm),
      m_N(static_cast<int>(residuals.size())),
      m_nPC(nPC),
      m_maskIdx(maskIdx),
      m_residSum(residuals.sum()),
      m_afModels(nullptr),
      m_genoToFlat(nullptr),
      m_XtX_inv_Xt(&XtX_inv_Xt),
      m_sqrt_XtX_inv_diag(&sqrt_XtX_inv_diag),
      m_AFVec(m_N),
      m_WVec(m_N),
      m_R_new(m_N),
      m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
      m_chunkGenoIndices(nullptr)
{
}

// ── Without GRM (SPAmix) ───────────────────────────────────────────

// Pre-computed AF, no GRM
SPAmixPlusMethod::SPAmixPlusMethod(
    const Eigen::VectorXd &residuals,
    const Eigen::VectorXd &resid2,
    const Eigen::MatrixXd &onePlusPCs,
    const OutlierData &outlier,
    double spaCutoff,
    const std::vector<AFModel> &afModels,
    const std::vector<uint32_t> &genoToFlat,
    int maskIdx
)
    : m_resid(residuals),
      m_resid2(resid2),
      m_onePlusPCs(onePlusPCs),
      m_outlier(outlier),
      m_spaCutoff(spaCutoff),
      m_hasGRM(false),
      m_grm(nullptr),
      m_N(static_cast<int>(residuals.size())),
      m_nPC(static_cast<int>(onePlusPCs.cols()) - 1),
      m_maskIdx(maskIdx),
      m_residSum(residuals.sum()),
      m_afModels(&afModels),
      m_genoToFlat(&genoToFlat),
      m_XtX_inv_Xt(nullptr),
      m_sqrt_XtX_inv_diag(nullptr),
      m_AFVec(m_N),
      m_WVec(m_N),
      m_R_new(0),
      m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
      m_chunkGenoIndices(nullptr)
{
}

// On-the-fly AF, no GRM
SPAmixPlusMethod::SPAmixPlusMethod(
    const Eigen::VectorXd &residuals,
    const Eigen::VectorXd &resid2,
    const Eigen::MatrixXd &onePlusPCs,
    const OutlierData &outlier,
    double spaCutoff,
    const Eigen::MatrixXd &XtX_inv_Xt,
    const Eigen::VectorXd &sqrt_XtX_inv_diag,
    int nPC,
    int maskIdx
)
    : m_resid(residuals),
      m_resid2(resid2),
      m_onePlusPCs(onePlusPCs),
      m_outlier(outlier),
      m_spaCutoff(spaCutoff),
      m_hasGRM(false),
      m_grm(nullptr),
      m_N(static_cast<int>(residuals.size())),
      m_nPC(nPC),
      m_maskIdx(maskIdx),
      m_residSum(residuals.sum()),
      m_afModels(nullptr),
      m_genoToFlat(nullptr),
      m_XtX_inv_Xt(&XtX_inv_Xt),
      m_sqrt_XtX_inv_diag(&sqrt_XtX_inv_diag),
      m_AFVec(m_N),
      m_WVec(m_N),
      m_R_new(0),
      m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
      m_chunkGenoIndices(nullptr)
{
}

std::unique_ptr<MethodBase> SPAmixPlusMethod::clone() const {
    std::unique_ptr<SPAmixPlusMethod> m;
    if (m_hasGRM) {
        if (m_XtX_inv_Xt) {
            m = std::make_unique<SPAmixPlusMethod>(
                m_resid,
                m_resid2,
                m_onePlusPCs,
                m_outlier,
                m_spaCutoff,
                *m_grm,
                *m_XtX_inv_Xt,
                *m_sqrt_XtX_inv_diag,
                m_nPC,
                m_maskIdx
            );
        } else {
            m = std::make_unique<SPAmixPlusMethod>(
                m_resid,
                m_resid2,
                m_onePlusPCs,
                m_outlier,
                m_spaCutoff,
                *m_grm,
                *m_afModels,
                *m_genoToFlat,
                m_maskIdx
            );
        }
    } else if (m_XtX_inv_Xt) {
        m = std::make_unique<SPAmixPlusMethod>(
            m_resid,
            m_resid2,
            m_onePlusPCs,
            m_outlier,
            m_spaCutoff,
            *m_XtX_inv_Xt,
            *m_sqrt_XtX_inv_diag,
            m_nPC,
            m_maskIdx
        );
    } else {
        m = std::make_unique<SPAmixPlusMethod>(
            m_resid,
            m_resid2,
            m_onePlusPCs,
            m_outlier,
            m_spaCutoff,
            *m_afModels,
            *m_genoToFlat,
            m_maskIdx
        );
    }

    // Bind a per-thread AFVec cache shared between SPAmix clones in the
    // same dedup mask group.  clone() is invoked from the worker thread,
    // so the thread_local map is per-worker; the shared_ptr's lifetime is
    // bounded by the lifetime of that worker thread.
    thread_local std::unordered_map<int, std::shared_ptr<SPAmixAFCache> > tlCacheMap;
    auto &slot = tlCacheMap[m_maskIdx];
    if (!slot) slot = std::make_shared<SPAmixAFCache>();
    m->m_afCache = slot;
    m->m_useAFCacheBatch = m_useAFCacheBatch;

    return m;
}

void SPAmixPlusMethod::prepareChunk(const std::vector<uint64_t> &gIndices) {
    m_chunkGenoIndices = &gIndices;
    // Invalidate any AFVec cache populated for an earlier chunk so the
    // first processScoreBatch of the new chunk repopulates it.
    if (m_afCache) m_afCache->lastChunkIdxs.clear();
}

// ======================================================================
// fillAFVecForMarker — populate per-individual AF for one marker
// ======================================================================

void SPAmixPlusMethod::fillAFVecForMarker(
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    int markerInChunkIdx,
    Eigen::Ref<Eigen::VectorXd> outAF
) {
    if (m_XtX_inv_Xt) {
        AFContext ctx{
            m_onePlusPCs,
            *m_XtX_inv_Xt,
            *m_sqrt_XtX_inv_diag,
            m_onePlusPCs.rightCols(m_nPC),
            m_N,
            m_nPC
        };
        computeAFVec(GVec, altFreq, ctx, outAF);
    } else {
        const uint64_t genoIdx = (*m_chunkGenoIndices)[markerInChunkIdx];
        const uint32_t flatIdx = (*m_genoToFlat)[genoIdx];
        const AFModel &model = (*m_afModels)[flatIdx];
        getAFVecFromModel(model, altFreq, m_onePlusPCs, m_N, outAF);
    }
}

// ======================================================================
// markerPvalFromAF — score test with optional GRM variance + SPA
//
//   rawScore  = Σ_i G_i · resid_i  (raw, uncentered; from fused GEMM or GEMV)
//   afVec     = per-individual ALT frequency
//   wVec      = per-individual W = 2·AF·(1−AF)
//
// Centering and variance use afVec / wVec; SPA uses the outlier split.
//
// spa_unify Stage 5.  The tail is now spamix_cgf::twoSidedSpa — the shared
// bracketed solver, the SIMD-dispatched spa_cgf::binomIndiv kernel and
// spa::bnTail's full guard set — in place of the two hand-written
// `spa::getProbSpaG` calls whose sum was returned unguarded.  Three
// behavioural consequences:
//
//   D5.  `fastGetRootK1` computed a convergence flag and returned it in
//        `RootResult::converge`; the caller read only `.root`
//        (common.cpp:100, :103, :123), so a saddlepoint that exhausted its
//        100-iteration budget produced an ordinary-looking finite p-value.
//        Worse, the non-finite exit at common.cpp:86-90 broke out with
//        `converge` still true, so the one genuinely divergent path was the
//        one reported as converged.  spa::Status is part of the return value
//        and reaches the user in SPA_STATUS.
//
//   D7.  The tolerance was a hard-coded `constexpr double tol = 0.001` tested
//        on the Newton STEP in zeta units (common.cpp:71, :97) — a thousand
//        times looser than the 1e-6 SPAGRM and SPAsqr use, with no CLI
//        exposure and no scaling.  It is now the solver's relative criterion
//        |K'(t) - s| <= rtol * sqrt(K''(t)) at rtol = 1e-6.  This is a real
//        numeric change and is the point of the stage.
//
//   Degenerate variance.  `VarS <= 0` returned p = 1.0 — a fabricated
//        perfectly-null p-value on a row whose Z, BETA and SE were all NaN.
//        L2 forbids substituting anything silently, so the row reports NA
//        with SPA_STATUS = 8 NA_NO_TEST: there is no statistic, hence nothing
//        for the D5 normal fallback to approximate.  Three markers of
//        examples/1kg reach it, all with ALT_FREQ > 0.99 and every
//        per-individual q̂ saturated at 1 by the AF model.
// ======================================================================

spa::TwoSided SPAmixPlusMethod::markerPvalFromAF(
    const Eigen::Ref<const Eigen::VectorXd> &afVec,
    const Eigen::Ref<const Eigen::VectorXd> &wVec,
    double rawScore,
    double &zScore,
    double &outVarS
) {
    // S_mean   = 2·Σ_i resid_i · AF_i
    // S_var_d  = Σ_i resid2_i · W_i             (diagonal variance)
    // When GRM is present, build R_new_i = resid_i · sqrt(W_i) and pass to
    // SparseGRM::spaVariance to obtain the GRM-based variance.
    double S_mean = 2.0 * m_resid.dot(afVec);
    double S_var_diag = m_resid2.dot(wVec);

    double VarS;
    if (m_hasGRM) {
        // R_new = resid · sqrt(W).  Use Eigen array ops so SIMD kicks in.
        m_R_new = m_resid.array() * wVec.array().sqrt();
        VarS = m_grm->spaVariance(m_R_new.data(), static_cast<uint32_t>(m_N));
    } else {
        VarS = S_var_diag;
    }

    if (VarS <= 0.0) {
        zScore = 0.0;
        outVarS = 0.0;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return spa::TwoSided{nan, nan, spa::Status::NaNoTest};
    }

    zScore = (rawScore - S_mean) / std::sqrt(VarS);
    outVarS = VarS;
    if (!std::isfinite(zScore)) {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        return spa::TwoSided{nan, nan, spa::Status::NaNoTest};
    }
    // spa::normalTwoSided routes through Boost's complement; for the standard
    // normal that is bit-for-bit `2 * pnorm(-|z|)`, since both reduce to
    // erfc(|z|/sqrt(2))/2 on the identical argument.  The branch therefore
    // reproduces the predecessor exactly and only gains LOG10P and the
    // NORMAL status.
    if (std::abs(zScore) <= m_spaCutoff) return spa::normalBranch(zScore);

    // ── SPA path with outlier / non-outlier split ──────────────────
    double S_spa, S_mean_spa;
    if (m_hasGRM) {
        double sqrt_ratio = std::sqrt(S_var_diag / VarS);
        S_spa = rawScore * sqrt_ratio;
        S_mean_spa = S_mean * sqrt_ratio;
    } else {
        S_spa = rawScore;
        S_mean_spa = S_mean;
    }

    // ── The Gaussian non-outlier block, by complement ──────────────
    //
    // The predecessor walked the NON-outlier set here — an indexed gather over
    // ~99 % of the cohort, per marker, inside the saddlepoint branch:
    //
    //     for (i < nNon) { af = afVec[posNonOutlier[i]];
    //                      mean += residNonOutlier[i]*af;
    //                      var  += resid2NonOutlier[i]*af*(1-af); }
    //
    // That loop is O(N) while everything else in the branch is O(nOutlier), and
    // on a 50000-subject synthetic cohort it measured at roughly twenty times
    // the cost of the entire CGF it was preparing — which is why replacing the
    // scalar CGF with the SIMD kernels moved the end-to-end saddlepoint time not
    // at all until this went with it.  01_findings.md does not record it.
    //
    // The two sums it produced are the complements of quantities already
    // computed above, over the whole cohort, by vectorized Eigen reductions:
    //
    //     S_mean     = 2 * Sum_all r_i * af_i
    //     S_var_diag =     Sum_all r_i^2 * w_i        w_i = 2 af_i (1 - af_i)
    //
    // so the non-outlier block is the total minus an O(nOutlier) sum, folded
    // into the gather that has to run anyway.  The outlier terms are formed from
    // `wVec[pos]` rather than from `2*af*(1-af)` recomputed, so that they are
    // bit-identical to the terms S_var_diag was accumulated from.
    //
    // Cancellation is bounded and mild.  Both subtractions take the total minus
    // the outliers' share f of it, so the relative error is amplified by
    // 1/(1-f); f is the IQR-outlier share of the variance, a few per cent in
    // practice, and even f = 0.9 leaves the error at 1e-13.  The variance is
    // clamped at zero because a negative one is not representable in the CGF and
    // could only arise when the true value is already below the noise floor of
    // the reduction, i.e. when the Gaussian block is contributing nothing.
    const int nOut = static_cast<int>(m_outlier.posOutlier.size());

    double outMean = 0.0, outVar = 0.0;
    for (int i = 0; i < nOut; ++i) {
        const uint32_t pos = m_outlier.posOutlier[i];
        const double af = afVec[pos];
        const double r = m_outlier.residOutlier[i];
        m_mafOutlier[i] = af;
        outMean += r * af;
        outVar += r * r * wVec[pos];
    }

    const double mean_nonOutlier = S_mean - 2.0 * outMean;
    double var_nonOutlier = S_var_diag - outVar;
    if (!(var_nonOutlier > 0.0)) var_nonOutlier = 0.0;

    const double absS = std::abs(S_spa - S_mean_spa);

    spamix_cgf::Context cgf;
    cgf.resid = m_outlier.residOutlier.data();
    cgf.af = m_mafOutlier.data();     // per-individual AF: binomIndiv
    cgf.nOutlier = nOut;
    cgf.mean = mean_nonOutlier;
    cgf.var = var_nonOutlier;

    // S_var_diag is the independence K''(0) and only sizes the initial
    // abscissa; the reflection point stays S_mean_spa, as before.
    return spamix_cgf::twoSidedSpa(cgf, S_mean_spa, absS, S_var_diag, zScore);
}

// ======================================================================
// pushResult — the seven output cells
// ======================================================================
//
// P LOG10P Z Z_Norm BETA SE SPA_STATUS.  Z is the p-consistent z (the normal
// quantile of P carrying Z_Norm's sign) and Z_Norm the raw score z; they differ
// in the tails, which is exactly where the saddlepoint does its work.

void SPAmixPlusMethod::pushResult(
    std::vector<double> &r,
    const spa::TwoSided &ts,
    double zScore,
    double varS
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double sqrtVarS = (varS > 0.0) ? std::sqrt(varS) : 0.0;
    const double zOut = (sqrtVarS > 0.0) ? zScore : nan;
    const double beta = (sqrtVarS > 0.0) ? zScore / sqrtVarS : nan;
    const double se   = (sqrtVarS > 0.0) ? 1.0 / sqrtVarS : nan;

    r.push_back(ts.p);
    r.push_back(ts.negLog10p);
    r.push_back(math::zFromPval(ts.p, zOut));
    r.push_back(zOut);
    r.push_back(beta);
    r.push_back(se);
    r.push_back(spamix_cgf::statusCode(ts.status));
}

// ======================================================================
// getResultVec — MethodBase scalar interface (kept for compatibility)
// ======================================================================

void SPAmixPlusMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int markerInChunkIdx,
    std::vector<double> &result
) {
    fillAFVecForMarker(GVec, altFreq, markerInChunkIdx, m_AFVec);
    m_WVec = 2.0 * m_AFVec.array() * (1.0 - m_AFVec.array());

    double rawScore = GVec.dot(m_resid);

    double zScore, VarS;
    const spa::TwoSided ts = markerPvalFromAF(m_AFVec, m_WVec, rawScore, zScore, VarS);
    pushResult(result, ts, zScore, VarS);
}

// ======================================================================
// ensureAFCacheFused — populate per-thread AFVec/W cache (pre-computed AF)
//
// The cache is shared by all clones with the same maskIdx on this thread.
// First call in a group fills the cache; subsequent calls (other phenos
// in the same group with identical chunkIdxs) reuse it for free.
// ======================================================================

void SPAmixPlusMethod::ensureAFCacheFused(
    const std::vector<int> &chunkIdxs,
    const std::vector<double> &altFreqs
) {
    auto &cache = *m_afCache;
    if (cache.lastChunkIdxs == chunkIdxs) return; // hot path: cache hit

    const int B = static_cast<int>(chunkIdxs.size());
    if (cache.afBatch.rows() != m_N || cache.afBatch.cols() < B) {
        cache.afBatch.resize(m_N, B);
        cache.wBatch.resize(m_N, B);
    }

    // m_chunkGenoIndices is set in prepareChunk; engine guarantees
    // chunkIdxs[b] falls within the current chunk.  Pre-computed AF mode
    // only depends on (model, altFreq) — independent of phenotype.
    for (int b = 0; b < B; ++b) {
        Eigen::Ref<Eigen::VectorXd> col = cache.afBatch.col(b);
        const uint64_t genoIdx = (*m_chunkGenoIndices)[chunkIdxs[b]];
        const uint32_t flatIdx = (*m_genoToFlat)[genoIdx];
        const AFModel &model = (*m_afModels)[flatIdx];
        getAFVecFromModel(model, altFreqs[b], m_onePlusPCs, m_N, col);
    }

    cache.wBatch.leftCols(B).array() =
        2.0 * cache.afBatch.leftCols(B).array() *
        (1.0 - cache.afBatch.leftCols(B).array());

    cache.lastChunkIdxs = chunkIdxs;
}

// ======================================================================
// processScoreBatch — fused-GEMM consumer (pre-computed AF mode)
// ======================================================================

void SPAmixPlusMethod::processScoreBatch(
    const Eigen::Ref<const Eigen::MatrixXd> &scores,
    const double *gSums,
    const double *gSumSqs,
    uint32_t nUsed,
    const std::vector<double> &altFreqs,
    const std::vector<int> &chunkIdxs,
    std::vector<std::vector<double> > &results
) {
    (void)gSums; (void)gSumSqs; (void)nUsed;

    const int B = static_cast<int>(scores.cols());
    results.resize(B);

    // Populate / reuse the per-thread AFVec batch for this maskIdx.
    ensureAFCacheFused(chunkIdxs, altFreqs);
    auto &cache = *m_afCache;

    for (int b = 0; b < B; ++b) {
        auto afCol = cache.afBatch.col(b);
        auto wCol  = cache.wBatch.col(b);

        const double rawScore = scores(0, b);
        double zScore, VarS;
        const spa::TwoSided ts = markerPvalFromAF(afCol, wCol, rawScore, zScore, VarS);

        auto &r = results[b];
        r.clear();
        r.reserve(7);
        pushResult(r, ts, zScore, VarS);
    }
}

// ======================================================================
// getResultBatch — non-fused batched path (on-the-fly AF mode)
//
// Used when supportsFusedGemm() = false (on-the-fly AF).  GBatch is the
// phenotype-local genotype matrix.  We still share AFVec across clones
// in the same maskIdx group via the per-thread cache, so the on-the-fly
// OLS / logistic fit runs once per marker per group instead of once per
// (marker × phenotype).
// ======================================================================

void SPAmixPlusMethod::getResultBatch(
    const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
    const std::vector<double> &altFreqs,
    const std::vector<int> &chunkIdxs,
    std::vector<std::vector<double> > &results
) {
    const int B = static_cast<int>(GBatch.cols());
    results.resize(B);

    if (m_useAFCacheBatch) {
        // On-the-fly AF: each marker's AFVec depends on its own G column, so
        // we cannot piggy-back on the simple genoIdx-keyed fill used in
        // ensureAFCache.  Instead, do an explicit per-marker fill into the
        // shared cache, gated by the cache's chunkIdxs key.
        auto &cache = *m_afCache;
        const bool needFill = (cache.lastChunkIdxs != chunkIdxs);
        if (needFill) {
            if (cache.afBatch.rows() != m_N || cache.afBatch.cols() < B) {
                cache.afBatch.resize(m_N, B);
                cache.wBatch.resize(m_N, B);
            }
            for (int b = 0; b < B; ++b) {
                // Eigen::Ref over the cached column lets fillAFVecForMarker
                // (and computeAFVec/getAFVecFromModel) write in place.
                Eigen::Ref<Eigen::VectorXd> afCol = cache.afBatch.col(b);
                fillAFVecForMarker(GBatch.col(b), altFreqs[b], chunkIdxs[b], afCol);
            }
            cache.wBatch.leftCols(B).array() =
                2.0 * cache.afBatch.leftCols(B).array() *
                (1.0 - cache.afBatch.leftCols(B).array());
            cache.lastChunkIdxs = chunkIdxs;
        }

        for (int b = 0; b < B; ++b) {
            auto afCol = cache.afBatch.col(b);
            auto wCol  = cache.wBatch.col(b);

            // Raw score = resid · G.  Done per phenotype since residuals differ.
            const double rawScore = GBatch.col(b).dot(m_resid);

            double zScore, VarS;
            const spa::TwoSided ts =
                markerPvalFromAF(afCol, wCol, rawScore, zScore, VarS);

            auto &r = results[b];
            r.clear();
            r.reserve(7);
            pushResult(r, ts, zScore, VarS);
        }
        return;
    }

    // Marker-by-marker path: no N × B cache.  Used when this phenotype is
    // the sole occupant of its maskIdx, so the cache would never deliver a
    // cross-call hit anyway.  Reuses the per-clone scratch m_AFVec / m_WVec
    // already sized to m_N at construction.
    for (int b = 0; b < B; ++b) {
        fillAFVecForMarker(GBatch.col(b), altFreqs[b], chunkIdxs[b], m_AFVec);
        m_WVec.array() = 2.0 * m_AFVec.array() * (1.0 - m_AFVec.array());

        const double rawScore = GBatch.col(b).dot(m_resid);

        double zScore, VarS;
        const spa::TwoSided ts =
            markerPvalFromAF(m_AFVec, m_WVec, rawScore, zScore, VarS);

        auto &r = results[b];
        r.clear();
        r.reserve(7);
        pushResult(r, ts, zScore, VarS);
    }
}

// ======================================================================
// fillUnionResiduals / fillResidualSums — fused-GEMM hooks
// ======================================================================

void SPAmixPlusMethod::fillUnionResiduals(
    Eigen::Ref<Eigen::MatrixXd> dest,
    const std::vector<uint32_t> &unionToLocal
) const {
    const uint32_t nUnion = static_cast<uint32_t>(unionToLocal.size());
    for (uint32_t i = 0; i < nUnion; ++i) {
        const uint32_t li = unionToLocal[i];
        if (li != UINT32_MAX) dest(i, 0) = m_resid[li];
    }
}

void SPAmixPlusMethod::fillResidualSums(double *dest) const {
    dest[0] = m_residSum;
}

// ======================================================================
// runSPAmixPlus — top-level orchestration
// ======================================================================

void runSPAmixPlus(
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &pcColNames,
    const std::string &phenoFile,
    const std::string &covarFile,
    const GenoSpec &geno,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &afFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierRatio,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    const std::vector<std::string> &covarNames,
    bool saveResid,
    bool longitudinal
) {
    // --longitudinal injects R_G (random-intercept residual) post-finalize, so
    // the GLM fit-path must not engage.
    const bool fitPath = !longitudinal && !phenoNameSpec.empty();
    nullmodel::RegressionModel regModel{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        regModel = nullmodel::parseRegressionModel(regressionModelStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
        infoMsg("SPAmix: fitting %s null model for %zu phenotype(s)",
                nullmodel::regressionModelName(regModel), phenoSpecs.size());
    }

    // ---- Load residual file and PC data ----
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);

    // --longitudinal: fit Y ~ X + (1|IID) on the long-format pheno file and
    // obtain R_G + the per-subject design before SubjectData is built.  Plain
    // SPAmix has no GRM (SPAmixPlus + --longitudinal is rejected at dispatch),
    // so grmIDs is normally empty.
    nsLongitudinal::LongResidResult Lr;
    if (longitudinal) {
        std::unordered_set<std::string> grmIDs;
        if (!spgrmGrabFile.empty() || !spgrmGctaFile.empty())
            grmIDs = SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, famIIDs);
        const auto outcomeNames = nsLongitudinal::splitOutcomeNames(phenoNameSpec);
        const auto kept = nsLongitudinal::buildKeptSet(keepFile, removeFile, famIIDs, grmIDs);
        infoMsg("SPAmix: fitting longitudinal Y ~ X + (1|IID) for %zu outcome(s)",
                outcomeNames.size());
        Lr = nsLongitudinal::computeLongitudinalResid(
            phenoFile, outcomeNames, covarNames, famIIDs, kept);
    }

    SubjectData sd(std::move(famIIDs));
    if (longitudinal) {
        sd.setKeepSubjects(
            std::unordered_set<std::string>(Lr.usedIIDs.begin(), Lr.usedIIDs.end()));
    } else if (fitPath) {
        // Load PC columns + null-model phenotype columns + null-model
        // covariates (when --covar is absent, --covar-name selects columns
        // from --pheno) in a single pass.  loadPhenoFile rejects duplicate
        // calls, so the column lists must be unioned here.
        std::vector<std::string> wanted = pcColNames;
        auto add = [&](const std::string &name) {
            if (name.empty()) return;
            if (std::find(wanted.begin(), wanted.end(), name) == wanted.end())
                wanted.push_back(name);
        };
        for (const auto &name : nullmodel::columnsNeeded(phenoSpecs)) add(name);
        if (covarFile.empty())
            for (const auto &name : covarNames) add(name);
        sd.loadPhenoFile(phenoFile, wanted);
    } else {
        sd.loadResidOne(phenoFile, residNames);
        if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile);
    }
    if (!covarFile.empty() && !longitudinal) {
        // Load the union of --pc-cols (needed by the per-individual AF model)
        // and --covar-name (needed by the null-model design).  Loading only
        // pcColNames would silently drop any --covar-name columns that are
        // present in the covar file but absent from --pc-cols, producing a
        // "column not found" error at the null-model fit.
        std::vector<std::string> covarLoadCols = pcColNames;
        for (const auto &name : covarNames) {
            if (std::find(covarLoadCols.begin(), covarLoadCols.end(), name) ==
                covarLoadCols.end())
                covarLoadCols.push_back(name);
        }
        sd.loadCovar(covarFile, covarLoadCols);
    }
    sd.setKeepRemove(keepFile, removeFile);
    if (!spgrmGrabFile.empty() || !spgrmGctaFile.empty())sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile,
                                                                                                      spgrmGctaFile, sd.
                                                                                                      famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();

    const int N = static_cast<int>(sd.nUsed());
    const int nPC = static_cast<int>(pcColNames.size());
    infoMsg("  %u subjects in union mask, %d PCs", sd.nUsed(), nPC);

    // ---- Inject longitudinal R_G + extract per-subject PCs ----
    // In longitudinal mode the PCs feeding the per-individual AF model come
    // from the long file's per-subject design (first row per IID).  --pc-cols
    // must be a subset of --covar-name (enforced at dispatch) so every PC is a
    // column of perSubjectX.
    Eigen::MatrixXd longPCs;
    if (longitudinal) {
        const auto used = sd.usedIIDs();
        std::vector<Eigen::VectorXd> rs;
        std::vector<std::string> ns;
        rs.reserve(Lr.R_G.size());
        ns.reserve(Lr.R_G.size());
        for (size_t k = 0; k < Lr.R_G.size(); ++k) {
            rs.push_back(nsLongitudinal::alignToUsed(Lr.usedIIDs, Lr.R_G[k], used));
            ns.push_back(Lr.names[k]);
        }
        sd.setResidualsFromFit(std::move(rs), std::move(ns));

        const Eigen::MatrixXd alignedX =
            nsLongitudinal::alignRowsToUsed(Lr.usedIIDs, Lr.perSubjectX, used);
        longPCs.resize(N, nPC);
        for (int j = 0; j < nPC; ++j) {
            auto it = std::find(Lr.xColNames.begin(), Lr.xColNames.end(), pcColNames[j]);
            if (it == Lr.xColNames.end())
                throw std::runtime_error(
                    "--longitudinal SPAmix: --pc-cols column '" + pcColNames[j] +
                    "' must also appear in --covar-name (PCs are sourced from the "
                    "long-format design)");
            longPCs.col(j) = alignedX.col(it - Lr.xColNames.begin());
        }
    }

    if (fitPath) {
        // --pc-cols and --covar-name are kept strictly separate in SPAmix /
        // SPAmixPlus:
        //   --pc-cols    → per-individual AF model (this file, line ~705)
        //                  and the AF model coefficients loaded/fit below.
        //   --covar-name → null-model design only (this block).
        // Users who want PCs adjusted in the null model must list them in
        // --covar-name explicitly.  This makes the two roles auditable and
        // prevents silent inclusion of PCs that the user only intended for
        // AF estimation.
        Eigen::MatrixXd covarUnion;
        if (covarNames.empty())
            covarUnion.resize(sd.nUsed(), 0);
        else
            covarUnion = sd.getColumns(covarNames);

        if (covarUnion.cols() > 0)
            infoMsg(
                "  Null-model design: intercept + %zu --covar-name covariate(s); "
                "--pc-cols is used only for the individual-AF model",
                covarNames.size()
            );
        else
            infoMsg(
                "  Null-model design: intercept only (no --covar-name supplied; "
                "--pc-cols is not used in the null model)"
            );

        nullmodel::EngineOptions eo;
        eo.nthreads = nthread;
        auto fits = nullmodel::fitAll(sd, phenoSpecs, regModel, covarUnion, eo);
        std::vector<Eigen::VectorXd> rs;
        std::vector<std::string> ns;
        rs.reserve(fits.size());
        ns.reserve(fits.size());
        for (auto &f : fits) {
            infoMsg("  Fitted '%s': %d subjects after NaN removal",
                    f.name.c_str(), f.nUsedRows);
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

    // ---- Design matrix: [1 | PCs] at union dimension ----
    Eigen::MatrixXd unionPCs = longitudinal ? longPCs : sd.getColumns(pcColNames);
    Eigen::MatrixXd unionOnePlusPCs(N, 1 + nPC);
    unionOnePlusPCs.col(0).setOnes();
    unionOnePlusPCs.rightCols(nPC) = unionPCs;

    // ---- Load sparse GRM at union dimension ----
    const bool hasGRM = !spgrmGrabFile.empty() || !spgrmGctaFile.empty();
    std::unique_ptr<SparseGRM> unionGrm;
    if (hasGRM) {
        infoMsg("Loading sparse GRM (raw)...");
        unionGrm =
            std::make_unique<SparseGRM>(SparseGRM::load(spgrmGrabFile, spgrmGctaFile, sd.usedIIDs(), sd.famIIDs()));
        infoMsg("  %u subjects, %zu entries", unionGrm->nSubjects(), unionGrm->nnz());
    }

    // ---- Load genotype data ----
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    const auto &markerInfo = genoData->markerInfo();
    const uint32_t nMarkers = static_cast<uint32_t>(markerInfo.size());

    // ---- genoIndex → flat marker index mapping ----
    std::vector<uint32_t> genoToFlat(genoData->nMarkers(), UINT32_MAX);
    std::vector<uint64_t> genoIndices(nMarkers);
    for (uint32_t fi = 0; fi < nMarkers; ++fi) {
        genoToFlat[markerInfo[fi].genoIndex] = fi;
        genoIndices[fi] = markerInfo[fi].genoIndex;
    }

    // ---- Load AF models (pre-computed) or prepare on-the-fly ----
    std::vector<AFModel> afModels;
    if (!afFile.empty()) {
        infoMsg("Loading pre-computed AF models: %s", afFile.c_str());
        afModels = loadAFModels(afFile, nPC, nMarkers, genoIndices);
        infoMsg("  %zu AF models loaded", afModels.size());
    }

    // ---- Build per-phenotype tasks ----
    auto phenoInfos = sd.buildPerColumnMasks();
    const int K = sd.residOneCols();
    if (K > 1) infoMsg("Multi-column residual file: %d phenotypes", K);

    const char *methodLabel = hasGRM ? "SPAmixP" : "SPAmix";

    // Per-phenotype data storage (must outlive PhenoTasks)
    std::vector<Eigen::VectorXd> pResid(K), pResid2(K);
    std::vector<OutlierData> pOutlier(K);

    // Deduped pools keyed by non-missingness pattern (unionToLocal).
    // Phenotypes sharing the same valid-subject set produce identical
    // PC matrices, OLS matrices, and GRM sub-matrices.
    //
    // CRITICAL: reserve(K) upfront so subsequent push_back calls cannot
    // reallocate the underlying buffer.  Each SPAmixPlusMethod stores a
    // raw reference / pointer into these pools at construction time
    // (m_resid, m_onePlusPCs, m_XtX_inv_Xt, m_sqrt_XtX_inv_diag, m_grm);
    // a single reallocation would invalidate all references handed to
    // previously-constructed methods, producing dangling reads inside
    // markerPvalFromAF (m_grm->m_nSubj returning garbage → spaVariance
    // size mismatch or SIGSEGV) or computeAFVec (m_XtX_inv_Xt deref into
    // freed memory).  K is an upper bound on the pool length, so a
    // single reserve(K) makes the buffers stable for the lifetime of
    // the tasks.
    std::vector<Eigen::MatrixXd> poolOnePlusPCs;
    std::vector<Eigen::MatrixXd> poolXtX_inv_Xt;
    std::vector<Eigen::VectorXd> poolSqrt_XtX_inv_diag;
    std::vector<SparseGRM> poolGrm;
    poolOnePlusPCs.reserve(static_cast<size_t>(K));
    poolXtX_inv_Xt.reserve(static_cast<size_t>(K));
    poolSqrt_XtX_inv_diag.reserve(static_cast<size_t>(K));
    poolGrm.reserve(static_cast<size_t>(K));
    std::vector<size_t> maskIdx(K);   // phenotype rc -> index into pools

    std::vector<PhenoTask> tasks(K);

    for (int rc = 0; rc < K; ++rc) {
        const auto &pi = phenoInfos[rc];

        // Extract per-phenotype residuals (always per-phenotype)
        pResid[rc] = (K > 1) ? extractPhenoVec(sd.residMatrix().col(rc), pi) : sd.residuals();
        pResid2[rc] = pResid[rc].array().square();
        pOutlier[rc] = detectOutliers(pResid[rc], outlierRatio);

        // Cache PC matrix, OLS matrices, and GRM by non-missingness pattern.
        size_t mIdx = poolOnePlusPCs.size(); // default: build new
        if (K > 1) {
            for (int j = 0; j < rc; ++j) {
                if (phenoInfos[j].unionToLocal == pi.unionToLocal) {
                    mIdx = maskIdx[j];
                    infoMsg("  Phenotype '%s': reusing PC/OLS/GRM from '%s'",
                            pi.name.c_str(), phenoInfos[j].name.c_str());
                    break;
                }
            }
        }
        if (mIdx == poolOnePlusPCs.size()) {
            // New non-missingness pattern — build all derived objects
            if (K > 1) {
                poolOnePlusPCs.push_back(extractPhenoMat(unionOnePlusPCs, pi));
            } else {
                poolOnePlusPCs.push_back(unionOnePlusPCs);
            }
            const auto &curPCs = poolOnePlusPCs.back();

            if (afFile.empty()) {
                Eigen::MatrixXd XtX = curPCs.transpose() * curPCs;
                Eigen::MatrixXd XtX_inv = XtX.ldlt().solve(Eigen::MatrixXd::Identity(1 + nPC, 1 + nPC));
                poolXtX_inv_Xt.push_back(XtX_inv * curPCs.transpose());
                poolSqrt_XtX_inv_diag.push_back(XtX_inv.diagonal().cwiseSqrt());
            }

            if (hasGRM && K > 1) {
                const auto &u2l = pi.unionToLocal;
                std::vector<SparseGRM::Entry> pEntries;
                for (const auto &e : unionGrm->entries()) {
                    uint32_t li = u2l[e.row], lj = u2l[e.col];
                    if (li != UINT32_MAX && lj != UINT32_MAX) pEntries.push_back({li, lj, e.value});
                }
                poolGrm.push_back(SparseGRM::fromEntries(pi.nUsed, std::move(pEntries)));
            }
        }
        maskIdx[rc] = mIdx;

        // Resolve pointers into deduped pools
        const auto &curPCs = poolOnePlusPCs[mIdx];

        SparseGRM *grmPtr = nullptr;
        if (hasGRM) {
            if (K > 1) {
                grmPtr = &poolGrm[mIdx];
            } else {
                grmPtr = unionGrm.get();
            }
        }

        // Create method
        std::unique_ptr<SPAmixPlusMethod> m;
        const int maskIdxArg = static_cast<int>(mIdx);
        if (hasGRM) {
            if (!afFile.empty()) {
                m = std::make_unique<SPAmixPlusMethod>(
                    pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                    spaCutoff, *grmPtr, afModels, genoToFlat, maskIdxArg);
            } else {
                m = std::make_unique<SPAmixPlusMethod>(
                    pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                    spaCutoff, *grmPtr,
                    poolXtX_inv_Xt[mIdx], poolSqrt_XtX_inv_diag[mIdx], nPC,
                    maskIdxArg);
            }
        } else {
            if (!afFile.empty()) {
                m = std::make_unique<SPAmixPlusMethod>(
                    pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                    spaCutoff, afModels, genoToFlat, maskIdxArg);
            } else {
                m = std::make_unique<SPAmixPlusMethod>(
                    pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                    spaCutoff,
                    poolXtX_inv_Xt[mIdx], poolSqrt_XtX_inv_diag[mIdx], nPC,
                    maskIdxArg);
            }
        }

        tasks[rc].phenoName = pi.name;
        tasks[rc].method = std::move(m);
        tasks[rc].unionToLocal = pi.unionToLocal;
        tasks[rc].nUsed = pi.nUsed;
        infoMsg("  Phenotype '%s': %u subjects", pi.name.c_str(), pi.nUsed);
    }
    if (K > 1)
        infoMsg("  %zu unique subject mask(s) for %d phenotypes", poolOnePlusPCs.size(), K);

    // Configure the on-the-fly AF cache per phenotype.  The N × B cache
    // amortizes AFVec computation across phenotypes that share a maskIdx
    // and therefore co-occupy a MissBatch.  When a maskIdx is occupied by
    // a single phenotype the cache delivers no cross-call reuse and only
    // inflates per-thread memory by ~2 · N · B · 8 bytes per (thread,
    // maskIdx).  We pass the share-count down so each method can skip the
    // cache when it would never be hit.
    std::vector<int> maskGroupSize(poolOnePlusPCs.size(), 0);
    for (int rc = 0; rc < K; ++rc) ++maskGroupSize[maskIdx[rc]];
    for (int rc = 0; rc < K; ++rc) {
        auto *spam = dynamic_cast<SPAmixPlusMethod *>(tasks[rc].method.get());
        if (spam) spam->setUseAFCacheBatch(maskGroupSize[maskIdx[rc]] >= 2);
    }

    infoMsg("Running %s marker tests (%d thread(s), %d phenotype(s))...", methodLabel, nthread, K);
    multiPhenoEngine(
        *genoData,
        tasks,
        outPrefix,
        methodLabel,
        compression,
        compressionLevel,
        nthread,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        hweCutoff
    );
}

// ══════════════════════════════════════════════════════════════════════
// runSPAmixPlusLoco — LOCO orchestration
// ══════════════════════════════════════════════════════════════════════
//
// The per-individual AF design (--pc-cols), the deduped OLS/GRM pools, and the
// AF-model cache are residual-independent and are built once.  Per chromosome
// the null model is refit with that chromosome's LOCO PGS appended as a
// covariate column, only the residual pools change, and each SPAmixPlusMethod is
// reconstructed.  Because SPAmixPlusMethod stores const references into the
// residual pools, the previous chromosome's methods are destroyed before the
// residual pools are overwritten.

void runSPAmixPlusLoco(
    const std::vector<std::string> &pcColNames,
    const std::string &phenoFile,
    const std::string &covarFile,
    const GenoSpec &geno,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &afFile,
    const std::string &predListFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierRatio,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec,
    const std::vector<std::string> &covarNames
) {
    if (phenoNameSpec.empty())
        throw std::runtime_error(
            "SPAmix-LOCO requires --pheno-name (an in-process null-model fit); "
            "precomputed residuals cannot be refit per chromosome");

    nullmodel::RegressionModel regModel =
        nullmodel::parseRegressionModel(regressionModelStr);
    std::vector<nullmodel::PhenoSpec> phenoSpecs =
        nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
    const int K = static_cast<int>(phenoSpecs.size());
    std::vector<std::string> specNames(K);
    for (int k = 0; k < K; ++k) specNames[k] = phenoSpecs[k].name;

    const bool hasGRM = !spgrmGrabFile.empty() || !spgrmGctaFile.empty();
    const char *methodLabel = hasGRM ? "SPAmixP" : "SPAmix";
    infoMsg("%s-LOCO: fitting %s null model for %d phenotype(s), "
            "per-chromosome LOCO PGS as covariate",
            methodLabel, nullmodel::regressionModelName(regModel), K);

    validatePredListPhenos(predListFile, specNames);

    // ---- Load pheno/PC/covar data (fit path only) ----
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    {
        std::vector<std::string> wanted = pcColNames;
        auto add = [&](const std::string &name) {
            if (name.empty()) return;
            if (std::find(wanted.begin(), wanted.end(), name) == wanted.end())
                wanted.push_back(name);
        };
        for (const auto &name : nullmodel::columnsNeeded(phenoSpecs)) add(name);
        if (covarFile.empty())
            for (const auto &name : covarNames) add(name);
        sd.loadPhenoFile(phenoFile, wanted);
    }
    if (!covarFile.empty()) {
        std::vector<std::string> covarLoadCols = pcColNames;
        for (const auto &name : covarNames)
            if (std::find(covarLoadCols.begin(), covarLoadCols.end(), name) ==
                covarLoadCols.end())
                covarLoadCols.push_back(name);
        sd.loadCovar(covarFile, covarLoadCols);
    }
    sd.setKeepRemove(keepFile, removeFile);
    if (hasGRM)
        sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();

    const int N = static_cast<int>(sd.nUsed());
    const int nPC = static_cast<int>(pcColNames.size());
    infoMsg("  %u subjects in union mask, %d PCs", sd.nUsed(), nPC);

    // ---- Base covariate matrix (--covar-name only, no PGS) ----
    Eigen::MatrixXd baseCovar;
    if (covarNames.empty()) baseCovar.resize(N, 0);
    else baseCovar = sd.getColumns(covarNames);

    // ---- One base fit to establish per-phenotype masks ----
    {
        nullmodel::EngineOptions eo;
        eo.nthreads = nthread;
        auto fits = nullmodel::fitAll(sd, phenoSpecs, regModel, baseCovar, eo);
        std::vector<Eigen::VectorXd> rs;
        std::vector<std::string> ns;
        rs.reserve(fits.size());
        ns.reserve(fits.size());
        for (auto &f : fits) {
            infoMsg("  Fitted '%s': %d subjects after NaN removal",
                    f.name.c_str(), f.nUsedRows);
            rs.push_back(std::move(f.residuals));
            ns.push_back(f.name);
        }
        sd.setResidualsFromFit(std::move(rs), std::move(ns));
    }
    auto phenoInfos = sd.buildPerColumnMasks();

    // Resolve each spec to a concrete model once (chromosome-invariant).
    std::vector<nullmodel::RegressionModel> specModel(K);
    for (int k = 0; k < K; ++k) {
        if (regModel != nullmodel::RegressionModel::Auto)
            specModel[k] = regModel;
        else if (nullmodel::isCoxSpec(phenoSpecs[k]))
            specModel[k] = nullmodel::RegressionModel::Cox;
        else
            specModel[k] = nullmodel::inferModelFromColumn(
                sd.getColumn(phenoSpecs[k].yColumn), phenoSpecs[k].yColumn,
                sd.usedIIDs()).model;
    }

    // ---- Individual-AF design [1 | PCs] (residual-independent) ----
    Eigen::MatrixXd unionPCs = sd.getColumns(pcColNames);
    Eigen::MatrixXd unionOnePlusPCs(N, 1 + nPC);
    unionOnePlusPCs.col(0).setOnes();
    unionOnePlusPCs.rightCols(nPC) = unionPCs;

    // ---- GRM / genotype / AF models (once) ----
    std::unique_ptr<SparseGRM> unionGrm;
    if (hasGRM) {
        infoMsg("Loading sparse GRM (raw)...");
        unionGrm = std::make_unique<SparseGRM>(
            SparseGRM::load(spgrmGrabFile, spgrmGctaFile, sd.usedIIDs(), sd.famIIDs()));
        infoMsg("  %u subjects, %zu entries", unionGrm->nSubjects(), unionGrm->nnz());
    }
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);
    const auto &markerInfo = genoData->markerInfo();
    const uint32_t nMarkers = static_cast<uint32_t>(markerInfo.size());
    std::vector<uint32_t> genoToFlat(genoData->nMarkers(), UINT32_MAX);
    std::vector<uint64_t> genoIndices(nMarkers);
    for (uint32_t fi = 0; fi < nMarkers; ++fi) {
        genoToFlat[markerInfo[fi].genoIndex] = fi;
        genoIndices[fi] = markerInfo[fi].genoIndex;
    }
    std::vector<AFModel> afModels;
    if (!afFile.empty()) {
        infoMsg("Loading pre-computed AF models: %s", afFile.c_str());
        afModels = loadAFModels(afFile, nPC, nMarkers, genoIndices);
        infoMsg("  %zu AF models loaded", afModels.size());
    }

    // ---- Deduped PC/OLS/GRM pools (residual-independent, built ONCE) ----
    std::vector<Eigen::MatrixXd> poolOnePlusPCs;
    std::vector<Eigen::MatrixXd> poolXtX_inv_Xt;
    std::vector<Eigen::VectorXd> poolSqrt_XtX_inv_diag;
    std::vector<SparseGRM> poolGrm;
    poolOnePlusPCs.reserve(static_cast<size_t>(K));
    poolXtX_inv_Xt.reserve(static_cast<size_t>(K));
    poolSqrt_XtX_inv_diag.reserve(static_cast<size_t>(K));
    poolGrm.reserve(static_cast<size_t>(K));
    std::vector<size_t> maskIdx(K);
    for (int rc = 0; rc < K; ++rc) {
        const auto &pi = phenoInfos[rc];
        size_t mIdx = poolOnePlusPCs.size();
        if (K > 1) {
            for (int j = 0; j < rc; ++j) {
                if (phenoInfos[j].unionToLocal == pi.unionToLocal) {
                    mIdx = maskIdx[j];
                    break;
                }
            }
        }
        if (mIdx == poolOnePlusPCs.size()) {
            if (K > 1) poolOnePlusPCs.push_back(extractPhenoMat(unionOnePlusPCs, pi));
            else poolOnePlusPCs.push_back(unionOnePlusPCs);
            const auto &curPCs = poolOnePlusPCs.back();
            if (afFile.empty()) {
                Eigen::MatrixXd XtX = curPCs.transpose() * curPCs;
                Eigen::MatrixXd XtX_inv =
                    XtX.ldlt().solve(Eigen::MatrixXd::Identity(1 + nPC, 1 + nPC));
                poolXtX_inv_Xt.push_back(XtX_inv * curPCs.transpose());
                poolSqrt_XtX_inv_diag.push_back(XtX_inv.diagonal().cwiseSqrt());
            }
            if (hasGRM && K > 1) {
                const auto &u2l = pi.unionToLocal;
                std::vector<SparseGRM::Entry> pEntries;
                for (const auto &e : unionGrm->entries()) {
                    uint32_t li = u2l[e.row], lj = u2l[e.col];
                    if (li != UINT32_MAX && lj != UINT32_MAX)
                        pEntries.push_back({li, lj, e.value});
                }
                poolGrm.push_back(SparseGRM::fromEntries(pi.nUsed, std::move(pEntries)));
            }
        }
        maskIdx[rc] = mIdx;
    }
    std::vector<int> maskGroupSize(poolOnePlusPCs.size(), 0);
    for (int rc = 0; rc < K; ++rc) ++maskGroupSize[maskIdx[rc]];

    // ---- Per-phenotype residual pools (overwritten per chromosome; methods
    //      hold const references into these) ----
    std::vector<Eigen::VectorXd> pResid(K), pResid2(K);
    std::vector<OutlierData> pOutlier(K);

    // ---- Non-missing masks + LOCO predictions ----
    std::vector<std::vector<bool> > nonMissing(K, std::vector<bool>(N, false));
    for (int k = 0; k < K; ++k)
        for (int i = 0; i < N; ++i)
            nonMissing[k][static_cast<size_t>(i)] =
                (phenoInfos[k].unionToLocal[static_cast<size_t>(i)] != UINT32_MAX);

    LocoData loco = LocoData::load(predListFile, specNames, sd.usedIIDs(), sd.famIIDs());
    auto locoChroms = loco.availableChromosomes();
    infoMsg("LOCO: %zu chromosomes available across all phenotypes", locoChroms.size());

    nullmodel::EngineOptions eo1;
    eo1.nthreads = 1;

    auto buildTasks = [&](const std::string &chr, std::vector<PhenoTask> &tasks) {
        tasks.resize(K);
        // Destroy previous chromosome's methods before overwriting the residual
        // pools they reference (lifetime discipline).
        for (auto &t : tasks) t.method.reset();

        for (int rc = 0; rc < K; ++rc) {
            const auto &pi = phenoInfos[rc];
            Eigen::MatrixXd covarUnion_k = appendLocoCovariate(
                loco, specNames[rc], chr, baseCovar, nonMissing[rc], "SPAmix-LOCO");
            auto fits1 = nullmodel::fitAll(
                sd, {phenoSpecs[rc]}, specModel[rc], covarUnion_k, eo1);
            pResid[rc] = extractPhenoVec(fits1[0].residuals, pi);
            pResid2[rc] = pResid[rc].array().square();
            pOutlier[rc] = detectOutliers(pResid[rc], outlierRatio);

            const size_t mIdx = maskIdx[rc];
            const auto &curPCs = poolOnePlusPCs[mIdx];
            SparseGRM *grmPtr = nullptr;
            if (hasGRM) grmPtr = (K > 1) ? &poolGrm[mIdx] : unionGrm.get();

            std::unique_ptr<SPAmixPlusMethod> m;
            const int maskIdxArg = static_cast<int>(mIdx);
            if (hasGRM) {
                if (!afFile.empty())
                    m = std::make_unique<SPAmixPlusMethod>(
                        pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                        spaCutoff, *grmPtr, afModels, genoToFlat, maskIdxArg);
                else
                    m = std::make_unique<SPAmixPlusMethod>(
                        pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                        spaCutoff, *grmPtr,
                        poolXtX_inv_Xt[mIdx], poolSqrt_XtX_inv_diag[mIdx], nPC,
                        maskIdxArg);
            } else {
                if (!afFile.empty())
                    m = std::make_unique<SPAmixPlusMethod>(
                        pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                        spaCutoff, afModels, genoToFlat, maskIdxArg);
                else
                    m = std::make_unique<SPAmixPlusMethod>(
                        pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                        spaCutoff,
                        poolXtX_inv_Xt[mIdx], poolSqrt_XtX_inv_diag[mIdx], nPC,
                        maskIdxArg);
            }
            m->setUseAFCacheBatch(maskGroupSize[mIdx] >= 2);
            tasks[rc].phenoName = pi.name;
            tasks[rc].method = std::move(m);
            tasks[rc].unionToLocal = pi.unionToLocal;
            tasks[rc].nUsed = pi.nUsed;
        }
    };

    infoMsg("%s-LOCO: starting LOCO association (%d phenotype(s), %zu chroms, %d threads)",
            methodLabel, K, locoChroms.size(), nthread);
    locoEngine(
        *genoData, locoChroms, specNames, buildTasks, outPrefix, methodLabel,
        compression, compressionLevel, nthread,
        missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);
}
