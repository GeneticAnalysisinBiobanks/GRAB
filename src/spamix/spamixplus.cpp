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

#include "geno_factory/geno_data.hpp"
#include "spamix/indiv_af.hpp"
#include "io/subject_data.hpp"
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
      m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size())),
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
      m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size())),
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
      m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size())),
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
      m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size())),
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
// ======================================================================

double SPAmixPlusMethod::markerPvalFromAF(
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
        return 1.0;
    }

    zScore = (rawScore - S_mean) / std::sqrt(VarS);
    outVarS = VarS;
    if (std::abs(zScore) < m_spaCutoff) {
        return 2.0 * math::pnorm(-std::abs(zScore));
    }

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

    const int nOut = static_cast<int>(m_outlier.posOutlier.size());
    const int nNon = static_cast<int>(m_outlier.posNonOutlier.size());

    for (int i = 0; i < nOut; ++i)
        m_mafOutlier[i] = afVec[m_outlier.posOutlier[i]];

    double mean_nonOutlier = 0.0, var_nonOutlier = 0.0;
    for (int i = 0; i < nNon; ++i) {
        const double af = afVec[m_outlier.posNonOutlier[i]];
        m_mafNonOutlier[i] = af;
        mean_nonOutlier += m_outlier.residNonOutlier[i] * af;
        var_nonOutlier += m_outlier.resid2NonOutlier[i] * af * (1.0 - af);
    }
    mean_nonOutlier *= 2.0;
    var_nonOutlier *= 2.0;

    double absS = std::abs(S_spa - S_mean_spa);

    double pval1 = spa::getProbSpaG(
        m_mafOutlier.data(),
        m_outlier.residOutlier.data(),
        nOut,
        absS + S_mean_spa,
        false,
        mean_nonOutlier,
        var_nonOutlier
    );

    double pval2 = spa::getProbSpaG(
        m_mafOutlier.data(),
        m_outlier.residOutlier.data(),
        nOut,
        -absS + S_mean_spa,
        true,
        mean_nonOutlier,
        var_nonOutlier
    );

    return pval1 + pval2;
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
    double pval = markerPvalFromAF(m_AFVec, m_WVec, rawScore, zScore, VarS);
    double sqrtVarS = (VarS > 0.0) ? std::sqrt(VarS) : 0.0;
    double beta = (sqrtVarS > 0.0) ? zScore / sqrtVarS : 0.0;
    double se   = (sqrtVarS > 0.0) ? 1.0 / sqrtVarS : 0.0;
    result.push_back(pval);
    result.push_back(beta);
    result.push_back(se);
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
        double pval = markerPvalFromAF(afCol, wCol, rawScore, zScore, VarS);
        double sqrtVarS = (VarS > 0.0) ? std::sqrt(VarS) : 0.0;
        double beta = (sqrtVarS > 0.0) ? zScore / sqrtVarS : 0.0;
        double se   = (sqrtVarS > 0.0) ? 1.0 / sqrtVarS : 0.0;

        auto &r = results[b];
        r.clear();
        r.reserve(3);
        r.push_back(pval);
        r.push_back(beta);
        r.push_back(se);
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
        double pval = markerPvalFromAF(afCol, wCol, rawScore, zScore, VarS);
        double sqrtVarS = (VarS > 0.0) ? std::sqrt(VarS) : 0.0;
        double beta = (sqrtVarS > 0.0) ? zScore / sqrtVarS : 0.0;
        double se   = (sqrtVarS > 0.0) ? 1.0 / sqrtVarS : 0.0;

        auto &r = results[b];
        r.clear();
        r.reserve(3);
        r.push_back(pval);
        r.push_back(beta);
        r.push_back(se);
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
    const std::string &traitTypeStr,
    const std::string &phenoNameSpec,
    const std::vector<std::string> &covarNames,
    bool saveResid
) {
    const bool fitPath = !phenoNameSpec.empty();
    nullmodel::TraitType traitT{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        traitT = nullmodel::parseTraitType(traitTypeStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(traitT, phenoNameSpec);
        infoMsg("SPAmix: fitting %s null model for %zu phenotype(s)",
                nullmodel::traitTypeName(traitT), phenoSpecs.size());
    }

    // ---- Load residual file and PC data ----
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (fitPath) {
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
    if (!covarFile.empty()) sd.loadCovar(covarFile, pcColNames);
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

    if (fitPath) {
        Eigen::MatrixXd covarUnion;
        if (!covarNames.empty()) {
            covarUnion = sd.getColumns(covarNames);
        } else {
            covarUnion.resize(sd.nUsed(), 0);
        }
        nullmodel::EngineOptions eo;
        eo.nthreads = nthread;
        auto fits = nullmodel::fitAll(sd, phenoSpecs, traitT, covarUnion, eo);
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
            nullmodel::writeResidualsFile(outPrefix + ".null.resid", sd, dumpFits);
        }
        sd.setResidualsFromFit(std::move(rs), std::move(ns));
    }

    // ---- Design matrix: [1 | PCs] at union dimension ----
    Eigen::MatrixXd unionPCs = sd.getColumns(pcColNames);
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
    std::vector<Eigen::MatrixXd> poolOnePlusPCs;
    std::vector<Eigen::MatrixXd> poolXtX_inv_Xt;
    std::vector<Eigen::VectorXd> poolSqrt_XtX_inv_diag;
    std::vector<SparseGRM> poolGrm;
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
