// spamixplus.cpp — Unified SPAmix / SPAmixPlus method implementation
//
// When a sparse GRM is provided, variance uses GRM-based covariance and
// the SPA tail is applied to variance-ratio-corrected S (SPAmixPlus).
// Without a GRM the diagonal variance Σ resid²·2·AF·(1−AF) is used (SPAmix).

#include "spamix/spamixplus.hpp"

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
#include <vector>

#include <zlib.h>

#include "geno_factory/geno_data.hpp"
#include "spamix/indiv_af.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

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
    const std::vector<uint32_t> &genoToFlat
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
      m_afModels(&afModels),
      m_genoToFlat(&genoToFlat),
      m_XtX_inv_Xt(nullptr),
      m_sqrt_XtX_inv_diag(nullptr),
      m_AFVec(m_N),
      m_R_new(m_N),
      m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
      m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size()))
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
    int nPC
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
      m_afModels(nullptr),
      m_genoToFlat(nullptr),
      m_XtX_inv_Xt(&XtX_inv_Xt),
      m_sqrt_XtX_inv_diag(&sqrt_XtX_inv_diag),
      m_AFVec(m_N),
      m_R_new(m_N),
      m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
      m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size()))
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
    const std::vector<uint32_t> &genoToFlat
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
      m_afModels(&afModels),
      m_genoToFlat(&genoToFlat),
      m_XtX_inv_Xt(nullptr),
      m_sqrt_XtX_inv_diag(nullptr),
      m_AFVec(m_N),
      m_R_new(0),
      m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
      m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size()))
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
    int nPC
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
      m_afModels(nullptr),
      m_genoToFlat(nullptr),
      m_XtX_inv_Xt(&XtX_inv_Xt),
      m_sqrt_XtX_inv_diag(&sqrt_XtX_inv_diag),
      m_AFVec(m_N),
      m_R_new(0),
      m_mafOutlier(static_cast<int>(outlier.posOutlier.size())),
      m_mafNonOutlier(static_cast<int>(outlier.posNonOutlier.size()))
{
}

std::unique_ptr<MethodBase> SPAmixPlusMethod::clone() const {
    if (m_hasGRM) {
        if (m_XtX_inv_Xt) {
            return std::make_unique<SPAmixPlusMethod>(
                m_resid,
                m_resid2,
                m_onePlusPCs,
                m_outlier,
                m_spaCutoff,
                *m_grm,
                *m_XtX_inv_Xt,
                *m_sqrt_XtX_inv_diag,
                m_nPC
            );
        }
        return std::make_unique<SPAmixPlusMethod>(
            m_resid,
            m_resid2,
            m_onePlusPCs,
            m_outlier,
            m_spaCutoff,
            *m_grm,
            *m_afModels,
            *m_genoToFlat
        );
    }
    // No GRM
    if (m_XtX_inv_Xt) {
        return std::make_unique<SPAmixPlusMethod>(
            m_resid,
            m_resid2,
            m_onePlusPCs,
            m_outlier,
            m_spaCutoff,
            *m_XtX_inv_Xt,
            *m_sqrt_XtX_inv_diag,
            m_nPC
        );
    }
    return std::make_unique<SPAmixPlusMethod>(
        m_resid,
        m_resid2,
        m_onePlusPCs,
        m_outlier,
        m_spaCutoff,
        *m_afModels,
        *m_genoToFlat
    );
}

void SPAmixPlusMethod::prepareChunk(const std::vector<uint64_t> &gIndices) {
    m_chunkGenoIndices = gIndices;
}

// ======================================================================
// getMarkerPval — score test with optional GRM variance + SPA
// ======================================================================

double SPAmixPlusMethod::getMarkerPval(
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    double &zScore,
    double &outVarS
) {

    // m_AFVec is pre-filled by caller (on-the-fly or from stored model)

    // Build score statistic S, mean, and diagonal variance in a single pass.
    // When using GRM we also build R_new for the GRM-based variance.
    double S = 0.0, S_mean = 0.0, S_var_diag = 0.0;
    for (int i = 0; i < m_N; ++i) {
        double af = m_AFVec[i];
        double gvar = 2.0 * af * (1.0 - af);
        S += GVec[i] * m_resid[i];
        S_mean += m_resid[i] * af;
        S_var_diag += m_resid2[i] * gvar;
        if (m_hasGRM) m_R_new[i] = m_resid[i] * std::sqrt(gvar);
    }
    S_mean *= 2.0;

    // Variance: GRM-based (SPAmixPlus) or diagonal (SPAmix)
    double VarS;
    if (m_hasGRM) {
        VarS = m_grm->spaVariance(m_R_new.data(), static_cast<uint32_t>(m_N));
    } else {
        VarS = S_var_diag;
    }

    if (VarS <= 0.0) {
        zScore = 0.0;
        outVarS = 0.0;
        return 1.0;
    }

    zScore = (S - S_mean) / std::sqrt(VarS);

    // Normal approximation when |z| < cutoff
    outVarS = VarS;
    if (std::abs(zScore) < m_spaCutoff) return 2.0 * math::pnorm(-std::abs(zScore));

    // ---- SPA with outlier / non-outlier split ----

    // When using GRM, apply variance-ratio correction to S so the SPA
    // operates on the diagonal-variance scale (needed by getProbSpaG).
    double S_spa, S_mean_spa;
    if (m_hasGRM) {
        double sqrt_ratio = std::sqrt(S_var_diag / VarS);
        S_spa = S * sqrt_ratio;
        S_mean_spa = S_mean * sqrt_ratio;
    } else {
        S_spa = S;
        S_mean_spa = S_mean;
    }

    const int nOut = static_cast<int>(m_outlier.posOutlier.size());
    const int nNon = static_cast<int>(m_outlier.posNonOutlier.size());

    // Gather outlier MAFs
    for (int i = 0; i < nOut; ++i)
        m_mafOutlier[i] = m_AFVec[m_outlier.posOutlier[i]];

    // Gather non-outlier MAFs and accumulate normal-approx terms in one pass
    double mean_nonOutlier = 0.0, var_nonOutlier = 0.0;
    for (int i = 0; i < nNon; ++i) {
        const double af = m_AFVec[m_outlier.posNonOutlier[i]];
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
// getResultVec — MethodBase interface (handles flip internally)
// ======================================================================

void SPAmixPlusMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int markerInChunkIdx,
    std::vector<double> &result
) {

    // GVec counts bim col5 (ALT) alleles; altFreq = freq(ALT).
    // No flip needed — plink always returns ALT counts.
    double maf = altFreq;

    // Fill m_AFVec: on-the-fly from genotype, or from pre-computed model
    if (m_XtX_inv_Xt) {
        AFContext ctx{m_onePlusPCs, *m_XtX_inv_Xt, *m_sqrt_XtX_inv_diag, m_onePlusPCs.rightCols(m_nPC), m_N, m_nPC};
        computeAFVec(GVec, maf, ctx, m_AFVec);
    } else {
        const uint64_t genoIdx = m_chunkGenoIndices[markerInChunkIdx];
        const uint32_t flatIdx = (*m_genoToFlat)[genoIdx];
        const AFModel &model = (*m_afModels)[flatIdx];
        getAFVecFromModel(model, maf, m_onePlusPCs, m_N, m_AFVec);
    }

    double zScore, VarS;
    double pval = getMarkerPval(GVec, maf, zScore, VarS);
    double sqrtVarS = (VarS > 0.0) ? std::sqrt(VarS) : 0.0;
    double beta = (sqrtVarS > 0.0) ? zScore / sqrtVarS : 0.0;
    result.push_back(pval);
    result.push_back(zScore);
    result.push_back(beta);
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
    const std::string &removeFile
) {
    // ---- Load residual file and PC data ----
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    sd.loadResidOne(phenoFile, residNames);
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile);
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
        if (hasGRM) {
            if (!afFile.empty()) {
                m = std::make_unique<SPAmixPlusMethod>(
                    pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                    spaCutoff, *grmPtr, afModels, genoToFlat);
            } else {
                m = std::make_unique<SPAmixPlusMethod>(
                    pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                    spaCutoff, *grmPtr,
                    poolXtX_inv_Xt[mIdx], poolSqrt_XtX_inv_diag[mIdx], nPC);
            }
        } else {
            if (!afFile.empty()) {
                m = std::make_unique<SPAmixPlusMethod>(
                    pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                    spaCutoff, afModels, genoToFlat);
            } else {
                m = std::make_unique<SPAmixPlusMethod>(
                    pResid[rc], pResid2[rc], curPCs, pOutlier[rc],
                    spaCutoff,
                    poolXtX_inv_Xt[mIdx], poolSqrt_XtX_inv_diag[mIdx], nPC);
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
