// spagxe.cpp — SPAGxE_CCT G×E test implementation.
//
// Phase 1: base SPAGxE, SPA-only, unrelated population.  Per phenotype a single
// genotype-independent residual R is fit once (trait ~ X + E).  Per variant:
//   1. marginal genetic screen  S_G = Σ G_i R_i ,  Var = 2q(1−q) Σ R_i²
//      (always the normal approximation; routes Branch A vs Branch B).
//   2a. Branch A  (p_marg > ε, the common case): the λ-orthogonalised G×E score
//       S_GxE = Σ G_i w_i ,  w_i = (E_i − λ) R_i ,  λ = Σ E_i R_i² / Σ R_i² ,
//       tested by the SPA / normal hybrid.  w and its outlier partition are
//       marker-independent → precomputed once.
//   2b. Branch B  (p_marg ≤ ε, rare): project the marginal effect out of the
//       residual, R̃ = R − [1,G]([1,G]ᵀ[1,G])⁻¹[1,G]ᵀ R, and test
//       Ŝ_GxE = Σ G_i E_i R̃_i by SPA on the per-marker weight E∘R̃.
//
// The SPA (spa::getProbSpaG) uses the IQR outlier / Gaussian non-outlier split
// (util/outlier.hpp), so it deliberately differs from the reference R package's
// full-sum SPA and its p=0-on-nonconvergence; a degenerate saddlepoint yields
// NaN here, never 0 (GRAB2 convention).  Reference math:
// tmp/SPAGxECCT/.../R/SPAGxECCT.R  (SPAGxE_CCT_one_SNP, SPA_G_one_SNP_homo_new).

#include "spagxe/spagxe.hpp"

#include "geno_factory/geno_data.hpp" // parseGenoIIDs, makeGenoData, GenoMeta
#include "io/subject_data.hpp"        // SubjectData, extractPhenoVec, PerPhenoInfo
#include "spamix/common.hpp"          // spa::getProbSpaG
#include "util/logging.hpp"           // infoMsg, warnMsg
#include "util/math_helper.hpp"       // math::pnorm, math::zFromPval
#include "util/null_model.hpp"        // nullmodel::fitAll and friends

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Two-sided SPA p-value for the genotype score  S = Σ G_i w_i  under the
// homogeneous binomial law  G_i ~ Bin(2, q), reusing spa::getProbSpaG with a
// precomputed outlier / non-outlier partition of the weight vector w.  Mirrors
// the reference SPA_G_one_SNP_homo_new: the normal approximation when
// |z| ≤ spaCutoff, otherwise the reflection-about-mean two-sided saddlepoint.
//   part     — IQR partition of w  (residOutlier = w[outlier],
//              resid2NonOutlier = w²[nonOutlier])
//   wSum     — Σ w  (over all i)          wSq — Σ w²  (over all i)
//   rawScore — Σ G_i w_i
// Sets zScore (raw score z) and outVar = Var(S); returns the p-value (NaN when
// the variance is non-positive or the saddlepoint fails to converge).
double spaScorePval(
    const OutlierData &part,
    double wSum,
    double wSq,
    double rawScore,
    double q,
    double spaCutoff,
    double &zScore,
    double &outVar
) {
    const double var = 2.0 * q * (1.0 - q) * wSq;
    if (!(var > 0.0)) { // monomorphic / degenerate weight
        zScore = 0.0;
        outVar = 0.0;
        return kNaN;
    }
    const double sMean = 2.0 * q * wSum; // retrospective mean E[S] = 2q·Σw
    zScore = (rawScore - sMean) / std::sqrt(var);
    outVar = var;
    if (std::abs(zScore) <= spaCutoff)
        return 2.0 * math::pnorm(-std::abs(zScore));

    // ── Saddlepoint tail with outlier / non-outlier split ──────────────
    const int nOut = static_cast<int>(part.posOutlier.size());
    const double meanNon = 2.0 * q * part.residNonOutlier.sum();
    const double varNon = 2.0 * q * (1.0 - q) * part.resid2NonOutlier.sum();
    std::vector<double> mafOut(static_cast<size_t>(nOut), q);
    const double absS = std::abs(rawScore - sMean);
    // Reflect about the fitted mean: upper = max(S, 2·sMean−S),
    // lower = min(S, 2·sMean−S)  (== sMean ± absS).
    const double pUpper = spa::getProbSpaG(mafOut.data(), part.residOutlier.data(),
                                           nOut, sMean + absS, false, meanNon, varNon);
    const double pLower = spa::getProbSpaG(mafOut.data(), part.residOutlier.data(),
                                           nOut, sMean - absS, true, meanNon, varNon);
    return pUpper + pLower;
}

} // namespace

// ══════════════════════════════════════════════════════════════════════
// SPAGxEMethod
// ══════════════════════════════════════════════════════════════════════

SPAGxEMethod::SPAGxEMethod(
    Eigen::VectorXd resid,
    std::vector<std::string> envNames,
    std::vector<Eigen::VectorXd> envVecs,
    double marginalCutoff,
    double spaCutoff,
    double outlierRatio
)
    : m_resid(std::move(resid)),
      m_envNames(std::move(envNames)),
      m_marginalCutoff(marginalCutoff),
      m_spaCutoff(spaCutoff),
      m_outlierRatio(outlierRatio)
{
    const Eigen::Index N = m_resid.size();
    m_residSum = m_resid.sum();
    m_residSq = m_resid.squaredNorm();
    const int nEnv = static_cast<int>(m_envNames.size());

    // m_W = [ R | w_1 | … | w_nEnv ]; one GEMM against a genotype batch then
    // yields the marginal score (col 0) and every Branch-A G×E score.
    m_W.resize(N, 1 + nEnv);
    m_W.col(0) = m_resid;

    const Eigen::ArrayXd r2 = m_resid.array().square();
    m_envs.reserve(nEnv);
    for (int e = 0; e < nEnv; ++e) {
        EnvData ed;
        ed.name = m_envNames[e];
        ed.E = std::move(envVecs[e]);
        // λ = Σ E R² / Σ R²  (variance-weighted regression of E on the constant).
        ed.lambda = (m_residSq > 0.0)
                        ? (ed.E.array() * r2).sum() / m_residSq
                        : 0.0;
        ed.w = (ed.E.array() - ed.lambda) * m_resid.array();
        ed.wSum = ed.w.sum();
        ed.wSq = ed.w.squaredNorm();
        ed.wOutlier = detectOutliers(ed.w, m_outlierRatio);
        m_W.col(1 + e) = ed.w;
        m_envs.push_back(std::move(ed));
    }
}

std::unique_ptr<MethodBase> SPAGxEMethod::clone() const {
    return std::make_unique<SPAGxEMethod>(*this);
}

std::string SPAGxEMethod::getHeaderColumns() const {
    // Leading marginal block P_G Z_G BETA_G SE_G (normal-approx, so its Z is
    // already p-consistent and needs no Z_Norm), then a 5-wide G×E block per
    // environment (SPA score test → p-consistent Z and raw-score Z_Norm differ
    // in the tails).
    std::string h = "\tP_G\tZ_G\tBETA_G\tSE_G";
    for (const auto &n : m_envNames)
        h += "\tP_Gx" + n + "\tZ_Gx" + n + "\tZ_Norm_Gx" + n + "\tBETA_Gx" + n +
             "\tSE_Gx" + n;
    return h;
}

void SPAGxEMethod::evalMarker(
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    const Eigen::VectorXd &rawScores,
    std::vector<double> &out
) {
    const int nEnv = static_cast<int>(m_envNames.size());
    out.assign(static_cast<size_t>(4 + 5 * nEnv), kNaN);

    const double q = altFreq;

    // ── Marginal genetic block (always the normal approximation) ────────
    // S_G = Σ G_i R_i (uncentered, matching the reference S1; E[S_G]=2q·ΣR≈0).
    const double varSG = m_residSq * 2.0 * q * (1.0 - q);
    double pMarg = 1.0; // degenerate variance ⇒ take Branch A (no projection)
    if (varSG > 0.0) {
        const double sG = rawScores[0];
        const double sdG = std::sqrt(varSG);
        const double zG = sG / sdG;
        pMarg = 2.0 * math::pnorm(-std::abs(zG));
        out[0] = pMarg;       // P_G
        out[1] = zG;          // Z_G (normal ⇒ p-consistent)
        out[2] = sG / varSG;  // BETA_G
        out[3] = 1.0 / sdG;   // SE_G
    }

    const bool branchB = (pMarg <= m_marginalCutoff);

    // Branch B: genotype-adjusted residual R̃ = R − α − β·G, projecting the
    // marginal genetic effect out.  [α, β]ᵀ = ([1,G]ᵀ[1,G])⁻¹ [1,G]ᵀ R.
    // R̃ is env-independent, so it is formed once per marker.
    Eigen::VectorXd R0;
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
        }
        // det ≤ 0 (monomorphic G): leave R0 empty ⇒ per-env NaN below.
    }

    for (int e = 0; e < nEnv; ++e) {
        const int base = 4 + 5 * e;
        double z = 0.0, var = 0.0, p, score;
        if (!branchB) {
            // Branch A: precomputed λ-orthogonalised weight.
            const EnvData &ed = m_envs[e];
            score = rawScores[1 + e]; // Σ G_i w_{e,i}
            p = spaScorePval(ed.wOutlier, ed.wSum, ed.wSq, score, q, m_spaCutoff,
                             z, var);
        } else {
            if (R0.size() == 0) continue; // degenerate G ⇒ leave NaN
            // Branch B per-marker weight u = E ∘ R̃; Ŝ_GxE = Σ G_i u_i.
            const Eigen::VectorXd u = m_envs[e].E.array() * R0.array();
            score = GVec.dot(u);
            const OutlierData ub = detectOutliers(u, m_outlierRatio);
            p = spaScorePval(ub, u.sum(), u.squaredNorm(), score, q, m_spaCutoff,
                             z, var);
        }
        if (var > 0.0) {
            out[base + 0] = p;                       // P_Gx<E>
            out[base + 1] = math::zFromPval(p, z);   // Z_Gx<E> (p-consistent)
            out[base + 2] = z;                       // Z_Norm_Gx<E> (raw score z)
            out[base + 3] = score / var;             // BETA_Gx<E>
            out[base + 4] = 1.0 / std::sqrt(var);    // SE_Gx<E>
        }
    }
}

void SPAGxEMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int /*markerInChunkIdx*/,
    std::vector<double> &result
) {
    const Eigen::VectorXd rawScores = m_W.transpose() * GVec;
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
    // Fused: marginal + every Branch-A G×E score in one GEMM.
    //   scores(b, 0)     = Σ G_b·R
    //   scores(b, 1 + e) = Σ G_b·w_e
    const Eigen::MatrixXd scores = GBatch.transpose() * m_W; // B × (1+nEnv)
    for (int b = 0; b < B; ++b) {
        const Eigen::VectorXd rs = scores.row(b).transpose();
        evalMarker(GBatch.col(b), altFreqs[b], rs, results[b]);
    }
}

// ══════════════════════════════════════════════════════════════════════
// runSPAGxE — full workflow (mirrors runSPAGRM, no sparse GRM in Phase 1)
// ══════════════════════════════════════════════════════════════════════

void runSPAGxE(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &envNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
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
    // Phase 1 is the base (unrelated) test; the sparse-GRM (+) variant is
    // wired in Phase 2.  Accept but ignore a supplied GRM for now.
    if (!spgrmGrabFile.empty() || !spgrmGctaFile.empty())
        warnMsg("SPAGxE: the sparse-GRM (+) variant is not yet enabled in this "
                "build; running the base unrelated test (GRM ignored).");
    (void)pairwiseIBDFile;

    const bool fitPath = !phenoNameSpec.empty();
    nullmodel::RegressionModel regModel{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        regModel = nullmodel::parseRegressionModel(regressionModelStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
        infoMsg("SPAGxE: fitting %s null model for %zu phenotype(s)",
                nullmodel::regressionModelName(regModel), phenoSpecs.size());
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
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();

    // Environment columns must enter the null model in fit mode (so R ⟂ E),
    // exactly as SAGELD requires --sageld-x ⊆ --covar-name.
    if (fitPath) {
        for (const auto &en : envNames)
            if (std::find(covarNames.begin(), covarNames.end(), en) == covarNames.end())
                throw std::runtime_error(
                    "SPAGxE: environment '" + en +
                    "' must also appear in --covar-name (it enters the "
                    "genotype-independent null model).");
    }
    // Drop subjects with a missing environment value (a no-op in fit mode,
    // where E is a covariate; necessary in residual mode, where E is only
    // used to form λ and is not otherwise NA-filtered).
    sd.dropNaInColumns(envNames);
    const uint32_t N = sd.nUsed();
    infoMsg("  %u subjects in union mask", N);

    if (fitPath) {
        Eigen::MatrixXd covarUnion;
        if (!covarNames.empty()) covarUnion = sd.getColumns(covarNames);
        else covarUnion.resize(N, 0);
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

        tasks[rc].phenoName = pi.name;
        tasks[rc].method = std::make_unique<SPAGxEMethod>(
            std::move(phenoResid), envNames, std::move(envVecs), marginalCutoff,
            spaCutoff, outlierIqrRatio);
        tasks[rc].unionToLocal = pi.unionToLocal;
        tasks[rc].nUsed = pi.nUsed;
        infoMsg("  Phenotype '%s': %u subjects", pi.name.c_str(), pi.nUsed);
    }

    infoMsg("Running SPAGxE marker tests (%d thread(s), %d phenotype(s))...",
            nthreads, K);
    multiPhenoEngine(*genoData, tasks, outPrefix, "SPAGxE", compression,
                     compressionLevel, nthreads, missingCutoff, minMafCutoff,
                     minMacCutoff, hweCutoff);
}
