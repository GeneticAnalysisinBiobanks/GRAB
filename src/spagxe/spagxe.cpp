// spagxe.cpp — SPAGxE_CCT G×E test implementation.
//
// Base (unrelated) and SPAGxE+ (sparse-GRM) share one code path: the GRM Φ
// enters only through the retrospective score variance as a quadratic form
// xᵀΦx = SparseGRM::spaVariance(x) = 2·Σstored − Σx², and the base path is the
// Φ = identity special case (spaVariance → Σx²).  Per phenotype a single
// genotype-independent residual R is fit once (trait ~ X + E).  Per variant:
//   1. marginal screen  S_G = Σ G_i R_i ,  Var = 2q(1−q)·RᵀΦR  (always normal;
//      routes Branch A vs Branch B on p_marg vs ε).
//   2a. Branch A (p_marg > ε): λ_GRM = RᵀΦR_E / RᵀΦR, w = (E − λ_GRM)R,
//       S_GxE = Σ G_i w_i, Var = 2q(1−q)·wᵀΦw, SPA / normal hybrid.
//   2b. Branch B (p_marg ≤ ε): R̃ = R − [1,G]([1,G]ᵀ[1,G])⁻¹[1,G]ᵀR (unweighted
//       even for +), Ŝ_GxE = Σ G_i E_i R̃_i, Var = 2q(1−q)·(E∘R̃)ᵀΦ(E∘R̃).
// The tail SPA rescales the statistic by √(indepVar/grmVar) (SAIGE variance
// ratio, the paper's non-mix convention) and applies the independence-CGF SPA
// via spa::getProbSpaG with the IQR outlier / Gaussian non-outlier split.  A
// degenerate saddlepoint yields NaN, never 0 (GRAB2 convention).  Reference:
// tmp/SPAGxECCT/.../R  (SPAGxE_CCT_one_SNP, SPA_G_one_SNP_homo_new,
// SPAGxE_Plus_one_SNP, SPAGxE_Plus_Nullmodel).

#include "spagxe/spagxe.hpp"

#include "geno_factory/geno_data.hpp" // parseGenoIIDs, makeGenoData, GenoMeta
#include "io/sparse_grm.hpp"          // SparseGRM
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
// homogeneous binomial law  G_i ~ Bin(2, q).  The GRM enters only through
// grmQuad = wᵀΦw (== wSq = Σw² in the base / unrelated case), which sets the
// score variance and — in the tail — a SAIGE-style variance-ratio rescaling of
// the statistic before the independence-CGF SPA (spa::getProbSpaG with the IQR
// outlier / Gaussian non-outlier split).  Mirrors SPA_G_one_SNP_homo_new /
// SPAGxE_Plus_one_SNP: normal approximation when |z| ≤ spaCutoff, otherwise the
// reflection-about-mean two-sided saddlepoint.
//   part     — IQR partition of w  (residOutlier = w[outlier],
//              resid2NonOutlier = w²[nonOutlier])
//   wSum     — Σ w      wSq — Σ w²      grmQuad — wᵀΦw   (all over every i)
//   rawScore — Σ G_i w_i
// Sets zScore (raw score z) and outVar = Var(S); returns the p-value (NaN when
// the variance is non-positive or the saddlepoint fails to converge).
double spaScorePval(
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
        return kNaN;
    }
    const double sMean = 2.0 * q * wSum; // retrospective mean E[S] = 2q·Σw
    zScore = (rawScore - sMean) / std::sqrt(scoreVar);
    outVar = scoreVar;
    if (std::abs(zScore) <= spaCutoff)
        return 2.0 * math::pnorm(-std::abs(zScore));

    // ── Saddlepoint tail: variance-ratio rescale + outlier / non-outlier split.
    // sqrtRatio = √(indepVar / grmVar) = √(wSq / grmQuad) (the 2q(1−q) cancels);
    // == 1 in the base case, so this reduces to the plain independence SPA.
    const double sqrtRatio = (grmQuad > 0.0) ? std::sqrt(wSq / grmQuad) : 1.0;
    const double absDev = std::abs(rawScore - sMean) * sqrtRatio;
    const int nOut = static_cast<int>(part.posOutlier.size());
    const double meanNon = 2.0 * q * part.residNonOutlier.sum();
    const double varNon = f * part.resid2NonOutlier.sum();
    std::vector<double> mafOut(static_cast<size_t>(nOut), q);
    // Reflect about the fitted mean: upper = sMean + |Δ|, lower = sMean − |Δ|.
    const double pUpper = spa::getProbSpaG(mafOut.data(), part.residOutlier.data(),
                                           nOut, sMean + absDev, false, meanNon, varNon);
    const double pLower = spa::getProbSpaG(mafOut.data(), part.residOutlier.data(),
                                           nOut, sMean - absDev, true, meanNon, varNon);
    return pUpper + pLower;
}

} // namespace

// ══════════════════════════════════════════════════════════════════════
// SPAGxEMethod
// ══════════════════════════════════════════════════════════════════════

double SPAGxEMethod::phiQuad(const Eigen::VectorXd &x) const {
    if (m_grm)
        return m_grm->spaVariance(x.data(), static_cast<uint32_t>(x.size()));
    return x.squaredNorm();
}

SPAGxEMethod::SPAGxEMethod(
    Eigen::VectorXd resid,
    std::vector<std::string> envNames,
    std::vector<Eigen::VectorXd> envVecs,
    std::shared_ptr<const SparseGRM> grm,
    double marginalCutoff,
    double spaCutoff,
    double outlierRatio
)
    : m_resid(std::move(resid)),
      m_envNames(std::move(envNames)),
      m_grm(std::move(grm)),
      m_marginalCutoff(marginalCutoff),
      m_spaCutoff(spaCutoff),
      m_outlierRatio(outlierRatio)
{
    const Eigen::Index N = m_resid.size();
    m_residSum = m_resid.sum();
    m_R_phiQuad = phiQuad(m_resid); // RᵀΦR (= ΣR² in the base case)
    const int nEnv = static_cast<int>(m_envNames.size());

    // m_W = [ R | w_1 | … | w_nEnv ]; one GEMM against a genotype batch then
    // yields the marginal score (col 0) and every Branch-A G×E score.
    m_W.resize(N, 1 + nEnv);
    m_W.col(0) = m_resid;

    m_envs.reserve(nEnv);
    for (int e = 0; e < nEnv; ++e) {
        EnvData ed;
        ed.name = m_envNames[e];
        ed.E = std::move(envVecs[e]);

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
    // S_G = Σ G_i R_i (uncentered, matching the reference S1; E[S_G]=2q·ΣR≈0),
    // variance 2q(1−q)·RᵀΦR.
    const double varSG = m_R_phiQuad * 2.0 * q * (1.0 - q);
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
    // R̃ is env-independent, so it is formed once per marker.  (Unweighted even
    // for +; the GRM re-enters only through the variance — paper eq. 3.)
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
            // Branch A: precomputed λ-orthogonalised weight and its Φ-form.
            const EnvData &ed = m_envs[e];
            score = rawScores[1 + e]; // Σ G_i w_{e,i}
            p = spaScorePval(ed.wOutlier, ed.wSum, ed.wSq, ed.wPhiQuad, score, q,
                             m_spaCutoff, z, var);
        } else {
            if (R0.size() == 0) continue; // degenerate G ⇒ leave NaN
            // Branch B per-marker weight u = E ∘ R̃; Ŝ_GxE = Σ G_i u_i.
            const Eigen::VectorXd u = m_envs[e].E.array() * R0.array();
            score = GVec.dot(u);
            const OutlierData ub = detectOutliers(u, m_outlierRatio);
            p = spaScorePval(ub, u.sum(), u.squaredNorm(), phiQuad(u), score, q,
                             m_spaCutoff, z, var);
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
// runSPAGxE — full workflow (mirrors runSPAGRM)
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
    const bool hasGrm = !spgrmGrabFile.empty() || !spgrmGctaFile.empty();

    const bool fitPath = !phenoNameSpec.empty();
    nullmodel::RegressionModel regModel{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        regModel = nullmodel::parseRegressionModel(regressionModelStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
        infoMsg("SPAGxE%s: fitting %s null model for %zu phenotype(s)",
                hasGrm ? "+" : "", nullmodel::regressionModelName(regModel),
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

        tasks[rc].phenoName = pi.name;
        tasks[rc].method = std::make_unique<SPAGxEMethod>(
            std::move(phenoResid), envNames, std::move(envVecs), phenoGrm,
            marginalCutoff, spaCutoff, outlierIqrRatio);
        tasks[rc].unionToLocal = pi.unionToLocal;
        tasks[rc].nUsed = pi.nUsed;
        infoMsg("  Phenotype '%s': %u subjects", pi.name.c_str(), pi.nUsed);
    }

    infoMsg("Running SPAGxE%s marker tests (%d thread(s), %d phenotype(s))...",
            hasGrm ? "+" : "", nthreads, K);
    multiPhenoEngine(*genoData, tasks, outPrefix, "SPAGxE", compression,
                     compressionLevel, nthreads, missingCutoff, minMafCutoff,
                     minMacCutoff, hweCutoff);
}
