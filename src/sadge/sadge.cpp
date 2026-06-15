// sadge.cpp — SADGE method implementation + runSADGE workflow.
#include "sadge/sadge.hpp"

#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "sadge/impute_kernel.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/null_model.hpp"

#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace sadge {

namespace {
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Two-sided normal-approximation z and p from score S, mean E, variance V.
inline void scoreZP(double S, double E, double V, double &z, double &p) {
    if (V > 0.0) {
        z = (S - E) / std::sqrt(V);
        p = 2.0 * math::pnorm(-std::abs(z));
    } else {
        z = kNaN;
        p = kNaN;
    }
}
} // namespace

SADGEMethod::SADGEMethod(
    std::shared_ptr<const SADGEShared> shared,
    std::shared_ptr<const PerPhenoData> pheno
)
    : m_shared(std::move(shared)),
      m_pheno(std::move(pheno))
{
}

std::unique_ptr<MethodBase> SADGEMethod::clone() const {
    // m_shared / m_pheno are immutable and shared; the per-thread worker
    // context is acquired lazily and never copied.
    return std::make_unique<SADGEMethod>(*this);
}

std::string SADGEMethod::getHeaderColumns() const {
    return "\tG_I.Z\tG_I.P\tG_I.S\tG_I.VarS"
           "\tG_D.Z\tG_D.P\tG_D.S\tG_D.VarS"
           "\tG.Z\tG.P\tG.S\tG.VarS"
           "\tSADGE_MAF";
}

SADGEWorkerCtx &SADGEMethod::workerCtx() {
    // One context per worker thread, shared by the K phenotype-clones running
    // on that thread (keyed by the shared object so concurrent SADGE runs do
    // not collide).  G_par is phenotype-independent, so the chunk is decoded
    // and imputed once per thread, not once per phenotype.
    static thread_local std::unordered_map<const SADGEShared *,
                                           std::shared_ptr<SADGEWorkerCtx> >
        tls;
    auto &slot = tls[m_shared.get()];
    if (!slot) {
        slot = std::make_shared<SADGEWorkerCtx>();
        slot->cursor = m_shared->parGeno->makeCursor();
        slot->Gp.resize(m_shared->nPedUsed);
        slot->missScratch.reserve(m_shared->nPedUsed / 8 + 1);
    }
    return *slot;
}

void SADGEMethod::prepareChunk(const std::vector<uint64_t> &gIdx) {
    if (gIdx.empty()) return;
    SADGEWorkerCtx &ctx = workerCtx();
    // prepareChunk is called once per clone; the K clones share ctx, so only
    // reset when the chunk genuinely changes.
    if (!ctx.chunkGIdx.empty() && ctx.chunkGIdx.size() == gIdx.size() &&
        ctx.chunkGIdx.front() == gIdx.front())
        return;
    ctx.chunkGIdx = gIdx;
    ctx.cache.clear();
    ctx.decodedUpTo = -1;
    ctx.cursor->beginSequentialBlock(gIdx.front()); // resets to file data start
}

const MarkerCache &SADGEMethod::ensureMarker(SADGEWorkerCtx &ctx, int chunkRel) {
    // Decode ascending up to chunkRel (the BGEN cursor only advances forward;
    // contiguous ascending decode keeps the cursor in sync and lets later,
    // possibly out-of-order per-phenotype passing subsets hit the cache).
    while (ctx.decodedUpTo < chunkRel) {
        const int c = ++ctx.decodedUpTo;
        const uint64_t g = ctx.chunkGIdx[static_cast<size_t>(c)];
        double altFreq, altCounts, missingRate, hweP, maf, mac;
        ctx.cursor->getGenotypes(g, ctx.Gp, altFreq, altCounts, missingRate,
                                 hweP, maf, mac, ctx.missScratch);
        // Mean-impute any missing dosages (none for imputed BGEN dosages).
        for (uint32_t idx : ctx.missScratch)
            ctx.Gp[idx] = 2.0 * altFreq;

        const double mu = altFreq; // pedigree-set allele frequency (matches R)
        PerMarkerCoeffs co;
        precomputeCoeffs(mu, co);

        MarkerCache mc;
        mc.mu = mu;
        mc.H = co.H;
        mc.I = co.I;
        mc.J = co.J;
        const auto &NS = m_shared->nonSing;
        mc.gImp.resize(static_cast<Eigen::Index>(NS.size()));
        const double *Gp = ctx.Gp.data();
        for (size_t j = 0; j < NS.size(); ++j) {
            const NonSingKernel &k = NS[j];
            const double gSelf = Gp[k.pedSelf];
            const double gCo = (k.pedCoSib >= 0) ? Gp[k.pedCoSib] : 0.0;
            const double gP1 = (k.pedPar1 >= 0) ? Gp[k.pedPar1] : 0.0;
            const double gP2 = (k.pedPar2 >= 0) ? Gp[k.pedPar2] : 0.0;
            mc.gImp[static_cast<Eigen::Index>(j)] =
                imputeGpar(k.fs1, gSelf, gCo, gP1, gP2, co);
        }
        ctx.cache.emplace(c, std::move(mc));
    }
    return ctx.cache.at(chunkRel);
}

void SADGEMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd>,
    double,
    int,
    std::vector<double> &
) {
    // SADGE is fused-only and always runs through multiPhenoEngine's fused
    // path (processScoreBatch); the per-marker GEMV path is never used.
    throw std::logic_error("SADGE: getResultVec called on a fused-only method");
}

void SADGEMethod::fillUnionResiduals(
    Eigen::Ref<Eigen::MatrixXd> dest,
    const std::vector<uint32_t> & /*unionToLocal*/
) const {
    // dest is pre-zeroed, nUnion × 2.  residUnion / residSing are already in
    // union order with 0 at absent subjects, so copy directly.
    dest.col(0) = m_pheno->residUnion;
    dest.col(1) = m_pheno->residSing;
}

void SADGEMethod::fillResidualSums(double *dest) const {
    dest[0] = m_pheno->sumR;
    dest[1] = m_pheno->sumRSing;
}

void SADGEMethod::processScoreBatch(
    const Eigen::Ref<const Eigen::MatrixXd> &scores,
    const double * /*gSums*/,
    const double * /*gSumSqs*/,
    uint32_t /*nUsed*/,
    const std::vector<double> & /*altFreqs*/,
    const std::vector<int> &chunkIdxs,
    std::vector<std::vector<double> > &results
) {
    const int B = static_cast<int>(chunkIdxs.size());
    results.resize(B);
    if (B == 0) return;

    SADGEWorkerCtx &ctx = workerCtx();
    // Decode (ascending) every marker up to the largest needed index; passing
    // subsets are ascending so chunkIdxs.back() is the maximum.
    ensureMarker(ctx, chunkIdxs.back());
    // Bound the cache to a few windows' worth of markers (the current window's
    // markers are all >= decodedUpTo - windowSpan, so this never evicts a
    // marker still needed by another phenotype of the same window).
    if (!ctx.cache.empty()) {
        const int lo = ctx.decodedUpTo - 256;
        if (lo > 0)
            for (auto it = ctx.cache.begin(); it != ctx.cache.end();)
                it = (it->first < lo) ? ctx.cache.erase(it) : std::next(it);
    }

    const PerPhenoData &ph = *m_pheno;
    const double sumR = ph.sumR, Rsq = ph.Rsq, Rpr = ph.Rpr;

    for (int b = 0; b < B; ++b) {
        const MarkerCache &mc = ctx.cache.at(chunkIdxs[b]);
        const double mu = mc.mu;

        const double S_G = scores(0, b);                          // Gᵀ R
        const double singG = scores(1, b);                        // Σ_singleton G·R
        const double nonSingScore = mc.gImp.dot(ph.residNonSingComp);
        const double S_par = nonSingScore + singG + 2.0 * mu * ph.sumRSing;

        const double varG = 2.0 * mu * (1.0 - mu);
        const double covGGsib = mu * (1.0 - mu);
        const double covGGpar = 2.0 * mu * (1.0 - mu);
        const double rho = interpRho(ph.rhoGrid, mu);

        std::array<double, kNStruct> vgp;
        varGenoPar(mu, mc.H, mc.I, mc.J, vgp);
        double var_par_sq = 0.0, var_par_pr = 0.0;
        for (int fs = 0; fs < kNStruct; ++fs)
            var_par_sq += ph.sumSquare[fs] * vgp[fs];
        // 2-sib structures only for the cross (product) term.
        for (int fs = static_cast<int>(FamStruct::s0par2sib);
             fs <= static_cast<int>(FamStruct::s2par2sib); ++fs)
            var_par_pr += ph.sumProd[fs] * vgp[fs];

        // ── G (total) ──
        const double E_G = 2.0 * mu * sumR;
        const double Var_G = Rsq * varG + 2.0 * Rpr * covGGsib;
        double zG, pG;
        scoreZP(S_G, E_G, Var_G, zG, pG);

        // ── G_D (direct) ──
        const double hr = 0.5 + rho;
        const double S_D = S_G - hr * S_par;
        const double E_D = -rho * 4.0 * mu * sumR;
        const double vssD = hr * hr * var_par_sq + (varG - 2.0 * hr * covGGpar) * Rsq;
        const double vspD = hr * hr * var_par_pr + (covGGsib - 2.0 * hr * covGGpar) * Rpr;
        const double Var_D = vssD + 2.0 * vspD;
        double zD, pD;
        scoreZP(S_D, E_D, Var_D, zD, pD);

        // ── G_I (indirect / nurture) ──
        const double S_I = S_par - S_G;
        const double E_I = 2.0 * mu * sumR;
        const double vssI = var_par_sq + (varG - 2.0 * covGGpar) * Rsq;
        const double vspI = var_par_pr + (covGGsib - 2.0 * covGGpar) * Rpr;
        const double Var_I = vssI + 2.0 * vspI;
        double zI, pI;
        scoreZP(S_I, E_I, Var_I, zI, pI);

        std::vector<double> &out = results[b];
        out.clear();
        out.reserve(13);
        out.push_back(zI);
        out.push_back(pI);
        out.push_back(S_I);
        out.push_back(Var_I);
        out.push_back(zD);
        out.push_back(pD);
        out.push_back(S_D);
        out.push_back(Var_D);
        out.push_back(zG);
        out.push_back(pG);
        out.push_back(S_G);
        out.push_back(Var_G);
        out.push_back(mu);
    }
}

} // namespace sadge

// ══════════════════════════════════════════════════════════════════════
// runSADGE
// ══════════════════════════════════════════════════════════════════════

void runSADGE(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &sadgeFamFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
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
    bool saveResid,
    uint64_t seed
) {
    const bool fitPath = !phenoNameSpec.empty();
    nullmodel::RegressionModel regModel{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        regModel = nullmodel::parseRegressionModel(regressionModelStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
        infoMsg("SADGE: fitting %s null model for %zu phenotype(s)",
                nullmodel::regressionModelName(regModel), phenoSpecs.size());
    }

    auto genoIIDs = parseGenoIIDs(geno);

    infoMsg("SADGE: loading pedigree %s", sadgeFamFile.c_str());
    sadge::SampleInfo si = sadge::buildSampleInfo(sadgeFamFile, genoIIDs);
    infoMsg("SADGE: %zu genotyped siblings, %zu pedigree members (siblings+parents)",
            si.siblings.size(), si.pedigreeIIDs.size());

    // ── Subject data: restrict the analysis union to siblings ─────────────
    SubjectData sd(std::move(genoIIDs));
    if (fitPath) {
        sd.loadPhenoFile(phenoFile, nullmodel::columnsNeeded(phenoSpecs));
    } else {
        sd.loadResidOne(phenoFile, residNames);
    }
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    {
        std::unordered_set<std::string> sibSet;
        sibSet.reserve(si.siblings.size() * 2);
        for (const auto &s : si.siblings)
            sibSet.insert(s.iid);
        sd.setKeepSubjects(std::move(sibSet)); // union ⊆ genotyped siblings
    }
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();
    const uint32_t N = sd.nUsed();
    infoMsg("SADGE: %u sibling subjects in analysis union", N);

    if (fitPath) {
        Eigen::MatrixXd covarUnion;
        if (!covarNames.empty()) {
            covarUnion = sd.getColumns(covarNames);
        } else if (sd.hasCovar()) {
            covarUnion = sd.covar();
        } else {
            covarUnion.resize(N, 0);
        }
        nullmodel::EngineOptions eo;
        eo.nthreads = nthreads;
        eo.seed = seed;
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
            nullmodel::writeResidualsFile(outPrefix + ".null.resid", sd, dumpFits);
        }
        sd.setResidualsFromFit(std::move(rs), std::move(ns));
    }

    // ── SADGE's own pedigree-set genotype meta (siblings + parents) ───────
    const auto &genoIIDsRef = sd.famIIDs();
    const uint32_t nFam = sd.nFam();
    std::vector<uint64_t> pedMask((nFam + 63) / 64, 0);
    uint32_t nPed = 0;
    {
        std::unordered_set<std::string> pedSet(
            si.pedigreeIIDs.begin(), si.pedigreeIIDs.end());
        for (uint32_t f = 0; f < nFam; ++f) {
            if (pedSet.count(genoIIDsRef[f])) {
                pedMask[f / 64] |= (1ULL << (f % 64));
                ++nPed;
            }
        }
    }
    auto parGenoU = makeGenoData(geno, pedMask, nFam, nPed, nSnpPerChunk);
    auto engineGeno = makeGenoData(geno, sd.usedMask(), nFam, N, nSnpPerChunk);
    if (parGenoU->nMarkers() != engineGeno->nMarkers())
        throw std::runtime_error(
            "SADGE: pedigree and analysis genotype metas disagree on marker count");
    infoMsg("SADGE: %u markers, %u sibling + %u pedigree subjects",
            engineGeno->nMarkers(), N, nPed);

    // ── Shared read-only state ────────────────────────────────────────────
    auto shared = std::make_shared<sadge::SADGEShared>();
    shared->parGeno = std::shared_ptr<const GenoMeta>(std::move(parGenoU));
    shared->nUnion = N;
    shared->nPedUsed = nPed;

    const auto usedIIDs = sd.usedIIDs();
    std::vector<int> unionToNonSing(N, -1);
    for (uint32_t u = 0; u < N; ++u) {
        auto it = si.sibIndexByIID.find(usedIIDs[u]);
        if (it == si.sibIndexByIID.end())
            throw std::runtime_error(
                "SADGE: union subject not found among siblings: " + usedIIDs[u]);
        const sadge::SibInfo &s = si.siblings[it->second];
        if (s.fs1 != sadge::FamStruct::s0par1sib) {
            unionToNonSing[u] = static_cast<int>(shared->nonSing.size());
            shared->nonSing.push_back(
                {u, s.fs1, s.pedSelf, s.pedCoSib, s.pedPar1, s.pedPar2});
        }
    }
    const int nNonSing = static_cast<int>(shared->nonSing.size());

    // ── Family grouping over the union (for per-phenotype precond) ─────────
    std::vector<sadge::FamilyGroup> familyGroups;
    {
        std::unordered_map<long, int> famIdxOf;
        for (uint32_t u = 0; u < N; ++u) {
            const sadge::SibInfo &s = si.siblings[si.sibIndexByIID.at(usedIIDs[u])];
            auto it = famIdxOf.find(s.fid);
            if (it == famIdxOf.end()) {
                famIdxOf.emplace(s.fid, static_cast<int>(familyGroups.size()));
                sadge::FamilyGroup g;
                g.genoNPar = s.genoNPar;
                g.n = 1;
                g.u = {static_cast<int>(u), -1};
                familyGroups.push_back(g);
            } else {
                sadge::FamilyGroup &g = familyGroups[it->second];
                if (g.n < 2) g.u[g.n] = static_cast<int>(u);
                ++g.n;
            }
        }
    }

    // ── Per-phenotype tasks ───────────────────────────────────────────────
    auto phenoInfos = sd.buildPerColumnMasks();
    const int K = sd.residOneCols();
    if (K > 1) infoMsg("SADGE: %d phenotypes", K);

    std::vector<PhenoTask> tasks(K);
    for (int rc = 0; rc < K; ++rc) {
        const auto &pi = phenoInfos[rc];
        Eigen::VectorXd residRaw;
        if (K == 1)
            residRaw = sd.residuals();
        else
            residRaw = sd.residMatrix().col(rc);
        auto ph = std::make_shared<sadge::PerPhenoData>(
            sadge::buildPerPhenoData(residRaw, familyGroups, N, unionToNonSing, nNonSing));
        tasks[rc].phenoName = pi.name;
        tasks[rc].method = std::make_unique<sadge::SADGEMethod>(shared, std::move(ph));
        tasks[rc].unionToLocal = pi.unionToLocal;
        tasks[rc].nUsed = pi.nUsed;
        infoMsg("  Phenotype '%s': %u subjects", pi.name.c_str(), pi.nUsed);
    }

    infoMsg("SADGE: running marker tests (%d thread(s), %d phenotype(s))...", nthreads, K);
    multiPhenoEngine(*engineGeno, tasks, outPrefix, "SADGE", compression,
                     compressionLevel, nthreads, missingCutoff, minMafCutoff,
                     minMacCutoff, hweCutoff);
}
