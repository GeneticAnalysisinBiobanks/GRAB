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
    // not collide).  G_par is phenotype-independent, so each window is imputed
    // once per thread (from the engine's shared genotype window), not once per
    // phenotype.
    static thread_local std::unordered_map<const SADGEShared *,
                                           std::shared_ptr<SADGEWorkerCtx> >
        tls;
    auto &slot = tls[m_shared.get()];
    if (!slot)
        slot = std::make_shared<SADGEWorkerCtx>();
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
    ctx.imputedWindowStart = -1;
}

void SADGEMethod::imputeMarker(
    MarkerCache &mc, const Eigen::Ref<const Eigen::VectorXd> &gCol) const {
    // The engine has already mean-imputed missing dosages to the pedigree-set
    // mean (2·AF_ped) before the GEMM, so the column mean / 2 is the pedigree
    // allele frequency (matches the R reference's mu).
    const double mu = gCol.sum() / (2.0 * static_cast<double>(m_shared->nUnion));
    PerMarkerCoeffs co;
    precomputeCoeffs(mu, co);

    mc.mu = mu;
    mc.H = co.H;
    mc.I = co.I;
    mc.J = co.J;
    const auto &NS = m_shared->nonSing;
    mc.gImp.resize(static_cast<Eigen::Index>(NS.size()));
    const double *Gp = gCol.data();
    for (size_t j = 0; j < NS.size(); ++j) {
        const NonSingKernel &k = NS[j];
        const double gSelf = Gp[k.pedSelf];
        const double gCo = (k.pedCoSib >= 0) ? Gp[k.pedCoSib] : 0.0;
        const double gP1 = (k.pedPar1 >= 0) ? Gp[k.pedPar1] : 0.0;
        const double gP2 = (k.pedPar2 >= 0) ? Gp[k.pedPar2] : 0.0;
        mc.gImp[static_cast<Eigen::Index>(j)] =
            imputeGpar(k.fs1, gSelf, gCo, gP1, gP2, co);
    }
}

void SADGEMethod::onFusedGenoWindow(
    const Eigen::Ref<const Eigen::MatrixXd> &Gwin, int windowStart) {
    SADGEWorkerCtx &ctx = workerCtx();
    // The K phenotype-clones on this thread all receive the identical window;
    // impute it only once (G_par is phenotype-independent).
    if (ctx.imputedWindowStart == windowStart) return;
    ctx.imputedWindowStart = windowStart;

    const int wlen = static_cast<int>(Gwin.cols());
    for (int bi = 0; bi < wlen; ++bi) {
        const int chunkRel = windowStart + bi;
        MarkerCache mc;
        imputeMarker(mc, Gwin.col(bi));
        ctx.cache[chunkRel] = std::move(mc);
        if (chunkRel > ctx.decodedUpTo) ctx.decodedUpTo = chunkRel;
    }
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
    // The genotype window for these markers has already been imputed by
    // onFusedGenoWindow (Phase 2b) before this call, so the cache is populated.
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

    // ── Engine genotype meta = the pedigree set (siblings + parents) ──────
    // SADGE decodes the genotype file ONCE, over the pedigree superset.  The
    // fused GEMM rides on these rows: parents (and non-phenotyped siblings)
    // carry zero residual and zero mask, so S_G / the singleton partial are
    // unchanged, while the same decoded window feeds the parental-genotype
    // imputation via onFusedGenoWindow.  Engine row order == pedigree-decode
    // order == genotype-file order, so si.pedRowByIID rows (pedSelf, ...) are
    // valid engine-union rows with no remapping.
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
    auto engineGeno = makeGenoData(geno, pedMask, nFam, nPed, nSnpPerChunk);
    infoMsg("SADGE: %u markers, %u sibling (union) + %u pedigree subjects",
            engineGeno->nMarkers(), N, nPed);

    // ── Shared read-only state ────────────────────────────────────────────
    auto shared = std::make_shared<sadge::SADGEShared>();
    shared->nUnion = nPed; // engine union = pedigree-set rows

    // Build the fs1-non-singleton kernels and family grouping over the
    // phenotyped union siblings (as before), but in PEDIGREE-row coordinates
    // (s.pedSelf), so they index directly into the engine union / GBatch_union.
    const auto usedIIDs = sd.usedIIDs();
    std::vector<int> unionRowToPed(N);          // union sibling row → pedigree row
    std::vector<int> unionToNonSing(nPed, -1);  // pedigree row → compact non-sing idx
    std::vector<sadge::FamilyGroup> familyGroups;
    {
        std::unordered_map<long, int> famIdxOf;
        for (uint32_t u = 0; u < N; ++u) {
            auto it = si.sibIndexByIID.find(usedIIDs[u]);
            if (it == si.sibIndexByIID.end())
                throw std::runtime_error(
                    "SADGE: union subject not found among siblings: " + usedIIDs[u]);
            const sadge::SibInfo &s = si.siblings[it->second];
            const int row = s.pedSelf; // engine union row of this sibling
            unionRowToPed[u] = row;

            if (s.fs1 != sadge::FamStruct::s0par1sib) {
                unionToNonSing[row] = static_cast<int>(shared->nonSing.size());
                shared->nonSing.push_back(
                    {static_cast<uint32_t>(row), s.fs1,
                     s.pedSelf, s.pedCoSib, s.pedPar1, s.pedPar2});
            }

            auto fit = famIdxOf.find(s.fid);
            if (fit == famIdxOf.end()) {
                famIdxOf.emplace(s.fid, static_cast<int>(familyGroups.size()));
                sadge::FamilyGroup g;
                g.genoNPar = s.genoNPar;
                g.n = 1;
                g.u = {row, -1};
                familyGroups.push_back(g);
            } else {
                sadge::FamilyGroup &g = familyGroups[fit->second];
                if (g.n < 2) g.u[g.n] = row;
                ++g.n;
            }
        }
    }
    const int nNonSing = static_cast<int>(shared->nonSing.size());

    // ── Per-phenotype tasks ───────────────────────────────────────────────
    auto phenoInfos = sd.buildPerColumnMasks(); // used for phenotype names
    const int K = sd.residOneCols();
    if (K > 1) infoMsg("SADGE: %d phenotypes", K);

    const double kNaNd = std::numeric_limits<double>::quiet_NaN();
    std::vector<PhenoTask> tasks(K);
    for (int rc = 0; rc < K; ++rc) {
        // Sibling-union residuals (NaN where this phenotype is absent).
        const Eigen::VectorXd residSib =
            (K == 1) ? sd.residuals() : Eigen::VectorXd(sd.residMatrix().col(rc));

        // Scatter into pedigree-row coordinates; NaN at parents and at
        // siblings absent from this phenotype (→ mask 0, residual 0).
        Eigen::VectorXd residRaw = Eigen::VectorXd::Constant(nPed, kNaNd);
        for (uint32_t u = 0; u < N; ++u)
            residRaw[unionRowToPed[u]] = residSib[u];

        auto ph = std::make_shared<sadge::PerPhenoData>(
            sadge::buildPerPhenoData(residRaw, familyGroups, nPed, unionToNonSing, nNonSing));

        // unionToLocal: mark phenotype-present rows (non-NaN residual).  SADGE
        // is fused-only, so only the present/absent distinction (the AugResid
        // mask) and nUsed (the altFreq denominator) matter; the local values
        // are never gathered.
        std::vector<uint32_t> u2l(nPed, UINT32_MAX);
        uint32_t loc = 0;
        for (uint32_t r = 0; r < nPed; ++r)
            if (!std::isnan(residRaw[r])) u2l[r] = loc++;

        tasks[rc].phenoName = phenoInfos[rc].name;
        tasks[rc].method = std::make_unique<sadge::SADGEMethod>(shared, std::move(ph));
        tasks[rc].unionToLocal = std::move(u2l);
        tasks[rc].nUsed = loc;
        infoMsg("  Phenotype '%s': %u subjects", phenoInfos[rc].name.c_str(), loc);
    }

    infoMsg("SADGE: running marker tests (%d thread(s), %d phenotype(s))...", nthreads, K);
    multiPhenoEngine(*engineGeno, tasks, outPrefix, "SADGE", compression,
                     compressionLevel, nthreads, missingCutoff, minMafCutoff,
                     minMacCutoff, hweCutoff);
}
