// spagrm.cpp — nsSPAGRM free functions, SPAGRMClass, and runSPAGRM

#include "spagrm/spagrm.hpp"
#include "engine/loco.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "spagrm/grm_null.hpp"
#include "spagrm/longitudinal_resid.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/null_model.hpp"
#include "util/text_scanner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>

namespace nsSPAGRM {

// ══════════════════════════════════════════════════════════════════════
// Two-sided saddlepoint p-value — on the shared solver and tail kernel
// ══════════════════════════════════════════════════════════════════════
//
// spa_unify Stage 4.  What used to live here was `mgf` (three MGF vectors plus
// a quotient temporary), a private Family-B Newton iteration `fastGetRoot`, and
// `getProbSpa`, a private and completely unguarded copy of the
// Barndorff-Nielsen modified signed root.  All three are gone.  The root finder
// is spa::solveSaddlepoint and the tail is spa::bnTail / spa::bnTailLog; what
// stays is the CGF, which is SPAGRM's own (tier 3, spagrm_cgf.hpp).
//
// Six behavioural consequences, in decreasing order of importance:
//
//   D6 + D1 jointly.  The old kernel computed
//
//       w = sgn(zeta)·sqrt(2·(zeta·Score − CGF0))
//       v = zeta·sqrt(CGF2)
//       u = w + log(v/w)/w
//
//       with NO isfinite(zeta) check, no zeta·Score − CGF0 ≥ 0 check, no
//       CGF2 > 0 check, no w ≠ 0 check and no v/w > 0 check.  01_findings.md
//       describes D1's observable failure as the guard `k2_total <= 0` firing,
//       but that guard belongs to WtCoxG: grep confirms the ONLY test of CGF2
//       anywhere in this file was `std::isinf(CGF2)` inside the root finder.  A
//       non-positive CGF2 therefore reached `std::sqrt` bare and produced NaN,
//       which math::pnorm forwards silently (math_helper.hpp:70).  D1 and D6
//       are one defect jointly here, not two independent ones, and the joint
//       repair is the shared guard set in spa::bnTail plus the mean-centred
//       K'' in spagrm_cgf.hpp.
//
//       Nor is the damage merely "a conservative p = 2.0" as D6 states.  With
//       CGF2 == 0 exactly, v = 0, log(v/w) = −inf and u = ∓inf, math::pnorm's
//       non-finite branch returns 1.0 on BOTH tails, and the sum 2.0 reached
//       the P column unclamped.  math::zFromPval(2.0, ·) then clamps inside
//       qnorm and reports |Z| ≈ 8.2, so the row was internally inconsistent —
//       P = 2 alongside a genome-wide-significant Z — rather than conservative.
//       spa::combineTails clamps into [0, 1] and, more to the point, reports NA
//       with a named status instead.
//
//   D5.  Neither `fastGetRoot` nor `getProbSpa` had a convergence flag at all
//       (their signatures return a bare double), so exhausting the 49-iteration
//       budget was indistinguishable from converging.  spa::Saddle carries a
//       Status that is part of the return value and cannot be dropped, and it
//       reaches the user in the SPA_STATUS column with P set to NA.
//
//   D7.  The upper tail was called with a bare literal 1e-4 and the lower tail
//       with the configured `s.tol` (1e-6), so the two probabilities being
//       summed were not computed to the same accuracy and `--spagrm-tol` did
//       not do what its name implies.  Both tails now use `s.tol`, and the
//       criterion is relative (|K'(t) − s| ≤ tol·sqrt(K''(t))) rather than an
//       absolute test on the pending step.  This is a real numeric change and
//       is the point of the stage.
//
//   P6.  The CGF is now evaluated in ONE pass per term class.  `getProbSpa`
//       walked the workspace three times per tail and `fastGetRoot` four times
//       per Newton iteration; MgfWorkspace's MGF0 / MGF1 / MGF2 / temp vectors
//       (8 × nOutlier doubles per clone) are deleted outright.
//
//   Bracketing and the terminal cumulants.  The old finder maintained no
//       bracket — its two "bisection" loops only halved the step — and
//       evaluated the terminal K'' at the post-step abscissa only because it
//       re-ran `mgf` after the loop.  The solver expands a bracket until the
//       residual changes sign, bisects whenever Newton misbehaves, and
//       evaluates the terminal cumulants at the root it returns.
//
// ── The lower-tail initial abscissa ──────────────────────────────────
//
// The two Family-B implementations did not agree here, which neither the audit
// nor its cross-check flagged: SPAsqr used `zeta2 = -zeta1`, minus the
// data-dependent upper-tail guess, while SPAGRM used
// `zeta2 = -std::abs(s.zeta)`, i.e. minus the fixed configured constant 0.01.
// Unification has to choose one.  This stage takes SPAsqr's convention,
// `init_lower = -init_upper = -min(|Score_adj|/Score_var, 1.2)`, for two
// reasons:
//
//   * It is an estimate of the quantity it initializes.  The first-order
//     saddlepoint is zeta ≈ s/K''(0), and Score_var is the nominal K''(0), so
//     |Score_adj|/Score_var IS that estimate; 0.01 estimates nothing.  The
//     lower-tail problem is the upper-tail problem at s → −s, so its
//     first-order root is the negation.
//   * It makes the two tails cost the same.  Starting the lower tail at −0.01
//     when the root is near −1 forced the bracket expansion to walk outward
//     from almost the origin every time.
//
// With a bracketing solver the choice cannot change the converged root, only
// the work needed to reach it, so this is an efficiency and symmetry decision
// rather than a numeric one.  SPAGRMClass::SharedData::zeta and the
// `zeta` constructor parameter are deleted with it; nsGRMNull::ZETA_DEFAULT is
// no longer referenced.

spa::TwoSided twoSidedSpa(
    const spagrm_cgf::Context &cgf,
    double absScore,
    double initZeta,
    double tol,
    double zNorm
) {
    double p[2], logp[2];
    spa::Status st[2];

    for (int k = 0; k < 2; ++k) {
        const bool lowerTail = (k == 1);
        const double s = lowerTail ? -absScore : absScore;

        spa::SolveOpts opt;
        opt.init = lowerTail ? -initZeta : initZeta;
        opt.scoreSign = s;   // only the sign is read
        opt.rtol = tol;      // D7: the configured tolerance reaches both tails

        const spa::Saddle sd = spa::solveSaddlepoint(
            s,
            [&](double t) { return spagrm_cgf::k12(t, cgf, s); },
            [&](double t) { return spagrm_cgf::kFull(t, cgf, s); },
            opt);

        spa::Status stLin = spa::Status::SpaOk;
        spa::Status stLog = spa::Status::SpaOk;
        p[k]    = spa::bnTail(sd.zeta, s, sd.K0, sd.K2, lowerTail, stLin);
        logp[k] = spa::bnTailLog(sd.zeta, s, sd.K0, sd.K2, lowerTail, stLog);
        st[k] = spa::worseStatus(sd.status, spa::worseStatus(stLin, stLog));
    }

    // spa::assemble owns the splice and the D5 fallback rule: P from the
    // linear assembly, -log10(P) from the log-domain one, and on a saddlepoint
    // failure in either tail both replaced by the two-sided normal tail at
    // zNorm under a status naming the substitution.
    return spa::assemble(p[0], p[1], logp[0], logp[1], st[0], st[1], zNorm);
}

} // namespace nsSPAGRM

// ══════════════════════════════════════════════════════════════════════
// SPAGRMClass
// ══════════════════════════════════════════════════════════════════════

SPAGRMClass::SPAGRMClass(
    Eigen::VectorXd resid,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double R_GRM_R_TwoSubjOutlier,
    double R_GRM_R,
    std::vector<double> MAF_interval,
    nsSPAGRM::FamilyData fam,
    double SPA_Cutoff,
    double tol
)
{
    auto sd = std::make_shared<SharedData>();
    sd->resid = std::move(resid);
    sd->resid_unrelated_outliers = std::move(fam.resid_unrelated_outliers);
    sd->sum_unrelated_outliers2 = sd->resid_unrelated_outliers.squaredNorm();
    sd->sum_R_nonOutlier = sum_R_nonOutlier;
    sd->R_GRM_R_nonOutlier = R_GRM_R_nonOutlier;
    sd->R_GRM_R_TwoSubjOutlier = R_GRM_R_TwoSubjOutlier;
    sd->R_GRM_R = R_GRM_R;
    sd->resid_sum = sd->resid.sum();
    sd->MAF_interval = std::move(MAF_interval);
    sd->TwoSubj_resid_list = std::move(fam.twoSubj_resid);
    sd->TwoSubj_rho_list = std::move(fam.twoSubj_rho);
    sd->ThreeSubj_standS_list = std::move(fam.threeSubj_standS);
    sd->ThreeSubj_CLT_list = std::move(fam.threeSubj_CLT);
    sd->SPA_Cutoff = SPA_Cutoff;
    sd->tol = tol;
    m_shared = std::move(sd);
    rebuildScratch();
}

SPAGRMClass::SPAGRMClass(const SPAGRMClass &o)
    : m_shared(o.m_shared)
{
    rebuildScratch();
}

void SPAGRMClass::rebuildScratch() {
    const auto &s = *m_shared;
    const int n3 = static_cast<int>(s.ThreeSubj_standS_list.size());
    m_threeSubj_scratch.resize(n3);
    size_t widest = 0;
    for (int i = 0; i < n3; ++i) {
        m_threeSubj_scratch[i].stand_S = s.ThreeSubj_standS_list[i];
        m_threeSubj_scratch[i].arr_prob.resize(s.ThreeSubj_standS_list[i].size());
        widest = std::max(widest, s.ThreeSubj_standS_list[i].size());
    }
    // One class-3 family's tilted weights, so the mean-centred K'' needs no
    // second exponential pass.  At most 3^MAX_NUM_IN_FAM = 243 doubles.
    m_cgfScratch.assign(widest, 0.0);
}

// ══════════════════════════════════════════════════════════════════════
// Per-marker two-sided p-value
// ══════════════════════════════════════════════════════════════════════
//
// The GVec-taking overload `SPAGRMClass::getMarkerPval` is deleted.  It was a
// ~70-line near-verbatim copy of this function differing only in computing
// `Score` itself, and grep over src/ and tests/ found no call site: the three
// SPAGRMMethod entry points and both SAGELD entry points all go through
// getMarkerPvalFromScore.  It also carried its own copies of the D6 and D7
// defects, so keeping it would have left two of the seven copies of D1's
// cancelling block alive in a function nothing calls.

spa::TwoSided SPAGRMClass::getMarkerPvalFromScore(
    double Score,
    double altFreq,
    double &zScore,
    double *outScoreVar
) {
    const auto &s = *m_shared;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double MAF = std::min(altFreq, 1.0 - altFreq);
    const double G_var = 2.0 * MAF * (1.0 - MAF);
    const double Score_var = G_var * s.R_GRM_R;

    if (outScoreVar) *outScoreVar = Score_var;

    // Monomorphic or degenerate marker — no statistic exists, so there is no
    // saddlepoint to attempt and nothing to report but the reason.
    if (Score_var <= 0.0 || MAF <= 0.0) {
        zScore = 0.0;
        return spa::TwoSided{nan, nan, spa::Status::NaNoTest};
    }

    zScore = Score / std::sqrt(Score_var);

    if (!std::isfinite(zScore))
        return spa::TwoSided{nan, nan, spa::Status::NaNoTest};

    // Normal approximation below the SPA cutoff.  spa::normalBranch routes the
    // tail through Boost's complement, which is the same double the old
    // `2 * pnorm(|z|, upper)` produced, and adds the log-domain -log10(p).
    if (std::abs(zScore) <= s.SPA_Cutoff) return spa::normalBranch(zScore);

    int order2 =
        static_cast<int>(std::lower_bound(s.MAF_interval.begin(), s.MAF_interval.end(), MAF) - s.MAF_interval.begin());
    const int nBins = static_cast<int>(s.MAF_interval.size());
    if (order2 <= 0) order2 = 1;
    if (order2 >= nBins) order2 = nBins - 1;
    const int order1 = order2 - 1;
    const double MAF_ratio = (s.MAF_interval[order2] - MAF) / (s.MAF_interval[order2] - s.MAF_interval[order1]);
    const double one_minus_mr = 1.0 - MAF_ratio;

    // Refresh each class-3 family's probability table at this marker's MAF, and
    // accumulate its untilted variance for the EmpVar rescaling below.  This is
    // the one part of the CGF input that is per-marker rather than per-tail.
    double Var_ThreeOutlier = 0.0;
    const int n3 = static_cast<int>(m_threeSubj_scratch.size());
    for (int i = 0; i < n3; ++i) {
        const double *c1 = s.ThreeSubj_CLT_list[i].col(order1).data();
        const double *c2 = s.ThreeSubj_CLT_list[i].col(order2).data();
        const double *ss = m_threeSubj_scratch[i].stand_S.data();
        double *ap = m_threeSubj_scratch[i].arr_prob.data();
        const size_t sz = m_threeSubj_scratch[i].stand_S.size();
        double s1 = 0.0, s2 = 0.0;
        for (size_t k = 0; k < sz; ++k) {
            ap[k] = MAF_ratio * c1[k] + one_minus_mr * c2[k];
            const double tmp = ss[k] * ap[k];
            s1 += tmp;
            s2 += ss[k] * tmp;
        }
        Var_ThreeOutlier += s2 - s1 * s1;
    }

    const double EmpVar =
        G_var * (s.R_GRM_R_nonOutlier + s.sum_unrelated_outliers2 + s.R_GRM_R_TwoSubjOutlier) + Var_ThreeOutlier;
    // EmpVar is a quadratic form over a thresholded sparse GRM plus a tabulated
    // variance and is NOT guaranteed non-negative; a negative EmpVar makes
    // Score_adj NaN.  That is not silently absorbed: the solver's first
    // evaluation is non-finite, so it returns Status::NaNoTest and the row
    // reports NA with a reason rather than an unexplained NaN.
    const double Score_adj = Score * std::sqrt(EmpVar / Score_var);

    // First-order saddlepoint estimate zeta ≈ s/K''(0), capped; the lower tail
    // starts at its negation (see the convention note above twoSidedSpa).
    const double initZeta = std::min(std::abs(Score_adj) / Score_var, 1.2);

    spagrm_cgf::Context cgf;
    cgf.outlierResid = s.resid_unrelated_outliers.data();
    cgf.nOutlier = static_cast<int>(s.resid_unrelated_outliers.size());
    cgf.twoResid = s.TwoSubj_resid_list.data();
    cgf.twoRho = s.TwoSubj_rho_list.data();
    cgf.nTwo = static_cast<int>(s.TwoSubj_resid_list.size());
    cgf.three = m_threeSubj_scratch.data();
    cgf.nThree = n3;
    cgf.scratch = m_cgfScratch.data();
    cgf.MAF = MAF;
    cgf.mean = 2.0 * MAF * s.sum_R_nonOutlier;
    cgf.var = 2.0 * MAF * (1.0 - MAF) * s.R_GRM_R_nonOutlier;

    return nsSPAGRM::twoSidedSpa(cgf, std::abs(Score_adj), initZeta, s.tol, zScore);
}

// ══════════════════════════════════════════════════════════════════════
// Helper — build edges, components, family entries from GRM entries
// ══════════════════════════════════════════════════════════════════════
namespace {

struct GRMTopology {
    std::unordered_set<uint32_t> singletonSet;
    std::vector<std::vector<uint32_t> > families;
    std::vector<std::vector<SparseGRM::Entry> > familyEntries;
};

GRMTopology buildTopology(
    uint32_t N,
    const std::vector<SparseGRM::Entry> &entries
) {
    std::vector<std::pair<uint32_t, uint32_t> > edges;
    {
        std::unordered_set<uint64_t> seen;
        for (const auto &e : entries) {
            if (e.row == e.col) continue;
            uint32_t lo = std::min(e.row, e.col);
            uint32_t hi = std::max(e.row, e.col);
            uint64_t key = (static_cast<uint64_t>(lo) << 32) | hi;
            if (seen.insert(key).second) edges.push_back({lo, hi});
        }
    }
    auto components = nsGRMNull::getComponents(N, edges);

    GRMTopology topo;
    std::vector<std::vector<uint32_t> > singletons;
    for (auto &comp : components) {
        if (comp.size() == 1)
            singletons.push_back(std::move(comp));
        else
            topo.families.push_back(std::move(comp));
    }
    for (const auto &s : singletons)
        topo.singletonSet.insert(s[0]);

    std::unordered_map<uint32_t, size_t> subjToFamily;
    for (size_t fi = 0; fi < topo.families.size(); ++fi)
        for (uint32_t idx : topo.families[fi])
            subjToFamily[idx] = fi;

    topo.familyEntries.resize(topo.families.size());
    for (const auto &e : entries) {
        auto it = subjToFamily.find(e.row);
        if (it != subjToFamily.end()) topo.familyEntries[it->second].push_back(e);
    }
    return topo;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// runSPAGRM — entry point
// ══════════════════════════════════════════════════════════════════════

void runSPAGRM(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    bool controlOutlier,
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
        infoMsg("SPAGRM: fitting %s null model for %zu phenotype(s)",
                nullmodel::regressionModelName(regModel), phenoSpecs.size());
    }

    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);

    // --longitudinal: fit Y ~ X + (1|IID) on the long-format pheno file and
    // obtain R_G before SubjectData is built.  The kept set is GRM ∩ keep −
    // remove, matching SAGELD pheno-mode.
    nsLongitudinal::LongResidResult Lr;
    std::unordered_set<std::string> grmIDs;
    if (longitudinal) {
        grmIDs = SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, famIIDs);
        const auto outcomeNames = nsLongitudinal::splitOutcomeNames(phenoNameSpec);
        const auto kept = nsLongitudinal::buildKeptSet(keepFile, removeFile, famIIDs, grmIDs);
        infoMsg("SPAGRM: fitting longitudinal Y ~ X + (1|IID) for %zu outcome(s)",
                outcomeNames.size());
        Lr = nsLongitudinal::computeLongitudinalResid(
            phenoFile, outcomeNames, covarNames, famIIDs, kept);
    }

    SubjectData sd(std::move(famIIDs));
    if (longitudinal) {
        sd.setKeepSubjects(
            std::unordered_set<std::string>(Lr.usedIIDs.begin(), Lr.usedIIDs.end()));
    } else if (fitPath) {
        sd.loadPhenoFile(phenoFile, nullmodel::columnsNeeded(phenoSpecs));
    } else {
        sd.loadResidOne(phenoFile, residNames);
    }
    if (!covarFile.empty() && !longitudinal) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    if (longitudinal)
        sd.setGrmSubjects(std::move(grmIDs));
    else
        sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();
    const uint32_t N = sd.nUsed();
    infoMsg("  %u subjects in union mask", N);

    // ---- Inject the longitudinal residual R_G as the per-subject residual ----
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
    }

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

    auto subjIDs = sd.usedIIDs();
    auto subjIdMap = text::buildIIDMap(subjIDs);

    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, subjIDs, sd.famIIDs());
    infoMsg("Sparse GRM: %zu entries (diagonal + off-diag)", grm.nnz());

    infoMsg("Loading pairwise IBD: %s", pairwiseIBDFile.c_str());
    auto ibdEntries = nsGRMNull::loadIndexedIBD(pairwiseIBDFile, subjIdMap);
    infoMsg("Loaded %zu IBD records", ibdEntries.size());

    const auto &allEntries = grm.entries();
    const auto &grmDiag = grm.diagonal();

    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);
    infoMsg("Genotype data: %u markers, %u subjects", genoData->nMarkers(), genoData->nSubjUsed());

    // ---- Build per-phenotype tasks ----
    auto phenoInfos = sd.buildPerColumnMasks();
    const int K = sd.residOneCols();
    if (K > 1) infoMsg("Multi-column residual file: %d phenotypes", K);

    std::vector<PhenoTask> tasks(K);

    if (K == 1) {
        // Single phenotype: use union structures directly
        auto topo = buildTopology(N, allEntries);
        auto ibdPairMap = nsGRMNull::buildIBDPairMap(ibdEntries);
        infoMsg("Building null model for '%s'...", phenoInfos[0].name.c_str());
        SPAGRMClass sg = nsGRMNull::buildSPAGRMNullModel(sd.residuals(), N, topo.singletonSet, grmDiag, topo.families,
                                                         topo.familyEntries, allEntries, ibdEntries, ibdPairMap,
                                                         spaCutoff, minMafCutoff, minMacCutoff,
                                                         outlierIqrRatio, controlOutlier, nthreads);
        tasks[0].phenoName = phenoInfos[0].name;
        tasks[0].method = std::make_unique<SPAGRMMethod>(std::move(sg));
        tasks[0].unionToLocal = phenoInfos[0].unionToLocal;
        tasks[0].nUsed = phenoInfos[0].nUsed;
        infoMsg("  Phenotype '%s': %u subjects", phenoInfos[0].name.c_str(), phenoInfos[0].nUsed);
    } else {
        for (int rc = 0; rc < K; ++rc) {
            const auto &pi = phenoInfos[rc];
            const auto &u2l = pi.unionToLocal;
            const uint32_t Np = pi.nUsed;

            // Re-index GRM entries to per-phenotype dense indices
            std::vector<SparseGRM::Entry> pEntries;
            std::vector<double> pDiag(Np, 1.0);
            for (const auto &e : allEntries) {
                uint32_t li = u2l[e.row], lj = u2l[e.col];
                if (li != UINT32_MAX && lj != UINT32_MAX) {
                    pEntries.push_back({li, lj, e.value});
                    if (li == lj) pDiag[li] = e.value;
                }
            }

            auto topo = buildTopology(Np, pEntries);

            // Re-index IBD entries
            std::vector<nsGRMNull::IndexedIBD> pIbd;
            for (const auto &ibd : ibdEntries) {
                uint32_t li = u2l[ibd.idx1], lj = u2l[ibd.idx2];
                if (li != UINT32_MAX && lj != UINT32_MAX) pIbd.push_back({li, lj, ibd.pa, ibd.pb, ibd.pc});
            }
            auto pIbdMap = nsGRMNull::buildIBDPairMap(pIbd);

            Eigen::VectorXd phenoResid = extractPhenoVec(sd.residMatrix().col(rc), pi);

            infoMsg("Building null model for '%s' (%u subjects)...", pi.name.c_str(), Np);
            SPAGRMClass sg = nsGRMNull::buildSPAGRMNullModel(phenoResid, Np, topo.singletonSet, pDiag, topo.families,
                                                             topo.familyEntries, pEntries, pIbd, pIbdMap, spaCutoff,
                                                             minMafCutoff, minMacCutoff,
                                                             outlierIqrRatio, controlOutlier, nthreads);

            tasks[rc].phenoName = pi.name;
            tasks[rc].method = std::make_unique<SPAGRMMethod>(std::move(sg));
            tasks[rc].unionToLocal = pi.unionToLocal;
            tasks[rc].nUsed = pi.nUsed;
            infoMsg("  Phenotype '%s': %u subjects", pi.name.c_str(), pi.nUsed);
        }
    }

    infoMsg("Running SPAGRM marker tests (%d thread(s), %d phenotype(s))...", nthreads, K);
    multiPhenoEngine(*genoData, tasks, outPrefix, "SPAGRM", compression, compressionLevel, nthreads,
                     missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);
}

// ══════════════════════════════════════════════════════════════════════
// runSPAGRMLoco — LOCO orchestration
// ══════════════════════════════════════════════════════════════════════
//
// The sparse-GRM topology (families / singletons / re-indexed entries) and the
// IBD structures are residual-independent, so they are built ONCE.  Per
// chromosome the null model is refit with that chromosome's LOCO PGS appended as
// a covariate column, and buildSPAGRMNullModel is re-run to rebuild R_GRM_R and
// the outlier / family / Chow-Liu tables from the refreshed residuals.
// SPAGRMMethod owns its data (constructed by move), so there is no reference-
// lifetime discipline to observe here.

void runSPAGRMLoco(
    const std::string &phenoFile,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &predListFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    bool controlOutlier,
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
    const std::vector<std::string> &covarNames
) {
    if (phenoNameSpec.empty())
        throw std::runtime_error(
            "SPAGRM-LOCO requires --pheno-name (an in-process null-model fit); "
            "precomputed residuals cannot be refit per chromosome");

    nullmodel::RegressionModel regModel =
        nullmodel::parseRegressionModel(regressionModelStr);
    std::vector<nullmodel::PhenoSpec> phenoSpecs =
        nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
    const int K = static_cast<int>(phenoSpecs.size());
    std::vector<std::string> specNames(K);
    for (int k = 0; k < K; ++k) specNames[k] = phenoSpecs[k].name;

    infoMsg("SPAGRM-LOCO: fitting %s null model for %d phenotype(s), "
            "per-chromosome LOCO PGS as covariate",
            nullmodel::regressionModelName(regModel), K);

    validatePredListPhenos(predListFile, specNames);

    // ---- Load pheno/covar/GRM data ----
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    sd.loadPhenoFile(phenoFile, nullmodel::columnsNeeded(phenoSpecs));
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();
    const uint32_t N = sd.nUsed();
    infoMsg("  %u subjects in union mask", N);

    // ---- Base covariate matrix (no PGS) ----
    Eigen::MatrixXd baseCovar;
    if (!covarNames.empty()) baseCovar = sd.getColumns(covarNames);
    else if (sd.hasCovar()) baseCovar = sd.covar();
    else baseCovar.resize(N, 0);

    // ---- One base fit to establish per-phenotype masks ----
    {
        nullmodel::EngineOptions eo;
        eo.nthreads = nthreads;
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
        if (regModel != nullmodel::RegressionModel::Auto) {
            specModel[k] = regModel;
        } else if (nullmodel::isCoxSpec(phenoSpecs[k])) {
            specModel[k] = nullmodel::RegressionModel::Cox;
        } else {
            auto info = nullmodel::inferModelFromColumn(
                sd.getColumn(phenoSpecs[k].yColumn), phenoSpecs[k].yColumn,
                sd.usedIIDs());
            specModel[k] = info.model;
        }
    }

    // ---- Load GRM, IBD, genotype data ----
    auto subjIDs = sd.usedIIDs();
    auto subjIdMap = text::buildIIDMap(subjIDs);
    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, subjIDs, sd.famIIDs());
    infoMsg("Sparse GRM: %zu entries (diagonal + off-diag)", grm.nnz());
    infoMsg("Loading pairwise IBD: %s", pairwiseIBDFile.c_str());
    auto ibdEntries = nsGRMNull::loadIndexedIBD(pairwiseIBDFile, subjIdMap);
    infoMsg("Loaded %zu IBD records", ibdEntries.size());
    const auto &allEntries = grm.entries();
    const auto &grmDiag = grm.diagonal();
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    // ---- Per-phenotype residual-independent GRM/IBD structures (built once) ----
    struct PhenoNull {
        uint32_t Np = 0;
        GRMTopology topo;
        std::vector<double> diag;
        std::vector<SparseGRM::Entry> entries;
        std::vector<nsGRMNull::IndexedIBD> ibd;
        std::unordered_map<uint64_t, uint32_t> ibdMap;
    };
    std::vector<PhenoNull> pn(K);
    for (int rc = 0; rc < K; ++rc) {
        const auto &pi = phenoInfos[rc];
        if (K == 1) {
            pn[rc].Np = N;
            pn[rc].entries = allEntries;
            pn[rc].diag = grmDiag;
            pn[rc].topo = buildTopology(N, allEntries);
            pn[rc].ibd = ibdEntries;
            pn[rc].ibdMap = nsGRMNull::buildIBDPairMap(ibdEntries);
        } else {
            const auto &u2l = pi.unionToLocal;
            const uint32_t Np = pi.nUsed;
            std::vector<SparseGRM::Entry> pEntries;
            std::vector<double> pDiag(Np, 1.0);
            for (const auto &e : allEntries) {
                uint32_t li = u2l[e.row], lj = u2l[e.col];
                if (li != UINT32_MAX && lj != UINT32_MAX) {
                    pEntries.push_back({li, lj, e.value});
                    if (li == lj) pDiag[li] = e.value;
                }
            }
            std::vector<nsGRMNull::IndexedIBD> pIbd;
            for (const auto &ibd : ibdEntries) {
                uint32_t li = u2l[ibd.idx1], lj = u2l[ibd.idx2];
                if (li != UINT32_MAX && lj != UINT32_MAX)
                    pIbd.push_back({li, lj, ibd.pa, ibd.pb, ibd.pc});
            }
            pn[rc].Np = Np;
            pn[rc].topo = buildTopology(Np, pEntries);
            pn[rc].entries = std::move(pEntries);
            pn[rc].diag = std::move(pDiag);
            pn[rc].ibdMap = nsGRMNull::buildIBDPairMap(pIbd);
            pn[rc].ibd = std::move(pIbd);
        }
    }

    // ---- Per-phenotype non-missing masks (union space, chromosome-invariant) ----
    std::vector<std::vector<bool> > nonMissing(K, std::vector<bool>(N, false));
    for (int k = 0; k < K; ++k)
        for (uint32_t i = 0; i < N; ++i)
            nonMissing[k][i] = (phenoInfos[k].unionToLocal[i] != UINT32_MAX);

    // ---- Load LOCO predictions ----
    LocoData loco = LocoData::load(predListFile, specNames, sd.usedIIDs(), sd.famIIDs());
    auto locoChroms = loco.availableChromosomes();
    infoMsg("LOCO: %zu chromosomes available across all phenotypes", locoChroms.size());

    nullmodel::EngineOptions eo1;
    eo1.nthreads = 1;

    auto buildTasks = [&](const std::string &chr, std::vector<PhenoTask> &tasks) {
        tasks.resize(K);
        for (int rc = 0; rc < K; ++rc) {
            Eigen::MatrixXd covarUnion_k = appendLocoCovariate(
                loco, specNames[rc], chr, baseCovar, nonMissing[rc], "SPAGRM-LOCO");
            auto fits1 = nullmodel::fitAll(
                sd, {phenoSpecs[rc]}, specModel[rc], covarUnion_k, eo1);
            Eigen::VectorXd phenoResid =
                extractPhenoVec(fits1[0].residuals, phenoInfos[rc]);

            SPAGRMClass sg = nsGRMNull::buildSPAGRMNullModel(
                phenoResid, pn[rc].Np, pn[rc].topo.singletonSet, pn[rc].diag,
                pn[rc].topo.families, pn[rc].topo.familyEntries, pn[rc].entries,
                pn[rc].ibd, pn[rc].ibdMap, spaCutoff, minMafCutoff, minMacCutoff,
                outlierIqrRatio, controlOutlier, nthreads);

            tasks[rc].phenoName = phenoInfos[rc].name;
            tasks[rc].method = std::make_unique<SPAGRMMethod>(std::move(sg));
            tasks[rc].unionToLocal = phenoInfos[rc].unionToLocal;
            tasks[rc].nUsed = phenoInfos[rc].nUsed;
        }
    };

    infoMsg("SPAGRM-LOCO: starting LOCO association (%d phenotype(s), %zu chroms, %d threads)",
            K, locoChroms.size(), nthreads);
    locoEngine(
        *genoData, locoChroms, specNames, buildTasks, outPrefix, "SPAGRM",
        compression, compressionLevel, nthreads,
        missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);
}
