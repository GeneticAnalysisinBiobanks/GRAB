// spacox.cpp — SPACox full implementation

#include "spacox/spacox.hpp"
#include "engine/loco.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "spagrm/longitudinal_resid.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/null_model.hpp"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>

// ======================================================================
// Worker-local SPA scratch
// ======================================================================
//
// The SPA branch of getMarkerPvalCore needs a handful of transient
// buffers.  These were previously per-clone fields on SPACoxMethod, which
// meant the multi-phenotype engine paid (K phenotypes × T worker threads)
// times their footprint.  Promoting them to a translation-unit-local
// thread_local struct lets all K phenotype-clones in one worker share
// one set of buffers — the engine spawns at most T concurrent threads,
// so total scratch footprint is bounded by T × (≈ 20N bytes) regardless
// of K.  Resize() is idempotent when the requested length matches the
// existing capacity, so phenotype switching inside a thread is O(1)
// after the first call.
//
// `weights` (P1) replaces the former length-N adjGNorm array on the
// stage-1 path: for hard-called input it holds at most four (weight,
// multiplicity) pairs, so the stage-1 CGF no longer materializes a
// per-subject vector at all.  `adjGNorm` survives for stage 2, whose
// covariate-adjusted weights are a projection and genuinely per-subject,
// and `nzSet` survives for the projection itself, where it is rebuilt on
// entry to stage 2 rather than on every marker.
namespace {
struct SPACoxScratch {
    spacox_cgf::GenoWeights weights;
    Eigen::VectorXd adjGNorm;
    Eigen::VectorXd adjGVec;
    std::vector<uint32_t> nzSet;
};

thread_local SPACoxScratch tlScratch;
}  // namespace

// ======================================================================
// DesignMatrix
// ======================================================================

DesignMatrix::DesignMatrix(const Eigen::MatrixXd &X)
    : m_X(X)
{
    const int nCols = static_cast<int>(m_X.cols());
    m_tX = m_X.transpose();
    Eigen::MatrixXd XtX = m_tX * m_X;
    m_XinvXX = m_X * XtX.ldlt().solve(Eigen::MatrixXd::Identity(nCols, nCols));
}

void DesignMatrix::adjustGenotype(
    const double *G,
    const uint32_t *nzIdx,
    int nNz,
    Eigen::Ref<Eigen::VectorXd> adjG
) const {

    const int p = nCols();
    const int N = nRows();

    Eigen::VectorXd tX_g = Eigen::VectorXd::Zero(p);
    for (int k = 0; k < nNz; ++k) {
        uint32_t j = nzIdx[k];
        double gj = G[j];
        tX_g.noalias() += gj * m_tX.col(j);
    }

    Eigen::Map<const Eigen::VectorXd> gVec(G, N);
    adjG.noalias() = gVec - m_XinvXX * tX_g;
}

// ======================================================================
// SPACoxMethod — construction & clone
// ======================================================================

SPACoxMethod::SPACoxMethod(
    const Eigen::VectorXd &residuals,
    double varResid,
    const CumulantTable &cumul,
    const DesignMatrix &design,
    double pvalCovAdjCut,
    double spaCutoff
)
    : m_resid(residuals),
      m_varResid(varResid),
      m_cumul(cumul),
      m_design(design),
      m_N(static_cast<int>(residuals.size())),
      m_pvalCovAdjCut(pvalCovAdjCut),
      m_spaCutoff(spaCutoff)
{
}

std::unique_ptr<MethodBase> SPACoxMethod::clone() const {
    // The clone shares with the master:
    //   * m_resid   (const reference to N-vector of residuals, owned by runSPACox)
    //   * m_cumul   (const reference to L=10000 CGF table, owned by runSPACox)
    //   * m_design  (const reference to N×p design / projection matrices,
    //                deduplicated across phenotypes by missingness pattern
    //                inside runSPACox)
    //   * the SPA scratch (weights / adjGNorm / adjGVec / nzSet) is
    //     thread_local in this translation unit, so K phenotype-clones in one
    //     worker share a single set of buffers.
    // The clone only owns its small scalar metadata (m_varResid, m_N,
    // m_pvalCovAdjCut, m_spaCutoff) and the references themselves — the
    // entire clone is therefore on the order of one cache line.
    return std::make_unique<SPACoxMethod>(m_resid, m_varResid, m_cumul, m_design, m_pvalCovAdjCut, m_spaCutoff);
}

// ======================================================================
// Two-sided saddlepoint p-value — on the shared solver and tail kernel
// ======================================================================
//
// spa_unify Stage 3.  What used to live here was a private Family-A Newton
// iteration (fastGetRootK1) plus a private copy of the Barndorff-Nielsen
// modified signed root.  Both are gone; the root finder is
// spa::solveSaddlepoint and the tail is spa::bnTailLog.  What
// stays is the CGF, which is SPACox's own (tier 3, see spacox_cgf.hpp).
//
// Four behavioural consequences, in decreasing order of importance:
//
//   D5.  The old finder computed a `converge` flag, returned it, and no
//        caller ever read it — `grep -n converge src/spacox/*` found the
//        field written and never consulted.  A saddlepoint that exhausted
//        its 100 iterations therefore produced an ordinary-looking finite
//        p-value with no NaN, no flag and no warning.  spa::Saddle carries a
//        Status that is part of the return value and cannot be dropped, and
//        it now reaches the user in the SPA_STATUS column with P set to NA.
//
//   The stale-K'' pairing.  The old code evaluated K at the returned root
//        but reused `rr.K2` from the finder's last completed iteration.  On
//        the converged path those coincide, because Family A breaks on
//        |diffX| < tol *before* applying the step; on the maxiter-exhausted
//        path they do not, so the guard `k2val <= 0.0` was testing curvature
//        at a different abscissa than the K it was paired with.  The solver
//        evaluates the terminal cumulants at the root it returns, so the
//        mismatch cannot recur.  (This was never an independent defect —
//        it was only reachable on the same path D5 was silent about.)
//
//   Tolerance.  The old break test was |K'/K''| < 1e-3, absolute and
//        applied to the *step*, so the returned abscissa carried an unapplied
//        correction of up to 1e-3.  The solver's test is relative and applied
//        to the residual after the step: |K'(t) - s| <= 1e-6 * sqrt(K''(t)).
//        Since w is stationary at the root and only v = zeta*sqrt(K'') carries
//        first-order sensitivity, the p-value effect is ~4e-4 in -log10 p, but
//        it is a real change and is the dominant source of the Stage 3 diff.
//
//   Bracketing.  The old finder had no bracket and its damping heuristic
//        could stall; the solver expands a bracket until the residual changes
//        sign and bisects whenever Newton misbehaves, so non-convergence is
//        now reported rather than silently converged-to-something.
//
// The initial abscissa is s itself.  SPACox standardizes the weights so that
// sum a_i^2 = 1/varResid and hence K''(0) ~ 1, which makes the first-order
// saddlepoint estimate zeta ~ s/K''(0) = s.  The sign constraint
// sgn(zeta) == sgn(s) is passed as well: K'(0) is the (near-zero) mean of the
// tilted score, so the root lies on the same side as s, and restricting the
// search to that half-line is the one safeguard the original Family B finders
// had that was worth keeping.

namespace {

// Two-sided assembly shared by both stages.  `cgf(t, s)` must return
// spa::K012{K(t), K'(t) - s, K''(t)}.
template <class Cgf>
spa::Result spaTwoSided(const Cgf &cgf, double absZ, double zNorm) {
    double logp[2];
    spa::Status st[2];

    for (int k = 0; k < 2; ++k) {
        const bool lowerTail = (k == 1);
        const double s = lowerTail ? -absZ : absZ;

        spa::SolveOpts opt;
        opt.init = s;
        opt.scoreSign = s;   // only the sign is read

        const spa::Saddle sd = spa::solveSaddlepoint(
            s, [&](double t) { return cgf(t, s); }, opt);

        spa::Status stLog = spa::Status::SpaOk;
        logp[k] = spa::bnTailLog(sd.zeta, s, sd.K0, sd.K2, lowerTail, stLog);
        st[k] = spa::worseStatus(sd.status, stLog);
    }

    // spa::assemble owns the two-sided combination and the D5 fallback rule:
    // on a saddlepoint failure in either tail the result is replaced by the
    // two-sided normal tail at zNorm, under a status naming the substitution.
    // One tail evaluation per side since log10p_unify Stage 3 deleted the
    // parallel linear tail this loop used to run beside `bnTailLog`.
    return spa::assemble(logp[0], logp[1], st[0], st[1], zNorm);
}

}  // namespace

spa::Result SPACoxMethod::getProbSpaBucketed(
    const spacox_cgf::GenoWeights &gw,
    double absZ,
    double zNorm
) const {
    const spacox_cgf::TableView tv(m_cumul);
    return spaTwoSided(
        [&](double t, double s) { return spacox_cgf::evalBucketed(tv, t, gw, s); },
        absZ, zNorm);
}

spa::Result SPACoxMethod::getProbSpaDense(
    const double *w,
    int n,
    double absZ,
    double zNorm
) const {
    const spacox_cgf::TableView tv(m_cumul);
    return spaTwoSided(
        [&](double t, double s) { return spacox_cgf::evalDense(tv, t, w, n, s); },
        absZ, zNorm);
}

// ======================================================================
// Per-marker p-value (two-stage SPA)
// ======================================================================

spa::Result SPACoxMethod::getMarkerPvalCore(
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    double S,
    double &zScore,
    double &outScoreVar
) {
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // S is pre-computed by the caller as GVec · m_resid.  Inside
    // getResultBatch this is folded over B markers in one GEMV; inside
    // getResultVec it is computed as a single dot product.  Either way
    // the operation is bit-identical to the column-by-column dot used
    // by the original implementation, so the variance loop below is
    // also left in its original per-marker form to preserve bit-equal
    // outputs across the refactor.
    const double twoMAF = 2.0 * altFreq;
    const double *gp = GVec.data();

    // VarS = varResid * Σ(g_i - 2·MAF)²
    //
    // Left as a per-subject reduction on purpose.  The bucket counts P1
    // introduces below could produce Σ(g_i - 2·MAF)² in O(1) as
    // n0·d0² + n1·d1² + n2·d2², but this loop runs for *every* marker while
    // the bucket build runs only inside the SPA branch, and replacing it
    // would perturb Var(S) — hence Z — for all markers including the ~95 %
    // that never reach the saddlepoint.  The loop is a single vectorizable
    // FMA reduction; there is nothing to win here and a whole-output
    // perturbation to lose.
    double sumAdjG2 = 0.0;
    for (int i = 0; i < m_N; ++i) {
        double adj = gp[i] - twoMAF;
        sumAdjG2 += adj * adj;
    }
    double VarS = m_varResid * sumAdjG2;
    outScoreVar = VarS;
    if (VarS <= 0.0) {
        zScore = 0.0;
        return spa::Result{nan, spa::Status::NaNoTest};
    }
    zScore = S / std::sqrt(VarS);

    // Normal approximation if |z| below SPA cutoff.  spa::normalBranch routes
    // the tail through Boost's complement, which is the same double the old
    // `2 * pnorm(-|z|)` produced, and adds the log-domain -log10(p).
    if (std::abs(zScore) <= m_spaCutoff) return spa::normalBranch(zScore);

    // ---- Stage 1: unadjusted SPA over the distinct genotype weights ----
    //
    // P1.  a_i = (g_i - 2·MAF)/sqrt(Var S) is a function of g_i alone, so for
    // hard-called input it takes at most four values.  buildGenoWeights
    // collapses them and keeps a leftover list for genuine dosages, which
    // makes each Newton iteration O(4) instead of O(nNonZero) — the previous
    // code already exploited this for g == 0 (the N0 term) and then walked
    // every non-zero subject individually, one std::atan each.
    double sqrtVarS = std::sqrt(VarS);
    spacox_cgf::buildGenoWeights(gp, m_N, twoMAF, sqrtVarS, tlScratch.weights);

    const double absZ1 = std::abs(zScore);
    spa::Result ts = getProbSpaBucketed(tlScratch.weights, absZ1, zScore);

    // A stage-1 failure (the p-value is NaN) falls through to stage 2 exactly
    // as it did before, when `NaN > m_pvalCovAdjCut` evaluated false: the
    // covariate-adjusted statistic is the better one, so it is worth trying.
    //
    // The comparison stays on the linear scale because --covar-p-threshold is
    // an INPUT parameter and decision D8 keeps input thresholds in linear p.
    // `std::pow(10, -L)` propagates NaN, so the fall-through above is
    // preserved; this is the site the deleted `spa::pFromNegLog10P` served
    // that is not an output column, and it is the reason the conversion is
    // written here rather than shared (log10p_unify Stage 8).
    if (std::pow(10.0, -ts.negLog10p) > m_pvalCovAdjCut) return ts;

    // ---- Stage 2: covariate-adjusted SPA ----
    //
    // adjG is the genotype projected onto the orthogonal complement of the
    // design matrix, so it carries no bucket structure and every subject must
    // be visited.  This runs only for markers whose stage-1 p is at or below
    // --covar-p-threshold (5e-5), so the dense path costs almost nothing over
    // a full scan.  nzSet is rebuilt here rather than on every marker, since
    // the projection is its only consumer.
    tlScratch.nzSet.clear();
    if (static_cast<int>(tlScratch.nzSet.capacity()) < m_N)
        tlScratch.nzSet.reserve(m_N);
    for (int i = 0; i < m_N; ++i)
        if (gp[i] != 0.0) tlScratch.nzSet.push_back(static_cast<uint32_t>(i));
    const int nNz = static_cast<int>(tlScratch.nzSet.size());

    tlScratch.adjGVec.resize(m_N);
    m_design.adjustGenotype(gp, tlScratch.nzSet.data(), nNz, tlScratch.adjGVec);

    VarS = m_varResid * tlScratch.adjGVec.squaredNorm();
    outScoreVar = VarS;
    if (VarS <= 0.0) {
        zScore = 0.0;
        return spa::Result{nan, spa::Status::NaNoTest};
    }
    zScore = S / std::sqrt(VarS);
    sqrtVarS = std::sqrt(VarS);

    tlScratch.adjGNorm.resize(m_N);
    double *anp = tlScratch.adjGNorm.data();
    const double *avp = tlScratch.adjGVec.data();
    for (int i = 0; i < m_N; ++i)
        anp[i] = avp[i] / sqrtVarS;

    return getProbSpaDense(anp, m_N, std::abs(zScore), zScore);
}

// ======================================================================
// Result assembly
// ======================================================================

void SPACoxMethod::pushResult(
    std::vector<double> &out,
    const spa::Result &ts,
    double S,
    double scoreVar
) {
    // LOG10P is the sole p-value column since log10p_unify Stage 8 (D1);
    // a consumer that needs the linear p takes 10^(-LOG10P).
    out.push_back(ts.negLog10p);  // LOG10P
    if (scoreVar > 0.0) {
        const double sd = std::sqrt(scoreVar);
        const double zNorm = S / sd;
        // Z from L, not from P: the linear inversion saturated at
        // |Z| = 37.0470962993612 for every L >= 299.698970 (Stage 4).
        out.push_back(math::zFromNegLog10P(ts.negLog10p, zNorm));  // Z (p-consistent)
        out.push_back(zNorm);                        // Z_Norm (raw score z)
        out.push_back(S / scoreVar);                 // BETA
        out.push_back(1.0 / sd);                     // SE
    } else {
        const double nan = std::numeric_limits<double>::quiet_NaN();
        out.push_back(nan);  // Z
        out.push_back(nan);  // Z_Norm
        out.push_back(nan);  // BETA
        out.push_back(nan);  // SE
    }
    out.push_back(static_cast<double>(static_cast<uint8_t>(ts.status)));
}

// ======================================================================
// MethodBase::getResultVec — called by marker engine per marker
// ======================================================================

void SPACoxMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double altFreq,
    int /*markerInChunkIdx*/,
    std::vector<double> &result
) {

    // altFreq = freq(bim col5 = ALT); GVec counts bim col5 alleles.
    //
    // ── Deliberate divergence from GRAB 0.2.4 (R) on altFreq > 0.5 ──
    //
    // GRAB 0.2.4 invokes Main.cpp::imputeGenoAndFlip which replaces
    //   GVec ← 2 − GVec        whenever altFreq > 0.5,
    // and then calls SPACoxClass::getMarkerPval with the original altFreq.
    // Inside that call the variance is built from
    //   adjGVec_R = GVec_flipped − 2·altFreq_original
    // whose sample mean is 2·(1 − altFreq) − 2·altFreq = 2 − 4·altFreq ≠ 0.
    // That is, the centring point used in the variance formula does not
    // equal the sample mean of the genotype vector actually entering the
    // score sum.  Letting Δ = 4·altFreq − 2,
    //   Σ adjGVec_R²  =  Σ adjGVec_correct²  +  N · Δ²
    // and the reported V̂_R = varResid · Σ adjGVec_R² over-estimates the
    // permutation-null variance by a factor 1 + N·Δ²/Σ adjGVec_correct².
    // Under H₀: β = 0 the R statistic
    //   Z_R = S_R / √V̂_R   ~  N(0, V_true / V̂_R)
    // is sub-Gaussian rather than N(0,1); the p-value 2·Φ(−|Z_R|) is
    // systematically conservative (too large) for altFreq > 0.5 markers,
    // and the score test loses its flip-invariance.  The defect depends
    // only on the genotype centring point, so it is present in both the
    // `time-to-event` and `Residual` paths of GRAB 0.2.4, regardless of
    // whether Σ Rᵢ = 0.
    //
    // We do not reproduce that defect here: GVec is consumed as-is, the
    // variance is built from (G − 2·altFreq)² whose centre 2·altFreq
    // equals the sample mean of G, and the resulting Z is exactly the
    // permutation-null score statistic — flip-invariant in the sense
    // that swapping REF/ALT in the input file negates S, leaves Var(S)
    // unchanged, and therefore leaves |Z| and the two-sided p-value
    // unchanged.  For markers with altFreq ≤ 0.5 the output matches
    // GRAB 0.2.4 to six-significant-digit print precision; for markers
    // with altFreq > 0.5 the output differs from GRAB 0.2.4 by exactly
    // the variance-inflation factor above, and represents the
    // statistically correct value rather than the R one.
    //
    // See dev-notes/methods/spacox.md §6 for the full derivation.
    const double S = GVec.dot(m_resid);
    double zScore, scoreVar;
    const spa::Result ts = getMarkerPvalCore(GVec, altFreq, S, zScore, scoreVar);
    pushResult(result, ts, S, scoreVar);
}

// ======================================================================
// MethodBase::getResultBatch — batched per-marker analysis
// ======================================================================
//
// The marker engine calls this with up to `preferredBatchSize()` markers
// per invocation.  The base-class default loops getResultVec() and copies
// each genotype column into a fresh Eigen::VectorXd; for SPACox that
// copy is N doubles per marker (~22000 × 2504 × 8 B ≈ 420 MB of memory
// traffic over a full 1KG chromosome) and the per-marker dot product is
// a BLAS-1 call.  Overriding here lets us:
//
//   1. fold B per-marker dot products into ONE BLAS-2 GEMV
//        scores = GBatchᵀ · m_resid                (B × 1, output)
//      Internally Eigen computes this column-by-column as the same
//      dot operation, so each entry scores[b] is bit-identical to the
//      original GBatch.col(b).dot(m_resid).
//
//   2. pass GBatch.col(b) into the SPA core as Eigen::Ref<const ...>
//      (no copy, no allocation) — the per-marker loop in
//      getMarkerPvalCore still walks gp[i] = GBatch(i, b) directly.
//
// Multi-phenotype fused-GEMM across the K phenotypes (i.e. one BLAS-3
// scoreMat = GBatchᵀ · residMat) is NOT done here because that would
// require the engine to route SPACox through its Phase-3 fused path,
// and Phase-3 hands the method only scalar score summaries — without
// access to the per-marker genotype vector required by the SPA branch
// (|Z| ≥ SPA_Cutoff).  See the comment on getResultVec above and
// dev-notes/methods/spacox.md §6 for the broader context.
void SPACoxMethod::getResultBatch(
    const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
    const std::vector<double> &altFreqs,
    const std::vector<int> & /*chunkIdxs*/,
    std::vector<std::vector<double> > &results
) {
    const int B = static_cast<int>(GBatch.cols());
    results.resize(B);

    // Batched score: one Eigen GEMV instead of B BLAS-1 dot products.
    const Eigen::VectorXd scores = GBatch.transpose() * m_resid;

    for (int b = 0; b < B; ++b) {
        auto &r = results[b];
        r.clear();
        r.reserve(6);
        double zScore, scoreVar;
        const spa::Result ts =
            getMarkerPvalCore(GBatch.col(b), altFreqs[b], scores[b], zScore, scoreVar);
        pushResult(r, ts, scores[b], scoreVar);
    }
}

// ======================================================================
// runSPACox — orchestration
// ======================================================================

void runSPACox(
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &covarNames,
    const std::string &phenoFile,
    const std::string &covarFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double pvalCovAdjCut,
    double spaCutoff,
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
    bool saveResid,
    bool longitudinal
) {
    // ---- Decide path: residual passthrough vs in-process null-model fit ----
    // --longitudinal pre-empts both: the per-subject residual R_G is computed
    // from a random-intercept LMM on the long-format pheno file and injected
    // after finalize (see below), so the GLM fit-path must not engage.
    const bool fitPath = !longitudinal && !phenoNameSpec.empty();
    nullmodel::RegressionModel regModel{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        regModel = nullmodel::parseRegressionModel(regressionModelStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
        infoMsg("SPACox: fitting %s null model for %zu phenotype(s)",
                nullmodel::regressionModelName(regModel), phenoSpecs.size());
    }

    // ---- Load resid/pheno file and covariate data ----
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);

    // --longitudinal: fit Y ~ X + (1|IID) on the long-format pheno file and
    // obtain the per-subject main-effect residual R_G before SubjectData is
    // built (parseLongPheno needs the .fam IID list).  SPACox has no GRM, so
    // the kept set is just (keep − remove).
    nsLongitudinal::LongResidResult Lr;
    if (longitudinal) {
        const auto outcomeNames = nsLongitudinal::splitOutcomeNames(phenoNameSpec);
        const auto kept = nsLongitudinal::buildKeptSet(keepFile, removeFile, famIIDs, {});
        infoMsg("SPACox: fitting longitudinal Y ~ X + (1|IID) for %zu outcome(s)",
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
        if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile);
    }
    if (!covarFile.empty() && !longitudinal) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();
    infoMsg("  %u subjects in union mask", sd.nUsed());

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

    // ---- If fit-path, run the null-model engine and inject residuals ----
    if (fitPath) {
        Eigen::MatrixXd covarUnion;
        if (!covarNames.empty()) {
            covarUnion = sd.getColumns(covarNames);
        } else if (sd.hasCovar()) {
            covarUnion = sd.covar();
        } else {
            covarUnion.resize(sd.nUsed(), 0);
        }
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
            // Rebuild fits vector with already-moved-from residuals for the writer:
            // we cannot reuse fits[*].residuals (moved into rs).  Reconstruct
            // NullModelFit views from rs/ns for writing.
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

    // ---- Build design-matrix (intercept + covariates) at union dimension ----
    infoMsg("Assembling union-level design matrix X = [intercept | covariates]");
    Eigen::MatrixXd unionX;
    if (longitudinal) {
        // The covariate-adjustment design matches the LMM fixed effects X
        // (intercept + --covar-name covariates absorbed into R_G); residualising
        // the genotype on the same X is the standard two-stage SPA procedure.
        unionX = nsLongitudinal::alignRowsToUsed(Lr.usedIIDs, Lr.perSubjectX, sd.usedIIDs());
    } else if (!covarNames.empty()) {
        auto cov = sd.getColumns(covarNames);
        unionX.resize(sd.nUsed(), cov.cols() + 1);
        unionX.col(0).setOnes();
        unionX.rightCols(cov.cols()) = cov;
    } else if (sd.hasCovar()) {
        const auto &cov = sd.covar();
        unionX.resize(sd.nUsed(), cov.cols() + 1);
        unionX.col(0).setOnes();
        unionX.rightCols(cov.cols()) = cov;
    } else {
        unionX.resize(sd.nUsed(), 1);
        unionX.col(0).setOnes();
    }

    // ---- Load genotype data (union mask) ----
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    // ---- Build per-phenotype tasks ----
    auto phenoInfos = sd.buildPerColumnMasks();
    const int K = sd.residOneCols();
    if (K > 1)
        infoMsg("Pre-computing per-phenotype design matrix X(X'X)^{-1} and SPA cumulant table (%d phenotypes, sequential):", K);
    else
        infoMsg("Pre-computing design matrix X(X'X)^{-1} and SPA cumulant table:");

    // Per-phenotype data (must outlive tasks — SPACoxMethod stores references)
    std::vector<Eigen::VectorXd> pResid(K);
    std::vector<CumulantTable> pCumul(K);
    std::vector<DesignMatrix> pDesign;
    pDesign.reserve(K);
    std::vector<size_t> designIdx(K);

    std::vector<PhenoTask> tasks(K);
    for (int rc = 0; rc < K; ++rc) {
        const auto &pi = phenoInfos[rc];
        const auto rcStart = std::chrono::steady_clock::now();

        // Extract per-phenotype residuals (always per-phenotype — residuals differ)
        pResid[rc] = (K > 1) ? extractPhenoVec(sd.residMatrix().col(rc), pi) : sd.residuals();

        // Cache DesignMatrix by non-missingness pattern (unionToLocal).
        // Phenotypes sharing the same valid-subject set produce identical
        // covariate matrices, so we deduplicate the expensive (X'X)^{-1}.
        size_t dIdx = pDesign.size(); // default: build new
        std::string reusedFrom;
        if (K > 1) {
            for (int j = 0; j < rc; ++j) {
                if (phenoInfos[j].unionToLocal == pi.unionToLocal) {
                    dIdx = designIdx[j];
                    reusedFrom = phenoInfos[j].name;
                    break;
                }
            }
        }
        if (dIdx == pDesign.size()) {
            Eigen::MatrixXd phenoX = (K > 1) ? extractPhenoMat(unionX, pi) : unionX;
            pDesign.emplace_back(phenoX);
        }
        designIdx[rc] = dIdx;

        // Build cumulant table from per-phenotype residuals (not cacheable).
        // The 10 000 grid points are independent, so the build itself is
        // threaded (P5); the per-phenotype loop stays sequential because each
        // table is already wide enough to saturate the pool.
        pCumul[rc] = spacox_cgf::buildCumulantTable(pResid[rc], nthread);
        double meanR = pResid[rc].mean();
        double varR = (pResid[rc].array() - meanR).square().mean();
        double N = static_cast<double>(pResid[rc].size());
        double varResid = varR * N / (N - 1.0);

        tasks[rc].phenoName = pi.name;
        tasks[rc].method = std::make_unique<SPACoxMethod>(
            pResid[rc], varResid, pCumul[rc], pDesign[designIdx[rc]], pvalCovAdjCut, spaCutoff);
        tasks[rc].unionToLocal = pi.unionToLocal;
        tasks[rc].nUsed = pi.nUsed;

        const double rcSec = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - rcStart).count();
        if (!reusedFrom.empty())
            infoMsg("  '%s': cumulant table built (%.1fs); reuses design matrix from '%s'",
                    pi.name.c_str(), rcSec, reusedFrom.c_str());
        else
            infoMsg("  '%s': design matrix + cumulant table built (%.1fs)",
                    pi.name.c_str(), rcSec);
    }
    if (K > 1)
        infoMsg("  %zu unique design matrix(es) for %d phenotypes", pDesign.size(), K);

    infoMsg("Running SPACox marker tests (%d thread(s), %d phenotype(s))...", nthread, K);
    multiPhenoEngine(
        *genoData,
        tasks,
        outPrefix,
        "SPACox",
        compression,
        compressionLevel,
        nthread,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        hweCutoff
    );
}

// ======================================================================
// runSPACoxLoco — LOCO orchestration
// ======================================================================
//
// For each chromosome the null model is refit with that chromosome's LOCO PGS
// appended as one estimated covariate column, and both the residuals AND the
// stage-2 covariate-projection DesignMatrix are rebuilt so the PGS is projected
// out of the genotype consistently.  Because each phenotype has its own PGS
// column, the refit is per-phenotype (fitAll with a single spec).

void runSPACoxLoco(
    const std::vector<std::string> &covarNames,
    const std::string &phenoFile,
    const std::string &covarFile,
    const GenoSpec &geno,
    const std::string &predListFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double pvalCovAdjCut,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &regressionModelStr,
    const std::string &phenoNameSpec
) {
    if (phenoNameSpec.empty())
        throw std::runtime_error(
            "SPACox-LOCO requires --pheno-name (an in-process null-model fit); "
            "precomputed residuals cannot be refit per chromosome");

    nullmodel::RegressionModel regModel =
        nullmodel::parseRegressionModel(regressionModelStr);
    std::vector<nullmodel::PhenoSpec> phenoSpecs =
        nullmodel::parsePhenoSpecList(regModel, phenoNameSpec);
    const int K = static_cast<int>(phenoSpecs.size());
    std::vector<std::string> specNames(K);
    for (int k = 0; k < K; ++k) specNames[k] = phenoSpecs[k].name;

    infoMsg("SPACox-LOCO: fitting %s null model for %d phenotype(s), "
            "per-chromosome LOCO PGS as covariate",
            nullmodel::regressionModelName(regModel), K);

    // Early validation: every phenotype must have an entry in the pred.list.
    validatePredListPhenos(predListFile, specNames);

    // ---- Load pheno/covar data ----
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    sd.loadPhenoFile(phenoFile, nullmodel::columnsNeeded(phenoSpecs));
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();
    infoMsg("  %u subjects in union mask", sd.nUsed());

    // ---- Base covariate matrix (no PGS), used as the augmentation base ----
    Eigen::MatrixXd baseCovar;
    if (!covarNames.empty()) baseCovar = sd.getColumns(covarNames);
    else if (sd.hasCovar()) baseCovar = sd.covar();
    else baseCovar.resize(sd.nUsed(), 0);

    // ---- One base fit to establish per-phenotype masks (residual values are
    //      discarded — they are recomputed per chromosome). ----
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

    // Resolve each spec to a concrete model once (inference is from the
    // chromosome-invariant phenotype column), so the per-chromosome refits do
    // not re-emit an "Auto-detected" log line for every chromosome.
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

    // ---- Union-level design base [intercept | covariates] (PGS appended per chr) ----
    const Eigen::Index nUnion = sd.nUsed();
    Eigen::MatrixXd unionXBase(nUnion, baseCovar.cols() + 1);
    unionXBase.col(0).setOnes();
    if (baseCovar.cols() > 0) unionXBase.rightCols(baseCovar.cols()) = baseCovar;

    // ---- Per-phenotype non-missing masks (union space, chromosome-invariant) ----
    std::vector<std::vector<bool> > nonMissing(K, std::vector<bool>(nUnion, false));
    for (int k = 0; k < K; ++k)
        for (Eigen::Index i = 0; i < nUnion; ++i)
            nonMissing[k][static_cast<size_t>(i)] =
                (phenoInfos[k].unionToLocal[static_cast<size_t>(i)] != UINT32_MAX);

    // ---- Load genotype data + LOCO predictions ----
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);
    LocoData loco = LocoData::load(predListFile, specNames, sd.usedIIDs(), sd.famIIDs());
    auto locoChroms = loco.availableChromosomes();
    infoMsg("LOCO: %zu chromosomes available across all phenotypes", locoChroms.size());

    // ---- Persistent per-phenotype storage (must outlive tasks — SPACoxMethod
    //      stores const references into these).  reserve(K) keeps addresses
    //      stable across clear()/emplace_back cycles. ----
    std::vector<Eigen::VectorXd> pResid(K);
    std::vector<CumulantTable> pCumul(K);
    std::vector<DesignMatrix> pDesign;
    pDesign.reserve(K);

    nullmodel::EngineOptions eo1;
    eo1.nthreads = 1;

    auto buildTasks = [&](const std::string &chr, std::vector<PhenoTask> &tasks) {
        tasks.resize(K);
        // Destroy the previous chromosome's methods BEFORE overwriting the pools
        // they reference (lifetime discipline).
        for (auto &t : tasks) t.method.reset();
        pDesign.clear();

        for (int k = 0; k < K; ++k) {
            Eigen::MatrixXd covarUnion_k = appendLocoCovariate(
                loco, specNames[k], chr, baseCovar, nonMissing[k], "SPACox-LOCO");

            auto fits1 = nullmodel::fitAll(
                sd, {phenoSpecs[k]}, specModel[k], covarUnion_k, eo1);
            pResid[k] = extractPhenoVec(fits1[0].residuals, phenoInfos[k]);

            // Augmented design [intercept | covariates | PGS_chr]; the PGS is the
            // last column of covarUnion_k.
            Eigen::MatrixXd unionX_k(nUnion, unionXBase.cols() + 1);
            unionX_k.leftCols(unionXBase.cols()) = unionXBase;
            unionX_k.col(unionXBase.cols()) =
                covarUnion_k.col(covarUnion_k.cols() - 1);
            Eigen::MatrixXd phenoX = extractPhenoMat(unionX_k, phenoInfos[k]);
            pDesign.emplace_back(phenoX);

            pCumul[k] = spacox_cgf::buildCumulantTable(pResid[k], nthread);
            double meanR = pResid[k].mean();
            double varR = (pResid[k].array() - meanR).square().mean();
            double N = static_cast<double>(pResid[k].size());
            double varResid = varR * N / (N - 1.0);

            tasks[k].phenoName = phenoInfos[k].name;
            tasks[k].method = std::make_unique<SPACoxMethod>(
                pResid[k], varResid, pCumul[k], pDesign[k], pvalCovAdjCut, spaCutoff);
            tasks[k].unionToLocal = phenoInfos[k].unionToLocal;
            tasks[k].nUsed = phenoInfos[k].nUsed;
        }
    };

    infoMsg("SPACox-LOCO: starting LOCO association (%d phenotype(s), %zu chroms, %d threads)",
            K, locoChroms.size(), nthread);
    locoEngine(
        *genoData,
        locoChroms,
        specNames,
        buildTasks,
        outPrefix,
        "SPACox",
        compression,
        compressionLevel,
        nthread,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        hweCutoff
    );
}
