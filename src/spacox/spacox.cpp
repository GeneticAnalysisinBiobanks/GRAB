// spacox.cpp — SPACox full implementation

#include "spacox/spacox.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/null_model.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ======================================================================
// Worker-local SPA scratch
// ======================================================================
//
// The SPA branch of getMarkerPvalCore needs three transient buffers
// (adjGNorm of length N, adjGVec of length N, nzSet of capacity ≤ N).
// These were previously per-clone fields on SPACoxMethod, which meant
// the multi-phenotype engine paid (K phenotypes × T worker threads)
// times their footprint.  Promoting them to a translation-unit-local
// thread_local struct lets all K phenotype-clones in one worker share
// one set of buffers — the engine spawns at most T concurrent threads,
// so total scratch footprint is bounded by T × (≈ 20N bytes) regardless
// of K.  Resize() is idempotent when the requested length matches the
// existing capacity, so phenotype switching inside a thread is O(1)
// after the first call.
namespace {
struct SPACoxScratch {
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
// Build empirical CGF interpolation table
// ======================================================================

CumulantTable buildCumulantTable(const Eigen::VectorXd &residuals) {
    constexpr double rangeMax = 100.0;
    constexpr int L = 10000;

    // Cauchy-quantile spacing: idx0 = tan(pi*(k/(L+1) - 0.5))
    // then rescale so max(idx1) = rangeMax.
    Eigen::VectorXd xGrid(L);
    double gridScale;
    {
        double maxAbs = 0.0;
        for (int k = 0; k < L; ++k) {
            double p = static_cast<double>(k + 1) / static_cast<double>(L + 1);
            xGrid[k] = std::tan(M_PI * (p - 0.5));
            double a = std::fabs(xGrid[k]);
            if (a > maxAbs) maxAbs = a;
        }
        gridScale = rangeMax / maxAbs;
        xGrid *= gridScale;
    }

    const int N = static_cast<int>(residuals.size());
    const double *rp = residuals.data();

    Eigen::VectorXd yK0(L), yK1(L), yK2(L);

    for (int i = 0; i < L; ++i) {
        double t = xGrid[i];
        double M0 = 0.0, M1 = 0.0, M2 = 0.0;
        for (int j = 0; j < N; ++j) {
            double r = rp[j];
            double e = std::exp(r * t);
            M0 += e;
            double re = r * e;
            M1 += re;
            M2 += r * re;
        }
        double invN = 1.0 / static_cast<double>(N);
        M0 *= invN;
        M1 *= invN;
        M2 *= invN;
        yK0[i] = std::log(M0);
        yK1[i] = M1 / M0;
        yK2[i] = (M0 * M2 - M1 * M1) / (M0 * M0);
    }

    // Pre-compute slopes for piecewise-linear interpolation.
    auto computeSlopes = [](const Eigen::VectorXd &x, const Eigen::VectorXd &y) -> Eigen::VectorXd {
        const int n = static_cast<int>(x.size());
        Eigen::VectorXd s(n - 1);
        for (int i = 0; i < n - 1; ++i)
            s[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
        return s;
    };

    CumulantTable ct;
    ct.xGrid = std::move(xGrid);
    ct.nGrid = L;
    ct.invScale = 1.0 / gridScale;
    ct.Lp1 = static_cast<double>(L + 1);
    ct.yK0 = std::move(yK0);
    ct.slopeK0 = computeSlopes(ct.xGrid, ct.yK0);
    ct.yK1 = std::move(yK1);
    ct.slopeK1 = computeSlopes(ct.xGrid, ct.yK1);
    ct.yK2 = std::move(yK2);
    ct.slopeK2 = computeSlopes(ct.xGrid, ct.yK2);
    return ct;
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
    //   * the SPA scratch (adjGNorm / adjGVec / nzSet) is thread_local in
    //     this translation unit, so K phenotype-clones in one worker share
    //     a single set of buffers.
    // The clone only owns its small scalar metadata (m_varResid, m_N,
    // m_pvalCovAdjCut, m_spaCutoff) and the references themselves — the
    // entire clone is therefore on the order of one cache line.
    return std::make_unique<SPACoxMethod>(m_resid, m_varResid, m_cumul, m_design, m_pvalCovAdjCut, m_spaCutoff);
}

// ======================================================================
// CGF interpolation — O(1) Cauchy-inverse index lookup
// ======================================================================

int SPACoxMethod::interpIdx(double v) const {
    // Inverse of Cauchy-quantile grid: idx = (atan(v/scale)/π + 0.5)*(L+1) - 1
    double fIdx = (std::atan(v * m_cumul.invScale) * (1.0 / M_PI) + 0.5) * m_cumul.Lp1 - 1.0;
    int lo = static_cast<int>(fIdx);
    if (lo < 0) lo = 0;
    if (lo >= m_cumul.nGrid - 1) lo = m_cumul.nGrid - 2;
    return lo;
}

double SPACoxMethod::interp(
    const double *yp,
    const double *sp,
    int lo,
    double v
) const {
    return yp[lo] + (v - m_cumul.xGrid.data()[lo]) * sp[lo];
}

double SPACoxMethod::interpK0(double v) const {
    const int n = m_cumul.nGrid;
    if (v <= m_cumul.xGrid.data()[0]) return m_cumul.yK0.data()[0];
    if (v >= m_cumul.xGrid.data()[n - 1]) return m_cumul.yK0.data()[n - 1];
    int lo = interpIdx(v);
    return interp(m_cumul.yK0.data(), m_cumul.slopeK0.data(), lo, v);
}

double SPACoxMethod::interpK1(double v) const {
    const int n = m_cumul.nGrid;
    if (v <= m_cumul.xGrid.data()[0]) return m_cumul.yK1.data()[0];
    if (v >= m_cumul.xGrid.data()[n - 1]) return m_cumul.yK1.data()[n - 1];
    int lo = interpIdx(v);
    return interp(m_cumul.yK1.data(), m_cumul.slopeK1.data(), lo, v);
}

double SPACoxMethod::interpK2(double v) const {
    const int n = m_cumul.nGrid;
    if (v <= m_cumul.xGrid.data()[0]) return m_cumul.yK2.data()[0];
    if (v >= m_cumul.xGrid.data()[n - 1]) return m_cumul.yK2.data()[n - 1];
    int lo = interpIdx(v);
    return interp(m_cumul.yK2.data(), m_cumul.slopeK2.data(), lo, v);
}

// ======================================================================
// Cumulant evaluation — fused K1+K2 for Newton-Raphson
// ======================================================================

double SPACoxMethod::evalK0(
    double t,
    int N0,
    double adjG0,
    const double *adjG,
    const uint32_t *idx,
    int n
) const {
    double sum = N0 * interpK0(t * adjG0);
    if (idx) {
        for (int k = 0; k < n; ++k)
            sum += interpK0(t * adjG[idx[k]]);
    } else {
        for (int k = 0; k < n; ++k)
            sum += interpK0(t * adjG[k]);
    }
    return sum;
}

std::pair<double, double> SPACoxMethod::evalK1K2(
    double t,
    int N0,
    double adjG0,
    const double *adjG,
    const uint32_t *idx,
    int n,
    double q2
) const {

    const double *xp = m_cumul.xGrid.data();
    const double *y1 = m_cumul.yK1.data();
    const double *s1 = m_cumul.slopeK1.data();
    const double *y2 = m_cumul.yK2.data();
    const double *s2 = m_cumul.slopeK2.data();
    const int nG = m_cumul.nGrid;
    const double invS = m_cumul.invScale;
    const double Lp1 = m_cumul.Lp1;
    const double invPi = 1.0 / M_PI;

    // Inline interpolation: compute bin once, read K1 and K2 from same bin
    auto interpBoth = [&](double v) -> std::pair<double, double> {
        if (v <= xp[0]) return {y1[0], y2[0]};
        if (v >= xp[nG - 1]) return {y1[nG - 1], y2[nG - 1]};
        double fIdx = (std::atan(v * invS) * invPi + 0.5) * Lp1 - 1.0;
        int lo = static_cast<int>(fIdx);
        if (lo < 0) lo = 0;
        if (lo >= nG - 1) lo = nG - 2;
        double dx = v - xp[lo];
        return {y1[lo] + dx * s1[lo], y2[lo] + dx * s2[lo]};
    };

    double sumK1 = 0.0, sumK2 = 0.0;

    // Contribution from N0 zero-genotype subjects
    {
        double v0 = t * adjG0;
        auto [k1, k2] = interpBoth(v0);
        sumK1 = N0 * adjG0 * k1;
        sumK2 = N0 * adjG0 * adjG0 * k2;
    }

    if (idx) {
        for (int k = 0; k < n; ++k) {
            double a = adjG[idx[k]];
            auto [k1, k2] = interpBoth(t * a);
            sumK1 += a * k1;
            sumK2 += a * a * k2;
        }
    } else {
        for (int k = 0; k < n; ++k) {
            double a = adjG[k];
            auto [k1, k2] = interpBoth(t * a);
            sumK1 += a * k1;
            sumK2 += a * a * k2;
        }
    }

    return {sumK1 - q2, sumK2};
}

// ======================================================================
// Newton-Raphson root-finding on K1(t) - q2 = 0
// ======================================================================

SPACoxMethod::RootResult SPACoxMethod::fastGetRootK1(
    double initX,
    int N0,
    double adjG0,
    const double *adjG,
    const uint32_t *idx,
    int n,
    double q2
) const {

    double x = initX, oldX;
    double K1val = 0.0, K2val = 0.0, oldK1;
    double diffX = std::numeric_limits<double>::infinity(), oldDiffX;
    bool converge = true;
    constexpr double tol = 0.001;
    constexpr int maxiter = 100;

    for (int iter = 0; iter < maxiter; ++iter) {
        oldX = x;
        oldDiffX = diffX;
        oldK1 = K1val;

        auto [k1, k2] = evalK1K2(x, N0, adjG0, adjG, idx, n, q2);
        K1val = k1;
        K2val = k2;

        diffX = -K1val / K2val;

        if (!std::isfinite(K1val)) {
            x = std::numeric_limits<double>::infinity();
            K2val = 0.0;
            break;
        }

        if ((K1val > 0) != (oldK1 > 0)) {
            while (std::abs(diffX) > std::abs(oldDiffX) - tol)
                diffX *= 0.5;
        }

        if (std::abs(diffX) < tol) break;
        x = oldX + diffX;

        if (iter == maxiter - 1) converge = false;
    }
    return {x, converge, K2val};
}

// ======================================================================
// SPA tail probability (Lugannani-Rice)
// ======================================================================

double SPACoxMethod::getProbSpa(
    double adjG0,
    const double *adjG,
    const uint32_t *idx,
    int n,
    int N0,
    double q2,
    bool lowerTail
) const {

    double initX = (q2 > 0) ? 3.0 : -3.0;

    RootResult rr = fastGetRootK1(initX, N0, adjG0, adjG, idx, n, q2);
    double zeta = rr.root;

    double k0val = evalK0(zeta, N0, adjG0, adjG, idx, n);
    double k2val = rr.K2; // already computed at converged root
    double temp1 = zeta * q2 - k0val;

    if (!std::isfinite(zeta) || temp1 < 0.0 || k2val <= 0.0) return std::numeric_limits<double>::quiet_NaN();

    double w = std::copysign(std::sqrt(2.0 * temp1), zeta);
    double v = zeta * std::sqrt(k2val);

    if (w == 0.0 || v == 0.0 || (v / w) <= 0.0) return std::numeric_limits<double>::quiet_NaN();

    double sign = lowerTail ? 1.0 : -1.0;
    return math::pnorm(sign * (w + (1.0 / w) * std::log(v / w)));
}

// ======================================================================
// Per-marker p-value (two-stage SPA)
// ======================================================================

double SPACoxMethod::getMarkerPvalCore(
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    double S,
    double &zScore,
    double &outScoreVar
) {

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
    double sumAdjG2 = 0.0;
    for (int i = 0; i < m_N; ++i) {
        double adj = gp[i] - twoMAF;
        sumAdjG2 += adj * adj;
    }
    double VarS = m_varResid * sumAdjG2;
    outScoreVar = VarS;
    if (VarS <= 0.0) {
        zScore = 0.0;
        return std::numeric_limits<double>::quiet_NaN();
    }
    zScore = S / std::sqrt(VarS);

    // Normal approximation if |z| below SPA cutoff
    if (std::abs(zScore) < m_spaCutoff) return 2.0 * math::pnorm(-std::abs(zScore));

    // ---- Stage 1: unadjusted SPA with non-zero index optimisation ----
    double sqrtVarS = std::sqrt(VarS);
    double adjG0 = -twoMAF / sqrtVarS;

    // Acquire (and resize if necessary) the thread-local SPA scratch.
    // resize() is a no-op when the requested length already matches.
    tlScratch.adjGNorm.resize(m_N);
    tlScratch.nzSet.clear();
    if (static_cast<int>(tlScratch.nzSet.capacity()) < m_N)
        tlScratch.nzSet.reserve(m_N);

    double *anp = tlScratch.adjGNorm.data();
    for (int i = 0; i < m_N; ++i) {
        anp[i] = (gp[i] - twoMAF) / sqrtVarS;
        if (gp[i] != 0.0) tlScratch.nzSet.push_back(static_cast<uint32_t>(i));
    }
    int nNz = static_cast<int>(tlScratch.nzSet.size());
    int N0 = m_N - nNz;

    double absZ = std::abs(zScore);
    double pval = getProbSpa(adjG0, anp, tlScratch.nzSet.data(), nNz, N0, absZ, false) +
                  getProbSpa(adjG0, anp, tlScratch.nzSet.data(), nNz, N0, -absZ, true);

    if (pval > m_pvalCovAdjCut) return pval;

    // ---- Stage 2: covariate-adjusted SPA ----
    tlScratch.adjGVec.resize(m_N);
    m_design.adjustGenotype(gp, tlScratch.nzSet.data(), nNz, tlScratch.adjGVec);

    VarS = m_varResid * tlScratch.adjGVec.squaredNorm();
    outScoreVar = VarS;
    if (VarS <= 0.0) {
        zScore = 0.0;
        return std::numeric_limits<double>::quiet_NaN();
    }
    zScore = S / std::sqrt(VarS);
    sqrtVarS = std::sqrt(VarS);

    const double *avp = tlScratch.adjGVec.data();
    for (int i = 0; i < m_N; ++i)
        anp[i] = avp[i] / sqrtVarS;

    // Full-vector SPA (no sparsity shortcut after adjustment)
    absZ = std::abs(zScore);
    pval = getProbSpa(0.0, anp, nullptr, m_N, 0, absZ, false) + getProbSpa(0.0, anp, nullptr, m_N, 0, -absZ, true);

    return pval;
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
    // See docs/methods/spacox.md §6 for the full derivation.
    const double S = GVec.dot(m_resid);
    double zScore, scoreVar;
    double pval = getMarkerPvalCore(GVec, altFreq, S, zScore, scoreVar);

    result.push_back(pval);
    if (scoreVar > 0.0) {
        result.push_back(S / scoreVar);
        result.push_back(1.0 / std::sqrt(scoreVar));
    } else {
        result.push_back(std::numeric_limits<double>::quiet_NaN());
        result.push_back(std::numeric_limits<double>::quiet_NaN());
    }
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
// docs/methods/spacox.md §6 for the broader context.
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
        r.reserve(3);
        double zScore, scoreVar;
        double pval =
            getMarkerPvalCore(GBatch.col(b), altFreqs[b], scores[b], zScore, scoreVar);
        r.push_back(pval);
        if (scoreVar > 0.0) {
            r.push_back(scores[b] / scoreVar);
            r.push_back(1.0 / std::sqrt(scoreVar));
        } else {
            r.push_back(std::numeric_limits<double>::quiet_NaN());
            r.push_back(std::numeric_limits<double>::quiet_NaN());
        }
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
    const std::string &traitTypeStr,
    const std::string &phenoNameSpec,
    bool saveResid
) {
    // ---- Decide path: residual passthrough vs in-process null-model fit ----
    const bool fitPath = !phenoNameSpec.empty();
    nullmodel::TraitType traitT{};
    std::vector<nullmodel::PhenoSpec> phenoSpecs;
    if (fitPath) {
        traitT = nullmodel::parseTraitType(traitTypeStr);
        phenoSpecs = nullmodel::parsePhenoSpecList(traitT, phenoNameSpec);
        infoMsg("SPACox: fitting %s null model for %zu phenotype(s)",
                nullmodel::traitTypeName(traitT), phenoSpecs.size());
    }

    // ---- Load resid/pheno file and covariate data ----
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
    infoMsg("  %u subjects in union mask", sd.nUsed());

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
            // Rebuild fits vector with already-moved-from residuals for the writer:
            // we cannot reuse fits[*].residuals (moved into rs).  Reconstruct
            // NullModelFit views from rs/ns for writing.
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

    // ---- Build design-matrix (intercept + covariates) at union dimension ----
    infoMsg("Building design matrix projection...");
    Eigen::MatrixXd unionX;
    if (!covarNames.empty()) {
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
    if (K > 1) infoMsg("Multi-column residual file: %d phenotypes", K);

    // Per-phenotype data (must outlive tasks — SPACoxMethod stores references)
    std::vector<Eigen::VectorXd> pResid(K);
    std::vector<CumulantTable> pCumul(K);
    std::vector<DesignMatrix> pDesign;
    pDesign.reserve(K);
    std::vector<size_t> designIdx(K);

    std::vector<PhenoTask> tasks(K);
    for (int rc = 0; rc < K; ++rc) {
        const auto &pi = phenoInfos[rc];

        // Extract per-phenotype residuals (always per-phenotype — residuals differ)
        pResid[rc] = (K > 1) ? extractPhenoVec(sd.residMatrix().col(rc), pi) : sd.residuals();

        // Cache DesignMatrix by non-missingness pattern (unionToLocal).
        // Phenotypes sharing the same valid-subject set produce identical
        // covariate matrices, so we deduplicate the expensive (X'X)^{-1}.
        size_t dIdx = pDesign.size(); // default: build new
        if (K > 1) {
            for (int j = 0; j < rc; ++j) {
                if (phenoInfos[j].unionToLocal == pi.unionToLocal) {
                    dIdx = designIdx[j];
                    infoMsg("  Phenotype '%s': reusing design matrix from '%s'",
                            pi.name.c_str(), phenoInfos[j].name.c_str());
                    break;
                }
            }
        }
        if (dIdx == pDesign.size()) {
            Eigen::MatrixXd phenoX = (K > 1) ? extractPhenoMat(unionX, pi) : unionX;
            pDesign.emplace_back(phenoX);
        }
        designIdx[rc] = dIdx;

        // Build cumulant table from per-phenotype residuals (not cacheable)
        pCumul[rc] = buildCumulantTable(pResid[rc]);
        double meanR = pResid[rc].mean();
        double varR = (pResid[rc].array() - meanR).square().mean();
        double N = static_cast<double>(pResid[rc].size());
        double varResid = varR * N / (N - 1.0);

        tasks[rc].phenoName = pi.name;
        tasks[rc].method = std::make_unique<SPACoxMethod>(
            pResid[rc], varResid, pCumul[rc], pDesign[designIdx[rc]], pvalCovAdjCut, spaCutoff);
        tasks[rc].unionToLocal = pi.unionToLocal;
        tasks[rc].nUsed = pi.nUsed;
        infoMsg("  Phenotype '%s': %u subjects", pi.name.c_str(), pi.nUsed);
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
