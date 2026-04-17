// spacox.cpp — SPACox full implementation

#include "spacox/spacox.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
      m_spaCutoff(spaCutoff),
      m_adjGNorm(residuals.size()),
      m_adjGVec(residuals.size())
{
    m_nzSet.reserve(m_N);
}

std::unique_ptr<MethodBase> SPACoxMethod::clone() const {
    // Shared const refs stay the same; scratch buffers are freshly allocated.
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

double SPACoxMethod::getMarkerPval(
    const Eigen::Ref<Eigen::VectorXd> &GVec,
    double MAF,
    double &zScore
) {

    // Score statistic S = G' * resid
    double S = GVec.dot(m_resid);
    double twoMAF = 2.0 * MAF;
    const double *gp = GVec.data();

    // VarS = varResid * Σ(g_i - 2·MAF)²
    double sumAdjG2 = 0.0;
    for (int i = 0; i < m_N; ++i) {
        double adj = gp[i] - twoMAF;
        sumAdjG2 += adj * adj;
    }
    double VarS = m_varResid * sumAdjG2;
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

    double *anp = m_adjGNorm.data();
    m_nzSet.clear();
    for (int i = 0; i < m_N; ++i) {
        anp[i] = (gp[i] - twoMAF) / sqrtVarS;
        if (gp[i] != 0.0) m_nzSet.push_back(static_cast<uint32_t>(i));
    }
    int nNz = static_cast<int>(m_nzSet.size());
    int N0 = m_N - nNz;

    double absZ = std::abs(zScore);
    double pval = getProbSpa(adjG0, anp, m_nzSet.data(), nNz, N0, absZ, false) +
                  getProbSpa(adjG0, anp, m_nzSet.data(), nNz, N0, -absZ, true);

    if (pval > m_pvalCovAdjCut) return pval;

    // ---- Stage 2: covariate-adjusted SPA ----
    m_design.adjustGenotype(gp, m_nzSet.data(), nNz, m_adjGVec);

    VarS = m_varResid * m_adjGVec.squaredNorm();
    if (VarS <= 0.0) {
        zScore = 0.0;
        return std::numeric_limits<double>::quiet_NaN();
    }
    zScore = S / std::sqrt(VarS);
    sqrtVarS = std::sqrt(VarS);

    const double *avp = m_adjGVec.data();
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
    double zScore;
    double pval = getMarkerPval(GVec, altFreq, zScore);

    result.push_back(pval);
    result.push_back(zScore);
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
    const std::string &removeFile
) {

    // ---- Load resid file and covariate data ----
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    sd.loadResidOne(phenoFile, residNames);
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();
    infoMsg("  %u subjects in union mask", sd.nUsed());

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
