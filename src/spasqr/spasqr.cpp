// spasqr.cpp — SPAsqr: SPA-squared multi-tau marker association (pure C++17 / Eigen / Boost)

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "spasqr/spasqr.hpp"
#include "engine/loco.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "spasqr/conquer.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <Eigen/Dense>
#include <immintrin.h>
#include "util/simd_dispatch.hpp"
#include "util/simd_math.hpp"

// ══════════════════════════════════════════════════════════════════════
// SPAsqr-specific scalar SPA functions
//
// SPAsqr always uses empty FamilyData (no TwoSubj/ThreeSubj), so the
// MGF computation reduces to a loop over unrelated outlier residuals
// only.  We replace the vector-based MgfWorkspace (8 × nOutlier doubles
// per clone) with scalar accumulators (O(1) extra memory per clone),
// eliminating the dominant per-clone memory cost.
// ══════════════════════════════════════════════════════════════════════

namespace {

// ── Scalar CGF result from outlier loop ────────────────────────────
struct ScalarCgfResult {
    double logMgf0Sum;    // Σ log(MGF0_i)  — for CGF0
    double cgf1Outlier;   // Σ (MGF1_i / MGF0_i)  — for CGF1
    double cgf2Outlier;   // Σ (MGF2_i / MGF0_i) - Σ (MGF1_i/MGF0_i)^2  — for CGF2
};

// Single-pass scalar computation of CGF contributions from unrelated
// outliers.  Equivalent to nsSPAGRM::mgf() for the outlier-only case,
// but uses O(1) memory instead of O(nOutlier).
//
// Per-outlier MGF for diploid genotype model (same as spagrm.cpp):
//   alpha_i = (1-MAF) + MAF * exp(t * r_i)
//   MGF0_i  = alpha_i^2
//   MGF1_i  = 2 * alpha_i * alpha1_i     where alpha1_i = MAF * r_i * exp(t*r_i)
//   MGF2_i  = 2 * (alpha1_i^2 + alpha_i * alpha2_i)  where alpha2_i = r_i * alpha1_i
inline ScalarCgfResult scalarOutlierCgf(
    double t,
    const double *outlierResid,
    int nOutlier,
    double MAF
) {
    double logMgf0Sum = 0.0;
    double tempSum    = 0.0;   // Σ MGF1_i / MGF0_i
    double s1         = 0.0;   // Σ MGF2_i / MGF0_i
    double s2         = 0.0;   // Σ (MGF1_i / MGF0_i)^2
    const double oneMmaf = 1.0 - MAF;

    for (int i = 0; i < nOutlier; ++i) {
        const double r      = outlierResid[i];
        const double lambda = std::exp(t * r);
        const double alpha  = oneMmaf + MAF * lambda;
        const double alpha1 = MAF * r * lambda;
        const double alpha2 = r * alpha1;

        // MGF0 = alpha^2,  MGF1 = 2*alpha*alpha1,  MGF2 = 2*(alpha1^2 + alpha*alpha2)
        const double mgf0   = alpha * alpha;
        const double mgf1   = 2.0 * alpha * alpha1;
        const double mgf2   = 2.0 * (alpha1 * alpha1 + alpha * alpha2);

        logMgf0Sum += std::log(mgf0);
        const double ratio  = mgf1 / mgf0;   // MGF1/MGF0
        tempSum += ratio;
        s1      += mgf2 / mgf0;              // MGF2/MGF0
        s2      += ratio * ratio;             // (MGF1/MGF0)^2
    }

    return {logMgf0Sum, tempSum, s1 - s2};
}

// ── SIMD-vectorized outlier CGF (C1): processes outliers in SIMD lanes ──

__attribute__((target("avx2,avx512f,avx512vl,fma")))
static ScalarCgfResult simdOutlierCgf_avx512(
    double t,
    const double *outlierResid,
    int nOutlier,
    double MAF
) {
    const __m512d vt       = _mm512_set1_pd(t);
    const __m512d vMAF     = _mm512_set1_pd(MAF);
    const __m512d vOneMmaf = _mm512_set1_pd(1.0 - MAF);
    const __m512d vTwo     = _mm512_set1_pd(2.0);

    __m512d vLogMgf0 = _mm512_setzero_pd();
    __m512d vTempSum = _mm512_setzero_pd();
    __m512d vS1      = _mm512_setzero_pd();
    __m512d vS2      = _mm512_setzero_pd();

    int i = 0;
    const int n8 = nOutlier & ~7;
    for (; i < n8; i += 8) {
        const __m512d vr      = _mm512_loadu_pd(outlierResid + i);
        const __m512d vlambda = avx512_exp_pd(_mm512_mul_pd(vt, vr));
        const __m512d valpha  = _mm512_fmadd_pd(vMAF, vlambda, vOneMmaf);
        const __m512d valpha1 = _mm512_mul_pd(vMAF, _mm512_mul_pd(vr, vlambda));
        const __m512d valpha2 = _mm512_mul_pd(vr, valpha1);

        const __m512d vmgf0 = _mm512_mul_pd(valpha, valpha);
        const __m512d vmgf1 = _mm512_mul_pd(vTwo, _mm512_mul_pd(valpha, valpha1));
        const __m512d vmgf2 = _mm512_mul_pd(vTwo,
                                            _mm512_fmadd_pd(valpha1, valpha1, _mm512_mul_pd(valpha, valpha2)));

        vLogMgf0 = _mm512_add_pd(vLogMgf0, avx512_log_pd(vmgf0));
        const __m512d vratio = _mm512_div_pd(vmgf1, vmgf0);
        vTempSum = _mm512_add_pd(vTempSum, vratio);
        vS1      = _mm512_add_pd(vS1, _mm512_div_pd(vmgf2, vmgf0));
        vS2      = _mm512_fmadd_pd(vratio, vratio, vS2);
    }

    // Masked tail: r=0 for inactive lanes → lambda=1, alpha=1, mgf0=1,
    // log(1)=0, mgf1=0, ratio=0 — all zero contributions.
    if (i < nOutlier) {
        const __mmask8 mask = static_cast<__mmask8>((1u << (nOutlier - i)) - 1u);
        const __m512d vr      = _mm512_maskz_loadu_pd(mask, outlierResid + i);
        const __m512d vlambda = avx512_exp_pd(_mm512_mul_pd(vt, vr));
        const __m512d valpha  = _mm512_fmadd_pd(vMAF, vlambda, vOneMmaf);
        const __m512d valpha1 = _mm512_mul_pd(vMAF, _mm512_mul_pd(vr, vlambda));
        const __m512d valpha2 = _mm512_mul_pd(vr, valpha1);

        const __m512d vmgf0 = _mm512_mul_pd(valpha, valpha);
        const __m512d vmgf1 = _mm512_mul_pd(vTwo, _mm512_mul_pd(valpha, valpha1));
        const __m512d vmgf2 = _mm512_mul_pd(vTwo,
                                            _mm512_fmadd_pd(valpha1, valpha1, _mm512_mul_pd(valpha, valpha2)));

        vLogMgf0 = _mm512_add_pd(vLogMgf0, avx512_log_pd(vmgf0));
        const __m512d vratio = _mm512_div_pd(vmgf1, vmgf0);
        vTempSum = _mm512_add_pd(vTempSum, vratio);
        vS1      = _mm512_add_pd(vS1, _mm512_div_pd(vmgf2, vmgf0));
        vS2      = _mm512_fmadd_pd(vratio, vratio, vS2);
    }

    return {_mm512_reduce_add_pd(vLogMgf0),
            _mm512_reduce_add_pd(vTempSum),
            _mm512_reduce_add_pd(vS1) - _mm512_reduce_add_pd(vS2)};
}

__attribute__((target("avx2,fma")))
static ScalarCgfResult simdOutlierCgf_avx2(
    double t,
    const double *outlierResid,
    int nOutlier,
    double MAF
) {
    const __m256d vt       = _mm256_set1_pd(t);
    const __m256d vMAF     = _mm256_set1_pd(MAF);
    const __m256d vOneMmaf = _mm256_set1_pd(1.0 - MAF);
    const __m256d vTwo     = _mm256_set1_pd(2.0);

    __m256d vLogMgf0 = _mm256_setzero_pd();
    __m256d vTempSum = _mm256_setzero_pd();
    __m256d vS1      = _mm256_setzero_pd();
    __m256d vS2      = _mm256_setzero_pd();

    int i = 0;
    const int n4 = nOutlier & ~3;
    for (; i < n4; i += 4) {
        const __m256d vr      = _mm256_loadu_pd(outlierResid + i);
        const __m256d vlambda = avx2_exp_pd(_mm256_mul_pd(vt, vr));
        const __m256d valpha  = _mm256_fmadd_pd(vMAF, vlambda, vOneMmaf);
        const __m256d valpha1 = _mm256_mul_pd(vMAF, _mm256_mul_pd(vr, vlambda));
        const __m256d valpha2 = _mm256_mul_pd(vr, valpha1);

        const __m256d vmgf0 = _mm256_mul_pd(valpha, valpha);
        const __m256d vmgf1 = _mm256_mul_pd(vTwo, _mm256_mul_pd(valpha, valpha1));
        const __m256d vmgf2 = _mm256_mul_pd(vTwo,
                                            _mm256_fmadd_pd(valpha1, valpha1, _mm256_mul_pd(valpha, valpha2)));

        vLogMgf0 = _mm256_add_pd(vLogMgf0, avx2_log_pd(vmgf0));
        const __m256d vratio = _mm256_div_pd(vmgf1, vmgf0);
        vTempSum = _mm256_add_pd(vTempSum, vratio);
        vS1      = _mm256_add_pd(vS1, _mm256_div_pd(vmgf2, vmgf0));
        vS2      = _mm256_fmadd_pd(vratio, vratio, vS2);
    }

    // Horizontal reduction of AVX2 accumulators
    auto hsum = [](__m256d v) -> double {
        __m128d lo = _mm256_castpd256_pd128(v);
        __m128d hi = _mm256_extractf128_pd(v, 1);
        lo = _mm_add_pd(lo, hi);
        return _mm_cvtsd_f64(lo) + _mm_cvtsd_f64(_mm_unpackhi_pd(lo, lo));
    };
    double logMgf0Sum = hsum(vLogMgf0);
    double tempSum    = hsum(vTempSum);
    double s1         = hsum(vS1);
    double s2         = hsum(vS2);

    // Scalar tail (1-3 remaining outliers)
    const double oneMmaf = 1.0 - MAF;
    for (; i < nOutlier; ++i) {
        const double r      = outlierResid[i];
        const double lambda = std::exp(t * r);
        const double alpha  = oneMmaf + MAF * lambda;
        const double alpha1 = MAF * r * lambda;
        const double alpha2 = r * alpha1;

        const double mgf0   = alpha * alpha;
        const double mgf1   = 2.0 * alpha * alpha1;
        const double mgf2   = 2.0 * (alpha1 * alpha1 + alpha * alpha2);

        logMgf0Sum += std::log(mgf0);
        const double ratio  = mgf1 / mgf0;
        tempSum += ratio;
        s1      += mgf2 / mgf0;
        s2      += ratio * ratio;
    }

    return {logMgf0Sum, tempSum, s1 - s2};
}

// ── Dispatch: pick fastest available CGF implementation ─────────────
using OutlierCgfFn = ScalarCgfResult (*)(double, const double *, int, double);

static OutlierCgfFn pickOutlierCgfFn() {
    switch (simdLevel()) {
    case SimdLevel::AVX512: return simdOutlierCgf_avx512;
    case SimdLevel::AVX2:   return simdOutlierCgf_avx2;
    default:                return scalarOutlierCgf;
    }
}

static const OutlierCgfFn outlierCgf = pickOutlierCgfFn();

inline double signOf(double x) {
    return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);
}

// Scalar Newton-Raphson root finder (equivalent to nsSPAGRM::fastGetRoot
// for outlier-only case).  Returns converged t; writes final CGF result.
inline double scalarFastGetRoot(
    const double *outlierResid,
    int nOutlier,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score,
    double MAF,
    double init_t,
    double tol,
    ScalarCgfResult &finalCgf,
    int maxiter = 50
) {
    double t = init_t;
    double CGF1 = 0.0, CGF2 = 0.0;
    double diff_t = std::numeric_limits<double>::infinity();

    const double mean = 2.0 * MAF * sum_R_nonOutlier;
    const double var  = 2.0 * MAF * (1.0 - MAF) * R_GRM_R_nonOutlier;

    ScalarCgfResult cgf{};

    for (int iter = 1; iter < maxiter; ++iter) {
        const double old_t      = t;
        const double old_diff_t = diff_t;
        const double old_CGF1   = CGF1;

        cgf = outlierCgf(t, outlierResid, nOutlier, MAF);

        CGF1 = cgf.cgf1Outlier + mean + var * t - Score;
        CGF2 = cgf.cgf2Outlier + var;

        diff_t = -CGF1 / CGF2;

        if (std::isnan(diff_t) || std::isinf(CGF2)) {
            t = t / 2.0;
            diff_t = std::min(std::abs(t), 1.0) * signOf(Score);
            continue;
        }

        if (std::isnan(old_CGF1) || (signOf(old_CGF1) != 0 && signOf(CGF1) != signOf(old_CGF1))) {
            if (std::abs(diff_t) < tol) {
                t = old_t + diff_t;
                break;
            } else {
                for (int bisect = 0; bisect < 200 && std::abs(old_diff_t) > tol &&
                     std::abs(diff_t) > std::abs(old_diff_t) - tol; ++bisect)
                    diff_t /= 2.0;
                t = old_t + diff_t;
                continue;
            }
        }

        if (signOf(Score) != signOf(old_t + diff_t) && (signOf(old_CGF1) == 0 || signOf(CGF1) == signOf(old_CGF1))) {
            for (int bisect = 0; bisect < 200 && signOf(Score) != signOf(old_t + diff_t); ++bisect)
                diff_t /= 2.0;
            t = old_t + diff_t;
            continue;
        }

        t = old_t + diff_t;
        if (std::abs(diff_t) < tol) break;
    }

    // Final evaluation at converged t.
    finalCgf = outlierCgf(t, outlierResid, nOutlier, MAF);
    return t;
}

// Scalar SPA tail probability (equivalent to nsSPAGRM::getProbSpa for
// outlier-only case).
inline double scalarGetProbSpa(
    const double *outlierResid,
    int nOutlier,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score,
    double MAF,
    bool lower_tail,
    double zeta,
    double tol
) {
    ScalarCgfResult cgf{};
    zeta = scalarFastGetRoot(outlierResid, nOutlier, sum_R_nonOutlier, R_GRM_R_nonOutlier,
                             Score, MAF, zeta, tol, cgf);

    const double mean = 2.0 * MAF * sum_R_nonOutlier;
    const double var  = 2.0 * MAF * (1.0 - MAF) * R_GRM_R_nonOutlier;

    const double CGF0 = cgf.logMgf0Sum + mean * zeta + 0.5 * var * zeta * zeta;
    const double CGF2 = cgf.cgf2Outlier + var;

    const double w_val = signOf(zeta) * std::sqrt(2.0 * (zeta * Score - CGF0));
    const double v_val = zeta * std::sqrt(CGF2);
    const double u     = w_val + (1.0 / w_val) * std::log(v_val / w_val);

    return math::pnorm(u, 0.0, 1.0, lower_tail);
}

// ══════════════════════════════════════════════════════════════════════
// Per-tau SPA data for SPAsqr (shared across clones via shared_ptr)
// ══════════════════════════════════════════════════════════════════════

struct SPAsqrPerTau {
    Eigen::VectorXd outlierResid;        // residuals of outlier subjects
    double sum_unrelated_outliers2;      // outlierResid.squaredNorm()
    double sum_R_nonOutlier;
    double R_GRM_R_nonOutlier;
    double R_GRM_R;                      // full variance term (all subjects)
};

struct SPAsqrSPAShared {
    std::vector<SPAsqrPerTau> perTau;    // one per tau
    double SPA_Cutoff;
    double zeta;                         // initial guess for Newton-Raphson
    double tol;
};

// ── Per-tau p-value (replaces SPAGRMClass::getMarkerPvalFromScore) ──
inline double scalarGetPvalFromScore(
    double Score,
    double altFreq,
    double &zScore,
    const SPAsqrPerTau &tau,
    const SPAsqrSPAShared &spa
) {
    const double MAF       = std::min(altFreq, 1.0 - altFreq);
    const double G_var     = 2.0 * MAF * (1.0 - MAF);
    const double Score_var = G_var * tau.R_GRM_R;

    if (Score_var <= 0.0 || MAF <= 0.0) {
        zScore = 0.0;
        return std::numeric_limits<double>::quiet_NaN();
    }

    zScore = Score / std::sqrt(Score_var);

    if (!std::isfinite(zScore))
        return std::numeric_limits<double>::quiet_NaN();

    if (std::abs(zScore) <= spa.SPA_Cutoff)
        return 2.0 * math::pnorm(std::abs(zScore), 0.0, 1.0, false);

    // EmpVar: SPAsqr has no TwoSubj/ThreeSubj, so simplifies to:
    //   G_var * (R_GRM_R_nonOutlier + sum_unrelated_outliers2)
    const double EmpVar    = G_var * (tau.R_GRM_R_nonOutlier + tau.sum_unrelated_outliers2);
    const double Score_adj = Score * std::sqrt(EmpVar / Score_var);

    double zeta1 = std::abs(Score_adj) / Score_var;
    zeta1 = std::min(zeta1, 1.2);
    const double zeta2 = -zeta1;

    const double *oResid  = tau.outlierResid.data();
    const int nOutlier = static_cast<int>(tau.outlierResid.size());

    const double pval1 = scalarGetProbSpa(oResid, nOutlier, tau.sum_R_nonOutlier, tau.R_GRM_R_nonOutlier,
                                          std::abs(Score_adj), MAF, false, zeta1, 1e-4);
    const double pval2 = scalarGetProbSpa(oResid, nOutlier, tau.sum_R_nonOutlier, tau.R_GRM_R_nonOutlier,
                                          -std::abs(Score_adj), MAF, true, zeta2, spa.tol);
    return pval1 + pval2;
}

// ── C3: Fused upper+lower SPA — shared setup, two Newtons ──────────
// Computes the SPA p-value for a marker that failed the normal-approx
// cutoff.  Both tails share the MAF/variance setup and outlier data
// pointer extraction; only the Newton root and saddlepoint differ.
inline double fusedGetPvalSpa(
    double Score,
    double altFreq,
    const SPAsqrPerTau &tau,
    const SPAsqrSPAShared &spa
) {
    const double MAF       = std::min(altFreq, 1.0 - altFreq);
    const double G_var     = 2.0 * MAF * (1.0 - MAF);
    const double Score_var = G_var * tau.R_GRM_R;

    const double EmpVar    = G_var * (tau.R_GRM_R_nonOutlier + tau.sum_unrelated_outliers2);
    const double Score_adj = Score * std::sqrt(EmpVar / Score_var);
    const double absAdj    = std::abs(Score_adj);

    double zeta1 = std::min(absAdj / Score_var, 1.2);
    const double zeta2 = -zeta1;

    const double *oResid   = tau.outlierResid.data();
    const int nOutlier = static_cast<int>(tau.outlierResid.size());
    const double sumNO    = tau.sum_R_nonOutlier;
    const double rgrNO    = tau.R_GRM_R_nonOutlier;

    // Mean and variance of non-outlier contribution (shared)
    const double mean = 2.0 * MAF * sumNO;
    const double var  = 2.0 * MAF * (1.0 - MAF) * rgrNO;

    // ── Upper tail Newton ──────────────────────────────────────────
    ScalarCgfResult cgfU{};
    const double tU = scalarFastGetRoot(oResid, nOutlier, sumNO, rgrNO,
                                        absAdj, MAF, zeta1, 1e-4, cgfU);
    const double CGF0_U = cgfU.logMgf0Sum + mean * tU + 0.5 * var * tU * tU;
    const double CGF2_U = cgfU.cgf2Outlier + var;
    const double wU = signOf(tU) * std::sqrt(2.0 * (tU * absAdj - CGF0_U));
    const double vU = tU * std::sqrt(CGF2_U);
    const double uU = wU + (1.0 / wU) * std::log(vU / wU);
    const double pval1 = math::pnorm(uU, 0.0, 1.0, false);

    // ── Lower tail Newton ──────────────────────────────────────────
    ScalarCgfResult cgfL{};
    const double ScoreL = -absAdj;
    const double tL = scalarFastGetRoot(oResid, nOutlier, sumNO, rgrNO,
                                        ScoreL, MAF, zeta2, spa.tol, cgfL);
    const double CGF0_L = cgfL.logMgf0Sum + mean * tL + 0.5 * var * tL * tL;
    const double CGF2_L = cgfL.cgf2Outlier + var;
    const double wL = signOf(tL) * std::sqrt(2.0 * (tL * ScoreL - CGF0_L));
    const double vL = tL * std::sqrt(CGF2_L);
    const double uL = wL + (1.0 / wL) * std::log(vL / wL);
    const double pval2 = math::pnorm(uL, 0.0, 1.0, true);

    return pval1 + pval2;
}

// ══════════════════════════════════════════════════════════════════════
// SPAsqrMethod — MethodBase adapter with zero per-clone heap allocation
//
// All per-tau SPA data is shared via shared_ptr (read-only).
// Cloning copies two shared_ptr's — no MgfWorkspace vectors allocated.
// ══════════════════════════════════════════════════════════════════════

class SPAsqrMethod : public MethodBase {
  public:
    SPAsqrMethod(
        int ntaus,
        std::shared_ptr<const SPAsqrSPAShared> spaShared,
        Eigen::MatrixXd residMat,
        Eigen::VectorXd residSums,
        std::vector<std::string> tauLabels
    )
        : m_ntaus(ntaus),
          m_spaShared(std::move(spaShared)),
          m_tauLabels(std::move(tauLabels)),
          m_hasLabels(true)
    {
        auto sd = std::make_shared<SharedMethodData>();
        sd->residMat  = std::move(residMat);
        sd->residSums = std::move(residSums);
        m_methodShared = std::move(sd);
    }

    std::unique_ptr<MethodBase> clone() const override {
        // Both shared_ptr's are read-only — no per-clone scratch at all.
        return std::make_unique<SPAsqrMethod>(*this);
    }

    int resultSize() const override {
        return 2 * m_ntaus + 1;
    }

    std::string getHeaderColumns() const override {
        std::ostringstream oss;
        if (m_hasLabels) {
            oss << "\tP_CCT";
            for (int i = 0; i < m_ntaus; ++i)
                oss << "\tP_" << m_tauLabels[i];
            for (int i = 0; i < m_ntaus; ++i)
                oss << "\tZ_" << m_tauLabels[i];
        } else {
            oss << "\tP_CCT";
            for (int i = 1; i <= m_ntaus; ++i)
                oss << "\tP_tau" << i;
            for (int i = 1; i <= m_ntaus; ++i)
                oss << "\tZ_tau" << i;
        }
        return oss.str();
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int /*markerInChunkIdx*/,
        std::vector<double> &result
    ) override {
        result.clear();
        result.reserve(2 * m_ntaus + 1);

        const double gMean = GVec.mean();
        Eigen::VectorXd scores = m_methodShared->residMat.transpose() * GVec;
        for (int i = 0; i < m_ntaus; ++i)
            scores[i] -= gMean * m_methodShared->residSums[i];

        processOneMarker(scores.data(), altFreq, result);
    }

    // ── Batched analysis: B markers at once ────────────────────────────
    void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        std::vector<std::vector<double> > &results
    ) override {
        const int B = static_cast<int>(GBatch.cols());
        results.resize(B);

        Eigen::MatrixXd scoreMatrix;
        scoreMatrix.noalias() = m_methodShared->residMat.transpose() * GBatch;

        const Eigen::VectorXd gMeans = GBatch.colwise().mean();
        scoreMatrix.noalias() -= m_methodShared->residSums * gMeans.transpose();

        for (int b = 0; b < B; ++b)
            processOneMarker(scoreMatrix.col(b).data(), altFreqs[b], results[b]);
    }

    int preferredBatchSize() const override {
        return std::min(std::max(4, 2 * m_ntaus), 16);
    }

    // ── Fused union-level GEMM interface ───────────────────────────────
    bool supportsFusedGemm() const override {
        return true;
    }

    int fusedGemmColumns() const override {
        return m_ntaus;
    }

    void fillUnionResiduals(
        Eigen::Ref<Eigen::MatrixXd> dest,
        const std::vector<uint32_t> &unionToLocal
    ) const override {
        const auto &residMat = m_methodShared->residMat;
        const uint32_t nUnion = static_cast<uint32_t>(unionToLocal.size());
        for (uint32_t i = 0; i < nUnion; ++i) {
            const uint32_t li = unionToLocal[i];
            if (li != UINT32_MAX)
                dest.row(i) = residMat.row(li);
        }
    }

    void fillResidualSums(double *dest) const override {
        const auto &rs = m_methodShared->residSums;
        std::copy(rs.data(), rs.data() + m_ntaus, dest);
    }

    void processScoreBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &scores,
        const double *gSums,
        uint32_t nUsed,
        const std::vector<double> &altFreqs,
        std::vector<std::vector<double> > &results
    ) override {
        const int B = static_cast<int>(scores.cols());
        results.resize(B);

        // ── Center scores ──────────────────────────────────────────
        m_centeredBuf = scores;
        const double invN = 1.0 / static_cast<double>(nUsed);
        for (int b = 0; b < B; ++b) {
            const double gMean = gSums[b] * invN;
            for (int t = 0; t < m_ntaus; ++t)
                m_centeredBuf(t, b) -= gMean * m_methodShared->residSums[t];
        }

        // ── C2: Tau-first SPA computation ──────────────────────────
        // Process all B markers for one tau before moving to the next.
        // This keeps the tau's outlierResid data hot in cache.
        // Layout: pBuf[b * ntaus + t], zBuf[b * ntaus + t]
        double zBuf[20 * 16];   // B(≤16) × ntaus(≤20)
        double pBuf[20 * 16];

        for (int t = 0; t < m_ntaus; ++t) {
            const SPAsqrPerTau &tau = m_spaShared->perTau[t];
            const double cutoff = m_spaShared->SPA_Cutoff;

            for (int b = 0; b < B; ++b) {
                const double Score   = m_centeredBuf(t, b);
                const double altFreq = altFreqs[b];
                const double MAF     = std::min(altFreq, 1.0 - altFreq);
                const double G_var   = 2.0 * MAF * (1.0 - MAF);
                const double Svar    = G_var * tau.R_GRM_R;

                double z, p;
                if (Svar <= 0.0 || MAF <= 0.0) {
                    z = 0.0;
                    p = std::numeric_limits<double>::quiet_NaN();
                } else {
                    z = Score / std::sqrt(Svar);
                    if (!std::isfinite(z)) {
                        p = std::numeric_limits<double>::quiet_NaN();
                    } else if (std::abs(z) <= cutoff) {
                        // Fast path: normal approximation
                        p = 2.0 * math::pnorm(std::abs(z), 0.0, 1.0, false);
                    } else {
                        // Slow path: C3 fused upper+lower SPA
                        p = fusedGetPvalSpa(Score, altFreq, tau, *m_spaShared);
                    }
                }
                zBuf[b * m_ntaus + t] = z;
                pBuf[b * m_ntaus + t] = p;
            }
        }

        // ── Assemble results with CCT ──────────────────────────────
        for (int b = 0; b < B; ++b) {
            auto &result = results[b];
            result.clear();
            result.reserve(2 * m_ntaus + 1);

            const double *pvals   = pBuf + b * m_ntaus;
            const double *zScores = zBuf + b * m_ntaus;

            // CCT (Cauchy combination test)
            double valid_p[20];
            int nValid = 0;
            for (int i = 0; i < m_ntaus; ++i)
                if (!std::isnan(pvals[i])) valid_p[nValid++] = pvals[i];

            double pCCT = std::numeric_limits<double>::quiet_NaN();
            if (nValid > 0) {
                bool hasZero = false;
                double tStat = 0.0;
                for (int vi = 0; vi < nValid; ++vi) {
                    double p = valid_p[vi];
                    if (p <= 0.0) { hasZero = true; break; }
                    double pc = (p >= 1.0) ? 0.999 : p;
                    tStat += (pc < 1e-15) ? (1.0 / (pc * M_PI))
                                          : std::tan((0.5 - pc) * M_PI);
                }
                if (hasZero) {
                    pCCT = 0.0;
                } else {
                    tStat /= static_cast<double>(nValid);
                    pCCT = (tStat > 1e15) ? (1.0 / tStat) / M_PI
                                           : 0.5 - std::atan(tStat) / M_PI;
                }
            }

            result.push_back(pCCT);
            for (int i = 0; i < m_ntaus; ++i)
                result.push_back(pvals[i]);
            for (int i = 0; i < m_ntaus; ++i)
                result.push_back(zScores[i]);
        }
    }

  private:
    int m_ntaus;
    std::shared_ptr<const SPAsqrSPAShared> m_spaShared;
    std::vector<std::string> m_tauLabels;
    bool m_hasLabels;

    struct SharedMethodData {
        Eigen::MatrixXd residMat;   // N × ntaus
        Eigen::VectorXd residSums;  // ntaus
    };

    std::shared_ptr<const SharedMethodData> m_methodShared;
    Eigen::MatrixXd m_centeredBuf;  // reused across processScoreBatch calls

    void processOneMarker(
        const double *centeredScores,
        double altFreq,
        std::vector<double> &result
    ) {
        result.clear();
        result.reserve(2 * m_ntaus + 1);

        // Stack arrays — ntaus is always small (≤20).
        double zScores[20];
        double pvals[20];

        for (int i = 0; i < m_ntaus; ++i) {
            double z;
            double p = scalarGetPvalFromScore(centeredScores[i], altFreq, z,
                                              m_spaShared->perTau[i], *m_spaShared);
            zScores[i] = z;
            pvals[i] = p;
        }

        // CCT (Cauchy combination test) p-value
        double valid_p[20];
        int nValid = 0;
        for (int i = 0; i < m_ntaus; ++i)
            if (!std::isnan(pvals[i])) valid_p[nValid++] = pvals[i];

        double pCCT = std::numeric_limits<double>::quiet_NaN();
        if (nValid > 0) {
            bool hasZero = false;
            double tStat = 0.0;
            for (int vi = 0; vi < nValid; ++vi) {
                double p = valid_p[vi];
                if (p <= 0.0) {
                    hasZero = true;
                    break;
                }
                double pc = (p >= 1.0) ? 0.999 : p;
                tStat += (pc < 1e-15) ? (1.0 / (pc * M_PI))
                                      : std::tan((0.5 - pc) * M_PI);
            }
            if (hasZero) {
                pCCT = 0.0;
            } else {
                tStat /= static_cast<double>(nValid);
                pCCT = (tStat > 1e15) ? (1.0 / tStat) / M_PI : 0.5 - std::atan(tStat) / M_PI;
            }
        }

        result.push_back(pCCT);
        for (int i = 0; i < m_ntaus; ++i)
            result.push_back(pvals[i]);
        for (int i = 0; i < m_ntaus; ++i)
            result.push_back(zScores[i]);
    }

};

// ══════════════════════════════════════════════════════════════════════
// Outlier detection (IQR-based, per column)
// ══════════════════════════════════════════════════════════════════════

// Returns an N × K boolean matrix (as std::vector<std::vector<bool>>).
// outlierIqrRatio  = multiplier for IQR (default 1.5)
// outlierAbsBound  = absolute clamp for cutoffs (default 0.55)
struct OutlierInfo {
    // per-column: indices of outlier subjects
    std::vector<std::vector<int> > outlierIdx;
    // per-column: boolean mask (size N)
    std::vector<std::vector<bool> > isOutlier;
};

OutlierInfo detectOutliers(
    const Eigen::MatrixXd &ResidMat,
    double outlierIqrRatio,
    double outlierAbsBound
) {
    const Eigen::Index N = ResidMat.rows();
    const Eigen::Index K = ResidMat.cols();

    OutlierInfo info;
    info.outlierIdx.resize(K);
    info.isOutlier.resize(K);

    // Scratch for sorting
    std::vector<double> scratch(N);

    for (Eigen::Index col = 0; col < K; ++col) {
        // Copy column for sorting
        for (Eigen::Index i = 0; i < N; ++i)
            scratch[i] = ResidMat(i, col);

        std::sort(scratch.begin(), scratch.end());

        // Q1, Q3 via linear interpolation (same as R type=7)
        auto quantile = [&](double prob) -> double {
            const double idx = prob * (N - 1);
            const Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
            const Eigen::Index hi = std::min(lo + 1, N - 1);
            const double frac = idx - lo;
            return scratch[lo] * (1.0 - frac) + scratch[hi] * frac;
        };

        const double Q1 = quantile(0.25);
        const double Q3 = quantile(0.75);
        const double IQR = Q3 - Q1;

        double cutLo = Q1 - outlierIqrRatio * IQR;
        double cutHi = Q3 + outlierIqrRatio * IQR;
        cutLo = std::max(cutLo, -outlierAbsBound);
        cutHi = std::min(cutHi, outlierAbsBound);

        info.isOutlier[col].resize(N, false);
        int nOutlier = 0;
        for (Eigen::Index i = 0; i < N; ++i) {
            const double v = ResidMat(i, col);
            if (v < cutLo || v > cutHi) {
                info.isOutlier[col][i] = true;
                info.outlierIdx[col].push_back(static_cast<int>(i));
                ++nOutlier;
            }
        }
    }
    return info;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// Shared pipeline: outlier detection → GRM → SPAGRM → marker engine
// ══════════════════════════════════════════════════════════════════════

// GRMEntry is now declared in spasqr.hpp

// Load GRM entries from disk and convert to flat GRMEntry vector.
std::vector<GRMEntry> loadGrmEntries(
    const std::vector<std::string> &subjOrder,
    const std::vector<std::string> &famIIDs,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile
) {
    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, subjOrder, famIIDs);
    infoMsg("Sparse GRM: %zu entries used", grm.nnz());
    std::vector<GRMEntry> entries;
    entries.reserve(grm.nnz());
    for (const auto &e : grm.entries())
        entries.push_back({e.row, e.col, e.value, (e.row == e.col) ? 1.0 : 2.0});
    return entries;
}

// Build SPAsqrMethod from a pre-computed residual matrix and pre-loaded GRM entries.
std::unique_ptr<MethodBase> buildSPAsqrMethod(
    Eigen::MatrixXd &ResidMat,
    const std::vector<GRMEntry> &grmEntries,
    uint32_t nUsed,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    double /*minMafCutoff*/,
    double /*minMacCutoff*/,
    std::vector<std::string> tauLabels,
    std::vector<double> *outlierRatiosOut
)
{
    const Eigen::Index N = ResidMat.rows();
    const Eigen::Index K = ResidMat.cols();
    const int ntaus = static_cast<int>(K);

    // ── 3. Outlier detection ───────────────────────────────────────────
    OutlierInfo outlierInfo = detectOutliers(ResidMat, outlierIqrRatio, outlierAbsBound);

    // Populate outlier ratios for caller to format
    if (outlierRatiosOut) {
        outlierRatiosOut->resize(K);
        for (Eigen::Index c = 0; c < K; ++c)
            (*outlierRatiosOut)[c] = static_cast<double>(outlierInfo.outlierIdx[c].size()) / N;
    }

    // ── 4. Compute per-column variance terms + build SPAsqrSPAShared ──
    auto spaShared = std::make_shared<SPAsqrSPAShared>();
    spaShared->perTau.resize(ntaus);
    spaShared->SPA_Cutoff = spaCutoff;
    spaShared->zeta       = 0.01;
    spaShared->tol        = 1e-6;

    for (int col = 0; col < ntaus; ++col) {
        const auto &isOut = outlierInfo.isOutlier[col];
        auto &pt = spaShared->perTau[col];

        // R_GRM_R and R_GRM_R_nonOutlier
        double rgrm_r = 0.0;
        double rgrm_r_no = 0.0;
        for (const auto &e : grmEntries) {
            const double contrib = e.factor * e.value * ResidMat(e.row, col) * ResidMat(e.col, col);
            rgrm_r += contrib;
            if (!isOut[e.row] && !isOut[e.col]) rgrm_r_no += contrib;
        }
        pt.R_GRM_R            = rgrm_r;
        pt.R_GRM_R_nonOutlier = rgrm_r_no;

        // sum_R_nonOutlier and outlier residual values
        double sumNO = 0.0;
        std::vector<double> outVals;
        outVals.reserve(outlierInfo.outlierIdx[col].size());
        for (uint32_t i = 0; i < nUsed; ++i) {
            if (!isOut[i]) sumNO += ResidMat(i, col);
            else outVals.push_back(ResidMat(i, col));
        }
        pt.sum_R_nonOutlier = sumNO;
        pt.outlierResid =
            Eigen::Map<Eigen::VectorXd>(outVals.data(), static_cast<Eigen::Index>(outVals.size()));
        pt.sum_unrelated_outliers2 = pt.outlierResid.squaredNorm();
    }

    // ── 5. Build residMat / residSums for fused dot products ──────────
    Eigen::VectorXd residSums(ntaus);
    for (int t = 0; t < ntaus; ++t)
        residSums[t] = ResidMat.col(t).sum();

    // ── 6. Build method ──────────────────────────────────────────────
    std::unique_ptr<SPAsqrMethod> method = std::make_unique<SPAsqrMethod>(
        ntaus,
        std::move(spaShared),
        Eigen::MatrixXd(ResidMat),   // copy — caller may reuse ResidMat
        std::move(residSums),
        std::move(tauLabels)
    );

    return method;
}

// ══════════════════════════════════════════════════════════════════════
// Re-index GRM entries from union-dense space to pheno-dense space
// ══════════════════════════════════════════════════════════════════════

static std::vector<GRMEntry> reindexGrm(
    const std::vector<GRMEntry> &unionGrm,
    const std::vector<uint32_t> &unionToLocal, // union index → pheno index (UINT32_MAX = absent)
    uint32_t nUnion
) {
    std::vector<GRMEntry> out;
    out.reserve(unionGrm.size());
    for (const auto &e : unionGrm) {
        if (e.row >= nUnion || e.col >= nUnion) continue;
        uint32_t lr = unionToLocal[e.row];
        uint32_t lc = unionToLocal[e.col];
        if (lr == UINT32_MAX || lc == UINT32_MAX) continue;
        out.push_back({lr, lc, e.value, e.factor});
    }
    return out;
}

// ══════════════════════════════════════════════════════════════════════
// runSPAsqr — multi-phenotype entry point with parallel conquer fits
// ══════════════════════════════════════════════════════════════════════

void runSPAsqr(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<double> &taus,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    double spasqrTol,
    double spasqrH,
    double spasqrHScale,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const int K = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());

    // ── 1. Load phenotype/covariate data (union mask) ───────────────
    // Union = subjects with genotype ∩ GRM ∩ keep/remove.  Per-phenotype
    // NA filtering is deferred — each phenotype uses its own non-missing
    // subset of the union.
    infoMsg("SPAsqr: Loading phenotype and covariate data (%d phenotypes, %d taus)", K, ntaus);
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile, phenoNames);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();

    const uint32_t nUnion = sd.nUsed();
    const Eigen::Index N = static_cast<Eigen::Index>(nUnion);

    // ── 2. Extract union-space covariates ───────────────────────────
    Eigen::MatrixXd unionX = covarNames.empty()
        ? (sd.hasCovar() ? Eigen::MatrixXd(sd.covar()) : Eigen::MatrixXd(N, 0))
        : sd.getColumns(covarNames);
    const int nCov = static_cast<int>(unionX.cols());

    // ── 3. Per-phenotype: build non-missing mask, extract Y/X ───────
    struct PhenoWork {
        std::vector<uint32_t> unionToLocal; // size nUnion; UINT32_MAX = absent
        uint32_t nk;                        // non-missing count
        Eigen::VectorXd Y;                  // nk
        Eigen::MatrixXd X;                  // nk × nCov
        double h;                           // bandwidth
        Eigen::MatrixXd ResidMat;           // nk × ntaus (filled by conquer)
    };

    std::vector<PhenoWork> pw(K);

    for (int k = 0; k < K; ++k) {
        Eigen::VectorXd fullY = sd.getColumn(phenoNames[k]);
        pw[k].unionToLocal.resize(nUnion, UINT32_MAX);
        uint32_t localIdx = 0;
        for (uint32_t i = 0; i < nUnion; ++i) {
            if (!std::isnan(fullY[i])) {
                pw[k].unionToLocal[i] = localIdx++;
            }
        }
        pw[k].nk = localIdx;
        if (pw[k].nk == 0)
            throw std::runtime_error("SPAsqr: phenotype '" + phenoNames[k] + "' has no non-missing subjects");

        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        pw[k].Y.resize(Nk);
        pw[k].X.resize(Nk, nCov);
        for (uint32_t i = 0; i < nUnion; ++i) {
            uint32_t li = pw[k].unionToLocal[i];
            if (li == UINT32_MAX) continue;
            pw[k].Y[li] = fullY[i];
            if (nCov > 0) pw[k].X.row(li) = unionX.row(i);
        }

        // Per-phenotype bandwidth
        if (spasqrH >= 0.0) {
            pw[k].h = spasqrH;
        } else {
            std::vector<double> ysorted(Nk);
            Eigen::VectorXd::Map(ysorted.data(), Nk) = pw[k].Y;
            std::sort(ysorted.begin(), ysorted.end());
            auto quantile = [&](double prob) -> double {
                double idx = prob * (Nk - 1);
                Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
                Eigen::Index hi = std::min(lo + 1, Nk - 1);
                double frac = idx - lo;
                return ysorted[lo] * (1.0 - frac) + ysorted[hi] * frac;
            };
            double iqr = quantile(0.75) - quantile(0.25);
            double scale = (spasqrHScale >= 0.0) ? spasqrHScale : 3.0;
            pw[k].h = iqr / scale;
            if (pw[k].h <= 0.0)
                pw[k].h = std::max(std::pow((std::log(Nk) + nCov) / static_cast<double>(Nk), 0.4), 0.05);
        }

        pw[k].ResidMat.resize(Nk, ntaus);
    }

    // Log N/bandwidth table
    {
        infoMsg("Sample size and smooth bandwidth per phenotype:");
        size_t nameW = 4;
        for (int k = 0; k < K; ++k)
            nameW = std::max(nameW, phenoNames[k].size());
        nameW += 2;
        char row[256];
        std::snprintf(row, sizeof(row), "    %-*s %10s %12s\n", static_cast<int>(nameW), "", "N", "bandwidth");
        fprintf(stderr, "%s", row);
        for (int k = 0; k < K; ++k) {
            std::snprintf(row, sizeof(row), "    %-*s %10u %12.6f\n",
                          static_cast<int>(nameW), phenoNames[k].c_str(), pw[k].nk, pw[k].h);
            fprintf(stderr, "%s", row);
        }
    }

    // ── 4. Parallel conquer fits: K × ntaus ─────────────────────────
    const int totalFits = K * ntaus;
    const int nWorkers = std::min(nthreads, totalFits);
    infoMsg("SPAsqr: Running %d conquer fits with %d threads", totalFits, nWorkers);

    std::atomic<int> nextFit{0};
    std::vector<std::string> fitErrors(totalFits);

    auto fitWorker = [&]() {
        for (;;) {
            int idx = nextFit.fetch_add(1, std::memory_order_relaxed);
            if (idx >= totalFits) break;
            int k = idx / ntaus;
            int t = idx % ntaus;
            const double h1 = 1.0 / pw[k].h;

            try {
                Eigen::VectorXd resid;
                Eigen::VectorXd beta = conquer::smqrGauss(pw[k].X, pw[k].Y, taus[t], pw[k].h, &resid, spasqrTol);
                infoMsg("[%s] tau=%.4f intercept=%.6f", phenoNames[k].c_str(), taus[t], beta(0));

                const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
                for (Eigen::Index i = 0; i < Nk; ++i)
                    pw[k].ResidMat(i, t) = taus[t] - math::pnorm(-resid(i) * h1);
            } catch (const std::exception &ex) {
                fitErrors[idx] = ex.what();
            }
        }
    };

    {
        std::vector<std::thread> threads;
        threads.reserve(nWorkers - 1);
        for (int t = 0; t < nWorkers - 1; ++t)
            threads.emplace_back(fitWorker);
        fitWorker();
        for (auto &th : threads)
            th.join();
    }

    for (int idx = 0; idx < totalFits; ++idx) {
        if (!fitErrors[idx].empty()) {
            int k = idx / ntaus;
            int t = idx % ntaus;
            throw std::runtime_error("SPAsqr: conquer failed for phenotype '" +
                                     phenoNames[k] + "' tau=" + std::to_string(taus[t]) + ": " + fitErrors[idx]);
        }
    }

    // ── 5. Build tau labels ────────────────────────────────────────────
    std::vector<std::string> tauLabels;
    tauLabels.reserve(ntaus);
    for (double tau : taus) {
        std::ostringstream oss;
        oss << "tau" << tau;
        tauLabels.push_back(oss.str());
    }

    // ── 6. Load genotype data and GRM once (shared, union space) ────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    std::vector<GRMEntry> unionGrm = loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), spgrmGrabFile, spgrmGctaFile);

    // ── 7. Build per-phenotype SPAsqrMethod + PhenoTask ─────────────
    std::vector<PhenoTask> tasks(K);
    std::vector<std::vector<double> > allOutlierRatios(K);

    for (int k = 0; k < K; ++k) {
        infoMsg("[%s] Building SPAsqr method (%d taus, %u subjects)",
                phenoNames[k].c_str(), ntaus, pw[k].nk);

        // Re-index GRM to pheno-dense space
        auto phenoGrm = reindexGrm(unionGrm, pw[k].unionToLocal, nUnion);

        auto method = buildSPAsqrMethod(
            pw[k].ResidMat,
            phenoGrm,
            pw[k].nk,
            spaCutoff,
            outlierIqrRatio,
            outlierAbsBound,
            minMafCutoff,
            minMacCutoff,
            tauLabels,
            &allOutlierRatios[k]
        );

        tasks[k].phenoName = phenoNames[k];
        tasks[k].method = std::move(method);
        tasks[k].unionToLocal = pw[k].unionToLocal;
        tasks[k].nUsed = pw[k].nk;
    }

    // Print outlier ratio table line-by-line
    {
        size_t nameW = 4;
        for (int k = 0; k < K; ++k)
            nameW = std::max(nameW, phenoNames[k].size());
        nameW += 2;

        infoMsg("Outlier ratios (IQR=%.2f, bound=%.2f):", outlierIqrRatio, outlierAbsBound);

        std::ostringstream hdr;
        hdr << "    " << std::setw(static_cast<int>(nameW)) << std::left << "";
        for (const auto &tl : tauLabels)
            hdr << std::setw(10) << std::right << tl;
        fprintf(stderr, "%s\n", hdr.str().c_str());

        for (int k = 0; k < K; ++k) {
            std::ostringstream row;
            row << "    " << std::setw(static_cast<int>(nameW)) << std::left << phenoNames[k];
            for (double r : allOutlierRatios[k])
                row << std::setw(10) << std::right << std::fixed << std::setprecision(4) << r;
            fprintf(stderr, "%s\n", row.str().c_str());
        }
    }

    // Free per-phenotype work data
    pw.clear();

    // ── 8. Run multi-phenotype engine ───────────────────────────────
    infoMsg("SPAsqr: Starting multi-phenotype association (%d phenotypes, %d taus, %d threads)", K, ntaus, nthreads);
    multiPhenoEngine(
        *genoData, tasks, outPrefix, "SPAsqr", compression, compressionLevel,
        nthreads, missingCutoff, minMafCutoff, minMacCutoff, hweCutoff
    );
}

// ══════════════════════════════════════════════════════════════════════
// runSPAsqrLoco — LOCO entry point
//
// Matches the serial workflow (2.serial_loco.sh):
//   - conquer fits are chromosome-independent (same Y, same X)
//   - per-phenotype bandwidth h
//   - SPAsqrMethod built once per phenotype, cloned per chromosome
//   - locoEngine iterates chromosomes, testing each chromosome's markers
// ══════════════════════════════════════════════════════════════════════

void runSPAsqrLoco(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<double> &taus,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const GenoSpec &geno,
    const std::string &predListFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierIqrRatio,
    double outlierAbsBound,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    double spasqrTol,
    double spasqrH,
    double spasqrHScale,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const int K = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());

    // ── 1. Load phenotype/covariate data (union mask) ───────────────
    // Union = subjects with genotype ∩ GRM ∩ keep/remove.  Per-phenotype
    // NA filtering is deferred — each phenotype uses its own non-missing
    // subset of the union.
    infoMsg("SPAsqr-LOCO: Loading phenotype and covariate data (%d phenotypes, %d taus)", K, ntaus);
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile, phenoNames);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();

    const uint32_t nUnion = sd.nUsed();
    const Eigen::Index N = static_cast<Eigen::Index>(nUnion);

    // ── 2. Extract union-space covariates ───────────────────────────
    Eigen::MatrixXd unionX = covarNames.empty()
        ? (sd.hasCovar() ? Eigen::MatrixXd(sd.covar()) : Eigen::MatrixXd(N, 0))
        : sd.getColumns(covarNames);
    const int nCov = static_cast<int>(unionX.cols());

    // ── 3. Per-phenotype: build non-missing mask, extract Y/baseX ──
    struct PhenoWork {
        std::vector<uint32_t> unionToLocal; // size nUnion; UINT32_MAX = absent
        uint32_t nk;                        // non-missing count
        Eigen::VectorXd Y;                  // nk
        Eigen::MatrixXd baseX;              // nk × nCov
        double h;                           // bandwidth
    };

    std::vector<PhenoWork> pw(K);

    for (int k = 0; k < K; ++k) {
        Eigen::VectorXd fullY = sd.getColumn(phenoNames[k]);
        pw[k].unionToLocal.resize(nUnion, UINT32_MAX);
        uint32_t localIdx = 0;
        for (uint32_t i = 0; i < nUnion; ++i) {
            if (!std::isnan(fullY[i])) {
                pw[k].unionToLocal[i] = localIdx++;
            }
        }
        pw[k].nk = localIdx;
        if (pw[k].nk == 0)
            throw std::runtime_error("SPAsqr-LOCO: phenotype '" + phenoNames[k] + "' has no non-missing subjects");

        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        pw[k].Y.resize(Nk);
        pw[k].baseX.resize(Nk, nCov);
        for (uint32_t i = 0; i < nUnion; ++i) {
            uint32_t li = pw[k].unionToLocal[i];
            if (li == UINT32_MAX) continue;
            pw[k].Y[li] = fullY[i];
            if (nCov > 0) pw[k].baseX.row(li) = unionX.row(i);
        }

        // Per-phenotype bandwidth
        if (spasqrH >= 0.0) {
            pw[k].h = spasqrH;
        } else {
            std::vector<double> ysorted(Nk);
            Eigen::VectorXd::Map(ysorted.data(), Nk) = pw[k].Y;
            std::sort(ysorted.begin(), ysorted.end());
            auto quantile = [&](double prob) -> double {
                double idx = prob * (Nk - 1);
                Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
                Eigen::Index hi = std::min(lo + 1, Nk - 1);
                double frac = idx - lo;
                return ysorted[lo] * (1.0 - frac) + ysorted[hi] * frac;
            };
            double iqr = quantile(0.75) - quantile(0.25);
            double scale = (spasqrHScale >= 0.0) ? spasqrHScale : 3.0;
            pw[k].h = iqr / scale;
            if (pw[k].h <= 0.0)
                pw[k].h = std::max(std::pow((std::log(Nk) + nCov) / static_cast<double>(Nk), 0.4), 0.05);
        }
    }

    // Log N/bandwidth table
    {
        infoMsg("Sample size and smooth bandwidth per phenotype:");
        size_t nameW = 4;
        for (int k = 0; k < K; ++k)
            nameW = std::max(nameW, phenoNames[k].size());
        nameW += 2;
        char row[256];
        std::snprintf(row, sizeof(row), "    %-*s %10s %12s\n", static_cast<int>(nameW), "", "N", "bandwidth");
        fprintf(stderr, "%s", row);
        for (int k = 0; k < K; ++k) {
            std::snprintf(row, sizeof(row), "    %-*s %10u %12.6f\n",
                          static_cast<int>(nameW), phenoNames[k].c_str(), pw[k].nk, pw[k].h);
            fprintf(stderr, "%s", row);
        }
    }

    // ── 4. Load LOCO predictions (union space) ─────────────────────
    LocoData loco = LocoData::load(predListFile, phenoNames, sd.usedIIDs(), sd.famIIDs());
    auto locoChroms = loco.availableChromosomes();
    infoMsg("LOCO: %zu chromosomes available across all phenotypes", locoChroms.size());

    // ── 5. Build tau labels ─────────────────────────────────────────
    std::vector<std::string> tauLabels;
    tauLabels.reserve(ntaus);
    for (double tau : taus) {
        std::ostringstream oss;
        oss << "tau" << tau;
        tauLabels.push_back(oss.str());
    }

    // ── 6. Load genotype data and GRM once (union space) ────────────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    std::vector<GRMEntry> unionGrm = loadGrmEntries(sd.usedIIDs(), sd.famIIDs(), spgrmGrabFile, spgrmGctaFile);

    // Pre-compute per-phenotype re-indexed GRM (shared across chromosomes)
    std::vector<std::vector<GRMEntry> > phenoGrms(K);
    for (int k = 0; k < K; ++k)
        phenoGrms[k] = reindexGrm(unionGrm, pw[k].unionToLocal, nUnion);

    // ── 7. Build LocoTaskBuilder callback ──────────────────────────
    // For each chromosome, augment base covariates with the LOCO column
    // (mapped to pheno-dense space), run K × ntaus conquer fits, and
    // build K SPAsqrMethods.
    auto buildTasks = [&](const std::string &chr, std::vector<PhenoTask> &tasks) {
        tasks.resize(K);

        // Augment base covariates with LOCO column per phenotype (pheno-dense)
        std::vector<Eigen::MatrixXd> X_augs(K);
        for (int k = 0; k < K; ++k) {
            const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
            X_augs[k].resize(Nk, nCov + 1);
            X_augs[k].leftCols(nCov) = pw[k].baseX;
            // Extract LOCO scores from union to pheno-dense space
            const auto &locoVec = loco.scores.at(phenoNames[k]).at(chr);
            for (uint32_t i = 0; i < nUnion; ++i) {
                uint32_t li = pw[k].unionToLocal[i];
                if (li != UINT32_MAX)
                    X_augs[k](li, nCov) = locoVec[i];
            }
        }

        // Parallel conquer fits: K × ntaus
        std::vector<Eigen::MatrixXd> ResidMats(K);
        for (int k = 0; k < K; ++k)
            ResidMats[k].resize(static_cast<Eigen::Index>(pw[k].nk), ntaus);

        const int totalFits = K * ntaus;
        const int nWorkers = std::min(nthreads, totalFits);
        infoMsg("SPAsqr-LOCO chr%s: Running %d conquer fits with %d threads",
                chr.c_str(), totalFits, nWorkers);

        std::atomic<int> nextFit{0};
        std::vector<std::string> fitErrors(totalFits);

        auto fitWorker = [&]() {
            for (;;) {
                int idx = nextFit.fetch_add(1, std::memory_order_relaxed);
                if (idx >= totalFits) break;
                int k = idx / ntaus;
                int t = idx % ntaus;
                const double h1 = 1.0 / pw[k].h;

                try {
                    Eigen::VectorXd resid;
                    conquer::smqrGauss(X_augs[k], pw[k].Y, taus[t], pw[k].h, &resid, spasqrTol);

                    const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
                    for (Eigen::Index i = 0; i < Nk; ++i)
                        ResidMats[k](i, t) = taus[t] - math::pnorm(-resid(i) * h1);
                } catch (const std::exception &ex) {
                    fitErrors[idx] = ex.what();
                }
            }
        };

        {
            std::vector<std::thread> threads;
            threads.reserve(nWorkers - 1);
            for (int t = 0; t < nWorkers - 1; ++t)
                threads.emplace_back(fitWorker);
            fitWorker();
            for (auto &th : threads)
                th.join();
        }

        for (int idx = 0; idx < totalFits; ++idx) {
            if (!fitErrors[idx].empty()) {
                int k = idx / ntaus;
                int t = idx % ntaus;
                throw std::runtime_error("SPAsqr-LOCO chr" + chr + ": conquer failed for phenotype '" +
                                         phenoNames[k] + "' tau=" + std::to_string(taus[t]) + ": " + fitErrors[idx]);
            }
        }

        // Build SPAsqrMethod for each phenotype (pheno-dense space)
        std::vector<std::vector<double> > allOutlierRatios(K);
        for (int k = 0; k < K; ++k) {
            auto method = buildSPAsqrMethod(
                ResidMats[k],
                phenoGrms[k],
                pw[k].nk,
                spaCutoff,
                outlierIqrRatio,
                outlierAbsBound,
                minMafCutoff,
                minMacCutoff,
                tauLabels,
                &allOutlierRatios[k]
            );

            tasks[k].phenoName = phenoNames[k];
            tasks[k].method = std::move(method);
            tasks[k].unionToLocal = pw[k].unionToLocal;
            tasks[k].nUsed = pw[k].nk;
        }

        // Print outlier ratio table for this chromosome line-by-line
        {
            size_t nameW = 4;
            for (int k = 0; k < K; ++k)
                nameW = std::max(nameW, phenoNames[k].size());
            nameW += 2;

            infoMsg("chr%s outlier ratios (IQR=%.2f, bound=%.2f):",
                    chr.c_str(), outlierIqrRatio, outlierAbsBound);

            std::ostringstream hdr;
            hdr << "    " << std::setw(static_cast<int>(nameW)) << std::left << "";
            for (const auto &tl : tauLabels)
                hdr << std::setw(10) << std::right << tl;
            fprintf(stderr, "%s\n", hdr.str().c_str());

            for (int k = 0; k < K; ++k) {
                std::ostringstream row;
                row << "    " << std::setw(static_cast<int>(nameW)) << std::left << phenoNames[k];
                for (double r : allOutlierRatios[k])
                    row << std::setw(10) << std::right << std::fixed << std::setprecision(4) << r;
                fprintf(stderr, "%s\n", row.str().c_str());
            }
        }
    };

    // ── 8. Run LOCO engine ─────────────────────────────────────────
    const int nChroms = static_cast<int>(locoChroms.size());
    infoMsg("SPAsqr-LOCO: Starting LOCO association (%d phenotypes, %d taus, %d chroms, %d threads)",
            K, ntaus, nChroms, nthreads);
    locoEngine(
        *genoData,
        locoChroms,
        phenoNames,
        buildTasks,
        outPrefix,
        "SPAsqr",
        compression,
        compressionLevel,
        nthreads,
        missingCutoff,
        minMafCutoff,
        minMacCutoff,
        hweCutoff
    );
}
