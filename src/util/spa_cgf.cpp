// spa_cgf.cpp — Cancellation-free binomial CGF kernels: scalar / AVX2 / AVX-512
//
// See src/util/spa_cgf.hpp for the algebra, the K'' derivation, the k12/kFull
// split, and the degenerate-case contract.  This file carries only the loop
// structure and the SIMD dispatch.
//
// A .cpp is required (unlike a header-only kernel) because the
// __attribute__((target(...))) variants and the `static const` dispatch pointers
// must have exactly one definition in the program.
//
// Structure, repeated identically for each of the three variants:
//
//     <variant>K12_scalar        plain loop, compiles on any architecture
//     <variant>K12_avx2          4 x double, scalar tail
//     <variant>K12_avx512        8 x double, masked tail
//     pick<Variant>K12Fn()       resolves via simdLevel()
//     static const ... K12Fn     resolved once, at static-init
//     <variant>K12               public entry point: endpoint hoist + dispatch
//     <variant>KFull             public entry point: K12 dispatch + scalar K
//
// Masked / padded lanes contribute exactly zero: both K' and K'' carry a factor
// r (and, for binomHapcount, a factor h), and inactive lanes load zero.

#include "util/spa_cgf.hpp"

#include <cmath>

#if defined(__x86_64__) || defined(_M_X64)
#  include <immintrin.h>
#endif
#include "util/simd_dispatch.hpp"
#include "util/simd_math.hpp"

namespace spa_cgf {

namespace {

// ══════════════════════════════════════════════════════════════════════
// Shared scalar loop body
// ══════════════════════════════════════════════════════════════════════
//
// The three scalar variants differ only in where p and h come from, so they all
// accumulate through this one step.  `guardAlpha` is a template parameter rather
// than a runtime flag so the compiler drops the compare entirely for the two
// variants whose endpoints are hoisted out of the loop.
template <bool guardAlpha>
inline void accumulateK12(
    double t,
    double r,
    double p,
    double h,
    double &K1,
    double &K2
) {
    const double u  = 1.0 - p;
    const double tr = clampTr(t * r);
    const double e  = p * std::exp(tr);

    double a  = u + e;
    double en = e;
    if (guardAlpha && a == 0.0) { a = 1.0; en = 1.0; }

    const double inva = 1.0 / a;
    const double pi   = en * inva;    // pi   = e/alpha
    const double omp  = u  * inva;    // 1-pi = u/alpha, never a subtraction
    const double hr   = h * r;

    K1 += hr * pi;
    K2 += hr * r * pi * omp;
}

#if defined(__x86_64__) || defined(_M_X64)

// Horizontal sum of an AVX2 accumulator.  Follows src/spasqr/spasqr.cpp.
__attribute__((target("avx2,fma")))
inline double hsum256(__m256d v) {
    __m128d lo = _mm256_castpd256_pd128(v);
    __m128d hi = _mm256_extractf128_pd(v, 1);
    lo = _mm_add_pd(lo, hi);
    return _mm_cvtsd_f64(lo) + _mm_cvtsd_f64(_mm_unpackhi_pd(lo, lo));
}

#endif  // x86_64

}  // namespace

// ══════════════════════════════════════════════════════════════════════
// binomUniform — G_i ~ Binomial(2, p), one p for all subjects
// ══════════════════════════════════════════════════════════════════════
//
// h == 2 throughout, folded into the constant 2.0 rather than carried as a
// vector.  The public entry point resolves p == 0 and p == 1 before dispatching,
// so these loops may assume 0 < p < 1, hence u > 0 and alpha >= u > 0: no
// per-lane alpha == 0 guard is compiled in.

namespace tier {

Cgf12 binomUniformK12_scalar(double t, const double *resid, int n, double p) noexcept {
    double K1 = 0.0, K2 = 0.0;
    for (int i = 0; i < n; ++i)
        accumulateK12<false>(t, resid[i], p, 2.0, K1, K2);
    return Cgf12{K1, K2};
}

#if defined(__x86_64__) || defined(_M_X64)

__attribute__((target("avx2,avx512f,avx512vl,fma")))
Cgf12 binomUniformK12_avx512(double t, const double *resid, int n, double p) noexcept {
    const __m512d vt    = _mm512_set1_pd(t);
    const __m512d vp    = _mm512_set1_pd(p);
    const __m512d vu    = _mm512_set1_pd(1.0 - p);
    const __m512d vTwo  = _mm512_set1_pd(2.0);
    const __m512d vHi   = _mm512_set1_pd(kTrClamp);
    const __m512d vLo   = _mm512_set1_pd(-kTrClamp);

    __m512d vK1 = _mm512_setzero_pd();
    __m512d vK2 = _mm512_setzero_pd();

    int i = 0;
    const int n8 = n & ~7;
    for (; i < n8; i += 8) {
        const __m512d vr  = _mm512_loadu_pd(resid + i);
        const __m512d vtr = _mm512_max_pd(_mm512_min_pd(_mm512_mul_pd(vt, vr), vHi), vLo);
        const __m512d ve  = _mm512_mul_pd(vp, avx512_exp_pd(vtr));
        const __m512d va  = _mm512_add_pd(vu, ve);

        const __m512d vinva = _mm512_div_pd(_mm512_set1_pd(1.0), va);
        const __m512d vpi   = _mm512_mul_pd(ve, vinva);
        const __m512d vomp  = _mm512_mul_pd(vu, vinva);
        const __m512d vhr   = _mm512_mul_pd(vTwo, vr);

        vK1 = _mm512_add_pd(vK1, _mm512_mul_pd(vhr, vpi));
        vK2 = _mm512_add_pd(vK2,
                            _mm512_mul_pd(_mm512_mul_pd(vhr, vr),
                                          _mm512_mul_pd(vpi, vomp)));
    }

    // Masked tail: inactive lanes load r = 0, so tr = 0, pi = p, omp = u, and
    // both accumulated terms carry the factor r == 0 — exactly zero contribution.
    if (i < n) {
        const __mmask8 mask = static_cast<__mmask8>((1u << (n - i)) - 1u);
        const __m512d vr  = _mm512_maskz_loadu_pd(mask, resid + i);
        const __m512d vtr = _mm512_max_pd(_mm512_min_pd(_mm512_mul_pd(vt, vr), vHi), vLo);
        const __m512d ve  = _mm512_mul_pd(vp, avx512_exp_pd(vtr));
        const __m512d va  = _mm512_add_pd(vu, ve);

        const __m512d vinva = _mm512_div_pd(_mm512_set1_pd(1.0), va);
        const __m512d vpi   = _mm512_mul_pd(ve, vinva);
        const __m512d vomp  = _mm512_mul_pd(vu, vinva);
        const __m512d vhr   = _mm512_mul_pd(vTwo, vr);

        vK1 = _mm512_add_pd(vK1, _mm512_mul_pd(vhr, vpi));
        vK2 = _mm512_add_pd(vK2,
                            _mm512_mul_pd(_mm512_mul_pd(vhr, vr),
                                          _mm512_mul_pd(vpi, vomp)));
    }

    return Cgf12{_mm512_reduce_add_pd(vK1), _mm512_reduce_add_pd(vK2)};
}

__attribute__((target("avx2,fma")))
Cgf12 binomUniformK12_avx2(double t, const double *resid, int n, double p) noexcept {
    const __m256d vt   = _mm256_set1_pd(t);
    const __m256d vp   = _mm256_set1_pd(p);
    const __m256d vu   = _mm256_set1_pd(1.0 - p);
    const __m256d vTwo = _mm256_set1_pd(2.0);
    const __m256d vOne = _mm256_set1_pd(1.0);
    const __m256d vHi  = _mm256_set1_pd(kTrClamp);
    const __m256d vLo  = _mm256_set1_pd(-kTrClamp);

    __m256d vK1 = _mm256_setzero_pd();
    __m256d vK2 = _mm256_setzero_pd();

    int i = 0;
    const int n4 = n & ~3;
    for (; i < n4; i += 4) {
        const __m256d vr  = _mm256_loadu_pd(resid + i);
        const __m256d vtr = _mm256_max_pd(_mm256_min_pd(_mm256_mul_pd(vt, vr), vHi), vLo);
        const __m256d ve  = _mm256_mul_pd(vp, avx2_exp_pd(vtr));
        const __m256d va  = _mm256_add_pd(vu, ve);

        const __m256d vinva = _mm256_div_pd(vOne, va);
        const __m256d vpi   = _mm256_mul_pd(ve, vinva);
        const __m256d vomp  = _mm256_mul_pd(vu, vinva);
        const __m256d vhr   = _mm256_mul_pd(vTwo, vr);

        vK1 = _mm256_add_pd(vK1, _mm256_mul_pd(vhr, vpi));
        vK2 = _mm256_add_pd(vK2,
                            _mm256_mul_pd(_mm256_mul_pd(vhr, vr),
                                          _mm256_mul_pd(vpi, vomp)));
    }

    double K1 = hsum256(vK1);
    double K2 = hsum256(vK2);

    // Scalar tail (1-3 subjects); AVX2 has no masked load.
    for (; i < n; ++i)
        accumulateK12<false>(t, resid[i], p, 2.0, K1, K2);

    return Cgf12{K1, K2};
}

#endif  // x86_64

}  // namespace tier

namespace {

using UniformK12Fn = Cgf12 (*)(double, const double *, int, double) noexcept;

UniformK12Fn pickUniformK12Fn() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return tier::binomUniformK12_avx512;
    case SimdLevel::AVX2:   return tier::binomUniformK12_avx2;
    default: break;
    }
#endif
    return tier::binomUniformK12_scalar;
}

const UniformK12Fn uniformK12Fn = pickUniformK12Fn();

}  // namespace

Cgf12 binomUniformK12(double t, const double *resid, int n, double p) noexcept {
    if (n <= 0) return Cgf12{0.0, 0.0};

    // Endpoint hoist.  p == 1: G_i == 2 almost surely, so 2·r_i is deterministic
    // and K'_i = 2·r_i, K''_i = 0.  p == 0: G_i == 0 almost surely, all cumulants
    // vanish.  Resolving both here lets the loops omit the alpha == 0 guard.
    if (p >= 1.0) {
        double sumR = 0.0;
        for (int i = 0; i < n; ++i) sumR += resid[i];
        return Cgf12{2.0 * sumR, 0.0};
    }
    if (p <= 0.0) return Cgf12{0.0, 0.0};

    return uniformK12Fn(t, resid, n, p);
}

Cgf012 binomUniformKFull(double t, const double *resid, int n, double p) noexcept {
    const Cgf12 d = binomUniformK12(t, resid, n, p);

    // K accumulates scalar: expm1/log1p have no vectorized form and, per N1, do
    // not need one — this runs once per tail, not once per Newton iteration.
    // Multiplication by 2.0 is exact, so factoring it out of the sum is
    // bit-identical to folding it into each term.
    double sumLogAlpha = 0.0;
    for (int i = 0; i < n; ++i) sumLogAlpha += logAlpha(t * resid[i], p);

    return Cgf012{2.0 * sumLogAlpha, d.K1, d.K2};
}

// ══════════════════════════════════════════════════════════════════════
// binomIndiv — G_i ~ Binomial(2, af[i]), per-individual allele frequency
// ══════════════════════════════════════════════════════════════════════
//
// af varies per lane, so u = 1 - af[i] is a vector and the endpoints cannot be
// hoisted.  af[i] == 1 combined with an underflowed lambda makes alpha vanish,
// so the alpha == 0 guard is applied per lane: one compare and two blends per
// vector, negligible beside the exp and the divide.

namespace tier {

Cgf12 binomIndivK12_scalar(
    double t,
    const double *resid,
    const double *af,
    int n
) noexcept {
    double K1 = 0.0, K2 = 0.0;
    for (int i = 0; i < n; ++i)
        accumulateK12<true>(t, resid[i], af[i], 2.0, K1, K2);
    return Cgf12{K1, K2};
}

#if defined(__x86_64__) || defined(_M_X64)

__attribute__((target("avx2,avx512f,avx512vl,fma")))
Cgf12 binomIndivK12_avx512(
    double t,
    const double *resid,
    const double *af,
    int n
) noexcept {
    const __m512d vt    = _mm512_set1_pd(t);
    const __m512d vOne  = _mm512_set1_pd(1.0);
    const __m512d vTwo  = _mm512_set1_pd(2.0);
    const __m512d vZero = _mm512_setzero_pd();
    const __m512d vHi   = _mm512_set1_pd(kTrClamp);
    const __m512d vLo   = _mm512_set1_pd(-kTrClamp);

    __m512d vK1 = _mm512_setzero_pd();
    __m512d vK2 = _mm512_setzero_pd();

    int i = 0;
    const int n8 = n & ~7;
    for (; i < n8; i += 8) {
        const __m512d vr = _mm512_loadu_pd(resid + i);
        const __m512d vp = _mm512_loadu_pd(af + i);
        const __m512d vu = _mm512_sub_pd(vOne, vp);

        const __m512d vtr = _mm512_max_pd(_mm512_min_pd(_mm512_mul_pd(vt, vr), vHi), vLo);
        __m512d ve = _mm512_mul_pd(vp, avx512_exp_pd(vtr));
        __m512d va = _mm512_add_pd(vu, ve);

        // alpha == 0 requires p == 1 (u == 0) and an underflowed lambda; the law
        // is then G == 2 a.s., so pi = 1 and pi·(1-pi) = 0.  (e, alpha) <- (1, 1)
        // with u == 0 reproduces both exactly.
        const __mmask8 deg = _mm512_cmp_pd_mask(va, vZero, _CMP_EQ_OQ);
        ve = _mm512_mask_blend_pd(deg, ve, vOne);
        va = _mm512_mask_blend_pd(deg, va, vOne);

        const __m512d vinva = _mm512_div_pd(vOne, va);
        const __m512d vpi   = _mm512_mul_pd(ve, vinva);
        const __m512d vomp  = _mm512_mul_pd(vu, vinva);
        const __m512d vhr   = _mm512_mul_pd(vTwo, vr);

        vK1 = _mm512_add_pd(vK1, _mm512_mul_pd(vhr, vpi));
        vK2 = _mm512_add_pd(vK2,
                            _mm512_mul_pd(_mm512_mul_pd(vhr, vr),
                                          _mm512_mul_pd(vpi, vomp)));
    }

    // Masked tail: r = 0 and af = 0 in inactive lanes, so alpha = 1 and both
    // accumulated terms carry the factor r == 0.
    if (i < n) {
        const __mmask8 mask = static_cast<__mmask8>((1u << (n - i)) - 1u);
        const __m512d vr = _mm512_maskz_loadu_pd(mask, resid + i);
        const __m512d vp = _mm512_maskz_loadu_pd(mask, af + i);
        const __m512d vu = _mm512_sub_pd(vOne, vp);

        const __m512d vtr = _mm512_max_pd(_mm512_min_pd(_mm512_mul_pd(vt, vr), vHi), vLo);
        __m512d ve = _mm512_mul_pd(vp, avx512_exp_pd(vtr));
        __m512d va = _mm512_add_pd(vu, ve);

        const __mmask8 deg = _mm512_cmp_pd_mask(va, vZero, _CMP_EQ_OQ);
        ve = _mm512_mask_blend_pd(deg, ve, vOne);
        va = _mm512_mask_blend_pd(deg, va, vOne);

        const __m512d vinva = _mm512_div_pd(vOne, va);
        const __m512d vpi   = _mm512_mul_pd(ve, vinva);
        const __m512d vomp  = _mm512_mul_pd(vu, vinva);
        const __m512d vhr   = _mm512_mul_pd(vTwo, vr);

        vK1 = _mm512_add_pd(vK1, _mm512_mul_pd(vhr, vpi));
        vK2 = _mm512_add_pd(vK2,
                            _mm512_mul_pd(_mm512_mul_pd(vhr, vr),
                                          _mm512_mul_pd(vpi, vomp)));
    }

    return Cgf12{_mm512_reduce_add_pd(vK1), _mm512_reduce_add_pd(vK2)};
}

__attribute__((target("avx2,fma")))
Cgf12 binomIndivK12_avx2(
    double t,
    const double *resid,
    const double *af,
    int n
) noexcept {
    const __m256d vt    = _mm256_set1_pd(t);
    const __m256d vOne  = _mm256_set1_pd(1.0);
    const __m256d vTwo  = _mm256_set1_pd(2.0);
    const __m256d vZero = _mm256_setzero_pd();
    const __m256d vHi   = _mm256_set1_pd(kTrClamp);
    const __m256d vLo   = _mm256_set1_pd(-kTrClamp);

    __m256d vK1 = _mm256_setzero_pd();
    __m256d vK2 = _mm256_setzero_pd();

    int i = 0;
    const int n4 = n & ~3;
    for (; i < n4; i += 4) {
        const __m256d vr = _mm256_loadu_pd(resid + i);
        const __m256d vp = _mm256_loadu_pd(af + i);
        const __m256d vu = _mm256_sub_pd(vOne, vp);

        const __m256d vtr = _mm256_max_pd(_mm256_min_pd(_mm256_mul_pd(vt, vr), vHi), vLo);
        __m256d ve = _mm256_mul_pd(vp, avx2_exp_pd(vtr));
        __m256d va = _mm256_add_pd(vu, ve);

        const __m256d deg = _mm256_cmp_pd(va, vZero, _CMP_EQ_OQ);
        ve = _mm256_blendv_pd(ve, vOne, deg);
        va = _mm256_blendv_pd(va, vOne, deg);

        const __m256d vinva = _mm256_div_pd(vOne, va);
        const __m256d vpi   = _mm256_mul_pd(ve, vinva);
        const __m256d vomp  = _mm256_mul_pd(vu, vinva);
        const __m256d vhr   = _mm256_mul_pd(vTwo, vr);

        vK1 = _mm256_add_pd(vK1, _mm256_mul_pd(vhr, vpi));
        vK2 = _mm256_add_pd(vK2,
                            _mm256_mul_pd(_mm256_mul_pd(vhr, vr),
                                          _mm256_mul_pd(vpi, vomp)));
    }

    double K1 = hsum256(vK1);
    double K2 = hsum256(vK2);

    for (; i < n; ++i)
        accumulateK12<true>(t, resid[i], af[i], 2.0, K1, K2);

    return Cgf12{K1, K2};
}

#endif  // x86_64

}  // namespace tier

namespace {

using IndivK12Fn = Cgf12 (*)(double, const double *, const double *, int) noexcept;

IndivK12Fn pickIndivK12Fn() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return tier::binomIndivK12_avx512;
    case SimdLevel::AVX2:   return tier::binomIndivK12_avx2;
    default: break;
    }
#endif
    return tier::binomIndivK12_scalar;
}

const IndivK12Fn indivK12Fn = pickIndivK12Fn();

}  // namespace

Cgf12 binomIndivK12(double t, const double *resid, const double *af, int n) noexcept {
    if (n <= 0) return Cgf12{0.0, 0.0};
    return indivK12Fn(t, resid, af, n);
}

Cgf012 binomIndivKFull(
    double t,
    const double *resid,
    const double *af,
    int n
) noexcept {
    const Cgf12 d = binomIndivK12(t, resid, af, n);

    double sumLogAlpha = 0.0;
    for (int i = 0; i < n; ++i) sumLogAlpha += logAlpha(t * resid[i], af[i]);

    return Cgf012{2.0 * sumLogAlpha, d.K1, d.K2};
}

// ══════════════════════════════════════════════════════════════════════
// binomHapcount — G_i ~ Binomial(hap[i], q), variable trial count
// ══════════════════════════════════════════════════════════════════════
//
// q is shared, so its endpoints hoist exactly as in binomUniform and the loops
// may assume 0 < q < 1.  hap varies per lane and is carried as a vector; it
// needs no guard of its own, every cumulant being proportional to it — hap == 0
// yields exactly zero rather than a special case.

namespace tier {

Cgf12 binomHapcountK12_scalar(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q
) noexcept {
    double K1 = 0.0, K2 = 0.0;
    for (int i = 0; i < n; ++i)
        accumulateK12<false>(t, resid[i], q, hap[i], K1, K2);
    return Cgf12{K1, K2};
}

#if defined(__x86_64__) || defined(_M_X64)

__attribute__((target("avx2,avx512f,avx512vl,fma")))
Cgf12 binomHapcountK12_avx512(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q
) noexcept {
    const __m512d vt   = _mm512_set1_pd(t);
    const __m512d vq   = _mm512_set1_pd(q);
    const __m512d vu   = _mm512_set1_pd(1.0 - q);
    const __m512d vOne = _mm512_set1_pd(1.0);
    const __m512d vHi  = _mm512_set1_pd(kTrClamp);
    const __m512d vLo  = _mm512_set1_pd(-kTrClamp);

    __m512d vK1 = _mm512_setzero_pd();
    __m512d vK2 = _mm512_setzero_pd();

    int i = 0;
    const int n8 = n & ~7;
    for (; i < n8; i += 8) {
        const __m512d vr = _mm512_loadu_pd(resid + i);
        const __m512d vh = _mm512_loadu_pd(hap + i);

        const __m512d vtr = _mm512_max_pd(_mm512_min_pd(_mm512_mul_pd(vt, vr), vHi), vLo);
        const __m512d ve  = _mm512_mul_pd(vq, avx512_exp_pd(vtr));
        const __m512d va  = _mm512_add_pd(vu, ve);

        const __m512d vinva = _mm512_div_pd(vOne, va);
        const __m512d vpi   = _mm512_mul_pd(ve, vinva);
        const __m512d vomp  = _mm512_mul_pd(vu, vinva);
        const __m512d vhr   = _mm512_mul_pd(vh, vr);

        vK1 = _mm512_add_pd(vK1, _mm512_mul_pd(vhr, vpi));
        vK2 = _mm512_add_pd(vK2,
                            _mm512_mul_pd(_mm512_mul_pd(vhr, vr),
                                          _mm512_mul_pd(vpi, vomp)));
    }

    // Masked tail: r = 0 and hap = 0 in inactive lanes — zero twice over.
    if (i < n) {
        const __mmask8 mask = static_cast<__mmask8>((1u << (n - i)) - 1u);
        const __m512d vr = _mm512_maskz_loadu_pd(mask, resid + i);
        const __m512d vh = _mm512_maskz_loadu_pd(mask, hap + i);

        const __m512d vtr = _mm512_max_pd(_mm512_min_pd(_mm512_mul_pd(vt, vr), vHi), vLo);
        const __m512d ve  = _mm512_mul_pd(vq, avx512_exp_pd(vtr));
        const __m512d va  = _mm512_add_pd(vu, ve);

        const __m512d vinva = _mm512_div_pd(vOne, va);
        const __m512d vpi   = _mm512_mul_pd(ve, vinva);
        const __m512d vomp  = _mm512_mul_pd(vu, vinva);
        const __m512d vhr   = _mm512_mul_pd(vh, vr);

        vK1 = _mm512_add_pd(vK1, _mm512_mul_pd(vhr, vpi));
        vK2 = _mm512_add_pd(vK2,
                            _mm512_mul_pd(_mm512_mul_pd(vhr, vr),
                                          _mm512_mul_pd(vpi, vomp)));
    }

    return Cgf12{_mm512_reduce_add_pd(vK1), _mm512_reduce_add_pd(vK2)};
}

__attribute__((target("avx2,fma")))
Cgf12 binomHapcountK12_avx2(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q
) noexcept {
    const __m256d vt   = _mm256_set1_pd(t);
    const __m256d vq   = _mm256_set1_pd(q);
    const __m256d vu   = _mm256_set1_pd(1.0 - q);
    const __m256d vOne = _mm256_set1_pd(1.0);
    const __m256d vHi  = _mm256_set1_pd(kTrClamp);
    const __m256d vLo  = _mm256_set1_pd(-kTrClamp);

    __m256d vK1 = _mm256_setzero_pd();
    __m256d vK2 = _mm256_setzero_pd();

    int i = 0;
    const int n4 = n & ~3;
    for (; i < n4; i += 4) {
        const __m256d vr = _mm256_loadu_pd(resid + i);
        const __m256d vh = _mm256_loadu_pd(hap + i);

        const __m256d vtr = _mm256_max_pd(_mm256_min_pd(_mm256_mul_pd(vt, vr), vHi), vLo);
        const __m256d ve  = _mm256_mul_pd(vq, avx2_exp_pd(vtr));
        const __m256d va  = _mm256_add_pd(vu, ve);

        const __m256d vinva = _mm256_div_pd(vOne, va);
        const __m256d vpi   = _mm256_mul_pd(ve, vinva);
        const __m256d vomp  = _mm256_mul_pd(vu, vinva);
        const __m256d vhr   = _mm256_mul_pd(vh, vr);

        vK1 = _mm256_add_pd(vK1, _mm256_mul_pd(vhr, vpi));
        vK2 = _mm256_add_pd(vK2,
                            _mm256_mul_pd(_mm256_mul_pd(vhr, vr),
                                          _mm256_mul_pd(vpi, vomp)));
    }

    double K1 = hsum256(vK1);
    double K2 = hsum256(vK2);

    for (; i < n; ++i)
        accumulateK12<false>(t, resid[i], q, hap[i], K1, K2);

    return Cgf12{K1, K2};
}

#endif  // x86_64

}  // namespace tier

namespace {

using HapcountK12Fn =
    Cgf12 (*)(double, const double *, const double *, int, double) noexcept;

HapcountK12Fn pickHapcountK12Fn() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return tier::binomHapcountK12_avx512;
    case SimdLevel::AVX2:   return tier::binomHapcountK12_avx2;
    default: break;
    }
#endif
    return tier::binomHapcountK12_scalar;
}

const HapcountK12Fn hapcountK12Fn = pickHapcountK12Fn();

}  // namespace

Cgf12 binomHapcountK12(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q
) noexcept {
    if (n <= 0) return Cgf12{0.0, 0.0};

    // Endpoint hoist.  q == 1: G_i == hap[i] almost surely, so K'_i = h_i·r_i and
    // K''_i = 0.  q == 0: G_i == 0 almost surely.  Both are exact, and resolving
    // them here is what lets the loops omit the alpha == 0 guard.
    if (q >= 1.0) {
        double sumHR = 0.0;
        for (int i = 0; i < n; ++i) sumHR += hap[i] * resid[i];
        return Cgf12{sumHR, 0.0};
    }
    if (q <= 0.0) return Cgf12{0.0, 0.0};

    return hapcountK12Fn(t, resid, hap, n, q);
}

Cgf012 binomHapcountKFull(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q
) noexcept {
    const Cgf12 d = binomHapcountK12(t, resid, hap, n, q);

    // hap varies per subject, so h cannot be factored out of the sum here.
    double K0 = 0.0;
    for (int i = 0; i < n; ++i) K0 += hap[i] * logAlpha(t * resid[i], q);

    return Cgf012{K0, d.K1, d.K2};
}

// ══════════════════════════════════════════════════════════════════════
// Tier reporting
// ══════════════════════════════════════════════════════════════════════

int simdLevelValue() {
    return static_cast<int>(simdLevel());
}

}  // namespace spa_cgf
