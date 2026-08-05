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
//     <variant>K0_scalar         the terminal K, same three tiers
//     <variant>K0_avx2
//     <variant>K0_avx512
//     pick<Variant>K12Fn()       resolves via simdLevel()
//     pick<Variant>K0Fn()
//     static const ... Fn        resolved once, at static-init
//     <variant>K12               public entry point: endpoint hoist + dispatch
//     <variant>KFull             public entry point: K12 dispatch + K0 dispatch
//     <variant>KFullExact        public entry point: K12 dispatch + scalar
//                                log1p/expm1 K, for tests and for any consumer
//                                that needs K relatively accurate
//
// Masked / padded lanes contribute exactly zero in both halves: K' and K'' carry
// a factor r (and, for binomHapcount, a factor h), and in the K half a lane with
// r = 0 has alpha = u + p = 1 exactly and log(1) = 0 exactly.  See the SIMD
// section of spa_cgf.hpp.

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

// ── The vector K body, once per width ──
//
// Returns log(alpha) for a vector of arguments, with the past-clamp linear
// excess added back.  `raw` is the unclamped t*r, `vu` and `vp` the (possibly
// broadcast) 1-p and p.  See spa_cgf.hpp's `logAlphaFast` for the three
// properties this relies on.
//
// `guardEndpoints` covers the per-lane p == 0 and p == 1 cases that only
// binomIndiv can present; for the two variants whose p is a hoisted scalar in
// (0, 1) the compiler drops both masks.

__attribute__((target("avx2,avx512f,avx512vl,fma")))
inline __m512d logAlphaVec512(__m512d raw, __m512d vu, __m512d vp,
                              bool guardEndpoints) {
    const __m512d vHi = _mm512_set1_pd(kTrClamp);
    const __m512d vLo = _mm512_set1_pd(-kTrClamp);
    const __m512d vZ  = _mm512_setzero_pd();

    const __m512d tr = _mm512_max_pd(_mm512_min_pd(raw, vHi), vLo);
    const __m512d a  = _mm512_add_pd(vu, _mm512_mul_pd(vp, avx512_exp_pd(tr)));

    const __m512d excess = _mm512_sub_pd(raw, tr);
    __m512d add = _mm512_max_pd(excess, vZ);
    if (guardEndpoints) {
        // p == 0: alpha is identically 1, so no excess applies at any x.
        add = _mm512_mask_blend_pd(_mm512_cmp_pd_mask(vp, vZ, _CMP_EQ_OQ), add, vZ);
        // u == 0 (p == 1): log(alpha) is x itself, so the excess applies on both
        // sides of the clamp, not only above it.
        add = _mm512_mask_blend_pd(_mm512_cmp_pd_mask(vu, vZ, _CMP_EQ_OQ),
                                   add, excess);
    }
    return _mm512_add_pd(avx512_log_pd(a), add);
}

__attribute__((target("avx2,fma")))
inline __m256d logAlphaVec256(__m256d raw, __m256d vu, __m256d vp,
                              bool guardEndpoints) {
    const __m256d vHi = _mm256_set1_pd(kTrClamp);
    const __m256d vLo = _mm256_set1_pd(-kTrClamp);
    const __m256d vZ  = _mm256_setzero_pd();

    const __m256d tr = _mm256_max_pd(_mm256_min_pd(raw, vHi), vLo);
    const __m256d a  = _mm256_add_pd(vu, _mm256_mul_pd(vp, avx2_exp_pd(tr)));

    const __m256d excess = _mm256_sub_pd(raw, tr);
    __m256d add = _mm256_max_pd(excess, vZ);
    if (guardEndpoints) {
        add = _mm256_blendv_pd(add, vZ, _mm256_cmp_pd(vp, vZ, _CMP_EQ_OQ));
        add = _mm256_blendv_pd(add, excess, _mm256_cmp_pd(vu, vZ, _CMP_EQ_OQ));
    }
    return _mm256_add_pd(avx2_log_pd(a), add);
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

// ── the terminal K, sum_i log(alpha_i); the factor h == 2 is the caller's ──

double binomUniformK0_scalar(double t, const double *resid, int n, double p) noexcept {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += logAlphaFast(t * resid[i], p);
    return s;
}

#if defined(__x86_64__) || defined(_M_X64)

__attribute__((target("avx2,avx512f,avx512vl,fma")))
double binomUniformK0_avx512(double t, const double *resid, int n, double p) noexcept {
    const __m512d vt = _mm512_set1_pd(t);
    const __m512d vp = _mm512_set1_pd(p);
    const __m512d vu = _mm512_set1_pd(1.0 - p);

    __m512d acc = _mm512_setzero_pd();

    int i = 0;
    const int n8 = n & ~7;
    for (; i < n8; i += 8) {
        const __m512d vr = _mm512_loadu_pd(resid + i);
        acc = _mm512_add_pd(
            acc, logAlphaVec512(_mm512_mul_pd(vt, vr), vu, vp, /*guard=*/false));
    }

    // Masked tail: inactive lanes load r = 0, so alpha = u + p = 1 exactly and
    // log(1) = 0 exactly — the lane contributes nothing.
    if (i < n) {
        const __mmask8 mask = static_cast<__mmask8>((1u << (n - i)) - 1u);
        const __m512d vr = _mm512_maskz_loadu_pd(mask, resid + i);
        acc = _mm512_add_pd(
            acc, logAlphaVec512(_mm512_mul_pd(vt, vr), vu, vp, /*guard=*/false));
    }

    return _mm512_reduce_add_pd(acc);
}

__attribute__((target("avx2,fma")))
double binomUniformK0_avx2(double t, const double *resid, int n, double p) noexcept {
    const __m256d vt = _mm256_set1_pd(t);
    const __m256d vp = _mm256_set1_pd(p);
    const __m256d vu = _mm256_set1_pd(1.0 - p);

    __m256d acc = _mm256_setzero_pd();

    int i = 0;
    const int n4 = n & ~3;
    for (; i < n4; i += 4) {
        const __m256d vr = _mm256_loadu_pd(resid + i);
        acc = _mm256_add_pd(
            acc, logAlphaVec256(_mm256_mul_pd(vt, vr), vu, vp, /*guard=*/false));
    }

    double s = hsum256(acc);
    for (; i < n; ++i) s += logAlphaFast(t * resid[i], p);
    return s;
}

#endif  // x86_64

}  // namespace tier

namespace {

using UniformK12Fn = Cgf12 (*)(double, const double *, int, double) noexcept;
using UniformK0Fn  = double (*)(double, const double *, int, double) noexcept;

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

UniformK0Fn pickUniformK0Fn() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return tier::binomUniformK0_avx512;
    case SimdLevel::AVX2:   return tier::binomUniformK0_avx2;
    default: break;
    }
#endif
    return tier::binomUniformK0_scalar;
}

const UniformK12Fn uniformK12Fn = pickUniformK12Fn();
const UniformK0Fn  uniformK0Fn  = pickUniformK0Fn();

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

// The p endpoints, resolved in closed form so that the K reduction may assume
// 0 < p < 1 exactly as the k12 reduction does.  p == 1: G_i == 2 a.s., so
// r_i*G_i is deterministic and K = 2*t*sum(r_i).  p == 0: G_i == 0 a.s. and every
// cumulant vanishes.  Both are exact, and both agree with `logAlpha`.
namespace {

bool uniformKEndpoint(double t, const double *resid, int n, double p, double *k0) {
    if (p >= 1.0) {
        double sumR = 0.0;
        for (int i = 0; i < n; ++i) sumR += resid[i];
        *k0 = 2.0 * t * sumR;
        return true;
    }
    if (p <= 0.0) { *k0 = 0.0; return true; }
    return false;
}

}  // namespace

Cgf012 binomUniformKFull(double t, const double *resid, int n, double p) noexcept {
    const Cgf12 d = binomUniformK12(t, resid, n, p);
    if (n <= 0) return Cgf012{0.0, d.K1, d.K2};

    double k0;
    if (uniformKEndpoint(t, resid, n, p, &k0)) return Cgf012{k0, d.K1, d.K2};

    // Multiplication by 2.0 is exact, so factoring it out of the sum is
    // bit-identical to folding it into each term.
    return Cgf012{2.0 * uniformK0Fn(t, resid, n, p), d.K1, d.K2};
}

Cgf012 binomUniformKFullExact(
    double t,
    const double *resid,
    int n,
    double p
) noexcept {
    const Cgf12 d = binomUniformK12(t, resid, n, p);

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

// ── the terminal K.  af varies per lane, so the p == 0 and p == 1 endpoints
// ── cannot be hoisted and the vector body carries both masks.

double binomIndivK0_scalar(
    double t,
    const double *resid,
    const double *af,
    int n
) noexcept {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += logAlphaFast(t * resid[i], af[i]);
    return s;
}

#if defined(__x86_64__) || defined(_M_X64)

__attribute__((target("avx2,avx512f,avx512vl,fma")))
double binomIndivK0_avx512(
    double t,
    const double *resid,
    const double *af,
    int n
) noexcept {
    const __m512d vt   = _mm512_set1_pd(t);
    const __m512d vOne = _mm512_set1_pd(1.0);

    __m512d acc = _mm512_setzero_pd();

    int i = 0;
    const int n8 = n & ~7;
    for (; i < n8; i += 8) {
        const __m512d vr = _mm512_loadu_pd(resid + i);
        const __m512d vp = _mm512_loadu_pd(af + i);
        acc = _mm512_add_pd(acc,
                            logAlphaVec512(_mm512_mul_pd(vt, vr),
                                           _mm512_sub_pd(vOne, vp), vp,
                                           /*guard=*/true));
    }

    // Masked tail: r = 0 and af = 0, so alpha = 1 - 0 + 0 = 1, log(1) = 0, and
    // the excess is zero because the argument is zero.  The p == 0 mask makes
    // that hold for any argument, not only this one.
    if (i < n) {
        const __mmask8 mask = static_cast<__mmask8>((1u << (n - i)) - 1u);
        const __m512d vr = _mm512_maskz_loadu_pd(mask, resid + i);
        const __m512d vp = _mm512_maskz_loadu_pd(mask, af + i);
        acc = _mm512_add_pd(acc,
                            logAlphaVec512(_mm512_mul_pd(vt, vr),
                                           _mm512_sub_pd(vOne, vp), vp,
                                           /*guard=*/true));
    }

    return _mm512_reduce_add_pd(acc);
}

__attribute__((target("avx2,fma")))
double binomIndivK0_avx2(
    double t,
    const double *resid,
    const double *af,
    int n
) noexcept {
    const __m256d vt   = _mm256_set1_pd(t);
    const __m256d vOne = _mm256_set1_pd(1.0);

    __m256d acc = _mm256_setzero_pd();

    int i = 0;
    const int n4 = n & ~3;
    for (; i < n4; i += 4) {
        const __m256d vr = _mm256_loadu_pd(resid + i);
        const __m256d vp = _mm256_loadu_pd(af + i);
        acc = _mm256_add_pd(acc,
                            logAlphaVec256(_mm256_mul_pd(vt, vr),
                                           _mm256_sub_pd(vOne, vp), vp,
                                           /*guard=*/true));
    }

    double s = hsum256(acc);
    for (; i < n; ++i) s += logAlphaFast(t * resid[i], af[i]);
    return s;
}

#endif  // x86_64

}  // namespace tier

namespace {

using IndivK12Fn = Cgf12 (*)(double, const double *, const double *, int) noexcept;
using IndivK0Fn  = double (*)(double, const double *, const double *, int) noexcept;

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

IndivK0Fn pickIndivK0Fn() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return tier::binomIndivK0_avx512;
    case SimdLevel::AVX2:   return tier::binomIndivK0_avx2;
    default: break;
    }
#endif
    return tier::binomIndivK0_scalar;
}

const IndivK12Fn indivK12Fn = pickIndivK12Fn();
const IndivK0Fn  indivK0Fn  = pickIndivK0Fn();

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
    if (n <= 0) return Cgf012{0.0, d.K1, d.K2};
    return Cgf012{2.0 * indivK0Fn(t, resid, af, n), d.K1, d.K2};
}

Cgf012 binomIndivKFullExact(
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

// ── the terminal K, sum_i hap_i * log(alpha_i).  hap varies per subject, so
// ── unlike the two diploid variants the weight cannot be factored out.

double binomHapcountK0_scalar(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q
) noexcept {
    double s = 0.0;
    for (int i = 0; i < n; ++i) s += hap[i] * logAlphaFast(t * resid[i], q);
    return s;
}

#if defined(__x86_64__) || defined(_M_X64)

__attribute__((target("avx2,avx512f,avx512vl,fma")))
double binomHapcountK0_avx512(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q
) noexcept {
    const __m512d vt = _mm512_set1_pd(t);
    const __m512d vq = _mm512_set1_pd(q);
    const __m512d vu = _mm512_set1_pd(1.0 - q);

    __m512d acc = _mm512_setzero_pd();

    int i = 0;
    const int n8 = n & ~7;
    for (; i < n8; i += 8) {
        const __m512d vr = _mm512_loadu_pd(resid + i);
        const __m512d vh = _mm512_loadu_pd(hap + i);
        acc = _mm512_add_pd(
            acc, _mm512_mul_pd(vh, logAlphaVec512(_mm512_mul_pd(vt, vr), vu, vq,
                                                  /*guard=*/false)));
    }

    // Masked tail: r = 0 gives log(alpha) = 0 and hap = 0 zeroes the lane again.
    if (i < n) {
        const __mmask8 mask = static_cast<__mmask8>((1u << (n - i)) - 1u);
        const __m512d vr = _mm512_maskz_loadu_pd(mask, resid + i);
        const __m512d vh = _mm512_maskz_loadu_pd(mask, hap + i);
        acc = _mm512_add_pd(
            acc, _mm512_mul_pd(vh, logAlphaVec512(_mm512_mul_pd(vt, vr), vu, vq,
                                                  /*guard=*/false)));
    }

    return _mm512_reduce_add_pd(acc);
}

__attribute__((target("avx2,fma")))
double binomHapcountK0_avx2(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q
) noexcept {
    const __m256d vt = _mm256_set1_pd(t);
    const __m256d vq = _mm256_set1_pd(q);
    const __m256d vu = _mm256_set1_pd(1.0 - q);

    __m256d acc = _mm256_setzero_pd();

    int i = 0;
    const int n4 = n & ~3;
    for (; i < n4; i += 4) {
        const __m256d vr = _mm256_loadu_pd(resid + i);
        const __m256d vh = _mm256_loadu_pd(hap + i);
        acc = _mm256_add_pd(
            acc, _mm256_mul_pd(vh, logAlphaVec256(_mm256_mul_pd(vt, vr), vu, vq,
                                                  /*guard=*/false)));
    }

    double s = hsum256(acc);
    for (; i < n; ++i) s += hap[i] * logAlphaFast(t * resid[i], q);
    return s;
}

#endif  // x86_64

}  // namespace tier

namespace {

using HapcountK12Fn =
    Cgf12 (*)(double, const double *, const double *, int, double) noexcept;
using HapcountK0Fn =
    double (*)(double, const double *, const double *, int, double) noexcept;

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

HapcountK0Fn pickHapcountK0Fn() {
#if defined(__x86_64__) || defined(_M_X64)
    switch (simdLevel()) {
    case SimdLevel::AVX512: return tier::binomHapcountK0_avx512;
    case SimdLevel::AVX2:   return tier::binomHapcountK0_avx2;
    default: break;
    }
#endif
    return tier::binomHapcountK0_scalar;
}

const HapcountK12Fn hapcountK12Fn = pickHapcountK12Fn();
const HapcountK0Fn  hapcountK0Fn  = pickHapcountK0Fn();

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

// The q endpoints, exactly as for binomUniform.  q == 1: G_i == hap_i a.s., so
// K = t*sum(hap_i r_i).  q == 0: every cumulant vanishes.
namespace {

bool hapcountKEndpoint(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q,
    double *k0
) {
    if (q >= 1.0) {
        double sumHR = 0.0;
        for (int i = 0; i < n; ++i) sumHR += hap[i] * resid[i];
        *k0 = t * sumHR;
        return true;
    }
    if (q <= 0.0) { *k0 = 0.0; return true; }
    return false;
}

}  // namespace

Cgf012 binomHapcountKFull(
    double t,
    const double *resid,
    const double *hap,
    int n,
    double q
) noexcept {
    const Cgf12 d = binomHapcountK12(t, resid, hap, n, q);
    if (n <= 0) return Cgf012{0.0, d.K1, d.K2};

    double k0;
    if (hapcountKEndpoint(t, resid, hap, n, q, &k0))
        return Cgf012{k0, d.K1, d.K2};

    return Cgf012{hapcountK0Fn(t, resid, hap, n, q), d.K1, d.K2};
}

Cgf012 binomHapcountKFullExact(
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
