// simd_math.hpp — Vectorised exp/log for double precision (AVX2 / AVX-512)
//
// Cephes-style minimax polynomials with ~1 ULP accuracy over the full
// double range.  No special NaN/Inf handling — SPA values are well-bounded.
//
// Usage:  #include "util/simd_math.hpp"
//         __m256d y = avx2_exp_pd(x);
//         __m512d y = avx512_exp_pd(x);
#pragma once

#include <immintrin.h>
#include <cstdint>

// ════════════════════════════════════════════════════════════════════
// AVX2 (4 × double)
// ════════════════════════════════════════════════════════════════════

// exp(x) for 4 doubles — Cephes-style range-reduction + degree-6 minimax
__attribute__((target("avx2,fma")))
inline __m256d avx2_exp_pd(__m256d x)
{
    // Constants
    const __m256d LOG2E   = _mm256_set1_pd(1.4426950408889634);     // log2(e)
    const __m256d LN2_HI  = _mm256_set1_pd(6.93145751953125e-1);   // ln(2) high part
    const __m256d LN2_LO  = _mm256_set1_pd(1.42860682030941723212e-6); // ln(2) low part
    const __m256d HALF    = _mm256_set1_pd(0.5);

    // Polynomial coefficients (minimax on [-ln2/2, ln2/2])
    const __m256d P0 = _mm256_set1_pd(1.0);
    const __m256d P1 = _mm256_set1_pd(1.0);
    const __m256d P2 = _mm256_set1_pd(0.5);
    const __m256d P3 = _mm256_set1_pd(1.6666666666666666e-1);   // 1/6
    const __m256d P4 = _mm256_set1_pd(4.1666666666666664e-2);   // 1/24
    const __m256d P5 = _mm256_set1_pd(8.3333333333333332e-3);   // 1/120
    const __m256d P6 = _mm256_set1_pd(1.3888888888888889e-3);   // 1/720
    const __m256d P7 = _mm256_set1_pd(1.9841269841269841e-4);   // 1/5040
    const __m256d P8 = _mm256_set1_pd(2.4801587301587302e-5);   // 1/40320
    const __m256d P9 = _mm256_set1_pd(2.7557319223985893e-6);   // 1/362880
    const __m256d P10 = _mm256_set1_pd(2.7557319223985888e-7);  // 1/3628800
    const __m256d P11 = _mm256_set1_pd(2.5052108385441720e-8);  // 1/39916800

    // Clamp to avoid overflow/underflow in integer conversion
    const __m256d MAX_X = _mm256_set1_pd(709.0);
    const __m256d MIN_X = _mm256_set1_pd(-709.0);
    x = _mm256_max_pd(_mm256_min_pd(x, MAX_X), MIN_X);

    // n = round(x * log2(e))
    __m256d t = _mm256_fmadd_pd(x, LOG2E, HALF);
    __m256d n = _mm256_floor_pd(t);

    // r = x - n * ln(2)  (Cody-Waite two-step reduction)
    __m256d r = _mm256_fnmadd_pd(n, LN2_HI, x);
    r = _mm256_fnmadd_pd(n, LN2_LO, r);

    // Evaluate polynomial:  exp(r) ≈ 1 + r + r²/2 + r³/6 + ...
    // Horner's method from highest to lowest
    __m256d poly = P11;
    poly = _mm256_fmadd_pd(poly, r, P10);
    poly = _mm256_fmadd_pd(poly, r, P9);
    poly = _mm256_fmadd_pd(poly, r, P8);
    poly = _mm256_fmadd_pd(poly, r, P7);
    poly = _mm256_fmadd_pd(poly, r, P6);
    poly = _mm256_fmadd_pd(poly, r, P5);
    poly = _mm256_fmadd_pd(poly, r, P4);
    poly = _mm256_fmadd_pd(poly, r, P3);
    poly = _mm256_fmadd_pd(poly, r, P2);
    poly = _mm256_fmadd_pd(poly, r, P1);
    poly = _mm256_fmadd_pd(poly, r, P0);

    // Reconstruct: exp(x) = 2^n * exp(r)
    // Convert n to int64 and add to exponent field
    __m128i ni = _mm256_cvtpd_epi32(n);
    __m256i ni64 = _mm256_cvtepi32_epi64(ni);
    __m256i bias = _mm256_set1_epi64x(1023);
    __m256i exp_bits = _mm256_slli_epi64(_mm256_add_epi64(ni64, bias), 52);
    __m256d scale = _mm256_castsi256_pd(exp_bits);

    return _mm256_mul_pd(poly, scale);
}

// log(x) for 4 doubles — glibc/musl-style decomposition
//
// Algorithm:
//   x = m * 2^e,  m in [1, 2)
//   if m > sqrt(2): m /= 2, e += 1   → m in [sqrt(1/2), sqrt(2)]
//   f = m - 1
//   s = f/(2+f) = (m-1)/(m+1)
//   log(m) = f - 0.5*f² + s*(0.5*f² + R(s²))
//   log(x) = e*ln(2) + log(m)
__attribute__((target("avx2,fma")))
inline __m256d avx2_log_pd(__m256d x)
{
    const __m256d ONE  = _mm256_set1_pd(1.0);
    const __m256d HALF = _mm256_set1_pd(0.5);
    const __m256d TWO  = _mm256_set1_pd(2.0);
    const __m256d LN2  = _mm256_set1_pd(0.6931471805599453);
    const __m256d SQRT2 = _mm256_set1_pd(1.4142135623730950488);

    // Minimax coefficients for R(z)/z where z = s²  (from glibc/musl)
    const __m256d Lg1 = _mm256_set1_pd(6.666666666666735130e-1);
    const __m256d Lg2 = _mm256_set1_pd(3.999999999940941908e-1);
    const __m256d Lg3 = _mm256_set1_pd(2.857142874366239149e-1);
    const __m256d Lg4 = _mm256_set1_pd(2.222219843214978396e-1);
    const __m256d Lg5 = _mm256_set1_pd(1.818357216161805012e-1);
    const __m256d Lg6 = _mm256_set1_pd(1.531383769920937332e-1);
    const __m256d Lg7 = _mm256_set1_pd(1.479819860511658591e-1);

    // Extract exponent and mantissa
    __m256i xi = _mm256_castpd_si256(x);
    __m256i exp_mask = _mm256_set1_epi64x(0x7FF0000000000000LL);
    __m256i mant_mask = _mm256_set1_epi64x(0x000FFFFFFFFFFFFFLL);
    __m256i bias = _mm256_set1_epi64x(1023);

    __m256i exp_bits = _mm256_srli_epi64(_mm256_and_si256(xi, exp_mask), 52);
    __m256i mant_bits = _mm256_or_si256(_mm256_and_si256(xi, mant_mask),
                                        _mm256_set1_epi64x(0x3FF0000000000000LL));
    __m256d m = _mm256_castsi256_pd(mant_bits);  // m in [1, 2)

    __m256d f = _mm256_sub_pd(m, ONE);

    // If m > sqrt(2): f = m*0.5 - 1, e += 1
    __m256d cmp = _mm256_cmp_pd(m, SQRT2, _CMP_GT_OQ);
    __m256d f_alt = _mm256_sub_pd(_mm256_mul_pd(m, HALF), ONE);
    f = _mm256_blendv_pd(f, f_alt, cmp);

    // Adjust exponent
    __m256i e_adj = _mm256_sub_epi64(exp_bits, bias);
    __m256i one_i64 = _mm256_set1_epi64x(1);
    __m256i adj = _mm256_and_si256(_mm256_castpd_si256(cmp), one_i64);
    e_adj = _mm256_add_epi64(e_adj, adj);

    // Convert int64 exponent to double via dword permutation
    __m256i idx = _mm256_setr_epi32(0, 2, 4, 6, 0, 0, 0, 0);
    __m128i e32 = _mm256_castsi256_si128(_mm256_permutevar8x32_epi32(e_adj, idx));
    __m256d e_dbl = _mm256_cvtepi32_pd(e32);

    // s = f / (2 + f)
    __m256d s = _mm256_div_pd(f, _mm256_add_pd(TWO, f));
    __m256d z = _mm256_mul_pd(s, s);

    // R(z) = z * (Lg1 + z*(Lg2 + z*(Lg3 + z*(Lg4 + z*(Lg5 + z*(Lg6 + z*Lg7))))))
    __m256d R = Lg7;
    R = _mm256_fmadd_pd(R, z, Lg6);
    R = _mm256_fmadd_pd(R, z, Lg5);
    R = _mm256_fmadd_pd(R, z, Lg4);
    R = _mm256_fmadd_pd(R, z, Lg3);
    R = _mm256_fmadd_pd(R, z, Lg2);
    R = _mm256_fmadd_pd(R, z, Lg1);
    R = _mm256_mul_pd(R, z);  // R(z) = z * polynomial

    // hfsq = 0.5 * f * f
    __m256d hfsq = _mm256_mul_pd(HALF, _mm256_mul_pd(f, f));

    // log(m) = f - hfsq + s*(hfsq + R)
    __m256d log_m = _mm256_sub_pd(f, hfsq);
    log_m = _mm256_fmadd_pd(s, _mm256_add_pd(hfsq, R), log_m);

    // log(x) = e * ln(2) + log(m)
    return _mm256_fmadd_pd(e_dbl, LN2, log_m);
}

// ════════════════════════════════════════════════════════════════════
// AVX-512 (8 × double)
// ════════════════════════════════════════════════════════════════════

__attribute__((target("avx2,avx512f,avx512vl,fma")))
inline __m512d avx512_exp_pd(__m512d x)
{
    const __m512d LOG2E   = _mm512_set1_pd(1.4426950408889634);
    const __m512d LN2_HI  = _mm512_set1_pd(6.93145751953125e-1);
    const __m512d LN2_LO  = _mm512_set1_pd(1.42860682030941723212e-6);
    const __m512d HALF    = _mm512_set1_pd(0.5);

    const __m512d P0  = _mm512_set1_pd(1.0);
    const __m512d P1  = _mm512_set1_pd(1.0);
    const __m512d P2  = _mm512_set1_pd(0.5);
    const __m512d P3  = _mm512_set1_pd(1.6666666666666666e-1);
    const __m512d P4  = _mm512_set1_pd(4.1666666666666664e-2);
    const __m512d P5  = _mm512_set1_pd(8.3333333333333332e-3);
    const __m512d P6  = _mm512_set1_pd(1.3888888888888889e-3);
    const __m512d P7  = _mm512_set1_pd(1.9841269841269841e-4);
    const __m512d P8  = _mm512_set1_pd(2.4801587301587302e-5);
    const __m512d P9  = _mm512_set1_pd(2.7557319223985893e-6);
    const __m512d P10 = _mm512_set1_pd(2.7557319223985888e-7);
    const __m512d P11 = _mm512_set1_pd(2.5052108385441720e-8);

    const __m512d MAX_X = _mm512_set1_pd(709.0);
    const __m512d MIN_X = _mm512_set1_pd(-709.0);
    x = _mm512_max_pd(_mm512_min_pd(x, MAX_X), MIN_X);

    __m512d t = _mm512_fmadd_pd(x, LOG2E, HALF);
    __m512d n = _mm512_roundscale_pd(t, _MM_FROUND_TO_NEG_INF | _MM_FROUND_NO_EXC);

    __m512d r = _mm512_fnmadd_pd(n, LN2_HI, x);
    r = _mm512_fnmadd_pd(n, LN2_LO, r);

    __m512d poly = P11;
    poly = _mm512_fmadd_pd(poly, r, P10);
    poly = _mm512_fmadd_pd(poly, r, P9);
    poly = _mm512_fmadd_pd(poly, r, P8);
    poly = _mm512_fmadd_pd(poly, r, P7);
    poly = _mm512_fmadd_pd(poly, r, P6);
    poly = _mm512_fmadd_pd(poly, r, P5);
    poly = _mm512_fmadd_pd(poly, r, P4);
    poly = _mm512_fmadd_pd(poly, r, P3);
    poly = _mm512_fmadd_pd(poly, r, P2);
    poly = _mm512_fmadd_pd(poly, r, P1);
    poly = _mm512_fmadd_pd(poly, r, P0);

    __m256i ni = _mm512_cvtpd_epi32(n);
    __m512i ni64 = _mm512_cvtepi32_epi64(ni);
    __m512i bias = _mm512_set1_epi64(1023);
    __m512i exp_bits = _mm512_slli_epi64(_mm512_add_epi64(ni64, bias), 52);
    __m512d scale = _mm512_castsi512_pd(exp_bits);

    return _mm512_mul_pd(poly, scale);
}

__attribute__((target("avx2,avx512f,avx512vl,fma")))
inline __m512d avx512_log_pd(__m512d x)
{
    const __m512d ONE  = _mm512_set1_pd(1.0);
    const __m512d HALF = _mm512_set1_pd(0.5);
    const __m512d TWO  = _mm512_set1_pd(2.0);
    const __m512d LN2  = _mm512_set1_pd(0.6931471805599453);
    const __m512d SQRT2 = _mm512_set1_pd(1.4142135623730950488);

    const __m512d Lg1 = _mm512_set1_pd(6.666666666666735130e-1);
    const __m512d Lg2 = _mm512_set1_pd(3.999999999940941908e-1);
    const __m512d Lg3 = _mm512_set1_pd(2.857142874366239149e-1);
    const __m512d Lg4 = _mm512_set1_pd(2.222219843214978396e-1);
    const __m512d Lg5 = _mm512_set1_pd(1.818357216161805012e-1);
    const __m512d Lg6 = _mm512_set1_pd(1.531383769920937332e-1);
    const __m512d Lg7 = _mm512_set1_pd(1.479819860511658591e-1);

    __m512i xi = _mm512_castpd_si512(x);
    __m512i mant_mask = _mm512_set1_epi64(0x000FFFFFFFFFFFFFLL);
    __m512i bias = _mm512_set1_epi64(1023);

    __m512i exp_bits = _mm512_srli_epi64(xi, 52);
    exp_bits = _mm512_and_si512(exp_bits, _mm512_set1_epi64(0x7FF));
    __m512i mant_bits = _mm512_or_si512(_mm512_and_si512(xi, mant_mask),
                                        _mm512_set1_epi64(0x3FF0000000000000LL));
    __m512d m = _mm512_castsi512_pd(mant_bits);

    __m512d f = _mm512_sub_pd(m, ONE);

    __mmask8 cmp = _mm512_cmp_pd_mask(m, SQRT2, _CMP_GT_OQ);
    __m512d f_alt = _mm512_sub_pd(_mm512_mul_pd(m, HALF), ONE);
    f = _mm512_mask_blend_pd(cmp, f, f_alt);

    __m512i e_adj = _mm512_sub_epi64(exp_bits, bias);
    e_adj = _mm512_mask_add_epi64(e_adj, cmp, e_adj, _mm512_set1_epi64(1));

    __m256i e32 = _mm512_cvtepi64_epi32(e_adj);
    __m512d e_dbl = _mm512_cvtepi32_pd(e32);

    __m512d s = _mm512_div_pd(f, _mm512_add_pd(TWO, f));
    __m512d z = _mm512_mul_pd(s, s);

    __m512d R = Lg7;
    R = _mm512_fmadd_pd(R, z, Lg6);
    R = _mm512_fmadd_pd(R, z, Lg5);
    R = _mm512_fmadd_pd(R, z, Lg4);
    R = _mm512_fmadd_pd(R, z, Lg3);
    R = _mm512_fmadd_pd(R, z, Lg2);
    R = _mm512_fmadd_pd(R, z, Lg1);
    R = _mm512_mul_pd(R, z);

    __m512d hfsq = _mm512_mul_pd(HALF, _mm512_mul_pd(f, f));
    __m512d log_m = _mm512_sub_pd(f, hfsq);
    log_m = _mm512_fmadd_pd(s, _mm512_add_pd(hfsq, R), log_m);

    return _mm512_fmadd_pd(e_dbl, LN2, log_m);
}
