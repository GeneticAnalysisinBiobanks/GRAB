// simd_dispatch.hpp — Runtime SIMD tier detection
//
// Detects AVX-512 / AVX2 / scalar at process startup via __builtin_cpu_supports().
// Used by hot kernels (phi variance, IBD accumulation) to dispatch to the
// widest available SIMD implementation at runtime, enabling a single binary
// that exploits AVX-512 when present and falls back gracefully.
#pragma once

enum class SimdLevel : int { Scalar = 0, AVX2 = 1, AVX512 = 2 };

inline SimdLevel simdLevel() {
    static const SimdLevel level = [] {
        __builtin_cpu_init();
        if (__builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512vl"))
            return SimdLevel::AVX512;
        if (__builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma"))
            return SimdLevel::AVX2;
        return SimdLevel::Scalar;
    }();
    return level;
}
