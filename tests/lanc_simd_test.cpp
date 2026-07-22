// lanc_simd_test.cpp — P2 gate: byte-identity of the SIMD unpack kernels.
//
// Verifies that lanc_simd::unpackNibbles and lanc_simd::unpack2bit produce
// output BYTE-IDENTICAL to their scalar reference for every SIMD tier the host
// supports, across a range of n that exercises small values, exact multiples
// of the lane count, and multiples +/- 1/2/3 (tails and partial 2-bit bytes).
//
// The scalar variant is the reference.  Each host-supported variant (AVX2 if
// simdLevel() >= AVX2, AVX-512 if >= AVX512) and the dispatched entry point are
// compared against it.  Variants whose ISA the host does not support are never
// called, avoiding SIGILL.  A 64-byte guard region past each output buffer is
// checked to catch any over-write by a SIMD store.
//
// Deterministic PRNG (splitmix64, fixed seed): the environment forbids
// time()/random_device.
//
// Build (standalone; mirrors the P1 round-trip test):
//   g++ -std=c++17 -O2 -Wall -Wextra -Isrc -Ithird_party/eigen-5.0.0
//       -Ithird_party/zstd-1.5.7/lib tests/lanc_simd_test.cpp
//       src/localplus/lanc_io.cpp src/geno_factory/variant_filter.cpp
//       build/zstd/common/*.o build/zstd/compress/*.o build/zstd/decompress/*.o
//       -lstdc++fs -o <bin>
//
// Prints the levels exercised and a total check count; exits non-zero on any
// mismatch.

#include "localplus/lanc_io.hpp"

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

// ── failure tracking ────────────────────────────────────────────────────
static long g_checks = 0;
static long g_failures = 0;

// ── deterministic PRNG (splitmix64, fixed seed) ─────────────────────────
struct Splitmix64 {
    uint64_t state;
    explicit Splitmix64(uint64_t seed) : state(seed) {}
    uint64_t next() {
        uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
        return z ^ (z >> 31);
    }
};

using UnpackFn = void (*)(const uint8_t *, int, uint8_t *, uint8_t *);

// Run one variant into guard-padded output buffers and assert byte-identity
// with the scalar reference (refA/refB, length n) plus an intact guard tail.
static void runOne(
    const char *label,
    UnpackFn f,
    int n,
    const std::vector<uint8_t> &in,
    const std::vector<uint8_t> &refA,
    const std::vector<uint8_t> &refB
) {
    const uint8_t GUARD = 0xA5;
    const int PAD = 64;
    std::vector<uint8_t> A(static_cast<size_t>(n) + PAD, GUARD);
    std::vector<uint8_t> B(static_cast<size_t>(n) + PAD, GUARD);

    f(in.data(), n, A.data(), B.data());

    ++g_checks;
    for (int i = 0; i < n; ++i) {
        if (A[i] != refA[i] || B[i] != refB[i]) {
            ++g_failures;
            std::fprintf(stderr, "  FAIL %s n=%d i=%d: got(%u,%u) exp(%u,%u)\n",
                         label, n, i, A[i], B[i], refA[i], refB[i]);
            return;
        }
    }
    for (int i = n; i < n + PAD; ++i) {
        if (A[i] != GUARD || B[i] != GUARD) {
            ++g_failures;
            std::fprintf(stderr, "  FAIL %s n=%d: SIMD store overwrote output guard byte %d\n",
                         label, n, i);
            return;
        }
    }
}

int main() {
    const int lvl = lanc_simd::simdLevelValue();
    const char *lvlName = (lvl >= 2) ? "AVX512" : (lvl >= 1) ? "AVX2" : "Scalar";
    std::printf("[lanc_simd] runtime SIMD level: %d (%s)\n", lvl, lvlName);

    // n values: small, exact multiples of the lane counts (32/64), and +/-1/2/3
    // around them, to exercise scalar tails and partial trailing 2-bit bytes.
    const int ns[] = {0,   1,   2,   3,   4,   5,   7,   8,   15,  16,  17, 31,
                      32,  33,  63,  64,  65,  100, 255, 256, 257, 1000, 4095};

    Splitmix64 rng(0x0123456789ABCDEFULL);

    bool dispRun = false, avx2Run = false, avx512Run = false;

    for (int n : ns) {
        // ── Kernel 1: unpackNibbles — input is n bytes (1 byte/individual) ──
        {
            std::vector<uint8_t> in(static_cast<size_t>(n));
            for (auto &b : in) b = static_cast<uint8_t>(rng.next());

            std::vector<uint8_t> refA(static_cast<size_t>(n));
            std::vector<uint8_t> refB(static_cast<size_t>(n));
            lanc_simd::unpackNibbles_scalar(in.data(), n, refA.data(), refB.data());

            runOne("nibbles/dispatch", lanc_simd::unpackNibbles, n, in, refA, refB);
            dispRun = true;
#if defined(__x86_64__) || defined(_M_X64)
            if (lvl >= 1) {
                runOne("nibbles/avx2", lanc_simd::unpackNibbles_avx2, n, in, refA, refB);
                avx2Run = true;
            }
            if (lvl >= 2) {
                runOne("nibbles/avx512", lanc_simd::unpackNibbles_avx512, n, in, refA, refB);
                avx512Run = true;
            }
#endif
        }

        // ── Kernel 2: unpack2bit — input is ceil(n/4) bytes (2 bits/ind) ──
        {
            const int inBytes = (n + 3) / 4;
            std::vector<uint8_t> in(static_cast<size_t>(inBytes));
            for (auto &b : in) b = static_cast<uint8_t>(rng.next());

            std::vector<uint8_t> ref0(static_cast<size_t>(n));
            std::vector<uint8_t> ref1(static_cast<size_t>(n));
            lanc_simd::unpack2bit_scalar(in.data(), n, ref0.data(), ref1.data());

            runOne("2bit/dispatch", lanc_simd::unpack2bit, n, in, ref0, ref1);
#if defined(__x86_64__) || defined(_M_X64)
            if (lvl >= 1)
                runOne("2bit/avx2", lanc_simd::unpack2bit_avx2, n, in, ref0, ref1);
            if (lvl >= 2)
                runOne("2bit/avx512", lanc_simd::unpack2bit_avx512, n, in, ref0, ref1);
#endif
        }
    }

    std::printf("[lanc_simd] levels exercised: scalar(reference)%s%s%s\n",
                dispRun ? " + dispatch" : "",
                avx2Run ? " + AVX2" : "",
                avx512Run ? " + AVX512" : "");
    std::printf("[lanc_simd] %ld checks run, %ld failures\n", g_checks, g_failures);

    if (g_failures == 0) {
        std::printf("PASS: all %ld checks passed\n", g_checks);
        return 0;
    }
    std::printf("FAIL: %ld of %ld checks failed\n", g_failures, g_checks);
    return 1;
}
