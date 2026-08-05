// tinytest.hpp — Minimal header-only test harness for GRAB.
//
// GRAB vendors no test framework and must not acquire one: CLAUDE.md forbids
// adding third-party libraries that are not vendored under third_party/, and a
// unit-test dependency would have to be carried on all three platforms for no
// runtime benefit.  This header supplies the smallest facility the SPA test
// suite needs — test registration, scalar and floating-point assertions, and a
// main() that reports counts and exits non-zero on failure.
//
// Usage:
//     #include "tinytest.hpp"
//     TEST(kernel_matches_reference) {
//         CHECK_REL(newForm(t, p), reference(t, p), 4e-16);
//     }
//     TINYTEST_MAIN
//
// Assertions record a failure and continue within the test body; a test that
// records at least one failure is reported as FAIL.  REQUIRE_* variants abort
// the current test body instead, for cases where continuing would be
// meaningless (a null pointer, a size mismatch before an element-wise loop).

#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <string>
#include <vector>

namespace tinytest {

// ──────────────────────────────────────────────────────────────────────
// Registry
// ──────────────────────────────────────────────────────────────────────

struct Case {
    const char *name;
    void (*fn)();
};

inline std::vector<Case> &registry() {
    static std::vector<Case> r;
    return r;
}

// Per-test failure counter.  Reset by the runner before each case.
inline int &failuresInCase() {
    static int n = 0;
    return n;
}

// Thrown by REQUIRE_* to abandon the remainder of a test body.
struct Abort {};

struct Registrar {
    Registrar(const char *name, void (*fn)()) {
        registry().push_back(Case{name, fn});
    }
};

// ──────────────────────────────────────────────────────────────────────
// Failure reporting
// ──────────────────────────────────────────────────────────────────────

inline void reportFailure(
    const char *file,
    int line,
    const char *expr,
    const std::string &detail
) {
    ++failuresInCase();
    std::fprintf(stderr, "    FAIL %s:%d\n      %s\n", file, line, expr);
    if (!detail.empty())
        std::fprintf(stderr, "      %s\n", detail.c_str());
}

// Relative difference |a-b| / max(|a|,|b|), with exact equality (including
// both-zero and both-infinite-same-sign) reported as 0.
inline double relDiff(double a, double b) {
    if (a == b) return 0.0;
    if (std::isnan(a) || std::isnan(b)) return std::numeric_limits<double>::infinity();
    const double m = std::fmax(std::fabs(a), std::fabs(b));
    if (m == 0.0) return 0.0;
    return std::fabs(a - b) / m;
}

inline std::string fmt2(double a, double b, double got, double want) {
    char buf[256];
    std::snprintf(buf, sizeof(buf),
                  "lhs=%.17g  rhs=%.17g  reldiff=%.3g  tol=%.3g", a, b, got, want);
    return std::string(buf);
}

// ──────────────────────────────────────────────────────────────────────
// Runner
// ──────────────────────────────────────────────────────────────────────

inline int run(int argc, char **argv) {
    const char *filter = (argc > 1) ? argv[1] : nullptr;
    int nRun = 0, nFail = 0;

    for (const Case &c : registry()) {
        if (filter && std::string(c.name).find(filter) == std::string::npos)
            continue;
        ++nRun;
        failuresInCase() = 0;
        std::fprintf(stderr, "  [ RUN  ] %s\n", c.name);
        try {
            c.fn();
        } catch (const Abort &) {
            // REQUIRE_* already recorded the failure.
        } catch (const std::exception &e) {
            reportFailure(__FILE__, __LINE__, "uncaught std::exception", e.what());
        } catch (...) {
            reportFailure(__FILE__, __LINE__, "uncaught unknown exception", "");
        }
        if (failuresInCase() == 0) {
            std::fprintf(stderr, "  [  OK  ] %s\n", c.name);
        } else {
            std::fprintf(stderr, "  [ FAIL ] %s  (%d assertion failure(s))\n",
                         c.name, failuresInCase());
            ++nFail;
        }
    }

    std::fprintf(stderr, "\n%d test(s) run, %d passed, %d failed\n",
                 nRun, nRun - nFail, nFail);
    return nFail == 0 ? 0 : 1;
}

} // namespace tinytest

// ──────────────────────────────────────────────────────────────────────
// Macros
// ──────────────────────────────────────────────────────────────────────

#define TINYTEST_CAT2(a, b) a##b
#define TINYTEST_CAT(a, b) TINYTEST_CAT2(a, b)

#define TEST(name)                                                            \
    static void TINYTEST_CAT(tinytest_fn_, name)();                           \
    static ::tinytest::Registrar TINYTEST_CAT(tinytest_reg_, name)(           \
        #name, &TINYTEST_CAT(tinytest_fn_, name));                            \
    static void TINYTEST_CAT(tinytest_fn_, name)()

#define CHECK(cond)                                                           \
    do {                                                                      \
        if (!(cond)) ::tinytest::reportFailure(__FILE__, __LINE__, #cond, ""); \
    } while (0)

#define CHECK_MSG(cond, msg)                                                  \
    do {                                                                      \
        if (!(cond))                                                          \
            ::tinytest::reportFailure(__FILE__, __LINE__, #cond, (msg));      \
    } while (0)

// Relative-tolerance floating-point comparison.  Prefer this to CHECK_NEAR
// for quantities whose magnitude varies across the tested domain.
#define CHECK_REL(a, b, tol)                                                  \
    do {                                                                      \
        const double tt_a = (a), tt_b = (b), tt_t = (tol);                    \
        const double tt_d = ::tinytest::relDiff(tt_a, tt_b);                  \
        if (!(tt_d <= tt_t))                                                  \
            ::tinytest::reportFailure(__FILE__, __LINE__,                     \
                                      #a " ~= " #b,                           \
                                      ::tinytest::fmt2(tt_a, tt_b, tt_d, tt_t)); \
    } while (0)

// Mixed absolute/relative comparison: passes when
//     |a − b| <= atol + rtol·max(|a|,|b|).
// This is the correct criterion for a quantity that passes through zero
// within the tested domain — a cumulant generating function at t = 0, for
// instance — where a pure relative tolerance is unsatisfiable by construction.
#define CHECK_CLOSE(a, b, rtol, atol)                                         \
    do {                                                                      \
        const double tt_a = (a), tt_b = (b);                                  \
        const double tt_r = (rtol), tt_at = (atol);                           \
        const double tt_m = std::fmax(std::fabs(tt_a), std::fabs(tt_b));      \
        const double tt_e = std::fabs(tt_a - tt_b);                           \
        const double tt_lim = tt_at + tt_r * tt_m;                            \
        if (!(tt_e <= tt_lim))                                                \
            ::tinytest::reportFailure(__FILE__, __LINE__,                     \
                                      #a " ~= " #b,                           \
                                      ::tinytest::fmt2(tt_a, tt_b, tt_e, tt_lim)); \
    } while (0)

// Absolute-tolerance comparison.  Use for quantities with a natural zero.
#define CHECK_NEAR(a, b, tol)                                                 \
    do {                                                                      \
        const double tt_a = (a), tt_b = (b), tt_t = (tol);                    \
        const double tt_d = std::fabs(tt_a - tt_b);                           \
        if (!(tt_d <= tt_t))                                                  \
            ::tinytest::reportFailure(__FILE__, __LINE__,                     \
                                      #a " ~= " #b,                           \
                                      ::tinytest::fmt2(tt_a, tt_b, tt_d, tt_t)); \
    } while (0)

#define REQUIRE(cond)                                                         \
    do {                                                                      \
        if (!(cond)) {                                                        \
            ::tinytest::reportFailure(__FILE__, __LINE__, #cond, "(required)"); \
            throw ::tinytest::Abort{};                                        \
        }                                                                     \
    } while (0)

#define TINYTEST_MAIN                                                         \
    int main(int argc, char **argv) { return ::tinytest::run(argc, argv); }
