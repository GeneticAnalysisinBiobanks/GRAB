// spa_solver_test.cpp — Tests for src/util/spa.hpp.
//
// Stage 1 of the saddlepoint unification (dev-notes/methods/spa_unify).  The
// header under test has no callers in src/ yet, so this binary is the only
// thing that exercises it and is correspondingly the whole gate.
//
// What is proved here:
//
//   § 1  Closed-form saddlepoints.  Gaussian, Poisson and single-subject
//        binomial CGFs all have an analytic root; the solver must recover it
//        to ~1e-14 relative.  These are the only cases where "correct" is
//        knowable independently of the solver, so they come first.
//   § 2  Property test over random problems: the achieved residual satisfies
//        the *stated relative* criterion, and K0/K2 are the cumulants at the
//        returned zeta.  The second half is the test that would have caught
//        SPACox's stale-K'' hole, where the finder returned a post-step
//        abscissa together with the pre-step K''.
//   § 3  Every Status value constructed and returned.
//   § 4  The bracket invariant, checked at every iteration through the
//        instrumentation hook.
//   § 5  bnTail / bnTailLog agreement, and log-domain survival past the point
//        where the linear scale reaches zero.
//   § 6  bnTail's treatment of (K0, K2) as opaque scalars (WtCoxG's D4 pair).
//   § 7  Two-sided assembly: no half-sized p-values, no NaN becoming 1.0, no
//        DBL_MIN substitution.
//   § 8  The normal short-circuit.
//
// Every random draw is seeded from a fixed constant so the binary is
// reproducible run to run.
//
// Builds with no extra objects: spa.hpp reaches only math::pnorm, which is
// inline in util/math_helper.hpp.  `make test` therefore picks this file up
// from the tests/*_test.cpp wildcard with no Makefile entry.

#include "tinytest.hpp"

#include "util/spa.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <random>
#include <string>
#include <vector>

namespace {

const double kInf = std::numeric_limits<double>::infinity();
const double kNaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kEps = std::numeric_limits<double>::epsilon();
constexpr double kLn10 = 2.30258509299404568402;
constexpr double kLn2  = 0.69314718055994530942;

std::string statusStr(spa::Status s) { return std::string(spa::statusName(s)); }

// Solver options that demand the tightest root the arithmetic allows: no
// relative slack, no bracket-collapse slack, so the loop runs until the
// residual reaches the representation-error floor of K'(t) - s and of the
// abscissa itself, or the bracket is down to adjacent doubles.  With rtol = 0
// the only surviving terms of the convergence criterion are 4*eps*|s| and
// 4*eps*K''(t)*|t|, which is exactly what "as tight as double precision
// allows" means; `Converged` is still a statement about the residual, so these
// options ask for the most it can be.
spa::SolveOpts tightOpts(double init = 0.0) {
    spa::SolveOpts o;
    o.init = init;
    o.rtol = 0.0;
    o.bracketRtol = 0.0;
    o.maxIter = 300;
    return o;
}

// ──────────────────────────────────────────────────────────────────────
// Reference CGFs with analytic saddlepoints
// ──────────────────────────────────────────────────────────────────────

// S ~ Normal(mu, var).  K(t) = mu*t + var*t^2/2, so K'(t) = mu + var*t and
// the root of K'(zeta) = s is zeta = (s - mu)/var exactly.
struct GaussianCgf {
    double mu, var, s;
    spa::K12  k12(double t) const { return spa::K12{mu + var * t - s, var}; }
    spa::K012 kFull(double t) const {
        return spa::K012{mu * t + 0.5 * var * t * t, mu + var * t - s, var};
    }
    double root() const { return (s - mu) / var; }
};

// S ~ Poisson(lambda).  K(t) = lambda*(e^t - 1), K' = K'' = lambda*e^t, so
// zeta = log(s/lambda).
struct PoissonCgf {
    double lambda, s;
    spa::K12 k12(double t) const {
        const double e = lambda * std::exp(t);
        return spa::K12{e - s, e};
    }
    spa::K012 kFull(double t) const {
        const double e = lambda * std::exp(t);
        return spa::K012{lambda * std::expm1(t), e - s, e};
    }
    double root() const { return std::log(s / lambda); }
};

// One subject contributing r * G with G ~ Binomial(2, p).  This is the atom
// out of which all four of GRAB's binomial CGFs are summed.
//
//   K(t)   = 2*log(1 - p + p*e^{t*r})
//   K'(t)  = 2*r*p*e^{t*r} / (1 - p + p*e^{t*r})
//   K''(t) = 2*r^2*p*e^{t*r}*(1-p) / (1 - p + p*e^{t*r})^2
//
// Writing x = p*e^{zeta*r} / (1 - p + p*e^{zeta*r}) we have K' = 2*r*x, so
// K'(zeta) = s gives x = s/(2r) and
//
//   zeta = log( x*(1-p) / (p*(1-x)) ) / r,     valid for 0 < x < 1.
//
// K'' is written in the cancellation-free form 2*r^2*e*u/a^2 and K through
// log1p/expm1 (spa_unify D1 and N1), so the reference is not itself carrying
// the defects Stage 2 removes.
struct BinomOneCgf {
    double p, r, s;
    spa::K12 k12(double t) const {
        const double lam = std::exp(t * r);
        const double u = 1.0 - p, e = p * lam, a = u + e;
        return spa::K12{2.0 * r * e / a - s, 2.0 * r * r * e * u / (a * a)};
    }
    spa::K012 kFull(double t) const {
        const double lam = std::exp(t * r);
        const double u = 1.0 - p, e = p * lam, a = u + e;
        return spa::K012{2.0 * std::log1p(p * std::expm1(t * r)),
                         2.0 * r * e / a - s,
                         2.0 * r * r * e * u / (a * a)};
    }
    double root() const {
        const double x = s / (2.0 * r);
        return std::log(x * (1.0 - p) / (p * (1.0 - x))) / r;
    }
};

// A many-subject binomial CGF plus a Gaussian block, which is the shape every
// production site actually has: an outlier set evaluated exactly and a
// non-outlier remainder folded in as mean + var*t.
struct MixedCgf {
    std::vector<double> resid, maf;
    double mean = 0.0, var = 0.0, s = 0.0;

    spa::K12 k12(double t) const {
        double k1 = mean + var * t - s, k2 = var;
        for (size_t i = 0; i < resid.size(); ++i) {
            const double r = resid[i], p = maf[i];
            const double lam = std::exp(t * r);
            const double u = 1.0 - p, e = p * lam, a = u + e;
            k1 += 2.0 * r * e / a;
            k2 += 2.0 * r * r * e * u / (a * a);
        }
        return spa::K12{k1, k2};
    }
    spa::K012 kFull(double t) const {
        const spa::K12 c = k12(t);
        double k0 = mean * t + 0.5 * var * t * t;
        for (size_t i = 0; i < resid.size(); ++i)
            k0 += 2.0 * std::log1p(maf[i] * std::expm1(t * resid[i]));
        return spa::K012{k0, c.k1, c.k2};
    }
};

// Records every (iteration, iterate, bracket) the solver visits.
struct Recorder {
    struct Row { int iter; double t, k1, k2, lo, hi; };
    std::vector<Row> rows;
    void operator()(int iter, double t, double k1, double k2, double lo, double hi) {
        rows.push_back(Row{iter, t, k1, k2, lo, hi});
    }
};

} // namespace

// ══════════════════════════════════════════════════════════════════════
// § 1  Closed-form saddlepoints
// ══════════════════════════════════════════════════════════════════════

TEST(gaussian_closed_form_root) {
    const double mus[]  = {0.0, 0.3, -2.5, 100.0};
    const double vars[] = {1.0, 2.0, 0.25, 1e4};
    const double ss[]   = {1.5, -1.5, 40.0, -0.75};

    for (double mu : mus)
        for (double var : vars)
            for (double s : ss) {
                const GaussianCgf c{mu, var, s};
                const spa::Saddle r = spa::solveSaddlepoint(
                    s, [&](double t) { return c.k12(t); },
                    [&](double t) { return c.kFull(t); }, tightOpts());
                CHECK_MSG(r.status == spa::Status::SpaOk, statusStr(r.status));
                CHECK_CLOSE(r.zeta, c.root(), 1e-14, 1e-300);
                // Cumulants at the root, in closed form.
                CHECK_REL(r.K2, var, 0.0);
                CHECK_CLOSE(r.K0, mu * r.zeta + 0.5 * var * r.zeta * r.zeta,
                            1e-14, 1e-300);
            }
}

TEST(poisson_closed_form_root) {
    const double lambdas[] = {0.5, 1.0, 7.0, 1000.0};
    // s/lambda spans three orders of magnitude either side of 1, so both the
    // large-positive and large-negative root regimes are covered.
    const double ratios[] = {1e-3, 0.1, 0.5, 2.0, 10.0, 1e3};

    for (double lam : lambdas)
        for (double ratio : ratios) {
            const double s = lam * ratio;
            const PoissonCgf c{lam, s};
            const spa::Saddle r = spa::solveSaddlepoint(
                s, [&](double t) { return c.k12(t); },
                [&](double t) { return c.kFull(t); }, tightOpts());
            CHECK_MSG(r.status == spa::Status::SpaOk, statusStr(r.status));
            CHECK_CLOSE(r.zeta, c.root(), 1e-14, 1e-14);
            // K'(zeta) = lambda*e^zeta = s at the root, and K'' equals K'.
            CHECK_CLOSE(r.K2, s, 1e-13, 0.0);
        }
}

TEST(binomial_single_subject_closed_form_root) {
    const double ps[] = {0.01, 0.1, 0.3, 0.5, 0.9};
    const double rs[] = {0.25, 1.0, 3.0, -1.0, -2.5};
    // Fractions of the admissible range for s, i.e. x = s/(2r) in (0, 1).
    const double fracs[] = {0.05, 0.25, 0.5, 0.75, 0.95};

    int nCases = 0;
    for (double p : ps)
        for (double r0 : rs)
            for (double f : fracs) {
                const double s = f * 2.0 * r0;
                const BinomOneCgf c{p, r0, s};
                const spa::Saddle r = spa::solveSaddlepoint(
                    s, [&](double t) { return c.k12(t); },
                    [&](double t) { return c.kFull(t); }, tightOpts());
                CHECK_MSG(r.status == spa::Status::SpaOk,
                          statusStr(r.status) + " p=" + std::to_string(p) +
                              " r=" + std::to_string(r0) + " f=" + std::to_string(f));
                CHECK_CLOSE(r.zeta, c.root(), 1e-13, 1e-13);
                ++nCases;
            }
    CHECK(nCases == 125);
}

// The initial guess must not change the answer, only the work done.  Every
// production site supplies a per-marker guess of its own devising.
//
// The guess list deliberately includes the ROOT ITSELF, and the root the solver
// returned on a previous call.  Until the bracket-expansion repair those two
// entries failed, and they failed *because* they were good: the first outward
// probe is sized at the bare Newton step |K'(init) - s| / K''(init), so a guess
// already within half an ulp of the root produced a probe that rounded back onto
// `init`, the expansion loop broke on its first pass with no bracket, and the
// solve returned MaxIter — i.e. NA at the caller — while holding the exact root.
// A solver whose failure probability rises with the quality of its initial guess
// is not usable by any warm start or cache, so the property is pinned here.
TEST(root_independent_of_initial_guess) {
    // s = 3 (ratio 1) is excluded: the root is then zero and every guess is
    // trivially good.  s = 3000 and s = 0.003 are the two entries that failed
    // before the repair; whether a given problem falls inside the half-ulp
    // window depends on where fl(log(s/lambda)) lands relative to the true root,
    // so the sweep carries several ratios rather than one.
    for (double s : {0.003, 0.3, 6.0, 30.0, 300.0, 3000.0}) {
        const PoissonCgf c{3.0, s};

        // A cold solve first, so its returned root can be fed back in below.
        const spa::Saddle cold = spa::solveSaddlepoint(
            c.s, [&](double t) { return c.k12(t); },
            [&](double t) { return c.kFull(t); }, tightOpts(0.0));
        CHECK_MSG(cold.status == spa::Status::SpaOk, statusStr(cold.status));

        const double inits[] = {-20.0, -1.0, 0.0, 0.5, 2.0, 5.0, 50.0,
                                c.root(),      // the exact root
                                cold.zeta};    // a previously returned root
        double first = kNaN;
        for (double init : inits) {
            const spa::Saddle r = spa::solveSaddlepoint(
                c.s, [&](double t) { return c.k12(t); },
                [&](double t) { return c.kFull(t); }, tightOpts(init));
            CHECK_MSG(r.status == spa::Status::SpaOk,
                      statusStr(r.status) + " s=" + std::to_string(s) +
                          " init=" + std::to_string(init));
            CHECK_CLOSE(r.zeta, c.root(), 1e-14, 1e-14);
            if (std::isnan(first)) first = r.zeta;
            else CHECK_CLOSE(r.zeta, first, 1e-14, 1e-14);
        }
    }
}

// The same property, swept rather than sampled: for every closed-form problem
// in this file's grids, re-solving with `init` set to the analytic root and with
// `init` set to the root the solver just returned must converge and return that
// root.  Compiled against the pre-repair util/spa.hpp this case records 177
// assertion failures — 29 of the 197 problems below return MAXITER while holding
// the exact root — and against the repaired header, none.
TEST(solving_from_the_root_itself_converges) {
    int nCases = 0;

    // Both entry points must agree that `init` is already the answer.
    auto check = [&](const spa::Saddle &r, double root, const std::string &tag) {
        CHECK_MSG(r.status == spa::Status::SpaOk, statusStr(r.status) + " " + tag);
        CHECK_MSG(std::fabs(r.zeta - root) <=
                      1e-13 * std::fmax(1.0, std::fabs(root)),
                  "root moved, " + tag);
        // A degenerate bracket is still a bracket: the root is pinned to a
        // single representable double and nothing narrower exists.
        CHECK_MSG(r.bracketed, "not bracketed, " + tag);
        CHECK_MSG(r.zeta >= r.lo && r.zeta <= r.hi, "outside bracket, " + tag);
        // The terminal cumulants must be the ones at the returned root, on this
        // path as on every other.
        CHECK_MSG(r.K2 > 0.0, "K2 not positive, " + tag);
        ++nCases;
    };

    // ── Poisson: 16 of these 24 collapsed before the repair ──
    for (double lam : {0.5, 1.0, 3.0, 7.0, 1000.0})
        for (double ratio : {1e-3, 0.1, 0.5, 2.0, 10.0, 1e3}) {
            const PoissonCgf c{lam, lam * ratio};
            const spa::Saddle cold = spa::solveSaddlepoint(
                c.s, [&](double t) { return c.k12(t); },
                [&](double t) { return c.kFull(t); }, tightOpts(0.0));
            CHECK(cold.status == spa::Status::SpaOk);
            for (double init : {c.root(), cold.zeta}) {
                const spa::Saddle r = spa::solveSaddlepoint(
                    c.s, [&](double t) { return c.k12(t); },
                    [&](double t) { return c.kFull(t); }, tightOpts(init));
                check(r, c.root(),
                      "poisson lam=" + std::to_string(lam) + " ratio=" +
                          std::to_string(ratio));
            }
        }

    // ── Single-subject binomial: 13 of these 125 collapsed before ──
    for (double p : {0.01, 0.1, 0.3, 0.5, 0.9})
        for (double r0 : {0.25, 1.0, 3.0, -1.0, -2.5})
            for (double f : {0.05, 0.25, 0.5, 0.75, 0.95}) {
                const BinomOneCgf c{p, r0, f * 2.0 * r0};
                const spa::Saddle cold = spa::solveSaddlepoint(
                    c.s, [&](double t) { return c.k12(t); },
                    [&](double t) { return c.kFull(t); }, tightOpts(0.0));
                CHECK(cold.status == spa::Status::SpaOk);
                for (double init : {c.root(), cold.zeta}) {
                    const spa::Saddle r = spa::solveSaddlepoint(
                        c.s, [&](double t) { return c.k12(t); },
                        [&](double t) { return c.kFull(t); }, tightOpts(init));
                    check(r, c.root(),
                          "binom p=" + std::to_string(p) + " r=" +
                              std::to_string(r0) + " f=" + std::to_string(f));
                }
            }

    // ── Gaussian, including the two doubles adjacent to the exact root ──
    // The root is representable exactly here, so `init = root` takes the
    // zero-residual fast path; its neighbours are what probe the ulp-scale
    // regime on a CGF whose K'' is constant.
    for (double var : {1.0, 2.0, 0.25, 1e4})
        for (double s : {1.5, -1.5, 40.0, -0.75}) {
            const GaussianCgf c{0.3, var, s};
            for (int d = -1; d <= 1; ++d) {
                const double init =
                    (d == 0) ? c.root() : std::nextafter(c.root(), d * kInf);
                const spa::Saddle r = spa::solveSaddlepoint(
                    c.s, [&](double t) { return c.k12(t); },
                    [&](double t) { return c.kFull(t); }, tightOpts(init));
                check(r, c.root(),
                      "gauss var=" + std::to_string(var) + " s=" +
                          std::to_string(s) + " d=" + std::to_string(d));
            }
        }

    CHECK(nCases == 2 * 30 + 2 * 125 + 48);
}

// The warm start on the shape production actually has: solve, then re-solve
// from the returned root.  This is the reachability argument for the repair —
// the collapse is unreachable from the initializers the six methods ship today
// (0 occurrences in the 8288 solves examples/baseline.sh performs), but any
// per-marker cache or LOCO-style reuse turns it on immediately: against the
// pre-repair header this case reports 85 failures out of 991 re-solves.
TEST(warm_start_from_the_returned_root_never_fails) {
    std::mt19937_64 rng(20260731ULL);
    std::uniform_real_distribution<double> umaf(5e-4, 0.5);
    std::uniform_real_distribution<double> uresid(-3.0, 3.0);
    std::uniform_int_distribution<int> unout(1, 60);
    std::uniform_real_distribution<double> uvarScale(-3.0, 3.0);
    std::uniform_real_distribution<double> uz(2.0, 12.0);

    const spa::SolveOpts base;   // stated defaults, as production uses
    int nWarm = 0, nWarmFail = 0, nDrift = 0;

    for (int trial = 0; trial < 500; ++trial) {
        MixedCgf c;
        const int nOut = unout(rng);
        c.resid.resize(nOut);
        c.maf.resize(nOut);
        double sumR2 = 0.0;
        for (int i = 0; i < nOut; ++i) {
            c.resid[i] = uresid(rng);
            c.maf[i] = umaf(rng);
            sumR2 += c.resid[i] * c.resid[i];
        }
        c.var = std::pow(10.0, uvarScale(rng)) * sumR2;

        const spa::K012 at0 = c.kFull(0.0);
        const double sd = std::sqrt(at0.k2);

        for (int sgn = -1; sgn <= 1; sgn += 2) {
            c.s = at0.k1 + sgn * uz(rng) * sd;

            const spa::Saddle cold = spa::solveSaddlepoint(
                c.s, [&](double t) { return c.k12(t); },
                [&](double t) { return c.kFull(t); }, base);
            if (cold.status != spa::Status::SpaOk) continue;

            spa::SolveOpts w = base;
            w.init = cold.zeta;
            const spa::Saddle warm = spa::solveSaddlepoint(
                c.s, [&](double t) { return c.k12(t); },
                [&](double t) { return c.kFull(t); }, w);
            ++nWarm;
            if (!spa::statusIsUsable(warm.status)) ++nWarmFail;
            else if (!(std::fabs(warm.zeta - cold.zeta) <=
                       1e-6 * std::fmax(1.0, std::fabs(cold.zeta))))
                ++nDrift;
        }
    }

    CHECK(nWarm > 900);
    CHECK_MSG(nWarmFail == 0,
              "warm start failed on " + std::to_string(nWarmFail) + " of " +
                  std::to_string(nWarm) + " re-solves");
    CHECK_MSG(nDrift == 0,
              "warm start moved the root on " + std::to_string(nDrift) + " of " +
                  std::to_string(nWarm));
}

// ══════════════════════════════════════════════════════════════════════
// § 2  Property test over random problems
// ══════════════════════════════════════════════════════════════════════

TEST(property_residual_meets_stated_relative_tolerance) {
    std::mt19937_64 rng(20260730ULL);
    std::uniform_real_distribution<double> umaf(0.01, 0.5);
    std::uniform_real_distribution<double> uresid(-3.0, 3.0);
    std::uniform_int_distribution<int> unout(1, 40);
    std::uniform_real_distribution<double> uvar(1.0, 50.0);
    std::uniform_real_distribution<double> uz(2.0, 6.0);

    const spa::SolveOpts opt;  // stated defaults: rtol 1e-6, bracketRtol 4 eps
    int nConverged = 0;

    for (int trial = 0; trial < 400; ++trial) {
        MixedCgf c;
        const int nOut = unout(rng);
        c.resid.resize(nOut);
        c.maf.resize(nOut);
        for (int i = 0; i < nOut; ++i) {
            c.resid[i] = uresid(rng);
            c.maf[i] = umaf(rng);
        }
        c.var = uvar(rng);

        // K'(0) is the mean of S; place s a few standard deviations away, on
        // the scale the SPA branch is actually entered at.  A strictly
        // positive Gaussian block makes K' unbounded in both directions, so a
        // root exists for every draw.
        const spa::K012 at0 = c.kFull(0.0);
        const double sd = std::sqrt(at0.k2);
        const double sign = (trial % 2 == 0) ? 1.0 : -1.0;
        c.s = at0.k1 + sign * uz(rng) * sd;

        const spa::Saddle r = spa::solveSaddlepoint(
            c.s, [&](double t) { return c.k12(t); },
            [&](double t) { return c.kFull(t); }, opt);

        CHECK_MSG(r.status == spa::Status::SpaOk,
                  statusStr(r.status) + " trial " + std::to_string(trial));
        if (r.status != spa::Status::SpaOk) continue;
        ++nConverged;

        // (a) the achieved residual meets the criterion the solver documents.
        const double tol =
            std::fmax(opt.rtol * std::sqrt(r.K2), 4.0 * kEps * std::fabs(c.s));
        CHECK_MSG(std::fabs(r.residual) <= tol,
                  "residual " + std::to_string(r.residual) + " tol " +
                      std::to_string(tol));

        // (b) K0 and K2 are the cumulants AT zeta.  Re-evaluated here from an
        // independent call rather than trusted, which is the check that
        // catches a finder returning a post-step abscissa with a pre-step K''.
        const spa::K012 atRoot = c.kFull(r.zeta);
        CHECK_REL(r.K0, atRoot.k0, 0.0);
        CHECK_REL(r.K2, atRoot.k2, 0.0);
        CHECK_REL(r.residual, atRoot.k1, 0.0);

        // (c) the returned root lies inside the final bracket, and the bracket
        // straddles the root.
        CHECK(r.bracketed);
        CHECK(r.zeta >= r.lo && r.zeta <= r.hi);

        // (d) sgn(zeta) agrees with sgn(s - K'(0)), which is what makes
        // sgn(zeta) a usable proxy for sgn(w) in the tail kernel.
        CHECK_MSG(spa::signOf(r.zeta) == sign,
                  "zeta " + std::to_string(r.zeta) + " sign " + std::to_string(sign));

        // (e) K2 > 0 at the root: a genuine CGF is convex, and the tail kernel
        // depends on it.
        CHECK(r.K2 > 0.0);
    }
    CHECK(nConverged == 400);
}

// ══════════════════════════════════════════════════════════════════════
// § 2a  Scale equivariance
// ══════════════════════════════════════════════════════════════════════
//
// Multiplying every residual by a positive constant c is a pure change of
// units.  The statistic becomes c*S, so s becomes c*s, K(t) becomes K(c*t),
// K'(t) becomes c*K'(c*t) and K''(t) becomes c^2*K''(c*t); the saddlepoint
// equation K'(zeta) = s is solved by zeta/c; and
//
//     w^2 = 2*(zeta*s - K(zeta)),   v = zeta*sqrt(K''(zeta)),   r* = w + log(v/w)/w
//
// are each invariant.  The reported p-value therefore MUST NOT MOVE.  Nothing
// about this is approximate: it is an identity, and a solver that violates it
// is reading a threshold in the wrong units.
//
// This is the property the pre-repair header failed.  Its bracket-collapse test
// was `hi - lo <= stepTol * max(1, |t|)` with stepTol = 1e-8 absolute, and it
// reported Converged there with no residual test; since |zeta| ~ |z|/sqrt(Var S)
// shrinks as the residual scale grows, a large enough c drove the root below the
// tolerance and the solver accepted an abscissa it had never located.  The
// coarse bracket probe carried the same absolute `max(1, |init|)` floor.
// Measured on this population before the repair: 179 of 2000 problems moved P at
// c = 2^10, 121 of 2000 changed SPA_STATUS at c = 2^30, and the worst |dlog10 P|
// was 0.298 at c = 2^20 and infinite (a p-value replaced by NA) at c = 2^30.
//
// The check is BIT-IDENTITY, not closeness, and it is legitimate to demand it:
// for c a power of two the rescaling is exact in binary floating point — every
// product, quotient, sum and square root above scales exactly — so every
// comparison the solver makes has the same outcome at both scales and the whole
// computation is the same computation.  Anything less than bit-identity would
// mean some decision was still being made in the wrong units.
namespace {

struct ScaleOutcome {
    double p, negLog10p, zetaTimesScale;
    spa::Status status;
};

// The two-sided saddlepoint every method assembles, at a chosen residual scale.
ScaleOutcome twoSidedAtScale(const MixedCgf &base, double z, double scale) {
    MixedCgf c = base;
    for (double &r : c.resid) r *= scale;
    c.mean *= scale;
    c.var *= scale * scale;

    const spa::K012 at0 = c.kFull(0.0);
    const double sd = std::sqrt(at0.k2);

    double p[2], logp[2];
    spa::Status st[2];
    double zetaUpper = 0.0;

    for (int k = 0; k < 2; ++k) {
        const bool lowerTail = (k == 1);
        const double sign = lowerTail ? -1.0 : 1.0;
        c.s = at0.k1 + sign * z * sd;

        spa::SolveOpts opt;
        // The uncapped first-order saddlepoint |Score| / Var.  Every method in
        // the tree passes this capped at 1.2; the cap is a tier-3 constant in
        // the abscissa and is therefore not equivariant itself, so it is left
        // off here — what is under test is the solver, not the initializer.
        opt.init = sign * std::fabs(c.s - at0.k1) / at0.k2;

        const spa::Saddle r = spa::solveSaddlepoint(
            c.s, [&](double t) { return c.k12(t); },
            [&](double t) { return c.kFull(t); }, opt);
        if (k == 0) zetaUpper = r.zeta * scale;

        spa::Status stLin = spa::Status::SpaOk;
        spa::Status stLog = spa::Status::SpaOk;
        p[k] = spa::bnTail(r.zeta, c.s, r.K0, r.K2, lowerTail, stLin);
        logp[k] = spa::bnTailLog(r.zeta, c.s, r.K0, r.K2, lowerTail, stLog);
        st[k] = spa::worseStatus(r.status, spa::worseStatus(stLin, stLog));
    }

    const spa::TwoSided lin = spa::combineTails(p[0], p[1], st[0], st[1]);
    const spa::TwoSided lg = spa::combineTailsLog(logp[0], logp[1], st[0], st[1]);
    return ScaleOutcome{lin.p, lg.negLog10p, zetaUpper, lin.status};
}

bool bitIdentical(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return true;
    return std::memcmp(&a, &b, sizeof(double)) == 0;
}

// The production-shaped population used by § 2 and § 4a: a handful of binomial
// outlier terms plus the analytic Gaussian block for the non-outliers.
MixedCgf drawProductionShapedCgf(std::mt19937_64 &rng, double &zOut) {
    std::uniform_real_distribution<double> umaf(0.005, 0.5);
    std::uniform_real_distribution<double> uresid(-3.0, 3.0);
    std::uniform_int_distribution<int> unout(1, 60);
    std::uniform_real_distribution<double> udil(1.0, 200.0);
    std::uniform_real_distribution<double> uz(2.0, 8.0);

    MixedCgf c;
    const int nOut = unout(rng);
    c.resid.resize(nOut);
    c.maf.resize(nOut);
    double sumR2 = 0.0;
    for (int i = 0; i < nOut; ++i) {
        c.resid[i] = uresid(rng);
        c.maf[i] = umaf(rng);
        sumR2 += c.resid[i] * c.resid[i];
    }
    c.var = 2.0 * 0.25 * sumR2 * udil(rng);
    zOut = uz(rng);
    return c;
}

} // namespace

TEST(scale_equivariance_is_exact_under_a_power_of_two_rescaling) {
    // 2^40 is about 1.1e12, well past the residual scale any unstandardized
    // phenotype or --resid-name vector reaches, and short of the exponent range
    // where K'' itself would overflow.
    const double scales[] = {1024.0, 1048576.0, 1073741824.0, 1099511627776.0};
    const char *names[] = {"2^10", "2^20", "2^30", "2^40"};

    std::mt19937_64 rng(20260731ULL);
    int nCase = 0, nP = 0, nLog = 0, nZeta = 0, nStatus = 0;

    for (int trial = 0; trial < 600; ++trial) {
        double z = 0.0;
        const MixedCgf c = drawProductionShapedCgf(rng, z);
        const ScaleOutcome ref = twoSidedAtScale(c, z, 1.0);

        for (int k = 0; k < 4; ++k) {
            const ScaleOutcome got = twoSidedAtScale(c, z, scales[k]);
            ++nCase;
            if (!bitIdentical(got.p, ref.p)) {
                if (nP == 0)
                    CHECK_MSG(false, std::string("P moved at ") + names[k] +
                                         ": " + std::to_string(ref.p) + " -> " +
                                         std::to_string(got.p));
                ++nP;
            }
            if (!bitIdentical(got.negLog10p, ref.negLog10p)) ++nLog;
            if (!bitIdentical(got.zetaTimesScale, ref.zetaTimesScale)) ++nZeta;
            if (got.status != ref.status) ++nStatus;
        }
    }

    CHECK(nCase == 2400);
    CHECK_MSG(nP == 0, "P moved on " + std::to_string(nP) + " of " +
                           std::to_string(nCase) + " rescalings");
    CHECK_MSG(nLog == 0, "LOG10P moved on " + std::to_string(nLog) + " of " +
                             std::to_string(nCase));
    CHECK_MSG(nZeta == 0, "zeta*scale moved on " + std::to_string(nZeta) +
                              " of " + std::to_string(nCase));
    CHECK_MSG(nStatus == 0, "SPA_STATUS moved on " + std::to_string(nStatus) +
                                " of " + std::to_string(nCase));
}

// A decimal rescaling is not exact in binary, so the inputs themselves differ in
// their last bits and bit-identity is not available.  What must hold is the
// weaker but still binding statement: the answer agrees to within the tolerance
// the solver states.  rtol = 1e-6 on the residual pins zeta to ~rtol/|z| in
// relative terms, and -log10 P is stationary in zeta at the root (d(w^2)/dzeta =
// 2*(s - K'(zeta)) = 0 there), so the induced movement is far smaller again.
TEST(scale_equivariance_holds_under_a_decimal_rescaling) {
    const double scales[] = {1e3, 1e6, 1e8, 1e9};

    std::mt19937_64 rng(20260731ULL);
    double worst = 0.0;
    int nStatus = 0, nCase = 0;

    for (int trial = 0; trial < 600; ++trial) {
        double z = 0.0;
        const MixedCgf c = drawProductionShapedCgf(rng, z);
        const ScaleOutcome ref = twoSidedAtScale(c, z, 1.0);

        for (double sc : scales) {
            const ScaleOutcome got = twoSidedAtScale(c, z, sc);
            ++nCase;
            if (got.status != ref.status) ++nStatus;
            CHECK(std::isfinite(got.negLog10p) == std::isfinite(ref.negLog10p));
            if (std::isfinite(got.negLog10p) && std::isfinite(ref.negLog10p)) {
                const double d = std::fabs(got.negLog10p - ref.negLog10p);
                if (d > worst) worst = d;
            }
        }
    }

    std::fprintf(stderr,
                 "\n      decimal rescaling over %d cases: worst |dlog10P| "
                 "%.4g, status changes %d\n\n", nCase, worst, nStatus);
    CHECK_MSG(nStatus == 0,
              "SPA_STATUS moved on " + std::to_string(nStatus) + " cases");
    // Before the repair this reached 0.245 at 1e6 on the same population.
    CHECK_MSG(worst <= 1e-5, "worst |dlog10P| " + std::to_string(worst));
}

// The same population, solved through the single-callable overload, must give
// the same root — and must actually exercise the terminal reuse.
TEST(single_callable_overload_agrees_and_reuses_terminal) {
    std::mt19937_64 rng(31415926ULL);
    std::uniform_real_distribution<double> umaf(0.01, 0.5);
    std::uniform_real_distribution<double> uresid(-3.0, 3.0);
    std::uniform_real_distribution<double> uz(2.0, 6.0);

    int nReused = 0, nTrials = 0;
    for (int trial = 0; trial < 120; ++trial) {
        MixedCgf c;
        c.resid.resize(12);
        c.maf.resize(12);
        for (int i = 0; i < 12; ++i) {
            c.resid[i] = uresid(rng);
            c.maf[i] = umaf(rng);
        }
        c.var = 20.0;
        const spa::K012 at0 = c.kFull(0.0);
        c.s = at0.k1 + ((trial % 2) ? -1.0 : 1.0) * uz(rng) * std::sqrt(at0.k2);

        const spa::Saddle a = spa::solveSaddlepoint(
            c.s, [&](double t) { return c.k12(t); },
            [&](double t) { return c.kFull(t); }, spa::SolveOpts{});
        const spa::Saddle b = spa::solveSaddlepoint(
            c.s, [&](double t) { return c.kFull(t); }, spa::SolveOpts{});

        CHECK(a.status == b.status);
        CHECK_REL(a.zeta, b.zeta, 0.0);
        CHECK_REL(a.K0, b.K0, 0.0);
        CHECK_REL(a.K2, b.K2, 0.0);
        // The two-callable form cannot reuse: its cheap callable never
        // produced K0, so the terminal pass is unavoidable.
        CHECK(!a.reusedTerminal);

        if (b.status != spa::Status::SpaOk) continue;
        ++nTrials;
        if (b.reusedTerminal) ++nReused;
    }
    // Every converged solve ends at the abscissa of its last loop evaluation,
    // so the single-callable form reuses on every one of them.  This is where
    // spa_unify P3 is realized; the two-callable form pays one terminal pass
    // instead of the logarithm on every iteration (P2), which is the cheaper
    // trade for the binomial kernels.
    CHECK(nTrials > 100);
    CHECK_MSG(nReused == nTrials,
              "reused " + std::to_string(nReused) + " of " + std::to_string(nTrials));
}

// ══════════════════════════════════════════════════════════════════════
// § 3  Every Status branch
// ══════════════════════════════════════════════════════════════════════

TEST(status_names_and_classification) {
    CHECK(statusStr(spa::Status::SpaOk) == "SPA_OK");
    CHECK(statusStr(spa::Status::Normal) == "NORMAL");
    CHECK(statusStr(spa::Status::SpaWSingular) == "SPA_W_SINGULAR");
    CHECK(statusStr(spa::Status::FallbackMaxIter) == "FALLBACK_MAXITER");
    CHECK(statusStr(spa::Status::FallbackGuardTemp) == "FALLBACK_GUARD_TEMP");
    CHECK(statusStr(spa::Status::FallbackGuardCurv) == "FALLBACK_GUARD_CURV");
    CHECK(statusStr(spa::Status::FallbackNonFinite) == "FALLBACK_NONFINITE");
    CHECK(statusStr(spa::Status::NaPostFail) == "NA_POST_FAIL");
    CHECK(statusStr(spa::Status::NaNoTest) == "NA_NO_TEST");

    // The numeric encoding is the output column and is a contract, not an
    // implementation detail: SPA_STATUS <= 2 means LOG10P is trustworthy,
    // 3..6 that it is the substituted normal tail, >= 7 that it is NA.
    CHECK(static_cast<int>(spa::Status::SpaOk) == 0);
    CHECK(static_cast<int>(spa::Status::Normal) == 1);
    CHECK(static_cast<int>(spa::Status::SpaWSingular) == 2);
    CHECK(static_cast<int>(spa::Status::FallbackMaxIter) == 3);
    CHECK(static_cast<int>(spa::Status::FallbackGuardTemp) == 4);
    CHECK(static_cast<int>(spa::Status::FallbackGuardCurv) == 5);
    CHECK(static_cast<int>(spa::Status::FallbackNonFinite) == 6);
    CHECK(static_cast<int>(spa::Status::NaPostFail) == 7);
    CHECK(static_cast<int>(spa::Status::NaNoTest) == 8);

    const spa::Status all[9] = {
        spa::Status::SpaOk, spa::Status::Normal, spa::Status::SpaWSingular,
        spa::Status::FallbackMaxIter, spa::Status::FallbackGuardTemp,
        spa::Status::FallbackGuardCurv, spa::Status::FallbackNonFinite,
        spa::Status::NaPostFail, spa::Status::NaNoTest};

    // The three predicates partition the enumeration and agree with the
    // numeric blocks.  Exactly one of them holds for every value.
    for (int i = 0; i < 9; ++i) {
        const int u = spa::statusIsUsable(all[i]) ? 1 : 0;
        const int f = spa::statusIsFallback(all[i]) ? 1 : 0;
        const int n = spa::statusIsNA(all[i]) ? 1 : 0;
        CHECK(u + f + n == 1);
        CHECK(u == (i <= 2 ? 1 : 0));
        CHECK(f == (i >= 3 && i <= 6 ? 1 : 0));
        CHECK(n == (i >= 7 ? 1 : 0));
    }

    // SpaWSingular is a DEGRADED SUCCESS: below kWSingularity the r*
    // correction is not computable, the tail is evaluated as Phi(+/-w), and
    // LOG10P is a number rather than NA (spa_unify N3).
    CHECK(spa::statusIsUsable(spa::Status::SpaWSingular));

    // Severity is monotone in the encoding, with SpaOk and Normal tied so
    // that a genuine saddlepoint tail is not displaced by a normal label.
    CHECK(spa::statusSeverity(spa::Status::SpaOk) ==
          spa::statusSeverity(spa::Status::Normal));
    for (int i = 1; i < 9; ++i)
        CHECK(spa::statusSeverity(all[i]) >= spa::statusSeverity(all[i - 1]));
    for (int i = 2; i < 9; ++i)
        CHECK(spa::statusSeverity(all[i]) > spa::statusSeverity(all[i - 1]));

    // Every fallback and every NA outranks every usable value, so a tail that
    // merely degraded can never mask a tail that failed.
    for (int i = 3; i < 9; ++i) {
        CHECK(spa::statusSeverity(all[i]) >
              spa::statusSeverity(spa::Status::SpaWSingular));
        CHECK(spa::worseStatus(spa::Status::SpaWSingular, all[i]) == all[i]);
        CHECK(spa::worseStatus(all[i], spa::Status::SpaWSingular) == all[i]);
    }
    // ...but it outranks a plain success, so it is what gets reported.
    CHECK(spa::worseStatus(spa::Status::SpaOk, spa::Status::SpaWSingular) ==
          spa::Status::SpaWSingular);
    CHECK(spa::worseStatus(spa::Status::Normal, spa::Status::SpaWSingular) ==
          spa::Status::SpaWSingular);
}

// The D5 fallback rule, stated on spa::assemble itself.
TEST(assemble_fallback_and_na_blocks) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double z = 4.25;
    const double pNorm = spa::normalTwoSided(z);
    const double lNorm = -spa::normalTwoSidedLog(z) / std::log(10.0);

    // Both tails usable: the saddlepoint result passes through untouched.
    {
        const double lu = std::log(0.03), ll = std::log(0.02);
        const spa::TwoSided r = spa::assemble(0.03, 0.02, lu, ll,
                                              spa::Status::SpaOk,
                                              spa::Status::SpaOk, z);
        CHECK(r.status == spa::Status::SpaOk);
        CHECK(std::fabs(r.p - 0.05) < 1e-15);
        CHECK(std::fabs(r.negLog10p - (-std::log10(0.05))) < 1e-12);
    }

    // One failed tail falls the WHOLE two-sided p back, and to the normal
    // tail exactly — bit-for-bit 2*Phi(-|z|), not a mixture of the two.
    for (spa::Status f : {spa::Status::FallbackMaxIter,
                          spa::Status::FallbackGuardTemp,
                          spa::Status::FallbackGuardCurv,
                          spa::Status::FallbackNonFinite}) {
        const spa::TwoSided r = spa::assemble(0.03, nan, std::log(0.03), nan,
                                              spa::Status::SpaOk, f, z);
        CHECK(r.status == f);
        CHECK(r.p == pNorm);
        CHECK(std::fabs(r.negLog10p - lNorm) <= 0.0);
    }

    // The NA block never falls back, whatever z is.
    for (spa::Status n : {spa::Status::NaPostFail, spa::Status::NaNoTest}) {
        const spa::TwoSided r = spa::assemble(0.03, nan, std::log(0.03), nan,
                                              spa::Status::SpaOk, n, z);
        CHECK(r.status == n);
        CHECK(std::isnan(r.p));
        CHECK(std::isnan(r.negLog10p));
    }

    // A non-finite z is "no test exists", not a saddlepoint failure, and it
    // outranks a fallback status: there is nothing to fall back to.
    {
        const spa::TwoSided r = spa::assemble(0.03, nan, std::log(0.03), nan,
                                              spa::Status::SpaOk,
                                              spa::Status::FallbackGuardCurv, nan);
        CHECK(r.status == spa::Status::NaNoTest);
        CHECK(std::isnan(r.p));
        CHECK(std::isnan(r.negLog10p));
    }
}

TEST(status_converged) {
    const GaussianCgf c{0.0, 1.0, 3.0};
    const spa::Saddle r = spa::solveSaddlepoint(
        c.s, [&](double t) { return c.k12(t); },
        [&](double t) { return c.kFull(t); }, spa::SolveOpts{});
    CHECK(r.status == spa::Status::SpaOk);
    CHECK(statusStr(r.status) == "SPA_OK");
}

// Three independent routes to FallbackMaxIter, because they mean different things.
TEST(status_maxiter_budget_exhausted) {
    // Poisson with the root 13.8 away from the origin: the exponential K'
    // means Newton cannot arrive in two iterations however well safeguarded.
    const PoissonCgf c{1.0, 1e6};
    spa::SolveOpts o;
    o.maxIter = 2;
    const spa::Saddle r = spa::solveSaddlepoint(
        c.s, [&](double t) { return c.k12(t); },
        [&](double t) { return c.kFull(t); }, o);
    CHECK(r.status == spa::Status::FallbackMaxIter);
    CHECK(r.bracketed);
    CHECK(r.iters == 2);
    // Even on the failure path the cumulants belong to the returned abscissa.
    // This is the structural close of the stale-K'' hole: SPACox's Family A
    // finder advanced x on its last iteration and then reported the K''
    // computed before that step.
    const spa::K012 atRoot = c.kFull(r.zeta);
    CHECK_REL(r.K0, atRoot.k0, 0.0);
    CHECK_REL(r.K2, atRoot.k2, 0.0);
    CHECK_REL(r.residual, atRoot.k1, 0.0);

    // The same problem with the full budget converges, so the MaxIter above
    // really is a budget effect and not a broken bracket.
    const spa::Saddle ok = spa::solveSaddlepoint(
        c.s, [&](double t) { return c.k12(t); },
        [&](double t) { return c.kFull(t); }, spa::SolveOpts{});
    CHECK_MSG(ok.status == spa::Status::SpaOk,
              statusStr(ok.status) + " iters " + std::to_string(ok.iters));
    CHECK_CLOSE(ok.zeta, c.root(), 1e-6, 1e-6);
}

TEST(status_maxiter_no_bracket_exists) {
    // K'(t) = tanh(t) is bounded in (-1, 1), so K'(zeta) = 2 has no solution
    // and the residual never changes sign.  This is the real shape of the
    // failure for a purely-outlier binomial CGF whose s falls outside the
    // reachable range of K': there is no saddlepoint to find, and reporting a
    // number would be reporting a fiction.
    const double s = 2.0;
    // log cosh, written so it does not overflow at the far end of the
    // expansion sweep.
    auto logCosh = [](double t) {
        const double a = std::fabs(t);
        return a + std::log1p(std::exp(-2.0 * a)) - 0.69314718055994530942;
    };
    auto k12 = [&](double t) {
        const double th = std::tanh(t);
        return spa::K12{th - s, 1.0 - th * th};
    };
    auto kFull = [&](double t) {
        const double th = std::tanh(t);
        return spa::K012{logCosh(t), th - s, 1.0 - th * th};
    };
    const spa::Saddle r = spa::solveSaddlepoint(s, k12, kFull, spa::SolveOpts{});
    CHECK_MSG(r.status == spa::Status::FallbackMaxIter, statusStr(r.status));
    CHECK(!r.bracketed);
}

TEST(status_maxiter_sign_constraint_excludes_the_root) {
    // The Gaussian root is at +3, but the caller asserts sgn(zeta) < 0.  The
    // admissible half-line contains no sign change, so no probability is
    // produced -- which is the honest outcome, not a root on the wrong side.
    const GaussianCgf c{0.0, 1.0, 3.0};
    spa::SolveOpts o;
    o.scoreSign = -1.0;
    const spa::Saddle r = spa::solveSaddlepoint(
        c.s, [&](double t) { return c.k12(t); },
        [&](double t) { return c.kFull(t); }, o);
    CHECK_MSG(r.status == spa::Status::FallbackMaxIter, statusStr(r.status));
    CHECK(!r.bracketed);
    CHECK(r.zeta <= 0.0);
}

TEST(sign_constraint_confines_the_root_when_it_is_satisfiable) {
    const GaussianCgf pos{0.0, 1.0, 3.0};
    spa::SolveOpts o;
    o.scoreSign = 1.0;
    const spa::Saddle r = spa::solveSaddlepoint(
        pos.s, [&](double t) { return pos.k12(t); },
        [&](double t) { return pos.kFull(t); }, o);
    CHECK(r.status == spa::Status::SpaOk);
    CHECK(r.zeta > 0.0);
    CHECK_CLOSE(r.zeta, 3.0, 1e-6, 1e-6);
    CHECK(r.lo >= 0.0);

    const GaussianCgf neg{0.0, 1.0, -3.0};
    o.scoreSign = -1.0;
    const spa::Saddle r2 = spa::solveSaddlepoint(
        neg.s, [&](double t) { return neg.k12(t); },
        [&](double t) { return neg.kFull(t); }, o);
    CHECK(r2.status == spa::Status::SpaOk);
    CHECK(r2.zeta < 0.0);
    CHECK_CLOSE(r2.zeta, -3.0, 1e-6, 1e-6);
    // The bracket never left the admissible half-line.
    CHECK(r2.hi <= 0.0);
}

TEST(status_guardcurv_from_solver) {
    // A monotone residual with a persistently negative "K''".  The bracket
    // logic uses only the sign of the residual, so the root is still found by
    // pure bisection -- which is the point: the K'' > 0 test happens before
    // every division, so a bad curvature costs the quality of the step, not a
    // division by a negative number.  Five of the six original finders divided
    // by K'' unchecked.
    auto k12 = [](double t) { return spa::K12{t - 0.5, -1.0}; };
    auto kFull = [](double t) {
        return spa::K012{0.5 * t * t - 0.5 * t, t - 0.5, -1.0};
    };
    const spa::Saddle r = spa::solveSaddlepoint(0.0, k12, kFull, spa::SolveOpts{});
    CHECK_MSG(r.status == spa::Status::FallbackGuardCurv, statusStr(r.status));
    CHECK(r.bracketed);
    CHECK(r.K2 == -1.0);
    // Bisection still located the root despite Newton being disabled.
    CHECK_CLOSE(r.zeta, 0.5, 1e-6, 1e-6);
}

TEST(status_guardcurv_from_tail) {
    spa::Status st = spa::Status::SpaOk;
    CHECK(std::isnan(spa::bnTail(1.0, 1.0, 0.0, -1.0, false, st)));
    CHECK(st == spa::Status::FallbackGuardCurv);

    st = spa::Status::SpaOk;
    CHECK(std::isnan(spa::bnTail(1.0, 1.0, 0.0, 0.0, false, st)));
    CHECK(st == spa::Status::FallbackGuardCurv);

    st = spa::Status::SpaOk;
    CHECK(std::isnan(spa::bnTailLog(1.0, 1.0, 0.0, -1.0, false, st)));
    CHECK(st == spa::Status::FallbackGuardCurv);
}

TEST(status_guardtemp) {
    // zeta*s - K0 = 1*1 - 2 = -1 < 0, so w is not real.
    spa::Status st = spa::Status::SpaOk;
    CHECK(std::isnan(spa::bnTail(1.0, 1.0, 2.0, 1.0, false, st)));
    CHECK(st == spa::Status::FallbackGuardTemp);

    st = spa::Status::SpaOk;
    CHECK(std::isnan(spa::bnTailLog(1.0, 1.0, 2.0, 1.0, true, st)));
    CHECK(st == spa::Status::FallbackGuardTemp);

    // Also reached when zeta*s overflows negative.
    st = spa::Status::SpaOk;
    CHECK(std::isnan(spa::bnTail(-1e300, 1e300, 0.0, 1.0, false, st)));
    CHECK(st == spa::Status::FallbackGuardTemp);
}

TEST(status_guardw) {
    // The threshold is pinned here so that moving it is a deliberate act.  See
    // the derivation at spa::kWSingularity: it sits above the abscissa at which
    // the terminal K's absolute error makes r* catastrophic (2.8e-4 at
    // n = 2e5, 4.7e-4 at n = 1e6) and a decade below the |w| ~ 1e-2 that the
    // CLI-enforced minimum spaCutoff of 0.01 makes reachable.
    CHECK(spa::kWSingularity == 1e-3);

    // GuardW is a degraded SUCCESS: the tail is Phi(+/-w), not NaN.  w == 0
    // happens two ways -- temp == 0 exactly with zeta != 0, which is
    // spa_unify D6's surviving statement, and zeta == 0.  In SPAGRM and SPAsqr
    // the first was unguarded, giving v/w = +/-inf, r* = +/-inf and an emitted
    // p-value of exactly 0.0 or 1.0 depending on the tail.  Phi(0) = 0.5 in
    // both tails is the right answer: the statistic sits at the mean.
    spa::Status st = spa::Status::SpaOk;
    CHECK(spa::bnTail(2.0, 1.5, 3.0, 1.0, false, st) == 0.5);  // 2*1.5-3 == 0
    CHECK(st == spa::Status::SpaWSingular);

    st = spa::Status::SpaOk;
    CHECK(spa::bnTail(0.0, 1.0, 0.0, 1.0, false, st) == 0.5);   // zeta == 0
    CHECK(st == spa::Status::SpaWSingular);

    st = spa::Status::SpaOk;
    CHECK(spa::bnTail(-0.0, 1.0, 0.0, 1.0, true, st) == 0.5);   // zeta == -0
    CHECK(st == spa::Status::SpaWSingular);

    st = spa::Status::SpaOk;
    CHECK_REL(spa::bnTailLog(2.0, 1.5, 3.0, 1.0, true, st), -kLn2, 1e-15);
    CHECK(st == spa::Status::SpaWSingular);

    // Two-sided, that is a p of exactly 1 with a status that says why.
    spa::Status su = spa::Status::SpaOk, sl = spa::Status::SpaOk;
    const double pu = spa::bnTail(0.0, 1.0, 0.0, 1.0, false, su);
    const double pl = spa::bnTail(0.0, 1.0, 0.0, 1.0, true, sl);
    const spa::TwoSided ts = spa::combineTails(pu, pl, su, sl);
    CHECK(ts.status == spa::Status::SpaWSingular);
    CHECK(spa::statusIsUsable(ts.status));
    CHECK(ts.p == 1.0);
    CHECK(ts.negLog10p == 0.0);

    // Inside the band the returned probability is Phi(+/-w) — the same phiTail
    // call serves both branches, so the only difference from an exact match is
    // whether the compiler contracts `zeta*s - K0` into an FMA in one
    // translation unit and not the other.
    const double zeta = 0.5 * spa::kWSingularity;   // Gaussian CGF, sigma == 1
    const double s = zeta, K0 = 0.5 * zeta * zeta, K2 = 1.0;
    const double w = std::sqrt(2.0 * (zeta * s - K0));
    CHECK(w < spa::kWSingularity);
    for (bool lower : {false, true}) {
        st = spa::Status::SpaOk;
        CHECK_REL(spa::bnTail(zeta, s, K0, K2, lower, st),
                  math::pnorm(w, 0.0, 1.0, lower), 1e-14);
        CHECK(st == spa::Status::SpaWSingular);
    }

    // Immediately outside it the modified root is formed again.
    const double zOut = 2.0 * spa::kWSingularity;
    st = spa::Status::FallbackMaxIter;
    (void)spa::bnTail(zOut, zOut, 0.5 * zOut * zOut, 1.0, false, st);
    CHECK(st == spa::Status::SpaOk);
}

// The property the guard exists for, stated as a sensitivity measurement.
//
// A guard applied to r* AFTER computing it would be useless: dr*/dK ~ 1/w^3,
// so an absolute perturbation of K far below every budget the CGF kernels
// defend is already an O(1) perturbation of r* at the threshold, and no
// post-hoc inspection of the result can tell that apart from a real root.  The
// only defence is not to form r* at all, and the observable consequence is that
// inside the band the answer responds to dK like w does (dK/|w|) rather than
// like r* does (dK/|w|^3).
TEST(guardw_fires_before_r_star_is_formed) {
    // Gaussian CGF with sigma == 1: s = zeta, K0 = zeta^2/2, K2 = 1, so
    // w = zeta exactly and v/w = 1, i.e. the removable singularity itself.
    const double K2 = 1.0;
    const double dK = 1e-11;   // the terminal K's absolute error at n ~ 1e6

    auto rStarTailAsWritten = [&](double zeta, double s, double k0, bool lower) {
        const double w = spa::signOf(zeta) * std::sqrt(2.0 * (zeta * s - k0));
        const double v = zeta * std::sqrt(K2);
        return math::pnorm(w + std::log(v / w) / w, 0.0, 1.0, lower);
    };

    {   // Inside the band.
        const double zeta = 0.5 * spa::kWSingularity;   // w = 5e-4
        const double s = zeta, K0 = 0.5 * zeta * zeta;

        spa::Status s0 = spa::Status::SpaOk, s1 = spa::Status::SpaOk;
        const double p0 = spa::bnTail(zeta, s, K0, K2, false, s0);
        const double p1 = spa::bnTail(zeta, s, K0 + dK, K2, false, s1);
        CHECK(s0 == spa::Status::SpaWSingular);
        CHECK(s1 == spa::Status::SpaWSingular);

        // Through w alone: |dp| ~ phi(0)*dK/|w| = 0.4*1e-11/5e-4 = 8e-9.
        const double moved = std::fabs(p1 - p0);
        CHECK_MSG(moved < 1e-7, "guarded path moved " + std::to_string(moved));

        // What r* would have done with the identical perturbation:
        // dr* = dK/|w|^3 = 1e-11/1.25e-10 = 0.08, so |dp| ~ 0.03.
        const double q0 = rStarTailAsWritten(zeta, s, K0, false);
        const double q1 = rStarTailAsWritten(zeta, s, K0 + dK, false);
        const double wouldMove = std::fabs(q1 - q0);
        CHECK_MSG(wouldMove > 1e-3,
                  "r* path moved only " + std::to_string(wouldMove));
        // Five orders of magnitude between the two responses is the whole
        // content of the guard.
        CHECK(wouldMove > 1e4 * moved);
    }

    {   // Outside the band, where r* IS formed: the residual sensitivity is
        // real but bounded, and this records how much of it survives.
        const double zeta = 4.0 * spa::kWSingularity;   // w = 4e-3
        const double s = zeta, K0 = 0.5 * zeta * zeta;
        spa::Status s0 = spa::Status::FallbackMaxIter, s1 = spa::Status::FallbackMaxIter;
        const double p0 = spa::bnTail(zeta, s, K0, K2, false, s0);
        const double p1 = spa::bnTail(zeta, s, K0 + dK, K2, false, s1);
        CHECK(s0 == spa::Status::SpaOk);
        CHECK(s1 == spa::Status::SpaOk);
        // dK/|w|^3 = 1e-11/6.4e-8 = 1.6e-4, so |dp| ~ 6e-5, and the two-sided
        // p here is 1 to within 1e-2, hence |dlog10 P| ~ 3e-5.
        CHECK(std::fabs(p1 - p0) < 1e-3);
    }
}

// 01_findings.md N3 states that at small |w| "the dropped term is O(zeta)".
// It is not: expanding about the mean gives r* -> w + rho3/6 with
// rho3 = K'''(0)/K''(0)^(3/2), a non-zero constant.  Poisson(lambda) has
// rho3 = 1/sqrt(lambda) and is the cleanest witness.  The distinction matters
// because it is what bounds the cost of the Phi(+/-w) degradation: O(rho3),
// not O(zeta), and therefore not something that vanishes at the threshold.
TEST(dropped_correction_tends_to_the_skewness_constant_not_to_zero) {
    for (double lambda : {0.5, 1.0, 4.0, 100.0}) {
        const double want = 1.0 / (6.0 * std::sqrt(lambda));
        double last = 0.0;
        for (double zeta : {1e-2, 1e-3, 1e-4}) {
            // Poisson CGF: K = lambda*(e^t - 1), K' = K'' = lambda*e^t.
            const double s = lambda * std::exp(zeta);
            const double K0 = lambda * std::expm1(zeta);
            const double K2 = s;
            const double w = std::sqrt(2.0 * (zeta * s - K0));
            const double v = zeta * std::sqrt(K2);
            last = std::log(v / w) / w;
        }
        // Converged to rho3/6 rather than to zero.  The residual at
        // zeta = 1e-4 is O(zeta) and measures 0.5% of the limit.
        CHECK_MSG(std::fabs(last - want) < 0.02 * want,
                  "lambda=" + std::to_string(lambda) + " limit=" +
                      std::to_string(last) + " want=" + std::to_string(want));
        CHECK(std::fabs(last) > 0.5 * want);   // emphatically not O(zeta)
    }
}

TEST(status_nonfinite_from_tail) {
    const double bad[] = {kNaN, kInf, -kInf};
    for (double b : bad) {
        spa::Status st = spa::Status::SpaOk;
        CHECK(std::isnan(spa::bnTail(b, 1.0, 0.0, 1.0, false, st)));
        CHECK(st == spa::Status::FallbackNonFinite);

        st = spa::Status::SpaOk;
        CHECK(std::isnan(spa::bnTail(1.0, b, 0.0, 1.0, false, st)));
        CHECK(st == spa::Status::FallbackNonFinite);

        st = spa::Status::SpaOk;
        CHECK(std::isnan(spa::bnTail(1.0, 1.0, b, 1.0, false, st)));
        CHECK(st == spa::Status::FallbackNonFinite);

        st = spa::Status::SpaOk;
        CHECK(std::isnan(spa::bnTail(1.0, 1.0, 0.0, b, false, st)));
        CHECK(st == spa::Status::FallbackNonFinite);

        st = spa::Status::SpaOk;
        CHECK(std::isnan(spa::bnTailLog(b, 1.0, 0.0, 1.0, false, st)));
        CHECK(st == spa::Status::FallbackNonFinite);
    }
}

TEST(status_nonfinite_from_solver_everywhere) {
    // A CGF that is non-finite at every abscissa: the retreat toward the
    // origin finds nothing usable, so no bracket is even attempted.
    auto k12 = [](double) { return spa::K12{kNaN, kNaN}; };
    auto kFull = [](double) { return spa::K012{kNaN, kNaN, kNaN}; };
    spa::SolveOpts o;
    o.init = 1.0;
    const spa::Saddle r = spa::solveSaddlepoint(0.0, k12, kFull, o);
    CHECK_MSG(r.status == spa::Status::FallbackNonFinite, statusStr(r.status));
    CHECK(!r.bracketed);
    CHECK(std::isnan(r.K2));
}

TEST(status_nonfinite_from_solver_interior) {
    // Finite and monotone at the bracket endpoints, non-finite in the middle,
    // where both the Newton step and the bisection retry land.
    //
    // The width of the non-finite interval is chosen so that the expansion
    // schedule brackets around it: from init = 0 the Newton probe at 0.5 is
    // non-finite, the halving retreat lands at 0.25, and the coarse jump of
    // max(|t|, |Newton step|) = 0.5 from there reaches 0.75, giving the bracket
    // [0.25, 0.75] with the hole strictly inside it.  A hole WIDER than the
    // coarse distance cannot be stepped over at all — the halving retreat
    // creeps up to its edge, as the comment at that retreat says — and the
    // solver then reports MaxIter for want of a sign change rather than
    // NonFinite.  That is a property of the expansion, not of this status, and
    // no CGF in GRAB can produce it: spa_cgf clamps t*r, so its cumulants are
    // finite at every finite abscissa and the non-finite region, where one
    // exists at all, is a half-line rather than an interval.
    auto res = [](double t) {
        if (t > 0.3 && t < 0.7) return kNaN;
        return t - 0.5;
    };
    auto k12 = [&](double t) { return spa::K12{res(t), 1.0}; };
    auto kFull = [&](double t) {
        return spa::K012{0.5 * t * t - 0.5 * t, res(t), 1.0};
    };
    const spa::Saddle r = spa::solveSaddlepoint(0.0, k12, kFull, spa::SolveOpts{});
    CHECK_MSG(r.status == spa::Status::FallbackNonFinite, statusStr(r.status));
    CHECK(r.bracketed);
}

TEST(status_normalbranch) {
    const spa::TwoSided b = spa::normalBranch(1.5);
    CHECK(b.status == spa::Status::Normal);
    CHECK(statusStr(b.status) == "NORMAL");
    CHECK(spa::statusIsUsable(b.status));

    // It survives the two-sided assembly, so a method that took the normal
    // branch on both sides still reports NORMAL rather than OK.
    const spa::TwoSided both =
        spa::combineTails(0.05, 0.05, spa::Status::Normal,
                          spa::Status::Normal);
    CHECK(both.status == spa::Status::Normal);
    const spa::TwoSided bothLog =
        spa::combineTailsLog(-3.0, -3.0, spa::Status::Normal,
                             spa::Status::Normal);
    CHECK(bothLog.status == spa::Status::Normal);
}

// ══════════════════════════════════════════════════════════════════════
// § 4  The bracket invariant
// ══════════════════════════════════════════════════════════════════════

TEST(bracket_invariant_holds_at_every_iteration) {
    struct Case { const char *name; int kind; double a, b, c; };
    // kind 0 = Gaussian(mu=a, var=b, s=c), 1 = Poisson(lambda=a, s=c),
    // 2 = BinomOne(p=a, r=b, s=c)
    const Case cases[] = {
        {"gauss", 0, 0.0, 4.0, 1.4},
        {"gauss-neg", 0, 0.5, 0.4, -0.3},
        {"poisson-far", 1, 1.0, 0.0, 1e6},
        {"poisson-near", 1, 4.0, 0.0, 1.0},
        {"poisson-tiny", 1, 1000.0, 0.0, 1.0},
        {"binom", 2, 0.3, 1.0, 0.9},
        {"binom-neg", 2, 0.05, -2.0, -3.6},
    };

    int totalRows = 0;
    for (const Case &cs : cases) {
        Recorder rec;
        spa::Saddle r{};
        double trueRoot = kNaN;

        if (cs.kind == 0) {
            const GaussianCgf g{cs.a, cs.b, cs.c};
            trueRoot = g.root();
            r = spa::solveSaddlepointTraced(
                g.s, [&](double t) { return g.k12(t); },
                [&](double t) { return g.kFull(t); }, tightOpts(), rec);
        } else if (cs.kind == 1) {
            const PoissonCgf pz{cs.a, cs.c};
            trueRoot = pz.root();
            r = spa::solveSaddlepointTraced(
                pz.s, [&](double t) { return pz.k12(t); },
                [&](double t) { return pz.kFull(t); }, tightOpts(), rec);
        } else {
            const BinomOneCgf bn{cs.a, cs.b, cs.c};
            trueRoot = bn.root();
            r = spa::solveSaddlepointTraced(
                bn.s, [&](double t) { return bn.k12(t); },
                [&](double t) { return bn.kFull(t); }, tightOpts(), rec);
        }

        CHECK_MSG(r.status == spa::Status::SpaOk,
                  std::string(cs.name) + ": " + statusStr(r.status));
        // The trace fires once per loop iteration, so the recorded count is
        // the reported count.
        CHECK_MSG(static_cast<int>(rec.rows.size()) == r.iters,
                  std::string(cs.name) + ": rows " +
                      std::to_string(rec.rows.size()) + " iters " +
                      std::to_string(r.iters));
        totalRows += static_cast<int>(rec.rows.size());

        // The bracket is compared against the closed-form root with a few-ULP
        // slack, because the two are not the same number.  The solver brackets
        // a root of the *evaluated* residual, and once the bracket has
        // collapsed to adjacent doubles -- which tightOpts drives it to -- the
        // residual is exactly zero at more than one representable abscissa.
        // In the "binom" case below the bracket collapses to a single double
        // that sits 0.78 ULP from the closed form, and the residual is exactly
        // zero at both.  Demanding exact containment would be demanding that
        // the closed-form expression and the CGF round identically.
        const double slack = 8.0 * kEps * std::fmax(1.0, std::fabs(trueRoot));

        double prevWidth = kInf;
        for (const Recorder::Row &row : rec.rows) {
            const std::string where =
                std::string(cs.name) + " iter " + std::to_string(row.iter);
            // (a) the iterate is inside its own bracket
            CHECK_MSG(row.t >= row.lo && row.t <= row.hi, where + ": t outside [lo, hi]");
            // (b) the true root is inside the bracket
            CHECK_MSG(trueRoot >= row.lo - slack && trueRoot <= row.hi + slack,
                      where + ": root " + std::to_string(trueRoot) + " outside [" +
                          std::to_string(row.lo) + ", " + std::to_string(row.hi) + "]");
            // (c) the bracket never widens
            const double width = row.hi - row.lo;
            CHECK_MSG(width <= prevWidth, where + ": bracket widened");
            prevWidth = width;
        }
        // (d) the final bracket still straddles the root, and the returned
        //     root is inside it.
        CHECK(r.bracketed);
        CHECK_MSG(trueRoot >= r.lo - slack && trueRoot <= r.hi + slack,
                  std::string(cs.name) + ": final bracket lost the root");
        CHECK(r.zeta >= r.lo && r.zeta <= r.hi);
    }
    // The sweep has to do real work somewhere, or the invariant above is
    // vacuously true.
    CHECK_MSG(totalRows > 30, "only " + std::to_string(totalRows) + " iterations traced");
}

TEST(bisection_recovers_when_newton_is_sabotaged) {
    // A correct monotone residual paired with a K'' a thousand times too
    // small, so every Newton step overshoots the bracket by a factor of a
    // thousand and must be rejected.  Convergence has to come from bisection.
    auto k12 = [](double t) { return spa::K12{t - 0.375, 1e-3}; };
    auto kFull = [](double t) {
        return spa::K012{0.5 * t * t - 0.375 * t, t - 0.375, 1e-3};
    };
    Recorder rec;
    const spa::Saddle r =
        spa::solveSaddlepointTraced(0.0, k12, kFull, tightOpts(), rec);
    CHECK_MSG(r.status == spa::Status::SpaOk, statusStr(r.status));
    CHECK_CLOSE(r.zeta, 0.375, 1e-14, 1e-15);
    CHECK(rec.rows.size() > 20);  // it really did bisect its way there
    for (const Recorder::Row &row : rec.rows) {
        CHECK(row.t >= row.lo && row.t <= row.hi);
        CHECK(0.375 >= row.lo && 0.375 <= row.hi);
    }
}

TEST(forced_bisection_bounds_the_iteration_count) {
    // Newton on an exponential K' creeps toward the root by a fixed decrement
    // per step while reducing |residual| every time, so neither of the two
    // bisection triggers in the design table fires.  Without the
    // bracket-stalled safeguard this needs ~475 iterations; with it, well
    // under the default budget of 100.
    const PoissonCgf c{1.0, 1e6};
    const spa::Saddle r = spa::solveSaddlepoint(
        c.s, [&](double t) { return c.k12(t); },
        [&](double t) { return c.kFull(t); }, spa::SolveOpts{});
    CHECK_MSG(r.status == spa::Status::SpaOk,
              statusStr(r.status) + " iters " + std::to_string(r.iters));
    CHECK_MSG(r.iters < 60, "iters " + std::to_string(r.iters));
    CHECK_CLOSE(r.zeta, c.root(), 1e-6, 1e-6);
}

// ══════════════════════════════════════════════════════════════════════
// § 4a  The evaluation budget
// ══════════════════════════════════════════════════════════════════════
//
// The CGF evaluation is the entire cost of the saddlepoint branch: one call is
// an O(nOutlier) SIMD reduction, and the solver makes several per tail per
// marker per phenotype.  A schedule change that keeps every guarantee and every
// answer while doubling the number of evaluations is therefore a serious
// regression that no other test in this file would notice — the roots, the
// statuses and the brackets would all still be right.
//
// The counts below are pinned against the population the SPA branch actually
// runs on, and they are pinned as MEANS with a ceiling rather than exactly,
// because the individual counts depend on the last bits of the residual.  The
// two phases are separated: `nEvalBracket` measures how many evaluations went
// into locating the bracket from the caller's initial guess, and
// `nEval - nEvalBracket` how many the safeguarded-Newton loop needed.  The
// terminal `evalFull` is in neither, being one call per solve by construction.
//
// Recorded history, so a future change can see which direction it moved.  The
// production rows were measured by instrumenting the solver and running the two
// migrated methods over the 50000-subject x 20000-marker synthetic cohort; this
// test's own population is deliberately broader (nOutlier 1 to 60, MAF 0.005 to
// 0.5, the Gaussian block scaled over two orders of magnitude), so its loop
// count is higher than production's and the two are not directly comparable.
//
//     schedule                       population    bracket   loop   total
//     ---------------------------    -----------   -------   ----   -----
//     Stage 1: first probe at        SPAsqr           2.00  12.71   14.71
//     max(coarse, Newton),           SPAGRM           2.00   9.98   11.98
//     forced bisection when the
//     BRACKET stopped halving
//
//     Stage 4 rework: first probe    SPAsqr           2.78   0.88    3.66
//     at Newton, forced bisection    SPAGRM           2.91   1.03    3.94
//     when the STEP stopped          this test        2.88   1.77    4.65
//     contracting
//
// The ceilings below are set about 25 % above this test's measured means: loose
// enough not to fire on an arithmetic reordering, tight enough that a return to
// the Stage 1 schedule fails loudly, since that schedule spends ten or more loop
// evaluations against a ceiling of two.
TEST(evaluation_budget_per_tail_is_pinned) {
    std::mt19937_64 rng(20260731ULL);
    std::uniform_real_distribution<double> umaf(0.005, 0.5);
    std::uniform_real_distribution<double> uresid(-3.0, 3.0);
    std::uniform_int_distribution<int> unout(1, 60);
    std::uniform_real_distribution<double> udil(1.0, 200.0);
    std::uniform_real_distribution<double> uz(2.0, 8.0);

    long long totEval = 0, totBrk = 0, totIter = 0;
    int nSolve = 0, worstEval = 0;

    for (int trial = 0; trial < 600; ++trial) {
        MixedCgf c;
        const int nOut = unout(rng);
        c.resid.resize(nOut);
        c.maf.resize(nOut);
        double sumR2 = 0.0;
        for (int i = 0; i < nOut; ++i) {
            c.resid[i] = uresid(rng);
            c.maf[i] = umaf(rng);
            sumR2 += c.resid[i] * c.resid[i];
        }
        // The analytic non-outlier Gaussian block every production site adds.
        c.var = 2.0 * 0.25 * sumR2 * udil(rng);

        const spa::K012 at0 = c.kFull(0.0);
        const double sd = std::sqrt(at0.k2);

        for (int side = 0; side < 2; ++side) {
            const double sign = side ? -1.0 : 1.0;
            MixedCgf cc = c;
            cc.s = at0.k1 + sign * uz(rng) * sd;

            // Exactly what SPAsqr and SPAGRM pass: the first-order saddlepoint
            // estimate s/K''(0), capped at 1.2, signed to match the score, and
            // the sign constraint asserted.
            spa::SolveOpts opt;
            opt.init = sign * std::fmin(std::fabs(cc.s - at0.k1) / at0.k2, 1.2);
            opt.scoreSign = sign;

            const spa::Saddle r = spa::solveSaddlepoint(
                cc.s, [&](double t) { return cc.k12(t); },
                [&](double t) { return cc.kFull(t); }, opt);

            CHECK_MSG(r.status == spa::Status::SpaOk,
                      statusStr(r.status) + " trial " + std::to_string(trial));
            if (r.status != spa::Status::SpaOk) continue;

            // The split is a partition of nEval, and the bracket phase always
            // costs at least the one evaluation at `init`.
            CHECK(r.nEvalBracket >= 1);
            CHECK(r.nEvalBracket <= r.nEval);

            ++nSolve;
            totEval += r.nEval;
            totBrk  += r.nEvalBracket;
            totIter += r.iters;
            if (r.nEval > worstEval) worstEval = r.nEval;
        }
    }

    CHECK(nSolve == 1200);
    const double mEval = static_cast<double>(totEval) / nSolve;
    const double mBrk  = static_cast<double>(totBrk) / nSolve;
    const double mLoop = mEval - mBrk;

    std::fprintf(stderr,
                 "\n      evaluations per tail over %d solves: bracket %.3f, "
                 "loop %.3f, total %.3f (worst single solve %d)\n"
                 "      terminal evalFull: 1 per solve by construction\n\n",
                 nSolve, mBrk, mLoop, mEval, worstEval);

    CHECK_MSG(mBrk <= 3.4, "bracket evaluations per tail " + std::to_string(mBrk));
    CHECK_MSG(mLoop <= 2.0, "loop evaluations per tail " + std::to_string(mLoop));
    CHECK_MSG(mEval <= 5.0, "total evaluations per tail " + std::to_string(mEval));
    // No single solve may blow up: the safeguard bounds the iteration count.
    CHECK_MSG(worstEval <= 40, "worst single solve " + std::to_string(worstEval));
}

// ══════════════════════════════════════════════════════════════════════
// § 5  bnTail and bnTailLog
// ══════════════════════════════════════════════════════════════════════

// Construct (zeta, s, K0, K2) producing a chosen w, so the tail can be probed
// anywhere in the range including far past the underflow floor.  For a
// Gaussian CGF with mu = 0: s = var*zeta and K0 = var*zeta^2/2, hence
//   temp = zeta*s - K0 = var*zeta^2/2,
//   w    = sgn(zeta)*sqrt(var)*|zeta| = sqrt(var)*zeta,
//   v    = zeta*sqrt(var) = w,
// so log(v/w) = 0 and r* = w exactly.  The saddlepoint tail of a Gaussian is
// the normal tail, which is the strongest available pin on the sign
// conventions.
struct GaussTailPoint {
    double zeta, s, K0, K2;
};
static GaussTailPoint gaussTailPoint(double w, double var) {
    const double zeta = w / std::sqrt(var);
    return GaussTailPoint{zeta, var * zeta, 0.5 * var * zeta * zeta, var};
}

TEST(bntail_gaussian_case_reduces_to_phi) {
    for (double w = 0.25; w <= 10.0; w += 0.25) {
        const GaussTailPoint g = gaussTailPoint(w, 2.5);
        spa::Status st = spa::Status::FallbackMaxIter;

        // temp and v are two algebraically-equal but numerically distinct
        // routes to the same quantity, so v/w differs from 1 by a few ULP and
        // r* differs from w by ~eps/w.  d(log p)/dw = -w, so the induced
        // relative error in p grows like w.  The bound admits that and
        // nothing more.
        const double tol = 1e-15 * (1.0 + w * w);

        const double up = spa::bnTail(g.zeta, g.s, g.K0, g.K2, false, st);
        CHECK(st == spa::Status::SpaOk);
        CHECK_REL(up, math::pnorm(w, 0.0, 1.0, false), tol);

        const GaussTailPoint gn = gaussTailPoint(-w, 2.5);
        const double lo = spa::bnTail(gn.zeta, gn.s, gn.K0, gn.K2, true, st);
        CHECK(st == spa::Status::SpaOk);
        CHECK_REL(lo, math::pnorm(-w, 0.0, 1.0, true), tol);

        // The upper tail at +w and the lower tail at -w are the same number.
        CHECK_REL(up, lo, 1e-15);
    }
}

TEST(bntaillog_agrees_with_log_bntail) {
    int nCompared = 0;
    for (double w = 0.5; w <= 37.0; w += 0.125) {
        const GaussTailPoint g = gaussTailPoint(w, 1.75);
        spa::Status stL = spa::Status::FallbackMaxIter, stP = spa::Status::FallbackMaxIter;
        const double p = spa::bnTail(g.zeta, g.s, g.K0, g.K2, false, stP);
        const double lp = spa::bnTailLog(g.zeta, g.s, g.K0, g.K2, false, stL);
        CHECK(stP == spa::Status::SpaOk);
        CHECK(stL == spa::Status::SpaOk);
        if (p > 0.0 && std::isfinite(std::log(p))) {
            CHECK_MSG(std::fabs(lp - std::log(p)) <=
                          1e-13 * std::fmax(1.0, std::fabs(std::log(p))),
                      "w=" + std::to_string(w) + " log(p)=" +
                          std::to_string(std::log(p)) + " lp=" + std::to_string(lp));
            ++nCompared;
        }
        // The lower tail must behave identically under the sign flip.
        const GaussTailPoint gn = gaussTailPoint(-w, 1.75);
        const double p2 = spa::bnTail(gn.zeta, gn.s, gn.K0, gn.K2, true, stP);
        const double lp2 = spa::bnTailLog(gn.zeta, gn.s, gn.K0, gn.K2, true, stL);
        CHECK_REL(p2, p, 1e-15);
        CHECK_NEAR(lp2, lp, 1e-14 * std::fmax(1.0, std::fabs(lp)));
    }
    CHECK_MSG(nCompared > 250, "only " + std::to_string(nCompared) + " comparisons");
}

TEST(bntaillog_survives_past_linear_underflow) {
    // Phi(-x) becomes denormal near x = 38.0 and flushes to zero near 38.5.
    // The log-domain tail must keep going: an association with true
    // p = 1e-400 has to remain reportable in LOG10P.
    bool sawLinearZero = false;
    for (double w = 36.0; w <= 200.0; w += 2.0) {
        const GaussTailPoint g = gaussTailPoint(w, 1.0);
        spa::Status st = spa::Status::FallbackMaxIter;
        const double p = spa::bnTail(g.zeta, g.s, g.K0, g.K2, false, st);
        const double lp = spa::bnTailLog(g.zeta, g.s, g.K0, g.K2, false, st);
        CHECK(st == spa::Status::SpaOk);
        CHECK_MSG(std::isfinite(lp), "bnTailLog non-finite at w=" + std::to_string(w));
        // log Phi(-w) = -w^2/2 - log(w*sqrt(2*pi)) + O(1/w^2); check that the
        // leading term is present and that nothing else dominates it.
        CHECK_MSG(lp < -0.49 * w * w && lp > -0.51 * w * w - 10.0,
                  "w=" + std::to_string(w) + " lp=" + std::to_string(lp));
        if (p == 0.0) sawLinearZero = true;
    }
    CHECK_MSG(sawLinearZero, "linear bnTail never underflowed to zero");

    // p = 1e-400 has log p = -920.9.  Invert log Phi(-w) = that by bisection
    // on w, then confirm the log tail reproduces it while the linear tail has
    // nothing left at all.
    const double targetLog = -400.0 * kLn10;
    double blo = 1.0, bhi = 200.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (blo + bhi);
        const GaussTailPoint g = gaussTailPoint(mid, 1.0);
        spa::Status st;
        const double lp = spa::bnTailLog(g.zeta, g.s, g.K0, g.K2, false, st);
        if (lp > targetLog) blo = mid; else bhi = mid;
    }
    const GaussTailPoint gg = gaussTailPoint(0.5 * (blo + bhi), 1.0);
    spa::Status st;
    const double lp = spa::bnTailLog(gg.zeta, gg.s, gg.K0, gg.K2, false, st);
    CHECK_REL(lp, targetLog, 1e-12);
    CHECK(spa::bnTail(gg.zeta, gg.s, gg.K0, gg.K2, false, st) == 0.0);
    // Expressed the way the LOG10P column will carry it.
    CHECK_CLOSE(-lp / kLn10, 400.0, 1e-12, 0.0);
}

TEST(bntaillog_matches_an_independent_higher_precision_reference) {
    // Sweeps across x = -37, where bnTailLog switches from log(Phi) to the
    // Mills-ratio asymptotic expansion.  Both branches are checked against the
    // same expansion evaluated in long double, which pins each of them
    // absolutely and therefore pins the continuity of the pair.  A missing
    // log(sqrt(2*pi)) would show as 0.92, a missing series as 7e-4, a
    // sign-flipped series as 1.5e-3 -- all far above the bound.
    int nBelow = 0, nAbove = 0;
    for (double w = 34.0; w <= 40.0; w += 0.05) {
        const GaussTailPoint g = gaussTailPoint(w, 1.0);
        spa::Status st;
        const double lp = spa::bnTailLog(g.zeta, g.s, g.K0, g.K2, false, st);
        CHECK(st == spa::Status::SpaOk);

        const long double x = static_cast<long double>(w);
        const long double y = 1.0L / (x * x);
        const long double series =
            y * (-1.0L + y * (3.0L + y * (-15.0L + y * (105.0L +
                 y * (-945.0L + y * 10395.0L)))));
        const long double ref = -0.5L * x * x - std::log(x) -
                                0.9189385332046727417803297364056L +
                                std::log1p(series);
        CHECK_MSG(std::fabs(lp - static_cast<double>(ref)) <= 1e-10,
                  "w=" + std::to_string(w) + " lp=" + std::to_string(lp) +
                      " ref=" + std::to_string(static_cast<double>(ref)));
        if (w < 37.0) ++nBelow; else ++nAbove;
    }
    CHECK(nBelow > 30);
    CHECK(nAbove > 30);
}

// ══════════════════════════════════════════════════════════════════════
// § 6  bnTail treats K0 and K2 as opaque scalars
// ══════════════════════════════════════════════════════════════════════

TEST(bntail_does_not_reconcile_an_inconsistent_k0_k2_pair) {
    // WtCoxG builds k0_total from var_n_K01 and k2_total from var_n_K2, which
    // are different variances (wtcoxg.cpp:430-435 documents this as
    // intentional; spa_unify D4 rules it out of scope pending a maintainer
    // decision).  The kernel must consume the pair exactly as handed over: a
    // kernel that "corrected" K2 to the second derivative of K0 would silently
    // change WtCoxG's statistic rather than compute it.
    const double zeta = 0.8, s = 4.0, K0 = 1.2;
    for (double K2 : {0.5, 1.0, 2.0, 7.5, 1e4}) {
        spa::Status st = spa::Status::FallbackMaxIter;
        const double got = spa::bnTail(zeta, s, K0, K2, false, st);
        CHECK(st == spa::Status::SpaOk);
        const double temp = zeta * s - K0;
        const double w = std::sqrt(2.0 * temp);   // zeta > 0
        const double v = zeta * std::sqrt(K2);
        const double rStar = w + std::log(v / w) / w;
        CHECK_REL(got, math::pnorm(rStar, 0.0, 1.0, false), 0.0);
    }
    // And the answers really are distinct, so nothing was normalized away.
    spa::Status st;
    CHECK(spa::bnTail(zeta, s, K0, 0.5, false, st) !=
          spa::bnTail(zeta, s, K0, 7.5, false, st));
}

TEST(bntail_matches_the_production_kernel_on_random_inputs) {
    // Reproduction of the arithmetic at wtcoxg.cpp:502-514 and
    // spacox.cpp:385-400, so the migrations in Stages 3 to 7 have a fixed
    // point to land on.
    std::mt19937_64 rng(2718281828ULL);
    std::uniform_real_distribution<double> uz(-4.0, 4.0);
    std::uniform_real_distribution<double> uk(0.1, 20.0);

    int nChecked = 0;
    for (int i = 0; i < 500; ++i) {
        const double zeta = uz(rng);
        if (zeta == 0.0) continue;
        const double K2 = uk(rng);
        const double s = uk(rng) * ((zeta > 0.0) ? 1.0 : -1.0);
        const double K0 = 0.4 * zeta * s;  // keeps temp = 0.6*zeta*s > 0

        const double temp1 = zeta * s - K0;
        if (temp1 <= 0.0) continue;
        const double w = ((zeta >= 0) ? 1.0 : -1.0) * std::sqrt(2.0 * temp1);
        const double v = zeta * std::sqrt(K2);
        // Skip the band in which bnTail deliberately does NOT form r*.
        if (std::fabs(w) <= spa::kWSingularity || v == 0.0 || (v / w) <= 0.0)
            continue;

        for (bool lower : {false, true}) {
            const double want = math::pnorm(w + (1.0 / w) * std::log(v / w), 0.0,
                                            1.0, lower);
            spa::Status st = spa::Status::FallbackMaxIter;
            const double got = spa::bnTail(zeta, s, K0, K2, lower, st);
            CHECK(st == spa::Status::SpaOk);
            // Not bit-identical, and deliberately so: the production sites
            // write (1/w)*log(v/w), which rounds twice, while bnTail writes
            // log(v/w)/w, which rounds once and is the more accurate of the
            // two.  The residual difference is a few ULP in r*, amplified into
            // the p-value by |d log p / d r*| = |r*|.  Measured worst case
            // across this sweep is 9.4e-15 relative.
            CHECK_REL(got, want, 1e-13);
        }
        ++nChecked;
    }
    CHECK_MSG(nChecked > 400, "only " + std::to_string(nChecked) + " checked");
}

// ══════════════════════════════════════════════════════════════════════
// § 7  Two-sided assembly
// ══════════════════════════════════════════════════════════════════════

TEST(combine_tails_sums_and_reports_negative_log10) {
    const spa::TwoSided r = spa::combineTails(1e-8, 3e-9, spa::Status::SpaOk,
                                              spa::Status::SpaOk);
    CHECK(r.status == spa::Status::SpaOk);
    CHECK_REL(r.p, 1.3e-8, 1e-15);
    CHECK_REL(r.negLog10p, -std::log10(1.3e-8), 0.0);
    CHECK(r.negLog10p > 0.0);

    // p == 1 gives +0, not -0.
    const spa::TwoSided one = spa::combineTails(0.5, 0.5, spa::Status::SpaOk,
                                                spa::Status::SpaOk);
    CHECK(one.p == 1.0);
    CHECK(one.negLog10p == 0.0);
    CHECK(!std::signbit(one.negLog10p));
}

TEST(combine_tails_never_reports_half_a_p_value) {
    // spamixlocalp.cpp:1044-1054 added each tail only if that tail's root
    // converged and produced NaN only when both failed, so a marker with one
    // good tail was reported at approximately half its correct p-value, with
    // no diagnostic.
    // GuardW is deliberately absent: it is a degraded success that always
    // returns a probability, so the (GuardW, NaN) pair below is not a state
    // bnTail can produce.  Its own behaviour is asserted after the loop.
    const spa::Status failures[] = {
        spa::Status::FallbackMaxIter, spa::Status::FallbackGuardTemp, spa::Status::FallbackGuardCurv,
        spa::Status::FallbackNonFinite,
    };
    for (spa::Status f : failures) {
        const spa::TwoSided a =
            spa::combineTails(1e-8, kNaN, spa::Status::SpaOk, f);
        CHECK(std::isnan(a.p));
        CHECK(std::isnan(a.negLog10p));
        CHECK(a.status == f);

        const spa::TwoSided b =
            spa::combineTails(kNaN, 1e-8, f, spa::Status::SpaOk);
        CHECK(std::isnan(b.p));
        CHECK(b.status == f);

        const spa::TwoSided c =
            spa::combineTailsLog(-18.4, kNaN, spa::Status::SpaOk, f);
        CHECK(std::isnan(c.p));
        CHECK(c.status == f);
    }

    // One degraded tail and one ordinary tail: a real p, and a status that
    // says the modified root was dropped on one side.
    const spa::TwoSided g = spa::combineTails(
        1e-8, 0.5, spa::Status::SpaOk, spa::Status::SpaWSingular);
    CHECK(g.status == spa::Status::SpaWSingular);
    CHECK(spa::statusIsUsable(g.status));
    CHECK_REL(g.p, 0.50000001, 1e-15);

    // ...but a genuine failure on the other side still wins.
    const spa::TwoSided h = spa::combineTails(
        kNaN, 0.5, spa::Status::FallbackMaxIter, spa::Status::SpaWSingular);
    CHECK(std::isnan(h.p));
    CHECK(h.status == spa::Status::FallbackMaxIter);
}

TEST(combine_tails_never_turns_nan_into_one) {
    // wtcoxg.cpp:618 wrote std::min(1.0, pval1 + pval2); std::min(a, NaN)
    // returns a, so every SPA failure surfaced as a perfectly null marker and
    // was then inverted through qchisq into an inflated var_S -- a different
    // statistic, not a conservative one.
    CHECK(std::min(1.0, kNaN) == 1.0);  // the behavior being fixed

    const spa::TwoSided r =
        spa::combineTails(kNaN, 0.25, spa::Status::SpaOk, spa::Status::SpaOk);
    CHECK(std::isnan(r.p));
    CHECK(r.status == spa::Status::FallbackNonFinite);

    const spa::TwoSided r2 =
        spa::combineTails(0.25, kInf, spa::Status::SpaOk, spa::Status::SpaOk);
    CHECK(std::isnan(r2.p));
    CHECK(r2.status == spa::Status::FallbackNonFinite);
}

TEST(combine_tails_clamps_into_the_unit_interval_without_flooring_at_dbl_min) {
    // Over-unity sums are clamped down...
    const spa::TwoSided hi = spa::combineTails(0.9, 0.9, spa::Status::SpaOk,
                                               spa::Status::SpaOk);
    CHECK(hi.p == 1.0);
    CHECK(hi.status == spa::Status::SpaOk);
    CHECK(hi.negLog10p == 0.0);

    // ...but a p that underflowed to zero is reported as zero.  Substituting
    // DBL_MIN, as spamixlocalp.cpp:987 does, manufactures a
    // genome-wide-significant hit out of a numerical failure.
    const spa::TwoSided lo = spa::combineTails(0.0, 0.0, spa::Status::SpaOk,
                                               spa::Status::SpaOk);
    CHECK(lo.p == 0.0);
    CHECK(lo.p != std::numeric_limits<double>::min());
    CHECK(lo.negLog10p == kInf);
    CHECK(lo.status == spa::Status::SpaOk);
}

TEST(combine_tails_log_is_log_sum_exp) {
    // Agreement with the closed form wherever both are representable.
    for (double e1 = -1.0; e1 >= -900.0; e1 -= 37.0)
        for (double e2 = -1.0; e2 >= -900.0; e2 -= 53.0) {
            const spa::TwoSided lg = spa::combineTailsLog(
                e1, e2, spa::Status::SpaOk, spa::Status::SpaOk);
            const double m = std::fmax(e1, e2);
            const double want =
                -(m + std::log1p(std::exp(std::fmin(e1, e2) - m))) / kLn10;
            CHECK_REL(lg.negLog10p, want, 1e-14);
            CHECK(lg.status == spa::Status::SpaOk);
        }

    // Two tails that both underflow on the linear scale still combine.
    const spa::TwoSided deep = spa::combineTailsLog(
        -920.0, -921.0, spa::Status::SpaOk, spa::Status::SpaOk);
    CHECK(deep.p == 0.0);
    CHECK(std::isfinite(deep.negLog10p));
    CHECK_REL(deep.negLog10p,
              -(-920.0 + std::log1p(std::exp(-1.0))) / kLn10, 1e-14);

    // Equal tails: log(2*e^x) = x + log 2.
    const spa::TwoSided eq = spa::combineTailsLog(
        -1000.0, -1000.0, spa::Status::SpaOk, spa::Status::SpaOk);
    CHECK_REL(eq.negLog10p, (1000.0 - 0.69314718055994530942) / kLn10, 1e-14);

    // Both tails exactly zero probability.
    const spa::TwoSided zero = spa::combineTailsLog(
        -kInf, -kInf, spa::Status::SpaOk, spa::Status::SpaOk);
    CHECK(zero.p == 0.0);
    CHECK(zero.negLog10p == kInf);

    // Over-unity is clamped the same way the linear form clamps it.
    const spa::TwoSided over = spa::combineTailsLog(
        -0.1, -0.1, spa::Status::SpaOk, spa::Status::SpaOk);
    CHECK(over.p == 1.0);
    CHECK(over.negLog10p == 0.0);
    CHECK(!std::signbit(over.negLog10p));
}

// ══════════════════════════════════════════════════════════════════════
// § 7b  math::pnormLog and the log domain end to end  (spa_unify N2)
// ══════════════════════════════════════════════════════════════════════

// Property 1 of the Stage 8 gate: pnormLog matches log(pnorm) to 1e-13
// wherever the latter is representable, in both tails and with a non-unit
// scale.  Above the -37 asymptote the two are the same computation, so what
// this really pins is that the standardization, the tail selection and the
// branch point are all right; below it, the next test takes over.
TEST(pnormlog_matches_log_pnorm_where_pnorm_is_representable) {
    int nChecked = 0;
    for (double x = -37.0; x <= 8.0; x += 0.03125) {
        for (bool lower : {true, false}) {
            const double p = math::pnorm(x, 0.0, 1.0, lower);
            if (!(p > 0.0) || !std::isfinite(std::log(p))) continue;
            const double lp = math::pnormLog(x, 0.0, 1.0, lower);
            CHECK_MSG(std::fabs(lp - std::log(p)) <=
                          1e-13 * std::fmax(1.0, std::fabs(std::log(p))),
                      "x=" + std::to_string(x) + " lower=" +
                          std::to_string(lower) + " lp=" + std::to_string(lp));
            ++nChecked;
        }
    }
    CHECK(nChecked > 2000);

    // Non-unit mean and scale, which is how WtCoxG calls it.
    for (double sd : {0.03, 1.0, 17.0}) {
        for (double mu : {-2.0, 0.0, 5.0}) {
            for (double x : {-4.0, -0.5, 0.0, 1.25, 6.0}) {
                const double xs = mu + sd * x;
                for (bool lower : {true, false}) {
                    const double p = math::pnorm(xs, mu, sd, lower);
                    const double lp = math::pnormLog(xs, mu, sd, lower);
                    CHECK_NEAR(lp, std::log(p),
                               1e-13 * std::fmax(1.0, std::fabs(std::log(p))));
                }
            }
        }
    }

    // Edges.  An infinite abscissa in the vanishing tail is -inf, not the
    // log(0 + 1e-300) = -690.8 that the removed `log_p` flag substituted:
    // a manufactured floor is exactly what decision L3 forbids.
    const double inf = std::numeric_limits<double>::infinity();
    CHECK(math::pnormLog(-inf, 0.0, 1.0, true) == -inf);
    CHECK(math::pnormLog(inf, 0.0, 1.0, false) == -inf);
    CHECK(math::pnormLog(inf, 0.0, 1.0, true) == 0.0);
    CHECK(math::pnormLog(-inf, 0.0, 1.0, false) == 0.0);
    CHECK(std::isnan(math::pnormLog(kNaN)));
}

// The half of pnormLog that has no linear-scale counterpart, checked against
// the expansion evaluated in long double.  The sweep straddles -37 so the two
// branches are pinned against one reference and the join is therefore pinned
// too, and it continues to -700 where Phi is 1e-106000 and no double could
// hold it.
TEST(pnormlog_is_accurate_far_past_the_linear_underflow) {
    for (double x = -30.0; x >= -700.0; x -= 0.5) {
        const long double t = static_cast<long double>(x);
        const long double y = 1.0L / (t * t);
        const long double series =
            y * (-1.0L + y * (3.0L + y * (-15.0L + y * (105.0L +
                 y * (-945.0L + y * 10395.0L)))));
        const long double ref = -0.5L * t * t - std::log(-t) -
                                0.9189385332046727417803297364056L +
                                std::log1p(series);
        const double lp = math::pnormLog(x, 0.0, 1.0, true);
        CHECK_MSG(std::isfinite(lp), "non-finite at x=" + std::to_string(x));
        CHECK_MSG(std::fabs(lp - static_cast<double>(ref)) <=
                      1e-13 * std::fabs(static_cast<double>(ref)),
                  "x=" + std::to_string(x) + " lp=" + std::to_string(lp));
        // The upper tail is the mirror image, exactly.
        CHECK(math::pnormLog(-x, 0.0, 1.0, false) == lp);
    }
    // Where the linear scale has nothing left at all.
    CHECK(math::pnorm(-40.0, 0.0, 1.0, true) == 0.0);
    CHECK(std::isfinite(math::pnormLog(-40.0, 0.0, 1.0, true)));
}

// Property 2 of the gate, stated on the output surface rather than on bnTail:
// a two-sided p of 1e-400 survives into the LOG10P column, and P reports the
// honest 0 rather than a fabricated floor.
TEST(a_p_of_1e_minus_400_is_recovered_in_log10p) {
    // Two-sided, so each tail carries p/2.
    const double logHalfTarget = -400.0 * kLn10 - kLn2;
    double blo = 1.0, bhi = 200.0;
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (blo + bhi);
        if (math::pnormLog(-mid, 0.0, 1.0, true) > logHalfTarget) blo = mid;
        else bhi = mid;
    }
    const double w = 0.5 * (blo + bhi);
    const double lTail = math::pnormLog(-w, 0.0, 1.0, true);

    const spa::TwoSided lg = spa::combineTailsLog(
        lTail, lTail, spa::Status::SpaOk, spa::Status::SpaOk);
    CHECK(lg.status == spa::Status::SpaOk);
    CHECK_CLOSE(lg.negLog10p, 400.0, 1e-12, 0.0);
    CHECK(lg.p == 0.0);           // honest zero, not DBL_MIN

    // The linear assembly of the same two tails has nothing left: both tails
    // underflowed, so it reports p = 0 with negLog10p = +inf.  That is the
    // defect L3 exists to route around, and the reason LOG10P is a column.
    const spa::TwoSided lin = spa::combineTails(
        std::exp(lTail), std::exp(lTail), spa::Status::SpaOk,
        spa::Status::SpaOk);
    CHECK(lin.p == 0.0);
    CHECK(std::isinf(lin.negLog10p));

    // log-sum-exp, not a naive sum of exponentials: the two tails are equal,
    // so their combination must be exactly one ln 2 above either of them.
    CHECK_NEAR(lg.negLog10p * -kLn10, lTail + kLn2, 1e-11);
}

// Property 3 of the gate: LOG10P and -log10(P) agree to 1e-9 wherever
// P > 1e-300.  This is checked at FULL double precision here because it cannot
// be checked end to end on the fixture -- LOG10P is printed to six significant
// figures, so a file-level comparison resolves no better than ~1e-6 relative.
TEST(log10p_agrees_with_minus_log10_p_at_full_precision) {
    double worst = 0.0, worstAt = 0.0;
    int nChecked = 0;
    for (double w = 0.05; w <= 26.0; w += 0.005) {
        // Two independent tails, deliberately unequal so the assembly is a
        // genuine sum rather than a doubling.
        const double lu = math::pnormLog(-w, 0.0, 1.0, true);
        const double ll = math::pnormLog(-1.37 * w - 0.2, 0.0, 1.0, true);
        const double pu = math::pnorm(-w, 0.0, 1.0, true);
        const double pl = math::pnorm(-1.37 * w - 0.2, 0.0, 1.0, true);

        const spa::TwoSided lin = spa::combineTails(
            pu, pl, spa::Status::SpaOk, spa::Status::SpaOk);
        const spa::TwoSided lg = spa::combineTailsLog(
            lu, ll, spa::Status::SpaOk, spa::Status::SpaOk);
        if (!(lin.p > 1e-300)) continue;

        const double d = std::fabs(lg.negLog10p - lin.negLog10p) /
                         std::fmax(1.0, lin.negLog10p);
        if (d > worst) { worst = d; worstAt = w; }
        CHECK_MSG(d <= 1e-9, "w=" + std::to_string(w) + " rel diff " +
                                 std::to_string(d));
        // The two assemblies must also agree on p itself while it is
        // representable.
        CHECK_REL(lg.p, lin.p, 1e-12);
        ++nChecked;
    }
    CHECK(nChecked > 5000);
    CHECK_MSG(worst < 1e-9, "worst " + std::to_string(worst) + " at w=" +
                                std::to_string(worstAt));
}

TEST(status_severity_ordering_is_total_and_deterministic) {
    const spa::Status all[] = {
        spa::Status::SpaOk, spa::Status::Normal, spa::Status::FallbackMaxIter,
        spa::Status::SpaWSingular, spa::Status::FallbackGuardTemp,
        spa::Status::FallbackGuardCurv, spa::Status::FallbackNonFinite,
        spa::Status::NaPostFail, spa::Status::NaNoTest,
    };
    for (spa::Status a : all)
        for (spa::Status b : all) {
            // Symmetric, and never invents a status neither side reported.
            CHECK(spa::worseStatus(a, b) == spa::worseStatus(b, a));
            const spa::Status w = spa::worseStatus(a, b);
            CHECK(w == a || w == b);
            CHECK(spa::statusSeverity(w) ==
                  std::max(spa::statusSeverity(a), spa::statusSeverity(b)));
            // A failure on either side always wins over a success.
            if (!spa::statusIsUsable(a) || !spa::statusIsUsable(b))
                CHECK(!spa::statusIsUsable(w));
        }
}

// ══════════════════════════════════════════════════════════════════════
// § 8  The normal short-circuit
// ══════════════════════════════════════════════════════════════════════

TEST(normal_two_sided_and_its_log) {
    CHECK_REL(spa::normalTwoSided(0.0), 1.0, 0.0);
    CHECK_REL(spa::normalTwoSided(1.959963984540054), 0.05, 1e-12);
    CHECK_REL(spa::normalTwoSided(-1.959963984540054), 0.05, 1e-12);
    // 2*Phi(-5) = 2 * 2.866515718791939e-7, an independently tabulated value.
    CHECK_REL(spa::normalTwoSided(5.0), 5.733031437583878e-7, 1e-12);

    for (double z = 0.5; z <= 37.0; z += 0.5) {
        const double p = spa::normalTwoSided(z);
        const double lp = spa::normalTwoSidedLog(z);
        CHECK_MSG(std::fabs(lp - std::log(p)) <= 1e-13 * std::fmax(1.0, std::fabs(lp)),
                  "z=" + std::to_string(z));
    }

    // Past the linear floor the log form keeps working.
    for (double z = 40.0; z <= 300.0; z += 20.0) {
        const double lp = spa::normalTwoSidedLog(z);
        CHECK(std::isfinite(lp));
        CHECK(lp < -0.49 * z * z);
        CHECK(spa::normalTwoSided(z) == 0.0);
    }

    CHECK(std::isnan(spa::normalTwoSidedLog(kNaN)));
    CHECK(std::isnan(spa::normalBranch(kNaN).p));
    CHECK(spa::normalBranch(kNaN).status == spa::Status::NaNoTest);
}

TEST(normal_branch_bundle_is_self_consistent) {
    for (double z : {0.0, 0.5, 2.0, 6.0, 25.0}) {
        const spa::TwoSided b = spa::normalBranch(z);
        CHECK(b.status == spa::Status::Normal);
        CHECK_REL(b.p, spa::normalTwoSided(z), 0.0);
        if (b.p > 0.0)
            CHECK_MSG(std::fabs(b.negLog10p + std::log10(b.p)) <=
                          1e-13 * std::fmax(1.0, b.negLog10p),
                      "z=" + std::to_string(z));
    }
    CHECK(spa::normalBranch(0.0).p == 1.0);
    CHECK(spa::normalBranch(0.0).negLog10p == 0.0);
    CHECK(!std::signbit(spa::normalBranch(0.0).negLog10p));
}

TEST(sign_of_is_zero_at_zero) {
    // The distinguishing property against std::copysign and the
    // `zeta >= 0 ? 1 : -1` spelling, both of which return +/-1 at zero and so
    // produce a signed sqrt(0) that is then divided by.
    CHECK(spa::signOf(0.0) == 0.0);
    CHECK(spa::signOf(-0.0) == 0.0);
    CHECK(spa::signOf(1e-300) == 1.0);
    CHECK(spa::signOf(-1e-300) == -1.0);
    CHECK(spa::signOf(kInf) == 1.0);
    CHECK(spa::signOf(-kInf) == -1.0);
    CHECK(spa::signOf(kNaN) == 0.0);
}

// ══════════════════════════════════════════════════════════════════════
// End-to-end: solver + tail + assembly on a realistic problem
// ══════════════════════════════════════════════════════════════════════

TEST(end_to_end_two_sided_spa_p_value) {
    // The shape every production site has: an outlier set plus a Gaussian
    // remainder, upper and lower tails solved from mirrored scores, then
    // combined.  Checked properties: both roots take the sign of their score,
    // the p-value is in (0, 1] and strictly decreasing in |z|, and the
    // log-domain assembly agrees with the linear one wherever the latter is
    // representable.
    std::mt19937_64 rng(161803398ULL);
    std::uniform_real_distribution<double> umaf(0.005, 0.05);
    std::normal_distribution<double> nresid(0.0, 1.0);

    MixedCgf c;
    c.resid.resize(30);
    c.maf.resize(30);
    for (int i = 0; i < 30; ++i) {
        c.resid[i] = nresid(rng) * 4.0;
        c.maf[i] = umaf(rng);
    }
    c.var = 5.0;

    const spa::K012 at0 = c.kFull(0.0);
    const double mean0 = at0.k1;   // K'(0), since c.s is still zero
    const double sd = std::sqrt(at0.k2);

    double prevP = kInf;
    int nChecked = 0;
    for (double z = 2.5; z <= 8.0; z += 0.5) {
        const double sUpper = mean0 + z * sd;
        const double sLower = mean0 - z * sd;

        MixedCgf cu = c; cu.s = sUpper;
        MixedCgf cl = c; cl.s = sLower;

        spa::SolveOpts ou; ou.init = 0.01;  ou.scoreSign = 1.0;
        spa::SolveOpts ol; ol.init = -0.01; ol.scoreSign = -1.0;

        const spa::Saddle ru = spa::solveSaddlepoint(
            sUpper, [&](double t) { return cu.k12(t); },
            [&](double t) { return cu.kFull(t); }, ou);
        const spa::Saddle rl = spa::solveSaddlepoint(
            sLower, [&](double t) { return cl.k12(t); },
            [&](double t) { return cl.kFull(t); }, ol);

        CHECK_MSG(ru.status == spa::Status::SpaOk, "upper " + statusStr(ru.status));
        CHECK_MSG(rl.status == spa::Status::SpaOk, "lower " + statusStr(rl.status));
        if (ru.status != spa::Status::SpaOk || rl.status != spa::Status::SpaOk)
            continue;
        CHECK(ru.zeta > 0.0);
        CHECK(rl.zeta < 0.0);

        spa::Status su, sl;
        const double pu = spa::bnTail(ru.zeta, sUpper, ru.K0, ru.K2, false, su);
        const double pl = spa::bnTail(rl.zeta, sLower, rl.K0, rl.K2, true, sl);
        const double lu = spa::bnTailLog(ru.zeta, sUpper, ru.K0, ru.K2, false, su);
        const double ll = spa::bnTailLog(rl.zeta, sLower, rl.K0, rl.K2, true, sl);

        CHECK_MSG(su == spa::Status::SpaOk, "tail upper " + statusStr(su));
        CHECK_MSG(sl == spa::Status::SpaOk, "tail lower " + statusStr(sl));
        if (su != spa::Status::SpaOk || sl != spa::Status::SpaOk) continue;

        const spa::TwoSided lin = spa::combineTails(pu, pl, su, sl);
        const spa::TwoSided lg = spa::combineTailsLog(lu, ll, su, sl);

        CHECK(lin.status == spa::Status::SpaOk);
        CHECK(lin.p > 0.0 && lin.p <= 1.0);
        CHECK_REL(lg.p, lin.p, 1e-13);
        CHECK_NEAR(lg.negLog10p, lin.negLog10p, 1e-11);

        // Strictly decreasing in |z|, since both tails are.
        CHECK_MSG(lin.p < prevP, "z=" + std::to_string(z) + " p not decreasing");
        prevP = lin.p;
        ++nChecked;
    }
    CHECK(nChecked >= 11);
}

TINYTEST_MAIN
