// spa.hpp — Unified saddlepoint solver, tail kernel and two-sided assembly.
//
// GRAB carried six independent saddlepoint implementations serving nine
// methods.  All six computed the same tail formula — the Barndorff-Nielsen
// modified signed root
//
//     zeta solves K'(zeta) = s
//     w  = sgn(zeta) * sqrt(2 * (zeta*s - K(zeta)))
//     v  = zeta * sqrt(K''(zeta))
//     r* = w + log(v/w) / w
//     p  = Phi(r*)              (lower tail)
//        = 1 - Phi(r*)          (upper tail)
//
// — but none computed it the same way.  What differed was the guard set (two
// sites had none), the root finder (three incompatible families), the
// tolerance conventions, and the failure semantics.  This header is tier 1 of
// the three-tier split described in dev-notes/methods/spa_unify/02_design.md:
// everything shared by all six sites, and nothing else.
//
//   Tier 1  solver, BN tail, guards, normal short-circuit, two-sided          <- here
//           assembly, status
//   Tier 2  binomial CGF kernels (uniform-MAF, individual-MAF, hapcount)      util/spa_cgf
//   Tier 3  SPAGRM family blocks, SPACox empirical table, score prep,
//           variance-ratio rescaling, initializers, per-method tolerances     per-method
//
// THE BOUNDARY RULE — shared code never touches a subject.
//
// Anything that iterates over subjects, or decides what a subject
// contributes, stays in the owning method and reaches the solver as an
// inlined template callable.  This is not stylistic.  The callback fires once
// per Newton iteration, not once per marker: worst case is 2 tails * ~100
// iterations per p-value, times markers * phenotypes * (for LEAF) clusters,
// which is 1e8 to 1e9 invocations on a realistic run.  A std::function there
// is an indirect branch that also blocks inlining of the CGF body, defeating
// the SIMD kernels entirely.  Hence every entry point below is either a
// function template on the callable type, or a scalar-argument leaf function.
//
// COST STRUCTURE.  The Newton loop needs only K' and K''; K itself is needed
// once per tail, at the converged root.  Since K' and K'' can be evaluated
// without a logarithm while K cannot (see spa_unify P2 and N1), the solver
// takes two callables:
//
//     eval(t)     -> K12 {k1, k2}       cheap; the Newton-loop payload
//     evalFull(t) -> K012 {k0, k1, k2}  terminal only; adds the logarithm
//
// `evalFull` is invoked exactly once, at the returned root.  A single-callable
// overload is provided for methods whose CGF produces K0 for free; there the
// terminal evaluation is skipped outright when the loop's last abscissa
// already equals the returned root.
//
// FAILURE SEMANTICS.  Status is part of every return value and is never an
// out-parameter a caller can forget to read.  Three of the six original sites
// computed a convergence flag and dropped it on the floor; two never computed
// one at all.  A saddlepoint that never converged then produced an
// ordinary-looking finite p-value with no NaN, no flag and no warning.  Here a
// degenerate input yields NaN in the probability *and* a Status naming the
// reason, which the calling method reports in its SPA_STATUS column.

#pragma once

#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>

#include "util/math_helper.hpp"

namespace spa {

// ──────────────────────────────────────────────────────────────────────
// § 1  Status
// ──────────────────────────────────────────────────────────────────────

enum class Status : uint8_t {
    Converged = 0,   // step and residual both within tolerance; no guard fired
    MaxIter,         // iteration or bracket-expansion budget exhausted
    GuardTemp,       // zeta*s - K(zeta) < 0, so w is not real
    GuardCurv,       // K''(zeta) <= 0, so v is not real
    GuardW,          // |w| at or below the removable-singularity threshold
    NonFinite,       // zeta, a cumulant, or r* left the reals
    NormalBranch,    // |z| <= spaCutoff; the saddlepoint was never attempted
};

// Token emitted in the SPA_STATUS output column.
inline const char *statusName(Status s) noexcept {
    switch (s) {
        case Status::Converged:    return "OK";
        case Status::MaxIter:      return "MAXITER";
        case Status::GuardTemp:    return "GUARD_TEMP";
        case Status::GuardCurv:    return "GUARD_CURV";
        case Status::GuardW:       return "GUARD_W";
        case Status::NonFinite:    return "NONFINITE";
        case Status::NormalBranch: return "NORMAL";
    }
    return "NONFINITE";
}

// True when the status means "no usable probability was produced", i.e. the
// method must report NA.  Converged and NormalBranch are the two successes.
inline bool statusIsFailure(Status s) noexcept {
    return s != Status::Converged && s != Status::NormalBranch;
}

// Ranking used when two tails disagree.  A specific guard outranks a bare
// non-convergence because it names the quantity that went wrong; NonFinite
// outranks everything because it is the least diagnosable.
inline int statusSeverity(Status s) noexcept {
    switch (s) {
        case Status::Converged:    return 0;
        case Status::NormalBranch: return 0;
        case Status::MaxIter:      return 1;
        case Status::GuardW:       return 2;
        case Status::GuardTemp:    return 3;
        case Status::GuardCurv:    return 4;
        case Status::NonFinite:    return 5;
    }
    return 5;
}

// The status a two-sided result inherits from its two tails.
inline Status worseStatus(Status a, Status b) noexcept {
    const int sa = statusSeverity(a), sb = statusSeverity(b);
    if (sa != sb) return (sa > sb) ? a : b;
    if (a == b) return a;
    // Equal severity but different value: the only such pair is
    // Converged / NormalBranch.  A genuine saddlepoint tail dominates a
    // normal-branch label, since the result did go through the SPA.
    return Status::Converged;
}

// ──────────────────────────────────────────────────────────────────────
// § 2  Small shared scalars
// ──────────────────────────────────────────────────────────────────────

// The three spellings of sgn(zeta) found across the original six sites
// (signOf, std::copysign, `zeta >= 0 ? 1 : -1`) are not equivalent at zero:
// copysign and the ternary both return +/-1 for +/-0.0, while this form
// returns 0.  Zero is the correct answer — it makes w exactly zero and routes
// the degenerate case to GuardW rather than emitting a signed square root of
// zero and dividing by it.
inline double signOf(double x) noexcept {
    return (x > 0.0) ? 1.0 : ((x < 0.0) ? -1.0 : 0.0);
}

namespace detail {

constexpr double kLn10      = 2.30258509299404568402;
constexpr double kLn2       = 0.69314718055994530942;
constexpr double kLogSqrt2Pi = 0.91893853320467274178;  // log(sqrt(2*pi))

// Below this abscissa the linear-scale normal CDF has lost all relative
// accuracy (Phi(-38) is denormal, Phi(-38.5) flushes to zero), so log Phi is
// taken from the asymptotic expansion instead.  -37 leaves Phi(-37) ~ 5.7e-300
// comfortably normal, so the two branches overlap in a regime where both are
// accurate.
constexpr double kLogPhiAsymptote = -37.0;

inline double quietNaN() noexcept {
    return std::numeric_limits<double>::quiet_NaN();
}

// log Phi(x).
//
// For x > kLogPhiAsymptote this is log(Phi(x)) evaluated through math::pnorm,
// which is accurate because Phi(x) is a normal double there and its relative
// error of a few ULP becomes an absolute error of a few ULP in the logarithm.
//
// For x <= kLogPhiAsymptote the linear scale has underflowed, and the Mills
// ratio asymptotic expansion is used instead:
//
//     Phi(x) = phi(x)/(-x) * (1 - 1/x^2 + 3/x^4 - 15/x^6 + 105/x^8 - ...)
//     log Phi(x) = -x^2/2 - log(-x) - log(sqrt(2*pi))
//                  + log1p(-1/x^2 + 3/x^4 - 15/x^6 + ...)
//
// At |x| = 37 the six retained terms leave a truncation error of order 1e-15
// in the log1p argument, hence 1e-15 absolute in a logarithm whose magnitude
// is ~690 — a relative accuracy of 1e-18.  The series is asymptotic, not
// convergent, which is exactly why it is used only far out in the tail.
//
// STAGE 8 (spa_unify N2).  This helper is the local stand-in for
// math::pnormLog.  When that lands in util/math_helper.hpp — where WtCoxG's
// two existing log_p=true consumers also need it — this body becomes a
// forward to it and the expansion is deleted from here.
inline double logPhi(double x) noexcept {
    if (std::isnan(x)) return quietNaN();
    if (x > kLogPhiAsymptote)
        return math::pnorm(x, 0.0, 1.0, /*lower_tail=*/true, /*log_p=*/true);
    if (x == -std::numeric_limits<double>::infinity())
        return -std::numeric_limits<double>::infinity();
    const double y = 1.0 / (x * x);
    const double series =
        y * (-1.0 + y * (3.0 + y * (-15.0 + y * (105.0 + y * (-945.0 + y * 10395.0)))));
    return -0.5 * x * x - std::log(-x) - kLogSqrt2Pi + std::log1p(series);
}

// The single site at which the tail probability meets the normal CDF.  Both
// bnTail and bnTailLog route through here so that the linear and log domains
// cannot drift apart — the class of defect this whole header exists to remove.
inline double phiTail(double x, bool lowerTail) noexcept {
    return math::pnorm(x, 0.0, 1.0, lowerTail);
}

inline double phiTailLog(double x, bool lowerTail) noexcept {
    return logPhi(lowerTail ? x : -x);
}

} // namespace detail

// ──────────────────────────────────────────────────────────────────────
// § 3  Cumulant bundles exchanged with the caller's CGF
// ──────────────────────────────────────────────────────────────────────

// Newton-loop payload.  k1 is the *residual* K'(t) - s, not K'(t): the caller
// subtracts s inside its own accumulation so that the association order stays
// under the owning method's control.  k2 is K''(t).
struct K12 {
    double k1;
    double k2;
};

// Terminal payload.  k0 is K(t); k1 and k2 are as in K12.
struct K012 {
    double k0;
    double k1;
    double k2;
};

// ──────────────────────────────────────────────────────────────────────
// § 4  Solver
// ──────────────────────────────────────────────────────────────────────

struct Saddle {
    double zeta;      // the returned root
    double K0;        // K(zeta)   — always evaluated AT zeta, never at a
    double K2;        // K''(zeta)   stale abscissa
    double residual;  // K'(zeta) - s, as achieved
    double lo;        // final bracket, containing the true root whenever
    double hi;        //   `bracketed` is true
    int    iters;     // safeguarded-Newton iterations executed
    int    nEval;     // cheap-callable evaluations (bracketing + loop)
    int    nEvalBracket;   // of those, the ones spent locating the bracket,
                           //   counting the evaluation at `init`.  Reported
                           //   separately because the two phases answer
                           //   different questions: nEvalBracket measures the
                           //   quality of the caller's initial guess, and
                           //   nEval - nEvalBracket the quality of the Newton
                           //   schedule.  The terminal evalFull is in neither.
    bool   bracketed; // a sign change in the residual was located
    bool   reusedTerminal;  // terminal abscissa equalled the loop's last one
    Status status;
};

struct SolveOpts {
    // Initial abscissa.  Methods that carry a per-marker guess (SPAGRM's
    // |Score|/Var, capped at 1.2) pass it here; 0 is always admissible
    // because K'(0) is the CGF mean and is finite for every CGF in GRAB.
    double init = 0.0;

    // Relative residual tolerance.  Convergence is
    //     |K'(t) - s| <= max(rtol * sqrt(K''(t)), 4 * eps * |s|).
    // Scaling by sqrt(K'') makes the criterion dimensionless: K' carries the
    // units of s and K'' the units of s^2.  Every tolerance in the original
    // six sites was absolute and unscaled, so the same numeric value meant
    // different things at different sample sizes (spa_unify D7).  The second
    // term is the representation-error floor of the difference K'(t) - s; it
    // is inactive by roughly nine orders of magnitude in the regime the SPA
    // branch is entered in, and exists only so that a pathological problem
    // stops instead of burning the whole iteration budget below the noise.
    double rtol = 1e-6;

    // Bracket-collapse tolerance.  When hi - lo falls to
    // stepTol * max(1, |t|) the root is pinned to that width and no further
    // refinement is possible, so the result is accepted.
    double stepTol = 1e-8;

    int maxIter = 100;

    // When non-zero, restrict the search to sgn(zeta) == sgn(scoreSign).
    // This is what makes sgn(zeta) a valid proxy for sgn(w), and is the one
    // genuinely useful safeguard the original Family B finders had.
    double scoreSign = 0.0;

    // Bracket expansion.  The first outward probe is the Newton step at
    // `init`, |K'(init) - s| / K''(init), which for a residual that is convex
    // on the side being searched lands beyond the root and brackets it in one
    // evaluation, with a bracket only as wide as the first-order estimate.
    //
    // `bracketStep` is the fallback distance, used when K''(init) is not
    // usable, and the distance the schedule jumps to on the first probe that
    // lands on the same side of the root: the probe distance becomes
    //     max(bracketGrow * previous, bracketStep * max(1, |init|))
    // once, and grows by bracketGrow thereafter.  Falling back to a coarse
    // constant rather than doubling a possibly-tiny Newton step is what bounds
    // the worst-case expansion cost at one evaluation more than a schedule that
    // started coarse in the first place.
    double bracketStep = 1.0;
    double bracketGrow = 2.0;
    int    maxExpand   = 64;
};

// Instrumentation hook.  An empty struct with an inline no-op call operator
// costs nothing: the compiler eliminates it entirely.  tests/spa_solver_test
// substitutes a recorder to verify the bracket invariant at every iteration,
// which cannot be checked from the return value alone.
struct NullTrace {
    void operator()(int, double, double, double, double, double) const noexcept {}
};

namespace detail {

// Wraps a single K012-producing callable so it can serve both the Newton loop
// and the terminal evaluation, caching the most recent result so that a
// terminal abscissa coinciding with the loop's last one costs nothing.  This
// is the reuse spa_unify P3 asks for, placed where it belongs: only the
// solver knows whether the last evaluation matches the returned root.
template <class EvalFull>
struct CachedFull {
    EvalFull &fn;
    double last = quietNaN();
    K012  value{quietNaN(), quietNaN(), quietNaN()};
    bool  valid = false;

    explicit CachedFull(EvalFull &f) : fn(f) {}

    K12 cheap(double t) {
        value = fn(t);
        last = t;
        valid = true;
        return K12{value.k1, value.k2};
    }

    K012 full(double t) {
        if (valid && t == last) return value;
        value = fn(t);
        last = t;
        valid = true;
        return value;
    }
};

// Two distinct callables: the cheap one never produces K0, so the terminal
// evaluation is unconditional.
template <class Eval, class EvalFull>
struct SplitFull {
    Eval &cheapFn;
    EvalFull &fullFn;

    SplitFull(Eval &c, EvalFull &f) : cheapFn(c), fullFn(f) {}

    K12  cheap(double t) { return cheapFn(t); }
    K012 full(double t) { return fullFn(t); }
};

// ──────────────────────────────────────────────────────────────────
// The solver proper.
//
// Dynamically bracketed safeguarded Newton.  Properties, and where each came
// from (spa_unify 02_design, "The root finder"):
//
//   * a bracket [lo, hi] expanded outward from `init` until the residual
//     changes sign.  None of the three original families maintained one;
//     Family A's "bisection" only damped the step and never established a
//     sign change, so it could not guarantee convergence.  The first probe is
//     placed at the Newton point rather than a constant distance away, so the
//     bracket is as tight as the caller's first-order guess allows.
//   * bisection whenever the Newton step leaves the bracket, produces a
//     non-finite residual, or fails to reduce |residual|; plus a forced
//     bisection whenever the Newton step itself has stopped contracting, which
//     the design table does not list but which is what actually bounds the
//     iteration count (see the comment at the safeguard).
//   * optional enforcement of sgn(zeta) == sgn(Score), implemented as a
//     restriction of the search interval to a closed half-line rather than
//     Family B's halve-until-the-sign-agrees loop.
//   * apply the step, then test.  Family A tested |step| < tol and broke
//     *before* applying, returning a root carrying an unapplied correction of
//     up to 1e-3.
//   * K'' > 0 checked before every division by it.  Five of the six original
//     finders divided by K'' unchecked.
//   * relative convergence on the residual, scaled by sqrt(K'').
//   * terminal cumulants evaluated at the returned root.
//
// The returned root is the last accepted iterate, which is always an endpoint
// of the final bracket and therefore always inside it.
// ──────────────────────────────────────────────────────────────────
template <class Evaluator, class Trace>
Saddle solveCore(double s, Evaluator &ev, const SolveOpts &opt, Trace &trace) {
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = quietNaN();

    Saddle out{};
    out.zeta = nan;
    out.K0 = nan;
    out.K2 = nan;
    out.residual = nan;
    out.lo = nan;
    out.hi = nan;
    out.iters = 0;
    out.nEval = 0;
    out.nEvalBracket = 0;
    out.bracketed = false;
    out.reusedTerminal = false;
    out.status = Status::MaxIter;

    // The residual noise floor: below 4 eps |s| the difference K'(t) - s is
    // dominated by the rounding error of its own subtraction.
    const double noiseFloor =
        4.0 * std::numeric_limits<double>::epsilon() * std::fabs(s);

    // Half-line implied by the sign constraint.
    double L = -inf, U = inf;
    if (opt.scoreSign > 0.0) L = 0.0;
    else if (opt.scoreSign < 0.0) U = 0.0;

    double t = opt.init;
    if (!std::isfinite(t)) t = 0.0;
    if (t < L) t = L;
    if (t > U) t = U;

    // ── Initial evaluation, retreating toward the origin if it fails ──
    // K'(0) is the CGF mean and is finite for every CGF in GRAB, so halving
    // toward zero always reaches a usable abscissa.
    K12 e = ev.cheap(t);
    ++out.nEval;
    for (int k = 0; k < opt.maxExpand && !std::isfinite(e.k1); ++k) {
        const double tn = 0.5 * t;
        if (tn == t) break;
        t = tn;
        e = ev.cheap(t);
        ++out.nEval;
    }

    double bestT = t, bestF = std::fabs(e.k1);
    if (!std::isfinite(bestF)) bestF = inf;

    double lo = nan, hi = nan;

    if (!std::isfinite(e.k1)) {
        // Nothing usable anywhere along the retreat path.
        out.status = Status::NonFinite;
    } else if (e.k1 == 0.0) {
        lo = hi = t;
        out.bracketed = true;
        out.status = Status::Converged;
    } else {
        // ── Bracket expansion ──
        // K' is non-decreasing (K'' >= 0 for any CGF), so the residual
        // K'(t) - s is non-decreasing: a negative residual puts the root to
        // the right, a positive one to the left.  This direction test uses
        // only the sign of the residual, so it survives a K'' that is noisy
        // or wrongly signed.
        const double dir = (e.k1 < 0.0) ? 1.0 : -1.0;
        // `a` is the near endpoint of the probe interval and `ea` its
        // cumulants, retained so that the bracketing exit can choose the
        // better of the two endpoints to start the Newton loop from.
        double a = t;
        K12 ea = e;

        // ── Sizing the first probe ──
        //
        // The first probe is placed at the Newton point, |K'(t) - s| / K''(t)
        // away from `init`.  For a residual that is convex on the side being
        // searched the tangent line lies below it, so the Newton point falls
        // beyond the root and ONE probe brackets it — with a bracket only as
        // wide as the first-order estimate, rather than as wide as an arbitrary
        // constant.
        //
        // This is where the excess evaluations were.  The previous schedule took
        // the LARGER of `bracketStep * max(1, |init|)` (1.0 by default) and the
        // Newton step, so for the roots the SPA branch actually visits — zeta of
        // order 1e-3 to 1, initialized from the caller's s/K''(0) estimate — the
        // first probe landed a full unit away, bracketed on the first try, and
        // handed the Newton loop a bracket three orders of magnitude wider than
        // the distance to the root.  The loop then spent its budget bisecting
        // that width.  Measured on the 50000 x 20000 synthetic cohort: SPAsqr
        // 14.71 cheap evaluations per tail, of which 12.71 were loop iterations
        // and 7.47 of those were bisections.
        //
        // `bracketStep` is retained as the fallback for an unusable K'' and as
        // the SECOND probe distance: if the Newton point does not bracket, the
        // schedule falls back to the coarse constant immediately, so the
        // worst-case expansion cost is one evaluation more than before rather
        // than a long geometric walk up from a tiny first step.
        double step = 0.0;
        if (e.k2 > 0.0 && std::isfinite(e.k2)) {
            const double newton = std::fabs(e.k1 / e.k2);
            if (std::isfinite(newton) && newton > 0.0) step = newton;
        }
        const double coarse =
            std::fmax(opt.bracketStep * std::fmax(1.0, std::fabs(t)), step);
        if (!(step > 0.0)) step = coarse;
        if (!(step > 0.0)) step = opt.bracketStep;
        // Consumed by the first probe that lands on the same side of the root:
        // see the `step` update at the bottom of the expansion loop.  It counts
        // usable probes only, not probes that fell in a non-finite region — the
        // step-halving retreat those trigger must not be mistaken for progress,
        // or the schedule creeps up to the edge of the non-finite region and
        // never steps over it.
        bool coarseUnused = true;

        // ── The first probe must land on a different double ──
        //
        // The bare Newton step, unlike the coarse constant it replaced, can be
        // arbitrarily small.  When it falls below half an ulp of `init` the
        // probe `init + dir*step` rounds straight back onto `init`, the loop
        // below breaks on its first pass with `bracketed` still false, and the
        // solve returns MaxIter after a single evaluation — while `zeta` holds
        // the correct root to the last bit.  Since statusIsFailure(MaxIter) is
        // true, the caller then reports NA for a marker whose saddlepoint was
        // solved exactly.  The trigger is |K'(init) - s| / K''(init) <
        // ulp(init)/2, i.e. "init already lies within half an ulp of the root",
        // so the failure probability would be monotone INCREASING in the quality
        // of the caller's initial guess.  That is the opposite of what a
        // bracketing solver must offer, and it is exactly the shape any warm
        // start or per-marker cache would take: re-solving each converged
        // problem with init set to the root it just returned failed on 8641 of
        // 79859 production-shaped problems (10.8 %) before this guard and on
        // none of them after.  `warm_start_from_the_returned_root_never_fails`
        // in tests/spa_solver_test.cpp is that experiment, at test scale.
        //
        // Two situations hide under that one condition, and they get different
        // answers:
        //
        //   * The residual at `init` already meets the convergence criterion.
        //     Then `init` is an answer by the solver's own stated contract, and
        //     under a locally accurate K'' the root lies within half an ulp of
        //     it, so no bracket narrower than the degenerate [init, init] is
        //     representable.  Report Converged with that bracket — the same
        //     statement the `e.k1 == 0.0` fast path above already makes, with
        //     the tolerance in place of exact zero.
        //   * The residual is still above tolerance.  That needs
        //     sqrt(K''(init)) * |init| * eps > rtol, so it takes an
        //     astronomically ill-scaled problem, but `init` is then genuinely
        //     not an answer and the expansion must proceed.  Fall back to the
        //     coarse distance, which is at least max(1, |init|) and therefore
        //     always moves.
        //
        // Gating BOTH branches on "the probe would not move", rather than
        // testing the tolerance unconditionally on entry, is deliberate and is
        // NOT merely conservatism.  Measured over examples/baseline.sh, 830 of
        // 8288 production solves (10.0 %, all of them SPAGRM) already satisfy
        // the convergence criterion at `init`; an unconditional early return
        // would take every one of those out of the Newton loop and move its
        // p-value in the last bits.  The collapse condition is reached 0 times
        // in the same 8288, so gating on it makes this repair provably inert on
        // every input the regression suite covers and a repair only where the
        // predecessor manufactured an NA.
        if (t + dir * step == t) {
            const double tol0 = std::fmax(
                (e.k2 > 0.0 && std::isfinite(e.k2)) ? opt.rtol * std::sqrt(e.k2)
                                                    : 0.0,
                noiseFloor);
            if (std::fabs(e.k1) <= tol0) {
                lo = hi = t;
                out.bracketed = true;
                out.status = Status::Converged;
            } else {
                step = coarse;
            }
        }

        for (int k = 0; k < opt.maxExpand && !out.bracketed; ++k) {
            double b = a + dir * step;
            if (!std::isfinite(b)) break;
            if (b < L) b = L;
            if (b > U) b = U;
            // Either the sign constraint's endpoint has been reached, or the
            // non-finite retreat below has halved `step` down past the local
            // ulp.  Both mean the admissible half-line is exhausted.  The first
            // probe cannot arrive here: the guard above floors it at a distance
            // that moves.
            if (b == a) break;

            const K12 eb = ev.cheap(b);
            ++out.nEval;

            if (!std::isfinite(eb.k1)) {
                // Probed outside the numerically representable region; come
                // back in rather than abandoning the expansion.
                step *= 0.5;
                if (!(step > 0.0)) break;
                continue;
            }

            if (std::fabs(eb.k1) < bestF) { bestF = std::fabs(eb.k1); bestT = b; }

            if (eb.k1 == 0.0) {
                t = b;
                e = eb;
                lo = hi = b;
                out.bracketed = true;
                out.status = Status::Converged;
                break;
            }

            if ((eb.k1 > 0.0) != (ea.k1 > 0.0)) {
                if (dir > 0.0) { lo = a; hi = b; }
                else           { lo = b; hi = a; }
                // Enter the loop at whichever endpoint has the smaller
                // |residual|.  Both are already evaluated, so the choice is
                // free, and starting from the better iterate is what lets
                // Newton converge without a safeguard step: the probe point b
                // overshoots the root by construction, so when the overshoot is
                // the larger error the correct starting point is a.
                if (std::fabs(ea.k1) < std::fabs(eb.k1)) {
                    t = a;
                    e = ea;
                } else {
                    t = b;
                    e = eb;
                }
                out.bracketed = true;
                break;
            }

            a = b;
            ea = eb;
            t = b;
            e = eb;
            // First usable probe that failed to bracket: jump straight to the
            // coarse constant rather than doubling a possibly-tiny Newton step.
            // Afterwards, grow geometrically as before.  This is what keeps the
            // worst case at one evaluation more than the predecessor's rather
            // than a long geometric walk up from a tiny first step.
            step *= opt.bracketGrow;
            if (coarseUnused) {
                coarseUnused = false;
                if (coarse > step) step = coarse;
            }
        }

        if (!out.bracketed) {
            // No sign change inside the admissible interval.  Either the
            // problem has no root there (K' is bounded for a
            // purely-outlier binomial CGF, so an |s| outside its range is
            // genuinely unreachable) or the sign constraint excluded it.
            // Both are reported as a budget exhaustion: no probability is
            // produced either way.
            out.status = Status::MaxIter;
            t = bestT;
        }
    }

    out.nEvalBracket = out.nEval;

    // ── Safeguarded Newton inside the bracket ──
    if (out.bracketed && out.status != Status::Converged) {
        Status st = Status::MaxIter;

        // Progress reference for the step-stall safeguard below.  Seeded with
        // the full bracket width so the first iteration is unconditionally
        // eligible for a Newton step, and thereafter carrying the length of the
        // last step taken.
        double prevStep = hi - lo;

        for (int iter = 0; iter < opt.maxIter; ++iter) {
            ++out.iters;
            trace(iter, t, e.k1, e.k2, lo, hi);

            if (lo == hi) { st = Status::Converged; break; }

            // ── The step-stall safeguard ──
            //
            // A Newton step is accepted only when it is at most half the length
            // of the previous step.  This is the classical safeguard (Press et
            // al., `rtsafe`), and it is what bounds the iteration count: either
            // the steps contract geometrically by a factor of two or better — in
            // which case the iterate reaches the root in O(log2(width / tol))
            // steps on its own — or a bisection is forced and the BRACKET
            // contracts by two, giving the same bound.  Newton on a convex CGF
            // converges quadratically, so its steps contract far faster than
            // required and the safeguard costs nothing on ordinary input.  On
            // the Poisson CGF, where Newton creeps toward the root by a fixed
            // decrement of one per step (t_{n+1} = t_n - 1 + s e^{-t_n}) while
            // reducing |residual| every time, the steps do not contract at all
            // and the safeguard fires — which is the case Stage 1 introduced a
            // safeguard for.
            //
            // STAGE 4 REWORK.  The predecessor tested the BRACKET width instead:
            // it forced a bisection whenever the bracket had not halved since the
            // previous iteration.  That trigger is far too eager, because Newton
            // from a bracketed start approaches the root monotonically from one
            // side, so it moves only ONE endpoint and the width shrinks slowly
            // even while the iterate converges quadratically.  A root at 0.5 in
            // [0, 1] reached by the Newton sequence 0.9, 0.7, 0.55, 0.501 gives
            // width ratios 0.78, 0.79, 0.91 — never 0.5 — so a bisection was
            // forced on essentially every other iteration and the solver
            // converged linearly on problems where Newton alone converges
            // quadratically.  Measured on the 50000 x 20000 synthetic cohort:
            // 6.23 of SPAsqr's 12.71 loop iterations per tail were forced
            // bisections.
            bool bisected = true;
            double tn = 0.0;
            if (e.k2 > 0.0 && std::isfinite(e.k2)) {
                const double cand = t - e.k1 / e.k2;
                const double len = std::fabs(cand - t);
                if (std::isfinite(cand) && cand > lo && cand < hi &&
                    len <= 0.5 * prevStep) {
                    tn = cand;
                    bisected = false;
                    prevStep = len;
                }
            }
            if (bisected) {
                const double half = 0.5 * (hi - lo);
                tn = lo + half;
                prevStep = half;
            }
            if (!(tn > lo && tn < hi)) {
                // The midpoint itself is not strictly interior: the bracket
                // is down to adjacent representable doubles.
                st = Status::Converged;
                break;
            }

            K12 en = ev.cheap(tn);
            ++out.nEval;

            if (!std::isfinite(en.k1)) {
                if (!bisected) {
                    const double half = 0.5 * (hi - lo);
                    tn = lo + half;
                    prevStep = half;
                    en = ev.cheap(tn);
                    ++out.nEval;
                    bisected = true;
                }
                if (!std::isfinite(en.k1)) {
                    t = tn;
                    e = en;
                    st = Status::NonFinite;
                    break;
                }
            }

            if (!bisected && !(std::fabs(en.k1) < std::fabs(e.k1))) {
                // The Newton step did not reduce |residual|.  It still
                // tightens the bracket, so fold it in and then bisect the
                // smaller interval.
                if (en.k1 > 0.0) hi = tn; else lo = tn;
                const double half = 0.5 * (hi - lo);
                const double tb = lo + half;
                if (tb > lo && tb < hi) {
                    const K12 eb = ev.cheap(tb);
                    ++out.nEval;
                    if (std::isfinite(eb.k1)) {
                        tn = tb;
                        en = eb;
                        prevStep = half;
                    }
                }
            }

            // Accept the point and tighten the bracket with it.
            if (en.k1 == 0.0)      lo = hi = tn;
            else if (en.k1 > 0.0)  hi = tn;
            else                   lo = tn;

            t = tn;
            e = en;

            // Test after applying the step, never before.
            const double tol = std::fmax(
                (e.k2 > 0.0 && std::isfinite(e.k2)) ? opt.rtol * std::sqrt(e.k2) : 0.0,
                noiseFloor);
            if (std::fabs(e.k1) <= tol) { st = Status::Converged; break; }
            if (hi - lo <= opt.stepTol * std::fmax(1.0, std::fabs(t))) {
                st = Status::Converged;
                break;
            }
        }
        out.status = st;
    }

    // ── Terminal cumulants, at the returned root ──
    // Unconditional, including on the failure paths, so that the Saddle never
    // reports a cumulant taken from a different abscissa than `zeta`.  This is
    // the structural close of SPACox's stale-K'' hole: its Family A finder
    // returned the post-step abscissa when the iteration budget ran out while
    // handing back the pre-step K''.
    out.zeta = t;
    out.lo = lo;
    out.hi = hi;
    out.reusedTerminal = ev.willReuse(t);

    const K012 term = ev.full(t);
    out.K0 = term.k0;
    out.K2 = term.k2;
    out.residual = term.k1;

    // The terminal guards refine a *successful* solve only.  When the root was
    // never located, the reason the root was never located is the more
    // informative report: a solver that walked out to t = 1e19 looking for a
    // sign change that does not exist will find K'' = 0 there, and reporting
    // GUARD_CURV would name the symptom instead of the cause.  Either way the
    // caller reports NA, so nothing is lost by preferring the cause.
    if (out.status == Status::Converged) {
        if (!std::isfinite(t) || !std::isfinite(term.k0) ||
            !std::isfinite(term.k1) || !std::isfinite(term.k2)) {
            out.status = Status::NonFinite;
        } else if (!(term.k2 > 0.0)) {
            out.status = Status::GuardCurv;
        }
    }

    return out;
}

} // namespace detail

// ──────────────────────────────────────────────────────────────────────
// § 4a  Solver entry points
// ──────────────────────────────────────────────────────────────────────
//
// eval(t)     -> spa::K12  {K'(t) - s, K''(t)}
// evalFull(t) -> spa::K012 {K(t), K'(t) - s, K''(t)}
//
// Both must be inlinable callables — a lambda, a function object, a function
// pointer.  Never std::function: see the boundary rule at the top of this
// header.

// Two-callable form, with instrumentation.  This is the cost-optimal shape:
// the loop never computes a logarithm, and `evalFull` runs exactly once.
template <class Eval, class EvalFull, class Trace>
Saddle solveSaddlepointTraced(
    double s,
    Eval &&eval,
    EvalFull &&evalFull,
    const SolveOpts &opt,
    Trace &&trace
) {
    using EvalT = typename std::remove_reference<Eval>::type;
    using FullT = typename std::remove_reference<EvalFull>::type;

    struct Wrapper {
        detail::SplitFull<EvalT, FullT> inner;
        K12  cheap(double t) { return inner.cheap(t); }
        K012 full(double t) { return inner.full(t); }
        bool willReuse(double) const noexcept { return false; }
    } w{detail::SplitFull<EvalT, FullT>(eval, evalFull)};

    return detail::solveCore(s, w, opt, trace);
}

template <class Eval, class EvalFull>
Saddle solveSaddlepoint(
    double s,
    Eval &&eval,
    EvalFull &&evalFull,
    const SolveOpts &opt
) {
    NullTrace tr;
    return solveSaddlepointTraced(s, static_cast<Eval &&>(eval),
                                  static_cast<EvalFull &&>(evalFull), opt, tr);
}

// Single-callable form, for a CGF that yields K0 at no extra cost (SPACox's
// empirical table, or any caller that would rather keep one code path).  The
// terminal evaluation is skipped when the loop's last abscissa already equals
// the returned root.
template <class EvalFull, class Trace>
Saddle solveSaddlepointTraced(
    double s,
    EvalFull &&evalFull,
    const SolveOpts &opt,
    Trace &&trace
) {
    using FullT = typename std::remove_reference<EvalFull>::type;

    struct Wrapper {
        detail::CachedFull<FullT> inner;
        K12  cheap(double t) { return inner.cheap(t); }
        K012 full(double t) { return inner.full(t); }
        bool willReuse(double t) const noexcept {
            return inner.valid && t == inner.last;
        }
    } w{detail::CachedFull<FullT>(evalFull)};

    return detail::solveCore(s, w, opt, trace);
}

template <class EvalFull>
Saddle solveSaddlepoint(double s, EvalFull &&evalFull, const SolveOpts &opt) {
    NullTrace tr;
    return solveSaddlepointTraced(s, static_cast<EvalFull &&>(evalFull), opt, tr);
}

// ──────────────────────────────────────────────────────────────────────
// § 5  Barndorff-Nielsen modified signed root
// ──────────────────────────────────────────────────────────────────────

// The removable singularity of r* = w + log(v/w)/w at w -> 0.
//
// Writing rho = v/w - 1, both v and w carry ~1 ULP relative error, so rho
// carries an *absolute* error of ~2 eps and rho/w an absolute error of
// ~2 eps/|w|.  Below |w| ~ 1e-5 the correction term is noise rather than a
// correction, and the correct degradation is Phi(+/-w): the dropped term is
// O(zeta) there whereas the computed term is O(eps/zeta) noise.
//
// STAGE 1 sets the threshold to zero, reproducing the `w == 0.0` guard that
// spamix, spacox and wtcoxg already carried and supplying it to spagrm and
// spasqr, which had none.  STAGE 8 (spa_unify N3) raises it to ~1e-4 and adds
// the Phi(+/-w) fallback in place of the NaN.  Keeping it as a named constant
// is what makes that a local change.
//
// WARNING TO STAGE 8 — raising this constant is not numerically neutral above
// the threshold, and the interaction is easy to miss because it lives in a
// different file.  util/spa_cgf.hpp's terminal-K section justifies the cheap,
// vectorized spelling of K on the grounds that only the ABSOLUTE error in K
// propagates, since dw = -dK/w and the |z| > spaCutoff entry gate keeps |w| of
// order 2 or more (min |w| measured at 1.779 over the production grid).  That
// argument is scale-dependent.  Perturbing K alone moves r* by
//
//     dr*/dK = -(1/w) * d/dw [ w + log(v/w)/w ]
//            = -(1/w) * [ 1 - (1 + log(v/w))/w^2 ]  ~  1/w^3   as w -> 0,
//
// because log(v/w) -> 0 at the removable singularity.  Immediately above a
// threshold of 1e-4 that is an amplification of 1e12: an absolute error of 1e-11
// in K, comfortably inside every budget spa_cgf defends, becomes an error of
// order 10 in r*, and 1e-10 becomes order 100.  So the region 1e-4 < |w| < ~1e-2
// is not a region in which r* may simply be evaluated as written.  Stage 8 must
// choose the threshold and the width of the Phi(+/-w) degradation band together,
// with dK/w^3 in view, rather than treating the guard as a pure NaN-avoidance
// device.  This constant is 0.0 today, so nothing here is presently exposed.
constexpr double kWSingularity = 0.0;

namespace detail {

// Guards and r*, shared by the linear and log tails so the two cannot drift.
// Returns r*, or NaN with `st` naming the guard that fired.
//
// K0 and K2 are treated as opaque scalars.  In particular this routine never
// assumes K2 is the second derivative of K0: WtCoxG legitimately supplies an
// inconsistent pair — its k0_total is built from var_n_K01 (a batch-effect
// Gaussian variance) while its k2_total is built from var_n_K2 (a
// finite-reference-panel correction), documented as intentional at
// wtcoxg.cpp:430-435 — and a kernel that silently "corrected" that would
// change WtCoxG's statistic rather than compute it.
inline double modifiedSignedRoot(
    double zeta,
    double s,
    double K0,
    double K2,
    Status &st
) noexcept {
    if (!std::isfinite(zeta) || !std::isfinite(s) || !std::isfinite(K0) ||
        !std::isfinite(K2)) {
        st = Status::NonFinite;
        return quietNaN();
    }
    if (!(K2 > 0.0)) {
        st = Status::GuardCurv;
        return quietNaN();
    }

    const double temp = zeta * s - K0;
    if (!(temp >= 0.0)) {
        st = Status::GuardTemp;
        return quietNaN();
    }

    const double w = signOf(zeta) * std::sqrt(2.0 * temp);
    if (!(std::fabs(w) > kWSingularity)) {
        st = Status::GuardW;
        return quietNaN();
    }

    const double v = zeta * std::sqrt(K2);
    const double vOverW = v / w;
    if (!(vOverW > 0.0) || !std::isfinite(vOverW)) {
        // Unreachable given the guards above — v and w both carry sgn(zeta),
        // so their ratio is positive whenever K2 > 0 and temp > 0.  Retained
        // as a belt: the ratio is the one place an inconsistent (K0, K2) pair
        // could in principle misbehave.
        st = Status::NonFinite;
        return quietNaN();
    }

    const double rStar = w + std::log(vOverW) / w;
    if (!std::isfinite(rStar)) {
        st = Status::NonFinite;
        return quietNaN();
    }

    st = Status::Converged;
    return rStar;
}

} // namespace detail

// Tail probability from the modified signed root.
//   lowerTail == true   -> Phi(r*)      approximates P(S <= s)
//   lowerTail == false  -> 1 - Phi(r*)  approximates P(S >= s)
//
// Scalars rather than a struct by deliberate choice: under the SysV AMD64 ABI
// four doubles pass in xmm0-xmm3 and the bool in edi, so the call is
// register-only.  A 32-byte struct by value is classified MEMORY and spills.
//
// Returns NaN and sets `st` on every degenerate input; sets
// Status::Converged (statusName "OK") when no guard fired.
inline double bnTail(
    double zeta,
    double s,
    double K0,
    double K2,
    bool lowerTail,
    Status &st
) noexcept {
    const double rStar = detail::modifiedSignedRoot(zeta, s, K0, K2, st);
    if (!std::isfinite(rStar)) return detail::quietNaN();
    return detail::phiTail(rStar, lowerTail);
}

// Natural logarithm of the same tail probability.
//
// Every original site returned pval1 + pval2 on the linear scale, where
// Phi(-x) becomes denormal at x ~ 38.0 and flushes to zero at x ~ 38.5, so an
// association with true p = 1e-400 was reported as 0.  This variant keeps the
// extreme tail (spa_unify N2, L3), and feeds the LOG10P output column.
inline double bnTailLog(
    double zeta,
    double s,
    double K0,
    double K2,
    bool lowerTail,
    Status &st
) noexcept {
    const double rStar = detail::modifiedSignedRoot(zeta, s, K0, K2, st);
    if (!std::isfinite(rStar)) return detail::quietNaN();
    return detail::phiTailLog(rStar, lowerTail);
}

// ──────────────────────────────────────────────────────────────────────
// § 6  Two-sided assembly
// ──────────────────────────────────────────────────────────────────────

// `negLog10p` is -log10(p), i.e. the quantity reported in the LOG10P column
// and therefore non-negative.  The field is named for the sign it carries so
// that a caller cannot get it backwards; the design sketch called this
// `log10p`, which reads as log10(p) and invites exactly that error.
struct TwoSided {
    double p;
    double negLog10p;
    Status status;
};

// Linear-scale assembly.  A single clamp policy for all six methods:
//
//   * failure in *either* tail yields NaN in p, never a half-sized p-value.
//     SPAmixLocalPlus added each tail only if that tail's root converged and
//     produced NaN only when both failed, so a marker whose upper tail
//     converged and whose lower tail did not was reported at approximately
//     half its correct p-value with no diagnostic.
//   * p is clamped into [0, 1].  It is *not* floored at
//     DBL_MIN: substituting 2.2e-308 for an underflowed tail manufactures a
//     genome-wide-significant hit.  A p that underflows to zero is reported
//     as zero with negLog10p = +inf; combineTailsLog is the way to keep the
//     magnitude.
//   * NaN never becomes 1.0.  WtCoxG wrote std::min(1.0, pval1 + pval2),
//     and std::min(1.0, NaN) returns 1.0, so every SPA failure surfaced as a
//     perfectly null marker.
inline TwoSided combineTails(
    double pUpper,
    double pLower,
    Status sUpper,
    Status sLower
) noexcept {
    const double nan = detail::quietNaN();
    const Status st = worseStatus(sUpper, sLower);
    if (statusIsFailure(st)) return TwoSided{nan, nan, st};
    if (!std::isfinite(pUpper) || !std::isfinite(pLower))
        return TwoSided{nan, nan, Status::NonFinite};

    double p = pUpper + pLower;
    if (p > 1.0) p = 1.0;
    if (p < 0.0) p = 0.0;
    // -log10(1.0) is -0.0, which prints as "-0"; normalize the sign of zero.
    double neg = -std::log10(p);
    if (neg == 0.0) neg = 0.0;
    return TwoSided{p, neg, st};
}

// Log-domain assembly, via log-sum-exp so that the sum of two underflowed
// tails does not re-underflow.  Inputs are natural logarithms, as returned by
// bnTailLog.
inline TwoSided combineTailsLog(
    double logPUpper,
    double logPLower,
    Status sUpper,
    Status sLower
) noexcept {
    const double nan = detail::quietNaN();
    const double inf = std::numeric_limits<double>::infinity();
    const Status st = worseStatus(sUpper, sLower);
    if (statusIsFailure(st)) return TwoSided{nan, nan, st};
    if (std::isnan(logPUpper) || std::isnan(logPLower))
        return TwoSided{nan, nan, Status::NonFinite};

    const double hi = std::fmax(logPUpper, logPLower);
    const double lo = std::fmin(logPUpper, logPLower);
    double lse;
    if (hi == -inf)      lse = -inf;               // both tails exactly zero
    else if (hi == inf)  lse = inf;
    else                 lse = hi + std::log1p(std::exp(lo - hi));

    if (lse > 0.0) lse = 0.0;                      // p clamped at 1
    double p = std::exp(lse);
    if (p > 1.0) p = 1.0;
    double neg = -lse / detail::kLn10;
    if (neg == 0.0) neg = 0.0;
    return TwoSided{p, neg, st};
}

// ──────────────────────────────────────────────────────────────────────
// § 7  The normal short-circuit
// ──────────────────────────────────────────────────────────────────────

// Two-sided normal p-value, for the |z| <= spaCutoff branch where the
// saddlepoint is never attempted.  Routed through Boost's complement rather
// than 1 - Phi(|z|), which floors at zero for |z| > 8.3.
inline double normalTwoSided(double z) noexcept {
    return 2.0 * math::pnorm(std::fabs(z), 0.0, 1.0, /*lower_tail=*/false);
}

// log of the same, valid to |z| well past the linear-scale underflow.
inline double normalTwoSidedLog(double z) noexcept {
    if (std::isnan(z)) return detail::quietNaN();
    return detail::kLn2 + detail::logPhi(-std::fabs(z));
}

// The complete normal-branch result, so that a method emitting P, LOG10P and
// SPA_STATUS has one call for the non-SPA path.
inline TwoSided normalBranch(double z) noexcept {
    const double nan = detail::quietNaN();
    if (std::isnan(z)) return TwoSided{nan, nan, Status::NonFinite};
    double p = normalTwoSided(z);
    if (p > 1.0) p = 1.0;
    if (p < 0.0) p = 0.0;
    const double lg = normalTwoSidedLog(z);
    double neg = -lg / detail::kLn10;
    if (neg == 0.0) neg = 0.0;
    return TwoSided{p, neg, Status::NormalBranch};
}

} // namespace spa
