// spagrm_ibd_test.cpp — the SPAGRM relatedness inputs, end to end
//
// SPAGRM's family CGF is built from two files it does not produce itself in
// the same run: a sparse GRM and a pairwise-IBD table.  Everything the CGF
// does with them presupposes that (pa, pb, pc) is a probability vector, and
// nothing downstream can tell you when it is not:
//
//   * nsGRMNull::buildChowLiuTree mixes three conditional joint tables, each
//     with marginal p0, using pa, pb and pc as weights.  The family table
//     therefore sums to the product of (pa+pb+pc) over the spanning-tree
//     edges, and spagrm_cgf's class-3 block contributes log of that to K.  A
//     constant in K leaves zeta, K' and K'' exactly unchanged and moves only
//     w = sgn(zeta)*sqrt(2*(zeta*s - K(zeta))): every SPA-branch p-value is
//     inflated by a factor approaching exp(K(0)), and once K(0) exceeds the
//     large-deviation rate the marker is reported NA with GUARD_TEMP.
//     GUARD_TEMP cannot arise from a genuine CGF — zeta*s - K(zeta) is the
//     Legendre transform and is non-negative by construction — so it is a
//     proof of bad input arriving far too late and attached to the wrong
//     subsystem.
//   * the two-subject block factorizes the pair through the roots of
//     x^2 - 2*(pa + pb/2)*x + pa.  Those leave [0, 1] as soon as
//     pa + pb/2 > 1, the four-point mixture acquires negative point masses,
//     and K'' can go non-positive.
//
// The suite pins the three places that must hold: the estimator produces an
// admissible vector for every input, the loader refuses one that is not, and
// the Chow-Liu table it feeds is a probability distribution with the right
// marginals.  It also pins the quadratic-form unification that removed
// SparseGRM::spaVariance.
//
// Links: grm_null.o spagrm.o io/sparse_grm.o util/{spa_cgf,math_helper,
// text_stream}.o plus the compression objects text_stream needs — see the
// TESTOBJS_spagrm_ibd_test list in the Makefile.

#include "tinytest.hpp"

#include "io/sparse_grm.hpp"
#include "spagrm/grm_null.hpp"
#include "spagrm/ibd.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

// ──────────────────────────────────────────────────────────────────────
// Shared helpers
// ──────────────────────────────────────────────────────────────────────

namespace {

// The GRM off-diagonals a real cohort produces, plus the inadmissible ones
// the bundled 1kg fixture actually contains (its 173 offending pairs span
// 1.0805 to 1.6592) and the degenerate corners the clamp arithmetic has to
// survive.
const std::vector<double> kGrmGrid = {
    -0.5, -1e-12, 0.0, 1e-12, 1e-9, 1e-6, 1e-5, 0.05, 0.125, 0.2,
    0.25, 0.35,   0.5, 0.7,   0.9,  0.99, 1.0,  1.0805, 1.3, 1.6592257, 3.0
};

// rho2 is a weighted mean of non-negative per-marker terms, so it is >= 0 in
// practice; values above 1 and the pathological 0 are included because the
// estimator must not depend on the moment estimate landing anywhere sensible.
const std::vector<double> kPcGrid = {
    0.0, 0.05, 0.1, 0.25, 0.4, 0.5, 0.6, 0.75, 0.9, 0.99, 1.0, 1.5, 4.0
};

// The per-allele sharing probabilities the two-subject CGF factorizes
// through: the roots of x^2 - 2*Rho*x + pa with Rho = pa + pb/2, formed
// exactly as nsGRMNull::buildSPAGRMNullModel forms them.
struct Roots {
    double lo, hi, disc;
};

Roots perAlleleRoots(const ibd::Triple &t) {
    const double rho = t.pa + 0.5 * t.pb;
    const double disc = rho * rho - t.pa;
    const double mid = std::sqrt(std::max(disc, 0.0));
    return Roots{rho - mid, rho + mid, disc};
}

// Write `body` to `path` and return the path, so a test can hand a file to
// loadIndexedIBD.  `make test` runs from the repository root, where the
// Makefile has already made tmp/ an order-only prerequisite of every test
// binary — but the binary is also runnable by hand from anywhere, and a
// missing directory would then surface as six unrelated assertion failures
// rather than as the one real problem.  Create it, and let REQUIRE report a
// genuinely unwritable location.
std::string writeTempFile(const char *name, const std::string &body) {
    std::error_code ec;
    std::filesystem::create_directories("tmp", ec);
    const std::string path = std::string("tmp/") + name;
    std::FILE *f = std::fopen(path.c_str(), "w");
    if (!f) return {};
    std::fwrite(body.data(), 1, body.size(), f);
    std::fclose(f);
    return path;
}

const std::unordered_map<std::string, uint32_t> &twoSubjectIdMap() {
    static const std::unordered_map<std::string, uint32_t> m = {{"A", 0}, {"B", 1}};
    return m;
}

bool loadThrows(const std::string &path) {
    try {
        nsGRMNull::loadIndexedIBD(path, twoSubjectIdMap());
    } catch (const std::exception &) {
        return true;
    }
    return false;
}

} // namespace

// ──────────────────────────────────────────────────────────────────────
// The estimator:  ibd::deriveIBD
// ──────────────────────────────────────────────────────────────────────

// Eq. (17)'s third equation.  This is the property the whole file exists to
// protect, and it must hold for every input, not merely for plausible ones.
TEST(derive_ibd_sums_to_one_everywhere) {
    for (double v : kGrmGrid) {
        for (double pc : kPcGrid) {
            const ibd::Triple t = ibd::deriveIBD(v, pc);
            CHECK_NEAR(t.pa + t.pb + t.pc, 1.0, 1e-15);
        }
    }
}

// Non-negativity, which Eq. (17) does not give for free and which the three
// admissibility constraints are there to supply.
TEST(derive_ibd_is_a_probability_vector_everywhere) {
    for (double v : kGrmGrid) {
        for (double pc : kPcGrid) {
            const ibd::Triple t = ibd::deriveIBD(v, pc);
            CHECK(t.pa >= -1e-15 && t.pa <= 1.0 + 1e-15);
            CHECK(t.pb >= -1e-15 && t.pb <= 1.0 + 1e-15);
            CHECK(t.pc >= -1e-15 && t.pc <= 1.0 + 1e-15);
        }
    }
}

// Constraint (ii) is exactly the condition that the two-subject
// factorization has real roots, and (iii) exactly the condition that the
// larger one does not exceed 1.  With all three imposed, the discriminant is
// non-negative on its own — buildSPAGRMNullModel's std::max(..., 0.0) guard
// never has to fire — and both roots are probabilities.
TEST(derive_ibd_per_allele_roots_are_probabilities) {
    for (double v : kGrmGrid) {
        for (double pc : kPcGrid) {
            const Roots r = perAlleleRoots(ibd::deriveIBD(v, pc));
            CHECK(r.disc >= -1e-15);
            CHECK(r.lo >= -1e-15 && r.lo <= 1.0 + 1e-15);
            CHECK(r.hi >= -1e-15 && r.hi <= 1.0 + 1e-15);
        }
    }
}

// A negative root or a root above 1 makes tj = (1-rho)*MAF*(1-MAF) negative
// in spagrm_cgf's four-point mixture, i.e. a signed measure rather than a
// distribution.  Stated separately from the test above because this is the
// consequence, and it is what a reader of a GUARD_CURV row needs to know.
TEST(derive_ibd_four_point_weights_are_non_negative) {
    for (double maf : {0.001, 0.01, 0.1, 0.3, 0.5}) {
        const double maf01 = maf * (1.0 - maf);
        for (double v : kGrmGrid) {
            for (double pc : kPcGrid) {
                const Roots r = perAlleleRoots(ibd::deriveIBD(v, pc));
                for (double rho : {r.lo, r.hi}) {
                    const double tj = (1.0 - rho) * maf01;
                    CHECK(tj >= -1e-18);              // weight on R1 and on R2
                    CHECK(maf - tj >= -1e-18);        // weight on R1+R2
                    CHECK(1.0 - maf - tj >= -1e-18);  // weight on 0
                }
            }
        }
    }
}

// The textbook relationships, each entered through the moment estimates the
// estimator actually receives: rho1 = 2*phi from the GRM and rho2 = k0 from
// the IBS0 statistic.
TEST(derive_ibd_reproduces_known_relationships) {
    struct Case {
        const char *name;
        double v, pcRaw, pa, pb, pc;
    };
    const Case cases[] = {
        // name                 rho1   rho2   k2    k1    k0
        {"unrelated",           0.00,  1.00,  0.00, 0.00, 1.00},
        {"first cousins",       0.125, 0.75,  0.00, 0.25, 0.75},
        {"half sibs",           0.25,  0.50,  0.00, 0.50, 0.50},
        {"parent-offspring",    0.50,  0.00,  0.00, 1.00, 0.00},
        {"full sibs",           0.50,  0.25,  0.25, 0.50, 0.25},
        {"monozygotic",         1.00,  0.00,  1.00, 0.00, 0.00},
    };
    for (const Case &c : cases) {
        const ibd::Triple t = ibd::deriveIBD(c.v, c.pcRaw);
        CHECK_NEAR(t.pa, c.pa, 1e-9);
        CHECK_NEAR(t.pb, c.pb, 1e-9);
        CHECK_NEAR(t.pc, c.pc, 1e-9);
    }
}

// Constraint (iii).  The 173 offending pairs of the bundled fixture all land
// here.  The predecessor's repair branch produced (v, 0, 0) — pa above 1 and
// a vector summing to v; the cap produces the duplicate vector instead.
TEST(derive_ibd_caps_above_the_monozygotic_bound) {
    for (double v : {1.0 + 1e-12, 1.0805, 1.3, 1.6592257, 3.0, 1e6}) {
        for (double pc : kPcGrid) {
            const ibd::Triple t = ibd::deriveIBD(v, pc);
            CHECK_NEAR(t.pa, 1.0, 1e-9);
            CHECK_NEAR(t.pb, 0.0, 1e-9);
            CHECK_NEAR(t.pc, 0.0, 1e-9);
            // Specifically NOT the predecessor's answer.
            CHECK(t.pa <= 1.0 + 1e-15);
        }
    }
}

// std::clamp is UNDEFINED when hi < lo, so the band's endpoints must be
// ordered for every representable v — not merely in exact arithmetic.  The
// predecessor's (1-v)^2 - 1e-10 inverted them for every v below 1e-5.  Simply
// dropping the epsilon is not enough either: evaluating (1-v)*(1-v) and 1-2v
// independently rounds each, so the difference is v^2 + O(2^-53) and the
// rounding term dominates below v ~ 1.3e-8.  A sweep of 1e-30 <= v < 1 finds
// 211 490 abscissae where the independent form inverts.  This pins the
// ordering directly, at the same abscissae, rather than inferring it from the
// output of a call that would already have been undefined.
TEST(derive_ibd_clamp_interval_is_ordered_for_every_v) {
    long inversions = 0;
    double worstV = 0.0;
    for (double v = 1e-30; v < 1.0; v *= 1.0000231) {
        const double lower = 1.0 - 2.0 * v;
        const double upper = lower + v * v;   // the shipped spelling
        if (upper < lower) { ++inversions; worstV = v; }
    }
    CHECK(inversions == 0);
    (void)worstV;

    // The independently-evaluated spelling really does invert, so the choice
    // above is load-bearing and cannot be "simplified" back.
    long naiveInversions = 0;
    for (double v = 1e-30; v < 1.0; v *= 1.0000231) {
        if ((1.0 - v) * (1.0 - v) < 1.0 - 2.0 * v) ++naiveInversions;
    }
    CHECK(naiveInversions > 0);

    // And the two spellings agree wherever both are well ordered.
    double worstRel = 0.0;
    for (double v = 1e-6; v < 1.0; v *= 1.000131) {
        const double a = (1.0 - v) * (1.0 - v);
        const double b = (1.0 - 2.0 * v) + v * v;
        const double d = std::abs(a - b) / a;
        if (d > worstRel) worstRel = d;
    }
    CHECK(worstRel < 1e-15);
}

// The admissible band for rho2 is [1-2v, (1-v)^2], whose width is exactly
// v^2 and which therefore collapses to the single point rho2 = 1 at v = 0.
TEST(derive_ibd_survives_the_degenerate_clamp_interval) {
    for (double v : {0.0, 1e-12, 1e-9, 1e-7, 1e-6, 9.9e-6}) {
        for (double pc : kPcGrid) {
            const ibd::Triple t = ibd::deriveIBD(v, pc);
            CHECK_NEAR(t.pa + t.pb + t.pc, 1.0, 1e-15);
            CHECK(t.pa >= 0.0 && t.pb >= 0.0 && t.pc >= 0.0);
        }
    }
    // At v = 0 the band is the single point rho2 = 1, so every rho2 maps to
    // the unrelated vector exactly.
    for (double pc : kPcGrid) {
        const ibd::Triple t = ibd::deriveIBD(0.0, pc);
        CHECK(t.pa == 0.0);
        CHECK(t.pb == 0.0);
        CHECK(t.pc == 1.0);
    }
}

// Below the clamps the estimator must be Eq. (17) itself, unmodified: the
// constraints shape the input, they do not replace the solution.
TEST(derive_ibd_is_equation_17_in_the_interior) {
    // rho2 strictly inside (1-2v, (1-v)^2).  That band has width exactly
    // v^2, so the interior cases have to be chosen against it rather than
    // guessed: for v = 0.125 it is (0.750, 0.765625), and for v = 0.25 it is
    // (0.500, 0.562500).
    struct Case { double v, pc; };
    const Case cases[] = {{0.125, 0.76}, {0.25, 0.55}, {0.4, 0.30}, {0.5, 0.10}};
    for (const Case &c : cases) {
        const ibd::Triple t = ibd::deriveIBD(c.v, c.pc);
        CHECK_NEAR(t.pc, c.pc, 1e-15);                        // delta0 = rho2
        CHECK_NEAR(t.pb, 2.0 * (1.0 - c.v - c.pc), 1e-15);    // delta1
        CHECK_NEAR(t.pa, 2.0 * c.v + c.pc - 1.0, 1e-15);      // delta2
        CHECK_NEAR(t.pa + 0.5 * t.pb, c.v, 1e-15);            // rho1 recovered
    }
}

// ──────────────────────────────────────────────────────────────────────
// The loader:  nsGRMNull::loadIndexedIBD
// ──────────────────────────────────────────────────────────────────────

TEST(load_indexed_ibd_accepts_a_valid_file) {
    const std::string path = writeTempFile(
        "spagrm_ibd_ok.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tB\t0.25\t0.5\t0.25\n");
    REQUIRE(!path.empty());
    const auto rows = nsGRMNull::loadIndexedIBD(path, twoSubjectIdMap());
    REQUIRE(rows.size() == 1);
    CHECK_NEAR(rows[0].pa, 0.25, 1e-15);
    CHECK_NEAR(rows[0].pb, 0.50, 1e-15);
    CHECK_NEAR(rows[0].pc, 0.25, 1e-15);
}

// The exact row shape the predecessor emitted for all 173 fixture pairs.
TEST(load_indexed_ibd_rejects_the_predecessors_repaired_row) {
    const std::string path = writeTempFile(
        "spagrm_ibd_sum.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tB\t1.6592256999999999\t0\t0\n");
    REQUIRE(!path.empty());
    CHECK(loadThrows(path));
}

TEST(load_indexed_ibd_rejects_a_negative_component) {
    const std::string path = writeTempFile(
        "spagrm_ibd_neg.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tB\t1.4\t-0.4\t0\n");
    REQUIRE(!path.empty());
    CHECK(loadThrows(path));
}

// strtod parses "nan" and "inf" happily, and every consumer would propagate
// them into a quadratic form or a CGF without comment.
TEST(load_indexed_ibd_rejects_non_finite_values) {
    for (const char *tok : {"nan", "inf", "-inf"}) {
        const std::string path = writeTempFile(
            "spagrm_ibd_nonfinite.txt",
            std::string("ID1\tID2\tpa\tpb\tpc\nA\tB\t0.25\t0.5\t") + tok + "\n");
        REQUIRE(!path.empty());
        CHECK(loadThrows(path));
    }
}

// Condition (ii) of the admissible set.  Non-negativity and unit sum do NOT
// imply it, so a triple can pass every other check and still be
// undecomposable into two independent allele-sharing indicators.  Without
// this check buildSPAGRMNullModel's std::max(Rho*Rho - pa, 0.0) silently
// projects it onto a different distribution — same kinship, different pa.
TEST(load_indexed_ibd_rejects_an_undecomposable_triple) {
    // pb^2 = 0.01 against 4*pa*pc = 4*0.45*0.45 = 0.81.  Non-negative, sums
    // to one, and inadmissible: no (rho+, rho-) in [0,1]^2 produces it.
    const std::string path = writeTempFile(
        "spagrm_ibd_undecomposable.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tB\t0.45\t0.10\t0.45\n");
    REQUIRE(!path.empty());
    CHECK(loadThrows(path));
}

// The textbook relationships sit ON the boundary pb^2 = 4*pa*pc, so the
// inequality must not be strict and the tolerance must not reject them.
TEST(load_indexed_ibd_accepts_the_boundary_relationships) {
    struct Case { const char *name, *row; };
    const Case cases[] = {
        {"full sibs",        "A\tB\t0.25\t0.5\t0.25\n"},
        {"parent-offspring", "A\tB\t0\t1\t0\n"},
        {"half sibs",        "A\tB\t0\t0.5\t0.5\n"},
        {"first cousins",    "A\tB\t0\t0.25\t0.75\n"},
        {"monozygotic",      "A\tB\t1\t0\t0\n"},
        {"unrelated",        "A\tB\t0\t0\t1\n"},
    };
    for (const Case &c : cases) {
        const std::string path = writeTempFile(
            "spagrm_ibd_boundary.txt",
            std::string("ID1\tID2\tpa\tpb\tpc\n") + c.row);
        REQUIRE(!path.empty());
        CHECK_MSG(!loadThrows(path), c.name);
    }
    // And everything ibd::deriveIBD can produce satisfies it by construction,
    // which is what makes the loader's check a statement about foreign files.
    for (double v : kGrmGrid) {
        for (double pc : kPcGrid) {
            const ibd::Triple t = ibd::deriveIBD(v, pc);
            CHECK(t.pb * t.pb >= 4.0 * t.pa * t.pc - 1e-12);
        }
    }
}

// strtod returns 0 without consuming input on a non-numeric field, so a
// truncated or comma-delimited row used to be accepted as (0, 0, 0).  It must
// be named as a parse failure, not as a normalization failure.
TEST(load_indexed_ibd_rejects_a_malformed_row) {
    const char *bodies[] = {
        "ID1\tID2\tpa\tpb\tpc\nA\tB\n",                  // no numeric fields
        "ID1\tID2\tpa\tpb\tpc\nA\tB\t0.25\n",            // truncated
        "ID1\tID2\tpa\tpb\tpc\nA\tB\t0.25\t0.5\n",       // truncated
        "ID1\tID2\tpa\tpb\tpc\nA\tB\t0.25,0.5,0.25\n",   // comma-delimited
    };
    for (const char *b : bodies) {
        const std::string path = writeTempFile("spagrm_ibd_malformed.txt", b);
        REQUIRE(!path.empty());
        CHECK(loadThrows(path));
    }
}

// A row whose subjects are outside the analysis set is dropped, not
// validated: it never reaches a CGF, and rejecting the run because of it
// would make an unrelated subset filter fail on data it does not use.
TEST(load_indexed_ibd_skips_unknown_subjects_without_validating_them) {
    const std::string path = writeTempFile(
        "spagrm_ibd_unknown.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tZ\t9\t9\t9\n"
        "A\tB\t0.25\t0.5\t0.25\n");
    REQUIRE(!path.empty());
    const auto rows = nsGRMNull::loadIndexedIBD(path, twoSubjectIdMap());
    CHECK(rows.size() == 1);
}

// The precision policy: exact input passes through, printed-precision input is
// renormalized, and anything the predecessor produced is rejected.  The two
// populations are separated by more than two orders of magnitude — the
// smallest deviation the old repair branch emitted was 0.0805.
TEST(load_indexed_ibd_applies_the_two_tier_precision_policy) {
    // Exact: accepted unchanged.
    const std::string exact = writeTempFile(
        "spagrm_ibd_tol_exact.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tB\t0.25\t0.5\t0.2500000001\n");
    REQUIRE(!exact.empty());
    CHECK(!loadThrows(exact));

    // A four-decimal third-party table: the paper permits other tools to
    // supply these, and a printed table cannot sum to one exactly.  Accepted
    // and renormalized, so what reaches the CGF is a distribution.
    const std::string printed = writeTempFile(
        "spagrm_ibd_tol_printed.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tB\t0.2500\t0.5001\t0.2500\n");
    REQUIRE(!printed.empty());
    CHECK(!loadThrows(printed));
    const auto rows = nsGRMNull::loadIndexedIBD(printed, twoSubjectIdMap());
    REQUIRE(rows.size() == 1);
    CHECK_NEAR(rows[0].pa + rows[0].pb + rows[0].pc, 1.0, 1e-15);

    // A component printed as a small negative — round-off around an estimate
    // that sat on the boundary — is floored to zero, and the sum is still
    // repaired to exactly one.  The floor must be applied BEFORE the
    // renormalization, or it reintroduces the deficit it just removed.
    const std::string negRound = writeTempFile(
        "spagrm_ibd_tol_neground.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tB\t-0.00002\t0.99998\t0.00002\n");
    REQUIRE(!negRound.empty());
    CHECK(!loadThrows(negRound));
    const auto nr = nsGRMNull::loadIndexedIBD(negRound, twoSubjectIdMap());
    REQUIRE(nr.size() == 1);
    CHECK(nr[0].pa >= 0.0 && nr[0].pb >= 0.0 && nr[0].pc >= 0.0);
    CHECK_NEAR(nr[0].pa + nr[0].pb + nr[0].pc, 1.0, 1e-15);

    // Beyond printed precision: rejected.
    const std::string bad = writeTempFile(
        "spagrm_ibd_tol_bad.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tB\t0.25\t0.5\t0.26\n");
    REQUIRE(!bad.empty());
    CHECK(loadThrows(bad));

    // The smallest deviation the predecessor's repair branch ever produced on
    // the bundled fixture, 1.0805 - 1.  It must not survive renormalization.
    const std::string predecessor = writeTempFile(
        "spagrm_ibd_tol_pred.txt",
        "ID1\tID2\tpa\tpb\tpc\n"
        "A\tB\t1.0804676\t0\t0\n");
    REQUIRE(!predecessor.empty());
    CHECK(loadThrows(predecessor));
}

// ──────────────────────────────────────────────────────────────────────
// The table:  nsGRMNull::buildChowLiuTree
// ──────────────────────────────────────────────────────────────────────

namespace {

std::vector<nsGRMNull::IndexedIBD> allPairs(
    const std::vector<uint32_t> &fam,
    const ibd::Triple &t
) {
    std::vector<nsGRMNull::IndexedIBD> out;
    for (size_t i = 0; i < fam.size(); ++i)
        for (size_t j = i + 1; j < fam.size(); ++j)
            out.push_back({fam[i], fam[j], t.pa, t.pb, t.pc});
    return out;
}

const std::vector<double> &mafGrid() {
    static const std::vector<double> g = {1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
    return g;
}

} // namespace

// K(0) = 0 for the class-3 block is exactly the statement that every column
// of the table sums to one.  The table is built over the whole MAF grid
// because spagrm.cpp interpolates between two adjacent columns per marker,
// so a column that is wrong anywhere is wrong for a band of MAFs.
TEST(chow_liu_columns_are_probability_distributions) {
    for (int n : {3, 4, 5}) {
        std::vector<uint32_t> fam;
        for (int i = 0; i < n; ++i) fam.push_back(static_cast<uint32_t>(i));

        // Relationships spanning the admissible set, plus the vector the
        // monozygotic cap produces, which is the degenerate corner.
        const std::vector<std::pair<double, double> > rel = {
            {0.0, 1.0}, {0.125, 0.75}, {0.25, 0.5}, {0.5, 0.25}, {0.5, 0.0}, {1.0, 0.0}
        };
        for (const auto &r : rel) {
            const ibd::Triple t = ibd::deriveIBD(r.first, r.second);
            const Eigen::MatrixXd clt =
                nsGRMNull::buildChowLiuTree(n, allPairs(fam, t), fam, mafGrid(), 1);
            REQUIRE(clt.cols() == static_cast<Eigen::Index>(mafGrid().size()));
            for (Eigen::Index c = 0; c < clt.cols(); ++c) {
                CHECK_NEAR(clt.col(c).sum(), 1.0, 1e-12);
                CHECK(clt.col(c).minCoeff() >= -1e-18);
            }
        }
    }
}

// Normalization on its own would hide a wrong marginal, so the marginals are
// pinned separately.  Each of pa_arr, pb_arr and pc_arr has both marginals
// equal to p0, so every single-subject marginal of the family table must be
// Hardy-Weinberg at that bin's MAF regardless of the relationship.
TEST(chow_liu_single_subject_marginals_are_hardy_weinberg) {
    const int n = 4;
    std::vector<uint32_t> fam = {0, 1, 2, 3};
    const ibd::Triple t = ibd::deriveIBD(0.5, 0.25); // full sibs
    const Eigen::MatrixXd clt =
        nsGRMNull::buildChowLiuTree(n, allPairs(fam, t), fam, mafGrid(), 1);

    int stride[5] = {1, 3, 9, 27, 81};
    const Eigen::Index arrSize = clt.rows();
    for (Eigen::Index c = 0; c < clt.cols(); ++c) {
        const double mu = mafGrid()[static_cast<size_t>(c)];
        const double q = 1.0 - mu;
        const double p0[3] = {q * q, 2.0 * mu * q, mu * mu};
        for (int subj = 0; subj < n; ++subj) {
            double marg[3] = {0.0, 0.0, 0.0};
            for (Eigen::Index idx = 0; idx < arrSize; ++idx)
                marg[(idx / stride[subj]) % 3] += clt(idx, c);
            for (int g = 0; g < 3; ++g)
                CHECK_NEAR(marg[g], p0[g], 1e-12);
        }
    }
}

// The unrelated vector must give the independence table exactly, which is
// the one case with a closed form to check the machinery against.
TEST(chow_liu_unrelated_family_is_the_product_of_marginals) {
    const int n = 3;
    std::vector<uint32_t> fam = {0, 1, 2};
    const ibd::Triple t = ibd::deriveIBD(0.0, 1.0); // (0, 0, 1)
    const Eigen::MatrixXd clt =
        nsGRMNull::buildChowLiuTree(n, allPairs(fam, t), fam, mafGrid(), 1);

    for (Eigen::Index c = 0; c < clt.cols(); ++c) {
        const double mu = mafGrid()[static_cast<size_t>(c)];
        const double q = 1.0 - mu;
        const double p0[3] = {q * q, 2.0 * mu * q, mu * mu};
        for (Eigen::Index idx = 0; idx < clt.rows(); ++idx) {
            const double want = p0[idx % 3] * p0[(idx / 3) % 3] * p0[(idx / 9) % 3];
            CHECK_REL(clt(idx, c), want, 1e-12);
        }
    }
}

// The build is parallelized over MAF bins; the result must not depend on how
// many workers ran, per CLAUDE.md's reproducibility constraint.
TEST(chow_liu_is_independent_of_thread_count) {
    const int n = 5;
    std::vector<uint32_t> fam = {0, 1, 2, 3, 4};
    const ibd::Triple t = ibd::deriveIBD(0.25, 0.5);
    const Eigen::MatrixXd a = nsGRMNull::buildChowLiuTree(n, allPairs(fam, t), fam, mafGrid(), 1);
    const Eigen::MatrixXd b = nsGRMNull::buildChowLiuTree(n, allPairs(fam, t), fam, mafGrid(), 8);
    REQUIRE(a.rows() == b.rows() && a.cols() == b.cols());
    for (Eigen::Index c = 0; c < a.cols(); ++c)
        for (Eigen::Index r = 0; r < a.rows(); ++r)
            CHECK(a(r, c) == b(r, c));
}

// ──────────────────────────────────────────────────────────────────────
// The quadratic form:  SparseGRM::quadForm
// ──────────────────────────────────────────────────────────────────────

// SPAGxE+ and SPAmixPlus used to reach the same quantity through
// SparseGRM::spaVariance, which evaluated 2*sum_stored(v*Ri*Rj) - sum(Ri^2)
// = R'GR + sum_i (G_ii - 1) Ri^2.  That is R'GR only for a unit diagonal,
// and the bundled fixture's diagonal spans 0.4597 to 3.8087.  WtCoxG already
// used quadForm, so the two disagreed.  spaVariance is gone; this pins the
// survivor against a dense evaluation and pins the size of what changed.
TEST(quad_form_matches_a_dense_evaluation_with_non_unit_diagonal) {
    // Lower triangle + diagonal, three subjects, deliberately non-unit
    // diagonal entries of the kind a structure-inflated GRM produces.
    std::vector<SparseGRM::Entry> e = {
        {0, 0, 0.46}, {1, 0, 0.32}, {1, 1, 1.90}, {2, 1, 0.21}, {2, 2, 3.81}
    };
    const SparseGRM g = SparseGRM::fromEntries(3, e);

    const double R[3] = {1.7, -0.9, 1.3};

    double dense[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    for (const auto &en : e) {
        dense[en.row][en.col] = en.value;
        dense[en.col][en.row] = en.value;
    }
    double want = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            want += R[i] * dense[i][j] * R[j];

    CHECK_REL(g.quadForm(R, 3), want, 1e-14);

    // And the deleted spelling really did differ, by exactly the diagonal
    // excess — recorded so the change cannot be undone as a no-op.
    double stored = 0.0;
    for (const auto &en : e) stored += en.value * R[en.row] * R[en.col];
    double dotRR = 0.0;
    for (int i = 0; i < 3; ++i) dotRR += R[i] * R[i];
    const double oldValue = 2.0 * stored - dotRR;

    double excess = 0.0;
    for (int i = 0; i < 3; ++i) excess += (dense[i][i] - 1.0) * R[i] * R[i];
    CHECK_REL(oldValue - want, excess, 1e-13);
    CHECK(std::abs(oldValue - want) > 1.0);
}

// The identity GRM is the case the two spellings agree on, and it is the
// base (unrelated) path SPAGxE takes when no GRM is supplied.
TEST(quad_form_reduces_to_the_squared_norm_on_the_identity) {
    std::vector<SparseGRM::Entry> e = {{0, 0, 1.0}, {1, 1, 1.0}, {2, 2, 1.0}};
    const SparseGRM g = SparseGRM::fromEntries(3, e);
    const double R[3] = {1.7, -0.9, 1.3};
    const double want = R[0] * R[0] + R[1] * R[1] + R[2] * R[2];
    CHECK_REL(g.quadForm(R, 3), want, 1e-15);
}

// A non-finite entry must be named as such.  The diagonal case is the one
// that matters: buildDiagonal uses NaN as its "no entry" sentinel, so unless
// the range check runs FIRST a literal `nan` diagonal is misreported as a
// missing diagonal entry.
TEST(sparse_grm_rejects_non_finite_entries) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();

    struct Case { const char *name; std::vector<SparseGRM::Entry> e; };
    const Case cases[] = {
        {"NaN on the diagonal",     {{0, 0, nan}, {1, 1, 1.0}}},
        {"NaN off the diagonal",    {{0, 0, 1.0}, {1, 0, nan}, {1, 1, 1.0}}},
        {"infinity on the diagonal",{{0, 0, inf}, {1, 1, 1.0}}},
    };
    for (const Case &c : cases) {
        bool threw = false;
        std::string what;
        try {
            SparseGRM::fromEntries(2, c.e);
        } catch (const std::exception &ex) {
            threw = true;
            what = ex.what();
        }
        CHECK_MSG(threw, c.name);
        // Named as non-finite, not as a missing diagonal entry.
        CHECK_MSG(what.find("non-finite") != std::string::npos, c.name);
    }
}

TINYTEST_MAIN
