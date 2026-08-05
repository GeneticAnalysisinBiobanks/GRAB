#!/usr/bin/env python3
"""make_fallback.py — a NULL cohort on which the saddlepoint actually fails.

Why this script exists
----------------------
Decision D5 of dev-notes/methods/log10p_unify/00_decisions.md makes the
saddlepoint's failure modes report the two-sided normal tail instead of NA,
under SPA_STATUS codes 3..6, and requires --help and CLAUDE.md to carry a
QUANTIFIED warning about that substitution.  Quantifying it needs rows with
SPA_STATUS in [3,6], and tests/make_synthetic.py produces exactly zero of
them: on its 50 000 x 20 000 cohort all 260 000 (marker, column) cells come
back 0 (SPA_OK) or 1 (NORMAL).  The only fallbacks the repository could
observe were 14 markers of the bundled examples/1kg fixture -- far too few for
a rate, a KS test, or a MAC binning.

This generator supplies the missing input.  It is a separate script rather
than another switch on make_synthetic.py because its genotype construction is
different in kind (a simulated pedigree rather than independent draws) and
because make_synthetic.py's PRNG stream is pinned by the byte-reproducibility
of bench_data/syn.*, bench_data/null.* and bench_data/witness.*.

What actually makes the saddlepoint fail (measured, not assumed)
---------------------------------------------------------------
A sweep over SPACox, SPAmix, SPAsqr, SPAGxE, SPAGxEmix, WtCoxG, LEAF and
SPAGRM, at --spa-z-threshold values down to 0.05 so that essentially every
marker takes the saddlepoint branch, found that on independently-drawn
genotypes NONE of the eight adapters ever fails -- 0 cells in [3,6] out of
tens of thousands.  Exactly one configuration in the whole tree fails, and it
is the one the bundled fixture already exhibits:

    SPAGRM, on a residual vector with many IQR outliers (a logistic or Cox
    null model), when those outliers sit in GRM families of three or more
    members, so that SPAGRM's CLASS-3 tabulated family CGF
    (src/spagrm/spagrm_cgf.hpp, the `c.three` loop) is populated.

Two controls pin it down.  Replacing the 1kg GRM by its own diagonal -- which
turns every outlier into a class-1 "unrelated outlier" and builds no family
table -- takes the failure count from 846/997 to 0/997 with the residuals
unchanged.  Conversely, keeping the 1kg phenotype, GRM and IBD but swapping in
an independently-drawn synthetic marker panel keeps the failure rate at
492/600, so the failure is a property of the NULL MODEL, not of the markers.

The failure is `zeta*s - K(zeta) < 0` (FALLBACK_GUARD_TEMP) and it is
confined to SMALL |Z|: on the bundled fixture at n = 2504 every failing row
has |Z_Norm| <= 1.52 and every succeeding row has |Z_Norm| >= 1.29.  That
signature is what a CGF with K(0) = delta > 0 produces, since the Legendre
transform then satisfies  sup_t (t*s - K(t)) >= -K(0),  i.e. it is negative
for every s within about sqrt(2*delta) standard deviations of the mean.  This
script does not diagnose the class-3 table further; it reproduces the regime
so that the fallback's behaviour can be measured.

Three further controls narrow it to the GRM's TOPOLOGY:

  * A simulated pedigree -- founders plus Mendelian offspring, families of
    four, 60 % of subjects related, 158 class-3 family models built -- does
    NOT fail, at any --spa-z-threshold.  Correctly-specified relatedness is
    therefore not sufficient.
  * A synthetic cohort wearing the BUNDLED 1kg GRM's topology (one dense
    component of 1165 of the 2504 subjects, 39 799 off-diagonal entries,
    values 0.2 .. 1.659) fails on 2 370 / 3 000 markers.  Since the cohort is
    otherwise entirely synthetic, the 1kg genotypes play no part.
  * The same cohort with that GRM stripped to its diagonal does not fail.

So the failing configuration is a large, densely connected declared-relatedness
component -- relatedness the genotypes do not exhibit, which is exactly what
the bundled fixture declares -- combined with an outlier-rich residual vector.
The generator reproduces that configuration directly (--grm dense) and offers
the correctly-specified pedigree (--grm pedigree) as the control that does not
fail.

A cohort built with --grm dense is null for the ASSOCIATION test being
measured -- the markers are independent of the phenotype -- while its declared
relatedness is deliberately not the relatedness the genotypes carry.  That
mis-specification is a property of the variance model, not of the null, and it
is the property under study: it is the only configuration in this repository
in which the D5 fallback is reachable.

Null by construction
--------------------
Founder haplotypes and Mendelian transmissions are drawn from a PRNG stream
that never sees the phenotype table, and the phenotype table is a function of
covariates and noise alone.  Genotype and phenotype are therefore independent,
exactly as in make_synthetic.py.  Relatives share genotypes but not phenotype
effects; that is a correctly-specified null for a score test, and it is the
configuration SPAGRM's family CGF exists to model.

Output files (layouts identical to make_synthetic.py, so the same command
lines work):

    <out>.bed/.bim/.fam        PLINK1 triple, variant-major, no missingness
    <out>.pheno                #IID MALE PC1..PC4 Quantitative Skew Time Event
                               Binary
    <out>.grm.sp/.grm.id       sparse GRM: pedigree kinship, 0.5 within each
                               nuclear family
    <out>.ref_pop<k>.afreq     plink2 --freq layout for WtCoxG / LEAF
    <out>.clusters.txt         --leaf-cluster-file layout
    <out>.design.txt           per-marker designed founder frequency and
                               realised MAC, so an analysis script can bin on
                               the design without parsing an output table

Reproducible: every stream is keyed off --seed, so identical arguments give
byte-identical files.

Usage (the invocation the Stage 0c ledger records):

    tests/make_fallback.py --n 20000 --m 40000 --seed 20260801 \
        --case-frac 0.05 --event-frac 0.05 --family-size 4 --rel-frac 0.6 \
        --min-maf 2e-4 --out bench_data/fallback
"""

import argparse
import os
import sys

import numpy as np


def write_fam(path, iids):
    with open(path, "w") as fh:
        for iid in iids:
            fh.write("0\t%s\t0\t0\t0\t-9\n" % iid)


def write_bim(path, m, rng):
    """22 autosomes, positions strictly increasing within each chromosome.

    Returns the (CHROM, ID) pairs in .bim order; loadRefAfFile/matchMarkers
    (src/wtcoxg/wtcoxg.cpp) key the reference panels on them.
    """
    per = max(1, m // 22)
    keys = []
    with open(path, "w") as fh:
        written = 0
        for chrom in range(1, 23):
            k = per if chrom < 22 else m - written
            if k <= 0:
                break
            pos = np.sort(rng.choice(np.arange(10_000, 250_000_000), size=k,
                                     replace=False))
            for i in range(k):
                fh.write("%d\tfb%d_%d\t0\t%d\tA\tG\n" % (chrom, chrom, i, pos[i]))
                keys.append((str(chrom), "fb%d_%d" % (chrom, i)))
            written += k
    return keys


def build_pedigree(n, rel_frac, fam_size):
    """Assign subjects to nuclear families; return (father, mother) index arrays.

    The first `n_fam * fam_size` subjects occupy consecutive family blocks of
    `fam_size`: two founders followed by `fam_size - 2` full siblings.  The
    remainder are unrelated founders.  Consecutive blocks keep the sparse GRM
    banded, which is what a real pedigree-derived GRM looks like once subjects
    are family-ordered.

    A founder is marked by parent index -1.
    """
    fam_size = max(3, int(fam_size))
    n_in_fam = int(n * rel_frac)
    n_fam = n_in_fam // fam_size
    father = np.full(n, -1, dtype=np.int64)
    mother = np.full(n, -1, dtype=np.int64)
    for f in range(n_fam):
        st = f * fam_size
        for c in range(2, fam_size):
            father[st + c] = st
            mother[st + c] = st + 1
    return father, mother, n_fam, fam_size


def write_bed(path, n, freqs, father, mother, rng, chunk=256):
    """Variant-major PLINK1 .bed drawn from the pedigree.

    PLINK1 two-bit codes, four samples per byte, least significant pair first:
        00  homozygous ALT (dosage 2)   01  missing
        10  heterozygous (dosage 1)     11  homozygous REF (dosage 0)

    Founders receive two independent Bernoulli(p_j) haplotypes; each offspring
    receives one uniformly chosen haplotype from each parent.  Offspring are
    always at a higher index than their parents, so a single forward pass
    suffices.

    No missing calls are written: missingness would blur the MAC binning the
    G4 report rests on, and make_synthetic.py already covers that path.

    Returns the realised per-marker minor allele count.
    """
    m = freqs.size
    nbytes = (n + 3) // 4
    code_of = np.array([0b11, 0b10, 0b00], dtype=np.uint8)
    kids = np.nonzero(father >= 0)[0]
    kf, km = father[kids], mother[kids]
    macs = np.empty(m, dtype=np.int64)

    with open(path, "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        for start in range(0, m, chunk):
            k = min(chunk, m - start)
            p = freqs[start:start + k][:, None]
            # Two haplotypes per subject; offspring haplotypes are overwritten
            # by a transmitted copy below, which is why they are drawn too
            # (keeping the draw shape constant keeps the stream reproducible
            # independently of the pedigree).
            h1 = (rng.random((k, n)) < p).astype(np.uint8)
            h2 = (rng.random((k, n)) < p).astype(np.uint8)
            if kids.size:
                pick_f = rng.integers(0, 2, size=(k, kids.size))
                pick_m = rng.integers(0, 2, size=(k, kids.size))
                fa = np.where(pick_f == 0, h1[:, kf], h2[:, kf])
                mo = np.where(pick_m == 0, h1[:, km], h2[:, km])
                h1[:, kids] = fa
                h2[:, kids] = mo
            g = (h1 + h2).astype(np.uint8)

            cnt = g.sum(axis=1, dtype=np.int64)
            macs[start:start + k] = np.minimum(cnt, 2 * n - cnt)

            codes = code_of[g]
            pad = nbytes * 4 - n
            if pad:
                codes = np.concatenate(
                    [codes, np.zeros((k, pad), dtype=np.uint8)], axis=1)
            c = codes.reshape(k, nbytes, 4)
            packed = (c[:, :, 0] | (c[:, :, 1] << 2) |
                      (c[:, :, 2] << 4) | (c[:, :, 3] << 6)).astype(np.uint8)
            fh.write(packed.tobytes())

    return macs


def _solve_intercept(lin, target, link):
    """Bisection for the intercept b with E[link(lin + b)] == target.

    `link` maps the linear predictor to a probability and is monotone in b, so
    bisection over a wide bracket converges unconditionally; 200 iterations
    take it to the last bits of a double, making the result a deterministic
    function of `lin` and `target` alone.
    """
    lo, hi = -60.0, 60.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if link(lin + mid).mean() > target:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def write_pheno(path, iids, rng, case_frac, event_frac, skew_df, frac_grid):
    """Covariates plus four traits, none of them a function of any genotype.

    Quantitative  Gaussian residuals -- the calibrated control trait, whose
                  null model produces no IQR outliers and therefore no family
                  models, hence no saddlepoint failures.  It is present so
                  that every report has a within-cohort control column.
    Skew          heavy-tailed, right-skewed residuals (exponentiated
                  Student-t with --skew-df degrees of freedom).
    Time/Event    Cox, event fraction driven to --event-frac.
    Binary        logistic, case fraction driven to --case-frac.

    B<ppp>        one extra logistic trait per entry of --case-frac-grid, at
                  case fraction ppp/1000.  Whether SPAGRM's class-3 CGF fails,
                  and up to which |Z|, turns out to depend erratically on the
                  particular outlier-family configuration a null model
                  produces -- at n = 2504 with the 1kg GRM topology, case
                  fractions 0.012, 0.108 and 0.249 fail on ~87 % of markers
                  while 0.048 and 0.498 do not fail at all.  A grid of traits
                  in one file therefore samples that space in a single run and
                  supplies both failing columns and within-cohort controls.

    The thresholded traits are what produce a residual vector with many IQR
    outliers, which is the first of the two conditions SPAGRM's class-3 family
    CGF needs.
    """
    n = len(iids)
    male = rng.integers(0, 2, size=n)
    pcs = rng.normal(0.0, 0.02, size=(n, 4))

    lin = 0.30 * male + 3.0 * pcs[:, 0] - 2.0 * pcs[:, 1]
    quant = lin + rng.normal(0.0, 1.0, size=n)

    t_noise = rng.standard_t(skew_df, size=n)
    skew = lin + np.exp(np.clip(t_noise, -20.0, 6.0))

    def expit(x):
        return 1.0 / (1.0 + np.exp(-x))

    b0 = _solve_intercept(lin, case_frac, expit)
    binary = (rng.random(n) < expit(lin + b0)).astype(int)

    t_cens = rng.uniform(0.0, 5.0, size=n)
    haz_lin = 0.5 * lin

    def ev_link(x):
        return 1.0 - np.exp(-np.exp(x) * t_cens)

    c0 = _solve_intercept(haz_lin, event_frac, ev_link)
    t_event = rng.exponential(1.0 / np.exp(haz_lin + c0))
    time = np.minimum(t_event, t_cens)
    event = (t_event <= t_cens).astype(int)

    grid = {}
    for f in frac_grid:
        name = "B%03d" % int(round(f * 1000))
        grid[name] = (rng.random(n) < expit(lin + _solve_intercept(lin, f, expit))
                      ).astype(int)
    gnames = sorted(grid)

    with open(path, "w") as fh:
        fh.write("#IID\tMALE\tPC1\tPC2\tPC3\tPC4\tQuantitative\tSkew"
                 "\tTime\tEvent\tBinary"
                 + "".join("\t" + g for g in gnames) + "\n")
        for i, iid in enumerate(iids):
            fh.write("%s\t%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%d\t%d%s\n"
                     % (iid, male[i], pcs[i, 0], pcs[i, 1], pcs[i, 2],
                        pcs[i, 3], quant[i], skew[i], time[i], event[i],
                        binary[i],
                        "".join("\t%d" % grid[g][i] for g in gnames)))

    return {"case_frac": float(np.mean(binary)),
            "event_frac": float(np.mean(event)),
            "grid": {g: float(grid[g].mean()) for g in gnames},
            "pc1": pcs[:, 0]}


def write_sparse_grm(prefix, iids, father, mother, fam_size, rng):
    """Pedigree kinship as a sparse GRM, in examples/1kg.grm.{sp,id} layout.

    'i<TAB>j<TAB>value' over zero-based indices into the .id file, lower
    triangle plus diagonal.  Within a nuclear family every parent-offspring
    and every sibling pair gets 0.5 plus a small jitter, which is the relatedness
    the genotypes were actually simulated with; the two founders of a family
    are not related to each other and get no entry.

    Family blocks of three or more are the point of this writer:
    make_synthetic.py emits sibling-like PAIRS only, so SPAGRM's class-3
    tabulated family CGF is never populated and its saddlepoint never fails.
    """
    n = len(iids)
    with open(prefix + ".grm.id", "w") as fh:
        for iid in iids:
            fh.write("0\t%s\n" % iid)

    n_pairs = 0
    with open(prefix + ".grm.sp", "w") as fh:
        for i in range(n):
            fh.write("%d\t%d\t%.8g\n" % (i, i, 1.0 + rng.normal(0.0, 0.01)))
        kids = np.nonzero(father >= 0)[0]
        seen = set()
        for c in kids:
            for par in (int(father[c]), int(mother[c])):
                key = (int(c), par)
                if key not in seen:
                    seen.add(key)
                    fh.write("%d\t%d\t%.8g\n"
                             % (c, par, 0.5 + rng.normal(0.0, 0.02)))
                    n_pairs += 1
        # Sibling pairs: same (father, mother).
        for a in range(kids.size):
            for b in range(a):
                ca, cb = int(kids[a]), int(kids[b])
                if father[ca] == father[cb] and mother[ca] == mother[cb]:
                    fh.write("%d\t%d\t%.8g\n"
                             % (max(ca, cb), min(ca, cb),
                                0.5 + rng.normal(0.0, 0.02)))
                    n_pairs += 1
    return n_pairs


def write_dense_grm(prefix, iids, dense_frac, degree, rng):
    """A 1kg-shaped sparse GRM: one large, densely connected component.

    The bundled examples/1kg.grm.sp declares 39 799 off-diagonal entries over
    a single component of 1165 of its 2504 subjects (mean degree 68, values
    0.2 to 1.659) while the 1kg genotypes themselves carry no such
    relatedness.  That combination -- a large dense DECLARED component over
    genotypes that do not exhibit it -- is the configuration in which SPAGRM's
    class-3 family CGF fails; see the module docstring for the three controls
    that isolate it.  This writer reproduces its shape:

        the first  dense_frac * n  subjects form one component,
        each is joined to `degree` uniformly chosen partners inside it,
        entry value  0.2 + Exponential(0.05), matching the bundled file's
        0.2 floor and long right tail.

    Returns the number of off-diagonal entries written.
    """
    n = len(iids)
    with open(prefix + ".grm.id", "w") as fh:
        for iid in iids:
            fh.write("0\t%s\n" % iid)

    nd = int(n * dense_frac)
    pairs = set()
    if nd >= 2 and degree > 0:
        for i in range(nd):
            partners = rng.choice(nd, size=min(degree, nd - 1), replace=False)
            for j in partners:
                j = int(j)
                if j == i:
                    continue
                pairs.add((max(i, j), min(i, j)))
    ordered = sorted(pairs)
    vals = 0.2 + rng.exponential(0.05, size=len(ordered))

    # The bundled fixture's DIAGONAL is the second thing that distinguishes it
    # from a textbook GRM: 0.460 .. 3.809 with median 0.744, i.e. mostly BELOW
    # one.  A diagonal below one shrinks R' Phi R relative to sum(R^2) and
    # therefore moves SPAGRM's EmpVar / Score_var ratio, which is the factor
    # its saddlepoint target is rescaled by (spagrm.cpp: Score_adj).  A
    # diagonal pinned at one does not reproduce the failure; this one does.
    diag = 0.45 + rng.exponential(0.30, size=n)
    with open(prefix + ".grm.sp", "w") as fh:
        for i in range(n):
            fh.write("%d\t%d\t%.8g\n" % (i, i, diag[i]))
        for (a, b), v in zip(ordered, vals):
            fh.write("%d\t%d\t%.8g\n" % (a, b, v))
    return len(ordered), nd


def write_template_grm(prefix, iids, template):
    """Copy an existing sparse GRM's .sp verbatim; write a matching .id.

    The bundled examples/1kg.grm.sp is the only sparse GRM in the repository
    on which SPAGRM's saddlepoint is known to fail, and the parametric
    imitation in write_dense_grm does not reproduce it: matching the component
    size, the degree, the 0.2 off-diagonal floor and even the sub-unit
    diagonal is not enough.  Rather than keep guessing at which further moment
    of that file matters, this mode uses the file itself and leaves the
    remaining question -- exactly which feature of it drives K(0) away from
    zero -- open and stated.

    Only the .id is rewritten, so the cohort's own IIDs are used and nothing
    of the 1kg genotypes or phenotypes enters.  --n must equal the template's
    subject count.
    """
    n = len(iids)
    src = template
    idsrc = template[:-3] + ".id" if template.endswith(".sp") else template + ".id"
    if not os.path.exists(idsrc):
        sys.exit("--grm-template: expected a companion %s" % idsrc)
    with open(idsrc) as fh:
        ntpl = sum(1 for _ in fh)
    if ntpl != n:
        sys.exit("--grm-template has %d subjects, --n is %d" % (ntpl, n))
    with open(src) as fin, open(prefix + ".grm.sp", "w") as fout:
        nline = 0
        for line in fin:
            fout.write(line)
            nline += 1
    with open(prefix + ".grm.id", "w") as fh:
        for iid in iids:
            fh.write("0\t%s\n" % iid)
    return nline


def write_ref_afreq(path, keys, freqs, obs_ct, rng):
    """An independent reference panel: Binomial(obs_ct, p_j) / obs_ct.

    Identical construction and identical null argument to
    make_synthetic.py's write_ref_afreq at --ref-batch-shift 0.  p_j is the
    DESIGNED founder frequency; on a pedigree cohort the realised study
    frequency deviates from it by drift within families, which is a real
    property of any family cohort rather than a defect of the fixture, and it
    inflates WtCoxG's batch test slightly.  Any use of this cohort as a
    batch-test calibration gate must account for that; the fallback
    measurement this script exists for does not use the batch test.
    """
    af = rng.binomial(int(obs_ct), np.asarray(freqs, dtype=float)) / float(obs_ct)
    with open(path, "w") as fh:
        fh.write("#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n")
        for (chrom, mid), a in zip(keys, af):
            fh.write("%s\t%s\tG\tA\t%.10g\t%d\n" % (chrom, mid, a, obs_ct))


def write_cluster_file(path, iids, pc1, n_clusters):
    """--leaf-cluster-file: equal-sized PC1 quantile bins (see make_synthetic)."""
    order = np.argsort(pc1, kind="stable")
    labels = np.empty(len(iids), dtype=int)
    edges = np.linspace(0, len(iids), n_clusters + 1).astype(int)
    for c in range(n_clusters):
        labels[order[edges[c]:edges[c + 1]]] = c + 1
    with open(path, "w") as fh:
        fh.write("#IID\tcluster\n")
        for iid, lab in zip(iids, labels):
            fh.write("%s\t%d\n" % (iid, lab))
    return labels


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20000, help="subjects")
    ap.add_argument("--m", type=int, default=40000, help="markers")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--seed", type=int, default=20260801)
    ap.add_argument("--min-maf", type=float, default=2e-4,
                    help="lower end of the log-uniform founder MAF spectrum")
    ap.add_argument("--case-frac", type=float, default=0.05,
                    help="Binary trait case fraction")
    ap.add_argument("--event-frac", type=float, default=0.05,
                    help="Time/Event trait event fraction")
    ap.add_argument("--skew-df", type=float, default=3.0,
                    help="degrees of freedom of the Skew trait's t noise")
    ap.add_argument("--case-frac-grid", default="0.01,0.02,0.05,0.10,0.15,"
                                                "0.20,0.25,0.30,0.40",
                    help="extra logistic traits B<ppp> at these case "
                         "fractions; empty for none")
    ap.add_argument("--grm", choices=("template", "dense", "pedigree", "none"),
                    default="template",
                    help="template: copy --grm-template's .sp verbatim and "
                         "write a matching .id (the configuration in which "
                         "the saddlepoint demonstrably fails).  dense: a "
                         "parametric imitation of it.  pedigree: real "
                         "Mendelian families, the control that does not fail")
    ap.add_argument("--grm-template", default="examples/1kg.grm.sp",
                    help="--grm template: sparse GRM whose .sp is copied; "
                         "--n must equal its subject count")
    ap.add_argument("--dense-frac", type=float, default=0.465,
                    help="--grm dense: fraction of subjects in the component "
                         "(1165/2504 in the bundled fixture)")
    ap.add_argument("--dense-degree", type=int, default=64,
                    help="--grm dense: partners drawn per subject")
    ap.add_argument("--rel-frac", type=float, default=0.6,
                    help="--grm pedigree: fraction of subjects in families")
    ap.add_argument("--family-size", type=int, default=4,
                    help="--grm pedigree: members per family (2 founders + "
                         "the rest full siblings)")
    ap.add_argument("--ref-panels", type=int, default=2)
    ap.add_argument("--ref-obs-ct", default="20000,16000")
    ap.add_argument("--n-clusters", type=int, default=3)
    args = ap.parse_args()

    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)

    iids = ["FB%07d" % i for i in range(args.n)]

    # Four independent streams.  Keeping the genotype stream separate from the
    # phenotype stream is what makes "the genotype draw never sees the
    # phenotype" a property of the code rather than of statement order.
    bim_rng = np.random.default_rng([args.seed, 0])
    geno_rng = np.random.default_rng([args.seed, 1])
    pheno_rng = np.random.default_rng([args.seed, 2])
    grm_rng = np.random.default_rng([args.seed, 3])

    lo, hi = np.log(args.min_maf), np.log(0.5)
    freqs = np.exp(geno_rng.uniform(lo, hi, size=args.m))

    # --grm dense declares relatedness the genotypes do not carry, exactly as
    # the bundled fixture does, so the genotype draw stays independent there.
    if args.grm == "pedigree":
        father, mother, n_fam, fam_size = build_pedigree(
            args.n, args.rel_frac, args.family_size)
    else:
        father = np.full(args.n, -1, dtype=np.int64)
        mother = np.full(args.n, -1, dtype=np.int64)
        n_fam, fam_size = 0, 0

    grid = [float(x) for x in args.case_frac_grid.split(",") if x.strip()]

    sys.stderr.write(
        "writing %s.{bed,bim,fam,pheno} — %d subjects x %d markers, grm=%s\n"
        % (args.out, args.n, args.m, args.grm))

    write_fam(args.out + ".fam", iids)
    keys = write_bim(args.out + ".bim", args.m, bim_rng)
    macs = write_bed(args.out + ".bed", args.n, freqs, father, mother, geno_rng)
    stats = write_pheno(args.out + ".pheno", iids, pheno_rng,
                        args.case_frac, args.event_frac, args.skew_df, grid)
    if args.grm == "pedigree":
        n_pairs = write_sparse_grm(args.out, iids, father, mother, fam_size,
                                   grm_rng)
        sys.stderr.write("  %s.grm.sp — %d families of %d, %d related pairs\n"
                         % (args.out, n_fam, fam_size, n_pairs))
    elif args.grm == "template":
        nline = write_template_grm(args.out, iids, args.grm_template)
        sys.stderr.write("  %s.grm.sp — %d lines copied from %s\n"
                         % (args.out, nline, args.grm_template))
    elif args.grm == "dense":
        n_pairs, nd = write_dense_grm(args.out, iids, args.dense_frac,
                                      args.dense_degree, grm_rng)
        sys.stderr.write("  %s.grm.sp — dense component of %d, %d off-diagonals\n"
                         % (args.out, nd, n_pairs))

    with open(args.out + ".design.txt", "w") as fh:
        fh.write("#ID\tAF_FOUNDER\tMAC_REALISED\n")
        for (_, mid), f, mac in zip(keys, freqs, macs):
            fh.write("%s\t%.10g\t%d\n" % (mid, f, mac))

    if args.ref_panels > 0:
        obs = [int(float(x)) for x in args.ref_obs_ct.split(",") if x.strip()]
        while len(obs) < args.ref_panels:
            obs.append(obs[-1])
        for k in range(args.ref_panels):
            ref_rng = np.random.default_rng([args.seed, 101, k])
            path = "%s.ref_pop%d.afreq" % (args.out, k + 1)
            write_ref_afreq(path, keys, freqs, obs[k], ref_rng)
            sys.stderr.write("  %s — OBS_CT %d\n" % (path, obs[k]))

    if args.n_clusters > 0:
        path = args.out + ".clusters.txt"
        labels = write_cluster_file(path, iids, stats["pc1"], args.n_clusters)
        sizes = ",".join(str(int((labels == c + 1).sum()))
                         for c in range(args.n_clusters))
        sys.stderr.write("  %s — %d clusters, sizes %s\n"
                         % (path, args.n_clusters, sizes))

    total = sum(os.path.getsize(args.out + ext)
                for ext in (".bed", ".bim", ".fam", ".pheno"))
    sys.stderr.write("done — %.1f MiB\n" % (total / (1 << 20)))
    sys.stderr.write("case fraction (Binary) %.6f — pass it as --prevalence\n"
                     % stats["case_frac"])
    sys.stderr.write("event fraction (Event) %.6f\n" % stats["event_frac"])
    if stats["grid"]:
        sys.stderr.write("grid case fractions: %s\n"
                         % ", ".join("%s=%.5f" % (k, v)
                                     for k, v in sorted(stats["grid"].items())))
    sys.stderr.write("realised MAC: min %d  median %d  max %d\n"
                     % (macs.min(), int(np.median(macs)), macs.max()))
    sys.stderr.write(
        "NULL by construction: genotypes are drawn from a PRNG stream that "
        "never sees the phenotype table.\n")


if __name__ == "__main__":
    main()
