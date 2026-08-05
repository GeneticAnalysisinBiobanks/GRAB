#!/usr/bin/env python3
"""make_synthetic.py — Generate a synthetic cohort for end-to-end benchmarking.

The bundled examples/1kg.* fixture holds 2504 subjects, which is adequate for
correctness regression but far too small to measure the saddlepoint rework's
effect: at that size the per-marker cost is dominated by I/O and engine
overhead rather than by the CGF kernel and root finder under change.  This
script writes a PLINK1 BED/BIM/FAM triple, a phenotype table matching the
schema of examples/1kg.pheno, and a sparse GRM in the format GRAB's
--sparse-grm reader expects, at whatever scale the benchmark needs.

The data are pure noise with no genotype-phenotype association by
construction.  That is deliberate: under the null the p-value distribution is
uniform, so the same output doubles as a calibration check (lambda_GC should
sit near 1), and the fraction of markers taking the saddlepoint branch is
governed only by --spa-z-threshold rather than by a planted signal.

Usage:
    tests/make_synthetic.py --n 50000 --m 20000 --out bench_data/syn
    build/grab2 --method spacox --pfile ... (see tests/bench_e2e.sh)

Reproducible: the PRNG is seeded from --seed (default 20260730), so two runs
with identical arguments produce byte-identical files.

--------------------------------------------------------------------------
WtCoxG / LEAF inputs
--------------------------------------------------------------------------

WtCoxG and LEAF cannot run on a genotype/phenotype/GRM quartet alone: WtCoxG
needs an external reference allele frequency (--ref-af) and LEAF needs one
reference file per ancestral population plus a cluster assignment.  Both are
written here so that the two methods the saddlepoint fails on most often are
also the two the synthetic null cohort can measure.

  <out>.ref_pop<k>.afreq   plink2 --freq layout, one row per .bim record, in
                           .bim order.  Header
                               #CHROM  ID  REF  ALT  ALT_FREQS  OBS_CT
                           matched to the .bim by (CHROM, ID);  REF/ALT are
                           the .bim col6/col5 alleles, so no strand flip is
                           exercised and ALT_FREQS is the frequency of the
                           allele GRAB's cursor counts.
  <out>.clusters.txt       --leaf-cluster-file layout: header `#IID cluster`,
                           cluster in {1, ..., --n-clusters}.

Null by construction, for the batch-effect test as well.  WtCoxG's batch test
is

    z_bat = (weighted study MAF - AF_ref) / sqrt(v),
    v     = (1/(2N) + 1/OBS_CT) * p (1 - p)      (uniform weights),

so z_bat is standard normal exactly when AF_ref is an INDEPENDENT sample of
OBS_CT alleles from the same allele frequency the genotypes were drawn from.
Each panel is therefore Binomial(OBS_CT, p_j) / OBS_CT with the same per-
marker p_j the .bed was generated from, drawn from a dedicated PRNG stream so
that adding or removing panels cannot perturb the .bed.  Two consequences are
worth stating because they are properties of the fixture, not of the method:

  * The panels are mutually independent, which is what LEAF's summix step
    assumes when it forms  AF_ref = sum_k pi_k AF_k  and reports an effective
    allele count  1 / sum_k (pi_k^2 / OBS_CT_k).  With panels drawn from one
    common p_j the individual pi_k are not identified -- any convex
    combination estimates the same p_j -- but the constrained least squares
    converges on the inverse-variance weights, and both AF_ref and the
    effective allele count are then exactly the quantities the batch test's
    null distribution requires.  Give the panels DIFFERENT OBS_CT (the
    default does) so this is visible rather than degenerate.
  * --prevalence must be set to the cohort's own case fraction.  WtCoxG
    weights controls by (1-K)/K / ((1-er)/er) with er the sample case
    fraction; a randomly ascertained cohort like this one has er = K, the
    weights collapse to 1, and the variance formula above is exact.  Passing
    a prevalence that disagrees with the realised case fraction re-weights a
    sample that was never ascertained on the outcome and miscalibrates the
    batch test.  The realised case fraction is printed at the end of the run.

--ref-batch-shift introduces a controlled, non-null batch effect: the panel
frequency is drawn from expit(logit(p_j) + delta_j) with delta_j ~ N(0, s^2).
It is 0 by default and must stay 0 for any cohort used as a calibration gate.

--------------------------------------------------------------------------
Planted signal (NOT null -- witness fixtures only)
--------------------------------------------------------------------------

--signal-markers / --signal-beta plant a genotype effect on the Time and
Binary traits.  A cohort generated with either of them is NOT a null cohort
and must never be used for a calibration gate.  Their purpose is to reach
statistics large enough to underflow the tail assemblies the log10p_unify
project is replacing, which no null cohort of any size can do: under the null
the largest |Z| over 1e6 markers is about 5, and the +Inf that
wtcoxg_cond::conditionalP can emit needs |Z| in the tens.
"""

import argparse
import os
import sys

import numpy as np


def write_fam(path, iids):
    with open(path, "w") as fh:
        for iid in iids:
            # FID IID PAT MAT SEX PHENO — GRAB keys on IID; the rest are padding.
            fh.write("0\t%s\t0\t0\t0\t-9\n" % iid)


def write_bim(path, m, rng):
    """22 autosomes, positions strictly increasing within each chromosome.

    Returns the (CHROM, ID) pairs in .bim order.  The reference-AF writer
    needs them because loadRefAfFile/matchMarkers key on (CHROM, ID); nothing
    about the returned value changes what is written or what the PRNG
    consumes.
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
                fh.write("%d\tsyn%d_%d\t0\t%d\tA\tG\n"
                         % (chrom, chrom, i, pos[i]))
                keys.append((str(chrom), "syn%d_%d" % (chrom, i)))
            written += k
    return keys


def write_bed(path, n, m, mafs, miss_rate, rng, chunk=512, keep_rows=()):
    """Write a variant-major PLINK1 .bed.

    PLINK1 two-bit codes, four samples per byte, least significant pair first:
        00  homozygous for the .bim A1 allele   (ALT dosage 2)
        01  missing
        10  heterozygous                        (ALT dosage 1)
        11  homozygous for the .bim A2 allele   (ALT dosage 0)

    `keep_rows` is a set of marker indices whose ALT dosages are returned to
    the caller (as {index: (n,) uint8, 3 = missing}) so that a planted-signal
    phenotype can be built from genotypes that were already written.  It does
    not alter what is drawn or written: a cohort generated with a planted
    signal has a .bed byte-identical to the same cohort generated without one,
    and only its phenotype table differs.
    """
    nbytes = (n + 3) // 4
    # Map ALT dosage 0,1,2 and missing(3) to the PLINK code.
    code_of = np.array([0b11, 0b10, 0b00, 0b01], dtype=np.uint8)
    keep_rows = set(int(i) for i in keep_rows)
    kept = {}

    with open(path, "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))       # magic + variant-major
        for start in range(0, m, chunk):
            k = min(chunk, m - start)
            # Genotypes for k variants at once: (k, n) ALT dosages.
            p = mafs[start:start + k][:, None]
            g = (rng.random((k, n)) < p).astype(np.uint8) + \
                (rng.random((k, n)) < p).astype(np.uint8)
            if miss_rate > 0:
                g = np.where(rng.random((k, n)) < miss_rate, 3, g)

            for j in keep_rows:
                if start <= j < start + k:
                    kept[j] = g[j - start].copy()

            codes = code_of[g]                     # (k, n) uint8 in 0..3

            # Pack four samples per byte, low pair first.
            pad = nbytes * 4 - n
            if pad:
                # Padding bits must be zero; code 00 is "hom A1", but PLINK
                # ignores bits past the sample count, so any value is legal.
                codes = np.concatenate(
                    [codes, np.zeros((k, pad), dtype=np.uint8)], axis=1)
            c = codes.reshape(k, nbytes, 4)
            packed = (c[:, :, 0] |
                      (c[:, :, 1] << 2) |
                      (c[:, :, 2] << 4) |
                      (c[:, :, 3] << 6)).astype(np.uint8)
            fh.write(packed.tobytes())

    return kept


def write_pheno(path, iids, rng, signal=None):
    """Phenotype table with the same columns examples/1kg.pheno provides.

    `signal`, when given, is a length-n additive genetic contribution on the
    linear-predictor scale.  It enters the Time and Binary traits only, and
    leaves the Quantitative trait null so that one calibrated trait survives
    on a witness fixture.  It consumes no randomness of its own, so the PRNG
    stream -- and therefore the .bed, .fam, .bim and GRM -- is identical with
    and without it.
    """
    n = len(iids)
    male = rng.integers(0, 2, size=n)
    pcs = rng.normal(0.0, 0.02, size=(n, 4))

    # No genotype effect anywhere: every trait is a function of covariates and
    # noise only, so the null holds exactly and lambda_GC is interpretable.
    lin = 0.30 * male + 3.0 * pcs[:, 0] - 2.0 * pcs[:, 1]
    quant = lin + rng.normal(0.0, 1.0, size=n)

    lin_g = lin if signal is None else lin + signal

    prob = 1.0 / (1.0 + np.exp(-(lin_g - 1.2)))
    binary = (rng.random(n) < prob).astype(int)

    # Exponential survival with administrative censoring at t = 5.
    rate = np.exp(0.5 * lin_g - 1.0)
    t_event = rng.exponential(1.0 / rate)
    t_cens = rng.uniform(0.0, 5.0, size=n)
    time = np.minimum(t_event, t_cens)
    event = (t_event <= t_cens).astype(int)

    with open(path, "w") as fh:
        fh.write("#IID\tMALE\tPC1\tPC2\tPC3\tPC4\tQuantitative\tTime\tEvent\tBinary\n")
        for i, iid in enumerate(iids):
            fh.write("%s\t%d\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%.8g\t%d\t%d\n"
                     % (iid, male[i], pcs[i, 0], pcs[i, 1], pcs[i, 2], pcs[i, 3],
                        quant[i], time[i], event[i], binary[i]))

    # PC1 is what the cluster file bins on and what --pc-cols hands LEAF's
    # K-means; return it so the two agree by construction.
    return {"case_frac": float(np.mean(binary)),
            "event_frac": float(np.mean(event)),
            "pc1": pcs[:, 0]}


def write_sparse_grm(prefix, iids, rng, rel_frac=0.02):
    """Sparse GRM: unit-ish diagonal plus a small number of related pairs.

    Format mirrors examples/1kg.grm.{sp,id}: the .sp file holds
    'i<TAB>j<TAB>value' with zero-based indices into the .id file, which in
    turn holds 'FID<TAB>IID' lines.  Only the lower triangle plus the diagonal
    is stored, matching the bundled fixture.
    """
    n = len(iids)
    with open(prefix + ".grm.id", "w") as fh:
        for iid in iids:
            fh.write("0\t%s\n" % iid)

    n_pairs = int(n * rel_frac / 2)
    # Sibling-like pairs among adjacent indices keep the matrix well-banded,
    # which is what a real pedigree-derived sparse GRM looks like after the
    # subjects have been ordered by family.
    left = rng.choice(np.arange(0, n - 1, 2), size=min(n_pairs, (n - 1) // 2),
                      replace=False)
    with open(prefix + ".grm.sp", "w") as fh:
        for i in range(n):
            fh.write("%d\t%d\t%.8g\n" % (i, i, 1.0 + rng.normal(0.0, 0.01)))
        for i in left:
            fh.write("%d\t%d\t%.8g\n" % (i + 1, i, 0.5 + rng.normal(0.0, 0.02)))


def write_ref_afreq(path, keys, mafs, obs_ct, rng, batch_shift=0.0):
    """One plink2 --freq (.afreq) file: an independent reference panel.

    `keys` are the (CHROM, ID) pairs in .bim order and `mafs` the per-marker
    ALT frequency the .bed was generated from.  The panel frequency is

        AF_ref_j = Binomial(obs_ct, q_j) / obs_ct,
        q_j      = p_j                                  (batch_shift = 0)
        q_j      = expit(logit(p_j) + N(0, batch_shift^2))   otherwise,

    so at batch_shift = 0 the panel is an unbiased sample of obs_ct alleles
    from the same allele frequency as the study genotypes, with variance
    p(1-p)/obs_ct -- exactly the `1/OBS_CT` term of WtCoxG's batch-test
    variance.  REF/ALT are written as the .bim col6/col5 alleles (G/A), the
    orientation matchMarkers() treats as unflipped.

    OBS_CT is an ALLELE count, as in plink2 --freq: a panel of `obs_ct/2`
    diploid individuals.
    """
    q = np.asarray(mafs, dtype=float)
    if batch_shift > 0.0:
        logit = np.log(q / (1.0 - q)) + rng.normal(0.0, batch_shift, size=q.size)
        q = 1.0 / (1.0 + np.exp(-logit))
    af = rng.binomial(int(obs_ct), q) / float(obs_ct)

    with open(path, "w") as fh:
        fh.write("#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n")
        for (chrom, mid), a in zip(keys, af):
            fh.write("%s\t%s\tG\tA\t%.10g\t%d\n" % (chrom, mid, a, obs_ct))


def write_cluster_file(path, iids, pc1, n_clusters):
    """A --leaf-cluster-file: equal-sized PC1 quantile bins.

    LEAF's own path runs K-means on --pc-cols, so binning PC1 puts the file
    and the K-means it replaces on the same footing.  Because PC1 is pure
    noise in this cohort the bins are a random partition of the subjects, and
    every cluster therefore has the same allele frequencies and the same null.
    Quantile bins rather than a fresh random draw so the assignment is a
    deterministic function of the phenotype file: it can be recomputed by any
    reader without knowing the PRNG state.
    """
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
    ap.add_argument("--n", type=int, default=50000, help="subjects")
    ap.add_argument("--m", type=int, default=20000, help="markers")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--miss-rate", type=float, default=0.002)
    ap.add_argument("--min-maf", type=float, default=0.005)
    ap.add_argument("--no-grm", action="store_true")
    ap.add_argument("--ref-panels", type=int, default=2,
                    help="number of <out>.ref_pop<k>.afreq files (0 = none)")
    ap.add_argument("--ref-obs-ct", default="20000,16000",
                    help="comma-separated ALLELE count per reference panel; "
                         "the last value is reused if fewer than --ref-panels")
    ap.add_argument("--ref-batch-shift", type=float, default=0.0,
                    help="SD of a per-marker logit-scale batch offset applied "
                         "to the reference panels; 0 = null (default)")
    ap.add_argument("--n-clusters", type=int, default=3,
                    help="clusters in <out>.clusters.txt (0 = no file)")
    ap.add_argument("--signal-markers", type=int, default=0,
                    help="plant a genotype effect on this many markers; "
                         "NON-NULL, witness fixtures only")
    ap.add_argument("--signal-beta", type=float, default=0.0,
                    help="per-marker effect on the standardized genotype")
    args = ap.parse_args()

    if args.signal_markers > 0 and args.signal_beta == 0.0:
        sys.exit("--signal-markers requires a non-zero --signal-beta")
    if args.signal_markers > args.m:
        sys.exit("--signal-markers exceeds --m")

    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    iids = ["SYN%07d" % i for i in range(args.n)]

    # Log-uniform MAF spectrum between --min-maf and 0.5: this puts most
    # markers at low frequency, as in a real cohort, which matters because the
    # saddlepoint branch fires disproportionately at low MAF.
    lo, hi = np.log(args.min_maf), np.log(0.5)
    mafs = np.exp(rng.uniform(lo, hi, size=args.m))

    # Planted markers are spread evenly through the .bim so that every
    # chromosome and every chunk of the marker engine carries some.
    sig_idx = (np.linspace(0, args.m - 1, args.signal_markers).astype(int)
               if args.signal_markers > 0 else np.empty(0, dtype=int))

    sys.stderr.write("writing %s.{bed,bim,fam,pheno} — %d subjects x %d markers\n"
                     % (args.out, args.n, args.m))
    write_fam(args.out + ".fam", iids)
    keys = write_bim(args.out + ".bim", args.m, rng)
    kept = write_bed(args.out + ".bed", args.n, args.m, mafs, args.miss_rate,
                     rng, keep_rows=sig_idx)

    signal = None
    if args.signal_markers > 0:
        signal = np.zeros(args.n)
        for j in sig_idx:
            g = kept[int(j)].astype(float)
            g[g == 3] = 2.0 * mafs[j]              # mean-impute the missing
            sd = np.sqrt(2.0 * mafs[j] * (1.0 - mafs[j]))
            signal += args.signal_beta * (g - 2.0 * mafs[j]) / sd

    stats = write_pheno(args.out + ".pheno", iids, rng, signal=signal)
    if not args.no_grm:
        write_sparse_grm(args.out, iids, rng)

    # Everything below draws from PRNG streams of its own, keyed off --seed,
    # so the four files above are byte-identical whatever is asked for here.
    if args.ref_panels > 0:
        obs = [int(float(x)) for x in args.ref_obs_ct.split(",") if x.strip()]
        if not obs:
            sys.exit("--ref-obs-ct is empty")
        while len(obs) < args.ref_panels:
            obs.append(obs[-1])
        for k in range(args.ref_panels):
            ref_rng = np.random.default_rng([args.seed, 101, k])
            path = "%s.ref_pop%d.afreq" % (args.out, k + 1)
            write_ref_afreq(path, keys, mafs, obs[k], ref_rng,
                            batch_shift=args.ref_batch_shift)
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
    if args.signal_markers > 0 or args.ref_batch_shift > 0.0:
        sys.stderr.write(
            "WARNING: this cohort is NOT null "
            "(--signal-markers %d, --ref-batch-shift %g); "
            "do not use it for a calibration gate\n"
            % (args.signal_markers, args.ref_batch_shift))


if __name__ == "__main__":
    main()
