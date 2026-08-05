#!/usr/bin/env python3
"""Build a pure-null cohort for SPAmixLocalPlus, from the fixtures in examples/.

`tests/make_synthetic.py` serves the six SPA methods that read a PLINK/PGEN
genotype file.  SPAmixLocalPlus does not: it reads a `.lanc` local-ancestry
binary, whose statistic is per ancestry and whose genotype law has a per-subject
trial count.  Nothing in the repository produced a null cohort for it, so
spa_unify Stage 7 — whose gate is a null calibration rather than a regression
baseline, SPAmixLocalPlus being unpinned — needed one.  This is it.

═════════════════════════════════════════════════════════════════════════
What the null has to be, and why the shipped fixture is not it
═════════════════════════════════════════════════════════════════════════

SPAmixLocalPlus tests S_k = sum_i R_i G_{i,k}, where G_{i,k} is subject i's
count of alternate alleles on the haplotypes RFMix assigned to ancestry k, and
h_{i,k} in {0, 1, 2} is how many haplotypes that is.  The saddlepoint models

    G_{i,k} ~ Binomial(h_{i,k}, q_k)   independently over i,

with one allele frequency q_k per (marker, ancestry).  A calibration run must
therefore supply genotypes that actually obey that law; otherwise it measures
the data, not the method.

`examples/hgdp_lanc.*` does not.  Its 925 query samples are the full HGDP panel,
including MID, OCE and AMR individuals that RFMix force-assigns to the nearest
of four continental reference ancestries, so within one ancestry class the
allele frequency still varies across subjects.  The Wahlund excess that
produces makes sum_i (G_i - q h_i)^2 exceed q(1-q) sum_i h_i, and the score is
inflated before the saddlepoint is ever reached: measured on the shipped
fixture with i.i.d. N(0,1) residuals, sd(z) = 1.0195 and lambda_GC(z^2) = 1.042.
That is a property of the cohort and would fail a [0.95, 1.05] gate for reasons
having nothing to do with the saddlepoint.

The construction here keeps the ANCESTRY plane and replaces the ALLELE plane.
Local ancestry comes verbatim from `examples/rfmix.chr{20,21,22}.msp.tsv`, so
the hapcount vectors the method sees are the real ones, with their real
imbalance across ancestries; the alleles are redrawn as independent Bernoulli
variates at a per-(marker, ancestry) frequency.  The result satisfies the model
exactly, and conditional on the residual vector the p-values are i.i.d.
Uniform(0,1), which is the frame in which a saddlepoint approximation is stated.

═════════════════════════════════════════════════════════════════════════
Outputs, and the commands that consume them
═════════════════════════════════════════════════════════════════════════

Written under OUTDIR:

    synlanc.chr{20,21,22}.vcf   phased query VCFs, simulated alleles
    hgdp.grm.id / hgdp.grm.sp   sparse GRM over the same subjects.  Diagonal
                                only by default, so the phi off-diagonal block
                                is empty, Var(S) = q(1-q) sum_i h_i R_i^2
                                exactly, and the variance-ratio rescaling
                                inside spaLocalPval is inert (varRatio == 1).
    hgdp.resid                  --resid-name input, NPHENO columns of i.i.d.
                                N(0,1).  Residuals are supplied directly rather
                                than fitted, so the null is exact rather than
                                asymptotic.

Then, with $G = build/grab2 and $O = OUTDIR:

    $G --make-lanc --vcf $O/synlanc. --rfmix-msp examples/rfmix. \
       --out $O/synlanc
    $G --cal-phi --lanc $O/synlanc --sp-grm-plink2 $O/hgdp.grm.sp \
       --out $O/synlanc
    $G --method SPAmixLocalPlus --lanc $O/synlanc \
       --admix-phi $O/synlanc.phi --pheno $O/hgdp.resid \
       --resid-name R01,R02,R03,R04,R05,R06 --out $O/cal

Every (marker x ancestry x phenotype) cell of $O/cal.R*.LocalP is one null
test.  At the defaults above that is 17200 x 4 x 6 = 412800 cells, of which
about 91 % survive the --maf / --mac filters; the standard error of lambda_GC
is then about 0.004.  Add `--spa-z-threshold 0.01` to route every cell through
the saddlepoint instead of only the |z| > 2 tail, which is the stronger test of
the solver.

Requires numpy, as tests/make_synthetic.py does.  Neither script is part of the
build; the binary remains the deliverable.
"""

import argparse
import math
import os

import numpy as np


def read_msp(path):
    """Parse one RFMix .msp.tsv into (sample IDs, [(chrom, spos, epos, calls)]).

    `calls` is one int8 per haplotype, in the file's column order, i.e.
    calls[2*s + h] for sample s haplotype h, with -1 for unassigned.
    """
    with open(path) as fh:
        fh.readline()                       # "#Subpopulation order/codes: ..."
        hdr = fh.readline().rstrip("\n").split("\t")
        cols = hdr[6:]
        samples = [cols[j].rsplit(".", 1)[0] for j in range(0, len(cols), 2)]
        wins = []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f or not f[0]:
                continue
            wins.append((f[0], int(f[1]), int(f[2]),
                         np.asarray(f[6:], dtype=np.int8)))
    return samples, wins


def write_grm(base, iids, rng, rel_frac):
    """Sparse GRM in the format examples/1kg.grm.{sp,id} uses: zero-based
    indices into the .id file, lower triangle plus diagonal."""
    n = len(iids)
    with open(base + ".grm.id", "w") as fh:
        for iid in iids:
            fh.write("0\t%s\n" % iid)
    with open(base + ".grm.sp", "w") as fh:
        for i in range(n):
            fh.write("%d\t%d\t%.8g\n" % (i, i, 1.0))
        if rel_frac > 0.0:
            npair = int(n * rel_frac / 2)
            left = rng.choice(np.arange(0, n - 1, 2),
                              size=min(npair, (n - 1) // 2), replace=False)
            for i in sorted(left):
                fh.write("%d\t%d\t%.8g\n" % (i + 1, i, 0.5))


def write_resid(base, iids, rng, npheno):
    r = rng.standard_normal((len(iids), npheno))
    names = ["R%02d" % (p + 1) for p in range(npheno)]
    with open(base + ".resid", "w") as fh:
        fh.write("#IID\t" + "\t".join(names) + "\n")
        for i, iid in enumerate(iids):
            fh.write(iid + "\t" + "\t".join("%.10g" % v for v in r[i]) + "\n")
    return names


def write_vcf(path, chrom, samples, wins, rng, per_window, k, min_af, max_af):
    """One phased VCF whose alleles obey Binomial(h, q) given the RFMix calls.

    Markers are placed strictly inside each window, `per_window` of them evenly
    spaced, so every marker has a defined ancestry call for every haplotype and
    none falls in a gap between windows.
    """
    nhap = 2 * len(samples)
    maxpos = max(w[2] for w in wins)
    written = 0
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##contig=<ID=%s,length=%d>\n" % (chrom, maxpos + 1000))
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,'
                 'Description="Genotype">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                 + "\t".join(samples) + "\n")
        for (wc, spos, epos, calls) in wins:
            if wc != chrom or calls.size != nhap:
                raise ValueError("msp window does not match %s / %d haplotypes"
                                 % (chrom, nhap))
            lo, hi = spos + 5, epos - 5
            if hi <= lo:
                continue
            step = max(1, (hi - lo) // per_window)
            positions = range(lo, min(hi, lo + step * per_window), step)
            # An unassigned haplotype (-1) is given ancestry 0's frequency; the
            # decoder counts it under no ancestry at all, so the choice is
            # immaterial and only has to be in range.
            anc = np.where(calls < 0, 0, calls).astype(np.int64)
            for pos in positions:
                q = np.exp(rng.uniform(math.log(min_af), math.log(max_af),
                                       size=k))
                a = (rng.random(nhap) < q[anc]).astype(np.int8)
                gt = np.char.add(np.char.add(a[0::2].astype("U1"), "|"),
                                 a[1::2].astype("U1"))
                fh.write("%s\t%d\t%s:%d\tA\tG\t.\t.\t.\tGT\t"
                         % (chrom, pos, chrom, pos))
                fh.write("\t".join(gt.tolist()))
                fh.write("\n")
                written += 1
    return written


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir")
    ap.add_argument("--msp-prefix", default="examples/rfmix.",
                    help="prefix of the RFMix .msp.tsv files")
    ap.add_argument("--chroms", default="chr20,chr21,chr22")
    ap.add_argument("--per-window", type=int, default=100,
                    help="markers per RFMix window")
    ap.add_argument("--k", type=int, default=4, help="ancestries in the msp")
    ap.add_argument("--npheno", type=int, default=20)
    ap.add_argument("--seed", type=int, default=20260801)
    ap.add_argument("--rel-frac", type=float, default=0.0,
                    help="fraction of subjects placed in related pairs; 0 keeps "
                         "the phi off-diagonal block empty")
    ap.add_argument("--min-af", type=float, default=0.02)
    ap.add_argument("--max-af", type=float, default=0.50)
    ap.add_argument("--prefix", default="synlanc.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    chroms = args.chroms.split(",")
    samples = None
    total = 0
    for chrom in chroms:
        s, wins = read_msp("%s%s.msp.tsv" % (args.msp_prefix, chrom))
        if samples is None:
            samples = s
        elif s != samples:
            raise ValueError("sample lists differ between msp files")
        out = os.path.join(args.outdir, "%s%s.vcf" % (args.prefix, chrom))
        total += write_vcf(out, chrom, samples, wins, rng, args.per_window,
                           args.k, args.min_af, args.max_af)
        print("%s: %d windows -> %s" % (chrom, len(wins), out))

    base = os.path.join(args.outdir, "hgdp")
    write_grm(base, samples, rng, args.rel_frac)
    names = write_resid(base, samples, rng, args.npheno)

    print("subjects=%d markers=%d phenotypes=%d" % (len(samples), total,
                                                    args.npheno))
    print("resid columns: " + ",".join(names))


if __name__ == "__main__":
    main()
