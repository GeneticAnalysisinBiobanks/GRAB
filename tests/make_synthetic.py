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
    """22 autosomes, positions strictly increasing within each chromosome."""
    per = max(1, m // 22)
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
            written += k


def write_bed(path, n, m, mafs, miss_rate, rng, chunk=512):
    """Write a variant-major PLINK1 .bed.

    PLINK1 two-bit codes, four samples per byte, least significant pair first:
        00  homozygous for the .bim A1 allele   (ALT dosage 2)
        01  missing
        10  heterozygous                        (ALT dosage 1)
        11  homozygous for the .bim A2 allele   (ALT dosage 0)
    """
    nbytes = (n + 3) // 4
    # Map ALT dosage 0,1,2 and missing(3) to the PLINK code.
    code_of = np.array([0b11, 0b10, 0b00, 0b01], dtype=np.uint8)

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


def write_pheno(path, iids, rng):
    """Phenotype table with the same columns examples/1kg.pheno provides."""
    n = len(iids)
    male = rng.integers(0, 2, size=n)
    pcs = rng.normal(0.0, 0.02, size=(n, 4))

    # No genotype effect anywhere: every trait is a function of covariates and
    # noise only, so the null holds exactly and lambda_GC is interpretable.
    lin = 0.30 * male + 3.0 * pcs[:, 0] - 2.0 * pcs[:, 1]
    quant = lin + rng.normal(0.0, 1.0, size=n)

    prob = 1.0 / (1.0 + np.exp(-(lin - 1.2)))
    binary = (rng.random(n) < prob).astype(int)

    # Exponential survival with administrative censoring at t = 5.
    rate = np.exp(0.5 * lin - 1.0)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=50000, help="subjects")
    ap.add_argument("--m", type=int, default=20000, help="markers")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--seed", type=int, default=20260730)
    ap.add_argument("--miss-rate", type=float, default=0.002)
    ap.add_argument("--min-maf", type=float, default=0.005)
    ap.add_argument("--no-grm", action="store_true")
    args = ap.parse_args()

    outdir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(outdir, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    iids = ["SYN%07d" % i for i in range(args.n)]

    # Log-uniform MAF spectrum between --min-maf and 0.5: this puts most
    # markers at low frequency, as in a real cohort, which matters because the
    # saddlepoint branch fires disproportionately at low MAF.
    lo, hi = np.log(args.min_maf), np.log(0.5)
    mafs = np.exp(rng.uniform(lo, hi, size=args.m))

    sys.stderr.write("writing %s.{bed,bim,fam,pheno} — %d subjects x %d markers\n"
                     % (args.out, args.n, args.m))
    write_fam(args.out + ".fam", iids)
    write_bim(args.out + ".bim", args.m, rng)
    write_bed(args.out + ".bed", args.n, args.m, mafs, args.miss_rate, rng)
    write_pheno(args.out + ".pheno", iids, rng)
    if not args.no_grm:
        write_sparse_grm(args.out, iids, rng)

    total = sum(os.path.getsize(args.out + ext)
                for ext in (".bed", ".bim", ".fam", ".pheno"))
    sys.stderr.write("done — %.1f MiB\n" % (total / (1 << 20)))


if __name__ == "__main__":
    main()
