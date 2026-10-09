#!/usr/bin/env python3
"""fallback_report.py — G1/G2/G4 for the D5 normal fallback.

Reads GRAB output tables and reports, for every LOG10P column that has a
matching SPA_STATUS column:

  G1  the number and rate of rows with SPA_STATUS in [3,6], overall and by
      MAC bin, the failure-mode split (3 FALLBACK_MAXITER, 4 GUARD_TEMP,
      5 GUARD_CURV, 6 NONFINITE), and a Kolmogorov-Smirnov test of the
      fallback rows' p-values against U(0,1);
  G2  the enrichment of fallback rows above -log10 p thresholds, relative to
      all rows of the same table, at 7.301 and at a ladder of lower
      thresholds (7.301 is unreachable on any null cohort of feasible size:
      the expected count is M * 5e-8);
  G4  the MAC binning that answers "is the fallback concentrated at low MAC".

and, because an enrichment on its own cannot say WHY, a two-way decomposition
(section D below):

  selection   the fallback rows are a selected subset of markers.  Comparing
              P(L > t) computed with the SAME estimator -- the normal tail at
              the printed Z_Norm -- on the fallback rows and on all rows
              isolates that selection, because the estimator is held fixed.
  estimator   on rows where the saddlepoint SUCCEEDED, both estimators are
              available: the reported LOG10P is the saddlepoint's, and the
              normal tail at the printed Z_Norm is the value the fallback
              would have reported.  Their difference, stratified by MAC and
              |Z_Norm|, is the estimator's bias with selection held fixed.
              Applying the matched-stratum median to the fallback rows gives
              the saddlepoint value they would have had.

Print resolution is a real limit and is respected: GRAB prints six
significant digits (plink2::dtoa_g, decision D3), so Z_Norm read back from a
table carries about 5e-7 relative error and every consistency check below is
made against the INTERVAL the printed Z_Norm represents, never to more digits
than were printed.

Usage:
    tests/fallback_report.py OUT.spagrm.B200.SPAGRM [more tables ...]
    tests/fallback_report.py --bins 0,10,50,200,1000,5000 FILE...
"""

import argparse
import collections
import gzip
import math
import subprocess
import sys

try:
    from scipy import stats as _st
except ImportError:                                     # pragma: no cover
    _st = None

LN10 = math.log(10.0)


def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".zst"):
        return subprocess.Popen(["zstd", "-dc", path], stdout=subprocess.PIPE,
                                text=True).stdout
    return open(path)


def read_table(path):
    fh = _open(path)
    hdr = fh.readline().rstrip("\n").split("\t")
    cols = {h: [] for h in hdr}
    for line in fh:
        t = line.rstrip("\n").split("\t")
        for h, v in zip(hdr, t):
            cols[h].append(v)
    return hdr, cols


def fnum(s):
    """Parse a printed cell; NA/nan/inf all become None."""
    if s in ("NA", "nan", "NaN", ""):
        return None
    try:
        v = float(s)
    except ValueError:
        return None
    return v if math.isfinite(v) else None


def norm_neglog10(z):
    """-log10(2 * Phi(-|z|)), in the log domain so it survives large |z|.

    math.erfc underflows at |z| ~ 27, so the asymptotic expansion takes over
    there; the crossover value agrees to 12 significant digits, which is far
    finer than the six digits any of this is compared against.
    """
    a = abs(z)
    if a < 26.0:
        return -math.log10(math.erfc(a / math.sqrt(2.0)))
    # log(2*Phi(-a)) = -a^2/2 - log(a) - log(2pi)/2 + log(1 - 1/a^2 + 3/a^4 ...)
    inv = 1.0 / (a * a)
    series = 1.0 - inv * (1.0 - 3.0 * inv * (1.0 - 5.0 * inv))
    ln = (-0.5 * a * a - math.log(a) - 0.5 * math.log(2.0 * math.pi)
          + math.log(series))
    return -ln / LN10


def pair_columns(hdr):
    """Match every SPA_STATUS* column to its LOG10P* companion.

    The two names differ only by the token SPA_STATUS vs LOG10P, which is the
    rule tests/regress.py's C3 check uses; both are applied to underscore
    segments so that prefixed and suffixed spellings pair up.
    """
    out = []
    for h in hdr:
        if "SPA_STATUS" not in h:
            continue
        cand = h.replace("SPA_STATUS", "LOG10P")
        if cand in hdr:
            out.append((cand, h))
    return out


def quantiles(xs, qs=(0, 5, 25, 50, 75, 95, 100)):
    if not xs:
        return {q: float("nan") for q in qs}
    s = sorted(xs)
    out = {}
    for q in qs:
        i = min(len(s) - 1, max(0, int(round(q / 100.0 * (len(s) - 1)))))
        out[q] = s[i]
    return out


def ks_uniform(ps):
    """KS statistic and p-value of `ps` against U(0,1)."""
    if not ps:
        return float("nan"), float("nan")
    if _st is not None:
        d, p = _st.kstest(ps, "uniform")
        return float(d), float(p)
    s = sorted(ps)
    n = len(s)
    d = max(max((i + 1) / n - v, v - i / n) for i, v in enumerate(s))
    # Kolmogorov asymptotic tail
    lam = (math.sqrt(n) + 0.12 + 0.11 / math.sqrt(n)) * d
    p = 2.0 * sum((-1) ** (k - 1) * math.exp(-2.0 * k * k * lam * lam)
                  for k in range(1, 101))
    return d, min(1.0, max(0.0, p))


STATUS_NAME = {0: "SPA_OK", 1: "NORMAL", 2: "SPA_W_SINGULAR",
               3: "FALLBACK_MAXITER", 4: "FALLBACK_GUARD_TEMP",
               5: "FALLBACK_GUARD_CURV", 6: "FALLBACK_NONFINITE",
               7: "NA_POST_FAIL", 8: "NA_NO_TEST"}


def report(path, bins, thresholds):
    hdr, cols = read_table(path)
    pairs = pair_columns(hdr)
    if not pairs:
        return
    macs = [fnum(x) for x in cols.get("MAC", [])]
    name = path.split("/")[-1]

    for lcol, scol in pairs:
        st = [fnum(x) for x in cols[scol]]
        lg = [fnum(x) for x in cols[lcol]]
        zn_name = scol.replace("SPA_STATUS", "Z_Norm")
        zn = [fnum(x) for x in cols[zn_name]] if zn_name in cols else [None] * len(st)
        n = len(st)

        cnt = collections.Counter(int(s) for s in st if s is not None)
        fb = [i for i, s in enumerate(st) if s is not None and 3 <= int(s) <= 6]
        spa_ok = [i for i, s in enumerate(st) if s is not None and int(s) == 0]
        have_l = [i for i in range(n) if lg[i] is not None]

        print("=" * 78)
        print("%s   column %s" % (name, lcol))
        print("  rows %d   LOG10P present %d   NA %d"
              % (n, len(have_l), n - len(have_l)))
        print("  SPA_STATUS: " + "  ".join(
            "%d %s=%d" % (k, STATUS_NAME.get(k, "?"), v)
            for k, v in sorted(cnt.items())))

        # ── calibration of the column as a whole (invariant 1) ────────
        if have_l:
            ps = sorted(10.0 ** (-lg[i]) for i in have_l)
            d0, p0 = ks_uniform(ps)
            med = ps[len(ps) // 2]
            # lambda_GC from the median of the 1-df chi-square implied by p
            if _st is not None:
                lam = _st.chi2.isf(med, 1) / _st.chi2.isf(0.5, 1)
            else:
                lam = float("nan")
            print("  cal LOG10P over all rows: lambda_GC = %.4f  KS D = %.4f "
                  " KS p = %.3e  #{L>7.301} = %d"
                  % (lam, d0, p0, sum(1 for i in have_l if lg[i] > 7.301)))

        # ── G1 ────────────────────────────────────────────────────────
        rate = len(fb) / n if n else float("nan")
        print("  G1  fallback rows %d / %d = %.4e" % (len(fb), n, rate))
        if not fb:
            continue
        modes = collections.Counter(int(st[i]) for i in fb)
        print("      failure modes: " + "  ".join(
            "%d %s=%d" % (k, STATUS_NAME[k], v) for k, v in sorted(modes.items())))

        fb_l = [lg[i] for i in fb if lg[i] is not None]
        fb_p = [10.0 ** (-v) for v in fb_l]
        d, kp = ks_uniform(fb_p)
        print("      KS of fallback p-values vs U(0,1): D = %.4f  p = %.3e"
              % (d, kp))
        print("      (a rejection here is EXPECTED and is not miscalibration:"
              " the fallback set is")
        print("       selected on |Z|, so its p-values cannot be uniform"
              " whatever estimator is used)")
        qs = quantiles(fb_l)
        print("      fallback LOG10P quantiles  " +
              "  ".join("q%d=%.4g" % (q, qs[q]) for q in sorted(qs)))
        qz = quantiles([abs(zn[i]) for i in fb if zn[i] is not None])
        print("      fallback |Z_Norm| quantiles " +
              "  ".join("q%d=%.4g" % (q, qz[q]) for q in sorted(qz)))

        # consistency of the fallback value with the normal tail, at the
        # resolution the table was printed with
        bad = 0
        worst = 0.0
        for i in fb:
            if zn[i] is None or lg[i] is None:
                continue
            z = abs(zn[i])
            # Six significant digits -> half-ulp of the printed decimal, on
            # BOTH columns: Z_Norm was rounded before we read it, and LOG10P
            # was rounded after it was computed, so the admissible interval is
            # the image of the Z_Norm interval widened by LOG10P's own ulp.
            half = 0.5 * 10.0 ** (math.floor(math.log10(z)) - 5) if z > 0 else 0.0
            # -log10 p increases with |z|, so the printed interval maps to
            # [f(z - half), f(z + half)].
            lo = norm_neglog10(max(0.0, z - half))
            hi = norm_neglog10(z + half)
            lh = (0.5 * 10.0 ** (math.floor(math.log10(abs(lg[i]))) - 5)
                  if lg[i] != 0.0 else 0.0)
            if not (lo - lh <= lg[i] <= hi + lh):
                bad += 1
            ref = norm_neglog10(z)
            if ref > 0:
                worst = max(worst, abs(lg[i] - ref) / ref)
        print("      fallback LOG10P vs -log10(2*Phi(-|Z_Norm|)): "
              "%d / %d outside the printed interval, max rel %.3e"
              % (bad, len(fb), worst))

        # ── G2 ────────────────────────────────────────────────────────
        print("  G2  enrichment above -log10 p thresholds"
              " (fallback rate vs all-row rate)")
        for t in thresholds:
            a = sum(1 for i in have_l if lg[i] > t)
            b = sum(1 for i in fb if lg[i] is not None and lg[i] > t)
            pa = a / len(have_l) if have_l else float("nan")
            pb = b / len(fb) if fb else float("nan")
            enr = (pb / pa) if pa > 0 else float("inf") if pb > 0 else float("nan")
            print("      L > %-6.3f  all %6d/%6d = %.3e   fallback %5d/%5d "
                  "= %.3e   enrichment %s"
                  % (t, a, len(have_l), pa, b, len(fb), pb,
                     ("%.3f" % enr) if enr == enr else "n/a (0/0)"))
        print("      max LOG10P: all %.4g   fallback %.4g"
              % (max(lg[i] for i in have_l), max(fb_l)))

        # ── G4 ────────────────────────────────────────────────────────
        if macs:
            print("  G4  MAC bins")
            print("      %-16s %8s %8s %10s" % ("bin", "n_all", "n_fb", "rate"))
            edges = list(bins) + [float("inf")]
            for k in range(len(edges) - 1):
                lo, hi = edges[k], edges[k + 1]
                idx = [i for i in range(n)
                       if macs[i] is not None and lo <= macs[i] < hi]
                nf = sum(1 for i in idx if st[i] is not None
                         and 3 <= int(st[i]) <= 6)
                r = nf / len(idx) if idx else float("nan")
                print("      [%-6g,%-7g) %8d %8d %10.3e"
                      % (lo, hi, len(idx), nf, r))

        # ── D  selection vs estimator ─────────────────────────────────
        print("  D   decomposition of the enrichment")
        # (a) selection: same estimator (the normal tail) on both row sets
        all_spa = [i for i in range(n)
                   if st[i] is not None and int(st[i]) in (0, 2)
                   or (st[i] is not None and 3 <= int(st[i]) <= 6)]
        for t in thresholds:
            na = sum(1 for i in have_l
                     if zn[i] is not None and norm_neglog10(zn[i]) > t)
            nb = sum(1 for i in fb
                     if zn[i] is not None and norm_neglog10(zn[i]) > t)
            pa = na / len(have_l) if have_l else float("nan")
            pb = nb / len(fb) if fb else float("nan")
            enr = (pb / pa) if pa > 0 else float("nan")
            print("      selection  L_normal > %-6.3f  all %.3e  fallback %.3e"
                  "  ratio %s"
                  % (t, pa, pb, ("%.3f" % enr) if enr == enr else "n/a"))
        # (b) estimator: on rows where the saddlepoint succeeded, both values
        deltas = []
        for i in spa_ok:
            if lg[i] is None or zn[i] is None:
                continue
            deltas.append(lg[i] - norm_neglog10(zn[i]))
        if deltas:
            qd = quantiles(deltas)
            print("      estimator  on %d SPA_OK rows, "
                  "LOG10P_saddlepoint - LOG10P_normal:" % len(deltas))
            print("                 " +
                  "  ".join("q%d=%+.4g" % (q, qd[q]) for q in sorted(qd)))
            print("                 (negative => the saddlepoint is LESS "
                  "significant than the normal,")
            print("                  i.e. the normal fallback would have been"
                  " anti-conservative)")
            if macs:
                print("                 by MAC bin (median delta, n):")
                edges = list(bins) + [float("inf")]
                for k in range(len(edges) - 1):
                    lo_, hi_ = edges[k], edges[k + 1]
                    dd = [lg[i] - norm_neglog10(zn[i]) for i in spa_ok
                          if lg[i] is not None and zn[i] is not None
                          and macs[i] is not None and lo_ <= macs[i] < hi_]
                    if dd:
                        print("                   [%-6g,%-7g) median %+.4g  n=%d"
                              % (lo_, hi_, sorted(dd)[len(dd) // 2], len(dd)))
        # (c) matched: apply the SPA_OK stratum median delta to fallback rows.
        # Strata are tried from most to least specific so that a fallback row
        # is never dropped for want of an exactly-matching SPA_OK row; the
        # coarsest level is the global median, which is the estimate the
        # column-level quantiles above already report.
        if deltas and macs:
            fine = collections.defaultdict(list)
            coarse = collections.defaultdict(list)

            def macbin(i):
                m = macs[i] if macs[i] is not None else -1
                mb = 0
                for e in bins:
                    if m >= e:
                        mb = e
                return mb

            def zbin(i):
                return round(abs(zn[i]) * 2) / 2 if zn[i] is not None else -1

            for i in spa_ok:
                if lg[i] is None or zn[i] is None:
                    continue
                dv = lg[i] - norm_neglog10(zn[i])
                fine[(macbin(i), zbin(i))].append(dv)
                coarse[zbin(i)].append(dv)
            gmed = sorted(deltas)[len(deltas) // 2]
            matched, unmatched = [], 0
            for i in fb:
                if lg[i] is None or zn[i] is None:
                    continue
                d_ = fine.get((macbin(i), zbin(i))) or coarse.get(zbin(i))
                med = sorted(d_)[len(d_) // 2] if d_ else gmed
                if not d_:
                    unmatched += 1
                matched.append(lg[i] + med)
            if matched:
                qm = quantiles(matched)
                print("      matched    imputed saddlepoint LOG10P for %d "
                      "fallback rows (%d fell back to the global median):"
                      % (len(matched), unmatched))
                print("                 " +
                      "  ".join("q%d=%.4g" % (q, qm[q]) for q in sorted(qm)))
                for t in thresholds:
                    b1 = sum(1 for v in fb_l if v > t)
                    b2 = sum(1 for v in matched if v > t)
                    print("                 L > %-6.3f  reported(normal) %d"
                          "   imputed(saddlepoint) %d" % (t, b1, b2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--bins", default="0,5,10,20,50,100,200,500,1000,2000",
                    help="MAC bin lower edges")
    ap.add_argument("--thresholds", default="1,2,3,4,5,7.301",
                    help="-log10 p thresholds for the G2 ladder")
    args = ap.parse_args()
    bins = [float(x) for x in args.bins.split(",") if x.strip()]
    ths = [float(x) for x in args.thresholds.split(",") if x.strip()]
    for f in args.files:
        report(f, bins, ths)


if __name__ == "__main__":
    main()
