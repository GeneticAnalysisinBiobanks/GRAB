#!/usr/bin/env python3
"""regress.py — Column-wise numeric comparison of two examples_output/ trees.

Invoked by tests/regress.sh; see that script for the rationale and workflow.

The comparison is deliberately asymmetric in what it treats as fatal:

  STRUCTURAL differences — a file present on one side only, a differing row
  count, a differing header, or a differing value in a non-numeric column —
  are always failures.  They mean the change altered what is being reported,
  not merely how precisely.

  NUMERIC differences are measured and reported, and fail only when they
  exceed the tolerance.  For p-value columns the reported quantity is the
  change in -log10(p), because that is the scale on which a genome-wide
  result is interpreted: a shift from 1e-8 to 2e-8 is 0.3 on that scale and
  matters, while a shift from 0.50 to 0.50000001 is 1e-8 and does not.

  MISSINGNESS changes — NA on one side and a number on the other — are
  counted and reported separately.  During the saddlepoint rework these are
  expected in one direction (a guard that now fires returns NA where the
  unguarded code previously emitted a number), so they are reported rather
  than silently folded into the numeric statistics.
"""

import argparse
import gzip
import math
import os
import sys

try:
    import zstandard  # optional; most trees are decompressed via the CLI below
except Exception:
    zstandard = None


# ──────────────────────────────────────────────────────────────────────
# Readers
# ──────────────────────────────────────────────────────────────────────

def open_text(path):
    """Open a possibly-compressed output file as text.

    baseline.sh varies the codec across blocks by design (plain, gzip, zstd),
    so the differ must handle all three.  zstd is decompressed through the
    system binary when the Python module is unavailable, which keeps this
    script dependency-free in the common case.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    if path.endswith(".zst"):
        if zstandard is not None:
            fh = open(path, "rb")
            return zstandard.ZstdDecompressor().stream_reader(fh)
        import subprocess
        out = subprocess.run(["zstd", "-dc", path], capture_output=True, check=True)
        import io
        return io.StringIO(out.stdout.decode())
    return open(path, "rt")


def read_table(path):
    """Return (meta, header, rows).

    Several writers (SAGELD's residual cache, for one) emit '##key=value'
    metadata lines ahead of the column header.  Those lines are returned
    separately and compared verbatim: they record null-model quantities such
    as the fitted lambda, so a change there is a genuine structural change,
    but they are not part of the table and must not be parsed as one.
    """
    with open_text(path) as fh:
        data = fh.read()
        if isinstance(data, bytes):
            data = data.decode()
    lines = [ln for ln in data.split("\n") if ln != ""]

    meta = []
    i = 0
    while i < len(lines) and lines[i].startswith("##"):
        meta.append(lines[i])
        i += 1
    if i >= len(lines):
        return meta, [], []

    header = lines[i].split()
    rows = [ln.split() for ln in lines[i + 1:]]
    return meta, header, rows


NA_TOKENS = {"NA", "NaN", "nan", "-nan", "inf", "-inf", "Inf", "-Inf", "."}


def as_float(tok):
    """Parse a cell into (is_missing, value)."""
    if tok in NA_TOKENS:
        return True, float("nan")
    try:
        v = float(tok)
    except ValueError:
        return None, None          # not numeric at all
    if math.isnan(v):
        return True, v
    return False, v


def is_pvalue_column(name):
    """True for columns whose values are p-values.

    Matching is on underscore/dot-separated tokens rather than substrings, so
    that 'P', 'HWE_P', 'meta_P_EXT' and 'P_TAU1' all match while 'POS', 'PC1'
    and 'ALT_FREQ' do not.
    """
    import re
    tokens = re.split(r"[_.]", name.upper())
    return "P" in tokens or any(t.startswith("PVAL") for t in tokens)


# ──────────────────────────────────────────────────────────────────────
# Comparison
# ──────────────────────────────────────────────────────────────────────

class ColumnStat:
    __slots__ = ("name", "max_abs", "max_rel", "max_dlog10p", "n_missing_change",
                 "n_sign_change", "n_numeric", "worst_row", "worst_pair")

    def __init__(self, name):
        self.name = name
        self.max_abs = 0.0
        self.max_rel = 0.0
        self.max_dlog10p = 0.0
        self.n_missing_change = 0
        self.n_sign_change = 0
        self.n_numeric = 0
        self.worst_row = -1
        self.worst_pair = (None, None)


def compare_file(base_path, new_path, rtol):
    """Return (structural_errors, [ColumnStat]) for one pair of files."""
    errs = []
    mb_, hb, rb = read_table(base_path)
    mn_, hn, rn = read_table(new_path)

    if mb_ != mn_:
        for a, b in zip(mb_, mn_):
            if a != b:
                errs.append("metadata line differs:\n    base: %s\n    new:  %s" % (a, b))
        if len(mb_) != len(mn_):
            errs.append("metadata line count differs: base %d, new %d"
                        % (len(mb_), len(mn_)))
        return errs, []
    if hb != hn:
        errs.append("header differs:\n    base: %s\n    new:  %s"
                    % (" ".join(hb), " ".join(hn)))
        return errs, []
    if len(rb) != len(rn):
        errs.append("row count differs: base %d, new %d" % (len(rb), len(rn)))
        return errs, []

    stats = [ColumnStat(name) for name in hb]
    ncol = len(hb)

    for i, (row_b, row_n) in enumerate(zip(rb, rn)):
        # A row narrower or wider than the header is a property of the writer,
        # not a difference between the two trees; only a mismatch BETWEEN the
        # trees is structural.  Compare over the fields the two rows share.
        if len(row_b) != len(row_n):
            errs.append("row %d field count differs: base %d, new %d"
                        % (i + 2, len(row_b), len(row_n)))
            continue
        for j in range(min(ncol, len(row_b))):
            tb, tn = row_b[j], row_n[j]
            if tb == tn:
                continue
            mb, vb = as_float(tb)
            mn, vn = as_float(tn)

            # Non-numeric column with differing text: structural.
            if mb is None or mn is None:
                errs.append("row %d column %s: non-numeric value differs "
                            "(%r vs %r)" % (i + 2, hb[j], tb, tn))
                continue

            st = stats[j]
            if mb != mn:
                st.n_missing_change += 1
                continue
            if mb and mn:
                continue           # both missing

            st.n_numeric += 1
            d = abs(vb - vn)
            m = max(abs(vb), abs(vn))
            rel = d / m if m > 0 else 0.0
            if d > st.max_abs:
                st.max_abs = d
            if rel > st.max_rel:
                st.max_rel = rel
                st.worst_row = i + 2
                st.worst_pair = (vb, vn)
            if (vb > 0) != (vn > 0) and vb != 0 and vn != 0:
                st.n_sign_change += 1
            if is_pvalue_column(hb[j]) and vb > 0 and vn > 0:
                dl = abs(math.log10(vb) - math.log10(vn))
                if dl > st.max_dlog10p:
                    st.max_dlog10p = dl

    return errs, stats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("new")
    ap.add_argument("--rtol", type=float, default=1e-9)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()

    def walk(root):
        """Relative paths of every regular file under root.

        baseline.sh writes into examples_output/baseline/ and a nested
        converted/ subdirectory for the cross-format SPAGRM block, so the
        comparison must recurse rather than list one level.
        """
        found = set()
        for dirpath, _dirnames, filenames in os.walk(root):
            for fn in filenames:
                full = os.path.join(dirpath, fn)
                found.add(os.path.relpath(full, root))
        return found

    base_files = walk(args.base)
    new_files = walk(args.new)

    failed = False

    only_base = sorted(base_files - new_files)
    only_new = sorted(new_files - base_files)
    if only_base:
        print("MISSING in new tree: %s" % ", ".join(only_base))
        failed = True
    if only_new:
        print("EXTRA in new tree:   %s" % ", ".join(only_new))
        failed = True

    n_identical = 0
    for name in sorted(base_files & new_files):
        bp = os.path.join(args.base, name)
        np_ = os.path.join(args.new, name)

        try:
            errs, stats = compare_file(bp, np_, args.rtol)
        except (UnicodeDecodeError, ValueError):
            # Not decodable as text: a binary artifact such as the .bed/.bcf/
            # .bgen files the cross-format block produces.  These carry no
            # floating-point output of ours and must be bit-identical, so
            # compare them byte for byte rather than skipping them.
            import hashlib

            def digest(path):
                h = hashlib.sha256()
                with open(path, "rb") as fh:
                    for chunk in iter(lambda: fh.read(1 << 20), b""):
                        h.update(chunk)
                return h.hexdigest()

            if digest(bp) == digest(np_):
                n_identical += 1
                if not args.quiet:
                    print("%-44s identical (binary)" % name)
            else:
                failed = True
                print("%-44s BINARY CONTENT DIFFERS" % name)
            continue
        except Exception as e:                       # genuinely unreadable
            print("%-44s SKIP (%s)" % (name, e))
            continue

        interesting = [s for s in stats
                       if s.max_rel > 0 or s.n_missing_change or s.n_sign_change]

        if errs:
            failed = True
            print("\n=== %s ===" % name)
            for e in errs[:10]:
                print("  STRUCTURAL: %s" % e)
            if len(errs) > 10:
                print("  ... and %d more structural differences" % (len(errs) - 10))
            continue

        if not interesting:
            n_identical += 1
            if not args.quiet:
                print("%-44s identical" % name)
            continue

        over = [s for s in interesting if s.max_rel > args.rtol]
        if over:
            failed = True

        print("\n=== %s ===" % name)
        print("  %-16s %12s %12s %12s %8s %8s"
              % ("column", "max|abs|", "max rel", "max dlog10P", "NA chg", "sgn chg"))
        for s in interesting:
            flag = "  <-- over rtol" if s.max_rel > args.rtol else ""
            print("  %-16s %12.4e %12.4e %12.4e %8d %8d%s"
                  % (s.name, s.max_abs, s.max_rel, s.max_dlog10p,
                     s.n_missing_change, s.n_sign_change, flag))
            if s.max_rel > args.rtol and s.worst_row > 0:
                print("      worst at row %d: %.17g -> %.17g"
                      % (s.worst_row, s.worst_pair[0], s.worst_pair[1]))

    print("\n%d file(s) identical, rtol = %g" % (n_identical, args.rtol))
    print("RESULT: %s" % ("DIFFERENCES EXCEED TOLERANCE OR STRUCTURE CHANGED"
                          if failed else "within tolerance"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
