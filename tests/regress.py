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

  INVARIANT violations on the new tree — the C1/C2/C3 checks of
  dev-notes/methods/log10p_unify/04_validation.md §2 — are always failures.
  They are properties of one tree rather than of a difference between two, so
  they are reported in their own section at the end and fire whether or not
  anything changed.

Two p-value representations coexist while the log10p_unify project is in
flight, and this script must read both:

  * a LINEAR p-value column ('P', 'HWE_P', 'meta_P_EXT', 'P_tau0.5', ...),
    whose -log10 scale value is -log10(v); and
  * a LOG10P column ('LOG10P', 'LOG10P_EXT', 'cl3_LOG10P_NOEXT',
    'LOG10P_tau0.5', 'anc2_LOG10P', 'LOG10P_CCT', ...), which already *is*
    -log10(p).

'max dlog10P' is the change on the -log10(p) scale regardless of which
representation a column uses.  Taking log10 of a LOG10P value, which is what
this script did before the log10p_unify project, would report a change of
1e-9 for a marker that moved from 300 to 301 and would report exactly zero
once the linear P column is deleted in Stage 8 — a clean-looking report that
measures nothing.

Note that '+Inf' and '-Inf' are treated as missing by the numeric comparison
(see NA_TOKENS).  That conflation is deliberate for the difference report, and
the C1 invariant below is what makes an infinity in a LOG10P column a failure
rather than a silently-counted change in missingness.
"""

import argparse
import gzip
import math
import os
import re
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


# ──────────────────────────────────────────────────────────────────────
# Column classification
# ──────────────────────────────────────────────────────────────────────
#
# Matching is on underscore/dot-separated tokens rather than substrings, so
# that 'P', 'HWE_P', 'meta_P_EXT' and 'P_tau0.5' all match as linear p-values
# while 'POS', 'PC1', 'PI_BAT' and 'ALT_FREQ' do not, and so that 'LOG10P',
# 'LOG10P_EXT', 'cl3_LOG10P_NOEXT', 'LOG10P_tau0.5', 'anc2_LOG10P',
# 'LOG10P_Gx<E>', 'LOG10P_CCT', 'LOG10P_BAT' and 'LOG10P_HWE' all match as
# -log10(p) columns.

def _tokens(name):
    return re.split(r"[_.]", name.upper())


def is_log10p_column(name):
    """True for columns whose values are -log10(p) (decision D2's spelling)."""
    return "LOG10P" in _tokens(name)


def is_linear_pvalue_column(name):
    """True for columns whose values are p-values on the linear scale."""
    toks = _tokens(name)
    if "LOG10P" in toks:
        return False
    return "P" in toks or any(t.startswith("PVAL") for t in toks)


def is_pvalue_column(name):
    """True for any column carrying a p-value, in either representation."""
    return is_log10p_column(name) or is_linear_pvalue_column(name)


def neg_log10_p(name, value):
    """The -log10(p) value of a cell, or None when it is not defined.

    A LOG10P column already holds -log10(p) and is returned unchanged; a
    linear p-value column is converted.  A linear p of exactly 0 has no
    -log10 and returns None (it is counted separately as an exact zero); a
    non-finite value in either representation returns None.
    """
    if value is None or not math.isfinite(value):
        return None
    if is_log10p_column(name):
        return value
    if value > 0.0:
        return -math.log10(value)
    return None


# ──────────────────────────────────────────────────────────────────────
# Explicit column pairing (--pair)
# ──────────────────────────────────────────────────────────────────────
#
# Stage 8 of the log10p_unify project deletes the linear 'P' column, so a
# marker's p-value lives in 'P' on the base side and in 'LOG10P' on the new
# side and the two are no longer comparable by name.  '--pair P:LOG10P'
# declares the correspondence once, as a rule over underscore-separated name
# segments rather than a per-column enumeration, so that a single spec covers
# every prefixed and suffixed variant:
#
#     P            -> LOG10P
#     P_EXT        -> LOG10P_EXT
#     meta_P_NOEXT -> meta_LOG10P_NOEXT
#     cl3_P_BAT    -> cl3_LOG10P_BAT
#     P_tau0.5     -> LOG10P_tau0.5
#     anc2_P       -> anc2_LOG10P
#     P_CCT        -> LOG10P_CCT
#     P_Wald_GxMALE-> LOG10P_Wald_GxMALE
#
# A spec whose left-hand side matches a whole column name is also honored
# verbatim, which is how a rename that is not a segment substitution — say
# '--pair HWE_P:LOG10P_HWE' for decision D7 — is expressed.

def paired_name(name, spec):
    """Apply one 'BASE:NEW' pairing rule to a base column name.

    Returns the corresponding new-tree column name, or None when the rule
    does not apply to this column.
    """
    lhs, rhs = spec
    if name == lhs:
        return rhs
    segs = name.split("_")
    lseg = lhs.split("_")
    n = len(lseg)
    for i in range(len(segs) - n + 1):
        if segs[i:i + n] == lseg:
            return "_".join(segs[:i] + rhs.split("_") + segs[i + n:])
    return None


# ──────────────────────────────────────────────────────────────────────
# Invariant checks C1-C3 (04_validation.md §2, invariant 4)
# ──────────────────────────────────────────────────────────────────────
#
# These are properties of a single tree, checked on the NEW tree, and they
# fail the run.  They exist because the failure mode this branch is removing
# is a quantity that a reader would mistake for a measurement: an infinity
# that came from -log10 of an underflowed ratio, a negative -log10(p), or a
# p-value present on a row whose status says no p-value exists.
#
#   C1  no LOG10P-family column holds +Inf or -Inf
#   C2  no LOG10P-family column holds a negative value
#   C3  SPA_STATUS and its matching LOG10P agree about whether a p exists
#
# C3 depends on which SPA_STATUS partition is in force.  Two are selectable:
#
#   legacy7  the seven-value enumeration in force before Stage 2.  Statuses
#            0 (OK), 4 (GUARD_W) and 6 (NORMAL) carry a p-value; 1 (MAXITER),
#            2 (GUARD_TEMP), 3 (GUARD_CURV) and 5 (NONFINITE) are NA.
#   nine     the nine-value partition of decision D4.  Statuses 0-6 carry a
#            p-value (3-6 being the named normal-approximation fallback);
#            7 (NA_POST_FAIL) and 8 (NA_NO_TEST) are NA.
#
# The default is 'nine', which is what the tree has emitted since log10p_unify
# Stage 2 re-partitioned spa::Status.  'legacy7' is kept only so that an
# examples_output tree produced by a pre-Stage-2 binary can still be checked;
# reading a Stage-2 tree under 'legacy7' reports every correctly substituted
# fallback row (SPA_STATUS 3-6 with a numeric LOG10P) as a C3 violation.

STATUS_SEMANTICS = {
    # name: (set of statuses that must carry a finite non-negative LOG10P,
    #        set of statuses whose LOG10P must be NA)
    "legacy7": ({0, 4, 6}, {1, 2, 3, 5}),
    "nine":    ({0, 1, 2, 3, 4, 5, 6}, {7, 8}),
}


def log10p_partner_of_status(name):
    """The LOG10P column that a SPA_STATUS column qualifies, by name.

    'SPA_STATUS' -> 'LOG10P', 'SPA_STATUS_EXT' -> 'LOG10P_EXT',
    'cl3_SPA_STATUS_NOEXT' -> 'cl3_LOG10P_NOEXT',
    'SPA_STATUS_tau0.5' -> 'LOG10P_tau0.5', 'anc2_SPA_STATUS' -> 'anc2_LOG10P',
    'meta_SPA_STATUS_EXT' -> 'meta_LOG10P_EXT'.
    """
    if "SPA_STATUS" not in name:
        return None
    return name.replace("SPA_STATUS", "LOG10P", 1)


class Violation:
    __slots__ = ("check", "path", "column", "count", "example")

    def __init__(self, check, path, column, count, example):
        self.check = check
        self.path = path
        self.column = column
        self.count = count
        self.example = example


def scan_invariants(root, semantics, quiet):
    """Run C1/C2/C3 over every parseable table under root.

    Returns (violations, unpaired_status_columns).
    """
    violations = []
    unpaired = set()
    carries_p, must_be_na = STATUS_SEMANTICS[semantics]

    for dirpath, _dirnames, filenames in sorted(os.walk(root)):
        for fn in sorted(filenames):
            full = os.path.join(dirpath, fn)
            rel = os.path.relpath(full, root)
            try:
                _meta, header, rows = read_table(full)
            except Exception:
                continue                       # binary artifact or unreadable
            if not header:
                continue

            for j, name in enumerate(header):
                if not is_log10p_column(name):
                    continue
                n_inf = n_neg = 0
                ex_inf = ex_neg = None
                for i, row in enumerate(rows):
                    if j >= len(row):
                        continue
                    tok = row[j]
                    if tok in NA_TOKENS and tok not in ("inf", "-inf", "Inf", "-Inf"):
                        continue
                    try:
                        v = float(tok)
                    except ValueError:
                        continue
                    if math.isinf(v):
                        n_inf += 1
                        if ex_inf is None:
                            ex_inf = "row %d: %s" % (i + 2, tok)
                    elif not math.isnan(v) and v < 0.0:
                        n_neg += 1
                        if ex_neg is None:
                            ex_neg = "row %d: %s" % (i + 2, tok)
                if n_inf:
                    violations.append(Violation("C1", rel, name, n_inf, ex_inf))
                if n_neg:
                    violations.append(Violation("C2", rel, name, n_neg, ex_neg))

            for j, name in enumerate(header):
                partner = log10p_partner_of_status(name)
                if partner is None:
                    continue
                if partner not in header:
                    unpaired.add("%s: %s (want %s)" % (rel, name, partner))
                    continue
                jp = header.index(partner)
                n_bad = 0
                example = None
                for i, row in enumerate(rows):
                    if j >= len(row) or jp >= len(row):
                        continue
                    st_tok, p_tok = row[j], row[jp]
                    try:
                        st = int(float(st_tok))
                    except ValueError:
                        continue
                    p_missing, p_val = as_float(p_tok)
                    if p_missing is None:
                        continue
                    bad = False
                    if st in must_be_na and not p_missing:
                        bad = True
                    elif st in carries_p:
                        if p_missing or not math.isfinite(p_val) or p_val < 0.0:
                            bad = True
                    if bad:
                        n_bad += 1
                        if example is None:
                            example = "row %d: %s=%s, %s=%s" % (
                                i + 2, name, st_tok, partner, p_tok)
                if n_bad:
                    violations.append(
                        Violation("C3", rel, "%s / %s" % (name, partner),
                                  n_bad, example))

    return violations, sorted(unpaired)


# ──────────────────────────────────────────────────────────────────────
# Comparison
# ──────────────────────────────────────────────────────────────────────

class ColumnStat:
    __slots__ = ("name", "max_abs", "max_rel", "max_dlog10p", "n_missing_change",
                 "n_na_gained", "n_na_lost", "n_sign_change", "n_numeric",
                 "n_text_change", "text_trans", "n_zero_base", "n_zero_new",
                 "worst_row", "worst_pair", "worst_dlog_row", "worst_dlog_pair",
                 "diagnostic")

    def __init__(self, name, diagnostic=False):
        self.name = name
        self.max_abs = 0.0
        self.max_rel = 0.0
        self.max_dlog10p = 0.0
        self.n_missing_change = 0
        self.n_na_gained = 0        # a number in base, NA in new
        self.n_na_lost = 0          # NA in base, a number in new
        self.n_sign_change = 0
        self.n_numeric = 0
        self.n_text_change = 0
        self.text_trans = {}        # "old -> new" -> count
        self.n_zero_base = 0
        self.n_zero_new = 0
        self.worst_row = -1
        self.worst_pair = (None, None)
        self.worst_dlog_row = -1
        self.worst_dlog_pair = (None, None)
        # A diagnostic column is an explicit --pair comparison across two
        # differently-named columns.  Its absolute and relative differences
        # are meaningless (the two sides are on different scales), so it is
        # reported for its max dlog10P and never graded against --rtol.
        self.diagnostic = diagnostic

    def interesting(self):
        # A column that carries an exactly-zero p-value in BOTH trees has not
        # changed; only a change in that count is a signal.  (HWE_P legitimately
        # holds a zero in several fixtures.)
        if self.diagnostic:
            return True
        return (self.max_rel > 0 or self.n_missing_change or self.n_sign_change
                or self.n_text_change or self.n_zero_base != self.n_zero_new)


def compare_file(base_path, new_path, rtol, pair_specs):
    """Return (structural_errors, [ColumnStat]) for one pair of files.

    A differing header is reported as a structural difference AND the columns the
    two trees share are still compared.  The saddlepoint rework adds LOG10P and
    SPA_STATUS to every SPA method one stage at a time, and each of those stages
    is graded on a maximum change in -log10(P); returning no statistics the
    moment a column is added would withhold exactly the number the acceptance
    gate is stated in.  Nothing is relaxed by this: an added or removed column is
    still structural and still fails.  A --pair rule adds a comparison across
    two differently-named columns; it does not suppress that structural finding,
    it only restores the ability to state how far the numbers moved.
    """
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

    # entries: (label, name_base, name_new, index_base, index_new, diagnostic)
    if hb != hn:
        added = [c for c in hn if c not in hb]
        removed = [c for c in hb if c not in hn]
        shared_b = [c for c in hb if c in hn]
        shared_n = [c for c in hn if c in hb]
        if added:
            errs.append("columns added: %s" % ", ".join(added))
        if removed:
            errs.append("columns removed: %s" % ", ".join(removed))
        if shared_b != shared_n:
            errs.append("shared columns reordered:\n    base: %s\n    new:  %s"
                        % (" ".join(shared_b), " ".join(shared_n)))
        if not added and not removed and shared_b == shared_n:
            errs.append("header differs:\n    base: %s\n    new:  %s"
                        % (" ".join(hb), " ".join(hn)))
        entries = [(c, c, c, hb.index(c), hn.index(c), False) for c in shared_b]
    else:
        entries = [(c, c, c, j, j, False) for j, c in enumerate(hb)]

    # Explicit pairings apply to base columns the new tree no longer has.
    for cb in hb:
        if cb in hn:
            continue
        for spec in pair_specs:
            cn = paired_name(cb, spec)
            if cn is not None and cn in hn:
                entries.append(("%s->%s" % (cb, cn), cb, cn,
                                hb.index(cb), hn.index(cn), True))
                break

    if not entries:
        return errs, []

    if len(rb) != len(rn):
        errs.append("row count differs: base %d, new %d" % (len(rb), len(rn)))
        return errs, []

    stats = [ColumnStat(label, diag) for label, _nb, _nn, _jb, _jn, diag in entries]

    for i, (row_b, row_n) in enumerate(zip(rb, rn)):
        # A row narrower or wider than the header is a property of the writer,
        # not a difference between the two trees; only a mismatch BETWEEN the
        # trees is structural.  Skip that check when the headers differ, since
        # the row widths legitimately differ then.
        if len(row_b) != len(row_n) and hb == hn:
            errs.append("row %d field count differs: base %d, new %d"
                        % (i + 2, len(row_b), len(row_n)))
            continue
        for k, (_label, nb, nn, jb, jn, diag) in enumerate(entries):
            if jb >= len(row_b) or jn >= len(row_n):
                continue
            tb, tn = row_b[jb], row_n[jn]
            st = stats[k]

            mb, vb = as_float(tb)
            mn, vn = as_float(tn)

            # Count exact-zero p-values whether or not the cell changed.  An
            # emitted p of exactly 0.0 is a defect in its own right and log10
            # cannot report it, so it is invisible in max dlog10P.  A LOG10P of
            # exactly 0 is p = 1 and is not a defect, so only linear p-value
            # columns are counted here.
            if is_linear_pvalue_column(nb) and mb is False and vb == 0.0:
                st.n_zero_base += 1
            if is_linear_pvalue_column(nn) and mn is False and vn == 0.0:
                st.n_zero_new += 1

            if tb == tn and not diag:
                continue

            # Non-numeric column with differing text: structural, aggregated
            # into a transition histogram rather than one error line per row, so
            # that a column such as SPA_STATUS can actually be enumerated.
            if mb is None or mn is None:
                if not diag:
                    st.n_text_change += 1
                    key = "%s -> %s" % (tb, tn)
                    st.text_trans[key] = st.text_trans.get(key, 0) + 1
                continue

            if mb != mn:
                st.n_missing_change += 1
                if mn:
                    st.n_na_gained += 1
                else:
                    st.n_na_lost += 1
                continue
            if mb and mn:
                continue           # both missing

            st.n_numeric += 1

            # The change on the -log10(p) scale, computed from whichever
            # representation each side uses.  For two LOG10P columns this is
            # |v_base - v_new|; for two linear P columns it is
            # |log10(v_base) - log10(v_new)|, as before; for a --pair across
            # the two representations the base side is converted and the new
            # side is taken as it stands.
            if is_pvalue_column(nb) and is_pvalue_column(nn):
                lb = neg_log10_p(nb, vb)
                ln = neg_log10_p(nn, vn)
                if lb is not None and ln is not None:
                    dl = abs(lb - ln)
                    if dl > st.max_dlog10p:
                        st.max_dlog10p = dl
                        st.worst_dlog_row = i + 2
                        st.worst_dlog_pair = (vb, vn)

            if diag:
                continue           # abs/rel/sign are meaningless across scales

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

    return errs, stats


def parse_pair(spec):
    parts = spec.split(":")
    if len(parts) != 2 or not parts[0] or not parts[1]:
        raise argparse.ArgumentTypeError(
            "--pair expects BASE_COLUMN:NEW_COLUMN, for example P:LOG10P")
    return (parts[0], parts[1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("new")
    ap.add_argument("--rtol", type=float, default=1e-9)
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--pair", type=parse_pair, action="append", default=[],
                    metavar="BASE:NEW",
                    help="compare a base-tree column against a differently "
                         "named new-tree column, as a rule over underscore-"
                         "separated name segments (P:LOG10P also pairs P_EXT "
                         "with LOG10P_EXT and cl3_P_BAT with cl3_LOG10P_BAT). "
                         "May be given more than once.")
    ap.add_argument("--spa-status-semantics", choices=sorted(STATUS_SEMANTICS),
                    default="nine",
                    help="which SPA_STATUS partition the C3 invariant is "
                         "checked against (default: nine, the partition in "
                         "force since log10p_unify Stage 2; legacy7 is the "
                         "seven-value enumeration that preceded it)")
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
            errs, stats = compare_file(bp, np_, args.rtol, args.pair)
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

        interesting = [s for s in stats if s.interesting()]

        if errs:
            failed = True
        if not errs and not interesting:
            n_identical += 1
            if not args.quiet:
                print("%-44s identical" % name)
            continue

        over = [s for s in interesting
                if not s.diagnostic and s.max_rel > args.rtol]
        if over or any(s.n_text_change for s in interesting):
            failed = True

        print("\n=== %s ===" % name)
        for e in errs[:10]:
            print("  STRUCTURAL: %s" % e)
        if len(errs) > 10:
            print("  ... and %d more structural differences" % (len(errs) - 10))

        if interesting:
            print("  %-24s %12s %12s %12s %9s %8s %8s"
                  % ("column", "max|abs|", "max rel", "max dlog10P",
                     "NA gained", "NA lost", "sgn chg"))
            for s in interesting:
                flag = ("  <-- over rtol"
                        if (not s.diagnostic and s.max_rel > args.rtol) else "")
                if s.diagnostic:
                    flag = "  (paired, dlog10P only)"
                print("  %-24s %12.4e %12.4e %12.4e %9d %8d %8d%s"
                      % (s.name, s.max_abs, s.max_rel, s.max_dlog10p,
                         s.n_na_gained, s.n_na_lost, s.n_sign_change, flag))
                if not s.diagnostic and s.max_rel > args.rtol and s.worst_row > 0:
                    print("      worst rel at row %d: %.17g -> %.17g"
                          % (s.worst_row, s.worst_pair[0], s.worst_pair[1]))
                if s.max_dlog10p > 0 and s.worst_dlog_row > 0:
                    print("      worst dlog10P at row %d: %.6g -> %.6g"
                          % (s.worst_dlog_row, s.worst_dlog_pair[0],
                             s.worst_dlog_pair[1]))
                if s.n_zero_base or s.n_zero_new:
                    print("      p == 0 exactly: base %d, new %d"
                          % (s.n_zero_base, s.n_zero_new))
                if s.n_text_change:
                    trans = sorted(s.text_trans.items(), key=lambda kv: -kv[1])
                    print("      TEXT CHANGED in %d row(s): %s"
                          % (s.n_text_change,
                             "; ".join("%s (x%d)" % (k, v) for k, v in trans[:8])))
                    if len(trans) > 8:
                        print("      ... and %d more distinct transitions"
                              % (len(trans) - 8))

    print("\n%d file(s) identical, rtol = %g" % (n_identical, args.rtol))

    # ── Invariants on the new tree (04_validation.md §2, C1-C3) ──────────
    violations, unpaired = scan_invariants(args.new, args.spa_status_semantics,
                                           args.quiet)
    print("\n=== invariant checks on %s (SPA_STATUS semantics: %s) ==="
          % (args.new, args.spa_status_semantics))
    for check, desc in (("C1", "LOG10P-family column holds +Inf or -Inf"),
                        ("C2", "LOG10P-family column holds a negative value"),
                        ("C3", "SPA_STATUS disagrees with its LOG10P about "
                               "whether a p-value exists")):
        hits = [v for v in violations if v.check == check]
        total = sum(v.count for v in hits)
        print("  %s  %-62s %s" % (check, desc,
                                  "PASS" if not hits else
                                  "FAIL (%d cell(s) in %d column(s))"
                                  % (total, len(hits))))
        for v in hits:
            print("        %s  %s  x%d   e.g. %s"
                  % (v.path, v.column, v.count, v.example))
    if unpaired:
        print("  SPA_STATUS columns with no matching LOG10P column "
              "(C3 not checked for these):")
        for u in unpaired:
            print("        %s" % u)
    if violations:
        failed = True

    print("\nRESULT: %s" % ("DIFFERENCES EXCEED TOLERANCE, STRUCTURE CHANGED, "
                            "OR AN INVARIANT WAS VIOLATED"
                            if failed else "within tolerance"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
