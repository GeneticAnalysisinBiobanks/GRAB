#!/usr/bin/env bash
# regress.sh — Tolerance-aware comparison of two examples_output/ trees.
#
# The saddlepoint rework changes results by design: the cancellation-free CGF,
# the unified root finder, and the repaired guards all move the last digits of
# every p-value, and some markers change status outright.  A byte-for-byte
# `cmp` is therefore the wrong instrument — it reports "different" for every
# file and tells you nothing about whether the difference is a rounding shift
# or a regression.
#
# This script instead reports, per file and per numeric column:
#
#   * the maximum absolute and relative difference,
#   * the number of cells where one side is NA/NaN and the other is not,
#   * the number of cells whose sign differs,
#   * for any column whose name looks like a p-value, the largest change in
#     -log10(p), which is the scale on which a GWAS result is actually read.
#     Both representations are understood: a linear P column is converted,
#     a LOG10P column is differenced as it stands.
#
# Non-numeric columns (IDs, alleles, chromosome labels) are compared exactly;
# any difference there is a structural change and is reported as such.
#
# It also runs the C1-C3 invariant checks of
# dev-notes/methods/log10p_unify/04_validation.md §2 over the NEW tree — no
# infinity and no negative value in a LOG10P column, and agreement between
# SPA_STATUS and its LOG10P about whether a p-value exists — and fails on a
# violation.  Those are properties of one tree, so they fire whether or not
# anything changed between the two.
#
# Binary artifacts (.bed / .bcf / .bgen, from the cross-format SPAGRM block)
# are compared byte for byte, with ONE documented exception.  `plink2 --export
# bcf` stamps its VCF header with '##fileDate=YYYYMMDD', so two runs on
# different days differ on an unmodified tree; for a BCF the comparison
# therefore decompresses the BGZF container and normalizes exactly eight ASCII
# digits after a literal '##fileDate=' inside the file's own declared
# header-text region.  Every other byte — the magic, l_text, every other header
# line, every variant and genotype record — still has to match, and the
# replacement is length preserving so a header whose length changed still
# fails.  `python3 tests/regress.py --self-test` proves both halves of that on
# constructed inputs and is run below before every comparison.
#
# That is one of three per-artifact exceptions in this repository, and they are
# not interchangeable:
#   baseline/converted/1kg.log  plink2's log: a wall-clock timestamp AND a
#                               memory estimate, not reproducible at all,
#                               discounted BY HAND (CLAUDE.md).
#   baseline/fit.ibd.zst        thread-order non-determinism in
#                               src/spagrm/ibd.cpp; excluded from the
#                               REPRODUCIBILITY gate by ruling R4, and compared
#                               normally between trees.
#   baseline/converted/1kg.bcf  the ##fileDate normalization above.
#
# Usage:
#     tests/regress.sh BASE_DIR NEW_DIR [--rtol 1e-9] [--quiet]
#                      [--pair BASE_COL:NEW_COL]...
#                      [--spa-status-semantics legacy7|nine]
#
# --pair restores comparability across a column rename: '--pair P:LOG10P'
# compares the base tree's P against the new tree's LOG10P and reports the
# change on the -log10(p) scale, and being a rule over underscore-separated
# name segments it covers P_EXT/LOG10P_EXT, cl3_P_BAT/cl3_LOG10P_BAT and the
# rest without per-column enumeration.  It does not suppress the structural
# finding that a column was removed.
#
# Exit status is 0 when every difference is within tolerance, no structural
# change was found, and no invariant was violated; 1 otherwise.  The intended
# workflow between stages is:
#
#     make -j && bash examples/baseline.sh && mv examples_output examples_output.base
#     # ... make a change ...
#     make -j && bash examples/baseline.sh
#     tests/regress.sh examples_output.base examples_output

set -uo pipefail

BASE=""
NEW=""
RTOL="1e-9"
QUIET=0
EXTRA=()

while [ $# -gt 0 ]; do
    case "$1" in
        --rtol)  RTOL="$2"; shift 2 ;;
        --quiet) QUIET=1; shift ;;
        --pair)  EXTRA+=(--pair "$2"); shift 2 ;;
        --spa-status-semantics)
                 EXTRA+=(--spa-status-semantics "$2"); shift 2 ;;
        --self-test)
            exec "${PYTHON:-python3}" "$(dirname "$0")/regress.py" --self-test ;;
        -h|--help)
            sed -n '2,70p' "$0"; exit 0 ;;
        *)
            if   [ -z "$BASE" ]; then BASE="$1"
            elif [ -z "$NEW"  ]; then NEW="$1"
            else echo "unexpected argument: $1" >&2; exit 2
            fi
            shift ;;
    esac
done

if [ -z "$BASE" ] || [ -z "$NEW" ]; then
    echo "usage: $0 BASE_DIR NEW_DIR [--rtol RTOL] [--quiet]" \
         "[--pair BASE_COL:NEW_COL]... [--spa-status-semantics legacy7|nine]" >&2
    exit 2
fi
if [ ! -d "$BASE" ] || [ ! -d "$NEW" ]; then
    echo "error: both arguments must be existing directories" >&2
    exit 2
fi

# Locate a python3 for the comparison core.  The repository ships no Python
# dependency; this is a developer-only regression aid, not part of the build.
PY="${PYTHON:-python3}"
if ! command -v "$PY" >/dev/null 2>&1; then
    echo "error: python3 not found; set PYTHON=/path/to/python3" >&2
    exit 2
fi

# Re-prove the one relaxation in the comparison (the BCF ##fileDate
# normalization) on constructed inputs before using it.  Milliseconds, and it
# means the exception is checked on every gate rather than trusted.
echo "=== regress.py self-test ==="
if ! "$PY" "$(dirname "$0")/regress.py" --self-test; then
    echo "error: regress.py self-test failed; the comparison is not trustworthy" >&2
    exit 2
fi
echo

exec "$PY" "$(dirname "$0")/regress.py" "$BASE" "$NEW" --rtol "$RTOL" \
    $([ "$QUIET" -eq 1 ] && echo --quiet) ${EXTRA[@]+"${EXTRA[@]}"}
