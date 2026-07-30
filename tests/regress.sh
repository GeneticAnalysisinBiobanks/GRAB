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
#
# Non-numeric columns (IDs, alleles, chromosome labels) are compared exactly;
# any difference there is a structural change and is reported as such.
#
# Usage:
#     tests/regress.sh BASE_DIR NEW_DIR [--rtol 1e-9] [--quiet]
#
# Exit status is 0 when every difference is within tolerance and no structural
# change was found, 1 otherwise.  The intended workflow between stages is:
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

while [ $# -gt 0 ]; do
    case "$1" in
        --rtol)  RTOL="$2"; shift 2 ;;
        --quiet) QUIET=1; shift ;;
        -h|--help)
            sed -n '2,32p' "$0"; exit 0 ;;
        *)
            if   [ -z "$BASE" ]; then BASE="$1"
            elif [ -z "$NEW"  ]; then NEW="$1"
            else echo "unexpected argument: $1" >&2; exit 2
            fi
            shift ;;
    esac
done

if [ -z "$BASE" ] || [ -z "$NEW" ]; then
    echo "usage: $0 BASE_DIR NEW_DIR [--rtol RTOL] [--quiet]" >&2
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

exec "$PY" "$(dirname "$0")/regress.py" "$BASE" "$NEW" --rtol "$RTOL" \
    $([ "$QUIET" -eq 1 ] && echo --quiet)
