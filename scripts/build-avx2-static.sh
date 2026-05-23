#!/usr/bin/env bash
# build-avx2-static.sh
#
# Build a portable GRAB binary that:
#   1. Caps the instruction set at AVX2 (no AVX-512), so it runs on any
#      x86-64-v3 CPU — Intel Haswell (2013) / AMD Excavator (2015) and later.
#   2. Statically links libstdc++ and libgcc, so the binary does not depend
#      on the target host's libstdc++.so.6 / libgcc_s.so.1 versions.
#
# Output: build/grab2
#
# Linux/Windows g++ only.  Apple clang on macOS does not accept
# -static-libstdc++, so this script refuses to run there.

set -euo pipefail

# Resolve to project root regardless of where the script is invoked from.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"

OS="$(uname -s)"
case "${OS}" in
    Linux|MINGW*|MSYS*) ;;
    Darwin)
        echo "Error: macOS / Apple clang does not support -static-libstdc++." >&2
        echo "       Use the plain 'make' target on macOS." >&2
        exit 1
        ;;
    *)
        echo "Warning: untested platform '${OS}'; proceeding anyway." >&2
        ;;
esac

# Verify that a static libstdc++ archive is available; otherwise the link
# step silently falls back to dynamic linking and the output is not portable.
GXX_BIN="${CXX:-g++}"
LIBSTDCXX_A="$("${GXX_BIN}" -print-file-name=libstdc++.a 2>/dev/null || true)"
if [[ "${LIBSTDCXX_A}" == "libstdc++.a" || ! -f "${LIBSTDCXX_A}" ]]; then
    echo "Error: libstdc++.a not found via '${GXX_BIN} -print-file-name=libstdc++.a'." >&2
    echo "       Install the static libstdc++ archive:" >&2
    echo "         Debian/Ubuntu : apt-get install libstdc++-<ver>-dev (static archive included)" >&2
    echo "         Rocky/RHEL    : dnf install libstdc++-static" >&2
    echo "         Arch          : already shipped with gcc" >&2
    exit 1
fi

JOBS="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"

echo "==> Cleaning previous build artifacts"
make clean >/dev/null

echo "==> Building with:"
echo "     GRAB_MARCH = -march=x86-64-v3      (AVX2 baseline, no AVX-512)"
echo "     STATIC_LIBS = -static-libstdc++ -static-libgcc"
echo "     parallel jobs = ${JOBS}"

make -j"${JOBS}" \
    GRAB_MARCH="-march=x86-64-v3" \
    STATIC_LIBS="-static-libstdc++ -static-libgcc"

BIN="build/grab2"
if [[ ! -x "${BIN}" ]]; then
    echo "Error: build/grab2 was not produced." >&2
    exit 1
fi

echo
echo "==> Built ${BIN}"
ls -lh "${BIN}"

# Sanity check 1: ensure no AVX-512 instructions made it into the binary.
# AVX-512 uses the EVEX prefix; objdump -d -M intel emits a 'zmm*' register
# operand for every AVX-512 vector instruction.  Grep counts any reference.
if command -v objdump >/dev/null 2>&1; then
    AVX512_HITS="$(objdump -d -M intel "${BIN}" 2>/dev/null | grep -cE '\bzmm[0-9]+\b' || true)"
    echo "==> AVX-512 instruction references in disassembly: ${AVX512_HITS}"
    if [[ "${AVX512_HITS}" -gt 0 ]]; then
        echo "    (kernels guarded by __attribute__((target(\"avx2,avx512f,...\")))"
        echo "     are skipped at runtime via simdLevel(); the binary still runs"
        echo "     on AVX2-only hardware.)"
    fi
fi

# Sanity check 2: confirm static linkage of libstdc++ / libgcc on Linux.
if [[ "${OS}" == "Linux" ]] && command -v ldd >/dev/null 2>&1; then
    echo "==> Shared library dependencies:"
    ldd "${BIN}" || true
    if ldd "${BIN}" 2>/dev/null | grep -qE 'libstdc\+\+|libgcc_s'; then
        echo "Warning: libstdc++ or libgcc_s appears in ldd output — static" >&2
        echo "         linkage may have silently failed." >&2
        exit 1
    fi
fi

echo
echo "==> Done."
