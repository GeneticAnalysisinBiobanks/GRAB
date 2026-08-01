# ── Toolchain ─────────────────────────────────────────────────────────────────
# Use ?= so external build environments (the GitHub Actions release workflow,
# manylinux2014 containers, CI-injected sysroot-aware compilers) can override
# the source-compile defaults.
CXX      ?= g++
CC       ?= gcc

# ── Version ───────────────────────────────────────────────────────────────────
# Single source of truth for the GRAB version string.  The file `VERSION` at
# the repo root is read once and the value is injected as the C++ macro
# GRAB_VERSION during compilation of src/.  When releasing, update VERSION,
# commit, and tag (e.g. `git tag v$(cat VERSION)`); the GitHub Actions
# release workflow consumes the same file via this Makefile, so the binary's
# `--version` output, the README, and the release artifact filenames stay
# aligned.
GRAB_VERSION ?= $(strip $(shell cat VERSION 2>/dev/null || echo 0.0.0+unknown))

# ── Platform detection ────────────────────────────────────────────────────────
# UNAME_S is empty on plain Windows cmd; MSYS2/MinGW sets it to MINGW* or MSYS*.
UNAME_S := $(shell uname -s 2>/dev/null)

ifeq ($(filter MINGW% MSYS%,$(UNAME_S)),)
  ifeq ($(filter Darwin%,$(UNAME_S)),)
    # ── Linux ─────────────────────────────────────────────────────────────────
    PLATFORM  := linux
    EXE       :=
    PLATFORM_LDLIBS :=
    # -pipe: pass data between gcc stages via pipes (no temp files).
    PLATFORM_FLAGS  := -pipe
    SHELL     := /bin/bash
  else
    # ── macOS ─────────────────────────────────────────────────────────────────
    PLATFORM  := macos
    EXE       :=
    PLATFORM_LDLIBS :=
    PLATFORM_FLAGS  := -pipe
    SHELL     := /bin/bash
  endif
else
  # ── Windows (MSYS2 / MinGW / Rtools) ─────────────────────────────────────
  PLATFORM  := windows
  EXE       := .exe
  # ws2_32: Winsock2 (htslib hfile.c socket I/O)
  # regex:  POSIX regex (htslib hts_expr.c).  Static-link chain on
  # MSYS2/MinGW:
  #   libregex.a → tre_regcomp / tre_regfree / tre_regerror / tre_regexec
  #   libtre.a   → libintl_gettext (i18n of tre's error strings)
  #   libintl.a  → libiconv_open / libiconv_close / libiconv
  # In dynamic-mode builds the libregex DLL has the entire chain
  # baked in, so the explicit `-ltre -lintl -liconv` block is a
  # no-op; in static-mode builds each .a is required.  Wrapping the
  # block in -Bstatic/-Bdynamic forces these three to be picked up
  # from .a archives regardless of the surrounding global link mode
  # so the resulting .exe is fully self-contained.
  PLATFORM_LDLIBS := -lws2_32 -lregex -Wl,-Bstatic -ltre -lintl -liconv -Wl,-Bdynamic
  PLATFORM_FLAGS  := -pipe
  SHELL     := /usr/bin/bash
  # Redirect linker temp dir to project-local tmp/ (Windows path format).
  WIN_TMP   := $(shell cygpath -w $(CURDIR)/tmp)
  export TMP  := $(WIN_TMP)
  export TEMP := $(WIN_TMP)
endif

# ── Architecture detection ────────────────────────────────────────────────────
# IS_X86 = 1 on x86_64; otherwise 0 (covers arm64/aarch64 macOS M-series, etc.).
UNAME_M := $(shell uname -m 2>/dev/null)
ifeq ($(UNAME_M),x86_64)
  IS_X86 := 1
else
  IS_X86 := 0
endif

# ── AVX2 detection (x86 only) ─────────────────────────────────────────────────
# Apple clang on ARM warns "-mavx2 unused" instead of erroring, so the probe
# would falsely succeed.  Restrict the probe to x86 outright.
ifeq ($(IS_X86),1)
  AVX2_OK := $(shell TMP_AVX=$$(mktemp -u 2>/dev/null || echo __avx2_test__); \
    echo 'int main(){return 0;}' | $(CXX) -x c++ -mavx2 - -o "$$TMP_AVX" 2>/dev/null \
    && { rm -f "$$TMP_AVX"; echo 1; } || rm -f "$$TMP_AVX")
  ifeq ($(AVX2_OK),1)
    SIMD_FLAGS := -mavx2 -mbmi -mbmi2 -mlzcnt -mfma
  else
    SIMD_FLAGS :=
  endif
else
  SIMD_FLAGS :=
endif

# ── Architecture baseline for GRAB sources (runtime SIMD dispatch) ────────────
# GRAB's own SIMD kernels use __attribute__((target(...))) to generate
# AVX-512 / AVX2 / scalar variants and resolve at startup via
# __builtin_cpu_supports().  The baseline march must NOT include AVX2 so
# that the scalar codepaths remain runnable on older x86-64 hardware.
# Third-party code (pgenlib, bgen, htslib, …) keeps compile-time SIMD_FLAGS.
ifeq ($(IS_X86),1)
  # Default: native arch for best SIMD on the build host (AVX-512 on Zen 5,
  # Intel Sapphire Rapids, etc.).  This is the right choice for users who
  # compile GRAB on the same machine that will run it.
  #
  # For portable binary distribution (GitHub Actions release artifacts,
  # docker images, shared cluster software trees) override with:
  #     make GRAB_MARCH=-march=x86-64-v3
  # which pins the baseline ISA at AVX2/FMA/BMI2 (Haswell, 2013+).  GRAB's
  # own hand-written SPAsqr-style kernels keep their AVX-512 variants via
  # __attribute__((target(...))), so the runtime dispatcher in
  # simd_dispatch.hpp still picks the AVX-512 path on capable hosts; only
  # Eigen-expressed and compiler-auto-vectorized code is capped at AVX2.
  GRAB_MARCH ?= -march=native
else
  # ARM (arm64 / aarch64): rely on default tuning + scalar fallbacks.
  GRAB_MARCH ?=
endif

# ── Compiler & linker flags ───────────────────────────────────────────────────
# Naming convention:
#   TP_CXXFLAGS / TP_CFLAGS  : internal defaults for third-party (zlib, zstd,
#                              libdeflate, pgenlib, bgen, htslib).  Include
#                              SIMD_FLAGS so the C/C++ deps get AVX2 codegen.
#   GRAB_CXXFLAGS            : internal defaults for src/*.cpp — no compile-
#                              time SIMD ISA selected; the runtime dispatcher
#                              in simd_dispatch.hpp resolves AVX2/AVX-512
#                              variants via __attribute__((target(...))).
#   LINK_FLAGS               : internal defaults for the final link step.
#
# The standard $(CPPFLAGS) / $(CXXFLAGS) / $(CFLAGS) / $(LDFLAGS) names are
# left for the build environment (external sysroots, hardening flags such
# as -fstack-protector-strong / -D_FORTIFY_SOURCE=2 / -fdebug-prefix-map=…,
# and rpath/link-search additions) to fill in.  They are appended after the
# internal flag sets so external overrides win for any conflicting option.
TP_CXXFLAGS := -std=c++17 -O3 -DNDEBUG $(PLATFORM_FLAGS) \
            -ffunction-sections -fdata-sections \
            -funroll-loops $(SIMD_FLAGS) \
            -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare
GRAB_CXXFLAGS := -std=c++17 -O3 -DNDEBUG $(GRAB_MARCH) $(PLATFORM_FLAGS) \
            -ffunction-sections -fdata-sections \
            -funroll-loops \
            -DGRAB_VERSION=\"$(GRAB_VERSION)\" \
            -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare \
            -Wno-maybe-uninitialized
TP_CFLAGS   := -O3 -DNDEBUG $(PLATFORM_FLAGS) -ffunction-sections -fdata-sections $(SIMD_FLAGS)
# libstdc++/libgcc linkage: dynamic by default (works everywhere; requires the
# system's libstdc++.so.6/libgcc_s.so.1 at runtime — both are part of any
# Linux/Windows g++ install).  Override with STATIC_LIBS="-static-libstdc++
# -static-libgcc" if you have the libstdc++-static archive installed and want
# a binary portable to older systems.  Apple clang on macOS doesn't accept
# these flags, so the override is Linux/Windows-only.
STATIC_LIBS ?=
ifeq ($(PLATFORM),macos)
  LINK_FLAGS := -Wl,-dead_strip -lpthread $(PLATFORM_LDLIBS) $(STATIC_LIBS)
else
  # -lstdc++fs is required by libstdc++ before GCC 9 to resolve
  # std::filesystem symbols. Newer toolchains expose an empty stub library
  # so the flag is harmless there.
  LINK_FLAGS := -Wl,--gc-sections -lpthread -lstdc++fs $(PLATFORM_LDLIBS) $(STATIC_LIBS)
endif

# ── Directories ───────────────────────────────────────────────────────────────
SRC_DIR      := src
BUILD_DIR    := build
ZLIB_DIR     := third_party/zlib-1.3.2
ZSTD_DIR     := third_party/zstd-1.5.7/lib
PGENLIB_DIR  := third_party/plink2-a.6.33
DEFLATE_DIR  := $(PGENLIB_DIR)/libdeflate
BGEN_DIR     := third_party/bgen-1.2.0
HTSLIB_DIR   := third_party/htslib-1.23.1
SQLITE_DIR   := third_party/sqlite-3.53.2

# ── Output binary ─────────────────────────────────────────────────────────────
BIN := $(BUILD_DIR)/grab2$(EXE)

# ── Include paths ─────────────────────────────────────────────────────────────
# Boost.Math and Eigen are header-only — no -l flags required.
INCLUDES := \
    -I$(SRC_DIR) \
    -Ithird_party/eigen-5.0.0 \
    -Ithird_party/boost-1.90.0 \
    -I$(ZLIB_DIR) \
    -I$(ZSTD_DIR) \
    -I$(PGENLIB_DIR)/include \
    -I$(BGEN_DIR)/genfile/include \
    -I$(HTSLIB_DIR) \
    -I$(HTSLIB_DIR)/htscodecs \
    -I$(SQLITE_DIR) \
    -I$(DEFLATE_DIR)

# ── Source discovery ──────────────────────────────────────────────────────────
# Pure-Make recursive wildcard: avoids calling the shell 'find' command,
# which can resolve to Windows FIND.EXE instead of the Unix utility.
rwildcard = $(foreach d,$(wildcard $(addsuffix /*,$(1))),\
                $(call rwildcard,$(d),$(2)) $(filter $(subst *,%,$(2)),$(d)))

# ── GRAB application sources (C++) ────────────────────────────────────────────
SRCS := $(call rwildcard,$(SRC_DIR),*.cpp)
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))
DEPS := $(OBJS:.o=.d)

# ── zlib (C) ──────────────────────────────────────────────────────────────────
ZLIB_SRCS := $(wildcard $(ZLIB_DIR)/*.c)
ZLIB_OBJS := $(patsubst $(ZLIB_DIR)/%.c, $(BUILD_DIR)/zlib/%.o, $(ZLIB_SRCS))

# ── zstd (C) — shared by pgenlib and bgen ─────────────────────────────────────
ZSTD_SRCS := $(wildcard $(ZSTD_DIR)/common/*.c) \
             $(wildcard $(ZSTD_DIR)/compress/*.c) \
             $(wildcard $(ZSTD_DIR)/decompress/*.c)
ZSTD_OBJS := $(patsubst $(ZSTD_DIR)/%.c, $(BUILD_DIR)/zstd/%.o, $(ZSTD_SRCS))

# ── libdeflate (C) — used by pgenlib ──────────────────────────────────────────
DEFLATE_SRCS := $(wildcard $(DEFLATE_DIR)/lib/*.c)
# Per-architecture CPU-feature detection.  Each *_cpu_features.c is guarded
# internally by arch #ifdefs, so compiling both on any target is harmless (the
# non-matching one becomes an empty translation unit); the matching one supplies
# the runtime-dispatch symbols (libdeflate_{x86,arm}_cpu_features) referenced by
# adler32.c / crc32.c.  Omitting the ARM file breaks the aarch64 link.
ifneq ($(wildcard $(DEFLATE_DIR)/lib/x86/x86_cpu_features.c),)
  DEFLATE_SRCS += $(DEFLATE_DIR)/lib/x86/x86_cpu_features.c
endif
ifneq ($(wildcard $(DEFLATE_DIR)/lib/arm/arm_cpu_features.c),)
  DEFLATE_SRCS += $(DEFLATE_DIR)/lib/arm/arm_cpu_features.c
endif
DEFLATE_OBJS := $(patsubst $(DEFLATE_DIR)/lib/%.c, $(BUILD_DIR)/libdeflate/%.o, $(DEFLATE_SRCS))

# ── pgenlib (C++ .cc files + SFMT.c) ─────────────────────────────────────────
PGEN_CC_SRCS := $(wildcard $(PGENLIB_DIR)/include/*.cc)
PGEN_CC_OBJS := $(patsubst $(PGENLIB_DIR)/include/%.cc, $(BUILD_DIR)/pgenlib/%.o, $(PGEN_CC_SRCS))
PGEN_C_SRCS  := $(PGENLIB_DIR)/include/SFMT.c
PGEN_C_OBJS  := $(BUILD_DIR)/pgenlib/SFMT.o
PGEN_OBJS    := $(PGEN_CC_OBJS) $(PGEN_C_OBJS)

# ── bgen (C++) — exclude View.cpp (needs boost::filesystem + IndexQuery) ─────
BGEN_SRCS := $(BGEN_DIR)/src/bgen.cpp \
             $(BGEN_DIR)/src/zlib.cpp \
             $(BGEN_DIR)/src/MissingValue.cpp
BGEN_OBJS := $(patsubst $(BGEN_DIR)/src/%.cpp, $(BUILD_DIR)/bgen/%.o, $(BGEN_SRCS))

# ── sqlite (C) — read-only .bgi (bgenix) index access for BGEN ─────────────────
SQLITE_SRCS := $(SQLITE_DIR)/sqlite3.c
SQLITE_OBJS := $(BUILD_DIR)/sqlite/sqlite3.o

# ── htslib (C) ────────────────────────────────────────────────────────────────
HTSLIB_SRCS    := $(wildcard $(HTSLIB_DIR)/*.c)
HTSCODEC_SRCS  := $(wildcard $(HTSLIB_DIR)/htscodecs/htscodecs/*.c)
# os/rand.c is #included by hts_os.c on Windows — not compiled separately.
HTSLIB_OBJS    := $(patsubst $(HTSLIB_DIR)/%.c, $(BUILD_DIR)/htslib/%.o, $(HTSLIB_SRCS))
HTSCODEC_OBJS  := $(patsubst $(HTSLIB_DIR)/htscodecs/htscodecs/%.c, \
                     $(BUILD_DIR)/htslib/htscodecs/%.o, $(HTSCODEC_SRCS))

# ── All objects ───────────────────────────────────────────────────────────────
ALL_OBJS := $(OBJS) $(ZLIB_OBJS) $(ZSTD_OBJS) $(DEFLATE_OBJS) \
            $(PGEN_OBJS) $(BGEN_OBJS) \
            $(HTSLIB_OBJS) $(HTSCODEC_OBJS) $(SQLITE_OBJS)

# ── Default target ────────────────────────────────────────────────────────────
.PHONY: all clean run

all: $(BIN)

# ── Link ──────────────────────────────────────────────────────────────────────
$(BIN): $(ALL_OBJS) | $(BUILD_DIR) tmp
	$(CXX) $(TP_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $^ -o $@ $(LINK_FLAGS) $(LDLIBS)

# ── Directory creation ────────────────────────────────────────────────────────
$(BUILD_DIR) tmp:
	mkdir -p $@

# ── Compile GRAB (C++) ────────────────────────────────────────────────────────
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | tmp
	@mkdir -p $(@D)
	$(CXX) $(GRAB_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -MMD -MP $(INCLUDES) -c $< -o $@

-include $(DEPS)

# ── Compile zlib (C) ──────────────────────────────────────────────────────────
$(BUILD_DIR)/zlib/%.o: $(ZLIB_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(TP_CFLAGS) $(CPPFLAGS) $(CFLAGS) -I$(ZLIB_DIR) -c $< -o $@

# ── Compile zstd (C) ──────────────────────────────────────────────────────────
$(BUILD_DIR)/zstd/%.o: $(ZSTD_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(TP_CFLAGS) $(CPPFLAGS) $(CFLAGS) -DZSTD_DISABLE_ASM -I$(ZSTD_DIR) -I$(ZSTD_DIR)/common -c $< -o $@

# ── Compile libdeflate (C) ────────────────────────────────────────────────────
$(BUILD_DIR)/libdeflate/%.o: $(DEFLATE_DIR)/lib/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(TP_CFLAGS) $(CPPFLAGS) $(CFLAGS) -I$(DEFLATE_DIR) -c $< -o $@

# ── Compile pgenlib (.cc → C++) ───────────────────────────────────────────────
#  IGNORE_BUNDLED_ZSTD: use our vendored zstd via -I instead of ../zstd/ path.
PGEN_CXXFLAGS := $(TP_CXXFLAGS) -DIGNORE_BUNDLED_ZSTD \
    -Wno-sign-compare -Wno-unused-function -Wno-missing-field-initializers \
    -Wno-maybe-uninitialized

$(BUILD_DIR)/pgenlib/%.o: $(PGENLIB_DIR)/include/%.cc | tmp
	@mkdir -p $(@D)
	$(CXX) $(PGEN_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -I$(PGENLIB_DIR)/include -I$(ZSTD_DIR) \
	    -I$(DEFLATE_DIR) -c $< -o $@

# SFMT.c is plain C
$(BUILD_DIR)/pgenlib/SFMT.o: $(PGENLIB_DIR)/include/SFMT.c | tmp
	@mkdir -p $(@D)
	$(CC) $(TP_CFLAGS) $(CPPFLAGS) $(CFLAGS) -I$(PGENLIB_DIR)/include -c $< -o $@

# ── Compile bgen (C++) ────────────────────────────────────────────────────────
BGEN_CXXFLAGS := $(TP_CXXFLAGS) -Wno-sign-compare -Wno-unused-variable

$(BUILD_DIR)/bgen/%.o: $(BGEN_DIR)/src/%.cpp | tmp
	@mkdir -p $(@D)
	$(CXX) $(BGEN_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) -I$(BGEN_DIR)/genfile/include \
	    -I$(ZSTD_DIR) -I$(ZLIB_DIR) -c $< -o $@

# ── Compile sqlite (C) — read-only .bgi index; lean, no extension loading ─────
# SQLITE_OMIT_LOAD_EXTENSION avoids a libdl dependency; THREADSAFE=1 keeps the
# one-shot index read (main thread, BgenData construction) safe in the otherwise
# multi-threaded process.  No -march=native (uses third-party SIMD_FLAGS only).
# -Wno-stringop-overread silences a GCC false positive in sqlite3ColumnSetColl
# (the inlined sqlite3Strlen30 -> strlen call): GCC cannot prove the collation-
# name pointer is non-null at that site and warns about a "region of size 0".
SQLITE_CFLAGS := $(TP_CFLAGS) \
    -DSQLITE_OMIT_LOAD_EXTENSION -DSQLITE_THREADSAFE=1 \
    -DSQLITE_DEFAULT_MEMSTATUS=0 -DSQLITE_OMIT_DEPRECATED \
    -Wno-stringop-overread

$(BUILD_DIR)/sqlite/%.o: $(SQLITE_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(SQLITE_CFLAGS) $(CPPFLAGS) $(CFLAGS) -I$(SQLITE_DIR) -c $< -o $@

# ── Compile htslib (C) ────────────────────────────────────────────────────────
HTSLIB_CFLAGS := $(TP_CFLAGS) -DHAVE_CONFIG_H \
    -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htscodecs -I$(ZLIB_DIR) \
    -Wno-sign-compare

$(BUILD_DIR)/htslib/%.o: $(HTSLIB_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(HTSLIB_CFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(BUILD_DIR)/htslib/htscodecs/%.o: $(HTSLIB_DIR)/htscodecs/htscodecs/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(HTSLIB_CFLAGS) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

# ── Install ───────────────────────────────────────────────────────────────────
# Standard PREFIX / DESTDIR semantics so the GitHub Actions release workflow
# and ad-hoc local installations (`make install PREFIX=$HOME/.local`) can
# stage the binary into a chroot via DESTDIR while keeping the runtime path
# under PREFIX.  Both `install -d` and `install -m` are supported by GNU
# coreutils and BSD `install` (macOS), avoiding the GNU-only `-D` shortcut.
PREFIX  ?= /usr/local
DESTDIR ?=

.PHONY: install
install: $(BIN)
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 0755 $(BIN) $(DESTDIR)$(PREFIX)/bin/grab2$(EXE)

# ── Dist (release packaging) ──────────────────────────────────────────────────
# Produces dist/grab2-$(GRAB_VERSION)-$(PLATFORM)-$(UNAME_M).tar.gz, a
# self-contained archive holding the stripped binary alongside README.md,
# LICENSE, and VERSION.  Intended for ad-hoc local testing of the same
# packaging recipe that .github/workflows/release.yml runs in CI; the
# GitHub Actions workflow performs the equivalent steps explicitly so
# it can handle the Windows .zip vs Unix .tar.gz difference.
#
# `strip` removes debug / symbol-table sections to shrink the binary by
# roughly 3-4x.  STRIP is overridable so cross-compile / sysroot setups
# can substitute the target-matched tool (e.g. aarch64-linux-gnu-strip).
STRIP   ?= strip
DIST_DIR := dist
DIST_PKG := grab2-$(GRAB_VERSION)-$(PLATFORM)-$(UNAME_M)

.PHONY: dist
dist: $(BIN)
	@rm -rf $(DIST_DIR)/$(DIST_PKG) $(DIST_DIR)/$(DIST_PKG).tar.gz
	@mkdir -p $(DIST_DIR)/$(DIST_PKG)
	cp $(BIN) $(DIST_DIR)/$(DIST_PKG)/
	$(STRIP) $(DIST_DIR)/$(DIST_PKG)/grab2$(EXE) || true
	cp README.md LICENSE VERSION $(DIST_DIR)/$(DIST_PKG)/
	cd $(DIST_DIR) && tar czf $(DIST_PKG).tar.gz $(DIST_PKG)
	@echo "Created $(DIST_DIR)/$(DIST_PKG).tar.gz"

# ── Tests and benchmarks ──────────────────────────────────────────────────────
# Developer-only targets.  `tests/` is deliberately outside the rwildcard over
# $(SRC_DIR), so nothing here is linked into the shipped binary and the
# "binary is the deliverable" property is preserved.
#
# Each test/benchmark is a self-contained translation unit with its own main().
# Those that exercise only header-level code (the SPA kernels, the reference
# CGF) link against nothing; those that need compiled GRAB objects list them
# explicitly in TEST_<name>_OBJS below.  There is no test-runner library and no
# third-party framework: tests/tinytest.hpp is ~150 lines of header.
#
#   make test    build and run every test, stop at the first failing binary
#   make bench   build and run every benchmark
#
# Tests are compiled with the same GRAB_CXXFLAGS as src/ so that the SIMD
# dispatch, -march setting and optimization level under test match the shipped
# configuration.  Assertions are wanted, so -DNDEBUG is filtered out.
TEST_DIR       := tests
TEST_BUILD_DIR := $(BUILD_DIR)/tests

TEST_SRCS  := $(wildcard $(TEST_DIR)/*_test.cpp)
TEST_BINS  := $(patsubst $(TEST_DIR)/%.cpp, $(TEST_BUILD_DIR)/%$(EXE), $(TEST_SRCS))

BENCH_SRCS := $(wildcard $(TEST_DIR)/bench_*.cpp)
BENCH_BINS := $(patsubst $(TEST_DIR)/%.cpp, $(TEST_BUILD_DIR)/%$(EXE), $(BENCH_SRCS))

TEST_CXXFLAGS := $(filter-out -DNDEBUG,$(GRAB_CXXFLAGS))

# Per-test object requirements.  A test that exercises only header-level code
# declares nothing.  The lists below replace the hand-written g++ command lines
# previously carried in each test's header comment, so `make test` is now the
# single source of truth for how a test links.
# The SPA CGF kernels live in a .cpp rather than a header: the
# __attribute__((target(...))) SIMD variants and the static const dispatch
# pointers must have exactly one definition.  Both the test and the benchmark
# therefore link the object rather than including the implementation.
TESTOBJS_spa_cgf_test  := $(BUILD_DIR)/util/spa_cgf.o
TESTOBJS_bench_spa_cgf := $(BUILD_DIR)/util/spa_cgf.o
# bench_spa_tail times the tail path against the linear tail log10p_unify
# Stage 3 deleted; spa.hpp is header-only but reaches math_helper's pnorm and
# pnormLog, which are out of line.
TESTOBJS_bench_spa_tail := $(BUILD_DIR)/util/math_helper.o
# bench_hwe times the deleted linear HWE test against plink2's HweLnP, which
# reaches it through the geno_factory wrapper.
TESTOBJS_bench_hwe := \
    $(BUILD_DIR)/geno_factory/hwe.o \
    $(BUILD_DIR)/pgenlib/plink2_stats.o \
    $(BUILD_DIR)/pgenlib/plink2_base.o \
    $(BUILD_DIR)/pgenlib/SFMT.o
# spagrm_cgf.hpp is header-only but delegates its class-1 term to the dispatched
# binomial kernels, so it needs the same object.
TESTOBJS_spagrm_cgf_test := $(BUILD_DIR)/util/spa_cgf.o
# spamix_cgf.hpp is likewise header-only over the dispatched binomial kernels.
TESTOBJS_spamix_cgf_test := $(BUILD_DIR)/util/spa_cgf.o
# spamixlocalp_cgf.hpp is header-only over spa_cgf's hapcount variant.
TESTOBJS_spamixlocalp_cgf_test := $(BUILD_DIR)/util/spa_cgf.o
# wtcoxg_cgf.hpp is header-only over the same kernels; its bivariate-normal
# tests additionally need math_helper.o for pmvnorm2dHalfRect / bvnCdf.
TESTOBJS_wtcoxg_cgf_test := \
    $(BUILD_DIR)/util/spa_cgf.o \
    $(BUILD_DIR)/util/math_helper.o
# log10p_test covers the log-domain distribution tier of the log10p_unify
# project (zFromNegLog10P, chisq1FromNegLog10P, cauchyCombineLog10, ptLog,
# pmvnorm2dHalfRectLog), all of which live in math_helper.cpp.  The object is
# listed from Stage 0, when the suite is still a skeleton, so that Stage 1 adds
# assertions to a suite that already links rather than changing both at once.
# Stage 9 additionally cross-checks LOG10P_HWE, so the suite links the HWE
# wrapper and plink2's statistics object that supplies HweLnP behind it.
TESTOBJS_log10p_test := \
    $(BUILD_DIR)/util/math_helper.o \
    $(BUILD_DIR)/geno_factory/hwe.o \
    $(BUILD_DIR)/pgenlib/plink2_stats.o \
    $(BUILD_DIR)/pgenlib/plink2_base.o \
    $(BUILD_DIR)/pgenlib/SFMT.o

# spagrm_ibd_test exercises the two relatedness inputs SPAGRM does not build
# in the same run — the pairwise-IBD table and the sparse GRM — so it links
# the real loader and the real Chow-Liu builder rather than a copy of either.
# grm_null.o pulls in SPAGRMClass (spagrm.o, and through it the saddlepoint
# kernels) and the compressed-text reader (text_stream.o and the three
# compression libraries).
TESTOBJS_spagrm_ibd_test := \
    $(BUILD_DIR)/spagrm/grm_null.o \
    $(BUILD_DIR)/spagrm/spagrm.o \
    $(BUILD_DIR)/io/sparse_grm.o \
    $(BUILD_DIR)/util/spa_cgf.o \
    $(BUILD_DIR)/util/math_helper.o \
    $(BUILD_DIR)/util/text_stream.o \
    $(ZSTD_OBJS) $(ZLIB_OBJS) $(DEFLATE_OBJS)

TESTOBJS_lanc_simd_test := \
    $(BUILD_DIR)/localplus/lanc_io.o \
    $(BUILD_DIR)/geno_factory/variant_filter.o \
    $(ZSTD_OBJS)

TESTOBJS_lanc_roundtrip_test := $(TESTOBJS_lanc_simd_test)

TESTOBJS_lanc_convert_rfmix_smoke_test := \
    $(BUILD_DIR)/localplus/lanc_io.o \
    $(BUILD_DIR)/localplus/lanc_convert_rfmix.o \
    $(BUILD_DIR)/geno_factory/variant_filter.o \
    $(BUILD_DIR)/io/subject_filter.o \
    $(HTSLIB_OBJS) $(HTSCODEC_OBJS) \
    $(ZSTD_OBJS) $(ZLIB_OBJS) $(DEFLATE_OBJS)

# Second expansion lets the pattern rule name a per-target variable derived
# from the stem, so each test pulls in exactly the objects it needs.
.SECONDEXPANSION:

$(TEST_BUILD_DIR)/%$(EXE): $(TEST_DIR)/%.cpp $$(TESTOBJS_$$*) | tmp
	@mkdir -p $(@D)
	$(CXX) $(TEST_CXXFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -MMD -MP \
	    -I$(TEST_DIR) $(INCLUDES) $< $(TESTOBJS_$*) -o $@ $(LINK_FLAGS) $(LDLIBS)

# Header dependencies for tests and benchmarks.  Without these, editing a
# header under tests/ (spa_reference.hpp, tinytest.hpp) would leave a stale
# binary in place and a test run would silently exercise the previous code.
TEST_DEPS := $(TEST_BINS:$(EXE)=.d) $(BENCH_BINS:$(EXE)=.d)
-include $(TEST_DEPS)

.PHONY: test bench
test: $(TEST_BINS)
	@fail=0; \
	for t in $(TEST_BINS); do \
	    echo "── $$t ─────────────────────────────────────────────"; \
	    "$$t" || fail=1; \
	    echo; \
	done; \
	if [ $$fail -ne 0 ]; then echo "TEST SUITE FAILED"; exit 1; fi; \
	echo "ALL TESTS PASSED"

bench: $(BENCH_BINS)
	@for b in $(BENCH_BINS); do \
	    echo "── $$b ─────────────────────────────────────────────"; \
	    "$$b"; \
	    echo; \
	done

# ── Helpers ───────────────────────────────────────────────────────────────────
.PHONY: run
run: all
	$(BIN)

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(BUILD_DIR) $(DIST_DIR)
