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
  # regex:  POSIX regex (htslib hts_expr.c).  MSYS2's libregex.a
  # references tre_regcomp/tre_regexec/etc., which are provided by
  # the separate libtre package.  In dynamic links the libregex DLL
  # has tre baked in, so the explicit -ltre is a no-op; in static
  # links it is required to satisfy the unresolved tre_* symbols.
  # Wrapping -ltre in -Bstatic/-Bdynamic forces it to be picked up
  # from libtre.a regardless of the surrounding global link mode,
  # which keeps the resulting .exe self-contained.
  PLATFORM_LDLIBS := -lws2_32 -lregex -Wl,-Bstatic -ltre -Wl,-Bdynamic
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
    -I$(HTSLIB_DIR)/htscodecs

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
# x86 CPU feature detection (only on x86/x86_64)
ifneq ($(wildcard $(DEFLATE_DIR)/lib/x86/x86_cpu_features.c),)
  DEFLATE_SRCS += $(DEFLATE_DIR)/lib/x86/x86_cpu_features.c
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
            $(HTSLIB_OBJS) $(HTSCODEC_OBJS)

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

# ── Helpers ───────────────────────────────────────────────────────────────────
.PHONY: run
run: all
	$(BIN)

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(BUILD_DIR) $(DIST_DIR)
