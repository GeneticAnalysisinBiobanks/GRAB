# ── Toolchain ─────────────────────────────────────────────────────────────────
CXX      := g++
CC       := gcc

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
  # regex:  POSIX regex  (htslib hts_expr.c)
  PLATFORM_LDLIBS := -lws2_32 -lregex
  PLATFORM_FLAGS  := -pipe
  SHELL     := /usr/bin/bash
  # Redirect linker temp dir to project-local tmp/ (Windows path format).
  WIN_TMP   := $(shell cygpath -w $(CURDIR)/tmp)
  export TMP  := $(WIN_TMP)
  export TEMP := $(WIN_TMP)
endif

# ── AVX2 detection ────────────────────────────────────────────────────────────
# Try to compile a small AVX2 test; if the compiler supports it, enable -mavx2.
# Use a temp file as output target for portability (MSYS2/Windows lacks /dev/null for ld).
AVX2_OK := $(shell TMP_AVX=$$(mktemp -u 2>/dev/null || echo __avx2_test__); \
  echo 'int main(){return 0;}' | $(CXX) -x c++ -mavx2 - -o "$$TMP_AVX" 2>/dev/null \
  && { rm -f "$$TMP_AVX"; echo 1; } || rm -f "$$TMP_AVX")
ifeq ($(AVX2_OK),1)
  SIMD_FLAGS := -mavx2 -mbmi -mbmi2 -mlzcnt -mfma
else
  SIMD_FLAGS :=
endif

# ── Architecture baseline for GRAB sources (runtime SIMD dispatch) ────────────
# GRAB's own SIMD kernels use __attribute__((target(...))) to generate
# AVX-512 / AVX2 / scalar variants and resolve at startup via
# __builtin_cpu_supports().  The baseline march must NOT include AVX2 so
# that the scalar codepaths remain runnable on older x86-64 hardware.
# Third-party code (pgenlib, bgen, htslib, …) keeps compile-time SIMD_FLAGS.
UNAME_M := $(shell uname -m 2>/dev/null)
ifeq ($(UNAME_M),x86_64)
  # x86-64-v2 = SSE4.2, SSSE3, POPCNT, CX16 (Nehalem 2008+)
  GRAB_MARCH := -march=x86-64-v2
else
  GRAB_MARCH :=
endif

# ── Compiler & linker flags ───────────────────────────────────────────────────
# CXXFLAGS: used for third-party C++ (pgenlib, bgen) — includes SIMD_FLAGS.
CXXFLAGS := -std=c++17 -O3 -DNDEBUG $(PLATFORM_FLAGS) \
            -ffunction-sections -fdata-sections \
            -funroll-loops $(SIMD_FLAGS) \
            -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare
# GRAB_CXXFLAGS: used for src/*.cpp — no compile-time SIMD; uses target attrs.
GRAB_CXXFLAGS := -std=c++17 -O3 -DNDEBUG $(GRAB_MARCH) $(PLATFORM_FLAGS) \
            -ffunction-sections -fdata-sections \
            -funroll-loops \
            -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare
CFLAGS   := -O3 -DNDEBUG $(PLATFORM_FLAGS) -ffunction-sections -fdata-sections $(SIMD_FLAGS)
STATIC_LIBS := -static-libstdc++ -static-libgcc
ifeq ($(PLATFORM),macos)
  LDFLAGS := -Wl,-dead_strip -lpthread $(PLATFORM_LDLIBS) $(STATIC_LIBS)
else
  LDFLAGS := -Wl,--gc-sections -lpthread $(PLATFORM_LDLIBS) $(STATIC_LIBS)
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
BIN := $(BUILD_DIR)/grab$(EXE)

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
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# ── Directory creation ────────────────────────────────────────────────────────
$(BUILD_DIR) tmp:
	mkdir -p $@

# ── Compile GRAB (C++) ────────────────────────────────────────────────────────
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | tmp
	@mkdir -p $(@D)
	$(CXX) $(GRAB_CXXFLAGS) -MMD -MP $(INCLUDES) -c $< -o $@

-include $(DEPS)

# ── Compile zlib (C) ──────────────────────────────────────────────────────────
$(BUILD_DIR)/zlib/%.o: $(ZLIB_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -I$(ZLIB_DIR) -c $< -o $@

# ── Compile zstd (C) ──────────────────────────────────────────────────────────
$(BUILD_DIR)/zstd/%.o: $(ZSTD_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -DZSTD_DISABLE_ASM -I$(ZSTD_DIR) -I$(ZSTD_DIR)/common -c $< -o $@

# ── Compile libdeflate (C) ────────────────────────────────────────────────────
$(BUILD_DIR)/libdeflate/%.o: $(DEFLATE_DIR)/lib/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -I$(DEFLATE_DIR) -c $< -o $@

# ── Compile pgenlib (.cc → C++) ───────────────────────────────────────────────
#  IGNORE_BUNDLED_ZSTD: use our vendored zstd via -I instead of ../zstd/ path.
PGEN_CXXFLAGS := $(CXXFLAGS) -DIGNORE_BUNDLED_ZSTD \
    -Wno-sign-compare -Wno-unused-function -Wno-missing-field-initializers \
    -Wno-maybe-uninitialized

$(BUILD_DIR)/pgenlib/%.o: $(PGENLIB_DIR)/include/%.cc | tmp
	@mkdir -p $(@D)
	$(CXX) $(PGEN_CXXFLAGS) -I$(PGENLIB_DIR)/include -I$(ZSTD_DIR) \
	    -I$(DEFLATE_DIR) -c $< -o $@

# SFMT.c is plain C
$(BUILD_DIR)/pgenlib/SFMT.o: $(PGENLIB_DIR)/include/SFMT.c | tmp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -I$(PGENLIB_DIR)/include -c $< -o $@

# ── Compile bgen (C++) ────────────────────────────────────────────────────────
BGEN_CXXFLAGS := $(CXXFLAGS) -Wno-sign-compare -Wno-unused-variable

$(BUILD_DIR)/bgen/%.o: $(BGEN_DIR)/src/%.cpp | tmp
	@mkdir -p $(@D)
	$(CXX) $(BGEN_CXXFLAGS) -I$(BGEN_DIR)/genfile/include \
	    -I$(ZSTD_DIR) -I$(ZLIB_DIR) -c $< -o $@

# ── Compile htslib (C) ────────────────────────────────────────────────────────
HTSLIB_CFLAGS := $(CFLAGS) -DHAVE_CONFIG_H \
    -I$(HTSLIB_DIR) -I$(HTSLIB_DIR)/htscodecs -I$(ZLIB_DIR) \
    -Wno-sign-compare

$(BUILD_DIR)/htslib/%.o: $(HTSLIB_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(HTSLIB_CFLAGS) -c $< -o $@

$(BUILD_DIR)/htslib/htscodecs/%.o: $(HTSLIB_DIR)/htscodecs/htscodecs/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(HTSLIB_CFLAGS) -c $< -o $@

# ── Helpers ───────────────────────────────────────────────────────────────────
.PHONY: run
run: all
	$(BIN)

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(BUILD_DIR)
