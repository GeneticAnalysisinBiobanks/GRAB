# ── Toolchain ─────────────────────────────────────────────────────────────────
# -pipe: pass data between gcc compilation stages via pipes (no temp files).
# TMP/TEMP: redirect the linker's temp dir to project-local tmp/ using the
# Windows path format that MinGW ld expects (cygpath -w converts POSIX→Win32).
SHELL    := /usr/bin/bash
WIN_TMP  := $(shell cygpath -w $(CURDIR)/tmp)
export TMP  := $(WIN_TMP)
export TEMP := $(WIN_TMP)
CXX      := g++
CC       := gcc
CXXFLAGS := -std=c++17 -O3 -march=native -DNDEBUG -pipe \
            -ffunction-sections -fdata-sections \
            -funroll-loops \
            -Wall -Wextra -Wno-unused-parameter
CFLAGS   := -O3 -DNDEBUG -pipe -ffunction-sections -fdata-sections
LDFLAGS  := -Wl,--gc-sections

# ── Directories ───────────────────────────────────────────────────────────────
SRC_DIR   := src
BUILD_DIR := build
ZLIB_DIR  := third_party/zlib-1.3.2

# ── Output binary ─────────────────────────────────────────────────────────────
BIN := $(BUILD_DIR)/grab.exe

# ── Include paths ─────────────────────────────────────────────────────────────
# Boost.Math and Eigen are header-only — no -l flags required.
# zlib is compiled from source (C files in third_party/zlib-1.3.2).
INCLUDES := \
    -I$(SRC_DIR) \
    -Ithird_party/eigen-5.0.0 \
    -Ithird_party/boost-1.90.0 \
    -I$(ZLIB_DIR)

# ── Source discovery ──────────────────────────────────────────────────────────
# Pure-Make recursive wildcard: avoids calling the shell 'find' command,
# which can resolve to Windows FIND.EXE instead of the Unix utility.
rwildcard = $(foreach d,$(wildcard $(addsuffix /*,$(1))),\
                $(call rwildcard,$(d),$(2)) $(filter $(subst *,%,$(2)),$(d)))

SRCS := $(call rwildcard,$(SRC_DIR),*.cpp)

# zlib C sources (compiled separately with CC, not CXX).
ZLIB_SRCS := $(wildcard $(ZLIB_DIR)/*.c)
ZLIB_OBJS := $(patsubst $(ZLIB_DIR)/%.c, $(BUILD_DIR)/zlib/%.o, $(ZLIB_SRCS))

# Preserve the directory structure inside build/ so parallel builds don't clash.
OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(SRCS))

# Auto-generated dependency files (one per .o) for incremental header tracking.
DEPS := $(OBJS:.o=.d)

# ── Default target ────────────────────────────────────────────────────────────
.PHONY: all clean run

all: $(BIN)

# ── Link ──────────────────────────────────────────────────────────────────────
$(BIN): $(OBJS) $(ZLIB_OBJS) | $(BUILD_DIR) tmp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

# ── Directory creation ────────────────────────────────────────────────────────
# Order-only prerequisite: created once, never triggers re-linking.
$(BUILD_DIR) tmp:
	mkdir -p $@

# ── Compile ───────────────────────────────────────────────────────────────────
# -MMD -MP writes dependency rules into build/<rel>.d alongside each .o file,
# so edits to any included header trigger exactly the right recompilations.
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | tmp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -MMD -MP $(INCLUDES) -c $< -o $@

# Pull in the generated dependency rules (the leading dash suppresses errors
# on the first build when no .d files exist yet).
-include $(DEPS)

# ── Compile zlib (C) ──────────────────────────────────────────────────────────
$(BUILD_DIR)/zlib/%.o: $(ZLIB_DIR)/%.c | tmp
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -I$(ZLIB_DIR) -c $< -o $@

# ── Helpers ───────────────────────────────────────────────────────────────────
.PHONY: run
run: all
	$(BIN)

# ── Clean ─────────────────────────────────────────────────────────────────────
clean:
	rm -rf $(BUILD_DIR)
