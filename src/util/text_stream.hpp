// text_stream.hpp — Unified compressed/plain text I/O
//
// TextWriter / TextReader auto-detect compression from file extension:
//   .gz  → gzip  (zlib)
//   .zst → zstd
//   else → plain text (FILE*)
//
// Both are move-only, RAII wrappers.  The write/read APIs mirror what
// marker.cpp and the utility modes need: bulk string writes and
// line-by-line reads.
#pragma once

#include <cstddef>
#include <cstdio>
#include <string>
#include <string_view>

// ══════════════════════════════════════════════════════════════════════
// TextWriter
// ══════════════════════════════════════════════════════════════════════

class TextWriter {
  public:
    enum class Mode { Plain, Gzip, Zstd };

    /// Open `path` for writing.  Compression inferred from extension.
    explicit TextWriter(const std::string &path);

    /// Open `path` with explicit mode and compression level.
    /// Level 0 means library default.
    TextWriter(
        const std::string &path,
        Mode mode,
        int level
    );

    ~TextWriter();

    TextWriter(const TextWriter &) = delete;
    TextWriter &operator=(const TextWriter &) = delete;

    TextWriter(TextWriter &&o) noexcept;
    TextWriter &operator=(TextWriter &&o) noexcept;

    /// Write raw bytes.
    void write(
        const char *data,
        size_t len
    );

    /// Convenience: write a string_view (also handles std::string and const char*).
    void write(std::string_view sv) {
        write(sv.data(), sv.size());
    }

    /// Flush and close.  Called automatically by destructor.
    void close();

    Mode mode() const {
        return m_mode;
    }

    /// Detect mode from file extension (static helper).
    static Mode inferMode(const std::string &path);

    /// Convert compression string ("gz", "zst", or empty) to Mode.
    static Mode modeFromString(const std::string &comp);

    /// Build output path: prefix.pheno.method[.gz|.zst]
    static std::string buildOutputPath(
        const std::string &prefix,
        const std::string &phenoName,
        const std::string &methodName,
        const std::string &compression
    );

  private:
    void cleanup() noexcept;

    Mode m_mode = Mode::Plain;
    FILE *m_fp = nullptr;   // plain
    void *m_gz = nullptr;   // gzFile
    void *m_zctx = nullptr; // ZSTD_CCtx*
    bool m_closed = false;
};

// ══════════════════════════════════════════════════════════════════════
// TextReader
// ══════════════════════════════════════════════════════════════════════

class TextReader {
  public:
    enum class Mode { Plain, Gzip, Zstd };

    /// Open `path` for reading.  Decompression inferred from extension.
    explicit TextReader(const std::string &path);

    ~TextReader();

    TextReader(const TextReader &) = delete;
    TextReader &operator=(const TextReader &) = delete;

    TextReader(TextReader &&o) noexcept;
    TextReader &operator=(TextReader &&o) noexcept;

    /// Read one line (without trailing newline).  Returns false at EOF.
    bool getline(std::string &line);

    /// Close the file.  Called automatically by destructor.
    void close();

    Mode mode() const {
        return m_mode;
    }

    static Mode inferMode(const std::string &path);

  private:
    void cleanup() noexcept;

    Mode m_mode = Mode::Plain;
    FILE *m_fp = nullptr;
    void *m_gz = nullptr;    // gzFile
    void *m_zdctx = nullptr; // ZSTD_DCtx*
    FILE *m_zfp = nullptr;   // underlying FILE* for zstd input
    bool m_closed = false;

    // Buffered decompression state for zstd
    static constexpr size_t kZstdInBufSize = 256u * 1024u;
    static constexpr size_t kZstdOutBufSize = 256u * 1024u;
    char *m_zInBuf = nullptr;
    char *m_zOutBuf = nullptr;
    size_t m_zInPos = 0;
    size_t m_zInSize = 0;
    size_t m_zOutPos = 0;
    size_t m_zOutSize = 0;
    bool m_zEof = false; // underlying file exhausted

    // Fill decompressed output buffer; returns false when done
    bool zstdFillBuf();

};
