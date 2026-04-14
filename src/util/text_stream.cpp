// text_stream.cpp — TextWriter / TextReader implementation

#include "util/text_stream.hpp"

#include <algorithm>
#include <cstring>
#include <stdexcept>
#include <string>

#include <zlib.h>
#include <zstd.h>

// ══════════════════════════════════════════════════════════════════════
//  Helper: extension-based mode detection
// ══════════════════════════════════════════════════════════════════════

static bool endsWith(
    const std::string &s,
    const char *suffix
) {
    size_t n = std::strlen(suffix);
    return s.size() >= n && s.compare(s.size() - n, n, suffix) == 0;
}

TextWriter::Mode TextWriter::inferMode(const std::string &path) {
    if (endsWith(path, ".gz")) return Mode::Gzip;
    if (endsWith(path, ".zst")) return Mode::Zstd;
    return Mode::Plain;
}

TextWriter::Mode TextWriter::modeFromString(const std::string &comp) {
    if (comp == "gz") return Mode::Gzip;
    if (comp == "zst") return Mode::Zstd;
    return Mode::Plain;
}

std::string TextWriter::buildOutputPath(
    const std::string &prefix,
    const std::string &phenoName,
    const std::string &methodName,
    const std::string &compression
) {
    std::string path = prefix + "." + phenoName + "." + methodName;
    if (compression == "gz")
        path += ".gz";
    else if (compression == "zst")
        path += ".zst";
    return path;
}

TextReader::Mode TextReader::inferMode(const std::string &path) {
    if (endsWith(path, ".gz")) return Mode::Gzip;
    if (endsWith(path, ".zst")) return Mode::Zstd;
    return Mode::Plain;
}

// ══════════════════════════════════════════════════════════════════════
//  TextWriter
// ══════════════════════════════════════════════════════════════════════

TextWriter::TextWriter(const std::string &path)
    : TextWriter(path, inferMode(path), 0)
{
}

TextWriter::TextWriter(
    const std::string &path,
    Mode mode,
    int level
)
    : m_mode(mode)
{
    switch (m_mode) {
    case Mode::Gzip: {
        // gzopen mode: "wb" + optional level digit
        char gzMode[8] = "wb";
        if (level > 0 && level <= 9) {
            gzMode[2] = static_cast<char>('0' + level);
            gzMode[3] = '\0';
        }
        gzFile gz = gzopen(path.c_str(), gzMode);
        if (!gz) throw std::runtime_error("Cannot open gzip output: " + path);
        gzbuffer(gz, 256u * 1024u);
        m_gz = static_cast<void *>(gz);
        break;
    }
    case Mode::Zstd: {
        FILE *fp = std::fopen(path.c_str(), "wb");
        if (!fp) throw std::runtime_error("Cannot open zstd output: " + path);
        m_fp = fp;
        ZSTD_CCtx *cctx = ZSTD_createCCtx();
        if (!cctx) {
            std::fclose(fp);
            m_fp = nullptr;
            throw std::runtime_error("ZSTD_createCCtx failed");
        }
        if (level > 0) ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, level);
        m_zctx = static_cast<void *>(cctx);
        break;
    }
    default: {
        FILE *fp = std::fopen(path.c_str(), "w");
        if (!fp) throw std::runtime_error("Cannot open output: " + path);
        // 256 KB user-space buffer
        std::setvbuf(fp, nullptr, _IOFBF, 256u * 1024u);
        m_fp = fp;
        break;
    }
    }
}

TextWriter::~TextWriter()
{
    cleanup();
}

TextWriter::TextWriter(TextWriter &&o) noexcept
    : m_mode(o.m_mode),
      m_fp(o.m_fp),
      m_gz(o.m_gz),
      m_zctx(o.m_zctx),
      m_closed(o.m_closed)
{
    o.m_fp = nullptr;
    o.m_gz = nullptr;
    o.m_zctx = nullptr;
    o.m_closed = true;
}

TextWriter &TextWriter::operator=(TextWriter &&o) noexcept {
    if (this != &o) {
        cleanup();
        m_mode = o.m_mode;
        m_fp = o.m_fp;
        m_gz = o.m_gz;
        m_zctx = o.m_zctx;
        m_closed = o.m_closed;
        o.m_fp = nullptr;
        o.m_gz = nullptr;
        o.m_zctx = nullptr;
        o.m_closed = true;
    }
    return *this;
}

void TextWriter::write(
    const char *data,
    size_t len
) {
    if (m_closed || len == 0) return;
    switch (m_mode) {
    case Mode::Gzip:
        gzwrite(static_cast<gzFile>(m_gz), data, static_cast<unsigned>(len));
        break;
    case Mode::Zstd: {
        ZSTD_CCtx *cctx = static_cast<ZSTD_CCtx *>(m_zctx);
        // Compress in streaming fashion, writing to m_fp
        ZSTD_inBuffer in = {data, len, 0};
        char outBuf[64u * 1024u];
        while (in.pos < in.size) {
            ZSTD_outBuffer out = {outBuf, sizeof(outBuf), 0};
            ZSTD_compressStream2(cctx, &out, &in, ZSTD_e_continue);
            if (out.pos > 0) std::fwrite(outBuf, 1, out.pos, m_fp);
        }
        break;
    }
    default:
        std::fwrite(data, 1, len, m_fp);
        break;
    }
}

void TextWriter::close() {
    if (m_closed) return;
    m_closed = true;

    if (m_mode == Mode::Zstd && m_zctx && m_fp) {
        // Flush remaining compressed data
        ZSTD_CCtx *cctx = static_cast<ZSTD_CCtx *>(m_zctx);
        char outBuf[64u * 1024u];
        ZSTD_inBuffer in = {nullptr, 0, 0};
        size_t remaining;
        do {
            ZSTD_outBuffer out = {outBuf, sizeof(outBuf), 0};
            remaining = ZSTD_compressStream2(cctx, &out, &in, ZSTD_e_end);
            if (out.pos > 0) std::fwrite(outBuf, 1, out.pos, m_fp);
        } while (remaining > 0);
    }

    if (m_zctx) {
        ZSTD_freeCCtx(static_cast<ZSTD_CCtx *>(m_zctx));
        m_zctx = nullptr;
    }
    if (m_gz) {
        gzclose(static_cast<gzFile>(m_gz));
        m_gz = nullptr;
    }
    if (m_fp) {
        std::fclose(m_fp);
        m_fp = nullptr;
    }
}

void TextWriter::cleanup() noexcept {
    if (!m_closed) {
        m_closed = true;
        // Best-effort flush on destruction
        if (m_mode == Mode::Zstd && m_zctx && m_fp) {
            ZSTD_CCtx *cctx = static_cast<ZSTD_CCtx *>(m_zctx);
            char outBuf[64u * 1024u];
            ZSTD_inBuffer in = {nullptr, 0, 0};
            size_t remaining;
            do {
                ZSTD_outBuffer out = {outBuf, sizeof(outBuf), 0};
                remaining = ZSTD_compressStream2(cctx, &out, &in, ZSTD_e_end);
                if (out.pos > 0) std::fwrite(outBuf, 1, out.pos, m_fp);
            } while (remaining > 0);
        }
        if (m_zctx) {
            ZSTD_freeCCtx(static_cast<ZSTD_CCtx *>(m_zctx));
            m_zctx = nullptr;
        }
        if (m_gz) {
            gzclose(static_cast<gzFile>(m_gz));
            m_gz = nullptr;
        }
        if (m_fp) {
            std::fclose(m_fp);
            m_fp = nullptr;
        }
    }
}

// ══════════════════════════════════════════════════════════════════════
//  TextReader
// ══════════════════════════════════════════════════════════════════════

TextReader::TextReader(const std::string &path)
    : m_mode(inferMode(path))
{
    switch (m_mode) {
    case Mode::Gzip: {
        gzFile gz = gzopen(path.c_str(), "rb");
        if (!gz) throw std::runtime_error("Cannot open gzip input: " + path);
        gzbuffer(gz, 256u * 1024u);
        m_gz = static_cast<void *>(gz);
        break;
    }
    case Mode::Zstd: {
        FILE *fp = std::fopen(path.c_str(), "rb");
        if (!fp) throw std::runtime_error("Cannot open zstd input: " + path);
        m_zfp = fp;
        ZSTD_DCtx *dctx = ZSTD_createDCtx();
        if (!dctx) {
            std::fclose(fp);
            m_zfp = nullptr;
            throw std::runtime_error("ZSTD_createDCtx failed");
        }
        m_zdctx = static_cast<void *>(dctx);
        m_zInBuf = new char[kZstdInBufSize];
        m_zOutBuf = new char[kZstdOutBufSize];
        break;
    }
    default: {
        FILE *fp = std::fopen(path.c_str(), "r");
        if (!fp) throw std::runtime_error("Cannot open input: " + path);
        m_fp = fp;
        break;
    }
    }
}

TextReader::~TextReader()
{
    cleanup();
}

TextReader::TextReader(TextReader &&o) noexcept
    : m_mode(o.m_mode),
      m_fp(o.m_fp),
      m_gz(o.m_gz),
      m_zdctx(o.m_zdctx),
      m_zfp(o.m_zfp),
      m_closed(o.m_closed),
      m_zInBuf(o.m_zInBuf),
      m_zOutBuf(o.m_zOutBuf),
      m_zInPos(o.m_zInPos),
      m_zInSize(o.m_zInSize),
      m_zOutPos(o.m_zOutPos),
      m_zOutSize(o.m_zOutSize),
      m_zEof(o.m_zEof)
{
    o.m_fp = nullptr;
    o.m_gz = nullptr;
    o.m_zdctx = nullptr;
    o.m_zfp = nullptr;
    o.m_zInBuf = nullptr;
    o.m_zOutBuf = nullptr;
    o.m_closed = true;
}

TextReader &TextReader::operator=(TextReader &&o) noexcept {
    if (this != &o) {
        cleanup();
        m_mode = o.m_mode;
        m_fp = o.m_fp;
        m_gz = o.m_gz;
        m_zdctx = o.m_zdctx;
        m_zfp = o.m_zfp;
        m_closed = o.m_closed;
        m_zInBuf = o.m_zInBuf;
        m_zOutBuf = o.m_zOutBuf;
        m_zInPos = o.m_zInPos;
        m_zInSize = o.m_zInSize;
        m_zOutPos = o.m_zOutPos;
        m_zOutSize = o.m_zOutSize;
        m_zEof = o.m_zEof;
        o.m_fp = nullptr;
        o.m_gz = nullptr;
        o.m_zdctx = nullptr;
        o.m_zfp = nullptr;
        o.m_zInBuf = nullptr;
        o.m_zOutBuf = nullptr;
        o.m_closed = true;
    }
    return *this;
}

bool TextReader::getline(std::string &line) {
    line.clear();
    if (m_closed) return false;

    switch (m_mode) {
    case Mode::Gzip: {
        gzFile gz = static_cast<gzFile>(m_gz);
        char buf[8192];
        while (gzgets(gz, buf, sizeof(buf))) {
            size_t n = std::strlen(buf);
            if (n > 0 && buf[n - 1] == '\n') {
                line.append(buf, n - 1);
                // Strip \r
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            line.append(buf, n);
        }
        // EOF or error — return what we have
        if (!line.empty()) {
            if (line.back() == '\r') line.pop_back();
            return true;
        }
        return false;
    }
    case Mode::Zstd: {
        while (true) {
            // Scan decompressed buffer for newline — batch append between newlines
            const char *scanStart = m_zOutBuf + m_zOutPos;
            const char *scanEnd   = m_zOutBuf + m_zOutSize;
            const char *nl = static_cast<const char *>(
                std::memchr(scanStart, '\n', static_cast<size_t>(scanEnd - scanStart)));
            if (nl) {
                line.append(scanStart, static_cast<size_t>(nl - scanStart));
                m_zOutPos = static_cast<size_t>(nl - m_zOutBuf) + 1;
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            // No newline in buffer — append remainder and refill
            line.append(scanStart, static_cast<size_t>(scanEnd - scanStart));
            m_zOutPos = m_zOutSize;
            if (!zstdFillBuf()) {
                // No more data
                if (!line.empty()) {
                    if (line.back() == '\r') line.pop_back();
                    return true;
                }
                return false;
            }
        }
    }
    default: {
        // Plain text: use fgets-based loop for portability (MSVC lacks POSIX getline)
        char tmp[8192];
        while (std::fgets(tmp, sizeof(tmp), m_fp)) {
            size_t n = std::strlen(tmp);
            if (n > 0 && tmp[n - 1] == '\n') {
                line.append(tmp, n - 1);
                if (!line.empty() && line.back() == '\r') line.pop_back();
                return true;
            }
            line.append(tmp, n);
        }
        if (!line.empty()) {
            if (line.back() == '\r') line.pop_back();
            return true;
        }
        return false;
    }
    }
}

bool TextReader::zstdFillBuf() {
    m_zOutPos = 0;
    m_zOutSize = 0;
    if (m_zEof && m_zInPos >= m_zInSize) return false;

    ZSTD_DCtx *dctx = static_cast<ZSTD_DCtx *>(m_zdctx);

    // Refill compressed input if needed
    if (m_zInPos >= m_zInSize && !m_zEof) {
        m_zInSize = std::fread(m_zInBuf, 1, kZstdInBufSize, m_zfp);
        m_zInPos = 0;
        if (m_zInSize == 0) m_zEof = true;
    }

    if (m_zInPos >= m_zInSize) return false;

    ZSTD_inBuffer in = {m_zInBuf, m_zInSize, m_zInPos};
    ZSTD_outBuffer out = {m_zOutBuf, kZstdOutBufSize, 0};
    ZSTD_decompressStream(dctx, &out, &in);
    m_zInPos = in.pos;
    m_zOutSize = out.pos;
    return m_zOutSize > 0;
}

void TextReader::close() {
    if (m_closed) return;
    m_closed = true;
    if (m_zdctx) {
        ZSTD_freeDCtx(static_cast<ZSTD_DCtx *>(m_zdctx));
        m_zdctx = nullptr;
    }
    delete[] m_zInBuf;
    m_zInBuf = nullptr;
    delete[] m_zOutBuf;
    m_zOutBuf = nullptr;
    if (m_zfp) {
        std::fclose(m_zfp);
        m_zfp = nullptr;
    }
    if (m_gz) {
        gzclose(static_cast<gzFile>(m_gz));
        m_gz = nullptr;
    }
    if (m_fp) {
        std::fclose(m_fp);
        m_fp = nullptr;
    }
}

void TextReader::cleanup() noexcept {
    if (!m_closed) {
        m_closed = true;
        if (m_zdctx) {
            ZSTD_freeDCtx(static_cast<ZSTD_DCtx *>(m_zdctx));
            m_zdctx = nullptr;
        }
        delete[] m_zInBuf;
        m_zInBuf = nullptr;
        delete[] m_zOutBuf;
        m_zOutBuf = nullptr;
        if (m_zfp) {
            std::fclose(m_zfp);
            m_zfp = nullptr;
        }
        if (m_gz) {
            gzclose(static_cast<gzFile>(m_gz));
            m_gz = nullptr;
        }
        if (m_fp) {
            std::fclose(m_fp);
            m_fp = nullptr;
        }
    }
}
