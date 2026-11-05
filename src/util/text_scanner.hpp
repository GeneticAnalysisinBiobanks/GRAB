// text_scanner.hpp — Shared text-parsing primitives
//
// Provides:
//   § 1  skipLine()       — strip trailing \r, return true if blank or '#' comment
//   § 2  TokenScanner     — zero-alloc token iteration (whitespace- or tab-delimited)
//   § 3  buildIIDMap()    — IID string → uint32_t index map
#pragma once

#include <cstdint>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace text {

// Strip trailing \r from `line`.  Returns true if the line should be skipped
// (blank or starts with '#').
inline bool skipLine(std::string &line) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    return line.empty() || line[0] == '#';
}

// Field-delimiter mode for TokenScanner.
//
//   Whitespace — runs of spaces and/or tabs delimit a field; consecutive
//                delimiters collapse, leading/trailing whitespace is trimmed,
//                and no empty field can be produced.  This is the historical
//                default and is used everywhere except the pheno/covar/resid
//                loaders.
//   Tab        — a single '\t' delimits a field; spaces are field *content*,
//                not delimiters.  This makes tab-separated files (`.tsv`) with
//                space-bearing values (e.g. `DRUG: CAPTOPRIL`) tokenise to the
//                correct column count.  Consecutive tabs and trimming are NOT
//                collapsed in the same way: each tab is exactly one boundary.
enum class Delim { Whitespace, Tab };

// Auto-detect the field delimiter from a header line: if the line contains a
// tab character it is treated as tab-separated, otherwise as whitespace-
// separated.  The rule is robust because a whitespace-delimited header never
// contains a tab, and a genuine TSV header always does; mixed files do not
// occur in practice.  Callers detect once on the header and reuse the result
// for every data row.
inline Delim detectDelim(const std::string &headerLine) {
    return headerLine.find('\t') != std::string::npos ? Delim::Tab
                                                       : Delim::Whitespace;
}

// Lightweight cursor for delimited tokenisation of a std::string.
// Designed to replace the repeated skipWS/nextTok lambda pattern.
//
// Usage:
//   TokenScanner tok(line);                 // whitespace-delimited (default)
//   TokenScanner tok(line, Delim::Tab);     // tab-delimited (.tsv)
//   std::string id1 = tok.next();
//   std::string id2 = tok.next();
//   tok.skipWS();
//   double val = std::strtod(tok.pos(), &endPtr);
struct TokenScanner {
    const char *p;
    const char *end;
    Delim delim;

    explicit TokenScanner(const std::string &line, Delim d = Delim::Whitespace)
        : p(line.c_str()),
          end(p + line.size()),
          delim(d)
    {
    }

    // Skip inter-field whitespace.  In Tab mode this is a no-op: spaces are
    // field content and tabs are consumed by next()/nextView() as the field
    // boundary, so there is nothing to skip between fields.
    void skipWS() {
        if (delim == Delim::Tab) return;
        while (p < end && (*p == ' ' || *p == '\t'))
            ++p;
    }

    const char *pos() const {
        return p;
    }

    bool atEnd() const {
        return p >= end;
    }

    std::string next() {
        return std::string(nextView());
    }

    // Zero-copy token as string_view (valid until the underlying string is destroyed).
    std::string_view nextView() {
        if (delim == Delim::Tab) {
            const char *s = p;
            while (p < end && *p != '\t')
                ++p;
            std::string_view tok{s, static_cast<size_t>(p - s)};
            if (p < end) ++p;   // consume the single '\t' field separator
            return tok;
        }
        skipWS();
        const char *s = p;
        while (p < end && *p != ' ' && *p != '\t')
            ++p;
        return {s, static_cast<size_t>(p - s)};
    }

    // Count remaining tokens without advancing.
    int countRemaining() const {
        const char *q = p;
        if (delim == Delim::Tab) {
            if (q >= end) return 0;
            int n = 1;
            while (q < end) { if (*q == '\t') ++n; ++q; }
            return n;
        }
        int n = 0;
        while (q < end) {
            while (q < end && (*q == ' ' || *q == '\t')) ++q;
            if (q >= end) break;
            ++n;
            while (q < end && *q != ' ' && *q != '\t') ++q;
        }
        return n;
    }

    // Count all tokens (call on a fresh scanner).
    int countAll() const {
        return countRemaining();
    }

};

// Build an IID → index map from a vector of subject IDs.
inline std::unordered_map<std::string, uint32_t> buildIIDMap(const std::vector<std::string> &iids) {
    std::unordered_map<std::string, uint32_t> m;
    m.reserve(iids.size());
    for (uint32_t i = 0; i < static_cast<uint32_t>(iids.size()); ++i)
        m.emplace(iids[i], i);
    return m;
}

} // namespace text
