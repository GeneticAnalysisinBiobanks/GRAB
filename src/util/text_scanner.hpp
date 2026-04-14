// text_scanner.hpp — Shared text-parsing primitives
//
// Provides:
//   § 1  skipLine()       — strip trailing \r, return true if blank or '#' comment
//   § 2  TokenScanner     — zero-alloc whitespace-delimited token iteration
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

// Lightweight cursor for whitespace-delimited tokenisation of a std::string.
// Designed to replace the repeated skipWS/nextTok lambda pattern.
//
// Usage:
//   TokenScanner tok(line);
//   std::string id1 = tok.next();
//   std::string id2 = tok.next();
//   tok.skipWS();
//   double val = std::strtod(tok.pos(), &endPtr);
struct TokenScanner {
    const char *p;
    const char *end;

    explicit TokenScanner(const std::string &line)
        : p(line.c_str()),
          end(p + line.size())
    {
    }

    void skipWS() {
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
        skipWS();
        const char *s = p;
        while (p < end && *p != ' ' && *p != '\t')
            ++p;
        return std::string(s, p);
    }

    // Zero-copy token as string_view (valid until the underlying string is destroyed).
    std::string_view nextView() {
        skipWS();
        const char *s = p;
        while (p < end && *p != ' ' && *p != '\t')
            ++p;
        return {s, static_cast<size_t>(p - s)};
    }

    // Count remaining whitespace-delimited tokens without advancing.
    int countRemaining() const {
        const char *q = p;
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
