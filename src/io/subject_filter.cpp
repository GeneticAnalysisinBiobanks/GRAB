// subject_filter.cpp — PLINK2-compatible --keep / --remove file parser

#include "io/subject_filter.hpp"
#include "util/logging.hpp"
#include "util/text_scanner.hpp"

#include <fstream>
#include <stdexcept>

// ── Internal helpers ──────────────────────────────────────────────────────────

static std::vector<std::string> splitWS(const std::string &s) {
    text::TokenScanner ts(s);
    std::vector<std::string> tokens;
    while (!ts.atEnd()) {
        auto sv = ts.nextView();
        if (sv.empty()) break;
        tokens.emplace_back(sv);
    }
    return tokens;
}

// Determine which 0-based column index is the IID given the first line.
// Returns the IID column index (0-based) and whether the first line is a header.
static void detectFormat(
    const std::string &firstLine,
    int &iidCol,
    bool &hasHeader
) {
    if (!firstLine.empty() && firstLine[0] == '#') {
        // Header line
        hasHeader = true;
        auto tokens = splitWS(firstLine.substr(1)); // strip leading '#'
        if (!tokens.empty() && tokens[0] == "FID") {
            // "#FID IID [SID ...]"
            iidCol = 1;
        } else if (!tokens.empty() && tokens[0] == "IID") {
            // "#IID [SID ...]"
            iidCol = 0;
        } else {
            throw std::runtime_error("subject_filter: unrecognised header '" + firstLine +
                                     "'.\n"
                                     "Expected '#FID IID ...', '#FID IID SID ...', '#IID', or '#IID SID ...'");
        }
    } else {
        hasHeader = false;
        auto tokens = splitWS(firstLine);
        if (tokens.size() == 1) {
            iidCol = 0; // single-column: IID only
        } else {
            iidCol = 1; // multi-column: first col is FID, second is IID (plink1)
        }
    }
}

// ── Public API ────────────────────────────────────────────────────────────────

std::unordered_set<std::string> parseSubjectIDFile(const std::string &path) {
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("subject_filter: cannot open '" + path + "'");

    std::unordered_set<std::string> result;
    std::string line;

    // Read first non-empty line to detect format
    while (std::getline(ifs, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) break;
    }
    if (line.empty()) return result; // empty file

    int iidCol = 1;
    bool hasHeader = false;
    detectFormat(line, iidCol, hasHeader);

    // Process first line if it is a data line
    if (!hasHeader) {
        auto tokens = splitWS(line);
        if (iidCol < (int)tokens.size()) result.insert(tokens[(size_t)iidCol]);
    }

    // Remaining lines
    while (std::getline(ifs, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        auto tokens = splitWS(line);
        if (iidCol < (int)tokens.size()) result.insert(tokens[(size_t)iidCol]);
    }

    return result;
}

std::vector<uint32_t> buildKeptIndices(
    const std::vector<std::string> &sampleIDs,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const uint32_t N = static_cast<uint32_t>(sampleIDs.size());
    std::vector<uint32_t> kept;
    kept.reserve(N);

    // Start with all indices and apply filters
    std::vector<bool> include(N, true);

    if (!keepFile.empty()) {
        auto keepSet = parseSubjectIDFile(keepFile);
        uint32_t nMatched = 0;
        for (uint32_t i = 0; i < N; ++i) {
            if (keepSet.count(sampleIDs[i])) {
                ++nMatched;
            } else {
                include[i] = false;
            }
        }
        uint32_t nUnmatched = static_cast<uint32_t>(keepSet.size()) - nMatched;
        infoMsg("--keep: %u/%zu subjects in keep file matched genotype (%u not found)", nMatched, keepSet.size(),
                nUnmatched);
    }

    if (!removeFile.empty()) {
        auto removeSet = parseSubjectIDFile(removeFile);
        uint32_t nRemoved = 0;
        for (uint32_t i = 0; i < N; ++i) {
            if (include[i] && removeSet.count(sampleIDs[i])) {
                include[i] = false;
                ++nRemoved;
            }
        }
        infoMsg("--remove: %u subjects removed from analysis", nRemoved);
    }

    for (uint32_t i = 0; i < N; ++i)
        if (include[i]) kept.push_back(i);

    if (!keepFile.empty() || !removeFile.empty())
        infoMsg("Subject filter: %u / %u subjects retained", (uint32_t)kept.size(), N);

    return kept;
}
