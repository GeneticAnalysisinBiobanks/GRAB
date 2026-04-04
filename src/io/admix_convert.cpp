// admix_convert.cpp — Text-to-ABED conversion

#include "io/admix_convert.hpp"
#include "io/admix.hpp"
#include "util/logging.hpp"

#include <zlib.h>
#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

// ── Gzip line reader ─────────────────────────────────────────────────

static bool gzReadLine(gzFile gz, std::string& line) {
    line.clear();
    char buf[65536];
    while (gzgets(gz, buf, sizeof(buf)) != nullptr) {
        size_t len = std::strlen(buf);
        line.append(buf, len);
        if (len > 0 && buf[len - 1] == '\n') {
            while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
                line.pop_back();
            return true;
        }
    }
    return !line.empty();
}

static std::vector<std::string> splitTabs(const std::string& s) {
    std::vector<std::string> tokens;
    size_t start = 0;
    while (start <= s.size()) {
        auto pos = s.find('\t', start);
        if (pos == std::string::npos) { tokens.push_back(s.substr(start)); break; }
        tokens.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return tokens;
}


// ── Main conversion ──────────────────────────────────────────────────

void convertTextToAbed(
    const std::string& textPrefix,
    const std::vector<int>& ancIdx,
    const std::string& outPrefix
) {
    if (ancIdx.empty())
        throw std::runtime_error("No ancestry indices provided");

    int K = static_cast<int>(ancIdx.size());

    // Open all 2K files (dosage + hapcount for each ancestry)
    struct AncFile { gzFile dosGz; gzFile hapGz; };
    std::vector<AncFile> files(K);
    for (int k = 0; k < K; ++k) {
        std::string dosFile = textPrefix + std::to_string(ancIdx[k]) + "_Dosage.txt.gz";
        std::string hapFile = textPrefix + std::to_string(ancIdx[k]) + "_HapCount.txt.gz";

        files[k].dosGz = gzopen(dosFile.c_str(), "rb");
        if (!files[k].dosGz)
            throw std::runtime_error("Cannot open dosage file: " + dosFile);

        files[k].hapGz = gzopen(hapFile.c_str(), "rb");
        if (!files[k].hapGz) {
            gzclose(files[k].dosGz);
            throw std::runtime_error("Cannot open hapcount file: " + hapFile);
        }
    }

    // Read header from first dosage file to get sample IDs and validate
    std::string headerLine;
    if (!gzReadLine(files[0].dosGz, headerLine))
        throw std::runtime_error("Empty dosage file for ancestry " + std::to_string(ancIdx[0]));

    auto headerTokens = splitTabs(headerLine);
    if (headerTokens.size() < 6)
        throw std::runtime_error("Dosage file header has fewer than 6 columns");

    // First 5 columns: CHROM POS ID REF ALT
    uint32_t nSamples = static_cast<uint32_t>(headerTokens.size() - 5);
    std::vector<std::string> sampleIDs(headerTokens.begin() + 5, headerTokens.end());

    infoMsg("Converting admix text to ABED: K=%d, N=%u samples", K, nSamples);

    // Skip headers of remaining files
    std::string skipLine;
    if (!gzReadLine(files[0].hapGz, skipLine))
        throw std::runtime_error("Empty hapcount file for ancestry " + std::to_string(ancIdx[0]));
    for (int k = 1; k < K; ++k) {
        if (!gzReadLine(files[k].dosGz, skipLine))
            throw std::runtime_error("Empty dosage file for ancestry " + std::to_string(ancIdx[k]));
        if (!gzReadLine(files[k].hapGz, skipLine))
            throw std::runtime_error("Empty hapcount file for ancestry " + std::to_string(ancIdx[k]));
    }

    // First pass: count markers (read first dosage file line-by-line)
    // Actually, let's do a single pass: read all files in parallel, writing as we go.
    // We need to count markers first for the header, or buffer.
    // Strategy: count lines in first dosage file, then rewind all.

    uint32_t nMarkers = 0;
    {
        gzFile countGz = gzopen(
            (textPrefix + std::to_string(ancIdx[0]) + "_Dosage.txt.gz").c_str(), "rb");
        if (!countGz) throw std::runtime_error("Cannot reopen for counting");
        std::string cntLine;
        gzReadLine(countGz, cntLine);  // skip header
        while (gzReadLine(countGz, cntLine)) ++nMarkers;
        gzclose(countGz);
    }
    infoMsg("Detected %u markers", nMarkers);

    // Rewind all file handles (close and reopen)
    for (int k = 0; k < K; ++k) {
        gzclose(files[k].dosGz);
        gzclose(files[k].hapGz);
        std::string dosFile = textPrefix + std::to_string(ancIdx[k]) + "_Dosage.txt.gz";
        std::string hapFile = textPrefix + std::to_string(ancIdx[k]) + "_HapCount.txt.gz";
        files[k].dosGz = gzopen(dosFile.c_str(), "rb");
        files[k].hapGz = gzopen(hapFile.c_str(), "rb");
        if (!files[k].dosGz || !files[k].hapGz)
            throw std::runtime_error("Cannot reopen files for ancestry " + std::to_string(ancIdx[k]));
        // Skip headers
        gzReadLine(files[k].dosGz, skipLine);
        gzReadLine(files[k].hapGz, skipLine);
    }

    // Create ABED writer (flag set at end if no missing data found)
    AbedWriter writer(outPrefix + ".abed", static_cast<uint16_t>(K), nSamples, nMarkers,
                      ABED_FLAG_NO_MISSING);

    // Open .bim for writing
    std::ofstream bimOut(outPrefix + ".bim");
    if (!bimOut) throw std::runtime_error("Cannot create " + outPrefix + ".bim");

    // Open .fam for writing
    std::ofstream famOut(outPrefix + ".fam");
    if (!famOut) throw std::runtime_error("Cannot create " + outPrefix + ".fam");
    for (const auto& sid : sampleIDs)
        famOut << sid << "\t" << sid << "\t0\t0\t0\t-9\n";
    famOut.close();

    uint64_t bytesPerTrack = (static_cast<uint64_t>(nSamples) + 3) / 4;
    std::vector<int8_t> values(nSamples);
    std::string dosLine, hapLine;
    bool hasMissing = false;

    for (uint32_t m = 0; m < nMarkers; ++m) {
        // For each ancestry, write dosage track then hapcount track
        for (int k = 0; k < K; ++k) {
            // Read dosage line
            if (!gzReadLine(files[k].dosGz, dosLine))
                throw std::runtime_error("Premature EOF in dosage file for ancestry "
                                         + std::to_string(ancIdx[k]) + " at marker " + std::to_string(m));
            if (!gzReadLine(files[k].hapGz, hapLine))
                throw std::runtime_error("Premature EOF in hapcount file for ancestry "
                                         + std::to_string(ancIdx[k]) + " at marker " + std::to_string(m));

            auto dosToks = splitTabs(dosLine);
            auto hapToks = splitTabs(hapLine);

            // Write .bim from first ancestry, first marker fields
            if (k == 0) {
                if (dosToks.size() < 5 + nSamples)
                    throw std::runtime_error("Dosage line too short at marker " + std::to_string(m));
                // BIM: CHROM ID CM POS REF ALT
                bimOut << dosToks[0] << "\t" << dosToks[2] << "\t0\t"
                       << dosToks[1] << "\t" << dosToks[3] << "\t" << dosToks[4] << "\n";
            }

            // Encode dosage track
            for (uint32_t s = 0; s < nSamples; ++s) {
                const auto& v = dosToks[5 + s];
                if (v.empty() || v == "NA" || v == "." || v == "-9") {
                    values[s] = -1; hasMissing = true;
                } else {
                    int iv = std::stoi(v);
                    values[s] = static_cast<int8_t>(iv >= 0 && iv <= 2 ? iv : -1);
                    if (values[s] < 0) hasMissing = true;
                }
            }
            auto packed = AbedWriter::encode(values);
            writer.writeTrack(packed.data(), bytesPerTrack);

            // Encode hapcount track
            for (uint32_t s = 0; s < nSamples; ++s) {
                const auto& v = hapToks[5 + s];
                if (v.empty() || v == "NA" || v == "." || v == "-9") {
                    values[s] = -1; hasMissing = true;
                } else {
                    int iv = std::stoi(v);
                    values[s] = static_cast<int8_t>(iv >= 0 && iv <= 2 ? iv : -1);
                    if (values[s] < 0) hasMissing = true;
                }
            }
            packed = AbedWriter::encode(values);
            writer.writeTrack(packed.data(), bytesPerTrack);
        }

        if ((m + 1) % 10000 == 0)
            infoMsg("  Converted %u / %u markers", m + 1, nMarkers);
    }

    if (hasMissing)
        writer.clearFlags(ABED_FLAG_NO_MISSING);

    writer.close();
    bimOut.close();

    for (int k = 0; k < K; ++k) {
        gzclose(files[k].dosGz);
        gzclose(files[k].hapGz);
    }

    infoMsg("Conversion complete: %s.abed (%u markers x %u samples x %d ancestries)",
            outPrefix.c_str(), nMarkers, nSamples, K);
}


// ── convertExtractTractsToAbed ───────────────────────────────────────
//
// Converts extract_tracts_fast_pgzip.py output to .abed.
//
// Supports two naming conventions (auto-detected):
//   New:  {prefix}.anc{k}.dosage.txt[.gz]  and  {prefix}.anc{k}.hapcount.txt[.gz]  (k=0,1,...)
//   Old:  {prefix}{k}_Dosage.txt.gz         and  {prefix}{k}_HapCount.txt.gz         (k=1,2,...)

static std::string findExtractTractsFile(const std::string& prefix, int k,
                                          const char* type, std::string& foundPath) {
    // Try .txt.gz first, then .txt
    for (const char* ext : {".txt.gz", ".txt"}) {
        std::string path = prefix + ".anc" + std::to_string(k) + "." + type + ext;
        // Try opening (gzopen handles both plain and gzip transparently)
        gzFile gz = gzopen(path.c_str(), "rb");
        if (gz) { gzclose(gz); foundPath = path; return ext; }
    }
    foundPath.clear();
    return "";
}

// Old naming: {prefix}{k}_Dosage.txt.gz / {prefix}{k}_HapCount.txt.gz (k = 1-based)
static bool findOldStyleFile(const std::string& prefix, int k,
                              bool isDosage, std::string& foundPath) {
    const char* tag   = isDosage ? "_Dosage" : "_HapCount";
    for (const char* ext : {".txt.gz", ".txt"}) {
        std::string path = prefix + std::to_string(k) + tag + ext;
        gzFile gz = gzopen(path.c_str(), "rb");
        if (gz) { gzclose(gz); foundPath = path; return true; }
    }
    foundPath.clear();
    return false;
}

void convertExtractTractsToAbed(
    const std::string& textPrefix,
    const std::string& outPrefix
) {
    // Auto-detect naming convention and K
    bool oldNaming = false;
    int K = 0;

    // Try new naming first (k=0,1,...)
    while (true) {
        std::string found;
        if (findExtractTractsFile(textPrefix, K, "dosage", found).empty()) break;
        ++K;
    }

    // Fall back to old naming (k=1,2,...) if new naming finds nothing
    if (K == 0) {
        int kOld = 1;
        while (true) {
            std::string found;
            if (!findOldStyleFile(textPrefix, kOld, true, found)) break;
            ++kOld;
        }
        K = kOld - 1;
        if (K > 0) oldNaming = true;
    }

    if (K == 0)
        throw std::runtime_error(
            "No ancestry files found with prefix '" + textPrefix + "'.\n"
            "Expected: " + textPrefix + ".anc0.dosage.txt[.gz]  (new naming)\n"
            "      or: " + textPrefix + "1_Dosage.txt[.gz]       (old naming)");

    infoMsg("Detected %d ancestries from prefix '%s' (%s naming)",
            K, textPrefix.c_str(), oldNaming ? "old" : "new");

    // Helper: resolve a path given naming convention + index (0-based for new, 1-based for old)
    auto resolvePath = [&](int k, bool isDosage, std::string& outPath) {
        if (oldNaming) {
            findOldStyleFile(textPrefix, k + 1, isDosage, outPath);
        } else {
            findExtractTractsFile(textPrefix, k, isDosage ? "dosage" : "hapcount", outPath);
        }
    };

    struct AncFile { gzFile dosGz; gzFile hapGz; };
    std::vector<AncFile> files(K);

    for (int k = 0; k < K; ++k) {
        std::string dosPath, hapPath;
        resolvePath(k, true,  dosPath);
        resolvePath(k, false, hapPath);

        if (dosPath.empty())
            throw std::runtime_error("Cannot find dosage file for anc" + std::to_string(k));
        if (hapPath.empty())
            throw std::runtime_error("Cannot find hapcount file for anc" + std::to_string(k));

        files[k].dosGz = gzopen(dosPath.c_str(), "rb");
        files[k].hapGz = gzopen(hapPath.c_str(), "rb");
        if (!files[k].dosGz)
            throw std::runtime_error("Cannot open: " + dosPath);
        if (!files[k].hapGz) {
            gzclose(files[k].dosGz);
            throw std::runtime_error("Cannot open: " + hapPath);
        }
    }

    // Read header from first dosage file
    std::string headerLine, skipLine;
    if (!gzReadLine(files[0].dosGz, headerLine))
        throw std::runtime_error("Empty dosage file for anc0");

    auto headerTokens = splitTabs(headerLine);
    if (headerTokens.size() < 6)
        throw std::runtime_error("Dosage file header has fewer than 6 columns");

    uint32_t nSamples = static_cast<uint32_t>(headerTokens.size() - 5);
    std::vector<std::string> sampleIDs(headerTokens.begin() + 5, headerTokens.end());
    infoMsg("Converting: K=%d, N=%u samples", K, nSamples);

    // Skip headers on remaining file handles
    if (!gzReadLine(files[0].hapGz, skipLine))
        throw std::runtime_error("Empty hapcount file for anc0");
    for (int k = 1; k < K; ++k) {
        if (!gzReadLine(files[k].dosGz, skipLine) || !gzReadLine(files[k].hapGz, skipLine))
            throw std::runtime_error("Empty file for anc" + std::to_string(k));
    }

    // Count markers by scanning first dosage file
    uint32_t nMarkers = 0;
    {
        std::string dosPath;
        resolvePath(0, true, dosPath);
        gzFile countGz = gzopen(dosPath.c_str(), "rb");
        if (!countGz) throw std::runtime_error("Cannot reopen for counting");
        std::string cntLine;
        gzReadLine(countGz, cntLine);  // skip header
        while (gzReadLine(countGz, cntLine)) ++nMarkers;
        gzclose(countGz);
    }
    infoMsg("Detected %u markers", nMarkers);

    // Rewind all files
    for (int k = 0; k < K; ++k) {
        gzclose(files[k].dosGz);
        gzclose(files[k].hapGz);
        std::string dosPath, hapPath;
        resolvePath(k, true,  dosPath);
        resolvePath(k, false, hapPath);
        files[k].dosGz = gzopen(dosPath.c_str(), "rb");
        files[k].hapGz = gzopen(hapPath.c_str(), "rb");
        if (!files[k].dosGz || !files[k].hapGz)
            throw std::runtime_error("Cannot reopen files for anc" + std::to_string(k));
        gzReadLine(files[k].dosGz, skipLine);
        gzReadLine(files[k].hapGz, skipLine);
    }

    // Create ABED writer + output files
    AbedWriter writer(outPrefix + ".abed", static_cast<uint16_t>(K), nSamples, nMarkers,
                      ABED_FLAG_NO_MISSING);

    std::ofstream bimOut(outPrefix + ".bim");
    if (!bimOut) throw std::runtime_error("Cannot create " + outPrefix + ".bim");

    std::ofstream famOut(outPrefix + ".fam");
    if (!famOut) throw std::runtime_error("Cannot create " + outPrefix + ".fam");
    for (const auto& sid : sampleIDs)
        famOut << sid << "\t" << sid << "\t0\t0\t0\t-9\n";
    famOut.close();

    uint64_t bytesPerTrack = (static_cast<uint64_t>(nSamples) + 3) / 4;
    std::vector<int8_t> values(nSamples);
    std::string dosLine, hapLine;
    bool hasMissing = false;

    auto encodeCol = [&](const std::vector<std::string>& toks, uint32_t off) {
        for (uint32_t s = 0; s < nSamples; ++s) {
            const auto& v = toks[off + s];
            if (v.empty() || v == "NA" || v == "." || v == "-9") {
                values[s] = -1; hasMissing = true;
            } else {
                int iv = std::stoi(v);
                values[s] = static_cast<int8_t>(iv >= 0 && iv <= 2 ? iv : -1);
                if (values[s] < 0) hasMissing = true;
            }
        }
        return AbedWriter::encode(values);
    };

    for (uint32_t m = 0; m < nMarkers; ++m) {
        for (int k = 0; k < K; ++k) {
            if (!gzReadLine(files[k].dosGz, dosLine))
                throw std::runtime_error("Premature EOF at marker " + std::to_string(m) + " anc" + std::to_string(k));
            if (!gzReadLine(files[k].hapGz, hapLine))
                throw std::runtime_error("Premature EOF (hapcount) at marker " + std::to_string(m) + " anc" + std::to_string(k));

            auto dosToks = splitTabs(dosLine);
            auto hapToks = splitTabs(hapLine);

            if (k == 0) {
                if (dosToks.size() < 5 + nSamples)
                    throw std::runtime_error("Dosage line too short at marker " + std::to_string(m));
                bimOut << dosToks[0] << "\t" << dosToks[2] << "\t0\t"
                       << dosToks[1] << "\t" << dosToks[3] << "\t" << dosToks[4] << "\n";
            }

            auto packed = encodeCol(dosToks, 5); writer.writeTrack(packed.data(), bytesPerTrack);
            packed       = encodeCol(hapToks, 5); writer.writeTrack(packed.data(), bytesPerTrack);
        }

        if ((m + 1) % 10000 == 0)
            infoMsg("  Converted %u / %u markers", m + 1, nMarkers);
    }

    if (hasMissing)
        writer.clearFlags(ABED_FLAG_NO_MISSING);

    writer.close();
    bimOut.close();
    for (int k = 0; k < K; ++k) { gzclose(files[k].dosGz); gzclose(files[k].hapGz); }

    infoMsg("Conversion complete: %s.abed (%u markers x %u samples x %d ancestries)",
            outPrefix.c_str(), nMarkers, nSamples, K);
}
