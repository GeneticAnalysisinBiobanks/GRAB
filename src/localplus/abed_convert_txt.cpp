// abed_convert_txt.cpp — Text-to-ABED conversion

#include "localplus/abed_convert_txt.hpp"
#include "localplus/abed_io.hpp"
#include "io/subject_filter.hpp"
#include "util/logging.hpp"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <zlib.h>

// ── Gzip line reader ─────────────────────────────────────────────────

static bool gzReadLine(
    gzFile gz,
    std::string &line
) {
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

static std::vector<std::string> splitTabs(const std::string &s) {
    std::vector<std::string> tokens;
    size_t start = 0;
    while (start <= s.size()) {
        auto pos = s.find('\t', start);
        if (pos == std::string::npos) {
            tokens.push_back(s.substr(start));
            break;
        }
        tokens.push_back(s.substr(start, pos - start));
        start = pos + 1;
    }
    return tokens;
}

// ── convertTextToAbed ───────────────────────────────────────────────
//
// Converts extract_tracts_fast_pgzip.py output to .abed.
// File naming: {prefix}.anc{k}.dosage[.gz] and {prefix}.anc{k}.hapcount[.gz]  (k=0,1,...)

static std::string findExtractTractsFile(
    const std::string &prefix,
    int k,
    const char *type,
    std::string &foundPath
) {
    // Try .gz first, then plain
    for (const char *ext : {".gz", ""}) {
        std::string path = prefix + ".anc" + std::to_string(k) + "." + type + ext;
        gzFile gz = gzopen(path.c_str(), "rb");
        if (gz) {
            gzclose(gz);
            foundPath = path;
            return ext;
        }
    }
    foundPath.clear();
    return "";
}

void convertTextToAbed(
    const std::string &textPrefix,
    const std::string &outPrefix,
    const std::string &keepFile,
    const std::string &removeFile,
    int nthreads
) {
    if (nthreads < 1) nthreads = 1;
    // Auto-detect K by scanning for anc0, anc1, ...
    int K = 0;
    while (true) {
        std::string found;
        if (findExtractTractsFile(textPrefix, K, "dosage", found).empty()) break;
        ++K;
    }

    if (K == 0)
        throw std::runtime_error("No ancestry files found with prefix '" + textPrefix +
                                 "'.\n"
                                 "Expected: " +
                                 textPrefix + ".anc0.dosage[.gz]");

    infoMsg("Detected %d ancestries from prefix '%s'", K, textPrefix.c_str());

    // Resolve a path for a given ancestry index and type
    auto resolvePath = [&](int k, bool isDosage, std::string &outPath) {
        findExtractTractsFile(textPrefix, k, isDosage ? "dosage" : "hapcount", outPath);
    };

    struct AncFile {
        gzFile dosGz;
        gzFile hapGz;
    };

    std::vector<AncFile> files(K);

    for (int k = 0; k < K; ++k) {
        std::string dosPath, hapPath;
        resolvePath(k, true, dosPath);
        resolvePath(k, false, hapPath);

        if (dosPath.empty()) throw std::runtime_error("Cannot find dosage file for anc" + std::to_string(k));
        if (hapPath.empty()) throw std::runtime_error("Cannot find hapcount file for anc" + std::to_string(k));

        files[k].dosGz = gzopen(dosPath.c_str(), "rb");
        files[k].hapGz = gzopen(hapPath.c_str(), "rb");
        if (!files[k].dosGz) throw std::runtime_error("Cannot open: " + dosPath);
        if (!files[k].hapGz) {
            gzclose(files[k].dosGz);
            throw std::runtime_error("Cannot open: " + hapPath);
        }
    }

    // Read header from first dosage file
    std::string headerLine, skipLine;
    if (!gzReadLine(files[0].dosGz, headerLine)) throw std::runtime_error("Empty dosage file for anc0");

    auto headerTokens = splitTabs(headerLine);
    if (headerTokens.size() < 6) throw std::runtime_error("Dosage file header has fewer than 6 columns");

    uint32_t nSamples = static_cast<uint32_t>(headerTokens.size() - 5);
    std::vector<std::string> sampleIDs(headerTokens.begin() + 5, headerTokens.end());
    infoMsg("Converting: K=%d, N=%u samples", K, nSamples);

    // Skip headers on remaining file handles
    if (!gzReadLine(files[0].hapGz, skipLine)) throw std::runtime_error("Empty hapcount file for anc0");
    for (int k = 1; k < K; ++k) {
        if (!gzReadLine(files[k].dosGz, skipLine) || !gzReadLine(files[k].hapGz, skipLine))
            throw std::runtime_error("Empty file for anc" + std::to_string(k));
    }

    // Apply --keep / --remove subject filtering
    auto keptIndices = buildKeptIndices(sampleIDs, keepFile, removeFile);
    uint32_t nKept = static_cast<uint32_t>(keptIndices.size());
    if (nKept == 0) throw std::runtime_error("convertTextToAbed: no subjects remain after --keep/--remove filters");

    // Build filtered sample ID list
    std::vector<std::string> keptSampleIDs(nKept);
    for (uint32_t j = 0; j < nKept; ++j)
        keptSampleIDs[j] = sampleIDs[keptIndices[j]];

    // Create ABED writer + output files  (no nMarkers needed in v2 header)
    // Text input may contain missing values, so default to noMissing=false;
    // we track hasMissing and the flag is set conservatively.
    AbedWriter writer(outPrefix + ".abed", static_cast<uint8_t>(K), nKept,
                      /*noMissing=*/ false, nthreads);

    std::ofstream bimOut(outPrefix + ".bim");
    if (!bimOut) throw std::runtime_error("Cannot create " + outPrefix + ".bim");

    std::ofstream famOut(outPrefix + ".fam");
    if (!famOut) throw std::runtime_error("Cannot create " + outPrefix + ".fam");
    for (const auto &sid : keptSampleIDs)
        famOut << sid << "\t" << sid << "\t0\t0\t0\t-9\n";
    famOut.close();

    const uint64_t bytesPerTrack = (static_cast<uint64_t>(nKept) + 3) / 4;
    std::vector<int8_t> values(nKept);
    std::string dosLine, hapLine;
    bool hasMissing = false;

    // Fast column encoder: scan tab-separated line directly, no string allocation.
    // Reads columns at sorted keptIndices[], skips others via forward scan.
    // nHeaderCols: number of non-sample columns before the sample data starts.
    auto fastEncodeKept = [&](const std::string &line, uint32_t nHeaderCols) {
        const char *p = line.c_str();
        // Skip the header columns (chrom, pos, id, ref, alt = 5 cols)
        for (uint32_t f = 0; f < nHeaderCols; ++f) {
            while (*p && *p != '\t')
                ++p;
            if (*p == '\t') ++p;
        }
        // p is now at the start of sample column 0.
        // keptIndices[] is sorted: scan left-to-right collecting requested columns.
        uint32_t col = 0;
        for (uint32_t ki = 0; ki < nKept; ++ki) {
            uint32_t target = keptIndices[ki];
            // Advance col to target by skipping (target - col) fields
            while (col < target) {
                while (*p && *p != '\t')
                    ++p;
                if (*p == '\t') ++p;
                ++col;
            }
            // Parse the value in [p, tab/end)
            const char *vStart = p;
            while (*p && *p != '\t')
                ++p;
            ptrdiff_t vLen = p - vStart;
            int8_t v;
            if (vLen == 0 || (vLen == 1 && *vStart == '.') || (vLen == 2 && vStart[0] == 'N' && vStart[1] == 'A') ||
                (vLen == 3 && vStart[0] == 'N' && vStart[1] == 'A' && vStart[2] == 'N') ||
                (vLen == 2 && vStart[0] == '-' && vStart[1] == '9')) {
                v = -1;
                hasMissing = true;
            } else {
                char *ep;
                long iv = std::strtol(vStart, &ep, 10);
                v = static_cast<int8_t>(iv >= 0 && iv <= 2 ? iv : -1);
                if (v < 0) hasMissing = true;
            }
            values[ki] = v;
            // p is now at '\t' or end; col still == target.
        }
        return AbedWriter::encode(values);
    };

    // Fast 5-field BIM extractor: scan the first 5 tab fields without string alloc.
    auto writeBimLine = [&](const std::string &line) {
        const char *p = line.c_str();
        auto nextField = [&]() -> std::string_view {
            const char *s = p;
            while (*p && *p != '\t')
                ++p;
            std::string_view sv(s, static_cast<size_t>(p - s));
            if (*p == '\t') ++p;
            return sv;
        };
        auto chrom = nextField(); // col 0
        auto pos = nextField();   // col 1
        auto id = nextField();    // col 2
        auto ref = nextField();   // col 3
        auto alt = nextField();   // col 4
        // BIM columns: CHROM ID CM POS REF ALT
        bimOut << chrom << '\t' << id << "\t0\t" << pos << '\t' << ref << '\t' << alt << '\n';
    };

    if (nthreads > 1) infoMsg("Note: --threads=%d has no effect on the text conversion path (I/O-bound)", nthreads);

    // Stream all markers until EOF — no pre-count needed.
    uint32_t nMarkers = 0;
    while (true) {
        // Read one line from each ancestry (dosage + hapcount)
        bool eof = false;
        for (int k = 0; k < K; ++k) {
            bool gotDos = gzReadLine(files[k].dosGz, dosLine);
            bool gotHap = gzReadLine(files[k].hapGz, hapLine);
            if (k == 0 && !gotDos) {
                eof = true;
                break;
            }
            if (!gotDos || !gotHap)
                throw std::runtime_error("Premature EOF at marker " + std::to_string(nMarkers) + " anc" +
                                         std::to_string(k));

            if (k == 0) writeBimLine(dosLine);

            auto packed = fastEncodeKept(dosLine, 5);
            writer.writeTrack(packed.data(), bytesPerTrack);
            packed = fastEncodeKept(hapLine, 5);
            writer.writeTrack(packed.data(), bytesPerTrack);
        }
        if (eof) break;
        ++nMarkers;

        if (nMarkers % 10000 == 0) infoMsg("  Converted %u markers", nMarkers);
    }

    writer.close();
    bimOut.close();
    for (int k = 0; k < K; ++k) {
        gzclose(files[k].dosGz);
        gzclose(files[k].hapGz);
    }

    infoMsg("Conversion complete: %s.abed (%u markers x %u samples x %d ancestries)", outPrefix.c_str(), nMarkers,
            nKept, K);
}
