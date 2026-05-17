// make_ldak_predlist.cpp — Build a grab --pred-list from LDAK Step 1 output.

#include "util/make_ldak_predlist.hpp"
#include "util/logging.hpp"

#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

namespace fs = std::filesystem;

struct LdakHit {
    fs::path singlePath;                 // <prefix>.step1.loco.prs (if any)
    std::map<int, fs::path> phenoPaths;  // 1-based idx → <prefix>.step1.phenoN.loco.prs
};

// Scan CWD and group LDAK Step 1 files by prefix. If `prefixFilter` is non-
// empty, only consider files whose captured prefix exactly matches it.
std::map<std::string, LdakHit> scanCwd(const std::string &prefixFilter) {
    std::map<std::string, LdakHit> hits;
    static const std::regex multiPat(
        R"((.*)\.step1\.pheno([0-9]+)\.loco\.prs)");
    static const std::regex singlePat(
        R"((.*)\.step1\.loco\.prs)");

    std::error_code ec;
    auto it = fs::directory_iterator(fs::current_path(), ec);
    if (ec)
        throw std::runtime_error("Cannot scan current directory '" +
                                 fs::current_path().string() +
                                 "' for LDAK Step 1 output: " + ec.message());

    for (const auto &entry : it) {
        if (!entry.is_regular_file(ec)) continue;
        const std::string name = entry.path().filename().string();

        std::smatch m;
        if (std::regex_match(name, m, multiPat)) {
            const std::string prefix = m[1].str();
            if (!prefixFilter.empty() && prefix != prefixFilter) continue;
            hits[prefix].phenoPaths[std::stoi(m[2].str())] = entry.path();
        } else if (std::regex_match(name, m, singlePat)) {
            const std::string prefix = m[1].str();
            if (!prefixFilter.empty() && prefix != prefixFilter) continue;
            hits[prefix].singlePath = entry.path();
        }
    }
    return hits;
}

// Tokenize a header line on whitespace.
std::vector<std::string> tokens(const std::string &line) {
    std::vector<std::string> out;
    const char *p   = line.c_str();
    const char *end = p + line.size();
    while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t')) ++p;
        if (p >= end) break;
        const char *s = p;
        while (p < end && *p != ' ' && *p != '\t') ++p;
        out.emplace_back(s, p - s);
    }
    return out;
}

} // namespace

void runMakeLdakPredlist(
    const std::string &phenoFile,
    const std::string &outPath,
    const std::string &ldakPrefix
) {
    if (phenoFile.empty())
        throw std::runtime_error("runMakeLdakPredlist: --pheno is required");
    if (outPath.empty())
        throw std::runtime_error("runMakeLdakPredlist: --out is required");

    // ── Read pheno header ────────────────────────────────────────────
    std::ifstream ifs(phenoFile);
    if (!ifs) throw std::runtime_error("Cannot open --pheno file: " + phenoFile);
    std::string line;
    if (!std::getline(ifs, line))
        throw std::runtime_error("Empty --pheno file: " + phenoFile);
    if (!line.empty() && line.back() == '\r') line.pop_back();

    const std::vector<std::string> header = tokens(line);
    // Accept either grab-style "IID Y1 ..." or PLINK/REGENIE-style "FID IID Y1 ...".
    std::size_t firstYCol = 0;
    if (header.size() >= 2 && header[0] == "IID") {
        firstYCol = 1;
    } else if (header.size() >= 3 && header[0] == "FID" && header[1] == "IID") {
        firstYCol = 2;
    } else {
        throw std::runtime_error(
            "--pheno header must start with 'IID' or 'FID IID' followed by "
            "Y columns; got: " + line);
    }
    std::vector<std::string> phenoNames(header.begin() + firstYCol, header.end());
    if (phenoNames.empty())
        throw std::runtime_error("--pheno has no Y columns after IID: " + phenoFile);
    // Duplicate Y names would silently collapse to one entry in the downstream
    // pred-list parser (last-write-wins on predMap[name]), so reject them up
    // front instead of producing a pred-list that maps two phenos to one PRS.
    {
        std::unordered_set<std::string> seen;
        for (const auto &name : phenoNames) {
            if (!seen.insert(name).second)
                throw std::runtime_error("--pheno header has duplicate Y column "
                                         "name: '" + name + "' in " + phenoFile);
        }
    }
    const std::size_t K = phenoNames.size();
    infoMsg("LDAK pred-list: %zu phenotype(s) in %s — %s",
            K, phenoFile.c_str(),
            [&] {
                std::string s;
                for (std::size_t i = 0; i < K; ++i) {
                    if (i) s += ',';
                    s += phenoNames[i];
                }
                return s;
            }().c_str());

    // ── Scan CWD ─────────────────────────────────────────────────────
    auto hits = scanCwd(ldakPrefix);
    if (hits.empty()) {
        std::string msg = "No LDAK Step 1 output (*.step1.loco.prs / "
                          "*.step1.phenoN.loco.prs) found in current directory: "
                        + fs::current_path().string();
        if (!ldakPrefix.empty())
            msg += "  (--prefix filter: '" + ldakPrefix + "')";
        throw std::runtime_error(msg);
    }

    // ── Match: prefix that covers all K phenoNames by count + position ──
    std::vector<std::pair<std::string, std::vector<fs::path>>> candidates;
    for (const auto &[prefix, hit] : hits) {
        std::vector<fs::path> paths;
        if (K == 1) {
            if (!hit.singlePath.empty())
                paths.push_back(hit.singlePath);
            else if (hit.phenoPaths.size() == 1 && hit.phenoPaths.count(1))
                paths.push_back(hit.phenoPaths.at(1));
        } else {
            bool full = true;
            for (std::size_t i = 1; i <= K; ++i) {
                auto it = hit.phenoPaths.find(static_cast<int>(i));
                if (it == hit.phenoPaths.end()) { full = false; break; }
                paths.push_back(it->second);
            }
            if (!full) paths.clear();
        }
        if (!paths.empty()) candidates.emplace_back(prefix, std::move(paths));
    }

    if (candidates.empty()) {
        std::ostringstream msg;
        msg << "No LDAK prefix in CWD covers all " << K << " phenotype(s). "
            << "Detected prefixes:";
        for (const auto &[prefix, hit] : hits) {
            msg << " '" << prefix << "'(single="
                << (hit.singlePath.empty() ? "no" : "yes")
                << ", phenoN={";
            bool first = true;
            for (const auto &[idx, _] : hit.phenoPaths) {
                if (!first) msg << ',';
                msg << idx;
                first = false;
            }
            msg << "})";
        }
        throw std::runtime_error(msg.str());
    }

    if (candidates.size() > 1) {
        std::string list;
        for (const auto &[prefix, _] : candidates) {
            if (!list.empty()) list += ", ";
            list += prefix;
        }
        throw std::runtime_error(
            "Ambiguous: multiple LDAK prefixes in CWD cover all " +
            std::to_string(K) + " phenotype(s) — " + list +
            ". Remove the extras or move the wanted prefix to a separate directory.");
    }

    const auto &[prefix, paths] = candidates.front();

    // ── Write pred-list ──────────────────────────────────────────────
    std::ofstream ofs(outPath);
    if (!ofs) throw std::runtime_error("Cannot write --out file: " + outPath);
    infoMsg("LDAK pred-list: matched prefix '%s' → writing %s",
            prefix.c_str(), outPath.c_str());
    for (std::size_t i = 0; i < K; ++i) {
        const std::string abs = fs::absolute(paths[i]).string();
        // The pred-list parser splits each line on whitespace, so a path with
        // an embedded space/tab gets silently truncated at the first split.
        // Refuse to write a broken file — surface the problem to the user.
        if (abs.find_first_of(" \t") != std::string::npos)
            throw std::runtime_error("LDAK pred-list: LOCO file path contains "
                                     "whitespace, which the pred-list parser "
                                     "would split: '" + abs + "'. Move the "
                                     "LDAK output (or working directory) to a "
                                     "path without spaces/tabs and re-run.");
        ofs << phenoNames[i] << '\t' << abs << '\n';
        infoMsg("  %s → %s", phenoNames[i].c_str(), abs.c_str());
    }
    ofs.close();
}
