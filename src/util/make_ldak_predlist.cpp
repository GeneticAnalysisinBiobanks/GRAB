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
#include <vector>

namespace {

namespace fs = std::filesystem;

struct LdakHit {
    fs::path singlePath;                 // <prefix>.step1.loco.prs (if any)
    std::map<int, fs::path> phenoPaths;  // 1-based idx → <prefix>.step1.phenoN.loco.prs
};

// Scan CWD and group LDAK Step 1 files by prefix.
std::map<std::string, LdakHit> scanCwd() {
    std::map<std::string, LdakHit> hits;
    static const std::regex multiPat(
        R"((.*)\.step1\.pheno([0-9]+)\.loco\.prs)");
    static const std::regex singlePat(
        R"((.*)\.step1\.loco\.prs)");

    std::error_code ec;
    auto it = fs::directory_iterator(fs::current_path(), ec);
    if (ec) return hits;

    for (const auto &entry : it) {
        if (!entry.is_regular_file(ec)) continue;
        const std::string name = entry.path().filename().string();

        std::smatch m;
        if (std::regex_match(name, m, multiPat)) {
            hits[m[1].str()].phenoPaths[std::stoi(m[2].str())] = entry.path();
        } else if (std::regex_match(name, m, singlePat)) {
            hits[m[1].str()].singlePath = entry.path();
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
    const std::string &outPath
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
    if (header.size() < 3 || header[0] != "FID" || header[1] != "IID")
        throw std::runtime_error(
            "--pheno file must start with header 'FID IID Y1 ...': " + phenoFile);

    std::vector<std::string> phenoNames(header.begin() + 2, header.end());
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
    auto hits = scanCwd();
    if (hits.empty())
        throw std::runtime_error(
            "No LDAK Step 1 output (*.step1.loco.prs / *.step1.phenoN.loco.prs) "
            "found in current directory: " + fs::current_path().string());

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
        ofs << phenoNames[i] << '\t' << abs << '\n';
        infoMsg("  %s → %s", phenoNames[i].c_str(), abs.c_str());
    }
    ofs.close();
}
