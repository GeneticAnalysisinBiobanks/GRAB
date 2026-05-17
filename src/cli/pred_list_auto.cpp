// pred_list_auto.cpp — Auto-discover LDAK LOCO pred-list in the CWD.

#include "cli/pred_list_auto.hpp"
#include "util/logging.hpp"

#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <string>
#include <vector>

namespace cli {

namespace {

namespace fs = std::filesystem;

struct PrefixHit {
    fs::path singlePath;                 // <prefix>.step1.loco.prs (if any)
    std::map<int, fs::path> phenoPaths;  // 1-based idx → <prefix>.step1.phenoN.loco.prs
};

// Scan cwd, group by prefix.
std::map<std::string, PrefixHit> scan() {
    std::map<std::string, PrefixHit> hits;
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
        // Try multi-pheno first (more specific).
        if (std::regex_match(name, m, multiPat)) {
            const std::string prefix = m[1].str();
            const int idx = std::stoi(m[2].str());
            hits[prefix].phenoPaths[idx] = entry.path();
        } else if (std::regex_match(name, m, singlePat)) {
            const std::string prefix = m[1].str();
            hits[prefix].singlePath = entry.path();
        }
    }
    return hits;
}

} // namespace

std::string autoBuildPredList(const std::vector<std::string> &phenoNames) {
    const std::size_t K = phenoNames.size();
    if (K == 0) return "";

    auto hits = scan();
    if (hits.empty()) return "";

    // Collect prefixes that satisfy the requested coverage:
    //   K == 1 → either a single-pheno file OR a single phenoN file with N=1
    //   K  > 1 → phenoN files with N = 1..K all present
    std::vector<std::pair<std::string, std::vector<fs::path>>> candidates;

    for (const auto &[prefix, hit] : hits) {
        std::vector<fs::path> paths;

        if (K == 1) {
            if (!hit.singlePath.empty()) {
                paths.push_back(hit.singlePath);
            } else if (hit.phenoPaths.size() == 1 && hit.phenoPaths.count(1)) {
                paths.push_back(hit.phenoPaths.at(1));
            }
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

    if (candidates.empty()) return "";

    if (candidates.size() > 1) {
        std::string list;
        for (const auto &[prefix, _] : candidates) {
            if (!list.empty()) list += ", ";
            list += prefix;
        }
        warnMsg(
            "Auto pred-list: multiple matching LOCO prefixes in CWD (%s). "
            "Pass --pred-list explicitly to disambiguate; skipping LOCO.",
            list.c_str());
        return "";
    }

    // Single unambiguous match — write a temp pred-list.
    const auto &[prefix, paths] = candidates.front();
    const fs::path tmpPath = fs::temp_directory_path() / "grab_auto_predlist.txt";

    std::ofstream ofs(tmpPath);
    if (!ofs) {
        warnMsg("Auto pred-list: cannot write temp file %s; skipping LOCO.",
                tmpPath.string().c_str());
        return "";
    }

    infoMsg("Auto pred-list: matched LDAK prefix '%s' for %zu phenotype(s):",
            prefix.c_str(), K);
    for (std::size_t i = 0; i < K; ++i) {
        const std::string abs = fs::absolute(paths[i]).string();
        ofs << phenoNames[i] << '\t' << abs << '\n';
        infoMsg("  %s → %s", phenoNames[i].c_str(), abs.c_str());
    }
    ofs.close();
    return tmpPath.string();
}

} // namespace cli
