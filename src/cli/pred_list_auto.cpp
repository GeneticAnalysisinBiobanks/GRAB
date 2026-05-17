// pred_list_auto.cpp — Auto-discover LOCO pred-list in the CWD.
//
// Two file shapes are recognized:
//   (1) REGENIE-style native pred-list:  <prefix>_pred.list
//         A complete pred-list already in the format grab --pred-list expects
//         (phenoName <TAB> locoPath per row). Returned directly, no rewriting.
//   (2) LDAK-KVIK Step 1 LOCO output:
//         <prefix>.step1.loco.prs           (single-pheno mode)
//         <prefix>.step1.phenoN.loco.prs    (--mpheno ALL: per-trait, N = 1, 2, …)
//         Synthesized into a temp pred-list with positional name binding.
//
// Match order: REGENIE first (it's the more "direct" form, no synthesis), then
// LDAK. If both LDAK and REGENIE produce a matching candidate in the same CWD,
// that is reported as ambiguous (the user explicitly has to disambiguate).

#include "cli/pred_list_auto.hpp"
#include "util/logging.hpp"

#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace cli {

namespace {

namespace fs = std::filesystem;

// ── REGENIE _pred.list parse ──────────────────────────────────────────
// Reads "phenoName <TAB-or-space> locoPath" lines, returns the set of
// pheno-names present. Empty set on parse error / unreadable file.
std::set<std::string> regeniePredListPhenos(const fs::path &p) {
    std::set<std::string> names;
    std::ifstream ifs(p);
    if (!ifs) return names;
    std::string line;
    while (std::getline(ifs, line)) {
        std::istringstream iss(line);
        std::string ph;
        iss >> ph;
        if (!ph.empty()) names.insert(ph);
    }
    return names;
}

// ── LDAK output grouping ──────────────────────────────────────────────
struct LdakPrefixHit {
    fs::path singlePath;                 // <prefix>.step1.loco.prs (if any)
    std::map<int, fs::path> phenoPaths;  // 1-based idx → <prefix>.step1.phenoN.loco.prs
};

struct Scan {
    std::vector<fs::path>             regenieFiles;  // *_pred.list paths
    std::map<std::string, LdakPrefixHit> ldakByPrefix;
};

Scan scanCwd() {
    Scan out;
    static const std::regex ldakMultiPat(
        R"((.*)\.step1\.pheno([0-9]+)\.loco\.prs)");
    static const std::regex ldakSinglePat(
        R"((.*)\.step1\.loco\.prs)");
    static const std::regex regenieListPat(
        R"((.*)_pred\.list)");

    std::error_code ec;
    auto it = fs::directory_iterator(fs::current_path(), ec);
    if (ec) return out;

    for (const auto &entry : it) {
        if (!entry.is_regular_file(ec)) continue;
        const std::string name = entry.path().filename().string();

        std::smatch m;
        if (std::regex_match(name, m, ldakMultiPat)) {
            const std::string prefix = m[1].str();
            const int idx = std::stoi(m[2].str());
            out.ldakByPrefix[prefix].phenoPaths[idx] = entry.path();
        } else if (std::regex_match(name, m, ldakSinglePat)) {
            const std::string prefix = m[1].str();
            out.ldakByPrefix[prefix].singlePath = entry.path();
        } else if (std::regex_match(name, m, regenieListPat)) {
            out.regenieFiles.push_back(entry.path());
        }
    }
    return out;
}

struct Candidate {
    std::string source;            // "REGENIE" or "LDAK"
    std::string prefixOrPath;      // for logging
    fs::path    regeniePath;       // non-empty if source == "REGENIE"
    std::vector<fs::path> ldakPaths; // non-empty if source == "LDAK"
};

} // namespace

std::string autoBuildPredList(
    const std::vector<std::string> &phenoNames,
    bool strict
) {
    const std::size_t K = phenoNames.size();
    if (K == 0) return "";

    const Scan sc = scanCwd();
    if (sc.regenieFiles.empty() && sc.ldakByPrefix.empty()) return "";

    std::vector<Candidate> candidates;

    // ── REGENIE: a *_pred.list whose pheno-name set covers every requested
    //            phenoName ─────────────────────────────────────────────
    for (const auto &p : sc.regenieFiles) {
        const auto names = regeniePredListPhenos(p);
        bool covers = true;
        for (const auto &pn : phenoNames)
            if (!names.count(pn)) { covers = false; break; }
        if (covers) {
            candidates.push_back({"REGENIE", p.filename().string(), p, {}});
        }
    }

    // ── LDAK: a prefix whose .step1.loco.prs (or phenoN files) covers
    //         every requested phenoName positionally ─────────────────
    for (const auto &[prefix, hit] : sc.ldakByPrefix) {
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
        if (!paths.empty())
            candidates.push_back({"LDAK", prefix, {}, std::move(paths)});
    }

    if (candidates.empty()) return "";

    if (candidates.size() > 1) {
        std::string list;
        for (const auto &c : candidates) {
            if (!list.empty()) list += ", ";
            list += c.source + ":" + c.prefixOrPath;
        }
        if (strict) {
            throw std::runtime_error(
                "--pred-list auto: multiple matching LOCO candidates in CWD (" +
                list + "). Pass --pred-list <FILE> explicitly to disambiguate.");
        }
        warnMsg(
            "Auto pred-list: multiple matching LOCO candidates in CWD (%s). "
            "Pass --pred-list explicitly to disambiguate; skipping LOCO.",
            list.c_str());
        return "";
    }

    // Single unambiguous candidate.
    const Candidate &c = candidates.front();

    if (c.source == "REGENIE") {
        const std::string abs = fs::absolute(c.regeniePath).string();
        infoMsg("Auto pred-list: matched REGENIE pred-list '%s' for %zu phenotype(s)",
                c.prefixOrPath.c_str(), K);
        return abs;
    }

    // LDAK — synthesize a temp pred-list.
    const fs::path tmpPath = fs::temp_directory_path() / "grab_auto_predlist.txt";
    std::ofstream ofs(tmpPath);
    if (!ofs) {
        warnMsg("Auto pred-list: cannot write temp file %s; skipping LOCO.",
                tmpPath.string().c_str());
        return "";
    }

    infoMsg("Auto pred-list: matched LDAK prefix '%s' for %zu phenotype(s):",
            c.prefixOrPath.c_str(), K);
    for (std::size_t i = 0; i < K; ++i) {
        const std::string abs = fs::absolute(c.ldakPaths[i]).string();
        ofs << phenoNames[i] << '\t' << abs << '\n';
        infoMsg("  %s → %s", phenoNames[i].c_str(), abs.c_str());
    }
    ofs.close();
    return tmpPath.string();
}

} // namespace cli
