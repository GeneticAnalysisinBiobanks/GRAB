// variant_filter.cpp — Shared --extract / --exclude logic for all readers

#include "geno_factory/variant_filter.hpp"

#include "util/text_scanner.hpp"

#include <fstream>
#include <stdexcept>

namespace geno_factory {

std::unordered_set<std::string> loadVariantIdSet(const std::string &path) {
    std::ifstream in(path);
    if (!in.is_open()) throw std::runtime_error("Cannot open variant ID file: " + path);

    std::unordered_set<std::string> ids;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        text::TokenScanner ts(line);
        auto tok = ts.nextView();
        if (tok.empty()) continue;
        ids.emplace(tok);
    }
    return ids;
}

void filterMarkersByIds(
    std::vector<GenoMeta::MarkerInfo> &markers,
    const std::string &extractFile,
    const std::string &excludeFile
) {
    const bool anyExtract = !extractFile.empty();
    const bool anyExclude = !excludeFile.empty();
    if (!anyExtract && !anyExclude) return;

    std::unordered_set<std::string> includeSet, excludeSet;
    if (anyExtract) includeSet = loadVariantIdSet(extractFile);
    if (anyExclude) excludeSet = loadVariantIdSet(excludeFile);

    std::vector<GenoMeta::MarkerInfo> filtered;
    filtered.reserve(markers.size());
    for (auto &m : markers) {
        if (anyExtract && includeSet.find(m.id) == includeSet.end()) continue;
        if (anyExclude && excludeSet.find(m.id) != excludeSet.end()) continue;
        filtered.push_back(std::move(m));
    }
    markers = std::move(filtered);

    if (markers.empty())
        throw std::runtime_error("No markers remain after --extract/--exclude");
}

} // namespace geno_factory
