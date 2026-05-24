// variant_filter.hpp — Shared --extract / --exclude logic for all readers
//
// Provides format-agnostic helpers for variant ID filtering, used by the
// plink / pgen / bgen / vcf readers (and the admixture .abed reader) so
// that --extract and --exclude behave identically regardless of genotype
// format.  Pure C++17, no external dependencies.
#pragma once

#include "geno_factory/geno_data.hpp"

#include <string>
#include <unordered_set>
#include <vector>

namespace geno_factory {

// Read a single-column ID file (one ID per line; lines starting with '#'
// and blank lines are ignored).  Returns the set of IDs.  Throws on
// missing or unreadable file.
std::unordered_set<std::string> loadVariantIdSet(const std::string &path);

// Apply --extract / --exclude to the markers vector, mutating it in place.
//
//   extractFile empty ⇒ no include filter applied.
//   excludeFile empty ⇒ no exclude filter applied.
//   Both empty        ⇒ immediate return (no work, no allocation).
//
// When at least one filter is active and the result is empty, throws
// std::runtime_error so the caller surfaces a clear diagnostic instead
// of silently producing zero variants.  IDs present in the include /
// exclude file that do not match any marker are silently ignored, mirroring
// plink2 and grab's --keep / --remove behaviour.
void filterMarkersByIds(
    std::vector<GenoMeta::MarkerInfo> &markers,
    const std::string &extractFile,
    const std::string &excludeFile
);

} // namespace geno_factory
