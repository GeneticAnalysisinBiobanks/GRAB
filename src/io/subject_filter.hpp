// subject_filter.hpp — PLINK2-compatible --keep / --remove file parser
#pragma once

#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

// Parse a PLINK2-style sample ID file and return the set of IIDs it contains.
//
// Supported formats:
//   Header "#FID IID [SID ...]"   → IID is column 2
//   Header "#IID [SID ...]"       → IID is column 1
//   No header, 1 column           → column 1 treated as IID
//   No header, 2+ columns         → column 2 treated as IID (PLINK 1.x FID/IID)
//
// When FID is undefined it is treated as '0'.
std::unordered_set<std::string> parseSubjectIDFile(const std::string &path);

// Given an ordered list of sample IDs from the genotype file and optional
// keep/remove filter files, return the 0-based indices (into sampleIDs) of
// subjects that survive the filters.
//
//   keepFile   — empty string means keep everyone
//   removeFile — empty string means remove nobody
//
// Subjects in keepFile that are not present in sampleIDs are silently ignored
// (with an info message).  The returned vector preserves the original order.
std::vector<uint32_t>buildKeptIndices(
    const std::vector<std::string> &sampleIDs,
    const std::string &keepFile,
    const
    std::string &removeFile
);
