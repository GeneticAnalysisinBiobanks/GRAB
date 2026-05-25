// int_pheno.hpp — Inverse-normal-transform utility for phenotype files.
//
// Reads a whitespace-separated phenotype file whose header follows the
// standard GRAB subject-data conventions:
//     #IID col1 col2 ...
//     IID  col1 col2 ...
//     #FID IID  col1 col2 ...
//     FID  IID  col1 col2 ...
// applies the inverse normal transform (Blom plotting position,
// average-rank ties) to the selected trait columns on their own
// non-missing scope, and writes the result to outPath.
//
// When traitNames is non-empty, only the listed columns are transformed
// and written; remaining trait columns are dropped from the output.
// When traitNames is empty, every trait column in the input is
// transformed and written.
//
// The output preserves the input header's ID-column layout (single
// "#IID" column vs. "FID IID" pair).  Entries that were missing in the
// input stay missing in the output (written as NA).
//
// Missing-value tokens recognized (same as the rest of GRAB): empty, ".",
// "NA", "na", "NaN", "nan", "-".
#pragma once

#include <string>
#include <vector>

void runIntPheno(
    const std::string &phenoFile,
    const std::string &outPath,
    const std::vector<std::string> &traitNames = {}
);
