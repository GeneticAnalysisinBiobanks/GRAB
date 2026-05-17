// int_pheno.hpp — Inverse-normal-transform utility for phenotype files.
//
// Reads a whitespace-separated phenotype file
//     FID  IID  Y1  Y2  ...
//     fam1 sample1 11.8 9.1
//     ...
// applies the inverse normal transform (Blom plotting position, average-rank
// ties) to every trait column on its own non-missing scope, and writes the
// result to outPath.
// FID/IID and the original column names are preserved; entries that were
// missing in the input stay missing in the output (written as NA).
//
// Missing-value tokens recognized (same as the rest of GRAB): empty, ".",
// "NA", "na", "NaN", "nan", "-".
#pragma once

#include <string>

void runIntPheno(
    const std::string &phenoFile,
    const std::string &outPath
);
