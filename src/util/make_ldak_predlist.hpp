// make_ldak_predlist.hpp — Build a grab --pred-list file from LDAK-KVIK Step 1
// output discovered in the current working directory.
//
// Scans the CWD for either pattern:
//   <prefix>.step1.loco.prs           (single-pheno mode)
//   <prefix>.step1.phenoN.loco.prs    (--mpheno ALL: per-trait, N = 1, 2, …)
//
// Matches against the Y columns of `phenoFile` (every column after FID and
// IID) by count + position, then writes a pred-list file
//   <phenoName_1>\t<abs_loco_path_1>
//   <phenoName_2>\t<abs_loco_path_2>
//   ...
// to `outPath` (e.g. `<--out>.txt`). Throws on:
//   - pheno file unreadable or missing FID/IID header
//   - 0 matching LDAK prefixes in CWD
//   - 2+ matching LDAK prefixes in CWD (ambiguous)
//   - matched prefix doesn't cover all requested pheno columns
#pragma once

#include <string>

void runMakeLdakPredlist(
    const std::string &phenoFile,
    const std::string &outPath
);
