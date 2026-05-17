// pred_list_auto.hpp — Auto-discover a LOCO pred-list in the current directory
// when the user omits --pred-list.
//
// Scans the working directory for LDAK-KVIK Step 1 output:
//   - <prefix>.step1.loco.prs            (single-pheno mode)
//   - <prefix>.step1.phenoN.loco.prs     (--mpheno ALL: per-trait, N = 1, 2, …)
//
// Matches the discovered files against the user's --pheno-name list by count
// (and, for the multi-pheno case, position). On a clean match, writes a
// temporary pred-list file and returns its path. Otherwise returns "" — the
// caller should silently fall back to "no LOCO" (the existing default when
// --pred-list is omitted).
#pragma once

#include <string>
#include <vector>

namespace cli {

// Returns:
//   - non-empty string  → path to a temp pred-list file matching `phenoNames`
//   - ""                → no match (zero / partial / ambiguous coverage)
//                         Caller falls back to no-LOCO behavior.
std::string autoBuildPredList(const std::vector<std::string> &phenoNames);

} // namespace cli
