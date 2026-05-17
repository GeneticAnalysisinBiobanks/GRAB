// pred_list_auto.hpp — Auto-discover a LOCO pred-list in the current directory
// when the user omits --pred-list (or passes the sentinel "auto").
//
// Scans the working directory for both supported LOCO sources:
//   (1) REGENIE native pred-list:    <prefix>_pred.list
//         Already in the grab --pred-list format (phenoName <TAB> locoPath).
//         Used directly — its path is returned as-is.
//   (2) LDAK-KVIK Step 1 outputs:
//         <prefix>.step1.loco.prs           (single-pheno mode)
//         <prefix>.step1.phenoN.loco.prs    (--mpheno ALL, positional)
//         Synthesized into a temp pred-list with positional binding.
//
// If exactly one source provides full coverage of the requested phenoNames →
// that source is used. Two+ candidates → ambiguous (strict: throw; lenient:
// warn + ""). Zero candidates → return "" (caller falls back to no-LOCO).
#pragma once

#include <string>
#include <vector>

namespace cli {

// Returns:
//   - non-empty string  → path to a temp pred-list file matching `phenoNames`
//   - ""                → no match (zero coverage). Caller falls back to
//                         no-LOCO behavior in both strict and lenient modes.
//
// `strict = true` (--pred-list auto explicit sentinel): ambiguous matches
//   (multiple matching prefixes in CWD) throw std::runtime_error so the
//   user is forced to disambiguate.
// `strict = false` (--pred-list omitted): ambiguous matches log a warning
//   and return "" — quietly fall back to no-LOCO.
std::string autoBuildPredList(
    const std::vector<std::string> &phenoNames,
    bool strict = false
);

} // namespace cli
