// longitudinal_resid.hpp — shared --longitudinal driver for SPACox / SPAmix /
//                          SPAGRM.
//
// These three single-residual marginal-association methods can accept a
// long-format (repeated-measures) phenotype directly.  This module fits the
// random-intercept null model  Y ~ X + (1 | IID)  per outcome and aggregates
// the per-row BLUP residuals to one value per subject,
//   R_G_i = Σ_j r_ij,
// which is then handed to each method's existing residual-mode marker test
// (the main genetic effect only — there is no environment / random-slope term
// and hence no --envir-name and no G×E).
//
// The fitting reuses nsSAGELDFit::parseLongPheno / fitRandomInterceptWithCovar
// / aggregatePerIID; SAGELD's validated random-slope path is not touched.
#pragma once

#include <Eigen/Dense>
#include <string>
#include <unordered_set>
#include <vector>

namespace nsLongitudinal {

struct LongResidResult {
    std::vector<std::string>     usedIIDs;     // .fam-ordered fit subjects
    std::vector<Eigen::VectorXd> R_G;          // one per outcome, length usedIIDs.size()
    std::vector<std::string>     names;        // outcome names (= phenoNames)
    Eigen::MatrixXd              perSubjectX;  // usedIIDs.size() × p, first row per IID
    std::vector<std::string>     xColNames;    // length p: "intercept", covarNames...
};

// Split a comma-separated --pheno-name string into outcome column names.
// Throws if empty or if a token carries a ':' (the Cox TIME:EVENT syntax is
// not valid for a longitudinal quantitative outcome).
std::vector<std::string> splitOutcomeNames(const std::string &commaList);

// Build the kept-subject filter  GRM ∩ (keep) − (remove)  in .fam order, the
// same intersection SAGELD pheno-mode passes to parseLongPheno.  An empty
// grmSubjects skips the GRM step; empty keep/remove paths skip those filters.
std::unordered_set<std::string> buildKeptSet(
    const std::string &keepFile,
    const std::string &removeFile,
    const std::vector<std::string> &famIIDs,
    const std::unordered_set<std::string> &grmSubjects);

// Fit Y ~ X + (1 | IID) per outcome on a long-format phenotype file and
// aggregate per-IID residuals R_G.  keptSubjects (when non-empty) restricts
// the fit sample (pass buildKeptSet's result).
LongResidResult computeLongitudinalResid(
    const std::string &phenoFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<std::string> &famIIDs,
    const std::unordered_set<std::string> &keptSubjects);

// Reorder a per-fit-subject vector onto a target IID order (e.g.
// SubjectData::usedIIDs()).  Throws if a target IID is absent from srcIIDs.
Eigen::VectorXd alignToUsed(
    const std::vector<std::string> &srcIIDs,
    const Eigen::VectorXd &srcVals,
    const std::vector<std::string> &targetIIDs);

// Row-wise counterpart of alignToUsed for a per-fit-subject design matrix.
Eigen::MatrixXd alignRowsToUsed(
    const std::vector<std::string> &srcIIDs,
    const Eigen::MatrixXd &srcRows,
    const std::vector<std::string> &targetIIDs);

} // namespace nsLongitudinal
