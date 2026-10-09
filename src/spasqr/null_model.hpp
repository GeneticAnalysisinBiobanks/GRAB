// null_model.hpp — Null-model setup shared by the three SPAsqr entry points
//
// runSPAsqr (score), runSPAsqrLoco (score + LOCO offset) and runSPAsqrWald all
// begin the same way: take the union subject set, split it per phenotype on
// that phenotype's non-missing rows, transform the response, and derive the
// smoothing bandwidth.  They diverge only afterwards, in what they fit and how
// often.  Keeping the common prologue here means "what the analysis set is"
// has one definition instead of three that drift apart.
//
// Module-internal: included only by spasqr.cpp and spasqr_wald.cpp.
#pragma once

#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace spasqr_null {

// ── Phenotype transform ────────────────────────────────────────────────

inline void applyPhenoTransform(
    Eigen::VectorXd &Y,
    const std::string &mode
) {
    if (mode == "raw") return;
    if (mode == "int") {
        Y = math::inverseRankNormal(Y);
        return;
    }
    if (mode == "standardize") {
        const Eigen::Index n = Y.size();
        if (n <= 1) return;
        const double mean = Y.mean();
        const double ssq  = (Y.array() - mean).square().sum();
        const double sd   = std::sqrt(ssq / static_cast<double>(n - 1));
        if (sd > 0.0) Y = (Y.array() - mean) / sd;
        return;
    }
    throw std::runtime_error("applyPhenoTransform: unknown mode '" + mode + "'");
}

// ── Smoothing bandwidth ────────────────────────────────────────────────

// h = IQR(Y) / scale, with R's type-7 quantile interpolation.
//
// A non-positive IQR means the response is degenerate (over half its mass on a
// single value), which would leave the smoothed check function with no scale at
// all; the fallback is the usual n^(-2/5)-flavoured rule, floored at 0.05.
inline double iqrBandwidth(
    const Eigen::VectorXd &Y,
    double scale,
    int nCov
) {
    const Eigen::Index n = Y.size();
    if (n < 4) return 1.0;

    std::vector<double> v(static_cast<size_t>(n));
    Eigen::VectorXd::Map(v.data(), n) = Y;
    std::sort(v.begin(), v.end());
    auto quantile = [&](double prob) {
        const double idx = prob * (static_cast<double>(n) - 1);
        const Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
        const Eigen::Index hi = std::min(lo + 1, n - 1);
        const double frac = idx - lo;
        return v[lo] * (1.0 - frac) + v[hi] * frac;
    };

    const double h = (quantile(0.75) - quantile(0.25)) / scale;
    if (h > 0.0) return h;
    return std::max(std::pow((std::log(static_cast<double>(n)) + nCov) /
                             static_cast<double>(n), 0.4),
                    0.05);
}

// ── Per-phenotype analysis set ─────────────────────────────────────────

// One phenotype's slice of the union: the subjects it is non-missing on, its
// transformed response, and the matching covariate rows.  `h` is left to the
// caller — score mode derives it once from Y, LOCO mode per chromosome from
// the LOCO-adjusted response, and Wald mode does both depending on --pred-list.
struct PhenoWork {
    std::vector<uint32_t> unionToLocal; // size nUnion; UINT32_MAX = absent
    uint32_t nk = 0;                    // non-missing count
    Eigen::VectorXd Y;                  // nk, transformed
    Eigen::MatrixXd X;                  // nk × nCov
    double h = 0.0;                     // smoothing bandwidth (caller-filled)
};

// Split the union into per-phenotype workspaces.  Union membership is decided
// once (genotype ∩ GRM ∩ keep/remove); NA filtering is per phenotype, so each
// phenotype gets its own contiguous re-indexing of the union.
//
// The transform runs on the per-phenotype non-missing scope, and everything
// downstream — bandwidth, QMME fit, and in LOCO mode the per-chromosome y_adj —
// operates on the transformed Y.  A LOCO PRS therefore has to be on the same
// scale as the transform the user asked for.
//
// `label` only selects the prefix of the error message.
inline std::vector<PhenoWork> buildPhenoWorkspaces(
    SubjectData &sd,
    const Eigen::MatrixXd &unionX,
    const std::vector<std::string> &phenoNames,
    const std::string &phenoTransform,
    const char *label
) {
    const uint32_t nUnion = sd.nUsed();
    const int nCov = static_cast<int>(unionX.cols());
    const int K = static_cast<int>(phenoNames.size());
    std::vector<PhenoWork> pw(K);

    for (int k = 0; k < K; ++k) {
        Eigen::VectorXd fullY = sd.getColumn(phenoNames[k]);
        pw[k].unionToLocal.assign(nUnion, UINT32_MAX);
        uint32_t localIdx = 0;
        for (uint32_t i = 0; i < nUnion; ++i)
            if (!std::isnan(fullY[i])) pw[k].unionToLocal[i] = localIdx++;

        pw[k].nk = localIdx;
        if (pw[k].nk == 0)
            throw std::runtime_error(std::string(label) + ": phenotype '" +
                                     phenoNames[k] + "' has no non-missing subjects");

        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        pw[k].Y.resize(Nk);
        pw[k].X.resize(Nk, nCov);
        for (uint32_t i = 0; i < nUnion; ++i) {
            const uint32_t li = pw[k].unionToLocal[i];
            if (li == UINT32_MAX) continue;
            pw[k].Y[li] = fullY[i];
            if (nCov > 0) pw[k].X.row(li) = unionX.row(i);
        }

        applyPhenoTransform(pw[k].Y, phenoTransform);
    }
    return pw;
}

// ── Reporting ──────────────────────────────────────────────────────────

// "<pheno>  N  <bandwidth>" table on stderr.  LOCO mode refits the bandwidth
// per chromosome, so it labels this column as the global reference value.
inline void logPhenoTable(
    const std::vector<std::string> &phenoNames,
    const std::vector<PhenoWork> &pw,
    const char *bwLabel
) {
    infoMsg("Sample size and smooth bandwidth per phenotype:");
    size_t nameW = 4;
    for (const auto &n : phenoNames) nameW = std::max(nameW, n.size());
    nameW += 2;

    char row[256];
    std::snprintf(row, sizeof(row), "    %-*s %10s %14s\n",
                  static_cast<int>(nameW), "", "N", bwLabel);
    fprintf(stderr, "%s", row);
    for (size_t k = 0; k < phenoNames.size(); ++k) {
        std::snprintf(row, sizeof(row), "    %-*s %10u %14.6f\n",
                      static_cast<int>(nameW), phenoNames[k].c_str(), pw[k].nk, pw[k].h);
        fprintf(stderr, "%s", row);
    }
}

// Per-phenotype × per-tau outlier fractions, one row per phenotype.  Score mode
// prints this once, LOCO mode once per chromosome, hence `title`.
inline void logOutlierTable(
    const std::string &title,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &tauLabels,
    const std::vector<std::vector<double> > &ratios,
    double outlierIqrRatio,
    double outlierAbsBound
) {
    size_t nameW = 4;
    for (const auto &n : phenoNames) nameW = std::max(nameW, n.size());
    nameW += 2;

    infoMsg("%s (IQR=%.2f, bound=%.2f):", title.c_str(),
            outlierIqrRatio, outlierAbsBound);

    std::ostringstream hdr;
    hdr << "    " << std::setw(static_cast<int>(nameW)) << std::left << "";
    for (const auto &tl : tauLabels)
        hdr << std::setw(10) << std::right << tl;
    fprintf(stderr, "%s\n", hdr.str().c_str());

    for (size_t k = 0; k < phenoNames.size(); ++k) {
        std::ostringstream row;
        row << "    " << std::setw(static_cast<int>(nameW)) << std::left << phenoNames[k];
        for (double r : ratios[k])
            row << std::setw(10) << std::right << std::fixed << std::setprecision(4) << r;
        fprintf(stderr, "%s\n", row.str().c_str());
    }
}

// ── LOCO offset ────────────────────────────────────────────────────────

// Scatter one chromosome's LOCO PGS from union space into this phenotype's
// dense space.  Score and Wald mode both need exactly this, including the
// hard failure below: LDAK / Regenie Step 1 is supposed to emit a score for
// every analysed subject, and the parsers leave NaN where one is missing.  A
// NaN here would propagate silently through y_adj into the QMME fit, so it is
// worth a stop with an actionable message instead.  `label` prefixes it.
inline Eigen::VectorXd locoDense(
    const Eigen::VectorXd &unionScores,
    const PhenoWork &pk,
    const std::string &phenoName,
    const std::string &chr,
    const char *label
) {
    const Eigen::Index Nk = static_cast<Eigen::Index>(pk.nk);
    Eigen::VectorXd dense(Nk);
    for (uint32_t i = 0; i < pk.unionToLocal.size(); ++i) {
        const uint32_t li = pk.unionToLocal[i];
        if (li != UINT32_MAX) dense[li] = unionScores[i];
    }
    if (!dense.allFinite()) {
        const Eigen::Index nBad = Nk - dense.array().isFinite().count();
        throw std::runtime_error(
            std::string(label) + ": LOCO file for phenotype '" + phenoName +
            "' chr " + chr + " is missing " + std::to_string(nBad) +
            " subject(s) that have non-missing Y. The LOCO PGS file must "
            "contain every subject in the --pheno analysis set. Re-run "
            "LDAK / Regenie Step 1 on the same sample set, or remove those "
            "subjects from --pheno.");
    }
    return dense;
}

inline std::vector<std::string> makeTauLabels(const std::vector<double> &taus) {
    std::vector<std::string> labels;
    labels.reserve(taus.size());
    for (double tau : taus) {
        std::ostringstream oss;
        oss << "tau" << tau;
        labels.push_back(oss.str());
    }
    return labels;
}

} // namespace spasqr_null
