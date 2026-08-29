// spasqr.hpp — SPAsqr: SPA-squared (multi-tau quantile regression)
//
// Pure C++17 / Eigen / Boost port of ref_code/src/mtSPAsqr.h + SPAsqr.R.
//
// Three run modes, one settings bundle.  Which mode the CLI picks:
//
//   --spasqr-mode wald                → runSPAsqrWald
//   --spasqr-mode score + --pred-list → runSPAsqrLoco
//   --spasqr-mode score               → runSPAsqr
//
// Score-mode workflow:
//   1. Load phenotype/covariates, run QMME quantile regression per tau
//   2. Build smoothed residual matrix: R(i,t) = tau - Phi(-resid(i)/h)
//   3. Detect outliers per column (IQR-based, configurable)
//   4. Load sparse GRM and compute variance terms per column
//   5. Per-marker: score test per tau, then CCT across taus
//
// LOCO mode repeats 1-4 per chromosome against y - loco_chr.  Wald mode skips
// the null model entirely and refits the full model per (marker, tau).
#pragma once

#include "geno_factory/geno_data.hpp"
#include <string>
#include <vector>

// Everything the three entry points need, as named fields rather than a
// positional argument list.
//
// The previous signatures ran to 27 parameters each, with long runs of
// same-typed arguments -- three doubles for the QC cutoffs, three strings for
// keep / remove / transform -- so a transposed pair compiled silently and
// showed up only as a wrong answer.  Named initialisation makes that a typo the
// compiler can see, and it keeps the three modes' settings visibly one thing
// instead of three lists that have to be kept in sync by eye.
//
// Fields a given mode does not use are simply ignored: wald reads no GRM and no
// outlier settings, and plain score mode reads no --pred-list.
struct SPAsqrConfig {
    // ── Inputs ─────────────────────────────────────────────────────────
    std::string phenoFile;
    std::string covarFile;
    std::vector<std::string> phenoNames;
    std::vector<std::string> covarNames;
    std::vector<double> taus;
    GenoSpec geno;
    std::string spgrmGrabFile;      // --sp-grm-grab    (score/LOCO)
    std::string spgrmGctaFile;      // --sp-grm-plink2  (score/LOCO)
    std::string predListFile;       // --pred-list      (LOCO; empty ⇒ no LOCO)
    std::string keepFile;
    std::string removeFile;

    // ── Output ─────────────────────────────────────────────────────────
    std::string outPrefix;
    std::string compression;        // "" | "gz" | "zst"
    int compressionLevel = 3;

    // ── Model ──────────────────────────────────────────────────────────
    std::string phenoTransform = "int";   // "raw" | "int" | "standardize"
    double spaCutoff        = 2.0;        // |z| above which the SPA runs
    double outlierIqrRatio  = 1.5;
    double outlierAbsBound  = 0.55;
    double spasqrTol        = 1e-6;
    double spasqrH          = -1.0;       // -1 ⇒ IQR-based auto
    double spasqrHScale     = -1.0;       // -1 ⇒ 3 (score/LOCO), 5 (wald)

    // ── Marker QC ──────────────────────────────────────────────────────
    double missingCutoff = 0.1;
    double minMafCutoff  = 1e-5;
    double minMacCutoff  = 10.0;
    double hweCutoff     = 0.0;

    // ── Execution ──────────────────────────────────────────────────────
    int nthreads     = 1;
    int nSnpPerChunk = 8192;              // wald: 8192 sentinel ⇒ auto-shrink
};

// Score mode: one null-model fit, reused for every marker.
void runSPAsqr(const SPAsqrConfig &cfg);

// Score mode with a per-chromosome LOCO offset: the null model is refitted per
// chromosome against y - loco_chr, so the residuals (and their CGF tables) are
// chromosome-specific.  Requires cfg.predListFile.
void runSPAsqrLoco(const SPAsqrConfig &cfg);

// Wald mode: per-marker × per-τ full-model refit + M-estimation sandwich
// variance.  Emits one marker per line:
//   CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC LOG10P_HWE
//   LOG10P_CCT LOG10P_tau{val}... Z_tau{val}... BETA_tau{val}... SE_tau{val}...
// Response is Y_transformed, or Y_transformed - loco_chr when cfg.predListFile
// is set.  No GRM is used — point estimation, not a score test.
void runSPAsqrWald(const SPAsqrConfig &cfg);
