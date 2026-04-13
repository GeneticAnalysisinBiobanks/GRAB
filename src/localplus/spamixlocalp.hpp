// spamixlocalp.hpp — SPAmixLocalPlus: local-ancestry-specific GWAS
//
// Phase 1: Phi estimation — streaming ancestry-specific kinship from admix .abed
// Phase 2: Per-ancestry GWAS — score test with SPA tail, CCT meta-analysis
//
// Uses .abed binary format for local ancestry dosage/hapcount data,
// and SparseGRM for related-pair structure.
#pragma once

#include "localplus/abed_io.hpp"
#include "spamix/common.hpp"
#include "io/sparse_grm.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <string>
#include <vector>

// ======================================================================
// PhiMatrices — ancestry-specific kinship for 4 haplotype scenarios
// ======================================================================

// Sparse COO storage for phi values, indexed into the "used subject" order.
struct PhiEntry {
    uint32_t i; // row index (used-subject space)
    uint32_t j; // col index (used-subject space)
    double value;
};

// Four scenarios based on (h_i, h_j) haplotype counts from a specific ancestry:
//   A: h_i=2, h_j=2   multiplier = 4 * q(1-q)
//   B: h_i=2, h_j=1   multiplier = 2 * q(1-q)
//   C: h_i=1, h_j=2   multiplier = 2 * q(1-q)
//   D: h_i=1, h_j=1   multiplier = 1 * q(1-q)
struct PhiMatrices {
    std::vector<PhiEntry> A; // (2,2) pairs
    std::vector<PhiEntry> B; // (2,1) pairs
    std::vector<PhiEntry> C; // (1,2) pairs
    std::vector<PhiEntry> D; // (1,1) pairs
};

// Pre-computed phi entries with R[i]*R[j]*phi*multiplier baked in.
// Eliminates random R[] lookups in the per-marker variance hot path.
//
// SoA (Structure-of-Arrays) layout for cache-friendly sequential scan
// and AVX2 vectorized processing.  All four scenarios (A/B/C/D) are
// merged into a single flat array so the hot loop makes ONE pass.
struct RprodSoA {
    std::vector<uint32_t> idx_i;    // subject index i
    std::vector<uint32_t> idx_j;    // subject index j
    std::vector<double> rprod;      // multiplier * phi * R[i] * R[j]
    std::vector<uint8_t> target_hi; // required hInt[i] value (1 or 2)
    std::vector<uint8_t> target_hj; // required hInt[j] value (1 or 2)

    size_t size() const {
        return rprod.size();
    }

    void reserve(size_t n) {
        idx_i.reserve(n);
        idx_j.reserve(n);
        rprod.reserve(n);
        target_hi.reserve(n);
        target_hj.reserve(n);
    }

    void push_back(
        uint32_t i,
        uint32_t j,
        double rp,
        uint8_t th_i,
        uint8_t th_j
    ) {
        idx_i.push_back(i);
        idx_j.push_back(j);
        rprod.push_back(rp);
        target_hi.push_back(th_i);
        target_hj.push_back(th_j);
    }

};

// Build unified SoA from PhiMatrices and residual vector (once per phenotype).
RprodSoA buildRprodSoA(
    const PhiMatrices &phi,
    const Eigen::VectorXd &R
);

// Compute off-diagonal variance from SoA phi entries and hapcount array.
// Uses AVX2 when available, scalar fallback otherwise.
double computeVarOffSoA(
    const RprodSoA &rp,
    const uint32_t *hInt,
    uint32_t nUsed
);

// Batch size for amortizing phi scan across multiple markers.
// Scanning phi is memory-bandwidth-bound; processing B markers per scan
// reduces DRAM traffic by factor B.
static constexpr int PHI_BATCH = 8;

// Batch off-diagonal variance: scans phi entries ONCE for PHI_BATCH markers.
//   hIntSM:   subject-major uint32 hapcount [nUsed * PHI_BATCH]
//             layout: hIntSM[subj * PHI_BATCH + batchIdx]
//   batchLen: actual marker count in this batch (≤ PHI_BATCH)
//   varOff:   output array [batchLen]
void computeVarOffSoABatch(
    const RprodSoA &rp,
    const uint32_t *hIntSM,
    uint32_t nUsed,
    int batchLen,
    double *varOff
);

// ======================================================================
// Phi estimation
// ======================================================================

// Estimate phi matrices for one ancestry from admix binary data + sparse GRM.
//
// For each GRM pair (i,j), accumulates:
//   phi_{ij}^{scenario} = (1/M) sum_s { (g_is - h_is*q_s)(g_js - h_js*q_s) / (h_is*h_js*q_s*(1-q_s)) }
// where the sum is over SNPs where both individuals match the scenario's (h_i, h_j) pattern.
//
//   admixData:  the admix binary data backend
//   grm:        sparse GRM (determines which pairs to compute)
//   ancIdx:     which ancestry to estimate (0-based)
//   MAF cutoff is hardcoded to 0.01.
PhiMatrices estimatePhiOneAncestry(
    const AdmixData &admixData,
    const SparseGRM &grm,
    int ancIdx,
    int nthreads = 1
);

// ======================================================================
// SPAmixLocalPlus variance computation
// ======================================================================

// Compute GRM-based variance for the score statistic, using phi matrices.
//   Var(S) = sum_scenario sum_{(i,j) in phi_scenario} w_scenario * q(1-q) * phi_ij * R_i * R_j
//            + sum_i R_i^2 * h_i * q(1-q)
//
// where w_scenario is the haplotype-count multiplier (4, 2, 2, or 1).
double computePhiVariance(
    const Eigen::VectorXd &R,
    const Eigen::VectorXd &hapcount,
    double q,
    const PhiMatrices &phi
);

// ======================================================================
// SPA p-value with outlier split
// ======================================================================

// Compute SPA p-value for the local-ancestry score test.
//   S:        score statistic = sum(dosage * R)
//   sMean:    pre-computed mean of S = q * hapcount.dot(R)
//   varDiag:  diagonal-only variance = q(1-q) * sum(R_i^2 * h_i)
//   R:        residual vector
//   hapcount: hapcount vector for this ancestry at this marker
//   q:        allele frequency
//   varS:     pre-computed variance (from computePhiVariance)
//   outlier:  outlier positions
//   spaCutoff: threshold for switching from normal to SPA
//
// Returns: {pval_spa, pval_normal}
std::pair<double, double> spaLocalPval(
    double S,
    double sMean,
    double varDiag,
    const Eigen::VectorXd &R,
    const Eigen::VectorXd &hapcount,
    double q,
    double varS,
    const OutlierData &outlier,
    double spaCutoff
);

// ======================================================================
// Full pipeline entry points
// ======================================================================

// Estimate phi matrices for all ancestries and write single wide file.
//   admixPrefix:  prefix for .abed/.bim/.fam
//   grmGrabFile / grmGctaFile: sparse GRM (exactly one non-empty)
//   phiOutputFile: output path for wide phi file
//   extractFile / excludeFile: for marker filtering
void runPhiEstimation(
    const std::string &admixPrefix,
    const std::string &grmGrabFile,
    const std::string &grmGctaFile,
    const std::string &phiOutputFile,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    const std::string &extractFile = {},
    const std::string &excludeFile = {},
    int nthreads = 1
);

// Run per-ancestry GWAS with SPAmixLocalPlus.
//   phenoFile:   phenotype file (columns selected by residNames)
//   residNames:  column names to use as residuals from phenoFile
//   admixPrefix:  prefix for .abed/.bim/.fam
//   admixPhiFile: pre-computed wide phi file
//   outPrefix:   output prefix for per-phenotype GWAS results
//   spaCutoff, outlierRatio, nthread, nSnpPerChunk: analysis params
void runSPAmixLocalPlus(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &admixPrefix,
    const std::string &admixPhiFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    double outlierRatio,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    const std::string &extractFile = {},
    const std::string &excludeFile = {}

);
