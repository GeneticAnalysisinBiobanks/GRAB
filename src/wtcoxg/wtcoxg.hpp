// wtcoxg.hpp — WtCoxG: Weighted Cox-type G-test (Li et al., 2025)
//
// Phases:
//   1. Data loading & marker matching  (resid-file, bfile, ref-af-file)
//   2. Batch-effect testing & parameter estimation (TPR, sigma2, w.ext, var ratios)
//   3. Marker-level SPA tests via markerEngine (MethodBase interface)
#pragma once

#include "geno_factory/geno_data.hpp"
#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <unordered_map>
#include <vector>

#include "engine/marker.hpp"

class GenoMeta;
class SparseGRM;

// ======================================================================
// Per-marker reference info  (populated during Phase 2)
// ======================================================================

struct WtCoxGRefInfo {
    double AF_ref = std::numeric_limits<double>::quiet_NaN();
    double obs_ct = std::numeric_limits<double>::quiet_NaN(); // allele count
    double TPR = std::numeric_limits<double>::quiet_NaN();
    double sigma2 = std::numeric_limits<double>::quiet_NaN();
    double pvalue_bat = std::numeric_limits<double>::quiet_NaN();
    double w_ext = std::numeric_limits<double>::quiet_NaN();
    double var_ratio_w0 = 1.0;
    double var_ratio_int = 1.0;
    double var_ratio_ext = 1.0;
};

// ======================================================================
// WtCoxGMethod — implements MethodBase for marker engine
// ======================================================================

class WtCoxGMethod : public MethodBase {
  public:
// Construct from pre-computed null-model quantities.
//   R:           martingale residuals (nSubj)
//   w:           sampling weights     (nSubj)
//   cutoff:      batch-effect p-value threshold (e.g. 0.05)
//   SPA_Cutoff:  z-score threshold to switch from normal to SPA
//   refMap:      genoIndex → WtCoxGRefInfo (from Phase 2)
    WtCoxGMethod(
        Eigen::VectorXd R,
        Eigen::VectorXd w,
        double cutoff,
        double SPA_Cutoff,
        std::shared_ptr<const std::unordered_map<uint64_t, WtCoxGRefInfo> > refMap
    );

// ---- MethodBase interface ----
    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 5;
    }

    std::string getHeaderColumns() const override;

    void prepareChunk(const std::vector<uint64_t> &gIndices) override;

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result
    ) override;

// For LEAF: compute ext/noext results with raw scores for meta-analysis.
    struct DualResult {
        double p_ext, p_noext;
        double score_ext, score_noext;
    };

    DualResult computeDual(
        Eigen::Ref<Eigen::VectorXd> GVec,
        int markerInChunkIdx
    );

// Access per-chunk ref info (for LEAF meta-analysis).
    const WtCoxGRefInfo &chunkRefInfoAt(int idx) const {
        return m_chunkRefInfo[idx];
    }

  private:
    struct WtResult {
        double pval;
        double score;
        double zscore;
    };

// Core SPA test for one marker with external adjustment.
    WtResult wtCoxGTest(
        const Eigen::Ref<const Eigen::VectorXd> &g,
        double p_bat,
        double TPR,
        double sigma2,
        double b,
        double var_ratio_int,
        double var_ratio_w0,
        double var_ratio_w1,
        double var_ratio0,
        double var_ratio1,
        double mu_ext,
        double obs_ct,
        double p_cut
    ) const;

// Members (const after construction, except per-chunk scratch)
    Eigen::VectorXd m_R;
    Eigen::VectorXd m_w;
    Eigen::VectorXd m_w1; // w / (2 * sum(w))
    double m_meanR;
    double m_sumR;
    double m_cutoff;
    double m_SPA_Cutoff;
    std::shared_ptr<const std::unordered_map<uint64_t, WtCoxGRefInfo> > m_refMap;

// Per-chunk scratch (rebuilt in prepareChunk)
    std::vector<WtCoxGRefInfo> m_chunkRefInfo;
    std::array<double, 2> m_scoreArr;
    std::array<double, 2> m_zScoreArr;
};

// ======================================================================
// Phase 1 — Data loading & marker matching
// ======================================================================

// Parsed plink2 .afreq record (one per variant)
// Columns: #CHROM  ID  REF  ALT  [PROVISIONAL_REF?]  ALT_FREQS  OBS_CT
struct RefAfRecord {
    std::string chrom;
    std::string id;
    std::string ref_allele; // REF column (plink2 reference allele)
    std::string alt_allele; // ALT column (plink2 alternate allele)
    double alt_freq;        // ALT_FREQS column
    double obs_ct;          // OBS_CT column — total allele number
};

// Parse a plink2 --freq .afreq file.  Header line (starting with '#')
// is used to detect column positions for CHROM, ID, REF, ALT, ALT_FREQS, OBS_CT.
// Two-column numeric fallback: if there is no header and each line has exactly
// two numeric values, they are treated as (ALT_FREQS, OBS_CT) in .bim order.
// When isNumericFallback is non-null it is set to true in that case.
std::vector<RefAfRecord> loadRefAfFile(
    const std::string &filename,
    bool *isNumericFallback = nullptr
);

// Match bim markers to ref-af records by (chr, bp, a1, a2) — exact only.
// Returns a map: genoIndex → {AF_ref, obs_ct, ...} for matched markers.
// mu0, mu1, n0, n1 per marker are computed from the genotype data.
struct MatchedMarkerInfo {
    uint64_t genoIndex;
    double AF_ref;
    double obs_ct; // total allele number (used directly in formulas)
    double mu0;    // control allele freq
    double mu1;    // case allele freq
    double n0;     // effective control count
    double n1;     // effective case count
    double mu_int; // internal MAF = 0.5*mu0 + 0.5*mu1 (folded ≤0.5)
};

// Sequential assignment for the two-column numeric fallback: rows are
// assumed to be in .bim order, AF_ref = ALT_FREQS directly (col 5 freq).
// Throws if row count != bim marker count.
std::vector<MatchedMarkerInfo> matchMarkersNumeric(
    const GenoMeta &plinkData,
    const std::vector<RefAfRecord> &refAf
);

// Match bim markers to ref-af records by (CHROM, ID) with allele orientation:
//   If afreq (ALT,REF) matches bim (col5,col6) -> AF_ref = ALT_FREQS
//   If afreq (REF,ALT) matches bim (col5,col6) -> AF_ref = 1-ALT_FREQS
//   Otherwise the marker is dropped.
// obs_ct is passed through directly (allele count).
std::vector<MatchedMarkerInfo> matchMarkers(
    const GenoMeta &plinkData,
    const std::vector<RefAfRecord> &refAf
);

// Scan genotypes to compute per-marker case/control allele frequencies.
// Populates mu0, mu1, n0, n1, mu_int in each MatchedMarkerInfo.
void computeMarkerStats(
    std::vector<MatchedMarkerInfo> &matched,
    const GenoMeta &plinkData,
    const Eigen::VectorXd &indicator
);

// ======================================================================
// Phase 2 — Batch-effect testing and parameter estimation
// ======================================================================

// Full batch-effect QC pipeline:
//   1. Compute batch p-values per marker
//   2. Compute variance ratios from sparse GRM (if provided)
//   3. Estimate TPR, sigma2 per MAF group (Nelder-Mead)
//   4. Compute optimal w.ext per MAF group (Brent 1D)
//   5. Compute var_ratio_ext per MAF group from sparse GRM (if provided)
//
// Populates refInfo for each matched marker and returns a shared map.
std::shared_ptr<std::unordered_map<uint64_t, WtCoxGRefInfo> >testBatchEffects(
    const std::vector<MatchedMarkerInfo> &matched,
    const Eigen::VectorXd &residuals,
    const Eigen::VectorXd &weights,
    const Eigen::VectorXd &indicator,
    const SparseGRM *grm,
    double refPrevalence,
    double cutoff
);

// ======================================================================
// Top-level orchestration — called from main()
// ======================================================================

// --pheno path: compute regression residuals internally
void runWtCoxGPheno(
    const std::string &phenoFile,
    const std::string &covarFile,                               // empty = no separate covar file
    const std::vector<std::string> &covarNames,                 // empty = intercept only
    const std::vector<std::string> &phenoNames,                 // selected phenotype columns
    const GenoSpec &geno,
    const std::string &refAfFile,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {}

);

// Multi-phenotype entry point: loads data/ref-AF/GRM once, parallelizes
// P null-model fits and batch-effect tests with min(nthreads, P) workers,
// then runs a single multiPhenoEngine call for all P phenotypes.
// Each phenoSpec is "TIME:EVENT" (survival) or "COLNAME" (binary).
void runWtCoxG(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    const std::vector<std::string> &phenoSpecs,
    const GenoSpec &geno,
    const std::string &refAfFile,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {}
);
