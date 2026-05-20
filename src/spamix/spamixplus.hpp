// spamixplus.hpp — Unified SPAmix / SPAmixPlus implementation
//
// When a sparse GRM is provided the variance uses GRM-based covariance
// (SPAmixPlus mode); otherwise the diagonal  Σ resid²·2·AF·(1−AF) is
// used (SPAmix mode).  Passing an identity GRM to SPAmixPlus produces
// results identical to SPAmix.
//
// Output columns: [Pvalue, zScore, BETA]
#pragma once

#include "geno_factory/geno_data.hpp"
#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

#include "engine/marker.hpp"
#include "spamix/common.hpp"
#include "spamix/indiv_af.hpp"
#include "io/sparse_grm.hpp"

// ======================================================================
// SPAmixAFCache — per-thread shared cache of per-marker AFVec / W
//
// Phenotypes with identical missingness patterns (same maskIdx) within a
// fused-GEMM group reuse the same AFVec across all batch markers.  The
// engine processes one group's full window in a single processScoreBatch
// call sequence, so a per-thread cache keyed by (maskIdx, chunkIdxs) is
// sufficient: the first phenotype in a group fills the cache, the rest
// reuse it.  Recomputed when the chunkIdxs vector changes.
// ======================================================================

struct SPAmixAFCache {
    std::vector<int> lastChunkIdxs;
    Eigen::MatrixXd  afBatch;   // N × B (one AFVec column per marker)
    Eigen::MatrixXd  wBatch;    // N × B,  W = 2 · AF · (1 − AF)
};

// ======================================================================
// SPAmixPlusMethod — MethodBase implementation (unified SPAmix/SPAmixPlus)
// ======================================================================

class SPAmixPlusMethod : public MethodBase {
  public:
// ── With GRM (SPAmixPlus) ──────────────────────────────────────

// Pre-computed AF + GRM
    SPAmixPlusMethod(
        const Eigen::VectorXd &residuals,
        const Eigen::VectorXd &resid2,
        const Eigen::MatrixXd &onePlusPCs,
        const OutlierData &outlier,
        double spaCutoff,
        const SparseGRM &grm,
        const std::vector<AFModel> &afModels,
        const std::vector<uint32_t> &genoToFlat,
        int maskIdx = 0
    );

// On-the-fly AF + GRM
    SPAmixPlusMethod(
        const Eigen::VectorXd &residuals,
        const Eigen::VectorXd &resid2,
        const Eigen::MatrixXd &onePlusPCs,
        const OutlierData &outlier,
        double spaCutoff,
        const SparseGRM &grm,
        const Eigen::MatrixXd &XtX_inv_Xt,
        const Eigen::VectorXd &sqrt_XtX_inv_diag,
        int nPC,
        int maskIdx = 0
    );

// ── Without GRM (SPAmix) ──────────────────────────────────────

// Pre-computed AF, no GRM
    SPAmixPlusMethod(
        const Eigen::VectorXd &residuals,
        const Eigen::VectorXd &resid2,
        const Eigen::MatrixXd &onePlusPCs,
        const OutlierData &outlier,
        double spaCutoff,
        const std::vector<AFModel> &afModels,
        const std::vector<uint32_t> &genoToFlat,
        int maskIdx = 0
    );

// On-the-fly AF, no GRM
    SPAmixPlusMethod(
        const Eigen::VectorXd &residuals,
        const Eigen::VectorXd &resid2,
        const Eigen::MatrixXd &onePlusPCs,
        const OutlierData &outlier,
        double spaCutoff,
        const Eigen::MatrixXd &XtX_inv_Xt,
        const Eigen::VectorXd &sqrt_XtX_inv_diag,
        int nPC,
        int maskIdx = 0
    );

    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 4;
    }

    std::string getHeaderColumns() const override {
        return "\tP\tZ\tBETA\tSE";
    }

    void prepareChunk(const std::vector<uint64_t> &gIndices) override;

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result
    ) override;

    void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        const std::vector<int> &chunkIdxs,
        std::vector<std::vector<double> > &results
    ) override;

    int preferredBatchSize() const override {
        return 16;
    }

// ── Fused-GEMM interface ──────────────────────────────────────
// Enabled only in pre-computed AF mode.  On-the-fly AF needs a
// phenotype-local genotype vector to fit the AF model, which the
// fused path does not currently expose, so on-the-fly stays on
// the non-fused getResultBatch path.

    bool supportsFusedGemm() const override {
        return m_afModels != nullptr;
    }

    int fusedGemmColumns() const override {
        return 1;
    }

    void fillUnionResiduals(
        Eigen::Ref<Eigen::MatrixXd> dest,
        const std::vector<uint32_t> &unionToLocal
    ) const override;

    void fillResidualSums(double *dest) const override;

    void processScoreBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &scores,
        const double *gSums,
        const double *gSumSqs,
        uint32_t nUsed,
        const std::vector<double> &altFreqs,
        const std::vector<int> &chunkIdxs,
        std::vector<std::vector<double> > &results
    ) override;

  private:
    double markerPvalFromAF(
        const Eigen::Ref<const Eigen::VectorXd> &afVec,
        const Eigen::Ref<const Eigen::VectorXd> &wVec,
        double rawScore,
        double &zScore,
        double &outVarS
    );

    void fillAFVecForMarker(
        const Eigen::Ref<const Eigen::VectorXd> &GVec,
        double altFreq,
        int markerInChunkIdx,
        Eigen::Ref<Eigen::VectorXd> outAF
    );

    void ensureAFCacheFused(
        const std::vector<int> &chunkIdxs,
        const std::vector<double> &altFreqs
    );

// Read-only shared data (stable references — owner outlives all clones)
    const Eigen::VectorXd &m_resid;
    const Eigen::VectorXd &m_resid2;
    const Eigen::MatrixXd &m_onePlusPCs;
    const OutlierData &m_outlier;
    double m_spaCutoff;
    bool m_hasGRM;
    const SparseGRM *m_grm; // nullptr when !m_hasGRM
    int m_N;
    int m_nPC;
    int m_maskIdx; // identifies the dedup mask group (shared by phenotypes
                   // with identical missingness pattern)

    double m_residSum;

// Pre-computed AF (non-null in pre-computed mode)
    const std::vector<AFModel> *m_afModels;
    const std::vector<uint32_t> *m_genoToFlat;

// On-the-fly AF (non-null in on-the-fly mode)
    const Eigen::MatrixXd *m_XtX_inv_Xt;
    const Eigen::VectorXd *m_sqrt_XtX_inv_diag;

// Per-thread scratch (mutable, freshly allocated in clone)
    Eigen::VectorXd m_AFVec;
    Eigen::VectorXd m_WVec;  // W = 2 AF (1 - AF) scratch for batched path
    Eigen::VectorXd m_R_new; // only used when m_hasGRM
    Eigen::VectorXd m_mafOutlier;
    Eigen::VectorXd m_mafNonOutlier;

// Per-thread AFVec cache, shared between clones with the same maskIdx.
// Populated lazily by clone() through a thread_local registry.
    std::shared_ptr<SPAmixAFCache> m_afCache;

// Chunk state.  Pointer into the engine-owned vector (alive for the chunk).
    const std::vector<uint64_t> *m_chunkGenoIndices;
};

// ======================================================================
// Orchestration — unified for SPAmix and SPAmixPlus
// ======================================================================

// spgrmGrabFile + spgrmGctaFile both empty → SPAmix (diagonal variance).
void runSPAmixPlus(
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &pcColNames,
    const std::string &phenoFile,
    const std::string &covarFile,
    const GenoSpec &geno,
    const std::string &spgrmGrabFile,                // empty → no GRM
    const std::string &spgrmGctaFile,                // empty → no GRM
    const std::string &afFile,                       // empty → compute on-the-fly
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
    const std::string &regressionModelStr = {},
    const std::string &phenoNameSpec = {},
    const std::vector<std::string> &covarNames = {},
    bool saveResid = false,
    uint64_t seed = 0                                      // ordinal surrogate-residual RNG seed
);
