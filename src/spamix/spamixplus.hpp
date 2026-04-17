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
        const std::vector<uint32_t> &genoToFlat
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
        int nPC
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
        const std::vector<uint32_t> &genoToFlat
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
        int nPC
    );

    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 3;
    }

    std::string getHeaderColumns() const override {
        return "\tP\tZ\tBETA";
    }

    void prepareChunk(const std::vector<uint64_t> &gIndices) override;

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result
    ) override;

  private:
    double getMarkerPval(
        const Eigen::Ref<const Eigen::VectorXd> &GVec,
        double altFreq,
        double &zScore,
        double &outVarS
    )
    ;

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

// Pre-computed AF (non-null in pre-computed mode)
    const std::vector<AFModel> *m_afModels;
    const std::vector<uint32_t> *m_genoToFlat;

// On-the-fly AF (non-null in on-the-fly mode)
    const Eigen::MatrixXd *m_XtX_inv_Xt;
    const Eigen::VectorXd *m_sqrt_XtX_inv_diag;

// Per-thread scratch (mutable, freshly allocated in clone)
    Eigen::VectorXd m_AFVec;
    Eigen::VectorXd m_R_new; // only used when m_hasGRM
    Eigen::VectorXd m_mafOutlier;
    Eigen::VectorXd m_mafNonOutlier;

// Chunk state
    std::vector<uint64_t> m_chunkGenoIndices;
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
    const std::string &removeFile = {}

);
