// spagrm.hpp — Saddlepoint approximation with GRM (pure C++17 / Eigen / Boost)
//
// Translated from mtSPAGRM.h / mtSPAGRM.cpp (RcppArmadillo) to Eigen.
// Provides the nsSPAGRM namespace (mgf, fastGetRoot, getProbSpa) and
// the SPAGRMClass marker-level evaluator used by SPAsqr and SPAGRM.
#pragma once

#include <Eigen/Dense>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "engine/marker.hpp"          // MethodBase
#include "geno_factory/geno_data.hpp" // GenoSpec

namespace nsSPAGRM {

// Per-marker updated data for three-or-more-subject families.
struct UpdatedThreeSubj {
    std::vector<double> stand_S;
    std::vector<double> arr_prob;
};

// Family data for one tau (or one residual column).
struct FamilyData {
    Eigen::VectorXd resid_unrelated_outliers;
    std::vector<std::array<double, 2> > twoSubj_resid;
    std::vector<std::vector<double> > twoSubj_rho;
    std::vector<std::vector<double> > threeSubj_standS;
    std::vector<Eigen::MatrixXd> threeSubj_CLT;
};

// Total number of elements in the mgf output vectors.
inline size_t mgfOutputSize(
    size_t n_unrelated,
    const std::vector<std::vector<double> > &TwoSubj_rho,
    size_t nThreeSubj
)
{
    size_t sz = n_unrelated;
    for (const auto &r : TwoSubj_rho)
        sz += r.size();
    sz += nThreeSubj;
    return sz;
}

// Pre-allocated scratch workspace (one per SPAGRMClass instance).
struct MgfWorkspace {
    Eigen::VectorXd MGF0, MGF1, MGF2;
    Eigen::VectorXd temp;
    // Unrelated-outlier intermediates
    Eigen::VectorXd ul_lambda, ul_alpha, ul_alpha_1, ul_alpha_2;

    MgfWorkspace() = default;
    MgfWorkspace(
        Eigen::Index mgfSz,
        Eigen::Index n_unrelated
    )
        : MGF0(mgfSz),
          MGF1(mgfSz),
          MGF2(mgfSz),
          temp(mgfSz),
          ul_lambda(n_unrelated),
          ul_alpha(n_unrelated),
          ul_alpha_1(n_unrelated),
          ul_alpha_2(n_unrelated)
    {
    }

};

// MGF and its first two derivatives → ws.MGF0/MGF1/MGF2.
void mgf(
    double t,
    const Eigen::VectorXd &resid_unrelated_outliers,
    const std::vector<std::array<double, 2> > &TwoSubj_resid,
    const std::vector<std::vector<double> > &TwoSubj_rho,
    const std::vector<UpdatedThreeSubj> &threeSubj,
    double MAF,
    MgfWorkspace &ws
);

// Newton-Raphson root finder for the CGF equation.
double fastGetRoot(
    const Eigen::VectorXd &resid_unrelated_outliers,
    const std::vector<std::array<double, 2> > &TwoSubj_resid,
    const std::vector<std::vector<double> > &TwoSubj_rho,
    const std::vector<UpdatedThreeSubj> &threeSubj,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score,
    double MAF,
    double init_t,
    double tol,
    MgfWorkspace &ws,
    int maxiter = 50
);

// Saddlepoint approximation tail probability.
double getProbSpa(
    const Eigen::VectorXd &resid_unrelated_outliers,
    const std::vector<std::array<double, 2> > &TwoSubj_resid,
    const std::vector<std::vector<double> > &TwoSubj_rho,
    const std::vector<UpdatedThreeSubj> &threeSubj,
    double sum_R_nonOutlier,
    double R_GRM_R_nonOutlier,
    double Score,
    double MAF,
    bool lower_tail,
    double zeta,
    double tol,
    MgfWorkspace &ws
);

} // namespace nsSPAGRM

// ══════════════════════════════════════════════════════════════════════
// SPAGRMClass — per-column (per-tau) marker evaluator
// ══════════════════════════════════════════════════════════════════════

class SPAGRMClass {
  public:
    SPAGRMClass(
        Eigen::VectorXd resid,
        double sum_R_nonOutlier,
        double R_GRM_R_nonOutlier,
        double R_GRM_R_TwoSubjOutlier,
        double R_GRM_R,
        std::vector<double> MAF_interval,
        nsSPAGRM::FamilyData fam,
        double SPA_Cutoff,
        double zeta,
        double tol
    );

    // Deep copy (including scratch state) for per-thread isolation.
    SPAGRMClass(const SPAGRMClass &o);
    SPAGRMClass &operator=(const SPAGRMClass &) = delete;

    double getMarkerPval(
        const Eigen::VectorXd &GVec,
        double altFreq,
        double &zScore,
        double gMean = std::numeric_limits<double>::quiet_NaN()
    );

    // Fast path: caller pre-computes Score = GVec.dot(m_resid) - gMean * m_resid_sum
    // via a fused multi-tau matrix multiply, avoiding redundant GVec reads.
    double getMarkerPvalFromScore(
        double Score,
        double altFreq,
        double &zScore
    );

    const Eigen::VectorXd &resid() const {
        return m_shared->resid;
    }

    double residSum() const {
        return m_shared->resid_sum;
    }

    // Read-only data shared across clones via shared_ptr.
    struct SharedData {
        Eigen::VectorXd resid;
        Eigen::VectorXd resid_unrelated_outliers;
        double sum_unrelated_outliers2;
        double sum_R_nonOutlier;
        double R_GRM_R_nonOutlier;
        double R_GRM_R_TwoSubjOutlier;
        double R_GRM_R;
        double resid_sum;
        std::vector<double> MAF_interval;
        std::vector<std::array<double, 2> > TwoSubj_resid_list;
        std::vector<std::vector<double> > TwoSubj_rho_list;
        std::vector<std::vector<double> > ThreeSubj_standS_list;
        std::vector<Eigen::MatrixXd> ThreeSubj_CLT_list;
        double SPA_Cutoff;
        double zeta;
        double tol;
    };

    const std::shared_ptr<const SharedData> &sharedData() const {
        return m_shared;
    }

  private:
    std::shared_ptr<const SharedData> m_shared;

    // Per-thread mutable scratch (rebuilt on copy).
    nsSPAGRM::MgfWorkspace m_workspace;
    std::vector<nsSPAGRM::UpdatedThreeSubj> m_threeSubj_scratch;

    void rebuildScratch();

};

// ══════════════════════════════════════════════════════════════════════
// SPAGRMMethod — MethodBase adapter wrapping a single SPAGRMClass
// ══════════════════════════════════════════════════════════════════════

class SPAGRMMethod : public MethodBase {
  public:
    explicit SPAGRMMethod(SPAGRMClass spagrm)
        : m_spagrm(std::move(spagrm))
    {
    }

    std::unique_ptr<MethodBase> clone() const override {
        return std::make_unique<SPAGRMMethod>(*this);
    }

    int resultSize() const override {
        return 2;
    }

    std::string getHeaderColumns() const override {
        return "\tP\tZ";
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int /*markerInChunkIdx*/,
        std::vector<double> &result
    ) override {
        result.clear();
        double z;
        double p = m_spagrm.getMarkerPval(GVec, altFreq, z);
        result.push_back(p);
        result.push_back(z);
    }

    // Batch analysis: fuse B dot products into one matrix-vector multiply.
    // Score[b] = GBatch.col(b).dot(resid) - gMean[b] * resid_sum
    // becomes: scores = GBatch^T * resid  (B×1, one Eigen matvec)
    void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        std::vector<std::vector<double> > &results
    ) override {
        const int B = static_cast<int>(GBatch.cols());
        results.resize(B);

        // Fused: B dot products in one matrix-vector multiply.
        // GBatch is N × B, resid is N × 1 → scores is B × 1
        const Eigen::VectorXd &resid = m_spagrm.resid();
        const double resid_sum = m_spagrm.residSum();

        Eigen::VectorXd scores;
        scores.noalias() = GBatch.transpose() * resid;

        // gMeans (length B) and mean adjustment
        const Eigen::VectorXd gMeans = GBatch.colwise().mean();
        scores.array() -= gMeans.array() * resid_sum;

        // Per-marker SPA (not batchable)
        for (int b = 0; b < B; ++b) {
            results[b].clear();
            double z;
            double p = m_spagrm.getMarkerPvalFromScore(scores[b], altFreqs[b], z);
            results[b].push_back(p);
            results[b].push_back(z);
        }
    }

    int preferredBatchSize() const override {
        return 8;
    }

    // ── Fused union-level GEMM interface ───────────────────────────────
    bool supportsFusedGemm() const override {
        return true;
    }

    int fusedGemmColumns() const override {
        return 1;
    }

    void fillUnionResiduals(
        Eigen::Ref<Eigen::MatrixXd> dest,
        const std::vector<uint32_t> &unionToLocal
    ) const override {
        // dest is pre-zeroed, N_union × 1.
        const Eigen::VectorXd &r = m_spagrm.resid();
        const uint32_t nUnion = static_cast<uint32_t>(unionToLocal.size());
        for (uint32_t i = 0; i < nUnion; ++i) {
            const uint32_t li = unionToLocal[i];
            if (li != UINT32_MAX)
                dest(i, 0) = r[li];
        }
    }

    void fillResidualSums(double *dest) const override {
        dest[0] = m_spagrm.residSum();
    }

    void processScoreBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &scores,
        const double *gSums,
        uint32_t nUsed,
        const std::vector<double> &altFreqs,
        std::vector<std::vector<double> > &results
    ) override {
        const int B = static_cast<int>(scores.cols());
        results.resize(B);
        const double residSum = m_spagrm.residSum();
        const double invN = 1.0 / static_cast<double>(nUsed);

        for (int b = 0; b < B; ++b) {
            const double gMean = gSums[b] * invN;
            const double centeredScore = scores(0, b) - gMean * residSum;
            results[b].clear();
            double z;
            double p = m_spagrm.getMarkerPvalFromScore(centeredScore, altFreqs[b], z);
            results[b].push_back(p);
            results[b].push_back(z);
        }
    }

  private:
    SPAGRMClass m_spagrm;
};

// ══════════════════════════════════════════════════════════════════════
// runSPAGRM — full workflow entry point
// ══════════════════════════════════════════════════════════════════════

void runSPAGRM(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
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
