// sadge.hpp — SADGE family-based direct/indirect genetic-effect score test.
//
// SADGE decomposes a SNP's association into a direct genetic effect (G_D), an
// indirect / parental-nurture effect (G_I), and a standard total effect (G),
// using sibling + parent family structure.  It is implemented as a fused
// MethodBase that rides on multiPhenoEngine: the shared fused GEMM supplies
// S_G = Gᵀ·R (column 0) and the singleton partial Σ_singleton G·R (column 1),
// while SADGE imputes the parental genotype G_par via its OWN genotype cursor
// over the pedigree (siblings+parents) subject set and computes S_par itself.
//
// Two stages of the R reference are fused with no intermediate file:
//   stage 1 (impute G_par)  — example_sadge_zouyu/impute_func.R
//   stage 2 (score test)    — example_sadge_zouyu/score_func_vec.R
#pragma once

#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "sadge/precond.hpp"
#include "sadge/sample_info.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace sadge {

// Per non-singleton (fs1) union sibling: where to find its family members in
// SADGE's pedigree decode, so the imputation kernel can be evaluated.
struct NonSingKernel {
    uint32_t unionIdx;
    FamStruct fs1;
    int32_t pedSelf, pedCoSib, pedPar1, pedPar2;
};

// Read-only data shared by every method clone across all threads.
struct SADGEShared {
    std::shared_ptr<const GenoMeta> parGeno; // SADGE's own pedigree-set meta
    uint32_t nUnion = 0;
    uint32_t nPedUsed = 0;
    std::vector<NonSingKernel> nonSing;      // fs1-non-singleton union siblings
};

// One cached marker: the imputed parental genotype for the non-singleton
// siblings (compact), plus the per-marker coefficients.
struct MarkerCache {
    double mu = 0.0, H = 0.0, I = 0.0, J = 0.0;
    Eigen::VectorXd gImp; // length = nonSing.size()
};

// Per-thread mutable state shared by the K phenotype-clones running on that
// worker thread (so the chunk is decoded/imputed once, not K times).
struct SADGEWorkerCtx {
    std::unique_ptr<GenoCursor> cursor;
    std::vector<uint64_t> chunkGIdx;          // current chunk's marker indices
    int decodedUpTo = -1;                     // highest chunk-relative idx decoded
    std::unordered_map<int, MarkerCache> cache; // key = chunk-relative idx
    Eigen::VectorXd Gp;                        // pedigree decode scratch (nPedUsed)
    std::vector<uint32_t> missScratch;
};

class SADGEMethod : public MethodBase {
  public:
    SADGEMethod(
        std::shared_ptr<const SADGEShared> shared,
        std::shared_ptr<const PerPhenoData> pheno);

    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 13; // G_I/G_D/G × {Z,P,S,VarS} + SADGE_MAF
    }

    std::string getHeaderColumns() const override;

    void prepareChunk(const std::vector<uint64_t> &gIndices) override;

    // Fused-only method; getResultVec is never used by multiPhenoEngine.
    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result) override;

    bool supportsFusedGemm() const override {
        return true;
    }

    int fusedGemmColumns() const override {
        return 2; // col 0 = residUnion (S_G), col 1 = residSing (singleton Σ G·R)
    }

    void fillUnionResiduals(
        Eigen::Ref<Eigen::MatrixXd> dest,
        const std::vector<uint32_t> &unionToLocal) const override;

    void fillResidualSums(double *dest) const override;

    void processScoreBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &scores,
        const double *gSums,
        const double *gSumSqs,
        uint32_t nUsed,
        const std::vector<double> &altFreqs,
        const std::vector<int> &chunkIdxs,
        std::vector<std::vector<double> > &results) override;

  private:
    SADGEWorkerCtx &workerCtx();
    const MarkerCache &ensureMarker(SADGEWorkerCtx &c, int chunkRel);

    std::shared_ptr<const SADGEShared> m_shared;
    std::shared_ptr<const PerPhenoData> m_pheno;
};

} // namespace sadge

// ══════════════════════════════════════════════════════════════════════
// runSADGE — full workflow entry point (mirrors runSPAGRM)
// ══════════════════════════════════════════════════════════════════════

void runSADGE(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &sadgeFamFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {},
    const std::string &regressionModelStr = {},
    const std::string &phenoNameSpec = {},
    const std::string &covarFile = {},
    const std::vector<std::string> &covarNames = {},
    bool saveResid = false,
    uint64_t seed = 0);
