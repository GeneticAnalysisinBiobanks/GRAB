// marker.hpp — Thread-pool marker-level association engine (pure C++17 / Eigen)
#pragma once

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

class GenoMeta;

// ======================================================================
// MethodBase — abstract interface for statistical methods.
//
// Each concrete method (SPAGRM, SPACox, WtCoxG, LEAF, …) derives from
// this and implements clone(), resultSize(), getHeaderColumns(), and
// getResultVec().  The engine clones one copy per worker thread.
// ======================================================================

class MethodBase {
  public:
    virtual ~MethodBase() = default;

    // Deep-copy for per-thread isolation (mutable scratch state, etc.).
    virtual std::unique_ptr<MethodBase> clone() const = 0;

    // Number of result columns after the 9 meta columns.
    virtual int resultSize() const = 0;

    // Tab-separated column headers for the result columns (leading tab included).
    virtual std::string getHeaderColumns() const = 0;

    // Called once per chunk before processing markers.  Default: no-op.
    virtual void prepareChunk(const std::vector<uint64_t> & /*gIndices*/) {
    }

    // Per-marker analysis.  The engine handles QC and imputation before
    // calling this.  `result` is cleared before each call.
    //
    //   GVec             — imputed genotype vector (counts bim col5 = ALT allele)
    //   altFreq          — ALT allele frequency (= freq of bim col5)
    //   markerInChunkIdx — 0-based index within current chunk
    //   result           — output: method-specific result values
    virtual void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result
    ) = 0;

    // Batch analysis: process B markers at once.  Default loops getResultVec().
    //
    //   GBatch    — imputed genotype matrix (N × B), one column per marker
    //   altFreqs  — ALT allele frequencies, one per marker (length B)
    //   results   — output: results[b] = method-specific result values for marker b
    //
    // Override in methods where the score computation is memory-bound and
    // can be converted from GEMV to GEMM (e.g. SPAsqr, SPAGRM).
    virtual void getResultBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &GBatch,
        const std::vector<double> &altFreqs,
        std::vector<std::vector<double> > &results
    ) {
        const int B = static_cast<int>(GBatch.cols());
        results.resize(B);
        for (int b = 0; b < B; ++b) {
            results[b].clear();
            // const_cast safe: getResultVec doesn't modify GVec for most methods.
            // For methods that do modify GVec, they must override getResultBatch.
            Eigen::VectorXd col = GBatch.col(b);
            getResultVec(col, altFreqs[b], b, results[b]);
        }
    }

    // Suggested marker batch size for this method (0 = use engine default).
    virtual int preferredBatchSize() const {
        return 0;
    }

    // ── Fused union-level GEMM interface ───────────────────────────────
    // Methods that compute Score = resid^T × G (or residMat^T × G for
    // multi-tau) can participate in a single fused GEMM across ALL
    // phenotypes, eliminating per-phenotype extraction entirely.
    //
    // Flow: engine builds AugResid (N_union × totalCols) from all fuseable
    // methods, does ONE GEMM per window → distributes score slices.

    // True if this method can participate in the fused GEMM.
    virtual bool supportsFusedGemm() const {
        return false;
    }

    // Number of residual columns (ntaus for SPAsqr, 1 for SPAGRM).
    virtual int fusedGemmColumns() const {
        return 0;
    }

    // Scatter per-phenotype residuals into union-level columns of dest.
    // dest is pre-zeroed (N_union × fusedGemmColumns()).
    // unionToLocal[i] maps union index → phenotype index (UINT32_MAX = absent).
    virtual void fillUnionResiduals(
        Eigen::Ref<Eigen::MatrixXd> dest,
        const std::vector<uint32_t> &unionToLocal
    ) const {
        (void)dest; (void)unionToLocal;
    }

    // Copy residual sums into dest (length = fusedGemmColumns()).
    virtual void fillResidualSums(double *dest) const {
        (void)dest;
    }

    // Process pre-computed raw scores (from the fused GEMM).
    //   scores   — fusedGemmColumns() × B matrix of raw scores
    //   gSums    — per-phenotype genotype sums, length B
    //   nUsed    — this phenotype's sample count (for gMean = gSum / nUsed)
    //   altFreqs — ALT allele frequencies, length B (pre-computed)
    //   results  — output: results[b] = method-specific result values
    virtual void processScoreBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &scores,
        const double *gSums,
        uint32_t nUsed,
        const std::vector<double> &altFreqs,
        std::vector<std::vector<double> > &results
    ) {
        (void)scores; (void)gSums; (void)nUsed;
        (void)altFreqs; (void)results;
    }

};

// ======================================================================
// Engine entry point
// ======================================================================

void markerEngine(
    const GenoMeta &genoData,
    const MethodBase &method,
    const std::string &outputFile,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff
);

// ======================================================================
// PhenoTask — one per phenotype for multiPhenoEngine
// ======================================================================

struct PhenoTask {
    std::string phenoName;
    std::unique_ptr<MethodBase> method; // single-phenotype method
    std::vector<uint32_t> unionToLocal; // union-dense index → pheno-dense index (UINT32_MAX = absent)
    uint32_t nUsed;                     // this phenotype's sample count
};

class TextWriter;

// Process a range of chunks [chunkStart, chunkEnd) through K phenotype tasks,
// writing results to pre-opened writers.  Writers must already have headers.
// Used internally by multiPhenoEngine and locoEngine.
void multiPhenoEngineRange(
    const GenoMeta &genoData,
    std::vector<PhenoTask> &tasks,
    const std::vector<std::string> &naSuffixes,
    size_t chunkStart,
    size_t chunkEnd,
    std::vector<TextWriter> &writers,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff
);



// Per-phenotype independent GWAS engine.
//
// Decodes genotypes once via the union mask, then for each phenotype:
//   - extracts per-phenotype genotype vector using unionToLocal
//   - computes per-phenotype allele stats (altFreq, MAC, missRate, HWE)
//   - applies QC filters per phenotype
//   - runs the method and writes to PREFIX.PHENO.METHOD[.gz|.zst]
//
// Threading: chunk-level work-stealing, K phenotypes sequential per chunk.
void multiPhenoEngine(
    const GenoMeta &genoData,
    std::vector<PhenoTask> &tasks,
    const std::string &outPrefix,
    const std::string &methodName,
    const std::string &compression,
    int compressionLevel,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff
);



// ======================================================================
// MultiMethod — runs N inner methods in one pass, interleaves results
//
// Column layout: for each residual name × each suffix
//   e.g. residNames=["R1","R2"], suffixes=["_P","_Z"]
//        → R1_P  R1_Z  R2_P  R2_Z
// ======================================================================

class MultiMethod : public MethodBase {
  public:
    MultiMethod(
        std::vector<std::unique_ptr<MethodBase> > methods,
        std::vector<std::string> residNames,
        std::vector<std::string> suffixes
    );

    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override;

    std::string getHeaderColumns() const override;

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
        std::vector<std::vector<double> > &results
    ) override;

    int preferredBatchSize() const override;

  private:
    std::vector<std::unique_ptr<MethodBase> > m_methods;
    std::vector<std::string> m_residNames;
    std::vector<std::string> m_suffixes;
};
