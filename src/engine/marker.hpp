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

    // Expand internal residuals from per-phenotype-dense (Np) to
    // union-dense (nUnion) space, filling absent entries with 0.
    // After this call, getResultVec() expects GVec of size nUnion.
    virtual void padToUnionSpace(
        const uint32_t * /*unionToLocal*/,
        uint32_t                                                             /*nUnion*/
    ) {
    }

    // Return true if this method supports the impute engine
    // (i.e. padToUnionSpace has been implemented and called).
    virtual bool supportsImputeEngine() const {
        return false;
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

// Impute-mode multi-phenotype engine: all phenotypes share one union GVec.
//
// Assumes every task's method has been padToUnionSpace'd so that residuals
// are union-dense (missing subjects have residual 0).  Decodes genotypes
// once, computes union-level stats once, and calls each method with the
// same GVec.  Eliminates K per-phenotype GVec buffers and extraction loops.
void imputeMultiPhenoEngineRange(
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

// Impute-mode wrapper: pads tasks to union space, then delegates to
// imputeMultiPhenoEngineRange.  Same interface as multiPhenoEngine.
void imputeMultiPhenoEngine(
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

  private:
    std::vector<std::unique_ptr<MethodBase> > m_methods;
    std::vector<std::string> m_residNames;
    std::vector<std::string> m_suffixes;
};
