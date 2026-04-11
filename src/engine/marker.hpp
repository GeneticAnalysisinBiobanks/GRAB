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
    virtual void prepareChunk(const std::vector<uint64_t> & /*gIndices*/) {}

    // Per-marker analysis.  The engine handles QC and imputation before
    // calling this.  `result` is cleared before each call.
    //
    //   GVec             — imputed genotype vector (counts bim col5 = ALT allele)
    //   altFreq          — ALT allele frequency (= freq of bim col5)
    //   markerInChunkIdx — 0-based index within current chunk
    //   result           — output: method-specific result values
    virtual void getResultVec(Eigen::Ref<Eigen::VectorXd> GVec,
                              double altFreq,
                              int markerInChunkIdx,
                              std::vector<double> &result) = 0;
};

// ======================================================================
// Engine entry point
// ======================================================================

void markerEngine(const GenoMeta &genoData,
                  const MethodBase &method,
                  const std::string &outputFile,
                  int nthreads,
                  double missingCutoff,
                  double minMafCutoff,
                  double minMacCutoff,
                  double hweCutoff);

// ======================================================================
// PhenoTask — one per phenotype for multiPhenoEngine
// ======================================================================

struct PhenoTask {
    std::string phenoName;
    std::unique_ptr<MethodBase> method; // single-phenotype method
    std::vector<uint32_t> unionToLocal; // union-dense index → pheno-dense index (UINT32_MAX = absent)
    uint32_t nUsed;                     // this phenotype's sample count
};

// Per-phenotype independent GWAS engine.
//
// Decodes genotypes once via the union mask, then for each phenotype:
//   - extracts per-phenotype genotype vector using unionToLocal
//   - computes per-phenotype allele stats (altFreq, MAC, missRate, HWE)
//   - applies QC filters per phenotype
//   - runs the method and writes to PREFIX.PHENO.METHOD[.gz|.zst]
//
// Threading: chunk-level work-stealing, K phenotypes sequential per chunk.
void multiPhenoEngine(const GenoMeta &genoData,
                      std::vector<PhenoTask> &tasks,
                      const std::string &outPrefix,
                      const std::string &methodName,
                      const std::string &compression,
                      int compressionLevel,
                      int nthreads,
                      double missingCutoff,
                      double minMafCutoff,
                      double minMacCutoff,
                      double hweCutoff);

// ======================================================================
// MultiMethod — runs N inner methods in one pass, interleaves results
//
// Column layout: for each residual name × each suffix
//   e.g. residNames=["R1","R2"], suffixes=["_P","_Z"]
//        → R1_P  R1_Z  R2_P  R2_Z
// ======================================================================

class MultiMethod : public MethodBase {
  public:
    MultiMethod(std::vector<std::unique_ptr<MethodBase>> methods,
                std::vector<std::string> residNames,
                std::vector<std::string> suffixes);

    std::unique_ptr<MethodBase> clone() const override;
    int resultSize() const override;
    std::string getHeaderColumns() const override;
    void prepareChunk(const std::vector<uint64_t> &gIndices) override;
    void getResultVec(Eigen::Ref<Eigen::VectorXd> GVec,
                      double altFreq,
                      int markerInChunkIdx,
                      std::vector<double> &result) override;

  private:
    std::vector<std::unique_ptr<MethodBase>> m_methods;
    std::vector<std::string> m_residNames;
    std::vector<std::string> m_suffixes;
};

// ======================================================================
// buildResidNames — resolve residual column names for output headers
//
// If the --null-resid file had a header, colNames carries those names;
// otherwise colNames is empty and we generate "R1", "R2", ...
// ======================================================================

inline std::vector<std::string> buildResidNames(const std::vector<std::string> &colNames, int nRC) {
    if (!colNames.empty()) return colNames;
    std::vector<std::string> names;
    names.reserve(nRC);
    for (int i = 0; i < nRC; ++i)
        names.push_back("R" + std::to_string(i + 1));
    return names;
}
