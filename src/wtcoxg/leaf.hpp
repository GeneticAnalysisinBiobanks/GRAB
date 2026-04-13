// leaf.hpp — LEAF: cluster-stratified WtCoxG with ancestry-matched reference AFs
//
// Workflow:
//   1. Load per-cluster resid files + per-cluster ref-af files (plink2 .afreq)
//   2. Load PLINK data with combined subjects
//   3. Match markers by (CHROM, ID) with allele flip detection
//   4. Per cluster: compute internal AF → summix → ancestry-matched AF_ref/obs_ct
//   5. Per cluster: batch-effect testing → refInfoMap
//   6. Create LEAFMethod (N WtCoxGMethod objects) → markerEngine → meta-analysis
#pragma once

#include "geno_factory/geno_data.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "engine/marker.hpp"

class WtCoxGMethod;
struct WtCoxGRefInfo;

// ======================================================================
// Per-population reference AF matched against .bim
// ======================================================================

struct PopMatchedAF {
    uint64_t genoIndex;
    double af;     // allele frequency aligned to .bim col5 orientation
    double obs_ct; // total allele number
};

// Load one plink2 .afreq file and match against PlinkData .bim markers
// by (CHROM, ID) with allele flip detection.  Returns matched entries.
std::vector<PopMatchedAF> loadAndMatchRefAf(
    const class GenoMeta &plinkData,
    const std::string &afreqFile
);

// ======================================================================
// K-means clustering
// ======================================================================

// K-means clustering (Lloyd's algorithm with K-means++ initialisation).
// Equivalent to R's kmeans(X, centers=k, nstart=nstart)$cluster.
//   X:       (n × p) feature matrix (e.g. principal components)
//   k:       number of clusters
//   nstart:  number of random restarts (best result kept)
// Returns integer cluster labels in {1, …, k} (length n).
Eigen::VectorXi kmeansCluster(
    const Eigen::Ref<const Eigen::MatrixXd> &X,
    int k,
    int nstart = 25,
    uint64_t seed = 0,
    int nthreads = 1
);

// ======================================================================
// Summix — ancestry proportion estimation
// ======================================================================

// Estimate ancestry proportions by minimising ||D*p − o||²
// subject to p ≥ 0, Σp = 1.
// Uses exact active-set enumeration (K ≤ 8 populations, up to 2^K subsets)
// with KKT / Lagrange-multiplier equality-constrained least squares.
//   observedAF: per-marker internal AF   (nSNP)
//   refAF:      reference pop AFs        (nSNP × nPop, nPop ≤ 8)
// Returns: proportion vector (nPop), sums to 1.
Eigen::VectorXd summixEstimate(
    const Eigen::VectorXd &observedAF,
    const Eigen::MatrixXd &refAF
);

// ======================================================================
// LEAFMethod — MethodBase implementation for cluster-stratified WtCoxG
// ======================================================================

class LEAFMethod : public MethodBase {
  public:
// Constructs from per-cluster WtCoxGMethod objects and index arrays.
//   clusterMethods:  one WtCoxGMethod per cluster (ownership transferred)
//   clusterIndices:  per-cluster indices into the full subject vector
    LEAFMethod(
        std::vector<std::unique_ptr<WtCoxGMethod> > clusterMethods,
        std::vector<std::vector<uint32_t> > clusterIndices
    );

// ---- MethodBase interface ----
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
    int m_nCluster;
    std::vector<std::unique_ptr<WtCoxGMethod> > m_clusterMethods;
    std::vector<std::vector<uint32_t> > m_clusterIndices;
    std::vector<Eigen::VectorXd> m_clusterGVec; // pre-allocated buffers
};

// ======================================================================
// Top-level orchestration — called from main()
// ======================================================================

// --pheno path: K-means clustering + per-cluster regression
void runLEAFPheno(
    const std::string &phenoFile,
    const std::string &covarFile,                             // empty = no separate covar file
    const std::vector<std::string> &covarNames,               // covariate columns for regression
    const std::vector<std::string> &phenoNames,               // selected phenotype columns
    const std::vector<std::string> &pcColNames,               // PC columns for K-means
    int nClusters,                                            // K for K-means
    uint64_t seed,                                            // RNG seed (0 = random)
    const GenoSpec &geno,
    const std::vector<std::string> &refAfFiles,               // one .afreq per ref pop
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

// Multi-phenotype entry point: shared K-means + geno scan + summix + GRM,
// parallel per-phenotype regression and batch-effect testing,
// single multiPhenoEngine call.
void runLEAF(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &covarNames,
    const std::vector<std::string> &phenoSpecs,               // "TIME:EVENT" or "COL"
    const std::vector<std::string> &pcColNames,
    int nClusters,
    uint64_t seed,
    const GenoSpec &geno,
    const std::vector<std::string> &refAfFiles,
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
