// grm_null.hpp — Shared GRM null-model building utilities
//
// Extracted from geno_prob.cpp so that both SPAGRM and SAGELD (and any
// future GRM-based method) can reuse the same topology decomposition,
// outlier detection, family splitting, and Chow-Liu tree construction.
#pragma once

#include "spagrm/spagrm.hpp"
#include "io/sparse_grm.hpp"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

namespace nsGRMNull {

// ══════════════════════════════════════════════════════════════════════
// Constants (matching R defaults — not exposed as CLI flags)
// ══════════════════════════════════════════════════════════════════════
constexpr double MAX_QUANTILE = 0.75;
constexpr double MIN_QUANTILE = 0.25;
constexpr double INIT_OUTLIER_RATIO = 1.5;
constexpr bool CONTROL_OUTLIER = true;
constexpr int MAX_NUM_IN_FAM = 5;
constexpr double ZETA_DEFAULT = 0.01;
constexpr double TOL_DEFAULT = 1e-6;

/// Build MAF grid from QC cutoffs: starts at min(mafCutoff, macCutoff/(2*n)),
/// uses half-decade steps (1-3-10-30-...) up to 0.5.
std::vector<double> buildMafInterval(
    double mafCutoff,
    double macCutoff,
    uint32_t nUsed
);

// ══════════════════════════════════════════════════════════════════════
// Indexed IBD entry (dense subject indices, no strings)
// ══════════════════════════════════════════════════════════════════════
struct IndexedIBD {
    uint32_t idx1, idx2;
    double pa, pb, pc;
};

/// Load pairwise IBD file into dense-indexed entries.
std::vector<IndexedIBD> loadIndexedIBD(
    const std::string &filename,
    const std::unordered_map<std::string, uint32_t> &idMap
);

/// Build hash map for O(1) IBD pair lookup: (min_idx << 32 | max_idx) → index.
std::unordered_map<uint64_t, uint32_t> buildIBDPairMap(const std::vector<IndexedIBD> &ibdEntries);

// ══════════════════════════════════════════════════════════════════════
// Union-Find for connected components
// ══════════════════════════════════════════════════════════════════════
class UnionFind {
  public:
    explicit UnionFind(uint32_t n);
    uint32_t find(uint32_t x);

    void unite(
        uint32_t a,
        uint32_t b
    );

    uint32_t componentSize(uint32_t x) {
        return size_[find(x)];
    }

  private:
    std::vector<uint32_t> parent_;
    std::vector<uint32_t> rank_;
    std::vector<uint32_t> size_;
};

/// Return connected components as vectors of node indices.
std::vector<std::vector<uint32_t> > getComponents(
    uint32_t nNodes,
    const std::vector<std::pair<uint32_t, uint32_t> > &edges
);

// ══════════════════════════════════════════════════════════════════════
// Quantile (R type=7 linear interpolation)
// ══════════════════════════════════════════════════════════════════════
double quantile_r7(
    std::vector<double> &sorted,
    double prob
);

// ══════════════════════════════════════════════════════════════════════
// Quadratic form R' * blockGRM * R for a (sub-)family.
// ══════════════════════════════════════════════════════════════════════
double familyQuadForm(
    const std::vector<uint32_t> &famMembers,
    const std::vector<SparseGRM::Entry> &entries,
    const Eigen::VectorXd &resid
);

// ══════════════════════════════════════════════════════════════════════
// Prim's MST on small dense graphs (used for Chow-Liu tree)
// ══════════════════════════════════════════════════════════════════════
struct MSTEdge {
    int from, to;
};

std::vector<MSTEdge> primMST(
    int N,
    const std::vector<std::vector<double> > &weight
);

// ══════════════════════════════════════════════════════════════════════
// Chow-Liu tree builder — 3^N × nMAF probability matrix
// ══════════════════════════════════════════════════════════════════════
Eigen::MatrixXd buildChowLiuTree(
    int N,
    const std::vector<IndexedIBD> &familyIBD,
    const std::vector<uint32_t> &famMembers,
    const std::vector<double> &maf_interval
);

// ══════════════════════════════════════════════════════════════════════
// Build stand.S array for a family of size N
// ══════════════════════════════════════════════════════════════════════
std::vector<double> buildStandS(
    int N,
    const std::vector<double> &resid
);

// ══════════════════════════════════════════════════════════════════════
// buildSPAGRMNullModel — per-column null model construction
//
// Shared by SPAGRM and SAGELD.  Takes a residual vector and everything
// loaded by the run* entry point, returns a ready-to-use SPAGRMClass.
// ══════════════════════════════════════════════════════════════════════
SPAGRMClass buildSPAGRMNullModel(
    const Eigen::VectorXd &Resid,
    uint32_t N,
    const std::unordered_set<uint32_t> &singletonSet,
    const std::vector<double> &grmDiag,
    const std::vector<std::vector<uint32_t> > &families,
    const std::vector<std::vector<SparseGRM::Entry> > &familyEntries,
    const std::vector<SparseGRM::Entry> &allGrmEntries,
    const std::vector<IndexedIBD> &ibdEntries,
    const std::unordered_map<uint64_t, uint32_t> &ibdPairMap,
    double spaCutoff,
    double minMafCutoff,
    double minMacCutoff,
    int nthreads = 1
);

} // namespace nsGRMNull
