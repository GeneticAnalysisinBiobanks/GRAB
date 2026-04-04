// geno_prob.cpp — SPAGRM null model + marker runner (pure C++17 / Eigen / Boost)
//
// Translates SPAGRM.NullModel() + GRAB.mtMarker() from R to C++.
// See geno_prob.hpp for the public API.

#include "spagrm/geno_prob.hpp"
#include "spagrm/spagrm.hpp"
#include "engine/marker.hpp"
#include "io/geno_data.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/text_scanner.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>


namespace {

// ══════════════════════════════════════════════════════════════════════
// Constants (matching R defaults — not exposed as CLI flags)
// ══════════════════════════════════════════════════════════════════════
constexpr double MAX_QUANTILE    = 0.75;
constexpr double MIN_QUANTILE    = 0.25;
constexpr double INIT_OUTLIER_RATIO = 1.5;
constexpr bool   CONTROL_OUTLIER = true;
constexpr int    MAX_NUM_IN_FAM  = 5;
static const std::vector<double> MAF_INTERVAL =
    {0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5};
constexpr double ZETA_DEFAULT = 0.01;
constexpr double TOL_DEFAULT  = 1e-6;





// ══════════════════════════════════════════════════════════════════════
// Indexed IBD entry (dense subject indices, no strings)
// ══════════════════════════════════════════════════════════════════════
struct IndexedIBD {
  uint32_t idx1, idx2;
  double pa, pb, pc;
};

std::vector<IndexedIBD> loadIndexedIBD(
    const std::string& filename,
    const std::unordered_map<std::string, uint32_t>& idMap)
{
  std::ifstream ifs(filename);
  if (!ifs) throw std::runtime_error("Cannot open IBD file: " + filename);
  std::vector<IndexedIBD> out;
  std::string line;
  while (std::getline(ifs, line)) {
    if (text::skipLine(line)) continue;
    if (line.size() >= 3 && line[0] == 'I' && line[1] == 'D' && line[2] == '1') continue;
    if (line.size() >= 4 && line[0] == '#' && line[1] == 'I' && line[2] == 'D' && line[3] == '1') continue;
    text::TokenScanner tok(line);
    std::string id1 = tok.next();
    std::string id2 = tok.next();
    auto it1 = idMap.find(id1);
    auto it2 = idMap.find(id2);
    if (it1 == idMap.end() || it2 == idMap.end()) continue;
    char* endPtr;
    tok.skipWS(); double pa = std::strtod(tok.pos(), &endPtr); tok.p = endPtr;
    tok.skipWS(); double pb = std::strtod(tok.pos(), &endPtr); tok.p = endPtr;
    tok.skipWS(); double pc = std::strtod(tok.pos(), &endPtr);
    out.push_back({it1->second, it2->second, pa, pb, pc});
  }
  return out;
}

// Build hash map for O(1) IBD pair lookup: (min_idx << 32 | max_idx) → index
std::unordered_map<uint64_t, uint32_t> buildIBDPairMap(
    const std::vector<IndexedIBD>& ibdEntries)
{
  std::unordered_map<uint64_t, uint32_t> m;
  m.reserve(ibdEntries.size());
  for (uint32_t i = 0; i < static_cast<uint32_t>(ibdEntries.size()); ++i) {
    uint32_t lo = std::min(ibdEntries[i].idx1, ibdEntries[i].idx2);
    uint32_t hi = std::max(ibdEntries[i].idx1, ibdEntries[i].idx2);
    m[(static_cast<uint64_t>(lo) << 32) | hi] = i;
  }
  return m;
}


// ══════════════════════════════════════════════════════════════════════
// Pairwise IBD file loader
// ══════════════════════════════════════════════════════════════════════


// ══════════════════════════════════════════════════════════════════════
// Union-Find for connected components
// ══════════════════════════════════════════════════════════════════════
class UnionFind {
public:
  explicit UnionFind(uint32_t n) : parent_(n), rank_(n, 0) {
    std::iota(parent_.begin(), parent_.end(), 0u);
  }
  uint32_t find(uint32_t x) {
    while (parent_[x] != x) {
      parent_[x] = parent_[parent_[x]]; // path halving
      x = parent_[x];
    }
    return x;
  }
  void unite(uint32_t a, uint32_t b) {
    a = find(a); b = find(b);
    if (a == b) return;
    if (rank_[a] < rank_[b]) std::swap(a, b);
    parent_[b] = a;
    if (rank_[a] == rank_[b]) ++rank_[a];
  }
private:
  std::vector<uint32_t> parent_;
  std::vector<uint32_t> rank_;
};

// Return connected components as vectors of node indices.
std::vector<std::vector<uint32_t>> getComponents(
    uint32_t nNodes,
    const std::vector<std::pair<uint32_t, uint32_t>>& edges)
{
  UnionFind uf(nNodes);
  for (const auto& [a, b] : edges)
    uf.unite(a, b);

  std::unordered_map<uint32_t, std::vector<uint32_t>> compMap;
  for (uint32_t i = 0; i < nNodes; ++i)
    compMap[uf.find(i)].push_back(i);

  std::vector<std::vector<uint32_t>> result;
  result.reserve(compMap.size());
  for (auto& [root, members] : compMap)
    result.push_back(std::move(members));
  return result;
}


// ══════════════════════════════════════════════════════════════════════
// Quantile (R type=7 linear interpolation)
// ══════════════════════════════════════════════════════════════════════
double quantile_r7(std::vector<double>& sorted, double prob) {
  const size_t n = sorted.size();
  const double idx = prob * static_cast<double>(n - 1);
  const size_t lo = static_cast<size_t>(std::floor(idx));
  const size_t hi = std::min(lo + 1, n - 1);
  const double frac = idx - static_cast<double>(lo);
  return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}


// ══════════════════════════════════════════════════════════════════════
// Compute quadratic form R' * blockGRM * R for a (sub-)family.
// `entries` should be pre-filtered to the family's connected component.
// Members of the sub-family are given by `famMembers`; entries not
// touching those members are skipped.
// ══════════════════════════════════════════════════════════════════════
double familyQuadForm(
    const std::vector<uint32_t>& famMembers,
    const std::vector<SparseGRM::Entry>& entries,
    const Eigen::VectorXd& resid)
{
  std::unordered_set<uint32_t> famSet(famMembers.begin(), famMembers.end());

  // Sum the quadratic form:  sum_{entries in family} factor * value * R[i] * R[j]
  // factor=1 for diagonal, factor=2 for off-diagonal (since GRM stores lower-tri only)
  double result = 0.0;
  for (const auto& e : entries) {
    if (!famSet.count(e.row) || !famSet.count(e.col)) continue;
    double factor = (e.row == e.col) ? 1.0 : 2.0;
    result += factor * e.value * resid[e.row] * resid[e.col];
  }
  return result;
}


// ══════════════════════════════════════════════════════════════════════
// Prim's MST algorithm on small dense graphs (used for Chow-Liu tree)
// Returns edges as (parentIdx, childIdx) pairs.
// ══════════════════════════════════════════════════════════════════════
struct MSTEdge { int from, to; };

std::vector<MSTEdge> primMST(
    int N,
    const std::vector<std::vector<double>>& weight)  // weight[i][j], maximise
{
  std::vector<MSTEdge> edges;
  edges.reserve(N - 1);
  std::vector<bool> inTree(N, false);
  // key[i] = max weight edge connecting i to the tree
  std::vector<double> key(N, -std::numeric_limits<double>::infinity());
  std::vector<int> parent(N, -1);

  key[0] = 0;
  for (int iter = 0; iter < N; ++iter) {
    // pick vertex with max key not in tree
    int u = -1;
    double best = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < N; ++i) {
      if (!inTree[i] && key[i] > best) { best = key[i]; u = i; }
    }
    inTree[u] = true;
    if (parent[u] >= 0)
      edges.push_back({parent[u], u});

    for (int v = 0; v < N; ++v) {
      if (!inTree[v] && weight[u][v] > key[v]) {
        key[v] = weight[u][v];
        parent[v] = u;
      }
    }
  }
  return edges;
}


// ══════════════════════════════════════════════════════════════════════
// Chow-Liu tree builder — matches chow.liu.tree() in R
//
// For each MAF in MAF_interval, builds a joint genotype probability
// array of size 3^N, stored column-major (flattened).
//
// Returns matrix (3^N × nMAF).
// ══════════════════════════════════════════════════════════════════════
Eigen::MatrixXd buildChowLiuTree(
    int N,
    const std::vector<IndexedIBD>& familyIBD,
    const std::vector<uint32_t>& famMembers,
    const std::vector<double>& maf_interval)
{
  // Build global → local index
  std::unordered_map<uint32_t, int> globalToLocal;
  for (int i = 0; i < N; ++i)
    globalToLocal[famMembers[i]] = i;

  // Number of IBD entries (directed edges in the family)
  const int nIBD = static_cast<int>(familyIBD.size());

  // total array size = 3^N
  int arrSize = 1;
  for (int i = 0; i < N; ++i) arrSize *= 3;

  const int nMAF = static_cast<int>(maf_interval.size());
  Eigen::MatrixXd CLT(arrSize, nMAF);

  // Powers of 3 for indexing: index = sum_i geno[i] * stride[i]
  std::vector<int> stride(N);
  stride[0] = 1;
  for (int i = 1; i < N; ++i)
    stride[i] = stride[i - 1] * 3;

  // Dense entropy matrix for Prim's MST
  std::vector<std::vector<double>> entropyMat(N, std::vector<double>(N, 0.0));

  for (int mi = 0; mi < nMAF; ++mi) {
    const double mu  = maf_interval[mi];
    const double omu = 1.0 - mu;

    // Baseline genotype probabilities: P(G=0), P(G=1), P(G=2)
    const double p0[3] = {omu * omu, 2.0 * mu * omu, mu * mu};

    // IBD state joint distributions (3x3 row-major)
    // pa (allele2 = IBD2), pb (allele1 = IBD1), pc (allele0 = IBD0)
    double pa_arr[9] = {omu*omu, 0, 0,  0, 2*mu*omu, 0,  0, 0, mu*mu};

    double pb_arr[9] = {
      omu*omu*omu,     mu*omu*omu,       0,
      mu*omu*omu,      mu*omu,           mu*mu*omu,
      0,               mu*mu*omu,        mu*mu*mu
    };

    double pc_arr[9] = {
      omu*omu*omu*omu,         2*mu*omu*omu*omu,       mu*mu*omu*omu,
      2*mu*omu*omu*omu,        4*mu*mu*omu*omu,        2*mu*mu*mu*omu,
      mu*mu*omu*omu,           2*mu*mu*mu*omu,          mu*mu*mu*mu
    };

    // Compute entropy for each IBD pair
    for (int j = 0; j < nIBD; ++j) {
      auto it1 = globalToLocal.find(familyIBD[j].idx1);
      auto it2 = globalToLocal.find(familyIBD[j].idx2);
      if (it1 == globalToLocal.end() || it2 == globalToLocal.end()) continue;

      double entropy = 0.0;
      for (int k = 0; k < 9; ++k) {
        double pro = familyIBD[j].pa * pa_arr[k]
                   + familyIBD[j].pb * pb_arr[k]
                   + familyIBD[j].pc * pc_arr[k];
        if (pro > 0.0 && pc_arr[k] > 0.0)
          entropy += pro * std::log(pro / pc_arr[k]);
      }
      entropyMat[it1->second][it2->second] = entropy;
      entropyMat[it2->second][it1->second] = entropy;
    }

    // Prim's maximum spanning tree
    std::vector<MSTEdge> mstEdges = primMST(N, entropyMat);

    // Retrieve IBD values for MST edges
    struct MSTEdgeIBD { int idx1, idx2; double pa, pb, pc; };
    std::vector<MSTEdgeIBD> mstIBD;
    mstIBD.reserve(mstEdges.size());
    for (const auto& e : mstEdges) {
      // Find the IBD entry matching this edge
      double epa = 0, epb = 0, epc = 1;
      for (const auto& ibd : familyIBD) {
        auto it1 = globalToLocal.find(ibd.idx1);
        auto it2 = globalToLocal.find(ibd.idx2);
        if (it1 == globalToLocal.end() || it2 == globalToLocal.end()) continue;
        if ((it1->second == e.from && it2->second == e.to) ||
            (it1->second == e.to   && it2->second == e.from)) {
          epa = ibd.pa; epb = ibd.pb; epc = ibd.pc;
          break;
        }
      }
      mstIBD.push_back({e.from, e.to, epa, epb, epc});
    }

    // Count how many times each node appears as non-root in the tree
    // (nodes appearing in multiple edges need division by marginal)
    std::vector<int> edgeCount(N, 0);
    for (const auto& e : mstEdges) {
      edgeCount[e.from]++;
      edgeCount[e.to]++;
    }
    // The "vec" in R code: duplicated entries among concatenated from/to
    // These are the internal (non-leaf) nodes that need marginal division
    std::vector<int> divideNodes;
    for (int i = 0; i < N; ++i) {
      // A node needs division (N-2) times if it appears in multiple edges
      // R code: vec <- c(mst.IBD$idxID1, mst.IBD$idxID2); vec <- vec[duplicated(vec)]
      if (edgeCount[i] > 1) {
        for (int d = 0; d < edgeCount[i] - 1; ++d)
          divideNodes.push_back(i);
      }
    }

    // Build probability array
    std::vector<double> arr(arrSize, 1.0);

    // For each MST edge, multiply in the pairwise conditional
    for (const auto& e : mstIBD) {
      double pro[9];
      for (int k = 0; k < 9; ++k)
        pro[k] = e.pa * pa_arr[k] + e.pb * pb_arr[k] + e.pc * pc_arr[k];

      // Multiply into the N-dimensional array
      // For each cell in the array, multiply by pro[g_idx1 * 3 + g_idx2]
      for (int idx = 0; idx < arrSize; ++idx) {
        int g1 = (idx / stride[e.idx1]) % 3;
        int g2 = (idx / stride[e.idx2]) % 3;
        arr[idx] *= pro[g1 * 3 + g2];
      }
    }

    // Divide by marginal for internal nodes
    for (int node : divideNodes) {
      for (int idx = 0; idx < arrSize; ++idx) {
        int g = (idx / stride[node]) % 3;
        arr[idx] /= p0[g];
      }
    }

    // Store column
    for (int idx = 0; idx < arrSize; ++idx)
      CLT(idx, mi) = arr[idx];
  }

  return CLT;
}


// ══════════════════════════════════════════════════════════════════════
// Build the stand.S array for a family of size N
//
// In R:
//   arr.index[[n]] is a list of N arrays, each of size 3^N.
//   arr.index[[n]][[i]] has value residuals[i] at the appropriate
//   genotype of subject i, and 1 elsewhere.
//   stand.S = array(rowSums(mapply(function(x,y) x*y, arr.index, Resid)), rep(3,N))
//
// This computes the equivalent: for each cell in the 3^N array,
//   stand.S[cell] = sum_i ( Resid[i] * geno_of_i_in_cell )
// where geno_of_i_in_cell is the genotype (0,1,2) of subject i
// encoded in the cell index.
// ══════════════════════════════════════════════════════════════════════
std::vector<double> buildStandS(int N, const std::vector<double>& resid) {
  int arrSize = 1;
  for (int i = 0; i < N; ++i) arrSize *= 3;

  std::vector<int> stride(N);
  stride[0] = 1;
  for (int i = 1; i < N; ++i)
    stride[i] = stride[i - 1] * 3;

  std::vector<double> standS(arrSize, 0.0);
  for (int idx = 0; idx < arrSize; ++idx) {
    double val = 0.0;
    for (int i = 0; i < N; ++i) {
      int g = (idx / stride[i]) % 3;
      val += static_cast<double>(g) * resid[i];
    }
    standS[idx] = val;
  }
  return standS;
}


// ══════════════════════════════════════════════════════════════════════
// SPAGRMMethod — MethodBase adapter wrapping a single SPAGRMClass
// ══════════════════════════════════════════════════════════════════════
class SPAGRMMethod : public MethodBase {
public:
  explicit SPAGRMMethod(SPAGRMClass spagrm) : m_spagrm(std::move(spagrm)) {}

  std::unique_ptr<MethodBase> clone() const override {
    return std::make_unique<SPAGRMMethod>(*this);
  }

  int resultSize() const override { return 2; }

  std::string getHeaderColumns() const override {
    return "\tSPAGRM_P\tSPAGRM_Z";
  }

  void getResultVec(
      Eigen::Ref<Eigen::VectorXd> GVec,
      double altFreq,
      int /*markerInChunkIdx*/,
      std::vector<double>& result) override
  {
    result.clear();
    double z;
    double p = m_spagrm.getMarkerPval(GVec, altFreq, z);
    result.push_back(p);
    result.push_back(z);
  }

private:
  SPAGRMClass m_spagrm;
};


} // anonymous namespace


// ══════════════════════════════════════════════════════════════════════
// buildSPAGRMNullModel — per-column null model construction
//
// Builds outlier detection, family R'GRM R terms, Chow-Liu trees, and
// returns a ready-to-use SPAGRMClass.  Called once per residual column
// in the multi-column loop.
// ══════════════════════════════════════════════════════════════════════

static SPAGRMClass buildSPAGRMNullModel(
    const Eigen::VectorXd& Resid,
    uint32_t N,
    const std::unordered_set<uint32_t>& singletonSet,
    const std::vector<double>& grmDiag,
    const std::vector<std::vector<uint32_t>>& families,
    const std::vector<std::vector<SparseGRM::Entry>>& familyEntries,
    const std::vector<SparseGRM::Entry>& allGrmEntries,
    const std::vector<IndexedIBD>& ibdEntries,
    const std::unordered_map<uint64_t, uint32_t>& ibdPairMap,
    double spaCutoff)
{
  // ── Outlier detection (IQR) ──────────────────────────────────────
  std::vector<double> sortedResid(N);
  for (uint32_t i = 0; i < N; ++i) sortedResid[i] = Resid[i];
  std::sort(sortedResid.begin(), sortedResid.end());

  const double Q1 = quantile_r7(sortedResid, MIN_QUANTILE);
  const double Q3 = quantile_r7(sortedResid, MAX_QUANTILE);
  const double IQR = Q3 - Q1;

  double outlierRatio = INIT_OUTLIER_RATIO;
  double cutLo = Q1 - outlierRatio * IQR;
  double cutHi = Q3 + outlierRatio * IQR;

  std::vector<bool> isOutlier(N, false);
  auto recomputeOutliers = [&]() {
    cutLo = Q1 - outlierRatio * IQR;
    cutHi = Q3 + outlierRatio * IQR;
    int cnt = 0;
    for (uint32_t i = 0; i < N; ++i) {
      isOutlier[i] = (Resid[i] < cutLo || Resid[i] > cutHi);
      if (isOutlier[i]) ++cnt;
    }
    return cnt;
  };

  int nOutlier = recomputeOutliers();

  if (CONTROL_OUTLIER) {
    while (nOutlier == 0) {
      outlierRatio *= 0.8;
      nOutlier = recomputeOutliers();
    }
    while (static_cast<double>(nOutlier) / N > 0.05) {
      outlierRatio += 0.5;
      nOutlier = recomputeOutliers();
    }
  }

  // ── Accumulate global variance terms for singletons ──────────────
  double R_GRM_R = 0.0;
  double sum_R_nonOutlier = 0.0;
  double R_GRM_R_nonOutlier = 0.0;

  for (uint32_t idx : singletonSet) {
    double diagVal = grmDiag[idx];
    double contrib = diagVal * Resid[idx] * Resid[idx];
    R_GRM_R += contrib;
    if (!isOutlier[idx])
      R_GRM_R_nonOutlier += contrib;
  }
  for (uint32_t idx : singletonSet) {
    if (!isOutlier[idx])
      sum_R_nonOutlier += Resid[idx];
  }

  std::vector<double> unrelatedOutlierResids;
  for (uint32_t idx : singletonSet) {
    if (isOutlier[idx])
      unrelatedOutlierResids.push_back(Resid[idx]);
  }

  double R_GRM_R_TwoSubjOutlier = 0.0;

  std::vector<std::array<double, 2>> twoSubj_resid_list;
  std::vector<std::vector<double>> twoSubj_rho_list;
  std::vector<std::vector<double>> threeSubj_standS_list;
  std::vector<Eigen::MatrixXd> threeSubj_CLT_list;

  // ── Process families ─────────────────────────────────────────────
  std::vector<std::vector<uint32_t>> outlierFamilies;

  for (size_t fi = 0; fi < families.size(); ++fi) {
    const auto& fam = families[fi];

    bool hasOutlier = false;
    for (uint32_t idx : fam)
      if (isOutlier[idx]) { hasOutlier = true; break; }

    double famQuad = familyQuadForm(fam, familyEntries[fi], Resid);
    R_GRM_R += famQuad;

    if (!hasOutlier) {
      double famSum = 0.0;
      for (uint32_t idx : fam) famSum += Resid[idx];
      sum_R_nonOutlier += famSum;
      R_GRM_R_nonOutlier += famQuad;
      continue;
    }

    if (static_cast<int>(fam.size()) <= MAX_NUM_IN_FAM) {
      outlierFamilies.push_back(fam);
      continue;
    }

    // ── Greedy family splitting ──────────────────────────────────
    std::unordered_set<uint32_t> famSet(fam.begin(), fam.end());

    struct OffDiagEntry {
      uint32_t row, col;
      double value, cov;
    };
    std::vector<OffDiagEntry> famEdges;
    for (const auto& e : familyEntries[fi]) {
      if (e.row == e.col) continue;
      if (e.row < e.col)
        famEdges.push_back({e.row, e.col, e.value,
                           std::abs(e.value * Resid[e.row] * Resid[e.col])});
    }
    std::sort(famEdges.begin(), famEdges.end(),
              [](const OffDiagEntry& a, const OffDiagEntry& b) {
                return a.cov < b.cov;
              });

    std::vector<bool> edgeRemoved(famEdges.size(), false);
    int removeUpTo = -1;
    for (size_t j = 0; j < famEdges.size(); ++j) {
      edgeRemoved[j] = true;
      std::vector<std::pair<uint32_t, uint32_t>> remainingEdges;
      for (size_t k = 0; k < famEdges.size(); ++k) {
        if (!edgeRemoved[k])
          remainingEdges.push_back({famEdges[k].row, famEdges[k].col});
      }
      auto subComps = getComponents(N, remainingEdges);
      int maxComp = 0;
      for (const auto& sc : subComps) {
        int cnt = 0;
        for (uint32_t idx : sc) if (famSet.count(idx)) ++cnt;
        maxComp = std::max(maxComp, cnt);
      }
      if (maxComp <= MAX_NUM_IN_FAM) {
        removeUpTo = static_cast<int>(j);
        break;
      }
    }

    if (removeUpTo >= 0) {
      for (int j = removeUpTo; j >= 0; --j) {
        edgeRemoved[j] = false;
        std::vector<std::pair<uint32_t, uint32_t>> remainingEdges;
        for (size_t k = 0; k < famEdges.size(); ++k) {
          if (!edgeRemoved[k])
            remainingEdges.push_back({famEdges[k].row, famEdges[k].col});
        }
        auto subComps = getComponents(N, remainingEdges);
        int maxComp = 0;
        for (const auto& sc : subComps) {
          int cnt = 0;
          for (uint32_t idx : sc) if (famSet.count(idx)) ++cnt;
          maxComp = std::max(maxComp, cnt);
        }
        if (maxComp > MAX_NUM_IN_FAM) edgeRemoved[j] = true;
      }
    }

    std::vector<std::pair<uint32_t, uint32_t>> finalEdges;
    for (size_t k = 0; k < famEdges.size(); ++k) {
      if (!edgeRemoved[k])
        finalEdges.push_back({famEdges[k].row, famEdges[k].col});
    }
    auto subComps = getComponents(N, finalEdges);

    for (const auto& sc : subComps) {
      std::vector<uint32_t> subFam;
      for (uint32_t idx : sc) if (famSet.count(idx)) subFam.push_back(idx);
      if (subFam.empty()) continue;

      bool subHasOutlier = false;
      for (uint32_t idx : subFam) if (isOutlier[idx]) { subHasOutlier = true; break; }

      double subQuad = familyQuadForm(subFam, familyEntries[fi], Resid);

      if (!subHasOutlier) {
        double subSum = 0.0;
        for (uint32_t idx : subFam) subSum += Resid[idx];
        sum_R_nonOutlier += subSum;
        R_GRM_R_nonOutlier += subQuad;
      } else {
        outlierFamilies.push_back(std::move(subFam));
      }
    }
  }

  // ── Build Chow-Liu trees for outlier families ────────────────────
  for (const auto& fam : outlierFamilies) {
    const int n1 = static_cast<int>(fam.size());

    if (n1 == 1) {
      unrelatedOutlierResids.push_back(Resid[fam[0]]);
      continue;
    }

    if (n1 == 2) {
      double R1 = Resid[fam[0]], R2 = Resid[fam[1]];
      double pairQuad = familyQuadForm(fam, allGrmEntries, Resid);
      R_GRM_R_TwoSubjOutlier += pairQuad;

      uint32_t lo = std::min(fam[0], fam[1]);
      uint32_t hi = std::max(fam[0], fam[1]);
      uint64_t key = (static_cast<uint64_t>(lo) << 32) | hi;
      double pa = 0, pb = 0;
      auto ibdIt = ibdPairMap.find(key);
      if (ibdIt != ibdPairMap.end()) {
        pa = ibdEntries[ibdIt->second].pa;
        pb = ibdEntries[ibdIt->second].pb;
      }

      double Rho = pa + 0.5 * pb;
      double midterm = std::sqrt(std::max(Rho * Rho - pa, 0.0));

      twoSubj_resid_list.push_back({R1, R2});
      twoSubj_rho_list.push_back({Rho + midterm, Rho - midterm});
      continue;
    }

    // Three-or-more subject outlier family
    std::vector<double>     famResid(n1);
    for (int i = 0; i < n1; ++i) {
      famResid[i] = Resid[fam[i]];
    }

    std::unordered_set<uint32_t> famIdxSet(fam.begin(), fam.end());
    std::vector<IndexedIBD> famIBD;
    for (const auto& ibd : ibdEntries) {
      if (famIdxSet.count(ibd.idx1) && famIdxSet.count(ibd.idx2))
        famIBD.push_back(ibd);
    }

    Eigen::MatrixXd CLT = buildChowLiuTree(n1, famIBD, fam, MAF_INTERVAL);
    std::vector<double> standS = buildStandS(n1, famResid);

    threeSubj_standS_list.push_back(std::move(standS));
    threeSubj_CLT_list.push_back(std::move(CLT));
  }

  // ── Assemble SPAGRMClass ─────────────────────────────────────────
  Eigen::VectorXd residOutliers = Eigen::Map<Eigen::VectorXd>(
      unrelatedOutlierResids.data(),
      static_cast<Eigen::Index>(unrelatedOutlierResids.size()));

  nsSPAGRM::FamilyData fd;
  fd.resid_unrelated_outliers = std::move(residOutliers);
  fd.twoSubj_resid = std::move(twoSubj_resid_list);
  fd.twoSubj_rho   = std::move(twoSubj_rho_list);
  fd.threeSubj_standS = std::move(threeSubj_standS_list);
  fd.threeSubj_CLT    = std::move(threeSubj_CLT_list);

  return SPAGRMClass(
      Resid,
      sum_R_nonOutlier,
      R_GRM_R_nonOutlier,
      R_GRM_R_TwoSubjOutlier,
      R_GRM_R,
      MAF_INTERVAL,
      std::move(fd),
      spaCutoff,
      ZETA_DEFAULT,
      TOL_DEFAULT);
}


// ══════════════════════════════════════════════════════════════════════
// runSPAGRM — entry point
// ══════════════════════════════════════════════════════════════════════

void runSPAGRM(
    const std::string& residFile,
    const std::string& spgrmGrabFile,
    const std::string& spgrmGctaFile,
    const std::string& pairwiseIBDFile,
    const GenoSpec& geno,
    const std::string& outputFile,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff)
{
  // ══════════════════════════════════════════════════════════════════
  // 1. Load residual file (2 columns: SubjID, Resid)
  // ══════════════════════════════════════════════════════════════════
  infoMsg("Loading residual file: %s", residFile.c_str());
  auto famIIDs = parseGenoIIDs(geno);
  SubjectData sd(std::move(famIIDs));
  sd.loadResidOne(residFile);
  sd.finalize();
  const uint32_t N = sd.nUsed();
  infoMsg("Loaded %u subjects (intersected with .fam)", N);

  auto subjIDs = sd.usedIIDs();
  auto subjIdMap = text::buildIIDMap(subjIDs);

  // ══════════════════════════════════════════════════════════════════
  // 2. Load sparse GRM (unsymmetrised — lower-tri + diagonal)
  // ══════════════════════════════════════════════════════════════════
  SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile,
                                   subjIDs, sd.famIIDs());
  infoMsg("Sparse GRM: %zu entries (diagonal + off-diag)", grm.nnz());

  // ══════════════════════════════════════════════════════════════════
  // 3. Load pairwise IBD (indexed — no strings in hot paths)
  // ══════════════════════════════════════════════════════════════════
  infoMsg("Loading pairwise IBD: %s", pairwiseIBDFile.c_str());
  auto ibdEntries = loadIndexedIBD(pairwiseIBDFile, subjIdMap);
  auto ibdPairMap = buildIBDPairMap(ibdEntries);
  infoMsg("Loaded %zu IBD records", ibdEntries.size());

  // ══════════════════════════════════════════════════════════════════
  // 4. Build GRM topology (components, singletons, families)
  //    and pre-bucket entries by component
  // ══════════════════════════════════════════════════════════════════
  const auto& allEntries = grm.entries();

  // Extract unique off-diagonal edges for connected-component detection
  std::vector<std::pair<uint32_t, uint32_t>> edges;
  {
    std::unordered_set<uint64_t> seen;
    for (const auto& e : allEntries) {
      if (e.row == e.col) continue;
      uint32_t lo = std::min(e.row, e.col);
      uint32_t hi = std::max(e.row, e.col);
      uint64_t key = (static_cast<uint64_t>(lo) << 32) | hi;
      if (seen.insert(key).second)
        edges.push_back({lo, hi});
    }
  }
  auto components = getComponents(N, edges);
  infoMsg("Found %zu connected components", components.size());

  std::vector<std::vector<uint32_t>> singletons, families;
  for (auto& comp : components) {
    if (comp.size() == 1)
      singletons.push_back(std::move(comp));
    else
      families.push_back(std::move(comp));
  }
  infoMsg("Singletons: %zu, Families: %zu", singletons.size(), families.size());

  std::unordered_set<uint32_t> singletonSet;
  for (const auto& s : singletons) singletonSet.insert(s[0]);

  // Use cached GRM diagonal for singleton variance
  const auto& grmDiag = grm.diagonal();

  // Pre-bucket: per-family GRM entries
  std::unordered_map<uint32_t, size_t> subjToFamily;
  for (size_t fi = 0; fi < families.size(); ++fi)
    for (uint32_t idx : families[fi])
      subjToFamily[idx] = fi;

  std::vector<std::vector<SparseGRM::Entry>> familyEntries(families.size());
  for (const auto& e : allEntries) {
    auto it = subjToFamily.find(e.row);
    if (it != subjToFamily.end())
      familyEntries[it->second].push_back(e);
  }

  // ══════════════════════════════════════════════════════════════════
  // 5. Load PLINK data
  // ══════════════════════════════════════════════════════════════════
  auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(),
                               nSnpPerChunk);

  // ══════════════════════════════════════════════════════════════════
  // 6. Per-residual-column loop
  // ══════════════════════════════════════════════════════════════════
  const int nRC = sd.residOneCols();
  if (nRC > 1)
    infoMsg("Multi-column residual file: %d columns", nRC);

  for (int rc = 0; rc < nRC; ++rc) {
    Eigen::VectorXd colBuf;
    if (nRC > 1) colBuf = sd.residMatrix().col(rc);
    const Eigen::VectorXd& Resid = (nRC > 1) ? colBuf : sd.residuals();

    std::string outFile = (nRC == 1) ? outputFile
        : (outputFile + "." + std::to_string(rc + 1) + ".gz");
    if (nRC > 1)
      infoMsg("  Column %d/%d%s -> %s", rc + 1, nRC,
              (rc < static_cast<int>(sd.residColNames().size())
                   ? (" (" + sd.residColNames()[rc] + ")").c_str() : ""),
              outFile.c_str());

    SPAGRMClass spagrm = buildSPAGRMNullModel(
        Resid, N, singletonSet, grmDiag, families, familyEntries,
        allEntries, ibdEntries, ibdPairMap, spaCutoff);

    SPAGRMMethod method(std::move(spagrm));

    if (nRC == 1) infoMsg("Starting SPAGRM marker-level association (%d threads)", nthreads);
    markerEngine(*genoData, method, outFile,
                 nthreads, missingCutoff, minMafCutoff, minMacCutoff,
                 /*exactHwe=*/false);
  }
}
