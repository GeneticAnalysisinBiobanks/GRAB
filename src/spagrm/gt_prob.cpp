// gt_prob.cpp — SPAGRM null model + marker runner (pure C++17 / Eigen / Boost)
//
// Translates SPAGRM.NullModel() + GRAB.mtMarker() from R to C++.
// See gt_prob.hpp for the public API.

#include "spagrm/gt_prob.hpp"
#include "spagrm/spagrm.hpp"
#include "engine/marker.hpp"
#include "io/plink.hpp"
#include "io/sparse_grm.hpp"
#include "util/logging.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
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
// Simple 2-column residual file loader: SubjID  Resid
// Skips '#' header lines.  Returns parallel vectors.
// ══════════════════════════════════════════════════════════════════════
struct ResidEntry { std::string id; double resid; };

std::vector<ResidEntry> loadResid2Col(const std::string& filename) {
  std::ifstream ifs(filename);
  if (!ifs) throw std::runtime_error("Cannot open residual file: " + filename);
  std::vector<ResidEntry> out;
  std::string line;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    const char*       p   = line.c_str();
    const char* const end = p + line.size();
    auto skipWS  = [&]() { while (p < end && (*p == ' ' || *p == '\t')) ++p; };
    auto nextTok = [&]() -> std::string {
      skipWS();
      const char* s = p;
      while (p < end && *p != ' ' && *p != '\t') ++p;
      return std::string(s, p);
    };
    std::string id = nextTok();
    skipWS();
    if (id.empty() || p >= end) continue;
    char* endPtr;
    double val = std::strtod(p, &endPtr);
    if (endPtr == p)
      throw std::runtime_error("Bad residual value in " + filename + " for ID " + id);
    out.push_back({std::move(id), val});
  }
  return out;
}


// ══════════════════════════════════════════════════════════════════════
// GRM raw entry (with both directions kept separate)
// ══════════════════════════════════════════════════════════════════════
struct GRMEntry {
  uint32_t row, col;
  double value;
};

// Parse sparse GRM into indexed entries, keeping BOTH directions for
// off-diagonal (for connected-component detection and quadratic form).
// Subjects not in idMap are silently dropped.
std::vector<GRMEntry> parseGRM(
    const std::string& filename,
    const std::unordered_map<std::string, uint32_t>& idMap)
{
  std::ifstream ifs(filename);
  if (!ifs) throw std::runtime_error("Cannot open sparse GRM: " + filename);
  std::vector<GRMEntry> entries;
  std::string line;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    const char*       p   = line.c_str();
    const char* const end = p + line.size();
    auto skipWS  = [&]() { while (p < end && (*p == ' ' || *p == '\t')) ++p; };
    auto nextTok = [&]() -> std::string {
      skipWS();
      const char* s = p;
      while (p < end && *p != ' ' && *p != '\t') ++p;
      return std::string(s, p);
    };
    skipWS();
    if (p >= end) continue;
    std::string id1 = nextTok();
    std::string id2 = nextTok();
    skipWS();
    if (p >= end) continue;
    char* endPtr;
    double val = std::strtod(p, &endPtr);
    if (endPtr == p) continue;
    auto it1 = idMap.find(id1);
    auto it2 = idMap.find(id2);
    if (it1 == idMap.end() || it2 == idMap.end()) continue;
    entries.push_back({it1->second, it2->second, val});
  }
  return entries;
}


// ══════════════════════════════════════════════════════════════════════
// Pairwise IBD file loader
// ══════════════════════════════════════════════════════════════════════
struct IBDEntry {
  std::string id1, id2;
  double pa, pb, pc;
};

std::vector<IBDEntry> loadIBD(
    const std::string& filename,
    const std::unordered_set<std::string>& subjSet)
{
  std::ifstream ifs(filename);
  if (!ifs) throw std::runtime_error("Cannot open IBD file: " + filename);
  std::vector<IBDEntry> out;
  std::string line;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;
    // Skip header line "ID1 ID2 pa pb pc"
    if (line.substr(0, 3) == "ID1") continue;
    const char*       p   = line.c_str();
    const char* const end = p + line.size();
    auto skipWS  = [&]() { while (p < end && (*p == ' ' || *p == '\t')) ++p; };
    auto nextTok = [&]() -> std::string {
      skipWS();
      const char* s = p;
      while (p < end && *p != ' ' && *p != '\t') ++p;
      return std::string(s, p);
    };
    std::string id1 = nextTok();
    std::string id2 = nextTok();
    if (subjSet.count(id1) == 0 || subjSet.count(id2) == 0) continue;
    char* endPtr;
    skipWS(); double pa = std::strtod(p, &endPtr); p = endPtr;
    skipWS(); double pb = std::strtod(p, &endPtr); p = endPtr;
    skipWS(); double pc = std::strtod(p, &endPtr);
    out.push_back({std::move(id1), std::move(id2), pa, pb, pc});
  }
  return out;
}


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
// Compute quadratic form R' * blockGRM * R for a small family.
// blockGRM is implicitly constructed from grmEntries that fall within
// the family members.  Indices in grmEntries use global subject
// ordering; famMembers maps global→local within the family.
// ══════════════════════════════════════════════════════════════════════
double familyQuadForm(
    const std::vector<uint32_t>& famMembers,
    const std::vector<GRMEntry>& grmEntries,
    const Eigen::VectorXd& resid)
{
  // Build fast local lookup
  std::unordered_map<uint32_t, uint32_t> globalToLocal;
  globalToLocal.reserve(famMembers.size());
  for (uint32_t i = 0; i < static_cast<uint32_t>(famMembers.size()); ++i)
    globalToLocal[famMembers[i]] = i;

  // Collect local residuals
  const int n = static_cast<int>(famMembers.size());
  Eigen::VectorXd localR(n);
  for (int i = 0; i < n; ++i)
    localR[i] = resid[famMembers[i]];

  // Sum the quadratic form:  sum_{entries in family} factor * value * R[i] * R[j]
  // factor=1 for diagonal, factor=2 for off-diagonal (since GRM stores lower-tri only)
  double result = 0.0;
  for (const auto& e : grmEntries) {
    auto it1 = globalToLocal.find(e.row);
    auto it2 = globalToLocal.find(e.col);
    if (it1 == globalToLocal.end() || it2 == globalToLocal.end()) continue;
    double factor = (e.row == e.col) ? 1.0 : 2.0;
    result += factor * e.value * localR[it1->second] * localR[it2->second];
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
    const std::vector<IBDEntry>& familyIBD,
    const std::vector<std::string>& familyIDs,
    const std::vector<double>& maf_interval)
{
  // Build ID → local index
  std::unordered_map<std::string, int> idToIdx;
  for (int i = 0; i < N; ++i)
    idToIdx[familyIDs[i]] = i;

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
      auto it1 = idToIdx.find(familyIBD[j].id1);
      auto it2 = idToIdx.find(familyIBD[j].id2);
      if (it1 == idToIdx.end() || it2 == idToIdx.end()) continue;

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
        auto it1 = idToIdx.find(ibd.id1);
        auto it2 = idToIdx.find(ibd.id2);
        if (it1 == idToIdx.end() || it2 == idToIdx.end()) continue;
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

  int resultSize() const override { return 2; }  // zScore, Pvalue

  std::string getHeaderColumns() const override {
    return "\tzScore\tPvalue";
  }

  void getResultVec(
      Eigen::Ref<Eigen::VectorXd> GVec,
      double altFreq,
      int /*markerInChunkIdx*/,
      bool /*flipped*/,
      std::vector<double>& result) override
  {
    result.clear();
    double z;
    double p = m_spagrm.getMarkerPval(GVec, altFreq, z);
    result.push_back(z);
    result.push_back(p);
  }

private:
  SPAGRMClass m_spagrm;
};


} // anonymous namespace


// ══════════════════════════════════════════════════════════════════════
// runSPAGRM — entry point
// ══════════════════════════════════════════════════════════════════════

void runSPAGRM(
    const std::string& residFile,
    const std::string& sparseGrmFile,
    const std::string& pairwiseIBDFile,
    const std::string& bfilePrefix,
    const std::string& outputFile,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff)
{
  // ══════════════════════════════════════════════════════════════════
  // 1. Load residual matrix (2 columns: SubjID, Resid)
  // ══════════════════════════════════════════════════════════════════
  infoMsg("Loading residual file: %s", residFile.c_str());
  auto residEntries = loadResid2Col(residFile);
  const uint32_t N = static_cast<uint32_t>(residEntries.size());
  infoMsg("Loaded %u subjects from residual file", N);

  // Build subject list and residual vector
  std::vector<std::string> subjIDs(N);
  Eigen::VectorXd Resid(N);
  std::unordered_map<std::string, uint32_t> subjIdMap;
  subjIdMap.reserve(N);
  for (uint32_t i = 0; i < N; ++i) {
    subjIDs[i] = std::move(residEntries[i].id);
    Resid[i]   = residEntries[i].resid;
    subjIdMap.emplace(subjIDs[i], i);
  }
  std::unordered_set<std::string> subjSet(subjIDs.begin(), subjIDs.end());

  // ══════════════════════════════════════════════════════════════════
  // 2. Load sparse GRM
  // ══════════════════════════════════════════════════════════════════
  infoMsg("Loading sparse GRM: %s", sparseGrmFile.c_str());
  auto grmEntries = parseGRM(sparseGrmFile, subjIdMap);
  infoMsg("Sparse GRM: %zu entries (diagonal + off-diag)", grmEntries.size());

  // ══════════════════════════════════════════════════════════════════
  // 3. Load pairwise IBD
  // ══════════════════════════════════════════════════════════════════
  infoMsg("Loading pairwise IBD: %s", pairwiseIBDFile.c_str());
  auto ibdEntries = loadIBD(pairwiseIBDFile, subjSet);
  infoMsg("Loaded %zu IBD records", ibdEntries.size());

  // ══════════════════════════════════════════════════════════════════
  // 4. Outlier detection (IQR-based, matching R ControlOutlier logic)
  // ══════════════════════════════════════════════════════════════════
  infoMsg("Detecting outliers (IQR method)");

  // Sort for quantile computation
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
  infoMsg("Outlier cutoffs: [%.3f, %.3f], outliers: %d", cutLo, cutHi, nOutlier);

  if (CONTROL_OUTLIER) {
    // Shrink ratio until we get at least 1 outlier
    while (nOutlier == 0) {
      outlierRatio *= 0.8;
      nOutlier = recomputeOutliers();
      infoMsg("Adjusted cutoffs: [%.3f, %.3f], outliers: %d", cutLo, cutHi, nOutlier);
    }
    // Expand ratio if outliers exceed 5%
    while (static_cast<double>(nOutlier) / N > 0.05) {
      outlierRatio += 0.5;
      nOutlier = recomputeOutliers();
      infoMsg("Reducing outliers: cutoffs [%.3f, %.3f], outliers: %d (%.1f%%)",
              cutLo, cutHi, nOutlier, 100.0 * nOutlier / N);
    }
  }
  infoMsg("Final outlier count: %d / %u (%.1f%%)",
          nOutlier, N, 100.0 * nOutlier / N);

  // ══════════════════════════════════════════════════════════════════
  // 5. Compute |Cov| = |value * R[i] * R[j]| for off-diagonal GRM
  //    entries, used to weight edges for greedy family splitting.
  // ══════════════════════════════════════════════════════════════════

  // Separate diagonal and off-diagonal entries
  struct OffDiagEntry {
    uint32_t row, col;
    double value;
    double cov;  // |value * R[row] * R[col]|
  };
  std::vector<OffDiagEntry> offDiag;
  offDiag.reserve(grmEntries.size());
  for (const auto& e : grmEntries) {
    if (e.row != e.col) {
      offDiag.push_back({e.row, e.col, e.value,
                         std::abs(e.value * Resid[e.row] * Resid[e.col])});
    }
  }

  // Build edges for connected components (unique undirected)
  std::vector<std::pair<uint32_t, uint32_t>> edges;
  {
    std::unordered_set<uint64_t> seen;
    for (const auto& e : offDiag) {
      uint32_t lo = std::min(e.row, e.col);
      uint32_t hi = std::max(e.row, e.col);
      uint64_t key = (static_cast<uint64_t>(lo) << 32) | hi;
      if (seen.insert(key).second)
        edges.push_back({lo, hi});
    }
  }

  // ══════════════════════════════════════════════════════════════════
  // 6. Connected-component decomposition
  // ══════════════════════════════════════════════════════════════════
  auto components = getComponents(N, edges);
  infoMsg("Found %zu connected components", components.size());

  // Classify into singletons (unrelated) and families
  std::vector<std::vector<uint32_t>> singletons, families;
  for (auto& comp : components) {
    if (comp.size() == 1)
      singletons.push_back(std::move(comp));
    else
      families.push_back(std::move(comp));
  }
  infoMsg("Singletons: %zu, Families: %zu", singletons.size(), families.size());

  // ══════════════════════════════════════════════════════════════════
  // 7. Accumulate global variance terms
  // ══════════════════════════════════════════════════════════════════

  // For singletons (unrelated): count their GRM diagonal contribution
  std::unordered_set<uint32_t> singletonSet;
  for (const auto& s : singletons) singletonSet.insert(s[0]);

  std::unordered_set<uint32_t> singletonNonOutlierSet;
  for (uint32_t idx : singletonSet) {
    if (!isOutlier[idx]) singletonNonOutlierSet.insert(idx);
  }

  // R_GRM_R: sum over all entries where both subjects are singletons
  // (This is just the diagonal entries for singletons since they have no off-diag within-group.)
  double R_GRM_R = 0.0;
  double sum_R_nonOutlier = 0.0;
  double R_GRM_R_nonOutlier = 0.0;

  // Unrelated singletons: their GRM contribution is just the diagonal entry
  for (const auto& e : grmEntries) {
    if (e.row == e.col && singletonSet.count(e.row)) {
      double contrib = e.value * Resid[e.row] * Resid[e.col];
      R_GRM_R += contrib;
      if (!isOutlier[e.row])
        R_GRM_R_nonOutlier += contrib;
    }
  }
  for (uint32_t idx : singletonSet) {
    if (!isOutlier[idx])
      sum_R_nonOutlier += Resid[idx];
  }

  // Unrelated outlier residuals
  std::vector<double> unrelatedOutlierResids;
  for (uint32_t idx : singletonSet) {
    if (isOutlier[idx])
      unrelatedOutlierResids.push_back(Resid[idx]);
  }

  double R_GRM_R_TwoSubjOutlier = 0.0;

  // Lists for SPAGRMClass
  std::vector<std::array<double, 2>> twoSubj_resid_list;
  std::vector<std::vector<double>> twoSubj_rho_list;
  std::vector<std::vector<double>> threeSubj_standS_list;
  std::vector<Eigen::MatrixXd> threeSubj_CLT_list;

  // ══════════════════════════════════════════════════════════════════
  // 8. Process families
  // ══════════════════════════════════════════════════════════════════
  if (!families.empty()) {
    infoMsg("Processing %zu family groups", families.size());

    // Structures to hold outlier families after greedy splitting
    std::vector<std::vector<uint32_t>> outlierFamilies;

    for (size_t fi = 0; fi < families.size(); ++fi) {
      const auto& fam = families[fi];

      // Check if any member is an outlier
      bool hasOutlier = false;
      for (uint32_t idx : fam) {
        if (isOutlier[idx]) { hasOutlier = true; break; }
      }

      // Compute R' * GRM * R for this family
      double famQuad = familyQuadForm(fam, grmEntries, Resid);
      R_GRM_R += famQuad;

      if (!hasOutlier) {
        // No outliers: add to non-outlier sums
        double famSum = 0.0;
        for (uint32_t idx : fam) famSum += Resid[idx];
        sum_R_nonOutlier += famSum;
        R_GRM_R_nonOutlier += famQuad;
        continue;
      }

      // Family has outliers
      if (static_cast<int>(fam.size()) <= MAX_NUM_IN_FAM) {
        outlierFamilies.push_back(fam);
        continue;
      }

      // ── Greedy family splitting ──────────────────────────────────
      // Need to split this family so max component <= MAX_NUM_IN_FAM
      //
      // Step 1: Remove edges in ascending |Cov| order until max component <= K
      // Step 2: Add back removed edges in descending |Cov| if they don't violate

      // Collect edges within this family, sorted by Cov ascending
      std::unordered_set<uint32_t> famSet(fam.begin(), fam.end());
      std::vector<OffDiagEntry> famEdges;
      for (const auto& e : offDiag) {
        if (famSet.count(e.row) && famSet.count(e.col)) {
          // Keep only one direction (lower < upper)
          if (e.row < e.col)
            famEdges.push_back(e);
        }
      }
      std::sort(famEdges.begin(), famEdges.end(),
                [](const OffDiagEntry& a, const OffDiagEntry& b) {
                  return a.cov < b.cov;
                });

      // Remove edges one by one until max component <= MAX_NUM_IN_FAM
      std::vector<bool> edgeRemoved(famEdges.size(), false);
      int removeUpTo = -1;
      for (size_t j = 0; j < famEdges.size(); ++j) {
        edgeRemoved[j] = true;

        // Build remaining edges and check components
        std::vector<std::pair<uint32_t, uint32_t>> remainingEdges;
        for (size_t k = 0; k < famEdges.size(); ++k) {
          if (!edgeRemoved[k])
            remainingEdges.push_back({famEdges[k].row, famEdges[k].col});
        }
        auto subComps = getComponents(N, remainingEdges);
        // Check max component size (only among family members)
        int maxComp = 0;
        for (const auto& sc : subComps) {
          int cnt = 0;
          for (uint32_t idx : sc) {
            if (famSet.count(idx)) ++cnt;
          }
          maxComp = std::max(maxComp, cnt);
        }
        if (maxComp <= MAX_NUM_IN_FAM) {
          removeUpTo = static_cast<int>(j);
          break;
        }
      }

      // Step 2: Add back removed edges in descending Cov order
      if (removeUpTo >= 0) {
        for (int j = removeUpTo; j >= 0; --j) {
          edgeRemoved[j] = false;  // try adding back

          std::vector<std::pair<uint32_t, uint32_t>> remainingEdges;
          for (size_t k = 0; k < famEdges.size(); ++k) {
            if (!edgeRemoved[k])
              remainingEdges.push_back({famEdges[k].row, famEdges[k].col});
          }
          auto subComps = getComponents(N, remainingEdges);
          int maxComp = 0;
          for (const auto& sc : subComps) {
            int cnt = 0;
            for (uint32_t idx : sc) {
              if (famSet.count(idx)) ++cnt;
            }
            maxComp = std::max(maxComp, cnt);
          }
          if (maxComp > MAX_NUM_IN_FAM) {
            edgeRemoved[j] = true;  // can't add back
          }
        }
      }

      // Build final sub-components from remaining edges
      std::vector<std::pair<uint32_t, uint32_t>> finalEdges;
      for (size_t k = 0; k < famEdges.size(); ++k) {
        if (!edgeRemoved[k])
          finalEdges.push_back({famEdges[k].row, famEdges[k].col});
      }
      auto subComps = getComponents(N, finalEdges);

      // Process each sub-component
      for (const auto& sc : subComps) {
        // Filter to family members only
        std::vector<uint32_t> subFam;
        for (uint32_t idx : sc) {
          if (famSet.count(idx)) subFam.push_back(idx);
        }
        if (subFam.empty()) continue;

        bool subHasOutlier = false;
        for (uint32_t idx : subFam) {
          if (isOutlier[idx]) { subHasOutlier = true; break; }
        }

        double subQuad = familyQuadForm(subFam, grmEntries, Resid);

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

    // ══════════════════════════════════════════════════════════════════
    // 9. Build Chow-Liu trees for outlier families
    // ══════════════════════════════════════════════════════════════════
    infoMsg("Building Chow-Liu trees for %zu outlier families",
            outlierFamilies.size());

    for (const auto& fam : outlierFamilies) {
      const int n1 = static_cast<int>(fam.size());

      if (n1 == 1) {
        // Singleton outlier from family splitting
        unrelatedOutlierResids.push_back(Resid[fam[0]]);
        continue;
      }

      if (n1 == 2) {
        // Two-subject outlier family
        double R1 = Resid[fam[0]];
        double R2 = Resid[fam[1]];

        // Compute quad form for this pair
        double pairQuad = familyQuadForm(fam, grmEntries, Resid);
        R_GRM_R_TwoSubjOutlier += pairQuad;

        // Find the IBD entry for this pair
        const std::string& sid1 = subjIDs[fam[0]];
        const std::string& sid2 = subjIDs[fam[1]];
        double pa = 0, pb = 0;
        for (const auto& ibd : ibdEntries) {
          if ((ibd.id1 == sid1 && ibd.id2 == sid2) ||
              (ibd.id1 == sid2 && ibd.id2 == sid1)) {
            pa = ibd.pa; pb = ibd.pb;
            break;
          }
        }

        double Rho = pa + 0.5 * pb;
        double midterm = std::sqrt(std::max(Rho * Rho - pa, 0.0));

        twoSubj_resid_list.push_back({R1, R2});
        twoSubj_rho_list.push_back({Rho + midterm, Rho - midterm});
        continue;
      }

      // Three-or-more subject outlier family
      // Collect family IDs and residuals
      std::vector<std::string> famIDs(n1);
      std::vector<double>     famResid(n1);
      for (int i = 0; i < n1; ++i) {
        famIDs[i]   = subjIDs[fam[i]];
        famResid[i] = Resid[fam[i]];
      }

      // Collect IBD entries within this family
      std::unordered_set<std::string> famIDSet(famIDs.begin(), famIDs.end());
      std::vector<IBDEntry> famIBD;
      for (const auto& ibd : ibdEntries) {
        if (famIDSet.count(ibd.id1) && famIDSet.count(ibd.id2))
          famIBD.push_back(ibd);
      }

      // Build Chow-Liu tree probability matrix
      Eigen::MatrixXd CLT = buildChowLiuTree(n1, famIBD, famIDs, MAF_INTERVAL);

      // Build stand.S vector
      std::vector<double> standS = buildStandS(n1, famResid);

      threeSubj_standS_list.push_back(std::move(standS));
      threeSubj_CLT_list.push_back(std::move(CLT));
    }
  }

  infoMsg("SPAGRM null model summary:");
  infoMsg("  R_GRM_R = %.6f", R_GRM_R);
  infoMsg("  R_GRM_R_nonOutlier = %.6f", R_GRM_R_nonOutlier);
  infoMsg("  R_GRM_R_TwoSubjOutlier = %.6f", R_GRM_R_TwoSubjOutlier);
  infoMsg("  sum_R_nonOutlier = %.6f", sum_R_nonOutlier);
  infoMsg("  Unrelated outliers: %zu", unrelatedOutlierResids.size());
  infoMsg("  Two-subject families: %zu", twoSubj_resid_list.size());
  infoMsg("  Three+ subject families: %zu", threeSubj_standS_list.size());

  // ══════════════════════════════════════════════════════════════════
  // 10. Build SPAGRMClass
  // ══════════════════════════════════════════════════════════════════
  Eigen::VectorXd residOutliers = Eigen::Map<Eigen::VectorXd>(
      unrelatedOutlierResids.data(),
      static_cast<Eigen::Index>(unrelatedOutlierResids.size()));

  nsSPAGRM::FamilyData fam;
  fam.resid_unrelated_outliers = std::move(residOutliers);
  fam.twoSubj_resid = std::move(twoSubj_resid_list);
  fam.twoSubj_rho   = std::move(twoSubj_rho_list);
  fam.threeSubj_standS = std::move(threeSubj_standS_list);
  fam.threeSubj_CLT    = std::move(threeSubj_CLT_list);

  SPAGRMClass spagrm(
      Resid,
      sum_R_nonOutlier,
      R_GRM_R_nonOutlier,
      R_GRM_R_TwoSubjOutlier,
      R_GRM_R,
      MAF_INTERVAL,
      std::move(fam),
      spaCutoff,
      ZETA_DEFAULT,
      TOL_DEFAULT);

  // ══════════════════════════════════════════════════════════════════
  // 11. Load PLINK and run marker engine
  // ══════════════════════════════════════════════════════════════════
  infoMsg("Loading PLINK files: %s", bfilePrefix.c_str());
  PlinkData plinkData(
      bfilePrefix + ".bed",
      bfilePrefix + ".bim",
      bfilePrefix + ".fam",
      subjIDs,
      "ref-first",
      {}, {}, {}, {},
      nSnpPerChunk);

  SPAGRMMethod method(std::move(spagrm));

  infoMsg("Starting SPAGRM marker-level association (%d threads)", nthreads);
  markerEngine(plinkData, method, outputFile,
               nthreads, missingCutoff, minMafCutoff, minMacCutoff,
               /*exactHwe=*/false);
}
