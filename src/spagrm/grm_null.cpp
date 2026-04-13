// grm_null.cpp — Shared GRM null-model building utilities (implementation)
//
// See grm_null.hpp for the public API.

#include "spagrm/grm_null.hpp"
#include "spagrm/spagrm.hpp"
#include "io/sparse_grm.hpp"
#include "util/logging.hpp"
#include "util/text_scanner.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numeric>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

namespace nsGRMNull {

// ══════════════════════════════════════════════════════════════════════
// buildMafInterval — dynamic MAF grid from QC cutoffs
// ══════════════════════════════════════════════════════════════════════

std::vector<double> buildMafInterval(
    double mafCutoff,
    double macCutoff,
    uint32_t nUsed
) {
    // Effective minimum MAF that can pass QC
    double effMin = mafCutoff;
    if (nUsed > 0) {
        double macBasedMaf = macCutoff / (2.0 * nUsed);
        effMin = std::min(effMin, macBasedMaf);
    }
    // Clamp to a reasonable floor (avoid degenerate grids)
    if (effMin < 1e-7) effMin = 1e-7;
    if (effMin > 0.01) effMin = 0.01;

    // Build half-decade grid: effMin, 3*effMin, 10*effMin, 30*effMin, ...
    // up to 0.1, then linear steps 0.2, 0.3, 0.4, 0.5
    std::vector<double> grid;
    grid.push_back(effMin);
    double val = effMin;
    // Half-decade steps: multiply alternately by 3 and ~3.33 (i.e., 3, 10, 30, 100, ...)
    bool step3 = true;
    while (true) {
        val = step3 ? val * 3.0 : val * (10.0 / 3.0);
        step3 = !step3;
        if (val >= 0.1) break;
        grid.push_back(val);
    }
    // Fixed upper bins (common variants)
    for (double v : {0.1, 0.2, 0.3, 0.4, 0.5}) {
        if (grid.back() < v) grid.push_back(v);
    }
    return grid;
}

// ══════════════════════════════════════════════════════════════════════
// loadIndexedIBD
// ══════════════════════════════════════════════════════════════════════

std::vector<IndexedIBD> loadIndexedIBD(
    const std::string &filename,
    const std::unordered_map<std::string, uint32_t> &idMap
) {
    std::ifstream ifs(filename);
    if (!ifs) throw std::runtime_error("Cannot open IBD file: " + filename);
    std::vector<IndexedIBD> out;
    std::string line;
    // Validate required header: first non-comment line must start with ID1 or #ID1
    bool headerFound = false;
    while (std::getline(ifs, line)) {
        if (text::skipLine(line)) continue;
        // Check header prefix
        if ((line.size() >= 3 && line[0] == 'I' && line[1] == 'D' && line[2] == '1') ||
            (line.size() >= 4 && line[0] == '#' && line[1] == 'I' && line[2] == 'D' && line[3] == '1')) {
            headerFound = true;
            break;
        }
        throw std::runtime_error(filename + ": first line must be a header starting with ID1 or #ID1");
    }
    if (!headerFound) throw std::runtime_error(filename + ": empty IBD file (no header)");
    while (std::getline(ifs, line)) {
        if (text::skipLine(line)) continue;
        text::TokenScanner tok(line);
        std::string id1 = tok.next();
        std::string id2 = tok.next();
        auto it1 = idMap.find(id1);
        auto it2 = idMap.find(id2);
        if (it1 == idMap.end() || it2 == idMap.end()) continue;
        char *endPtr;
        tok.skipWS();
        double pa = std::strtod(tok.pos(), &endPtr);
        tok.p = endPtr;
        tok.skipWS();
        double pb = std::strtod(tok.pos(), &endPtr);
        tok.p = endPtr;
        tok.skipWS();
        double pc = std::strtod(tok.pos(), &endPtr);
        out.push_back({it1->second, it2->second, pa, pb, pc});
    }
    return out;
}

// ══════════════════════════════════════════════════════════════════════
// buildIBDPairMap
// ══════════════════════════════════════════════════════════════════════

std::unordered_map<uint64_t, uint32_t> buildIBDPairMap(const std::vector<IndexedIBD> &ibdEntries) {
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
// UnionFind
// ══════════════════════════════════════════════════════════════════════

UnionFind::UnionFind(uint32_t n)
    : parent_(n),
      rank_(n, 0),
      size_(n, 1)
{
    std::iota(parent_.begin(), parent_.end(), 0u);
}

uint32_t UnionFind::find(uint32_t x) {
    while (parent_[x] != x) {
        parent_[x] = parent_[parent_[x]]; // path halving
        x = parent_[x];
    }
    return x;
}

void UnionFind::unite(
    uint32_t a,
    uint32_t b
) {
    a = find(a);
    b = find(b);
    if (a == b) return;
    if (rank_[a] < rank_[b]) std::swap(a, b);
    parent_[b] = a;
    size_[a] += size_[b];
    if (rank_[a] == rank_[b]) ++rank_[a];
}

// ══════════════════════════════════════════════════════════════════════
// getComponents
// ══════════════════════════════════════════════════════════════════════

std::vector<std::vector<uint32_t> > getComponents(
    uint32_t nNodes,
    const std::vector<std::pair<uint32_t, uint32_t> > &edges
) {
    UnionFind uf(nNodes);
    for (const auto &[a, b] : edges)
        uf.unite(a, b);

    std::unordered_map<uint32_t, std::vector<uint32_t> > compMap;
    for (uint32_t i = 0; i < nNodes; ++i)
        compMap[uf.find(i)].push_back(i);

    std::vector<std::vector<uint32_t> > result;
    result.reserve(compMap.size());
    for (auto &[root, members] : compMap)
        result.push_back(std::move(members));
    return result;
}

// ══════════════════════════════════════════════════════════════════════
// quantile_r7
// ══════════════════════════════════════════════════════════════════════

double quantile_r7(
    std::vector<double> &sorted,
    double prob
) {
    const size_t n = sorted.size();
    const double idx = prob * static_cast<double>(n - 1);
    const size_t lo = static_cast<size_t>(std::floor(idx));
    const size_t hi = std::min(lo + 1, n - 1);
    const double frac = idx - static_cast<double>(lo);
    return sorted[lo] * (1.0 - frac) + sorted[hi] * frac;
}

// ══════════════════════════════════════════════════════════════════════
// familyQuadForm
// ══════════════════════════════════════════════════════════════════════

double familyQuadForm(
    const std::vector<uint32_t> &famMembers,
    const std::vector<SparseGRM::Entry> &entries,
    const Eigen::VectorXd &resid
) {
    std::unordered_set<uint32_t> famSet(famMembers.begin(), famMembers.end());

    double result = 0.0;
    for (const auto &e : entries) {
        if (!famSet.count(e.row) || !famSet.count(e.col)) continue;
        double factor = (e.row == e.col) ? 1.0 : 2.0;
        result += factor * e.value * resid[e.row] * resid[e.col];
    }
    return result;
}

// ══════════════════════════════════════════════════════════════════════
// primMST
// ══════════════════════════════════════════════════════════════════════

std::vector<MSTEdge> primMST(
    int N,
    const std::vector<std::vector<double> > &weight
) {
    std::vector<MSTEdge> edges;
    edges.reserve(N - 1);
    std::vector<bool> inTree(N, false);
    std::vector<double> key(N, -std::numeric_limits<double>::infinity());
    std::vector<int> parent(N, -1);

    key[0] = 0;
    for (int iter = 0; iter < N; ++iter) {
        int u = -1;
        double best = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < N; ++i) {
            if (!inTree[i] && key[i] > best) {
                best = key[i];
                u = i;
            }
        }
        inTree[u] = true;
        if (parent[u] >= 0) edges.push_back({parent[u], u});

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
// buildChowLiuTree
// ══════════════════════════════════════════════════════════════════════

Eigen::MatrixXd buildChowLiuTree(
    int N,
    const std::vector<IndexedIBD> &familyIBD,
    const std::vector<uint32_t> &famMembers,
    const std::vector<double> &maf_interval
) {
    std::unordered_map<uint32_t, int> globalToLocal;
    for (int i = 0; i < N; ++i)
        globalToLocal[famMembers[i]] = i;

    const int nIBD = static_cast<int>(familyIBD.size());

    int arrSize = 1;
    for (int i = 0; i < N; ++i)
        arrSize *= 3;

    const int nMAF = static_cast<int>(maf_interval.size());
    Eigen::MatrixXd CLT(arrSize, nMAF);

    std::vector<int> stride(N);
    stride[0] = 1;
    for (int i = 1; i < N; ++i)
        stride[i] = stride[i - 1] * 3;

    std::vector<std::vector<double> > entropyMat(N, std::vector<double>(N, 0.0));

    for (int mi = 0; mi < nMAF; ++mi) {
        const double mu = maf_interval[mi];
        const double omu = 1.0 - mu;

        const double p0[3] = {omu * omu, 2.0 * mu * omu, mu * mu};

        double pa_arr[9] = {omu * omu, 0, 0, 0, 2 * mu * omu, 0, 0, 0, mu * mu};

        double pb_arr[9] = {omu * omu * omu, mu * omu * omu, 0, mu * omu * omu, mu * omu, mu * mu * omu, 0,
                            mu * mu * omu,   mu * mu * mu};

        double pc_arr[9] = {omu * omu * omu * omu,    2 * mu * omu * omu * omu, mu * mu * omu * omu,
                            2 * mu * omu * omu * omu, 4 * mu * mu * omu * omu,  2 * mu * mu * mu * omu,
                            mu * mu * omu * omu,      2 * mu * mu * mu * omu,   mu * mu * mu * mu};

        for (int j = 0; j < nIBD; ++j) {
            auto it1 = globalToLocal.find(familyIBD[j].idx1);
            auto it2 = globalToLocal.find(familyIBD[j].idx2);
            if (it1 == globalToLocal.end() || it2 == globalToLocal.end()) continue;

            double entropy = 0.0;
            for (int k = 0; k < 9; ++k) {
                double pro = familyIBD[j].pa * pa_arr[k] + familyIBD[j].pb * pb_arr[k] + familyIBD[j].pc * pc_arr[k];
                if (pro > 0.0 && pc_arr[k] > 0.0) entropy += pro * std::log(pro / pc_arr[k]);
            }
            entropyMat[it1->second][it2->second] = entropy;
            entropyMat[it2->second][it1->second] = entropy;
        }

        std::vector<MSTEdge> mstEdges = primMST(N, entropyMat);

        struct MSTEdgeIBD {
            int idx1, idx2;
            double pa, pb, pc;
        };

        std::vector<MSTEdgeIBD> mstIBD;
        mstIBD.reserve(mstEdges.size());
        for (const auto &e : mstEdges) {
            double epa = 0, epb = 0, epc = 1;
            for (const auto &ibd : familyIBD) {
                auto it1 = globalToLocal.find(ibd.idx1);
                auto it2 = globalToLocal.find(ibd.idx2);
                if (it1 == globalToLocal.end() || it2 == globalToLocal.end()) continue;
                if ((it1->second == e.from && it2->second == e.to) || (it1->second == e.to && it2->second == e.from)) {
                    epa = ibd.pa;
                    epb = ibd.pb;
                    epc = ibd.pc;
                    break;
                }
            }
            mstIBD.push_back({e.from, e.to, epa, epb, epc});
        }

        std::vector<int> edgeCount(N, 0);
        for (const auto &e : mstEdges) {
            edgeCount[e.from]++;
            edgeCount[e.to]++;
        }
        std::vector<int> divideNodes;
        for (int i = 0; i < N; ++i) {
            if (edgeCount[i] > 1) {
                for (int d = 0; d < edgeCount[i] - 1; ++d)
                    divideNodes.push_back(i);
            }
        }

        std::vector<double> arr(arrSize, 1.0);

        for (const auto &e : mstIBD) {
            double pro[9];
            for (int k = 0; k < 9; ++k)
                pro[k] = e.pa * pa_arr[k] + e.pb * pb_arr[k] + e.pc * pc_arr[k];

            for (int idx = 0; idx < arrSize; ++idx) {
                int g1 = (idx / stride[e.idx1]) % 3;
                int g2 = (idx / stride[e.idx2]) % 3;
                arr[idx] *= pro[g1 * 3 + g2];
            }
        }

        for (int node : divideNodes) {
            for (int idx = 0; idx < arrSize; ++idx) {
                int g = (idx / stride[node]) % 3;
                arr[idx] /= p0[g];
            }
        }

        for (int idx = 0; idx < arrSize; ++idx)
            CLT(idx, mi) = arr[idx];
    }

    return CLT;
}

// ══════════════════════════════════════════════════════════════════════
// buildStandS
// ══════════════════════════════════════════════════════════════════════

std::vector<double> buildStandS(
    int N,
    const std::vector<double> &resid
) {
    int arrSize = 1;
    for (int i = 0; i < N; ++i)
        arrSize *= 3;

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
// buildSPAGRMNullModel
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
    int nthreads
) {
    // ── Build dynamic MAF grid from QC cutoffs ───────────────────────
    const std::vector<double> mafInterval = buildMafInterval(minMafCutoff, minMacCutoff, N);
    infoMsg(
        "  MAF interval grid (%zu bins): [%.2g .. %.2g]",
        mafInterval.size(),
        mafInterval.front(),
        mafInterval.back()
    );

    // ── Outlier detection (IQR) ──────────────────────────────────────
    std::vector<double> sortedResid(N);
    for (uint32_t i = 0; i < N; ++i)
        sortedResid[i] = Resid[i];
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

    infoMsg("  Outlier detection: %d outliers (%.1f%%) from %u subjects", nOutlier, 100.0 * nOutlier / N, N);

    // ── Accumulate global variance terms for singletons ──────────────
    double R_GRM_R = 0.0;
    double sum_R_nonOutlier = 0.0;
    double R_GRM_R_nonOutlier = 0.0;

    for (uint32_t idx : singletonSet) {
        double diagVal = grmDiag[idx];
        double contrib = diagVal * Resid[idx] * Resid[idx];
        R_GRM_R += contrib;
        if (!isOutlier[idx]) R_GRM_R_nonOutlier += contrib;
    }
    for (uint32_t idx : singletonSet) {
        if (!isOutlier[idx]) sum_R_nonOutlier += Resid[idx];
    }

    std::vector<double> unrelatedOutlierResids;
    for (uint32_t idx : singletonSet) {
        if (isOutlier[idx]) unrelatedOutlierResids.push_back(Resid[idx]);
    }

    double R_GRM_R_TwoSubjOutlier = 0.0;

    std::vector<std::array<double, 2> > twoSubj_resid_list;
    std::vector<std::vector<double> > twoSubj_rho_list;
    std::vector<std::vector<double> > threeSubj_standS_list;
    std::vector<Eigen::MatrixXd> threeSubj_CLT_list;

    // ── Process families (multi-threaded) ──────────────────────────────
    struct FamAccum {
        double r_grm_r = 0.0;
        double sum_r_no = 0.0;
        double r_grm_r_no = 0.0;
        size_t nSplit = 0;
        std::vector<std::vector<uint32_t> > outlierFams;
    };

    const int nWorkers = std::max(1, std::min(nthreads, static_cast<int>(families.size())));
    std::vector<FamAccum> accums(nWorkers);
    std::atomic<size_t> nextFam{0};

    auto familyWorker = [&](int tid) {
        auto &acc = accums[tid];
        while (true) {
            const size_t fi = nextFam.fetch_add(1, std::memory_order_relaxed);
            if (fi >= families.size()) break;

            const auto &fam = families[fi];

            bool hasOutlier = false;
            for (uint32_t idx : fam)
                if (isOutlier[idx]) {
                    hasOutlier = true;
                    break;
                }

            // Direct quad form — no hash set needed, entries belong to this family
            double famQuad = 0.0;
            for (const auto &e : familyEntries[fi]) {
                double factor = (e.row == e.col) ? 1.0 : 2.0;
                famQuad += factor * e.value * Resid[e.row] * Resid[e.col];
            }
            acc.r_grm_r += famQuad;

            if (!hasOutlier) {
                double famSum = 0.0;
                for (uint32_t idx : fam)
                    famSum += Resid[idx];
                acc.sum_r_no += famSum;
                acc.r_grm_r_no += famQuad;
                continue;
            }

            if (static_cast<int>(fam.size()) <= MAX_NUM_IN_FAM) {
                acc.outlierFams.push_back(fam);
                continue;
            }

            // ── Greedy family splitting (binary search + Union-Find) ─────
            ++acc.nSplit;
            const uint32_t famN = static_cast<uint32_t>(fam.size());
            std::unordered_map<uint32_t, uint32_t> g2l;
            g2l.reserve(famN);
            for (uint32_t i = 0; i < famN; ++i)
                g2l[fam[i]] = i;

            struct OffDiagEntry {
                uint32_t localRow, localCol;
                double value, cov;
            };

            std::vector<OffDiagEntry> famEdges;
            for (const auto &e : familyEntries[fi]) {
                if (e.row == e.col) continue;
                if (e.row < e.col) famEdges.push_back(
                        {g2l[e.row], g2l[e.col], e.value,
                         std::abs(e.value * Resid[e.row] * Resid[e.col])}
                );
            }
            std::sort(
                famEdges.begin(),
                famEdges.end(),
                [](const OffDiagEntry &a, const OffDiagEntry &b) {
                return a.cov < b.cov;
            }
            );

            const int nEdges = static_cast<int>(famEdges.size());

            auto maxCompKeeping = [&](int j) -> int {
                UnionFind uf(famN);
                for (int k = j + 1; k < nEdges; ++k)
                    uf.unite(famEdges[k].localRow, famEdges[k].localCol);
                uint32_t maxSz = 0;
                for (uint32_t v = 0; v < famN; ++v)
                    maxSz = std::max(maxSz, uf.componentSize(v));
                return static_cast<int>(maxSz);
            };

            int removeUpTo = -1;
            if (maxCompKeeping(-1) > MAX_NUM_IN_FAM) {
                int lo = 0, hi = nEdges - 1;
                while (lo <= hi) {
                    int mid = lo + (hi - lo) / 2;
                    if (maxCompKeeping(mid) <= MAX_NUM_IN_FAM) {
                        removeUpTo = mid;
                        hi = mid - 1;
                    } else {
                        lo = mid + 1;
                    }
                }
            }

            std::vector<bool> edgeKept(nEdges, false);
            {
                UnionFind uf(famN);
                for (int k = removeUpTo + 1; k < nEdges; ++k) {
                    uf.unite(famEdges[k].localRow, famEdges[k].localCol);
                    edgeKept[k] = true;
                }
                for (int k = removeUpTo; k >= 0; --k) {
                    uint32_t ra = uf.find(famEdges[k].localRow);
                    uint32_t rb = uf.find(famEdges[k].localCol);
                    if (ra == rb) {
                        edgeKept[k] = true;
                    } else if (uf.componentSize(ra) + uf.componentSize(rb) <=
                               static_cast<uint32_t>(MAX_NUM_IN_FAM)) {
                        uf.unite(ra, rb);
                        edgeKept[k] = true;
                    }
                }
            }

            {
                UnionFind uf(famN);
                for (int k = 0; k < nEdges; ++k)
                    if (edgeKept[k]) uf.unite(famEdges[k].localRow, famEdges[k].localCol);

                // Compute all sub-family quad forms in one pass — O(E) total
                std::unordered_map<uint32_t, double> compQuad;
                for (const auto &e : familyEntries[fi]) {
                    uint32_t lr = g2l[e.row], lc = g2l[e.col];
                    uint32_t rr = uf.find(lr), rc = uf.find(lc);
                    if (rr != rc) continue;                     // cross-component edge, removed
                    double factor = (e.row == e.col) ? 1.0 : 2.0;
                    compQuad[rr] += factor * e.value * Resid[e.row] * Resid[e.col];
                }

                std::unordered_map<uint32_t, std::vector<uint32_t> > compMap;
                for (uint32_t v = 0; v < famN; ++v)
                    compMap[uf.find(v)].push_back(v);

                for (auto &[root, localMembers] : compMap) {
                    std::vector<uint32_t> subFam;
                    subFam.reserve(localMembers.size());
                    for (uint32_t li : localMembers)
                        subFam.push_back(fam[li]);

                    bool subHasOutlier = false;
                    for (uint32_t idx : subFam)
                        if (isOutlier[idx]) {
                            subHasOutlier = true;
                            break;
                        }

                    if (!subHasOutlier) {
                        double subSum = 0.0;
                        for (uint32_t idx : subFam)
                            subSum += Resid[idx];
                        acc.sum_r_no += subSum;
                        acc.r_grm_r_no += compQuad[root];
                    } else {
                        acc.outlierFams.push_back(std::move(subFam));
                    }
                }
            }
        }
    };

    {
        std::vector<std::thread> threads;
        for (int t = 1; t < nWorkers; ++t)
            threads.emplace_back(familyWorker, t);
        familyWorker(0);
        for (auto &t : threads)
            t.join();
    }

    // Merge thread-local accumulators
    std::vector<std::vector<uint32_t> > outlierFamilies;
    size_t nSplitFamilies = 0;
    for (auto &acc : accums) {
        R_GRM_R += acc.r_grm_r;
        sum_R_nonOutlier += acc.sum_r_no;
        R_GRM_R_nonOutlier += acc.r_grm_r_no;
        nSplitFamilies += acc.nSplit;
        for (auto &fam : acc.outlierFams)
            outlierFamilies.push_back(std::move(fam));
    }

    infoMsg(
        "  %zu singletons, %zu families (%zu needed splitting)",
        singletonSet.size(),
        families.size(),
        nSplitFamilies
    );
    infoMsg("  %zu outlier families to process", outlierFamilies.size());

    // ── Build Chow-Liu trees for outlier families ────────────────────

    // Build GRM off-diagonal map for O(1) pair lookups
    std::unordered_map<uint64_t, double> grmOffDiagMap;
    grmOffDiagMap.reserve(allGrmEntries.size());
    for (const auto &e : allGrmEntries) {
        if (e.row != e.col) {
            uint32_t lo = std::min(e.row, e.col);
            uint32_t hi = std::max(e.row, e.col);
            grmOffDiagMap[(static_cast<uint64_t>(lo) << 32) | hi] = e.value;
        }
    }

    for (const auto &fam : outlierFamilies) {
        const int n1 = static_cast<int>(fam.size());

        if (n1 == 1) {
            unrelatedOutlierResids.push_back(Resid[fam[0]]);
            continue;
        }

        if (n1 == 2) {
            double R1 = Resid[fam[0]], R2 = Resid[fam[1]];

            // Direct quad form: diag[a]*R1² + diag[b]*R2² + 2*grm(a,b)*R1*R2
            uint32_t lo = std::min(fam[0], fam[1]);
            uint32_t hi = std::max(fam[0], fam[1]);
            uint64_t key = (static_cast<uint64_t>(lo) << 32) | hi;

            double offDiag = 0.0;
            auto grmIt = grmOffDiagMap.find(key);
            if (grmIt != grmOffDiagMap.end()) offDiag = grmIt->second;
            double pairQuad = grmDiag[fam[0]] * R1 * R1 + grmDiag[fam[1]] * R2 * R2 + 2.0 * offDiag * R1 * R2;
            R_GRM_R_TwoSubjOutlier += pairQuad;

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
        std::vector<double> famResid(n1);
        for (int i = 0; i < n1; ++i) {
            famResid[i] = Resid[fam[i]];
        }

        // Enumerate all pairs in the family for IBD lookup (at most C(5,2)=10)
        std::vector<IndexedIBD> famIBD;
        for (int i = 0; i < n1; ++i) {
            for (int j = i + 1; j < n1; ++j) {
                uint32_t lo = std::min(fam[i], fam[j]);
                uint32_t hi = std::max(fam[i], fam[j]);
                uint64_t key = (static_cast<uint64_t>(lo) << 32) | hi;
                auto it = ibdPairMap.find(key);
                if (it != ibdPairMap.end()) famIBD.push_back(ibdEntries[it->second]);
            }
        }

        Eigen::MatrixXd CLT = buildChowLiuTree(n1, famIBD, fam, mafInterval);
        std::vector<double> standS = buildStandS(n1, famResid);

        threeSubj_standS_list.push_back(std::move(standS));
        threeSubj_CLT_list.push_back(std::move(CLT));
    }

    // ── Assemble SPAGRMClass ─────────────────────────────────────────
    infoMsg(
        "  Built %zu two-subject, %zu multi-subject family models, "
        "%zu unrelated outliers",
        twoSubj_resid_list.size(),
        threeSubj_CLT_list.size(),
        unrelatedOutlierResids.size()
    );

    Eigen::VectorXd residOutliers = Eigen::Map<Eigen::VectorXd>(
        unrelatedOutlierResids.data(),
        static_cast<Eigen::Index>(unrelatedOutlierResids.size())
    );

    nsSPAGRM::FamilyData fd;
    fd.resid_unrelated_outliers = std::move(residOutliers);
    fd.twoSubj_resid = std::move(twoSubj_resid_list);
    fd.twoSubj_rho = std::move(twoSubj_rho_list);
    fd.threeSubj_standS = std::move(threeSubj_standS_list);
    fd.threeSubj_CLT = std::move(threeSubj_CLT_list);

    return SPAGRMClass(
        Resid,
        sum_R_nonOutlier,
        R_GRM_R_nonOutlier,
        R_GRM_R_TwoSubjOutlier,
        R_GRM_R,
        mafInterval,
        std::move(fd),
        spaCutoff,
        ZETA_DEFAULT,
        TOL_DEFAULT
    );
}

} // namespace nsGRMNull
