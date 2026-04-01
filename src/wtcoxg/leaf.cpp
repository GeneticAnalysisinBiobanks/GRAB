// leaf.cpp — LEAF full implementation

#include "wtcoxg/leaf.hpp"
#include "wtcoxg/wtcoxg.hpp"
#include "io/plink.hpp"
#include "io/subject_data.hpp"
#include "io/sparse_grm.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>
#include <cctype>

// ======================================================================
// loadAndMatchRefAf — load one plink2 .afreq and match vs PlinkData
// ======================================================================

std::vector<PopMatchedAF> loadAndMatchRefAf(
    const PlinkData& plinkData,
    const std::string& afreqFile) {

  auto records = loadRefAfFile(afreqFile);

  // Build lookup: "chrom:id" → index
  auto makeKey = [](const std::string& chr, const std::string& id) -> std::string {
    std::string k;
    k.reserve(chr.size() + 1 + id.size());
    k += chr; k += ':'; k += id;
    return k;
  };

  std::unordered_map<std::string, size_t> refMap;
  refMap.reserve(records.size());
  for (size_t i = 0; i < records.size(); ++i)
    refMap.emplace(makeKey(records[i].chrom, records[i].id), i);

  std::vector<PopMatchedAF> matched;
  matched.reserve(plinkData.markerInfo().size());
  for (const auto& mi : plinkData.markerInfo()) {
    auto key = makeKey(mi.chrom, mi.id);
    auto it = refMap.find(key);
    if (it == refMap.end()) continue;

    const auto& rec = records[it->second];
    PopMatchedAF pm;
    pm.genoIndex = mi.genoIndex;

    if (rec.alt_allele == mi.ref && rec.ref_allele == mi.alt) {
      pm.af = rec.alt_freq;
    } else if (rec.ref_allele == mi.ref && rec.alt_allele == mi.alt) {
      pm.af = 1.0 - rec.alt_freq;
    } else {
      continue;
    }
    pm.n_ref = rec.obs_ct / 2.0;
    matched.push_back(pm);
  }
  return matched;
}


// ======================================================================
// Summix — exact active-set enumeration with KKT / Lagrange multiplier
//
// For each non-empty subset S ⊆ {0..K-1} (K ≤ 6 → max 63 subsets),
// solve the equality-constrained LS:
//   min ||D_S * p_S - o||²   s.t.  1'p_S = 1
// via the (|S|+1) × (|S|+1) KKT system:
//   [2*D_S'D_S  1] [p_S  ]   [2*D_S'o]
//   [1'         0] [lambda] = [1      ]
// Keep the solution with all p_S ≥ 0 and lowest objective.
// ======================================================================

Eigen::VectorXd summixEstimate(
    const Eigen::VectorXd& observedAF,
    const Eigen::MatrixXd& refAF) {

  const int K = static_cast<int>(refAF.cols());
  if (K == 0)
    return Eigen::VectorXd();
  if (K > 6)
    throw std::runtime_error("summix: nPop > 6 not supported (got "
                             + std::to_string(K) + ")");

  // Filter out rows with NaN
  std::vector<Eigen::Index> validRows;
  validRows.reserve(observedAF.size());
  for (Eigen::Index i = 0; i < observedAF.size(); ++i) {
    if (std::isnan(observedAF[i])) continue;
    bool ok = true;
    for (int p = 0; p < K; ++p)
      if (std::isnan(refAF(i, p))) { ok = false; break; }
    if (ok) validRows.push_back(i);
  }

  if (validRows.empty())
    return Eigen::VectorXd::Constant(K, 1.0 / K);

  const int nValid = static_cast<int>(validRows.size());
  Eigen::VectorXd obs(nValid);
  Eigen::MatrixXd D(nValid, K);
  for (int i = 0; i < nValid; ++i) {
    obs[i] = observedAF[validRows[i]];
    D.row(i) = refAF.row(validRows[i]);
  }

  // Precompute full D'D and D'o
  Eigen::MatrixXd DtD = D.transpose() * D;   // K×K
  Eigen::VectorXd Dto = D.transpose() * obs;  // K
  double oTo = obs.squaredNorm();

  Eigen::VectorXd bestP = Eigen::VectorXd::Constant(K, 1.0 / K);
  double bestObj = std::numeric_limits<double>::infinity();

  // Enumerate all non-empty subsets of {0..K-1}
  const int nSubsets = (1 << K);
  for (int mask = 1; mask < nSubsets; ++mask) {
    // Count and collect active indices
    int s = 0;
    int idx[6];
    for (int j = 0; j < K; ++j)
      if (mask & (1 << j)) idx[s++] = j;

    // Build the (s+1) × (s+1) KKT system
    // [2*A   e] [p_S  ]   [2*b]
    // [e'   0] [lambda] = [1  ]
    // where A = D_S'D_S (s×s), b = D_S'o (s), e = ones(s)
    const int dim = s + 1;
    Eigen::MatrixXd KKT = Eigen::MatrixXd::Zero(dim, dim);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(dim);

    for (int i = 0; i < s; ++i) {
      for (int j = 0; j < s; ++j)
        KKT(i, j) = 2.0 * DtD(idx[i], idx[j]);
      KKT(i, s) = 1.0;   // e column
      KKT(s, i) = 1.0;   // e' row
      rhs[i] = 2.0 * Dto[idx[i]];
    }
    rhs[s] = 1.0;  // sum constraint

    // Solve
    Eigen::VectorXd sol = KKT.partialPivLu().solve(rhs);

    // Check all p_S >= 0
    bool feasible = true;
    for (int i = 0; i < s; ++i) {
      if (sol[i] < -1e-12) { feasible = false; break; }
    }
    if (!feasible) continue;

    // Clamp tiny negatives to 0
    for (int i = 0; i < s; ++i)
      if (sol[i] < 0.0) sol[i] = 0.0;

    // Compute objective: ||D_S p_S - o||²  = p_S' A p_S - 2 b' p_S + o'o
    Eigen::VectorXd pS = sol.head(s);
    double obj = 0.0;
    for (int i = 0; i < s; ++i)
      for (int j = 0; j < s; ++j)
        obj += pS[i] * DtD(idx[i], idx[j]) * pS[j];
    for (int i = 0; i < s; ++i)
      obj -= 2.0 * Dto[idx[i]] * pS[i];
    obj += oTo;

    if (obj < bestObj) {
      bestObj = obj;
      bestP.setZero();
      for (int i = 0; i < s; ++i)
        bestP[idx[i]] = pS[i];
    }
  }

  return bestP;
}


// ======================================================================
// LEAFMethod — MethodBase implementation
// ======================================================================

LEAFMethod::LEAFMethod(
    std::vector<std::unique_ptr<WtCoxGMethod>> clusterMethods,
    std::vector<std::vector<uint32_t>>         clusterIndices)
  : m_nCluster(static_cast<int>(clusterMethods.size())),
    m_clusterMethods(std::move(clusterMethods)),
    m_clusterIndices(std::move(clusterIndices))
{
  m_clusterGVec.resize(m_nCluster);
  for (int c = 0; c < m_nCluster; ++c)
    m_clusterGVec[c].resize(m_clusterIndices[c].size());
}

std::unique_ptr<MethodBase> LEAFMethod::clone() const {
  std::vector<std::unique_ptr<WtCoxGMethod>> cloned;
  cloned.reserve(m_nCluster);
  for (const auto& m : m_clusterMethods) {
    auto p = m->clone();  // returns unique_ptr<MethodBase>
    cloned.push_back(std::unique_ptr<WtCoxGMethod>(
        static_cast<WtCoxGMethod*>(p.release())));
  }
  return std::make_unique<LEAFMethod>(std::move(cloned), m_clusterIndices);
}

int LEAFMethod::resultSize() const {
  return 2 + 3 * m_nCluster;
}

std::string LEAFMethod::getHeaderColumns() const {
  std::ostringstream oss;
  oss << "\tmeta.p_ext\tmeta.p_noext";
  for (int i = 1; i <= m_nCluster; ++i)
    oss << "\tcl" << i << ".p_ext\tcl" << i << ".p_noext\tcl" << i << ".p_batch";
  return oss.str();
}

void LEAFMethod::prepareChunk(const std::vector<uint64_t>& gIndices) {
  for (auto& m : m_clusterMethods)
    m->prepareChunk(gIndices);
}

void LEAFMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> GVec,
    double /*altFreq*/, int markerInChunkIdx,
    bool /*flipped*/, std::vector<double>& result) {

  std::vector<double> pExt(m_nCluster), pNoext(m_nCluster);
  std::vector<double> sExt(m_nCluster), sNoext(m_nCluster);

  for (int c = 0; c < m_nCluster; ++c) {
    // Gather cluster genotypes from full GVec
    const auto& idx = m_clusterIndices[c];
    Eigen::VectorXd& gClu = m_clusterGVec[c];
    for (size_t k = 0; k < idx.size(); ++k)
      gClu[static_cast<Eigen::Index>(k)] = GVec[idx[k]];

    auto dr = m_clusterMethods[c]->computeDual(gClu, markerInChunkIdx);
    pExt[c]   = dr.p_ext;
    pNoext[c] = dr.p_noext;
    sExt[c]   = dr.score_ext;
    sNoext[c] = dr.score_noext;
  }

  // Fixed-effects meta-analysis: pool scores across clusters
  auto metaP = [](const std::vector<double>& scores,
                  const std::vector<double>& pvals) -> double {
    double sumScore = 0.0, sumVar = 0.0;
    for (size_t c = 0; c < scores.size(); ++c) {
      if (std::isnan(scores[c]) || std::isnan(pvals[c])
          || pvals[c] <= 0.0 || pvals[c] >= 1.0) continue;
      double chisq = math::qchisq(pvals[c], 1.0, false, false);
      if (chisq < 1e-30) chisq = 1e-30;
      double var = (scores[c] * scores[c]) / chisq;
      if (std::isnan(var)) continue;
      sumScore += scores[c];
      sumVar   += var;
    }
    if (sumVar <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    double z = sumScore / std::sqrt(sumVar);
    return std::erfc(std::fabs(z) / std::sqrt(2.0));  // two-sided p-value
  };

  result.push_back(metaP(sExt, pExt));
  result.push_back(metaP(sNoext, pNoext));
  for (int c = 0; c < m_nCluster; ++c) {
    result.push_back(pExt[c]);
    result.push_back(pNoext[c]);
    result.push_back(m_clusterMethods[c]->chunkRefInfoAt(markerInChunkIdx).pvalue_bat);
  }
}


// ======================================================================
// runLEAF — full orchestration
// ======================================================================

void runLEAF(
    const std::vector<std::string>& residFiles,
    const std::string& bfilePrefix,
    const std::vector<std::string>& refAfFiles,
    const std::string& spgrmSaigeFile,
    const std::string& spgrmGctaPrefix,
    const std::string& outputFile,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff) {

  const int nCluster = static_cast<int>(residFiles.size());

  // ---- Load per-cluster resid files via SubjectData ----
  infoMsg("Loading %d cluster resid files...", nCluster);
  auto famIIDs = parseFamIIDs(bfilePrefix + ".fam");
  const uint32_t nFamTotal = static_cast<uint32_t>(famIIDs.size());

  std::vector<SubjectData> clusterSD;
  clusterSD.reserve(nCluster);
  for (int c = 0; c < nCluster; ++c) {
    infoMsg("  Cluster %d: %s", c + 1, residFiles[c].c_str());
    auto fc = famIIDs;  // copy for each cluster
    clusterSD.emplace_back(std::move(fc));
    clusterSD.back().loadResidWtCoxG(residFiles[c]);
    clusterSD.back().finalize();
    infoMsg("    %u subjects loaded", clusterSD.back().nUsed());
  }

  // Build union bitmask (OR of all cluster masks) and cluster indices
  const size_t nMaskWords = (nFamTotal + 63) / 64;
  std::vector<uint64_t> unionMask(nMaskWords, 0);
  for (int c = 0; c < nCluster; ++c) {
    const auto& cmask = clusterSD[c].usedMask();
    for (size_t w = 0; w < nMaskWords; ++w)
      unionMask[w] |= cmask[w];
  }
  uint32_t nUnion = 0;
  for (size_t w = 0; w < nMaskWords; ++w)
    nUnion += static_cast<uint32_t>(__builtin_popcountll(unionMask[w]));

  // Build per-cluster indices into the union-dense genotype vector
  std::vector<std::vector<uint32_t>> clusterIndices(nCluster);
  for (int c = 0; c < nCluster; ++c)
    clusterIndices[c].reserve(clusterSD[c].nUsed());

  uint32_t denseIdx = 0;
  for (size_t w = 0; w < nMaskWords; ++w) {
    uint64_t bits = unionMask[w];
    while (bits) {
      int bit = __builtin_ctzll(bits);
      uint64_t bitVal = uint64_t(1) << bit;
      for (int c = 0; c < nCluster; ++c) {
        if (clusterSD[c].usedMask()[w] & bitVal) {
          clusterIndices[c].push_back(denseIdx);
          break;
        }
      }
      ++denseIdx;
      bits &= bits - 1;
    }
  }
  infoMsg("  Total subjects (union): %u", nUnion);

  // ---- Load PLINK data with union bitmask ----
  infoMsg("Loading PLINK data: %s", bfilePrefix.c_str());
  PlinkData plinkData(
      bfilePrefix + ".bed",
      bfilePrefix + ".bim",
      bfilePrefix + ".fam",
      unionMask,
      nFamTotal,
      nUnion,
      "ref-first",
      {}, {}, {}, {},
      nSnpPerChunk);
  infoMsg("  %u subjects matched, %u markers",
          plinkData.nSubjUsed(), plinkData.nMarkers());

  // ---- Load per-population ref-af files (plink2 .afreq) ----
  const int nPop = static_cast<int>(refAfFiles.size());
  infoMsg("Loading %d reference AF files...", nPop);
  std::vector<std::vector<PopMatchedAF>> popMatched(nPop);
  for (int p = 0; p < nPop; ++p) {
    infoMsg("  Pop %d: %s", p + 1, refAfFiles[p].c_str());
    popMatched[p] = loadAndMatchRefAf(plinkData, refAfFiles[p]);
    infoMsg("    %zu markers matched", popMatched[p].size());
  }

  // ---- Build intersection of markers present in all populations ----
  // Collect genoIndex → {popAF[nPop], popN[nPop]} for markers in all pops
  struct MultiPopEntry {
    double af[6];
    double n_ref[6];
  };
  std::unordered_map<uint64_t, MultiPopEntry> allPopMap;
  // Seed from pop 0
  if (nPop > 0) {
    for (const auto& pm : popMatched[0]) {
      MultiPopEntry e{};
      e.af[0] = pm.af;
      e.n_ref[0] = pm.n_ref;
      allPopMap.emplace(pm.genoIndex, e);
    }
  }
  // Intersect with remaining pops
  for (int p = 1; p < nPop; ++p) {
    std::unordered_map<uint64_t, MultiPopEntry> next;
    for (const auto& pm : popMatched[p]) {
      auto it = allPopMap.find(pm.genoIndex);
      if (it == allPopMap.end()) continue;
      auto entry = it->second;
      entry.af[p] = pm.af;
      entry.n_ref[p] = pm.n_ref;
      next.emplace(pm.genoIndex, entry);
    }
    allPopMap = std::move(next);
  }

  // Convert to sorted vector for deterministic order
  struct MatchedMultiPop {
    uint64_t genoIndex;
    double popAF[6];
    double popN[6];
  };
  std::vector<MatchedMultiPop> matchedMulti;
  matchedMulti.reserve(allPopMap.size());
  for (auto& [gi, e] : allPopMap) {
    MatchedMultiPop mm{};
    mm.genoIndex = gi;
    for (int p = 0; p < nPop; ++p) {
      mm.popAF[p] = e.af[p];
      mm.popN[p]  = e.n_ref[p];
    }
    matchedMulti.push_back(mm);
  }
  std::sort(matchedMulti.begin(), matchedMulti.end(),
            [](const auto& a, const auto& b) { return a.genoIndex < b.genoIndex; });
  const size_t nMatched = matchedMulti.size();
  infoMsg("  %zu markers present in all %d populations", nMatched, nPop);

  // ---- Per-cluster summix + AF synthesis ----
  infoMsg("Per-cluster ancestry estimation and AF synthesis...");

  // Scan genotypes once for all matched markers → full GVec
  // Then compute per-cluster internal AF
  PlinkCursor cursor(plinkData.bedFile(),
                     static_cast<uint32_t>(plinkData.nMarkers()),
                     plinkData.nSubjInFile(),
                     plinkData.usedMask(),
                     plinkData.nSubjUsed(),
                     plinkData.isAltFirst(),
                     plinkData.allUsed());

  const uint32_t nFull = plinkData.nSubjUsed();
  Eigen::VectorXd fullGVec(nFull);

  // Pre-position the cursor at the first matched marker so that
  // the block-buffered sequential reader is used for the scan.
  if (!matchedMulti.empty())
    cursor.beginSequentialBlock(matchedMulti.front().genoIndex);

  // Per-cluster internal AF matrix: nMatched × nCluster
  // Per-cluster case/control stats
  struct ClMarkerStat {
    double intAF;
    double mu0, mu1, n0, n1, mu_int;
  };
  std::vector<std::vector<ClMarkerStat>> clStats(nCluster,
      std::vector<ClMarkerStat>(nMatched));

  infoMsg("  Scanning %zu matched markers for per-cluster allele frequencies...", nMatched);
  for (size_t m = 0; m < nMatched; ++m) {
    // getGenotypesSimple: no QC stats, missing → NaN, no indexForMissing overhead
    cursor.getGenotypesSimple(matchedMulti[m].genoIndex, fullGVec);

    for (int c = 0; c < nCluster; ++c) {
      const auto& idx = clusterIndices[c];
      const auto& ind = clusterSD[c].indicator();
      double sum0 = 0, sum1 = 0, cnt0 = 0, cnt1 = 0;
      for (size_t k = 0; k < idx.size(); ++k) {
        double g = fullGVec[idx[k]];
        if (std::isnan(g)) continue;
        if (ind[static_cast<Eigen::Index>(k)] == 1.0) {
          sum1 += g; cnt1 += 1.0;
        } else {
          sum0 += g; cnt0 += 1.0;
        }
      }
      double total = cnt0 + cnt1;
      auto& st = clStats[c][m];
      st.intAF  = (total > 0) ? (sum0 + sum1) / (2.0 * total) : std::numeric_limits<double>::quiet_NaN();
      st.mu0    = (cnt0 > 0) ? sum0 / (2.0 * cnt0) : 0.0;
      st.mu1    = (cnt1 > 0) ? sum1 / (2.0 * cnt1) : 0.0;
      st.n0     = cnt0;
      st.n1     = cnt1;
      st.mu_int = st.intAF;
    }
  }

  // Per cluster: run summix → compute AF_ref, AN_ref
  // Build per-cluster MatchedMarkerInfo for testBatchEffects
  std::vector<std::vector<MatchedMarkerInfo>> clMatchedInfo(nCluster);

  for (int c = 0; c < nCluster; ++c) {
    infoMsg("  Cluster %d: running summix...", c + 1);

    // Build observed AF and reference AF matrix for summix
    Eigen::VectorXd obsAF(nMatched);
    Eigen::MatrixXd refMat(nMatched, nPop);
    for (size_t m = 0; m < nMatched; ++m) {
      obsAF[m] = clStats[c][m].intAF;
      for (int p = 0; p < nPop; ++p)
        refMat(m, p) = matchedMulti[m].popAF[p];
    }

    Eigen::VectorXd proportions = summixEstimate(obsAF, refMat);
    for (int p = 0; p < nPop; ++p)
      infoMsg("    pop%d: %.4f", p + 1, proportions[p]);

    // Synthesise ancestry-matched AF_ref and AN_ref
    clMatchedInfo[c].resize(nMatched);
    for (size_t m = 0; m < nMatched; ++m) {
      auto& mi = clMatchedInfo[c][m];
      mi.genoIndex = matchedMulti[m].genoIndex;

      // AF_ref = Σ(pi_k * AF_k)
      double af_ref = 0.0;
      for (int p = 0; p < nPop; ++p)
        af_ref += proportions[p] * matchedMulti[m].popAF[p];
      mi.AF_ref = af_ref;

      // N_ref = 1 / Σ(pi_k² / N_k)  (effective sample size for admixed ancestry)
      double denom = 0.0;
      for (int p = 0; p < nPop; ++p) {
        double N = matchedMulti[m].popN[p];
        if (N > 0 && proportions[p] > 0)
          denom += (proportions[p] * proportions[p]) / N;
      }
      mi.N_ref = (denom > 0) ? 1.0 / denom : 0.0;

      mi.mu0    = clStats[c][m].mu0;
      mi.mu1    = clStats[c][m].mu1;
      mi.n0     = clStats[c][m].n0;
      mi.n1     = clStats[c][m].n1;
      mi.mu_int = clStats[c][m].mu_int;
    }
  }

  // ---- Per-cluster batch-effect testing ----
  infoMsg("Per-cluster batch-effect testing...");

  std::vector<std::shared_ptr<std::unordered_map<uint64_t, WtCoxGRefInfo>>>
      clRefMaps(nCluster);

  for (int c = 0; c < nCluster; ++c) {
    infoMsg("  Cluster %d: batch-effect testing (%zu markers)...",
            c + 1, clMatchedInfo[c].size());

    std::unique_ptr<SparseGRM> grm;
    if (!spgrmGctaPrefix.empty()) {
      grm = std::make_unique<SparseGRM>(
          SparseGRM::fromGCTA(spgrmGctaPrefix, clusterSD[c].usedIIDs()));
      infoMsg("    Sparse GRM: %u subjects, %zu non-zeros",
              grm->nSubjects(), grm->nnz());
    } else if (!spgrmSaigeFile.empty()) {
      grm = std::make_unique<SparseGRM>(spgrmSaigeFile, clusterSD[c].usedIIDs());
      infoMsg("    Sparse GRM: %u subjects, %zu non-zeros",
              grm->nSubjects(), grm->nnz());
    }

    clRefMaps[c] = testBatchEffects(
        clMatchedInfo[c], clusterSD[c].residuals(), clusterSD[c].weights(),
        clusterSD[c].indicator(), grm.get(), refPrevalence, cutoff);
    infoMsg("    %zu markers retained", clRefMaps[c]->size());
  }

  // ---- Build LEAFMethod and run marker engine ----
  infoMsg("Building LEAF method (%d clusters)...", nCluster);

  std::vector<std::unique_ptr<WtCoxGMethod>> clMethods;
  clMethods.reserve(nCluster);
  for (int c = 0; c < nCluster; ++c) {
    clMethods.push_back(std::make_unique<WtCoxGMethod>(
        clusterSD[c].residuals(),
        clusterSD[c].weights(),
        cutoff,
        spaCutoff,
        clRefMaps[c]));
  }

  LEAFMethod method(std::move(clMethods), clusterIndices);

  infoMsg("Running LEAF marker engine (%d thread(s))...", nthread);
  markerEngine(plinkData, method, outputFile,
               nthread,
               missingCutoff,
               minMafCutoff,
               minMacCutoff,
               /*exactHwe=*/false);
}
