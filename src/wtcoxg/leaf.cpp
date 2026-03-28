// leaf.cpp — LEAF full implementation

#include "wtcoxg/leaf.hpp"
#include "wtcoxg/wtcoxg.hpp"
#include "io/plink.hpp"
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
// Helper utilities
// ======================================================================

namespace {

// Split a line on whitespace (tabs and spaces)
std::vector<std::string> splitWS(const std::string& line) {
  std::vector<std::string> tokens;
  std::istringstream iss(line);
  std::string tok;
  while (iss >> tok) tokens.push_back(std::move(tok));
  return tokens;
}

} // anon namespace


// ======================================================================
// Multi-population ref-af file parser (6+2N columns, no header)
// Format: chr  id  cm  bp  a1  a2  AF1  AN1  AF2  AN2  ...
// Same first 6 columns as .bim. Lines starting with '#' are skipped.
// ======================================================================

MultiPopRefAf loadMultiPopRefAfFile(const std::string& filename) {
  std::ifstream ifs(filename);
  if (!ifs) throw std::runtime_error("Cannot open ref-af file: " + filename);

  MultiPopRefAf result;
  result.records.reserve(100000);

  std::string line;
  bool nPopDetected = false;

  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;
    auto fields = splitWS(line);
    if (fields.size() < 8) continue;  // at least 6 bim + 1 AF/AN pair

    if (!nPopDetected) {
      int extra = static_cast<int>(fields.size()) - 6;
      if (extra < 2 || extra % 2 != 0)
        throw std::runtime_error("ref-af file must have 6+2N columns (got "
                                 + std::to_string(fields.size()) + "): " + filename);
      result.nPop = extra / 2;
      nPopDetected = true;
    }

    const int nPop = result.nPop;
    if (static_cast<int>(fields.size()) < 6 + 2 * nPop) continue;

    MultiPopRefAf::Record rec;
    rec.chrom = fields[0];
    rec.id    = fields[1];
    rec.cm    = std::strtod(fields[2].c_str(), nullptr);
    rec.pos   = static_cast<uint32_t>(std::strtoul(fields[3].c_str(), nullptr, 10));
    rec.a1    = fields[4];
    rec.a2    = fields[5];
    rec.popAF.resize(nPop);
    rec.popAN.resize(nPop);
    for (int p = 0; p < nPop; ++p) {
      rec.popAF[p] = std::strtod(fields[6 + 2 * p].c_str(), nullptr);
      rec.popAN[p] = std::strtod(fields[6 + 2 * p + 1].c_str(), nullptr);
    }
    result.records.push_back(std::move(rec));
  }
  return result;
}


// ======================================================================
// Marker matching — exact match on (chr, id, cm, bp, a1, a2), no flip
// ======================================================================

std::vector<MultiPopMatchedMarker> matchMultiPopMarkers(
    const PlinkData& plinkData,
    const MultiPopRefAf& refAf) {

  // Build ref lookup: key = "chr:bp:a1:a2"
  auto makeKey = [](const std::string& chr, uint32_t pos,
                    const std::string& a1, const std::string& a2) -> std::string {
    std::string k;
    k.reserve(chr.size() + 20 + a1.size() + a2.size());
    k += chr; k += ':';
    k += std::to_string(pos); k += ':';
    k += a1; k += ':'; k += a2;
    return k;
  };

  std::unordered_map<std::string, size_t> refMap;
  refMap.reserve(refAf.records.size());
  for (size_t i = 0; i < refAf.records.size(); ++i) {
    const auto& r = refAf.records[i];
    refMap.emplace(makeKey(r.chrom, r.pos, r.a1, r.a2), i);
  }

  const int nPop = refAf.nPop;
  std::vector<MultiPopMatchedMarker> matched;
  matched.reserve(plinkData.markerInfo().size());

  for (const auto& mi : plinkData.markerInfo()) {
    auto key = makeKey(mi.chrom, mi.pos, mi.ref, mi.alt);
    auto it = refMap.find(key);
    if (it == refMap.end()) continue;  // no exact match → drop

    const auto& ref = refAf.records[it->second];
    MultiPopMatchedMarker m;
    m.genoIndex = mi.genoIndex;
    m.popAF.resize(nPop);
    m.popAN.resize(nPop);
    for (int p = 0; p < nPop; ++p) {
      m.popAF[p] = ref.popAF[p];
      m.popAN[p] = ref.popAN[p];
    }
    matched.push_back(std::move(m));
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
    const std::string& refAfFile,
    const std::string& sparseGrmFile,
    const std::string& outputFile,
    double refPrevalence,
    double cutoff,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk) {

  const int nCluster = static_cast<int>(residFiles.size());

  // ---- Load per-cluster resid files ----
  infoMsg("Loading %d cluster resid files...", nCluster);
  std::vector<ResidData> clusterResid(nCluster);
  for (int c = 0; c < nCluster; ++c) {
    infoMsg("  Cluster %d: %s", c + 1, residFiles[c].c_str());
    clusterResid[c] = loadResidFile(residFiles[c]);
    infoMsg("    %zu subjects loaded", clusterResid[c].subjects.size());
  }

  // Build combined subject list + per-cluster index arrays
  std::vector<std::string> combinedSubjects;
  std::vector<std::vector<uint32_t>> clusterIndices(nCluster);
  for (int c = 0; c < nCluster; ++c) {
    clusterIndices[c].reserve(clusterResid[c].subjects.size());
    for (const auto& s : clusterResid[c].subjects) {
      clusterIndices[c].push_back(static_cast<uint32_t>(combinedSubjects.size()));
      combinedSubjects.push_back(s);
    }
  }
  infoMsg("  Total subjects: %zu", combinedSubjects.size());

  // ---- Load PLINK data with combined subjects ----
  infoMsg("Loading PLINK data: %s", bfilePrefix.c_str());
  PlinkData plinkData(
      bfilePrefix + ".bed",
      bfilePrefix + ".bim",
      bfilePrefix + ".fam",
      combinedSubjects,
      "ref-first",
      {}, {}, {}, {},
      nSnpPerChunk);
  infoMsg("  %u subjects matched, %u markers",
          plinkData.nSubjUsed(), plinkData.nMarkers());

  // ---- Load multi-population ref-af file ----
  infoMsg("Loading multi-population ref-af file: %s", refAfFile.c_str());
  auto refAf = loadMultiPopRefAfFile(refAfFile);
  infoMsg("  %zu records, %d populations", refAf.records.size(), refAf.nPop);

  // ---- Match markers ----
  infoMsg("Matching markers...");
  auto matchedMulti = matchMultiPopMarkers(plinkData, refAf);
  infoMsg("  %zu markers matched", matchedMulti.size());

  const int nPop = refAf.nPop;
  const size_t nMatched = matchedMulti.size();

  // ---- Per-cluster summix + AF synthesis ----
  infoMsg("Per-cluster ancestry estimation and AF synthesis...");

  // Scan genotypes once for all matched markers → full GVec
  // Then compute per-cluster internal AF
  PlinkCursor cursor(plinkData.bedFile(),
                     static_cast<uint32_t>(plinkData.nMarkers()),
                     plinkData.nSubjInFile(),
                     plinkData.samplePosMap(),
                     plinkData.isAltFirst(),
                     plinkData.isIdentityMap());

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
      const auto& ind = clusterResid[c].indicator;
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

      // AN_ref = 1 / Σ(pi_k² / AN_k)
      double denom = 0.0;
      for (int p = 0; p < nPop; ++p) {
        double an = matchedMulti[m].popAN[p];
        if (an > 0 && proportions[p] > 0)
          denom += (proportions[p] * proportions[p]) / an;
      }
      mi.AN_ref = (denom > 0) ? 1.0 / denom : 0.0;

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
    if (!sparseGrmFile.empty()) {
      grm = std::make_unique<SparseGRM>(sparseGrmFile, clusterResid[c].subjects);
      infoMsg("    Sparse GRM: %u subjects, %zu non-zeros",
              grm->nSubjects(), grm->nnz());
    }

    clRefMaps[c] = testBatchEffects(
        clMatchedInfo[c], clusterResid[c], grm.get(), refPrevalence, cutoff);
    infoMsg("    %zu markers retained", clRefMaps[c]->size());
  }

  // ---- Build LEAFMethod and run marker engine ----
  infoMsg("Building LEAF method (%d clusters)...", nCluster);

  std::vector<std::unique_ptr<WtCoxGMethod>> clMethods;
  clMethods.reserve(nCluster);
  for (int c = 0; c < nCluster; ++c) {
    clMethods.push_back(std::make_unique<WtCoxGMethod>(
        clusterResid[c].residuals,
        clusterResid[c].weights,
        cutoff,
        spaCutoff,
        clRefMaps[c]));
  }

  LEAFMethod method(std::move(clMethods), clusterIndices);

  infoMsg("Running LEAF marker engine (%d thread(s))...", nthread);
  markerEngine(plinkData, method, outputFile,
               nthread,
               /*missingCutoff=*/0.15,
               /*minMafMarker=*/0.0,
               /*minMacMarker=*/0.5,
               /*exactHwe=*/false);
}
