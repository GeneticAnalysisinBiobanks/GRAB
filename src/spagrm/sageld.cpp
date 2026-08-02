// sageld.cpp — SAGELD: G×E interaction marker test (pure C++17 / Eigen)
//
// Per marker, produces:
//   P_G, LOG10P_G           — main genetic effect (normal approximation)
//   P_Gx<E>, LOG10P_Gx<E>   — G×E interaction (SPA on combined residual per env)
//   Z_G, Z_Gx<E>            — corresponding z-scores
//   SPA_STATUS_G, SPA_STATUS_Gx<E>
//                           — outcome of the test that produced each p-value,
//                             as static_cast<uint8_t>(spa::Status) (D4)
//
// Combined residual per env:  Resid_combined = R_Gx<E> - λ * R_G
// λ estimated genome-wide from residual covariance.
//
// Two input modes (caller-selected; see sageld.hpp for the CLI surface):
//   Residual mode — pre-computed per-IID R_G, R_<E>, R_Gx<E> from a file.
//   Pheno mode    — long-format Y, X, E; fit  Y ~ X + (E | IID)  internally
//                   via EM-ML; aggregate BLUP residuals to per-IID columns.

#include "spagrm/sageld.hpp"
#include "spagrm/sageld_fit.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "spagrm/grm_null.hpp"
#include "spagrm/spagrm.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/text_scanner.hpp"
#include "util/text_stream.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

namespace {

// ══════════════════════════════════════════════════════════════════════
// EnvSpec — residual-file column layout (residual mode only)
// ══════════════════════════════════════════════════════════════════════

struct EnvSpec {
    std::string name;
    int colE;
    int colGxE;
};

std::vector<EnvSpec> parseEnvSpecs(
    int nRC,
    const std::vector<std::string> &colNames
) {
    if (nRC < 3 || (nRC - 1) % 2 != 0)
        throw std::runtime_error("SAGELD residual file: expected odd number of columns ≥ 3 "
                                 "(R_G  R_<E1>  R_Gx<E1>  [R_<E2>  R_Gx<E2>  ...]), got " +
                                 std::to_string(nRC));

    const int nEnv = (nRC - 1) / 2;
    std::vector<EnvSpec> envs(nEnv);

    if (colNames.empty()) {
        infoMsg("Residual file has no header; parsing by column order: "
                "R_G  R_<E1>  R_Gx<E1>  ...");
        for (int i = 0; i < nEnv; ++i)
            envs[i] = {"E" + std::to_string(i + 1), 1 + i * 2, 2 + i * 2};
        return envs;
    }

    if (static_cast<int>(colNames.size()) != nRC)
        throw std::runtime_error("SAGELD: column name count (" + std::to_string(colNames.size()) +
                                 ") does not match data column count (" + std::to_string(nRC) + ")");

    if (colNames[0] != "R_G")
        throw std::runtime_error("SAGELD: first residual column must be 'R_G', got '" + colNames[0] + "'");

    for (int i = 0; i < nEnv; ++i) {
        const int ce = 1 + i * 2;
        const int cg = 2 + i * 2;
        const std::string &eName = colNames[ce];
        const std::string &gxName = colNames[cg];
        const bool eOk = (eName.size() >= 3 && eName.substr(0, 2) == "R_" &&
                          !(eName.size() >= 4 && eName.substr(0, 4) == "R_Gx"));
        if (!eOk)
            throw std::runtime_error("SAGELD: column " + std::to_string(ce) + " ('" + eName +
                                     "') must match R_<E> (environment residual, e.g. R_AGE)");
        if (gxName.size() < 5 || gxName.substr(0, 4) != "R_Gx")
            throw std::runtime_error("SAGELD: column " + std::to_string(cg) + " ('" + gxName +
                                     "') must match R_Gx<E> (G×E residual, e.g. R_GxAGE)");
        const std::string envFromE = eName.substr(2);
        const std::string envFromGxE = gxName.substr(4);
        if (envFromE != envFromGxE)
            infoMsg("  Warning: env name mismatch '%s' vs '%s' in column pair %d/%d", eName.c_str(), gxName.c_str(), ce,
                    cg);
        envs[i] = {envFromE, ce, cg};
    }
    return envs;
}

// ══════════════════════════════════════════════════════════════════════
// GRMTopology — GRM + IBD topology, loaded once and shared across envs
// ══════════════════════════════════════════════════════════════════════

struct GRMTopology {
    std::vector<SparseGRM::Entry> allEntries;
    std::vector<double> grmDiag;
    std::unordered_set<uint32_t> singletonSet;
    std::vector<std::vector<uint32_t> > families;
    std::vector<std::vector<SparseGRM::Entry> > familyEntries;
    std::vector<nsGRMNull::IndexedIBD> ibdEntries;
    std::unordered_map<uint64_t, uint32_t> ibdPairMap;
};

GRMTopology loadGRMTopology(
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const std::vector<std::string> &subjIDs,
    const std::vector<std::string> &famIIDs,
    uint32_t N
) {
    GRMTopology topo;

    SparseGRM grm = SparseGRM::load(spgrmGrabFile, spgrmGctaFile, subjIDs, famIIDs);
    infoMsg("Sparse GRM: %zu entries (diagonal + off-diag)", grm.nnz());

    auto subjIdMap = text::buildIIDMap(subjIDs);
    infoMsg("Loading pairwise IBD: %s", pairwiseIBDFile.c_str());
    topo.ibdEntries = nsGRMNull::loadIndexedIBD(pairwiseIBDFile, subjIdMap);
    topo.ibdPairMap = nsGRMNull::buildIBDPairMap(topo.ibdEntries);
    infoMsg("Loaded %zu IBD records", topo.ibdEntries.size());

    topo.allEntries = grm.entries();
    topo.grmDiag = grm.diagonal();

    std::vector<std::pair<uint32_t, uint32_t> > edges;
    {
        std::unordered_set<uint64_t> seen;
        for (const auto &e : topo.allEntries) {
            if (e.row == e.col) continue;
            uint32_t lo = std::min(e.row, e.col);
            uint32_t hi = std::max(e.row, e.col);
            uint64_t key = (static_cast<uint64_t>(lo) << 32) | hi;
            if (seen.insert(key).second) edges.push_back({lo, hi});
        }
    }
    auto components = nsGRMNull::getComponents(N, edges);
    infoMsg("Found %zu connected components", components.size());

    std::vector<std::vector<uint32_t> > singletons;
    for (auto &comp : components) {
        if (comp.size() == 1)
            singletons.push_back(std::move(comp));
        else
            topo.families.push_back(std::move(comp));
    }
    infoMsg("Singletons: %zu, Families: %zu", singletons.size(), topo.families.size());
    for (const auto &s : singletons)
        topo.singletonSet.insert(s[0]);

    std::unordered_map<uint32_t, size_t> subjToFamily;
    for (size_t fi = 0; fi < topo.families.size(); ++fi)
        for (uint32_t idx : topo.families[fi])
            subjToFamily[idx] = fi;

    topo.familyEntries.resize(topo.families.size());
    for (const auto &e : topo.allEntries) {
        auto it = subjToFamily.find(e.row);
        if (it != subjToFamily.end()) topo.familyEntries[it->second].push_back(e);
    }
    return topo;
}

double computeRGRMR(
    const Eigen::VectorXd &R,
    const GRMTopology &topo
) {
    double acc = 0.0;
    for (const auto &e : topo.allEntries) {
        double factor = (e.row == e.col) ? 1.0 : 2.0;
        acc += factor * e.value * R[e.row] * R[e.col];
    }
    return acc;
}

// ══════════════════════════════════════════════════════════════════════
// SAGELDSplitLedgers — the _G / _GxE / _G_GxE decomposition of every
//   non-partition-fixed SPAGRM ledger, built ONCE over the fixed partition
//   (the mean-λ combined residual's partition).  A marker-specific λ_i then
//   recombines them:
//     linear   X_i = X_GxE − λ·X_G
//     bilinear X_i = X_GxE + λ²·X_G − λ·X_G_GxE       (X_G_GxE folds the 2)
//   so the exact getMarkerPvalFromScore SPA can run per marker without
//   re-detecting outliers or rebuilding families.  See
//   sageldMarkerLambdaNegLog10P.
// ══════════════════════════════════════════════════════════════════════
struct SAGELDSplitLedgers {
    // Full residual vectors — for the centered Score_i and the ctor's resid.
    Eigen::VectorXd resid_G, resid_GxE;

    // Unrelated-outlier residual sub-vectors (linear).
    Eigen::VectorXd ruo_G, ruo_GxE;

    // sum_R_nonOutlier (linear).
    double sumRnon_G = 0.0, sumRnon_GxE = 0.0;
    // R_GRM_R total (bilinear) — the Score-variance denominator.
    double RGRMR_G = 0.0, RGRMR_GxE = 0.0, RGRMR_G_GxE = 0.0;
    // R_GRM_R_nonOutlier (bilinear).
    double RGRMRnon_G = 0.0, RGRMRnon_GxE = 0.0, RGRMRnon_G_GxE = 0.0;
    // R_GRM_R_TwoSubjOutlier (bilinear).
    double RGRMRtwo_G = 0.0, RGRMRtwo_GxE = 0.0, RGRMRtwo_G_GxE = 0.0;

    // Two-subject residual pairs (linear), aligned with TwoSubj_rho_list.
    std::vector<std::array<double, 2> > two_G, two_GxE;
    // Three-subject standS vectors (linear), aligned with ThreeSubj_CLT_list.
    std::vector<std::vector<double> > three_G, three_GxE;

    // Partition-fixed data copied verbatim.
    std::vector<double> MAF_interval;
    std::vector<std::vector<double> > TwoSubj_rho_list;
    std::vector<Eigen::MatrixXd> ThreeSubj_CLT_list;
    double SPA_Cutoff = 0.0, tol = 0.0;
};

// Build the split ledgers over a fixed partition (from buildSPAGRMNullModel's
// outPartition on the mean-λ combined residual).
SAGELDSplitLedgers buildSAGELDSplitLedgers(
    const nsGRMNull::SPAGRMPartition &part,
    const Eigen::VectorXd &Resid_G,
    const Eigen::VectorXd &Resid_GxE,
    const GRMTopology &topo,
    double spaCutoff
) {
    SAGELDSplitLedgers L;
    L.resid_G = Resid_G;
    L.resid_GxE = Resid_GxE;

    // Bilinear quad accumulator over an entry list into (G, GxE, cross).
    auto addQuad = [&](const SparseGRM::Entry &e, double &qG, double &qGxE, double &qX) {
        const double a = ((e.row == e.col) ? 1.0 : 2.0) * e.value;
        const double rgR = Resid_G[e.row], rgC = Resid_G[e.col];
        const double reR = Resid_GxE[e.row], reC = Resid_GxE[e.col];
        qG   += a * rgR * rgC;
        qGxE += a * reR * reC;
        qX   += a * (reR * rgC + rgR * reC);
    };

    // ── R_GRM_R total over the full block-diagonal GRM ──
    for (const auto &e : topo.allEntries)
        addQuad(e, L.RGRMR_G, L.RGRMR_GxE, L.RGRMR_G_GxE);

    // ── Non-outlier singletons: diagonal quad + linear sum ──
    for (uint32_t idx : part.nonOutlierSingletons) {
        const double d = topo.grmDiag[idx];
        const double rg = Resid_G[idx], re = Resid_GxE[idx];
        L.RGRMRnon_G     += d * rg * rg;
        L.RGRMRnon_GxE   += d * re * re;
        L.RGRMRnon_G_GxE += d * 2.0 * re * rg;
        L.sumRnon_G   += rg;
        L.sumRnon_GxE += re;
    }
    // ── Non-outlier family/sub-family members: linear sum ──
    for (uint32_t idx : part.nonOutlierFamilyMembers) {
        L.sumRnon_G   += Resid_G[idx];
        L.sumRnon_GxE += Resid_GxE[idx];
    }
    // ── Non-outlier kept entries: quad ──
    for (const auto &e : part.nonOutlierFamilyEntries)
        addQuad(e, L.RGRMRnon_G, L.RGRMRnon_GxE, L.RGRMRnon_G_GxE);

    // ── Unrelated outliers (ordered) ──
    const int nuo = static_cast<int>(part.unrelatedOutliers.size());
    L.ruo_G.resize(nuo);
    L.ruo_GxE.resize(nuo);
    for (int k = 0; k < nuo; ++k) {
        const uint32_t idx = part.unrelatedOutliers[k];
        L.ruo_G[k]   = Resid_G[idx];
        L.ruo_GxE[k] = Resid_GxE[idx];
    }

    // ── Two-subject outlier families ──
    L.two_G.reserve(part.twoSubjFamilies.size());
    L.two_GxE.reserve(part.twoSubjFamilies.size());
    L.TwoSubj_rho_list.reserve(part.twoSubjFamilies.size());
    for (const auto &g : part.twoSubjFamilies) {
        const double r1g = Resid_G[g.a], r2g = Resid_G[g.b];
        const double r1e = Resid_GxE[g.a], r2e = Resid_GxE[g.b];
        L.RGRMRtwo_G   += g.diagA * r1g * r1g + g.diagB * r2g * r2g + 2.0 * g.offDiag * r1g * r2g;
        L.RGRMRtwo_GxE += g.diagA * r1e * r1e + g.diagB * r2e * r2e + 2.0 * g.offDiag * r1e * r2e;
        L.RGRMRtwo_G_GxE += g.diagA * 2.0 * r1e * r1g + g.diagB * 2.0 * r2e * r2g +
                            2.0 * g.offDiag * (r1e * r2g + r1g * r2e);
        L.two_G.push_back({r1g, r2g});
        L.two_GxE.push_back({r1e, r2e});
        L.TwoSubj_rho_list.push_back({g.rho[0], g.rho[1]});
    }

    // ── Three-or-more-subject outlier families ──
    L.three_G.reserve(part.threeSubjFamilies.size());
    L.three_GxE.reserve(part.threeSubjFamilies.size());
    L.ThreeSubj_CLT_list.reserve(part.threeSubjFamilies.size());
    for (const auto &g : part.threeSubjFamilies) {
        const int n1 = static_cast<int>(g.members.size());
        std::vector<double> rg(n1), re(n1);
        for (int i = 0; i < n1; ++i) {
            rg[i] = Resid_G[g.members[i]];
            re[i] = Resid_GxE[g.members[i]];
        }
        L.three_G.push_back(nsGRMNull::buildStandS(n1, rg));
        L.three_GxE.push_back(nsGRMNull::buildStandS(n1, re));
        L.ThreeSubj_CLT_list.push_back(g.CLT);
    }

    L.MAF_interval = part.mafInterval;
    L.SPA_Cutoff = spaCutoff;
    L.tol = nsGRMNull::TOL_DEFAULT;
    return L;
}

// Marker-specific λ p-value: recombine the split ledgers at λ_i, build a
// per-marker SPAGRMClass, and run the validated getMarkerPvalFromScore.  The
// score keeps GRAB2's genotype-mean centering (SG-1 unchanged); outScore
// returns the centered score so the caller forms BETA = Score / Var(S).
//
// Returns the tier's complete result — the magnitude L = −log10(P) and the
// saddlepoint Status — because both are output columns as of log10p_unify
// Stage 7.  It returned the linear p until Stage 4, then the bare magnitude
// until Stage 7; the P column SAGELD still emits is `spa::pFromNegLog10P` of
// the magnitude, and goes in Stage 8.
spa::Result sageldMarkerLambdaResult(
    const SAGELDSplitLedgers &L,
    double lambda_i,
    const Eigen::Ref<const Eigen::VectorXd> &GVec,
    double altFreq,
    double &zScore,
    double &outScore,
    double *outScoreVar
) {
    const double lam = lambda_i;
    const double lam2 = lam * lam;

    const double R_GRM_R_i         = L.RGRMR_GxE    + lam2 * L.RGRMR_G    - lam * L.RGRMR_G_GxE;
    const double R_GRM_R_nonOut_i  = L.RGRMRnon_GxE + lam2 * L.RGRMRnon_G - lam * L.RGRMRnon_G_GxE;
    const double R_GRM_R_TwoSubj_i = L.RGRMRtwo_GxE + lam2 * L.RGRMRtwo_G - lam * L.RGRMRtwo_G_GxE;
    const double sum_R_nonOut_i    = L.sumRnon_GxE  - lam * L.sumRnon_G;

    nsSPAGRM::FamilyData fd;
    fd.resid_unrelated_outliers = L.ruo_GxE - lam * L.ruo_G;
    const size_t n2 = L.two_G.size();
    fd.twoSubj_resid.resize(n2);
    for (size_t k = 0; k < n2; ++k) {
        fd.twoSubj_resid[k][0] = L.two_GxE[k][0] - lam * L.two_G[k][0];
        fd.twoSubj_resid[k][1] = L.two_GxE[k][1] - lam * L.two_G[k][1];
    }
    fd.twoSubj_rho = L.TwoSubj_rho_list;
    const size_t n3 = L.three_G.size();
    fd.threeSubj_standS.resize(n3);
    for (size_t j = 0; j < n3; ++j) {
        const size_t sz = L.three_G[j].size();
        fd.threeSubj_standS[j].resize(sz);
        for (size_t k = 0; k < sz; ++k)
            fd.threeSubj_standS[j][k] = L.three_GxE[j][k] - lam * L.three_G[j][k];
    }
    fd.threeSubj_CLT = L.ThreeSubj_CLT_list;

    Eigen::VectorXd resid_i = L.resid_GxE - lam * L.resid_G;
    const double gMean = GVec.mean();
    const double resid_i_sum = L.resid_GxE.sum() - lam * L.resid_G.sum();
    const double Score_i = GVec.dot(resid_i) - gMean * resid_i_sum;
    outScore = Score_i;

    SPAGRMClass sg(std::move(resid_i), sum_R_nonOut_i, R_GRM_R_nonOut_i,
                   R_GRM_R_TwoSubj_i, R_GRM_R_i, L.MAF_interval, std::move(fd),
                   L.SPA_Cutoff, L.tol);
    // getMarkerPvalFromScore has returned {−log10(P), Status} since spa_unify
    // Stage 4; log10p_unify Stage 7 stops discarding the second field at this
    // call site, which is the first of the six known gaps CLAUDE.md lists.
    return sg.getMarkerPvalFromScore(Score_i, altFreq, zScore, outScoreVar);
}

// ══════════════════════════════════════════════════════════════════════
// SAGELDMethod — MethodBase adapter
// ══════════════════════════════════════════════════════════════════════

class SAGELDMethod : public MethodBase {
  public:
    struct PerEnv {
        explicit PerEnv(SPAGRMClass sg) : spagrm_combined(std::move(sg)) {}

        SPAGRMClass spagrm_combined;

        // ── SG-3 marker-specific λ branch (empty/null when disabled) ──
        Eigen::VectorXd Resid_E;                                 // zScoreE numerator
        double R_GRM_R_E = 0.0;                                  // zScoreE denominator
        double cutoff = 0.0;                                     // qnorm(PvalueCutoff/2, upper)
        std::shared_ptr<const nsSAGELDFit::GallopCache> cache;   // markerLambda(GVec, *cache)
        std::shared_ptr<const SAGELDSplitLedgers> splitLedgers;  // recombination ledgers

        bool markerLambdaEnabled() const {
            return cache != nullptr && splitLedgers != nullptr && R_GRM_R_E > 0.0;
        }
    };

    SAGELDMethod(
        Eigen::VectorXd resid_G,
        double R_GRM_R_G,
        std::vector<PerEnv> envs,
        std::vector<std::string> envNames
    )
        : m_resid_G(std::move(resid_G)),
          m_R_GRM_R_G(R_GRM_R_G),
          m_envs(std::move(envs)),
          m_envNames(std::move(envNames))
    {
    }

    std::unique_ptr<MethodBase> clone() const override {
        return std::make_unique<SAGELDMethod>(*this);
    }

    int resultSize() const override {
        return 6 + 7 * static_cast<int>(m_envs.size());
    }

    // Schema: a 6-wide G main-effect block (normal-approx only, so its Z is
    // already p-consistent and needs no Z_Norm) followed by a 7-wide block per
    // env G×E test (SPA score test, so the p-consistent Z and the raw-score
    // Z_Norm diverge in the tails):
    //   col 0..5      :  P_G  LOG10P_G  Z_G  BETA_G  SE_G  SPA_STATUS_G
    //   col 6 + 7·e+0 :  P_Gx<Ee>
    //   col 6 + 7·e+1 :  LOG10P_Gx<Ee>
    //   col 6 + 7·e+2 :  Z_Gx<Ee>        (p-consistent: 2·pnorm(−|Z|)=P)
    //   col 6 + 7·e+3 :  Z_Norm_Gx<Ee>   (raw Score/√Var)
    //   col 6 + 7·e+4 :  BETA_Gx<Ee>
    //   col 6 + 7·e+5 :  SE_Gx<Ee>
    //   col 6 + 7·e+6 :  SPA_STATUS_Gx<Ee>
    //
    // LOG10P and SPA_STATUS are new in log10p_unify Stage 7.  They were not
    // absent for a reason: SPAGRMClass::getMarkerPvalFromScore — the very
    // routine that produces the G×E p-value, and the same one SPAGRM itself
    // calls — has always returned all three, and this file kept only the
    // magnitude.  CLAUDE.md listed the discard as the first of its known gaps
    // in the output-column contract, and a tier change could move SAGELD's
    // p-values with no status column present to diagnose it with.
    //
    // The G main effect is a plain two-sided normal test with no saddlepoint,
    // so its status is the constant `1 NORMAL` wherever Var(S_G) > 0 and
    // `8 NA_NO_TEST` where it is not (decision D4).  Its LOG10P is evaluated in
    // the log domain by spa::normalBranch rather than as −log10 of the linear
    // tail, which is the same rule every other block in the tree now follows:
    // 2·Φ(−|z|) is exactly zero beyond |z| ≈ 38.6, and −log10 of that is +∞.
    std::string getHeaderColumns() const override {
        std::string h = "\tP_G\tLOG10P_G\tZ_G\tBETA_G\tSE_G\tSPA_STATUS_G";
        for (const auto &n : m_envNames)
            h += "\tP_Gx" + n + "\tLOG10P_Gx" + n + "\tZ_Gx" + n +
                 "\tZ_Norm_Gx" + n + "\tBETA_Gx" + n + "\tSE_Gx" + n +
                 "\tSPA_STATUS_Gx" + n;
        return h;
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int,
        std::vector<double> &result
    ) override {
        const int nE = static_cast<int>(m_envs.size());
        result.assign(6 + 7 * nE, std::numeric_limits<double>::quiet_NaN());
        // Every block gets a status even when it produces nothing, so that a
        // row of NA states which of the two reasons it is (D4).
        const double stNoTest =
            static_cast<double>(static_cast<uint8_t>(spa::Status::NaNoTest));
        result[5] = stNoTest;
        for (int i = 0; i < nE; ++i) result[6 + 7 * i + 6] = stNoTest;

        const double MAF = std::min(altFreq, 1.0 - altFreq);
        const double G_var = 2.0 * MAF * (1.0 - MAF);
        const double gMean = GVec.mean();

        // ── G main effect (normal approximation; SG-1 centering unchanged) ──
        const double Score_G = GVec.dot(m_resid_G) - gMean * m_resid_G.sum();
        const double Var_G = G_var * m_R_GRM_R_G;
        if (Var_G > 0.0 && MAF > 0.0) {
            const double sdG = std::sqrt(Var_G);
            const double zG = Score_G / sdG;
            // The magnitude first, the linear p from it (Stage 7).  The tail
            // is the same 2·Φ(−|z|) that `2.0 * math::pnorm(|z|, upper)` gave,
            // evaluated in the log domain so that it keeps a magnitude past
            // |z| ≈ 38.6 where the linear form flushes to exactly zero.
            const spa::Result rG = spa::normalBranch(zG);
            writePLZBetaSe(result, 0, rG, zG, Score_G / Var_G, 1.0 / sdG);
        }

        for (int i = 0; i < nE; ++i) {
            auto &env = m_envs[i];
            double zGxE = 0.0;
            double scoreVar = 0.0;
            double score = 0.0;
            // The magnitude L = −log10(P) is what both branches below produce
            // and what the Z inversion consumes; the P column is derived from
            // it (Stage 4).  The Status travels with it since Stage 7.
            spa::Result rGxE{std::numeric_limits<double>::quiet_NaN(),
                             spa::Status::NaNoTest};

            // ── SG-3 screen: per-marker environment z-score ──
            // zScoreE = (Resid_E · G) / sqrt(2·MAF·(1−MAF)·R_GRM_R_E), matching
            // estimateLambdaPerMarker; |zScoreE| ≥ cutoff triggers a marker-
            // specific λ_i (only when the GALLOP cache + split ledgers exist).
            bool useMarkerLambda = false;
            if (env.markerLambdaEnabled()) {
                const double denomE = std::sqrt(G_var * env.R_GRM_R_E);
                if (denomE > 0.0) {
                    const double zScoreE = GVec.dot(env.Resid_E) / denomE;
                    useMarkerLambda = (std::abs(zScoreE) >= env.cutoff);
                }
            }

            if (useMarkerLambda) {
                const double lambda_i = nsSAGELDFit::markerLambda(GVec, *env.cache);
                if (std::isfinite(lambda_i)) {
                    rGxE = sageldMarkerLambdaResult(*env.splitLedgers, lambda_i, GVec,
                                                    altFreq, zGxE, score, &scoreVar);
                } else {
                    useMarkerLambda = false; // degenerate V → mean-λ fallback
                }
            }
            if (!useMarkerLambda) {
                // Mean-λ path (SG-1 centering unchanged).
                score = GVec.dot(env.spagrm_combined.resid()) -
                        gMean * env.spagrm_combined.residSum();
                rGxE = env.spagrm_combined.getMarkerPvalFromScore(
                    score, altFreq, zGxE, &scoreVar);
            }

            const int base = 6 + 7 * i;
            const double pGxE = spa::pFromNegLog10P(rGxE.negLog10p);
            result[base + 6] =
                static_cast<double>(static_cast<uint8_t>(rGxE.status));
            if (scoreVar > 0.0) {
                const double sdE = std::sqrt(scoreVar);
                // Z from L (Stage 4): the linear inversion saturated at
                // |Z| = 37.0470962993612 for every L >= 299.698970.
                writePLZZnormBetaSe(result, base, pGxE, rGxE.negLog10p,
                                    math::zFromNegLog10P(rGxE.negLog10p, zGxE), zGxE,
                                    score / scoreVar, 1.0 / sdE);
            } else {
                // P / LOG10P even when var = 0: the SPA may still return a
                // finite magnitude, and its status says which case this is.
                result[base + 0] = pGxE;
                result[base + 1] = rGxE.negLog10p;
            }
        }
    }

    // SAGELD does NOT participate in the fused union-level GEMM.  The SG-3
    // marker-specific λ branch needs the raw per-marker genotype vector — for
    // markerLambda() and the zScoreE screen — which the fused score-only path
    // does not expose, so every marker routes through getResultVec.  The score
    // dot-product is ≈1% of runtime (the method is saddlepoint-dominated), so
    // dropping the fused GEMM is negligible; supportsFusedGemm() therefore
    // returns false and the fillUnionResiduals / fillResidualSums /
    // processScoreBatch hooks inherit the no-op base defaults.
    int preferredBatchSize() const override {
        return 16;
    }

    bool supportsFusedGemm() const override {
        return false;
    }

  private:
    // Write (P, LOG10P, Z, BETA, SE, SPA_STATUS) into out at offset, for the
    // normal-approximation G main-effect block.  `r` is the whole tier result,
    // so the magnitude and the status can only ever be written together; P is
    // derived from the magnitude until Stage 8 deletes the column.
    static void writePLZBetaSe(
        std::vector<double> &out,
        size_t offset,
        const spa::Result &r, double z, double beta, double se
    ) {
        out[offset + 0] = spa::pFromNegLog10P(r.negLog10p);
        out[offset + 1] = r.negLog10p;
        out[offset + 2] = z;
        out[offset + 3] = beta;
        out[offset + 4] = se;
        out[offset + 5] = static_cast<double>(static_cast<uint8_t>(r.status));
    }

    // Write (P, LOG10P, Z, Z_Norm, BETA, SE) into out at offset, for the SPA
    // G×E blocks.  zP is the p-consistent z (2·pnorm(−|zP|)=P after SPA); zNorm
    // is the raw score z (Score/√Var).  SPA_STATUS is written by the caller,
    // because it is set on every path including the ones that reach neither
    // this helper nor a p-value at all.
    static void writePLZZnormBetaSe(
        std::vector<double> &out,
        size_t offset,
        double p, double negLog10p, double zP, double zNorm, double beta, double se
    ) {
        out[offset + 0] = p;
        out[offset + 1] = negLog10p;
        out[offset + 2] = zP;
        out[offset + 3] = zNorm;
        out[offset + 4] = beta;
        out[offset + 5] = se;
    }

    Eigen::VectorXd m_resid_G;
    double m_R_GRM_R_G;
    std::vector<PerEnv> m_envs;
    std::vector<std::string> m_envNames;
};

// ══════════════════════════════════════════════════════════════════════
// GALLOPMethod — MethodBase adapter for the exact GALLOP Wald test
//   (develop-R SAGELD.NullModel UsedMethod="GALLOP").
//
// Pheno-mode only.  Per marker it solves the 2×2 system from the precomputed
// null-model projection cache (nsSAGELDFit::GallopCache) to obtain exact effect
// estimates and standard errors for the genetic main effect (G) and the G×E
// interaction; no fixed residual, no λ, no SPA, no sparse-GRM correction.
//
// The cache is shared (shared_ptr<const>) so per-thread clone() does not copy
// the nSubj×2nx matrices or the LDLT.  The per-subject genotype the engine
// passes to getResultVec is already mean-imputed (matches develop-R
// imputeMethod="mean"); the intercept column in X absorbs the genotype mean.
//
// Output reuses the SAGELD (P, LOG10P, Z, BETA, SE, SPA_STATUS) sextuple
// schema, with Z = β/SE the Wald z and the two-sided normal tail 2·Φ(−|z|)
// evaluated in the log domain (spa::normalBranch), so LOG10P keeps a magnitude
// where the linear tail underflows.  SPA_STATUS is `1 NORMAL` on both blocks
// wherever an estimate exists: GALLOP runs no saddlepoint at all, which is one
// of the two cases decision D4 assigns that code to.  A marker whose 2×2 solve
// failed, or whose SE is not positive, has no statistic and takes
// `8 NA_NO_TEST`.  Single env per object (--envir-name is one column).
// ══════════════════════════════════════════════════════════════════════
class GALLOPMethod : public MethodBase {
  public:
    GALLOPMethod(
        std::shared_ptr<const nsSAGELDFit::GallopCache> cache,
        std::string envName
    )
        : m_cache(std::move(cache)), m_envName(std::move(envName))
    {
    }

    std::unique_ptr<MethodBase> clone() const override {
        return std::make_unique<GALLOPMethod>(*this);
    }

    int resultSize() const override { return 12; }

    int preferredBatchSize() const override {
        return 1;   // per-marker Wald; non-fused, so do not widen the window.
    }

    // col 0..5  : P_G  LOG10P_G  Z_G  BETA_G  SE_G  SPA_STATUS_G
    // col 6..11 : P_Gx<E>  LOG10P_Gx<E>  Z_Gx<E>  BETA_Gx<E>  SE_Gx<E>
    //             SPA_STATUS_Gx<E>
    std::string getHeaderColumns() const override {
        return "\tP_G\tLOG10P_G\tZ_G\tBETA_G\tSE_G\tSPA_STATUS_G"
               "\tP_Gx" + m_envName + "\tLOG10P_Gx" + m_envName +
               "\tZ_Gx" + m_envName + "\tBETA_Gx" + m_envName +
               "\tSE_Gx" + m_envName + "\tSPA_STATUS_Gx" + m_envName;
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double /*altFreq*/,
        int,
        std::vector<double> &result
    ) override {
        result.assign(12, std::numeric_limits<double>::quiet_NaN());
        const double stNoTest =
            static_cast<double>(static_cast<uint8_t>(spa::Status::NaNoTest));
        result[5] = stNoTest;
        result[11] = stNoTest;
        const auto g = nsSAGELDFit::markerGallop(GVec, *m_cache);
        if (!g.ok) return;   // the 2x2 solve failed: no estimate, status stays 8
        writeWald(result, 0, g.beta(0), g.se(0));   // G main effect
        writeWald(result, 6, g.beta(1), g.se(1));   // G×E interaction
    }

  private:
    // Wald (P, LOG10P, Z, BETA, SE, SPA_STATUS) from an effect estimate and its
    // standard error.  The tail is spa::normalBranch, which returns the
    // magnitude and the `1 NORMAL` status together; se <= 0 leaves the row's
    // pre-set `8 NA_NO_TEST`.
    static void writeWald(std::vector<double> &out, size_t off,
                          double beta, double se) {
        out[off + 3] = beta;
        out[off + 4] = se;
        if (se > 0.0) {
            const double z = beta / se;
            const spa::Result r = spa::normalBranch(z);
            out[off + 0] = spa::pFromNegLog10P(r.negLog10p);
            out[off + 1] = r.negLog10p;
            out[off + 2] = z;
            out[off + 5] =
                static_cast<double>(static_cast<uint8_t>(r.status));
        }
    }

    std::shared_ptr<const nsSAGELDFit::GallopCache> m_cache;
    std::string m_envName;
};

// ══════════════════════════════════════════════════════════════════════
// buildSAGELDArtifacts — given per-IID R_G + per-env R_GxE plus shared
//   topology, build a SAGELDMethod (no engine call).  Multiple environments
//   are built sequentially within one call; parallelism across phenotypes
//   is the caller's responsibility (see runSAGELDPhenoMode).
// ══════════════════════════════════════════════════════════════════════
std::unique_ptr<SAGELDMethod> buildSAGELDArtifacts(
    const Eigen::VectorXd &Resid_G,
    const std::vector<std::string> &envNames,
    const std::vector<Eigen::VectorXd> &Resid_GxE_list,
    const GRMTopology &topo,
    uint32_t N,
    double spaCutoff,
    double minMafCutoff,
    double minMacCutoff,
    int innerThreads,
    const std::string &phenoLabel = {},
    const std::vector<double> &lambdaOverride = {},
    // SG-3 marker-specific λ inputs (marker-λ enabled per env only when a
    // GallopCache is supplied and the matching Resid_E has full length N):
    const std::vector<Eigen::VectorXd> &Resid_E_list = {},
    double markerLambdaCutoff = 0.0,
    std::shared_ptr<const nsSAGELDFit::GallopCache> cache = nullptr
) {
    const int nEnv = static_cast<int>(envNames.size());
    if (static_cast<int>(Resid_GxE_list.size()) != nEnv)
        throw std::runtime_error("buildSAGELDArtifacts: envNames / Resid_GxE_list size mismatch");
    if (!lambdaOverride.empty() && static_cast<int>(lambdaOverride.size()) != nEnv)
        throw std::runtime_error("buildSAGELDArtifacts: lambdaOverride size mismatch");

    const double R_GRM_R_G = computeRGRMR(Resid_G, topo);

    std::vector<SAGELDMethod::PerEnv> envData;
    envData.reserve(nEnv);
    std::vector<std::string> envNamesOut;
    envNamesOut.reserve(nEnv);

    for (int i = 0; i < nEnv; ++i) {
        const auto &Resid_GxE = Resid_GxE_list[i];
        // Closed-form OLS λ is the residual-mode fallback (no genotype access).
        // In pheno mode, lambdaOverride supplies the per-marker variance-ratio
        // mean computed via the develop-R GALLOP cache; see
        // estimateLambdaPerMarker() and R/SAGELD.R:347-369.
        const double lambdaFallback =
            Resid_G.squaredNorm() > 0
                ? Resid_GxE.dot(Resid_G) / Resid_G.squaredNorm()
                : 0.0;
        const double lambda = lambdaOverride.empty() ? lambdaFallback : lambdaOverride[i];
        if (phenoLabel.empty())
            infoMsg("Env '%s': lambda = %.6f", envNames[i].c_str(), lambda);
        else
            infoMsg("[%s] Env '%s': lambda = %.6f",
                    phenoLabel.c_str(), envNames[i].c_str(), lambda);
        Eigen::VectorXd envCombined = Resid_GxE - lambda * Resid_G;

        // Marker-specific λ is available only when a GALLOP projection cache
        // is supplied (pheno mode, or residual mode with a serialized cache)
        // and this env carries a full-length Resid_E for the zScoreE screen.
        const bool wantMarkerLambda =
            (cache != nullptr) &&
            (i < static_cast<int>(Resid_E_list.size())) &&
            (Resid_E_list[i].size() == static_cast<Eigen::Index>(N));

        // The partition (outlier detection + family splits + Chow-Liu trees)
        // is driven by the mean-λ combined residual, so the marker-λ split
        // ledgers share exactly the combined model's partition.
        nsGRMNull::SPAGRMPartition part;
        SPAGRMClass spagrm_combined = nsGRMNull::buildSPAGRMNullModel(
            envCombined, N, topo.singletonSet, topo.grmDiag, topo.families,
            topo.familyEntries, topo.allEntries, topo.ibdEntries, topo.ibdPairMap,
            spaCutoff, minMafCutoff, minMacCutoff,
            nsGRMNull::INIT_OUTLIER_RATIO, nsGRMNull::CONTROL_OUTLIER, innerThreads,
            wantMarkerLambda ? &part : nullptr);

        SAGELDMethod::PerEnv pe{std::move(spagrm_combined)};
        if (wantMarkerLambda) {
            pe.Resid_E = Resid_E_list[i];
            pe.R_GRM_R_E = computeRGRMR(Resid_E_list[i], topo);
            pe.cutoff = markerLambdaCutoff;
            pe.cache = cache;
            pe.splitLedgers = std::make_shared<const SAGELDSplitLedgers>(
                buildSAGELDSplitLedgers(part, Resid_G, Resid_GxE, topo, spaCutoff));
        }
        envData.push_back(std::move(pe));
        envNamesOut.push_back(envNames[i]);
    }

    return std::make_unique<SAGELDMethod>(
        Resid_G, R_GRM_R_G, std::move(envData), std::move(envNamesOut));
}

// ══════════════════════════════════════════════════════════════════════
// estimateLambdaPerMarker — develop-R parity λ estimator.
//
// Replicates R/SAGELD.R:347-369:
//   • Sample up to maxMarkers markers (stride-based deterministic sample).
//   • For each marker:
//        mu       = altFreq
//        zScoreE  = R_E · g  /  sqrt(2 mu (1-mu) R_GRM_R_E)
//        accept iff 0.05 < mu < 0.95  and  |zScoreE| < zScoreECutoff
//   • λ_i = V[1,2] / V[1,1] via the GallopCache.
//   • λ̂ = mean(λ_i) over accepted markers; require ≥ 100 accepted.
//
// Returns NaN if fewer than 100 markers passed screening.  Caller falls
// back to the closed-form OLS λ in that case.
// ══════════════════════════════════════════════════════════════════════
double estimateLambdaPerMarker(
    const GenoMeta &genoData,
    uint32_t N,
    const Eigen::VectorXd &Resid_E,
    double R_GRM_R_E,
    double zScoreECutoff,
    const nsSAGELDFit::GallopCache &cache,
    int maxMarkers,
    const std::string &phenoLabel = {}
) {
    if (!(R_GRM_R_E > 0.0))
        return std::numeric_limits<double>::quiet_NaN();

    const uint64_t total = static_cast<uint64_t>(genoData.nMarkers());
    if (total == 0) return std::numeric_limits<double>::quiet_NaN();

    const uint64_t target = std::min(static_cast<uint64_t>(maxMarkers), total);

    auto cursor = genoData.makeCursor();
    Eigen::VectorXd geno(N);
    std::vector<uint32_t> missingIdx;

    std::vector<double> lambdaObs;
    lambdaObs.reserve(static_cast<size_t>(target));

    cursor->beginSequentialBlock(0);
    for (uint64_t k = 0; k < target; ++k) {
        // Even, genome-spanning systematic sample: markerIdx = floor(k·total/target).
        // Exact integer arithmetic (k < target ≤ 2000, so k·total stays well
        // within uint64) keeps it deterministic and platform-independent.
        // Unlike an integer stride = total/target, this does not collapse to the
        // first `target` markers when target < total < 2·target (e.g. 2000 <
        // total < 4000): floor(k·total/target) yields `target` distinct indices
        // spread across the whole [0, total) range for any total ≥ target, since
        // k ≤ target−1 ⇒ markerIdx = floor(k·total/target) < total.
        const uint64_t markerIdx = (k * total) / target;
        double altFreq = 0.0, altCounts = 0.0, missingRate = 0.0;
        double log10pHwe = 0.0, maf = 0.0, mac = 0.0;
        cursor->getGenotypes(markerIdx, geno, altFreq, altCounts, missingRate,
                             log10pHwe, maf, mac, missingIdx);

        // Mean impute missing entries to 2·AF (matches develop-R imputeMethod="mean")
        const double meanGeno = 2.0 * altFreq;
        for (auto idx : missingIdx) geno[idx] = meanGeno;

        // Screening: common variant (MAF > 0.05) and small marginal G-E z-score
        const double mu = altFreq;
        if (!(mu > 0.05 && mu < 0.95)) continue;
        const double denom = std::sqrt(2.0 * mu * (1.0 - mu) * R_GRM_R_E);
        if (!(denom > 0.0)) continue;
        const double zScoreE = Resid_E.dot(geno) / denom;
        if (!(std::abs(zScoreE) < zScoreECutoff)) continue;

        const double lambda = nsSAGELDFit::markerLambda(geno, cache);
        if (std::isfinite(lambda)) lambdaObs.push_back(lambda);
    }

    if (lambdaObs.size() < 100) {
        const char *prefix = phenoLabel.empty() ? "" : phenoLabel.c_str();
        if (phenoLabel.empty())
            infoMsg("WARNING: only %zu markers passed λ screening (< 100); "
                    "falling back to closed-form OLS λ",
                    lambdaObs.size());
        else
            infoMsg("[%s] WARNING: only %zu markers passed λ screening (< 100); "
                    "falling back to closed-form OLS λ",
                    prefix, lambdaObs.size());
        return std::numeric_limits<double>::quiet_NaN();
    }

    double sum = 0.0;
    for (double v : lambdaObs) sum += v;
    const double lambdaMean = sum / static_cast<double>(lambdaObs.size());
    if (phenoLabel.empty())
        infoMsg("λ estimate from %zu markers (per-marker mean): %.6f",
                lambdaObs.size(), lambdaMean);
    else
        infoMsg("[%s] λ estimate from %zu markers (per-marker mean): %.6f",
                phenoLabel.c_str(), lambdaObs.size(), lambdaMean);
    return lambdaMean;
}

// Scan an IID-format residual file for a `##sageld-lambda=<value>` metadata
// comment.  Returns NaN if absent or malformed.  Stops at the first
// non-comment line (typically the data header `#IID`).
double readLambdaFromResidFile(const std::string &filename) {
    std::string line;
    try {
        TextReader ifs(filename); // auto-detects .gz/.zst by extension
        while (ifs.getline(line)) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            if (line.size() < 2 || line[0] != '#' || line[1] != '#') break;
            constexpr const char *kKey = "##sageld-lambda=";
            const size_t kKeyLen = std::strlen(kKey);
            if (line.size() > kKeyLen && line.compare(0, kKeyLen, kKey) == 0) {
                const char *p = line.c_str() + kKeyLen;
                char *ep = nullptr;
                const double v = std::strtod(p, &ep);
                if (ep != p && std::isfinite(v)) return v;
                return std::numeric_limits<double>::quiet_NaN();
            }
        }
    } catch (const std::exception &) {
        // Unreadable file → no λ metadata; fall through to NaN.
    }
    return std::numeric_limits<double>::quiet_NaN();
}

// ══════════════════════════════════════════════════════════════════════
// loadSAGELDCacheFromResid — reconstruct the GALLOP projection cache that
//   --save-resid serialized into the residual file, so residual-input mode
//   can also run the marker-specific λ branch.  The per-subject cache rows
//   are keyed by IID and re-aligned to `usedIIDs` (the residual-mode subject
//   order, matching the engine's genotype decode order).
//
// Returns a null cache (→ mean-λ only) when the file lacks the cache block,
// or when the residual-mode used set does not exactly match the cache's
// subject set (a --keep/--remove/GRM subset would make the cached Q and the
// per-subject sums inconsistent).  `cutoff` is the zScoreE screen threshold
// recovered from ##sageld-pvalue-cutoff (default 0.001).
// ══════════════════════════════════════════════════════════════════════
struct ResidCacheResult {
    std::shared_ptr<const nsSAGELDFit::GallopCache> cache; // null if unavailable
    double cutoff = 0.0;
    bool sawCacheBlock = false; // true once a ##sageld-cache-* / Si_00.. column block
                                // is found, even if later rejected for a subject-set
                                // mismatch — lets the caller warn only when the block
                                // is genuinely absent (legacy / hand-made residual file)
};

ResidCacheResult loadSAGELDCacheFromResid(
    const std::string &filename,
    const std::vector<std::string> &usedIIDs
) {
    ResidCacheResult res;
    double pcut = 0.001; // PvalueCutoff default; header may override

    std::vector<std::string> lines;
    try {
        TextReader ifs(filename); // auto-detects .gz/.zst by extension
        std::string ln;
        while (ifs.getline(ln)) {
            if (!ln.empty() && ln.back() == '\r') ln.pop_back();
            lines.push_back(ln);
        }
    } catch (const std::exception &) {
        res.cutoff = math::qnorm(pcut / 2.0, 0.0, 1.0, false);
        return res;
    }

    double sig = 0.0;
    int nx = -1, nsubj = -1, colHdrLine = -1;
    std::vector<double> qflat;

    for (size_t li = 0; li < lines.size(); ++li) {
        const std::string &ln = lines[li];
        if (ln.empty()) continue;
        if (ln.size() >= 2 && ln[0] == '#' && ln[1] == '#') {
            auto after = [&](const char *k) -> const char * {
                const size_t kl = std::strlen(k);
                return (ln.size() > kl && ln.compare(0, kl, k) == 0) ? ln.c_str() + kl : nullptr;
            };
            if (const char *p = after("##sageld-cache-sig=")) sig = std::strtod(p, nullptr);
            else if (const char *p = after("##sageld-cache-nx=")) nx = std::atoi(p);
            else if (const char *p = after("##sageld-cache-nsubj=")) nsubj = std::atoi(p);
            else if (const char *p = after("##sageld-pvalue-cutoff=")) pcut = std::strtod(p, nullptr);
            else if (const char *p = after("##sageld-cache-Q=")) {
                const char *s = p;
                char *e = nullptr;
                while (true) {
                    const double v = std::strtod(s, &e);
                    if (e == s) break;
                    qflat.push_back(v);
                    s = e;
                }
            }
            continue;
        }
        if (ln[0] == '#') { colHdrLine = static_cast<int>(li); break; } // #IID column header
        break; // a data line before any header → malformed
    }
    res.cutoff = math::qnorm(pcut / 2.0, 0.0, 1.0, false);

    // No cache block → mean-λ only (not an error; legacy/hand-made file).
    if (colHdrLine < 0 || nx <= 0 || nsubj <= 0 ||
        static_cast<int>(qflat.size()) != nx * nx)
        return res;

    // Locate the cache columns by name in the column header.
    std::vector<std::string> colHeaders;
    {
        text::TokenScanner ts(lines[colHdrLine]);
        while (!ts.atEnd()) {
            auto sv = ts.nextView();
            if (sv.empty()) break;
            colHeaders.emplace_back(sv);
        }
    }
    std::unordered_map<std::string, int> colIdx;
    for (int c = 0; c < static_cast<int>(colHeaders.size()); ++c) colIdx[colHeaders[c]] = c;
    auto findCol = [&](const std::string &name) -> int {
        auto it = colIdx.find(name);
        return it == colIdx.end() ? -1 : it->second;
    };

    const int si00 = findCol("Si_00"), si01 = findCol("Si_01"), si11 = findCol("Si_11");
    const int st00 = findCol("StS_00"), st01 = findCol("StS_01"), st11 = findCol("StS_11");
    if (si00 < 0 || si01 < 0 || si11 < 0 || st00 < 0 || st01 < 0 || st11 < 0)
        return res; // no cache block
    std::vector<int> xtsIdx(2 * nx), atsIdx(2 * nx);
    int maxIdx = std::max({si00, si01, si11, st00, st01, st11});
    for (int j = 0; j < 2 * nx; ++j) {
        xtsIdx[j] = findCol("XTs_" + std::to_string(j));
        atsIdx[j] = findCol("AtS_" + std::to_string(j));
        if (xtsIdx[j] < 0 || atsIdx[j] < 0) return res;
        maxIdx = std::max({maxIdx, xtsIdx[j], atsIdx[j]});
    }

    // A genuine cache column block is present; any failure past this point is a
    // subject-set mismatch (warned individually below), not an absent block.
    res.sawCacheBlock = true;

    // Index data rows by IID.
    std::unordered_map<std::string, std::vector<std::string> > rowByIID;
    rowByIID.reserve(lines.size());
    for (size_t li = static_cast<size_t>(colHdrLine) + 1; li < lines.size(); ++li) {
        const std::string &ln = lines[li];
        if (ln.empty() || ln[0] == '#') continue;
        std::vector<std::string> toks;
        text::TokenScanner ts(ln);
        while (!ts.atEnd()) {
            auto sv = ts.nextView();
            if (sv.empty()) break;
            toks.emplace_back(sv);
        }
        if (!toks.empty()) rowByIID.emplace(toks[0], std::move(toks));
    }

    const int N = static_cast<int>(usedIIDs.size());
    if (N != nsubj) {
        warnMsg("SAGELD residual mode: GALLOP cache nsubj=%d != used subjects=%d "
                "(--keep/--remove/GRM subset); marker-specific λ disabled",
                nsubj, N);
        return res;
    }

    auto cache = std::make_shared<nsSAGELDFit::GallopCache>();
    cache->nx = nx;
    cache->nSubj = N;
    cache->sig = sig;
    cache->Q.resize(nx, nx);
    for (int r = 0; r < nx; ++r)
        for (int c = 0; c < nx; ++c)
            cache->Q(r, c) = qflat[static_cast<size_t>(r) * nx + c];
    cache->Qsolver.compute(cache->Q);
    cache->Si.resize(N);
    cache->StS.resize(N);
    cache->XTs.resize(N, 2 * nx);
    cache->AtS.resize(N, 2 * nx);
    // Tr0 / sig drive only the GALLOP Wald path (markerGallop), not markerLambda,
    // so Tr0 is left empty here.

    for (int i = 0; i < N; ++i) {
        auto it = rowByIID.find(usedIIDs[i]);
        if (it == rowByIID.end() || static_cast<int>(it->second.size()) <= maxIdx) {
            warnMsg("SAGELD residual mode: IID '%s' missing/short GALLOP cache row; "
                    "marker-specific λ disabled",
                    usedIIDs[i].c_str());
            res.cache = nullptr;
            return res;
        }
        const auto &t = it->second;
        auto val = [&](int c) -> double { return std::strtod(t[c].c_str(), nullptr); };
        Eigen::Matrix2d Si;
        Si << val(si00), val(si01), val(si01), val(si11);
        Eigen::Matrix2d StS;
        StS << val(st00), val(st01), val(st01), val(st11);
        cache->Si[i] = Si;
        cache->StS[i] = StS;
        for (int j = 0; j < 2 * nx; ++j) {
            cache->XTs(i, j) = val(xtsIdx[j]);
            cache->AtS(i, j) = val(atsIdx[j]);
        }
    }

    res.cache = std::move(cache);
    return res;
}

// ══════════════════════════════════════════════════════════════════════
// runSAGELDCoreSingle — residual-mode wrapper: build artifacts then call
//                       markerEngine for a single output file.
// ══════════════════════════════════════════════════════════════════════
void runSAGELDCoreSingle(
    const Eigen::VectorXd &Resid_G,
    const std::vector<std::string> &envNames,
    const std::vector<Eigen::VectorXd> &Resid_GxE_list,
    const GRMTopology &topo,
    uint32_t N,
    GenoMeta &genoData,
    const std::string &outputFile,
    double spaCutoff,
    int nthreads,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::vector<double> &lambdaOverride = {},
    const std::vector<Eigen::VectorXd> &Resid_E_list = {},
    double markerLambdaCutoff = 0.0,
    std::shared_ptr<const nsSAGELDFit::GallopCache> cache = nullptr
) {
    auto method = buildSAGELDArtifacts(Resid_G, envNames, Resid_GxE_list,
                                       topo, N, spaCutoff, minMafCutoff, minMacCutoff,
                                       /*innerThreads=*/nthreads, /*phenoLabel=*/{},
                                       lambdaOverride, Resid_E_list, markerLambdaCutoff, cache);
    infoMsg("Running SAGELD marker-level association on '%s'...", outputFile.c_str());
    markerEngine(genoData, *method, outputFile, nthreads,
                 missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);
}

// ══════════════════════════════════════════════════════════════════════
// runSAGELDResidualMode — original residual-input mode
// ══════════════════════════════════════════════════════════════════════
void runSAGELDResidualMode(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &outputFile,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile
) {
    infoMsg("Loading residual file: %s", phenoFile.c_str());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    sd.loadResidOne(phenoFile, residNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGrmSubjects(SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, sd.famIIDs()));
    sd.setGenoLabel(geno.flagLabel());
    sd.setGrmLabel(grmFlagLabel(spgrmGrabFile, spgrmGctaFile));
    sd.finalize();
    const uint32_t N = sd.nUsed();
    infoMsg("Loaded %u subjects (intersected with .fam)", N);

    const auto envs_spec = parseEnvSpecs(sd.residOneCols(), sd.residColNames());
    const int nEnv = static_cast<int>(envs_spec.size());
    infoMsg("SAGELD: %d environment(s):", nEnv);
    for (const auto &e : envs_spec)
        infoMsg("  %s (cols %d / %d)", e.name.c_str(), e.colE, e.colGxE);

    auto topo = loadGRMTopology(spgrmGrabFile, spgrmGctaFile, pairwiseIBDFile,
                                sd.usedIIDs(), sd.famIIDs(), N);
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    const Eigen::MatrixXd &residMat = sd.residMatrix();
    const Eigen::VectorXd Resid_G = residMat.col(0);

    std::vector<std::string> envNames(nEnv);
    std::vector<Eigen::VectorXd> Resid_GxE_list(nEnv);
    std::vector<Eigen::VectorXd> Resid_E_list(nEnv);
    for (int i = 0; i < nEnv; ++i) {
        envNames[i] = envs_spec[i].name;
        Resid_GxE_list[i] = residMat.col(envs_spec[i].colGxE);
        Resid_E_list[i] = residMat.col(envs_spec[i].colE); // R_<E> for the zScoreE screen
    }

    // Inherit per-marker λ from fit-mode if the residual file carries one
    // (written by --save-resid as a `##sageld-lambda=<value>` header).  If
    // the user-supplied residual file lacks this comment, the closed-form
    // OLS λ is used as fallback (legacy behaviour).
    const double lambdaFromFile = readLambdaFromResidFile(phenoFile);
    std::vector<double> lambdaOverride;
    if (std::isfinite(lambdaFromFile)) {
        lambdaOverride.assign(nEnv, lambdaFromFile);
        infoMsg("Residual file supplies λ = %.6f via ##sageld-lambda metadata",
                lambdaFromFile);
    }

    // SG-3: reconstruct the GALLOP cache serialized by --save-resid so the
    // marker-specific λ branch runs in residual mode too.  When absent (or the
    // used subject set does not match the cache), cache is null → mean-λ only.
    // The cache is applied only when nEnv == 1 (it is tied to the fitted env).
    ResidCacheResult rc = loadSAGELDCacheFromResid(phenoFile, sd.usedIIDs());
    std::vector<Eigen::VectorXd> markerLambdaResidE;
    std::shared_ptr<const nsSAGELDFit::GallopCache> markerLambdaCache;
    if (rc.cache && nEnv == 1) {
        markerLambdaResidE = Resid_E_list;
        markerLambdaCache = rc.cache;
        infoMsg("Residual file supplies a GALLOP cache (nx=%d); marker-specific "
                "λ branch enabled",
                rc.cache->nx);
    } else if (rc.cache && nEnv != 1) {
        infoMsg("SAGELD residual mode: %d envs but a single-env GALLOP cache; "
                "marker-specific λ disabled",
                nEnv);
    } else if (!rc.sawCacheBlock) {
        // Legacy / hand-made residual file with no serialized GALLOP cache: the
        // marker-specific λ branch cannot run, so fall back to genome-wide mean λ.
        // (A cache block that WAS present but rejected for a subject-set mismatch
        // is already reported by loadSAGELDCacheFromResid, so it is not re-warned.)
        warnMsg("SAGELD residual mode: residual file carries no GALLOP cache block "
                "(##sageld-cache-*); running genome-wide mean λ only. Regenerate the "
                "residual file with pheno-mode --save-resid to enable the "
                "marker-specific λ branch.");
    }

    runSAGELDCoreSingle(Resid_G, envNames, Resid_GxE_list, topo, N, *genoData,
                        outputFile, spaCutoff, nthreads,
                        missingCutoff, minMafCutoff, minMacCutoff, hweCutoff,
                        lambdaOverride, markerLambdaResidE, rc.cutoff, markerLambdaCache);
}

// ── Helper: load --keep / --remove file into a string set ──────────────
std::unordered_set<std::string> loadIIDsFromFile(const std::string &path) {
    std::unordered_set<std::string> out;
    if (path.empty()) return out;
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open IID file: " + path);
    std::string line;
    while (std::getline(ifs, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        text::TokenScanner ts(line);
        if (ts.atEnd()) continue;
        ts.nextView();                       // FID (or IID if single column)
        std::string iid;
        if (!ts.atEnd()) iid = ts.next();    // IID
        else { text::TokenScanner ts2(line); iid = ts2.next(); }
        if (!iid.empty()) out.insert(iid);
    }
    return out;
}

// ══════════════════════════════════════════════════════════════════════
// runSAGELDPhenoMode — long-format pheno + per-phenotype LMM fit, driven
//                      by multiPhenoEngine for single-pass genotype decode.
//
// Pipeline:
//   1. parseLongPheno  → shared subject set (NaN-drop is global; identity
//                        unionToLocal per task).
//   2. loadGRMTopology + makeGenoData (built once).
//   3. Build PhenoTask vector: per phenotype, fit Y ~ X + (E | IID) and
//      construct SAGELDMethod via buildSAGELDArtifacts.  Work-stealing
//      parallelism across phenotypes with nOuter = min(nthreads, nPheno);
//      remaining threads stay idle inside the per-phenotype build (each
//      buildSAGELDArtifacts call does its env loop sequentially).
//   4. multiPhenoEngine: one genotype decode pass writes per-phenotype
//      outputs <outPrefix>.<phenoName>.SAGELD[.gz|.zst].
//
// The auxiliary E ~ (1 | IID) fit (Resid_E) is intentionally dropped:
// the closed-form lambda inside buildSAGELDArtifacts only consumes R_G
// and R_GxE, so the extra fit was logged but never propagated downstream.
// ══════════════════════════════════════════════════════════════════════
void runSAGELDPhenoMode(
    const std::string &phenoFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<std::string> &envNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    bool saveResid,
    bool gallop,
    const std::string &keepFile,
    const std::string &removeFile
) {
    if (envNames.size() != 1)
        throw std::runtime_error("SAGELD pheno mode: --envir-name currently accepts a single env column");
    for (const auto &en : envNames)
        if (std::find(covarNames.begin(), covarNames.end(), en) == covarNames.end())
            throw std::runtime_error("SAGELD pheno mode: env '" + en +
                                     "' must also appear in --covar-name");

    auto famIIDs = parseGenoIIDs(geno);

    // Build the kept-subject filter: GRM ∩ (keep) − (remove).  GALLOP ignores
    // the sparse GRM entirely (relatedness is captured by the random effects),
    // so the GRM intersection is dropped: kept = (keep) − (remove).
    auto grmIDs = gallop ? std::unordered_set<std::string>{}
                         : SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, famIIDs);
    auto keepSet = loadIIDsFromFile(keepFile);
    auto removeSet = loadIIDsFromFile(removeFile);

    std::unordered_set<std::string> kept;
    kept.reserve(famIIDs.size());
    for (const auto &iid : famIIDs) {
        if (!keepSet.empty() && !keepSet.count(iid)) continue;
        if (!removeSet.empty() && removeSet.count(iid)) continue;
        if (!grmIDs.empty() && !grmIDs.count(iid)) continue;
        kept.insert(iid);
    }
    if (kept.empty())
        throw std::runtime_error(
            gallop ? "GALLOP: no subjects survive --keep / --remove intersection"
                   : "SAGELD pheno mode: no subjects survive --keep / --remove / GRM intersection");

    auto longData = nsSAGELDFit::parseLongPheno(
        phenoFile, phenoNames, covarNames, envNames[0], famIIDs, kept);

    const uint32_t N = static_cast<uint32_t>(longData.uniqueIIDs.size());
    if (N == 0) throw std::runtime_error("SAGELD pheno mode: no subjects in long-format file");

    // Build used-mask in .fam coordinates for genoData
    const uint32_t nFam = static_cast<uint32_t>(famIIDs.size());
    std::vector<uint64_t> usedMask((nFam + 63) / 64, 0ULL);
    for (uint32_t i = 0; i < N; ++i) {
        const uint32_t fi = longData.famSubjIdx[i];
        usedMask[fi / 64] |= (1ULL << (fi % 64));
    }

    // GALLOP needs no sparse-GRM topology (no relatedness quadratic forms).
    GRMTopology topo;
    if (!gallop)
        topo = loadGRMTopology(spgrmGrabFile, spgrmGctaFile, pairwiseIBDFile,
                               longData.uniqueIIDs, famIIDs, N);
    auto genoData = makeGenoData(geno, usedMask, nFam, N, nSnpPerChunk);

    const int nPheno = static_cast<int>(phenoNames.size());
    const int nEnv = static_cast<int>(envNames.size());
    const int nOuter = std::max(1, std::min(nthreads, nPheno));
    const int innerThreads = std::max(1, nthreads / nOuter);
    infoMsg("SAGELD pheno mode: %d phenotype(s) × %d env(s) on %u subjects "
            "(%d outer × %d inner threads)",
            nPheno, nEnv, N, nOuter, innerThreads);

    // All phenotypes share the same subject set (parseLongPheno drops a row
    // if any listed pheno / covar / env is NaN), so unionToLocal is identity.
    std::vector<uint32_t> identityMap(N);
    for (uint32_t i = 0; i < N; ++i) identityMap[i] = i;

    std::vector<PhenoTask> tasks(nPheno);

    // develop-R parity λ estimation parameters (R/SAGELD.R:103,301-305,371-376):
    //   PvalueCutoff = 0.001  →  zScoreECutoff = qnorm(0.0005, lower.tail=FALSE)
    //   marker sample size cap = 2000; require ≥ 100 to pass screening
    constexpr int kLambdaSampleSize = 2000;
    constexpr double kLambdaPvalueCutoff = 0.001;
    const double zScoreECutoff = math::qnorm(
        kLambdaPvalueCutoff / 2.0, 0.0, 1.0, /*lower_tail=*/false);

    auto buildOne = [&](int p) {
        const std::string &phenoName = phenoNames[p];

        Eigen::VectorXd y = longData.Y.col(p);
        auto fitMain = nsSAGELDFit::fitRandomSlopeML(longData, y);
        infoMsg("[%s] Y ~ X + (E | IID) converged in %d iters: "
                "sigma2=%.6g, D=[%.6g, %.6g; %.6g, %.6g]",
                phenoName.c_str(), fitMain.iterations, fitMain.sigma2,
                fitMain.D(0, 0), fitMain.D(0, 1), fitMain.D(1, 0), fitMain.D(1, 1));

        // ── GALLOP (exact Wald) branch ─────────────────────────────────
        // No residual aggregation, no λ, no E-fit, no GRM quadratic forms:
        // the per-marker Wald solve reads only the projection cache.
        if (gallop) {
            auto cache = std::make_shared<const nsSAGELDFit::GallopCache>(
                nsSAGELDFit::buildGallopCache(
                    longData, fitMain.D, fitMain.sigma2, fitMain.residPerRow));
            tasks[p].phenoName = phenoName;
            tasks[p].method = std::make_unique<GALLOPMethod>(std::move(cache), envNames[0]);
            tasks[p].unionToLocal = identityMap;
            tasks[p].nUsed = N;
            return;
        }

        Eigen::VectorXd Resid_G = nsSAGELDFit::aggregatePerIID(longData, fitMain.residPerRow);
        Eigen::VectorXd Resid_GxE = nsSAGELDFit::aggregateWeightedPerIID(
            longData, fitMain.residPerRow, longData.E);

        // E ~ 1 + (1|IID) is needed for the marginal G–E screen used by
        // estimateLambdaPerMarker; the same residuals are also written
        // when --save-resid is set.  Develop-R fits this with lme4::lmer;
        // we use the 1-D Brent path of fitRandomInterceptML (matches lme4
        // to machine precision).
        auto fitE = nsSAGELDFit::fitRandomInterceptML(longData, longData.E);
        Eigen::VectorXd Resid_E = nsSAGELDFit::aggregatePerIID(longData, fitE.residPerRow);

        // Per-marker λ via the GALLOP cache (develop-R parity).  When fewer
        // than 100 markers pass the screen, lambdaMean is NaN and the
        // downstream buildSAGELDArtifacts falls back to the closed-form OLS λ.
        const double R_GRM_R_E = computeRGRMR(Resid_E, topo);
        auto gallopCache = std::make_shared<const nsSAGELDFit::GallopCache>(
            nsSAGELDFit::buildGallopCache(
                longData, fitMain.D, fitMain.sigma2, fitMain.residPerRow));
        const double lambdaMean = estimateLambdaPerMarker(
            *genoData, N, Resid_E, R_GRM_R_E, zScoreECutoff,
            *gallopCache, kLambdaSampleSize, phenoName);

        if (saveResid) {
            // When --compression is set, append the codec extension so residual
            // mode's TextReader auto-detects it; otherwise keep a plain .resid.
            std::string residPath = outPrefix + "." + phenoName + ".SAGELD.resid";
            if (compression == "gz") residPath += ".gz";
            else if (compression == "zst") residPath += ".zst";
            std::ostringstream rf;
            // %.17g preserves a double exactly across the round-trip
            // (DBL_DECIMAL_DIG = 17); matches nullmodel::writeResidualsFile.
            char buf[32];
            auto put = [&](double v) { // header value, no leading tab
                std::snprintf(buf, sizeof(buf), "%.17g", v);
                rf << buf;
            };
            auto writeVal = [&](double v) { // data cell, leading tab, NA-aware
                rf << '\t';
                if (std::isnan(v)) {
                    rf << "NA";
                } else {
                    std::snprintf(buf, sizeof(buf), "%.17g", v);
                    rf << buf;
                }
            };

            // VCF-style metadata comment lines.  parseIIDFile skips `##`
            // lines transparently; SAGELD residual-mode parses the λ value
            // back so that fit-mode → residual-mode produces byte-identical
            // marker outputs on the mean-λ markers even when fit-mode used
            // per-marker λ.
            if (std::isfinite(lambdaMean)) {
                rf << "##sageld-lambda=";
                put(lambdaMean);
                rf << "\n";
            }

            // ── SG-3 GALLOP cache block ──
            // Serialize the per-subject projection cache so residual-mode can
            // ALSO run the marker-specific λ branch (markerLambda consumes only
            // nx / nSubj / Q / Si / StS / XTs / AtS — Tr0 and sig are GALLOP-Wald
            // only and are omitted; sig is kept in the header for documentation).
            const nsSAGELDFit::GallopCache &gc = *gallopCache;
            const int nx = gc.nx;
            rf << "##sageld-cache-sig=";
            put(gc.sig);
            rf << "\n##sageld-cache-nx=" << nx << "\n";
            rf << "##sageld-cache-nsubj=" << gc.nSubj << "\n";
            rf << "##sageld-pvalue-cutoff=";
            put(kLambdaPvalueCutoff);
            rf << "\n##sageld-cache-Q=";
            for (int r = 0; r < nx; ++r)
                for (int c = 0; c < nx; ++c) {
                    if (r != 0 || c != 0) rf << ' ';
                    put(gc.Q(r, c));
                }
            rf << "\n";

            // Column header: residuals, then the per-subject cache columns
            // (Si/StS symmetric → 3 each; XTs/AtS → 2·nx each).
            rf << "#IID\tR_G\tR_" << envNames[0] << "\tR_Gx" << envNames[0];
            rf << "\tSi_00\tSi_01\tSi_11\tStS_00\tStS_01\tStS_11";
            for (int j = 0; j < 2 * nx; ++j) rf << "\tXTs_" << j;
            for (int j = 0; j < 2 * nx; ++j) rf << "\tAtS_" << j;
            rf << "\n";

            for (uint32_t i = 0; i < N; ++i) {
                rf << longData.uniqueIIDs[i];
                writeVal(Resid_G[i]);
                writeVal(Resid_E[i]);
                writeVal(Resid_GxE[i]);
                writeVal(gc.Si[i](0, 0));
                writeVal(gc.Si[i](0, 1));
                writeVal(gc.Si[i](1, 1));
                writeVal(gc.StS[i](0, 0));
                writeVal(gc.StS[i](0, 1));
                writeVal(gc.StS[i](1, 1));
                for (int j = 0; j < 2 * nx; ++j) writeVal(gc.XTs(i, j));
                for (int j = 0; j < 2 * nx; ++j) writeVal(gc.AtS(i, j));
                rf << '\n';
            }
            {
                TextWriter tw(residPath, TextWriter::modeFromString(compression),
                              compressionLevel);
                const std::string body = rf.str();
                tw.write(body);
                tw.close();
            }
            infoMsg("[%s] Saved residuals + GALLOP cache: %s (nx=%d, sigma2_E=%.6g, tau2_E=%.6g)",
                    phenoName.c_str(), residPath.c_str(), nx,
                    fitE.sigma2, fitE.D(0, 0));
        }

        // Single env supported here; future multi-env would loop here.
        std::vector<Eigen::VectorXd> Resid_GxE_list;
        Resid_GxE_list.push_back(std::move(Resid_GxE));
        std::vector<std::string> singleEnvName{envNames[0]};
        std::vector<double> lambdaOverride;
        if (std::isfinite(lambdaMean)) lambdaOverride.push_back(lambdaMean);

        // SG-3: pass the E-residual, screen cutoff, and GALLOP cache so the
        // marker-specific λ branch is available at test time.
        std::vector<Eigen::VectorXd> residEList{Resid_E};

        auto method = buildSAGELDArtifacts(
            Resid_G, singleEnvName, Resid_GxE_list, topo, N,
            spaCutoff, minMafCutoff, minMacCutoff, innerThreads, phenoName,
            lambdaOverride, residEList, zScoreECutoff, gallopCache);

        tasks[p].phenoName = phenoName;
        tasks[p].method = std::move(method);
        tasks[p].unionToLocal = identityMap;
        tasks[p].nUsed = N;
    };

    if (nOuter == 1) {
        for (int p = 0; p < nPheno; ++p) buildOne(p);
    } else {
        std::atomic<int> nextP{0};
        std::exception_ptr workerErr;
        std::mutex errMu;
        auto worker = [&]() {
            try {
                while (true) {
                    int p = nextP.fetch_add(1, std::memory_order_relaxed);
                    if (p >= nPheno) break;
                    buildOne(p);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lk(errMu);
                if (!workerErr) workerErr = std::current_exception();
            }
        };
        std::vector<std::thread> workers;
        workers.reserve(nOuter - 1);
        for (int t = 1; t < nOuter; ++t) workers.emplace_back(worker);
        worker();
        for (auto &w : workers) w.join();
        if (workerErr) std::rethrow_exception(workerErr);
    }

    const char *methodLabel = gallop ? "GALLOP" : "SAGELD";
    infoMsg("Running %s marker tests via multiPhenoEngine (%d phenotype(s), %d threads)...",
            methodLabel, nPheno, nthreads);
    multiPhenoEngine(*genoData, tasks, outPrefix, methodLabel,
                     compression, compressionLevel, nthreads,
                     missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// runSAGELD — public entry point
// ══════════════════════════════════════════════════════════════════════
void runSAGELD(
    const std::string &phenoFile,
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<std::string> &envNames,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &pairwiseIBDFile,
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    bool saveResid,
    bool gallop,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const bool isResidMode = !residNames.empty();
    const bool isPhenoMode = !phenoNames.empty();

    if (isResidMode && isPhenoMode)
        throw std::runtime_error("SAGELD: --resid-name and --pheno-name are mutually exclusive");
    if (!isResidMode && !isPhenoMode)
        throw std::runtime_error("SAGELD: need either --resid-name (residual mode) or --pheno-name + --envir-name (pheno mode)");

    if (isResidMode) {
        if (gallop)
            throw std::runtime_error("SAGELD: --sageld-method gallop requires --pheno-name (GALLOP needs the in-memory null-model projection cache and cannot consume a residual file)");
        if (saveResid)
            throw std::runtime_error("SAGELD: --save-resid requires --pheno-name (residual-input mode has no null model to save)");
        std::string resOut = outPrefix + ".SAGELD";
        if (compression == "gz") resOut += ".gz";
        else if (compression == "zst") resOut += ".zst";

        runSAGELDResidualMode(phenoFile, residNames, spgrmGrabFile, spgrmGctaFile, pairwiseIBDFile,
                              geno, resOut, spaCutoff, nthreads, nSnpPerChunk,
                              missingCutoff, minMafCutoff, minMacCutoff, hweCutoff,
                              keepFile, removeFile);
        return;
    }

    // Pheno mode
    if (envNames.empty())
        throw std::runtime_error("SAGELD pheno mode: --envir-name is required");
    if (covarNames.empty())
        throw std::runtime_error("SAGELD pheno mode: --covar-name is required (must include the --envir-name variable)");
    if (gallop && saveResid)
        throw std::runtime_error("SAGELD: --save-resid is incompatible with --sageld-method gallop (GALLOP produces no residual vector)");
    runSAGELDPhenoMode(phenoFile, phenoNames, covarNames, envNames,
                       spgrmGrabFile, spgrmGctaFile, pairwiseIBDFile, geno,
                       outPrefix, compression, compressionLevel,
                       spaCutoff, nthreads, nSnpPerChunk,
                       missingCutoff, minMafCutoff, minMacCutoff, hweCutoff,
                       saveResid, gallop, keepFile, removeFile);
}
