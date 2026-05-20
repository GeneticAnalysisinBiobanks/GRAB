// sageld.cpp — SAGELD: G×E interaction marker test (pure C++17 / Eigen)
//
// Per marker, produces:
//   P_G       — main genetic effect (normal approximation)
//   P_Gx<E>   — G×E interaction (SPA on combined residual per env)
//   Z_G, Z_Gx<E> — corresponding z-scores
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
#include <fstream>
#include <limits>
#include <memory>
#include <mutex>
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
// SAGELDMethod — MethodBase adapter
// ══════════════════════════════════════════════════════════════════════

class SAGELDMethod : public MethodBase {
  public:
    struct PerEnv {
        SPAGRMClass spagrm_combined;
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
        return 4 + 4 * static_cast<int>(m_envs.size());
    }

    // Paired schema, one (P, Z, BETA, SE) quadruple per test (G main effect
    // followed by each env's G×E test):
    //   col 0..3      :  P_G  Z_G  BETA_G  SE_G
    //   col 4 + 4·e+0 :  P_Gx<Ee>
    //   col 4 + 4·e+1 :  Z_Gx<Ee>
    //   col 4 + 4·e+2 :  BETA_Gx<Ee>
    //   col 4 + 4·e+3 :  SE_Gx<Ee>
    std::string getHeaderColumns() const override {
        std::string h = "\tP_G\tZ_G\tBETA_G\tSE_G";
        for (const auto &n : m_envNames)
            h += "\tP_Gx" + n + "\tZ_Gx" + n + "\tBETA_Gx" + n + "\tSE_Gx" + n;
        return h;
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int,
        std::vector<double> &result
    ) override {
        const int nE = static_cast<int>(m_envs.size());
        result.assign(4 + 4 * nE, std::numeric_limits<double>::quiet_NaN());

        const double MAF = std::min(altFreq, 1.0 - altFreq);
        const double G_var = 2.0 * MAF * (1.0 - MAF);
        const double Score_G = GVec.dot(m_resid_G) - GVec.mean() * m_resid_G.sum();
        const double Var_G = G_var * m_R_GRM_R_G;
        if (Var_G > 0.0 && MAF > 0.0) {
            const double sdG = std::sqrt(Var_G);
            const double zG = Score_G / sdG;
            const double pG = 2.0 * math::pnorm(std::abs(zG), 0.0, 1.0, false);
            writePZBetaSe(result, 0, pG, zG, Score_G / Var_G, 1.0 / sdG);
        }

        for (int i = 0; i < nE; ++i) {
            double zGxE = 0.0;
            double scoreVar = 0.0;
            const double centeredScore =
                GVec.dot(m_envs[i].spagrm_combined.resid()) -
                GVec.mean() * m_envs[i].spagrm_combined.residSum();
            const double pGxE = m_envs[i].spagrm_combined
                                    .getMarkerPvalFromScore(centeredScore, altFreq, zGxE, &scoreVar);
            if (scoreVar > 0.0) {
                const double sdE = std::sqrt(scoreVar);
                writePZBetaSe(result, 4 + 4 * i, pGxE, zGxE, centeredScore / scoreVar, 1.0 / sdE);
            } else {
                result[4 + 4 * i + 0] = pGxE; // P even when var=0 (SPA may still return finite p)
            }
        }
    }

    // ── Fused union-level GEMM interface ───────────────────────────────
    // Column layout for the augmented residual matrix:
    //   col 0           : R_G                    (main genetic effect)
    //   col 1 + e (e<E) : combined residual per env (= R_Gx<E> − λ_e · R_G)
    // For B markers per chunk, multiPhenoEngine produces a (1+E) × B
    // score matrix via one fused GEMM; processScoreBatch then converts the
    // raw scores into per-marker (P_G, P_GxE_e, Z_G, Z_GxE_e) via the same
    // normal-approx (G) and SPA (GxE) reductions as getResultVec.

    int preferredBatchSize() const override {
        return 16;
    }

    bool supportsFusedGemm() const override {
        return true;
    }

    int fusedGemmColumns() const override {
        return 1 + static_cast<int>(m_envs.size());
    }

    void fillUnionResiduals(
        Eigen::Ref<Eigen::MatrixXd> dest,
        const std::vector<uint32_t> &unionToLocal
    ) const override {
        // dest is pre-zeroed, N_union × (1 + nEnv).
        const uint32_t nUnion = static_cast<uint32_t>(unionToLocal.size());
        const int nE = static_cast<int>(m_envs.size());
        for (uint32_t i = 0; i < nUnion; ++i) {
            const uint32_t li = unionToLocal[i];
            if (li == UINT32_MAX) continue;
            dest(i, 0) = m_resid_G[li];
            for (int e = 0; e < nE; ++e)
                dest(i, 1 + e) = m_envs[e].spagrm_combined.resid()[li];
        }
    }

    void fillResidualSums(double *dest) const override {
        dest[0] = m_resid_G.sum();
        const int nE = static_cast<int>(m_envs.size());
        for (int e = 0; e < nE; ++e)
            dest[1 + e] = m_envs[e].spagrm_combined.residSum();
    }

    void processScoreBatch(
        const Eigen::Ref<const Eigen::MatrixXd> &scores,
        const double *gSums,
        const double *gSumSqs,
        uint32_t nUsed,
        const std::vector<double> &altFreqs,
        const std::vector<int> & /*chunkIdxs*/,
        std::vector<std::vector<double> > &results
    ) override {
        (void)gSumSqs;
        const int B = static_cast<int>(scores.cols());
        const int nE = static_cast<int>(m_envs.size());
        results.resize(B);

        const double residG_sum = m_resid_G.sum();
        const double invN = 1.0 / static_cast<double>(nUsed);

        for (int b = 0; b < B; ++b) {
            auto &out = results[b];
            out.assign(4 + 4 * nE, std::numeric_limits<double>::quiet_NaN());

            const double altFreq = altFreqs[b];
            const double MAF = std::min(altFreq, 1.0 - altFreq);
            const double gMean = gSums[b] * invN;

            // P_G via the normal approximation, sharing the SPAGRM-style
            // score centering: Score_G = scores(0,b) - gMean · Σ R_G.
            const double Score_G = scores(0, b) - gMean * residG_sum;
            const double G_var = 2.0 * MAF * (1.0 - MAF);
            const double Var_G = G_var * m_R_GRM_R_G;
            if (Var_G > 0.0 && MAF > 0.0) {
                const double sdG = std::sqrt(Var_G);
                const double zG = Score_G / sdG;
                const double pG = 2.0 * math::pnorm(std::abs(zG), 0.0, 1.0, false);
                writePZBetaSe(out, 0, pG, zG, Score_G / Var_G, 1.0 / sdG);
            }

            // P_GxE per env via SPAGRMClass::getMarkerPvalFromScore on the
            // centred combined-residual score; scoreVar captured for BETA/SE.
            for (int e = 0; e < nE; ++e) {
                const double envResidSum = m_envs[e].spagrm_combined.residSum();
                const double centeredScore = scores(1 + e, b) - gMean * envResidSum;
                double zGxE = 0.0;
                double scoreVar = 0.0;
                const double pGxE = m_envs[e].spagrm_combined
                                        .getMarkerPvalFromScore(centeredScore, altFreq, zGxE, &scoreVar);
                if (scoreVar > 0.0) {
                    const double sdE = std::sqrt(scoreVar);
                    writePZBetaSe(out, 4 + 4 * e, pGxE, zGxE, centeredScore / scoreVar, 1.0 / sdE);
                } else {
                    out[4 + 4 * e + 0] = pGxE;
                }
            }
        }
    }

  private:
    // Write (P, Z, BETA, SE) into out at offset.  Caller passes pre-computed
    // values so the same helper serves both the normal-approx G main effect
    // (BETA = S/Var, SE = 1/sqrt(Var)) and the SPA G×E path (Z from SPA, BETA
    // / SE from the nominal score variance returned by getMarkerPvalFromScore).
    static void writePZBetaSe(
        std::vector<double> &out,
        size_t offset,
        double p, double z, double beta, double se
    ) {
        out[offset + 0] = p;
        out[offset + 1] = z;
        out[offset + 2] = beta;
        out[offset + 3] = se;
    }

    Eigen::VectorXd m_resid_G;
    double m_R_GRM_R_G;
    std::vector<PerEnv> m_envs;
    std::vector<std::string> m_envNames;
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
    const std::string &phenoLabel = {}
) {
    const int nEnv = static_cast<int>(envNames.size());
    if (static_cast<int>(Resid_GxE_list.size()) != nEnv)
        throw std::runtime_error("buildSAGELDArtifacts: envNames / Resid_GxE_list size mismatch");

    const double R_GRM_R_G = computeRGRMR(Resid_G, topo);

    std::vector<SAGELDMethod::PerEnv> envData;
    envData.reserve(nEnv);
    std::vector<std::string> envNamesOut;
    envNamesOut.reserve(nEnv);

    for (int i = 0; i < nEnv; ++i) {
        const auto &Resid_GxE = Resid_GxE_list[i];
        const double lambda = Resid_G.squaredNorm() > 0
                                  ? Resid_GxE.dot(Resid_G) / Resid_G.squaredNorm()
                                  : 0.0;
        if (phenoLabel.empty())
            infoMsg("Env '%s': genome-wide lambda = %.6f", envNames[i].c_str(), lambda);
        else
            infoMsg("[%s] Env '%s': genome-wide lambda = %.6f",
                    phenoLabel.c_str(), envNames[i].c_str(), lambda);
        Eigen::VectorXd envCombined = Resid_GxE - lambda * Resid_G;
        SPAGRMClass spagrm_combined = nsGRMNull::buildSPAGRMNullModel(
            envCombined, N, topo.singletonSet, topo.grmDiag, topo.families,
            topo.familyEntries, topo.allEntries, topo.ibdEntries, topo.ibdPairMap,
            spaCutoff, minMafCutoff, minMacCutoff,
            nsGRMNull::INIT_OUTLIER_RATIO, nsGRMNull::CONTROL_OUTLIER, innerThreads);
        envData.push_back(SAGELDMethod::PerEnv{std::move(spagrm_combined)});
        envNamesOut.push_back(envNames[i]);
    }

    return std::make_unique<SAGELDMethod>(
        Resid_G, R_GRM_R_G, std::move(envData), std::move(envNamesOut));
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
    double hweCutoff
) {
    auto method = buildSAGELDArtifacts(Resid_G, envNames, Resid_GxE_list,
                                       topo, N, spaCutoff, minMafCutoff, minMacCutoff,
                                       /*innerThreads=*/nthreads);
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
    for (int i = 0; i < nEnv; ++i) {
        envNames[i] = envs_spec[i].name;
        Resid_GxE_list[i] = residMat.col(envs_spec[i].colGxE);
    }

    runSAGELDCoreSingle(Resid_G, envNames, Resid_GxE_list, topo, N, *genoData,
                        outputFile, spaCutoff, nthreads,
                        missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);
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
    const std::string &keepFile,
    const std::string &removeFile
) {
    if (envNames.size() != 1)
        throw std::runtime_error("SAGELD pheno mode: --sageld-x currently accepts a single env column");
    for (const auto &en : envNames)
        if (std::find(covarNames.begin(), covarNames.end(), en) == covarNames.end())
            throw std::runtime_error("SAGELD pheno mode: env '" + en +
                                     "' must also appear in --covar-name");

    auto famIIDs = parseGenoIIDs(geno);

    // Build the kept-subject filter: GRM ∩ (keep) − (remove)
    auto grmIDs = SparseGRM::parseSubjectIDs(spgrmGrabFile, spgrmGctaFile, famIIDs);
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
        throw std::runtime_error("SAGELD pheno mode: no subjects survive --keep / --remove / GRM intersection");

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

    auto topo = loadGRMTopology(spgrmGrabFile, spgrmGctaFile, pairwiseIBDFile,
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

    auto buildOne = [&](int p) {
        const std::string &phenoName = phenoNames[p];

        Eigen::VectorXd y = longData.Y.col(p);
        auto fitMain = nsSAGELDFit::fitRandomSlopeML(longData, y);
        infoMsg("[%s] Y ~ X + (E | IID) converged in %d NM iters: "
                "sigma2=%.6g, D=[%.6g, %.6g; %.6g, %.6g]",
                phenoName.c_str(), fitMain.iterations, fitMain.sigma2,
                fitMain.D(0, 0), fitMain.D(0, 1), fitMain.D(1, 0), fitMain.D(1, 1));

        Eigen::VectorXd Resid_G = nsSAGELDFit::aggregatePerIID(longData, fitMain.residPerRow);
        Eigen::VectorXd Resid_GxE = nsSAGELDFit::aggregateWeightedPerIID(
            longData, fitMain.residPerRow, longData.E);

        // ── --save-resid: write per-pheno residual file ─────────────────
        // Format matches src/spagrm/sageld.cpp::parseEnvSpecs (residual-input
        // mode):  #IID  R_G  R_<E>  R_Gx<E>.  R_<E> is the BLUP residual of
        // E ~ 1 + (1|IID), computed only when --save-resid is set (the C++
        // closed-form lambda does not consume R_<E>; we fit only to produce
        // a file compatible with the lme4 reference convention).
        if (saveResid) {
            auto fitE = nsSAGELDFit::fitRandomInterceptML(longData, longData.E);
            Eigen::VectorXd Resid_E = nsSAGELDFit::aggregatePerIID(longData, fitE.residPerRow);
            const std::string residPath = outPrefix + "." + phenoName + ".SAGELD.resid.tsv";
            std::ofstream rf(residPath);
            if (!rf) throw std::runtime_error("Cannot write residual file: " + residPath);
            rf << "#IID\tR_G\tR_" << envNames[0] << "\tR_Gx" << envNames[0] << "\n";
            rf.setf(std::ios::scientific);
            rf.precision(15);
            for (uint32_t i = 0; i < N; ++i)
                rf << longData.uniqueIIDs[i] << "\t"
                   << Resid_G[i] << "\t" << Resid_E[i] << "\t" << Resid_GxE[i] << "\n";
            infoMsg("[%s] Saved residuals: %s (sigma2_E=%.6g, tau2_E=%.6g)",
                    phenoName.c_str(), residPath.c_str(),
                    fitE.sigma2, fitE.D(0, 0));
        }

        // Single env supported here; future multi-env would loop here.
        std::vector<Eigen::VectorXd> Resid_GxE_list;
        Resid_GxE_list.push_back(std::move(Resid_GxE));
        std::vector<std::string> singleEnvName{envNames[0]};

        auto method = buildSAGELDArtifacts(
            Resid_G, singleEnvName, Resid_GxE_list, topo, N,
            spaCutoff, minMafCutoff, minMacCutoff, innerThreads, phenoName);

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

    infoMsg("Running SAGELD marker tests via multiPhenoEngine (%d phenotype(s), %d threads)...",
            nPheno, nthreads);
    multiPhenoEngine(*genoData, tasks, outPrefix, "SAGELD",
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
    const std::string &keepFile,
    const std::string &removeFile
) {
    const bool isResidMode = !residNames.empty();
    const bool isPhenoMode = !phenoNames.empty();

    if (isResidMode && isPhenoMode)
        throw std::runtime_error("SAGELD: --resid-name and --pheno-name are mutually exclusive");
    if (!isResidMode && !isPhenoMode)
        throw std::runtime_error("SAGELD: need either --resid-name (residual mode) or --pheno-name + --sageld-x (pheno mode)");

    if (isResidMode) {
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
        throw std::runtime_error("SAGELD pheno mode: --sageld-x is required");
    if (covarNames.empty())
        throw std::runtime_error("SAGELD pheno mode: --covar-name is required (must include every --sageld-x variable)");
    runSAGELDPhenoMode(phenoFile, phenoNames, covarNames, envNames,
                       spgrmGrabFile, spgrmGctaFile, pairwiseIBDFile, geno,
                       outPrefix, compression, compressionLevel,
                       spaCutoff, nthreads, nSnpPerChunk,
                       missingCutoff, minMafCutoff, minMacCutoff, hweCutoff,
                       saveResid, keepFile, removeFile);
}
