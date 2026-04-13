// sageld.cpp — SAGELD: G×E interaction marker test (pure C++17 / Eigen)
//
// Per marker, produces:
//   P_G       — main genetic effect (normal approximation)
//   P_Gx<E>   — G×E interaction (SPA on combined residual per env)
//   Z_G, Z_Gx<E> — corresponding z-scores
//
// Combined residual per env:  Resid_combined = R_Gx<E> - λ * R_G
// λ estimated genome-wide from residual covariance.

#include "spagrm/sageld.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "spagrm/grm_null.hpp"
#include "spagrm/spagrm.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/text_scanner.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <Eigen/Dense>

namespace {

// ══════════════════════════════════════════════════════════════════════
// EnvSpec — describes one G×E environment's column indices + name
// ══════════════════════════════════════════════════════════════════════

struct EnvSpec {
    std::string name; // env name (from header or "E1","E2",...); used in output
    int colE;         // residMat column index of R_<E>
    int colGxE;       // residMat column index of R_Gx<E>
};

// Validate and extract environment specs from the residual file metadata.
//
// Expected column order (nRC = total residual cols, excl. IID):
//   col 0   : R_G
//   col 1   : R_<E1>
//   col 2   : R_Gx<E1>
//   col 3   : R_<E2>        (optional extra envs)
//   col 4   : R_Gx<E2>
//   ...
// nRC must be odd and ≥ 3.
//
// If colNames is empty (no header), parse positionally and name envs E1, E2, ...
std::vector<EnvSpec> parseEnvSpecs(
    int nRC,
    const std::vector<std::string> &colNames
) {
    if (nRC < 3 || (nRC - 1) % 2 != 0)
        throw std::runtime_error("SAGELD residual file: expected odd number of columns ≥ 3 "
                                 "(R_G  R_<E1>  R_Gx<E1>  [R_<E2>  R_Gx<E2>  ...]), "
                                 "got " +
                                 std::to_string(nRC));

    const int nEnv = (nRC - 1) / 2;
    std::vector<EnvSpec> envs(nEnv);

    if (colNames.empty()) {
        // No header — parse by position, name envs "E1", "E2", ...
        infoMsg("Residual file has no header; parsing by column order: "
                "R_G  R_<E1>  R_Gx<E1>  ...");
        for (int i = 0; i < nEnv; ++i)
            envs[i] = {"E" + std::to_string(i + 1), 1 + i * 2, 2 + i * 2};
        return envs;
    }

    // Header present — validate pattern
    if (static_cast<int>(colNames.size()) != nRC)
        throw std::runtime_error("SAGELD: column name count (" + std::to_string(colNames.size()) +
                                 ") does not match data column count (" + std::to_string(nRC) + ")");

    if (colNames[0] != "R_G")
        throw std::runtime_error("SAGELD: first residual column must be 'R_G', got '" + colNames[0] + "'");

    for (int i = 0; i < nEnv; ++i) {
        const int ce = 1 + i * 2; // R_<E>
        const int cg = 2 + i * 2; // R_Gx<E>

        const std::string &eName = colNames[ce];
        const std::string &gxName = colNames[cg];

        // R_<E>: starts with "R_" but NOT "R_Gx"
        const bool eOk =
            (eName.size() >= 3 && eName.substr(0, 2) == "R_" && !(eName.size() >= 4 && eName.substr(0, 4) == "R_Gx"));
        if (!eOk)
            throw std::runtime_error("SAGELD: column " + std::to_string(ce) + " ('" + eName +
                                     "') must match R_<E> (environment residual, e.g. R_AGE)");

        // R_Gx<E>: starts with "R_Gx"
        if (gxName.size() < 5 || gxName.substr(0, 4) != "R_Gx")
            throw std::runtime_error("SAGELD: column " + std::to_string(cg) + " ('" + gxName +
                                     "') must match R_Gx<E> (G×E residual, e.g. R_GxAGE)");

        // Extract env name from E column ("R_AGE" → "AGE")
        const std::string envFromE = eName.substr(2);    // after "R_"
        const std::string envFromGxE = gxName.substr(4); // after "R_Gx"

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

    // Build connected components for family/singleton split
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

// Compute R' * GRM * R using pre-loaded topology.
// Single pass over all entries — O(E), no per-family hash sets.
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
//
// Holds:
//   m_resid_G     — main-effect residual vector (shared across envs)
//   m_R_GRM_R_G   — R_G' * GRM * R_G (for Var_G normal approx)
//   m_envs        — per-env: SPAGRMClass for the combined residual
//   m_envNames    — env name strings for output column headers
//
// Output columns (resultSize = 2 + 2*nEnv):
//   P_G, P_Gx<E1>, ..., Z_G, Z_Gx<E1>, ...
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
        return 2 + 2 * static_cast<int>(m_envs.size());
    }

    // P_G  P_Gx<E1>  ...  Z_G  Z_Gx<E1>  ...
    std::string getHeaderColumns() const override {
        std::string h = "\tP_G";
        for (const auto &n : m_envNames)
            h += "\tP_Gx" + n;
        h += "\tZ_G";
        for (const auto &n : m_envNames)
            h += "\tZ_Gx" + n;
        return h;
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int /*markerInChunkIdx*/,
        std::vector<double> &result
    ) override {
        const int nE = static_cast<int>(m_envs.size());
        result.resize(2 + 2 * nE);

        const double MAF = std::min(altFreq, 1.0 - altFreq);
        const double G_var = 2.0 * MAF * (1.0 - MAF);
        const double Score_G = GVec.dot(m_resid_G) - GVec.mean() * m_resid_G.sum();
        const double Var_G = G_var * m_R_GRM_R_G;
        const double zG = Score_G / std::sqrt(Var_G);
        const double pG = 2.0 * math::pnorm(std::abs(zG), 0.0, 1.0, false);

        result[0] = pG;      // P_G
        result[1 + nE] = zG; // Z_G

        for (int i = 0; i < nE; ++i) {
            double zGxE = 0.0;
            double pGxE = m_envs[i].spagrm_combined.getMarkerPval(GVec, altFreq, zGxE);
            result[1 + i] = pGxE;      // P_Gx<Ei>
            result[2 + nE + i] = zGxE; // Z_Gx<Ei>
        }
    }

  private:
    Eigen::VectorXd m_resid_G;
    double m_R_GRM_R_G;
    std::vector<PerEnv> m_envs;
    std::vector<std::string> m_envNames;
};

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// runSAGELD — entry point
// ══════════════════════════════════════════════════════════════════════

void runSAGELD(
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
    // 1. Load residual file
    infoMsg("Loading pheno file: %s", phenoFile.c_str());
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

    // 2. Parse env specs (validate header / assign positional names)
    const auto envs_spec = parseEnvSpecs(sd.residOneCols(), sd.residColNames());
    const int nEnv = static_cast<int>(envs_spec.size());
    infoMsg("SAGELD: %d environment(s): ", nEnv);
    for (const auto &e : envs_spec)
        infoMsg("  %s (cols %d / %d)", e.name.c_str(), e.colE, e.colGxE);

    // 3. Load GRM topology + IBD
    auto topo = loadGRMTopology(spgrmGrabFile, spgrmGctaFile, pairwiseIBDFile, sd.usedIIDs(), sd.famIIDs(), N);

    // 4. Load genotype data
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);

    // 5. Extract R_G (col 0) and compute its GRM quadratic form once
    const Eigen::MatrixXd &residMat = sd.residMatrix();
    const Eigen::VectorXd Resid_G = residMat.col(0);
    const double R_GRM_R_G = computeRGRMR(Resid_G, topo);

    // 6. Per-env: compute lambda, build null model for combined residual
    std::vector<SAGELDMethod::PerEnv> envData;
    std::vector<std::string> envNames;
    envData.reserve(nEnv);
    envNames.reserve(nEnv);

    for (const auto &es : envs_spec) {
        infoMsg("Building null model for environment '%s'...", es.name.c_str());

        const Eigen::VectorXd Resid_GxE = residMat.col(es.colGxE);

        // Genome-wide lambda: λ = (R_GxE · R_G) / (R_G · R_G)
        const double lambda = Resid_G.squaredNorm() > 0 ? Resid_GxE.dot(Resid_G) / Resid_G.squaredNorm() : 0.0;
        infoMsg("  Genome-wide lambda = %.6f", lambda);

        const Eigen::VectorXd Resid_combined = Resid_GxE - lambda * Resid_G;

        SPAGRMClass spagrm_combined = nsGRMNull::buildSPAGRMNullModel(
            Resid_combined, N, topo.singletonSet, topo.grmDiag, topo.families, topo.familyEntries, topo.allEntries,
            topo.ibdEntries, topo.ibdPairMap, spaCutoff, minMafCutoff, minMacCutoff, nthreads);

        envData.push_back({std::move(spagrm_combined)});
        envNames.push_back(es.name);
    }

    // 7. Build method and run single marker-engine pass
    SAGELDMethod method(Resid_G, R_GRM_R_G, std::move(envData), std::move(envNames));

    infoMsg("Starting SAGELD marker-level association (%d threads, %d env(s))...", nthreads, nEnv);
    markerEngine(*genoData, method, outputFile, nthreads, missingCutoff, minMafCutoff, minMacCutoff, hweCutoff);
}
