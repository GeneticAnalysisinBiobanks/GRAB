// spagxe.cpp — SPAGxE_CCT G×E test implementation.
//
// Phase 0 scaffold: the MethodBase interface compiles and registers, but the
// per-marker Branch A / Branch B logic and the runSPAGxE orchestration are
// filled in by Phase 1 (base SPA-only) and Phase 2 (sparse-GRM +).

#include "spagxe/spagxe.hpp"

#include <limits>
#include <stdexcept>
#include <utility>

// ══════════════════════════════════════════════════════════════════════
// SPAGxEMethod
// ══════════════════════════════════════════════════════════════════════

SPAGxEMethod::SPAGxEMethod(
    Eigen::VectorXd resid,
    std::vector<std::string> envNames,
    double marginalCutoff,
    double spaCutoff
)
    : m_resid(std::move(resid)),
      m_envNames(std::move(envNames)),
      m_marginalCutoff(marginalCutoff),
      m_spaCutoff(spaCutoff)
{
}

std::unique_ptr<MethodBase> SPAGxEMethod::clone() const {
    return std::make_unique<SPAGxEMethod>(*this);
}

std::string SPAGxEMethod::getHeaderColumns() const {
    // Leading marginal block P_G Z_G BETA_G SE_G (normal-approx, so its Z is
    // already p-consistent and needs no Z_Norm), then a 5-wide G×E block per
    // environment (SPA score test → p-consistent Z and raw-score Z_Norm differ
    // in the tails).
    std::string h = "\tP_G\tZ_G\tBETA_G\tSE_G";
    for (const auto &n : m_envNames)
        h += "\tP_Gx" + n + "\tZ_Gx" + n + "\tZ_Norm_Gx" + n + "\tBETA_Gx" + n +
             "\tSE_Gx" + n;
    return h;
}

void SPAGxEMethod::getResultVec(
    Eigen::Ref<Eigen::VectorXd> /*GVec*/,
    double /*altFreq*/,
    int /*markerInChunkIdx*/,
    std::vector<double> &result
) {
    // Phase 0 scaffold: emit NaN placeholders with the correct column count.
    // Phase 1 replaces this with the marginal screen + Branch A / Branch B SPA.
    result.assign(resultSize(), std::numeric_limits<double>::quiet_NaN());
}

// ══════════════════════════════════════════════════════════════════════
// runSPAGxE
// ══════════════════════════════════════════════════════════════════════

void runSPAGxE(
    const std::string & /*phenoFile*/,
    const std::vector<std::string> & /*residNames*/,
    const std::vector<std::string> & /*envNames*/,
    const std::string & /*spgrmGrabFile*/,
    const std::string & /*spgrmGctaFile*/,
    const std::string & /*pairwiseIBDFile*/,
    const GenoSpec & /*geno*/,
    const std::string & /*outPrefix*/,
    const std::string & /*compression*/,
    int /*compressionLevel*/,
    double /*marginalCutoff*/,
    double /*spaCutoff*/,
    double /*outlierIqrRatio*/,
    int /*nthreads*/,
    int /*nSnpPerChunk*/,
    double /*missingCutoff*/,
    double /*minMafCutoff*/,
    double /*minMacCutoff*/,
    double /*hweCutoff*/,
    const std::string & /*keepFile*/,
    const std::string & /*removeFile*/,
    const std::string & /*regressionModelStr*/,
    const std::string & /*phenoNameSpec*/,
    const std::string & /*covarFile*/,
    const std::vector<std::string> & /*covarNames*/,
    bool /*saveResid*/
) {
    throw std::runtime_error(
        "SPAGxE is not yet implemented (Phase 0 scaffold): the base SPA-only "
        "marker test lands in Phase 1.");
}
