// spamixlocalp.cpp — SPAmixLocalPlus implementation
//
// Phi estimation + per-ancestry GWAS with streaming admix binary I/O.

#include "spamix/spamixlocalp.hpp"
#include "io/admix.hpp"
#include "io/sparse_grm.hpp"
#include "io/subject_data.hpp"
#include "engine/marker.hpp"
#include "spamix/common.hpp"
#include "util/logging.hpp"
#include "util/text_stream.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <condition_variable>
#include <cstdio>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>
#include <boost/math/distributions/normal.hpp>

static constexpr double MIN_P_VALUE = std::numeric_limits<double>::min();


// ======================================================================
// Precomputed phi products — bake R[i]*R[j]*phi*mult into a flat double
// ======================================================================

RprodPhi buildRprodPhi(const PhiMatrices& phi, const Eigen::VectorXd& R) {
    RprodPhi rp;
    auto build = [&](const std::vector<PhiEntry>& src,
                     std::vector<RprodEntry>& dst, double mult) {
        dst.resize(src.size());
        for (size_t p = 0; p < src.size(); ++p) {
            dst[p].i     = src[p].i;
            dst[p].j     = src[p].j;
            dst[p].rprod = mult * src[p].value * R[src[p].i] * R[src[p].j];
        }
    };
    build(phi.A, rp.A, 4.0);
    build(phi.B, rp.B, 2.0);
    build(phi.C, rp.C, 2.0);
    build(phi.D, rp.D, 1.0);
    return rp;
}


// ======================================================================
// Phi estimation — streaming through .abed, no full-matrix materialization
// ======================================================================

static constexpr double PHI_MAF_CUTOFF = 0.01;

PhiMatrices estimatePhiOneAncestry(
    const AdmixData& admixData,
    const SparseGRM& grm,
    int ancIdx,
    int nthreads
) {
    const double mafCutoff = PHI_MAF_CUTOFF;
    uint32_t nUsed = admixData.nSubjUsed();
    uint32_t nMarkers = admixData.nMarkers();
    if (nthreads < 1) nthreads = 1;

    // Build directed pair list from GRM entries (both directions for off-diagonal)
    struct DPair { uint32_t i, j; };
    std::vector<DPair> pairs;
    pairs.reserve(grm.nnz() * 2);
    for (const auto& e : grm.entries()) {
        if (e.row != e.col) {
            pairs.push_back({e.row, e.col});
            pairs.push_back({e.col, e.row});
        }
    }
    const size_t nPairs = pairs.size();

    // Per-thread accumulator structure
    struct PairAcc { double ratioSum; uint32_t validCount; };

    // Worker function: processes markers [mBegin, mEnd) using its own cursor
    // and accumulates into local vectors.
    struct ThreadResult {
        std::vector<PairAcc> accA, accB, accC, accD;
        uint32_t snpsPassed = 0;
    };

    auto workerFn = [&](uint32_t mBegin, uint32_t mEnd) -> ThreadResult {
        ThreadResult res;
        res.accA.assign(nPairs, {0.0, 0});
        res.accB.assign(nPairs, {0.0, 0});
        res.accC.assign(nPairs, {0.0, 0});
        res.accD.assign(nPairs, {0.0, 0});

        auto cursor = admixData.makeCursor();
        cursor->beginSequentialBlock(mBegin);
        Eigen::VectorXd dosage(nUsed), hapcount(nUsed);
        double mafHigh = 1.0 - mafCutoff;
        bool checkMissing = !admixData.hasNoMissing();

        for (uint32_t m = mBegin; m < mEnd; ++m) {
            double q = cursor->getAdmixGenotypes(m, ancIdx, dosage, hapcount);

            if (q <= mafCutoff || q >= mafHigh) continue;
            ++res.snpsPassed;

            double qTerm = q * (1.0 - q);

            for (size_t p = 0; p < nPairs; ++p) {
                uint32_t i = pairs[p].i;
                uint32_t j = pairs[p].j;
                if (i >= nUsed || j >= nUsed) continue;

                double hi = hapcount[i];
                double hj = hapcount[j];
                double gi = dosage[i];
                double gj = dosage[j];

                if (checkMissing) {
                    if (!std::isfinite(hi) || !std::isfinite(hj) ||
                        !std::isfinite(gi) || !std::isfinite(gj)) continue;
                }

                int hiInt = static_cast<int>(std::round(hi));
                int hjInt = static_cast<int>(std::round(hj));

                double denom = hi * hj * qTerm;
                if (std::abs(denom) < 1e-15) continue;

                double numer = (gi - hi * q) * (gj - hj * q);
                double ratio = numer / denom;

                if (hiInt == 2 && hjInt == 2) {
                    res.accA[p].ratioSum += ratio;
                    res.accA[p].validCount++;
                } else if (hiInt == 2 && hjInt == 1) {
                    res.accB[p].ratioSum += ratio;
                    res.accB[p].validCount++;
                } else if (hiInt == 1 && hjInt == 2) {
                    res.accC[p].ratioSum += ratio;
                    res.accC[p].validCount++;
                } else if (hiInt == 1 && hjInt == 1) {
                    res.accD[p].ratioSum += ratio;
                    res.accD[p].validCount++;
                }
            }

            if ((m + 1) % 10000 == 0)
                infoMsg("  Phi estimation: %u / %u markers (anc%d)", m + 1, nMarkers, ancIdx);
        }
        return res;
    };

    // Split markers across threads
    std::vector<ThreadResult> results(nthreads);
    if (nthreads == 1) {
        results[0] = workerFn(0, nMarkers);
    } else {
        std::vector<std::thread> workers;
        workers.reserve(nthreads - 1);
        uint32_t chunkSize = (nMarkers + nthreads - 1) / nthreads;
        for (int t = 0; t < nthreads - 1; ++t) {
            uint32_t mBegin = t * chunkSize;
            uint32_t mEnd   = std::min(mBegin + chunkSize, nMarkers);
            workers.emplace_back([&results, &workerFn, t, mBegin, mEnd]() {
                results[t] = workerFn(mBegin, mEnd);
            });
        }
        // Main thread handles last chunk
        {
            uint32_t mBegin = (nthreads - 1) * chunkSize;
            uint32_t mEnd   = nMarkers;
            results[nthreads - 1] = workerFn(mBegin, mEnd);
        }
        for (auto& w : workers) w.join();
    }

    // Merge accumulators
    std::vector<PairAcc> accA(nPairs, {0.0, 0});
    std::vector<PairAcc> accB(nPairs, {0.0, 0});
    std::vector<PairAcc> accC(nPairs, {0.0, 0});
    std::vector<PairAcc> accD(nPairs, {0.0, 0});
    uint32_t snpsPassed = 0;

    for (int t = 0; t < nthreads; ++t) {
        snpsPassed += results[t].snpsPassed;
        for (size_t p = 0; p < nPairs; ++p) {
            accA[p].ratioSum   += results[t].accA[p].ratioSum;
            accA[p].validCount += results[t].accA[p].validCount;
            accB[p].ratioSum   += results[t].accB[p].ratioSum;
            accB[p].validCount += results[t].accB[p].validCount;
            accC[p].ratioSum   += results[t].accC[p].ratioSum;
            accC[p].validCount += results[t].accC[p].validCount;
            accD[p].ratioSum   += results[t].accD[p].ratioSum;
            accD[p].validCount += results[t].accD[p].validCount;
        }
    }

    infoMsg("  Phi estimation: %u markers passed MAF filter (anc%d)", snpsPassed, ancIdx);

    // Build PhiMatrices from accumulators
    PhiMatrices result;
    auto buildEntries = [&](const std::vector<PairAcc>& acc, std::vector<PhiEntry>& out) {
        for (size_t p = 0; p < nPairs; ++p) {
            if (acc[p].validCount > 0) {
                double phi = acc[p].ratioSum / acc[p].validCount;
                out.push_back({pairs[p].i, pairs[p].j, phi});
            }
        }
    };
    buildEntries(accA, result.A);
    buildEntries(accB, result.B);
    buildEntries(accC, result.C);
    buildEntries(accD, result.D);

    infoMsg("  Phi entries: A=%zu, B=%zu, C=%zu, D=%zu",
            result.A.size(), result.B.size(), result.C.size(), result.D.size());

    return result;
}


// ======================================================================
// Variance computation with phi matrices
// ======================================================================

double computePhiVariance(
    const Eigen::VectorXd& R,
    const Eigen::VectorXd& hapcount,
    double q,
    const PhiMatrices& phi
) {
    double qTerm = q * (1.0 - q);
    double var = 0.0;

    // Off-diagonal terms (phi already has bidirectional pairs)
    // Scenario A: h_i=2, h_j=2 → multiplier = 4
    for (const auto& e : phi.A) {
        if (std::abs(hapcount[e.i] - 2.0) < 0.5 && std::abs(hapcount[e.j] - 2.0) < 0.5)
            var += 4.0 * qTerm * e.value * R[e.i] * R[e.j];
    }
    // Scenario B: h_i=2, h_j=1 → multiplier = 2
    for (const auto& e : phi.B) {
        if (std::abs(hapcount[e.i] - 2.0) < 0.5 && std::abs(hapcount[e.j] - 1.0) < 0.5)
            var += 2.0 * qTerm * e.value * R[e.i] * R[e.j];
    }
    // Scenario C: h_i=1, h_j=2 → multiplier = 2
    for (const auto& e : phi.C) {
        if (std::abs(hapcount[e.i] - 1.0) < 0.5 && std::abs(hapcount[e.j] - 2.0) < 0.5)
            var += 2.0 * qTerm * e.value * R[e.i] * R[e.j];
    }
    // Scenario D: h_i=1, h_j=1 → multiplier = 1
    for (const auto& e : phi.D) {
        if (std::abs(hapcount[e.i] - 1.0) < 0.5 && std::abs(hapcount[e.j] - 1.0) < 0.5)
            var += 1.0 * qTerm * e.value * R[e.i] * R[e.j];
    }

    // Diagonal: sum_i R_i^2 * h_i * q(1-q)
    for (int i = 0; i < static_cast<int>(R.size()); ++i) {
        var += R[i] * R[i] * hapcount[i] * qTerm;
    }

    return var;
}


// ======================================================================
// SPA p-value with outlier split (adapted from reference implementation)
// ======================================================================

namespace {

// Fused CGF derivatives — single exp() per call instead of 3 separate calls.
// K0 = h * log((1-p) + p*e^t),  K1 = h * p*e^t / base,  K2 = h * p*e^t*(1-p) / base^2
inline void kG012Local(double t, double maf, double h,
                       double& K0out, double& K1out, double& K2out) {
    double et   = std::exp(std::clamp(t, -700.0, 700.0));
    double base = (1.0 - maf) + maf * et;
    if (base > 1e-15) {
        double pe   = maf * et;         // p * e^t
        double q1p  = (1.0 - maf);      // (1-p)
        double bsq  = base * base;
        K0out = h * std::log(base);
        K1out = h * pe / base;
        K2out = h * pe * q1p / bsq;
    } else {
        K0out = -std::numeric_limits<double>::infinity();
        K1out = 0.0;
        K2out = 0.0;
    }
}

// Newton-Raphson root finding for K'(t) = s
struct RootResult { double root; bool converged; };

RootResult findRoot(double s,
                    const double* rOut, const double* hOut, int nOut,
                    double q, double meanNorm, double varNorm
) {
    static constexpr double tol = 0.001;
    static constexpr int maxIter = 100;
    double initVals[] = {0.0, -1.0, 1.0, -2.0, 2.0};

    for (double init : initVals) {
        double t = init;
        bool conv = false;
        for (int iter = 0; iter < maxIter; ++iter) {
            double K1 = 0.0, K2 = 0.0;
            for (int i = 0; i < nOut; ++i) {
                double tR = std::clamp(t * rOut[i], -700.0, 700.0);
                double k0i, k1i, k2i;
                kG012Local(tR, q, hOut[i], k0i, k1i, k2i);
                K1 += rOut[i] * k1i;
                K2 += rOut[i] * rOut[i] * k2i;
            }
            double K1total = K1 + meanNorm + varNorm * t - s;
            double K2total = K2 + varNorm;

            if (std::abs(K1total) < tol) { conv = true; break; }
            if (K2total <= 1e-10) break;

            double dt = std::clamp(-K1total / K2total, -2.0, 2.0);
            t = std::clamp(t + dt, -20.0, 20.0);
            if (std::abs(dt) < tol) { conv = true; break; }
        }
        if (conv) return {t, true};
    }
    return {0.0, false};
}

// Lugannani-Rice tail probability
double lugannamiRicePval(double zeta, double s, double q,
                         const double* rOut, const double* hOut, int nOut,
                         double meanNorm, double varNorm, bool upperTail
) {
    double K0 = 0.0, K2 = 0.0;
    for (int i = 0; i < nOut; ++i) {
        double tR = std::clamp(zeta * rOut[i], -700.0, 700.0);
        double k0i, k1i, k2i;
        kG012Local(tR, q, hOut[i], k0i, k1i, k2i);
        K0 += k0i;
        K2 += rOut[i] * rOut[i] * k2i;
    }
    K0 += meanNorm * zeta + 0.5 * varNorm * zeta * zeta;
    K2 += varNorm;

    double temp = zeta * s - K0;
    if (temp <= 0 || K2 <= 0) {
        // Fallback to normal
        double z = std::abs(s) / std::sqrt(K2 > 0 ? K2 : 1.0);
        boost::math::normal_distribution<> norm;
        return boost::math::cdf(boost::math::complement(norm, z));
    }

    double w = std::copysign(std::sqrt(2.0 * temp), zeta);
    double v = zeta * std::sqrt(K2);

    if (std::abs(w) < 1e-12 || std::abs(v) < 1e-12 || !std::isfinite(w) || !std::isfinite(v))
        return MIN_P_VALUE;

    double lr_arg = w + (1.0 / w) * std::log(v / w);
    if (!std::isfinite(lr_arg)) return MIN_P_VALUE;

    boost::math::normal_distribution<> norm;
    if (upperTail) {
        return boost::math::cdf(boost::math::complement(norm, lr_arg));
    } else {
        return boost::math::cdf(norm, lr_arg);
    }
}

} // anonymous namespace

std::pair<double, double> spaLocalPval(
    double S,
    double sMean,
    double varDiag,
    const Eigen::VectorXd& R,
    const Eigen::VectorXd& hapcount,
    double q,
    double varS,
    const OutlierData& outlier,
    double spaCutoff
) {
    double z = (varS > 0.0) ? (S - sMean) / std::sqrt(varS) : 0.0;
    boost::math::normal_distribution<> norm;
    double pNorm = 2.0 * boost::math::cdf(boost::math::complement(norm, std::abs(z)));
    if (pNorm <= 0.0) pNorm = MIN_P_VALUE;

    if (std::abs(z) < spaCutoff || outlier.posOutlier.empty()) {
        return {pNorm, pNorm};
    }

    // varDiag already computed by caller
    double varRatio = (varS > 0.0) ? varDiag / varS : 1.0;
    double sNew = S * std::sqrt(varRatio);
    double sMeanNew = sMean * std::sqrt(varRatio);

    // Extract outlier data
    int nOut = static_cast<int>(outlier.posOutlier.size());
    std::vector<double> rOut(nOut), hOut(nOut);
    double meanNorm = 0.0, varNorm = 0.0;

    for (int i = 0; i < nOut; ++i) {
        uint32_t idx = outlier.posOutlier[i];
        rOut[i] = R[idx];
        hOut[i] = hapcount[idx];
    }
    for (uint32_t idx : outlier.posNonOutlier) {
        meanNorm += R[idx] * q * hapcount[idx];
        varNorm  += R[idx] * R[idx] * q * (1.0 - q) * hapcount[idx];
    }

    double sUpper = std::max(sNew, 2.0 * sMeanNew - sNew);
    double sLower = std::min(sNew, 2.0 * sMeanNew - sNew);

    auto rootUpper = findRoot(sUpper, rOut.data(), hOut.data(), nOut, q, meanNorm, varNorm);
    auto rootLower = findRoot(sLower, rOut.data(), hOut.data(), nOut, q, meanNorm, varNorm);

    double pval = 0.0;
    if (rootUpper.converged) {
        pval += lugannamiRicePval(rootUpper.root, sUpper, q, rOut.data(), hOut.data(), nOut,
                                  meanNorm, varNorm, true);
    }
    if (rootLower.converged) {
        pval += lugannamiRicePval(rootLower.root, sLower, q, rOut.data(), hOut.data(), nOut,
                                  meanNorm, varNorm, false);
    }

    if (!rootUpper.converged && !rootLower.converged)
        pval = std::numeric_limits<double>::quiet_NaN();

    if (std::isfinite(pval)) {
        if (pval <= 0.0) pval = MIN_P_VALUE;
        if (pval > 1.0) pval = 1.0;
    }

    return {pval, pNorm};
}


// ======================================================================
// Phi I/O — wide format (single file, all ancestries)
//
// Header: idx1\tidx2\tanc0_A\tanc0_B\tanc0_C\tanc0_D\tanc1_A\t...
// Indices are 0-based into .fam order.  No .id sidecar needed.
// ======================================================================

static void writePhiWide(const std::string& path,
                         const std::vector<PhiMatrices>& allPhi,
                         int K) {
    // Collect union of all (i,j) pairs across ancestries and scenarios
    struct PairKey { uint32_t i, j; };
    auto pairHash = [](const PairKey& p) {
        return std::hash<uint64_t>{}(uint64_t(p.i) << 32 | p.j);
    };
    auto pairEq = [](const PairKey& a, const PairKey& b) {
        return a.i == b.i && a.j == b.j;
    };
    std::unordered_map<PairKey, std::vector<double>, decltype(pairHash), decltype(pairEq)>
        rows(0, pairHash, pairEq);

    int nCols = K * 4;  // A,B,C,D per ancestry

    auto insertEntries = [&](const std::vector<PhiEntry>& entries, int colIdx) {
        for (const auto& e : entries) {
            auto [it, inserted] = rows.try_emplace({e.i, e.j},
                std::vector<double>(nCols, std::numeric_limits<double>::quiet_NaN()));
            it->second[colIdx] = e.value;
        }
    };

    for (int k = 0; k < K; ++k) {
        insertEntries(allPhi[k].A, k * 4 + 0);
        insertEntries(allPhi[k].B, k * 4 + 1);
        insertEntries(allPhi[k].C, k * 4 + 2);
        insertEntries(allPhi[k].D, k * 4 + 3);
    }

    // Write
    TextWriter out(path);

    // Header
    std::string hdr = "idx1\tidx2";
    const char* scenarioTag[] = {"_A", "_B", "_C", "_D"};
    for (int k = 0; k < K; ++k)
        for (int s = 0; s < 4; ++s)
            hdr += std::string("\tanc") + std::to_string(k) + scenarioTag[s];
    hdr += '\n';
    out.write(hdr);

    // Data rows
    char buf[64];
    for (const auto& [pair, vals] : rows) {
        std::string line;
        line.reserve(64 + 16 * nCols);
        std::snprintf(buf, sizeof(buf), "%u\t%u", pair.i, pair.j);
        line += buf;
        for (int c = 0; c < nCols; ++c) {
            if (std::isnan(vals[c]))
                line += "\tNA";
            else {
                std::snprintf(buf, sizeof(buf), "\t%.17g", vals[c]);
                line += buf;
            }
        }
        line += '\n';
        out.write(line);
    }

    infoMsg("Phi written: %zu pairs x %d columns -> %s",
            rows.size(), nCols, path.c_str());
}

static std::vector<PhiMatrices> readPhiWide(const std::string& path, int K) {
    TextReader reader(path);

    // Parse header to determine column mapping
    std::string header;
    if (!reader.getline(header))
        throw std::runtime_error("Cannot read phi file header: " + path);

    int nCols = K * 4;
    std::vector<PhiMatrices> result(K);

    std::string line;
    while (reader.getline(line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        uint32_t i, j;
        ss >> i >> j;

        for (int c = 0; c < nCols; ++c) {
            std::string tok;
            ss >> tok;
            if (tok == "NA" || tok.empty()) continue;
            double val = std::stod(tok);
            int k = c / 4;
            int s = c % 4;
            PhiEntry entry{i, j, val};
            switch (s) {
                case 0: result[k].A.push_back(entry); break;
                case 1: result[k].B.push_back(entry); break;
                case 2: result[k].C.push_back(entry); break;
                case 3: result[k].D.push_back(entry); break;
            }
        }
    }

    for (int k = 0; k < K; ++k)
        infoMsg("  anc%d phi: A=%zu B=%zu C=%zu D=%zu",
                k, result[k].A.size(), result[k].B.size(),
                result[k].C.size(), result[k].D.size());

    return result;
}


// ======================================================================
// runPhiEstimation — full pipeline (writes single wide phi file)
// ======================================================================

void runPhiEstimation(
    const std::string& admixPrefix,
    const std::string& grmGrabFile,
    const std::string& grmGctaFile,
    const std::string& phiOutputFile,
    const std::string& keepFile,
    const std::string& removeFile,
    const std::string& extractFile,
    const std::string& excludeFile,
    int nthreads
) {
    infoMsg("=== SPAmixLocalPlus Phi Estimation ===");

    // Load .fam IIDs (genotype subject list)
    auto famIIDs = parseFamIIDs(admixPrefix + ".fam");
    uint32_t nFam = static_cast<uint32_t>(famIIDs.size());
    infoMsg("Genotype (.fam): %u subjects", nFam);

    // Load GRM against full .fam order
    auto grm = SparseGRM::load(grmGrabFile, grmGctaFile, famIIDs, famIIDs);
    infoMsg("GRM loaded: %u subjects (dimension), %zu entries", grm.nSubjects(), grm.nnz());

    // Determine GRM subjects: those with a non-zero diagonal entry
    const auto& diag = grm.diagonal();
    uint32_t nGrm = 0;
    for (uint32_t i = 0; i < nFam; ++i)
        if (diag[i] != 0.0) ++nGrm;
    infoMsg("GRM: %u subjects with non-zero diagonal", nGrm);

    // Build intersection bitmask: only subjects present in both genotype and GRM
    uint32_t nMaskWords = (nFam + 63) / 64;
    std::vector<uint64_t> usedMask(nMaskWords, 0);
    uint32_t nUsed = 0;
    for (uint32_t i = 0; i < nFam; ++i) {
        if (diag[i] != 0.0) {
            usedMask[i / 64] |= (uint64_t(1) << (i % 64));
            ++nUsed;
        }
    }
    infoMsg("Intersection (genotype ∩ GRM): %u subjects will be analysed", nUsed);

    if (nUsed == 0)
        throw std::runtime_error("runPhiEstimation: intersection of genotype and GRM subjects is empty");

    // Load admix data using intersection mask
    AdmixData admixData(admixPrefix, usedMask, nFam, nUsed,
                        extractFile, excludeFile);
    int K = admixData.nAncestries();
    infoMsg("Ancestries: %d, Markers: %u, Samples: %u", K, admixData.nMarkers(), admixData.nSubjUsed());

    // Estimate phi for all ancestries
    std::vector<PhiMatrices> allPhi(K);
    for (int k = 0; k < K; ++k) {
        infoMsg("Estimating phi for anc%d (%d thread%s)...", k, nthreads, nthreads > 1 ? "s" : "");
        allPhi[k] = estimatePhiOneAncestry(admixData, grm, k, nthreads);
    }

    // Write single wide file
    writePhiWide(phiOutputFile, allPhi, K);

    infoMsg("Phi estimation complete -> %s", phiOutputFile.c_str());
}


// ======================================================================
// CCT (Cauchy Combination Test) — combine per-ancestry p-values
// ======================================================================

static double cauchyCombine(const double* pvals, int K) {
    int nValid = 0;
    double tStat = 0.0;
    for (int k = 0; k < K; ++k) {
        if (!std::isfinite(pvals[k]) || pvals[k] < 0.0) continue;
        if (pvals[k] <= 0.0) return 0.0;    // any exact zero ⇒ combined = 0
        double pc = (pvals[k] >= 1.0) ? 0.999 : pvals[k];
        tStat += std::tan((0.5 - pc) * M_PI);
        ++nValid;
    }
    if (nValid == 0) return std::numeric_limits<double>::quiet_NaN();
    tStat /= static_cast<double>(nValid);
    return (tStat > 1e15) ? (1.0 / tStat) / M_PI
                          : 0.5 - std::atan(tStat) / M_PI;
}


// ======================================================================
// Unified GWAS — process all K ancestries per marker, one output file
// ======================================================================

static void runUnifiedGWAS(
    const AdmixData& admixData,
    const Eigen::VectorXd& resid,
    const std::vector<PhiMatrices>& allPhi,
    const OutlierData& outlier,
    double spaCutoff,
    double mafCutoff,
    double macCutoff,
    const std::string& outputFile,
    const std::string& compression,
    int compressionLevel,
    int nthreads
) {
    const int K = admixData.nAncestries();
    uint32_t nUsed = admixData.nSubjUsed();
    const auto& markerInfo = admixData.markerInfo();

    // Build header: CHROM POS ID REF ALT  P_CCT  anc0_... anc1_...
    std::string header = "CHROM\tPOS\tID\tREF\tALT\tP_CCT";
    for (int k = 0; k < K; ++k) {
        std::string pfx = "\tanc" + std::to_string(k) + "_";
        header += pfx + "AltFreq";
        header += pfx + "MissingRate";
        header += pfx + "P";
        header += pfx + "Pnorm";
        header += pfx + "Stat";
        header += pfx + "Var";
        header += pfx + "zScore";
        header += pfx + "AltCounts";
        header += pfx + "BetaG";
    }
    header += '\n';

    TextWriter out(outputFile,
                   TextWriter::modeFromString(compression),
                   compressionLevel);
    out.write(header);

    // Parallel processing using chunks
    const auto& chunks = admixData.chunkIndices();
    std::atomic<size_t> nextChunk{0};
    std::atomic<size_t> chunksCompleted{0};
    size_t nChunks = chunks.size();
    uint32_t nMarkers = admixData.nMarkers();

    // Pre-compute R^2 and rprod arrays (once per phenotype, not per marker).
    // This bakes R[i]*R[j]*phi*multiplier into a flat double, eliminating
    // random R[] lookups in the per-marker variance hot path.
    Eigen::ArrayXd R2 = resid.array().square();
    std::vector<RprodPhi> rphi(K);
    for (int k = 0; k < K; ++k)
        rphi[k] = buildRprodPhi(allPhi[k], resid);

    const bool noMissing = admixData.hasNoMissing();

    struct PaddedFlag { alignas(64) int ready = 0; };
    std::vector<std::string> chunkOutput(nChunks);
    std::vector<PaddedFlag> chunkReady(nChunks);

    std::mutex writeMutex;
    std::condition_variable writeCv;
    bool stopWriter = false;

    // Writer thread
    std::thread writer([&]() {
        for (size_t ci = 0; ci < nChunks; ++ci) {
            std::unique_lock<std::mutex> lk(writeMutex);
            writeCv.wait(lk, [&] { return chunkReady[ci].ready || stopWriter; });
            if (chunkReady[ci].ready)
                out.write(chunkOutput[ci]);
        }
    });

    // Worker function
    auto workerFn = [&]() {
        auto cursor = admixData.makeCursor();
        Eigen::VectorXd dosage(nUsed), hapcount(nUsed);
        // Per-ancestry scratch for all-at-once decoding
        Eigen::MatrixXd dosMatrix(nUsed, K), hapMatrix(nUsed, K);
        char fmtBuf[64];
        std::vector<double> ancPvals(K);
        // Thread-local scratch: compact integer hapcount for branchless variance
        std::vector<uint8_t> hInt(nUsed);

        for (size_t ci = nextChunk.fetch_add(1); ci < nChunks; ci = nextChunk.fetch_add(1)) {
            const auto& gIndices = chunks[ci];
            // Estimate ~(5 + 9*K + 20) chars per column per line
            std::string buf;
            buf.reserve(gIndices.size() * (80 + 90 * K));

            cursor->beginSequentialBlock(gIndices.front());

            for (size_t mi = 0; mi < gIndices.size(); ++mi) {
                uint64_t localIdx = gIndices[mi];
                const auto& mInfo = markerInfo[localIdx];

                // Decode all K ancestries at once
                cursor->getAllAncestries(localIdx, dosMatrix, hapMatrix);

                // Marker meta prefix: CHROM POS ID REF ALT
                buf += mInfo.chrom;
                buf += '\t';
                std::snprintf(fmtBuf, sizeof(fmtBuf), "%u", mInfo.pos);
                buf += fmtBuf;
                buf += '\t';
                buf += mInfo.id;
                buf += '\t';
                buf += mInfo.ref;
                buf += '\t';
                buf += mInfo.alt;

                // Per-ancestry results (collected first, then CCT, then written)
                struct AncResult {
                    double maf, missRate, pSpa, pNorm, S, varS, z, dosSum, betaG;
                    bool pass;
                };
                std::vector<AncResult> ancRes(K);

                for (int k = 0; k < K; ++k) {
                    auto dosCol = dosMatrix.col(k);
                    auto hapCol = hapMatrix.col(k);

                    // QC: accumulate dosage/hapcount sums, handle missing
                    double hapSum = 0.0, dosSum = 0.0;
                    uint32_t nMissing = 0;
                    if (noMissing) {
                        // Fast path: no NaN values possible, skip isfinite checks
                        dosSum = dosCol.sum();
                        hapSum = hapCol.sum();
                    } else {
                        for (uint32_t s = 0; s < nUsed; ++s) {
                            if (!std::isfinite(dosCol[s]) || !std::isfinite(hapCol[s])) {
                                ++nMissing;
                                dosCol[s] = 0.0;
                                hapCol[s] = 0.0;
                            } else {
                                dosSum += dosCol[s];
                                hapSum += hapCol[s];
                            }
                        }
                    }

                    double missRate = static_cast<double>(nMissing) / nUsed;
                    double maf = (hapSum > 0) ? dosSum / hapSum : 0.0;
                    double mac = std::min(dosSum, hapSum - dosSum);

                    AncResult& ar = ancRes[k];
                    ar.maf = maf;
                    ar.missRate = missRate;
                    ar.dosSum = dosSum;

                    if (maf < mafCutoff || maf > (1.0 - mafCutoff) || mac < macCutoff) {
                        ar.pass = false;
                        ar.pSpa = std::numeric_limits<double>::quiet_NaN();
                        ancPvals[k] = std::numeric_limits<double>::quiet_NaN();
                        continue;
                    }
                    ar.pass = true;

                    double q = maf;
                    double qTerm = q * (1.0 - q);

                    // Score test — Eigen SIMD dot products (replaces scalar loop)
                    double S = dosCol.dot(resid);
                    double sMean = q * hapCol.dot(resid);

                    // ── Variance via precomputed rprod + compact hInt ──
                    // Convert hapcount to uint8 (values are exactly 0, 1, or 2)
                    if (noMissing) {
                        for (uint32_t s = 0; s < nUsed; ++s)
                            hInt[s] = static_cast<uint8_t>(hapCol[s]);
                    } else {
                        for (uint32_t s = 0; s < nUsed; ++s) {
                            double h = hapCol[s];
                            hInt[s] = std::isfinite(h) ? static_cast<uint8_t>(h) : 0;
                        }
                    }

                    // Off-diagonal phi variance — branchless accumulation.
                    // Replaces: abs(hapcount[i]-target)<0.5 branch + random R[] reads
                    // with: integer equality (branchless) + precomputed rprod (sequential).
                    const RprodPhi& rp = rphi[k];
                    double varOff = 0.0;
                    for (size_t p = 0; p < rp.A.size(); ++p) {
                        const auto& e = rp.A[p];
                        varOff += ((hInt[e.i] == 2) & (hInt[e.j] == 2)) * e.rprod;
                    }
                    for (size_t p = 0; p < rp.B.size(); ++p) {
                        const auto& e = rp.B[p];
                        varOff += ((hInt[e.i] == 2) & (hInt[e.j] == 1)) * e.rprod;
                    }
                    for (size_t p = 0; p < rp.C.size(); ++p) {
                        const auto& e = rp.C[p];
                        varOff += ((hInt[e.i] == 1) & (hInt[e.j] == 2)) * e.rprod;
                    }
                    for (size_t p = 0; p < rp.D.size(); ++p) {
                        const auto& e = rp.D[p];
                        varOff += ((hInt[e.i] == 1) & (hInt[e.j] == 1)) * e.rprod;
                    }

                    // Diagonal variance — vectorized dot product
                    double varDiag = qTerm * R2.matrix().dot(hapCol);
                    double varS = varOff * qTerm + varDiag;

                    double z = (varS > 0.0) ? (S - sMean) / std::sqrt(varS) : 0.0;
                    auto [pSpa, pNorm] = spaLocalPval(S, sMean, varDiag, resid, hapCol, q, varS, outlier, spaCutoff);
                    double betaG = (varS > 0.0) ? (S - sMean) / varS
                                                : std::numeric_limits<double>::quiet_NaN();

                    ar.S = S;
                    ar.varS = varS;
                    ar.z = z;
                    ar.pSpa = pSpa;
                    ar.pNorm = pNorm;
                    ar.betaG = betaG;
                    ancPvals[k] = pSpa;
                }

                // CCT across ancestries
                double pCCT = cauchyCombine(ancPvals.data(), K);

                // Write P_CCT
                buf += '\t';
                if (std::isfinite(pCCT)) {
                    std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", pCCT);
                    buf += fmtBuf;
                } else {
                    buf += "NA";
                }

                // Write per-ancestry columns
                for (int k = 0; k < K; ++k) {
                    const AncResult& ar = ancRes[k];
                    // AltFreq
                    buf += '\t';
                    std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", ar.maf);
                    buf += fmtBuf;
                    // MissingRate
                    buf += '\t';
                    std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", ar.missRate);
                    buf += fmtBuf;

                    if (!ar.pass) {
                        buf += "\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
                        continue;
                    }
                    // P
                    buf += '\t';
                    std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", ar.pSpa);
                    buf += fmtBuf;
                    // Pnorm
                    buf += '\t';
                    std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", ar.pNorm);
                    buf += fmtBuf;
                    // Stat
                    buf += '\t';
                    std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", ar.S);
                    buf += fmtBuf;
                    // Var
                    buf += '\t';
                    std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", ar.varS);
                    buf += fmtBuf;
                    // zScore
                    buf += '\t';
                    std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", ar.z);
                    buf += fmtBuf;
                    // AltCounts
                    buf += '\t';
                    std::snprintf(fmtBuf, sizeof(fmtBuf), "%.0f", ar.dosSum);
                    buf += fmtBuf;
                    // BetaG
                    buf += '\t';
                    if (std::isfinite(ar.betaG)) {
                        std::snprintf(fmtBuf, sizeof(fmtBuf), "%.6g", ar.betaG);
                        buf += fmtBuf;
                    } else {
                        buf += "NA";
                    }
                }
                buf += '\n';
            }

            {
                std::lock_guard<std::mutex> lk(writeMutex);
                chunkOutput[ci] = std::move(buf);
                chunkReady[ci].ready = 1;
            }
            writeCv.notify_all();

            // Progress logging: report at ~25%, 50%, 75%
            size_t done = chunksCompleted.fetch_add(1) + 1;
            if (nChunks >= 20) {
                size_t q1 = nChunks / 4, q2 = nChunks / 2, q3 = nChunks * 3 / 4;
                if (done == q1 || done == q2 || done == q3) {
                    uint32_t markersDone = static_cast<uint32_t>(
                        static_cast<uint64_t>(done) * nMarkers / nChunks);
                    infoMsg("    %u / %u markers (~%u%%)",
                            markersDone, nMarkers,
                            static_cast<unsigned>(done * 100 / nChunks));
                }
            }
        }
    };

    // Launch workers
    int nWorkers = std::max(1, nthreads);
    std::vector<std::thread> workers;
    workers.reserve(nWorkers);
    for (int t = 0; t < nWorkers; ++t)
        workers.emplace_back(workerFn);

    for (auto& w : workers) w.join();
    {
        std::lock_guard<std::mutex> lk(writeMutex);
        stopWriter = true;
    }
    writeCv.notify_all();
    writer.join();

    infoMsg("  %u markers processed -> %s", admixData.nMarkers(), outputFile.c_str());
}


// ======================================================================
// runSPAmixLocalPlus — main entry point
// ======================================================================

void runSPAmixLocalPlus(
    const std::string& residFile,
    const std::string& admixPrefix,
    const std::string& admixPhiFile,
    const std::string& outPrefix,
    const std::string& compression,
    int compressionLevel,
    double spaCutoff,
    double outlierRatio,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    const std::string& keepFile,
    const std::string& removeFile,
    const std::string& extractFile,
    const std::string& excludeFile
) {
    infoMsg("=== SPAmixLocalPlus GWAS ===");

    // Load subjects
    auto famIIDs = parseFamIIDs(admixPrefix + ".fam");
    SubjectData sd(famIIDs);
    sd.loadResidOne(residFile);
    sd.setKeepRemove(keepFile, removeFile);
    sd.finalize();

    uint32_t nUsed = sd.nUsed();
    infoMsg("Subjects: %u in .fam, %u used", sd.nFam(), nUsed);

    // Load admix data
    AdmixData admixData(admixPrefix, sd.usedMask(), sd.nFam(), nUsed,
                        extractFile, excludeFile, nSnpPerChunk);
    int K = admixData.nAncestries();
    infoMsg("Ancestries: %d, Markers: %u", K, admixData.nMarkers());

    // Load phi from wide file (all ancestries at once)
    auto allPhi = readPhiWide(admixPhiFile, K);

    // Per-residual-column loop
    const int nRC = sd.residOneCols();
    if (nRC > 1) infoMsg("Multi-column residual file: %d phenotypes", nRC);

    auto phenoInfos = sd.buildPerColumnMasks();

    for (int rc = 0; rc < nRC; ++rc) {
        const auto& pi = phenoInfos[rc];

        // Build union-dimension residual with 0 for missing subjects
        Eigen::VectorXd colResid;
        if (nRC > 1) {
            colResid = sd.residMatrix().col(rc);
            for (Eigen::Index s = 0; s < colResid.size(); ++s)
                if (std::isnan(colResid[s])) colResid[s] = 0.0;
        }
        const Eigen::VectorXd& resid = (nRC > 1) ? colResid : sd.residuals();

        std::string outFile = TextWriter::buildOutputPath(
            outPrefix, pi.name, "SPAmixLocalPlus", compression);

        infoMsg("  Phenotype '%s': %u subjects, %u markers, %d ancestries -> %s",
                pi.name.c_str(), pi.nUsed, admixData.nMarkers(), K, outFile.c_str());

        // Detect outliers per residual column
        OutlierData outlier = detectOutliers(resid, outlierRatio);

        runUnifiedGWAS(admixData, resid, allPhi, outlier,
                       spaCutoff, minMafCutoff, minMacCutoff,
                       outFile, compression, compressionLevel, nthread);
    }

    infoMsg("SPAmixLocalPlus GWAS complete.");
}
