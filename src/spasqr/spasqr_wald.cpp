// spasqr_wald.cpp — SPAsqr Wald mode (full-model effect size + SE)
//
// For every (marker, τ) pair, this mode refits the joint smoothed-QR model
// on Z = [1 | X | G] using QMME, then computes the M-estimation sandwich
//      V = A^{-1} B A^{-1} / n,
// where
//      A = (1/n) Σ K_h(-e_i)   Z_i Z_i^T   (smooth-QR Hessian)
//      B = (1/n) Σ R_i^2       Z_i Z_i^T   (R_i = τ - Φ(-e_i/h))
// Effect size β̂_G is the last entry of θ̂ and SE = sqrt(V[γγ]).
//
// Differences from score-mode runSPAsqr / runSPAsqrLoco:
//   - Refit per (marker, τ); no null-model caching.
//   - No GRM (point estimation, not score test).
//   - Output is long-format summary stats; one row per (marker, τ):
//       CHROM POS ID REF ALT MISS_RATE ALT_FREQ MAC HWE_P TAU BETA SE Z P
//   - --pred-list (optional): per-chromosome y_resp = Y_transformed - loco_chr,
//     same offset convention as runSPAsqrLoco (β=1, α=0).

#include "spasqr/spasqr.hpp"
#include "spasqr/qmme.hpp"

#include "engine/loco.hpp"
#include "engine/marker_impl.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "io/sparse_grm.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

// Phenotype pre-transform (mirrors helper in spasqr.cpp; small enough to copy).
void applyPhenoTransform(Eigen::VectorXd &Y, const std::string &mode) {
    if (mode == "raw") return;
    if (mode == "int") {
        Y = math::inverseRankNormal(Y);
        return;
    }
    if (mode == "standardize") {
        const Eigen::Index n = Y.size();
        if (n <= 1) return;
        const double mean = Y.mean();
        const double ssq  = (Y.array() - mean).square().sum();
        const double sd   = std::sqrt(ssq / static_cast<double>(n - 1));
        if (sd > 0.0) Y = (Y.array() - mean) / sd;
        return;
    }
    throw std::runtime_error("applyPhenoTransform: unknown mode '" + mode + "'");
}

// IQR-based bandwidth: h = IQR(Y) / scale.  Caller picks scale (Wald
// defaults to 10; score mode uses 3).  Falls back to a small positive
// constant if IQR == 0.
double iqrBandwidth(const Eigen::VectorXd &Y, double scale) {
    const Eigen::Index n = Y.size();
    if (n < 4) return 1.0;
    std::vector<double> v(static_cast<size_t>(n));
    Eigen::VectorXd::Map(v.data(), n) = Y;
    std::sort(v.begin(), v.end());
    auto q = [&](double prob) {
        const double idx = prob * (static_cast<double>(n) - 1);
        const Eigen::Index lo = static_cast<Eigen::Index>(std::floor(idx));
        const Eigen::Index hi = std::min(lo + 1, n - 1);
        const double frac = idx - lo;
        return v[lo] * (1.0 - frac) + v[hi] * frac;
    };
    const double iqr = q(0.75) - q(0.25);
    const double h   = iqr / scale;
    if (h > 0.0) return h;
    return std::max(std::pow(std::log(static_cast<double>(n)) / static_cast<double>(n), 0.4), 0.05);
}

// Compute Wald β̂_G + SE for one marker × one tau.  Inputs are pheno-dense:
// y_resp (response), X (covariates with no intercept; QMME prepends one),
// G (genotype, no NaN).
//
// `sandwichOk` reports whether the SE is numerically valid (LDLT(A) succeeded
// and V_γγ > 0).  `qmmeTolMet` reports whether QMME hit its strict gradient
// tolerance.  The two can differ at extreme τ — slow QMME convergence is
// usually fine for the sandwich since A^{-1} B A^{-1} doesn't need tight
// gradient zero, only an accurate residual vector.
struct WaldResult {
    bool sandwichOk;     // SE is finite and positive
    bool qmmeTolMet;     // diagnostic only — st.converged from QMME
    double beta;
    double se;
    int iters;
};

WaldResult fitWaldOne(
    const Eigen::VectorXd &y,           // pheno-dense response (n)
    const Eigen::MatrixXd &X,           // pheno-dense covariates (n × p), no intercept
    const Eigen::VectorXd &G,           // pheno-dense genotype (n), no NaN
    double tau,
    double h,
    double tol,
    int maxIter
) {
    const Eigen::Index n = y.size();
    const int p = static_cast<int>(X.cols());
    const int dim = p + 2;  // intercept + p covars + G

    // Build the joint design (no intercept; QMME prepends one).
    Eigen::MatrixXd XG(n, p + 1);
    if (p > 0) XG.leftCols(p) = X;
    XG.col(p) = G;

    qmme::SqrSolver solver(XG, /*delta*/ 1e-6);
    solver.prepareBandwidth(h);

    Eigen::VectorXd resid;
    conquer::ConquerStatus st;
    Eigen::VectorXd theta = solver.solve(y, tau, &resid, tol, maxIter, /*restart*/ 50, &st);

    WaldResult out{};
    out.iters = st.gaussIter;
    out.qmmeTolMet = st.converged;
    out.beta = theta(dim - 1);   // last entry: γ̂
    out.se   = std::numeric_limits<double>::quiet_NaN();
    out.sandwichOk = false;
    // Sandwich is computed regardless of QMME's strict gradient-norm
    // convergence flag.  At extreme τ, QMME often stops short of its
    // gradient cutoff but the residual vector is accurate enough for the
    // M-estimation variance.

    // Build Z_full = [1 | X | G] in original space, n × dim.
    Eigen::MatrixXd Z(n, dim);
    Z.col(0).setOnes();
    if (p > 0) Z.middleCols(1, p) = X;
    Z.col(dim - 1) = G;

    // Per-i kernel weight K_h(-e_i) = 1/(h√(2π)) · exp(-(e_i)²/(2h²)).
    const double inv_sqrt2pi_h = 1.0 / (h * std::sqrt(2.0 * M_PI));
    const double inv_2h2 = 1.0 / (2.0 * h * h);
    Eigen::ArrayXd K(n);
    Eigen::ArrayXd R(n);
    for (Eigen::Index i = 0; i < n; ++i) {
        const double ei = resid(i);
        K(i) = inv_sqrt2pi_h * std::exp(-ei * ei * inv_2h2);
        R(i) = tau - math::pnorm(-ei / h);
    }

    // A = (1/n) Z^T diag(K) Z; B = (1/n) Z^T diag(R²) Z.
    const double inv_n = 1.0 / static_cast<double>(n);
    Eigen::MatrixXd ZtKZ(dim, dim);
    Eigen::MatrixXd ZtR2Z(dim, dim);
    {
        Eigen::MatrixXd ZK = Z.array().colwise() * K;
        Eigen::MatrixXd ZR = Z.array().colwise() * (R * R);
        ZtKZ.noalias() = ZK.transpose() * Z * inv_n;
        ZtR2Z.noalias() = ZR.transpose() * Z * inv_n;
    }

    // V = A^{-1} B A^{-1} / n.   A and B are symmetric, but A^{-1} B is NOT
    // (in general).  So solve in two stages: M = A^{-1} B, then V = A^{-1} M^T
    // (= A^{-1} B A^{-T} = A^{-1} B A^{-1} since A is symmetric).
    Eigen::LDLT<Eigen::MatrixXd> ldlt(ZtKZ);
    if (ldlt.info() != Eigen::Success) return out;     // sandwichOk stays false
    Eigen::MatrixXd M = ldlt.solve(ZtR2Z);              // dim × dim
    Eigen::MatrixXd V = ldlt.solve(M.transpose()) * inv_n;
    const double v_gg = V(dim - 1, dim - 1);

    if (v_gg > 0.0) {
        out.se = std::sqrt(v_gg);
        out.sandwichOk = true;
    }
    return out;
}

} // namespace

void runSPAsqrWald(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<double> &taus,
    const GenoSpec &geno,
    const std::string &predListFile,
    const std::string &outPrefix,
    double spasqrTol,
    double spasqrH,
    double spasqrHScale,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile,
    const std::string &removeFile,
    const std::string &phenoTransform
) {
    const int K = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());
    const bool useLoco = !predListFile.empty();
    // Wald defaults to h-scale=10 (vs score-mode's 3): per-marker QMME refits
    // with G in the design benefit from a smaller bandwidth so the kernel
    // weight K_h(-e) better resolves the score density f(0) and the
    // sandwich-derived SE matches the Gaussian asymptotic limit.
    const double effHScale = (spasqrHScale >= 0.0) ? spasqrHScale : 10.0;
    // Wald refits per (marker, τ) — keep iter cap modest. The ε_grad tolerance
    // tracks the user's --spasqr-tol directly (not floored to 1e-9 like score
    // mode, which fits once and reuses) so a single bad fit can't hang the run.
    const double qmmeTol = spasqrTol;
    const int maxIter = 5000;

    infoMsg("SPAsqr (wald): pheno-transform = %s, %s, ntaus = %d",
            phenoTransform.c_str(),
            useLoco ? "with LOCO offset" : "no LOCO",
            ntaus);

    // ── 1. Subject filtering (genotype ∩ keep/remove ∩ pheno) ───────────
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile, phenoNames);
    if (!covarFile.empty()) sd.loadCovar(covarFile, covarNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();

    const uint32_t nUnion = sd.nUsed();
    const Eigen::Index N = static_cast<Eigen::Index>(nUnion);

    Eigen::MatrixXd unionX = covarNames.empty()
        ? (sd.hasCovar() ? Eigen::MatrixXd(sd.covar()) : Eigen::MatrixXd(N, 0))
        : sd.getColumns(covarNames);
    const int nCov = static_cast<int>(unionX.cols());

    // ── 2. Per-phenotype: build pheno-dense Y/X (skip NA in Y) ─────────
    struct PhenoWork {
        std::vector<uint32_t> unionToLocal;  // size nUnion; UINT32_MAX = absent
        uint32_t nk = 0;
        Eigen::VectorXd Y;                   // nk; transformed
        Eigen::MatrixXd X;                   // nk × nCov
        // For LOCO: per-chr y_resp = Y - loco_chr (built in §4)
        // For no-LOCO: y_resp_global = Y (single vector)
        std::unordered_map<std::string, Eigen::VectorXd> yRespByChr;
        std::unordered_map<std::string, double> hByChr;
        Eigen::VectorXd yRespNoLoco;
        double hNoLoco = 0.0;
    };
    std::vector<PhenoWork> pw(K);

    for (int k = 0; k < K; ++k) {
        Eigen::VectorXd fullY = sd.getColumn(phenoNames[k]);
        pw[k].unionToLocal.assign(nUnion, UINT32_MAX);
        uint32_t loc = 0;
        for (uint32_t i = 0; i < nUnion; ++i)
            if (!std::isnan(fullY[i])) pw[k].unionToLocal[i] = loc++;
        pw[k].nk = loc;
        if (pw[k].nk == 0)
            throw std::runtime_error("SPAsqr (wald): phenotype '" + phenoNames[k]
                                     + "' has no non-missing subjects");
        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        pw[k].Y.resize(Nk);
        pw[k].X.resize(Nk, nCov);
        for (uint32_t i = 0; i < nUnion; ++i) {
            const uint32_t li = pw[k].unionToLocal[i];
            if (li == UINT32_MAX) continue;
            pw[k].Y[li] = fullY[i];
            if (nCov > 0) pw[k].X.row(li) = unionX.row(i);
        }
        applyPhenoTransform(pw[k].Y, phenoTransform);
    }

    // ── 3. Load LOCO (optional) ──────────────────────────────────────────
    std::unique_ptr<LocoData> loco;
    std::unordered_set<std::string> locoChroms;
    if (useLoco) {
        loco = std::make_unique<LocoData>(LocoData::load(predListFile, phenoNames,
                                                          sd.usedIIDs(), sd.famIIDs()));
        locoChroms = loco->availableChromosomes();
        infoMsg("SPAsqr (wald): LOCO available for %zu chromosomes", locoChroms.size());
    }

    // ── 4. Pre-compute per-pheno y_resp + bandwidth h ──────────────────
    for (int k = 0; k < K; ++k) {
        const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);
        if (useLoco) {
            for (const auto &chr : locoChroms) {
                const auto &locoVec = loco->scores.at(phenoNames[k]).at(chr);
                Eigen::VectorXd loco_dense(Nk);
                for (uint32_t i = 0; i < nUnion; ++i) {
                    const uint32_t li = pw[k].unionToLocal[i];
                    if (li != UINT32_MAX) loco_dense[li] = locoVec[i];
                }
                // LDAK / Regenie Step 1 output should include a PGS for every
                // subject in the analysis set. If any non-missing-Y subject is
                // absent, the parser leaves NaN at that position. Hard-fail
                // instead of silently corrupting the per-marker QMME refit.
                if (!loco_dense.allFinite()) {
                    const Eigen::Index nBad = Nk - loco_dense.array().isFinite().count();
                    throw std::runtime_error(
                        "SPAsqr-LOCO (wald): LOCO file for phenotype '" + phenoNames[k] +
                        "' chr " + chr + " is missing " + std::to_string(nBad) +
                        " subject(s) that have non-missing Y. The LOCO PGS file "
                        "must contain every subject in the --pheno analysis set. "
                        "Re-run LDAK / Regenie Step 1 on the same sample set, or "
                        "remove those subjects from --pheno.");
                }
                Eigen::VectorXd yResp = pw[k].Y - loco_dense;
                const double h = (spasqrH >= 0.0) ? spasqrH : iqrBandwidth(yResp, effHScale);
                pw[k].hByChr[chr] = h;
                pw[k].yRespByChr.emplace(chr, std::move(yResp));
            }
            infoMsg("SPAsqr (wald): [%s] precomputed y_resp for %zu chromosomes (n=%u)",
                    phenoNames[k].c_str(), locoChroms.size(), pw[k].nk);
        } else {
            pw[k].yRespNoLoco = pw[k].Y;
            pw[k].hNoLoco = (spasqrH >= 0.0) ? spasqrH : iqrBandwidth(pw[k].Y, effHScale);
            infoMsg("SPAsqr (wald): [%s] n=%u, h=%.6f",
                    phenoNames[k].c_str(), pw[k].nk, pw[k].hNoLoco);
        }
    }

    // ── 5. Build genotype meta + cursor ─────────────────────────────────
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(),
                                 /*chunk*/ 1024);
    const auto &markerList = genoData->markerInfo();   // already filtered by --extract/--exclude/--chr
    const size_t nMarkers = markerList.size();
    infoMsg("SPAsqr (wald): %zu markers selected after --extract/--exclude/--chr", nMarkers);
    if (nMarkers == 0) {
        warnMsg("SPAsqr (wald): nothing to test (empty marker list).");
        return;
    }

    auto cursor = genoData->makeCursor();
    cursor->beginSequentialBlock(markerList.front().genoIndex);

    // ── 6. Open per-phenotype output files ──────────────────────────────
    std::vector<std::ofstream> outFiles(K);
    for (int k = 0; k < K; ++k) {
        const std::string path = outPrefix + "." + phenoNames[k] + ".SPAsqr";
        outFiles[k].open(path);
        if (!outFiles[k]) throw std::runtime_error("Cannot open output file: " + path);
        outFiles[k] << "CHROM\tPOS\tID\tREF\tALT\tMISS_RATE\tALT_FREQ\tMAC\tHWE_P"
                    << "\tTAU\tBETA\tSE\tZ\tP\n";
        infoMsg("SPAsqr (wald): writing %s", path.c_str());
    }

    // ── 7. Main loop: per marker × per phenotype × per tau ─────────────
    Eigen::VectorXd unionG(nUnion);
    std::vector<uint32_t> indexForMissing;
    indexForMissing.reserve(64);

    uint32_t emitted = 0;
    uint32_t qcSkipped = 0;
    uint32_t locoMissing = 0;
    uint32_t fitFailures = 0;

    for (size_t mi = 0; mi < nMarkers; ++mi) {
        const auto &mInfo = markerList[mi];
        double altFreq, altCounts, missingRate, hweP, maf, mac;
        cursor->getGenotypes(mInfo.genoIndex, unionG, altFreq, altCounts,
                             missingRate, hweP, maf, mac, indexForMissing);

        const std::string &chrStr = mInfo.chrom;
        const uint32_t pos        = mInfo.pos;
        const std::string &idStr  = mInfo.id;
        const std::string &refStr = mInfo.ref;
        const std::string &altStr = mInfo.alt;

        // QC at the union level
        const bool qcFail = std::isnan(altFreq) || missingRate > missingCutoff
                            || maf < minMafCutoff || mac < minMacCutoff
                            || (hweCutoff > 0.0 && !std::isnan(hweP) && hweP < hweCutoff);
        if (qcFail) { ++qcSkipped; continue; }

        for (int k = 0; k < K; ++k) {
            const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);

            // Pick y_resp / h for this chr × phenotype
            const Eigen::VectorXd *yRespPtr = nullptr;
            double h_use = 0.0;
            if (useLoco) {
                auto it = pw[k].yRespByChr.find(chrStr);
                if (it == pw[k].yRespByChr.end()) {
                    // Marker on a chromosome with no LOCO entry — skip.
                    ++locoMissing;
                    continue;
                }
                yRespPtr = &it->second;
                h_use = pw[k].hByChr.at(chrStr);
            } else {
                yRespPtr = &pw[k].yRespNoLoco;
                h_use = pw[k].hNoLoco;
            }

            // Map union → pheno-dense; impute pheno-dense missing with 2*altFreq.
            Eigen::VectorXd Gk(Nk);
            for (uint32_t i = 0; i < nUnion; ++i) {
                const uint32_t li = pw[k].unionToLocal[i];
                if (li == UINT32_MAX) continue;
                const double v = unionG(i);
                Gk[li] = std::isnan(v) ? 2.0 * altFreq : v;
            }

            for (int t = 0; t < ntaus; ++t) {
                const double tau = taus[t];
                WaldResult wr;
                try {
                    wr = fitWaldOne(*yRespPtr, pw[k].X, Gk, tau, h_use, qmmeTol, maxIter);
                } catch (const std::exception &ex) {
                    warnMsg("[%s] %s tau=%.3f fit threw: %s",
                            phenoNames[k].c_str(), idStr.c_str(), tau, ex.what());
                    wr = WaldResult{false, false, std::nan(""), std::nan(""), -1};
                }
                if (!wr.sandwichOk) ++fitFailures;
                if (mi < 5 || (mi + 1) % 100 == 0) {
                    const char *status = wr.sandwichOk
                        ? (wr.qmmeTolMet ? "ok" : "ok(qmme-loose)")
                        : "FAIL";
                    infoMsg("[%s] %s tau=%.3f β=%.4e SE=%.4e iters=%d %s",
                            phenoNames[k].c_str(), idStr.c_str(), tau,
                            wr.beta, wr.se, wr.iters, status);
                }
                // Z/P are valid whenever the sandwich SE is — QMME tol gating
                // here was overly strict and previously suppressed legitimate
                // p-values for fits that converged optimisation-wise but did
                // not meet QMME's gradient cutoff.
                const double Z = (wr.sandwichOk && std::isfinite(wr.beta))
                                 ? wr.beta / wr.se : std::nan("");
                const double P = std::isfinite(Z) ? 2.0 * (1.0 - math::pnorm(std::fabs(Z)))
                                                  : std::nan("");

                std::ostringstream row;
                row << chrStr << '\t' << pos << '\t' << idStr << '\t'
                    << refStr << '\t' << altStr << '\t'
                    << std::scientific << std::setprecision(6)
                    << missingRate << '\t' << altFreq << '\t' << mac << '\t' << hweP << '\t'
                    << std::fixed << std::setprecision(3) << tau << '\t'
                    << std::scientific << std::setprecision(6)
                    << wr.beta << '\t' << wr.se << '\t' << Z << '\t' << P << '\n';
                outFiles[k] << row.str();
                ++emitted;
            }
        }

        if ((mi + 1) % 200 == 0)
            infoMsg("SPAsqr (wald): processed %zu/%zu markers (emitted %u rows)",
                    mi + 1, nMarkers, emitted);
    }

    for (auto &of : outFiles) of.close();
    infoMsg("SPAsqr (wald): done. Emitted=%u, QC-skipped=%u, LOCO-missing=%u, fitFailures=%u",
            emitted, qcSkipped, locoMissing, fitFailures);
}
