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
// Threading model: per-marker QR refit is driven through MethodBase /
// multiPhenoEngine (no-LOCO) or locoEngine (LOCO) — identical to the
// score-mode dispatch.  Output is plink2-style one-marker-per-line wide
// format with P_CCT + P_tau* + Z_tau* + BETA_tau* + SE_tau* columns,
// written via TextWriter honoring --compression {gz, zst}.

#include "spasqr/spasqr.hpp"
#include "spasqr/qmme.hpp"

#include "engine/loco.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

// ── Phenotype pre-transform (mirrors helper in spasqr.cpp). ──────────
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

// ── IQR-based bandwidth.  Wald defaults to scale=10. ────────────────
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

// ── Per-(marker, τ) Wald refit + sandwich variance. ─────────────────
//
// y, X, G are all pheno-dense.  X has no intercept (QMME prepends one).
// G must be NaN-free (engine imputes before invocation).
struct WaldResult {
    bool sandwichOk;     // SE is finite and positive
    double beta;
    double se;
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
    // (in general).  So solve in two stages: M = A^{-1} B, then V = A^{-1} M^T.
    Eigen::LDLT<Eigen::MatrixXd> ldlt(ZtKZ);
    if (ldlt.info() != Eigen::Success) return out;     // sandwichOk stays false
    Eigen::MatrixXd M = ldlt.solve(ZtR2Z);
    Eigen::MatrixXd V = ldlt.solve(M.transpose()) * inv_n;
    const double v_gg = V(dim - 1, dim - 1);

    if (v_gg > 0.0) {
        out.se = std::sqrt(v_gg);
        out.sandwichOk = true;
    }
    return out;
}

// ══════════════════════════════════════════════════════════════════════
// SPAsqrWaldMethod — MethodBase adapter driving per-marker QR refits.
//
// Holds one phenotype's null-model data (Y_resp, X, h, tol, maxIter) via
// shared_ptr so clones cost only two pointer copies.  Each clone owns
// its own getResultVec stack; the qmme::SqrSolver is constructed per
// marker (the design [X | G] changes per marker, so the solver cannot
// be hoisted).
// ══════════════════════════════════════════════════════════════════════

class SPAsqrWaldMethod : public MethodBase {
  public:
    struct Shared {
        Eigen::VectorXd Y_resp;        // pheno-dense response (Nk)
        Eigen::MatrixXd X;             // pheno-dense covariates (Nk × p), no intercept
        std::vector<double> taus;
        std::vector<std::string> tauLabels;
        double h        = 0.0;
        double qmmeTol  = 1e-7;
        int    maxIter  = 5000;
    };

    explicit SPAsqrWaldMethod(std::shared_ptr<const Shared> sh)
        : m_shared(std::move(sh)) {
    }

    std::unique_ptr<MethodBase> clone() const override {
        return std::make_unique<SPAsqrWaldMethod>(*this);
    }

    int resultSize() const override {
        return 1 + 4 * static_cast<int>(m_shared->taus.size());
    }

    std::string getHeaderColumns() const override {
        std::ostringstream oss;
        const auto &labels = m_shared->tauLabels;
        oss << "\tP_CCT";
        for (const auto &lab : labels) oss << "\tP_"    << lab;
        for (const auto &lab : labels) oss << "\tZ_"    << lab;
        for (const auto &lab : labels) oss << "\tBETA_" << lab;
        for (const auto &lab : labels) oss << "\tSE_"   << lab;
        return oss.str();
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double /*altFreq*/,
        int /*markerInChunkIdx*/,
        std::vector<double> &result
    ) override {
        const auto &sh = *m_shared;
        const int ntaus = static_cast<int>(sh.taus.size());

        std::vector<double> betas(ntaus, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> ses  (ntaus, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> zs   (ntaus, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> ps   (ntaus, std::numeric_limits<double>::quiet_NaN());

        // GVec is pheno-dense, NaN-imputed by the engine.
        const Eigen::VectorXd G = GVec;

        for (int t = 0; t < ntaus; ++t) {
            WaldResult wr;
            try {
                wr = fitWaldOne(sh.Y_resp, sh.X, G, sh.taus[t], sh.h,
                                sh.qmmeTol, sh.maxIter);
            } catch (const std::exception &) {
                wr = WaldResult{false, std::nan(""), std::nan("")};
            }
            if (wr.sandwichOk && std::isfinite(wr.beta) && wr.se > 0.0) {
                const double z = wr.beta / wr.se;
                betas[t] = wr.beta;
                ses[t]   = wr.se;
                zs[t]    = z;
                ps[t]    = 2.0 * (1.0 - math::pnorm(std::fabs(z)));
            }
        }

        const double pCCT = math::cauchyCombine(ps.data(), static_cast<int>(ps.size()));

        result.clear();
        result.reserve(resultSize());
        result.push_back(pCCT);
        for (double p : ps)    result.push_back(p);
        for (double z : zs)    result.push_back(z);
        for (double b : betas) result.push_back(b);
        for (double s : ses)   result.push_back(s);
    }

    int preferredBatchSize() const override {
        return 1;   // per-marker QR refit; parallelism is at the chunk level.
    }

  private:
    std::shared_ptr<const Shared> m_shared;
};

// ── Phenotype work struct: subject filtering + transform happens here. ─
struct PhenoWork {
    std::vector<uint32_t> unionToLocal;   // size nUnion; UINT32_MAX = absent
    uint32_t nk = 0;
    Eigen::VectorXd Y;                    // nk; transformed
    Eigen::MatrixXd X;                    // nk × nCov
    // No-LOCO path:
    Eigen::VectorXd yRespNoLoco;
    double hNoLoco = 0.0;
};

std::vector<std::string> makeTauLabels(const std::vector<double> &taus) {
    std::vector<std::string> labels;
    labels.reserve(taus.size());
    for (double tau : taus) {
        std::ostringstream oss;
        oss << "tau" << tau;
        labels.push_back(oss.str());
    }
    return labels;
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
    const std::string &phenoTransform,
    int nthreads,
    int nSnpPerChunk,
    const std::string &compression,
    int compressionLevel
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

    // ── 4. Pre-compute per-pheno y_resp + bandwidth h (no-LOCO path only) ─
    if (!useLoco) {
        for (int k = 0; k < K; ++k) {
            pw[k].yRespNoLoco = pw[k].Y;
            pw[k].hNoLoco = (spasqrH >= 0.0) ? spasqrH : iqrBandwidth(pw[k].Y, effHScale);
            infoMsg("SPAsqr (wald): [%s] n=%u, h=%.6f",
                    phenoNames[k].c_str(), pw[k].nk, pw[k].hNoLoco);
        }
    }

    // ── 5. Tau labels (shared by all phenotypes) ────────────────────────
    const std::vector<std::string> tauLabels = makeTauLabels(taus);

    // ── 6. Build genotype meta with auto-shrunk chunk size ──────────────
    // When the user has not overridden --chunk-size (default 8192), shrink
    // it so chunk count ≥ 4·nthreads.  This keeps the work-stealing thread
    // pool fed even when wald is invoked against a small --extract subset.
    int effChunk = nSnpPerChunk;
    {
        // Probe the marker count cheaply via a first GenoMeta build.  Pre-1.0
        // makeGenoData is the only marker enumeration path; build twice if we
        // need to shrink (cost is tiny vs per-marker QR refits).
        auto probe = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(),
                                  /*chunk*/ 1);
        const size_t nMarkers = probe->markerInfo().size();
        if (nMarkers == 0) {
            warnMsg("SPAsqr (wald): nothing to test (empty marker list).");
            return;
        }
        if (effChunk == 8192) {  // CLI default sentinel — auto-tune.
            const int threads = std::max(1, nthreads);
            const size_t target = static_cast<size_t>(threads) * 4;
            const size_t autoChunk =
                std::max<size_t>(1, (nMarkers + target - 1) / target);
            effChunk = static_cast<int>(std::min<size_t>(autoChunk, 8192));
        }
        infoMsg("SPAsqr (wald): %zu markers, %d threads, chunk-size = %d",
                nMarkers, std::max(1, nthreads), effChunk);
    }
    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(),
                                 effChunk);

    // ── 7. Dispatch to multiPhenoEngine (no-LOCO) or locoEngine (LOCO) ─
    if (!useLoco) {
        std::vector<PhenoTask> tasks(K);
        for (int k = 0; k < K; ++k) {
            auto shared = std::make_shared<SPAsqrWaldMethod::Shared>();
            shared->Y_resp    = std::move(pw[k].yRespNoLoco);
            shared->X         = pw[k].X;
            shared->taus      = taus;
            shared->tauLabels = tauLabels;
            shared->h         = pw[k].hNoLoco;
            shared->qmmeTol   = qmmeTol;
            shared->maxIter   = maxIter;

            tasks[k].phenoName    = phenoNames[k];
            tasks[k].method       = std::make_unique<SPAsqrWaldMethod>(std::move(shared));
            tasks[k].unionToLocal = pw[k].unionToLocal;
            tasks[k].nUsed        = pw[k].nk;
        }

        infoMsg("SPAsqr (wald): starting association (%d phenotypes, %d taus, %d threads)",
                K, ntaus, std::max(1, nthreads));
        multiPhenoEngine(
            *genoData, tasks, outPrefix, "SPAsqr",
            compression, compressionLevel, nthreads,
            missingCutoff, minMafCutoff, minMacCutoff, hweCutoff
        );
        return;
    }

    // LOCO path: rebuild K tasks per chromosome with that chromosome's
    // y_resp and bandwidth h.
    auto buildTasks = [&](const std::string &chr, std::vector<PhenoTask> &tasks) {
        tasks.resize(K);
        for (int k = 0; k < K; ++k) {
            const Eigen::Index Nk = static_cast<Eigen::Index>(pw[k].nk);

            // Map LOCO scores from union → pheno-dense.
            const auto &locoVec = loco->scores.at(phenoNames[k]).at(chr);
            Eigen::VectorXd loco_dense(Nk);
            for (uint32_t i = 0; i < nUnion; ++i) {
                const uint32_t li = pw[k].unionToLocal[i];
                if (li != UINT32_MAX) loco_dense[li] = locoVec[i];
            }
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

            auto shared = std::make_shared<SPAsqrWaldMethod::Shared>();
            shared->Y_resp    = std::move(yResp);
            shared->X         = pw[k].X;
            shared->taus      = taus;
            shared->tauLabels = tauLabels;
            shared->h         = h;
            shared->qmmeTol   = qmmeTol;
            shared->maxIter   = maxIter;

            tasks[k].phenoName    = phenoNames[k];
            tasks[k].method       = std::make_unique<SPAsqrWaldMethod>(std::move(shared));
            tasks[k].unionToLocal = pw[k].unionToLocal;
            tasks[k].nUsed        = pw[k].nk;

            infoMsg("SPAsqr (wald) chr%s [%s]: n=%u, h=%.6f",
                    chr.c_str(), phenoNames[k].c_str(), pw[k].nk, h);
        }
    };

    infoMsg("SPAsqr (wald): starting LOCO association (%d phenotypes, %d taus, %zu chroms, %d threads)",
            K, ntaus, locoChroms.size(), std::max(1, nthreads));
    locoEngine(
        *genoData, locoChroms, phenoNames, buildTasks,
        outPrefix, "SPAsqr",
        compression, compressionLevel, nthreads,
        missingCutoff, minMafCutoff, minMacCutoff, hweCutoff
    );
}
