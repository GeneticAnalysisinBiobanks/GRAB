// spasqr_wald.cpp — SPAsqr Wald mode (full-model effect size + SE)
//
// For every (marker, τ) pair, this mode refits the joint smoothed-QR model
// on Z = [1 | X | G] by damped Newton (qmme::SqrSolver::solveNewton), started
// from the null-model fit at that τ with the genotype coefficient at 0, and
// converged to ||grad||_inf <= --spasqr-tol.  It then computes the
// M-estimation sandwich
//      V = A^{-1} B A^{-1} / n,
// where
//      A = (1/n) Σ K_h(-e_i)   Z_i Z_i^T   (smooth-QR Hessian)
//      B = (1/n) Σ R_i^2       Z_i Z_i^T   (R_i = τ - Φ(-e_i/h))
// Effect size β̂_G is the last entry of θ̂ and SE = sqrt(V[γγ]).
//
// Quantile-heterogeneity test (Koenker–Bassett, smoothed, robust).  Stacking
// the ntaus estimating equations gives the joint sandwich
//      Cov(θ̂_k, θ̂_l) = A_k^{-1} B_kl A_l^{-1} / n,
//      B_kl = (1/n) Σ ψ_ik ψ_il Z_i Z_i^T,
// evaluated at each tau's own full-model (with G) residuals.  Its [γγ] block
// V is ntaus × ntaus with diagonal SE_k^2; with c_ik = (Z_i^T A_k^{-1} e_γ) ψ_ik,
// V = C^T C / n^2.  Under H0: β(τ_1) = … = β(τ_K) (any common value),
//      β̄ = 1'V^{-1}β̂ / 1'V^{-1}1,   Q = (β̂ − β̄1)' V^{-1} (β̂ − β̄1) ~ χ²_{K−1},
// reported as BETA_HOM / SE_HOM / HET_Q / LOG10P_HET (df = ntaus − 1).  Everything is
// accumulated in double precision from the fits themselves, never from the
// rounded output columns: at h = IQR/3 the cross-tau correlation matrix has
// condition number ~1e9, so the test needs the fits converged well below the
// printed precision (see --spasqr-tol).
//
// Threading model: per-marker QR refit is driven through MethodBase /
// multiPhenoEngine (no-LOCO) or locoEngine (LOCO) — identical to the
// score-mode dispatch.  Output is plink2-style one-marker-per-line wide
// format with LOG10P_CCT + P_tau* + LOG10P_tau* + Z_tau* + BETA_tau* +
// SE_tau* + SPA_STATUS_tau* columns, written via TextWriter honoring
// --compression {gz, zst}.

#include "spasqr/spasqr.hpp"
#include "spasqr/null_model.hpp"
#include "spasqr/qmme.hpp"

#include "engine/loco.hpp"
#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/spa.hpp"       // spa::normalTwoSidedLog, spa::Status

#include <Eigen/Dense>
#include <boost/math/special_functions/gamma.hpp>

#include <atomic>

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

// ── Per-marker Wald refit (all τ) + sandwich variance. ──────────────
//
// y, X, G are all pheno-dense.  X has no intercept (the solver prepends one).
// G must be NaN-free (engine imputes before invocation).
//
// Per marker: one SqrSolver([X | G]) (caches the standardised design; its
// majoriser Cholesky backs the Newton fallback step), and every τ starts from
// that τ's null-model fit with the genotype coefficient at 0.  From there
// Newton needs 2-3 iterations to 1e-8.
struct WaldResult {
    bool sandwichOk;     // SE is finite and positive
    double beta;
    double se;
};

struct HetResult {
    bool ok = false;
    double betaHom = std::numeric_limits<double>::quiet_NaN();
    double seHom   = std::numeric_limits<double>::quiet_NaN();
    double q       = std::numeric_limits<double>::quiet_NaN();
    double log10p  = std::numeric_limits<double>::quiet_NaN();
};

// Solver work counters for the run log (fits, iterations, non-converged).
std::atomic<long long> g_nFits{0}, g_nIter{0}, g_nNotConv{0};

// −log10 of the χ²_df upper tail, finite past the double underflow of P.
double chisqUpperNegLog10(double q, int df) {
    if (!(q >= 0.0) || df < 1) return std::numeric_limits<double>::quiet_NaN();
    const double a = 0.5 * df, x = 0.5 * q;
    const double p = boost::math::gamma_q(a, x);
    if (p > 1e-280) return -std::log10(p);
    // Asymptotic series  Q(a,x) ~ x^{a-1} e^{-x} / Γ(a) · Σ_k (a-1)…(a-k)/x^k.
    double s = 1.0, t = 1.0;
    for (int k = 1; k < 30; ++k) {
        t *= (a - k) / x;
        s += t;
        if (std::fabs(t) < 1e-17 * std::fabs(s)) break;
    }
    const double lnp = (a - 1.0) * std::log(x) - x - std::lgamma(a) + std::log(s);
    return -lnp / math::kLn10;
}

void fitWaldAllTaus(
    const Eigen::VectorXd &y,                // pheno-dense response (n)
    const Eigen::MatrixXd &X,                // pheno-dense covariates (n × p), no intercept
    const Eigen::VectorXd &G,                // pheno-dense genotype (n), no NaN
    const std::vector<double> &taus,
    const std::vector<Eigen::VectorXd> &nullTheta,   // per tau, [intercept, β_X] (size p+1)
    double h,
    double tol,
    int maxIter,
    std::vector<WaldResult> &out,
    HetResult &het
) {
    const Eigen::Index n = y.size();
    const int p = static_cast<int>(X.cols());
    const int dim = p + 2;                   // intercept + p covars + G
    const int ntaus = static_cast<int>(taus.size());

    Eigen::MatrixXd XG(n, p + 1);
    if (p > 0) XG.leftCols(p) = X;
    XG.col(p) = G;

    qmme::SqrSolver solver(XG, /*delta*/ 1e-6);
    solver.prepareBandwidth(h);

    Eigen::MatrixXd Z(n, dim);
    Z.col(0).setOnes();
    if (p > 0) Z.middleCols(1, p) = X;
    Z.col(dim - 1) = G;

    const double inv_sqrt2pi_h = 1.0 / (h * std::sqrt(2.0 * M_PI));
    const double inv_2h2 = 1.0 / (2.0 * h * h);
    const double inv_n = 1.0 / static_cast<double>(n);

    out.assign(ntaus, WaldResult{false, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::quiet_NaN()});
    het = HetResult{};

    Eigen::ArrayXd K(n), R(n);
    Eigen::MatrixXd C(n, ntaus);             // c_ik = (Z_i^T A_k^{-1} e_γ) ψ_ik
    bool allOk = true;
    Eigen::VectorXd eG = Eigen::VectorXd::Zero(dim);
    eG(dim - 1) = 1.0;

    for (int t = 0; t < ntaus; ++t) {
        const double tau = taus[t];
        // Null-model warm start at this tau, genotype coefficient 0.
        Eigen::VectorXd init(dim);
        init.head(p + 1) = nullTheta[t];
        init(dim - 1) = 0.0;

        Eigen::VectorXd resid;
        qmme::SolverStatus st;
        Eigen::VectorXd theta = solver.solveNewton(y, tau, &resid, tol, maxIter, &st, &init);
        g_nFits.fetch_add(1, std::memory_order_relaxed);
        g_nIter.fetch_add(st.iter, std::memory_order_relaxed);
        if (!st.converged) g_nNotConv.fetch_add(1, std::memory_order_relaxed);

        WaldResult row{false, theta(dim - 1), std::numeric_limits<double>::quiet_NaN()};

        for (Eigen::Index i = 0; i < n; ++i) {
            const double ei = resid(i);
            K(i) = inv_sqrt2pi_h * std::exp(-ei * ei * inv_2h2);
            R(i) = tau - math::pnorm(-ei / h);
        }
        Eigen::MatrixXd ZK = Z.array().colwise() * K;
        Eigen::MatrixXd ZR = Z.array().colwise() * (R * R);
        Eigen::MatrixXd ZtKZ  = ZK.transpose() * Z * inv_n;
        Eigen::MatrixXd ZtR2Z = ZR.transpose() * Z * inv_n;

        Eigen::LDLT<Eigen::MatrixXd> ldlt(ZtKZ);
        if (ldlt.info() == Eigen::Success) {
            Eigen::MatrixXd M = ldlt.solve(ZtR2Z);
            Eigen::MatrixXd V = ldlt.solve(M.transpose()) * inv_n;
            const double v_gg = V(dim - 1, dim - 1);
            if (v_gg > 0.0) {
                row.se = std::sqrt(v_gg);
                row.sandwichOk = true;
            }
            const Eigen::VectorXd aK = ldlt.solve(eG);             // A_k^{-1} e_γ
            C.col(t) = ((Z * aK).array() * R).matrix();
        }
        if (!row.sandwichOk || !std::isfinite(row.beta)) allOk = false;
        out[t] = row;
    }

    if (!allOk || ntaus < 2) return;
    // V = C^T C / n^2  (diagonal = SE^2).
    const Eigen::MatrixXd V = (C.transpose() * C) * (inv_n * inv_n);
    Eigen::LDLT<Eigen::MatrixXd> vl(V);
    if (vl.info() != Eigen::Success || !vl.isPositive()) return;
    Eigen::VectorXd b(ntaus);
    for (int t = 0; t < ntaus; ++t) b(t) = out[t].beta;
    const Eigen::VectorXd one = Eigen::VectorXd::Ones(ntaus);
    const Eigen::VectorXd Vi1 = vl.solve(one);
    const double s11 = one.dot(Vi1);
    if (!(s11 > 0.0) || !Vi1.allFinite()) return;
    const double bbar = Vi1.dot(b) / s11;
    const Eigen::VectorXd e = b - bbar * one;
    const double q = e.dot(vl.solve(e));
    if (!std::isfinite(q)) return;
    het.ok = true;
    het.betaHom = bbar;
    het.seHom = 1.0 / std::sqrt(s11);
    het.q = std::max(q, 0.0);
    het.log10p = chisqUpperNegLog10(het.q, ntaus - 1);
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
        double tol      = 1e-8;
        int    maxIter  = 200;
        std::vector<Eigen::VectorXd> nullTheta;   // per tau: null fit [intercept, β_X]
    };

    explicit SPAsqrWaldMethod(std::shared_ptr<const Shared> sh)
        : m_shared(std::move(sh)) {
    }

    std::unique_ptr<MethodBase> clone() const override {
        return std::make_unique<SPAsqrWaldMethod>(*this);
    }

    int resultSize() const override {
        return 1 + 5 * static_cast<int>(m_shared->taus.size()) + 4;
    }

    // LOG10P_CCT, then five per-tau groups: LOG10P, Z, BETA, SE, SPA_STATUS.
    // The grouping (all taus of one quantity, then all taus of the next) is
    // score mode's, and LOG10P / SPA_STATUS are placed exactly where score mode
    // places them, so the two modes remain readable against each other.  The
    // linear P group left in log10p_unify Stage 8 (decision D1): LOG10P_tau is
    // the sole p-value per quantile level and a consumer that needs the linear
    // p takes 10^(-LOG10P_tau).
    //
    // SPA_STATUS_tau is `1 NORMAL` on every tau that produced a test: this leg
    // is a plain Wald z against the normal, so no saddlepoint is ever attempted
    // and decision D4 assigns exactly that code to a test which does not use
    // one.  A tau whose sandwich variance is not usable (the LDLT failed, or
    // V[γγ] <= 0) has no statistic at all and takes `8 NA_NO_TEST` with
    // LOG10P, Z, BETA and SE all NA.  There is no fallback in between, because
    // a variance that does not exist leaves nothing to fall back to.
    std::string getHeaderColumns() const override {
        std::ostringstream oss;
        const auto &labels = m_shared->tauLabels;
        oss << "\tLOG10P_CCT";
        for (const auto &lab : labels) oss << "\tLOG10P_"     << lab;
        for (const auto &lab : labels) oss << "\tZ_"          << lab;
        for (const auto &lab : labels) oss << "\tBETA_"       << lab;
        for (const auto &lab : labels) oss << "\tSE_"         << lab;
        for (const auto &lab : labels) oss << "\tSPA_STATUS_" << lab;
        oss << "\tBETA_HOM\tSE_HOM\tHET_Q\tLOG10P_HET";
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
        std::vector<double> Ls   (ntaus, std::numeric_limits<double>::quiet_NaN());
        std::vector<double> sts  (ntaus,
            static_cast<double>(static_cast<uint8_t>(spa::Status::NaNoTest)));

        // GVec is pheno-dense, NaN-imputed by the engine.
        const Eigen::VectorXd G = GVec;

        std::vector<WaldResult> rows;
        HetResult het;
        try {
            fitWaldAllTaus(sh.Y_resp, sh.X, G, sh.taus, sh.nullTheta, sh.h,
                           sh.tol, sh.maxIter, rows, het);
        } catch (const std::exception &) {
            rows.assign(ntaus, WaldResult{false, std::nan(""), std::nan("")});
            het = HetResult{};
        }

        for (int t = 0; t < ntaus; ++t) {
            const auto &wr = rows[t];
            if (wr.sandwichOk && std::isfinite(wr.beta) && wr.se > 0.0) {
                const double z = wr.beta / wr.se;
                betas[t] = wr.beta;
                ses[t]   = wr.se;
                zs[t]    = z;
                // The magnitude first, the linear p from it.  This leg is a
                // plain normal tail, so 2*pnorm(|z|, upper) through Boost's
                // complement was already accurate wherever it was representable
                // — but it stops being representable at |z| = 38.6, and it is
                // what the Cauchy combination is taken over.  Carrying L and
                // deriving P from it (log10p_unify Stage 5, the pattern Stage 3
                // established for the tier's own P column) removes that ceiling
                // from LOG10P_CCT.  Since Stage 7 the magnitude is also a
                // column of its own, and since Stage 8 it is the only one: the
                // derived P_tau group is gone.
                Ls[t]    = -spa::normalTwoSidedLog(z) / math::kLn10;
                sts[t]   = static_cast<double>(
                    static_cast<uint8_t>(spa::Status::Normal));
            }
        }

        // Over the magnitudes, not the linear p-values: the Cauchy statistic
        // 1/(pi*p) overflows for p <= 1e-308 and the linear routine returned
        // P_CCT = 0 there (01_numerics §2.1).
        const double lCCT = math::cauchyCombineLog10(Ls.data(), static_cast<int>(Ls.size()));

        result.clear();
        result.reserve(resultSize());
        result.push_back(lCCT);
        for (double L : Ls)    result.push_back(L);
        for (double z : zs)    result.push_back(z);
        for (double b : betas) result.push_back(b);
        for (double s : ses)   result.push_back(s);
        for (double s : sts)   result.push_back(s);
        result.push_back(het.betaHom);
        result.push_back(het.seHom);
        result.push_back(het.q);
        result.push_back(het.log10p);
    }

    int preferredBatchSize() const override {
        return 1;   // per-marker QR refit; parallelism is at the chunk level.
    }

  private:
    std::shared_ptr<const Shared> m_shared;
};

// ── Phenotype work struct: subject filtering + transform happens here. ─

// Null-model fit [1 | X] at every tau, used as the warm start of the per-marker
// refits (genotype coefficient 0).  Same solver and tolerance as the refits.
std::vector<Eigen::VectorXd> fitNullThetas(
    const Eigen::VectorXd &y, const Eigen::MatrixXd &X, const std::vector<double> &taus,
    double h, double tol, int maxIter
) {
    std::vector<Eigen::VectorXd> th(taus.size());
    if (X.cols() == 0) {
        // No covariates: the warm start is the empirical tau-quantile of y.
        for (size_t t = 0; t < taus.size(); ++t) {
            std::vector<double> v(y.data(), y.data() + y.size());
            std::sort(v.begin(), v.end());
            const double idx = taus[t] * (v.size() - 1);
            const size_t lo = static_cast<size_t>(idx);
            const size_t hi = std::min(lo + 1, v.size() - 1);
            th[t] = Eigen::VectorXd::Constant(1, v[lo] + (idx - lo) * (v[hi] - v[lo]));
        }
        return th;
    }
    qmme::SqrSolver s(X, 1e-6);
    s.prepareBandwidth(h);
    for (size_t t = 0; t < taus.size(); ++t) {
        th[t] = s.solveNewton(y, taus[t], nullptr, tol, maxIter);
    }
    return th;
}

} // namespace

using namespace spasqr_null;

void runSPAsqrWald(const SPAsqrConfig &cfg) {
    const auto &phenoNames = cfg.phenoNames;
    const auto &taus       = cfg.taus;
    const int K = static_cast<int>(phenoNames.size());
    const int ntaus = static_cast<int>(taus.size());
    const bool useLoco = !cfg.predListFile.empty();
    // Wald defaults to h-scale=5 (vs score-mode's 3): per-marker QMME refits
    // with G in the design benefit from a smaller bandwidth so the kernel
    // weight K_h(-e) better resolves the score density f(0) and the
    // sandwich-derived SE matches the Gaussian asymptotic limit.
    const double effHScale = (cfg.spasqrHScale >= 0.0) ? cfg.spasqrHScale : 5.0;
    // --spasqr-tol is the Newton ||grad||_inf tolerance of the null-model and
    // per-marker fits.  The heterogeneity test needs the refits converged well
    // below the printed precision of BETA (default 1e-8: |Δβ| < 2e-5 SE).
    const double tol = cfg.spasqrTol;
    const int maxIter = 200;   // Newton: 2-3 iterations from the null warm start
    g_nFits = 0; g_nIter = 0; g_nNotConv = 0;

    infoMsg("SPAsqr (wald): pheno-transform = %s, %s, ntaus = %d, Newton tol = %g",
            cfg.phenoTransform.c_str(),
            useLoco ? "with LOCO offset" : "no LOCO",
            ntaus, tol);
    auto reportSolver = [&]() {
        const long long nf = g_nFits.load();
        infoMsg("SPAsqr (wald): Newton refits = %lld, mean iterations = %.2f, not converged = %lld",
                nf, nf ? double(g_nIter.load()) / nf : 0.0,
                g_nNotConv.load());
    };

    // ── 1. Subject filtering (genotype ∩ keep/remove ∩ pheno) ───────────
    auto famIIDs = parseGenoIIDs(cfg.geno);
    SubjectData sd(std::move(famIIDs));
    if (!cfg.phenoFile.empty()) sd.loadPhenoFile(cfg.phenoFile, phenoNames);
    if (!cfg.covarFile.empty()) sd.loadCovar(cfg.covarFile, cfg.covarNames);
    sd.setKeepRemove(cfg.keepFile, cfg.removeFile);
    sd.setGenoLabel(cfg.geno.flagLabel());
    sd.finalize();

    const uint32_t nUnion = sd.nUsed();
    const Eigen::Index N = static_cast<Eigen::Index>(nUnion);

    Eigen::MatrixXd unionX = cfg.covarNames.empty()
        ? (sd.hasCovar() ? Eigen::MatrixXd(sd.covar()) : Eigen::MatrixXd(N, 0))
        : sd.getColumns(cfg.covarNames);
    const int nCov = static_cast<int>(unionX.cols());

    // ── 2. Per-phenotype: build pheno-dense Y/X (skip NA in Y) ─────────
    std::vector<PhenoWork> pw = buildPhenoWorkspaces(
        sd, unionX, phenoNames, cfg.phenoTransform, "SPAsqr (wald)");

    // ── 3. Load LOCO (optional) ──────────────────────────────────────────
    std::unique_ptr<LocoData> loco;
    std::unordered_set<std::string> locoChroms;
    if (useLoco) {
        loco = std::make_unique<LocoData>(LocoData::load(cfg.predListFile, phenoNames,
                                                          sd.usedIIDs(), sd.famIIDs()));
        locoChroms = loco->availableChromosomes();
        infoMsg("SPAsqr (wald): LOCO available for %zu chromosomes", locoChroms.size());
    }

    // ── 4. Per-pheno bandwidth (no-LOCO path only; with LOCO it is derived
    //      per chromosome from the LOCO-adjusted response). ──────────────
    if (!useLoco) {
        for (int k = 0; k < K; ++k) {
            pw[k].h = (cfg.spasqrH >= 0.0) ? cfg.spasqrH
                                       : iqrBandwidth(pw[k].Y, effHScale, nCov);
            infoMsg("SPAsqr (wald): [%s] n=%u, h=%.6f",
                    phenoNames[k].c_str(), pw[k].nk, pw[k].h);
        }
    }

    // ── 5. Tau labels (shared by all phenotypes) ────────────────────────
    const std::vector<std::string> tauLabels = makeTauLabels(taus);

    // ── 6. Build genotype meta with auto-shrunk chunk size ──────────────
    // When the user has not overridden --chunk-ksnp (default 8 ksnp = 8192),
    // shrink it so chunk count ≥ 4·cfg.nthreads.  This keeps the work-stealing
    // thread pool fed even when wald is invoked against a small --extract subset.
    int effChunk = cfg.nSnpPerChunk;
    {
        // Probe the marker count cheaply via a first GenoMeta build.  Pre-1.0
        // makeGenoData is the only marker enumeration path; build twice if we
        // need to shrink (cost is tiny vs per-marker QR refits).
        auto probe = makeGenoData(cfg.geno, sd.usedMask(), sd.nFam(), sd.nUsed(),
                                  /*chunk*/ 1);
        const size_t nMarkers = probe->markerInfo().size();
        if (nMarkers == 0) {
            warnMsg("SPAsqr (wald): nothing to test (empty marker list).");
            return;
        }
        if (effChunk == 8192) {  // CLI default sentinel — auto-tune.
            const int threads = std::max(1, cfg.nthreads);
            const size_t target = static_cast<size_t>(threads) * 4;
            const size_t autoChunk =
                std::max<size_t>(1, (nMarkers + target - 1) / target);
            effChunk = static_cast<int>(std::min<size_t>(autoChunk, 8192));
        }
        infoMsg("SPAsqr (wald): %zu markers, %d threads, chunk-size = %d",
                nMarkers, std::max(1, cfg.nthreads), effChunk);
    }
    auto genoData = makeGenoData(cfg.geno, sd.usedMask(), sd.nFam(), sd.nUsed(),
                                 effChunk);

    // ── 7. Dispatch to multiPhenoEngine (no-LOCO) or locoEngine (LOCO) ─
    if (!useLoco) {
        std::vector<PhenoTask> tasks(K);
        for (int k = 0; k < K; ++k) {
            auto shared = std::make_shared<SPAsqrWaldMethod::Shared>();
            shared->Y_resp    = pw[k].Y;
            shared->X         = pw[k].X;
            shared->taus      = taus;
            shared->tauLabels = tauLabels;
            shared->h         = pw[k].h;
            shared->tol       = tol;
            shared->maxIter   = maxIter;
            shared->nullTheta = fitNullThetas(pw[k].Y, pw[k].X, taus, pw[k].h, tol, maxIter);

            tasks[k].phenoName    = phenoNames[k];
            tasks[k].method       = std::make_unique<SPAsqrWaldMethod>(std::move(shared));
            tasks[k].unionToLocal = pw[k].unionToLocal;
            tasks[k].nUsed        = pw[k].nk;
        }

        infoMsg("SPAsqr (wald): starting association (%d phenotypes, %d taus, %d threads)",
                K, ntaus, std::max(1, cfg.nthreads));
        multiPhenoEngine(
            *genoData, tasks, cfg.outPrefix, "SPAsqr",
            cfg.compression, cfg.compressionLevel, cfg.nthreads,
            cfg.missingCutoff, cfg.minMafCutoff, cfg.minMacCutoff, cfg.hweCutoff
        );
        reportSolver();
        return;
    }

    // LOCO path: rebuild K tasks per chromosome with that chromosome's
    // y_resp and bandwidth h.
    auto buildTasks = [&](const std::string &chr, std::vector<PhenoTask> &tasks) {
        tasks.resize(K);
        for (int k = 0; k < K; ++k) {
            Eigen::VectorXd yResp = pw[k].Y - locoDense(
                loco->scores.at(phenoNames[k]).at(chr), pw[k], phenoNames[k], chr,
                "SPAsqr-LOCO (wald)");
            const double h = (cfg.spasqrH >= 0.0) ? cfg.spasqrH : iqrBandwidth(yResp, effHScale, nCov);

            auto shared = std::make_shared<SPAsqrWaldMethod::Shared>();
            shared->Y_resp    = std::move(yResp);
            shared->X         = pw[k].X;
            shared->taus      = taus;
            shared->tauLabels = tauLabels;
            shared->h         = h;
            shared->tol       = tol;
            shared->maxIter   = maxIter;
            shared->nullTheta = fitNullThetas(shared->Y_resp, pw[k].X, taus, h, tol, maxIter);

            tasks[k].phenoName    = phenoNames[k];
            tasks[k].method       = std::make_unique<SPAsqrWaldMethod>(std::move(shared));
            tasks[k].unionToLocal = pw[k].unionToLocal;
            tasks[k].nUsed        = pw[k].nk;

            infoMsg("SPAsqr (wald) chr%s [%s]: n=%u, h=%.6f",
                    chr.c_str(), phenoNames[k].c_str(), pw[k].nk, h);
        }
    };

    infoMsg("SPAsqr (wald): starting LOCO association (%d phenotypes, %d taus, %zu chroms, %d threads)",
            K, ntaus, locoChroms.size(), std::max(1, cfg.nthreads));
    locoEngine(
        *genoData, locoChroms, phenoNames, buildTasks,
        cfg.outPrefix, "SPAsqr",
        cfg.compression, cfg.compressionLevel, cfg.nthreads,
        cfg.missingCutoff, cfg.minMafCutoff, cfg.minMacCutoff, cfg.hweCutoff
    );
    reportSolver();
}
