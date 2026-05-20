// sageld_fit.cpp — Long-format LMM fitting for SAGELD's pheno-input mode.
//
// EM-ML updates for two model templates:
//   • fitRandomSlopeML:     y ~ X + (Z | IID),  Z = [1, E]   (q = 2)
//   • fitRandomInterceptML: y ~ 1 + (1 | IID)                (q = 1)
//
// Per-subject computations use the Woodbury form
//      V_i^{-1} = σ^{-2} (I - Z_i M_i^{-1} Z_i'),
//      M_i      = σ² D^{-1} + Z_i'Z_i,
// so the q × q inverse M_i^{-1} is the only matrix solve needed per
// subject; n_i × n_i operations on Z_i appear only as small dot products.
//
// Stopping rule: relative change in (β, σ², D) below `tol` and a hard cap
// at `maxIter`.  A final E-step is performed after convergence so the
// returned residuals correspond exactly to the returned (β̂, σ̂², D̂).

#include "spagrm/sageld_fit.hpp"

#include "util/logging.hpp"
#include "util/math_helper.hpp"
#include "util/text_scanner.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace nsSAGELDFit {

namespace {

constexpr double kMinSigma2 = 1e-12;
constexpr double kMinDDiag = 1e-12;

bool isValidHeaderName(const std::string &s) {
    if (s.empty()) return false;
    for (char c : s) {
        if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.'))
            return false;
    }
    return true;
}

bool isMissingToken(std::string_view sv) {
    return sv.empty() || sv == "." || sv == "NA" || sv == "na" || sv == "NaN" || sv == "nan" || sv == "-";
}

} // namespace

// ════════════════════════════════════════════════════════════════════════
// parseLongPheno — read long-format file, group by IID, restrict to .fam
// ════════════════════════════════════════════════════════════════════════
LongPhenoData parseLongPheno(
    const std::string &filename,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::string &envName,
    const std::vector<std::string> &famIIDs,
    const std::unordered_set<std::string> &keptSubjects
) {
    if (phenoNames.empty()) throw std::runtime_error("SAGELD: --pheno-name required for direct-phenotype mode");
    if (envName.empty()) throw std::runtime_error("SAGELD: --sageld-x required for direct-phenotype mode");

    std::ifstream ifs(filename, std::ios::in | std::ios::binary);
    if (!ifs) throw std::runtime_error("Cannot open file: " + filename);

    std::string line;
    while (std::getline(ifs, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty()) break;
    }
    if (line.empty()) throw std::runtime_error(filename + ": empty file, header required");

    std::vector<std::string> headers;
    {
        text::TokenScanner ts(line);
        while (!ts.atEnd()) {
            auto sv = ts.nextView();
            if (sv.empty()) break;
            headers.emplace_back(sv);
        }
    }
    if (headers.size() < 2) throw std::runtime_error(filename + ": header needs at least 2 columns");

    int iidCol = 0;
    int dataStart = 1;
    if (headers.size() >= 3 &&
        (headers[0] == "#FID" || headers[0] == "FID") && headers[1] == "IID") {
        iidCol = 1;
        dataStart = 2;
        infoMsg("%s: detected FID+IID header; using IID column, ignoring FID", filename.c_str());
    }

    for (size_t c = static_cast<size_t>(dataStart); c < headers.size(); ++c) {
        if (!isValidHeaderName(headers[c]))
            throw std::runtime_error(filename + ": invalid header name '" + headers[c] +
                                     "' (must match [0-9A-Za-z_\\-.]+ )");
    }

    // Resolve column indices for env / pheno / covar (in header coordinates)
    std::unordered_map<std::string, int> headerIdx;
    for (int c = dataStart; c < static_cast<int>(headers.size()); ++c)
        headerIdx[headers[c]] = c;

    auto resolve = [&](const std::string &name, const char *role) -> int {
        auto it = headerIdx.find(name);
        if (it == headerIdx.end())
            throw std::runtime_error(filename + ": " + role + " column '" + name + "' not found in header");
        return it->second;
    };

    const int envColIdx = resolve(envName, "--sageld-x");
    std::vector<int> phenoColIdx;
    phenoColIdx.reserve(phenoNames.size());
    for (const auto &pn : phenoNames) phenoColIdx.push_back(resolve(pn, "--pheno-name"));
    std::vector<int> covarColIdx;
    covarColIdx.reserve(covarNames.size());
    for (const auto &cn : covarNames) covarColIdx.push_back(resolve(cn, "--covar-name"));

    // .fam IID → index lookup
    std::unordered_map<std::string, uint32_t> famMap;
    famMap.reserve(famIIDs.size() * 2);
    for (uint32_t k = 0; k < famIIDs.size(); ++k) famMap.emplace(famIIDs[k], k);

    // ── Parse rows ────────────────────────────────────────────────────────
    // Buffer per-row (IID, data) and filter to subjects present in .fam.
    struct RowBuf {
        uint32_t famIdx;
        std::vector<double> X;     // 1 + p_cov  (intercept + covariates)
        std::vector<double> Y;     // q
        double E;
    };
    std::vector<RowBuf> rows;
    rows.reserve(65536);

    const int nExpectedToks = static_cast<int>(headers.size());
    uint32_t lineNo = 1;
    uint32_t droppedNotInFam = 0;
    uint32_t droppedNA = 0;
    std::unordered_set<std::string> seenSubjNotInFam;

    while (std::getline(ifs, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        text::TokenScanner ts(line);
        std::string iid;
        int tokIdx = 0;
        for (int c = 0; c <= iidCol; ++c) {
            if (ts.atEnd()) break;
            if (c == iidCol)
                iid = ts.next();
            else
                ts.nextView();
            ++tokIdx;
        }
        while (tokIdx < dataStart) {
            if (ts.atEnd()) break;
            ts.nextView();
            ++tokIdx;
        }

        // Collect all data tokens into a fixed-position buffer
        std::vector<double> tokVals(nExpectedToks, std::numeric_limits<double>::quiet_NaN());
        for (int ci = dataStart; ci < nExpectedToks; ++ci) {
            if (ts.atEnd()) throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                                                     ": expected " + std::to_string(nExpectedToks) +
                                                     " columns, got " + std::to_string(ci));
            ts.skipWS();
            const char *tp = ts.pos();
            auto sv = ts.nextView();
            if (isMissingToken(sv)) {
                tokVals[ci] = std::numeric_limits<double>::quiet_NaN();
            } else {
                char *ep;
                double v = std::strtod(tp, &ep);
                if (ep == tp) throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                                                       ": non-numeric value in column " + std::to_string(ci + 1));
                tokVals[ci] = v;
            }
        }
        ts.skipWS();
        if (!ts.atEnd()) throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                                                  ": expected " + std::to_string(nExpectedToks) + " columns, got more");

        auto famIt = famMap.find(iid);
        if (famIt == famMap.end()) {
            ++droppedNotInFam;
            seenSubjNotInFam.insert(iid);
            continue;
        }
        if (!keptSubjects.empty() && !keptSubjects.count(iid)) {
            continue;
        }

        // Materialize this row's RowBuf — only if E, all phenos, all covars non-NaN
        RowBuf rb;
        rb.famIdx = famIt->second;
        rb.E = tokVals[envColIdx];
        rb.Y.resize(phenoColIdx.size());
        for (size_t k = 0; k < phenoColIdx.size(); ++k) rb.Y[k] = tokVals[phenoColIdx[k]];
        rb.X.resize(1 + covarColIdx.size());
        rb.X[0] = 1.0;
        for (size_t k = 0; k < covarColIdx.size(); ++k) rb.X[1 + k] = tokVals[covarColIdx[k]];

        // Drop row if any required entry is NaN
        bool hasNa = std::isnan(rb.E);
        for (size_t k = 0; k < rb.X.size() && !hasNa; ++k) if (std::isnan(rb.X[k])) hasNa = true;
        for (size_t k = 0; k < rb.Y.size() && !hasNa; ++k) if (std::isnan(rb.Y[k])) hasNa = true;
        if (hasNa) { ++droppedNA; continue; }

        rows.push_back(std::move(rb));
    }

    if (rows.empty()) throw std::runtime_error(filename + ": no usable rows after filtering");
    if (droppedNotInFam > 0)
        infoMsg("  Dropped %u rows (%zu distinct IIDs) absent from .fam",
                droppedNotInFam, seenSubjNotInFam.size());
    if (droppedNA > 0)
        infoMsg("  Dropped %u rows with NA in pheno / covar / env", droppedNA);

    // ── Group rows by .fam index, build CSR layout ────────────────────────
    std::sort(rows.begin(), rows.end(),
              [](const RowBuf &a, const RowBuf &b) { return a.famIdx < b.famIdx; });

    LongPhenoData data;
    data.phenoNames = phenoNames;
    data.covarNames = covarNames;
    data.envName = envName;

    const int N = static_cast<int>(rows.size());
    const int p = static_cast<int>(rows[0].X.size());
    const int q = static_cast<int>(rows[0].Y.size());
    data.X.resize(N, p);
    data.Y.resize(N, q);
    data.E.resize(N);

    data.subjStart.push_back(0);
    uint32_t curFam = rows[0].famIdx;
    data.famSubjIdx.push_back(curFam);
    data.uniqueIIDs.push_back(famIIDs[curFam]);

    for (int r = 0; r < N; ++r) {
        if (rows[r].famIdx != curFam) {
            data.subjStart.push_back(static_cast<uint32_t>(r));
            curFam = rows[r].famIdx;
            data.famSubjIdx.push_back(curFam);
            data.uniqueIIDs.push_back(famIIDs[curFam]);
        }
        data.E[r] = rows[r].E;
        for (int k = 0; k < p; ++k) data.X(r, k) = rows[r].X[k];
        for (int k = 0; k < q; ++k) data.Y(r, k) = rows[r].Y[k];
    }
    data.subjStart.push_back(static_cast<uint32_t>(N));

    infoMsg("Long-format: %d rows × %zu subjects (mean %.1f rows/subj), %d covariates + intercept, %zu phenotype(s)",
            N, data.uniqueIIDs.size(),
            static_cast<double>(N) / static_cast<double>(data.uniqueIIDs.size()),
            p - 1, phenoNames.size());

    return data;
}

// ════════════════════════════════════════════════════════════════════════
// fitRandomSlopeML — Y ~ X + (1 + E | IID) via REML profile likelihood
//                    + Nelder-Mead simplex search (3-D Cholesky parameters).
//
// Reparameterise τ = D / σ² via lower-triangular Cholesky factor
//      L(θ) = [[exp(θ_0), 0],
//              [θ_1,      exp(θ_2)]],
//      τ    = L Lᵀ,
// so θ ∈ ℝ³ is unconstrained and τ is always SPD.
//
// REML profile loss (σ² profiled out):
//      −2 ℓ_REML(τ) = (N − p) log(PSS(τ)) + Σᵢ log|I + τ Zᵢ'Zᵢ| + log|X' V⁻¹ X|.
// Each subject contributes a 2×2 inverse Mᵢ = τ⁻¹ + Zᵢ'Zᵢ; the cost per
// loss evaluation matches one EM iteration.  Nelder-Mead typically
// converges in 50–200 evaluations, vs 5000+ for the previous EM-ML loop.
//
// Boundary behaviour: when the optimised Cholesky factor satisfies
// ‖L‖_F < 1e-6, D is forced to zero and the fit collapses to OLS.
// This mirrors lme4's behaviour at a singular variance component, and
// avoids simplex jitter inside the flat basin near the boundary.
// ════════════════════════════════════════════════════════════════════════

namespace {

// Per-subject scratch reused across τ evaluations.  Allocated once.
struct SubjectCache {
    uint32_t r0;
    int n_i;
    Eigen::Matrix2d ZtZ;       // 2×2 with [[n_i, ΣE_i], [ΣE_i, ΣE_i²]]
    Eigen::MatrixXd ZtX;       // 2×p
    Eigen::Vector2d Zty;       // 2
    double yty;                // ‖y_i‖²
    Eigen::MatrixXd XtX_i;     // p×p
    Eigen::VectorXd Xty_i;     // p
};

std::vector<SubjectCache> buildSubjectCache(
    const LongPhenoData &data,
    const Eigen::VectorXd &y
) {
    const int nSubj = static_cast<int>(data.uniqueIIDs.size());
    const int p = static_cast<int>(data.X.cols());
    std::vector<SubjectCache> caches(nSubj);
    for (int i = 0; i < nSubj; ++i) {
        const uint32_t r0 = data.subjStart[i];
        const uint32_t r1 = data.subjStart[i + 1];
        const int n_i = static_cast<int>(r1 - r0);
        caches[i].r0 = r0;
        caches[i].n_i = n_i;
        const auto Ei = data.E.segment(r0, n_i);
        const auto X_i = data.X.middleRows(r0, n_i);
        const auto y_i = y.segment(r0, n_i);
        const double sumE = Ei.sum();
        const double sumE2 = Ei.squaredNorm();
        caches[i].ZtZ << static_cast<double>(n_i), sumE,
                         sumE, sumE2;
        caches[i].ZtX.resize(2, p);
        caches[i].ZtX.row(0) = X_i.colwise().sum();
        caches[i].ZtX.row(1) = (Ei.transpose() * X_i).eval();
        caches[i].Zty(0) = y_i.sum();
        caches[i].Zty(1) = Ei.dot(y_i);
        caches[i].yty = y_i.squaredNorm();
        caches[i].XtX_i = X_i.transpose() * X_i;
        caches[i].Xty_i = X_i.transpose() * y_i;
    }
    return caches;
}

// Decode θ ∈ ℝ³ → SPD τ via L L^T with log-diagonal Cholesky.
// θ_0, θ_2 clamped to ±25 to keep exp() in finite range.
Eigen::Matrix2d unpackTau(const std::array<double, 3> &theta) {
    const double t0 = std::clamp(theta[0], -25.0, 25.0);
    const double t1 = theta[1];
    const double t2 = std::clamp(theta[2], -25.0, 25.0);
    const double a = std::exp(t0);
    const double c = std::exp(t2);
    Eigen::Matrix2d L;
    L << a, 0.0,
         t1, c;
    return L * L.transpose();
}

// Profile-likelihood evaluator at fixed τ.  Returns the intermediate
// quantities needed to (a) form the REML loss, and (b) reconstruct β̂,
// σ̂², residuals after optimisation.
struct ProfileResult {
    Eigen::VectorXd beta;            // β̂(τ)
    double PSS;                      // (y-Xβ̂)'(I+ZτZ')^{-1}(y-Xβ̂)
    double sumLogDetIPlusTauZtZ;     // Σ_i log|I + τ Z_i'Z_i|
    double logDetA;                  // log|X' V^{-1} X| (un-scaled by σ²)
    bool valid;
};

ProfileResult profileEval(
    const Eigen::Matrix2d &tau,
    const std::vector<SubjectCache> &cache,
    int p
) {
    ProfileResult res;
    res.valid = false;

    bool invertible = false;
    Eigen::Matrix2d tauInv;
    double tauDet = 0.0;
    tau.computeInverseAndDetWithCheck(tauInv, tauDet, invertible, 1e-30);
    if (!invertible || tauDet <= 0.0) return res;
    const double logTauDet = std::log(tauDet);
    const double nSubj = static_cast<double>(cache.size());

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(p, p);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(p);
    double sumLogDetM = 0.0;

    // First pass: build A, b, accumulate Σ log|M_i|.
    for (const auto &c : cache) {
        const Eigen::Matrix2d M = tauInv + c.ZtZ;
        const double detM = M.determinant();
        if (detM <= 0.0) return res;
        const Eigen::Matrix2d Minv = M.inverse();
        sumLogDetM += std::log(detM);

        A.noalias() += c.XtX_i;
        A.noalias() -= c.ZtX.transpose() * Minv * c.ZtX;
        b.noalias() += c.Xty_i;
        b.noalias() -= c.ZtX.transpose() * (Minv * c.Zty);
    }

    Eigen::LDLT<Eigen::MatrixXd> ldltA(A);
    if (ldltA.info() != Eigen::Success) return res;
    res.beta = ldltA.solve(b);
    // log|A| = Σ log(D_ii); for SPD A, all D_ii > 0.
    const Eigen::VectorXd Adiag = ldltA.vectorD();
    double logDetA = 0.0;
    for (int j = 0; j < Adiag.size(); ++j) {
        const double d = Adiag(j);
        if (d <= 0.0) return res;
        logDetA += std::log(d);
    }
    res.logDetA = logDetA;

    // Second pass: PSS = Σ_i [‖r_i‖² − (Z_i' r_i)' M_i^{-1} (Z_i' r_i)]
    // with r_i = y_i − X_i β̂; expand via cached XtX_i, Xty_i, ZtX, Zty.
    double PSS = 0.0;
    for (const auto &c : cache) {
        const Eigen::Matrix2d M = tauInv + c.ZtZ;
        const Eigen::Matrix2d Minv = M.inverse();
        // ‖r_i‖² = yty − 2 Xty_i' β + β' XtX_i β
        const double rti = c.yty - 2.0 * c.Xty_i.dot(res.beta) + res.beta.dot(c.XtX_i * res.beta);
        // Z_i'r_i = Zty_i − ZtX_i β
        const Eigen::Vector2d Ztr = c.Zty - c.ZtX * res.beta;
        const double quad = Ztr.dot(Minv * Ztr);
        PSS += rti - quad;
    }
    res.PSS = PSS;
    res.sumLogDetIPlusTauZtZ = nSubj * logTauDet + sumLogDetM;
    res.valid = true;
    return res;
}

} // anonymous namespace

LMMFit fitRandomSlopeML(
    const LongPhenoData &data,
    const Eigen::VectorXd &y,
    int maxIter,
    double tol
) {
    const int N = static_cast<int>(y.size());
    const int p = static_cast<int>(data.X.cols());

    if (N != static_cast<int>(data.X.rows())) throw std::runtime_error("fitRandomSlopeML: y size mismatch with X");
    if (N != static_cast<int>(data.E.size())) throw std::runtime_error("fitRandomSlopeML: y size mismatch with E");
    if (N <= p) throw std::runtime_error("fitRandomSlopeML: N (= " + std::to_string(N) +
                                         ") must exceed p (= " + std::to_string(p) + ")");

    // ── Build per-subject cache once (one pass over all rows). ──────────
    const auto cache = buildSubjectCache(data, y);

    // ── OLS reference (initial point + boundary fallback). ──────────────
    const Eigen::VectorXd betaOLS =
        (data.X.transpose() * data.X).ldlt().solve(data.X.transpose() * y);
    const Eigen::VectorXd e0 = y - data.X * betaOLS;
    const double sigma2_OLS =
        std::max(kMinSigma2, e0.squaredNorm() / static_cast<double>(N - p));

    // θ₀ = (½log(0.1 σ²₀), 0, ½log(0.1 σ²₀)) → D ≈ 0.1 σ²₀ I in original scale.
    // Since θ parameterises τ = D/σ², use ½log(0.1) which puts τ ≈ 0.1·I.
    const double thetaInit = 0.5 * std::log(0.1);
    const std::vector<double> init = {thetaInit, 0.0, thetaInit};

    // ── Nelder-Mead on the REML profile loss. ───────────────────────────
    auto lossFunc = [&](const std::vector<double> &vt) -> double {
        const std::array<double, 3> theta = {vt[0], vt[1], vt[2]};
        const Eigen::Matrix2d tau = unpackTau(theta);
        const auto pr = profileEval(tau, cache, p);
        if (!pr.valid || pr.PSS <= 0.0)
            return 1e30; // infeasible / numerically degenerate
        return static_cast<double>(N - p) * std::log(pr.PSS)
             + pr.sumLogDetIPlusTauZtZ
             + pr.logDetA;
    };
    const auto opt = math::nelderMead(lossFunc, init, tol, maxIter);

    // ── Decode optimum + boundary handling. ─────────────────────────────
    const std::array<double, 3> thetaHat = {opt.par[0], opt.par[1], opt.par[2]};
    const double t0 = std::clamp(thetaHat[0], -25.0, 25.0);
    const double t2 = std::clamp(thetaHat[2], -25.0, 25.0);
    const double a = std::exp(t0);
    const double b_chol = thetaHat[1];
    const double c_chol = std::exp(t2);
    const double L_normsq = a * a + b_chol * b_chol + c_chol * c_chol;
    const bool collapseToOLS = (std::sqrt(L_normsq) < 1e-6);

    LMMFit fit;
    fit.beta.resize(p);
    fit.D.resize(2, 2);
    fit.residPerRow.resize(N);
    fit.iterations = opt.niter;
    fit.logLik = -0.5 * opt.value; // up to constants dropped from REML profile

    if (collapseToOLS) {
        fit.beta = betaOLS;
        fit.D.setZero();
        fit.sigma2 = sigma2_OLS;
        fit.residPerRow = e0;
        return fit;
    }

    const Eigen::Matrix2d tauHat = unpackTau(thetaHat);
    const auto pr = profileEval(tauHat, cache, p);
    if (!pr.valid || pr.PSS <= 0.0) {
        // Numerical fall-through: should be rare.
        fit.beta = betaOLS;
        fit.D.setZero();
        fit.sigma2 = sigma2_OLS;
        fit.residPerRow = e0;
        return fit;
    }

    fit.beta = pr.beta;
    fit.sigma2 = std::max(kMinSigma2, pr.PSS / static_cast<double>(N - p));
    fit.D = fit.sigma2 * tauHat;
    fit.D(0, 0) = std::max(fit.D(0, 0), kMinDDiag);
    fit.D(1, 1) = std::max(fit.D(1, 1), kMinDDiag);

    // ── Final pass: per-row residuals at (β̂, τ̂). ───────────────────────
    const Eigen::Matrix2d tauInv = tauHat.inverse();
    for (const auto &c : cache) {
        const Eigen::Matrix2d M = tauInv + c.ZtZ;
        const Eigen::Matrix2d Minv = M.inverse();
        const auto Ei = data.E.segment(c.r0, c.n_i);
        const auto X_i = data.X.middleRows(c.r0, c.n_i);
        const auto y_i = y.segment(c.r0, c.n_i);
        Eigen::VectorXd r_i = y_i - X_i * fit.beta;
        Eigen::Vector2d Ztr;
        Ztr(0) = r_i.sum();
        Ztr(1) = Ei.dot(r_i);
        const Eigen::Vector2d bHat = Minv * Ztr;
        Eigen::VectorXd eps_i = r_i;
        eps_i.array() -= bHat(0);
        eps_i.array() -= bHat(1) * Ei.array();
        fit.residPerRow.segment(c.r0, c.n_i) = eps_i;
    }
    return fit;
}

// ════════════════════════════════════════════════════════════════════════
// fitRandomInterceptML — y ~ 1 + (1 | IID) via 1-D REML profile likelihood
//                        + Brent minimiser.  q = 1 specialisation of the
//                        random-slope path; intercept-only X.
//
// Parameterise τ = D₁₁/σ² as τ = exp(2θ), θ ∈ ℝ.  For each subject i with
// n_i observations, M_i = 1/τ + n_i (scalar) and (using Woodbury):
//   |I + τ Zᵢ'Zᵢ| = 1 + n_i τ           (scalar)
//   Z_i'r_i      = Σ_j r_ij             (scalar)
// GLS μ̂ has closed form once τ is fixed.  REML profile loss:
//   −2 ℓ_REML(τ) = (N − 1) log(PSS(τ)) + Σ_i log(1 + n_i τ) + log A(τ)
//   A(τ) = Σ_i n_i / (τ M_i)            (scalar, p=1 → log A is scalar log)
// 1-D Brent on θ converges in ≲ 60 evaluations, sub-millisecond.
// ════════════════════════════════════════════════════════════════════════
LMMFit fitRandomInterceptML(
    const LongPhenoData &data,
    const Eigen::VectorXd &y,
    int maxIter,
    double tol
) {
    const int N = static_cast<int>(y.size());
    const int nSubj = static_cast<int>(data.uniqueIIDs.size());

    if (N <= 1) throw std::runtime_error("fitRandomInterceptML: need ≥ 2 rows");

    // ── Per-subject sufficient stats (one pass over y). ─────────────────
    std::vector<int> n_i(nSubj);
    Eigen::VectorXd sumy(nSubj);   // Σ_j y_ij
    Eigen::VectorXd ssy(nSubj);    // Σ_j y_ij²
    for (int i = 0; i < nSubj; ++i) {
        const uint32_t r0 = data.subjStart[i];
        const uint32_t r1 = data.subjStart[i + 1];
        n_i[i] = static_cast<int>(r1 - r0);
        const auto y_i = y.segment(r0, n_i[i]);
        sumy(i) = y_i.sum();
        ssy(i)  = y_i.squaredNorm();
    }

    // ── REML profile loss as a function of θ where τ = exp(2θ). ─────────
    // Returns +∞-ish sentinel on infeasible (PSS ≤ 0 or A ≤ 0) so Brent
    // recovers gracefully without throwing.
    auto evalLoss = [&](double theta, double *outMu = nullptr,
                        double *outPSS = nullptr) -> double {
        const double tClamped = std::clamp(theta, -25.0, 25.0);
        const double tau = std::exp(2.0 * tClamped);

        // GLS μ̂: A · μ̂ = b; A = Σ n_i / (τ M_i), b = Σ sumy_i / (τ M_i).
        // Note M_i = 1/τ + n_i ⇒ τ M_i = 1 + n_i τ. So we can factor:
        //   A = Σ n_i / (1 + n_i τ),   b = Σ sumy_i / (1 + n_i τ).
        // No τ blow-up: (1 + n_i τ) stays finite for any finite τ.
        double A = 0.0;
        double b_acc = 0.0;
        double sumLog = 0.0;   // Σ_i log(1 + n_i τ)
        for (int i = 0; i < nSubj; ++i) {
            const double w = 1.0 / (1.0 + static_cast<double>(n_i[i]) * tau);
            A      += static_cast<double>(n_i[i]) * w;
            b_acc  += sumy(i) * w;
            sumLog += std::log(1.0 + static_cast<double>(n_i[i]) * tau);
        }
        if (!(A > 0.0)) return 1e30;
        const double mu = b_acc / A;
        if (outMu) *outMu = mu;

        // PSS = Σ_i [ssy_i − 2 sumy_i μ + n_i μ² − (sumy_i − n_i μ)² / M_i],
        // where M_i = 1/τ + n_i.  Using w = 1/(1 + n_i τ) = (1/τ)/M_i, so
        // (sumy_i − n_i μ)² / M_i = τ · w · (sumy_i − n_i μ)².
        double PSS = 0.0;
        for (int i = 0; i < nSubj; ++i) {
            const double n = static_cast<double>(n_i[i]);
            const double w = 1.0 / (1.0 + n * tau);
            const double sr = sumy(i) - n * mu;          // Σ_j r_ij
            PSS += ssy(i) - 2.0 * sumy(i) * mu + n * mu * mu - tau * w * sr * sr;
        }
        if (!(PSS > 0.0)) return 1e30;
        if (outPSS) *outPSS = PSS;

        return static_cast<double>(N - 1) * std::log(PSS) + sumLog + std::log(A);
    };

    // ── Brent search on θ ∈ [-25, 25]. ──────────────────────────────────
    const double thetaHat = math::brentMin(evalLoss, -25.0, 25.0, tol, maxIter);
    double muHat = 0.0;
    double PSS   = 0.0;
    evalLoss(thetaHat, &muHat, &PSS);
    const double tauHat = std::exp(2.0 * std::clamp(thetaHat, -25.0, 25.0));
    const double sigma2Hat = std::max(kMinSigma2, PSS / static_cast<double>(N - 1));

    // ── Final pass: BLUP residuals at (μ̂, τ̂). ─────────────────────────
    Eigen::VectorXd resid(N);
    for (int i = 0; i < nSubj; ++i) {
        const uint32_t r0 = data.subjStart[i];
        const int n = n_i[i];
        const auto y_i = y.segment(r0, n);
        const double sr = sumy(i) - static_cast<double>(n) * muHat;
        // M_i = 1/τ + n_i  ⇒  b̂_i = sr / M_i = τ · w · sr where w = 1/(1+nτ).
        const double bHat = tauHat * sr / (1.0 + static_cast<double>(n) * tauHat);
        Eigen::VectorXd eps_i = y_i.array() - muHat - bHat;
        resid.segment(r0, n) = eps_i;
    }

    LMMFit fit;
    fit.beta.resize(1);
    fit.beta(0) = muHat;
    fit.D.resize(1, 1);
    fit.D(0, 0) = sigma2Hat * tauHat;       // τ² in original (D) scale
    fit.sigma2 = sigma2Hat;
    fit.residPerRow = std::move(resid);
    fit.iterations = 0;                      // Brent doesn't expose iter count
    fit.logLik = 0.0;
    return fit;
}

// ════════════════════════════════════════════════════════════════════════
// aggregatePerIID / aggregateWeightedPerIID
// ════════════════════════════════════════════════════════════════════════
Eigen::VectorXd aggregatePerIID(
    const LongPhenoData &data,
    const Eigen::VectorXd &residPerRow
) {
    const int nSubj = static_cast<int>(data.uniqueIIDs.size());
    Eigen::VectorXd out(nSubj);
    for (int i = 0; i < nSubj; ++i) {
        const uint32_t r0 = data.subjStart[i];
        const uint32_t r1 = data.subjStart[i + 1];
        out(i) = residPerRow.segment(r0, r1 - r0).sum();
    }
    return out;
}

Eigen::VectorXd aggregateWeightedPerIID(
    const LongPhenoData &data,
    const Eigen::VectorXd &residPerRow,
    const Eigen::VectorXd &weightPerRow
) {
    if (residPerRow.size() != weightPerRow.size())
        throw std::runtime_error("aggregateWeightedPerIID: resid/weight size mismatch");
    const int nSubj = static_cast<int>(data.uniqueIIDs.size());
    Eigen::VectorXd out(nSubj);
    for (int i = 0; i < nSubj; ++i) {
        const uint32_t r0 = data.subjStart[i];
        const uint32_t r1 = data.subjStart[i + 1];
        out(i) = residPerRow.segment(r0, r1 - r0)
                     .cwiseProduct(weightPerRow.segment(r0, r1 - r0))
                     .sum();
    }
    return out;
}

} // namespace nsSAGELDFit
