// long_gee.cpp — marginal GEE null models for LoQus (see header).
#include "spasqr/long_gee.hpp"

#include "util/math_helper.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace longgee {

int Records::maxRecords() const {
    uint32_t m = 0;
    for (size_t i = 0; i + 1 < start.size(); ++i)
        m = std::max(m, start[i + 1] - start[i]);
    return static_cast<int>(m);
}

namespace {

// Outer iteration: at most this many solve / update-ρ rounds.  Observed in
// simulation (independence start, ψ-based ρ): 3–5.
constexpr int kMaxOuter = 50;
// Convergence is judged on dimensionless quantities only, so that rescaling
// the phenotype or a covariate cannot change the fit (scale equivariance):
// a change in fitted values relative to the residual scale (h for the
// quantile model, sd(y) for the linear one), and a change in ρ.
constexpr double kTolFit = 1e-8;
constexpr double kTolRho = 1e-8;
// Inner Newton (quantile model): ‖S‖∞ / N on the standardized design, where
// ψ is dimensionless, so the criterion is too.
constexpr int kMaxNewton = 100;
constexpr double kTolScore = 1e-10;
constexpr double kMinStep = 1e-4;

// Standardized design Z = [1 | (X − mean) / sd].  Working on Z keeps the
// Newton system well scaled and makes the score norm unit-free.  A constant
// covariate keeps sd = 1 (checkDesign rejects one before any fit).
//
// Z is never materialized.  At biobank scale it is a records × (p+1) matrix,
// and a copy per concurrent fit (plus a same-sized D·Z for the Hessian) was
// the dominant memory term: 60 MB per extra thread at 470k records × 21
// columns, ~2 GB per thread at 4M × 31.  Every product the fits need is
// formed from X instead, a block of kBlockRows rows at a time, centred and
// scaled on the fly exactly as Z's rows would be; X itself is shared
// read-only by all fits of a phenotype.
class StdDesign {
public:
    explicit StdDesign(const Records &rec) : rec_(rec), p_(rec.X.cols()), mx_(p_), sx_(p_) {
        const Eigen::Index n = rec.X.rows();
        for (Eigen::Index j = 0; j < p_; ++j) {
            const double m = rec.X.col(j).mean();
            const double v = (n > 1) ? (rec.X.col(j).array() - m).square().sum() / static_cast<double>(n - 1)
                                     : 0.0;
            mx_(j) = m;
            sx_(j) = (v > 0.0) ? std::sqrt(v) : 1.0;
        }
    }

    Eigen::Index cols() const { return p_ + 1; }

    // Z b
    Eigen::VectorXd times(const Eigen::VectorXd &b) const {
        Eigen::VectorXd out(rec_.X.rows());
        forBlocks([&](Eigen::Index a, const Eigen::MatrixXd &Zb) { out.segment(a, Zb.rows()) = Zb * b; });
        return out;
    }

    // Zᵀ v
    Eigen::VectorXd tTimes(const Eigen::VectorXd &v) const {
        Eigen::VectorXd out = Eigen::VectorXd::Zero(cols());
        forBlocks([&](Eigen::Index a, const Eigen::MatrixXd &Zb) {
            out.noalias() += Zb.transpose() * v.segment(a, Zb.rows());
        });
        return out;
    }

    // Σ_t w_t z_t z_tᵀ, or ZᵀZ when w is null
    Eigen::MatrixXd gram(const Eigen::VectorXd *w) const {
        Eigen::MatrixXd out = Eigen::MatrixXd::Zero(cols(), cols());
        forBlocks([&](Eigen::Index a, const Eigen::MatrixXd &Zb) {
            if (w) out.noalias() += Zb.transpose() * (Zb.array().colwise() * w->segment(a, Zb.rows()).array()).matrix();
            else   out.noalias() += Zb.transpose() * Zb;
        });
        return out;
    }

    // Σ_i c_i (Σ_{t∈i} z_t)(Σ_{t∈i} w_t z_t)ᵀ, with w ≡ 1 when null: the
    // between-record part of Σ_i Z_iᵀ R_i⁻¹ W_i Z_i for exchangeable R_i.
    Eigen::MatrixXd subjectCross(const Eigen::VectorXd &c, const Eigen::VectorXd *w) const {
        Eigen::MatrixXd out = Eigen::MatrixXd::Zero(cols(), cols());
        Eigen::MatrixXd Zi;
        Eigen::VectorXd sz(cols()), swz(cols());
        for (int i = 0; i < rec_.nSubj(); ++i) {
            const Eigen::Index a = rec_.start[i], m = rec_.start[i + 1] - a;
            rowsZ(a, m, Zi);
            sz = Zi.colwise().sum().transpose();
            swz = w ? (Zi.transpose() * w->segment(a, m)).eval() : sz;
            out.noalias() += c(i) * sz * swz.transpose();
        }
        return out;
    }

    // Original-space coefficients [b0, b_1..b_p] → standardized-space.
    Eigen::VectorXd toStd(const Eigen::VectorXd &bOrig) const {
        Eigen::VectorXd b(p_ + 1);
        b(0) = bOrig(0) + mx_.dot(bOrig.tail(p_));
        for (Eigen::Index j = 0; j < p_; ++j) b(j + 1) = bOrig(j + 1) * sx_(j);
        return b;
    }

private:
    static constexpr Eigen::Index kBlockRows = 4096;
    const Records &rec_;
    Eigen::Index p_;
    Eigen::RowVectorXd mx_, sx_;

    // Rows [a, a+m) of Z into out (m × (p+1)).
    void rowsZ(Eigen::Index a, Eigen::Index m, Eigen::MatrixXd &out) const {
        out.resize(m, p_ + 1);
        out.col(0).setOnes();
        out.rightCols(p_) = ((rec_.X.middleRows(a, m).rowwise() - mx_).array().rowwise() / sx_.array()).matrix();
    }

    template <typename F>
    void forBlocks(F &&f) const {
        const Eigen::Index n = rec_.X.rows();
        Eigen::MatrixXd Zb;
        for (Eigen::Index a = 0; a < n; a += kBlockRows) {
            rowsZ(a, std::min(kBlockRows, n - a), Zb);
            f(a, Zb);
        }
    }
};

Eigen::VectorXd perSubjectSum(const Records &rec, const Eigen::VectorXd &v) {
    const int nSubj = rec.nSubj();
    Eigen::VectorXd s(nSubj);
    for (int i = 0; i < nSubj; ++i)
        s(i) = v.segment(rec.start[i], rec.start[i + 1] - rec.start[i]).sum();
    return s;
}

Eigen::VectorXd recordCounts(const Records &rec) {
    const int nSubj = rec.nSubj();
    Eigen::VectorXd m(nSubj);
    for (int i = 0; i < nSubj; ++i) m(i) = static_cast<double>(rec.start[i + 1] - rec.start[i]);
    return m;
}

// Exchangeable moment estimator on the globally standardized residual:
//   ρ̂ = Σ_i Σ_{t≠s} z_it z_is / Σ_i m_i (m_i − 1).
// Returns NaN when it is undefined (no subject with ≥ 2 records, or a
// constant residual).
double momentRho(const Records &rec, const Eigen::VectorXd &r) {
    const Eigen::Index n = r.size();
    if (n < 2) return std::numeric_limits<double>::quiet_NaN();
    const double mean = r.mean();
    const double var = (r.array() - mean).square().sum() / static_cast<double>(n - 1);
    if (!(var > 0.0)) return std::numeric_limits<double>::quiet_NaN();
    const double inv = 1.0 / std::sqrt(var);
    double num = 0.0, den = 0.0;
    for (int i = 0; i < rec.nSubj(); ++i) {
        const uint32_t a = rec.start[i], b = rec.start[i + 1];
        const double m = static_cast<double>(b - a);
        if (m < 2.0) continue;
        double s = 0.0, q = 0.0;
        for (uint32_t t = a; t < b; ++t) {
            const double z = (r(t) - mean) * inv;
            s += z;
            q += z * z;
        }
        num += s * s - q;
        den += m * (m - 1.0);
    }
    if (den <= 0.0) return std::numeric_limits<double>::quiet_NaN();
    return num / den;
}

// ρ is admissible iff R_i is positive definite for every subject: ρ < 1 and
// 1 + (m_max − 1) ρ > 0.
bool admissible(double rho, int mMax) {
    return std::isfinite(rho) && rho < 1.0 && 1.0 + (mMax - 1) * rho > 0.0;
}

Eigen::VectorXd exchC(const Eigen::VectorXd &m, double rho) {
    return (rho / (1.0 + (m.array() - 1.0) * rho)).matrix();
}

Eigen::VectorXd exchWeights(const Eigen::VectorXd &subjSum, const Eigen::VectorXd &m, double rho) {
    return (subjSum.array() / (1.0 + (m.array() - 1.0) * rho)).matrix();
}

NullFit fallBack(NullFit fit, const Eigen::VectorXd &indepWeight, std::string why) {
    fit.weight = indepWeight;
    fit.rho = 0.0;
    fit.fellBack = true;
    fit.reason = std::move(why);
    return fit;
}

// ── Quantile score for fixed ρ ─────────────────────────────────────────
struct QScore {
    Eigen::VectorXd e, psi, S;
    double norm = 0.0;            // ‖S‖∞ / N
};

QScore quantileScore(const Records &rec, const StdDesign &d, const Eigen::VectorXd &b,
                     const Eigen::VectorXd &c, double rho, double tau, double h) {
    QScore q;
    q.e = rec.y - d.times(b);
    const Eigen::Index n = q.e.size();
    q.psi.resize(n);
    for (Eigen::Index t = 0; t < n; ++t) q.psi(t) = tau - math::pnorm(-q.e(t) / h);
    // Σ_i Z_iᵀ R_i⁻¹ ψ_i ∝ Zᵀ(ψ − u), u_t = c_i Σ_{s∈i} ψ_s for t ∈ i.
    const Eigen::VectorXd spsi = perSubjectSum(rec, q.psi);
    Eigen::VectorXd u(n);
    for (int i = 0; i < rec.nSubj(); ++i)
        u.segment(rec.start[i], rec.start[i + 1] - rec.start[i]).setConstant(c(i) * spsi(i));
    q.S = d.tTimes(q.psi - u) / (1.0 - rho);
    q.norm = q.S.cwiseAbs().maxCoeff() / static_cast<double>(n);
    return q;
}

// Newton on S(β) = Σ Z_iᵀ R_i⁻¹ ψ_i = 0 for fixed exchangeable ρ.  Returns
// false (with b untouched) when the Jacobian is singular or the iteration
// cannot reach the score tolerance.
bool newtonQuantile(const Records &rec, const StdDesign &d, Eigen::VectorXd &b,
                    const Eigen::VectorXd &m, double rho, double tau, double h) {
    const Eigen::VectorXd c = exchC(m, rho);
    Eigen::VectorXd cur = b;
    QScore q = quantileScore(rec, d, cur, c, rho, tau, h);
    for (int it = 0; it < kMaxNewton; ++it) {
        if (q.norm < kTolScore) break;
        // H = Σ Z_iᵀ R_i⁻¹ D_i Z_i,  D = φ(e/h)/h
        Eigen::VectorXd dd(q.e.size());
        for (Eigen::Index t = 0; t < dd.size(); ++t) dd(t) = math::dnorm(q.e(t) / h) / h;
        const Eigen::MatrixXd H = (d.gram(&dd) - d.subjectCross(c, &dd)) / (1.0 - rho);
        const Eigen::FullPivLU<Eigen::MatrixXd> lu(H);
        if (!lu.isInvertible()) return false;
        const Eigen::VectorXd step = lu.solve(q.S);
        if (!step.allFinite()) return false;

        double lam = 1.0;
        QScore next;
        for (;;) {
            next = quantileScore(rec, d, cur + lam * step, c, rho, tau, h);
            if (next.norm < q.norm || lam < kMinStep) break;
            lam *= 0.5;
        }
        const double moved = d.times(lam * step).cwiseAbs().maxCoeff() / h;
        cur += lam * step;
        q = std::move(next);
        if (moved < 1e-12) break;
    }
    if (!(q.norm < 1e-8)) return false;
    b = cur;
    return true;
}

} // namespace

void checkDesign(
    const Records &rec,
    const std::vector<std::string> &names,
    const std::string &context
) {
    const Eigen::Index n = rec.X.rows();
    const Eigen::Index p = rec.X.cols();
    // 1 − R² below this is taken as exact collinearity.  The Schur complement
    // is formed on a correlation matrix, so rounding sits near 1e-15 and a
    // genuinely (if strongly) collinear covariate with R² < 1 − 1e-10 passes.
    constexpr double kTolCollinear = 1e-10;
    for (Eigen::Index j = 0; j < p; ++j) {
        const auto x = rec.X.col(j);
        const double mean = x.mean();
        const double sd = std::sqrt((x.array() - mean).square().sum() / static_cast<double>(n));
        const double scale = x.cwiseAbs().maxCoeff();
        if (!(sd > 1e-12 * scale))
            throw std::runtime_error(context + ": covariate '" + names[j] +
                                     "' is constant over the analysed records (collinear with the"
                                     " intercept); the null-model design must have full rank");
    }
    // Correlation matrix of the covariates, from the standardized design's
    // Gram matrix (its intercept row and column are dropped; the centred
    // columns are orthogonal to the intercept).
    const Eigen::MatrixXd G = StdDesign(rec).gram(nullptr).bottomRightCorner(p, p);
    const Eigen::VectorXd gd = G.diagonal().cwiseSqrt();
    const Eigen::MatrixXd C = G.array() / (gd * gd.transpose()).array();
    std::vector<Eigen::Index> kept;
    for (Eigen::Index j = 0; j < p; ++j) {
        double resid = 1.0;                                // 1 − R² of column j on the kept ones
        if (!kept.empty()) {
            const Eigen::Index k = static_cast<Eigen::Index>(kept.size());
            Eigen::MatrixXd Ckk(k, k);
            Eigen::VectorXd ckj(k);
            for (Eigen::Index a = 0; a < k; ++a) {
                ckj(a) = C(kept[a], j);
                for (Eigen::Index b = 0; b < k; ++b) Ckk(a, b) = C(kept[a], kept[b]);
            }
            resid = 1.0 - ckj.dot(Ckk.ldlt().solve(ckj));
        }
        if (resid < kTolCollinear) {
            std::string before;
            for (Eigen::Index a : kept) before += (before.empty() ? "" : ", ") + names[a];
            throw std::runtime_error(context + ": covariate '" + names[j] +
                                     "' is a linear combination of the intercept and " + before +
                                     "; the null-model design must have full rank");
        }
        kept.push_back(j);
    }
}

bool hasResidualVariation(const Records &rec) {
    const StdDesign d(rec);
    // Center first so the check is invariant to an arbitrary outcome offset.
    const Eigen::VectorXd y = rec.y.array() - rec.y.mean();
    const Eigen::MatrixXd gram = d.gram(nullptr);
    const Eigen::VectorXd beta = gram.ldlt().solve(d.tTimes(y));
    const Eigen::VectorXd resid = y - d.times(beta);
    const double roundoff = 64.0 * std::numeric_limits<double>::epsilon() * y.norm();
    return !resid.allFinite() || resid.norm() > roundoff;
}

NullFit fitQuantile(
    const Records &rec,
    const qmme::SqrSolver &solver,
    double tau,
    double h,
    WorkingCorr corr,
    double qmmeTol
) {
    NullFit fit;
    Eigen::VectorXd resid;
    qmme::SolverStatus st;
    const Eigen::VectorXd bOrig = solver.solve(rec.y, tau, &resid, qmmeTol,
                                               /*maxIter*/ 50000, /*restartPeriod*/ 50, &st);
    fit.qmmeConverged = st.converged;

    const Eigen::Index n = resid.size();
    Eigen::VectorXd psi(n);
    for (Eigen::Index t = 0; t < n; ++t) psi(t) = tau - math::pnorm(-resid(t) / h);
    const Eigen::VectorXd indepWeight = perSubjectSum(rec, psi);
    fit.weight = indepWeight;
    if (corr == WorkingCorr::Independence) return fit;

    const int mMax = rec.maxRecords();
    if (mMax < 2)
        return fallBack(fit, indepWeight, "every subject has a single record");

    const StdDesign d(rec);
    const Eigen::VectorXd m = recordCounts(rec);
    Eigen::VectorXd b = d.toStd(bOrig);
    double rho = momentRho(rec, psi);

    for (int outer = 1; outer <= kMaxOuter; ++outer) {
        fit.outer = outer;
        if (!admissible(rho, mMax))
            return fallBack(fit, indepWeight,
                            "rho = " + std::to_string(rho) + " outside the admissible range");
        const Eigen::VectorXd bOld = b;
        if (!newtonQuantile(rec, d, b, m, rho, tau, h))
            return fallBack(fit, indepWeight, "Newton solve of the exchangeable GEE failed");
        const Eigen::VectorXd e = rec.y - d.times(b);
        for (Eigen::Index t = 0; t < n; ++t) psi(t) = tau - math::pnorm(-e(t) / h);
        const double rhoNew = momentRho(rec, psi);
        const double dFit = d.times(b - bOld).cwiseAbs().maxCoeff() / h;
        const bool done = outer > 1 && dFit < kTolFit && std::abs(rhoNew - rho) < kTolRho;
        if (done) {
            fit.rho = rho;
            fit.weight = exchWeights(perSubjectSum(rec, psi), m, rho);
            return fit;
        }
        rho = rhoNew;
    }
    return fallBack(fit, indepWeight,
                    "no convergence in " + std::to_string(kMaxOuter) + " outer rounds");
}

NullFit fitLinear(
    const Records &rec,
    WorkingCorr corr
) {
    NullFit fit;
    const StdDesign d(rec);
    const Eigen::Index n = rec.y.size();
    const Eigen::MatrixXd ZtZ = d.gram(nullptr);
    const Eigen::VectorXd Zty = d.tTimes(rec.y);
    const Eigen::LDLT<Eigen::MatrixXd> ols(ZtZ);
    if (ols.info() != Eigen::Success || !ols.isPositive())
        throw std::runtime_error("linear GEE: singular covariate design");
    Eigen::VectorXd b = ols.solve(Zty);
    Eigen::VectorXd e = rec.y - d.times(b);
    const Eigen::VectorXd indepWeight = perSubjectSum(rec, e);
    fit.weight = indepWeight;
    if (corr == WorkingCorr::Independence) return fit;

    const int mMax = rec.maxRecords();
    if (mMax < 2)
        return fallBack(fit, indepWeight, "every subject has a single record");

    const Eigen::VectorXd m = recordCounts(rec);
    const Eigen::VectorXd sy = perSubjectSum(rec, rec.y);
    const double ySd = std::sqrt((rec.y.array() - rec.y.mean()).square().sum() /
                                 static_cast<double>(std::max<Eigen::Index>(n - 1, 1)));
    double rho = momentRho(rec, e);

    for (int outer = 1; outer <= kMaxOuter; ++outer) {
        fit.outer = outer;
        if (!admissible(rho, mMax))
            return fallBack(fit, indepWeight,
                            "rho = " + std::to_string(rho) + " outside the admissible range");
        // GLS normal equations; the common 1/(1−ρ) factor cancels.
        const Eigen::VectorXd c = exchC(m, rho);
        const Eigen::MatrixXd A = ZtZ - d.subjectCross(c, nullptr);
        // Σ_i c_i SZ_i sy_i = Zᵀu with u_t = c_i sy_i for t ∈ i
        Eigen::VectorXd u(n);
        for (int i = 0; i < rec.nSubj(); ++i)
            u.segment(rec.start[i], rec.start[i + 1] - rec.start[i]).setConstant(c(i) * sy(i));
        const Eigen::VectorXd r = Zty - d.tTimes(u);
        const Eigen::LDLT<Eigen::MatrixXd> gls(A);
        if (gls.info() != Eigen::Success || !gls.isPositive())
            return fallBack(fit, indepWeight, "singular exchangeable GLS system");
        const Eigen::VectorXd bNew = gls.solve(r);
        const double dFit = d.times(bNew - b).cwiseAbs().maxCoeff() / ySd;
        b = bNew;
        e = rec.y - d.times(b);
        const double rhoNew = momentRho(rec, e);
        const bool done = outer > 1 && dFit < kTolFit && std::abs(rhoNew - rho) < kTolRho;
        if (done) {
            fit.rho = rho;
            fit.weight = exchWeights(perSubjectSum(rec, e), m, rho);
            return fit;
        }
        rho = rhoNew;
    }
    return fallBack(fit, indepWeight,
                    "no convergence in " + std::to_string(kMaxOuter) + " outer rounds");
}

} // namespace longgee
