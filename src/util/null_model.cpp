// null_model.cpp — Unified null-model fitting engine implementation
//
// See util/null_model.hpp for the public contract and CLI mapping.
#include "util/null_model.hpp"

#include "io/subject_data.hpp"
#include "util/logging.hpp"
#include "util/regression.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <exception>
#include <fstream>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace nullmodel {

namespace {

std::string toLower(std::string s) {
    for (auto &c : s)
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return s;
}

// Trim ASCII whitespace from both ends of a token.
std::string trim(const std::string &s) {
    size_t a = 0;
    while (a < s.size() && std::isspace(static_cast<unsigned char>(s[a])))
        ++a;
    size_t b = s.size();
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1])))
        --b;
    return s.substr(a, b - a);
}

std::vector<std::string> splitCommaRaw(const std::string &s) {
    std::vector<std::string> out;
    std::istringstream iss(s);
    std::string tok;
    while (std::getline(iss, tok, ',')) {
        std::string t = trim(tok);
        if (!t.empty()) out.push_back(std::move(t));
    }
    return out;
}

// Build keep-mask in union rows: true where every named column is non-NaN.
// Returns the kept-row union indices.
std::vector<Eigen::Index> buildKeepIndices(
    const std::vector<Eigen::VectorXd> &cols,
    const Eigen::MatrixXd &covarUnion
) {
    if (cols.empty())
        throw std::runtime_error("buildKeepIndices: no columns supplied");
    const Eigen::Index n = cols.front().size();
    std::vector<Eigen::Index> keep;
    keep.reserve(static_cast<size_t>(n));
    for (Eigen::Index i = 0; i < n; ++i) {
        bool ok = true;
        for (const auto &c : cols) {
            if (std::isnan(c[i])) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;
        for (Eigen::Index j = 0; j < covarUnion.cols(); ++j) {
            if (std::isnan(covarUnion(i, j))) {
                ok = false;
                break;
            }
        }
        if (ok) keep.push_back(i);
    }
    return keep;
}

Eigen::VectorXd subsetVecByIdx(
    const Eigen::VectorXd &v, const std::vector<Eigen::Index> &idx
) {
    Eigen::VectorXd out(static_cast<Eigen::Index>(idx.size()));
    for (size_t i = 0; i < idx.size(); ++i)
        out[static_cast<Eigen::Index>(i)] = v[idx[i]];
    return out;
}

Eigen::MatrixXd subsetCovarByIdx(
    const Eigen::MatrixXd &M, const std::vector<Eigen::Index> &idx
) {
    Eigen::MatrixXd out(static_cast<Eigen::Index>(idx.size()), M.cols());
    for (size_t i = 0; i < idx.size(); ++i)
        out.row(static_cast<Eigen::Index>(i)) = M.row(idx[i]);
    return out;
}

// Convert a kept-row double phenotype vector to a contiguous-coded integer
// vector for cumulativeLogitFit.  Rejects values that are not non-negative
// integers in [0, 1e6].  iidLookup supplies a human-readable subject ID for
// the error message; idx maps kept-row position back to its union-row index.
Eigen::VectorXi coerceOrdinal(
    const Eigen::VectorXd &y,
    const std::vector<Eigen::Index> &idx,
    const std::vector<std::string> &unionIIDs,
    const std::string &columnName,
    const std::string &phenoName
) {
    Eigen::VectorXi out(y.size());
    for (Eigen::Index i = 0; i < y.size(); ++i) {
        double v = y[i];
        if (!std::isfinite(v) || v < 0.0 || v > 1e6 || v != std::floor(v)) {
            std::string iid = (idx[static_cast<size_t>(i)] >= 0 &&
                               static_cast<size_t>(idx[static_cast<size_t>(i)]) < unionIIDs.size())
                                  ? unionIIDs[static_cast<size_t>(idx[static_cast<size_t>(i)])]
                                  : std::string("?");
            std::ostringstream msg;
            msg << "Null-model fit for '" << phenoName << "': ordinal column '"
                << columnName << "' must contain non-negative integer values"
                << " in [0, 1e6]; got " << v << " at IID=" << iid;
            throw std::runtime_error(msg.str());
        }
        out[i] = static_cast<int>(v);
    }
    return out;
}

// Build a design matrix with an intercept column of ones and the covariates
// (which may have zero columns).
Eigen::MatrixXd buildDesignWithIntercept(const Eigen::MatrixXd &covarKept) {
    const Eigen::Index n = covarKept.rows();
    const Eigen::Index p = covarKept.cols();
    Eigen::MatrixXd X(n, 1 + p);
    X.col(0).setOnes();
    if (p > 0)
        X.rightCols(p) = covarKept;
    return X;
}

// Apply binary recode to `y` in place: smaller distinct value → 0, larger → 1.
// NaN entries are preserved.
void applyBinaryRecode(Eigen::VectorXd &y, double v0, double v1) {
    for (Eigen::Index i = 0; i < y.size(); ++i) {
        if (std::isnan(y[i])) continue;
        if (y[i] == v0) y[i] = 0.0;
        else if (y[i] == v1) y[i] = 1.0;
    }
}

// Apply ordinal shift to `y` in place: subtract `base` from each non-NaN value.
void applyOrdinalShift(Eigen::VectorXd &y, double base) {
    for (Eigen::Index i = 0; i < y.size(); ++i)
        if (!std::isnan(y[i])) y[i] -= base;
}

// Fit one PhenoSpec.  When `mHint == Auto`, the regression model is inferred
// per spec (via inferModelFromColumn / inferCoxSurvival).  Otherwise the spec
// is validated against the explicit `mHint` and (for Logistic / Ordinal) a
// recode is still applied if the data uses non-{0,1} or non-0-based encoding.
NullModelFit fitOne(
    const PhenoSpec &spec,
    RegressionModel mHint,
    const SubjectData &sd,
    const Eigen::MatrixXd &covarUnion,
    const EngineOptions &opts,
    const std::vector<std::string> &unionIIDs
) {
    NullModelFit fit;
    fit.name = spec.name;
    const Eigen::Index nUnion = static_cast<Eigen::Index>(sd.nUsed());
    fit.residuals.setConstant(nUnion,
                              std::numeric_limits<double>::quiet_NaN());

    // ── Cox path (spec contains TIME:EVENT) ──────────────────────────────
    if (isCoxSpec(spec)) {
        if (mHint != RegressionModel::Auto && mHint != RegressionModel::Cox)
            throw std::runtime_error(
                "--regression-model " + std::string(regressionModelName(mHint)) +
                " specified but '" + spec.name +
                "' uses TIME:EVENT syntax (which requires cox)");

        Eigen::VectorXd time = sd.getColumn(spec.timeColumn);
        Eigen::VectorXd event = sd.getColumn(spec.eventColumn);
        inferCoxSurvival(time, event, spec.timeColumn, spec.eventColumn, unionIIDs);

        if (mHint == RegressionModel::Auto)
            infoMsg("  Auto-detected regression model for '%s': cox", spec.name.c_str());

        auto idx = buildKeepIndices({time, event}, covarUnion);
        if (idx.empty())
            throw std::runtime_error(
                "Null-model fit for '" + spec.name +
                "': no complete cases after NaN removal");
        Eigen::VectorXd timeK = subsetVecByIdx(time, idx);
        Eigen::VectorXd eventK = subsetVecByIdx(event, idx);
        Eigen::MatrixXd X = subsetCovarByIdx(covarUnion, idx);
        Eigen::VectorXd w = Eigen::VectorXd::Ones(static_cast<Eigen::Index>(idx.size()));
        Eigen::VectorXd r = regression::coxResiduals(timeK, eventK, X, w,
                                                    opts.coxTol, opts.coxMaxIter);
        for (size_t i = 0; i < idx.size(); ++i)
            fit.residuals[idx[i]] = r[static_cast<Eigen::Index>(i)];
        fit.nUsedRows = static_cast<int>(idx.size());
        return fit;
    }

    // ── Non-Cox path: model is Linear / Logistic / Ordinal ───────────────
    if (mHint == RegressionModel::Cox)
        throw std::runtime_error(
            "--regression-model cox requires TIME:EVENT syntax;"
            " '" + spec.name + "' is a single-column spec");

    Eigen::VectorXd y = sd.getColumn(spec.yColumn);

    InferredCol info = inferModelFromColumn(y, spec.yColumn, unionIIDs);
    RegressionModel m;
    if (mHint == RegressionModel::Auto) {
        m = info.model;
        infoMsg("  Auto-detected regression model for '%s': %s",
                spec.name.c_str(), regressionModelName(m));
    } else {
        if (info.model != mHint)
            throw std::runtime_error(
                "--regression-model " + std::string(regressionModelName(mHint)) +
                " specified for '" + spec.yColumn +
                "' but inference returned " + regressionModelName(info.model));
        m = mHint;
    }

    if (info.needRecode) {
        if (m == RegressionModel::Logistic) {
            double v0 = info.sortedDistinct[0], v1 = info.sortedDistinct[1];
            applyBinaryRecode(y, v0, v1);
            infoMsg("  Recoded logistic column '%s': {%g, %g} -> {0, 1}",
                    spec.yColumn.c_str(), v0, v1);
        } else if (m == RegressionModel::Ordinal) {
            double base = info.sortedDistinct.front();
            applyOrdinalShift(y, base);
            infoMsg("  Recoded ordinal column '%s': shift by -%g",
                    spec.yColumn.c_str(), base);
        }
    }

    auto idx = buildKeepIndices({y}, covarUnion);
    if (idx.empty())
        throw std::runtime_error(
            "Null-model fit for '" + spec.name +
            "': no complete cases after NaN removal");
    Eigen::VectorXd yK = subsetVecByIdx(y, idx);
    Eigen::MatrixXd covarK = subsetCovarByIdx(covarUnion, idx);
    Eigen::VectorXd w = Eigen::VectorXd::Ones(static_cast<Eigen::Index>(idx.size()));

    Eigen::VectorXd r;
    if (m == RegressionModel::Linear) {
        Eigen::MatrixXd X = buildDesignWithIntercept(covarK);
        r = regression::linearResiduals(yK, X, w, opts.linearTol, opts.linearMaxIter);
    } else if (m == RegressionModel::Logistic) {
        Eigen::MatrixXd X = buildDesignWithIntercept(covarK);
        r = regression::logisticResiduals(yK, X, w, opts.logisticTol, opts.logisticMaxIter);
    } else { // Ordinal
        Eigen::VectorXi yi = coerceOrdinal(yK, idx, unionIIDs, spec.yColumn, spec.name);
        Eigen::MatrixXd X = buildDesignWithIntercept(covarK);
        // Forward the base seed verbatim.  cumulativeLogitFit constructs a
        // local std::mt19937 per call, so two ordinal phenotypes fit in the
        // same process draw from independent RNG instances even when sharing
        // the same seed value.  Using a portable single seed (rather than a
        // C++-specific hash combiner) makes the surrogate residuals exactly
        // reproducible from any R / Python re-implementation that seeds its
        // RNG with the same uint64_t.
        auto res = regression::cumulativeLogitFit(yi, X, opts.ordinalTol,
                                                  opts.ordinalMaxIter, opts.seed);
        r = std::move(res.residuals);
    }

    if (r.size() != static_cast<Eigen::Index>(idx.size()))
        throw std::runtime_error(
            "Null-model fit for '" + spec.name +
            "': internal error - regression returned " + std::to_string(r.size()) +
            " residuals for " + std::to_string(idx.size()) + " kept rows");

    for (size_t i = 0; i < idx.size(); ++i)
        fit.residuals[idx[i]] = r[static_cast<Eigen::Index>(i)];
    fit.nUsedRows = static_cast<int>(idx.size());
    return fit;
}

} // anonymous namespace

// ══════════════════════════════════════════════════════════════════════
// Public API
// ══════════════════════════════════════════════════════════════════════

RegressionModel parseRegressionModel(const std::string &s) {
    std::string k = toLower(trim(s));
    if (k.empty()) return RegressionModel::Auto;  // empty string ↔ default 'auto'
    if (k == "auto") return RegressionModel::Auto;
    if (k == "linear") return RegressionModel::Linear;
    if (k == "logistic") return RegressionModel::Logistic;
    if (k == "cox") return RegressionModel::Cox;
    if (k == "ordinal") return RegressionModel::Ordinal;
    // Earlier spellings rejected with a migration hint pointing to the
    // canonical regression-family names.
    if (k == "quantitative")
        throw std::runtime_error(
            "--regression-model 'quantitative' was renamed to 'linear'.");
    if (k == "binary")
        throw std::runtime_error(
            "--regression-model 'binary' was renamed to 'logistic'.");
    if (k == "time-to-event")
        throw std::runtime_error(
            "--regression-model 'time-to-event' was renamed to 'cox'.");
    throw std::runtime_error(
        "Unknown --regression-model '" + s +
        "'; expected one of: auto, linear, logistic, cox, ordinal");
}

const char *regressionModelName(RegressionModel m) {
    switch (m) {
    case RegressionModel::Auto:     return "auto";
    case RegressionModel::Linear:   return "linear";
    case RegressionModel::Logistic: return "logistic";
    case RegressionModel::Cox:      return "cox";
    case RegressionModel::Ordinal:  return "ordinal";
    }
    return "unknown";
}

InferredCol inferModelFromColumn(
    const Eigen::VectorXd &y,
    const std::string &colName,
    const std::vector<std::string> &unionIIDs,
    int ordinalMaxLevels
) {
    // 1. Collect non-NaN values
    std::vector<double> v;
    v.reserve(static_cast<size_t>(y.size()));
    for (Eigen::Index i = 0; i < y.size(); ++i)
        if (!std::isnan(y[i])) v.push_back(y[i]);
    if (v.empty())
        throw std::runtime_error(
            "Column '" + colName + "': all values are NaN");

    // 2. Sort + unique → distinct values
    std::sort(v.begin(), v.end());
    std::vector<double> distinct;
    distinct.reserve(v.size());
    distinct.push_back(v.front());
    for (size_t i = 1; i < v.size(); ++i)
        if (v[i] != distinct.back()) distinct.push_back(v[i]);

    if (distinct.size() == 1) {
        std::ostringstream m;
        m << "Column '" << colName << "': constant value " << distinct[0]
          << " (no variation; cannot infer trait)";
        throw std::runtime_error(m.str());
    }

    InferredCol info;
    info.sortedDistinct = distinct;

    // 3. Two distinct → Binary (lenient recode to {0, 1})
    if (distinct.size() == 2) {
        info.model = RegressionModel::Logistic;
        info.needRecode = !(distinct[0] == 0.0 && distinct[1] == 1.0);
        return info;
    }

    // 4. All-integer test (NaNs already filtered out)
    bool allInt = true;
    for (double x : distinct) {
        if (x != std::floor(x) || !std::isfinite(x)) {
            allInt = false;
            break;
        }
    }

    // 5. Ordinal: ≤ ordinalMaxLevels distinct integers, contiguous range
    if (allInt && static_cast<int>(distinct.size()) <= ordinalMaxLevels) {
        double base = distinct.front();
        double top = distinct.back();
        bool contiguous = (top - base + 1.0 == static_cast<double>(distinct.size()));
        if (contiguous) {
            info.model = RegressionModel::Ordinal;
            info.needRecode = (base != 0.0);
            return info;
        }
        // Non-contiguous integer codes (e.g. {0, 2, 5}) are not ordinal-
        // shaped, so treat them as Quantitative.  Explicit --regression-model
        // ordinal will still be rejected downstream because inference
        // returns Quantitative here.
        (void)unionIIDs; // currently unused; reserved for per-IID error hints
    }

    // 6. Otherwise → Quantitative (no recode)
    info.model = RegressionModel::Linear;
    info.needRecode = false;
    return info;
}

InferredSurv inferCoxSurvival(
    const Eigen::VectorXd &time,
    const Eigen::VectorXd &event,
    const std::string &timeName,
    const std::string &eventName,
    const std::vector<std::string> &unionIIDs
) {
    if (time.size() != event.size())
        throw std::runtime_error(
            "inferCoxSurvival: time and event columns differ in length");

    auto iidAt = [&](Eigen::Index i) -> std::string {
        if (i >= 0 && static_cast<size_t>(i) < unionIIDs.size())
            return unionIIDs[static_cast<size_t>(i)];
        return "?";
    };

    // Validate time: every non-NaN value must be > 0
    for (Eigen::Index i = 0; i < time.size(); ++i) {
        double t = time[i];
        if (std::isnan(t)) continue;
        if (!(t > 0.0)) {
            std::ostringstream m;
            m << "time column '" << timeName << "': non-positive value " << t
              << " at IID=" << iidAt(i) << " (time-to-event requires t > 0)";
            throw std::runtime_error(m.str());
        }
    }

    // Validate event: strict {0, 1}.  Collect distinct non-NaN values for a
    // helpful error on degenerate columns.
    bool sawZero = false, sawOne = false;
    for (Eigen::Index i = 0; i < event.size(); ++i) {
        double e = event[i];
        if (std::isnan(e)) continue;
        if (e == 0.0) { sawZero = true; continue; }
        if (e == 1.0) { sawOne = true; continue; }
        std::ostringstream m;
        m << "event column '" << eventName << "': value " << e
          << " at IID=" << iidAt(i)
          << " is not 0 or 1; event must be a 0/1 indicator";
        throw std::runtime_error(m.str());
    }
    if (!sawZero && !sawOne)
        throw std::runtime_error(
            "event column '" + eventName + "': all NaN");
    if (sawZero && !sawOne)
        throw std::runtime_error(
            "event column '" + eventName + "': no observed events (all censored)");
    // sawOne && !sawZero is the all-events case; valid for Cox.

    InferredSurv s;
    s.model = RegressionModel::Cox;
    return s;
}

PhenoSpec parsePhenoSpecAuto(const std::string &token) {
    std::string tok = trim(token);
    if (tok.empty())
        throw std::runtime_error("--pheno-name: empty phenotype token");
    size_t colonCount = std::count(tok.begin(), tok.end(), ':');
    if (colonCount > 1)
        throw std::runtime_error(
            "Phenotype '" + tok +
            "' contains more than one ':' (Cox spec must be TIME:EVENT)");
    PhenoSpec p;
    if (colonCount == 1) {
        size_t pos = tok.find(':');
        std::string timeCol = trim(tok.substr(0, pos));
        std::string eventCol = trim(tok.substr(pos + 1));
        if (timeCol.empty() || eventCol.empty())
            throw std::runtime_error(
                "Cox phenotype '" + tok +
                "': both TIME and EVENT column names must be non-empty");
        p.name = timeCol + "_" + eventCol;
        p.timeColumn = std::move(timeCol);
        p.eventColumn = std::move(eventCol);
    } else {
        p.name = tok;
        p.yColumn = tok;
    }
    return p;
}

std::vector<PhenoSpec> parsePhenoSpecList(
    RegressionModel m, const std::string &commaList
) {
    auto tokens = splitCommaRaw(commaList);
    if (tokens.empty())
        throw std::runtime_error("--pheno-name is empty");

    if (m == RegressionModel::Auto) {
        std::vector<PhenoSpec> out;
        out.reserve(tokens.size());
        for (const auto &tok : tokens)
            out.push_back(parsePhenoSpecAuto(tok));
        return out;
    }

    std::vector<PhenoSpec> out;
    out.reserve(tokens.size());
    for (const auto &tok : tokens) {
        size_t colonCount = std::count(tok.begin(), tok.end(), ':');
        if (m == RegressionModel::Cox) {
            if (colonCount != 1)
                throw std::runtime_error(
                    "cox phenotype '" + tok +
                    "' must be specified as TIME:EVENT (got " +
                    std::to_string(colonCount) + " ':' separator(s))");
            size_t pos = tok.find(':');
            std::string timeCol = trim(tok.substr(0, pos));
            std::string eventCol = trim(tok.substr(pos + 1));
            if (timeCol.empty() || eventCol.empty())
                throw std::runtime_error(
                    "Cox phenotype '" + tok +
                    "': both TIME and EVENT column names must be non-empty");
            // Canonical name uses '_' instead of ':' so the residual column
            // header satisfies the strict-format regex [0-9A-Za-z_\-.]+ when
            // the engine writes PREFIX.null.resid for round-trip reload.
            PhenoSpec p;
            p.name = timeCol + "_" + eventCol;
            p.timeColumn = std::move(timeCol);
            p.eventColumn = std::move(eventCol);
            out.push_back(std::move(p));
        } else {
            if (colonCount != 0)
                throw std::runtime_error(
                    "Regression model '" + std::string(regressionModelName(m)) +
                    "' does not use TIME:EVENT syntax; got '" + tok + "'");
            PhenoSpec p;
            p.name = tok;
            p.yColumn = tok;
            out.push_back(std::move(p));
        }
    }
    return out;
}

std::vector<std::string> columnsNeeded(const std::vector<PhenoSpec> &specs) {
    std::set<std::string> seen;
    std::vector<std::string> out;
    for (const auto &s : specs) {
        auto add = [&](const std::string &name) {
            if (name.empty()) return;
            if (seen.insert(name).second) out.push_back(name);
        };
        add(s.yColumn);
        add(s.timeColumn);
        add(s.eventColumn);
    }
    return out;
}

std::vector<NullModelFit> fitAll(
    const SubjectData &sd,
    const std::vector<PhenoSpec> &specs,
    RegressionModel m,
    const Eigen::MatrixXd &covarUnion,
    const EngineOptions &opts
) {
    if (specs.empty())
        throw std::runtime_error("Null-model engine: empty phenotype list");
    if (covarUnion.rows() != static_cast<Eigen::Index>(sd.nUsed()))
        throw std::runtime_error(
            "Null-model engine: covariate matrix row count (" +
            std::to_string(covarUnion.rows()) + ") != sd.nUsed() (" +
            std::to_string(sd.nUsed()) + ")");

    const size_t K = specs.size();
    const std::vector<std::string> unionIIDs = sd.usedIIDs();

    std::vector<NullModelFit> out(K);
    std::vector<std::exception_ptr> errs(K);
    std::atomic<size_t> next(0);
    int nthreads = std::max(1, opts.nthreads);
    nthreads = std::min(nthreads, static_cast<int>(K));

    auto worker = [&]() {
        for (;;) {
            size_t i = next.fetch_add(1, std::memory_order_relaxed);
            if (i >= K) return;
            try {
                out[i] = fitOne(specs[i], m, sd, covarUnion, opts, unionIIDs);
            } catch (...) {
                errs[i] = std::current_exception();
            }
        }
    };

    if (nthreads <= 1) {
        worker();
    } else {
        std::vector<std::thread> pool;
        pool.reserve(static_cast<size_t>(nthreads));
        for (int i = 0; i < nthreads; ++i) pool.emplace_back(worker);
        for (auto &th : pool) th.join();
    }

    for (size_t i = 0; i < K; ++i) {
        if (errs[i]) std::rethrow_exception(errs[i]);
    }
    return out;
}

void writeResidualsFile(
    const std::string &outPath,
    const SubjectData &sd,
    const std::vector<NullModelFit> &fits
) {
    if (fits.empty())
        throw std::runtime_error("writeResidualsFile: no fits to write");
    const auto iids = sd.usedIIDs();
    const Eigen::Index N = static_cast<Eigen::Index>(iids.size());
    for (const auto &f : fits) {
        if (f.residuals.size() != N)
            throw std::runtime_error(
                "writeResidualsFile: residual length (" +
                std::to_string(f.residuals.size()) + ") != nUsed (" +
                std::to_string(N) + ") for '" + f.name + "'");
    }

    std::ofstream ofs(outPath);
    if (!ofs)
        throw std::runtime_error("writeResidualsFile: cannot open '" + outPath + "' for writing");

    ofs << "IID";
    for (const auto &f : fits) ofs << '\t' << f.name;
    ofs << '\n';

    char buf[32];
    for (Eigen::Index i = 0; i < N; ++i) {
        ofs << iids[static_cast<size_t>(i)];
        for (const auto &f : fits) {
            double v = f.residuals[i];
            ofs << '\t';
            if (std::isnan(v)) {
                ofs << "NA";
            } else {
                // %.17g preserves a double exactly across the round-trip
                // (DBL_DECIMAL_DIG = 17).  Without this, --save-resid +
                // --resid-name reload introduces ~1e-10 residual error
                // which propagates into ~1e-4 SPA p-value differences for
                // logistic/ordinal traits.
                std::snprintf(buf, sizeof(buf), "%.17g", v);
                ofs << buf;
            }
        }
        ofs << '\n';
    }
    if (!ofs)
        throw std::runtime_error("writeResidualsFile: write failed on '" + outPath + "'");
    infoMsg("Wrote fitted residuals to %s (%lld subjects, %zu phenotype(s))",
            outPath.c_str(), static_cast<long long>(N), fits.size());
}

} // namespace nullmodel
