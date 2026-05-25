// null_model.hpp — Unified null-model fitting engine
//
// Provides the upstream null-model fitting path consumed by SPACox / SPAGRM /
// SPAmix when the user supplies phenotype columns instead of pre-computed
// residuals.  Each PhenoSpec is fitted independently by the appropriate
// regression function in src/util/regression.hpp; multiple phenotypes are
// fitted in parallel via a std::thread / std::atomic pool that mirrors the
// marker-engine pattern in src/engine/marker.cpp.
//
// Output is a union-length residual vector per phenotype (NaN-padded at rows
// the per-pheno regression dropped because of missing inputs).  Downstream
// machinery (SubjectData::buildPerColumnMasks) handles the NaN-padded rows.
//
// CLI mapping:
//   --regression-model    auto | linear | logistic | cox | ordinal
//   --pheno-name          non-Cox:  Y1,Y2,...
//                         Cox:      TIME1:EVENT1,TIME2:EVENT2,...
//   --save-resid          writes PREFIX.null.resid in the strict format
//                         consumed by SubjectData::loadResidOne, so a
//                         subsequent invocation can reuse the residuals via
//                         --resid-name.
#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

class SubjectData;

namespace nullmodel {

// User-facing regression-model categories accepted by --regression-model.
// `Auto` triggers per-spec inference from the phenotype data (column
// distinct-value pattern plus colon-syntax for survival).  The names align
// with the regression family used downstream: linear (continuous response,
// least-squares), logistic (0/1 indicator), cox (proportional-hazards
// survival, TIME:EVENT spec), ordinal (cumulative-logit on integer codes).
enum class RegressionModel { Auto, Linear, Logistic, Cox, Ordinal };

// Case-insensitive parse: auto | linear | logistic | cox | ordinal.
// The earlier "quantitative / binary / time-to-event" spellings are rejected
// with a migration hint pointing to the canonical regression-model names.
RegressionModel parseRegressionModel(const std::string &s);

const char *regressionModelName(RegressionModel m);

// Auto-inference results for a single-column phenotype.  `needRecode == true`
// signals that the caller must map the original values to the regression's
// expected coding (logistic: smaller value → 0, larger → 1; ordinal: subtract
// min so codes start at 0).  `sortedDistinct` carries the original distinct
// values for both recode and log output.
struct InferredCol {
    RegressionModel model;
    bool needRecode;
    std::vector<double> sortedDistinct;
};

// Survival-pair validation result.  Always returns `Cox`; the event column
// is required to be strictly {0, 1}, so no recoding is performed.
struct InferredSurv {
    RegressionModel model; // always Cox
};

// Infer the regression model of a single-column phenotype from its values.
// Rules:
//   all NaN              → throw
//   constant             → throw
//   2 distinct           → Logistic (recode {v0,v1} → {0,1})
//   integer, ≤ 10 distinct, contiguous → Ordinal (shift to 0..J-1)
//   anything else        → Linear (no recode)
// Note: integer columns whose distinct values are non-contiguous (e.g.
// {0, 2, 5}) are treated as Linear.
InferredCol inferModelFromColumn(
    const Eigen::VectorXd &y,
    const std::string &colName,
    const std::vector<std::string> &unionIIDs,
    int ordinalMaxLevels = 10);

// Validate a survival pair (time, event) and confirm model = Cox.
// Throws on time ≤ 0, all-NaN event, single-value event ≠ {1}, or any event
// value outside {0, 1}.
InferredSurv inferCoxSurvival(
    const Eigen::VectorXd &time,
    const Eigen::VectorXd &event,
    const std::string &timeName,
    const std::string &eventName,
    const std::vector<std::string> &unionIIDs);

struct PhenoSpec {
    std::string name;        // canonical: "Y1" (non-Cox) or "T1_E1" (Cox)
    std::string yColumn;     // non-Cox: phenotype column header
    std::string timeColumn;  // Cox: time column header
    std::string eventColumn; // Cox: event column header
};

// True when the spec represents a survival (Cox) phenotype.
inline bool isCoxSpec(const PhenoSpec &s) {
    return !s.timeColumn.empty();
}

// Parse the raw --pheno-name string (commas separate phenotypes).
//   non-Cox models (linear / logistic / ordinal): each token is a single
//     column name and must not contain ':'.
//   Cox model: each token has the form "TIME:EVENT" (exactly one ':').
//   Auto:   per-token inference: ':' → Cox spec, otherwise → non-Cox spec.
std::vector<PhenoSpec> parsePhenoSpecList(
    RegressionModel m, const std::string &commaList);

// Convenience: parse a single token in Auto mode and return its PhenoSpec.
// Used by callers that already split on commas and want per-token inference.
PhenoSpec parsePhenoSpecAuto(const std::string &token);

// Distinct pheno-file column names that loadPhenoFile should request.
std::vector<std::string> columnsNeeded(const std::vector<PhenoSpec> &specs);

struct NullModelFit {
    std::string name;          // mirrors PhenoSpec::name
    Eigen::VectorXd residuals; // length = sd.nUsed(); NaN where regression dropped row
    int nUsedRows = 0;         // rows actually fitted (after NaN drop)
};

struct EngineOptions {
    int nthreads = 1;
    double linearTol = 1e-9;
    int linearMaxIter = 50;
    double logisticTol = 1e-9;
    int logisticMaxIter = 50;
    double coxTol = 1e-10;
    int coxMaxIter = 40;
    double ordinalTol = 1e-7;
    int ordinalMaxIter = 50;
    // RNG seed forwarded to cumulativeLogitFit for ordinal surrogate-residual
    // sampling.  seed == 0 → std::random_device (non-reproducible across runs).
    uint64_t seed = 0;
};

// Fit each PhenoSpec; up to opts.nthreads phenotypes in parallel.
// Preconditions:
//   sd.finalize() has been called; loadPhenoFile has loaded every column
//     listed in columnsNeeded(specs).
//   covarUnion has shape (sd.nUsed(), p) where p >= 0, and does NOT include
//     an intercept column.  Linear/logistic/ordinal prepend an intercept
//     internally; Cox does not.
// Throws std::runtime_error (with phenotype name) on per-pheno failure.
std::vector<NullModelFit> fitAll(
    const SubjectData &sd,
    const std::vector<PhenoSpec> &specs,
    RegressionModel m,
    const Eigen::MatrixXd &covarUnion,
    const EngineOptions &opts);

// Write the fitted residual matrix as a strict-format file consumable by
// SubjectData::loadResidOne via --pheno PREFIX.null.resid --resid-name ....
//   Header: IID<TAB>name1<TAB>name2...
//   Rows:   sd.usedIIDs() order; NaN entries written as "NA".
void writeResidualsFile(
    const std::string &outPath,
    const SubjectData &sd,
    const std::vector<NullModelFit> &fits);

} // namespace nullmodel
