// subject_data.hpp — Unified per-subject data: residuals, covariates, PCs
//
// Strict file format (all per-subject files):
//   - Exactly one header line (required)
//   - Header names must match [0-9A-Za-z_\-.]+
//   - First column = subject ID (any header name)
//   - All subsequent columns = data (tab or whitespace delimited)
//   - Missing values: NA, NaN, ., N/A, NULL (case-insensitive)
//
// Residuals:
//   loadResidOne — IID + 1..N residual columns (selected by --resid-name)
//
// Phenotype / covariate files:
//   pheno  — IID + named columns (selected by --pheno-name)
//   covar  — IID + named columns (selected by --covar-name)
//   Columns selected by name via getColumn() / getColumns() after finalize.
#pragma once

#include "io/subject_set.hpp"

#include <Eigen/Dense>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Parse .fam file → IIDs (column 2, 0-indexed).
std::vector<std::string> parseFamIIDs(const std::string &famFile);

// Per-phenotype mask info for multi-residual GWAS.
struct PerPhenoInfo {
    std::string name;   // column header name (e.g. "Pheno1")
    uint32_t nUsed = 0; // subjects non-NA in this column
    std::vector<uint32_t>
    unionToLocal;     // size = nUnionUsed; maps union-dense index → pheno-dense index (UINT32_MAX if absent)
};

// Extract dense per-phenotype vector from union-dimension vector.
inline Eigen::VectorXd extractPhenoVec(
    const Eigen::VectorXd &unionVec,
    const PerPhenoInfo &info
) {
    Eigen::VectorXd out(info.nUsed);
    for (Eigen::Index i = 0; i < unionVec.size(); ++i) {
        uint32_t li = info.unionToLocal[static_cast<size_t>(i)];
        if (li != UINT32_MAX) out[li] = unionVec[i];
    }
    return out;
}

// Extract dense per-phenotype matrix from union-dimension matrix.
inline Eigen::MatrixXd extractPhenoMat(
    const Eigen::MatrixXd &unionMat,
    const PerPhenoInfo &info
) {
    Eigen::MatrixXd out(info.nUsed, unionMat.cols());
    for (Eigen::Index i = 0; i < unionMat.rows(); ++i) {
        uint32_t li = info.unionToLocal[static_cast<size_t>(i)];
        if (li != UINT32_MAX) out.row(li) = unionMat.row(i);
    }
    return out;
}

class SubjectData {
  public:
// Construct with .fam IIDs (caller obtains from parseFamIIDs).
    explicit SubjectData(std::vector<std::string> famIIDs);

// Read data column names from a pheno/covar file header (lightweight,
// no data rows parsed).  Handles the FID+IID header convention.
    static std::vector<std::string> readColumnNames(const std::string &filename);

// ── Residual loaders (call exactly one before finalize) ────────────
    void loadResidOne(
        const std::string &filename,
        const std::vector<std::string> &residNames = {}

    );                                                                  // IID + selected resid cols

// ── Optional per-subject files ─────────────────────────────────────
    void loadPhenoFile(
        const std::string &filename,
        const std::vector<std::string> &neededNames = {}

    );                                                                   // strict format: IID + named columns

// Load covariate file.  neededNames selects columns by name (empty = all).
// When neededNames is empty, loads all data columns (columns 1..N).
    void loadCovar(
        const std::string &filename,
        const std::vector<std::string> &neededNames = {}

    );

// ── Subject filter (optional, call before finalize) ─────────────────
// Apply --keep / --remove subject-ID files (PLINK2-compatible format).
// Either path may be empty (no-op).  Must be called before finalize().
    void setKeepRemove(
        const std::string &keepFile,
        const std::string &removeFile
    );

// Set GRM subject IDs for intersection (GRM step).
// Call before finalize().  Empty set → skip GRM intersection.
    void setGrmSubjects(std::unordered_set<std::string> grmIDs);

// Set descriptive labels for the subject pipeline log.
    void setGenoLabel(const std::string &label) {
        m_subjectSet.setGenoLabel(label);
    }

    void setGrmLabel(const std::string &label) {
        m_subjectSet.setGrmLabel(label);
    }

// Access the underlying SubjectSet (valid after finalize).
    const SubjectSet &subjectSet() const {
        return m_subjectSet;
    }

// ── Build bitmask and reorder to .fam order ────────────────────────
// Intersects all loaded IID sets with .fam, builds bitmask, reorders
// dense data to .fam order.  Must be called exactly once, after all
// load* calls are complete.
    void finalize();

// ── Accessors (valid after finalize) ───────────────────────────────

    uint32_t nFam() const {
        return m_subjectSet.nFam();
    }

    uint32_t nUsed() const {
        return m_nUsed;
    }

    const std::vector<std::string> &famIIDs() const {
        return m_subjectSet.famIIDs();
    }

    const std::vector<uint64_t> &usedMask() const {
        return m_usedMask;
    }

    uint32_t nMaskWords() const {
        return static_cast<uint32_t>(m_usedMask.size());
    }

// Residuals — size nUsed, dense, .fam order among used subjects.
// Populated by residOne or residWtCoxG.
    const Eigen::VectorXd &residuals() const {
        return m_residuals;
    }

// Number of residual columns loaded by loadResidOne (1 when single-col).
    int residOneCols() const {
        return m_nResidOneCols;
    }

// Column names from header (empty for legacy headerless files).
    const std::vector<std::string> &residColNames() const {
        return m_residColNames;
    }

// Weights / indicator — only valid after residWtCoxG or setResidWeightIndicator.
    const Eigen::VectorXd &weights() const {
        return m_weights;
    }

    const Eigen::VectorXd &indicator() const {
        return m_indicator;
    }

// Residual matrix — nUsed × K, only valid after residSPAsqr.
    int residCols() const {
        return static_cast<int>(m_residMatrix.cols());
    }

    const Eigen::MatrixXd &residMatrix() const {
        return m_residMatrix;
    }

// Covariate matrix — nUsed × p (empty if not loaded).
    bool hasCovar() const {
        return m_covar.size() > 0;
    }

    int covarCols() const {
        return static_cast<int>(m_covar.cols());
    }

    const Eigen::MatrixXd &covar() const {
        return m_covar;
    }

    const std::vector<std::string> &covarColNames() const {
        return m_covarColNames;
    }

// ── Named column access (valid after finalize) ─────────────────────
// Searches covar columns first, then pheno columns.
// Throws if name not found.
    bool hasColumn(const std::string &name) const;

    Eigen::VectorXd getColumn(const std::string &name) const;

    Eigen::MatrixXd getColumns(const std::vector<std::string> &names) const;

// Zero-allocation: copy named columns directly into target matrix columns.
    void fillColumnsInto(
        const std::vector<std::string> &names,
        Eigen::Ref<Eigen::MatrixXd> target
    ) const;

// ── Post-finalize NA filtering ─────────────────────────────────────
// Drop subjects with NaN in any of the named pheno/covar columns.
// Updates bitmask, nUsed, and all dense arrays.  Call after finalize().
    void dropNaInColumns(const std::vector<std::string> &colNames);

// ── Post-finalize setters (for computed regression path) ───────────
// Sets residuals, weights, and indicator AFTER finalize, replacing
// whatever was loaded from file.
    void setResidWeightIndicator(
        Eigen::VectorXd resid,
        Eigen::VectorXd weight,
        Eigen::VectorXd ind
    );

// ── Direct initialization from a pre-computed bitmask ──────────────
// Bypasses the normal load/finalize pipeline.  Sets the used-mask,
// nUsed, and residuals/weights/indicator directly.  Useful when the
// caller has already computed per-cluster masks and regression output
// (e.g. LEAF pheno path).
    void initFromMask(
        std::vector<uint64_t> mask,
        uint32_t nUsed,
        Eigen::VectorXd resid,
        Eigen::VectorXd weight,
        Eigen::VectorXd ind
    );

// Used IIDs in .fam order (for SparseGRM construction, logging, etc.).
    std::vector<std::string> usedIIDs() const;

// ── Per-phenotype masks for multi-residual GWAS ────────────────────
// Returns one PerPhenoInfo per residual column.  The union mask
// (usedMask()) covers all subjects non-NA in at least one column.
// Each PerPhenoInfo::unionToLocal maps union-dense index to the
// phenotype's dense index (or UINT32_MAX if the subject is NA in
// that column).  Only valid after finalize() with ResidType::One.
    std::vector<PerPhenoInfo> buildPerColumnMasks() const;

  private:
// ── Raw parsed data (before finalize) ──────────────────────────────
    struct RawFile {
        std::vector<std::string> iids;
        int nCols = 0;
        std::vector<double> vals;          // row-major: vals[row * nCols + col]
        std::vector<std::string> colNames; // from header (always populated)
    };

    static RawFile parseIIDFile(
        const std::string &filename,
        int expectCols,
        const std::vector<std::string> &neededCols = {}
    );

// expectCols: -1 = auto-detect from header;
//              N = require exactly N data columns.

// ── State ──────────────────────────────────────────────────────────
    SubjectSet m_subjectSet;
    uint32_t m_nUsed = 0;
    bool m_finalized = false;

    std::vector<uint64_t> m_usedMask;
    std::vector<uint32_t> m_usedFamIndices; // dense index → .fam index

// Dense per-used-subject data (populated by finalize)
    Eigen::VectorXd m_residuals;
    Eigen::VectorXd m_weights;
    Eigen::VectorXd m_indicator;
    Eigen::MatrixXd m_residMatrix; // N × K
    Eigen::MatrixXd m_covar;       // N × p
    Eigen::MatrixXd m_phenoData;   // N × q

// Column name → column index maps (populated by finalize)
    std::vector<std::string> m_covarColNames;
    std::vector<std::string> m_phenoColNames;
    std::unordered_map<std::string, int> m_covarColMap;
    std::unordered_map<std::string, int> m_phenoColMap;

// Pre-finalize raw storage
    enum class ResidType { None, One };

    ResidType m_residType = ResidType::None;
    RawFile m_rawResid;
    RawFile m_rawCovar;
    RawFile m_rawPheno;
    bool m_hasRawResid = false;
    bool m_hasRawCovar = false;
    bool m_hasRawPheno = false;

    int m_nResidOneCols = 0;
    std::vector<std::string> m_residColNames;

// Covariate column means (stored by loadCovar for finalize mean-fill)
    std::vector<double> m_covarColMeans;
};
