// spacox.hpp — SPACox: Saddlepoint Approximation for Cox Model
//
// Port of ref_code/src/mtSPACox into the pure C++17 / Eigen framework.
//
// Workflow:
//   1. Pre-compute empirical CGF interpolation tables from residuals
//   2. Per-marker: score test → normal approx or SPA tail probability
//      Stage 1: unadjusted genotype
//      Stage 2: covariate-adjusted if p < pVal_covaAdj_cutoff
//   3. Output: [Pvalue, zScore]
#pragma once

#include <Eigen/Dense>
#include <string>
#include <vector>

#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"

// ======================================================================
// DesignMatrix — covariate projection (spacox-only)
// ======================================================================

class DesignMatrix {
  public:
    explicit DesignMatrix(const Eigen::MatrixXd &X);

    int nRows() const {
        return static_cast<int>(m_X.rows());
    }

    int nCols() const {
        return static_cast<int>(m_X.cols());
    }

    const Eigen::MatrixXd &X() const {
        return m_X;
    }

    const Eigen::MatrixXd &tX() const {
        return m_tX;
    }

    const Eigen::MatrixXd &XinvXX() const {
        return m_XinvXX;
    }

    void adjustGenotype(
        const double *G,
        const uint32_t *nzIdx,
        int nNz,
        Eigen::Ref<Eigen::VectorXd> adjG
    ) const;

  private:
    Eigen::MatrixXd m_X;      // N × p
    Eigen::MatrixXd m_tX;     // p × N
    Eigen::MatrixXd m_XinvXX; // N × p
};

// ======================================================================
// Empirical CGF interpolation table
// ======================================================================

struct CumulantTable {
    Eigen::VectorXd xGrid; // length L, strictly increasing
    int nGrid;
    double invScale; // 1/scale for O(1) bin lookup
    double Lp1;      // L + 1
    Eigen::VectorXd yK0, slopeK0;
    Eigen::VectorXd yK1, slopeK1;
    Eigen::VectorXd yK2, slopeK2;
};

// Build the CGF table from residuals using Cauchy-quantile spacing.
// Hardcoded: range = [-100, 100], length = 10000.
CumulantTable buildCumulantTable(const Eigen::VectorXd &residuals);

// ======================================================================
// SPACoxMethod — MethodBase implementation
// ======================================================================

class SPACoxMethod : public MethodBase {
  public:
// Construct from pre-computed CGF table, residuals, and design matrix.
// All large const objects are taken by const reference (shared, read-only).
    SPACoxMethod(
        const Eigen::VectorXd &residuals,
        double varResid,
        const CumulantTable &cumul,
        const DesignMatrix &design,
        double pvalCovAdjCut,
        double spaCutoff
    );

// ---- MethodBase interface ----
    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 2;
    }

    std::string getHeaderColumns() const override {
        return "\tP\tZ";
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result
    ) override;

  private:
// ---- CGF interpolation (O(1) Cauchy-inverse index) ----
    int interpIdx(double v) const;

    double interp(
        const double *yp,
        const double *sp,
        int lo,
        double v
    ) const;

    double interpK0(double v) const;

    double interpK1(double v) const;

    double interpK2(double v) const;

// ---- Cumulant evaluation (scalar loops, no heap alloc) ----
    double evalK0(
        double t,
        int N0,
        double adjG0,
        const double *adjG,
        const uint32_t *idx,
        int n
    ) const;

// Fused K1+K2 evaluation — single loop, one interp lookup per element
    std::pair<double, double>evalK1K2(
        double t,
        int N0,
        double adjG0,
        const double *adjG,
        const uint32_t *idx,
        int n,
        double q2
    ) const;

// ---- SPA root-finding & tail probability ----
    struct RootResult {
        double root;
        bool converge;
        double K2;
    };

    RootResult fastGetRootK1(
        double initX,
        int N0,
        double adjG0,
        const double *adjG,
        const uint32_t *idx,
        int n,
        double q2
    ) const;

    double getProbSpa(
        double adjG0,
        const double *adjG,
        const uint32_t *idx,
        int n,
        int N0,
        double q2,
        bool lowerTail
    ) const;

// ---- Per-marker score test ----
    double getMarkerPval(
        const Eigen::Ref<Eigen::VectorXd> &GVec,
        double MAF,
        double &zScore
    );

// Read-only shared data (references are stable because the owner outlives all clones)
    const Eigen::VectorXd &m_resid;
    double m_varResid;
    const CumulantTable &m_cumul;
    const DesignMatrix &m_design;
    int m_N;
    double m_pvalCovAdjCut;
    double m_spaCutoff;

// Per-thread scratch (mutable, cloned)
    Eigen::VectorXd m_adjGNorm;
    Eigen::VectorXd m_adjGVec;
    std::vector<uint32_t> m_nzSet;
};

// ======================================================================
// Orchestration
// ======================================================================

void runSPACox(
    const std::vector<std::string> &residNames,
    const std::vector<std::string> &covarNames,             // empty = no covariates
    const std::string &phenoFile,                            // pheno file (for residuals + covar columns)
    const std::string &covarFile,                          // covar file (for covar columns)
    const GenoSpec &geno,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double pvalCovAdjCut,
    double spaCutoff,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {}
);
