// polmm.hpp — POLMM: Proportional Odds Logistic Mixed Model
//
// All-in-one ordinal GWAS: fits null model internally, then runs
// marker-level association with SPA for ordinal phenotypes.
// Uses sparse GRM for random effects via PCG solver.
//
// Reference: Bi et al. (2021) "POLMM: A fast and accurate method
//   for detecting association in ordinal categorical data"
#pragma once

#include "engine/marker.hpp"
#include "geno_factory/geno_data.hpp"
#include "io/sparse_grm.hpp"

#include <Eigen/Dense>
#include <memory>
#include <string>
#include <vector>

// ======================================================================
// POLMMNullModel — fitted null model parameters and marker-test matrices
// ======================================================================

struct POLMMNullModel {
    int n; // number of subjects
    int J; // number of ordinal categories
    int p; // number of covariates (including intercept)

    // Null model parameters
    Eigen::VectorXd beta; // (p) fixed-effect coefficients
    Eigen::VectorXd eps;  // (J-1) cutpoints
    Eigen::VectorXd bVec; // (n) random effects
    double tau;           // variance component

    // Fitted probabilities and working matrices
    Eigen::MatrixXd muMat; // (n, J) category probabilities
    Eigen::MatrixXd iRMat; // (n, J-1) inverse working-correlation scales

    // Precomputed per-marker test matrices (computed once after null model)
    Eigen::MatrixXd XXR_Psi_RX_new; // (n, p) for covariate adjustment
    Eigen::MatrixXd XR_Psi_R_new;   // (p, n) for covariate adjustment
    Eigen::VectorXd RymuVec;        // (n) score vector constant
    Eigen::VectorXd RPsiR;          // (n) per-subject variance weight

    // MAC-binned variance ratios (Optimization 4)
    // Boundaries are upper-exclusive: bin i covers [macBounds[i], macBounds[i+1])
    // The last bin covers [macBounds.back(), +inf).
    std::vector<double> macBounds; // bin boundaries
    std::vector<double> varRatios; // one VR per bin
};

// ======================================================================
// POLMMMethod — MethodBase implementation for marker-level testing
// ======================================================================

class POLMMMethod : public MethodBase {
  public:
    POLMMMethod(
        const POLMMNullModel &null,
        double spaCutoff
    );

    std::unique_ptr<MethodBase> clone() const override;

    int resultSize() const override {
        return 4;
    }

    std::string getHeaderColumns() const override {
        return "\tPOLMM_P\tPOLMM_Z\tPOLMM_BETA\tPOLMM_SE";
    }

    void getResultVec(
        Eigen::Ref<Eigen::VectorXd> GVec,
        double altFreq,
        int markerInChunkIdx,
        std::vector<double> &result
    ) override;

  private:
// Shared (read-only) from null model
    int m_n, m_J;
    Eigen::MatrixXd m_XXR_Psi_RX_new; // (n, p)
    Eigen::MatrixXd m_XR_Psi_R_new;   // (p, n)
    Eigen::VectorXd m_RymuVec;        // (n)
    Eigen::VectorXd m_RPsiR;          // (n)
    double m_spaCutoff;

// MAC-binned variance ratios
    std::vector<double> m_macBounds;
    std::vector<double> m_varRatios;
    double lookupVarRatio(double mac) const;

// SPA data (per-subject ordinal distribution)
    Eigen::MatrixXd m_muMat; // (n, J) or (n, J-1)
    Eigen::MatrixXd m_iRMat; // (n, J-1)

// Per-thread scratch space
    Eigen::VectorXd m_adjG; // (n)
    Eigen::VectorXd m_tmpP; // (p)

// SPA internals
    double spaTest(
        double Stat,
        double VarW
    ) const;

    double spaTestBinary(
        double Stat,
        double VarW
    ) const;                                              // K=2 fast path

    double K0(
        double t,
        const double *cMat,
        const double *muSub,
        int nSub,
        int Jm1,
        double Ratio0,
        double m1
    ) const;

    double K1(
        double t,
        const double *cMat,
        const double *muSub,
        int nSub,
        int Jm1,
        double Ratio0,
        double m1
    ) const;

    double K2(
        double t,
        const double *cMat,
        const double *muSub,
        int nSub,
        int Jm1,
        double Ratio0
    ) const;

};

// ======================================================================
// Entry point
// ======================================================================

void runPOLMM(
    const std::string &phenoFile,
    const std::string &covarFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const GenoSpec &geno,
    const std::string &spgrmGrabFile,
    const std::string &spgrmGctaFile,
    const std::string &outPrefix,
    const std::string &compression,
    int compressionLevel,
    double spaCutoff,
    int nthreads,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {}

);
