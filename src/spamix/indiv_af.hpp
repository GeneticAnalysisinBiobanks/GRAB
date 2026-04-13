// indiv_af.hpp — Per-individual allele-frequency model:
//   compute, write (binary / text / gz), and reconstruct AF vectors.
//   Also: runSPAmixAF — pre-compute AF models for all markers.
//
// Binary format (.bin):
//   No header. Records indexed by raw .bim line (genoIndex).
//   filePos = genoIndex * recordSize
//   recordSize = sizeof(int8_t) + (1+nPC)*sizeof(double)
//   status int8_t: 0 = uniform, 1 = OLS, 2 = logistic
//   betas  (1+nPC) doubles: intercept first, then one per PC
//   Unprocessed markers: all zeros (status=0).
//
// Text / gz format:
//   Header: #STATUS\tBETA0\tBETA1\t...\n
//   One row per filtered marker, tab-separated, doubles as %.17g.
//
// Shared between SPAmixAF (pre-compute step) and SPAmix (on-the-fly).
#pragma once

#include "geno_factory/geno_data.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

// ======================================================================
// AFModel — result of fitting one marker
// ======================================================================

struct AFModel {
    int8_t status;         // 0=uniform, 1=OLS, 2=logistic
    Eigen::VectorXd betas; // length (1+nPC); all zeros for status 0
};

// ======================================================================
// AFContext — pre-computed OLS matrices for AF estimation
// ======================================================================

struct AFContext {
    const Eigen::MatrixXd &onePlusPCs;        // (N × 1+nPC) = [1 | PCs]
    const Eigen::MatrixXd &XtX_inv_Xt;        // (1+nPC × N) = (X'X)^{-1} X'
    const Eigen::VectorXd &sqrt_XtX_inv_diag; // (1+nPC) = sqrt(diag((X'X)^{-1}))
    Eigen::Ref<const Eigen::MatrixXd> PCs;    // (N × nPC) — accepts MatrixXd or block
    int N;
    int nPC;
};

// ======================================================================
// computeAFModel — fit model, return status + betas (for pre-compute storage)
//
// Cascade:
//   MAC <= 20               → status 0 (uniform altFreq)
//   OLS fit good            → status 1 (betas = full (1+nPC) OLS coefs)
//   OLS bad, no sig PCs     → status 0
//   logistic on sig PCs     → status 2 (betas mapped to full (1+nPC) vector)
// ======================================================================

AFModel computeAFModel(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    double altFreq,
    const AFContext &ctx
);

// ======================================================================
// computeAFVec — fused compute: directly writes per-individual AF into out
//
// Same cascade as computeAFModel, but writes AF values directly without
// materialising the intermediate betas vector.  Avoids redundant matrix-
// vector multiplies on the OLS path and uses the smaller sigPCs matrix
// for the logistic path.  Preferred for on-the-fly usage (SPAmix).
// ======================================================================

void computeAFVec(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    double altFreq,
    const AFContext &ctx,
    Eigen::Ref<Eigen::VectorXd> out
);

// ======================================================================
// getAFVecFromModel — reconstruct per-individual AF from a stored model
//
//  model    — AFModel as returned by computeAFModel / loaded from file
//  altFreq  — overall alt allele frequency (used for status 0)
//  onePlusPCs — (N × 1+nPC)
//  N        — sample size
//  out      — output vector (size N), caller-allocated
// ======================================================================

void getAFVecFromModel(
    const AFModel &model,
    double altFreq,
    const Eigen::MatrixXd &onePlusPCs,
    int N,
    Eigen::Ref<Eigen::VectorXd> out
);

// ======================================================================
// IndivAFWriter — write all marker AF models to disk
//
// Usage:
//   IndivAFWriter w(outputFile, nBimMarkers, nPC, markerInfo);
//   w.write(genoIndex, chr, bp, status, betas);  // called per marker
//   w.close();
// ======================================================================

// Forward-declared from plink.hpp — avoid circular include.
// Caller passes CHR / BP as strings / uint32_t directly.

class IndivAFWriter {
  public:
    enum class Mode { Binary, Text };

// Infer mode from outputFile extension: .bin → Binary, else Text.
// Text mode supports .gz, .zst, or plain via TextWriter.
    IndivAFWriter(
        const std::string &outputFile,
        uint64_t nBimMarkers,           // total .bim lines (binary pre-alloc)
        int nPC
    );

    ~IndivAFWriter();

// Write one record.  For Binary, genoIndex determines the file offset.
    void write(
        uint64_t genoIndex,
        int8_t status,
        const Eigen::VectorXd &betas
    );

    void close();

    Mode mode() const {
        return m_mode;
    }

  private:
    Mode m_mode;
    int m_nPC;
    long long m_recordSize; // valid only for Binary

// Binary output
    std::fstream m_binOut;

// Text output (TextWriter handles gz/zst/plain)
    std::unique_ptr<class TextWriter> m_writer;

    bool m_closed = false;
};

// ======================================================================
// IndivAFReader — load one record from a binary AF file by genoIndex
//
//   model.betas will have length (1+nPC).
// ======================================================================

class IndivAFReader {
  public:
    IndivAFReader(
        const std::string &binFile,
        int nPC
    );
    ~IndivAFReader();

// Seek to genoIndex and fill model.  Returns false if status = 0 (caller
// may still use the returned status-0 model for a uniform AF estimate).
    bool read(
        uint64_t genoIndex,
        AFModel &model
    );

  private:
    std::ifstream m_in;
    int m_nPC;
    long long m_recordSize;
};

// ======================================================================
// loadAFModels — load pre-computed AF models from file
//
// Dispatches by extension:
//   .bin  → binary random-access (seeks by genoIndices[fi])
//   .gz / other → tab-separated text, read sequentially
//
// genoIndices[fi] = markerInfo[fi].genoIndex (used only for .bin).
// Always returns exactly nMarkers models in flat marker order.
// ======================================================================

std::vector<AFModel>loadAFModels(
    const std::string &path,
    int nPC,
    uint32_t nMarkers,
    const std::vector<uint64_t> &
    genoIndices
);

// ======================================================================
// runSPAmixAF — pre-compute per-marker AF models and write to disk.
//
// The eigenvector file provides both subject IDs (first column) and PCs.
//
// Command:
//   build/grab --method SPAmixAF --eigen-vecs <PCs> --bfile <prefix>
//              --out <out[.bin|.gz|.txt]>
// ======================================================================

// Compute AF models for all markers in parallel, storing results in memory.
// genoToFlat maps raw genoIndex → flat marker index (UINT32_MAX = skip).
// Callers must also include "geno_factory/plink.hpp".
class GenoMeta;
std::vector<AFModel> computeAFModelsInMemory(
    const GenoMeta &plinkData,
    const AFContext &afCtx,
    const std::vector<uint32_t> &genoToFlat,
    int nthread
);

void runSPAmixAF(
    const std::vector<std::string> &pcColNames,
    const std::string &phenoFile,
    const std::string &covarFile,
    const GenoSpec &geno,
    const std::string &outputFile,
    int nthread,
    int nSnpPerChunk,
    double missingCutoff,
    double minMafCutoff,
    double minMacCutoff,
    double hweCutoff,
    const std::string &keepFile = {},
    const std::string &removeFile = {}

);
