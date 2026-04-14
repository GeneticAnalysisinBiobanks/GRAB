// indiv_af.cpp — Shared per-individual AF model computation, writing, and reading.
// See spamix/indiv_af.hpp for format documentation.

#include "spamix/indiv_af.hpp"

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "util/math_helper.hpp"
#include "util/text_scanner.hpp"
#include "util/text_stream.hpp"

// ======================================================================
// computeAFModel — returns status + betas for file storage
// ======================================================================

AFModel computeAFModel(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    double altFreq,
    const AFContext &ctx
) {
    constexpr double PCs_pval_cut = 0.05;
    constexpr double negRatioCut = 0.1;
    constexpr double MAC_cutoff = 20.0;

    const int N = ctx.N;
    const int nPC = ctx.nPC;

    AFModel model;
    model.status = 0;
    model.betas = Eigen::VectorXd::Zero(1 + nPC);

    const double MAC = altFreq * 2.0 * N;
    if (MAC <= MAC_cutoff) return model;

    // OLS: coef = (X'X)^{-1} X' g
    Eigen::VectorXd coef = ctx.XtX_inv_Xt * g;

    // Fitted genotype (single matvec — reused for negRatio check and rss)
    Eigen::VectorXd fitted = ctx.onePlusPCs * coef;

    // Check fitted AF = fitted_geno / 2 for out-of-bound
    int nError = 0;
    for (int i = 0; i < N; ++i) {
        double af = fitted[i] * 0.5;
        if (af < 0.0 || af > 1.0) ++nError;
    }

    if (static_cast<double>(nError) / N < negRatioCut) {
        model.status = 1;
        model.betas = coef;
        return model;
    }

    // Compute rss from fitted genotype (no redundant matvec)
    const double rss = (g - fitted).squaredNorm();
    const double df = static_cast<double>(N - nPC - 1);
    const double sqrtS2 = std::sqrt(rss / df);

    std::vector<int> sigIdx;
    for (int j = 0; j < nPC; ++j) {
        const double se = ctx.sqrt_XtX_inv_diag[j + 1] * sqrtS2;
        const double t = coef[j + 1] / se;
        if (2.0 * math::pt(-std::abs(t), df) < PCs_pval_cut) sigIdx.push_back(j);
    }

    if (sigIdx.empty()) return model;

    // Binarise genotype for logistic regression
    Eigen::VectorXd g0 = Eigen::VectorXd::Zero(N);
    double macAfter = 0.0;
    for (int i = 0; i < N; ++i) {
        if (g[i] > 0.5) {
            g0[i] = 1.0;
            macAfter += 1.0;
        }
    }

    if (macAfter <= MAC_cutoff) return model;

    const int nSig = static_cast<int>(sigIdx.size());
    Eigen::MatrixXd sigPCs(N, nSig);
    for (int j = 0; j < nSig; ++j)
        sigPCs.col(j) = ctx.PCs.col(sigIdx[j]);

    Eigen::VectorXd subBeta = math::logisticRegressionBeta(sigPCs, g0);

    model.status = 2;
    model.betas[0] = subBeta[0];
    for (int j = 0; j < nSig; ++j)
        model.betas[sigIdx[j] + 1] = subBeta[j + 1];

    return model;
}

// ======================================================================
// computeAFVec — fused: directly writes per-individual AF into out
// ======================================================================

void computeAFVec(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    double altFreq,
    const AFContext &ctx,
    Eigen::Ref<Eigen::VectorXd> out
) {
    constexpr double PCs_pval_cut = 0.05;
    constexpr double negRatioCut = 0.1;
    constexpr double MAC_cutoff = 20.0;

    const int N = ctx.N;
    const int nPC = ctx.nPC;

    const double MAC = altFreq * 2.0 * N;
    if (MAC <= MAC_cutoff) {
        out.setConstant(altFreq);
        return;
    }

    // OLS: coef = (X'X)^{-1} X' g
    Eigen::VectorXd coef = ctx.XtX_inv_Xt * g;

    // Fitted genotype — written into out (single matvec)
    out.noalias() = ctx.onePlusPCs * coef;

    // Check fitted AF = fitted_geno / 2 for out-of-bound
    int nError = 0;
    for (int i = 0; i < N; ++i) {
        double af = out[i] * 0.5;
        if (af < 0.0 || af > 1.0) ++nError;
    }

    if (static_cast<double>(nError) / N < negRatioCut) {
        // Status 1 — OLS: convert to AF and clamp
        for (int i = 0; i < N; ++i) {
            out[i] *= 0.5;
            if (out[i] < 0.0) out[i] = 0.0;
            if (out[i] > 1.0) out[i] = 1.0;
        }
        return;
    }

    // Compute rss from fitted genotype still in out (no redundant matvec)
    const double rss = (g - out).squaredNorm();
    const double df = static_cast<double>(N - nPC - 1);
    const double sqrtS2 = std::sqrt(rss / df);

    std::vector<int> sigIdx;
    for (int j = 0; j < nPC; ++j) {
        const double se = ctx.sqrt_XtX_inv_diag[j + 1] * sqrtS2;
        const double t = coef[j + 1] / se;
        if (2.0 * math::pt(-std::abs(t), df) < PCs_pval_cut) sigIdx.push_back(j);
    }

    if (sigIdx.empty()) {
        out.setConstant(altFreq);
        return;
    }

    // Binarise genotype for logistic regression
    Eigen::VectorXd g0 = Eigen::VectorXd::Zero(N);
    double macAfter = 0.0;
    for (int i = 0; i < N; ++i) {
        if (g[i] > 0.5) {
            g0[i] = 1.0;
            macAfter += 1.0;
        }
    }

    if (macAfter <= MAC_cutoff) {
        out.setConstant(altFreq);
        return;
    }

    const int nSig = static_cast<int>(sigIdx.size());
    Eigen::MatrixXd sigPCs(N, nSig);
    for (int j = 0; j < nSig; ++j)
        sigPCs.col(j) = ctx.PCs.col(sigIdx[j]);

    // Logistic regression → AF = 1 - sqrt(1 - sigmoid(eta))
    // Uses smaller sigPCs matrix instead of full onePlusPCs
    out = math::logisticRegression(sigPCs, g0);
}

// ======================================================================
// getAFVecFromModel
// ======================================================================

void getAFVecFromModel(
    const AFModel &model,
    double altFreq,
    const Eigen::MatrixXd &onePlusPCs,
    int N,
    Eigen::Ref<Eigen::VectorXd> out
) {
    switch (model.status) {
    case 0:
        out.setConstant(altFreq);
        break;

    case 1:
        // OLS: fitted = (onePlusPCs * betas) / 2, clamped to [0,1]
        out.noalias() = onePlusPCs * model.betas;
        out *= 0.5;
        for (int i = 0; i < N; ++i) {
            if (out[i] < 0.0) out[i] = 0.0;
            if (out[i] > 1.0) out[i] = 1.0;
        }
        break;

    case 2: {
        // Logistic: μ = sigmoid(onePlusPCs * betas); AF = 1 - sqrt(1-μ)
        out.noalias() = onePlusPCs * model.betas;
        for (int i = 0; i < N; ++i) {
            const double mu = 1.0 / (1.0 + std::exp(-out[i]));
            out[i] = 1.0 - std::sqrt(1.0 - mu);
        }
        break;
    }

    default:
        out.setConstant(altFreq);
        break;
    }
}

// ======================================================================
// IndivAFWriter — helpers
// ======================================================================

static IndivAFWriter::Mode inferMode(const std::string &path) {
    const auto n = path.size();
    if (n > 4 && path.compare(n - 4, 4, ".bin") == 0) return IndivAFWriter::Mode::Binary;
    return IndivAFWriter::Mode::Text;
}

// ======================================================================
// IndivAFWriter — constructor
// ======================================================================

IndivAFWriter::IndivAFWriter(
    const std::string &outputFile,
    uint64_t nBimMarkers,
    int nPC
)
    : m_mode(inferMode(outputFile)),
      m_nPC(nPC),
      m_recordSize(static_cast<long long>(sizeof(int8_t)) +
                   static_cast<long long>(1 + nPC) * static_cast<long long>(sizeof(double)))
{
    if (m_mode == Mode::Binary) {
        // Pre-allocate the file to nBimMarkers * recordSize zeros so that each
        // marker can be written at its genoIndex-based offset independently.
        const long long totalSize = static_cast<long long>(nBimMarkers) * m_recordSize;

        {
            std::ofstream pre(outputFile, std::ios::binary | std::ios::trunc);
            if (!pre.is_open()) throw std::runtime_error("Cannot create binary AF file: " + outputFile);
            constexpr std::streamsize CHUNK = 1 << 20;
            std::vector<char> zeros(static_cast<size_t>(std::min<long long>(CHUNK, totalSize > 0 ? totalSize : CHUNK)),
                                    0);
            for (long long rem = totalSize; rem > 0;) {
                const std::streamsize n =
                    static_cast<std::streamsize>(std::min<long long>(rem, static_cast<long long>(zeros.size())));
                pre.write(zeros.data(), n);
                rem -= n;
            }
        }

        m_binOut.open(outputFile, std::ios::binary | std::ios::in | std::ios::out);
        if (!m_binOut.is_open()) throw std::runtime_error("Cannot open binary AF file for random write: " + outputFile);

    } else {
        m_writer = std::make_unique<TextWriter>(outputFile);

        // Write TSV header
        std::string hdr = "#STATUS";
        for (int j = 0; j <= nPC; ++j)
            hdr += "\tBETA" + std::to_string(j);
        hdr += "\n";
        m_writer->write(hdr);
    }
}

// ======================================================================
// IndivAFWriter — write one record
// ======================================================================

void IndivAFWriter::write(
    uint64_t genoIndex,
    int8_t status,
    const Eigen::VectorXd &betas
) {
    if (m_mode == Mode::Binary) {
        const long long filePos = static_cast<long long>(genoIndex) * m_recordSize;
        m_binOut.seekp(filePos);
        m_binOut.write(reinterpret_cast<const char *>(&status), sizeof(int8_t));
        m_binOut.write(
            reinterpret_cast<const char *>(betas.data()),
            static_cast<std::streamsize>((1 + m_nPC) * sizeof(double))
        );

    } else {
        char buf[64];
        std::string line;
        line.reserve(64 + 26 * (1 + m_nPC));

        std::snprintf(buf, sizeof(buf), "%d", static_cast<int>(status));
        line += buf;
        for (int j = 0; j <= m_nPC; ++j) {
            line += '\t';
            std::snprintf(buf, sizeof(buf), "%.17g", betas[j]);
            line += buf;
        }
        line += '\n';

        m_writer->write(line);
    }
}

// ======================================================================
// IndivAFWriter — close
// ======================================================================

void IndivAFWriter::close() {
    if (m_closed) return;
    m_closed = true;

    if (m_mode == Mode::Binary && m_binOut.is_open())m_binOut.close();
    else if (m_writer)m_writer->close();
}

IndivAFWriter::~IndivAFWriter()
{
    close();
}

// ======================================================================
// IndivAFReader
// ======================================================================

IndivAFReader::IndivAFReader(
    const std::string &binFile,
    int nPC
)
    : m_nPC(nPC),
      m_recordSize(static_cast<long long>(sizeof(int8_t)) +
                   static_cast<long long>(1 + nPC) * static_cast<long long>(sizeof(double)))
{
    m_in.open(binFile, std::ios::binary);
    if (!m_in.is_open()) throw std::runtime_error("Cannot open binary AF file: " + binFile);
}

IndivAFReader::~IndivAFReader()
{
    m_in.close();
}

bool IndivAFReader::read(
    uint64_t genoIndex,
    AFModel &model
) {
    const long long filePos = static_cast<long long>(genoIndex) * m_recordSize;
    m_in.seekg(filePos);
    if (!m_in.good()) throw std::runtime_error("Binary AF seek failed at genoIndex=" + std::to_string(genoIndex));

    int8_t st = 0;
    m_in.read(reinterpret_cast<char *>(&st), sizeof(int8_t));
    model.status = st;
    model.betas.resize(1 + m_nPC);
    m_in.read(reinterpret_cast<char *>(model.betas.data()), static_cast<std::streamsize>((1 + m_nPC) * sizeof(double)));
    if (!m_in.good()) throw std::runtime_error("Binary AF read failed at genoIndex=" + std::to_string(genoIndex));
    return st != 0;
}

// ======================================================================
// ======================================================================
// loadAFModels — helpers (binary + text/gz) and dispatcher
// ======================================================================

static std::vector<AFModel>loadAFModelsBinary(
    const std::string &path,
    int nPC,
    uint32_t nMarkers,
    const std::vector<
        uint64_t> &genoIndices
) {
    std::vector<AFModel> models(nMarkers);
    IndivAFReader reader(path, nPC);
    for (uint32_t fi = 0; fi < nMarkers; ++fi)
        reader.read(genoIndices[fi], models[fi]);
    return models;
}

static std::vector<AFModel> loadAFModelsText(
    const std::string &path,
    int nPC,
    uint32_t nExpected
) {
    std::vector<AFModel> models;
    models.reserve(nExpected);

    uint32_t lineNo = 0;
    auto parseLine = [&](std::string &line) {
        ++lineNo;
        if (text::skipLine(line)) return;
        text::TokenScanner ts(line);
        auto sv = ts.nextView();
        if (sv.empty()) return;
        char *end;
        long status = std::strtol(sv.data(), &end, 10);
        if (end == sv.data()) throw std::runtime_error(path + " line " + std::to_string(lineNo) +
                                                       ": expected STATUS BETA0 ...");
        AFModel m;
        m.status = static_cast<int8_t>(status);
        m.betas.resize(1 + nPC);
        for (int j = 0; j <= nPC; ++j) {
            ts.skipWS();
            if (ts.atEnd()) throw std::runtime_error(
                          path + " line " + std::to_string(lineNo) + ": expected " +
                          std::to_string(1 + nPC) + " beta values"
            );
            const char *cur = ts.pos();
            m.betas[j] = std::strtod(cur, &end);
            if (end == cur) throw std::runtime_error(
                          path + " line " + std::to_string(lineNo) + ": expected " +
                          std::to_string(1 + nPC) + " beta values"
            );
            ts.p = end;
        }
        models.push_back(std::move(m));
    };

    TextReader reader(path);
    std::string line;
    while (reader.getline(line))
        parseLine(line);

    return models;
}

std::vector<AFModel>loadAFModels(
    const std::string &path,
    int nPC,
    uint32_t nMarkers,
    const std::vector<uint64_t> &
    genoIndices
) {
    const auto len = path.size();
    std::vector<AFModel> models;
    if (len > 4 && path.compare(len - 4, 4, ".bin") == 0)models = loadAFModelsBinary(path, nPC, nMarkers, genoIndices);
    else models = loadAFModelsText(path, nPC, nMarkers);

    if (models.size() != nMarkers)throw std::runtime_error(
                  "AF model count (" + std::to_string(models.size()) +
                  ") does not match marker count (" + std::to_string(nMarkers) + ")"
    );
    return models;
}

// ======================================================================
// runSPAmixAF — pre-compute per-marker AF models (top-level orchestration)
// ======================================================================

#include <atomic>
#include <chrono>
#include <ctime>
#include <exception>
#include <mutex>
#include <sstream>
#include <thread>

#include "geno_factory/geno_data.hpp"
#include "io/subject_data.hpp"
#include "util/logging.hpp"

namespace {

struct GenoQC {
    double altFreq = 0.0;
    double maf = 0.0;
    double mac = 0.0;
    int nValid = 0;
};

GenoQC computeGenoQC(
    const Eigen::Ref<const Eigen::VectorXd> &g,
    int N
) {
    GenoQC qc;
    double sum = 0.0;
    int nMiss = 0;
    for (int i = 0; i < N; ++i) {
        if (std::isnan(g[i]))++nMiss;
        else sum += g[i];
    }
    qc.nValid = N - nMiss;
    if (qc.nValid > 0) qc.altFreq = sum / (2.0 * qc.nValid);
    qc.maf = (qc.altFreq > 0.5) ? (1.0 - qc.altFreq) : qc.altFreq;
    qc.mac = qc.maf * 2.0 * qc.nValid;
    return qc;
}

void imputeMissing(
    Eigen::Ref<Eigen::VectorXd> g,
    int N,
    double altFreq
) {
    const double mean = 2.0 * altFreq;
    for (int i = 0; i < N; ++i)
        if (std::isnan(g[i])) g[i] = mean;
}

} // anonymous namespace

std::vector<AFModel> computeAFModelsInMemory(
    const GenoMeta &plinkData,
    const AFContext &afCtx,
    const std::vector<uint32_t> &genoToFlat,
    int nthread
) {
    const auto &markerInfo = plinkData.markerInfo();
    const auto &chunkIndices = plinkData.chunkIndices();
    const size_t nTotalChunks = chunkIndices.size();
    const uint32_t nMarkers = static_cast<uint32_t>(markerInfo.size());
    const int N = afCtx.N;

    std::vector<AFModel> models(nMarkers);
    for (auto &m : models) {
        m.status = 0;
        m.betas = Eigen::VectorXd::Zero(1 + afCtx.nPC);
    }

    std::atomic<size_t> nextChunk(0);
    std::exception_ptr workerError;
    std::mutex errorMutex;

    auto workerFn = [&]() {
        auto cursor = plinkData.makeCursor();
        Eigen::VectorXd GVec(N);

        while (true) {
            const size_t ci = nextChunk.fetch_add(1);
            if (ci >= nTotalChunks) break;

            const auto &chunk = chunkIndices[ci];
            if (!chunk.empty()) cursor->beginSequentialBlock(chunk.front());

            for (uint64_t gIdx : chunk) {
                const uint32_t flatIdx = genoToFlat[gIdx];
                if (flatIdx == UINT32_MAX) continue;

                cursor->getGenotypesSimple(gIdx, GVec);
                GenoQC qc = computeGenoQC(GVec, N);
                imputeMissing(GVec, N, qc.altFreq);

                try {
                    models[flatIdx] = computeAFModel(GVec, qc.altFreq, afCtx);
                } catch (...) {
                    std::lock_guard<std::mutex> lk(errorMutex);
                    if (!workerError) workerError = std::current_exception();
                    return;
                }
            }
            infoMsg("  AF chunk %zu/%zu done", ci + 1, nTotalChunks);
        }
    };

    const int nt = std::min(nthread, static_cast<int>(nTotalChunks));
    std::vector<std::thread> workers;
    workers.reserve(nt);
    for (int i = 0; i < nt; ++i)
        workers.emplace_back(workerFn);
    for (auto &w : workers)
        w.join();

    if (workerError) std::rethrow_exception(workerError);

    int cnt0 = 0, cnt1 = 0, cnt2 = 0;
    for (uint32_t i = 0; i < nMarkers; ++i) {
        if (models[i].status == 1)++cnt1;
        else if (models[i].status == 2)++cnt2;
        else ++cnt0;
    }
    infoMsg("AF models: status0=%d, status1=%d, status2=%d", cnt0, cnt1, cnt2);

    return models;
}

void runSPAmixAF(
    const std::vector<std::string> &pcColNames,
    const std::string &phenoFile,
    const std::string &covarFile,
    const GenoSpec &geno,
    const std::string &outputFile,
    int nthread,
    int nSnpPerChunk,
    double /*missingCutoff*/,
    double /*minMafCutoff*/,
    double /*minMacCutoff*/,
    double /*hweCutoff*/,
    const std::string &keepFile,
    const std::string &removeFile
) {
    const auto wallStart = std::chrono::steady_clock::now();
    const std::clock_t cpuStart = std::clock();

    // Load pheno / covar files and extract PC columns.
    infoMsg("Loading PC data for AF model (%zu columns)", pcColNames.size());
    auto famIIDs = parseGenoIIDs(geno);
    SubjectData sd(std::move(famIIDs));
    if (!phenoFile.empty()) sd.loadPhenoFile(phenoFile);
    if (!covarFile.empty()) sd.loadCovar(covarFile, pcColNames);
    sd.setKeepRemove(keepFile, removeFile);
    sd.setGenoLabel(geno.flagLabel());
    sd.finalize();

    const int N = static_cast<int>(sd.nUsed());
    const int nPC = static_cast<int>(pcColNames.size());
    Eigen::MatrixXd PCs = sd.getColumns(pcColNames);
    infoMsg("  %u subjects, %d PCs", sd.nUsed(), nPC);

    Eigen::MatrixXd onePlusPCs(N, 1 + nPC);
    onePlusPCs.col(0).setOnes();
    onePlusPCs.rightCols(nPC) = PCs;

    Eigen::MatrixXd XtX = onePlusPCs.transpose() * onePlusPCs;
    Eigen::MatrixXd XtX_inv = XtX.ldlt().solve(Eigen::MatrixXd::Identity(1 + nPC, 1 + nPC));
    Eigen::MatrixXd XtX_inv_Xt = XtX_inv * onePlusPCs.transpose();
    Eigen::VectorXd sqrt_XtX_inv_diag = XtX_inv.diagonal().cwiseSqrt();

    const AFContext afCtx{onePlusPCs, XtX_inv_Xt, sqrt_XtX_inv_diag, PCs, N, nPC};

    auto genoData = makeGenoData(geno, sd.usedMask(), sd.nFam(), sd.nUsed(), nSnpPerChunk);
    infoMsg("  %u subjects, %u markers", genoData->nSubjUsed(), genoData->nMarkers());

    const auto &markerInfo = genoData->markerInfo();
    const auto &chunkIndices = genoData->chunkIndices();
    const size_t nTotalChunks = chunkIndices.size();
    const uint32_t nMarkers = static_cast<uint32_t>(markerInfo.size());
    const uint32_t nBimMarkers = genoData->nMarkers();

    // Build genoIndex → flat index map.
    std::vector<uint32_t> genoToFlat(genoData->nMarkers(), UINT32_MAX);
    for (uint32_t fi = 0; fi < nMarkers; ++fi)
        genoToFlat[markerInfo[fi].genoIndex] = fi;

    infoMsg("Computing AF models (%d thread(s), %zu chunks)...", nthread, nTotalChunks);
    const std::vector<AFModel> allModels = computeAFModelsInMemory(*genoData, afCtx, genoToFlat, nthread);

    infoMsg("Writing output: %s", outputFile.c_str());
    IndivAFWriter writer(outputFile, nBimMarkers, nPC);
    for (uint32_t fi = 0; fi < nMarkers; ++fi) {
        const auto &mi = markerInfo[fi];
        writer.write(mi.genoIndex, allModels[fi].status, allModels[fi].betas);
    }
    writer.close();

    const auto wallEnd = std::chrono::steady_clock::now();
    const double wallSec = std::chrono::duration<double>(wallEnd - wallStart).count();
    const double cpuSec = static_cast<double>(std::clock() - cpuStart) / CLOCKS_PER_SEC;
    infoMsg("Wall time: %.1f seconds, CPU time: %.1f seconds", wallSec, cpuSec);
}
