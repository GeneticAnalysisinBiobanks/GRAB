// int_pheno.cpp — Inverse-normal-transform utility for phenotype files.

#include "util/int_pheno.hpp"
#include "util/logging.hpp"
#include "util/math_helper.hpp"

#include <Eigen/Core>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

static bool isMissingToken(const std::string &s) {
    return s.empty() || s == "." || s == "NA" || s == "na" ||
           s == "NaN" || s == "nan" || s == "-";
}

static std::vector<std::string> splitWhitespace(const std::string &line) {
    std::vector<std::string> toks;
    const char *p   = line.c_str();
    const char *end = p + line.size();
    while (p < end) {
        while (p < end && (*p == ' ' || *p == '\t')) ++p;
        if (p >= end) break;
        const char *start = p;
        while (p < end && *p != ' ' && *p != '\t') ++p;
        toks.emplace_back(start, p - start);
    }
    return toks;
}

void runIntPheno(
    const std::string &phenoFile,
    const std::string &outPath,
    const std::vector<std::string> &traitNames
) {
    if (phenoFile.empty()) throw std::runtime_error("runIntPheno: --pheno is required");
    if (outPath.empty())   throw std::runtime_error("runIntPheno: --out is required");

    std::ifstream ifs(phenoFile);
    if (!ifs) throw std::runtime_error("runIntPheno: cannot open " + phenoFile);

    // ── Header ─────────────────────────────────────────────────────────
    std::string line;
    if (!std::getline(ifs, line))
        throw std::runtime_error("runIntPheno: empty file " + phenoFile);
    if (!line.empty() && line.back() == '\r') line.pop_back();

    std::vector<std::string> header = splitWhitespace(line);
    if (header.size() < 2)
        throw std::runtime_error(
            "runIntPheno: header must have at least 2 columns (ID + data); got: " + line);

    // Detect ID-column layout (matches GRAB's parseIIDFile conventions).
    //   "#FID IID ..." / "FID IID ..."  → iidCol = 1, dataStart = 2
    //   "#IID ..."     / "IID ..."      → iidCol = 0, dataStart = 1
    //   anything else                   → iidCol = 0, dataStart = 1
    std::size_t iidCol   = 0;
    std::size_t dataStart = 1;
    if (header.size() >= 3 &&
        (header[0] == "#FID" || header[0] == "FID") &&
        header[1] == "IID")
    {
        iidCol    = 1;
        dataStart = 2;
    }

    const std::size_t nInputTraits = header.size() - dataStart;
    if (nInputTraits == 0)
        throw std::runtime_error(
            "runIntPheno: no trait columns after ID column(s); got: " + line);

    // ── Resolve which trait columns to transform & output ──────────────
    std::vector<std::size_t> keepCols;  // absolute column indices into header
    if (traitNames.empty()) {
        for (std::size_t c = dataStart; c < header.size(); ++c)
            keepCols.push_back(c);
    } else {
        std::unordered_set<std::string> wanted(traitNames.begin(), traitNames.end());
        for (std::size_t c = dataStart; c < header.size(); ++c) {
            if (wanted.count(header[c])) keepCols.push_back(c);
        }
        if (keepCols.size() != traitNames.size()) {
            std::ostringstream msg;
            msg << "runIntPheno: --pheno-name requested " << traitNames.size()
                << " columns but only " << keepCols.size()
                << " were found in " << phenoFile << " header";
            throw std::runtime_error(msg.str());
        }
    }

    const std::size_t nKept = keepCols.size();
    infoMsg("INT pheno: input = %s", phenoFile.c_str());
    infoMsg("INT pheno: traits (input) = %zu, (output) = %zu",
            nInputTraits, nKept);

    // ── Read rows ──────────────────────────────────────────────────────
    std::vector<std::string> fids, iids;
    iids.reserve(1u << 14);
    if (iidCol == 1) fids.reserve(1u << 14);

    std::vector<std::vector<double>> cols(nKept);
    for (auto &c : cols) c.reserve(1u << 14);

    std::size_t lineNo = 1;
    while (std::getline(ifs, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        std::vector<std::string> toks = splitWhitespace(line);
        if (toks.size() != header.size()) {
            std::ostringstream msg;
            msg << "runIntPheno: " << phenoFile << " line " << lineNo
                << ": expected " << header.size() << " fields, got " << toks.size();
            throw std::runtime_error(msg.str());
        }

        if (iidCol == 1) fids.push_back(std::move(toks[0]));
        iids.push_back(std::move(toks[iidCol]));

        for (std::size_t k = 0; k < nKept; ++k) {
            const std::string &s = toks[keepCols[k]];
            if (isMissingToken(s)) {
                cols[k].push_back(std::numeric_limits<double>::quiet_NaN());
                continue;
            }
            const char *cstr = s.c_str();
            char *endp = nullptr;
            double v = std::strtod(cstr, &endp);
            if (endp == cstr || *endp != '\0') {
                std::ostringstream msg;
                msg << "runIntPheno: " << phenoFile << " line " << lineNo
                    << ", column '" << header[keepCols[k]]
                    << "': non-numeric value '" << s << "'";
                throw std::runtime_error(msg.str());
            }
            cols[k].push_back(v);
        }
    }
    ifs.close();

    const std::size_t N = iids.size();
    infoMsg("INT pheno: rows = %zu", N);
    if (N == 0) throw std::runtime_error("runIntPheno: no data rows in " + phenoFile);

    // ── Apply INT per trait (non-missing scope) ────────────────────────
    for (std::size_t k = 0; k < nKept; ++k) {
        std::vector<Eigen::Index> okIdx;
        okIdx.reserve(N);
        for (std::size_t i = 0; i < N; ++i)
            if (!std::isnan(cols[k][i])) okIdx.push_back(static_cast<Eigen::Index>(i));

        const std::size_t m = okIdx.size();
        if (m < 2) {
            warnMsg("trait '%s': %zu non-missing value(s) — leaving as NA",
                    header[keepCols[k]].c_str(), m);
            for (std::size_t i = 0; i < N; ++i)
                cols[k][i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        Eigen::VectorXd y(m);
        for (std::size_t j = 0; j < m; ++j) y[static_cast<Eigen::Index>(j)] = cols[k][okIdx[j]];

        Eigen::VectorXd z = math::inverseRankNormal(y);

        for (std::size_t i = 0; i < N; ++i)
            cols[k][i] = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t j = 0; j < m; ++j)
            cols[k][okIdx[j]] = z[static_cast<Eigen::Index>(j)];

        infoMsg("trait '%s': n_obs = %zu / %zu  → INT applied",
                header[keepCols[k]].c_str(), m, N);
    }

    // ── Write output (preserves input ID-column layout) ────────────────
    std::ofstream ofs(outPath);
    if (!ofs) throw std::runtime_error("runIntPheno: cannot write " + outPath);

    if (iidCol == 1) ofs << header[0] << '\t' << header[1];
    else             ofs << header[0];
    for (std::size_t k = 0; k < nKept; ++k) ofs << '\t' << header[keepCols[k]];
    ofs << '\n';

    char buf[64];
    for (std::size_t i = 0; i < N; ++i) {
        if (iidCol == 1) ofs << fids[i] << '\t' << iids[i];
        else             ofs << iids[i];
        for (std::size_t k = 0; k < nKept; ++k) {
            const double v = cols[k][i];
            if (std::isnan(v)) {
                ofs << "\tNA";
            } else {
                std::snprintf(buf, sizeof(buf), "%.9g", v);
                ofs << '\t' << buf;
            }
        }
        ofs << '\n';
    }
    ofs.close();

    infoMsg("INT pheno: wrote %s", outPath.c_str());
}
