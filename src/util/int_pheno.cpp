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
    const std::string &outPath
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
    if (header.size() < 3 || header[0] != "FID" || header[1] != "IID")
        throw std::runtime_error(
            "runIntPheno: first two header columns must be 'FID' 'IID' followed "
            "by at least one trait column; got: " + line);

    const std::size_t nTraits = header.size() - 2;

    infoMsg("INT pheno: input = %s", phenoFile.c_str());
    infoMsg("INT pheno: traits = %zu", nTraits);

    // ── Read rows ──────────────────────────────────────────────────────
    std::vector<std::string> fids, iids;
    fids.reserve(1u << 14);
    iids.reserve(1u << 14);
    // cols[t][row] holds raw value (NaN = missing).
    std::vector<std::vector<double>> cols(nTraits);
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

        fids.push_back(std::move(toks[0]));
        iids.push_back(std::move(toks[1]));
        for (std::size_t t = 0; t < nTraits; ++t) {
            const std::string &s = toks[2 + t];
            if (isMissingToken(s)) {
                cols[t].push_back(std::numeric_limits<double>::quiet_NaN());
                continue;
            }
            const char *cstr = s.c_str();
            char *endp = nullptr;
            double v = std::strtod(cstr, &endp);
            if (endp == cstr || *endp != '\0') {
                std::ostringstream msg;
                msg << "runIntPheno: " << phenoFile << " line " << lineNo
                    << ", column '" << header[2 + t]
                    << "': non-numeric value '" << s << "'";
                throw std::runtime_error(msg.str());
            }
            cols[t].push_back(v);
        }
    }
    ifs.close();

    const std::size_t N = fids.size();
    infoMsg("INT pheno: rows = %zu", N);
    if (N == 0) throw std::runtime_error("runIntPheno: no data rows in " + phenoFile);

    // ── Apply INT per trait (non-missing scope) ────────────────────────
    for (std::size_t t = 0; t < nTraits; ++t) {
        std::vector<Eigen::Index> okIdx;
        okIdx.reserve(N);
        for (std::size_t i = 0; i < N; ++i)
            if (!std::isnan(cols[t][i])) okIdx.push_back(static_cast<Eigen::Index>(i));

        const std::size_t m = okIdx.size();
        if (m < 2) {
            warnMsg("trait '%s': %zu non-missing value(s) — leaving as NA",
                    header[2 + t].c_str(), m);
            for (std::size_t i = 0; i < N; ++i)
                cols[t][i] = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        Eigen::VectorXd y(m);
        for (std::size_t k = 0; k < m; ++k) y[static_cast<Eigen::Index>(k)] = cols[t][okIdx[k]];

        Eigen::VectorXd z = math::inverseRankNormal(y);

        for (std::size_t i = 0; i < N; ++i)
            cols[t][i] = std::numeric_limits<double>::quiet_NaN();
        for (std::size_t k = 0; k < m; ++k)
            cols[t][okIdx[k]] = z[static_cast<Eigen::Index>(k)];

        infoMsg("trait '%s': n_obs = %zu / %zu  → INT applied",
                header[2 + t].c_str(), m, N);
    }

    // ── Write output ───────────────────────────────────────────────────
    std::ofstream ofs(outPath);
    if (!ofs) throw std::runtime_error("runIntPheno: cannot write " + outPath);

    ofs << header[0];
    for (std::size_t c = 1; c < header.size(); ++c) ofs << '\t' << header[c];
    ofs << '\n';

    char buf[64];
    for (std::size_t i = 0; i < N; ++i) {
        ofs << fids[i] << '\t' << iids[i];
        for (std::size_t t = 0; t < nTraits; ++t) {
            const double v = cols[t][i];
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
