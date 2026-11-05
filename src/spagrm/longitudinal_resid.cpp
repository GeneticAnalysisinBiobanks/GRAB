// longitudinal_resid.cpp — see longitudinal_resid.hpp.

#include "spagrm/longitudinal_resid.hpp"

#include "spagrm/sageld_fit.hpp"
#include "util/logging.hpp"
#include "util/text_scanner.hpp"

#include <cctype>
#include <fstream>
#include <stdexcept>
#include <unordered_map>

namespace nsLongitudinal {

namespace {

// Trim ASCII whitespace from both ends.
std::string trim(const std::string &s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace(static_cast<unsigned char>(s[a]))) ++a;
    while (b > a && std::isspace(static_cast<unsigned char>(s[b - 1]))) --b;
    return s.substr(a, b - a);
}

// Load IIDs from a PLINK2-style --keep / --remove file (FID IID, or single
// IID column).  Lines beginning with '#' and blank lines are skipped.
std::unordered_set<std::string> loadIIDsFromFile(const std::string &path) {
    std::unordered_set<std::string> out;
    if (path.empty()) return out;
    std::ifstream ifs(path);
    if (!ifs) throw std::runtime_error("Cannot open IID file: " + path);
    std::string line;
    while (std::getline(ifs, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        text::TokenScanner ts(line);
        if (ts.atEnd()) continue;
        ts.nextView();                       // FID (or IID if single column)
        std::string iid;
        if (!ts.atEnd()) {
            iid = ts.next();                 // IID
        } else {
            text::TokenScanner ts2(line);
            iid = ts2.next();
        }
        if (!iid.empty()) out.insert(iid);
    }
    return out;
}

std::unordered_map<std::string, int> idMap(const std::vector<std::string> &ids) {
    std::unordered_map<std::string, int> m;
    m.reserve(ids.size() * 2);
    for (int i = 0; i < static_cast<int>(ids.size()); ++i) m.emplace(ids[i], i);
    return m;
}

} // namespace

std::vector<std::string> splitOutcomeNames(const std::string &commaList) {
    std::vector<std::string> out;
    size_t start = 0;
    while (start <= commaList.size()) {
        size_t comma = commaList.find(',', start);
        std::string tok = trim(commaList.substr(
            start, comma == std::string::npos ? std::string::npos : comma - start));
        if (!tok.empty()) {
            if (tok.find(':') != std::string::npos)
                throw std::runtime_error(
                    "--longitudinal: outcome '" + tok +
                    "' must be a single quantitative column (the Cox TIME:EVENT "
                    "syntax is not valid in longitudinal mode)");
            out.push_back(tok);
        }
        if (comma == std::string::npos) break;
        start = comma + 1;
    }
    if (out.empty())
        throw std::runtime_error("--longitudinal: --pheno-name lists no outcome columns");
    return out;
}

std::unordered_set<std::string> buildKeptSet(
    const std::string &keepFile,
    const std::string &removeFile,
    const std::vector<std::string> &famIIDs,
    const std::unordered_set<std::string> &grmSubjects
) {
    const auto keepSet = loadIIDsFromFile(keepFile);
    const auto removeSet = loadIIDsFromFile(removeFile);
    std::unordered_set<std::string> kept;
    kept.reserve(famIIDs.size());
    for (const auto &iid : famIIDs) {
        if (!keepSet.empty() && !keepSet.count(iid)) continue;
        if (!removeSet.empty() && removeSet.count(iid)) continue;
        if (!grmSubjects.empty() && !grmSubjects.count(iid)) continue;
        kept.insert(iid);
    }
    if (kept.empty())
        throw std::runtime_error(
            "--longitudinal: no subjects survive --keep / --remove / GRM intersection");
    return kept;
}

LongResidResult computeLongitudinalResid(
    const std::string &phenoFile,
    const std::vector<std::string> &phenoNames,
    const std::vector<std::string> &covarNames,
    const std::vector<std::string> &famIIDs,
    const std::unordered_set<std::string> &keptSubjects
) {
    using nsSAGELDFit::LMMFit;
    using nsSAGELDFit::LongPhenoData;

    LongPhenoData data = nsSAGELDFit::parseLongPheno(
        phenoFile, phenoNames, covarNames, /*envName=*/"", famIIDs, keptSubjects);

    const int nSubj = static_cast<int>(data.uniqueIIDs.size());
    if (nSubj == 0)
        throw std::runtime_error("--longitudinal: no subjects in long-format file");
    const int q = static_cast<int>(data.phenoNames.size());

    LongResidResult out;
    out.usedIIDs = data.uniqueIIDs;
    out.names = data.phenoNames;
    out.R_G.resize(q);
    for (int k = 0; k < q; ++k) {
        Eigen::VectorXd y = data.Y.col(k);
        LMMFit fit = nsSAGELDFit::fitRandomInterceptWithCovar(data, y);
        out.R_G[k] = nsSAGELDFit::aggregatePerIID(data, fit.residPerRow);
        infoMsg("  [%s] Y ~ X + (1 | IID): %d subjects, %lld rows; "
                "sigma2=%.6g, D11=%.6g",
                data.phenoNames[k].c_str(), nSubj,
                static_cast<long long>(y.size()), fit.sigma2, fit.D(0, 0));
    }

    // First row of each subject = the (time-invariant) per-subject design.
    const int p = static_cast<int>(data.X.cols());
    out.perSubjectX.resize(nSubj, p);
    for (int i = 0; i < nSubj; ++i)
        out.perSubjectX.row(i) = data.X.row(data.subjStart[i]);
    out.xColNames.reserve(p);
    out.xColNames.emplace_back("intercept");
    for (const auto &cn : data.covarNames) out.xColNames.push_back(cn);

    return out;
}

Eigen::VectorXd alignToUsed(
    const std::vector<std::string> &srcIIDs,
    const Eigen::VectorXd &srcVals,
    const std::vector<std::string> &targetIIDs
) {
    const auto m = idMap(srcIIDs);
    Eigen::VectorXd out(static_cast<Eigen::Index>(targetIIDs.size()));
    for (size_t i = 0; i < targetIIDs.size(); ++i) {
        auto it = m.find(targetIIDs[i]);
        if (it == m.end())
            throw std::runtime_error(
                "--longitudinal: used subject '" + targetIIDs[i] +
                "' missing from the fit set");
        out[static_cast<Eigen::Index>(i)] = srcVals[it->second];
    }
    return out;
}

Eigen::MatrixXd alignRowsToUsed(
    const std::vector<std::string> &srcIIDs,
    const Eigen::MatrixXd &srcRows,
    const std::vector<std::string> &targetIIDs
) {
    const auto m = idMap(srcIIDs);
    Eigen::MatrixXd out(static_cast<Eigen::Index>(targetIIDs.size()), srcRows.cols());
    for (size_t i = 0; i < targetIIDs.size(); ++i) {
        auto it = m.find(targetIIDs[i]);
        if (it == m.end())
            throw std::runtime_error(
                "--longitudinal: used subject '" + targetIIDs[i] +
                "' missing from the fit set");
        out.row(static_cast<Eigen::Index>(i)) = srcRows.row(it->second);
    }
    return out;
}

} // namespace nsLongitudinal
