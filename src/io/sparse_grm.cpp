// sparse_grm.cpp — SparseGRM parser and quadratic-form implementation

#include "io/sparse_grm.hpp"
#include "util/logging.hpp"
#include "util/text_scanner.hpp"

#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <unordered_set>

namespace {
// Derive .grm.id path from a .grm.sp file path.
std::string grmIdPath(const std::string &spFile) {
    const std::string suffix = ".grm.sp";
    if (spFile.size() > suffix.size() && spFile.compare(spFile.size() - suffix.size(), suffix.size(), suffix) == 0)
        return spFile.substr(0, spFile.size() - suffix.size()) + ".grm.id";
    // Fallback: strip last extension, add .grm.id
    auto dot = spFile.rfind('.');
    return (dot != std::string::npos ? spFile.substr(0, dot) : spFile) + ".grm.id";
}

} // anonymous namespace

SparseGRM::SparseGRM(
    const std::string &filename,
    const std::vector<std::string> &subjectOrder
)
{
    auto idMap = text::buildIIDMap(subjectOrder);
    m_nSubj = static_cast<uint32_t>(subjectOrder.size());

    // Parse file
    std::ifstream ifs(filename);
    if (!ifs) throw std::runtime_error("SparseGRM: cannot open " + filename);

    std::string line;
    std::vector<Entry> raw;
    raw.reserve(subjectOrder.size() * 4); // rough guess

    // Read header
    if (!std::getline(ifs, line))
        throw std::runtime_error(filename + ": empty GRM file");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.rfind("ID1", 0) != 0 && line.rfind("#ID1", 0) != 0 &&
        line.rfind("IID1", 0) != 0 && line.rfind("#IID1", 0) != 0)
        throw std::runtime_error(filename + ": header must start with 'ID1', '#ID1', 'IID1', or '#IID1'");

    uint32_t lineNo = 1;
    while (std::getline(ifs, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        text::TokenScanner tok(line);
        tok.skipWS();
        if (tok.atEnd()) continue;

        const std::string id1 = tok.next();
        const std::string id2 = tok.next();
        tok.skipWS();
        if (id2.empty() || tok.atEnd())
            throw std::runtime_error(filename + " line " + std::to_string(lineNo) + ": expected ID1 ID2 VALUE");
        char *endPtr;
        const double val = std::strtod(tok.pos(), &endPtr);
        if (endPtr == tok.pos())
            throw std::runtime_error(filename + " line " + std::to_string(lineNo) + ": invalid VALUE field");

        auto it1 = idMap.find(id1);
        auto it2 = idMap.find(id2);
        if (it1 == idMap.end() || it2 == idMap.end()) continue;

        const uint32_t r = it1->second;
        const uint32_t c = it2->second;
        raw.push_back({r, c, val});
    }

    // Sort by (row, col) for cache-friendly access in quadForm
    std::sort(raw.begin(), raw.end(),
              [](const Entry &a, const Entry &b) {
        return a.row < b.row || (a.row == b.row && a.col < b.col);
    });

    m_entries = std::move(raw);
    buildDiagonal();
}

// ══════════════════════════════════════════════════════════════════════
// readGctaIIDs — shared .grm.id parser
// ══════════════════════════════════════════════════════════════════════

std::vector<std::string> SparseGRM::readGctaIIDs(const std::string &spFile) {
    const std::string idFile = grmIdPath(spFile);
    std::ifstream idStream(idFile);
    if (!idStream) return {}; // file not found — caller handles fallback

    std::vector<std::string> iids;
    iids.reserve(65536);
    std::string line;
    while (std::getline(idStream, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        text::TokenScanner tok(line);
        /*FID*/ tok.next();
        std::string iid = tok.next();
        if (iid.empty()) throw std::runtime_error(idFile + ": missing IID on line " + std::to_string(iids.size() + 1));
        iids.push_back(std::move(iid));
    }
    return iids;
}

// ══════════════════════════════════════════════════════════════════════
// GCTA format: .grm.sp file + companion .grm.id
// ══════════════════════════════════════════════════════════════════════

SparseGRM SparseGRM::fromGCTA(
    const std::string &spFile,
    const std::vector<std::string> &subjectOrder,
    const std::vector<std::string> &famIIDs
) {
    SparseGRM grm;

    auto idMap = text::buildIIDMap(subjectOrder);
    grm.m_nSubj = static_cast<uint32_t>(subjectOrder.size());

    // ── Read .grm.id → 0-based file-index → IID ──────────────────────
    auto fileIIDs = readGctaIIDs(spFile);
    if (fileIIDs.empty()) {
        // .grm.id not found — fall back to .fam IID order
        if (famIIDs.empty())
            throw std::runtime_error("SparseGRM::fromGCTA: " + grmIdPath(spFile) +
                                     " not found and no .fam fallback provided");
        infoMsg("--sp-grm-plink2: %s not found; .grm.id assumed identical to .fam order", grmIdPath(spFile).c_str());
        fileIIDs = famIIDs;
    } else {
        // .grm.id found — validate overlap with subjectOrder (.fam)
        uint32_t nMatched = 0;
        for (const auto &iid : fileIIDs)
            if (idMap.count(iid)) ++nMatched;

        if (nMatched == 0 && !fileIIDs.empty())
            throw std::runtime_error("SparseGRM::fromGCTA: .grm.id (" + grmIdPath(spFile) + ") has " +
                                     std::to_string(fileIIDs.size()) + " subjects but none match .fam (" +
                                     std::to_string(subjectOrder.size()) + " subjects)");

        if (nMatched == subjectOrder.size() && fileIIDs.size() == subjectOrder.size()) {
            // Check if order is identical
            bool identical = true;
            for (size_t i = 0; i < fileIIDs.size(); ++i)
                if (fileIIDs[i] != subjectOrder[i]) {
                    identical = false;
                    break;
                }
            if (identical)
                infoMsg("--sp-grm-plink2: .grm.id identical to .fam (%u subjects)", nMatched);
            else
                infoMsg("--sp-grm-plink2: .grm.id has %u/%zu subjects matching .fam (reordered)", nMatched,
                        fileIIDs.size());
        } else {
            infoMsg("--sp-grm-plink2: .grm.id has %zu subjects, %u match .fam (%zu subjects)", fileIIDs.size(),
                    nMatched, subjectOrder.size());
        }
    }

    // Map file-index → canonical index (UINT32_MAX if not in subjectOrder)
    std::vector<uint32_t> fileToCanon(fileIIDs.size(), UINT32_MAX);
    for (uint32_t fi = 0; fi < static_cast<uint32_t>(fileIIDs.size()); ++fi) {
        auto it = idMap.find(fileIIDs[fi]);
        if (it != idMap.end()) fileToCanon[fi] = it->second;
    }

    // ── Read .grm.sp ──────────────────────────────────────────────────
    std::ifstream spStream(spFile);
    if (!spStream) throw std::runtime_error("SparseGRM::fromGCTA: cannot open " + spFile);

    std::vector<Entry> raw;
    raw.reserve(subjectOrder.size() * 4);

    const uint32_t nFileIDs = static_cast<uint32_t>(fileIIDs.size());
    std::string line;
    uint32_t lineNo = 0;
    while (std::getline(spStream, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        // Format: idx1  idx2  value  (whitespace-separated, 0-based)
        const char *p = line.c_str();
        char *endPtr;
        unsigned long idx1 = std::strtoul(p, &endPtr, 10);
        if (endPtr == p) throw std::runtime_error(spFile + " line " + std::to_string(lineNo) + ": invalid first index");
        p = endPtr;
        unsigned long idx2 = std::strtoul(p, &endPtr, 10);
        if (endPtr == p)
            throw std::runtime_error(spFile + " line " + std::to_string(lineNo) + ": invalid second index");
        p = endPtr;
        double val = std::strtod(p, &endPtr);
        if (endPtr == p) throw std::runtime_error(spFile + " line " + std::to_string(lineNo) + ": invalid value");

        if (idx1 >= nFileIDs || idx2 >= nFileIDs) continue;
        uint32_t r = fileToCanon[idx1];
        uint32_t c = fileToCanon[idx2];
        if (r == UINT32_MAX || c == UINT32_MAX) continue;

        raw.push_back({r, c, val});
    }

    std::sort(raw.begin(), raw.end(),
              [](const Entry &a, const Entry &b) {
        return a.row < b.row || (a.row == b.row && a.col < b.col);
    });

    grm.m_entries = std::move(raw);
    grm.buildDiagonal();
    return grm;
}

void SparseGRM::buildDiagonal() {
    m_diagonal.assign(m_nSubj, 0.0);
    for (const auto &e : m_entries) {
        if (e.row == e.col && e.row < m_nSubj) m_diagonal[e.row] = e.value;
    }
}

void SparseGRM::multiply(
    const double *x,
    double *result,
    uint32_t n
) const {
    if (n != m_nSubj) throw std::runtime_error("SparseGRM::multiply: size mismatch");
    std::fill(result, result + n, 0.0);
    for (const auto &e : m_entries) {
        result[e.row] += e.value * x[e.col];
        if (e.row != e.col) result[e.col] += e.value * x[e.row];
    }
}

double SparseGRM::quadForm(
    const double *x,
    uint32_t n
) const {
    if (n != m_nSubj) throw std::runtime_error("SparseGRM::quadForm: size mismatch");
    double sum = 0.0;
    for (const auto &e : m_entries) {
        double prod = e.value * x[e.row] * x[e.col];
        sum += prod;
        if (e.row != e.col) sum += prod; // off-diagonal: count both (i,j) and (j,i)
    }
    return sum;
}

double SparseGRM::bilinearForm(
    const double *x,
    const double *y,
    uint32_t n
) const {
    if (n != m_nSubj) throw std::runtime_error("SparseGRM::bilinearForm: size mismatch");
    double sum = 0.0;
    for (const auto &e : m_entries) {
        if (e.row == e.col)
            sum += e.value * x[e.row] * y[e.col];
        else
            sum += e.value * (x[e.row] * y[e.col] + x[e.col] * y[e.row]);
    }
    return sum;
}

double SparseGRM::spaVariance(
    const double *R,
    uint32_t n
) const {
    if (n != m_nSubj) throw std::runtime_error("SparseGRM::spaVariance: size mismatch");
    double covSum = 0.0;
    for (const auto &e : m_entries)
        covSum += e.value * R[e.row] * R[e.col];
    // dot(R, R)
    double dotRR = 0.0;
    for (uint32_t i = 0; i < n; ++i)
        dotRR += R[i] * R[i];
    return 2.0 * covSum - dotRR;
}

SparseGRM SparseGRM::load(
    const std::string &grabFile,
    const std::string &gctaFile,
    const std::vector<std::string> &subjectOrder,
    const std::vector<std::string> &famIIDs
) {
    if (!gctaFile.empty()) return fromGCTA(gctaFile, subjectOrder, famIIDs);
    return SparseGRM(grabFile, subjectOrder);
}

SparseGRM SparseGRM::fromEntries(
    uint32_t nSubj,
    std::vector<Entry> entries
) {
    SparseGRM g;
    g.m_nSubj = nSubj;
    g.m_entries = std::move(entries);
    g.buildDiagonal();
    return g;
}

// ══════════════════════════════════════════════════════════════════════
// Lightweight subject-ID scanner (no numeric data loaded)
// ══════════════════════════════════════════════════════════════════════

std::unordered_set<std::string> SparseGRM::parseSubjectIDs(
    const std::string &grabFile,
    const std::string &gctaFile,
    const std::vector<std::string> &famIIDs
) {
    std::unordered_set<std::string> ids;

    if (!gctaFile.empty()) {
        // ── GCTA format: read .grm.id ──────────────────────────────────
        auto fileIIDs = readGctaIIDs(gctaFile);
        if (fileIIDs.empty()) {
            // .grm.id not found — all genotype subjects assumed in GRM
            if (!famIIDs.empty()) {
                infoMsg("--sp-grm-plink2: .grm.id not found; assuming all .fam subjects are in GRM");
                ids.insert(famIIDs.begin(), famIIDs.end());
            }
        } else {
            ids.insert(fileIIDs.begin(), fileIIDs.end());
        }
        return ids;
    }

    // ── GRAB format: scan ID1/ID2 columns ─────────────────────────────
    std::ifstream ifs(grabFile);
    if (!ifs) throw std::runtime_error("SparseGRM::parseSubjectIDs: cannot open " + grabFile);

    std::string line;
    if (!std::getline(ifs, line))
        throw std::runtime_error(grabFile + ": empty GRM file");
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.rfind("ID1", 0) != 0 && line.rfind("#ID1", 0) != 0 &&
        line.rfind("IID1", 0) != 0 && line.rfind("#IID1", 0) != 0)
        throw std::runtime_error(grabFile + ": header must start with 'ID1', '#ID1', 'IID1', or '#IID1'");

    while (std::getline(ifs, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        text::TokenScanner tok(line);
        tok.skipWS();
        if (tok.atEnd()) continue;
        std::string id1 = tok.next();
        std::string id2 = tok.next();
        if (!id1.empty()) ids.insert(std::move(id1));
        if (!id2.empty()) ids.insert(std::move(id2));
    }
    return ids;
}
