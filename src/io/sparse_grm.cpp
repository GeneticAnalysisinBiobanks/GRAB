// sparse_grm.cpp — SparseGRM parser and quadratic-form implementation

#include "io/sparse_grm.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

SparseGRM::SparseGRM(const std::string& filename,
                     const std::vector<std::string>& subjectOrder,
                     bool symmetrize) {
  // Build subject → index map
  std::unordered_map<std::string, uint32_t> idMap;
  idMap.reserve(subjectOrder.size());
  for (uint32_t i = 0; i < subjectOrder.size(); ++i)
    idMap.emplace(subjectOrder[i], i);
  m_nSubj = static_cast<uint32_t>(subjectOrder.size());

  // Parse file
  std::ifstream ifs(filename);
  if (!ifs)
    throw std::runtime_error("SparseGRM: cannot open " + filename);

  std::string line;
  std::vector<Entry> raw;
  raw.reserve(subjectOrder.size() * 4);  // rough guess

  uint32_t lineNo = 0;
  while (std::getline(ifs, line)) {
    ++lineNo;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;

    // Whitespace-delimited parse: ID1  ID2  VALUE
    const char*       p   = line.c_str();
    const char* const end = p + line.size();
    auto skipWS  = [&]() { while (p < end && (*p == ' ' || *p == '\t')) ++p; };
    auto nextTok = [&]() -> std::string {
      skipWS();
      const char* s = p;
      while (p < end && *p != ' ' && *p != '\t') ++p;
      return std::string(s, p);
    };

    skipWS();
    if (p >= end) continue;  // blank / all-whitespace line

    const std::string id1 = nextTok();
    const std::string id2 = nextTok();
    skipWS();
    if (id2.empty() || p >= end)
      throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                               ": expected ID1 ID2 VALUE");
    char* endPtr;
    const double val = std::strtod(p, &endPtr);
    if (endPtr == p)
      throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                               ": invalid VALUE field");

    auto it1 = idMap.find(id1);
    auto it2 = idMap.find(id2);
    if (it1 == idMap.end() || it2 == idMap.end()) continue;

    const uint32_t r = it1->second;
    const uint32_t c = it2->second;
    raw.push_back({r, c, val});
    if (symmetrize && r != c)
      raw.push_back({c, r, val});  // symmetrise
  }

  // Sort by (row, col) for cache-friendly access in quadForm
  std::sort(raw.begin(), raw.end(), [](const Entry& a, const Entry& b) {
    return a.row < b.row || (a.row == b.row && a.col < b.col);
  });

  m_entries = std::move(raw);
}

// ══════════════════════════════════════════════════════════════════════
// readGctaIIDs — shared .grm.id parser
// ══════════════════════════════════════════════════════════════════════

std::vector<std::string> SparseGRM::readGctaIIDs(const std::string& prefix) {
  const std::string idFile = prefix + ".grm.id";
  std::ifstream idStream(idFile);
  if (!idStream)
    throw std::runtime_error("SparseGRM::readGctaIIDs: cannot open " + idFile);

  std::vector<std::string> iids;
  iids.reserve(65536);
  std::string line;
  while (std::getline(idStream, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    // .grm.id format: FID  IID  (whitespace-separated)
    const char*       p   = line.c_str();
    const char* const end = p + line.size();
    auto skipWS  = [&]() { while (p < end && (*p == ' ' || *p == '\t')) ++p; };
    auto nextTok = [&]() -> std::string {
      skipWS();
      const char* s = p;
      while (p < end && *p != ' ' && *p != '\t') ++p;
      return std::string(s, p);
    };
    /*FID*/ nextTok();
    std::string iid = nextTok();
    if (iid.empty())
      throw std::runtime_error(idFile + ": missing IID on line " +
                               std::to_string(iids.size() + 1));
    iids.push_back(std::move(iid));
  }
  return iids;
}

// ══════════════════════════════════════════════════════════════════════
// GCTA format: PREFIX.grm.id + PREFIX.grm.sp
// ══════════════════════════════════════════════════════════════════════

SparseGRM SparseGRM::fromGCTA(const std::string& prefix,
                              const std::vector<std::string>& subjectOrder,
                              bool symmetrize) {
  SparseGRM grm;

  // Build subject → canonical index map
  std::unordered_map<std::string, uint32_t> idMap;
  idMap.reserve(subjectOrder.size());
  for (uint32_t i = 0; i < static_cast<uint32_t>(subjectOrder.size()); ++i)
    idMap.emplace(subjectOrder[i], i);
  grm.m_nSubj = static_cast<uint32_t>(subjectOrder.size());

  // ── Read .grm.id → 0-based file-index → IID ──────────────────────
  auto fileIIDs = readGctaIIDs(prefix);

  // Map file-index → canonical index (UINT32_MAX if not in subjectOrder)
  std::vector<uint32_t> fileToCanon(fileIIDs.size(), UINT32_MAX);
  for (uint32_t fi = 0; fi < static_cast<uint32_t>(fileIIDs.size()); ++fi) {
    auto it = idMap.find(fileIIDs[fi]);
    if (it != idMap.end())
      fileToCanon[fi] = it->second;
  }

  // ── Read .grm.sp ──────────────────────────────────────────────────
  const std::string spFile = prefix + ".grm.sp";
  std::ifstream spStream(spFile);
  if (!spStream)
    throw std::runtime_error("SparseGRM::fromGCTA: cannot open " + spFile);

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
    const char* p = line.c_str();
    char* endPtr;
    unsigned long idx1 = std::strtoul(p, &endPtr, 10);
    if (endPtr == p)
      throw std::runtime_error(spFile + " line " + std::to_string(lineNo) +
                               ": invalid first index");
    p = endPtr;
    unsigned long idx2 = std::strtoul(p, &endPtr, 10);
    if (endPtr == p)
      throw std::runtime_error(spFile + " line " + std::to_string(lineNo) +
                               ": invalid second index");
    p = endPtr;
    double val = std::strtod(p, &endPtr);
    if (endPtr == p)
      throw std::runtime_error(spFile + " line " + std::to_string(lineNo) +
                               ": invalid value");

    if (idx1 >= nFileIDs || idx2 >= nFileIDs) continue;
    uint32_t r = fileToCanon[idx1];
    uint32_t c = fileToCanon[idx2];
    if (r == UINT32_MAX || c == UINT32_MAX) continue;

    raw.push_back({r, c, val});
    if (symmetrize && r != c)
      raw.push_back({c, r, val});
  }

  std::sort(raw.begin(), raw.end(), [](const Entry& a, const Entry& b) {
    return a.row < b.row || (a.row == b.row && a.col < b.col);
  });

  grm.m_entries = std::move(raw);
  return grm;
}

double SparseGRM::quadForm(const double* x, uint32_t n) const {
  if (n != m_nSubj)
    throw std::runtime_error("SparseGRM::quadForm: size mismatch");
  double sum = 0.0;
  for (const auto& e : m_entries)
    sum += e.value * x[e.row] * x[e.col];
  return sum;
}

double SparseGRM::bilinearForm(const double* x, const double* y, uint32_t n) const {
  if (n != m_nSubj)
    throw std::runtime_error("SparseGRM::bilinearForm: size mismatch");
  double sum = 0.0;
  for (const auto& e : m_entries)
    sum += e.value * x[e.row] * y[e.col];
  return sum;
}

double SparseGRM::spaVariance(const double* R, uint32_t n) const {
  if (n != m_nSubj)
    throw std::runtime_error("SparseGRM::spaVariance: size mismatch");
  double covSum = 0.0;
  for (const auto& e : m_entries)
    covSum += e.value * R[e.row] * R[e.col];
  // dot(R, R)
  double dotRR = 0.0;
  for (uint32_t i = 0; i < n; ++i)
    dotRR += R[i] * R[i];
  return 2.0 * covSum - dotRR;
}
