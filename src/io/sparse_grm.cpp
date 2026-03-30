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
