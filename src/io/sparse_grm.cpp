// sparse_grm.cpp — SparseGRM parser and quadratic-form implementation

#include "io/sparse_grm.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

SparseGRM::SparseGRM(const std::string& filename,
                     const std::vector<std::string>& subjectOrder) {
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

  while (std::getline(ifs, line)) {
    if (line.empty() || line[0] == '#') continue;

    // Fast tab-delimited parse: ID1 \t ID2 \t value
    const char* p = line.c_str();
    const char* t1 = p;
    while (*t1 && *t1 != '\t') ++t1;
    if (!*t1) continue;
    std::string id1(p, t1);

    const char* t2 = t1 + 1;
    while (*t2 && *t2 != '\t') ++t2;
    if (!*t2) continue;
    std::string id2(t1 + 1, t2);

    auto it1 = idMap.find(id1);
    auto it2 = idMap.find(id2);
    if (it1 == idMap.end() || it2 == idMap.end()) continue;

    double val;
    char* end;
    val = std::strtod(t2 + 1, &end);
    if (end == t2 + 1) continue;

    uint32_t r = it1->second;
    uint32_t c = it2->second;
    raw.push_back({r, c, val});
    if (r != c)
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
