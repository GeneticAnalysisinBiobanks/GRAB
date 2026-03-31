// resid_file.cpp — Residual file parser
//
// Format: whitespace-delimited, no header.
//   subject  residual  [weight]  [indicator]
// '#'-prefixed lines and blank lines are silently skipped.

#include "io/resid_file.hpp"

#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

ResidData loadResidFile(const std::string& filename) {
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  if (!ifs) throw std::runtime_error("Cannot open resid file: " + filename);

  ResidData rd;
  std::vector<double> v_res, v_wt, v_ind;
  rd.subjects.reserve(65536);
  v_res.reserve(65536);
  v_wt.reserve(65536);
  v_ind.reserve(65536);

  uint32_t lineNo = 0;
  std::string line;
  while (std::getline(ifs, line)) {
    ++lineNo;
    // Normalise Windows CRLF
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;

    const char*       p   = line.c_str();
    const char* const end = p + line.size();

    // Skip leading whitespace
    while (p < end && (*p == ' ' || *p == '\t')) ++p;
    if (p == end) continue;  // blank / all-whitespace line

    // Column 1: subject ID (required)
    const char* tokStart = p;
    while (p < end && *p != ' ' && *p != '\t') ++p;
    std::string subj(tokStart, p);

    // Column 2: residual (required)
    while (p < end && (*p == ' ' || *p == '\t')) ++p;
    if (p == end)
      throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                               ": missing RESID column");
    char* next;
    double res = std::strtod(p, &next);
    if (next == p)
      throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                               ": invalid RESID value");
    p = next;

    // Column 3: weight (optional, default 1.0)
    double wt = 1.0;
    while (p < end && (*p == ' ' || *p == '\t')) ++p;
    if (p < end) {
      double tmp = std::strtod(p, &next);
      if (next != p) { wt = tmp; p = next; }
    }

    // Column 4: indicator (optional, default 0.0)
    double ind = 0.0;
    while (p < end && (*p == ' ' || *p == '\t')) ++p;
    if (p < end) {
      double tmp = std::strtod(p, &next);
      if (next != p) ind = tmp;
    }

    rd.subjects.push_back(std::move(subj));
    v_res.push_back(res);
    v_wt.push_back(wt);
    v_ind.push_back(ind);
  }

  const Eigen::Index n = static_cast<Eigen::Index>(rd.subjects.size());
  rd.residuals = Eigen::Map<Eigen::VectorXd>(v_res.data(), n);
  rd.weights   = Eigen::Map<Eigen::VectorXd>(v_wt.data(),  n);
  rd.indicator = Eigen::Map<Eigen::VectorXd>(v_ind.data(), n);
  return rd;
}


// ══════════════════════════════════════════════════════════════════════
// Multi-column residual matrix loader
// ══════════════════════════════════════════════════════════════════════

ResidMatrixData loadResidMatrix(const std::string& filename) {
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  if (!ifs) throw std::runtime_error("Cannot open resid matrix file: " + filename);

  ResidMatrixData rmd;
  std::vector<std::vector<double>> rows;
  rmd.subjects.reserve(65536);
  rows.reserve(65536);

  int nCols = -1;   // determined from first data row
  uint32_t lineNo = 0;
  std::string line;

  while (std::getline(ifs, line)) {
    ++lineNo;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty()) continue;
    if (line[0] == '#') continue;   // skip header / comment

    const char*       p   = line.c_str();
    const char* const end = p + line.size();

    // Skip leading whitespace
    while (p < end && (*p == ' ' || *p == '\t')) ++p;
    if (p == end) continue;

    // Column 1: subject ID
    const char* tokStart = p;
    while (p < end && *p != ' ' && *p != '\t') ++p;
    std::string subj(tokStart, p);

    // Remaining columns: numeric residual values
    std::vector<double> vals;
    vals.reserve(nCols > 0 ? static_cast<size_t>(nCols) : 16);
    while (p < end) {
      while (p < end && (*p == ' ' || *p == '\t')) ++p;
      if (p >= end) break;
      char* next;
      double v = std::strtod(p, &next);
      if (next == p)
        throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                                 ": invalid numeric value");
      vals.push_back(v);
      p = next;
    }

    if (vals.empty())
      throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                               ": no residual columns found");

    if (nCols < 0)
      nCols = static_cast<int>(vals.size());
    else if (static_cast<int>(vals.size()) != nCols)
      throw std::runtime_error(filename + " line " + std::to_string(lineNo) +
                               ": expected " + std::to_string(nCols) +
                               " columns, got " + std::to_string(vals.size()));

    rmd.subjects.push_back(std::move(subj));
    rows.push_back(std::move(vals));
  }

  if (rmd.subjects.empty())
    throw std::runtime_error(filename + ": no data rows found");

  const Eigen::Index N = static_cast<Eigen::Index>(rmd.subjects.size());
  const Eigen::Index K = static_cast<Eigen::Index>(nCols);
  rmd.residuals.resize(N, K);
  for (Eigen::Index i = 0; i < N; ++i)
    for (Eigen::Index j = 0; j < K; ++j)
      rmd.residuals(i, j) = rows[static_cast<size_t>(i)][static_cast<size_t>(j)];

  return rmd;
}

