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
