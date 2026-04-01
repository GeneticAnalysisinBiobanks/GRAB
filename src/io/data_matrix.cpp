// data_matrix.cpp — Design file parser (with leading IID column), covariate projection, and generic matrix loader

#include "io/data_matrix.hpp"

#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

DesignMatrix::DesignMatrix(const Eigen::MatrixXd& X)
  : m_X(X)
{
  const int nCols = static_cast<int>(m_X.cols());
  // Pre-compute tX and XinvXX = X (X'X)^{-1}
  m_tX = m_X.transpose();                              // p × N
  Eigen::MatrixXd XtX = m_tX * m_X;                    // p × p
  m_XinvXX = m_X * XtX.ldlt().solve(
      Eigen::MatrixXd::Identity(nCols, nCols));         // N × p
}

void DesignMatrix::adjustGenotype(
    const double* G, const uint32_t* nzIdx, int nNz,
    Eigen::Ref<Eigen::VectorXd> adjG) const {

  const int p = nCols();
  const int N = nRows();

  // tX_g = tX * G  (p × 1), accumulated only over non-zero entries
  Eigen::VectorXd tX_g = Eigen::VectorXd::Zero(p);
  for (int k = 0; k < nNz; ++k) {
    uint32_t j = nzIdx[k];
    double gj = G[j];
    // m_tX is p × N  →  col(j) is the j-th subject's covariate column
    tX_g.noalias() += gj * m_tX.col(j);
  }

  // adjG = G - XinvXX * tX_g
  Eigen::Map<const Eigen::VectorXd> gVec(G, N);
  adjG.noalias() = gVec - m_XinvXX * tX_g;
}

// ======================================================================
// Lightweight numeric matrix loader
// ======================================================================

Eigen::MatrixXd loadNumericMatrix(const std::string& filename) {
  std::ifstream ifs(filename, std::ios::in | std::ios::binary);
  if (!ifs)
    throw std::runtime_error("Cannot open numeric matrix file: " + filename);

  std::vector<double> vals;
  vals.reserve(65536);
  int nCols = 0;
  int nRows = 0;

  std::string line;
  while (std::getline(ifs, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.empty() || line[0] == '#') continue;

    const char*       p   = line.c_str();
    const char* const end = p + line.size();
    int colsThisRow = 0;

    while (p < end) {
      while (p < end && (*p == ' ' || *p == '\t')) ++p;
      if (p >= end) break;
      char* next;
      double v = std::strtod(p, &next);
      if (next == p)
        throw std::runtime_error(
            "loadNumericMatrix: non-numeric token on row " + std::to_string(nRows + 1) +
            ": " + filename);
      vals.push_back(v);
      ++colsThisRow;
      p = next;
    }
    if (colsThisRow == 0) continue;  // blank / all-whitespace line

    if (nRows == 0) {
      nCols = colsThisRow;
    } else if (colsThisRow != nCols) {
      throw std::runtime_error(
          "loadNumericMatrix: row " + std::to_string(nRows + 1) +
          " has " + std::to_string(colsThisRow) +
          " columns (expected " + std::to_string(nCols) + ")");
    }
    ++nRows;
  }

  if (nRows == 0 || nCols == 0)
    throw std::runtime_error("loadNumericMatrix: empty file: " + filename);

  Eigen::MatrixXd M(nRows, nCols);
  for (int r = 0; r < nRows; ++r)
    for (int c = 0; c < nCols; ++c)
      M(r, c) = vals[r * nCols + c];

  return M;
}


